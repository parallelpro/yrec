!----------------------------------------------------------------------
! transport_checks
!----------------------------------------------------------------------
! The checks secular_transport makes after each iteration of its
! transport solve: on the new specific angular momentum, from which
! they also derive the new rotation curve (check_angular_momentum), and
! on the new abundances, from which they also derive the new mean
! molecular weight (check_composition).
module transport_checks
      use rotation_configuration, only: rotation_settings, last_mu_species
      use rotation_messages, only: rotation_report, rotation_note, message_length, &
           unwritten_line
      implicit none
      private
      public :: check_angular_momentum, check_composition, check_mixed_species

! check_angular_momentum's iteration history holds
! max_recorded_iterations (io's derive_controls rejects a rotating run
! with more). A negative interior abundance smaller than
! small_abundance_fraction of the surface value is repaired; a larger
! one cuts the sub-step, and the cut after max_substep_cuts stops the
! run.
      integer, parameter :: max_recorded_iterations = 16
      double precision, parameter :: small_abundance_fraction = 1.0d-5
      integer, parameter :: max_substep_cuts = 3
contains

! ---------------------------------------------------------------
! Checks the specific angular momentum j after an iteration of the
! transport solve and finds the rotation curve that goes with it (the
! legacy CHECKJ). In order:
!  1. reject_negative_j: a zone with j <= 0 stops the run (ierr = 1).
!  2. measure_convergence: the largest fractional change of j since the
!     previous iteration; the iteration has converged when its size is
!     at most rotation_settings%convergence_tolerance.
!  3. update_omega: omega and the moments of inertia for the new j
!     (omega_from_j).
!  4. report_largest_change: on the last iteration, or at convergence,
!     prints the largest fractional change of j over the sub-step and
!     the history of measure_convergence.
!
! Inputs:
!   log_density, log_radius, log_mass : log10 of rho, r and m at the
!      zones.
!   shell_mass : the mass of each zone (g).
!   specific_angular_momentum_start : j at the start of the sub-step.
!   iteration_number : this iteration, 1 to
!      rotation_settings%max_diffusion_iters.
!   am_transport_convective_flag : true where a zone is convective for
!      angular momentum transport (overshoot regions included).
!   num_zones : the number of zones.
! In/out:
!   specific_angular_momentum : j after this iteration (cm^2/s);
!      update_omega resets it to I omega/m for the omega it finds.
!   specific_angular_momentum_prev : j after the previous iteration;
!      set to j once compared with it.
!   cut_count : the number of sub-step cuts so far.
! Outputs:
!   omega : the angular velocity (rad/s).
!   moment_of_inertia, qiw : each zone's moment of inertia and its
!      derivative with respect to omega.
!   eta_squared, mean_radius : the shape quantities (eta_2 and r0) that
!      omega_from_j finds along with omega.
!   converged_flag : true if the iteration has converged.
!   redo_flag : false: a negative j stops the run rather than asking
!      for a redo (see reject_negative_j).
!   ierr : nonzero stops the run.
subroutine check_angular_momentum(s, log_density, specific_angular_momentum_prev, &
     specific_angular_momentum_start, log_radius, log_mass, shell_mass, &
     iteration_number, &
     am_transport_convective_flag, num_zones, eta_squared, &
     moment_of_inertia, specific_angular_momentum, cut_count, &
     converged_flag, redo_flag, omega, &
     qiw, mean_radius, ierr)
      use omega_profile, only: solid_body_omega, omega_from_j, &
           default_inertia_tolerance
      use star_info_lib, only: star_info, json
      use io_units
      implicit none
      type(star_info), intent(in) :: s

      double precision, intent(in) :: log_density(json)
      double precision, intent(inout) :: specific_angular_momentum_prev(json)
      double precision, intent(in) :: specific_angular_momentum_start(json)
      double precision, intent(in) :: log_radius(json), log_mass(json), &
           shell_mass(json)
      integer, intent(in) :: iteration_number
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      integer, intent(inout) :: cut_count
      logical, intent(inout) :: converged_flag
      logical, intent(out) :: redo_flag
      double precision, intent(inout) :: omega(json)
      double precision, intent(inout) :: qiw(json), mean_radius(json)
      integer, intent(out) :: ierr
! locals
! The history of measure_convergence, one entry per iteration: room for
! 16 iterations, and max_diffusion_iters is not checked against it.
      double precision :: max_delta_j_by_iter(max_recorded_iterations)
      integer :: max_delta_j_zone_by_iter(max_recorded_iterations)
      integer :: zone_index, scan_index
      double precision :: delta_j_fraction
      double precision :: inertia_tolerance, shape_tolerance
      double precision :: max_fractional_dj
      integer :: max_dj_zone
! a message's lines (rotation_messages)
      character(len=message_length) :: lines(8)

      ierr = 0

      converged_flag = .false.
      redo_flag = .false.
      call reject_negative_j()
      if (ierr /= 0 .or. redo_flag) return
      call measure_convergence()
      call update_omega()
      if(iteration_number.eq.rotation_settings%max_diffusion_iters.or.converged_flag) &
           call report_largest_change()

contains

! ---------------------------------------------------------------
! The first zone with j <= 0 stops the run (ierr = 1), as the F77
! CHECKJ in effect did: its cut-and-redo branch was unreachable
! (cut_count > 0 always held after the increment), so the stop
! message's "3 attempts" were never made.
subroutine reject_negative_j()
      do zone_index = 1,num_zones
         if(specific_angular_momentum(zone_index).le.0.0d0) then
            cut_count = cut_count + 1
            lines = unwritten_line
            write(lines,1000) zone_index
            call rotation_report(lines)
 1000       format(1x,39('>'),39('<')/5x,'ERROR IN SR CHECKJ'/ &
                    5x,'NEGATIVE J/M ENCOUNTERED IN ZONE',i5, &
                    ' AND 3 ATTEMPTS AT CUTTING TIMESTEP FAILED'/ &
                    'RUN STOPPED')
            ierr = 1
            return
         endif
      end do
end subroutine reject_negative_j

! ---------------------------------------------------------------
! The largest fractional change of j since the previous iteration,
!   max_delta_j_by_iter(it) = (j - j_prev)/j_prev at the zone where
!   |j - j_prev|/j_prev is largest (kept with its sign),
! and that zone, max_delta_j_zone_by_iter(it), for the printed history.
! The iteration has converged when |max_delta_j_by_iter(it)| <=
! convergence_tolerance. j_prev is set to j as the scan passes.
subroutine measure_convergence()
      max_delta_j_by_iter(iteration_number) = &
           (specific_angular_momentum(1)-specific_angular_momentum_prev(1))/ &
           specific_angular_momentum_prev(1)
      specific_angular_momentum_prev(1) = specific_angular_momentum(1)
      max_delta_j_zone_by_iter(iteration_number) = 1
      do zone_index = 2,num_zones
         delta_j_fraction = &
              (specific_angular_momentum(zone_index)- &
              specific_angular_momentum_prev(zone_index))/ &
              specific_angular_momentum_prev(zone_index)
         if(abs(delta_j_fraction).gt. &
              abs(max_delta_j_by_iter(iteration_number))) then
            max_delta_j_by_iter(iteration_number) = delta_j_fraction
            max_delta_j_zone_by_iter(iteration_number) = zone_index
         endif
         specific_angular_momentum_prev(zone_index) = &
              specific_angular_momentum(zone_index)
      end do
      if(abs(max_delta_j_by_iter(iteration_number)).le. &
           rotation_settings%convergence_tolerance) then
         converged_flag = .true.
      endif
end subroutine measure_convergence

! ---------------------------------------------------------------
! omega for the new j (omega_from_j: solid bodies in the convection
! zones, one zone at a time elsewhere). The moments of inertia depend
! on omega, so omega_from_j iterates. j itself is only good to about
! convergence_tolerance until the transport has converged, so before
! the last iteration and convergence the tolerances of that iteration
! (inertia_tolerance, otherwise zero) and of the shape (shape_tolerance,
! otherwise rotation_settings%shape_tolerance) are loosened to at least
! convergence_tolerance/100 for this call.
subroutine update_omega()
      inertia_tolerance = default_inertia_tolerance
      shape_tolerance = rotation_settings%shape_tolerance
      if(iteration_number.lt.rotation_settings%max_diffusion_iters.and..not.converged_flag)then
         inertia_tolerance = &
              max(rotation_settings%convergence_tolerance*1.0d-2,default_inertia_tolerance)
         shape_tolerance = max(rotation_settings%convergence_tolerance*1.0d-2, &
              rotation_settings%shape_tolerance)
      endif
      call omega_from_j(log_density,specific_angular_momentum,log_radius,log_mass, &
           shell_mass,am_transport_convective_flag,num_zones,eta_squared, &
           moment_of_inertia,omega,qiw,mean_radius, &
           s%disk_gate_age_gyr, &
           inertia_tolerance=inertia_tolerance, shape_tolerance=shape_tolerance)
end subroutine update_omega

! (The safeguard against reversals of the omega gradient, the F77
! CHECKJ's solid-body smoothing where omega had come to rise outward,
! is gone, 2026-10-28: its entry test, omega(i) - omega(i-1) > 1 rad/s,
! never held for a star, so it never acted, and at the 1e-15 rad/s of
! its other tests it would act on the map back's interpolation noise.
! docs/bug_sweep.md O15; git history keeps it.)


! ---------------------------------------------------------------
! Prints the largest fractional change of j over the sub-step,
! (j - j_start)/j_start, with its zone, and after it the history of
! measure_convergence (each iteration's largest change and its zone);
! nothing is printed when the solve had converged before this check.
subroutine report_largest_change()
      max_fractional_dj = 0.0d0
      max_dj_zone = 0
      do zone_index = 1,num_zones
         delta_j_fraction = &
              (specific_angular_momentum(zone_index)- &
              specific_angular_momentum_start(zone_index))/ &
              specific_angular_momentum_start(zone_index)
         if(abs(delta_j_fraction).gt.abs(max_fractional_dj)) then
            max_fractional_dj = delta_j_fraction
            max_dj_zone = zone_index
         endif
      end do
      lines = unwritten_line
      write(lines,160)max_fractional_dj,max_dj_zone, &
           (max_delta_j_by_iter(scan_index), &
           max_delta_j_zone_by_iter(scan_index),scan_index=1, &
           iteration_number)
      call rotation_note(lines)
  160       format(' MAX D(J/M)/(J/M)',1pe12.3,' AT PT.',i5, &
              ' BY ITERATION'/5(1x,e11.3,i4))
end subroutine report_largest_change

end subroutine check_angular_momentum

! ---------------------------------------------------------------
! Checks the abundances after an iteration of the rotational mixing and
! finds the mean molecular weight they imply (the legacy CHECKC).
!
! A mass fraction outside [0, 1] in an interior zone (2 to num_zones-1)
! is repaired if it is a small negative value, |X| < 1e-5 X(surface),
! and otherwise halves the sub-step and asks for a redo; the fourth cut
! of the sub-step (cut_count counts both checks' cuts) stops the run.
! Before the last iteration species 1-4 are checked (H1, He4 and He3;
! slot 3, Z, is skipped), on the last species 1-11.
! secular_transport mixes species 5 and up only after its iterations,
! so on the last iteration those are checked before their mixing.
!
! From the second iteration on, the mean molecular weight of a fully
! ionized gas, with A the atomic weights and the metals taken as C12
! for the ions and as 2 nucleons per electron:
!   1/mu_ion = amuenv + dX/A(H1) + dY/A(He4) + dZ/12 + dX3/A(He3),
!   1/mu_e   = X/A(H1) + 2 (X3/A(He3) + Y/A(He4)) + Z/2,
!   mu = 1/(1/mu_ion + 1/mu_e) = mu_ion mu_e/(mu_ion + mu_e),
! where dX, dY, dZ, dX3 are the differences from the envelope mixture
! and star%envelope_inverse_atomic_weight = sum X_i/A_i of that mixture.
!
! Inputs:
!   iteration_number : this iteration, 1 to
!      rotation_settings%max_diffusion_iters.
!   num_zones : the number of zones.
! In/out:
!   composition(15, json) : the mass fractions; small negative values
!      are repaired in place.
!   dt : the sub-step (s), halved on a cut.
!   cut_count : the number of sub-step cuts so far.
!   converged_flag : set false on a cut.
! Outputs:
!   redo_flag : true if the sub-step must be redone with the halved dt.
!   sub_physics%mu : the new mean molecular weight (from the second
!      iteration on).
!   ierr : nonzero stops the run.
subroutine check_composition(s, sub_physics, composition, iteration_number, num_zones, &
     dt, cut_count, converged_flag, redo_flag, ierr)
      use rotation_scratch_lib, only: substep_physics_state
      use star_info_lib, only: star_info, json, i_h1, i_he4, i_metals, i_he3, n_species_basic
      use io_units
      implicit none
      type(star_info), intent(in) :: s
      type(substep_physics_state), intent(inout) :: sub_physics

      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: iteration_number
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: dt
      integer, intent(inout) :: cut_count
      logical, intent(inout) :: converged_flag
      logical, intent(out) :: redo_flag
      integer, intent(out) :: ierr

! The atomic weights of H1, He4, the metals (as C12) and He3.
      double precision :: atomic_weight(4)
      data atomic_weight/1.007825d0,4.002603d0,12.0d0,3.01603d0/
! locals
      integer :: num_diffused_species, zone_index
      double precision :: delta_hydrogen, delta_helium, delta_metal, &
           delta_helium3
! amu_calc_temp holds 1/mu_ion and then 1/mu_e. Despite their names,
! ion_mean_weight_inverse and electron_mean_weight_inverse hold the
! reciprocals of those sums, mu_ion and mu_e.
      double precision :: amu_calc_temp, ion_mean_weight_inverse, &
           electron_mean_weight_inverse

      ierr = 0

! The species checked: 1-4 before the last iteration, whose mixing
! covers only those, and 1-11 on it.
      if(iteration_number.eq.rotation_settings%max_diffusion_iters)then
         num_diffused_species = n_species_basic
      else
         num_diffused_species = last_mu_species
      endif
      redo_flag = .false.
      call check_species_range(composition, 1, num_diffused_species, num_zones, &
           dt, cut_count, converged_flag, redo_flag, ierr)
      if (ierr /= 0 .or. redo_flag) return
! The mean molecular weight of a fully ionized gas (the equations in
! the header), from the second iteration on. star%envelope_inverse_atomic_weight is 1/mu_ion of
! the envelope mixture. Partial ionization is ignored; the original
! notes that massive stars need a correction for it.
      if(iteration_number.gt.1)then
         do zone_index = 1,num_zones
            delta_hydrogen = composition(i_h1,zone_index)- &
                 s%envelope_hydrogen_fraction
            delta_helium = composition(i_he4,zone_index)- &
                 s%envelope_helium_fraction
            delta_metal = composition(i_metals,zone_index)- &
                 s%envelope_metal_fraction
            delta_helium3 = composition(i_he3,zone_index)- &
                 s%envelope_he3_fraction
            amu_calc_temp = s%envelope_inverse_atomic_weight + delta_hydrogen/atomic_weight(1) + &
                 delta_helium/atomic_weight(2) + &
                 delta_metal/atomic_weight(3) + &
                 delta_helium3/atomic_weight(4)
            ion_mean_weight_inverse = 1.0d0/amu_calc_temp
            amu_calc_temp = composition(i_h1,zone_index)/atomic_weight(1)+ &
                 2.0d0*(composition(i_he3,zone_index)/atomic_weight(4) &
                 + composition(i_he4,zone_index)/atomic_weight(2)) + &
                 0.5d0*composition(i_metals,zone_index)
            electron_mean_weight_inverse = 1.0d0/amu_calc_temp
            sub_physics%mu(zone_index) = ion_mean_weight_inverse* &
                 electron_mean_weight_inverse/ &
                 (ion_mean_weight_inverse+electron_mean_weight_inverse)
         end do
      endif

      return
end subroutine check_composition

! ---------------------------------------------------------------
! The species species_begin..species_end checked after their mixing as
! check_composition checks slots 1-4 (secular_transport calls it for
! the species mixed once, after the iterations; docs/bug_sweep.md O14):
! a small negative interior abundance repaired, a larger one a halved
! sub-step and a redo, the fourth cut of the sub-step a stop (ierr = 1).
subroutine check_mixed_species(composition, species_begin, species_end, num_zones, &
     dt, cut_count, redo_flag, ierr)
      use star_info_lib, only: json
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: species_begin, species_end, num_zones
      double precision, intent(inout) :: dt
      integer, intent(inout) :: cut_count
      logical, intent(inout) :: redo_flag
      integer, intent(out) :: ierr
      logical :: converged_flag
      converged_flag = .true.
      call check_species_range(composition, species_begin, species_end, num_zones, &
           dt, cut_count, converged_flag, redo_flag, ierr)
end subroutine check_mixed_species

! ---------------------------------------------------------------
! Each species species_begin..species_end except Z (slot i_metals, not
! mixed as a unit), at the interior zones 2..num_zones-1: a mass
! fraction below 0 or above 1 is replaced by the mean of its
! neighbours (at least 0) when it is below zero by less than 1e-5 of
! the surface value; otherwise the sub-step is halved for a redo
! (redo_flag, converged_flag false), and the fourth cut of the sub-step
! stops the run (ierr = 1). Returns at the first cut.
subroutine check_species_range(composition, species_begin, species_end, num_zones, &
     dt, cut_count, converged_flag, redo_flag, ierr)
      use star_info_lib, only: json, i_metals
      use io_units
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: species_begin, species_end, num_zones
      double precision, intent(inout) :: dt
      integer, intent(inout) :: cut_count
      logical, intent(inout) :: converged_flag, redo_flag
      integer, intent(out) :: ierr
      integer :: species_index, zone_index
      character(len=message_length) :: lines(8)
      ierr = 0
      do species_index = species_begin,species_end
! Slot i_metals holds Z, which is not mixed as a unit.
         if(species_index.eq.i_metals)cycle
         do zone_index = 2,num_zones-1
            if(composition(species_index,zone_index).lt.0.0d0.or. &
                 composition(species_index,zone_index).gt.1.0d0)then
! A species mixed into a region that destroys it quickly can come out
! slightly negative from the interplay of mixing and burning. A value
! below 1e-5 of the surface abundance in size is replaced by the mean
! of its neighbours (at least 0) instead of cutting the sub-step.
               if(abs(composition(species_index,zone_index)).lt. &
                    small_abundance_fraction*composition(species_index,num_zones))then
                  composition(species_index,zone_index) = &
                       max(0.0d0,0.5d0* &
                       (composition(species_index,zone_index+1)+ &
                       composition(species_index,zone_index-1)))
                  cycle
               endif
! Otherwise halve the sub-step and redo it; the fourth cut of the
! sub-step stops the run instead.
               cut_count = cut_count + 1
               if(cut_count.gt.max_substep_cuts)then
                  lines = unwritten_line
                  write(lines,1010) species_index,zone_index, &
                       composition(species_index,zone_index)
                  call rotation_report(lines)
 1010 format(1x,39('>'),39('>')/' ERROR IN SR CHECKC'/ &
              ' ANOMALOUS COMP NUMBER',i2,' IN ZONE',i5,' ABUNDANCE ', &
              1pe12.3/' 3 ATTEMPTS AT TIMESTEP CUTTING FAILED'/ &
              'RUN STOPPED')
                  ! a nonzero ierr stops the run
                  ierr = 1
                  return
               else
                  redo_flag = .true.
                  converged_flag = .false.
                  dt = 0.5d0*dt
                  lines = unwritten_line
                  write(lines,1015)cut_count,species_index,zone_index, &
                       composition(species_index,zone_index)
                  call rotation_report(lines)
 1015 format(' ERROR IN SR CHECKC'/' TIMESTEP CUT NUMBER ',i2, &
              ' DUE TO ANOMALOUS COMP NUMBER',i2,' IN ZONE',i5, &
              ' ABUNDANCE ',1pe12.3)
                  return
               endif
            endif
         end do
      end do
end subroutine check_species_range

end module transport_checks

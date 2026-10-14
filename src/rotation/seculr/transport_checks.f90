!----------------------------------------------------------------------
! transport_checks
!----------------------------------------------------------------------
! The checks secular_transport makes after each iteration of its
! transport solve: on the new specific angular momentum, from which
! they also derive the new rotation curve (check_angular_momentum), and
! on the new abundances, from which they also derive the new mean
! molecular weight (check_composition).
module transport_checks
      use rotation_configuration, only: rotation_settings
      implicit none
      private
      public :: check_angular_momentum, check_composition

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
!  4. smooth_gradient_reversals: forces solid-body rotation where
!     d omega/dr has turned positive (in practice it never acts).
!  5. report_largest_change: on the last iteration, or at convergence,
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
!   already_converged_flag : true if the solve had converged before
!      this check (the advection-diffusion solve can say so); it marks
!      the iteration converged and suppresses the summary line.
! In/out:
!   specific_angular_momentum : j after this iteration (cm^2/s);
!      update_omega and smooth_gradient_reversals reset it to I omega/m
!      for the omega they find.
!   specific_angular_momentum_prev : j after the previous iteration;
!      set to j once compared with it.
!   cut_count : the number of sub-step cuts so far.
!   dt : the sub-step (s), which a cut would halve (see
!      reject_negative_j).
! Outputs:
!   omega : the angular velocity (rad/s).
!   moment_of_inertia, qiw : each zone's moment of inertia and its
!      derivative with respect to omega.
!   eta_squared, mean_radius : the shape quantities (eta_2 and r0) that
!      omega_from_j finds along with omega.
!   converged_flag : true if the iteration has converged.
!   redo_flag : true if the sub-step must be redone with a smaller dt;
!      in practice always false on return (see reject_negative_j and
!      smooth_gradient_reversals).
!   ierr : nonzero stops the run.
subroutine check_angular_momentum(s, log_density, specific_angular_momentum_prev, &
     specific_angular_momentum_start, log_radius, log_mass, shell_mass, &
     iteration_number, &
     am_transport_convective_flag, num_zones, dt, eta_squared, &
     moment_of_inertia, specific_angular_momentum, cut_count, &
     converged_flag, redo_flag, omega, &
     qiw, mean_radius, &
     already_converged_flag, ierr)
      use omega_profile, only: solid_body_omega, omega_from_j, &
           default_inertia_tolerance
      use rotation_scratch_lib

      use star_info_lib, only: star_info, json
      use io_units
      implicit none
      type(star_info), intent(inout) :: s

      double precision, intent(inout) :: log_density(json)
      double precision, intent(inout) :: specific_angular_momentum_prev(json)
      double precision, intent(in) :: specific_angular_momentum_start(json)
      double precision, intent(inout) :: log_radius(json), log_mass(json), &
           shell_mass(json)
      integer, intent(in) :: iteration_number
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: dt
      double precision, intent(inout) :: eta_squared(json), &
           moment_of_inertia(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      integer, intent(inout) :: cut_count
      logical, intent(inout) :: converged_flag
      logical, intent(out) :: redo_flag
      double precision, intent(inout) :: omega(json)
      double precision, intent(inout) :: qiw(json), mean_radius(json)
      logical, intent(in) :: already_converged_flag
      integer, intent(out) :: ierr
! locals
! The history of measure_convergence, one entry per iteration: room for
! 16 iterations, and max_diffusion_iters is not checked against it.
      double precision :: max_delta_j_by_iter(16)
      integer :: max_delta_j_zone_by_iter(16)
      integer :: zone_index, scan_index
      double precision :: delta_j_fraction
      double precision :: inertia_tolerance, shape_tolerance
      integer :: zone_bottom, zone_top
      double precision :: max_fractional_dj
      integer :: max_dj_zone

      ierr = 0

      converged_flag = .false.
      redo_flag = .false.
      if(already_converged_flag)then
         converged_flag = .true.
         redo_flag = .false.
      endif
      call reject_negative_j()
      if (ierr /= 0 .or. redo_flag) return
      call measure_convergence()
      call update_omega()
      call smooth_gradient_reversals()
      if(iteration_number.eq.rotation_settings%max_diffusion_iters.or.converged_flag) &
           call report_largest_change()

contains

! ---------------------------------------------------------------
! The first zone with j <= 0 stops the run (ierr = 1). The original
! halved dt and redid the sub-step, up to three times, before stopping;
! here cut_count > 0 always holds after the increment (it starts each
! sub-step at 0), so that branch (redo_flag, dt/2) is unreachable, and
! the stop message's "3 attempts" were not made.
subroutine reject_negative_j()
      do zone_index = 1,num_zones
         if(specific_angular_momentum(zone_index).le.0.0d0) then
            cut_count = cut_count + 1
            if(cut_count.gt.0)then
               write(6,1000) zone_index
               write(run_log_unit,1000) zone_index
 1000          format(1x,39('>'),39('<')/5x,'ERROR IN SR CHECKJ'/ &
                       5x,'NEGATIVE J/M ENCOUNTERED IN ZONE',i5, &
                       ' AND 3 ATTEMPTS AT CUTTING TIMESTEP FAILED'/ &
                       'RUN STOPPED')
               ! a nonzero ierr stops the run
               ierr = 1
               return
            else
               redo_flag = .true.
               dt = 0.5d0*dt
               write(6,1005)cut_count,zone_index
               write(run_log_unit,1005)cut_count,zone_index
 1005          format(5x,'ERROR IN SR CHECKJ'/5x,'TIMESTEP CUT,',1x, &
                       'NUMBER',i5,' DUE TO NEGATIVE J/M IN ZONE',i5)
               return
            endif
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

! ---------------------------------------------------------------
! A safeguard against reversals of the omega gradient, which the
! transport (a diffusion) should not create. Scanning inward from the
! surface: where omega rises outward across an interface,
! omega(i) - omega(i-1) > 0, but did not in the previous model
! (star%old_omega), the two zones, together with the rest of any
! convection zone either lies in, are forced to rotate as one solid
! body (solid_body_omega). The region then widens by a zone (and that
! zone's convection zone) at each edge where the new omega has made a
! fresh reversal, and the solid body is redone, until no edge moves.
! The entry test is omega(i) - omega(i-1) > 1 rad/s, inherited from the
! F77 checkj.f. Stellar rotation rates are far below 1 rad/s (the
! Sun's is 2.9e-6 rad/s), so the safeguard practically never fires.
! The test on the previous model and the edge tests use 1e-15 rad/s.
subroutine smooth_gradient_reversals()
! zone_index is the outer zone of the interface tested; zone_bottom the
! inner edge of the last region treated, below which the scan resumes.
      zone_index = num_zones
      zone_bottom = num_zones
      zone_scan: do
! omega rises outward across the interface (by more than 1 rad/s).
      if(omega(zone_index)-omega(zone_index-1).gt.1.0d0)then
! It did in the previous model too: leave it and move one zone in.
         if(s%old_omega(zone_index)-s%old_omega(zone_index-1).gt.1.0d-15)then
            zone_index = zone_bottom-1
            zone_bottom = zone_index
            if(zone_index.gt.1) cycle zone_scan
         else
! d omega/dr has changed sign, which signals an error in the transport:
! force the offending zones to solid-body rotation. zone_top is the
! upper one, extended to the top of its convection zone if it lies in
! one.
         zone_top = zone_index
         if(am_transport_convective_flag(zone_top) .and. zone_top.lt.num_zones) then
            do scan_index = zone_top + 1,num_zones
               if(.not.am_transport_convective_flag(scan_index)) exit
            end do
            zone_top = scan_index - 1
         endif
! zone_bottom is the lower one, extended likewise to the bottom of its
! convection zone.
         zone_bottom = zone_index - 1
         if(am_transport_convective_flag(zone_bottom) .and. zone_bottom.gt.1) then
            do scan_index = zone_bottom - 1,1,-1
               if(.not.am_transport_convective_flag(scan_index)) exit
            end do
            zone_bottom = scan_index + 1
         endif
! Solid-body rotation from zone_bottom to zone_top.
         call solid_body_omega(log_density,specific_angular_momentum,log_radius,log_mass, &
                    shell_mass,zone_bottom,zone_top,eta_squared, &
                    moment_of_inertia,omega,qiw,mean_radius,num_zones,s%disk_gate_age_gyr)
! Widen the region while the redistribution makes new reversals at its
! edges. redo_flag serves here as the flag that an edge moved; it ends
! false, so it never asks for a redo of the sub-step.
         do
         redo_flag = .false.
! A new reversal at the interface below zone_bottom: take in the zone
! below (and its convection zone).
         if(zone_bottom.gt.1) then
            if(omega(zone_bottom)-omega(zone_bottom-1).gt.1.0d-15)then
               if(s%old_omega(zone_bottom)-s%old_omega(zone_bottom-1).lt.1.0d-15)then
                  redo_flag = .true.
                  zone_bottom = zone_bottom - 1
                  if(am_transport_convective_flag(zone_bottom) .and. zone_bottom.gt.1) then
                     do scan_index = zone_bottom - 1,1,-1
                        if(.not.am_transport_convective_flag(scan_index)) exit
                     end do
                     zone_bottom = scan_index + 1
                  endif
               endif
            endif
         endif
! A new reversal at the interface above zone_top: take in the zone
! above (and its convection zone).
         if(zone_top.lt.num_zones) then
            if(omega(zone_top+1)-omega(zone_top).gt.1.0d-15)then
               if(s%old_omega(zone_top+1)-s%old_omega(zone_top).lt.1.0d-15)then
                  redo_flag = .true.
                  zone_top = zone_top+1
                  if(am_transport_convective_flag(zone_top) .and. zone_top.lt.num_zones) then
                     do scan_index = zone_top+1,num_zones
                        if(.not.am_transport_convective_flag(scan_index)) exit
                     end do
                     zone_top = scan_index - 1
                  endif
               endif
            endif
         endif
! Done when neither edge moved; otherwise redo the solid body over the
! widened region.
         if(.not. redo_flag) exit
            call solid_body_omega(log_density,specific_angular_momentum,log_radius,log_mass, &
                       shell_mass,zone_bottom,zone_top,eta_squared, &
                       moment_of_inertia,omega,qiw,mean_radius,num_zones,s%disk_gate_age_gyr)
         end do
         if(iteration_number.eq.rotation_settings%max_diffusion_iters.or.converged_flag) &
              write(*,120)zone_bottom,zone_top,iteration_number
  120    format(5x,'OMEGA GRADIENT REVERSAL BETWEEN ZONES ', &
                 i5,' AND ',i5,' ITERATION ',i5)
         endif
      endif
! Continue the scan below the region just treated.
      zone_index = zone_bottom - 1
      zone_bottom = zone_index
      if(zone_index.le.1) exit zone_scan
      end do zone_scan
end subroutine smooth_gradient_reversals

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
      if(.not.already_converged_flag)then
         write(*,160)max_fractional_dj,max_dj_zone, &
              (max_delta_j_by_iter(scan_index), &
              max_delta_j_zone_by_iter(scan_index),scan_index=1, &
              iteration_number)
  160       format(' MAX D(J/M)/(J/M)',1pe12.3,' AT PT.',i5, &
                 ' BY ITERATION'/5(1x,e11.3,i4))
      endif
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
! and star%amuenv = sum X_i/A_i of that mixture.
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
!   mid_scr%mu : the new mean molecular weight (from the second
!      iteration on).
!   ierr : nonzero stops the run.
subroutine check_composition(s, composition, iteration_number, num_zones, &
     dt, cut_count, converged_flag, redo_flag, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: star_info, json, i_h1, i_he4, i_metals, i_he3
      use io_units
      implicit none
      type(star_info), intent(inout) :: s

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
      integer :: num_diffused_species, species_index, zone_index
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
         num_diffused_species = 11
      else
         num_diffused_species = 4
      endif
      redo_flag = .false.
      do species_index = 1,num_diffused_species
! Slot i_metals holds Z, which is not mixed as a unit.
         if(species_index.eq.i_metals)cycle
         do zone_index = 2,num_zones-1
            if(composition(species_index,zone_index).lt.0.0d0.or. &
                 composition(species_index,zone_index).gt.1.0d0)then
! A species mixed into a region that destroys it quickly can come out
! slightly negative from the interplay of mixing and burning. A value
! below 1e-5 of the surface abundance in size is replaced by the mean
! of its neighbours (at least 0) instead of cutting the sub-step. (The
! end-zone branch is unreachable: the loop covers zones 2 to
! num_zones-1.)
               if(abs(composition(species_index,zone_index)).lt. &
                    1.0d-5*composition(species_index,num_zones))then
                  if(zone_index.eq.1.or.zone_index.eq.num_zones)then
                     composition(species_index,zone_index) = 0.0d0
                  else
                     composition(species_index,zone_index) = &
                          max(0.0d0,0.5d0* &
                          (composition(species_index,zone_index+1)+ &
                          composition(species_index,zone_index-1)))
                  endif
                  cycle
               endif
! Otherwise halve the sub-step and redo it; the fourth cut of the
! sub-step stops the run instead.
               cut_count = cut_count + 1
               if(cut_count.gt.3)then
                  write(6,1010) species_index,zone_index, &
                       composition(species_index,zone_index)
                  write(run_log_unit,1010) species_index,zone_index, &
                       composition(species_index,zone_index)
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
                  write(6,1015)cut_count,species_index,zone_index, &
                       composition(species_index,zone_index)
                  write(run_log_unit,1015)cut_count,species_index, &
                       zone_index,composition(species_index,zone_index)
 1015 format(' ERROR IN SR CHECKC'/' TIMESTEP CUT NUMBER ',i2, &
              ' DUE TO ANOMALOUS COMP NUMBER',i2,' IN ZONE',i5, &
              ' ABUNDANCE ',1pe12.3)
                  return
               endif
            endif
         end do
      end do
! The mean molecular weight of a fully ionized gas (the equations in
! the header), from the second iteration on. star%amuenv is 1/mu_ion of
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
            amu_calc_temp = s%amuenv + delta_hydrogen/atomic_weight(1) + &
                 delta_helium/atomic_weight(2) + &
                 delta_metal/atomic_weight(3) + &
                 delta_helium3/atomic_weight(4)
            ion_mean_weight_inverse = 1.0d0/amu_calc_temp
            amu_calc_temp = composition(i_h1,zone_index)/atomic_weight(1)+ &
                 2.0d0*(composition(i_he3,zone_index)/atomic_weight(4) &
                 + composition(i_he4,zone_index)/atomic_weight(2)) + &
                 0.5d0*composition(i_metals,zone_index)
            electron_mean_weight_inverse = 1.0d0/amu_calc_temp
            mid_scr%mu(zone_index) = ion_mean_weight_inverse* &
                 electron_mean_weight_inverse/ &
                 (ion_mean_weight_inverse+electron_mean_weight_inverse)
         end do
      endif

      return
end subroutine check_composition

end module transport_checks

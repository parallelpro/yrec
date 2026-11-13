!----------------------------------------------------------------------
! henyey_iteration (formerly crrect)
!----------------------------------------------------------------------
! YREC's Newton-Raphson relaxation driver: the Henyey structure solve
! for one iteration level, in two calls that evolve_step's solve_level
! makes in order, with the composition step between them at levels >= 3:
!
!   refresh_surface_bc   the surface boundary condition
!                        (compute_surface_boundary), when
!                        recompute_surface_bc asks for it. It runs
!                        before the composition step because it reads
!                        the surface composition.
!   henyey_iterate       the correction loop. Each iteration linearizes
!                        the structure equations at every shell
!                        (compute_henyey_coefficients), builds the outer
!                        boundary-condition rows from the envelope fit
!                        (star%envelope_fit_coeffs), solves for the
!                        P, T, R and L corrections (henyey_solve), tests
!                        them against the tolerances (physics%solver_tolerances),
!                        applies them, and with rotation updates the
!                        rotation curve and the rotational P and T
!                        factors (omega_from_j, rotation_shape_factors).
!
! The two share the Saha iteration counter of the eos calls (saha_iterations)
! and whether compute_surface_boundary integrated a new envelope
! (envelope_recomputed_flag). Legacy crrect did both, with the
! composition step in between.
!
! Nonzero ierr propagates a callee failure to evolve_step.
! Assumes log Teff in star%log_Teff and log L/Lsun in star%log_L_lsun.
module henyey_iteration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use henyey_elimination, only: henyey_solve
      use surface_boundary, only: compute_surface_boundary
      use structure_equations, only: compute_henyey_coefficients, &
           luminosity_channel_sum
      use omega_profile, only: omega_from_j
      use rotating_shape, only: rotation_shape_factors
      use star_info_lib, only: star_info, json
      implicit none
      private
      public :: refresh_surface_bc, henyey_iterate

contains

!----------------------------------------------------------------------
! refresh_surface_bc: start an iteration level. Sets log L
! (star%log_L_lsun), the log mass fraction inside the fitting point,
!   log_fit_point_mass_fraction = log m(nz) - log M,
! and resets the Saha iteration counter. When recompute_surface_bc is set
! (forced on for a new triangle, or for a reset one at level 2), it
! rebuilds the envelope fit with compute_surface_boundary.
!----------------------------------------------------------------------
!   star reads: fp_rot, ft_rot, log_mass, log_total_mass,
!     luminosity_lsun, nz, physics%rotation_active,
!     recompute_envelope_triangle, reset_triangle
!   star writes: log_fit_point_mass_fraction, log_L_lsun
subroutine refresh_surface_bc(s, iteration_level, recompute_surface_bc, &
     saha_iterations, envelope_recomputed_flag, ierr)
      use math_lib
      use star_info_lib, only: i_c12, i_c13, i_h1, i_he3, i_he4, i_metals, &
           i_n14, i_n15, i_o16, i_o17, i_o18
      implicit none
      type(star_info), intent(inout) :: s

      integer, intent(in) :: iteration_level
      logical, intent(inout) :: recompute_surface_bc
      integer, intent(out) :: saha_iterations
      logical, intent(out) :: envelope_recomputed_flag
      integer, intent(out) :: ierr

      real(dp) :: surface_pressure_rotation_factor, &
           surface_temperature_rotation_factor
      integer :: jerr

      ierr = 0
! False unless compute_surface_boundary integrates a new envelope.
      envelope_recomputed_flag = .false.
      s%log_L_lsun = log10(s%luminosity_lsun(s%nz))
      saha_iterations = 0
      s%log_fit_point_mass_fraction = s%log_mass(s%nz) - s%log_total_mass
      if (s%recompute_envelope_triangle .or. (s%reset_triangle .and. iteration_level == 2)) &
           recompute_surface_bc = .true.
! The surface rotation factors fp and ft for the envelope integrations
! (1 without rotation).
      if (s%physics%rotation_active.and.recompute_surface_bc) then
       surface_pressure_rotation_factor = s%fp_rot(s%nz)
       surface_temperature_rotation_factor = &
            s%ft_rot(s%nz)
      else
       surface_pressure_rotation_factor = 1.0d0
       surface_temperature_rotation_factor = 1.0d0
      end if
! The envelope integrations use the surface composition (X, Z and the
! Debye-Hueckel species of the outermost shell) and take log P(nz) as
! their guess at the fitting-point pressure (compute_surface_boundary
! reads them from the star).
      if (recompute_surface_bc) then
       call compute_surface_boundary(s, surface_pressure_rotation_factor, &
            surface_temperature_rotation_factor, saha_iterations, &
            envelope_recomputed_flag, jerr)
! An envelope or atmosphere integration failure (including
! envelope_termination, a negative code) comes back in jerr and goes up
! to the driver.
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
      end if
end subroutine refresh_surface_bc

!----------------------------------------------------------------------
! henyey_iterate: the Newton correction loop of one iteration level, at
! most max_iterations passes. Each pass linearizes the structure
! equations, accepts the model if the residuals are already small,
! closes the system with the surface rows, solves for the corrections,
! tests them and applies them, scaled by the correction factor.
! Inputs: delta_time, the time step in s; envelope_zone_index, from which
!   zone up the entropy term of compute_henyey_coefficients adds the
!   accretion shifts of T and P (see solve_level for its value).
! Outputs: converged, true once the model has converged (on entry, the
!   previous level's result); star%model_diverged_flag, true if a
!   correction exceeded its limit (false on every other return);
!   iterations_done, increased by one for each applied correction.
!----------------------------------------------------------------------
!   star reads: none
!   star writes: model_diverged_flag
subroutine henyey_iterate(s, delta_time, max_iterations, converged, &
     recompute_surface_bc, envelope_recomputed_flag, &
     saha_iterations, envelope_zone_index, iterations_done, ierr)
      implicit none
      type(star_info), intent(inout) :: s

      real(dp), intent(in) :: delta_time
      integer, intent(in) :: max_iterations
      logical, intent(inout) :: converged
      logical, intent(in) :: recompute_surface_bc, envelope_recomputed_flag
      integer, intent(inout) :: saha_iterations
      integer, intent(in) :: envelope_zone_index
      integer, intent(inout) :: iterations_done
      integer, intent(out) :: ierr

      integer :: iter
! The largest L correction (for the helium flash's correction factor) and
! the factor applied to the corrections.
      real(dp) :: luminosity_correction_max, correction_factor

      ierr = 0
! The divergence flag is set on every path: the early return below and
! an error return leave it false, as do the passes until
! test_henyey_convergence sets it.
      s%model_diverged_flag = .false.
! If the model converged at the previous level and the envelope triangle
! was checked without any new envelope, there is nothing to correct.
      if (recompute_surface_bc .and. .not.envelope_recomputed_flag .and. &
           converged) return
      converged = .false.
      do iter = 1,max_iterations
         call linearize_structure_equations(s, delta_time, saha_iterations, &
              envelope_zone_index, ierr)
         if (ierr /= 0) return
         call renormalize_luminosity_channels(s)
         call accept_small_residuals(s, iter, iterations_done, converged)
         if (converged) return
         call form_surface_bc_rows(s)
         call solve_henyey_corrections(s, ierr)
         if (ierr /= 0) return
         call find_largest_corrections(s)
         call test_henyey_convergence(s, converged, luminosity_correction_max)
         call set_henyey_correction_factor(s, luminosity_correction_max, converged, &
              iterations_done, correction_factor)
         call log_henyey_corrections(s, converged, iterations_done, correction_factor)
         if (s%model_diverged_flag) return
         call apply_henyey_corrections(s, correction_factor)
         call update_rotation_after_corrections(s, ierr)
         if (ierr /= 0) return
         iterations_done = iterations_done + 1
         if (converged) return
      end do
end subroutine henyey_iterate

! ---------------------------------------------------------------
! Solve for the corrections (henyey_solve); they replace the residuals
! in elim_rhs. A zero pivot of the surface system stops the run.
!   star reads: elim_coeff, nz, surface_bc
!   star writes: elim_rhs
subroutine solve_henyey_corrections(s, ierr)
      use core_messages, only: core_report, message_length, unwritten_line
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      character(len=message_length) :: lines(1)
      call henyey_solve(s%nz,s%elim_coeff,s%elim_rhs,s%surface_bc,ierr)
      if (ierr /= 0) then
         lines = unwritten_line
         lines(1) = ' henyey_solve: a zero pivot in the surface boundary system; run stopped'
         call core_report(lines)
      end if
end subroutine solve_henyey_corrections

! ---------------------------------------------------------------
! Zero the residual maxima and linearize the structure equations at
! every shell, with the forward elimination (compute_henyey_coefficients);
! a failure leaves ierr nonzero.
!   star reads: none
!   star writes: max_residual
subroutine linearize_structure_equations(s, delta_time, saha_iterations, envelope_zone_index, ierr)
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: delta_time
      integer, intent(inout) :: saha_iterations
      integer, intent(in) :: envelope_zone_index
      integer, intent(out) :: ierr
      integer :: jerr
      ierr = 0
       s%max_residual(1) = 0.0d0
       s%max_residual(2) = 0.0d0
       s%max_residual(3) = 0.0d0
       s%max_residual(4) = 0.0d0
       call compute_henyey_coefficients(s, delta_time, saha_iterations, &
            envelope_zone_index, jerr)
       if (jerr /= 0) then
       ! Propagate the failure to the caller.
          ierr = jerr
          return
       end if
end subroutine linearize_structure_equations

! ---------------------------------------------------------------
! Scale the luminosity channels so that they sum to the surface
! luminosity,
!   luminosity_breakdown_lsun(j) *= L(nz)/sum_j luminosity_breakdown_lsun(j),
! except during the helium flash, when the star is out of thermal
! equilibrium. The sum is structure_equations' luminosity_channel_sum,
! with the extra energy sources' channel when there are any (so that
! they are not folded into the others). This is the per-iteration twin
! of observables' renormalize_luminosity_breakdown; it also guards a
! nonpositive sum, and is kept separate on purpose
! (audit/readability-sweep-2026-09-03/reports/henyey-wave2.md item 6).
!   star reads: luminosity_lsun, nz, physics%helium_flash_active
!   star writes: luminosity_breakdown_lsun
subroutine renormalize_luminosity_channels(s)
      use energy_sources, only: num_energy_sources
      use star_info_lib, only: n_lum_channels
      type(star_info), intent(inout) :: s
      real(dp) :: total_luminosity_terms, temp
      integer :: j
       total_luminosity_terms = luminosity_channel_sum(s%luminosity_breakdown_lsun, &
            num_energy_sources() > 0)
       if (.not.s%physics%helium_flash_active .and. total_luminosity_terms > 0.0d0) &
            then
          temp = s%luminosity_lsun(s%nz)/total_luminosity_terms
          do j = 1,n_lum_channels
             s%luminosity_breakdown_lsun(j) = s%luminosity_breakdown_lsun(j)*temp
          end do
       end if
end subroutine renormalize_luminosity_channels

! ---------------------------------------------------------------
! From the second pass on, the model has converged, without a solve, if
! the largest residual |F_k| of every difference equation is within
! tolerance: P and R against htoler(5,1), T and L against htoler(5,2)
! (the L residual is absolute, in Lsun). The boundary conditions'
! residuals are not checked.
!   star reads: max_residual, physics%solver_tolerances
!   star writes: newton_iterations
subroutine accept_small_residuals(s, iter, iterations_done, converged)
      use run_messages, only: solver_diagnostics
      use core_messages, only: core_log, message_length, unwritten_line
      type(star_info), intent(inout) :: s
      integer, intent(in) :: iter, iterations_done
      logical, intent(inout) :: converged
      character(len=message_length) :: lines(2)
      integer :: j
       if (iter > 1) then
          if (s%max_residual(1) <= s%physics%solver_tolerances(5,1) .and. &
               s%max_residual(2) <= s%physics%solver_tolerances(5,2) .and. &
               s%max_residual(3) <= s%physics%solver_tolerances(5,1) .and. &
               s%max_residual(4) <= s%physics%solver_tolerances(5,2)) then
             s%newton_iterations = iterations_done + 1
             if (solver_diagnostics()) then
                lines = unwritten_line
                write(lines(1),20) (s%max_residual(j),j=1,4)
   20           format(' R.H.S. BELOW TOLERANCES--P',1PE9.2,'  T ',E9.2, &
              '  R ',E9.2,'  L ',E9.2)
                write(lines(2),75) iterations_done+1
                call core_log(lines)
             end if
             converged = .true.
             return
          end if
       end if
! A copy of log_henyey_corrections' format 75, since a statement label
! belongs to one procedure.
   75 format(10X,'MODEL CONVERGED AFTER ',I4,'  ITERATIONS')
end subroutine accept_small_residuals

! ---------------------------------------------------------------
! The outer boundary-condition rows from the envelope fit of
! compute_surface_boundary. At the fitting point, the outermost shell N,
!   log R = cf(1) log P + cf(2) log T + cf(3)
!   log L = cf(4) log P + cf(5) log T + cf(6)
! with cf = star%envelope_fit_coeffs. Linearized in henyey_solve's form,
! with s = star%surface_bc,
!   dlog R = s(3) - s(1) dlog P - s(2) dlog T
!   dL     = s(6) - s(4) dlog P - s(5) dlog T
! s(3) and s(6) carry the fits' residuals, and the L row is multiplied
! by dL/dlog L = ln10 L.
!   star reads: envelope_fit_coeffs, log_L_lsun, logP, logR, logT,
!     luminosity_lsun, nz
!   star writes: surface_bc
subroutine form_surface_bc_rows(s)
      use physical_constants, only: ln10
      type(star_info), intent(inout) :: s
      real(dp) :: temp
       s%surface_bc(1) = -s%envelope_fit_coeffs(1)
       s%surface_bc(2) = -s%envelope_fit_coeffs(2)
       s%surface_bc(3) = s%envelope_fit_coeffs(1)*s%logP(s%nz) + &
            s%envelope_fit_coeffs(2)*s%logT(s%nz) + &
            s%envelope_fit_coeffs(3) - s%logR(s%nz)
       temp = ln10*s%luminosity_lsun(s%nz)
       s%surface_bc(4) =-temp*s%envelope_fit_coeffs(4)
       s%surface_bc(5) =-temp*s%envelope_fit_coeffs(5)
       s%surface_bc(6) = temp*(s%envelope_fit_coeffs(4)*s%logP(s%nz)+ &
            s%envelope_fit_coeffs(5)*s%logT(s%nz)+ &
            s%envelope_fit_coeffs(6)-s%log_L_lsun)
end subroutine form_surface_bc_rows

! ---------------------------------------------------------------
! The largest |correction| of each variable and the shell where it
! occurs (ties go to the outer shell). For L it is the smaller of |dL|
! and |dL/L|, except at shell 1, which enters with |dL|.
!   star reads: elim_rhs, luminosity_lsun, nz
!   star writes: max_correction_index, max_residual
subroutine find_largest_corrections(s)
      type(star_info), intent(inout) :: s
      real(dp) :: test
      integer :: i, j
       do j = 1,4
          s%max_residual(j) = abs(s%elim_rhs(j,1))
          s%max_correction_index(j) = 1
       end do
       do i = 2,s%nz
          test = abs(s%elim_rhs(1,i))
          if (s%max_residual(1) <= test) then
             s%max_residual(1) = test
             s%max_correction_index(1) = i
          end if
          test = abs(s%elim_rhs(2,i))
          if (s%max_residual(2) <= test) then
             s%max_residual(2) = test
             s%max_correction_index(2) = i
          end if
          test = abs(s%elim_rhs(3,i))
          if (s%max_residual(3) <= test) then
             s%max_residual(3) = test
             s%max_correction_index(3) = i
          end if
          test = min(abs(s%elim_rhs(4,i)),abs(s%elim_rhs(4,i)/s%luminosity_lsun(i)))
          if (s%max_residual(4) <= test) then
             s%max_residual(4) = test
             s%max_correction_index(4) = i
          end if
       end do
end subroutine find_largest_corrections

! ---------------------------------------------------------------
! Converged if every largest correction is below its tolerance
! htoler(k,1); too large if any exceeds its limit htoler(k,2). The same
! test holds during the helium flash. Then keep the signed correction at
! each maximum, for the log line.
!   star reads: elim_rhs, max_correction_index,
!     physics%solver_tolerances
!   star writes: max_residual, model_diverged_flag
subroutine test_henyey_convergence(s, converged, luminosity_correction_max)
      type(star_info), intent(inout) :: s
      logical, intent(out) :: converged
      real(dp), intent(out) :: luminosity_correction_max
      integer :: j, max_correction_pos
! The largest L correction, for the helium flash's correction factor.
       luminosity_correction_max = s%max_residual(4)
       converged = s%max_residual(1) < s%physics%solver_tolerances(1,1) .and. &
            s%max_residual(2) < s%physics%solver_tolerances(2,1) &
         .and. s%max_residual(3) < s%physics%solver_tolerances(3,1) .and. &
         s%max_residual(4) < s%physics%solver_tolerances(4,1)
       s%model_diverged_flag = s%max_residual(1) > s%physics%solver_tolerances(1,2) .or. &
            s%max_residual(2) > s%physics%solver_tolerances(2,2) &
         .or. s%max_residual(3) > s%physics%solver_tolerances(3,2) .or. &
         s%max_residual(4) > s%physics%solver_tolerances(4,2)
       do j = 1,4
          max_correction_pos = s%max_correction_index(j)
          s%max_residual(j) = s%elim_rhs(j,max_correction_pos)
       end do
end subroutine test_henyey_convergence

! ---------------------------------------------------------------
! Ramp the under-relaxation factor and pick the correction factor (eased
! in during the helium flash); count the iterations if converged.
!   star reads: physics%helium_flash_active,
!     physics%solver_underrelax_increment,
!     physics%solver_underrelax_start
!   star writes: newton_correction_factor, newton_iterations
subroutine set_henyey_correction_factor(s, luminosity_correction_max, converged, iterations_done, correction_factor)
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: luminosity_correction_max
      logical, intent(in) :: converged
      integer, intent(in) :: iterations_done
      real(dp), intent(out) :: correction_factor
       if (s%physics%solver_underrelax_start > 0.0d0) s%newton_correction_factor = min(1.d0,s%newton_correction_factor+s%physics%solver_underrelax_increment)
! During the helium flash, apply 0.8 of the correction once the largest
! L correction is at most 0.5, and all of it once it is at most 5e-3.
       correction_factor = s%newton_correction_factor
       if (s%physics%helium_flash_active) then
          if (luminosity_correction_max <= 5.0d-1) correction_factor=8.0d-1
          if (luminosity_correction_max <= 5.0d-3) correction_factor=1.0d0
       end if
       if (converged) s%newton_iterations = iterations_done + 1
end subroutine set_henyey_correction_factor

! ---------------------------------------------------------------
! Log the corrections (solver diagnostics), and report corrections too
! large to apply.
!   star reads: luminosity_breakdown_lsun, max_correction_index,
!     max_residual, model_diverged_flag, physics%solver_tolerances
!   star writes: none
subroutine log_henyey_corrections(s, converged, iterations_done, correction_factor)
      use run_messages, only: solver_diagnostics
      use core_messages, only: core_report, core_log, message_length, unwritten_line
      use star_info_lib, only: i_lum_3alpha, i_lum_cno, i_lum_he_c, i_lum_pp1, &
           i_lum_pp2, i_lum_pp3
      type(star_info), intent(inout) :: s
      logical, intent(in) :: converged
      integer, intent(in) :: iterations_done
      real(dp), intent(in) :: correction_factor
      real(dp) :: hydrogen_burn_luminosity, helium_burn_luminosity
      character(len=message_length) :: lines(2)
      integer :: j
! The hydrogen- and helium-burning luminosities, for the log line.
       hydrogen_burn_luminosity = s%luminosity_breakdown_lsun(i_lum_pp1) + s%luminosity_breakdown_lsun(i_lum_pp2) &
            + s%luminosity_breakdown_lsun(i_lum_pp3) + s%luminosity_breakdown_lsun(i_lum_cno)
       helium_burn_luminosity= s%luminosity_breakdown_lsun(i_lum_3alpha) + s%luminosity_breakdown_lsun(i_lum_he_c)
! The correction trace is written only with solver diagnostics on.
       if (solver_diagnostics()) then
          lines = unwritten_line
          write (lines(1),60) converged,s%max_residual(4), &
               s%physics%solver_tolerances(4,1),s%max_correction_index(4)
   60       format (1X,'DEL-L/L  ',L2,1P2E12.4,5X,I5)
          write(lines(2),70)(s%max_correction_index(j),s%max_residual(j), &
               j=1,4),correction_factor,hydrogen_burn_luminosity, &
               helium_burn_luminosity,(s%luminosity_breakdown_lsun(j),j=6,7)
   70       format(' CORR',I5,'P',1PE9.2,I5,'T',E9.2,I5,'R',E9.2,I5,'L', &
            E9.2,'  F=',0PF5.3,'  E-HY',1PE10.3,' HE',E10.3,' NU',E10.3, &
            ' G',E10.3)
          call core_log(lines)
          if (converged) then
             lines = unwritten_line
             write(lines(1),75) iterations_done+1
   75          format(10X,'MODEL CONVERGED AFTER ',I4,'  ITERATIONS')
             call core_log(lines)
          end if
       end if
       if (s%model_diverged_flag) then
          lines = unwritten_line
          write(lines(1),80)
   80       format(1X,'-----CORRECTIONS EXCEEDED TOLERANCES')
          call core_report(lines)
            return
         end if
end subroutine log_henyey_corrections

! ---------------------------------------------------------------
! Apply the corrections, scaled by the correction factor f, at every
! shell: log P += f dlog P and log T += f dlog T (and the same to their
! changes over the step, log_pressure_delta and log_temperature_delta,
! which the entropy term reads), log R += f dlog R, L += f dL. Then the
! new log L, and log Teff from the envelope fit,
!   log Teff = cf(7) log P(nz) + cf(8) log T(nz) + cf(9).
!   star reads: elim_rhs, envelope_fit_coeffs, nz
!   star writes: log_L_lsun, log_pressure_delta, log_Teff,
!     log_temperature_delta, logP, logR, logT, luminosity_lsun
subroutine apply_henyey_corrections(s, correction_factor)
      use math_lib
      type(star_info), intent(inout) :: s
      real(dp), intent(in) :: correction_factor
      real(dp) :: temp
      integer :: i
       do i = 1,s%nz
          temp = correction_factor*s%elim_rhs(1,i)
          s%logP(i) = s%logP(i) + temp
          s%log_pressure_delta(i) = s%log_pressure_delta(i) + temp
          temp = correction_factor*s%elim_rhs(2,i)
          s%logT(i) = s%logT(i) + temp
          s%log_temperature_delta(i) = s%log_temperature_delta(i) + temp
          s%logR(i) = s%logR(i) + correction_factor*s%elim_rhs(3,i)
          s%luminosity_lsun(i) = s%luminosity_lsun(i) + &
               correction_factor*s%elim_rhs(4,i)
       end do
       s%log_L_lsun = log10(s%luminosity_lsun(s%nz))
       s%log_Teff = s%envelope_fit_coeffs(7)*s%logP(s%nz) + &
            s%envelope_fit_coeffs(8)*s%logT(s%nz) + &
            s%envelope_fit_coeffs(9)
end subroutine apply_henyey_corrections

! ---------------------------------------------------------------
! With rotation, the rotation curve, the P and T shape factors and the
! rotational kinetic energy of the corrected structure,
! E_rot = 0.5 omega j dm per shell.
!   star reads: am_transport_convective_flag, disk_gate_age_gyr, dm,
!     log_mass, logR, logRho, nz, physics%rotation_active
!   star writes: eta_squared, fp_rot, ft_rot, i_rot, j_rot,
!     kinetic_energy_rot, mean_gravity, mean_radius, omega
subroutine update_rotation_after_corrections(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      real(dp) :: shell_angular_momentum
      integer :: i
! di_domega of omega_from_j, which this caller does not use
      real(dp) :: qiw(json)
      ierr = 0
         if (s%physics%rotation_active) then
            call omega_from_j(s%logRho,s%j_rot,s%logR, &
                 s%log_mass,s%dm,s%am_transport_convective_flag, &
                 s%nz,s%eta_squared,s%i_rot,s%omega,qiw, &
                 s%mean_radius,s%disk_gate_age_gyr)
            call rotation_shape_factors(s%logRho,s%logR,s%log_mass,s%nz,s%omega, &
                 s%eta_squared,s%fp_rot, &
                 s%ft_rot,s%mean_gravity,s%mean_radius,ierr)
            if (ierr /= 0) return
            do i = 1,s%nz
               shell_angular_momentum = s%j_rot(i)*s%dm(i)
               s%kinetic_energy_rot(i) = 0.5d0*s%omega(i)*shell_angular_momentum
            end do
         end if
end subroutine update_rotation_after_corrections

end module henyey_iteration

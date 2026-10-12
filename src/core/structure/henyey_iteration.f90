!----------------------------------------------------------------------
! henyey_iteration (formerly crrect)
!----------------------------------------------------------------------
! YREC's Newton-Raphson relaxation driver: the Henyey structure solve
! for one iteration level, in two calls that evolve_step's solve_level
! makes in order, with the composition step between them at levels >= 3:
!
!   refresh_surface_bc   the surface boundary condition (surfbc), when
!                        recompute_surface_bc asks for it. It runs
!                        before the composition step because it reads
!                        the surface composition.
!   henyey_iterate       the correction loop. Each iteration calls
!                        henyey_coefficients to linearize the structure
!                        equations at every mesh point, builds the outer
!                        boundary-condition row from the envelope-fit
!                        coefficients (star%envelope_fit_coeffs, from
!                        surfbc), calls henyey_solve to back-solve the
!                        block-tridiagonal system for the P/T/R/L
!                        corrections, checks them against the
!                        convergence/divergence tolerances (ctrl%*),
!                        applies them, and (if rotation is active)
!                        updates the rotation curve and the rotational
!                        P/T factors via omega_from_j/rotation_shape_factors.
!
! The two share the Saha-state counter of the EOS calls (saha_state)
! and whether surfbc recomputed the envelope (envelope_recomputed_flag).
! (2026-10: one routine, crrect/henyey_iterate, used to do both and run
! the composition step itself in between.)
!
! Non-zero ierr propagates a callee failure to evolve_step.
! Assumes log Teff in star%log_Teff and log L/Lsun in star%log_L.
module henyey_iteration
      use henyey_elimination, only: henyey_solve
      use surface_boundary, only: compute_surface_boundary
      use structure_equations, only: compute_henyey_coefficients
      use omega_profile, only: omega_from_j
      use rotating_shape, only: rotation_shape_factors
      implicit none
      private
      public :: refresh_surface_bc, henyey_iterate

contains

!----------------------------------------------------------------------
! refresh_surface_bc: start an iteration level -- log L, the envelope
! mass, and (when recompute_surface_bc, forced on for a new or reset
! triangle) the surface boundary condition from surfbc.
!----------------------------------------------------------------------
subroutine refresh_surface_bc(iteration_level, start_new_triangle, &
     reset_triangle, recompute_surface_bc, tri_orientation, &
     stored_vertex_index, saha_state, envelope_recomputed_flag, ierr)
      use eos_lib, only: eos_set_debye_huckel_composition
      use star_info_lib, only: star, i_c12, i_c13, i_h1, i_he3, i_he4, &
           i_metals, i_n14, i_n15, i_o16, i_o17, i_o18
      implicit none

      integer, intent(in) :: iteration_level
      logical, intent(inout) :: start_new_triangle, reset_triangle
      logical, intent(inout) :: recompute_surface_bc
      double precision, intent(inout) :: tri_orientation
      integer, intent(inout) :: stored_vertex_index
      integer, intent(out) :: saha_state
      logical, intent(out) :: envelope_recomputed_flag
      integer, intent(out) :: ierr

      double precision :: hydrogen_fraction, metal_fraction, &
           log10_pressure_limit
      double precision :: surface_pressure_rotation_factor, &
           surface_temperature_rotation_factor
      integer :: jerr

      ierr = 0
! .false. unless surfbc runs and says otherwise
      envelope_recomputed_flag = .false.
      star%log_L = log10(star%luminosity_lsun(star%nz))
! RESET THE SAHA-STATE COUNTER (2026 W2: the envelope/atmosphere
! integrand call counters kenv/katm that used to be zeroed here were
! never read; removed)
      saha_state = 0
      star%senv = star%log_mass(star%nz) - star%log_total_mass
      if (start_new_triangle .or. (reset_triangle .and. iteration_level.eq.2)) &
           recompute_surface_bc = .true.
!  FIND NEW FP AND FT IF MODEL IS ROTATING
      if (star%job%rotation_active.and.recompute_surface_bc) then
       surface_pressure_rotation_factor = star%fp_rot(star%nz)
       surface_temperature_rotation_factor = &
            star%ft_rot(star%nz)
      else
       surface_pressure_rotation_factor = 1.0d0
       surface_temperature_rotation_factor = 1.0d0
      endif
! SET UP SURFACE BOUNDARY CONDITIONS-2ND AND 3RD LEVELS OF ITER ONLY
! FIND ENVELOPE MASS AND SET X AND Z TO ENVELOPE VALUES
      if (recompute_surface_bc) then
       hydrogen_fraction = star%xnew
       metal_fraction = star%znew
       log10_pressure_limit = star%logP(star%nz)
       call eos_set_debye_huckel_composition(star%xa(:,star%nz))
       call compute_surface_boundary(star%trial_log_temperature,star%trial_log_luminosity,star%envelope_fit_coeffs,star%fit_point_pressure,star%fit_point_temperature, &
            star%fit_point_radius,tri_orientation,stored_vertex_index, &
            star%stored_envelope_state, &
            start_new_triangle,reset_triangle,saha_state, &
            star%log_total_mass,star%log_L, &
            star%log_Teff,hydrogen_fraction,metal_fraction, &
            surface_pressure_rotation_factor, &
            surface_temperature_rotation_factor,envelope_recomputed_flag, &
            log10_pressure_limit,star%convective_flag,star%nz,jerr)
! 2026 numerics-gate opt-in: surfbc surfaces envelope/atmosphere
! integration failures (incl. numerics_termination, negative) via
! ierr; propagate to the driver.
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
      endif
end subroutine refresh_surface_bc

!----------------------------------------------------------------------
! henyey_iterate: the Newton correction loop of one iteration level.
! envelope_zone_index gates the accretion term of the entropy equation
! in henyey_coefficients (zones at or above it get the accreted-layer
! T and P changes); see solve_level for its value.
! INPUTS  delta_time = TIME STEP IN SECONDS
!         max_iterations = ITERATION NUMBER
! OUTPUTS converged = .T. IF MODEL HAS CONVERGED
!         corrections_too_large = .T. IF CORRECTIONS ARE TOO LARGE
!----------------------------------------------------------------------
subroutine henyey_iterate(delta_time, max_iterations, converged, &
     corrections_too_large, recompute_surface_bc, envelope_recomputed_flag, &
     saha_state, envelope_zone_index, dlnrho_dlnt, dlnrho_dlnp, &
     iterations_done, ierr)
      use star_info_lib, only: star, i_lum_3alpha, i_lum_cno, i_lum_grav, &
           i_lum_he_c, i_lum_neu, i_lum_pp1, i_lum_pp2, i_lum_pp3, n_lum_channels
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: delta_time
      integer, intent(in) :: max_iterations
      logical, intent(inout) :: converged
      logical, intent(out) :: corrections_too_large
      logical, intent(in) :: recompute_surface_bc, envelope_recomputed_flag
      integer, intent(inout) :: saha_state
      integer, intent(in) :: envelope_zone_index
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp
      integer, intent(inout) :: iterations_done
      integer, intent(out) :: ierr

      integer :: i, j, iter, max_correction_pos
! sign of the fcorr0 control: > 0 means the under-relaxation factor
! star%job%fcorr (started at |fcorr0| - fcorri by evolve_step) is
! ramped up by fcorri each iteration; <= 0 holds it fixed
      logical :: underrelax_ramp_active
      double precision :: temp, test, total_luminosity_terms, &
           luminosity_correction_max, correction_factor
      double precision :: hydrogen_burn_luminosity, helium_burn_luminosity
      double precision :: shell_angular_momentum
      integer :: jerr

      ierr = 0
      underrelax_ramp_active = star%ctrl%fcorr0.gt.0.0d0
!  IF MODEL CONVERGED ON THE PREVIOUS LEVEL OF ITERATION AND
!  THE TRIANGLE WAS CHECKED AND NOT FLIPPED,SKIP CORRECTION ROUTINE.
      if (recompute_surface_bc .and. .not.envelope_recomputed_flag .and. &
           converged) return
      converged = .false.
      do iter = 1,max_iterations
       call linearize_structure_equations()
       if (ierr /= 0) return
       call renormalize_luminosity_channels()
       call accept_small_residuals()
       if (converged) return
       call form_surface_bc_rows()
! DO BACK SOLUTION FOR CORRECTIONS
       call henyey_solve(star%nz,star%elim_coeff,star%elim_rhs,star%surface_bc)
       call find_largest_corrections()
       call test_henyey_convergence()
       call set_henyey_correction_factor()
       call log_henyey_corrections()
       if (corrections_too_large) return
       call apply_henyey_corrections()
       call update_rotation_after_corrections()
       if (ierr /= 0) return
       iterations_done = iterations_done + 1
       if (converged) return
      end do

      return

contains

! ---------------------------------------------------------------
! Zero the residual maxima and linearize the structure equations at
! every shell (henyey_coefficients); a failure leaves ierr nonzero.
subroutine linearize_structure_equations()
! DO HENYEY REDUCTION
       star%max_residual(1) = 0.0d0
       star%max_residual(2) = 0.0d0
       star%max_residual(3) = 0.0d0
       star%max_residual(4) = 0.0d0
! (2026 de-tramp: the star%-shaped slots are read/written by
! henyey_coefficients directly; only the per-call locals remain.)
       call compute_henyey_coefficients(delta_time,dlnrho_dlnt,dlnrho_dlnp,saha_state, &
            envelope_zone_index, jerr)
       if (jerr /= 0) then
       ! 2026 (phase five, step B): propagate instead of stopping
          ierr = jerr
          return
       end if
end subroutine linearize_structure_equations

! ---------------------------------------------------------------
! Scale the luminosity channels to sum to the surface luminosity
! (not during the helium flash).
subroutine renormalize_luminosity_channels()
! RENORMALIZE TLUMX-S
!CC   TAKE OUT RENORMALIZATION DURING HE FLASH (NON-THERMAL EQUALIBRIUM)
! (per-iteration twin of observables's renormalize_luminosity_breakdown;
! this one also guards total > 0 -- kept separate on purpose, see
! reports/henyey-wave2.md item 6)
       total_luminosity_terms = star%luminosity_breakdown(i_lum_pp1)+star%luminosity_breakdown(i_lum_pp2)+ &
            star%luminosity_breakdown(i_lum_pp3)+star%luminosity_breakdown(i_lum_cno)+star%luminosity_breakdown(i_lum_3alpha)+ &
            star%luminosity_breakdown(i_lum_neu)+star%luminosity_breakdown(i_lum_grav)+star%luminosity_breakdown(i_lum_he_c)
       if (.not.star%ctrl%helium_flash_active .and. total_luminosity_terms.gt.0.0d0) &
            then
          temp = star%luminosity_lsun(star%nz)/total_luminosity_terms
          do j = 1,n_lum_channels
             star%luminosity_breakdown(j) = star%luminosity_breakdown(j)*temp
          end do
       endif
end subroutine renormalize_luminosity_channels

! ---------------------------------------------------------------
! After the first iteration, the model has converged if the equations'
! residuals are all below their tolerances.
subroutine accept_small_residuals()
! CHECK ON SIGNIFICANCE OF R.H.S. EQUATIONS FOR P AND T
! N.B.  DOES NOT CHECK DIFFERENCES IN BOUNDARY EQUATIONS
       if (iter.gt.1) then
          if (star%max_residual(1).le.star%ctrl%htoler(5,1) .and. &
               star%max_residual(2).le.star%ctrl%htoler(5,2) .and. &
               star%max_residual(3).le.star%ctrl%htoler(5,1) .and. &
               star%max_residual(4).le. star%ctrl%htoler(5,2)) then
             star%newton_iterations = iterations_done + 1
             if (solver_diagnostics()) then
                write(run_log_unit,20) (star%max_residual(j),j=1,4)
   20           format(' R.H.S. BELOW TOLERANCES--P',1PE9.2,'  T ',E9.2, &
              '  R ',E9.2,'  L ',E9.2)
                write(run_log_unit,75) iterations_done+1
             end if
             converged = .true.
             return
          endif
       endif
! (a copy of log_henyey_corrections' format 75: a label belongs to
! one procedure)
   75 format(10X,'MODEL CONVERGED AFTER ',I4,'  ITERATIONS')
end subroutine accept_small_residuals

! ---------------------------------------------------------------
! The outer boundary-condition rows from the envelope-fit coefficients.
subroutine form_surface_bc_rows()
! SET UP HENYEY COEFFICIENTS FOR SURFACE
       star%surface_bc(1) = -star%envelope_fit_coeffs(1)
       star%surface_bc(2) = -star%envelope_fit_coeffs(2)
       star%surface_bc(3) = star%envelope_fit_coeffs(1)*star%logP(star%nz) + &
            star%envelope_fit_coeffs(2)*star%logT(star%nz) + &
            star%envelope_fit_coeffs(3) - star%logR(star%nz)
       temp = ln10*star%luminosity_lsun(star%nz)
       star%surface_bc(4) =-temp*star%envelope_fit_coeffs(4)
       star%surface_bc(5) =-temp*star%envelope_fit_coeffs(5)
       star%surface_bc(6) = temp*(star%envelope_fit_coeffs(4)*star%logP(star%nz)+ &
            star%envelope_fit_coeffs(5)*star%logT(star%nz)+ &
            star%envelope_fit_coeffs(6)-star%log_L)
end subroutine form_surface_bc_rows

! ---------------------------------------------------------------
! The largest P, T, R and L corrections and the shells where they occur.
subroutine find_largest_corrections()
! CHECK ON MAXIMUM CORRECTIONS
       do j = 1,4
          star%max_residual(j) = dabs(star%elim_rhs(j,1))
          star%max_correction_index(j) = 1
       end do
       do i = 2,star%nz
          test = dabs(star%elim_rhs(1,i))
          if (star%max_residual(1).le.test) then
             star%max_residual(1) = test
             star%max_correction_index(1) = i
          endif
          test = dabs(star%elim_rhs(2,i))
          if (star%max_residual(2).le.test) then
             star%max_residual(2) = test
             star%max_correction_index(2) = i
          endif
          test = dabs(star%elim_rhs(3,i))
          if (star%max_residual(3).le.test) then
             star%max_residual(3) = test
             star%max_correction_index(3) = i
          endif
          test = dmin1(dabs(star%elim_rhs(4,i)),dabs(star%elim_rhs(4,i)/star%luminosity_lsun(i)))
          if (star%max_residual(4).le.test) then
             star%max_residual(4) = test
             star%max_correction_index(4) = i
          endif
       end do
end subroutine find_largest_corrections

! ---------------------------------------------------------------
! Converged if every largest correction is below its tolerance, too
! large if any exceeds its limit; keep the signed largest corrections.
subroutine test_henyey_convergence()
!CC   HE FLASH -- OK FOR ALL
       luminosity_correction_max = star%max_residual(4)
! LFINI = T IF MAX CORRECTIONS LESS THAN CONVERGENCE CRITERIA SET IN
! HTOLER. LARGE = T IF MAX CORRECTIONS GREATER THAN LARGEST CORRECTIONS
! ALLOWED, ALSO SET IN HTOLER
       converged = star%max_residual(1).lt.star%ctrl%htoler(1,1) .and. &
            star%max_residual(2).lt.star%ctrl%htoler(2,1) &
         .and. star%max_residual(3).lt.star%ctrl%htoler(3,1) .and. &
         star%max_residual(4).lt.star%ctrl%htoler(4,1)
       corrections_too_large = star%max_residual(1).gt.star%ctrl%htoler(1,2) .or. &
            star%max_residual(2).gt.star%ctrl%htoler(2,2) &
         .or. star%max_residual(3).gt.star%ctrl%htoler(3,2) .or. &
         star%max_residual(4).gt.star%ctrl%htoler(4,2)
       do j = 1,4
          max_correction_pos = star%max_correction_index(j)
          star%max_residual(j) = star%elim_rhs(j,max_correction_pos)
       end do
end subroutine test_henyey_convergence

! ---------------------------------------------------------------
! Ramp the under-relaxation factor and pick the correction factor (eased
! in during the helium flash); count the iterations if converged.
subroutine set_henyey_correction_factor()
       if (underrelax_ramp_active) star%job%fcorr = dmin1(1.d0,star%job%fcorr+star%ctrl%fcorri)
! HE FLASH CHANGE
       correction_factor = star%job%fcorr
       if (star%ctrl%helium_flash_active) then
          if (luminosity_correction_max.le.5.0d-1) correction_factor=8.0d-1
          if (luminosity_correction_max.le.5.0d-3) correction_factor=1.0d0
       endif
       hydrogen_burn_luminosity = star%luminosity_breakdown(i_lum_pp1) + star%luminosity_breakdown(i_lum_pp2) &
            + star%luminosity_breakdown(i_lum_pp3) + star%luminosity_breakdown(i_lum_cno)
       helium_burn_luminosity= star%luminosity_breakdown(i_lum_3alpha) + star%luminosity_breakdown(i_lum_he_c)
       if (converged) star%newton_iterations = iterations_done + 1
end subroutine set_henyey_correction_factor

! ---------------------------------------------------------------
! Log the corrections (solver diagnostics), and report corrections too
! large to apply.
subroutine log_henyey_corrections()
! 2026 log redesign: the correction trace was gated on the retired
! LCORR control; it is solver forensics, behind the diagnostics flag.
       if (solver_diagnostics()) then
          write (run_log_unit,60) converged,star%max_residual(4), &
               star%ctrl%htoler(4,1),star%max_correction_index(4)
   60       format (1X,'DEL-L/L  ',L2,1P2E12.4,5X,I5)
          write(run_log_unit,70)(star%max_correction_index(j),star%max_residual(j), &
               j=1,4),correction_factor,hydrogen_burn_luminosity, &
               helium_burn_luminosity,(star%luminosity_breakdown(j),j=6,7)
   70       format(' CORR',I5,'P',1PE9.2,I5,'T',E9.2,I5,'R',E9.2,I5,'L', &
            E9.2,'  F=',0PF5.3,'  E-HY',1PE10.3,' HE',E10.3,' NU',E10.3, &
            ' G',E10.3)
          if (converged) then
             write(run_log_unit,75) iterations_done+1
   75          format(10X,'MODEL CONVERGED AFTER ',I4,'  ITERATIONS')
          endif
       endif
       if (corrections_too_large) then
          write(terminal_unit,80)
          write(run_log_unit,80)
   80       format(1X,'-----CORRECTIONS EXCEEDED TOLERANCES')
            return
         endif
end subroutine log_henyey_corrections

! ---------------------------------------------------------------
! Apply the scaled corrections to every shell; the new log L and
! log Teff.
subroutine apply_henyey_corrections()
! APPLY CORRECTIONS
       do i = 1,star%nz
          temp = correction_factor*star%elim_rhs(1,i)
          star%logP(i) = star%logP(i) + temp
          star%log_pressure_delta(i) = star%log_pressure_delta(i) + temp
          temp = correction_factor*star%elim_rhs(2,i)
          star%logT(i) = star%logT(i) + temp
          star%log_temperature_delta(i) = star%log_temperature_delta(i) + temp
          star%logR(i) = star%logR(i) + correction_factor*star%elim_rhs(3,i)
          star%luminosity_lsun(i) = star%luminosity_lsun(i) + &
               correction_factor*star%elim_rhs(4,i)
       end do
       star%log_L = log10(star%luminosity_lsun(star%nz))
       star%log_Teff = star%envelope_fit_coeffs(7)*star%logP(star%nz) + &
            star%envelope_fit_coeffs(8)*star%logT(star%nz) + &
            star%envelope_fit_coeffs(9)
end subroutine apply_henyey_corrections

! ---------------------------------------------------------------
! With rotation, the rotation curve, the P/T shape factors and the
! rotational kinetic energy of the corrected structure.
subroutine update_rotation_after_corrections()
         if (star%job%rotation_active) then
            call omega_from_j(star%logRho,star%j_rot,star%logR, &
                 star%log_mass,star%dm,star%am_transport_convective_flag, &
                 star%nz,star%eta_squared,star%i_rot,star%omega,star%qiw, &
                 star%mean_radius)
            call rotation_shape_factors(star%logRho,star%logR,star%log_mass,star%nz,star%omega, &
                 star%eta_squared,star%fp_rot, &
                 star%ft_rot,star%mean_gravity,star%mean_radius,ierr)
            if (ierr /= 0) return
            do i = 1,star%nz
               shell_angular_momentum = star%j_rot(i)*star%dm(i)
               star%kinetic_energy_rot(i) = 0.5d0*star%omega(i)*shell_angular_momentum
            end do
         endif
end subroutine update_rotation_after_corrections

end subroutine henyey_iterate

end module henyey_iteration

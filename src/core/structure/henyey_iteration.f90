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
!                        them against the tolerances (ctrl%htoler),
!                        applies them, and with rotation updates the
!                        rotation curve and the rotational P and T
!                        factors (omega_from_j, rotation_shape_factors).
!
! The two share the Saha-state counter of the eos calls (saha_state)
! and whether compute_surface_boundary integrated a new envelope
! (envelope_recomputed_flag). Legacy crrect did both, with the
! composition step in between.
!
! Nonzero ierr propagates a callee failure to evolve_step.
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
! refresh_surface_bc: start an iteration level. Sets log L
! (star%log_L), the log mass fraction inside the fitting point,
!   senv = log m(nz) - log M,
! and resets the Saha-state counter. When recompute_surface_bc is set
! (forced on for a new triangle, or for a reset one at level 2), it
! rebuilds the envelope fit with compute_surface_boundary.
!----------------------------------------------------------------------
subroutine refresh_surface_bc(iteration_level, start_new_triangle, &
     reset_triangle, recompute_surface_bc, tri_orientation, &
     stored_vertex_index, saha_state, envelope_recomputed_flag, ierr)
      use eos_lib, only: eos_set_debye_huckel_composition
      use math_lib
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
! False unless compute_surface_boundary integrates a new envelope.
      envelope_recomputed_flag = .false.
      star%log_L = log10(star%luminosity_lsun(star%nz))
      saha_state = 0
      star%senv = star%log_mass(star%nz) - star%log_total_mass
      if (start_new_triangle .or. (reset_triangle .and. iteration_level.eq.2)) &
           recompute_surface_bc = .true.
! The surface rotation factors fp and ft for the envelope integrations
! (1 without rotation).
      if (star%job%rotation_active.and.recompute_surface_bc) then
       surface_pressure_rotation_factor = star%fp_rot(star%nz)
       surface_temperature_rotation_factor = &
            star%ft_rot(star%nz)
      else
       surface_pressure_rotation_factor = 1.0d0
       surface_temperature_rotation_factor = 1.0d0
      endif
! The envelope integrations use the surface composition (X, Z and the
! Debye-Hueckel species of the outermost shell) and take log P(nz) as
! their guess at the fitting-point pressure.
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
! An envelope or atmosphere integration failure (including
! numerics_termination, a negative code) comes back in jerr and goes up
! to the driver.
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
      endif
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
!   previous level's result); corrections_too_large, true if a correction
!   exceeded its limit; iterations_done, increased by one for each
!   applied correction.
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
! The sign of the fcorr0 control: > 0 means the under-relaxation factor
! star%job%fcorr (started at |fcorr0| - fcorri by evolve_step) is ramped
! up by fcorri each iteration; <= 0 holds it fixed.
      logical :: underrelax_ramp_active
      double precision :: temp, test, total_luminosity_terms, &
           luminosity_correction_max, correction_factor
      double precision :: hydrogen_burn_luminosity, helium_burn_luminosity
      double precision :: shell_angular_momentum
      integer :: jerr

      ierr = 0
      underrelax_ramp_active = star%ctrl%fcorr0.gt.0.0d0
! If the model converged at the previous level and the envelope triangle
! was checked without any new envelope, there is nothing to correct.
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
! Solve for the corrections; they replace the residuals in elim_rhs.
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
! every shell, with the forward elimination (compute_henyey_coefficients);
! a failure leaves ierr nonzero.
subroutine linearize_structure_equations()
       star%max_residual(1) = 0.0d0
       star%max_residual(2) = 0.0d0
       star%max_residual(3) = 0.0d0
       star%max_residual(4) = 0.0d0
       call compute_henyey_coefficients(delta_time,dlnrho_dlnt,dlnrho_dlnp,saha_state, &
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
!   luminosity_breakdown(j) *= L(nz)/sum_j luminosity_breakdown(j),
! except during the helium flash, when the star is out of thermal
! equilibrium. This is the per-iteration twin of observables'
! renormalize_luminosity_breakdown; it also guards a nonpositive sum, and
! is kept separate on purpose (reports/henyey-wave2.md item 6).
subroutine renormalize_luminosity_channels()
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
! From the second pass on, the model has converged, without a solve, if
! the largest residual |F_k| of every difference equation is within
! tolerance: P and R against htoler(5,1), T and L against htoler(5,2)
! (the L residual is absolute, in Lsun). The boundary conditions'
! residuals are not checked.
subroutine accept_small_residuals()
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
subroutine form_surface_bc_rows()
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
! The largest |correction| of each variable and the shell where it
! occurs (ties go to the outer shell). For L it is the smaller of |dL|
! and |dL/L|, except at shell 1, which enters with |dL|.
subroutine find_largest_corrections()
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
! Converged if every largest correction is below its tolerance
! htoler(k,1); too large if any exceeds its limit htoler(k,2). The same
! test holds during the helium flash. Then keep the signed correction at
! each maximum, for the log line.
subroutine test_henyey_convergence()
! The largest L correction, for the helium flash's correction factor.
       luminosity_correction_max = star%max_residual(4)
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
! During the helium flash, apply 0.8 of the correction once the largest
! L correction is at most 0.5, and all of it once it is at most 5e-3.
       correction_factor = star%job%fcorr
       if (star%ctrl%helium_flash_active) then
          if (luminosity_correction_max.le.5.0d-1) correction_factor=8.0d-1
          if (luminosity_correction_max.le.5.0d-3) correction_factor=1.0d0
       endif
! The hydrogen- and helium-burning luminosities, for the log line.
       hydrogen_burn_luminosity = star%luminosity_breakdown(i_lum_pp1) + star%luminosity_breakdown(i_lum_pp2) &
            + star%luminosity_breakdown(i_lum_pp3) + star%luminosity_breakdown(i_lum_cno)
       helium_burn_luminosity= star%luminosity_breakdown(i_lum_3alpha) + star%luminosity_breakdown(i_lum_he_c)
       if (converged) star%newton_iterations = iterations_done + 1
end subroutine set_henyey_correction_factor

! ---------------------------------------------------------------
! Log the corrections (solver diagnostics), and report corrections too
! large to apply.
subroutine log_henyey_corrections()
! The correction trace is written only with solver diagnostics on.
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
! Apply the corrections, scaled by the correction factor f, at every
! shell: log P += f dlog P and log T += f dlog T (and the same to their
! changes over the step, log_pressure_delta and log_temperature_delta,
! which the entropy term reads), log R += f dlog R, L += f dL. Then the
! new log L, and log Teff from the envelope fit,
!   log Teff = cf(7) log P(nz) + cf(8) log T(nz) + cf(9).
subroutine apply_henyey_corrections()
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
! With rotation, the rotation curve, the P and T shape factors and the
! rotational kinetic energy of the corrected structure,
! E_rot = 0.5 omega j dm per shell.
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

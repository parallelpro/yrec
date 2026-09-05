!----------------------------------------------------------------------
! henyey_iterate (formerly crrect)
!----------------------------------------------------------------------
! YREC's Newton-Raphson relaxation driver: the Henyey structure solve
! for one time step. Before the correction loop it refreshes the
! surface boundary condition (surfbc; only when recompute_surface_bc
! asks for it) and, on the iteration levels that need it, mixes the
! convection zones (mix). Each correction iteration then calls
! henyey_coefficients to linearize the structure equations at every
! mesh point, builds the outer boundary-condition row from the
! envelope-fit coefficients (star%envelope_fit_coeffs, from surfbc),
! calls henyey_solve to back-solve the block-tridiagonal system for
! the P/T/R/L corrections, checks those corrections against the
! convergence/divergence tolerances (ctrl%*), applies them, and (if
! rotation is active) updates the rotation curve and the rotational
! P/T factors via omega_from_j/rotation_shape_factors.
!
! Non-zero ierr propagates a callee failure to evolve_step.
!
! INPUTS ASSUMES GIVEN LOG(TE) AS star%log_Teff
!        ASSUMES GIVEN LOG(L/LSUN) AS star%log_L
!        delta_time = TIME STEP IN SECONDS
!        max_iterations = ITERATION NUMBER
! OUTPUTS  converged = .T. IF MODEL HAS CONVERGED
!          corrections_too_large = .T. IF CORRECTIONS ARE TOO LARGE
subroutine henyey_iterate(delta_time, max_iterations, converged, &
     corrections_too_large, start_new_triangle, reset_triangle, &
     recompute_surface_bc, tri_orientation, stored_vertex_index, &
     dlnrho_dlnt, dlnrho_dlnp, iterations_done, iteration_level, ierr)
      use star_info_lib, only: star, i_c12, i_c13, i_h1, i_he3, i_he4, i_lum_3alpha, i_lum_cno, i_lum_grav, i_lum_he_c, i_lum_neu, i_lum_pp1, i_lum_pp2, i_lum_pp3, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18, &
           n_species_basic, n_species_extended, n_lum_channels, max_convective_zones
      use luout_lib
      use run_log_lib, only: solver_diagnostics
      use phys_const_lib
      use yale_eos_lib

      use math_lib
      implicit none

      double precision, intent(in) :: delta_time
      integer, intent(in) :: max_iterations
      logical, intent(inout) :: converged
      logical, intent(out) :: corrections_too_large
      logical, intent(inout) :: start_new_triangle, reset_triangle
      logical, intent(inout) :: recompute_surface_bc
      double precision, intent(inout) :: tri_orientation
      integer, intent(inout) :: stored_vertex_index
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp
      integer, intent(inout) :: iterations_done
      integer, intent(in) :: iteration_level

! --- locals ---
      integer :: kenv, katm, ksaha
      integer :: num_species, i, j, iter, max_correction_pos
      integer :: core_cz_edge, envelope_zone_index
      integer :: mixed_zone_bounds_no_overshoot(max_convective_zones,2)
      double precision :: hydrogen_fraction, metal_fraction, &
           log10_pressure_limit
      double precision :: surface_pressure_rotation_factor, &
           surface_temperature_rotation_factor
      logical :: envelope_recomputed_flag
      double precision :: timestep_years
      double precision :: temp, test, total_luminosity_terms, &
           luminosity_correction_max, correction_factor
      double precision :: hydrogen_burn_luminosity, helium_burn_luminosity
      double precision :: shell_angular_momentum
      ! 2026 (ROADMAP.md stage 3): library errors return here via ierr;
      ! this driver-side call site preserves the historical stop.
      integer :: jerr

      integer, intent(out) :: ierr

      ierr = 0

      if (max_iterations.le.0) return
      star%log_L = log10(star%luminosity_lsun(star%nz))
! ZERO COUNTERS
      kenv = 0
      katm = 0
      ksaha = 0
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
       if (use_debye_huckel_correction) then
          debye_huckel_x = star%xa(i_h1,star%nz)
          debye_huckel_y = star%xa(i_he4,star%nz)+star%xa(i_he3,star%nz)
          debye_huckel_z_total = star%xa(i_metals,star%nz)
          debye_huckel_z(1) = star%xa(i_c12,star%nz)+star%xa(i_c13,star%nz)
          debye_huckel_z(2) = star%xa(i_n14,star%nz)+star%xa(i_n15,star%nz)
          debye_huckel_z(3) = star%xa(i_o16,star%nz)+star%xa(i_o17,star%nz)+ &
               star%xa(i_o18,star%nz)
       end if
       call surfbc(star%trial_log_temperature,star%trial_log_luminosity,star%envelope_fit_coeffs,star%fit_point_pressure,star%fit_point_temperature, &
            star%fit_point_radius,tri_orientation,stored_vertex_index, &
            star%stored_envelope_state, &
            start_new_triangle,reset_triangle,ksaha,kenv,katm, &
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
! 7/91 ADD CALL TO MIX
      if (iteration_level.gt.2 .and. delta_time.gt.0.0d0) then
         num_species = n_species_basic
         if (star%job%use_extended_composition) num_species = n_species_extended
         do i = 1,star%nz
            do j = 1,num_species
               star%xa(j,i) = star%xa_start(j,i)
            end do
         end do
         call mix(delta_time, iteration_level, timestep_years, &
              core_cz_edge, envelope_zone_index, &
              mixed_zone_bounds_no_overshoot, jerr)
         if (jerr /= 0) then
         ! 2026 (phase five, step B): propagate instead of stopping
            ierr = jerr
            return
         end if
      endif
!  IF MODEL CONVERGED ON THE PREVIOUS LEVEL OF ITERATION AND
!  THE TRIANGLE WAS CHECKED AND NOT FLIPPED,SKIP CORRECTION ROUTINE.
      if (recompute_surface_bc .and. .not.envelope_recomputed_flag .and. &
           converged) return
      converged = .false.
      do iter = 1,max_iterations
! DO HENYEY REDUCTION
       star%max_residual(1) = 0.0d0
       star%max_residual(2) = 0.0d0
       star%max_residual(3) = 0.0d0
       star%max_residual(4) = 0.0d0
! (2026 de-tramp: the star%-shaped slots are read/written by
! henyey_coefficients directly; only the per-call locals remain.)
       call henyey_coefficients(delta_time,dlnrho_dlnt,dlnrho_dlnp,ksaha, &
            envelope_zone_index, jerr)
       if (jerr /= 0) then
       ! 2026 (phase five, step B): propagate instead of stopping
          ierr = jerr
          return
       end if
! RENORMALIZE TLUMX-S
!CC   TAKE OUT RENORMALIZATION DURING HE FLASH (NON-THERMAL EQUALIBRIUM)
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
! DO BACK SOLUTION FOR CORRECTIONS
       call henyey_solve(star%nz,star%elim_coeff,star%elim_rhs,star%surface_bc)
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
       if (star%ctrl%fcorr0.gt.0.0d0) star%job%fcorr = dmin1(1.d0,star%job%fcorr+star%ctrl%fcorri)
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
       iterations_done = iterations_done + 1
       if (converged) return
      end do

      return
end subroutine henyey_iterate

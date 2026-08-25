!----------------------------------------------------------------------
! crrect
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original crrect.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! JANUARY 28, 1991:
! REPLACED EGEN SO THAT IT FOLLOWS FULL CNO NOT CN IN EQUILIBRIUM AND
! FEBRUARY 11, 1991:
! REPLACED CONVEC SO THAT WOULD CALCULATE OVERSHOOT CORRECTLY.
! ADDED ROUTINES OVERSH AND HSUBP.
!
! This is YREC's Newton-Raphson relaxation driver: it orchestrates the
! Henyey structure solve for one time step. Each pass through the
! main DO 100 loop calls coefft to linearize the structure equations
! at every mesh point, builds the outer (surface) boundary-condition
! row from the envelope-fit coefficients (star%envelope_fit_coeffs, from
! surfbc), calls hsolve to back-solve the block-tridiagonal system for
! the P/T/R/L corrections, checks those corrections against the
! convergence/divergence tolerances in common/ctol/, applies them, and
! (if rotation is active) updates the rotation curve and rotational
! kinetic energy via getrot/fpft. surfbc (envelope/atmosphere fit) and
! mix (convective mixing/star%xa update) are called once near the
! top, ahead of the correction loop, on the iteration levels where
! they are needed.
!
! CROSS-CALLEE NAMING NOTE: several dummy arguments here are threaded
! into more than one already-converted callee, and those callees do
! not always agree on a name for the same physical data (this file's
! own dummy names were free to choose, per the project's incremental
! conversion order, but the callees' names were fixed by earlier
! batches). Judgment calls made below, all verified against the
! actual physics/usage in this file:
!   - BL is log10(L/Lsun) here (see BL = DLOG10(HL(M)) below), named
!     star%log_L. It is passed into surfbc.f90's parameter
!     named "luminosity_linear" -- that surfbc.f90 name is a misnomer
!     inherited from that file's own earlier conversion; out of scope
!     to fix here.
!   - HL is the linear luminosity (L/Lsun); DLOG10(HL(M)) above proves
!     it is not itself a log quantity. Named star%luminosity_lsun (matches
!     coefft.f90's slot name for it). mix.f90's slot name for this
!     same array is "log_luminosity", which is the same kind of
!     misnomer as the BL case above; also out of scope to fix here.
!   - HSTOT (total stellar mass, log10(M/Msun)) is named
!     star%log_total_mass (matches mix.f90's slot name); surfbc.f90's slot
!     name for the same value is "log10_star_mass".
!   - TEFFL (log10 Teff) is named star%log_Teff (matches mix.f90/
!     coefft.f90); surfbc.f90's slot name for the same value is
!     "log10_teff".
!   - DELTS (time step, seconds) is named delta_time (matches
!     coefft.f90); mix.f90's slot name for the same value is
!     "timestep".
!   - M (number of mesh points) is named star%nz (matches
!     coefft.f90/fpft.f90); mix.f90/getrot.f90 call the same count
!     "num_zones", hsolve.f90 calls it "num_shells", surfbc.f90 calls
!     it "zone_index" (there it is the single index M, not a count).
!   - HS1 is named star%m; coefft.f90's slot name for the same
!     array is "mass_weight_ln". Not read directly in this file
!     (only passed through to mix/coefft), so the more physically
!     transparent of the two established names was kept.
!   - FP/FT (rotational P/T correction factors) are named
!     star%fp_rot/star%ft_rot (matches
!     fpft.f90, which computes them); coefft.f90's slot names for the
!     same arrays are rotation_p_factor/rotation_t_factor.
!   - R0 is named star%mean_radius (matches getrot.f90, which computes it);
!     fpft.f90's slot name for the same array is "r0".
!   - ETA2 is named star%eta_squared (matches getrot.f90); fpft.f90's slot
!     name for the same array is "eta2".
!   - The Debye-Huckel hydrogen-fraction common/debhu/ member (XXDH in
!     the original) is named debye_huckel_x here to match coefft.f90 (this
!     file's closest sibling, called every iteration below); note
!     hsubp.f90 instead uses the physically-clearer "xxdh" for the
!     same slot -- see hsubp.f90's own comment on this pre-existing
!     inconsistency.
!
! INPUTS ASSUMES GIVEN LOG(TE) AS star%log_Teff
!        ASSUMES GIVEN LOG(L/LSUN) AS star%log_L
!        delta_time = TIME STEP IN SECONDS
!        max_iterations = ITERATION NUMBER
! OUTPUTS  converged = .T. IF MODEL HAS CONVERGED
!          corrections_too_large = .T. IF CORRECTIONS ARE TOO LARGE
subroutine crrect(delta_time, max_iterations, converged, &
     corrections_too_large, start_new_triangle, reset_triangle, &
     recompute_surface_bc, tri_orientation, stored_vertex_index, &
     in_atmosphere, want_derivatives, mixing_active, &
     conductive_opacity_flag, dlnrho_dlnt, dlnrho_dlnp, iterations_done, &
     iteration_level, ierr)
      use star_info_lib, only: star, i_c12, i_c13, i_h1, i_he3, i_he4, i_lum_3alpha, i_lum_cno, i_lum_grav, i_lum_he_c, i_lum_neu, i_lum_pp1, i_lum_pp2, i_lum_pp3, i_metals, i_n14, i_n15, i_o16, i_o17, i_o18, json
      use luout_lib
      use const_lib
      use yale_eos_lib

      implicit none

      double precision, intent(in) :: delta_time
      integer, intent(in) :: max_iterations
      logical, intent(inout) :: converged
      logical, intent(out) :: corrections_too_large
      logical, intent(inout) :: start_new_triangle, reset_triangle
      logical, intent(inout) :: recompute_surface_bc
      double precision, intent(inout) :: tri_orientation
      integer, intent(inout) :: stored_vertex_index
      logical, intent(out) :: in_atmosphere, want_derivatives, &
           mixing_active, conductive_opacity_flag
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp
      integer, intent(inout) :: iterations_done
      integer, intent(in) :: iteration_level

















! --- locals ---
      integer :: kenv, katm, ksaha
      integer :: num_species, i, j, iter, max_correction_pos
      integer :: core_cz_edge, envelope_zone_index
      integer :: mixed_zone_bounds_no_overshoot(12,2)
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
      star%log_L = dlog10(star%luminosity_lsun(star%nz))
! ZERO COUNTERS
      kenv = 0
      katm = 0
      ksaha = 0
      star%env_comp%senv = star%log_mass(star%nz) - star%log_total_mass
      if (start_new_triangle.or.reset_triangle .and.iteration_level.eq.2) &
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
       hydrogen_fraction = star%env_comp%xnew
       metal_fraction = star%env_comp%znew
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
!               LNEW,LRESET,LSBC,KSAHA,KENV,KATM,HSTOT,BL,  ! KC 2025-05-31
!               (recompute_surface_bc/LSBC removed from this call site)
            start_new_triangle,reset_triangle,ksaha,kenv,katm, &
            star%log_total_mass,star%log_L, &
            star%log_Teff,hydrogen_fraction,metal_fraction, &
            surface_pressure_rotation_factor, &
            surface_temperature_rotation_factor,envelope_recomputed_flag, &
            log10_pressure_limit,star%convective_flag,star%nz)
      endif
! 7/91 ADD CALL TO MIX
      if (iteration_level.gt.2 .and. delta_time.gt.0.0d0) then
         num_species = 11
         if (star%job%use_extended_composition) num_species = 15
         do i = 1,star%nz
            do j = 1,num_species
               star%xa(j,i) = star%prev%xa_start(j,i)
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
!      IF(LROT)THEN
!         CALL OVROT(HCOMP,HD,HP,HR,HS,HT,LC,M,LCZ,MRZONE,MXZONE,
!     *        NRZONE,NZONE)
!      ENDIF
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
       call coefft(delta_time,star%nz,star%logRho,star%elim_coeff,star%elim_rhs, &
            star%gravitational_luminosity,star%luminosity_lsun,star%max_residual, &
            star%logP,star%log_pressure_delta,star%logR,star%log_mass, &
            star%m,star%dm,star%logT,star%log_temperature_delta, &
            star%xa,star%convective_flag,star%luminosity_breakdown,in_atmosphere, &
            want_derivatives,mixing_active,conductive_opacity_flag, &
            dlnrho_dlnt,dlnrho_dlnp, &
!      *   KSAHA,MODEL,FP,FT,HKEROT,HKEROT0,JENV,TEFFL)  ! KC 2025-05-31
            ksaha,star%fp_rot,star%ft_rot, &
            star%kinetic_energy_rot,star%kinetic_energy_rot_old,envelope_zone_index, &
            star%log_Teff, jerr)
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
          do j = 1,8
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
             write(short_file_unit,20) (star%max_residual(j),j=1,4)
   20          format(' R.H.S. BELOW TOLERANCES--P',1PE9.2,'  T ',E9.2, &
             '  R ',E9.2,'  L ',E9.2)
             write(short_file_unit,75) iterations_done+1
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
       call hsolve(star%nz,star%elim_coeff,star%elim_rhs,star%surface_bc)
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
       if (star%ctrl%lcorr) then
          write (short_file_unit,60) converged,star%max_residual(4), &
               star%ctrl%htoler(4,1),star%max_correction_index(4)
   60       format (1X,'DEL-L/L  ',L2,1P2E12.4,5X,I5)
          write(short_file_unit,70)(star%max_correction_index(j),star%max_residual(j), &
               j=1,4),correction_factor,hydrogen_burn_luminosity, &
               helium_burn_luminosity,(star%luminosity_breakdown(j),j=6,7)
   70       format(' CORR',I5,'P',1PE9.2,I5,'T',E9.2,I5,'R',E9.2,I5,'L', &
            E9.2,'  F=',0PF5.3,'  E-HY',1PE10.3,' HE',E10.3,' NU',E10.3, &
            ' G',E10.3)
          if (converged) then
             write(short_file_unit,75) iterations_done+1
   75          format(10X,'MODEL CONVERGED AFTER ',I4,'  ITERATIONS')
          endif
       endif
       if (corrections_too_large) then
          write(iowr,80)
          write(short_file_unit,80)
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
       star%log_L = dlog10(star%luminosity_lsun(star%nz))
       star%log_Teff = star%envelope_fit_coeffs(7)*star%logP(star%nz) + &
            star%envelope_fit_coeffs(8)*star%logT(star%nz) + &
            star%envelope_fit_coeffs(9)
         if (star%job%rotation_active) then
            call getrot(star%logRho,star%j_rot,star%logR, &
                 star%log_mass,star%dm,star%am_transport_convective_flag, &
                 star%nz,star%eta_squared,star%i_rot,star%omega,star%qiw, &
                 star%mean_radius)
            call fpft(star%logRho,star%logR,star%log_mass,star%nz,star%omega, &
                 star%eta_squared,star%fp_rot, &
                 star%ft_rot,star%mean_gravity,star%mean_radius)
            do i = 1,star%nz
               shell_angular_momentum = star%j_rot(i)*star%dm(i)
               star%kinetic_energy_rot(i) = 0.5d0*star%omega(i)*shell_angular_momentum
            end do
         endif
       iterations_done = iterations_done + 1
       if (converged) return
      end do

      return
end subroutine crrect

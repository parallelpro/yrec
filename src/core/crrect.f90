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
! row from the envelope-fit coefficients (envelope_coeffs, from
! surfbc), calls hsolve to back-solve the block-tridiagonal system for
! the P/T/R/L corrections, checks those corrections against the
! convergence/divergence tolerances in common/ctol/, applies them, and
! (if rotation is active) updates the rotation curve and rotational
! kinetic energy via getrot/fpft. surfbc (envelope/atmosphere fit) and
! mix (convective mixing/composition update) are called once near the
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
!     log10_luminosity. It is passed into surfbc.f90's parameter
!     named "luminosity_linear" -- that surfbc.f90 name is a misnomer
!     inherited from that file's own earlier conversion; out of scope
!     to fix here.
!   - HL is the linear luminosity (L/Lsun); DLOG10(HL(M)) above proves
!     it is not itself a log quantity. Named luminosity_lsun (matches
!     coefft.f90's slot name for it). mix.f90's slot name for this
!     same array is "log_luminosity", which is the same kind of
!     misnomer as the BL case above; also out of scope to fix here.
!   - HSTOT (total stellar mass, log10(M/Msun)) is named
!     log_total_mass (matches mix.f90's slot name); surfbc.f90's slot
!     name for the same value is "log10_star_mass".
!   - TEFFL (log10 Teff) is named log_teff (matches mix.f90/
!     coefft.f90); surfbc.f90's slot name for the same value is
!     "log10_teff".
!   - DELTS (time step, seconds) is named delta_time (matches
!     coefft.f90); mix.f90's slot name for the same value is
!     "timestep".
!   - M (number of mesh points) is named num_points (matches
!     coefft.f90/fpft.f90); mix.f90/getrot.f90 call the same count
!     "num_zones", hsolve.f90 calls it "num_shells", surfbc.f90 calls
!     it "zone_index" (there it is the single index M, not a count).
!   - HS1 is named enclosed_mass; coefft.f90's slot name for the same
!     array is "mass_weight_ln". Not read directly in this file
!     (only passed through to mix/coefft), so the more physically
!     transparent of the two established names was kept.
!   - FP/FT (rotational P/T correction factors) are named
!     pressure_rotation_factor/temperature_rotation_factor (matches
!     fpft.f90, which computes them); coefft.f90's slot names for the
!     same arrays are rotation_p_factor/rotation_t_factor.
!   - R0 is named mean_radius (matches getrot.f90, which computes it);
!     fpft.f90's slot name for the same array is "r0".
!   - ETA2 is named eta_squared (matches getrot.f90); fpft.f90's slot
!     name for the same array is "eta2".
!   - The Debye-Huckel hydrogen-fraction common/debhu/ member (XXDH in
!     the original) is named xxdy here to match coefft.f90 (this
!     file's closest sibling, called every iteration below); note
!     hsubp.f90 instead uses the physically-clearer "xxdh" for the
!     same slot -- see hsubp.f90's own comment on this pre-existing
!     inconsistency.
!
! INPUTS ASSUMES GIVEN LOG(TE) AS log_teff
!        ASSUMES GIVEN LOG(L/LSUN) AS log10_luminosity
!        delta_time = TIME STEP IN SECONDS
!        max_iterations = ITERATION NUMBER
! OUTPUTS  converged = .T. IF MODEL HAS CONVERGED
!          corrections_too_large = .T. IF CORRECTIONS ARE TOO LARGE
subroutine crrect(delta_time, num_points, max_iterations, converged, &
     corrections_too_large, start_new_triangle, reset_triangle, &
     recompute_surface_bc, &
     tri_teffl, tri_logl, envelope_coeffs, vtx_logp, vtx_logt, vtx_logr, &
     tri_orientation, stored_vertex_index, stored_envelope_state, &
     log_total_mass, log10_luminosity, log_teff, &
     log_density, elim_coeff, elim_rhs, gravitational_luminosity, &
     luminosity_lsun, max_residual, log_pressure, log_pressure_delta, &
     log_radius, log_mass, enclosed_mass, shell_mass, surface_bc, &
     log_temperature, log_temperature_delta, composition, convective_flag, &
     max_correction_index, luminosity_terms, in_atmosphere, want_derivatives, &
     mixing_active, conductive_opacity_flag, dlnrho_dlnt, dlnrho_dlnp, &
     pressure_rotation_factor, temperature_rotation_factor, eta_squared, &
     omega, mean_radius, iterations_done, mean_gravity, moment_of_inertia, &
     specific_angular_momentum, iteration_level, am_transport_convective_flag, &
     mixed_zone_bounds, qiw, kinetic_energy_rot, kinetic_energy_rot_old)

      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: delta_time
      integer, intent(in) :: num_points
      integer, intent(in) :: max_iterations
      logical, intent(inout) :: converged
      logical, intent(out) :: corrections_too_large
      logical, intent(inout) :: start_new_triangle, reset_triangle
      logical, intent(inout) :: recompute_surface_bc
      double precision, intent(inout) :: tri_teffl(3), tri_logl(3), &
           envelope_coeffs(9)
      double precision, intent(inout) :: vtx_logp(3), vtx_logt(3), &
           vtx_logr(3)
      double precision, intent(inout) :: tri_orientation
      integer, intent(inout) :: stored_vertex_index
      double precision, intent(inout) :: stored_envelope_state(4)
      double precision, intent(in) :: log_total_mass
      double precision, intent(inout) :: log10_luminosity, log_teff
      double precision, intent(inout) :: log_density(json)
      double precision, intent(inout) :: elim_coeff(4,2,json), &
           elim_rhs(4,json)
      double precision, intent(inout) :: gravitational_luminosity(json)
      double precision, intent(inout) :: luminosity_lsun(json)
      double precision, intent(inout) :: max_residual(4)
      double precision, intent(inout) :: log_pressure(json), &
           log_pressure_delta(json), log_radius(json)
      double precision, intent(in) :: log_mass(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(in) :: shell_mass(json)
      double precision, intent(inout) :: surface_bc(6)
      double precision, intent(inout) :: log_temperature(json), &
           log_temperature_delta(json)
      double precision, intent(inout) :: composition(15,json)
      logical, intent(inout) :: convective_flag(json)
      integer, intent(out) :: max_correction_index(4)
      double precision, intent(inout) :: luminosity_terms(8)
      logical, intent(out) :: in_atmosphere, want_derivatives, &
           mixing_active, conductive_opacity_flag
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp
      double precision, intent(inout) :: pressure_rotation_factor(json), &
           temperature_rotation_factor(json)
      double precision, intent(inout) :: eta_squared(json), omega(json), &
           mean_radius(json)
      integer, intent(inout) :: iterations_done
      double precision, intent(out) :: mean_gravity(json)
      double precision, intent(inout) :: moment_of_inertia(json), &
           specific_angular_momentum(json)
      integer, intent(in) :: iteration_level
      logical, intent(in) :: am_transport_convective_flag(json)
      integer, intent(inout) :: mixed_zone_bounds(12,2)
      double precision, intent(inout) :: qiw(json)
      double precision, intent(inout) :: kinetic_energy_rot(json)
      double precision, intent(in) :: kinetic_energy_rot_old(json)

! common/luout/: short_file_unit/iowr are used here. Naming matches
! getopac.f90.
      integer :: ilast, idebug, itrack, short_file_unit, imilne, imodpt, &
           istor, iowr
      common/luout/ ilast, idebug, itrack, short_file_unit, imilne, &
           imodpt, istor, iowr

! common/ccout/: not used in this file. Naming matches coefft.f90.
      logical :: lstore, lstatm, lstenv, lstmod, lstphys, lstrot, lscrib, &
           lstch, lphhd
      common/ccout/ lstore, lstatm, lstenv, lstmod, lstphys, lstrot, &
           lscrib, lstch, lphhd

! common/ccout2/: lcorr is used here (LDEBUG/LMILNE/LTRACK/LSTPCH are
! unused placeholders). Naming matches hpoint.f90/meqos.f90.
      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common/ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch

! common/cenv/: not used in this file. Naming matches surfbc.f90.
      double precision :: tri_delta_teffl, tri_delta_logl, senv0_placeholder
      logical :: lsenv0_placeholder, lnew0_placeholder
      common/cenv/ tri_delta_teffl, tri_delta_logl, senv0_placeholder, &
           lsenv0_placeholder, lnew0_placeholder

! common/comp/: xnew/znew/senv are used here (the rest are unused
! placeholders). Naming matches hpoint.f90/getopac.f90.
      double precision :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, zenvm, amuenv, fxenv(12), xnew, znew, &
           stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv, xnew, znew, stotal, senv

! common/const/: not used in this file. Naming matches wrtout.f90.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/ solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude

! common/const1/: only ln10 (originally CLN) is used here, in the
! surface-BC luminosity-equation coefficients. Naming matches
! coefft.f90/eqburn.f90.
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
      common/const1/ ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! common/const2/: not used in this file. Naming matches coefft.f90.
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn
      common/const2/ gas_constant, radiation_constant_over_3, ca3l, csig, &
           csigl, cgl, cmkh, cmkhn

! common/ctol/: htoler/fcorr0/fcorri/fcorr are all used here for the
! Newton-Raphson convergence/divergence tolerances and correction
! damping. Naming matches hpoint.f90/mixgrid.f90.
      double precision :: htoler(5,2), fcorr0, fcorri, fcorr, &
           chi_grid_scale(12)
      integer :: niter1, niter2, niter3
      common/ctol/ htoler, fcorr0, fcorri, fcorr, chi_grid_scale, niter1, &
           niter2, niter3

! common/envgen/: not used in this file. Naming matches stitch.f90/
! wrtmod.f90.
      double precision :: atm_step_size, envelope_step_size
      logical :: envelope_generation_flag
      common/envgen/ atm_step_size, envelope_step_size, &
           envelope_generation_flag

! common/heflsh/: helium_flash_active (originally LKUTHE) is used here
! to gate the He-flash correction-damping logic. Naming matches
! wrtlst.f90/wrtout.f90.
      logical :: helium_flash_active
      common/heflsh/ helium_flash_active

! common/rot/: rotation_active (originally LROT) is used here to gate
! the rotational P/T factor and rotation-curve updates. Naming matches
! coefft.f90/fpft.f90.
      double precision :: wnew, walpcz, acfpft
      integer :: itfp1, itfp2
      logical :: rotation_active, instability_transport_active, lwnew
      common/rot/ wnew, walpcz, acfpft, itfp1, itfp2, rotation_active, &
           instability_transport_active, lwnew

! common/neweng/: not used in this file. Naming matches coefft.f90/
! mix.f90.
      integer :: niter4
      logical :: lnews, lsnu
      common/neweng/ niter4, lnews, lsnu

! common/oldmod/: old_composition (originally HCOMPP) is used here to
! restore the pre-Henyey composition ahead of calling mix. Naming
! matches mix.f90/hpoint.f90/dburn.f90.
      double precision :: old_pressure(json), old_temperature(json), &
           old_radius(json), old_luminosity(json), old_density(json), &
           old_composition(15,json), old_shell_mass(json), old_teff
      logical :: old_convective_flag(json), old_cz_flag(json)
      integer :: old_num_zones
      common/oldmod/ old_pressure, old_temperature, old_radius, &
           old_luminosity, old_density, old_composition, old_shell_mass, &
           old_convective_flag, old_cz_flag, old_teff, old_num_zones

! common/flag/: use_extended_composition (originally LEXCOM) is used
! here to pick the number of mixed species. Naming matches mix.f90/
! mixcz.f90.
      logical :: use_extended_composition
      common/flag/ use_extended_composition

! DBG 7/92 COMMON BLOCK ADDED TO COMPUTE DEBYE-HUCKEL CORRECTION.
! common/debhu/: cdh/ldh gate and select the correction; xxdy (H
! fraction)/yydh (He)/zzdh (metals)/zdh(1:3) (CNO groups) are set here
! ahead of calling surfbc. Naming matches coefft.f90 (xxdy there is an
! apparent pre-existing typo for the hydrogen fraction -- see the
! cross-callee naming note above).
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdy, yydh, zzdh, &
           dhnue(18)
      logical :: ldh
      common/debhu/ cdh, etadh0, etadh1, zdh, xxdy, yydh, zzdh, dhnue, ldh

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

      save

      if (max_iterations.le.0) return
      log10_luminosity = dlog10(luminosity_lsun(num_points))
! ZERO COUNTERS
      kenv = 0
      katm = 0
      ksaha = 0
      senv = log_mass(num_points) - log_total_mass
      if (start_new_triangle.or.reset_triangle .and.iteration_level.eq.2) &
           recompute_surface_bc = .true.
!  FIND NEW FP AND FT IF MODEL IS ROTATING
      if (rotation_active.and.recompute_surface_bc) then
       surface_pressure_rotation_factor = pressure_rotation_factor(num_points)
       surface_temperature_rotation_factor = &
            temperature_rotation_factor(num_points)
      else
       surface_pressure_rotation_factor = 1.0d0
       surface_temperature_rotation_factor = 1.0d0
      endif
! SET UP SURFACE BOUNDARY CONDITIONS-2ND AND 3RD LEVELS OF ITER ONLY
! FIND ENVELOPE MASS AND SET X AND Z TO ENVELOPE VALUES
      if (recompute_surface_bc) then
       hydrogen_fraction = xnew
       metal_fraction = znew
       log10_pressure_limit = log_pressure(num_points)
       if (ldh) then
          xxdy = composition(1,num_points)
          yydh = composition(2,num_points)+composition(4,num_points)
          zzdh = composition(3,num_points)
          zdh(1) = composition(5,num_points)+composition(6,num_points)
          zdh(2) = composition(7,num_points)+composition(8,num_points)
          zdh(3) = composition(9,num_points)+composition(10,num_points)+ &
               composition(11,num_points)
       end if
       call surfbc(tri_teffl,tri_logl,envelope_coeffs,vtx_logp,vtx_logt, &
            vtx_logr,tri_orientation,stored_vertex_index, &
            stored_envelope_state, &
!               LNEW,LRESET,LSBC,KSAHA,KENV,KATM,HSTOT,BL,  ! KC 2025-05-31
!               (recompute_surface_bc/LSBC removed from this call site)
            start_new_triangle,reset_triangle,ksaha,kenv,katm, &
            log_total_mass,log10_luminosity, &
            log_teff,hydrogen_fraction,metal_fraction, &
            surface_pressure_rotation_factor, &
            surface_temperature_rotation_factor,envelope_recomputed_flag, &
            log10_pressure_limit,convective_flag,num_points)
      endif
! 7/91 ADD CALL TO MIX
      if (iteration_level.gt.2 .and. delta_time.gt.0.0d0) then
         num_species = 11
         if (use_extended_composition) num_species = 15
         do 2 i = 1,num_points
            do 1 j = 1,num_species
               composition(j,i) = old_composition(j,i)
    1       continue
    2    continue
         call mix(delta_time,composition,log_density,luminosity_lsun, &
              log_pressure,log_radius,log_mass,enclosed_mass,shell_mass, &
              log_total_mass, &
!        *        HT,ITLVL,LC,M,QDT,QDP,DDAGE,JCORE,JENV,  ! KC 2025-05-31
              log_temperature,iteration_level,convective_flag,num_points, &
              timestep_years,core_cz_edge,envelope_zone_index, &
              mixed_zone_bounds,mixed_zone_bounds_no_overshoot,log_teff)
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
      do 100 iter = 1,max_iterations
! DO HENYEY REDUCTION
       max_residual(1) = 0.0d0
       max_residual(2) = 0.0d0
       max_residual(3) = 0.0d0
       max_residual(4) = 0.0d0
       call coefft(delta_time,num_points,log_density,elim_coeff,elim_rhs, &
            gravitational_luminosity,luminosity_lsun,max_residual, &
            log_pressure,log_pressure_delta,log_radius,log_mass, &
            enclosed_mass,shell_mass,log_temperature,log_temperature_delta, &
            composition,convective_flag,luminosity_terms,in_atmosphere, &
            want_derivatives,mixing_active,conductive_opacity_flag, &
            dlnrho_dlnt,dlnrho_dlnp, &
!      *   KSAHA,MODEL,FP,FT,HKEROT,HKEROT0,JENV,TEFFL)  ! KC 2025-05-31
            ksaha,pressure_rotation_factor,temperature_rotation_factor, &
            kinetic_energy_rot,kinetic_energy_rot_old,envelope_zone_index, &
            log_teff)
! RENORMALIZE TLUMX-S
!CC   TAKE OUT RENORMALIZATION DURING HE FLASH (NON-THERMAL EQUALIBRIUM)
       total_luminosity_terms = luminosity_terms(1)+luminosity_terms(2)+ &
            luminosity_terms(3)+luminosity_terms(4)+luminosity_terms(5)+ &
            luminosity_terms(6)+luminosity_terms(7)+luminosity_terms(8)
       if (.not.helium_flash_active .and. total_luminosity_terms.gt.0.0d0) &
            then
          temp = luminosity_lsun(num_points)/total_luminosity_terms
          do 10 j = 1,8
             luminosity_terms(j) = luminosity_terms(j)*temp
   10       continue
       endif
! CHECK ON SIGNIFICANCE OF R.H.S. EQUATIONS FOR P AND T
! N.B.  DOES NOT CHECK DIFFERENCES IN BOUNDARY EQUATIONS
       if (iter.gt.1) then
          if (max_residual(1).le.htoler(5,1) .and. &
               max_residual(2).le.htoler(5,2) .and. &
               max_residual(3).le.htoler(5,1) .and. &
               max_residual(4).le. htoler(5,2)) then
             write(short_file_unit,20) (max_residual(j),j=1,4)
   20          format(' R.H.S. BELOW TOLERANCES--P',1PE9.2,'  T ',E9.2, &
             '  R ',E9.2,'  L ',E9.2)
             write(short_file_unit,75) iterations_done+1
             converged = .true.
             return
          endif
       endif
! SET UP HENYEY COEFFICIENTS FOR SURFACE
       surface_bc(1) = -envelope_coeffs(1)
       surface_bc(2) = -envelope_coeffs(2)
       surface_bc(3) = envelope_coeffs(1)*log_pressure(num_points) + &
            envelope_coeffs(2)*log_temperature(num_points) + &
            envelope_coeffs(3) - log_radius(num_points)
       temp = ln10*luminosity_lsun(num_points)
       surface_bc(4) =-temp*envelope_coeffs(4)
       surface_bc(5) =-temp*envelope_coeffs(5)
       surface_bc(6) = temp*(envelope_coeffs(4)*log_pressure(num_points)+ &
            envelope_coeffs(5)*log_temperature(num_points)+ &
            envelope_coeffs(6)-log10_luminosity)
! DO BACK SOLUTION FOR CORRECTIONS
       call hsolve(num_points,elim_coeff,elim_rhs,surface_bc)
! CHECK ON MAXIMUM CORRECTIONS
       do 30 j = 1,4
          max_residual(j) = dabs(elim_rhs(j,1))
          max_correction_index(j) = 1
   30    continue
       do 40 i = 2,num_points
          test = dabs(elim_rhs(1,i))
          if (max_residual(1).le.test) then
             max_residual(1) = test
             max_correction_index(1) = i
          endif
          test = dabs(elim_rhs(2,i))
          if (max_residual(2).le.test) then
             max_residual(2) = test
             max_correction_index(2) = i
          endif
          test = dabs(elim_rhs(3,i))
          if (max_residual(3).le.test) then
             max_residual(3) = test
             max_correction_index(3) = i
          endif
          test = dmin1(dabs(elim_rhs(4,i)),dabs(elim_rhs(4,i)/luminosity_lsun(i)))
          if (max_residual(4).le.test) then
             max_residual(4) = test
             max_correction_index(4) = i
          endif
   40    continue
!CC   HE FLASH -- OK FOR ALL
       luminosity_correction_max = max_residual(4)
! LFINI = T IF MAX CORRECTIONS LESS THAN CONVERGENCE CRITERIA SET IN
! HTOLER. LARGE = T IF MAX CORRECTIONS GREATER THAN LARGEST CORRECTIONS
! ALLOWED, ALSO SET IN HTOLER
       converged = max_residual(1).lt.htoler(1,1) .and. &
            max_residual(2).lt.htoler(2,1) &
         .and. max_residual(3).lt.htoler(3,1) .and. &
         max_residual(4).lt.htoler(4,1)
       corrections_too_large = max_residual(1).gt.htoler(1,2) .or. &
            max_residual(2).gt.htoler(2,2) &
         .or. max_residual(3).gt.htoler(3,2) .or. &
         max_residual(4).gt.htoler(4,2)
       do 50 j = 1,4
          max_correction_pos = max_correction_index(j)
          max_residual(j) = elim_rhs(j,max_correction_pos)
   50    continue
       if (fcorr0.gt.0.0d0) fcorr = dmin1(1.d0,fcorr+fcorri)
! HE FLASH CHANGE
       correction_factor = fcorr
       if (helium_flash_active) then
          if (luminosity_correction_max.le.5.0d-1) correction_factor=8.0d-1
          if (luminosity_correction_max.le.5.0d-3) correction_factor=1.0d0
       endif
       hydrogen_burn_luminosity = luminosity_terms(1) + luminosity_terms(2) &
            + luminosity_terms(3) + luminosity_terms(4)
       helium_burn_luminosity= luminosity_terms(5) + luminosity_terms(8)
       if (lcorr) then
          write (short_file_unit,60) converged,max_residual(4), &
               htoler(4,1),max_correction_index(4)
   60       format (1X,'DEL-L/L  ',L2,1P2E12.4,5X,I5)
          write(short_file_unit,70)(max_correction_index(j),max_residual(j), &
               j=1,4),correction_factor,hydrogen_burn_luminosity, &
               helium_burn_luminosity,(luminosity_terms(j),j=6,7)
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
       do 90 i = 1,num_points
          temp = correction_factor*elim_rhs(1,i)
          log_pressure(i) = log_pressure(i) + temp
          log_pressure_delta(i) = log_pressure_delta(i) + temp
          temp = correction_factor*elim_rhs(2,i)
          log_temperature(i) = log_temperature(i) + temp
          log_temperature_delta(i) = log_temperature_delta(i) + temp
          log_radius(i) = log_radius(i) + correction_factor*elim_rhs(3,i)
          luminosity_lsun(i) = luminosity_lsun(i) + &
               correction_factor*elim_rhs(4,i)
   90    continue
       log10_luminosity = dlog10(luminosity_lsun(num_points))
       log_teff = envelope_coeffs(7)*log_pressure(num_points) + &
            envelope_coeffs(8)*log_temperature(num_points) + &
            envelope_coeffs(9)
         if (rotation_active) then
            call getrot(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,am_transport_convective_flag, &
                 num_points,eta_squared,moment_of_inertia,omega,qiw, &
                 mean_radius)
            call fpft(log_density,log_radius,log_mass,num_points,omega, &
                 eta_squared,pressure_rotation_factor, &
                 temperature_rotation_factor,mean_gravity,mean_radius)
            do i = 1,num_points
               shell_angular_momentum = specific_angular_momentum(i)*shell_mass(i)
               kinetic_energy_rot(i) = 0.5d0*omega(i)*shell_angular_momentum
            end do
         endif
       iterations_done = iterations_done + 1
       if (converged) return
  100 continue

      return
end subroutine crrect

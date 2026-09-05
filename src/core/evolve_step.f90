!----------------------------------------------------------------------
! evolve_step
!----------------------------------------------------------------------
! Added 2026 (phase five -- the embeddable engine; ROADMAP.md). One
! model advance, moved from run_yrec's model loop body: the
! divergence/retry loop (retry_step), the corrector ladder via
! henyey_iterate (solve_level), mixing, rotation
! (evolve_angular_momentum), rezoning, light-element burning, the
! timestep update, the stitched model + observables + output pass,
! and the per-model stop criteria. Exits:
!   step_status = step_continue        advance accepted, continue
!   step_status = step_kind_card_done  age/abundance stop reached
!   step_status = step_leave_run_loop  target radius crossed
!   (named constants from core/stop_conditions.f90)
!   ierr /= 0        error; run_yrec returns it to the CLI wrapper
!
! The blanket `save` is load-bearing exactly as in run_yrec: these
! locals were statics of program main and several carry state across
! model advances (the convergence bookkeeping, the envelope-CZ zone
! of the previous step); the reset block below (evolve_step_reset_pending)
! covers them alongside evolve_state.
subroutine evolve_step(model_iteration, step_status, ierr)

      use net_lib
      use star_info_lib, only: star, evolve_step_reset_pending, i_h1, i_h2, i_he4, &
           n_species_basic, n_species_extended
      use luout_lib
      use phys_const_lib
      use burn_lib
      use yrec_output, only: output_write_model
      use observables_lib, only: compute_observables
      use stitched_model_lib, only: build_stitched_model
      use stop_conditions
      implicit none

! the run index (former common/zramp/ NK) is star%job%nk (2026
! phase-A eviction), not an argument
      integer, intent(in) :: model_iteration
      integer, intent(out) :: step_status, ierr

      double precision :: delta_lum_step, delta_pressure_step, &
           delta_radius_step, delta_temp_step, &
           target_envelope_mass
      integer :: envelope_cz_zone_end, envelope_cz_zone_prev, &
           iteration_level, iterations_done, max_iterations, &
           num_mixed_zones, num_mixed_zones_no_overshoot, &
           num_radiative_zones, num_species, i, j, ii, itrot, jerr
      logical :: converged, evolve_model_flag, &
           new_atmosphere_fit_needed, recompute_surface_bc, &
           wind_loss_active
! compute_timestep's report that it switched the structure-based
! timestep limits off (log Tc > 7.1, L_grav < 0); the driver applies
! it to star%job (2026 W3, sidechan Deferred 4).
      logical :: disable_structure_dt_limits
! load-bearing: see header
      save
   ! INTENTIONAL: cross-step driver state; reset via evolve_step_reset_pending
      step_status = step_continue
      ierr = 0

! 2026 (phase five, step C): on a repeated run_yrec call, put this
! routine's SAVEd locals back to process-start state (static zero).
! Set by yrec_reset_lib's prologue; a fresh
! process never sees the flag, so call 1 is untouched.
      if (evolve_step_reset_pending) then
         delta_lum_step = 0.0d0
         delta_pressure_step = 0.0d0
         delta_radius_step = 0.0d0
         delta_temp_step = 0.0d0
         target_envelope_mass = 0.0d0
         envelope_cz_zone_end = 0
         envelope_cz_zone_prev = 0
         iteration_level = 0
         iterations_done = 0
         max_iterations = 0
         num_mixed_zones = 0
         num_mixed_zones_no_overshoot = 0
         num_radiative_zones = 0
         num_species = 0
         converged = .false.
         evolve_model_flag = .false.
         new_atmosphere_fit_needed = .false.
         recompute_surface_bc = .false.
         wind_loss_active = .false.
         evolve_step_reset_pending = .false.
      end if

! One model advance, as named phases (2026, core/ phase 3). Each
! phase is an internal subroutine below; error paths
! set ierr (checked after each call), divergence sets
! star%model_diverged_flag (host cycles the retry loop).
      call update_output_flags_for_step

! read_starting_model is re-run inside the loop for timestep cutting
! (the retry_step loop replaces the historical backward goto 15;
! divergence bailouts cycle it.)
      retry_step: do
         call reload_model_if_diverged
         if (ierr /= 0) return
         call advance_composition_and_age
         if (ierr /= 0) return
         call rezone_or_snapshot
         if (ierr /= 0) return
         call solve_structure
         if (ierr /= 0) return
         if (star%model_diverged_flag) cycle retry_step
         call converge_with_rotation
         if (ierr /= 0) return
         if (star%model_diverged_flag) cycle retry_step
         exit retry_step
      end do retry_step
! 2026 (bugsweep Batch 3): a non-finite converged structure is an
! error, not a model (see converged_model_is_nan in stop_conditions).
      if (converged_model_is_nan()) then
         ierr = 1
         return
      end if
! LOCATE THE HYDROGEN-BURNING SHELL AND THE BOUNDARIES OF THE CENTRAL
! AND SURFACE CONVECTION ZONES (IF APPLICABLE).
       call locate_shell_boundaries(star%xa,star%luminosity_lsun,star%convective_flag,star%nz, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%h_shell_zone_begin,star%h_shell_end_index,star%h_shell_midpoint_zone, &
              star%has_h_shell)
! PERFORM LIGHT ELEMENT BURNING (2026: extracted below)
         call burn_light_elements
! INCLUDED
         if (new_atmosphere_fit_needed) then
            call rebuild_envelope(target_envelope_mass,star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%m,star%dm, &
                            star%log_total_mass,star%logT,star%convective_flag,star%eta_squared,star%i_rot,star%j_rot,star%qiw,star%mean_radius, &
                            star%kinetic_energy_rot,star%log_L,star%total_angular_momentum,star%total_rotational_ke,star%log_Teff,star%nz,star%recompute_envelope_triangle,ierr)
            if (ierr /= 0) return
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
            call rotation_shape_factors(star%logRho,star%logR,star%log_mass,star%nz,star%omega,star%eta_squared,star%fp_rot,star%ft_rot,star%mean_gravity,star%mean_radius,ierr)
            if (ierr /= 0) return
            new_atmosphere_fit_needed = .false.
         endif
! DETERMINE TIMESTEP FOR NEXT MODEL
! compute_timestep ALSO LOCATES THE H-BURNING SHELL
       star%dt = dabs(star%dt)
       star%dt_saved = star%dt
       disable_structure_dt_limits = .false.
       call compute_timestep(star%dt,star%chosen_dt,star%nz,star%logRho,star%luminosity_lsun,star%m,star%dm,star%logT,star%xa,star%core_cz_top_index, &
              star%h_shell_midpoint_zone,star%luminosity_breakdown,star%dage,star%timestep_yr,star%job%nk,star%logP,star%logR,star%omega, &
              star%max_domega_frac,star%h_shell_zone_begin,star%log_Teff, disable_structure_dt_limits)
       if (disable_structure_dt_limits) star%job%use_structure_dt_limits = .false.
! IF EVOLVING TO A GIVEN AGE AND KIND CARD IS DONE, AVOID ZEROING OUT
! TIMESTEP WRITTEN TO MODEL (AS THIS MAKES CONTINUING A SEQUENCE AWKWARD.)
!     INSTEAD WRITE THE PREVIOUS MODEL TIMESTEP TO MODEL.
! ONLY IF A FIXED END AGE IS USED, NOT FOR OTHER STOPS
       if (star%job%end_age_stop_active(star%job%nk) .and. star%job%target_end_age(star%job%nk).gt.0.0D0) then
          if (reached_end_age(star%job%nk)) then
             star%dt = max(star%dt_saved,1.0D-3*star%dage*seconds_per_year)
             star%timestep_yr = star%dt/seconds_per_year
          else
             star%dt_saved = star%dt
          endif
       else
          star%dt_saved = star%dt
       endif
       if (star%job%rescale_kind(star%job%nk).ne.2) star%model_number = star%model_number+1
! 2026 stitched-model restructure: assemble the full converged star
! (interior + envelope + atmosphere) and materialize the stitched
! arrays BEFORE the observables pass. compute_observables' turnover
! theme walks them, star%pphot comes from this build, and the io
! writers below only copy -- nothing downstream integrates anything.
       call build_stitched_model
! 2026 (phase four, step 5): compute the per-model observables in
! the star layer (fills star%*, star%luminosity_breakdown
! renormalization, star% via gettau); wrtout below only
! reads. One theme subroutine per quantity family -- see
! core/observables_lib.f90.
       call compute_observables(ierr)
       if (ierr /= 0) return
! WRTOUT IS THE OUTPUT DRIVER ROUTINE
       call output_write_model()

! MHP 10/24 GENERALIZED STOP CONDITIONS
! (2026: one entry in stop_conditions -- end-age stop, the D/X/Y
! central-abundance table, and the star-calibration target-radius
! check, in that order)
      call check_stop_conditions(model_iteration, step_status)
      return
contains

! ---------------------------------------------------------------
! Per-step run-log bookkeeping: when the lrwsh control
! (rewind_short_file) is set, rewind the run log (run_log_unit,
! log_output_file / run.log) so it holds only the current model's
! solver diagnostics.
subroutine update_output_flags_for_step

          if (star%ctrl%rewind_short_file) then
             rewind(run_log_unit)
          endif
end subroutine update_output_flags_for_step

! ---------------------------------------------------------------
! If the previous attempt diverged, re-read the starting model with
! the cut timestep (starin) and re-check the rotation configuration.
! Sets ierr on failure. Always marks the punch (model-in-memory)
! flag on the way out.
subroutine reload_model_if_diverged
            if (star%model_diverged_flag) then
             call read_starting_model(star%timestep_yr, star%dt, star%chosen_dt, &
                  star%trial_sign_flag, star%ikut_flag, star%istore_flag, &
                  star%model_diverged_flag, star%recompute_envelope_triangle, star%job%nk, &
                  star%dlnrho_dlnp, star%dlnrho_dlnt, star%total_angular_momentum, &
                  star%total_rotational_ke, star%convective_velocity, &
                  star%job%mixture_weights, ierr)
             if (ierr /= 0) return
             call check_rotation_initialised(ierr)
             if (ierr /= 0) return
          endif
          star%punch_pending_flag = .true.
end subroutine reload_model_if_diverged

! ---------------------------------------------------------------
! The composition half of the step, when the model is to be aged:
! wind mass loss, store the start-of-step composition, burn + mix
! (mix pipeline), and advance the age. Sets ierr on failure.
subroutine advance_composition_and_age
! skip this section if model not to be aged
! MHP 7/98
! need to add logic to permit resacling + time evolution for
! pre-main sequence models
! (2026, same logic as the former if-ladder: rescale-only kind cards
! (rescale_kind = 2) skip aging EXCEPT while the center is still cool
! -- pre-main-sequence models rescale and age at the same time.)
            evolve_model_flag = star%model_number.ge.0 .and. &
                 (star%job%rescale_kind(star%job%nk).ne.2 .or. star%logT(1).lt.6.6D0)
            new_atmosphere_fit_needed = .false.
            if (evolve_model_flag) then
! ADD MASS LOSS CALCULATION
               call massloss(star%log_L,star%dt,star%xa,star%logRho,star%j_rot,star%logP,star%logR, &
                             star%log_mass,star%m,star%dm,star%log_total_mass,star%logT,star%envelope_cz_bottom_index,star%recompute_envelope_triangle, &
                             star%nz,star%omega,star%star_mass,star%log_Teff,target_envelope_mass,new_atmosphere_fit_needed, ierr)
               if (ierr /= 0) return
! STORE COMPOSITION MATRIX AT THE BEGINNING OF THE TIMESTEP.
               num_species = n_species_basic
               if (star%job%use_extended_composition) num_species=n_species_extended
               do i = 1,star%nz
                  do j = 1,num_species
                     star%xa_start(j,i) = star%xa(j,i)
                  end do
               end do
               iteration_level=1
! mixed_zone_bounds_no_overshoot stays an ARGUMENT of mix (not read as
! star% inside it) because henyey_iterate passes its own local array
! there -- storage deliberately separate from this one. The star copy
! is passed explicitly here.
               call mix(star%dt, iteration_level, star%timestep_yr, &
                    star%core_cz_top_index, star%envelope_cz_bottom_index, &
                    star%mixed_zone_bounds_no_overshoot, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
             star%timestep_yr = star%dt/seconds_per_year
             star%dage = star%dage + 1.0D-9*star%timestep_yr
            endif
end subroutine advance_composition_and_age

! ---------------------------------------------------------------
! Rezone the new model (hpoint), or -- during He-flash timestep
! cutting -- snapshot the previous structure instead; store the
! rotational KE distribution and the light-element burning
! start-of-step state (lirate88). Sets ierr on failure.
subroutine rezone_or_snapshot
      use math_lib
! rezone new model, except rezoning not performed for He flash calculations
          if (.not.star%ctrl%helium_flash_active) then
             call rezone(star%istore_flag, star%reset_triangle, star%h_shell_zone_begin, &
                  star%has_h_shell, star%total_angular_momentum, &
                  star%total_rotational_ke, ierr)
             if (ierr /= 0) return
! STORE NEW CZ BASE
               star%jcz = star%envelope_cz_bottom_index
            else
! save old model for PTIME
               do i=1, star%nz
                  star%logP_start(i) = star%logP(i)
                  star%logT_start(i) = star%logT(i)
                  star%logR_start(i) = star%logR(i)
                  star%luminosity_lsun_start(i) = star%luminosity_lsun(i)
                  star%logRho_start(i) = star%logRho(i)
               end do
! Teff and the zone count are saved as well
               star%log_Teff_start = star%log_Teff
               star%nz_start = star%nz
          endif
! store starting distribution of rotational kinetic energy.
            if (star%job%rotation_active) then
               do i = 1,star%nz
                  star%kinetic_energy_rot_old(i) = star%kinetic_energy_rot(i)
               end do
            endif
! changed for lithium burning with overshoot.
! store starting depth of C.Z. for light element burning.
            if (star%job%use_extended_composition) then
               star%cz_base_radius_prev = 0.0D0
               envelope_cz_zone_prev = star%envelope_cz_bottom_index
               star%pressure_scale_height_start = envelope_overshoot_depth()
! find burning rates at the beginning of the time step.
               call lirate88(star%xa,star%logRho,star%logT,star%nz,1)
            endif
end subroutine rezone_or_snapshot

! ---------------------------------------------------------------
! Corrector levels 1-3 (solve_level wrapper): level 1 on the cached
! envelope triangle, level 2 re-verifying the surface boundary,
! then the entropy-rate terms for the next initial guess, then
! level 3. Returns early with the divergence flag set for the host
! to cycle the retry loop; ierr on error.
subroutine solve_structure
! begin correction routines
! set flags for henyey_iterate:
! if ctrl%lnew0 = T, new envelope triangle calculated the 1st iteration
! (i.e. old triangle ignored)
! converged = T if model has converged
! star%model_diverged_flag = T if model has diverged
          if (star%ctrl%lnew0) star%recompute_envelope_triangle = .true.
! "do not age this model" is encoded as a NEGATIVE star%dt: the
! henyey pass (delta_time.gt.0 gates the entropy terms), the
! improved-first-guess test just below, and the rotation step
! (evolve_angular_momentum, which also does arithmetic on it) all
! read the sign; it is made positive again after convergence
! (star%dt = dabs(star%dt) before compute_timestep).
            if (.not.evolve_model_flag) star%dt = -dabs(star%dt)
! under-relaxation start: |fcorr0| - fcorri; the SIGN of the
! fcorr0 control selects whether henyey_iterate ramps fcorr back up
! by fcorri per iteration (fcorr0 > 0) or holds it fixed.
            star%job%fcorr = dabs(star%ctrl%fcorr0) - star%ctrl%fcorri
            iterations_done = 0
            star%model_diverged_flag = .false.
            converged = .false.
            if (.not.star%ctrl%improved_first_guess_flag .or. star%dt.le.0.0D0) then
               do i = 1,star%nz
! zero entropy terms
                  star%log_temperature_delta(i) = 0.0D0
                  star%log_pressure_delta(i) = 0.0D0
                  star%temperature_entropy_term(i) = 0.0D0
                  star%pressure_entropy_term(i) = 0.0D0
                  star%luminosity_entropy_term(i) = 0.0D0
                  star%radius_entropy_term(i) = 0.0D0
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
               end do
            else
! use the rate of change in the previous model to estimate the new
! run of structure variables.
               do i = 1,star%nz
                  delta_temp_step = star%temperature_entropy_term(i)*star%dt
                  delta_pressure_step = star%pressure_entropy_term(i)*star%dt
                  delta_lum_step = star%luminosity_lsun(i)*star%luminosity_entropy_term(i)*star%dt
                  delta_radius_step = star%radius_entropy_term(i)*star%dt
                  star%log_temperature_delta(i) = delta_temp_step
                  star%log_pressure_delta(i) = delta_pressure_step
                  star%logT(i) = star%logT(i) + delta_temp_step
                  star%logP(i) = star%logP(i) + delta_pressure_step
                  star%luminosity_lsun(i) = star%luminosity_lsun(i) + delta_lum_step
                  star%logR(i) = star%logR(i) + delta_radius_step
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
               end do
            endif

! FIRST LEVEL OF ITERATIONS
! USE ENVELOPE TRIANGLE OF THE PREVIOUS MODEL;
! FOR THE FIRST MODEL OF A RUN,THE TRIANGLE IS GENERATED HERE.
! (2026: the four-level ladder goes through the solve_level wrapper
! below; per-level differences are only (level, max iterations,
! recompute surface BC))
            call solve_level(1, star%ctrl%max_iter_level1, .false.)
            if (ierr /= 0) return
! SECOND LEVEL OF ITERATIONS
! CHECK ENVELOPE TRIANGLE BEFORE ITERATING FOR SOLUTION
            if (star%model_diverged_flag) return   ! (host cycles retry_step on the flag)
            call solve_level(2, star%ctrl%max_iter_level2, .true.)
            if (ierr /= 0) return
            if (star%model_diverged_flag) return   ! (host cycles retry_step on the flag)
! 7/91 STORE CHANGES IN THE STRUCTURE. THESE CHANGES ARE USED TO GET AN
! IMPROVED FIRST GUESS AT THE STRUCTURE FOR THE NEXT MODEL IF LNEWS=T.
            if (star%dt.gt.0.0D0) then
               do ii = 1,star%nz
                  star%temperature_entropy_term(ii)=star%log_temperature_delta(ii)/star%dt
                  star%pressure_entropy_term(ii)=star%log_pressure_delta(ii)/star%dt
                  star%luminosity_entropy_term(ii)=2.0D0*(star%luminosity_lsun(ii)-star%luminosity_lsun_start(ii))/(star%luminosity_lsun(ii)+star%luminosity_lsun_start(ii))/star%dt
                  star%radius_entropy_term(ii)=(star%logR(ii)-star%logR_start(ii))/star%dt
               end do
            endif
! THIRD LEVEL OF ITERATIONS
            call solve_level(3, star%ctrl%max_iter_level3, .false.)
            if (ierr /= 0) return
            if (star%model_diverged_flag) return   ! (host cycles retry_step on the flag)
end subroutine solve_structure

! ---------------------------------------------------------------
! The structure <-> rotation iteration (num_rotation_structure_iters passes): corrector
! level 4 INSIDE the loop, convergence check, convection-zone
! re-mix (mixcz), then the angular momentum update (getw) and the
! new shape factors (fpft). Start-of-timestep composition and
! angular momentum are restored on repeat passes so nothing is
! double-counted.
subroutine converge_with_rotation
            if (.not.star%job%rotation_active) then
               star%job%num_rotation_structure_iters = 1
            endif
! MHP 05/02
! IF THE CODE IS ITERATING BETWEEN THE STRUCTURE AND ROTATION
! SOLUTIONS, ENSURE THAT THE START-OF-TIMESTEP QUANTITIES
! HCOMPP (COMPOSITION) AND HJMSAV (ANGULAR MOMENTUM) ARE ONLY
! OVERWRITTEN ON THE LAST RUN THROUGH.
            if (star%job%num_rotation_structure_iters.gt.1) then
               do i = 1,star%nz
                  star%orig_specific_angular_momentum(i) = star%j_rot(i)
                  do j = 1,n_species_extended
                     star%orig_composition(j,i) = star%xa_start(j,i)
                  end do
               end do
            endif
            do itrot = 1, star%job%num_rotation_structure_iters
! MHP 05/02 RESTORE ORIGINAL "START OF TIMESTEP"
! VALUES FOR THE COMPOSITION MATRIX
               if (itrot.gt.1) then
                  do i = 1,star%nz
                     do j = 1,n_species_extended
                        star%xa_start(j,i) = star%orig_composition(j,i)
                     end do
                  end do
               endif
! 7/91 THE FOURTH LEVEL OF ITERATION REPEATS THE ITERATION BETWEEN THE
! MIXING AND THE STRUCTURE VARIABLES.  IT SHOULD NOT BE USED FOR MODELS
! WHERE SEMI-CONVECTION IS IMPORTANT (ITERATING BETWEEN THE BURNING AND
! STRUCTURE GENERATES OSCILLATIONS). IT SHOULD BE USED FOR HIGH-PRECISION
! WORK (E.G. SOLAR MODELS).
! Surface boundary conditions checked again since we've changed the
! star%xa (and hence the structure) of the model in ITLVL=3
! (to be implemented when I know the rest of it works!)
! (level 4 runs INSIDE the structure<->rotation iteration above,
! unlike levels 1-3 -- deliberate, see the itrot loop)
            call solve_level(4, star%ctrl%max_iter_level4, .false.)
            if (ierr /= 0) return
            if (.not.converged) then
! MODEL FAILED TO CONVERGE WITHIN(NITER1+NITER2+NITER3+NITER4)ITERATIONS
               star%model_diverged_flag = .true.
               return   ! (host cycles retry_step on the flag)
            endif

! MODEL HAS CONVERGED
! ENSURE CONVECTION ZONES ARE FULLY MIXED.
      call homogenize_convection_zones(star%xa,star%dm,star%convective_flag,star%nz)

! STORE TOTAL AGE FOR THE DISK-LOCKING GATE
            star%disk_gate_age_gyr = star%dage
            if (star%job%rotation_active) then
! RESTORE ORIGINAL START OF TIMESTEP VALUES
! TO THE ANGULAR MOMENTUM DISTRIBUTION
               if (itrot.gt.1) then
                  do i = 1,star%nz
                     star%j_rot(i) = star%orig_specific_angular_momentum(i)
                  end do
               endif
! (2026: the print_rotation_diagnostics pre-arm that lived here --
! MHP 9/94's end-of-card QUAD/PHIS terminal line, the flag's only
! effect -- is retired along with the flag.)
! FIND THE NEW RUN OF OMEGA
               wind_loss_active = star%job%use_wind_torque
               call evolve_angular_momentum(star%dt, star%max_domega_frac, wind_loss_active, &
                    envelope_cz_zone_prev, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
               call rotation_shape_factors(star%logRho,star%logR,star%log_mass,star%nz,star%omega,star%eta_squared,star%fp_rot,star%ft_rot,star%mean_gravity,star%mean_radius,ierr)
               if (ierr /= 0) return
            endif
            end do
end subroutine converge_with_rotation

! ---------------------------------------------------------------
! Li/Be/H2 burning over the step for non-rotating (or non-mixing)
! models: end-of-step convection-zone depth (convec), end-of-step
! rates (lirate88), then the burn (liburn).
subroutine burn_light_elements
      use math_lib
! PERFORM LIGHT ELEMENT BURNING
         if (star%job%use_extended_composition .and. star%model_number.ge.0 .and. star%dt.gt.0.0D0) then
! ONLY FOR MODELS WITHOUT ROTATION, OR WITHOUT ROTATIONAL MIXING.
            if (.not.star%job%rotation_active .or. .not.star%job%instability_transport_active) then
! FIND CONVECTION ZONE DEPTH AT THE END OF THE TIME STEP.
               call find_convection_zones(star%xa,star%logRho,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%nz,star%radiative_zone_bounds,star%mixed_zone_bounds, &
                            star%mixed_zone_bounds_no_overshoot,star%core_cz_top_index,star%envelope_cz_bottom_index,num_radiative_zones,num_mixed_zones,num_mixed_zones_no_overshoot, ierr)
               if (ierr /= 0) return
! CHANGED FOR LITHIUM BURNING WITH OVERSHOOT.
               envelope_cz_zone_end = star%envelope_cz_bottom_index
               star%pressure_scale_height_end = envelope_overshoot_depth()
! FIND BURNING RATES AT THE END OF THE TIME STEP.
               call lirate88(star%xa,star%logRho,star%logT,star%nz,2)
               call liburn(star%dt,star%xa,star%logR,star%m,star%dm,star%logT,envelope_cz_zone_end,envelope_cz_zone_prev,star%nz)
            endif
         endif
end subroutine burn_light_elements


! ---------------------------------------------------------------
! One level of the corrector ladder: same henyey_iterate call for every
! level; only the level number, its iteration budget, and whether
! the surface boundary condition is re-verified differ. All other
! arguments are the host's locals / star%evo state, via host
! association. ierr and the divergence flag are checked by the
! caller after each level.
subroutine solve_level(level, level_max_iterations, check_surface_bc)
      integer, intent(in) :: level, level_max_iterations
      logical, intent(in) :: check_surface_bc

      iteration_level = level
      max_iterations = level_max_iterations
      recompute_surface_bc = check_surface_bc
      call henyey_iterate(star%dt, max_iterations, converged, &
           star%model_diverged_flag, star%recompute_envelope_triangle, &
           star%reset_triangle, recompute_surface_bc, star%trial_sign_flag, &
           star%istore_flag, star%dlnrho_dlnt, star%dlnrho_dlnp, &
           iterations_done, iteration_level, ierr)
end subroutine solve_level

! ---------------------------------------------------------------
! Envelope-overshoot penetration depth (cm): overshoot_alpha_envelope
! times the pressure scale height at the base of the envelope
! convection zone, or zero when envelope overshoot is off. Same
! formula at the start and the end of the step (rezone_or_snapshot,
! burn_light_elements).
double precision function envelope_overshoot_depth()
      use math_lib
      if (star%job%envelope_overshoot_active) then
         envelope_overshoot_depth = star%ctrl%overshoot_alpha_envelope*exp(clndp*(star%logP(star%envelope_cz_bottom_index)+2.0D0*star%logR(star%envelope_cz_bottom_index) &
              -star%logRho(star%envelope_cz_bottom_index)-cgl-star%log_mass(star%envelope_cz_bottom_index)))
      else
         envelope_overshoot_depth = 0.0D0
      endif
end function envelope_overshoot_depth

end subroutine evolve_step

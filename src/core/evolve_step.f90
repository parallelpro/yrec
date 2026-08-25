!----------------------------------------------------------------------
! evolve_step
!----------------------------------------------------------------------
! Added 2026 (phase five -- the embeddable engine; ROADMAP.md). One
! model advance, moved verbatim from run_yrec's model loop body: the
! divergence/redo cycle (label 15), the corrector iteration via
! crrect, mixing, rotation (getw), rezoning (hpoint), the output
! diagnostics + wrtout, timestep update, and the per-model stop
! criteria. Exits:
!   step_status = step_continue        advance accepted, continue
!   step_status = step_kind_card_done  age/abundance stop reached
!   step_status = step_leave_run_loop  target radius crossed
!   (named constants from core/stop_conditions.f90)
!   ierr /= 0        error; run_yrec returns it to the CLI wrapper
! (Internally the historical goto 110/200 jumps now target the 810/
! 820 exit labels below -- the control flow is otherwise untouched.)
!
! The blanket `save` is load-bearing exactly as in run_yrec: these
! locals were statics of program main and several carry state across
! model advances (nao's data init, the convergence bookkeeping); the
! phase-C reset audit covers them alongside evolve_state.
subroutine evolve_step(model_iteration, step_status, ierr)

      use net_lib
      use star_info_lib, only: star, evolve_step_reset_pending, i_h1, i_h2, i_he4
      use luout_lib
      use const_lib
      use burn_lib
      use yrec_output, only: output_write_model
      use observables_lib, only: compute_observables
      use stop_conditions
      implicit none

! the run index (former common/zramp/ NK) is star%job%nk (2026
! phase-A eviction), not an argument
      integer, intent(in) :: model_iteration
      integer, intent(out) :: step_status, ierr

      double precision :: delta_lum_step, delta_pressure_step, &
           delta_radius_step, delta_temp_step, log_gravity, &
           target_envelope_mass, teff_kelvin_unused
      integer :: envelope_cz_zone_end, envelope_cz_zone_prev, &
           iteration_level, iterations_done, max_iterations, nao, &
           num_mixed_zones, num_mixed_zones_no_overshoot, &
           num_radiative_zones, num_species, i, j, ii, itrot, jerr
      logical :: conductive_opacity_flag, converged, &
           evolve_model_flag, in_atmosphere, mixing_active, &
           new_atmosphere_fit_needed, recompute_surface_bc, &
           use_correct_gradients, want_derivatives, wind_loss_active
      data nao/1/
! load-bearing: see header
      save
   ! INTENTIONAL: cross-step driver state; reset via evolve_step_reset_pending
      step_status = step_continue
      ierr = 0

! 2026 (phase five, step C): on a repeated run_yrec call, put this
! routine's SAVEd locals back to process-start state (static zero;
! nao's data value). Set by yrec_reset_lib's prologue; a fresh
! process never sees the flag, so call 1 is untouched.
      if (evolve_step_reset_pending) then
         delta_lum_step = 0.0d0
         delta_pressure_step = 0.0d0
         delta_radius_step = 0.0d0
         delta_temp_step = 0.0d0
         log_gravity = 0.0d0
         target_envelope_mass = 0.0d0
         teff_kelvin_unused = 0.0d0
         envelope_cz_zone_end = 0
         envelope_cz_zone_prev = 0
         iteration_level = 0
         iterations_done = 0
         max_iterations = 0
         num_mixed_zones = 0
         num_mixed_zones_no_overshoot = 0
         num_radiative_zones = 0
         num_species = 0
         nao = 1
         conductive_opacity_flag = .false.
         converged = .false.
         evolve_model_flag = .false.
         in_atmosphere = .false.
         mixing_active = .false.
         new_atmosphere_fit_needed = .false.
         recompute_surface_bc = .false.
         use_correct_gradients = .false.
         want_derivatives = .false.
         wind_loss_active = .false.
         evolve_step_reset_pending = .false.
      end if

! One model advance, as named phases (2026, core/ phase 3). Each
! phase is an internal subroutine below, moved verbatim; error paths
! set ierr (checked after each call), divergence sets
! star%evo%model_diverged_flag (host cycles the retry loop).
      call update_output_flags_for_step

! STARIN called here for timestep cutting
! (Restructured 2026: the backward goto 15 retry became the named
! retry_step loop; divergence bailouts cycle it.)
      retry_step: do
         call reload_model_if_diverged
         if (ierr /= 0) return
         call advance_composition_and_age
         if (ierr /= 0) return
         call rezone_or_snapshot
         if (ierr /= 0) return
         call solve_structure
         if (ierr /= 0) return
         if (star%evo%model_diverged_flag) cycle retry_step
         call converge_with_rotation
         if (ierr /= 0) return
         if (star%evo%model_diverged_flag) cycle retry_step
         exit retry_step
      end do retry_step
! LOCATE THE HYDROGEN-BURNING SHELL AND THE BOUNDARIES OF THE CENTRAL
! AND SURFACE CONVECTION ZONES (IF APPLICABLE).
       call findsh(star%xa,star%luminosity_lsun,star%convective_flag,star%nz, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%evo%h_shell_zone_begin,star%evo%h_shell_end_index,star%evo%h_shell_midpoint_zone, &
              star%evo%has_h_shell)
! PERFORM LIGHT ELEMENT BURNING (2026: extracted below)
         call burn_light_elements
! INCLUDED
         if (new_atmosphere_fit_needed) then
            call getnewenv(target_envelope_mass,star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%m,star%dm, &
!     *                     HSTOT,HT,LC,ETA2,HG,HI,HJM,QIW,R0,  ! KC 2025-05-31
                            star%log_total_mass,star%logT,star%convective_flag,star%eta_squared,star%i_rot,star%j_rot,star%qiw,star%mean_radius, &
                            star%kinetic_energy_rot,star%log_L,star%evo%total_angular_momentum,star%evo%total_rotational_ke,star%log_Teff,star%nz,star%evo%recompute_envelope_triangle)
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
            call fpft(star%logRho,star%logR,star%log_mass,star%nz,star%omega,star%eta_squared,star%fp_rot,star%ft_rot,star%mean_gravity,star%mean_radius)
            new_atmosphere_fit_needed = .false.
         endif
! DETERMINE TIMESTEP FOR NEXT MODEL
! HTIMER ALSO LOCATES THE H-BURNING SHELL
! JVS 04/14 added teffl to passed htimer variables
       star%evo%dt = dabs(star%evo%dt)
       star%evo%dt_saved = star%evo%dt
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *        JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(star%evo%dt,star%evo%hydrogen_dt,star%nz,star%logRho,star%luminosity_lsun,star%m,star%dm,star%logT,star%xa,star%core_cz_top_index, &
              star%evo%h_shell_midpoint_zone,star%luminosity_breakdown,star%run%dage,star%evo%timestep_yr,star%job%nk,star%logP,star%logR,star%omega, &
              star%evo%max_domega_frac,star%evo%h_shell_zone_begin,star%log_Teff)
! IF EVOLVING TO A GIVEN AGE AND KIND CARD IS DONE, AVOID ZEROING OUT
! TIMESTEP WRITTEN TO MODEL (AS THIS MAKES CONTINUING A SEQUENCE AWKWARD.)
!     INSTEAD WRITE THE PREVIOUS MODEL TIMESTEP TO MODEL.
! ONLY IF A FIXED END AGE IS USED, NOT FOR OTHER STOPS
       if (star%job%end_age_stop_active(star%job%nk) .and. star%job%target_end_age(star%job%nk).gt.0.0D0) then
          if (reached_end_age(star%job%nk)) then
             star%evo%dt = max(star%evo%dt_saved,1.0D-3*star%run%dage*seconds_per_year)
             star%evo%timestep_yr = star%evo%dt/seconds_per_year
          else
             star%evo%dt_saved = star%evo%dt
          endif
       else
          star%evo%dt_saved = star%evo%dt
       endif
       if (star%job%rescale_kind(star%job%nk).ne.2) star%model_number = star%model_number+1
! 2026 (phase four, step 5): compute the per-model observables in
! the star layer (fills star%run%*, star%luminosity_breakdown
! renormalization, star%turnover% via gettau); wrtout below only
! reads. One theme subroutine per quantity family -- see
! core/observables_lib.f90.
       call compute_observables(ierr)
       if (ierr /= 0) return
! WRTOUT IS THE OUTPUT DRIVER ROUTINE
       call output_write_model(star%evo%timestep_yr, log_gravity, star%evo%has_h_shell, &
            star%evo%h_shell_zone_begin, star%evo%h_shell_midpoint_zone, star%evo%h_shell_end_index, &
            star%evo%trial_sign_flag, star%evo%punch_pending_flag, star%evo%total_angular_momentum, &
            star%evo%total_rotational_ke)

! MHP 10/24 GENERALIZED STOP CONDITIONS
!     IF EVOLVING TO A GIVEN AGE AND AGE IS REACHED, KIND CARD IS DONE
!       IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)GOTO 110
       if (reached_end_age(star%job%nk)) then
          step_status = step_kind_card_done
          return
       end if
! MHP 10/24 CHECK ALL STOP CONDITIONS, EXIT IF ANY SATISFIED
! (2026: the D/X/Y checks are one table walk in stop_conditions)
         if (abundance_stop_triggered(star%job%nk)) then
! SET I/O FLAGS PROPERLY AND EXIT LOOP
            pulsation_output_active = star%evo%saved_pulse_output_flag
            star%run%sound_speed_output_active = .true.
            star%run%print_rotation_diagnostics = .true.
            step_status = step_kind_card_done
            return
         endif
! TEST IF MODEL IS NEAR DESIRED Teff AND L. IF NOT RESCALE AND TRY AGAIN.
         if (star%ctrl%calibrate_star_flag .and. .not. star_found_flag) then
            if (mod(star%job%nk,2).eq.0) then
! chkscal protocol: iteration 1 only primes its previous-model state
! (the value computed here is never read -- see chkscal.f90)
             if (model_iteration.eq.1) then
                teff_kelvin_unused = 10.0D0**star%log_Teff
             else
                call chkscal(star%log_L, star%log_Teff, star%run%dage, star%job%nk)
                if (just_passed_target_radius_flag) then
                   step_status = step_leave_run_loop
                   return
                end if
             end if
          endif
       endif

! END OF RUN

      step_status = step_continue
      return
contains

! ---------------------------------------------------------------
! Per-step output-flag choreography: optional short-file rewind,
! re-arming pulse output for the last model of the last run and for
! the JVS ages-of-interest brackets (ageout), pre-arming the final
! model of an age-stopped run, and the pulse path-length bookkeeping
! (pdist).
subroutine update_output_flags_for_step

! rewind ISHORT if LRWSH is true (keeps ISHORT small)
          if (star%ctrl%rewind_short_file) then
             rewind(short_file_unit)
          endif
! DBG PULSE:  if last model of last run then set LPULSE to LSAVPU
            if (model_iteration.eq.star%job%num_models(star%job%nk) .and. star%job%nk .eq. star%job%num_runs) then
                 pulsation_output_active = star%evo%saved_pulse_output_flag
            end if

! JVS 02/11: Also allow pulse output at particular ages along the way
!
!  If the step is bracketing an age of interest, turn on output. This will
! for the step before and step after the age in AGEOUT. Once the info has
! been printed out, AGEOUT is set to the next age.
!
! Turn on calcad:
      if (star%ctrl%acoustic_depth_output) then
            compute_acoustic_depth=.true.
      else
            compute_acoustic_depth = .false.
      endif
! If output has been turned on for a previous step, keep it on for the next
! step, but then turn it off.
      if (star%ctrl%acoustic_depth_output) then
            if (ageout_bracket_armed) then
                  print*, 'LJWRT on'
                  pulsation_output_active = star%evo%saved_pulse_output_flag
                  nao=nao+1
                  calcad_ageout_output_active =.true.
                  ageout_model_output_flag =.false.
                  ageout_bracket_armed=.false.
            else if (.not.ageout_bracket_armed) then
                  calcad_ageout_output_active=.false.
            endif
! If this is the step before one of the ages of interest, print everything out.
! Also, save model structure.
            if (nao.lt.6) then
                  if (star%run%dage+star%evo%timestep_yr/1.0D9-star%ctrl%output_ages_gyr(nao) .le. 0.0D0 .and. &
                  star%run%dage+2.0D0*star%evo%timestep_yr/1.0D9-star%ctrl%output_ages_gyr(nao) .ge. 0.0D0 .and. .not. ageout_bracket_armed) then
                        print*, 'AGEOUT reached'
                        pulsation_output_active = star%evo%saved_pulse_output_flag
                        calcad_ageout_output_active = .true.
                        ageout_model_output_flag = .true.
                        ageout_bracket_armed=.true.
                  endif
            endif
       endif


! JVS end
!
! DBG PULSE:  if endage reached then set LPULSE to LSAVPU
! MHP 10/24 GENERALIZE CHECK
         if (approaching_end_age(star%job%nk)) then
                 pulsation_output_active = star%evo%saved_pulse_output_flag
! MHP 7/96 compute sound speed for solar model
                 star%run%sound_speed_output_active = .true.
            end if


!FD echo LSOUND
!        print*,'MAIN LSOUND = ',LSOUND
!FD end
            if (star%ctrl%po_output_enabled) then
! MHP 8/25 changed to add file names as declared variables
             call pdist(star%evo%prev_log_l,star%evo%prev_log_teff,star%evo%prev_age,star%evo%path_length_sq,star%log_L,star%log_Teff,model_iteration,star%job%pulse_atm_path, &
             star%job%pulse_env_path,star%job%pulse_mod_path)
          endif
end subroutine update_output_flags_for_step

! ---------------------------------------------------------------
! If the previous attempt diverged, re-read the starting model with
! the cut timestep (starin) and re-check the rotation configuration.
! Sets ierr on failure. Always marks the punch (model-in-memory)
! flag on the way out.
subroutine reload_model_if_diverged
            if (star%evo%model_diverged_flag) then
!              CALL STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,  ! KC 2025-05-31
             call starin(star%evo%timestep_yr, star%evo%dt, star%evo%hydrogen_dt, &
                  star%evo%trial_sign_flag, star%evo%ikut_flag, star%evo%istore_flag, &
                  star%evo%model_diverged_flag, star%evo%recompute_envelope_triangle, star%job%nk, &
                  star%evo%dlnrho_dlnp, star%evo%dlnrho_dlnt, star%evo%total_angular_momentum, &
                  star%evo%total_rotational_ke, star%evo%convective_velocity, &
                  star%job%mixture_weights, ierr)
             if (ierr /= 0) return
             if ((star%omega(1) .eq. 0) .and. (rotation_active)) then
18               format('LROT set to TRUE, but OMEGA(1) = 0. Stopping.', &
                        ' Initialize rotation rates or set LROT to', &
                        ' FALSE.')
                 print 18
                 ! 2026 (phase five, step B): configuration error returns to the
                 ! CLI wrapper (which stops) instead of stopping here.
                 ierr = 1
                 return
             endif
          endif
          star%evo%punch_pending_flag = .true.
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
               call massloss(star%log_L,star%run%dage,star%evo%dt,star%xa,star%logRho,star%j_rot,star%logP,star%logR, &
                             star%log_mass,star%m,star%dm,star%log_total_mass,star%logT,star%envelope_cz_bottom_index,star%evo%recompute_envelope_triangle, &
                             star%nz,star%omega,star%star_mass,star%log_Teff,target_envelope_mass,new_atmosphere_fit_needed)
! STORE COMPOSITION MATRIX AT THE BEGINNING OF THE TIMESTEP.
               num_species = 11
               if (use_extended_composition) num_species=15
               do i = 1,star%nz
                  do j = 1,num_species
                     star%prev%xa_start(j,i) = star%xa(j,i)
                  end do
               end do
               iteration_level=1
! mixed_zone_bounds_no_overshoot stays an ARGUMENT of mix (not read as
! star% inside it) because crrect passes its own local array there --
! storage deliberately separate from main's. main passes the star copy
! explicitly.
               call mix(star%evo%dt, iteration_level, star%evo%timestep_yr, &
                    star%core_cz_top_index, star%envelope_cz_bottom_index, &
                    star%mixed_zone_bounds_no_overshoot, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
             star%evo%timestep_yr = star%evo%dt/seconds_per_year
             star%run%dage = star%run%dage + 1.0D-9*star%evo%timestep_yr
            endif
end subroutine advance_composition_and_age

! ---------------------------------------------------------------
! Rezone the new model (hpoint), or -- during He-flash timestep
! cutting -- snapshot the previous structure instead; store the
! rotational KE distribution and the light-element burning
! start-of-step state (lirate88). Sets ierr on failure.
subroutine rezone_or_snapshot
!***MHP 1/04 OPACITY TEST
! DBG 12/95 GET OPACITY
!*** END TEST
! rezone new model, except rezoning not performed for He flash calculations
          if (.not.star%ctrl%helium_flash_active) then
             call hpoint(star%evo%istore_flag, star%evo%reset_triangle, star%evo%h_shell_zone_begin, &
                  star%evo%has_h_shell, star%evo%total_angular_momentum, &
                  star%evo%total_rotational_ke, ierr)
             if (ierr /= 0) return
! STORE NEW CZ BASE
               star%light_burn%jcz = star%envelope_cz_bottom_index
            else
! save old model for PTIME
               do i=1, star%nz
                  star%prev%logP_start(i) = star%logP(i)
                  star%prev%logT_start(i) = star%logT(i)
                  star%prev%logR_start(i) = star%logR(i)
                  star%prev%luminosity_lsun_start(i) = star%luminosity_lsun(i)
                  star%prev%logRho_start(i) = star%logRho(i)
               end do
! JVS 04/14 Save Teff as well
               star%prev%log_Teff_start = star%log_Teff
!  JVS 05/25 Added model number to list of saved values
           star%prev%nz_start = star%nz

          endif
! store starting distribution of rotational kinetic energy.
            if (rotation_active) then
               do i = 1,star%nz
                  star%kinetic_energy_rot_old(i) = star%kinetic_energy_rot(i)
               end do
            endif
! changed for lithium burning with overshoot.
! store starting depth of C.Z. for light element burning.
            if (use_extended_composition) then
               star%light_burn%cz_base_radius_prev = 0.0D0
               envelope_cz_zone_prev = star%envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  star%light_burn%pressure_scale_height_start = star%ctrl%alphae*exp(clndp*(star%logP(star%envelope_cz_bottom_index)+2.0D0*star%logR(star%envelope_cz_bottom_index) &
                           -star%logRho(star%envelope_cz_bottom_index)-cgl-star%log_mass(star%envelope_cz_bottom_index)))
               else
                  star%light_burn%pressure_scale_height_start = 0.0D0
               endif
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
! set flags for CRRECT
! CRRECT checks surface boundary conditions in every iteration
! if LNEW0 = T, new envelope triangle calculated the 1st iteration
! (i.e. old triangle ignored)
! LFINI = T if model has converged
! LARGE = T if model has diverged
          if (star%ctrl%lnew0) star%evo%recompute_envelope_triangle = .true.
            if (.not.evolve_model_flag) star%evo%dt = -dabs(star%evo%dt)
            fcorr = dabs(star%ctrl%fcorr0) - star%ctrl%fcorri
            iterations_done = 0
            star%evo%model_diverged_flag = .false.
            converged = .false.
            if (.not.star%ctrl%lnews .or. star%evo%dt.le.0.0D0) then
               do i = 1,star%nz
! zero entropy terms
                  star%log_temperature_delta(i) = 0.0D0
                  star%log_pressure_delta(i) = 0.0D0
                  star%run%temperature_entropy_term(i) = 0.0D0
                  star%run%pressure_entropy_term(i) = 0.0D0
                  star%run%luminosity_entropy_term(i) = 0.0D0
                  star%run%radius_entropy_term(i) = 0.0D0
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
               end do
            else
! use the rate of change in the previous model to estimate the new
! run of structure variables.
               do i = 1,star%nz
                  delta_temp_step = star%run%temperature_entropy_term(i)*star%evo%dt
                  delta_pressure_step = star%run%pressure_entropy_term(i)*star%evo%dt
                  delta_lum_step = star%luminosity_lsun(i)*star%run%luminosity_entropy_term(i)*star%evo%dt
                  delta_radius_step = star%run%radius_entropy_term(i)*star%evo%dt
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
! CALL TO CRRECT - ADDED ITERATION LEVEL (2026: the four-level ladder
! goes through the solve_level wrapper below; per-level differences
! are only (level, max iterations, recompute surface BC))
            call solve_level(1, star%ctrl%niter1, .false.)
            if (ierr /= 0) return
! SECOND LEVEL OF ITERATIONS
! CHECK ENVELOPE TRIANGLE BEFORE ITERATING FOR SOLUTION
            if (star%evo%model_diverged_flag) return   ! (host cycles retry_step on the flag)
            call solve_level(2, star%ctrl%niter2, .true.)
            if (ierr /= 0) return
            if (star%evo%model_diverged_flag) return   ! (host cycles retry_step on the flag)
! 7/91 STORE CHANGES IN THE STRUCTURE. THESE CHANGES ARE USED TO GET AN
! IMPROVED FIRST GUESS AT THE STRUCTURE FOR THE NEXT MODEL IF LNEWS=T.
            if (star%evo%dt.gt.0.0D0) then
               do ii = 1,star%nz
                  star%run%temperature_entropy_term(ii)=star%log_temperature_delta(ii)/star%evo%dt
                  star%run%pressure_entropy_term(ii)=star%log_pressure_delta(ii)/star%evo%dt
                  star%run%luminosity_entropy_term(ii)=2.0D0*(star%luminosity_lsun(ii)-star%prev%luminosity_lsun_start(ii))/(star%luminosity_lsun(ii)+star%prev%luminosity_lsun_start(ii))/star%evo%dt
                  star%run%radius_entropy_term(ii)=(star%logR(ii)-star%prev%logR_start(ii))/star%evo%dt
               end do
            endif
! THIRD LEVEL OF ITERATIONS
            call solve_level(3, star%ctrl%niter3, .false.)
            if (ierr /= 0) return
            if (star%evo%model_diverged_flag) return   ! (host cycles retry_step on the flag)
end subroutine solve_structure

! ---------------------------------------------------------------
! The structure <-> rotation iteration (itdif1 passes): corrector
! level 4 INSIDE the loop, convergence check, convection-zone
! re-mix (mixcz), then the angular momentum update (getw) and the
! new shape factors (fpft). Start-of-timestep composition and
! angular momentum are restored on repeat passes so nothing is
! double-counted.
subroutine converge_with_rotation
            if (.not.rotation_active) then
               itdif1 = 1
            endif
! MHP 05/02
! IF THE CODE IS ITERATING BETWEEN THE STRUCTURE AND ROTATION
! SOLUTIONS, ENSURE THAT THE START-OF-TIMESTEP QUANTITIES
! HCOMPP (COMPOSITION) AND HJMSAV (ANGULAR MOMENTUM) ARE ONLY
! OVERWRITTEN ON THE LAST RUN THROUGH.
            if (itdif1.gt.1) then
               do i = 1,star%nz
                  star%run%orig_specific_angular_momentum(i) = star%j_rot(i)
                  do j = 1,15
                     star%run%orig_composition(j,i) = star%prev%xa_start(j,i)
                  end do
               end do
            endif
            do itrot = 1, itdif1
! MHP 05/02 RESTORE ORIGINAL "START OF TIMESTEP"
! VALUES FOR THE COMPOSITION MATRIX
               if (itrot.gt.1) then
                  do i = 1,star%nz
                     do j = 1,15
                        star%prev%xa_start(j,i) = star%run%orig_composition(j,i)
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
            call solve_level(4, star%ctrl%niter4, .false.)
            if (ierr /= 0) return
!  25         CONTINUE
            if (.not.converged) then
! MODEL FAILED TO CONVERGE WITHIN(NITER1+NITER2+NITER3+NITER4)ITERATIONS
               star%evo%model_diverged_flag = .true.
               return   ! (host cycles retry_step on the flag)
            endif

! MODEL HAS CONVERGED
! ENSURE CONVECTION ZONES ARE FULLY MIXED.
! MHP 02/12 INFER CONVECTIVE OVERTURN TIMESCALE (USED IN MDOT)
! JVS 02/12 CALL MIXCZ(HCOMP,HS2,LC,M)
! KC 2025-05-30 addressed warning messages from Makefile.legacy
! C G Somers 6/14, SET IMIX = .TRUE. SO THE CORRECT GRADS ARE USED.
!       IMIX = .TRUE.
!       CALL MIXCZ(HCOMP,HS2,HS1,LC,HR,HP,HD,HG,M,IMIX)
! G Somers 6/14, SET LIMIX = .TRUE. SO THE CORRECT GRADS ARE USED.
      use_correct_gradients = .true.
!       CALL MIXCZ(HCOMP,HS2,HS1,LC,HR,HP,HD,HG,M,LIMIX)  ! KC 2025-05-31
      call mixcz(star%xa,star%dm,star%convective_flag,star%nz)
! G Somers END

! MHP 9/94 STORE TOTAL AGE IN SAGE
            disk_lifetime = star%run%dage
            if (rotation_active) then
! RESTORE ORIGINAL START OF TIMESTEP VALUES
! TO THE ANGULAR MOMENTUM DISTRIBUTION
               if (itrot.gt.1) then
                  do i = 1,star%nz
                     star%j_rot(i) = star%run%orig_specific_angular_momentum(i)
                  end do
               endif
! MHP 9/94 ADDED FLAG TO TURN ON ROTATION OUTPUT WHEN END OF KIND
! CARD REACHED.
! MHP 10/24 GENERALIZE CHECK
               if (approaching_end_age(star%job%nk)) then
                  star%run%print_rotation_diagnostics = .true.
               else
                  star%run%print_rotation_diagnostics = .false.
               endif
! FIND THE NEW RUN OF OMEGA
! JENV0 ADDED TO SR CALL.
               wind_loss_active = ljdot0
               call getw(star%evo%dt, star%evo%max_domega_frac, wind_loss_active, &
                    envelope_cz_zone_prev, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
               call fpft(star%logRho,star%logR,star%log_mass,star%nz,star%omega,star%eta_squared,star%fp_rot,star%ft_rot,star%mean_gravity,star%mean_radius)
            endif
            end do
end subroutine converge_with_rotation

! ---------------------------------------------------------------
! Li/Be/H2 burning over the step for non-rotating (or non-mixing)
! models: end-of-step convection-zone depth (convec), end-of-step
! rates (lirate88), then the burn (liburn).
subroutine burn_light_elements
! PERFORM LIGHT ELEMENT BURNING
         if (use_extended_composition .and. star%model_number.ge.0 .and. star%evo%dt.gt.0.0D0) then
! ONLY FOR MODELS WITHOUT ROTATION, OR WITHOUT ROTATIONAL MIXING.
            if (.not.rotation_active .or. .not.instability_transport_active) then
! FIND CONVECTION ZONE DEPTH AT THE END OF THE TIME STEP.
               call convec(star%xa,star%logRho,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%nz,star%radiative_zone_bounds,star%mixed_zone_bounds, &
                            star%mixed_zone_bounds_no_overshoot,star%core_cz_top_index,star%envelope_cz_bottom_index,num_radiative_zones,num_mixed_zones,num_mixed_zones_no_overshoot)
! CHANGED FOR LITHIUM BURNING WITH OVERSHOOT.
               envelope_cz_zone_end = star%envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  star%light_burn%pressure_scale_height_end = star%ctrl%alphae*exp(clndp*(star%logP(star%envelope_cz_bottom_index)+2.0D0*star%logR(star%envelope_cz_bottom_index) &
                           -star%logRho(star%envelope_cz_bottom_index)-cgl-star%log_mass(star%envelope_cz_bottom_index)))
               else
                  star%light_burn%pressure_scale_height_end = 0.0D0
               endif
! FIND BURNING RATES AT THE END OF THE TIME STEP.
               call lirate88(star%xa,star%logRho,star%logT,star%nz,2)
!                CALL LIBURN(DELTS,HCOMP,HD,HR,HS1,HS2,HT,JENV1,JENV0,M)  ! KC 2025-05-31
               call liburn(star%evo%dt,star%xa,star%logR,star%m,star%dm,star%logT,envelope_cz_zone_end,envelope_cz_zone_prev,star%nz)
            endif
         endif
! MHP 07/02 RESTORE PRIOR FITTING POINT IF MASS ACCRETION IS BEING
end subroutine burn_light_elements


! ---------------------------------------------------------------
! One level of the corrector ladder: same crrect call for every
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
      call crrect(star%evo%dt, max_iterations, converged, &
           star%evo%model_diverged_flag, star%evo%recompute_envelope_triangle, &
           star%evo%reset_triangle, recompute_surface_bc, star%evo%trial_sign_flag, &
           star%evo%istore_flag, in_atmosphere, want_derivatives, &
           mixing_active, conductive_opacity_flag, star%evo%dlnrho_dlnt, &
           star%evo%dlnrho_dlnp, iterations_done, iteration_level, ierr)
end subroutine solve_level

end subroutine evolve_step

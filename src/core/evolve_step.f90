!----------------------------------------------------------------------
! evolve_step
!----------------------------------------------------------------------
! Added 2026 (phase five -- the embeddable engine; ROADMAP.md). One
! model advance, moved verbatim from run_yrec's model loop body: the
! divergence/redo cycle (label 15), the corrector iteration via
! crrect, mixing, rotation (getw), rezoning (hpoint), the output
! diagnostics + wrtout, timestep update, and the per-model stop
! criteria. Exits:
!   step_status = 0  advance accepted, continue the model loop
!   step_status = 1  run finished (age/abundance stop) -> label 110
!                    handling in run_yrec (punch bookkeeping)
!   step_status = 2  leave the run loop entirely (target radius
!                    passed) -> label 200 in run_yrec
!   ierr /= 0        error; run_yrec returns it to the CLI wrapper
! (Internally the historical goto 110/200 jumps now target the 810/
! 820 exit labels below -- the control flow is otherwise untouched.)
!
! The blanket `save` is load-bearing exactly as in run_yrec: these
! locals were statics of program main and several carry state across
! model advances (nao's data init, the convergence bookkeeping); the
! phase-C reset audit covers them alongside evolve_state.
subroutine evolve_step(model_iteration, step_status, ierr)

      use nuclear_lib
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use luout_lib
      use const_lib
      use star_info_lib, only: star
      use star_job_lib, only: job
      use evolve_state_lib, only: evo, evolve_step_reset_pending
      use burn_lib
      implicit none

! nk (the run index) is const_lib module state (former common/zramp/),
! not an argument
      integer, intent(in) :: model_iteration
      integer, intent(out) :: step_status, ierr

      double precision :: delta_lum_step, delta_pressure_step, &
           delta_radius_step, delta_temp_step, log_gravity, &
           target_envelope_mass, teff_kelvin_unused
      integer :: envelope_cz_zone_end, envelope_cz_zone_prev, &
           iteration_level, iterations_done, max_iterations, nao, &
           num_mixed_zones, num_mixed_zones_no_overshoot, &
           num_radiative_zones, num_species, i, j, ii, itrot, jerr
      logical :: conductive_opacity_flag, converged, end_kind_flag, &
           evolve_model_flag, in_atmosphere, mixing_active, &
           new_atmosphere_fit_needed, recompute_surface_bc, &
           use_correct_gradients, want_derivatives, wind_loss_active
      data nao/1/
! load-bearing: see header
      save

      step_status = 0
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
         end_kind_flag = .false.
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

! rewind ISHORT if LRWSH is true (keeps ISHORT small)
          if (lrwsh_placeholder) then
             rewind(short_file_unit)
          endif
! DBG PULSE:  if last model of last run then set LPULSE to LSAVPU
            if (model_iteration.eq.num_models(nk) .and. nk .eq. num_runs) then
                 pulsation_output_active = evo%saved_pulse_output_flag
            end if

! JVS 02/11: Also allow pulse output at particular ages along the way
!
!  If the step is bracketing an age of interest, turn on output. This will
! for the step before and step after the age in AGEOUT. Once the info has
! been printed out, AGEOUT is set to the next age.
!
! Turn on calcad:
      if (acoustic_depth_output) then
            compute_acoustic_depth=.true.
      else
            compute_acoustic_depth = .false.
      endif
! If output has been turned on for a previous step, keep it on for the next
! step, but then turn it off.
      if (acoustic_depth_output) then
            if (ljwrt_placeholder) then
                  print*, 'LJWRT on'
                  pulsation_output_active = evo%saved_pulse_output_flag
                  nao=nao+1
                  lclcd_placeholder =.true.
                  ljlast_placeholder =.false.
                  ljwrt_placeholder=.false.
            else if (.not.ljwrt_placeholder) then
                  lclcd_placeholder=.false.
            endif
! If this is the step before one of the ages of interest, print everything out.
! Also, save model structure.
            if (nao.lt.6) then
                  if (star%run%dage+evo%timestep_yr/1.0D9-ageout_placeholder(nao) .le. 0.0D0 .and. &
                  star%run%dage+2.0D0*evo%timestep_yr/1.0D9-ageout_placeholder(nao) .ge. 0.0D0 .and. .not. ljwrt_placeholder) then
                        print*, 'AGEOUT reached'
                        pulsation_output_active = evo%saved_pulse_output_flag
                        lclcd_placeholder = .true.
                        ljlast_placeholder = .true.
                        ljwrt_placeholder=.true.
                  endif
            endif
       endif


! JVS end
!
! DBG PULSE:  if endage reached then set LPULSE to LSAVPU
! MHP 10/24 GENERALIZE CHECK
         if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (abs(target_end_age(nk)-star%run%dage*1.0D9-evo%timestep_yr) .le. 1.0D0)) then
                 pulsation_output_active = evo%saved_pulse_output_flag
! MHP 7/96 compute sound speed for solar model
                 star%run%sound_speed_output_active = .true.
            end if


!FD echo LSOUND
!        print*,'MAIN LSOUND = ',LSOUND
!FD end
            if (po_output_enabled) then
! MHP 8/25 changed to add file names as declared variables
             call pdist(evo%prev_log_l,evo%prev_log_teff,evo%prev_age,evo%path_length_sq,star%log10_luminosity,star%log_teff,model_iteration,job%pulse_atm_path, &
             job%pulse_env_path,job%pulse_mod_path)
          endif

! STARIN called here for timestep cutting
   15       if (evo%model_diverged_flag) then
!              CALL STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,  ! KC 2025-05-31
             call starin(evo%timestep_yr, evo%delta_time, evo%hydrogen_dt, &
                  evo%trial_sign_flag, evo%ikut_flag, evo%istore_flag, &
                  evo%model_diverged_flag, evo%recompute_envelope_triangle, nk, &
                  evo%dlnrho_dlnp, evo%dlnrho_dlnt, evo%total_angular_momentum, &
                  evo%total_rotational_ke, evo%convective_velocity, &
                  job%mixture_weights, ierr)
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
          evo%punch_pending_flag = .true.

! skip this section if model not to be aged
! MHP 7/98
! need to add logic to permit resacling + time evolution for
! pre-main sequence models
            if (rescale_kind(nk).ne.2 .and. star%model_number.ge.0) then
               evolve_model_flag = .true.
            else if (star%model_number.ge.0 .and. star%log_temperature(1).lt.6.6D0) then
               evolve_model_flag = .true.
            else
               evolve_model_flag = .false.
            endif
            new_atmosphere_fit_needed = .false.
            if (evolve_model_flag) then
! ADD MASS LOSS CALCULATION
               call massloss(star%log10_luminosity,star%run%dage,evo%delta_time,star%composition,star%log_density,star%specific_angular_momentum,star%log_pressure,star%log_radius, &
                             star%log_mass,star%enclosed_mass,star%shell_mass,star%log_total_mass,star%log_temperature,star%envelope_cz_bottom_index,evo%recompute_envelope_triangle, &
                             star%num_zones,star%omega,star%total_mass_msun,star%log_teff,target_envelope_mass,new_atmosphere_fit_needed)
! STORE COMPOSITION MATRIX AT THE BEGINNING OF THE TIMESTEP.
               num_species = 11
               if (use_extended_composition) num_species=15
               do 33 i = 1,star%num_zones
                  do 32 j = 1,num_species
                     star%prev%old_composition(j,i) = star%composition(j,i)
   32             continue
   33          continue
               iteration_level=1
! mixed_zone_bounds_no_overshoot stays an ARGUMENT of mix (not read as
! star% inside it) because crrect passes its own local array there --
! storage deliberately separate from main's. main passes the star copy
! explicitly.
               call mix(evo%delta_time, iteration_level, evo%timestep_yr, &
                    star%core_cz_top_index, star%envelope_cz_bottom_index, &
                    star%mixed_zone_bounds_no_overshoot, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
             evo%timestep_yr = evo%delta_time/seconds_per_year
             star%run%dage = star%run%dage + 1.0D-9*evo%timestep_yr
            endif
!***MHP 1/04 OPACITY TEST
!      IDT = 15
!      DO JJJ = 1,4
!         IDD(JJJ) = 5
!      END DO
!      XXX = 0.7D0
!      ZZZ = 0.02D0
!      Do JJJ = 1,2000
!         READ(75,*)TL,DL,XX,ZZ,OO
!         IF(TT.GT.1.0D9)STOP911
!         TL = LOG10(TT)
!         DL = LOG10(DD)
! DBG 12/95 GET OPACITY
!         CALL GETOPAC(DL, TL,XXX,ZZZ, O, OL, QOD, QOT)
!         DIFF = (O-OO)/O
!         RL = DL - 3*TL +18.0D0
!         WRITE(76,1554)TL,RL,XXX,ZZZ,O,OO,DIFF
! 1554    format(4f11.6,3e20.10)
!      END DO
!*** END TEST
! rezone new model, except rezoning not performed for He flash calculations
          if (.not.helium_flash_active) then
             call hpoint(evo%istore_flag, evo%reset_triangle, evo%h_shell_zone_begin, &
                  evo%has_h_shell, evo%total_angular_momentum, &
                  evo%total_rotational_ke, ierr)
             if (ierr /= 0) return
! STORE NEW CZ BASE
               star%light_burn%jcz = star%envelope_cz_bottom_index
            else
! save old model for PTIME
               do i=1, star%num_zones
                  star%prev%old_pressure(i) = star%log_pressure(i)
                  star%prev%old_temperature(i) = star%log_temperature(i)
                  star%prev%old_radius(i) = star%log_radius(i)
                  star%prev%old_luminosity(i) = star%luminosity_lsun(i)
                  star%prev%old_density(i) = star%log_density(i)
               end do
! JVS 04/14 Save Teff as well
               star%prev%old_teff = star%log_teff
!  JVS 05/25 Added model number to list of saved values
           star%prev%old_num_zones = star%num_zones

          endif
! store starting distribution of rotational kinetic energy.
            if (rotation_active) then
               do i = 1,star%num_zones
                  star%kinetic_energy_rot_old(i) = star%kinetic_energy_rot(i)
               end do
            endif
! changed for lithium burning with overshoot.
! store starting depth of C.Z. for light element burning.
            if (use_extended_composition) then
               star%light_burn%cz_base_radius_prev = 0.0D0
               envelope_cz_zone_prev = star%envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  star%light_burn%pressure_scale_height_start = alphae*exp(clndp*(star%log_pressure(star%envelope_cz_bottom_index)+2.0D0*star%log_radius(star%envelope_cz_bottom_index) &
                           -star%log_density(star%envelope_cz_bottom_index)-cgl-star%log_mass(star%envelope_cz_bottom_index)))
               else
                  star%light_burn%pressure_scale_height_start = 0.0D0
               endif
! find burning rates at the beginning of the time step.
               call lirate88(star%composition,star%log_density,star%log_temperature,star%num_zones,1)
            endif
! begin correction routines
! set flags for CRRECT
! CRRECT checks surface boundary conditions in every iteration
! if LNEW0 = T, new envelope triangle calculated the 1st iteration
! (i.e. old triangle ignored)
! LFINI = T if model has converged
! LARGE = T if model has diverged
          if (lnew0) evo%recompute_envelope_triangle = .true.
            if (.not.evolve_model_flag) evo%delta_time = -dabs(evo%delta_time)
            fcorr = dabs(fcorr0) - fcorri
            iterations_done = 0
            evo%model_diverged_flag = .false.
            converged = .false.
            if (.not.lnews .or. evo%delta_time.le.0.0D0) then
               do 20 i = 1,star%num_zones
! zero entropy terms
                  star%log_temperature_delta(i) = 0.0D0
                  star%log_pressure_delta(i) = 0.0D0
                  star%run%temperature_entropy_term(i) = 0.0D0
                  star%run%pressure_entropy_term(i) = 0.0D0
                  star%run%luminosity_entropy_term(i) = 0.0D0
                  star%run%radius_entropy_term(i) = 0.0D0
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
 20            continue
            else
! use the rate of change in the previous model to estimate the new
! run of structure variables.
               do 30 i = 1,star%num_zones
                  delta_temp_step = star%run%temperature_entropy_term(i)*evo%delta_time
                  delta_pressure_step = star%run%pressure_entropy_term(i)*evo%delta_time
                  delta_lum_step = star%luminosity_lsun(i)*star%run%luminosity_entropy_term(i)*evo%delta_time
                  delta_radius_step = star%run%radius_entropy_term(i)*evo%delta_time
                  star%log_temperature_delta(i) = delta_temp_step
                  star%log_pressure_delta(i) = delta_pressure_step
                  star%log_temperature(i) = star%log_temperature(i) + delta_temp_step
                  star%log_pressure(i) = star%log_pressure(i) + delta_pressure_step
                  star%luminosity_lsun(i) = star%luminosity_lsun(i) + delta_lum_step
                  star%log_radius(i) = star%log_radius(i) + delta_radius_step
! zero gravitational energy terms.
                  star%gravitational_luminosity(i) = 0.0D0
 30            continue
            endif

! FIRST LEVEL OF ITERATIONS
! USE ENVELOPE TRIANGLE OF THE PREVIOUS MODEL;
! FOR THE FIRST MODEL OF A RUN,THE TRIANGLE IS GENERATED HERE.
            max_iterations = niter1
            recompute_surface_bc = .false.
! CALL TO CRRECT - ADDED ITERATION LEVEL
            iteration_level = 1
            call crrect(evo%delta_time, max_iterations, converged, &
                 evo%model_diverged_flag, evo%recompute_envelope_triangle, &
                 evo%reset_triangle, recompute_surface_bc, evo%trial_sign_flag, &
                 evo%istore_flag, in_atmosphere, want_derivatives, &
                 mixing_active, conductive_opacity_flag, evo%dlnrho_dlnt, &
                 evo%dlnrho_dlnp, iterations_done, iteration_level, ierr)
            if (ierr /= 0) return
! SECOND LEVEL OF ITERATIONS
! CHECK ENVELOPE TRIANGLE BEFORE ITERATING FOR SOLUTION
            if (evo%model_diverged_flag) goto 15
            recompute_surface_bc = .true.
            max_iterations = niter2
            iteration_level = 2
            call crrect(evo%delta_time, max_iterations, converged, &
                 evo%model_diverged_flag, evo%recompute_envelope_triangle, &
                 evo%reset_triangle, recompute_surface_bc, evo%trial_sign_flag, &
                 evo%istore_flag, in_atmosphere, want_derivatives, &
                 mixing_active, conductive_opacity_flag, evo%dlnrho_dlnt, &
                 evo%dlnrho_dlnp, iterations_done, iteration_level, ierr)
            if (ierr /= 0) return
            if (evo%model_diverged_flag) goto 15
! 7/91 STORE CHANGES IN THE STRUCTURE. THESE CHANGES ARE USED TO GET AN
! IMPROVED FIRST GUESS AT THE STRUCTURE FOR THE NEXT MODEL IF LNEWS=T.
            if (evo%delta_time.gt.0.0D0) then
               do 27 ii = 1,star%num_zones
                  star%run%temperature_entropy_term(ii)=star%log_temperature_delta(ii)/evo%delta_time
                  star%run%pressure_entropy_term(ii)=star%log_pressure_delta(ii)/evo%delta_time
                  star%run%luminosity_entropy_term(ii)=2.0D0*(star%luminosity_lsun(ii)-star%prev%old_luminosity(ii))/(star%luminosity_lsun(ii)+star%prev%old_luminosity(ii))/evo%delta_time
                  star%run%radius_entropy_term(ii)=(star%log_radius(ii)-star%prev%old_radius(ii))/evo%delta_time
 27            continue
            endif
! THIRD LEVEL OF ITERATIONS
            recompute_surface_bc = .false.
            max_iterations = niter3
            iteration_level = 3
            call crrect(evo%delta_time, max_iterations, converged, &
                 evo%model_diverged_flag, evo%recompute_envelope_triangle, &
                 evo%reset_triangle, recompute_surface_bc, evo%trial_sign_flag, &
                 evo%istore_flag, in_atmosphere, want_derivatives, &
                 mixing_active, conductive_opacity_flag, evo%dlnrho_dlnt, &
                 evo%dlnrho_dlnp, iterations_done, iteration_level, ierr)
            if (ierr /= 0) return
            if (evo%model_diverged_flag) goto 15
            if (.not.rotation_active) then
               itdif1 = 1
            endif
! MHP 05/02
! IF THE CODE IS ITERATING BETWEEN THE STRUCTURE AND ROTATION
! SOLUTIONS, ENSURE THAT THE START-OF-TIMESTEP QUANTITIES
! HCOMPP (COMPOSITION) AND HJMSAV (ANGULAR MOMENTUM) ARE ONLY
! OVERWRITTEN ON THE LAST RUN THROUGH.
            if (itdif1.gt.1) then
               do i = 1,star%num_zones
                  star%run%orig_specific_angular_momentum(i) = star%specific_angular_momentum(i)
                  do j = 1,15
                     star%run%orig_composition(j,i) = star%prev%old_composition(j,i)
                  end do
               end do
            endif
            do itrot = 1, itdif1
! MHP 05/02 RESTORE ORIGINAL "START OF TIMESTEP"
! VALUES FOR THE COMPOSITION MATRIX
               if (itrot.gt.1) then
                  do i = 1,star%num_zones
                     do j = 1,15
                        star%prev%old_composition(j,i) = star%run%orig_composition(j,i)
                     end do
                  end do
               endif
! 7/91 THE FOURTH LEVEL OF ITERATION REPEATS THE ITERATION BETWEEN THE
! MIXING AND THE STRUCTURE VARIABLES.  IT SHOULD NOT BE USED FOR MODELS
! WHERE SEMI-CONVECTION IS IMPORTANT (ITERATING BETWEEN THE BURNING AND
! STRUCTURE GENERATES OSCILLATIONS). IT SHOULD BE USED FOR HIGH-PRECISION
! WORK (E.G. SOLAR MODELS).
! Surface boundary conditions checked again since we've changed the
! star%composition (and hence the structure) of the model in ITLVL=3
! (to be implemented when I know the rest of it works!)
            max_iterations = niter4
            recompute_surface_bc=.false.
            iteration_level = 4
            call crrect(evo%delta_time, max_iterations, converged, &
                 evo%model_diverged_flag, evo%recompute_envelope_triangle, &
                 evo%reset_triangle, recompute_surface_bc, evo%trial_sign_flag, &
                 evo%istore_flag, in_atmosphere, want_derivatives, &
                 mixing_active, conductive_opacity_flag, evo%dlnrho_dlnt, &
                 evo%dlnrho_dlnp, iterations_done, iteration_level, ierr)
            if (ierr /= 0) return
!  25         CONTINUE
            if (.not.converged) then
! MODEL FAILED TO CONVERGE WITHIN(NITER1+NITER2+NITER3+NITER4)ITERATIONS
               evo%model_diverged_flag = .true.
               goto 15
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
      call mixcz(star%composition,star%shell_mass,star%convective_flag,star%num_zones)
! G Somers END

! MHP 9/94 STORE TOTAL AGE IN SAGE
            disk_lifetime = star%run%dage
            if (rotation_active) then
! RESTORE ORIGINAL START OF TIMESTEP VALUES
! TO THE ANGULAR MOMENTUM DISTRIBUTION
               if (itrot.gt.1) then
                  do i = 1,star%num_zones
                     star%specific_angular_momentum(i) = star%run%orig_specific_angular_momentum(i)
                  end do
               endif
! MHP 9/94 ADDED FLAG TO TURN ON ROTATION OUTPUT WHEN END OF KIND
! CARD REACHED.
! MHP 10/24 GENERALIZE CHECK
         if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (abs(target_end_age(nk)-star%run%dage*1.0D9-evo%timestep_yr) .le. 1.0D0)) then
!               IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)THEN
                  star%run%lprt0_placeholder = .true.
               else
                  star%run%lprt0_placeholder = .false.
               endif
! FIND THE NEW RUN OF OMEGA
! JENV0 ADDED TO SR CALL.
               wind_loss_active = ljdot0
               call getw(evo%delta_time, evo%max_domega_frac, wind_loss_active, &
                    envelope_cz_zone_prev, jerr)
               if (jerr /= 0) then
               ! 2026 (phase five, step B): propagate instead of stopping
                  ierr = jerr
                  return
               end if
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
               call fpft(star%log_density,star%log_radius,star%log_mass,star%num_zones,star%omega,star%eta_squared,star%pressure_rotation_factor,star%temperature_rotation_factor,star%mean_gravity,star%mean_radius)
            endif
            end do
! LOCATE THE HYDROGEN-BURNING SHELL AND THE BOUNDARIES OF THE CENTRAL
! AND SURFACE CONVECTION ZONES (IF APPLICABLE).
       call findsh(star%composition,star%luminosity_lsun,star%convective_flag,star%num_zones, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,evo%h_shell_zone_begin,evo%h_shell_end_index,evo%h_shell_midpoint_zone, &
              evo%has_h_shell)
! PERFORM LIGHT ELEMENT BURNING
         if (use_extended_composition .and. star%model_number.ge.0 .and. evo%delta_time.gt.0.0D0) then
! ONLY FOR MODELS WITHOUT ROTATION, OR WITHOUT ROTATIONAL MIXING.
            if (.not.rotation_active .or. .not.instability_transport_active) then
! FIND CONVECTION ZONE DEPTH AT THE END OF THE TIME STEP.
               call convec(star%composition,star%log_density,star%log_pressure,star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,star%num_zones,star%radiative_zone_bounds,star%mixed_zone_bounds, &
                            star%mixed_zone_bounds_no_overshoot,star%core_cz_top_index,star%envelope_cz_bottom_index,num_radiative_zones,num_mixed_zones,num_mixed_zones_no_overshoot)
! CHANGED FOR LITHIUM BURNING WITH OVERSHOOT.
               envelope_cz_zone_end = star%envelope_cz_bottom_index
               if (envelope_overshoot_active) then
                  star%light_burn%pressure_scale_height_end = alphae*exp(clndp*(star%log_pressure(star%envelope_cz_bottom_index)+2.0D0*star%log_radius(star%envelope_cz_bottom_index) &
                           -star%log_density(star%envelope_cz_bottom_index)-cgl-star%log_mass(star%envelope_cz_bottom_index)))
               else
                  star%light_burn%pressure_scale_height_end = 0.0D0
               endif
! FIND BURNING RATES AT THE END OF THE TIME STEP.
               call lirate88(star%composition,star%log_density,star%log_temperature,star%num_zones,2)
!                CALL LIBURN(DELTS,HCOMP,HD,HR,HS1,HS2,HT,JENV1,JENV0,M)  ! KC 2025-05-31
               call liburn(evo%delta_time,star%composition,star%log_radius,star%enclosed_mass,star%shell_mass,star%log_temperature,envelope_cz_zone_end,envelope_cz_zone_prev,star%num_zones)
            endif
         endif
! MHP 07/02 RESTORE PRIOR FITTING POINT IF MASS ACCRETION IS BEING
! INCLUDED
         if (new_atmosphere_fit_needed) then
            call getnewenv(target_envelope_mass,star%composition,star%log_density,star%luminosity_lsun,star%log_pressure,star%log_radius,star%log_mass,star%enclosed_mass,star%shell_mass, &
!     *                     HSTOT,HT,LC,ETA2,HG,HI,HJM,QIW,R0,  ! KC 2025-05-31
                            star%log_total_mass,star%log_temperature,star%convective_flag,star%eta_squared,star%moment_of_inertia,star%specific_angular_momentum,star%qiw,star%mean_radius, &
                            star%kinetic_energy_rot,star%log10_luminosity,evo%total_angular_momentum,evo%total_rotational_ke,star%log_teff,star%num_zones,evo%recompute_envelope_triangle)
! CALCULATE FP AND FT GIVEN OMEGA FOR THE NEW POINT DISTRIBUTION
            call fpft(star%log_density,star%log_radius,star%log_mass,star%num_zones,star%omega,star%eta_squared,star%pressure_rotation_factor,star%temperature_rotation_factor,star%mean_gravity,star%mean_radius)
            new_atmosphere_fit_needed = .false.
         endif
! DETERMINE TIMESTEP FOR NEXT MODEL
! HTIMER ALSO LOCATES THE H-BURNING SHELL
! JVS 04/14 added teffl to passed htimer variables
       evo%delta_time = dabs(evo%delta_time)
       evo%delta_time_saved = evo%delta_time
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *        JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(evo%delta_time,evo%hydrogen_dt,star%num_zones,star%log_density,star%luminosity_lsun,star%enclosed_mass,star%shell_mass,star%log_temperature,star%composition,star%core_cz_top_index, &
              evo%h_shell_midpoint_zone,star%luminosity_breakdown,star%run%dage,evo%timestep_yr,nk,star%log_pressure,star%log_radius,star%omega, &
              evo%max_domega_frac,evo%h_shell_zone_begin,star%log_teff)
! IF EVOLVING TO A GIVEN AGE AND KIND CARD IS DONE, AVOID ZEROING OUT
! TIMESTEP WRITTEN TO MODEL (AS THIS MAKES CONTINUING A SEQUENCE AWKWARD.)
!     INSTEAD WRITE THE PREVIOUS MODEL TIMESTEP TO MODEL.
! ONLY IF A FIXED END AGE IS USED, NOT FOR OTHER STOPS
       if (end_age_stop_active(nk) .and. target_end_age(nk).gt.0.0D0) then
          if (target_end_age(nk)-star%run%dage*1.0D9.le.1.0D0) then
             evo%delta_time = max(evo%delta_time_saved,1.0D-3*star%run%dage*seconds_per_year)
             evo%timestep_yr = evo%delta_time/seconds_per_year
          else
             evo%delta_time_saved = evo%delta_time
          endif
       else
          evo%delta_time_saved = evo%delta_time
       endif
       if (rescale_kind(nk).ne.2) star%model_number = star%model_number+1
! 2026 (phase four, step 5): compute the output diagnostics in the
! star layer (fills star%run%*, star%luminosity_breakdown
! renormalization, star%turnover% via gettau); wrtout below only reads.
       call update_output_diagnostics(ierr)
       if (ierr /= 0) return
! WRTOUT IS THE OUTPUT DRIVER ROUTINE
       call wrtout(evo%timestep_yr, log_gravity, evo%has_h_shell, &
            evo%h_shell_zone_begin, evo%h_shell_midpoint_zone, evo%h_shell_end_index, &
            evo%trial_sign_flag, evo%punch_pending_flag, evo%total_angular_momentum, &
            evo%total_rotational_ke)

! MHP 10/24 GENERALIZED STOP CONDITIONS
!     IF EVOLVING TO A GIVEN AGE AND AGE IS REACHED, KIND CARD IS DONE
!       IF(LENDAG(NK).AND.ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)GOTO 110
       if (end_age_stop_active(nk).and.target_end_age(nk).gt.0.0D0 .and. &
         (target_end_age(nk)-star%run%dage*1.0D9).le.1.0D0) goto 810
! MHP 10/24 CHECK ALL STOP CONDITIONS, EXIT IF ANY SATISFIED
         end_kind_flag = .false.
         if (end_age_stop_active(nk).and.central_deuterium_stop(nk).gt.0.0D0 .and. &
              star%composition(12,1).lt.central_deuterium_stop(nk)) then
            write(*,104)star%composition(12,1),central_deuterium_stop(nk)
 104        format('CENTRAL D ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
         if (end_age_stop_active(nk).and.central_hydrogen_stop(nk).gt.0.0D0 .and. &
              star%composition(1,1).lt.central_hydrogen_stop(nk)) then
            write(*,105)star%composition(1,1),central_hydrogen_stop(nk)
 105        format('CENTRAL X ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
         if (end_age_stop_active(nk).and.central_helium_stop(nk).gt.0.0D0 .and. &
              star%composition(2,1).lt.central_helium_stop(nk)) then
            write(*,106)star%composition(2,1),central_helium_stop(nk)
 106        format('CENTRAL Y ',E12.4,' BELOW STOP VALUE ',E12.4)
            end_kind_flag =.true.
         endif
! IF EXITING, SET I/O FLAGS PROPERLY AND EXIT LOOP
         if (end_kind_flag) then
            pulsation_output_active = evo%saved_pulse_output_flag
            star%run%sound_speed_output_active = .true.
            star%run%lprt0_placeholder = .true.
            goto 810
         endif
! TEST IF MODEL IS NEAR DESIRED Teff AND L. IF NOT RESCALE AND TRY AGAIN.
         if (calibrate_star_flag .and. .not. star_found_flag) then
            if (mod(nk,2).eq.0) then
             if (model_iteration.eq.1) then
                teff_kelvin_unused = 10.0D0**star%log_teff
             else
                call chkscal(star%log10_luminosity, star%log_teff, star%run%dage, nk)
                if (just_passed_target_radius_flag) goto 820
             end if
          endif
       endif

! END OF RUN

      step_status = 0
      return

! run finished (was: goto 110)
 810  step_status = 1
      return
! leave the run loop (was: goto 200)
 820  step_status = 2
      return
end subroutine evolve_step

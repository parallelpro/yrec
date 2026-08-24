!----------------------------------------------------------------------
! run_yrec
!----------------------------------------------------------------------
! Added 2026 (phase five, step A -- the embeddable engine; see
! ROADMAP.md "Next milestone"). This is the entire body of the former
! `program main`, moved verbatim: namelist/controls read (parmin),
! table setup (setups), and the Monte-Carlo/run loop containing the
! full evolution driver. `program main` is now a thin CLI wrapper.
!
! The blanket `save` below is load-bearing: program-unit variables are
! implicitly static, and this body's locals (including its `data`
! initializations and rescale bookkeeping) rely on that; as subroutine
! locals they would otherwise be automatic. With save they are exactly
! as static as before -- which also means run_yrec is NOT yet
! re-entrant (one call per process); phase C of this milestone (the
! SAVE/reset audit) is what will make repeated calls legal.
!
! ierr: reserved for phase B, when the driver-side stops in this body
! (and the `if (jerr /= 0) stop` seams from phase-three stage 3)
! become error returns. In phase A it is set to 0 and the historical
! stops remain exactly where they were.
!
subroutine run_yrec(ierr)

! the array size, i.e. max # of shells is specified in the parameter
! statement. it defines JSON. to change the array size do a global
! change on "JSON=2000" or whatever.
      use net_lib
      use star_info_lib, only: star, i_h1, i_he4, i_lum_grav, i_metals, json
      use yrec_output, only: output_run_header
      use luout_lib
      use const_lib
      use yrec_reset_lib, only: yrec_run_prologue
      use stop_conditions, only: step_kind_card_done, &
           step_leave_run_loop, disarm_satisfied_stops
      implicit none
      integer :: step_status
      integer, parameter :: numtt = 70
      integer, parameter :: numd = 19
      integer, parameter :: numx = 10
      integer, parameter :: numz = 13
      integer, parameter :: numxz = 126

! DBGLAOL - to save space make tables single precision
! MHP 8/25 Removed unused variables and added pass-through variables
! See the CROSS-CALLEE NAMING NOTE above for the file-path locals
! (declared further below, in the --- locals --- section, using
! setups.f90's/pdist.f90's own descriptive dummy-argument spelling).

!     MHP 10/24 FLAG FOR END OF RUN

! MHP 10/24 NEW VARIABLES FOR STOP CRITERIA ON CENTRAL ABUNDANCE are
! carried in common/sett/ above.

! latest values (Bahcall and Pinsonneault 1996)-actual values set in
! subroutine PARMIN
      double precision :: bp96_scale_factor(17)
      data bp96_scale_factor/0.9558,0.9690,0.9712,1.0,1.0,0.992,1.0,1.0, &
           1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.92088,0.1625/
! MHP 3/96 added data for base solar age, L
!       DATA SUNAGE,SUNL/4.57D09,3.844D33/  ! KC 2025-05-31
      double precision :: reference_solar_luminosity
      data reference_solar_luminosity/3.844D33/

! --- locals ---
! calibration card protocols (see the verdict block in the run loop)
      integer, parameter :: solar_calib_cards_per_cycle = 3
      integer, parameter :: star_calib_cards_per_cycle = 2
      logical :: runs_complete
      integer :: monte_carlo_run_number
      double precision :: age_scale_factor
      integer :: convergence_iterations
      double precision :: initial_x_guess, initial_alpha_guess
      logical :: saved_use_structure_dt_limits
      integer :: saved_atm_choice
      integer:: i, j
      integer :: wrtlst_unit, model_iteration
      double precision :: log_r_rsun, current_zx, surface_z_over_x
      character(len=256) :: alex06_table_path, allard_table_path, &
           atm_table_path, fermi_table_path, kurucz_table_path, &
           kurucz_table2_path, laol_table_path, laol_table2_path, &
           opal95_table_path, opal92_table_path
      character(len=256) :: zams_a_table_path, zams_b_table_path, &
           zams_c_table_path, centre1_table_path, centre2_table_path, &
           centre3_table_path, centre4_table_path, centre5_table_path
      character(len=256) :: opal92_table2_path, pure_z_table_path, &
           scv_h_table_path, scv_he_table_path, scv_z_table_path
      character(len=256) :: alex95_table_paths(7)
      double precision :: monte_helium_diffusion_fraction

      integer, intent(out) :: ierr
! load-bearing: see header
      save
   ! INTENTIONAL: run-level state across MC runs; reset via yrec_reset
!*******
! START
!*******
      ierr = 0

! 2026 (phase five, step C): fresh-process semantics for repeated
! calls -- see core/yrec_reset.f90.
      call yrec_run_prologue

      call setversion()

      iowr = 9
! LPUNCH is TRUE once first model is calculated
      star%evo%punch_pending_flag = .false.
! 2026 (phase five): controls read and setup are now star-layer
! routines operating on the star_job structure (state/star_job_lib).
      call read_controls(ierr)
      if (ierr /= 0) return
      call star_setup(ierr)
      if (ierr /= 0) return

      do monte_carlo_run_number = star%job%mc_run_start,star%job%mc_run_end
! for monte carlo run, input values of parameters being changed.
      if (lmonte) then
         cross_section_scale(1) = star%run%s11_rate(monte_carlo_run_number)*bp96_scale_factor(1)
         cross_section_scale(2) = star%run%s33_rate(monte_carlo_run_number)*bp96_scale_factor(2)
         cross_section_scale(3) = star%run%s34_rate(monte_carlo_run_number)*bp96_scale_factor(3)
         cross_section_scale(16) = star%run%s17_rate(monte_carlo_run_number)*bp96_scale_factor(16)
         monte_helium_diffusion_fraction = star%run%helium_fraction_param(monte_carlo_run_number)
         fgrz = star%run%diffusion_factor(monte_carlo_run_number)
         solar_luminosity_cgs = reference_solar_luminosity*star%run%luminosity_target(monte_carlo_run_number)
         log10_solar_luminosity = dlog10(solar_luminosity_cgs)
         ln_solar_luminosity = ln10/solar_luminosity_cgs
         age_scale_factor = star%run%age_target(monte_carlo_run_number)
! timestep and final age are altered in SR SETCAL; input #s should be
! scaled for a solar age of 4.57 Gyr
         target_end_age(2)=1.0D8
         target_end_age(3)=4.57D9
      else
         age_scale_factor = 1.0D0
      endif
! DBG PULSE: save LPULSE flag, set LPULSE to F except on last model of
! last run, then set LPULSE to saved value of LPULSE.
      star%evo%saved_pulse_output_flag = pulsation_output_active
! 02/11 JVS uncommented LPULSE=.FALSE.
!      LPULSE = .FALSE.
! MHP 1/93 add option to automatically calibrate solar model.
! MHP 3/96 added counter for # of iterations per converged model and
! starting estimate of ALPHA and X
      if (calibrate_solar_model) then
         call setcal(age_scale_factor)
         convergence_iterations = 1
         initial_x_guess = rescale_params(2,1)
         initial_alpha_guess = mixing_length_array(1)
         saved_use_structure_dt_limits = use_structure_dt_limits   ! save LPTIME for reuse during calibration
         saved_atm_choice  = atm_choice    ! save KTTAU for reuse during calibration
      else
         convergence_iterations = 0
      endif
! DBG 12/94 add option to automatically calculate a stellar model
! of specified Teff and L
      if (calibrate_star_flag) then
         call setscal
      endif

!**********
!     RUN THROUGH THE KIND CARDS IN ORDER
!**********
      run_loop: do nk = 1, num_runs
! Per-kind-card setup (2026, core/ phase 3): everything from the
! card's initial composition through the first timestep estimate is
! begin_kind_card below; config/read errors return through ierr.
         call begin_kind_card
         if (ierr /= 0) return

! for a given kind card, evolve NMODLS(NK) times
! if rescaling is being performed, NMODLS(NK) is the number of times
! the new model is being relaxed
       do model_iteration = 1,num_models(nk)
! 2026 (phase five): one model advance per iteration, extracted to
! core/evolve_step.f90 (see its header for the step_status contract).
       call evolve_step(model_iteration, step_status, ierr)
       if (ierr /= 0) return
       if (step_status == step_kind_card_done) exit
       if (step_status == step_leave_run_loop) cycle run_loop
       end do

! Store the card's last model if requested (2026: end_kind_card).
         call end_kind_card

! --- End-of-card calibration verdict ------------------------------
! The two calibration modes drive the run list with implicit card
! protocols, stated here once:
!
!  SOLAR calibration (calibrate_solar_model; setcal/chkcal): kind
!  cards run in TRIPLES of solar_calib_cards_per_cycle = 3 --
!  card 3k+1 rescales the seed model to the current (X, Z, alpha)
!  guess, card 3k+2 evolves it to 1e8 yr (settling), card 3k+3
!  evolves to the target solar age. After each completed triple
!  (mod(nk,3) == 0) chkcal tests (log L, log R[, log Z/X]) against
!  the Sun and either declares convergence or writes the
!  Newton-corrected (X, Z, alpha) into the NEXT triple's cards.
!  setcal pre-expands the run list to 16 such triples (48 cards);
!  the verdict below also caps the attempts at 15.
!
!  STAR calibration (calibrate_star_flag; setscal/chkscal): cards
!  run in PAIRS of star_calib_cards_per_cycle = 2 -- odd cards
!  rescale, even cards evolve. chkscal watches every model of an
!  even card for the target-radius crossing (that check lives in
!  evolve_step; a crossing leaves the run loop via
!  step_leave_run_loop) and, once the luminosity there matches too,
!  arms a final run stopped at the interpolated age
!  (star_found_flag) -- which the verdict below turns into the end
!  of the run list.
         call end_of_card_calibration(runs_complete)
         if (runs_complete) exit run_loop

! END RUN LOOP
      end do run_loop
! EXIT RUN LOOP

! FOR MONTE CARLO, REWIND OUTPUT FILES AND WRITE OUT SNU FLUXES AND
! MODEL PARAMETERS (legacy mode only; io/write_run_summaries.f90 --
! 2026, core/ phase 4: the last writer left in the driver moved to
! io/. surface_z_over_x is inout: the failed-convergence branch
! reports the value carried from a previous cycle, historical
! SAVE semantics preserved).
      call write_run_summaries(monte_carlo_run_number, &
           convergence_iterations, initial_x_guess, initial_alpha_guess, &
           log_r_rsun, surface_z_over_x)
      end do

! 2026 (phase five, step B): the normal end-of-job stop became this
! clean return (ierr stays 0); the CLI wrapper simply ends.
      return

contains

! ---------------------------------------------------------------
! End-of-kind-card calibration verdict (see the protocol comment
! at the call site). Sets runs_complete when the run list is done:
! solar calibration converged (or 15 attempts exhausted), or the
! star-calibration target has been hit on an even card. On a
! non-final solar triple, applies the between-cycle bookkeeping
! (iteration counter, pulse-file rewinds, LPTIME/KTTAU restores),
! ordering preserved exactly.
subroutine end_of_card_calibration(runs_complete)
      logical, intent(out) :: runs_complete

      runs_complete = .false.
! MHP 1/93 CHECK AUTOMATIC CALIBRATATION OF SOLAR MODEL.
!c MHP 5/96 changed solar calibration to perform solar models in 3 kind cards
         if (calibrate_solar_model) then
! JVS Turn off calcad - speeds things up
            compute_acoustic_depth=.false.
            if (mod(nk,solar_calib_cards_per_cycle).eq.0) then
               log_r_rsun = 0.5D0*(star%log_L+log10_solar_luminosity-c4pil-csigl-4.0D0*star%log_Teff)-log10_solar_radius
! MHP 06/13 Add solar Z/X to observables
               current_zx = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
               call chkcal(star%log_L,log_r_rsun,nk,current_zx)
!               CALL CHKCAL(BL,RLL,NK)
               use_structure_dt_limits = saved_use_structure_dt_limits  ! Restore LPTIME to original value for next cycle
               atm_choice  = saved_atm_choice    ! Restore KTTAU to original value for next cycle
               if (star%run%solar_calibration_active) then
                  runs_complete = .true.
                  return
               else
!c MHP 8/96 added counter for # of runs needed for calibration
                  convergence_iterations = convergence_iterations + 1
! MHP 6/97 STOP AFTER 10 ATTEMPTS AT CALIBRATION
!                  IF(ICONV.GE.11) GOTO 250
                  if (convergence_iterations.ge.15) then
                     runs_complete = .true.
                     return
                  end if
                  if (pulsation_output_active) then
! DBG 6/93 Need to delete pulse output because have not got ultimate
! model yet.
! MHP 8/25 Replaced delete file with rewind file. This is functionally the same and avoids the need to pass the character string for the file name from parmin.
                     rewind(opal_model_unit)
                     rewind(opal_envelope_unit)
                     rewind(opal_atm_unit)
!                     CLOSE(IOPMOD, STATUS='DELETE')
!                     CLOSE(IOPENV, STATUS='DELETE')
!                     CLOSE(IOPATM, STATUS='DELETE')
!                     OPEN(IOPMOD, FILE=FPMOD,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
!                     OPEN(IOPENV, FILE=FPENV,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
!                     OPEN(IOPATM, FILE=FPATM,STATUS='UNKNOWN',
!     *                    FORM='FORMATTED')
                  end if
               end if
            endif
         endif

! DBG 12/94 NO MORE RUNS NEEDED. HAVE CALIBRATED STELLAR MODEL
         if (calibrate_star_flag .and. star_found_flag.and.(mod(nk,star_calib_cards_per_cycle).eq.0)) then
            runs_complete = .true.
            return
         end if
end subroutine end_of_card_calibration


! ---------------------------------------------------------------
! Start one kind card: per-card composition/mixing-length settings,
! read or reuse the starting model (starin, incl. rescaling), check
! the rotation configuration, disarm already-satisfied stops, the
! optional neutrino table, output headers, the He-flash restore
! copy, locate shells/CZ edges (findsh), the first timestep
! estimate (htimer), and zero the entropy and light-element-rate
! state. Sets ierr on configuration or model-read errors.
subroutine begin_kind_card
         star%run%sound_speed_output_active = .false.
!         LPULSE=.FALSE.
         initial_envelope_x = initial_x_array(nk)
         initial_envelope_z = initial_z_array(nk)
         cmixl = mixing_length_array(nk)
       change_envelope_mass_flag = has_senv0_array(nk)
       requested_envelope_mass = senv0_array(nk)
       star%evo%reset_triangle = .false.
       star%evo%model_diverged_flag = .false.
! MHP 10/02 ZERO OUT INITIAL ANGULAR MOMENTUM
         star%evo%total_angular_momentum = 0.0D0
         star%evo%total_rotational_ke = 0.0D0
! read in the initial model here
! STARIN also calls RSCALE to perform rescaling if requested
!        CALL STARIN(BL,CFENV,DAGE,DDAGE,DELTS,DELTSH,DELTS0,ETA2,  ! KC 2025-05-31
       call starin(star%evo%timestep_yr, star%evo%dt, star%evo%hydrogen_dt, star%evo%trial_sign_flag, &
            star%evo%ikut_flag, star%evo%istore_flag, star%evo%model_diverged_flag, &
            star%evo%recompute_envelope_triangle, nk, star%evo%dlnrho_dlnp, star%evo%dlnrho_dlnt, &
            star%evo%total_angular_momentum, star%evo%total_rotational_ke, &
            star%evo%convective_velocity, star%job%mixture_weights, ierr)
       if (ierr /= 0) return

      if ((star%omega(1) .eq. 0) .and. (rotation_active)) then

1611      format('LROT set to TRUE, but OMEGA(1) = 0. Stopping.', &
                 ' Initialize rotation rates or set LROT to', &
                 ' FALSE.')
          print 1611
          ! 2026 (phase five, step B): configuration error returns to the
          ! CLI wrapper (which stops) instead of stopping here.
          ierr = 1
          return
      endif
!     MHP 10/24 CHECK STOP CONDITIONS AND DISABLE THEM IF THE STARTING VALUES ARE BELOW THE TARGET THRESHOLD
! (2026: one table walk in stop_conditions -- the hand-written D/X/Y
! triple that used to live here carried the disarm-the-wrong-stop bug
! fixed in phase 1.)
         call disarm_satisfied_stops(nk)
! Opt-in diagnostic (2026): the former LNUTAB per-zone neutrino
! table, off since 2004, is now the compute_neutrino_fluxes control
! (core/neutrino_flux_table.f90); it describes the starting model
! of each kind card.
      if (compute_neutrino_fluxes) then
         call neutrino_flux_table
      endif
! save mass in solar units
         pulsation_mass_msun=star%star_mass
! MHP 08/02 STORE STARTING CZ PROPERTIES
         star%light_burn%jcz = star%envelope_cz_bottom_index
         star%turnover%convective_turnover_timescale = 0.0D0
! write out headers of the appropriate output files
      call output_run_header(star%star_mass)
! DBG PULSE OUT 7/92
! initialize variables for calculating when to dump pulse output
         star%evo%prev_log_l = star%log_L
         star%evo%prev_log_teff = star%log_Teff
         star%evo%prev_age = star%run%dage
         star%evo%path_length_sq = 0.0D0

       if (helium_flash_active) then
! timestep cutting requires a model stored in logical unit ILAST
! or it will crash - so copy initial model to unit ILAST
          if (star%evo%punch_pending_flag) then
             wrtlst_unit = ilast
             call wrtlst(wrtlst_unit,star%xa,star%logRho,star%luminosity_lsun, &
                  star%logP,star%logR,star%log_mass,star%logT,star%convective_flag, &
                  star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure, &
                  star%fit_point_temperature,star%fit_point_radius,star%envelope_fit_coeffs, &
                  star%evo%trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index, &
                  star%envelope_cz_bottom_index,star%model_number,star%nz, &
                  star%star_mass,star%log_Teff,star%log_L,star%log_total_mass,star%run%dage, &
                  star%evo%timestep_yr,star%omega)
          endif
       endif

! locate the hydrogen-burning shell and the boundaries of the central
! and surface convection zones (if applicable).
         call findsh(star%xa,star%luminosity_lsun,star%convective_flag,star%nz, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%evo%h_shell_zone_begin, &
              star%evo%h_shell_end_index,star%evo%h_shell_midpoint_zone,star%evo%has_h_shell)
! determine timestep for model
! JVS 04/14 Added Teffl to passed variables
!        CALL HTIMER(DELTS,DELTSH,M,HD,HL,HS1,HS2,HT,LC,HCOMP,JCORE,
!      *               JXMID,TLUMX,DAGE,DDAGE,QDT,QDP,NK,HP,HR,OMEGA,  ! KC 2025-05-31
       call htimer(star%evo%dt,star%evo%hydrogen_dt,star%nz,star%logRho,star%luminosity_lsun, &
            star%m,star%dm,star%logT,star%xa,star%core_cz_top_index, &
            star%evo%h_shell_midpoint_zone,star%luminosity_breakdown,star%run%dage,star%evo%timestep_yr,nk, &
            star%logP,star%logR,star%omega,star%evo%max_domega_frac,star%evo%h_shell_zone_begin, &
            star%log_Teff)

       star%evo%dt_saved = star%evo%dt
! zero out entropy terms.
         do i = 1,star%nz
            star%run%temperature_entropy_term(i) = 0.0D0
            star%run%pressure_entropy_term(i) = 0.0D0
            star%run%luminosity_entropy_term(i) = 0.0D0
            star%run%radius_entropy_term(i) = 0.0D0
         end do

! zero out light element burning rates in the surface CZ.
         if (use_extended_composition) then
            star%light_burn%log_rate_li6_prev = 0.0D0
            star%light_burn%log_rate_li7_prev = 0.0D0
            star%light_burn%log_rate_be9_prev = 0.0D0
         endif
end subroutine begin_kind_card

! ---------------------------------------------------------------
! Finish one kind card: store the last model to the .store stream
! when configured (putstore) and clear the punch flag.
subroutine end_kind_card
! G Somers 11/14, CHANGE CALL TO PUTSTORE INSTEAD OF WRTLST.
! STORE LAST MODEL IN ISTOR IF LSTORE, LSTPCH, AND LPUNCH ARE .TRUE.
         if (lstore.and.lstpch.and.star%evo%punch_pending_flag) then
          call putstore(star%xa,star%logRho,star%luminosity_lsun,star%logP,star%logR,star%log_mass,star%logT,star%convective_flag,star%trial_log_temperature,star%trial_log_luminosity,star%fit_point_pressure,star%fit_point_temperature,star%fit_point_radius, &
                 star%envelope_fit_coeffs,star%evo%trial_sign_flag,star%luminosity_breakdown,star%core_cz_top_index,star%envelope_cz_bottom_index,star%model_number,star%nz,star%star_mass,star%log_Teff,star%log_L,star%log_total_mass, &
                 star%run%dage,star%evo%timestep_yr,star%omega,star%m,star%eta_squared,star%mean_radius,star%fp_rot,star%ft_rot,star%j_rot,star%i_rot)
            star%evo%punch_pending_flag = .false.
       endif
! 110  CONTINUE
! G Somers END
end subroutine end_kind_card

end subroutine run_yrec

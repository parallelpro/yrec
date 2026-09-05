!----------------------------------------------------------------------
! run_yrec
!----------------------------------------------------------------------
! Added 2026 (phase five -- the embeddable engine; see ROADMAP.md).
! The body of the former `program main`: namelist/controls read
! (read_parameters), table setup (setup_tables), and the Monte-Carlo /
! kind-card run loop around evolve_step. `program main` is a thin CLI
! wrapper; yrec_capi is the C entry.
!
! The blanket `save` below is load-bearing: program-unit variables are
! implicitly static, and this body's locals (the rescale and
! calibration bookkeeping) rely on that; as subroutine locals they
! would otherwise be automatic. Repeated calls are made legal by
! yrec_run_prologue (core/yrec_reset.f90), which puts the saved state
! back to process-start values.
!
! ierr: 0 on a normal end of job; positive on a configuration, read or
! numerical error (the CLI wrapper exits 1); negative for
! numerics_termination (the historical "diverged" stop, exit 0).
!
subroutine run_yrec(ierr)

      use net_lib
      use star_info_lib, only: star, i_h1, i_metals
      use yrec_output, only: output_run_header
      use luout_lib
      use phys_const_lib
      use yrec_reset_lib, only: yrec_run_prologue
      use monte_carlo_lib, only: apply_monte_carlo_parameters, &
           write_run_summaries
      use stop_conditions, only: step_kind_card_done, &
           step_leave_run_loop, init_stop_conditions, &
           check_rotation_initialised
      use run_log_lib, only: log_run_summary
      use controls_check_lib, only: warn_inconsistent_controls
      use observables_lib, only: log_r_surface_cm
      implicit none
      integer :: step_status

! --- locals ---
! calibration card protocols (see the verdict block in the run loop)
      integer, parameter :: solar_calib_cards_per_cycle = 3
      integer, parameter :: star_calib_cards_per_cycle = 2
      logical :: runs_complete
      integer :: monte_carlo_run_number
      double precision :: age_scale_factor
      integer :: convergence_iterations
! end-of-run wall-clock (terminal-only summary line)
      integer(kind=8) :: clock_start, clock_now, clock_rate
      logical :: saved_use_structure_dt_limits
! compute_timestep's report that it switched the structure-based
! timestep limits off; the driver applies it to star%job (2026 W3).
      logical :: disable_structure_dt_limits
      integer :: saved_atm_choice
      integer :: i
      integer :: kind_card, model_iteration
      double precision :: log_r_rsun, current_zx

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

! LPUNCH is TRUE once first model is calculated
      star%punch_pending_flag = .false.
! 2026 (phase five): controls read and setup are now star-layer
! routines operating on the star_job structure (state/star_job_lib).
      call read_controls(ierr)
      if (ierr /= 0) return
! a second opacity-table set at a different Z is needed under the
! Z ramp or Z diffusion (read by the kap table loaders from
! star_setup's kap_init on; formerly decided inside kap/setupopac)
      star%use_two_z_tables = star%ctrl%use_z_ramp .or. star%job%use_diffusion_z
      call star_setup(ierr)
      if (ierr /= 0) return
! sanity pass over the adopted controls: warn (terminal only) about
! silently-inert or self-cancelling combinations (2026 audit)
      call warn_inconsistent_controls

      do monte_carlo_run_number = star%job%mc_run_start,star%job%mc_run_end
! Per-run setup (2026): the two prologue blocks are contained
! subroutines below, mirroring begin_kind_card/end_kind_card.
      call apply_monte_carlo_parameters(monte_carlo_run_number, &
           age_scale_factor)
      call begin_calibration

!**********
!     RUN THROUGH THE KIND CARDS IN ORDER
!**********
! star%job%nk is the run list's cursor (formerly the module variable
! NK doubling as the DO index). A structure component cannot be a
! DO-variable, so a local drives the loop and star%job%nk shadows it;
! the exhaustion fix-up after the loop reproduces the historical DO
! semantics (nk = num_runs+1 when the list runs out, current card
! when the calibration verdict exits early) for the post-loop
! readers (write_run_summaries' nk-2 indexing).
      call system_clock(clock_start, clock_rate)
      runs_complete = .false.
      run_loop: do kind_card = 1, star%job%num_runs
         star%job%nk = kind_card
! Per-kind-card setup (2026, core/ phase 3): everything from the
! card's initial composition through the first timestep estimate is
! begin_kind_card below; config/read errors return through ierr.
         call begin_kind_card
         if (ierr /= 0) return

! for a given kind card, evolve NMODLS(NK) times
! if rescaling is being performed, NMODLS(NK) is the number of times
! the new model is being relaxed
       do model_iteration = 1,star%job%num_models(star%job%nk)
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
!  evolve_step; a crossing returns step_leave_run_loop, which
!  abandons the current card -- no end_kind_card, no verdict -- and
!  goes straight on to the next (rescale) card) and, once the
!  luminosity there matches too, arms a final run stopped at the
!  interpolated age (star_found_flag) -- which the verdict below
!  turns into the end of the run list.
         call end_of_card_calibration(runs_complete)
         if (runs_complete) exit run_loop

! END RUN LOOP
      end do run_loop
! EXIT RUN LOOP
      if (.not. runs_complete) star%job%nk = star%job%num_runs + 1

! End-of-run summary: why the run ended + the final model (2026,
! run-log verbosity item; wall-clock goes to the terminal only).
      call system_clock(clock_now)
      call log_run_summary(dble(clock_now-clock_start)/dble(clock_rate))

! FOR MONTE CARLO, REWIND OUTPUT FILES AND WRITE OUT SNU FLUXES AND
! MODEL PARAMETERS (legacy mode only; core/monte_carlo.f90).
      call write_run_summaries(monte_carlo_run_number, &
           convergence_iterations, log_r_rsun)
      end do

! 2026 (phase five, step B): the normal end-of-job stop became this
! clean return (ierr stays 0); the CLI wrapper simply ends.
      return

contains


! ---------------------------------------------------------------
! Arm the calibration protocols for this run, if configured: setcal
! expands the run list into solar-calibration triples (recording the
! starting X/alpha guesses and saving the LPTIME/KTTAU controls that
! chkcal's cycles restore); setscal expands it into star-calibration
! pairs. Also saves the pulse-output flag the calibration cycles
! toggle.
subroutine begin_calibration
! MHP 1/93 add option to automatically calibrate solar model.
! MHP 3/96 added counter for # of iterations per converged model
      if (star%ctrl%calibrate_solar_model) then
         call setup_solar_calibration(age_scale_factor)
         convergence_iterations = 1
         saved_use_structure_dt_limits = star%job%use_structure_dt_limits   ! save LPTIME for reuse during calibration
         saved_atm_choice  = star%job%atm_choice    ! save KTTAU for reuse during calibration
      else
         convergence_iterations = 0
      endif
! DBG 12/94 add option to automatically calculate a stellar model
! of specified Teff and L
      if (star%ctrl%calibrate_star_flag) then
         call setup_star_calibration
      endif
end subroutine begin_calibration

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
         if (star%ctrl%calibrate_solar_model) then
            if (mod(star%job%nk,solar_calib_cards_per_cycle).eq.0) then
               log_r_rsun = log_r_surface_cm(star%log_L,star%log_Teff)-star%log10_solar_radius
! MHP 06/13 Add solar Z/X to observables
               current_zx = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
               call check_solar_calibration(star%log_L,log_r_rsun,star%job%nk,current_zx)
               star%job%use_structure_dt_limits = saved_use_structure_dt_limits  ! Restore LPTIME to original value for next cycle
               star%job%atm_choice  = saved_atm_choice    ! Restore KTTAU to original value for next cycle
               if (star%solar_calibration_active) then
                  star%termination_reason = 'solar calibration converged'
                  runs_complete = .true.
                  return
               else
!c MHP 8/96 added counter for # of runs needed for calibration
                  convergence_iterations = convergence_iterations + 1
! STOP AFTER 15 ATTEMPTS AT CALIBRATION
                  if (convergence_iterations.ge.15) then
                     star%termination_reason = &
                          'solar calibration NOT converged after 15 cycles'
                     runs_complete = .true.
                     return
                  end if
! 2026 retire-legacy: the pulse-trio rewind (delete-and-redo on a
! non-ultimate model) went with the .pmod/.penv/.patm files.
               end if
            endif
         endif

! DBG 12/94 NO MORE RUNS NEEDED. HAVE CALIBRATED STELLAR MODEL
         if (star%ctrl%calibrate_star_flag .and. star%star_found_flag.and.(mod(star%job%nk,star_calib_cards_per_cycle).eq.0)) then
            star%termination_reason = 'star calibration converged'
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
         star%job%initial_envelope_x = star%job%initial_x_array(star%job%nk)
         star%job%initial_envelope_z = star%job%initial_z_array(star%job%nk)
         star%mixing_length_alpha = star%job%mixing_length_array(star%job%nk)
       star%job%change_envelope_mass_flag = star%job%has_senv0_array(star%job%nk)
       star%job%requested_envelope_mass = star%job%senv0_array(star%job%nk)
       star%reset_triangle = .false.
       star%model_diverged_flag = .false.
! MHP 10/02 ZERO OUT INITIAL ANGULAR MOMENTUM
         star%total_angular_momentum = 0.0D0
         star%total_rotational_ke = 0.0D0
! read in the initial model here
! read_starting_model also performs rescaling if requested
       call read_starting_model(star%timestep_yr, star%dt, star%chosen_dt, star%trial_sign_flag, &
            star%ikut_flag, star%istore_flag, star%model_diverged_flag, &
            star%recompute_envelope_triangle, star%job%nk, star%dlnrho_dlnp, star%dlnrho_dlnt, &
            star%total_angular_momentum, star%total_rotational_ke, &
            star%convective_velocity, star%job%mixture_weights, ierr)
       if (ierr /= 0) return

! 2026: a fresh start-model load restarts the model counter, so history
! numbering begins at set_initial_model_number (default 1) instead of
! whatever counter the model file stored (the birthline library files
! carry ~63). <= 0 keeps the stored counter (restart continuation).
! Placed after read_starting_model so the stored number still gates the
! evolved-model consistency warnings there.
       if (star%job%first_call_flag(star%job%nk) .and. &
           star%ctrl%set_initial_model_number > 0) then
          star%model_number = star%ctrl%set_initial_model_number - 1
       end if

      call check_rotation_initialised(ierr)
      if (ierr /= 0) return
!     MHP 10/24 CHECK STOP CONDITIONS AND DISABLE THEM IF THE STARTING VALUES ARE BELOW THE TARGET THRESHOLD
! (2026: one table walk in stop_conditions -- the hand-written D/X/Y
! triple that used to live here carried the disarm-the-wrong-stop bug
! fixed in phase 1.)
         call init_stop_conditions(star%job%nk)
! Opt-in diagnostic (2026): the former LNUTAB per-zone neutrino
! table, off since 2004, is now the compute_neutrino_fluxes control
! (core/neutrino_flux_table.f90); it describes the starting model
! of each kind card.
      if (star%ctrl%compute_neutrino_fluxes) then
         call neutrino_flux_table
      endif
! save mass in solar units
         star%pulsation_mass_msun=star%star_mass
! MHP 08/02 STORE STARTING CZ PROPERTIES
         star%jcz = star%envelope_cz_bottom_index
         star%convective_turnover_timescale = 0.0D0
! write out headers of the appropriate output files
      call output_run_header(star%star_mass)
! DBG PULSE OUT 7/92
! initialize variables for calculating when to dump pulse output

       if (star%ctrl%helium_flash_active) then
! timestep cutting requires a model stored in logical unit ILAST
! or it will crash - so copy initial model to unit ILAST
          if (star%punch_pending_flag) then
             call write_mod_model(last_model_unit)
          endif
       endif

! locate the hydrogen-burning shell and the boundaries of the central
! and surface convection zones (if applicable).
         call locate_shell_boundaries(star%xa,star%luminosity_lsun,star%convective_flag,star%nz, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%h_shell_zone_begin, &
              star%h_shell_end_index,star%h_shell_midpoint_zone,star%has_h_shell)
! determine timestep for model
       disable_structure_dt_limits = .false.
       call compute_timestep(star%dt,star%chosen_dt,star%nz,star%logRho,star%luminosity_lsun, &
            star%m,star%dm,star%logT,star%xa,star%core_cz_top_index, &
            star%h_shell_midpoint_zone,star%luminosity_breakdown,star%dage,star%timestep_yr,star%job%nk, &
            star%logP,star%logR,star%omega,star%max_domega_frac,star%h_shell_zone_begin, &
            star%log_Teff, disable_structure_dt_limits)
       if (disable_structure_dt_limits) star%job%use_structure_dt_limits = .false.

       star%dt_saved = star%dt
! zero out entropy terms.
         do i = 1,star%nz
            star%temperature_entropy_term(i) = 0.0D0
            star%pressure_entropy_term(i) = 0.0D0
            star%luminosity_entropy_term(i) = 0.0D0
            star%radius_entropy_term(i) = 0.0D0
         end do

! zero out light element burning rates in the surface CZ.
         if (star%job%use_extended_composition) then
            star%log_rate_li6_prev = 0.0D0
            star%log_rate_li7_prev = 0.0D0
            star%log_rate_be9_prev = 0.0D0
         endif
end subroutine begin_kind_card

! ---------------------------------------------------------------
! Finish one kind card: store the last model to the .store stream
! when configured (putstore) and clear the punch flag.
subroutine end_kind_card
      use run_log_lib, only: log_final_model_line
! 2026 retire-legacy: the end-of-card putstore call is deleted with
! the .store file (punch_pending_flag keeps its helium-flash reload
! role above).
! 2026 log redesign: the last converged model of a card always gets
! a progress line, whatever the terminal_interval phase.
      call log_final_model_line()
end subroutine end_kind_card

end subroutine run_yrec

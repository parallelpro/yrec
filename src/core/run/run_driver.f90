!----------------------------------------------------------------------
! run_driver
!----------------------------------------------------------------------
! The run itself: run_yrec reads the controls, sets up the tables, and
! goes through the run cards (inside the Monte-Carlo loop, and driving
! the calibration loops), evolving each card's model with evolve_step.
module run_driver
      use neutrino_flux_output, only: write_neutrino_flux_table
      use step_driver, only: evolve_step
      use run_setup, only: prepare_run
      use shell_boundaries, only: locate_shell_boundaries
      use starting_model, only: read_starting_model
      use calibration, only: check_solar_calibration, setup_solar_calibration, setup_star_calibration
      implicit none
      private
      public :: run_yrec

contains

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! run_yrec (formerly the body of program main)
!----------------------------------------------------------------------
! The whole run: read the controls (read_controls), set up the
! constants and tables (prepare_run), then, for each Monte-Carlo run,
! go through the kind cards in order, evolving each card's model with
! evolve_step and applying the calibration verdict after each card.
! program main is a thin command-line wrapper around this routine, and
! yrec_capi is its C entry point (the embeddable engine; see
! docs/ROADMAP.md, phase five).
!
! The blanket save below keeps the locals (the rescaling and
! calibration bookkeeping) static, as they were in the main program
! this body came from. Repeated calls are legal because
! yrec_run_prologue (core/run/yrec_reset.f90) first restores the star
! state to its process-start values.
!
! ierr: 0 on a normal end of job; positive on a configuration, read or
! numerical error (the command-line wrapper exits 1); negative on
! numerics_termination (the historical "solution diverged" stop, which
! exits 0).
!
subroutine run_yrec(ierr)
      use timestep_control, only: compute_timestep
      use yrec_version, only: set_version_strings

      use net_lib
      use burn_step, only: configure_network
      use star_info_lib, only: star, i_h1, i_metals
      use yrec_output, only: output_run_header
      use io_units
      use model_files, only: write_mod_model
      use physical_constants
      use yrec_reset, only: yrec_run_prologue
      use monte_carlo, only: apply_monte_carlo_parameters, &
           write_run_summaries
      use stop_conditions, only: step_kind_card_done, &
           step_leave_run_loop, init_stop_conditions, &
           check_rotation_initialised
      use run_log, only: log_run_summary
      use controls_check, only: warn_inconsistent_controls
      use convection_lib, only: set_convection_settings
      use rotation_configuration, only: set_rotation_settings
      use rotation_settings_setup, only: rotation_settings_from_star
      use controls_reader, only: read_controls
      use microphysics_settings, only: convection_settings_from_star
      implicit none
      integer :: step_status

! --- locals ---
! Cards per calibration cycle (see the protocol comment in the run
! loop).
      integer, parameter :: solar_calib_cards_per_cycle = 3
      integer, parameter :: star_calib_cards_per_cycle = 2
      logical :: runs_complete
      integer :: monte_carlo_run_number
      double precision :: age_scale_factor
      integer :: convergence_iterations
! Wall clock for the end-of-run summary line (terminal only).
      integer(kind=8) :: clock_start, clock_now, clock_rate
      logical :: saved_use_structure_dt_limits
! Set by compute_timestep when it switches the structure-based timestep
! limits off; the driver then applies the switch to star%job.
      logical :: disable_structure_dt_limits
      integer :: saved_atm_choice
      integer :: i
      integer :: kind_card, model_iteration
      double precision :: log_r_rsun, current_zx

      integer, intent(out) :: ierr
! Keep the locals static (see the header).
      save
      ierr = 0

! Restore the process-start state, so that a repeated call behaves like
! a fresh run (core/run/yrec_reset.f90).
      call yrec_run_prologue

      call set_version_strings()

! punch_pending_flag (legacy LPUNCH) becomes true once a model has been
! computed.
      star%punch_pending_flag = .false.
! Read the controls into star%ctrl and star%job.
      call read_controls(ierr)
      if (ierr /= 0) return
! A second opacity-table set at a different Z is needed under the Z
! ramp or with Z diffusion; the kap table loaders read this flag from
! kap_init (called in prepare_run) on.
      star%use_two_z_tables = star%ctrl%use_z_ramp .or. star%job%use_diffusion_z
      call prepare_run(ierr)
      if (ierr /= 0) return
! Warn (terminal only) about control combinations that are silently
! inert or cancel each other.
      call warn_inconsistent_controls

      do monte_carlo_run_number = star%job%mc_run_start,star%job%mc_run_end
! Per-run setup: this run's Monte-Carlo parameters, then the network
! settings (which take the Monte-Carlo rate factors), then the
! calibration protocols.
      call apply_monte_carlo_parameters(monte_carlo_run_number, &
           age_scale_factor)
      call configure_network()
      call begin_calibration

! Go through the kind cards in order. star%job%nk is the run list's
! cursor (legacy NK). A structure component cannot be a DO variable, so
! the local kind_card drives the loop and star%job%nk follows it. After
! the loop, nk is num_runs+1 when the list ran out, or the current card
! when the calibration verdict ended the run early, as the legacy DO
! loop left it; write_run_summaries indexes from it (nk-2).
      call system_clock(clock_start, clock_rate)
      runs_complete = .false.
      run_loop: do kind_card = 1, star%job%num_runs
         star%job%nk = kind_card
! Set up the card, from its initial composition to the first timestep
! estimate (begin_kind_card below); configuration and read errors
! return through ierr.
         call begin_kind_card
         if (ierr /= 0) return

! Evolve the card's model num_models(nk) times (legacy NMODLS). On a
! rescaling card, num_models is the number of times the rescaled model
! is relaxed.
       do model_iteration = 1,star%job%num_models(star%job%nk)
! One model step per iteration (core/step/step_driver.f90; its header
! gives the step_status contract).
       call evolve_step(model_iteration, step_status, ierr)
       if (ierr /= 0) return
       if (step_status == step_kind_card_done) exit
       if (step_status == step_leave_run_loop) cycle run_loop
       end do

! Finish the card (end_kind_card below).
         call end_kind_card

! --- End-of-card calibration verdict ------------------------------
! The two calibration modes drive the run list with implicit card
! protocols, stated here once.
!
!  Solar calibration (calibrate_solar_model; setup_solar_calibration and
!  check_solar_calibration): the kind cards run in triples of
!  solar_calib_cards_per_cycle = 3. Card 3k+1 rescales the seed model
!  to the current (X, Z, alpha) guess, card 3k+2 evolves it to 1e8 yr
!  (settling), and card 3k+3 evolves it to the target solar age. After
!  each completed triple (mod(nk,3) == 0), check_solar_calibration
!  tests log L, log R (and optionally log Z/X) against the Sun and
!  either declares convergence or writes the Newton-corrected
!  (X, Z, alpha) into the next triple's cards. setup_solar_calibration
!  expands the run list to 16 such triples (48 cards); the verdict
!  below gives up after 14 unconverged triples.
!
!  Star calibration (calibrate_star_flag; setup_star_calibration and
!  check_star_calibration): the cards run in pairs of
!  star_calib_cards_per_cycle = 2; odd cards rescale and even cards
!  evolve. check_star_calibration watches every model of an even card
!  for the crossing of the target radius. That check lives in
!  check_stop_conditions (called by evolve_step); a crossing returns
!  step_leave_run_loop, which abandons the current card (no
!  end_kind_card, no verdict) and goes straight on to the next
!  (rescale) card. Once the luminosity at the crossing matches too, it
!  arms a final run stopped at the interpolated age (star_found_flag),
!  which the verdict below turns into the end of the run list.
         call end_of_card_calibration(runs_complete)
         if (runs_complete) exit run_loop

      end do run_loop
! The list ran out: leave nk past the last card, as the legacy DO loop
! did.
      if (.not. runs_complete) star%job%nk = star%job%num_runs + 1

! End-of-run summary: why the run ended and the final model; the wall
! clock goes to the terminal only.
      call system_clock(clock_now)
      call log_run_summary(dble(clock_now-clock_start)/dble(clock_rate))

! Monte-Carlo and solar-calibration end-of-run bookkeeping: rewind the
! output units and write the Monte-Carlo summary
! (core/run/monte_carlo.f90).
      call write_run_summaries(monte_carlo_run_number, &
           convergence_iterations, log_r_rsun)
      end do

! Normal end of job: return with ierr = 0.
      return

contains


! ---------------------------------------------------------------
! Arm the calibration protocols for this run, if configured.
! setup_solar_calibration expands the run list into solar-calibration
! triples; the starting use_structure_dt_limits and atm_choice (legacy
! LPTIME, KTTAU), which the run may change, are saved so that each
! cycle starts from them again. setup_star_calibration expands the run
! list into star-calibration pairs.
subroutine begin_calibration
! Automatic solar calibration (MHP 1/93). convergence_iterations counts
! the calibration cycles: it starts at 1, and is 0 without calibration.
      if (star%ctrl%calibrate_solar_model) then
         call setup_solar_calibration(age_scale_factor)
         convergence_iterations = 1
         saved_use_structure_dt_limits = star%job%use_structure_dt_limits   ! restored after each cycle
         saved_atm_choice  = star%job%atm_choice    ! restored after each cycle
      else
         convergence_iterations = 0
      endif
! Automatic calibration (DBG 12/94) of a star to a target L and a
! target Teff or radius.
      if (star%ctrl%calibrate_star_flag) then
         call setup_star_calibration
      endif
end subroutine begin_calibration

! ---------------------------------------------------------------
! End-of-card calibration verdict (see the protocol comment at the call
! site). Sets runs_complete when the run list is done: the solar
! calibration has converged or given up, or the star calibration's
! final run has finished on an even card. After each solar triple it
! restores use_structure_dt_limits and atm_choice for the next cycle
! and, on a miss, advances the cycle counter.
subroutine end_of_card_calibration(runs_complete)
      logical, intent(out) :: runs_complete

      runs_complete = .false.
! Solar calibration: after each completed triple, compare the model
! with the Sun. log_r_rsun = log(R/Rsun), with R from
! L = 4 pi R^2 sigma Teff^4.
         if (star%ctrl%calibrate_solar_model) then
            if (mod(star%job%nk,solar_calib_cards_per_cycle).eq.0) then
               log_r_rsun = log10_radius_from_l_teff(star%log_L,star%log_Teff, &
                    star%log10_solar_luminosity)-star%log10_solar_radius
! The surface Z/X, for the optional Z/X calibration (MHP 06/13).
               current_zx = star%xa(i_metals,star%nz)/star%xa(i_h1,star%nz)
               call check_solar_calibration(star%log_L,log_r_rsun,star%job%nk,current_zx)
               star%job%use_structure_dt_limits = saved_use_structure_dt_limits  ! restore for the next cycle
               star%job%atm_choice  = saved_atm_choice    ! restore for the next cycle
               if (star%solar_calibration_active) then
                  star%termination_reason = 'solar calibration converged'
                  runs_complete = .true.
                  return
               else
! Not converged: count the cycle (MHP 8/96).
                  convergence_iterations = convergence_iterations + 1
! Give up when the counter reaches 15. It starts at 1, so this is
! after 14 unconverged triples.
                  if (convergence_iterations.ge.15) then
                     star%termination_reason = &
                          'solar calibration NOT converged after 14 cycles'
                     runs_complete = .true.
                     return
                  end if
               end if
            endif
         endif

! Star calibration: the final run to the interpolated age, armed by
! check_star_calibration, has finished on an even card, so no more runs
! are needed.
         if (star%ctrl%calibrate_star_flag .and. star%star_found_flag.and.(mod(star%job%nk,star_calib_cards_per_cycle).eq.0)) then
            star%termination_reason = 'star calibration converged'
            runs_complete = .true.
            return
         end if
end subroutine end_of_card_calibration


! ---------------------------------------------------------------
! Start one kind card: set the card's composition, mixing length and
! envelope-mass request; read the starting model or reuse the one in
! memory (read_starting_model, with any rescaling); check the rotation
! setup; disarm the stops the starting model has already passed; write
! the optional neutrino-flux table and the output headers; save the
! model for the helium flash's timestep cuts; locate the burning shell
! and the convection-zone edges; estimate the first timestep; and zero
! the entropy terms. Sets ierr on a configuration or model-read error.
subroutine begin_kind_card
         star%job%initial_envelope_x = star%job%initial_x_array(star%job%nk)
         star%job%initial_envelope_z = star%job%initial_z_array(star%job%nk)
         star%mixing_length_alpha = star%job%mixing_length_array(star%job%nk)
! Pass the card's mixing length to the convection domain's settings,
         call set_convection_settings(convection_settings_from_star())
! and refresh the rotation domain's settings.
         call set_rotation_settings(rotation_settings_from_star())
       star%job%change_envelope_mass_flag = star%job%has_senv0_array(star%job%nk)
       star%job%requested_envelope_mass = star%job%senv0_array(star%job%nk)
       star%reset_triangle = .false.
       star%model_diverged_flag = .false.
! Zero the total angular momentum and rotational kinetic energy before
! the model read (MHP 10/02).
         star%total_angular_momentum = 0.0D0
         star%total_rotational_ke = 0.0D0
! Read the starting model, or reuse the one in memory, and rescale it
! if the card asks for that.
       call read_starting_model(star%timestep_yr, star%dt, star%chosen_dt, star%trial_sign_flag, &
            star%ikut_flag, star%istore_flag, star%model_diverged_flag, &
            star%recompute_envelope_triangle, star%job%nk, star%dlnrho_dlnp, star%dlnrho_dlnt, &
            star%total_angular_momentum, star%total_rotational_ke, &
            star%convective_velocity, star%job%mixture_weights, ierr)
       if (ierr /= 0) return

! A fresh starting-model load restarts the model counter, so the
! history numbering begins at set_initial_model_number (default 1)
! rather than at the counter stored in the model file (the birthline
! library files carry about 63). A value <= 0 keeps the stored counter,
! for continuing a run. This comes after read_starting_model so that
! the stored number still gates the evolved-model consistency warnings
! there.
       if (star%job%first_call_flag(star%job%nk) .and. &
           star%ctrl%set_initial_model_number > 0) then
          star%model_number = star%ctrl%set_initial_model_number - 1
       end if

      call check_rotation_initialised(ierr)
      if (ierr /= 0) return
! Disarm the stop conditions whose thresholds the starting model has
! already passed (e.g. a central-abundance stop already below its
! target) (MHP 10/24).
         call init_stop_conditions(star%job%nk)
! Optional diagnostic: the per-shell neutrino-flux table of each card's
! starting model (compute_neutrino_fluxes, legacy LNUTAB;
! core/diagnostics/neutrino_flux_output.f90).
      if (star%ctrl%compute_neutrino_fluxes) then
         call write_neutrino_flux_table
      endif
! The stellar mass in solar units, for the pulsation output.
         star%pulsation_mass_msun=star%star_mass
! Record the starting base of the surface convection zone and zero the
! convective turnover timescale (MHP 08/02).
         star%jcz = star%envelope_cz_bottom_index
         star%convective_turnover_timescale = 0.0D0
! Write the headers of the output files.
      call output_run_header(star%star_mass)

       if (star%ctrl%helium_flash_active) then
! Through the helium flash a timestep cut reloads the model stored on
! last_model_unit, and fails if none is there; so once a model has been
! computed, store the card's starting model on that unit.
          if (star%punch_pending_flag) then
             call write_mod_model(last_model_unit)
          endif
       endif

! Locate the hydrogen-burning shell and the boundaries of the central
! and surface convection zones (if present).
         call locate_shell_boundaries(star%xa,star%luminosity_lsun,star%convective_flag,star%nz, &
              star%core_cz_top_index,star%envelope_cz_bottom_index,star%h_shell_zone_begin, &
              star%h_shell_end_index,star%h_shell_midpoint_zone,star%has_h_shell)
! Estimate the first timestep of the card.
       disable_structure_dt_limits = .false.
       call compute_timestep(star%dt,star%chosen_dt,star%nz,star%logRho,star%luminosity_lsun, &
            star%m,star%logT,star%xa,star%core_cz_top_index, &
            star%h_shell_midpoint_zone,star%luminosity_breakdown,star%dage,star%timestep_yr,star%job%nk, &
            star%logP,star%logR,star%omega,star%max_domega_frac,star%h_shell_zone_begin, &
            star%log_Teff, disable_structure_dt_limits)
       if (disable_structure_dt_limits) star%job%use_structure_dt_limits = .false.

       star%dt_saved = star%dt
! Zero the entropy terms.
         do i = 1,star%nz
            star%temperature_entropy_term(i) = 0.0D0
            star%pressure_entropy_term(i) = 0.0D0
            star%luminosity_entropy_term(i) = 0.0D0
            star%radius_entropy_term(i) = 0.0D0
         end do

end subroutine begin_kind_card

! ---------------------------------------------------------------
! Finish one kind card: write the progress line of the card's last
! converged model, whatever the terminal_interval phase.
subroutine end_kind_card
      use run_log, only: log_final_model_line
      call log_final_model_line()
end subroutine end_kind_card

end subroutine run_yrec

end module run_driver

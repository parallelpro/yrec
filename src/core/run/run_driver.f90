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
! this body came from. The run's star, s, is one of them, and every
! routine the run calls takes it as an argument. Repeated calls are legal
! because yrec_run_prologue (core/run/yrec_reset.f90) first resets s to
! its process-start state.
!
! ierr: 0 on a normal end of job; positive on a configuration, read or
! numerical error (the command-line wrapper exits 1); negative on
! envelope_termination (the historical "solution diverged" stop, which
! exits 0).
!
!   star reads: age_gyr, card_from_starting_model, card_initial_x,
!     card_initial_z, card_log_fit_point_mass_fraction,
!     card_mixing_length_alpha, card_moves_fit_point, card_num_models,
!     control%calibrate_solar_model, control%calibrate_star_flag,
!     control%set_initial_model_number, control%use_z_ramp,
!     convective_flag, first_mc_run, last_mc_run,
!     log10_solar_luminosity, log10_solar_radius, log_L_lsun, log_Teff,
!     logP, logR, logRho, logT, luminosity_breakdown_lsun,
!     luminosity_lsun, m, num_cards, nz, omega,
!     physics%compute_neutrino_fluxes, physics%diffuse_helium_active,
!     physics%helium_flash_active, physics%use_diffusion_z,
!     solar_calibration_converged, star_found_flag, star_mass, xa
!   star writes: atm_choice_in_use, convective_turnover_timescale,
!     core_cz_top_index, current_card, current_card_envelope_x,
!     current_card_envelope_z, dln_luminosity_dt, dlog_pressure_dt,
!     dlog_radius_dt, dlog_temperature_dt, dt, dt_saved,
!     envelope_cz_bottom_index, envelope_cz_bottom_index_start,
!     fit_point_move_requested, h_shell_bot_index, h_shell_mid_index,
!     h_shell_top_index, has_h_shell, max_domega_frac,
!     mixing_length_alpha, model_diverged_flag, model_number,
!     pulsation_mass_msun, punch_pending_flag,
!     recompute_envelope_triangle,
!     requested_log_fit_point_mass_fraction, reset_triangle,
!     settling_active, stored_vertex_index, structure_dt_limits_active,
!     termination_reason, timestep_yr, total_angular_momentum,
!     total_rotational_ke, triangle_orientation, use_two_z_tables
subroutine run_yrec(ierr)
      use timestep_control, only: compute_timestep
      use yrec_version, only: set_version_strings

      use settings_install, only: install_monte_carlo_settings, install_card_settings
      use star_info_lib, only: star_info, i_h1, i_metals
      use yrec_output, only: output_state, output_open, output_close, &
           output_run_header, output_final_model_line
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
      use controls_reader, only: read_controls
      use run_files, only: open_run_files
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
! limits off; the driver then applies the switch to
! star%structure_dt_limits_active.
      logical :: disable_structure_dt_limits
      integer :: saved_atm_choice
      integer :: i
      integer :: kind_card, model_iteration
      double precision :: log_r_rsun, current_zx

      integer, intent(out) :: ierr
! The run's star (static with the other locals, allocated on the first
! call); yrec_run_prologue resets it to its process-start state.
      type(star_info), allocatable :: s
! The job's output state (the history file, the profile columns and
! counter, the progress lines), set up by output_open after the read
! and passed to every step.
      type(output_state) :: output
! Keep the locals static (see the header).
      save
      ierr = 0

! Restore the process-start state, so that a repeated call behaves like
! a fresh run (core/run/yrec_reset.f90).
      call yrec_run_prologue(s)

      call set_version_strings()

! punch_pending_flag (legacy LPUNCH) becomes true once a model has been
! computed.
      s%punch_pending_flag = .false.
! Read the controls into star%control and star%physics, and arm the run's working copies of
! the controls it changes (the run plan, the latches).
      call read_controls(s, ierr)
      if (ierr /= 0) return
! the files the deck names: the starting and final models, the run log
! (which then receives the read's messages), the Monte Carlo and
! isochrone files
      call open_run_files(ierr)
      if (ierr /= 0) return
      call output_open(s, output, ierr)
      if (ierr /= 0) return
! A second opacity-table set at a different Z is needed under the Z
! ramp or with Z diffusion; the kap table loaders read this flag from
! kap_init (called in prepare_run) on.
      s%use_two_z_tables = s%control%use_z_ramp .or. s%physics%use_diffusion_z
      call prepare_run(s, ierr)
      if (ierr /= 0) return
! Warn (on the terminal and in the run log) about control combinations
! that are silently inert or cancel each other.
      call warn_inconsistent_controls(s)

      do monte_carlo_run_number = s%first_mc_run,s%last_mc_run
! Per-run setup: this run's Monte-Carlo parameters, then the network
! settings (which take the Monte-Carlo rate factors), then the
! calibration protocols.
      call apply_monte_carlo_parameters(s, monte_carlo_run_number, &
           age_scale_factor)
      call install_monte_carlo_settings(s, ierr)
      if (ierr /= 0) return
      call begin_calibration

! Go through the kind cards in order. star%current_card is the run list's
! cursor (legacy NK). A structure component cannot be a DO variable, so
! the local kind_card drives the loop and star%current_card follows it. After
! the loop, current_card is num_cards+1 when the list ran out, or the current card
! when the calibration verdict ended the run early, as the legacy DO
! loop left it; write_run_summaries indexes from it (current_card-2).
      call system_clock(clock_start, clock_rate)
      runs_complete = .false.
      run_loop: do kind_card = 1, s%num_cards
         s%current_card = kind_card
! Set up the card, from its initial composition to the first timestep
! estimate (begin_kind_card below); configuration and read errors
! return through ierr.
         call begin_kind_card
         if (ierr /= 0) return

! Evolve the card's model card_num_models(current_card) times (legacy
! NMODLS). On a rescaling card, it is the number of times the rescaled model
! is relaxed.
       do model_iteration = 1,s%card_num_models(s%current_card)
! One model step per iteration (core/step/step_driver.f90; its header
! gives the step_status contract).
       call evolve_step(s, output, model_iteration, step_status, ierr)
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
!  each completed triple (mod(current_card,3) == 0), check_solar_calibration
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
! The list ran out: leave current_card past the last card, as the legacy DO loop
! did.
      if (.not. runs_complete) s%current_card = s%num_cards + 1

! End-of-run summary: why the run ended and the final model; the wall
! clock goes to the terminal only.
      call system_clock(clock_now)
      call log_run_summary(s, dble(clock_now-clock_start)/dble(clock_rate))

! Monte-Carlo and solar-calibration end-of-run bookkeeping: rewind the
! output units and write the Monte-Carlo summary
! (core/run/monte_carlo.f90).
      call write_run_summaries(s, monte_carlo_run_number, &
           convergence_iterations, log_r_rsun)
      end do

! Normal end of job: the history file closed, and ierr = 0.
      call output_close(output)
      return

contains


! ---------------------------------------------------------------
! Arm the calibration protocols for this run, if configured.
! setup_solar_calibration expands the run list into solar-calibration
! triples; the starting structure_dt_limits_active and atm_choice_in_use
! (legacy LPTIME, KTTAU), which the run may change, are saved so that each
! cycle starts from them again. setup_star_calibration expands the run
! list into star-calibration pairs.
subroutine begin_calibration
! Automatic solar calibration (MHP 1/93). convergence_iterations counts
! the calibration cycles: it starts at 1, and is 0 without calibration.
      if (s%control%calibrate_solar_model) then
         call setup_solar_calibration(s, age_scale_factor)
         convergence_iterations = 1
         saved_use_structure_dt_limits = s%structure_dt_limits_active   ! restored after each cycle
         saved_atm_choice  = s%atm_choice_in_use    ! restored after each cycle
      else
         convergence_iterations = 0
      endif
! Automatic calibration (DBG 12/94) of a star to a target L and a
! target Teff or radius.
      if (s%control%calibrate_star_flag) then
         call setup_star_calibration(s)
      endif
end subroutine begin_calibration

! ---------------------------------------------------------------
! End-of-card calibration verdict (see the protocol comment at the call
! site). Sets runs_complete when the run list is done: the solar
! calibration has converged or given up, or the star calibration's
! final run has finished on an even card. After each solar triple it
! restores structure_dt_limits_active and atm_choice_in_use for the next cycle
! and, on a miss, advances the cycle counter.
subroutine end_of_card_calibration(runs_complete)
      logical, intent(out) :: runs_complete

      runs_complete = .false.
! Solar calibration: after each completed triple, compare the model
! with the Sun. log_r_rsun = log(R/Rsun), with R from
! L = 4 pi R^2 sigma Teff^4.
         if (s%control%calibrate_solar_model) then
            if (mod(s%current_card,solar_calib_cards_per_cycle).eq.0) then
               log_r_rsun = log10_radius_from_l_teff(s%log_L_lsun,s%log_Teff, &
                    s%log10_solar_luminosity)-s%log10_solar_radius
! The surface Z/X, for the optional Z/X calibration (MHP 06/13).
               current_zx = s%xa(i_metals,s%nz)/s%xa(i_h1,s%nz)
               call check_solar_calibration(s, s%log_L_lsun,log_r_rsun,s%current_card,current_zx)
               s%structure_dt_limits_active = saved_use_structure_dt_limits  ! restore for the next cycle
               s%atm_choice_in_use  = saved_atm_choice    ! restore for the next cycle
               if (s%solar_calibration_converged) then
                  s%termination_reason = 'solar calibration converged'
                  runs_complete = .true.
                  return
               else
! Not converged: count the cycle (MHP 8/96).
                  convergence_iterations = convergence_iterations + 1
! Give up when the counter reaches 15. It starts at 1, so this is
! after 14 unconverged triples.
                  if (convergence_iterations.ge.15) then
                     s%termination_reason = &
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
         if (s%control%calibrate_star_flag .and. s%star_found_flag.and.(mod(s%current_card,star_calib_cards_per_cycle).eq.0)) then
            s%termination_reason = 'star calibration converged'
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
! read_starting_model's and compute_timestep's outputs this caller
! does not use
      double precision :: chosen_dt, convective_velocity
      integer :: ikut_flag
         s%current_card_envelope_x = s%card_initial_x(s%current_card)
         s%current_card_envelope_z = s%card_initial_z(s%current_card)
         s%mixing_length_alpha = s%card_mixing_length_alpha(s%current_card)
! The card's settings records: convection's takes the card's mixing
! length; rotation's and element diffusion's are refreshed.
         call install_card_settings(s)
       s%fit_point_move_requested = s%card_moves_fit_point(s%current_card)
       s%requested_log_fit_point_mass_fraction = s%card_log_fit_point_mass_fraction(s%current_card)
       s%reset_triangle = .false.
       s%model_diverged_flag = .false.
! Arm helium settling for the card (the settling routines switch it off
! when the centre runs out of hydrogen).
       s%settling_active = s%physics%diffuse_helium_active
! Zero the total angular momentum and rotational kinetic energy before
! the model read (MHP 10/02).
         s%total_angular_momentum = 0.0D0
         s%total_rotational_ke = 0.0D0
! Read the starting model, or reuse the one in memory, and rescale it
! if the card asks for that.
       call read_starting_model(s, s%timestep_yr, s%dt, chosen_dt, s%triangle_orientation, &
            ikut_flag, s%stored_vertex_index, s%model_diverged_flag, &
            s%recompute_envelope_triangle, s%current_card, &
            s%total_angular_momentum, s%total_rotational_ke, &
            convective_velocity, ierr)
       if (ierr /= 0) return

! A fresh starting-model load restarts the model counter, so the
! history numbering begins at set_initial_model_number (default 1)
! rather than at the counter stored in the model file (the birthline
! library files carry about 63). A value <= 0 keeps the stored counter,
! for continuing a run. This comes after read_starting_model so that
! the stored number still gates the evolved-model consistency warnings
! there.
       if (s%card_from_starting_model(s%current_card) .and. &
           s%control%set_initial_model_number > 0) then
          s%model_number = s%control%set_initial_model_number - 1
       end if

      call check_rotation_initialised(s, ierr)
      if (ierr /= 0) return
! Disarm the stop conditions whose thresholds the starting model has
! already passed (e.g. a central-abundance stop already below its
! target) (MHP 10/24).
         call init_stop_conditions(s, s%current_card)
! Optional diagnostic: the per-shell neutrino-flux table of each card's
! starting model (compute_neutrino_fluxes, legacy LNUTAB;
! core/diagnostics/neutrino_flux_output.f90).
      if (s%physics%compute_neutrino_fluxes) then
         call write_neutrino_flux_table(s)
      endif
! The stellar mass in solar units, for the pulsation output.
         s%pulsation_mass_msun=s%star_mass
! Record the starting base of the surface convection zone and zero the
! convective turnover timescale (MHP 08/02).
         s%envelope_cz_bottom_index_start = s%envelope_cz_bottom_index
         s%convective_turnover_timescale = 0.0D0
! Write the headers of the output files.
      call output_run_header(s, s%star_mass)

       if (s%physics%helium_flash_active) then
! Through the helium flash a timestep cut reloads the model stored on
! last_model_unit, and fails if none is there; so once a model has been
! computed, store the card's starting model on that unit.
          if (s%punch_pending_flag) then
             call write_mod_model(s, last_model_unit, replace=.true.)
          endif
       endif

! Locate the hydrogen-burning shell and the boundaries of the central
! and surface convection zones (if present).
         call locate_shell_boundaries(s, s%xa,s%luminosity_lsun,s%convective_flag,s%nz, &
              s%core_cz_top_index,s%envelope_cz_bottom_index,s%h_shell_bot_index, &
              s%h_shell_top_index,s%h_shell_mid_index,s%has_h_shell)
! Estimate the first timestep of the card.
       disable_structure_dt_limits = .false.
       call compute_timestep(s, s%dt,chosen_dt,s%nz,s%logRho,s%luminosity_lsun, &
            s%m,s%logT,s%xa,s%core_cz_top_index, &
            s%h_shell_mid_index,s%luminosity_breakdown_lsun,s%age_gyr,s%timestep_yr,s%current_card, &
            s%logP,s%logR,s%omega,s%max_domega_frac,s%h_shell_bot_index, &
            s%log_Teff, disable_structure_dt_limits)
       if (disable_structure_dt_limits) s%structure_dt_limits_active = .false.

       s%dt_saved = s%dt
! Zero the entropy terms.
         do i = 1,s%nz
            s%dlog_temperature_dt(i) = 0.0D0
            s%dlog_pressure_dt(i) = 0.0D0
            s%dln_luminosity_dt(i) = 0.0D0
            s%dlog_radius_dt(i) = 0.0D0
         end do

end subroutine begin_kind_card

! ---------------------------------------------------------------
! Finish one kind card: write the progress line of the card's last
! converged model, whatever the terminal_interval phase.
subroutine end_kind_card
      call output_final_model_line(s, output)
end subroutine end_kind_card

end subroutine run_yrec

end module run_driver

!----------------------------------------------------------------------
! test_run_contracts
!----------------------------------------------------------------------
! ok/FAIL checks of run/'s contracts, each confirmed to fail under a
! mutation (test_run_contracts.mutations.toml):
!   the star calibration (calibration.f90) on a star set up here, with
!   models placed on either side of the target radius R*:
!     - a crossing before the last card sets up the next pair of cards
!       (a new X, or the final run to the age at R*);
!     - a crossing on the last card (card 50 = max_runs) sets up
!       nothing past it and ends with a verdict: not converged, or R*
!       and L* reached with no card left for the final run;
!     - a pass that starts above R* and contracts through it crosses
!       between its own models, not at its second model (M328), and
!       the next pass starts from its own first model;
!     - end_star_calibration_card: the final run's end is the
!       converged verdict, and a last card that ends without finding
!       the star is the not-converged one;
!   the solar calibration's check (check_solar_calibration) on a star
!   set up here:
!     - its Newton step solves the 2x2 system of the fixed partial
!       derivatives, and the next triple's cards take the new X and
!       alpha; within the tolerances it is converged and writes no card;
!       with the Z/X calibration, Z is rescaled by target/model;
!   the central-abundance stops (stop_conditions.f90):
!     - a stop the starting model has passed is disarmed
!       (init_stop_conditions), a later crossing ends the card, and
!       nothing stops a card whose stops are off;
!   a Monte-Carlo run's start (monte_carlo.f90's begin_monte_carlo_run)
!   on a star set up here:
!     - without the solar calibration, cards 2 and 3 keep the deck's
!       stop ages; with it, they get the calibration's 1e8 and 4.57e9 yr;
!     - the run's solar luminosity is the deck's times its sample's
!       luminosity factor;
!   the sample file's read (setup_monte_carlo_runs), on a scratch file:
!     - a file that ends before the last run asked for is an error
!       (ierr 1), and a complete one reads every run's samples;
!   the C entry (yrec_capi's yrec_run) and the re-entrancy prologue
!   (yrec_reset), called here as a C or Fortran host would:
!     - NULL pointers are read as empty paths (the default decks;
!       M331);
!     - a unit the host opened stays open across calls, and the run's
!       own units (the run log, the final model) are closed after one,
!       and so are the files it opened with newunit, their unit
!       variables back at -1 (M312);
!     - a namelist path of exactly 256 characters (the override slots)
!       runs, and a longer one is refused with status 1: cut to 256, it
!       would name the deck, which would then run.
!
! Usage: test_run_contracts DECK, with DECK an absolute path to a short
! deck (test_run_contracts.py passes the configuration matrix's
! cm_base.inlist), run from a directory with an empty output/
! and YREC_INPUT and YREC_START set. The 256-character paths are DECK
! with "./" components after its leading "/". Prints one line per check
! and exits with status 1 when one fails.
program test_run_contracts
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_loc, c_null_char, c_null_ptr
      use star_info_lib, only: star_info
      use control_slots, only: max_runs, i_target_x, i_target_z
      use composition_layout, only: i_h1, i_he4
      use physical_constants, only: c4pi, csig
      use calibration, only: setup_star_calibration, check_star_calibration, &
           end_star_calibration_card, check_solar_calibration, dlogl_dx_fixed, &
           dlogr_dx_fixed, dlogl_dalpha_fixed, dlogr_dalpha_fixed
      use stop_conditions, only: check_stop_conditions, init_stop_conditions, &
           step_continue, step_kind_card_done, step_leave_run_loop
      use monte_carlo, only: begin_monte_carlo_run, setup_monte_carlo_runs
      use yrec_capi, only: yrec_run
      use run_messages, only: run_log_unit
      use run_files, only: last_model_unit, dynamics_unit, starting_model_unit, &
           mc_summary_unit, mc_model_unit, isochrone_unit
      implicit none
      integer :: num_failed = 0
      character(len=1024) :: deck

      if (command_argument_count() /= 1) then
         write(*,'(a)') 'usage: test_run_contracts DECK (an absolute path)'
         stop 2
      end if
      call get_command_argument(1, deck)
      if (deck(1:1) /= '/' .or. len_trim(deck) > 200) then
         write(*,'(a)') 'test_run_contracts: DECK must be an absolute path of at most 200 characters'
         stop 2
      end if

      call star_calibration_checks()
      call star_calibration_pass_checks()
      call solar_calibration_checks()
      call abundance_stop_checks()
      call monte_carlo_checks()
      call c_entry_checks()

      if (num_failed > 0) then
         write(*,'(a,i0,a)') 'test_run_contracts: ', num_failed, ' checks FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_run_contracts: done'

contains

! ---------------------------------------------------------------
! The star calibration: a target of 1 Lsun at 1 Rsun, with a
! tolerance of 1e-3 Lsun, and cards 1 and 2 laid out over the 50 cards
! by setup_star_calibration.
subroutine star_calibration_checks()
      type(star_info), allocatable :: s
      logical :: runs_complete
      character(len=96) :: not_converged
      double precision :: x_before

      write(not_converged,'(a,i0,a)') 'star calibration NOT converged in ', max_runs, ' cards'

! A pass before the last card, with L off target, sets up the next
! pair's X by the secant step.
      call new_calibration_star(s)
      call cross_target_radius(s, 0.1d0, 48)
      call run_contract_check_that('star calibration: a crossing before the last card sets the next pair''s X', &
           s%card_initial_x(49) /= 0.70d0 .and. s%card_initial_x(50) == s%card_initial_x(49) &
           .and. s%card_rescale_targets(i_target_x, 49) == s%card_initial_x(49))
      call run_contract_check_that('star calibration: ... and leaves the end reason alone', &
           s%termination_reason == 'model budget exhausted')

! The same with L on target arms the final run to the age at R*.
      call new_calibration_star(s)
      call cross_target_radius(s, 0.0d0, 48)
      call run_contract_check_that('star calibration: L on target before the last card arms the final run', &
           s%star_found_flag .and. s%card_stops_active(49) .and. s%card_stops_active(50) &
           .and. abs(s%card_max_age_yr(50)/1.5d9 - 1.0d0) < 1.0d-9 &
           .and. s%termination_reason == 'model budget exhausted')

! On the last card nothing follows: the verdicts.
      call new_calibration_star(s)
      x_before = s%card_initial_x(max_runs)
      call cross_target_radius(s, 0.1d0, max_runs)
      call run_contract_check_that('star calibration: a crossing on the last card with L off target: not converged', &
           s%termination_reason == not_converged .and. s%just_passed_target_radius_flag &
           .and. .not. s%star_found_flag .and. s%card_initial_x(max_runs) == x_before)
      call new_calibration_star(s)
      call cross_target_radius(s, 0.0d0, max_runs)
      call run_contract_check_that('star calibration: a crossing on the last card with L on target: no final run', &
           s%termination_reason == &
           'star calibration: R* and L* reached on the last card, no card left for the final run' &
           .and. .not. s%star_found_flag)

! The end-of-card verdict.
      call new_calibration_star(s)
      s%star_found_flag = .true.
      s%current_card = max_runs
      call end_star_calibration_card(s, runs_complete)
      call run_contract_check_that('star calibration: the final run''s end is the converged verdict', &
           runs_complete .and. s%termination_reason == 'star calibration converged')
      s%current_card = max_runs - 1
      s%termination_reason = 'model budget exhausted'
      call end_star_calibration_card(s, runs_complete)
      call run_contract_check_that('star calibration: no verdict after the final run''s rescale card', &
           .not. runs_complete .and. s%termination_reason == 'model budget exhausted')
      s%star_found_flag = .false.
      s%current_card = max_runs - 2
      call end_star_calibration_card(s, runs_complete)
      call run_contract_check_that('star calibration: no verdict on an even card before the last', &
           .not. runs_complete .and. s%termination_reason == 'model budget exhausted')
      s%current_card = max_runs
      call end_star_calibration_card(s, runs_complete)
      call run_contract_check_that('star calibration: the last card ended without the star: not converged', &
           .not. runs_complete .and. s%termination_reason == not_converged)
end subroutine star_calibration_checks

! ---------------------------------------------------------------
! The solar calibration's check after a triple (check_solar_calibration)
! on a star whose first triple (cards 1 to 3) ran with X = 0.71,
! Z = 0.018 and alpha = 1.8 and ended at log L = 0.01, log R = -0.004:
! the Newton step solves its 2x2 system and goes into cards 4 to 6; a
! model within the tolerances is converged and writes no card; with the
! Z/X calibration, Z is rescaled by the ratio of the target to the
! model's Z/X, and Z/X must match as well.
subroutine solar_calibration_checks()
      type(star_info), allocatable :: s
      double precision, parameter :: log_l = 0.01d0, log_r = -0.004d0
      double precision :: dx, dalpha

      call new_solar_calibration_star(s)
      s%log_L_lsun = log_l
      call check_solar_calibration(s, log_r, 0.0250d0)
      dx = s%card_rescale_targets(i_target_x, 4) - 0.71d0
      dalpha = s%card_mixing_length_alpha(4) - 1.8d0
      call run_contract_check_that('solar calibration: the Newton step solves its 2x2 system for dX and dalpha', &
           abs(dlogl_dx_fixed*dx + dlogl_dalpha_fixed*dalpha + log_l) < 1.0d-15 .and. &
           abs(dlogr_dx_fixed*dx + dlogr_dalpha_fixed*dalpha + log_r) < 1.0d-15 .and. &
           .not. s%solar_calibration_converged .and. dx /= 0.0d0)
      call run_contract_check_that('solar calibration: the next triple''s three cards take the new X and alpha', &
           all(s%card_initial_x(4:6) == s%card_rescale_targets(i_target_x, 4)) .and. &
           all(s%card_mixing_length_alpha(4:6) == s%card_mixing_length_alpha(4)) .and. &
           s%card_mixing_length_alpha(7) == 0.0d0 .and. s%card_initial_z(4) == 0.0d0)

      call new_solar_calibration_star(s)
      s%log_L_lsun = 0.5d-5
      call check_solar_calibration(s, -0.5d-5, 0.0250d0)
      call run_contract_check_that('solar calibration: L and R within their tolerances: converged, no card written', &
           s%solar_calibration_converged .and. s%card_mixing_length_alpha(4) == 0.0d0 .and. &
           s%card_initial_x(4) == 0.0d0)

      call new_solar_calibration_star(s)
      s%control%calibrate_solar_zx = .true.
      s%control%target_solar_zx = 0.0245d0
      s%control%zx_tolerance = 1.0d-4
      s%log_L_lsun = 0.5d-5
      call check_solar_calibration(s, -0.5d-5, 0.0250d0)
      call run_contract_check_that('solar calibration: with Z/X off target, Z is rescaled by target/model', &
           .not. s%solar_calibration_converged .and. &
           abs(s%card_rescale_targets(i_target_z, 4)/(0.018d0*0.0245d0/0.0250d0) - 1.0d0) < 1.0d-15 &
           .and. all(s%card_initial_z(4:6) == s%card_rescale_targets(i_target_z, 4)))
end subroutine solar_calibration_checks

! A fresh star at the end of its first solar-calibration triple, with
! tolerances of 1e-5 on log L and log R and no Z/X calibration.
subroutine new_solar_calibration_star(s)
      type(star_info), allocatable, intent(inout) :: s
      if (allocated(s)) deallocate(s)
      allocate(s)
      s%control%luminosity_tolerance = 1.0d-5
      s%control%radius_tolerance = 1.0d-5
      s%control%calibrate_solar_zx = .false.
      s%current_card = 3
      s%card_rescale_targets(i_target_x, 1) = 0.71d0
      s%card_rescale_targets(i_target_z, 1) = 0.018d0
      s%card_mixing_length_alpha(1) = 1.8d0
end subroutine new_solar_calibration_star

! ---------------------------------------------------------------
! The central-abundance stops on a star set up here: card 1 stops at
! central X = 0.10 and Y = 0.20, and its starting model has X = 0.05
! (already past its stop) and Y = 0.60. init_stop_conditions disarms
! the X stop (its threshold negated), so the model does not end the
! card; Y falling below its stop does; nothing stops a card whose stops
! are off.
subroutine abundance_stop_checks()
      type(star_info), allocatable :: s
      integer :: status

      allocate(s)
      call no_other_stops(s)
      s%current_card = 1
      s%card_kind(1) = 1
      s%card_stops_active(1) = .true.
      s%card_central_hydrogen_stop(1) = 0.10d0
      s%card_central_helium_stop(1) = 0.20d0
      s%xa(i_h1,1) = 0.05d0
      s%xa(i_he4,1) = 0.60d0
      call init_stop_conditions(s)
      call check_stop_conditions(s, 1, status)
      call run_contract_check_that('abundance stops: a stop the starting model has passed is disarmed', &
           s%card_central_hydrogen_stop(1) == -0.10d0 .and. s%card_central_helium_stop(1) == 0.20d0 &
           .and. status == step_continue .and. s%termination_reason == 'model budget exhausted')
      s%xa(i_he4,1) = 0.15d0
      s%card_stops_active(1) = .false.
      call check_stop_conditions(s, 2, status)
      call run_contract_check_that('abundance stops: none while the card''s stops are off', &
           status == step_continue)
      s%card_stops_active(1) = .true.
      call check_stop_conditions(s, 2, status)
      call run_contract_check_that('abundance stops: central Y below its stop ends the card', &
           status == step_kind_card_done .and. s%termination_reason == 'central Y below stop value')
end subroutine abundance_stop_checks

! ---------------------------------------------------------------
! A pass of the star calibration through check_stop_conditions, as
! evolve_step calls it, on a track that starts above R* = 1 Rsun and
! contracts through it (M328): models at 1.2, 1.1 and 0.99 Rsun, aged
! 1, 2 and 3 Gyr. Only the third crosses, between the second and
! itself, at 3 - 0.01/0.11 Gyr. The next pass again starts above R*.
subroutine star_calibration_pass_checks()
      type(star_info), allocatable :: s
      integer :: status(3)
      double precision :: age_at_crossing

      call new_calibration_star(s)
      call no_other_stops(s)
      s%control%calibrate_star_flag = .true.
      s%current_card = 2
      call pass_model(s, 1, 1.2d0, 1.0d0, status(1))
      call pass_model(s, 2, 1.1d0, 2.0d0, status(2))
      call run_contract_check_that('star calibration: a pass that starts above R* does not cross at its second model', &
           all(status(1:2) == step_continue) .and. .not. s%just_passed_target_radius_flag)
      call pass_model(s, 3, 0.99d0, 3.0d0, status(3))
      age_at_crossing = 3.0d0 - 0.01d0/0.11d0
      call run_contract_check_that('star calibration: ... and crosses R* downward between its own models', &
           status(3) == step_leave_run_loop .and. s%just_passed_target_radius_flag .and. &
           abs(s%age_gyr_at_target_radius/age_at_crossing - 1.0d0) < 1.0d-12)
      s%current_card = 4
      call pass_model(s, 1, 1.3d0, 0.5d0, status(1))
      call pass_model(s, 2, 1.25d0, 0.6d0, status(2))
      call run_contract_check_that('star calibration: the next pass starts from its own first model', &
           all(status(1:2) == step_continue) .and. .not. s%just_passed_target_radius_flag)
end subroutine star_calibration_pass_checks

! One model of a pass: R/Rsun = r_rsun at log(L/Lsun) = 0.1 and the age
! age_gyr, as the model_iteration-th model of the current card.
subroutine pass_model(s, model_iteration, r_rsun, age_gyr, step_status)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: model_iteration
      double precision, intent(in) :: r_rsun, age_gyr
      integer, intent(out) :: step_status
      s%log_L_lsun = 0.1d0
      s%log_Teff = log_teff_at(s, 0.1d0, r_rsun)
      s%age_gyr = age_gyr
      call check_stop_conditions(s, model_iteration, step_status)
end subroutine pass_model

! The other stops of check_stop_conditions off, as a deck without them
! leaves them: no end age or abundance stop (the cards' stops are off),
! and the structure limits at their off values.
subroutine no_other_stops(s)
      type(star_info), intent(inout) :: s
      s%card_stops_active = .false.
      s%control%log_L_upper_limit = 1.0d99
      s%control%Teff_upper_limit = 1.0d99
      s%control%log_g_upper_limit = 1.0d99
      s%control%nu_max_upper_limit = 1.0d99
      s%control%log_L_lower_limit = -1.0d99
      s%control%Teff_lower_limit = -1.0d99
      s%control%log_g_lower_limit = -1.0d99
      s%control%nu_max_lower_limit = -1.0d99
end subroutine no_other_stops

! A fresh star with the target and cards 1 and 2 (X = 0.70), laid out by
! setup_star_calibration; the previous pass's X (card 45's rescale
! target) is 0.71 with log L(R*) = 0.2, for the secant step.
subroutine new_calibration_star(s)
      type(star_info), allocatable, intent(inout) :: s
      if (allocated(s)) deallocate(s)
      allocate(s)
      s%solar_luminosity_cgs = 3.828d33
      s%solar_radius_cgs = 6.957d10
      s%control%specify_teff_flag = .false.
      s%calibration_target_radius_rsun = 1.0d0
      s%control%target_luminosity_lsun = 1.0d0
      s%control%target_star_luminosity_tolerance = 1.0d-3
      s%card_initial_x(1) = 0.70d0
      s%card_rescale_targets(i_target_x, 1) = 0.70d0
      s%card_num_models(1) = 2
      s%card_num_models(2) = 100
      s%card_kind(1) = 3
      s%card_kind(2) = 1
      call setup_star_calibration(s)
      s%card_rescale_targets(i_target_x, max_runs - 5) = 0.71d0
      s%log_l_lsun_at_target_radius_prev_run = 0.2d0
end subroutine new_calibration_star

! Two models on run_index's card, at 0.99 and 1.01 Rsun, aged 1 and 2
! Gyr, with log(L/Lsun) = log_l: the second crosses R* = 1 Rsun (at
! 1.5 Gyr).
subroutine cross_target_radius(s, log_l, run_index)
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: log_l
      integer, intent(in) :: run_index
      s%current_card = run_index
      s%log_L_lsun = log_l
      s%log_Teff = log_teff_at(s, log_l, 0.99d0)
      s%age_gyr = 1.0d0
      call check_star_calibration(s)
      s%log_Teff = log_teff_at(s, log_l, 1.01d0)
      s%age_gyr = 2.0d0
      call check_star_calibration(s)
end subroutine cross_target_radius

! log Teff of a model of log(L/Lsun) = log_l and R/Rsun = r_rsun, from
! L = 4 pi R^2 sigma Teff^4.
double precision function log_teff_at(s, log_l, r_rsun)
      type(star_info), intent(in) :: s
      double precision, intent(in) :: log_l, r_rsun
      log_teff_at = 0.25d0*log10(10.0d0**log_l*s%solar_luminosity_cgs/ &
           (c4pi*csig*(r_rsun*s%solar_radius_cgs)**2))
end function log_teff_at

! ---------------------------------------------------------------
! The start of a Monte-Carlo run (begin_monte_carlo_run) on a star whose
! deck stops cards 2 and 3 at 5e6 and 3e9 yr, with run 1's samples all
! 1: without the solar calibration the cards keep the deck's ages; with
! it, the run sets the calibration's (setup_solar_calibration then lays
! out its cycles).
subroutine monte_carlo_checks()
      type(star_info), allocatable :: s
      double precision :: age_scale_factor
      integer :: ierr
      call new_monte_carlo_star(s)
      call begin_monte_carlo_run(s, 1, age_scale_factor)
      call run_contract_check_that('Monte Carlo: without the calibration, the deck''s stop ages', &
           s%card_max_age_yr(2) == 5.0d6 .and. s%card_max_age_yr(3) == 3.0d9)
      call new_monte_carlo_star(s)
      s%control%calibrate_solar_model = .true.
      call begin_monte_carlo_run(s, 1, age_scale_factor)
      call run_contract_check_that('Monte Carlo: with the calibration, 1e8 and 4.57e9 yr', &
           s%card_max_age_yr(2) == 1.0d8 .and. s%card_max_age_yr(3) == 4.57d9)
! A sample's luminosity factor of 1.01 on the deck's 3.828e33 erg/s.
      call new_monte_carlo_star(s)
      s%mc_samples%luminosity_factor(1) = 1.01d0
      call begin_monte_carlo_run(s, 1, age_scale_factor)
      call run_contract_check_that('Monte Carlo: L_sun is the deck''s times the luminosity factor', &
           s%solar_luminosity_cgs == 3.828d33*1.01d0 .and. &
           s%log10_solar_luminosity == log10(3.828d33*1.01d0))

! The sample file, with runs 1 to 2 asked for: one record, then two (run
! 2's luminosity factor 1.01).
      call new_monte_carlo_star(s)
      s%physics%mc_run_end = 2
      call open_sample_file(s, 1)
      call setup_monte_carlo_runs(s, ierr)
      call close_sample_file()
      call run_contract_check_that('Monte Carlo: a sample file that ends before the last run: ierr 1', &
           ierr == 1)
      call new_monte_carlo_star(s)
      s%physics%mc_run_end = 2
      call open_sample_file(s, 2)
      call setup_monte_carlo_runs(s, ierr)
      call close_sample_file()
      call run_contract_check_that('Monte Carlo: a complete sample file reads every run', &
           ierr == 0 .and. s%mc_samples%luminosity_factor(2) == 1.01d0)
end subroutine monte_carlo_checks

! A scratch sample file on run_files' dynamics_unit, rewound, with
! num_records records in the shipped mc_input.dat's layout (run k's
! luminosity factor 1 + (k-1)/100), and its name for the star's
! messages.
subroutine open_sample_file(s, num_records)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: num_records
      integer :: k
      s%control%dynamics_file = 'test_run_contracts samples (scratch)'
      open(newunit=dynamics_unit, status='scratch', form='formatted')
      do k = 1, num_records
         write(dynamics_unit,'(a,i2,1p7e10.3/1pe9.3)') 'run  ', k, 1.0d0, 1.0d0, 1.0d0, &
              1.0d0, 2.29d-2, 1.0d0, 1.0d0 + 0.01d0*(k - 1), 1.0d0
      end do
      rewind(dynamics_unit)
end subroutine open_sample_file

! Close it, and set the unit back to -1, as run_files' close_run_files
! does: the run's own close would otherwise find a stale number there.
subroutine close_sample_file()
      close(dynamics_unit)
      dynamics_unit = -1
end subroutine close_sample_file

! A fresh star for a Monte-Carlo run: three cards stopping at 1e6, 5e6
! and 3e9 yr, the deck's solar luminosity, and run 1's samples all 1.
subroutine new_monte_carlo_star(s)
      type(star_info), allocatable, intent(inout) :: s
      if (allocated(s)) deallocate(s)
      allocate(s)
      s%physics%monte_carlo_active = .true.
      s%physics%solar_luminosity_cgs = 3.828d33
      s%control%num_runs = 3
      s%control%target_end_age(1:3) = [1.0d6, 5.0d6, 3.0d9]
      s%mc_samples%s11_factor(1) = 1.0d0
      s%mc_samples%s33_factor(1) = 1.0d0
      s%mc_samples%s34_factor(1) = 1.0d0
      s%mc_samples%s17_factor(1) = 1.0d0
      s%mc_samples%settling_z_factor(1) = 1.0d0
      s%mc_samples%luminosity_factor(1) = 1.0d0
      s%mc_samples%age_factor(1) = 1.0d0
end subroutine new_monte_carlo_star

! ---------------------------------------------------------------
! The C entry and the prologue, with the host's own unit 42 open
! throughout: two calls on a missing deck (the second call's prologue is
! a re-entry), then the deck at exactly 256 characters, then 257.
subroutine c_entry_checks()
      integer, parameter :: host_unit = 42, slot_length = 256
      integer(c_int) :: status, status_null
      character(len=:), allocatable :: path_256, path_long
      character(kind=c_char), allocatable, target :: missing(:), deck_256(:), deck_long(:), empty(:)
      logical :: host_open, log_open, model_open

      open(unit=host_unit, status='scratch')
      missing = c_path('no_such_deck.inlist')

      status = yrec_run(c_loc(missing), c_loc(missing))
      inquire(unit=host_unit, opened=host_open)
      call run_contract_check_that('yrec_run: a missing deck returns a positive status', status > 0)
      call run_contract_check_that('yrec_run: the host''s unit stays open after a call', host_open)

      status = yrec_run(c_loc(missing), c_loc(missing))
      inquire(unit=host_unit, opened=host_open)
      call run_contract_check_that('yrec_run_prologue: the host''s unit stays open at a re-entry', host_open)

! NULL for "the default decks" (M331): read as empty paths, so the call
! returns what a call with empty strings returns.
      empty = c_path('')
      status = yrec_run(c_loc(empty), c_loc(empty))
      status_null = yrec_run(c_null_ptr, c_null_ptr)
      call run_contract_check_that('yrec_run: NULL paths are read as empty ones (the default decks)', &
           status_null == status)

      path_256 = padded_path(trim(deck), slot_length)
      path_long = path_256 // 'x'
      deck_256 = c_path(path_256)
      deck_long = c_path(path_long)
      status = yrec_run(c_loc(deck_256), c_loc(deck_256))
      inquire(unit=host_unit, opened=host_open)
      inquire(unit=run_log_unit, opened=log_open)
      inquire(unit=last_model_unit, opened=model_open)
      call run_contract_check_that('yrec_run: a path of exactly 256 characters runs (status 0)', status == 0)
      call run_contract_check_that('yrec_run: the run''s log and final-model units are closed after it', &
           .not. log_open .and. .not. model_open)
      call run_contract_check_that('yrec_run: ... and the host''s unit is still open', host_open)
! The files the run opened with newunit (M312): closed, their unit
! variables back at -1.
      call run_contract_check_that('yrec_run: the run''s newunit files are closed after it too', &
           all([starting_model_unit, dynamics_unit, mc_summary_unit, mc_model_unit, &
           isochrone_unit] == -1))

      status = yrec_run(c_loc(deck_long), c_loc(deck_256))
      call run_contract_check_that('yrec_run: a control path over 256 characters is refused (status 1)', &
           status == 1)
      status = yrec_run(c_loc(deck_256), c_loc(deck_long))
      call run_contract_check_that('yrec_run: a physics path over 256 characters is refused (status 1)', &
           status == 1)
      close(host_unit)
end subroutine c_entry_checks

! The path as a C string: its characters and a NUL.
function c_path(path) result(cstr)
      character(len=*), intent(in) :: path
      character(kind=c_char), allocatable :: cstr(:)
      integer :: i
      allocate(cstr(len(path)+1))
      do i = 1, len(path)
         cstr(i) = path(i:i)
      end do
      cstr(len(path)+1) = c_null_char
end function c_path

! The absolute path padded to length characters with "./" components
! after its leading "/" (and one extra "/" when the padding is odd): the
! same file.
function padded_path(path, length) result(padded)
      character(len=*), intent(in) :: path
      integer, intent(in) :: length
      character(len=:), allocatable :: padded
      integer :: pad
      pad = length - len(path)
      padded = '/' // repeat('./', pad/2) // repeat('/', mod(pad, 2)) // path(2:)
end function padded_path

! ---------------------------------------------------------------
! One check: its label and ok, or FAIL.
subroutine run_contract_check_that(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      character(len=96) :: padded
      padded = label
      if (condition) then
         write(*,'(a,a)') padded, '  ok'
      else
         write(*,'(a,a)') padded, '  FAIL'
         num_failed = num_failed + 1
      end if
end subroutine run_contract_check_that

end program test_run_contracts

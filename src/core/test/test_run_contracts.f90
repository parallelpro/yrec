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
!     - end_star_calibration_card: the final run's end is the
!       converged verdict, and a last card that ends without finding
!       the star is the not-converged one;
!   the C entry (yrec_capi's yrec_run) and the re-entrancy prologue
!   (yrec_reset), called here as a C or Fortran host would:
!     - a unit the host opened stays open across calls, and the run's
!       own units (the run log, the final model) are closed after one;
!     - a namelist path of exactly 256 characters (the override slots)
!       runs, and a longer one is refused with status 1: cut to 256, it
!       would name the deck, which would then run.
!
! Usage: test_run_contracts DECK, with DECK an absolute path to a short
! deck (test_run_contracts.py passes the configuration matrix's
! cm_env_overshoot.inlist), run from a directory with an empty output/
! and YREC_INPUT and YREC_START set. The 256-character paths are DECK
! with "./" components after its leading "/". Prints one line per check
! and exits with status 1 when one fails.
program test_run_contracts
      use, intrinsic :: iso_c_binding, only: c_char, c_int, c_loc, c_null_char
      use star_info_lib, only: star_info
      use control_slots, only: max_runs, i_target_x
      use physical_constants, only: c4pi, csig
      use calibration, only: setup_star_calibration, check_star_calibration, &
           end_star_calibration_card
      use yrec_capi, only: yrec_run
      use run_messages, only: run_log_unit
      use run_files, only: last_model_unit
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
! The C entry and the prologue, with the host's own unit 42 open
! throughout: two calls on a missing deck (the second call's prologue is
! a re-entry), then the deck at exactly 256 characters, then 257.
subroutine c_entry_checks()
      integer, parameter :: host_unit = 42, slot_length = 256
      integer(c_int) :: status
      character(len=:), allocatable :: path_256, path_long
      character(kind=c_char), allocatable, target :: missing(:), deck_256(:), deck_long(:)
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

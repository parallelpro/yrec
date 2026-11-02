!----------------------------------------------------------------------
! test_controls
!----------------------------------------------------------------------
! Reads one deck through read_controls, exactly as a run does, and
! prints every resulting control: each element of star%control and
! star%physics, of the buffer-only members and of the run copies
! (dump_controls, generated: controls_dump.inc), and the unit numbers.
! Doubles
! print at full precision, so two dumps agree exactly when the controls
! are bit-identical.
!   test_controls                  the pristine defaults (no deck read)
!   test_controls deck.inlist      a single-file inlist
!   test_controls deck.nml1 deck.nml2   a legacy CONTROL/PHYSICS pair
! With YREC_TEST_READ_TWICE set to a second deck (and, for a legacy pair,
! YREC_TEST_READ_TWICE_PHYSICS to its .nml2), the deck given is read, the
! run copies scribbled on as a run would, and then the second deck read
! into the same star: the dump is the second read's, which must equal
! that deck read alone.
! The output files a run opens go under $YREC_OUTPUT (test_controls.py
! points it at a scratch directory). test_controls.py runs this on every
! deck in examples/ and testsuite/.
program test_controls
      use star_info_lib, only: star_info, control_nml_override, physics_nml_override
      use run_messages, only: run_log_unit, terminal_unit
      use run_files, only: last_model_unit, open_run_files
      use controls_reader, only: read_controls
      implicit none
! the test's star (star_info_lib declares no instance); static, as it
! is too large for the stack
      type(star_info), save :: star
      character(len=256) :: deck1, deck2, second1, second2
      integer :: ierr, u, i1, nargs
      logical :: kept_first = .false.

      u = 6
      nargs = command_argument_count()
      if (nargs >= 1) then
         call get_command_argument(1, deck1)
         deck2 = deck1
         if (nargs >= 2) call get_command_argument(2, deck2)
         control_nml_override = deck1
         physics_nml_override = deck2
! read_controls writes its banner to standard output; the dump goes to
! unit 6 after a marker so the two can be told apart
         call read_controls(star, ierr)
! the run's files, as run_yrec opens them after the read
         if (ierr == 0) call open_run_files(star, ierr)
! YREC_TEST_READ_TWICE=<deck>: scribble on the run's working copies of
! the controls as a run would (a card's mixing length, the atmosphere in
! use, a Monte-Carlo cross-section scale), close the run's files, then
! read the second deck into the same star; the dump must be that deck's
! alone (nothing of the first read left), and the read must have armed
! the copies from their controls again
         call get_environment_variable("YREC_TEST_READ_TWICE", second1)
         if (len_trim(second1) > 0) then
            call get_environment_variable("YREC_TEST_READ_TWICE_PHYSICS", second2)
            if (len_trim(second2) == 0) second2 = second1
            control_nml_override = second1
            physics_nml_override = second2
            star%card_mixing_length_alpha(1) = 7.0d0
            star%atm_choice_in_use = 7
            star%cross_section_scale_in_use(1) = 7.0d0
! the run's files closed, as run_yrec's prologue closes them
            do i1 = 7, 99
               close(i1)
            end do
            call read_controls(star, ierr)
            if (ierr == 0) call open_run_files(star, ierr)
            kept_first = any(star%card_mixing_length_alpha /= star%control%mixing_length_array) &
                 .or. star%atm_choice_in_use /= star%physics%atm_choice .or. &
                 any(star%cross_section_scale_in_use /= star%physics%cross_section_scale)
         end if
         write(u,'(a)') "# controls"
! (in the dump, so that test_controls.py's comparison with one read sees it)
         if (kept_first) write(u,'(a)') "a run copy kept the first run's value"
         write(u,'(a,i0)') "ierr = ", ierr
      else
         write(u,'(a)') "# controls (pristine defaults)"
      end if
      call dump_controls(star, u)
      write(u,'(a," = ",i0)') 'last_model_unit', last_model_unit
      write(u,'(a," = ",i0)') 'run_log_unit', run_log_unit
      write(u,'(a," = ",i0)') 'terminal_unit', terminal_unit

contains

! every control (generated: the members the star does not hold, registry
! home buffer, from the read's buffer)
      include 'controls_dump.inc'

function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=12) :: buf
      write(buf,'(i0)') i
      s = trim(buf)
end function itoa

end program test_controls

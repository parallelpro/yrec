!----------------------------------------------------------------------
! test_controls
!----------------------------------------------------------------------
! Reads one deck through read_controls, exactly as a run does, and
! prints every resulting control: each element of star%ctrl
! (controls_dump.inc, generated from its declarations) and of the
! buffer-only members, the star fields read_controls sets directly and
! the unit numbers. Doubles
! print at full precision, so two dumps agree exactly when the controls
! are bit-identical.
!   test_controls                  the pristine defaults (no deck read)
!   test_controls deck.inlist      a single-file inlist
!   test_controls deck.nml1 deck.nml2   a legacy CONTROL/PHYSICS pair
! (with YREC_TEST_READ_TWICE=1, the deck read twice into the same star).
! The output files a run opens go under $YREC_OUTPUT (test_controls.py
! points it at a scratch directory). test_controls.py runs this on every
! deck in examples/ and testsuite/.
program test_controls
      use star_info_lib, only: star_info, control_nml_override, physics_nml_override
      use io_units, only: last_model_unit, run_log_unit, terminal_unit
      use controls_reader, only: read_controls
      use run_files, only: open_run_files
! the members the star does not hold (registry home buffer) are dumped
! from the read's buffer
      use controls_buffer
      implicit none
! the test's star (star_info_lib declares no instance); static, as it
! is too large for the stack
      type(star_info), save :: star
      character(len=256) :: deck1, deck2
      character(len=8) :: twice
      integer :: ierr, u, i1, i2, nargs

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
         if (ierr == 0) call open_run_files(ierr)
! YREC_TEST_READ_TWICE=1: scribble on the run's working copies of the
! controls as a run would (the Newton factor, the card's envelope X, a
! card's mixing length), close the run's files, then read the deck again
! into the same star; the dump must be a single read's, and the read
! must have armed the copies from their controls again
         call get_environment_variable("YREC_TEST_READ_TWICE", twice)
         if (twice == "1") then
            star%newton_correction_factor = 7.0d0
            star%current_card_envelope_x = 7.0d0
            star%card_mixing_length_alpha(1) = 7.0d0
! the run's files closed, as run_yrec's prologue closes them
            do i1 = 7, 99
               close(i1)
            end do
            call read_controls(star, ierr)
            if (ierr == 0) call open_run_files(ierr)
            if (star%newton_correction_factor /= star%ctrl%solver_underrelax_factor .or. &
                 star%current_card_envelope_x /= star%ctrl%initial_envelope_x .or. &
                 any(star%card_mixing_length_alpha /= star%ctrl%mixing_length_array)) &
                 write(u,'(a)') "a run copy kept the first run's value"
         end if
         write(u,'(a)') "# controls"
         write(u,'(a,i0)') "ierr = ", ierr
      else
         write(u,'(a)') "# controls (pristine defaults)"
      end if
      include 'controls_dump.inc'
      write(u,'(a," = ",i0)') 'star%atm_choice_initial', star%atm_choice_initial
      write(u,'(a," = ",l1)') 'star%use_ttau_relation', star%use_ttau_relation
      write(u,'(a," = ",es25.17e3)') 'star%envelope_store_log_temperature', star%envelope_store_log_temperature
      do i1 = 1, size(star%cross_section_scale)
         write(u,'(a," = ",es25.17e3)') 'star%cross_section_scale('//itoa(i1)//')', &
              star%cross_section_scale(i1)
      end do
      do i1 = 1, size(star%extra_cross_section_scale)
         write(u,'(a," = ",es25.17e3)') 'star%extra_cross_section_scale('//itoa(i1)//')', &
              star%extra_cross_section_scale(i1)
      end do
      write(u,'(a," = ",es25.17e3)') 'star%envelope_hydrogen_fraction', star%envelope_hydrogen_fraction
      write(u,'(a," = ",es25.17e3)') 'star%envelope_metal_fraction', star%envelope_metal_fraction
      write(u,'(a," = ",i0)') 'last_model_unit', last_model_unit
      write(u,'(a," = ",i0)') 'run_log_unit', run_log_unit
      write(u,'(a," = ",i0)') 'terminal_unit', terminal_unit

contains

function itoa(i) result(s)
      integer, intent(in) :: i
      character(len=:), allocatable :: s
      character(len=12) :: buf
      write(buf,'(i0)') i
      s = trim(buf)
end function itoa

end program test_controls

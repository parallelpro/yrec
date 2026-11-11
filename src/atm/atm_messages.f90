!----------------------------------------------------------------------
! atm_messages
!----------------------------------------------------------------------
! The atmosphere's messages, in one place (rotation_messages and
! eos_messages are the models). Each message is a procedure: the caller
! says what happened, and this module says it. The text and the
! destinations are the F77's and its port's, kept byte for byte (list-
! directed lines included, with their trailing blanks) until the B
! phase's message commit rewrites them.
!
! Destinations:
!   the terminal  standard output (what atm's write(*) meant); the
!                 Kurucz table's range message and the 2026 errors go to
!                 run_messages' terminal_unit instead (standard output
!                 once the run has assigned its units; test_atm points
!                 it at its log)
!   the run log   run_messages' run_log_unit
!
! The messages:
!   the surface-pressure tables (Kurucz, Castelli/Kurucz)
!     report_table_open_failure          the table cannot be opened: terminal
!     report_surface_pressure_outside    a point outside the table: terminal, run log
!     log_spline_failure                 a spline failed: run log
!   the Allard lookups
!     report_allard_teff_above           log Teff above the table: both
!     report_allard_gravity_outside      log g outside the table: both
!     report_allard_teff_below           log Teff below the table, an error: both
!     report_allard_gravity_outside_row  log g outside a stencil row's range: both
!     report_polint_failure              a polint failed: terminal
!     log_allard_lookup                  every lookup's result: run log (M235)
!   the Allard reader
!     report_allard_open_failure         the file cannot be opened: terminal
!     log_old_nextgen_file, log_allard_file_description
!                                        the file's layout: run log
!     report_allard_read_error           a record cannot be read: terminal
!     report_allard_unmatched            a record's Teff or log g is not in
!                                        the grid the first pass built: terminal
!     report_allard_invalid_file         neither layout: both
!     report_allard_reader_stopped       after each of the reader's errors: both
!   the Allard table's check
!     log_allard_bad_entry, log_allard_few_rows, log_allard_short_row
!                                        what is wrong: run log
!     report_allard_table_bad            the table is refused: both
!     log_allard_tables                  the table refused: run log
!   the star layer's switch between a table and the gray relation
!     note_choice_switch                 at the table's limit: terminal
!                                        (core/structure/surface_boundary's
!                                        note, kept as it was)
!   the 2026 checks
!     report_atm_error                   "ERROR: <text>": terminal, run log
!                                        (run_messages' report_error)
module atm_messages
      use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
      use run_messages, only: run_log_unit, terminal_unit, report_error
      use numerics_lib, only: numerics_status_text
      use atm_tables, only: allard_table
      use atm_choices, only: atm_kurucz, atm_allard
      use math_lib
      implicit none
      private
      public :: report_table_open_failure, report_surface_pressure_outside, &
           log_spline_failure
      public :: report_allard_teff_above, report_allard_gravity_outside, &
           report_allard_teff_below, report_allard_gravity_outside_row, &
           report_polint_failure, log_allard_lookup
      public :: report_allard_open_failure, log_old_nextgen_file, &
           log_allard_file_description, report_allard_read_error, report_allard_unmatched, &
           report_allard_invalid_file, report_allard_reader_stopped
      public :: log_allard_bad_entry, log_allard_few_rows, log_allard_short_row, &
           report_allard_table_bad, log_allard_tables
      public :: report_atm_error
      public :: note_choice_switch

contains

! ======================================================================
! The surface-pressure tables

subroutine report_table_open_failure(path)
      character(len=*), intent(in) :: path
      write(output_unit,*) 'atm_init: cannot open the atmosphere table ', trim(path)
end subroutine report_table_open_failure

subroutine report_surface_pressure_outside(log10_teff, log10_gravity)
      real(dp), intent(in) :: log10_teff, log10_gravity
      write(terminal_unit,911) log10_teff, log10_gravity
      write(run_log_unit,911) log10_teff, log10_gravity
  911 format(1X,'DESIRED ATMOSPHERE OUTSIDE TABLE RANGE'/ &
           ' LOG TEFF',F10.6,' LOG G',F10.6/' RUN STOPPED')
end subroutine report_surface_pressure_outside

subroutine log_spline_failure(status)
      integer, intent(in) :: status
      write(run_log_unit,*) 'splint: '//numerics_status_text(status)
end subroutine log_spline_failure

! ======================================================================
! The Allard lookups. The range messages: a blank line, a title and the
! limits with the value between them, list-directed, to the run log and
! the terminal.

subroutine report_allard_teff_above(log_teff_min, log_teff, log_teff_max)
      real(dp), intent(in) :: log_teff_min, log_teff, log_teff_max
      call write_range(run_log_unit, 'ALSURFP: TEFFL greater than TeffLmax: ', &
           '    TEFFLmin,TEFFL,TEFFLmax; ', log_teff_min, log_teff, log_teff_max)
      call write_range(output_unit, 'ALSURFP: TEFFL greater than TeffLmax: ', &
           '    TEFFLmin,TEFFL,TEFFLmax; ', log_teff_min, log_teff, log_teff_max)
end subroutine report_allard_teff_above

subroutine report_allard_gravity_outside(log_g_min, log_g, log_g_max)
      real(dp), intent(in) :: log_g_min, log_g, log_g_max
      call write_range(run_log_unit, 'ALSURFP: GL out of max range: ', &
           '    GLXmin,GL,GLXmax; ', log_g_min, log_g, log_g_max)
      call write_range(output_unit, 'ALSURFP: GL out of max range: ', &
           '    GLXmin,GL,GLXmax; ', log_g_min, log_g, log_g_max)
end subroutine report_allard_gravity_outside

subroutine report_allard_teff_below(log_teff_min, log_teff, log_teff_max)
      real(dp), intent(in) :: log_teff_min, log_teff, log_teff_max
      call write_range(run_log_unit, 'ALSURFP: TEFFL less than TEFFLmin: ', &
           '    TEFFLmin,TEFFL,TEFFLmax; ', log_teff_min, log_teff, log_teff_max)
      call write_range(output_unit, 'ALSURFP: TEFFL less than TEFFLmin: ', &
           '    TEFFLmin,TEFFL,TEFFLmax; ', log_teff_min, log_teff, log_teff_max)
      write(output_unit,*)
      write(output_unit,*) '******** ALSURFP: Program Terminated ********'
      write(output_unit,*)
      write(run_log_unit,*)
      write(run_log_unit,*) '******** ALSURF: Program Terminated ********'
      write(run_log_unit,*)
end subroutine report_allard_teff_below

subroutine report_allard_gravity_outside_row(log_g_min, log_g, log_g_max)
      real(dp), intent(in) :: log_g_min, log_g, log_g_max
      call write_range(run_log_unit, 'ALSURFP: GL out of extended range: ', &
           '    GLmin,GL,GLmax; ', log_g_min, log_g, log_g_max)
      call write_range(output_unit, 'ALSURFP: GL out of extended range: ', &
           '    GLmin,GL,GLmax; ', log_g_min, log_g, log_g_max)
end subroutine report_allard_gravity_outside_row

subroutine write_range(unit, title, labels, low, value, high)
      integer, intent(in) :: unit
      character(len=*), intent(in) :: title, labels
      real(dp), intent(in) :: low, value, high
      write(unit,*)
      write(unit,*) title
      write(unit,*) labels, low, value, high
end subroutine write_range

subroutine report_polint_failure(status)
      integer, intent(in) :: status
      write(output_unit,*) 'polint: '//numerics_status_text(status)
end subroutine report_polint_failure

! the result of a lookup: log P at (log Teff, log g), or at tau = 100
! with allard_use_tau100, and log T
subroutine log_allard_lookup(use_tau100, log_teff, log_g, log10_temperature, log10_pressure)
      logical, intent(in) :: use_tau100
      real(dp), intent(in) :: log_teff, log_g, log10_temperature, log10_pressure
      if (use_tau100) then
         write(run_log_unit,73)
   73    format('ALSURFP:  PL,TL at Tau=100 INTERPOLATED', ' FROM ALLARD TABULATED VALUES:')
      else
         write(run_log_unit,74)
   74    format('ALSURFP: PL,TL at T=TEFF,GL INTERPOLATED ', 'FROM ALLARD TABULATED VALUES:')
      end if
      write(run_log_unit,75) log_teff, log_g, log10_temperature, log10_pressure
   75 format(' ',10x,'LOG(Teffl) =',f10.5,', LOG(G) =',f10.5, &
           ', LOG(T) =',f10.5,', Log(P) =',f10.5)
end subroutine log_allard_lookup

! ======================================================================
! The Allard reader

subroutine report_allard_open_failure(path)
      character(len=*), intent(in) :: path
      write(output_unit,*)
      write(output_unit,*) 'ALFILEIN: Allard File OPEN Failure'
      write(output_unit,*) '     File name: ', trim(path)
end subroutine report_allard_open_failure

subroutine log_old_nextgen_file()
      write(run_log_unit,*) 'ALFileIn: File Description: 1999 NEXTGEN', ' (Old NEXTGEN)'
end subroutine log_old_nextgen_file

! the first 47 characters of the file's first record after its count
subroutine log_allard_file_description(description)
      character(len=*), intent(in) :: description
      write(run_log_unit,'(2a)') 'ALFilein: New Allard Atm: File Description: ', description
end subroutine log_allard_file_description

! text: the pass's message
subroutine report_allard_read_error(text)
      character(len=*), intent(in) :: text
      write(output_unit,*)
      write(output_unit,*) text
      write(output_unit,*)
end subroutine report_allard_read_error

! title, then "<label> <value> <tail>"
subroutine report_allard_unmatched(title, label, value, tail)
      character(len=*), intent(in) :: title, label, tail
      real(dp), intent(in) :: value
      write(output_unit,*) title
      write(output_unit,*) label, value, tail
end subroutine report_allard_unmatched

subroutine report_allard_invalid_file(first_record)
      character(len=*), intent(in) :: first_record
      write(run_log_unit,*) '*** Invalid Allard Atmosphere input file ***'
      write(run_log_unit,*) "First record was '", first_record, "'"
      write(output_unit,*) '*** Invalid Allard Atmosphere input file ***'
      write(output_unit,*) "First record was '", first_record, "'"
end subroutine report_allard_invalid_file

subroutine report_allard_reader_stopped()
      write(output_unit,*)
      write(output_unit,*)
      write(output_unit,*) '**************** PROGRAM ALFilein TERMINATED ', &
           '@ 9999 **************'
      write(output_unit,*)
      write(output_unit,*)
      write(output_unit,*)
      write(run_log_unit,*)
      write(run_log_unit,*) '**************** PROGRAM ALFilein TERMINATED ', &
           '@ 9999 *************'
      write(run_log_unit,*)
      write(run_log_unit,*)
end subroutine report_allard_reader_stopped

! ======================================================================
! The Allard table's check

subroutine log_allard_bad_entry(teff, log_g)
      real(dp), intent(in) :: teff, log_g
      write(run_log_unit,'(2a,f5.0,f7.2)') 'ALTABINIT: Bad input Allard Table: ', &
           'TEFF, GL: ', teff, log_g
end subroutine log_allard_bad_entry

subroutine log_allard_few_rows(num_teff)
      integer, intent(in) :: num_teff
      write(run_log_unit,'(2a,i3)') 'ALTABINIT: Bad input Allard Table: ', &
           'Less than 4 rows: nTeff = ', num_teff
end subroutine log_allard_few_rows

subroutine log_allard_short_row(row, num_valid, teff, log_g_min, log_g_max)
      integer, intent(in) :: row, num_valid
      real(dp), intent(in) :: teff, log_g_min, log_g_max
      write(run_log_unit,'(2a,2i4,2x,3f7.2)') 'ALTABINIT: Bad input Allard Table: ', &
           'Row with less that 4 elements: i,#,Teff,GLMin,GLmax: ', &
           row, num_valid, teff, log_g_min, log_g_max
end subroutine log_allard_short_row

subroutine report_allard_table_bad()
      write(output_unit,*)
      write(output_unit,*) '******** ALTABINIT: Program Terminated ********'
      write(output_unit,*)
      write(run_log_unit,*)
      write(run_log_unit,*) '******** ALTABINIT: Program Terminated ********'
      write(run_log_unit,*)
end subroutine report_allard_table_bad

! the table refused, with each row's index, Teff, first and last valid
! log g column and log g range, for each of its three quantities
subroutine log_allard_tables(tbl)
      type(allard_table), intent(in) :: tbl
      integer :: i, j
      write(run_log_unit,*)
      write(run_log_unit,*) ' **************** Allard Table  ***************'
      write(run_log_unit,*)
      write(run_log_unit,*) 'nTeff= ', tbl%num_teff, ' ,nGL=', tbl%num_g
      write(run_log_unit,*) 'nFeH= ', tbl%num_feh, ' ,nAlpha=', tbl%num_alpha
      write(run_log_unit,*)
      write(run_log_unit,*) 'TEFFLmin= ', tbl%log_teff_min, ' ,TEFFLmax= ', tbl%log_teff_max
      write(run_log_unit,*) 'GLXmin=', tbl%log_g_min, ' GLXmax=', tbl%log_g_max
      write(run_log_unit,'(A10,16X,A,20(5X,F6.2,5X))') ' PL @ Teff', 'GL ->', &
           (tbl%log_g(j), j=1,tbl%num_g)
      write(run_log_unit,*)
      write(run_log_unit,*) 'PL at Teff,GL'
      write(run_log_unit,*)
      do i = 1, tbl%num_teff
         call write_row(i, tbl%log10_pressure(i,1:tbl%num_g))
      end do
      write(run_log_unit,*)
      write(run_log_unit,*) 'PL at Tau=100'
      write(run_log_unit,*)
      do i = 1, tbl%num_teff
         call write_row(i, tbl%log10_pressure_tau100(i,1:tbl%num_g))
      end do
      write(run_log_unit,*)
      write(run_log_unit,*) 'TL at Tau=100'
      write(run_log_unit,*)
      do i = 1, tbl%num_teff
         call write_row(i, tbl%log10_temperature_tau100(i,1:tbl%num_g))
      end do
      write(run_log_unit,*)

contains

subroutine write_row(i, values)
      integer, intent(in) :: i
      real(dp), intent(in) :: values(:)
      write(run_log_unit,'(I5,F6.0,2I4,2F6.2,1P20E16.8)') i, pow(10.0d0, tbl%log_teff(i)), &
           tbl%first_column(i), tbl%last_column(i), &
           tbl%row_log_g_min(i), tbl%row_log_g_max(i), values
end subroutine write_row

end subroutine log_allard_tables

! ======================================================================
! The switch between a table and the gray relation: from_choice to
! to_choice at log10_teff; allard_log_teff_last, the Allard tables'
! switch

subroutine note_choice_switch(from_choice, to_choice, log10_teff, allard_log_teff_last)
      integer, intent(in) :: from_choice, to_choice
      real(dp), intent(in) :: log10_teff, allard_log_teff_last
      character(len=*), parameter :: kurucz_to_gray = &
           "('LOG TEFF OF ',F7.3,' ABOVE 3.95 - SWITCH','TO GRAY ATMOSPHERE BOUNDARY CONDITION')"
      character(len=*), parameter :: allard_to_gray = &
           "('LOG TEFF OF ',F7.3,' ABOVE Allard Table max ',F7.3," // &
           "'  - SWITCH TO GRAY ATMOSPHERE BOUNDARY CONDITION')"
      character(len=*), parameter :: gray_to_kurucz = &
           "('LOG TEFF OF ',F7.3,' BELOW 3.95 - SWITCH',' BACK TO KURUCZ ATMOSPHERE BOUNDARY CONDITION')"
      character(len=*), parameter :: gray_to_allard = &
           "('LOG TEFF OF ',F7.3,' below Allard Table max ',F7.3," // &
           "'  - SWITCH BACK TO ALLARD ATMOSPHERE BOUNDARY CONDITION')"

      if (from_choice == atm_kurucz) then
         write(output_unit,*)
         write(output_unit,kurucz_to_gray) log10_teff
         write(output_unit,*)
      else if (from_choice == atm_allard) then
         write(output_unit,*)
         write(output_unit,allard_to_gray) log10_teff, allard_log_teff_last
         write(output_unit,*)
      else if (to_choice == atm_kurucz) then
         write(output_unit,gray_to_kurucz) log10_teff
      else if (to_choice == atm_allard) then
         write(output_unit,gray_to_allard) log10_teff, allard_log_teff_last
      end if
end subroutine note_choice_switch

! ======================================================================
! The 2026 checks' errors

subroutine report_atm_error(text)
      character(len=*), intent(in) :: text
      call report_error(text)
end subroutine report_atm_error

end module atm_messages

!----------------------------------------------------------------------
! atm_messages
!----------------------------------------------------------------------
! The atmosphere's messages, in one place (rotation_messages and
! eos_messages are the models). Each message is a procedure: the caller
! says what happened, and this module says it, in one line of sentence
! case after the entry point that met it ("atm_init: ...",
! "atm_get_surface_pt: ...") or after "note: ". Two destinations:
!   the terminal and the run log  every error atm returns ierr for, a
!                 lookup outside the Allard tables (which the caller
!                 falls back from), and the switch between a table and
!                 the gray relation;
!   the run log   what the Allard reader found: the file's layout and
!                 description and, with a table it refuses, each reason
!                 and the table.
! The terminal is standard output, whether or not the run has assigned
! io's units (as kap's and eos's messages).
!
! The messages:
!   the surface-pressure tables (Kurucz, Castelli/Kurucz)
!     report_table_open_failure          the table cannot be opened
!     report_table_read_failure          a record cannot be read
!     report_surface_pressure_outside    a point outside the table
!     report_table_spline_failure        a spline failed
!   the Allard lookups
!     report_allard_teff_above           log Teff above the tables
!     report_allard_gravity_outside      log g outside the tables
!     report_allard_teff_below           log Teff below the tables (an error)
!     report_allard_gravity_outside_row  log g outside a stencil row's range
!     report_polint_failure              a polint failed
!   the Allard reader
!     report_allard_open_failure, report_allard_read_failure,
!     report_allard_invalid_file, report_allard_unmatched,
!     report_allard_composition_missing, report_allard_without_tau100
!                                        its errors
!     log_old_nextgen_file, log_allard_file_description
!                                        the file's layout: run log
!   the Allard table's check
!     report_allard_few_rows             fewer rows than a stencil needs
!     log_allard_missing_entry, log_allard_short_row
!                                        what is wrong: run log
!     report_allard_table_bad            the table is refused
!     log_allard_tables                  the table refused: run log
!   atm_lib's checks
!     report_unknown_choice, report_table_not_read, report_not_a_table
!   the star layer's switch between a table and the gray relation
!     note_choice_switch                 at the table's limit
module atm_messages
      use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
      use run_messages, only: run_log_unit
      use numerics_lib, only: numerics_status_text
      use atm_tables, only: allard_table
      use atm_choices, only: atm_kurucz, atm_allard, atm_choice_name, &
           kurucz_gray_switch_log_teff
      use math_lib
      implicit none
      private
      public :: report_table_open_failure, report_table_read_failure, &
           report_surface_pressure_outside, report_table_spline_failure
      public :: report_allard_teff_above, report_allard_gravity_outside, &
           report_allard_teff_below, report_allard_gravity_outside_row, &
           report_polint_failure
      public :: report_allard_open_failure, report_allard_read_failure, &
           report_allard_invalid_file, report_allard_unmatched, &
           report_allard_composition_missing, report_allard_without_tau100, &
           log_old_nextgen_file, log_allard_file_description
      public :: report_allard_few_rows, log_allard_missing_entry, log_allard_short_row, &
           report_allard_table_bad, log_allard_tables
      public :: report_unknown_choice, report_table_not_read, report_not_a_table
      public :: note_choice_switch

! long enough for a line with a table path of 256 characters
      integer, parameter :: line_length = 1024

contains

! ======================================================================
! The two destinations

! a line on the terminal and in the run log
subroutine atm_report(line)
      character(len=*), intent(in) :: line
      write(output_unit,'(a)') trim(line)
      write(run_log_unit,'(a)') trim(line)
end subroutine atm_report

! a line in the run log
subroutine atm_log(line)
      character(len=*), intent(in) :: line
      write(run_log_unit,'(a)') trim(line)
end subroutine atm_log

! x with the given decimals, without blanks
function fixed(x, decimals) result(text)
      real(dp), intent(in) :: x
      integer, intent(in) :: decimals
      character(len=:), allocatable :: text
      character(len=40) :: buffer
      character(len=16) :: form
      write(form,'(a,i0,a)') '(f40.', decimals, ')'
      write(buffer,form) x
      text = trim(adjustl(buffer))
end function fixed

! ======================================================================
! The surface-pressure tables

subroutine report_table_open_failure(path)
      character(len=*), intent(in) :: path
      call atm_report('atm_init: cannot open the atmosphere table ' // trim(path))
end subroutine report_table_open_failure

! message: the read's iomsg
subroutine report_table_read_failure(path, message)
      character(len=*), intent(in) :: path, message
      call atm_report('atm_init: cannot read the atmosphere table ' // trim(path) // ': ' // &
           trim(message))
end subroutine report_table_read_failure

! atm_choice: the table's (its name in the message)
subroutine report_surface_pressure_outside(atm_choice, log10_teff, log10_gravity)
      integer, intent(in) :: atm_choice
      real(dp), intent(in) :: log10_teff, log10_gravity
      call atm_report('atm_get_surface_pt: log Teff = ' // fixed(log10_teff, 6) // &
           ', log g = ' // fixed(log10_gravity, 6) // ' is outside the ' // &
           trim(atm_choice_name(atm_choice)))
end subroutine report_surface_pressure_outside

subroutine report_table_spline_failure(atm_choice, status)
      integer, intent(in) :: atm_choice, status
      call atm_report('atm_get_surface_pt: a spline in the ' // trim(atm_choice_name(atm_choice)) // &
           ' failed: ' // numerics_status_text(status))
end subroutine report_table_spline_failure

! ======================================================================
! The Allard lookups: the value, and the range it is outside (the range
! a lookup may use: the tables' first and last rows and the rows' log g
! ranges, extended).

subroutine report_allard_teff_above(log_teff_min, log_teff, log_teff_max)
      real(dp), intent(in) :: log_teff_min, log_teff, log_teff_max
      call atm_report('atm_get_surface_pt: log Teff = ' // fixed(log_teff, 6) // &
           ' is above the Allard tables (log Teff ' // fixed(log_teff_min, 6) // ' to ' // &
           fixed(log_teff_max, 6) // ')')
end subroutine report_allard_teff_above

subroutine report_allard_gravity_outside(log_g_min, log_g, log_g_max)
      real(dp), intent(in) :: log_g_min, log_g, log_g_max
      call atm_report('atm_get_surface_pt: log g = ' // fixed(log_g, 6) // &
           ' is outside the Allard tables (log g ' // fixed(log_g_min, 6) // ' to ' // &
           fixed(log_g_max, 6) // ')')
end subroutine report_allard_gravity_outside

subroutine report_allard_teff_below(log_teff_min, log_teff, log_teff_max)
      real(dp), intent(in) :: log_teff_min, log_teff, log_teff_max
      call atm_report('atm_get_surface_pt: log Teff = ' // fixed(log_teff, 6) // &
           ' is below the Allard tables (log Teff ' // fixed(log_teff_min, 6) // ' to ' // &
           fixed(log_teff_max, 6) // ')')
end subroutine report_allard_teff_below

! log_teff_row: the log Teff of the stencil's row
subroutine report_allard_gravity_outside_row(log_teff_row, log_g_min, log_g, log_g_max)
      real(dp), intent(in) :: log_teff_row, log_g_min, log_g, log_g_max
      call atm_report('atm_get_surface_pt: log g = ' // fixed(log_g, 6) // &
           ' is outside the Allard tables'' row at log Teff ' // fixed(log_teff_row, 6) // &
           ' (log g ' // fixed(log_g_min, 6) // ' to ' // fixed(log_g_max, 6) // ')')
end subroutine report_allard_gravity_outside_row

subroutine report_polint_failure(status)
      integer, intent(in) :: status
      call atm_report('atm_get_surface_pt: the interpolation in the Allard tables failed: ' // &
           numerics_status_text(status))
end subroutine report_polint_failure

! ======================================================================
! The Allard reader

subroutine report_allard_open_failure(path)
      character(len=*), intent(in) :: path
      call atm_report('atm_init: cannot open the Allard table ' // trim(path))
end subroutine report_allard_open_failure

! where: the part of the file (its header, the first pass, the second
! pass); message: the read's iomsg, if it has one
subroutine report_allard_read_failure(path, where, message)
      character(len=*), intent(in) :: path, where
      character(len=*), intent(in), optional :: message
      character(len=line_length) :: line
      line = 'atm_init: cannot read the Allard table ' // trim(path) // ' (' // where // ')'
      if (present(message)) line = trim(line) // ': ' // trim(message)
      call atm_report(line)
end subroutine report_allard_read_failure

subroutine report_allard_invalid_file(path, first_record)
      character(len=*), intent(in) :: path, first_record
      call atm_report('atm_init: ' // trim(path) // ' is in neither Allard layout (its first record: ''' // &
           trim(first_record) // ''')')
end subroutine report_allard_invalid_file

! quantity: 'Teff' or 'log g'; a record's value the first pass did not
! collect
subroutine report_allard_unmatched(path, quantity, value)
      character(len=*), intent(in) :: path, quantity
      real(dp), intent(in) :: value
      call atm_report('atm_init: the Allard table ' // trim(path) // ': a record''s ' // quantity // &
           ' = ' // fixed(value, 2) // ' is not in the grid of the first pass')
end subroutine report_allard_unmatched

! the deck's [Fe/H] and [alpha/Fe], and the file's values (sorted)
subroutine report_allard_composition_missing(path, target_feh, target_alpha, file_feh, &
     file_alpha)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: target_feh, target_alpha, file_feh(:), file_alpha(:)
      character(len=line_length) :: line
      integer :: k
      line = 'atm_init: allard_target_feh = ' // fixed(target_feh, 2) // &
           ' and allard_target_alpha = ' // fixed(target_alpha, 2) // ': the Allard file ' // &
           trim(path) // ' has no table there; its [Fe/H]:'
      do k = 1, size(file_feh)
         line = trim(line) // ' ' // fixed(file_feh(k), 2)
      end do
      line = trim(line) // '; its [alpha/Fe]:'
      do k = 1, size(file_alpha)
         line = trim(line) // ' ' // fixed(file_alpha(k), 2)
      end do
      call atm_report(line)
end subroutine report_allard_composition_missing

subroutine report_allard_without_tau100(path)
      character(len=*), intent(in) :: path
      call atm_report('atm_init: allard_use_tau100 asks for log P and log T at tau = 100, ' // &
           'which the old NextGen 1 file ' // trim(path) // ' does not have')
end subroutine report_allard_without_tau100

subroutine log_old_nextgen_file()
      call atm_log('atm_init: the Allard table is in the old NextGen 1 layout (1999)')
end subroutine log_old_nextgen_file

! the file's description: its first record after the count
subroutine log_allard_file_description(description)
      character(len=*), intent(in) :: description
      call atm_log('atm_init: the Allard table''s description: ' // trim(description))
end subroutine log_allard_file_description

! ======================================================================
! The Allard table's check

subroutine report_allard_few_rows(num_teff, min_rows)
      integer, intent(in) :: num_teff, min_rows
      character(len=line_length) :: line
      write(line,'(a,i0,a,i0)') 'atm_init: the Allard table has ', num_teff, &
           ' log Teff rows at the deck''s composition; its lookups need ', min_rows
      call atm_report(line)
end subroutine report_allard_few_rows

! a -999 between a row's first and last valid log g column
subroutine log_allard_missing_entry(teff, log_g)
      real(dp), intent(in) :: teff, log_g
      character(len=line_length) :: line
      write(line,'(a,i0,2a)') 'atm_init: the Allard table has no log P at Teff = ', nint(teff), &
           ' K, log g = ', fixed(log_g, 2) // ', inside its row''s valid log g'
      call atm_log(line)
end subroutine log_allard_missing_entry

subroutine log_allard_short_row(teff, num_valid, min_valid)
      real(dp), intent(in) :: teff
      integer, intent(in) :: num_valid, min_valid
      character(len=line_length) :: line
      write(line,'(a,i0,a,i0,a,i0)') 'atm_init: the Allard table''s row at Teff = ', nint(teff), &
           ' K has ', num_valid, ' valid log g columns; its lookups need ', min_valid
      call atm_log(line)
end subroutine log_allard_short_row

subroutine report_allard_table_bad()
      call atm_report('atm_init: the Allard table cannot be interpolated; ' // &
           'the run log lists why, and the table')
end subroutine report_allard_table_bad

! the table refused: its extents and ranges, then for each of its three
! quantities a row per log Teff (its index, Teff, first and last valid
! log g column, and log g range, then the values at each log g)
subroutine log_allard_tables(tbl)
      type(allard_table), intent(in) :: tbl
      character(len=line_length) :: line
      integer :: i, j
      write(line,'(a,i0,a,i0)') 'atm_init: the Allard table as read: log Teff rows ', &
           tbl%num_teff, ', log g columns ', tbl%num_g
      call atm_log(line)
      call atm_log('  the lookups'' range: log Teff ' // fixed(tbl%log_teff_min, 6) // ' to ' // &
           fixed(tbl%log_teff_max, 6) // ', log g ' // fixed(tbl%log_g_min, 6) // ' to ' // &
           fixed(tbl%log_g_max, 6))
      write(line,'(a31,20(5x,f6.2,5x))') '  row  Teff valid g range  g ->', &
           (tbl%log_g(j), j=1,tbl%num_g)
      call atm_log(line)
      call atm_log('  log P at T = Teff')
      do i = 1, tbl%num_teff
         call write_row(i, tbl%log10_pressure(i,1:tbl%num_g))
      end do
      call atm_log('  log P at tau = 100')
      do i = 1, tbl%num_teff
         call write_row(i, tbl%log10_pressure_tau100(i,1:tbl%num_g))
      end do
      call atm_log('  log T at tau = 100')
      do i = 1, tbl%num_teff
         call write_row(i, tbl%log10_temperature_tau100(i,1:tbl%num_g))
      end do

contains

subroutine write_row(i, values)
      integer, intent(in) :: i
      real(dp), intent(in) :: values(:)
      write(line,'(I5,F6.0,2I4,2F6.2,1P20E16.8)') i, pow(10.0d0, tbl%log_teff(i)), &
           tbl%first_column(i), tbl%last_column(i), &
           tbl%row_log_g_min(i), tbl%row_log_g_max(i), values
      call atm_log(line)
end subroutine write_row

end subroutine log_allard_tables

! ======================================================================
! atm_lib's checks

! atm_choice is none of the choices
subroutine report_unknown_choice(atm_choice)
      integer, intent(in) :: atm_choice
      character(len=line_length) :: line
      write(line,'(a,i0,a)') 'atm_init: atm_choice = ', atm_choice, ' is not an atmosphere choice' // &
           ' (0 Eddington gray, 1 Krishna-Swamy, 2 HSRA, 3 Kurucz, 4 Allard, 5 Castelli/Kurucz)'
      call atm_report(line)
end subroutine report_unknown_choice

! a lookup in a table atm_init did not read
subroutine report_table_not_read(atm_choice)
      integer, intent(in) :: atm_choice
      character(len=line_length) :: line
      write(line,'(a,i0,3a)') 'atm_get_surface_pt: the table of atm_choice ', atm_choice, &
           ' (', trim(atm_choice_name(atm_choice)), ') was not read; atm_init reads only the table ' // &
           'of the deck''s choice'
      call atm_report(line)
end subroutine report_table_not_read

subroutine report_not_a_table(atm_choice)
      integer, intent(in) :: atm_choice
      character(len=line_length) :: line
      write(line,'(a,i0,a)') 'atm_get_surface_pt: atm_choice ', atm_choice, &
           ' is not a table (3 Kurucz, 4 Allard, 5 Castelli/Kurucz)'
      call atm_report(line)
end subroutine report_not_a_table

! ======================================================================
! The switch between a table and the gray relation: from_choice to
! to_choice at log10_teff; allard_log_teff_last, the Allard tables'
! switch (their last log Teff)

subroutine note_choice_switch(from_choice, to_choice, log10_teff, allard_log_teff_last)
      integer, intent(in) :: from_choice, to_choice
      real(dp), intent(in) :: log10_teff, allard_log_teff_last
      character(len=*), parameter :: kurucz_limit = ', the Kurucz table''s limit', &
           allard_limit = ', the Allard tables'' last log Teff'

      if (from_choice == atm_kurucz) then
         call atm_report('note: log Teff = ' // fixed(log10_teff, 3) // ' is at or above ' // &
              fixed(kurucz_gray_switch_log_teff, 3) // kurucz_limit // &
              ': switching to a gray atmosphere')
      else if (from_choice == atm_allard) then
         call atm_report('note: log Teff = ' // fixed(log10_teff, 3) // ' is at or above ' // &
              fixed(allard_log_teff_last, 3) // allard_limit // &
              ': switching to a gray atmosphere')
      else if (to_choice == atm_kurucz) then
         call atm_report('note: log Teff = ' // fixed(log10_teff, 3) // ' is below ' // &
              fixed(kurucz_gray_switch_log_teff, 3) // kurucz_limit // &
              ': switching back to the Kurucz table')
      else if (to_choice == atm_allard) then
         call atm_report('note: log Teff = ' // fixed(log10_teff, 3) // ' is below ' // &
              fixed(allard_log_teff_last, 3) // allard_limit // &
              ': switching back to the Allard tables')
      end if
end subroutine note_choice_switch

!======================================================================
! Sources: the F77's SURFP, ALSURFP, ALFILEIN, ALTABINIT and SURFBC said
! most of these (in capitals, with their routines' names, and with the
! "Program Terminated" banners of a run that went on); the texts are
! 2026's, in sentence case after the entry point that met each.

end module atm_messages

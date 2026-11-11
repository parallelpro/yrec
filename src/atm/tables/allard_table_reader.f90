!----------------------------------------------------------------------
! allard_table_reader
!----------------------------------------------------------------------
! Reads the Allard model atmospheres (NextGen, AMES-Cond, AMES-Dusty)
! into an allard_table at atm_init, and checks them:
!     read_allard_tables  the file, in either of its two layouts; then
!                         find_table_ranges, the lookups' limits, and the
!                         table's check
! The two layouts, told apart by the first record:
!   the old NextGen 1 file (1999; '3.50' in columns 4-7): one composition,
!     old_nextgen_num_g log g columns and old_nextgen_num_teff log Teff
!     rows of log P at T = Teff, in fixed formats (no tau = 100 columns);
!   the Allard layout (2007; 'ALLARD' in columns 4-9): a count of header
!     records and the file's description, the rest of the header, then a
!     record per (Teff, log g, [Fe/H], [alpha/Fe]) with log P at T = Teff
!     and log P and log T at tau = 100. A first pass collects the Teff and
!     log g of the records at the deck's [Fe/H] and [alpha/Fe]
!     (allard_target_feh, allard_target_alpha), a second fills the table.
! An entry the file does not give stays allard_invalid (-999). Errors
! return ierr = 1 after their message (atm_messages).
module allard_table_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use atm_tables, only: allard_table, allard_max_teff
      use atm_messages, only: report_allard_open_failure, log_old_nextgen_file, &
           log_allard_file_description, report_allard_read_error, report_allard_unmatched, &
           report_allard_invalid_file, report_allard_reader_stopped, log_allard_bad_entry, &
           log_allard_few_rows, log_allard_short_row, report_allard_table_bad, &
           log_allard_tables, report_atm_error
      use math_lib
      implicit none
      private
      public :: read_allard_tables

! an entry the file does not give, and the test for a valid one
      real(dp), parameter :: allard_invalid = -999.0d0, valid_above = -998.0d0
! the old NextGen 1 layout: its mark, its extents and formats (log g,
! log Teff, then log P a log Teff row at a time)
      character(len=*), parameter :: old_nextgen_mark = '3.50'
      integer, parameter :: old_nextgen_num_teff = 54, old_nextgen_num_g = 5
      character(len=*), parameter :: old_log_g_format = '(5f7.2/)'
      character(len=*), parameter :: old_log_teff_format = '(10(1p5e16.8,/),1p4e16.8,/)'
      character(len=*), parameter :: old_pressure_format = '(1p5e16.8)'
! the Allard layout: its mark, the header's first record (the count of
! header records, the description), a record's grid point, and a whole
! record
      character(len=*), parameter :: allard_mark = 'ALLARD'
      character(len=*), parameter :: header_format = '(i2,x,a)'
      character(len=*), parameter :: grid_point_format = '(f6.0,3f6.2)'
      character(len=*), parameter :: record_format = '(f6.0,3f6.2,1p4d16.8)'
! the description's length in the run log's note
      integer, parameter :: description_length = 47
! A record is at the deck's composition when its [Fe/H] is within
! feh_tolerance of allard_target_feh and its [alpha/Fe] equals
! allard_target_alpha; the second pass stores it when both are within
! grid_tolerance, which also matches a record's Teff and log g to the
! grid (the F77's three tests).
      real(dp), parameter :: feh_tolerance = 1.0d-5, grid_tolerance = 1.0d-6
! A row's lookups may use log g down to low_gravity_widths of its first
! column width below its first valid column, and one width above its
! last (the F77's; its comments said one width at both ends, M233).
      real(dp), parameter :: low_gravity_widths = 4.0d0
! the fewest rows, and valid columns in a row, a 4-point stencil needs
      integer, parameter :: min_rows = 4, min_valid_columns = 4

contains

!----------------------------------------------------------------------
! read_allard_tables
!----------------------------------------------------------------------
! Reads the file at path into tbl (target_feh and target_alpha select the
! composition of the Allard layout), then finds the lookups' limits and
! checks the table. A file that cannot be opened or read, a record off
! the grid, a composition the file does not hold, or a table without a
! 4-point stencil returns ierr = 1 after its message.
subroutine read_allard_tables(path, target_feh, target_alpha, tbl, ierr)
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: target_feh, target_alpha
      type(allard_table), intent(out) :: tbl
      integer, intent(out) :: ierr

      character(len=35) :: first_record
      integer :: allard_unit, ios

      ierr = 0
      tbl%log_teff = allard_invalid
      tbl%log_g = allard_invalid
      tbl%log10_pressure = allard_invalid
      tbl%log10_pressure_tau100 = allard_invalid
      tbl%log10_temperature_tau100 = allard_invalid

      open(newunit=allard_unit, file=path, status='OLD', iostat=ios)
      if (ios /= 0) then
         call report_allard_open_failure(path)
         call report_allard_reader_stopped()
         ierr = 1
         return
      end if
      read(allard_unit, '(a)', iostat=ios) first_record
      rewind(allard_unit)
      if (ios /= 0) then
         call report_atm_error('cannot read the Allard table ' // trim(path))
         ierr = 1
      else if (first_record(4:7) == old_nextgen_mark) then
         call read_old_nextgen(allard_unit, path, tbl, ierr)
      else if (first_record(4:9) == allard_mark) then
         call read_allard_layout(allard_unit, path, target_feh, target_alpha, tbl, ierr)
      else
         call report_allard_invalid_file(first_record)
         call report_allard_reader_stopped()
         ierr = 1
      end if
      close(allard_unit)
      if (ierr /= 0) return

      call find_table_ranges(tbl, ierr)
      if (ierr /= 0) return
      tbl%loaded = .true.
end subroutine read_allard_tables

! The old NextGen 1 layout: the log g columns, the log Teff rows, then log
! P a row at a time. (Its tables have no tau = 100 columns; a run that
! asks for them gets -999 there: O19, retired/allard_tau100_check.f90.txt.)
subroutine read_old_nextgen(allard_unit, path, tbl, ierr)
      integer, intent(in) :: allard_unit
      character(len=*), intent(in) :: path
      type(allard_table), intent(inout) :: tbl
      integer, intent(inout) :: ierr

      integer :: row, ios
      character(len=256) :: message

      tbl%is_old_nextgen = .true.
      call log_old_nextgen_file()
      tbl%num_teff = old_nextgen_num_teff
      tbl%num_g = old_nextgen_num_g
      read(allard_unit, old_log_g_format, iostat=ios, iomsg=message) tbl%log_g(1:tbl%num_g)
      if (ios == 0) read(allard_unit, old_log_teff_format, iostat=ios, iomsg=message) &
           tbl%log_teff(1:tbl%num_teff)
      if (ios == 0) read(allard_unit, old_pressure_format, iostat=ios, iomsg=message) &
           (tbl%log10_pressure(row,1:tbl%num_g), row=1,tbl%num_teff)
      if (ios /= 0) then
         call report_atm_error('cannot read the Allard table ' // trim(path) // ': ' // &
              trim(message))
         ierr = 1
      end if
end subroutine read_old_nextgen

! The Allard layout: the header, the first pass (the grid at the deck's
! composition), the second (the table).
subroutine read_allard_layout(allard_unit, path, target_feh, target_alpha, tbl, ierr)
      integer, intent(in) :: allard_unit
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: target_feh, target_alpha
      type(allard_table), intent(inout) :: tbl
      integer, intent(inout) :: ierr

      character(len=256) :: description
      integer :: num_header_records, ios
      real(dp) :: teffs(allard_max_teff)

      tbl%is_old_nextgen = .false.
      read(allard_unit, header_format, iostat=ios) num_header_records, description
      if (ios /= 0) then
         call report_atm_error('cannot read the Allard table ' // trim(path))
         ierr = 1
         return
      end if
      call log_allard_file_description(description(1:description_length))
      call skip_records(allard_unit, num_header_records - 1)

      call collect_grid(allard_unit, path, target_feh, target_alpha, tbl, teffs, ierr)
      if (ierr /= 0) return

      rewind(allard_unit)
      call skip_records(allard_unit, num_header_records)
      call fill_table(allard_unit, target_feh, target_alpha, teffs, tbl, ierr)
end subroutine read_allard_layout

! skips n records (none for n < 1)
subroutine skip_records(allard_unit, n)
      integer, intent(in) :: allard_unit, n
      integer :: k
      do k = 1, n
         read(allard_unit, *)
      end do
end subroutine skip_records

! The first pass: the distinct Teff (teffs, sorted, not log) and log g
! (tbl%log_g, sorted) of the records at the deck's composition. A
! composition the file does not hold is an error that names the file's
! [Fe/H] and [alpha/Fe] values (M221).
subroutine collect_grid(allard_unit, path, target_feh, target_alpha, tbl, teffs, ierr)
      integer, intent(in) :: allard_unit
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: target_feh, target_alpha
      type(allard_table), intent(inout) :: tbl
      real(dp), intent(out) :: teffs(:)
      integer, intent(inout) :: ierr

! the [Fe/H] and [alpha/Fe] values of the file, for that message
      integer, parameter :: max_file_compositions = 50
      real(dp) :: file_feh(max_file_compositions), file_alpha(max_file_compositions)
      integer :: num_file_feh, num_file_alpha
      real(dp) :: teff, log_g, feh, alpha
      integer :: ios

      teffs = 0.0d0
      tbl%num_teff = 0
      tbl%num_g = 0
      num_file_feh = 0
      num_file_alpha = 0
      do
         read(allard_unit, grid_point_format, iostat=ios) teff, log_g, feh, alpha
         if (ios > 0) then
            call report_allard_read_error('AtTabInit Pass 1 File Error ')
            call report_allard_reader_stopped()
            ierr = 1
            return
         end if
         if (ios < 0) exit
         call add_distinct(feh, file_feh, num_file_feh)
         call add_distinct(alpha, file_alpha, num_file_alpha)
         if (abs(feh - target_feh) > feh_tolerance .or. abs(alpha - target_alpha) > 0.0d0) cycle
         call add_distinct(teff, teffs, tbl%num_teff)
         call add_distinct(log_g, tbl%log_g, tbl%num_g)
      end do
      if (tbl%num_teff == 0) then
         call report_missing_composition()
         ierr = 1
         return
      end if
! the table holds one composition (its count, for the table's printout)
      tbl%num_feh = 1
      tbl%num_alpha = 1
      call sort_shell(tbl%num_teff, teffs)
      call sort_shell(tbl%num_g, tbl%log_g)

contains

! the deck's [Fe/H] and [alpha/Fe], and the file's
subroutine report_missing_composition()
      character(len=1024) :: text
      character(len=16) :: item
      integer :: k
      call sort_shell(num_file_feh, file_feh)
      call sort_shell(num_file_alpha, file_alpha)
      write(text,'(a,f6.2,a,f6.2,3a)') 'allard_target_feh =', target_feh, &
           ' and allard_target_alpha =', target_alpha, &
           ': the Allard file ', trim(path), ' has no table there; its [Fe/H]:'
      do k = 1, num_file_feh
         write(item,'(f6.2)') file_feh(k)
         text = trim(text) // ' ' // adjustl(item)
      end do
      text = trim(text) // '; its [alpha/Fe]:'
      do k = 1, num_file_alpha
         write(item,'(f6.2)') file_alpha(k)
         text = trim(text) // ' ' // adjustl(item)
      end do
      call report_atm_error(trim(text))
end subroutine report_missing_composition

end subroutine collect_grid

! value appended to list(1:n) unless within grid_tolerance of an entry
! already there (or the list is full)
pure subroutine add_distinct(value, list, n)
      real(dp), intent(in) :: value
      real(dp), intent(inout) :: list(:)
      integer, intent(inout) :: n
      if (any(abs(list(1:n) - value) < grid_tolerance)) return
      if (n == size(list)) return
      n = n + 1
      list(n) = value
end subroutine add_distinct

! The second pass: each record at the deck's composition into its row
! (its Teff in teffs) and column (its log g in tbl%log_g); the rows' log
! Teff. A record whose Teff or log g the first pass did not collect is an
! error.
subroutine fill_table(allard_unit, target_feh, target_alpha, teffs, tbl, ierr)
      integer, intent(in) :: allard_unit
      real(dp), intent(in) :: target_feh, target_alpha
      real(dp), intent(in) :: teffs(:)
      type(allard_table), intent(inout) :: tbl
      integer, intent(inout) :: ierr

      real(dp) :: teff, log_g, feh, alpha, log_p, log_p_tau100, log_t_tau100
      integer :: row, column, ios

      do
         read(allard_unit, record_format, iostat=ios) teff, log_g, feh, alpha, &
              log_p, log_p_tau100, log_t_tau100
         if (ios > 0) then
            call report_allard_read_error('AtTabInit File Read error')
            call report_allard_reader_stopped()
            ierr = 1
            return
         end if
         if (ios < 0) exit
         if (abs(feh - target_feh) > feh_tolerance .or. abs(alpha - target_alpha) > 0.0d0) cycle

         do row = 1, tbl%num_teff
            if (abs(teff - teffs(row)) < grid_tolerance) exit
         end do
         if (row > tbl%num_teff) then
            call report_allard_unmatched('ALFilein: Impossible failure #1', 'Teff of ', teff, &
                 ' not in Teff table.')
            call report_allard_reader_stopped()
            ierr = 1
            return
         end if
         do column = 1, tbl%num_g
            if (abs(log_g - tbl%log_g(column)) < grid_tolerance) exit
         end do
         if (column > tbl%num_g) then
            call report_allard_unmatched('ALFilein: Impossible failure #2', 'GL of ', log_g, &
                 ' not in GL table.')
            call report_allard_reader_stopped()
            ierr = 1
            return
         end if

         if (abs(feh - target_feh) < grid_tolerance .and. &
              abs(alpha - target_alpha) < grid_tolerance) then
            tbl%log_teff(row) = log10(teff)
            tbl%log_g(column) = log_g
            tbl%log10_pressure(row,column) = log_p
            tbl%log10_pressure_tau100(row,column) = log_p_tau100
            tbl%log10_temperature_tau100(row,column) = log_t_tau100
         end if
      end do
end subroutine fill_table

!----------------------------------------------------------------------
! find_table_ranges
!----------------------------------------------------------------------
! The lookups' limits, then the table's check:
! 1. log Teff: the first and last rows, and one row's width beyond each
!    (log_teff_min, log_teff_max);
! 2. each row's first and last valid log g column;
! 3. each row's log g range: low_gravity_widths of its first column
!    width below its first valid column, one width above its last; the
!    table's, the widest of them (log_g_min, log_g_max);
! 4. the check: at least four rows (tested first: the limits need two),
!    no -999 between a row's first and last valid column, and four valid
!    columns in every row. A table that fails is printed to the run log,
!    and ierr = 1.
subroutine find_table_ranges(tbl, ierr)
      type(allard_table), intent(inout) :: tbl
      integer, intent(inout) :: ierr

      integer :: row, column, first, last, n
      logical :: table_is_bad

      n = tbl%num_teff
      if (n < min_rows) then
         call log_allard_few_rows(n)
         call report_allard_table_bad()
         ierr = 1
         return
      end if

! 1. log Teff
      tbl%log_teff_first = tbl%log_teff(1)
      tbl%log_teff_last = tbl%log_teff(n)
      tbl%log_teff_min = tbl%log_teff(1) - (tbl%log_teff(2) - tbl%log_teff(1))
      tbl%log_teff_max = tbl%log_teff(n) + (tbl%log_teff(n) - tbl%log_teff(n-1))

! 2. each row's valid columns
      do row = 1, n
         do column = 1, tbl%num_g
            if (tbl%log10_pressure(row,column) > valid_above) then
               tbl%first_column(row) = column
               exit
            end if
         end do
         do column = tbl%num_g, 1, -1
            if (tbl%log10_pressure(row,column) > valid_above) then
               tbl%last_column(row) = column
               exit
            end if
         end do
      end do

! 3. each row's log g range, and the table's
      tbl%log_g_min = 999.0d0
      tbl%log_g_max = -999.0d0
      do row = 1, n
         first = tbl%first_column(row)
         last = tbl%last_column(row)
         tbl%row_log_g_min(row) = tbl%log_g(first) - &
              low_gravity_widths*(tbl%log_g(first+1) - tbl%log_g(first))
         if (tbl%row_log_g_min(row) < tbl%log_g_min) tbl%log_g_min = tbl%row_log_g_min(row)
         tbl%row_log_g_max(row) = tbl%log_g(last) + (tbl%log_g(last) - tbl%log_g(last-1))
         if (tbl%row_log_g_max(row) > tbl%log_g_max) tbl%log_g_max = tbl%row_log_g_max(row)
      end do

! 4. the check
      table_is_bad = .false.
      do row = 1, n
         do column = tbl%first_column(row), tbl%last_column(row)
            if (tbl%log10_pressure(row,column) < valid_above) then
               table_is_bad = .true.
               call log_allard_bad_entry(pow(10.0d0, tbl%log_teff(row)), tbl%log_g(column))
            end if
         end do
      end do
      do row = 1, n
         first = tbl%first_column(row)
         last = tbl%last_column(row)
         if (last - first + 1 < min_valid_columns) then
            table_is_bad = .true.
            call log_allard_short_row(row, last - first + 1, pow(10.0d0, tbl%log_teff(row)), &
                 tbl%row_log_g_min(row), tbl%row_log_g_max(row))
         end if
      end do
      if (table_is_bad) then
         call report_allard_table_bad()
         call log_allard_tables(tbl)
         ierr = 1
      end if
end subroutine find_table_ranges

! Shell's sort (diminishing increments, Numerical Recipes' shell) of
! values(1:n) into ascending order.
pure subroutine sort_shell(n, values)
      integer, intent(in) :: n
      real(dp), intent(inout) :: values(:)
      integer :: increment, i, j
      real(dp) :: v

      if (n == 1) return
      increment = 1
      do
         increment = 3*increment + 1
         if (increment > n) exit
      end do
      do
         increment = increment/3
         do i = increment+1, n
            v = values(i)
            j = i
            do while (values(j-increment) > v)
               values(j) = values(j-increment)
               j = j - increment
               if (j <= increment) exit
            end do
            values(j) = v
         end do
         if (increment <= 1) exit
      end do
end subroutine sort_shell

!======================================================================
! Sources: the F77's ALFILEIN (the reader, both layouts) and ALTABINIT
! (the limits and the check), M. H. Pinsonneault's Allard atmospheres
! (6/97) and the 2007 NextGen layout (2008); Numerical Recipes' shell
! sort (Press et al., section 8.1). The model atmospheres: F. Allard,
! P. H. Hauschildt and collaborators (NextGen, AMES-Cond, AMES-Dusty).

end module allard_table_reader

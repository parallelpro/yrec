!----------------------------------------------------------------------
! opal95_reader
!----------------------------------------------------------------------
! Reading an OPAL 1995-format table file into an opal95_table_set (the
! family's only I/O; opal95_opacity builds the fixed-Z and surface tables
! from it and looks the opacity up).
!     read_opal95_tables   the file's 126 tables, checked against the set's
!                          X and Z grids
module opal95_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: n_opal95_d, n_opal95_t, n_opal95_tables, n_opal95_z, &
           n_opal95_x, opal95_missing_opacity, opal95_table_set
      use kap_messages, only: kap_report, message_length, unwritten_line
      use opacity_table_files, only: open_table
      implicit none
      private
      public :: read_opal95_tables

contains

!----------------------------------------------------------------------
! read_opal95_tables
!----------------------------------------------------------------------
! Reads the full set of OPAL95 opacity tables (every tabulated Z and X)
! into tbl, after the file's header, checking each table's X and Z
! against the set's grids. The tables come X by X, each X with its Z
! from 0 up (13 Z for the first eight X, 10 for X = 0.95, 12 for
! X = 1 - Z). A file that cannot be opened or read, a table whose X and
! Z differ from the grid's, or a file with fewer than n_opal95_tables
! tables is reported, with ierr 1.
subroutine read_opal95_tables(tbl, opal95_table_path, ierr)
      type(opal95_table_set), intent(inout) :: tbl
      character(len=*), intent(in) :: opal95_table_path
      integer, intent(out) :: ierr

      integer :: table_unit, ios
! a line of the file's header, read until the first TABLE line
      character(len=132) :: header_line
! the table being read: its number in the file, its X and Z indices and
! its slot in tbl%full_opacity; the X and Z its header states, and the X
! its slot is at
      integer :: table_number, ix, iz, slot
      real(dp) :: header_x, header_z, slot_x
      integer :: row, num_filled, num_z_at_x
      character(len=message_length) :: lines(1)

      call open_table(opal95_table_path, 'OPAL 1995', table_unit, ierr)
      if (ierr /= 0) return
! skip the file's header: up to the first table's header line
      do
         read(table_unit,'(a)',iostat=ios) header_line
         if (ios < 0) then
            call ends_inside_table(0)
            return
         else if (ios > 0) then
            call unreadable_header(0)
            return
         end if
         if (header_line(1:5) == 'TABLE') exit
      end do
      read(header_line,'(36X,F7.4,11X,F7.4)',iostat=ios) header_x, header_z
      if (ios /= 0) then
         call unreadable_header(0)
         return
      end if
      iz = 1
      ix = 1
      do table_number = 1, n_opal95_tables
!        the table's X and Z against its slot's
         slot = tbl%table_start_index(iz)+ix
         if (ix < n_opal95_x) then
            slot_x = tbl%grid_x(ix)
         else
            slot_x = 1.0d0 - tbl%grid_z(iz)
         end if
         if (tbl%grid_z(iz) /= header_z .or. slot_x /= header_x) then
            lines = unwritten_line
            write(lines,'(a,i0,a,2(f7.4,a),2(f7.4,a))') 'kap_init: table ', table_number, &
                 ' of the OPAL 1995 table ' // trim(opal95_table_path) // ' is at X =', header_x, &
                 ', Z =', header_z, ' where X =', slot_x, ', Z =', tbl%grid_z(iz), ' is expected'
            call report_and_close(lines)
            return
         end if
!        its log R grid (log R = log rho - 3 log T6)
         read(table_unit,'(///,4x,19f7.1,/)',iostat=ios) tbl%grid_logr
         if (ios < 0) then
            call ends_after_tables(table_number-1)
            return
         else if (ios > 0) then
            call unreadable_header(table_number)
            return
         end if
!        its rows: log T, then log kappa at each log R the row has (the
!        rows from log T 7.2 up end before the last log R; the cells
!        they lack are stored as opal95_missing_opacity)
         do row = 1,n_opal95_t
            num_filled = opal95_row_columns(row)
            read(table_unit,'(f4.2,19f7.3)',iostat=ios) tbl%grid_logt(row), &
                 tbl%full_opacity(slot,row,1:num_filled)
            if (ios < 0) then
               call ends_inside_table(table_number)
               return
            else if (ios > 0) then
               lines = unwritten_line
               write(lines,'(a,i0,2a)') 'kap_init: cannot read a row of table ', table_number, &
                    ' of the OPAL 1995 table ', trim(opal95_table_path)
               call report_and_close(lines)
               return
            end if
            tbl%full_opacity(slot,row,num_filled+1:) = opal95_missing_opacity
         end do
         if (table_number == n_opal95_tables) exit
!        the next table: the next Z of this X, or the first Z of the
!        next X (X = 0.95 has no table from Z = 0.06 up, X = 1 - Z none
!        at Z = 0.1); its header
         if (ix <= 8) then
            num_z_at_x = n_opal95_z
         else if (ix == 9) then
            num_z_at_x = n_opal95_z - 3
         else
            num_z_at_x = n_opal95_z - 1
         end if
         if (iz < num_z_at_x) then
            iz = iz + 1
         else
            iz = 1
            ix = ix + 1
         end if
         read(table_unit,'(/36x,f7.4,11x,f7.4)',iostat=ios) header_x, header_z
         if (ios < 0) then
            call ends_after_tables(table_number)
            return
         else if (ios > 0) then
            call unreadable_header(table_number+1)
            return
         end if
      end do
      close(table_unit)

contains

! the three ways the file can fail to read, each reported with the
! unit closed and ierr = 1
subroutine unreadable_header(number)
      integer, intent(in) :: number
      character(len=message_length) :: lines(1)

      lines = unwritten_line
      write(lines,'(a,i0,a)') 'kap_init: cannot read the header of table ', number, &
           ' of the OPAL 1995 table ' // trim(opal95_table_path)
      call report_and_close(lines)
end subroutine unreadable_header

subroutine ends_inside_table(number)
      integer, intent(in) :: number
      character(len=message_length) :: lines(1)

      lines = unwritten_line
      write(lines,'(a,i0,a)') 'kap_init: the OPAL 1995 table ' // trim(opal95_table_path) // &
           ' ends inside table ', number
      call report_and_close(lines)
end subroutine ends_inside_table

subroutine ends_after_tables(number)
      integer, intent(in) :: number
      character(len=message_length) :: lines(1)

      lines = unwritten_line
      write(lines,'(a,2(i0,a))') 'kap_init: the OPAL 1995 table ' // trim(opal95_table_path) // &
           ' ends after ', number, ' of its ', n_opal95_tables, ' tables'
      call report_and_close(lines)
end subroutine ends_after_tables

subroutine report_and_close(lines)
      character(len=*), intent(in) :: lines(:)

      call kap_report(lines)
      close(table_unit)
      ierr = 1
end subroutine report_and_close

end subroutine read_opal95_tables

! The number of log R columns a table's row has, of n_opal95_d (19):
! all of them up to the 57th row (log T 7.1), then fewer and fewer as
! the rows reach the tables' empty high-R, high-T corner (18 at log T
! 7.2, 14 at the last row, log T 8.7).
pure integer function opal95_row_columns(row) result(num_columns)
      integer, intent(in) :: row

      select case (row)
      case (:57)
         num_columns = n_opal95_d
      case (58)
         num_columns = n_opal95_d - 1
      case (59:60)
         num_columns = n_opal95_d - 2
      case (61:64)
         num_columns = n_opal95_d - 3
      case (65:69)
         num_columns = n_opal95_d - 4
      case default
         num_columns = n_opal95_d - 5
      end select
end function opal95_row_columns

end module opal95_reader

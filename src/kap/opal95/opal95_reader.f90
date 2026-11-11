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
      use kap_messages, only: kap_log, kap_report, kap_terminal, message_length, &
           unwritten_line
      use opacity_table_files, only: open_table
      implicit none
      private
      public :: read_opal95_tables

contains

!----------------------------------------------------------------------
! read_opal95_tables (formerly ll95tbl)
!----------------------------------------------------------------------
! Reads the full set of OPAL95 opacity tables (every tabulated Z and X)
! into tbl, after the file's header, checking each table's X and Z
! against the set's grids. A file that cannot be
! opened or read, a table whose X and Z differ from the grid's, or a
! file with fewer than n_opal95_tables tables is reported, with ierr 1.
subroutine read_opal95_tables(tbl, opal95_table_path, ierr)
      type(opal95_table_set), intent(inout) :: tbl
      integer, intent(out) :: ierr

! the table file's unit
      integer :: table_unit
      character(len=256), intent(in) :: opal95_table_path

! a line of the file's header, read until the first TABLE line
      character(len=132) :: header_line

      integer :: fmt_start, iz, ix, nn, n, i, j, nmax
      integer :: read_status
      real(dp) :: xx, zz, xxt
      character(len=message_length) :: lines(1)

! skip the header
      call open_table(opal95_table_path, 'OPAL 1995', table_unit, ierr)
      if (ierr /= 0) return
      fmt_start=1
! the table number the messages name (none until the first TABLE line)
      nn = 0

      do
         read(table_unit,1,end=9999) header_line
      if (header_line(fmt_start:fmt_start+4) == 'TABLE') exit
      end do
    1 format(a)
      read(header_line,'(36X,F7.4,11X,F7.4)',iostat=read_status) xx,zz
      if (read_status /= 0) goto 9997
      iz = 1
      ix = 1
      nn = 1
! one table per pass
      table_loop: do
      n = tbl%table_start_index(iz)+ix
      if (ix < n_opal95_x) then
         xxt = tbl%grid_x(ix)
      else
         xxt = 1.0d0 - tbl%grid_z(iz)
      end if
      if (tbl%grid_z(iz) /= zz .or. xxt /= xx) then
       lines = unwritten_line
       write(lines,*)' OPAL95: Z ERROR INCOMPATIBLE TABLE'
       call kap_log(lines)
         close(table_unit)
         ierr = 1
         return
      end if

!     the table's log R grid (log R = log rho - 3 log T6)
      read(table_unit,20,iostat=read_status) (tbl%grid_logr(i),i=1,n_opal95_d)
      if (read_status < 0) then
         lines = unwritten_line
         write(lines,*) 'read_opal95_tables: OPAL95 opacity table file ended after ', nn-1, &
              ' of ', n_opal95_tables, ' tables'
         call kap_terminal(lines)
         close(table_unit)
         ierr = 1
         return
      end if
      if (read_status > 0) goto 9997
   20 format(///,4x,19f7.1,/)
!     the table: log kappa in log T and log R; the rows from the 58th have
!     no data at the highest log R, stored as opal95_missing_opacity
      do i = 1,57
         read(table_unit,30,end=9999,err=9998) tbl%grid_logt(i), &
              (tbl%full_opacity(n,i,j),j=1,n_opal95_d)
      end do
   30 format(f4.2,19f7.3)
      read(table_unit,31,end=9999,err=9998) tbl%grid_logt(58), &
           (tbl%full_opacity(n,58,j),j=1,n_opal95_d-1)
   31 format(f4.2,18f7.3)

      tbl%full_opacity(n,58,n_opal95_d) = opal95_missing_opacity
      do i = 59,60
         read(table_unit,32,end=9999,err=9998) tbl%grid_logt(i), &
              (tbl%full_opacity(n,i,j),j=1,n_opal95_d-2)
   32    format(f4.2,17f7.3)
         do j = n_opal95_d-1,n_opal95_d
            tbl%full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      do i = 61,64
         read(table_unit,33,end=9999,err=9998) tbl%grid_logt(i), &
              (tbl%full_opacity(n,i,j),j=1,n_opal95_d-3)
   33    format(f4.2,16f7.3)
         do j = n_opal95_d-2,n_opal95_d
            tbl%full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      do i = 65,69
         read(table_unit,34,end=9999,err=9998) tbl%grid_logt(i), &
              (tbl%full_opacity(n,i,j),j=1,n_opal95_d-4)
   34    format(f4.2,15f7.3)
         do j = n_opal95_d-3,n_opal95_d
            tbl%full_opacity(n,i,j) = opal95_missing_opacity
         end do
      end do
      i = n_opal95_t
      read(table_unit,35,end=9999,err=9998) tbl%grid_logt(i), &
           (tbl%full_opacity(n,i,j),j=1,n_opal95_d-5)
   35 format(f4.2,14f7.3)
      do j = n_opal95_d-4,n_opal95_d
         tbl%full_opacity(n,i,j) = opal95_missing_opacity
      end do

!     Exit if correct number of tables read in.
      if (nn < n_opal95_tables) then
      nn = nn + 1
!     Need to account for fewer X values at higher Z.
      if (ix <= 8) then
         nmax = n_opal95_z
      else if (ix == 9) then
         nmax = n_opal95_z - 3
      else
         nmax = n_opal95_z - 1
      end if
!     Read in next metal abundance (at fixed X) or read in first metal
!     abundance at next X.
      if (iz < nmax) then
         iz = iz + 1
      else
         iz = 1
         ix = ix + 1
      end if

!     Return to read in next table (L10's gap: a file that ends after a
!     whole table is a short file, not a runtime error).
      read(table_unit,900,iostat=read_status) xx,zz
  900 format(/36x,f7.4,11x,f7.4)
      if (read_status < 0) then
         lines = unwritten_line
         write(lines,*) 'read_opal95_tables: OPAL95 opacity table file ended after ', nn-1, &
              ' of ', n_opal95_tables, ' tables'
         call kap_terminal(lines)
         close(table_unit)
         ierr = 1
         return
      end if
      if (read_status > 0) goto 9997
      cycle table_loop
      end if
      exit table_loop
      end do table_loop

      close(table_unit)
      return
!     A table header that does not read.
 9997 lines = unwritten_line
      write(lines,*) 'read_opal95_tables: malformed OPAL95 opacity table (header read failed)'
      call kap_terminal(lines)
      close(table_unit)
      ierr = 1
      return
!     A table row that does not read.
 9998 lines = unwritten_line
      write(lines,'(a,i0,2a)') 'kap_init: cannot read a row of table ', nn, &
           ' of the OPAL 1995 table ', trim(opal95_table_path)
      call kap_report(lines)
      close(table_unit)
      ierr = 1
      return
!     Premature end of file inside a table body (end=9999 above).
 9999 lines = unwritten_line
      write(lines,*) 'read_opal95_tables: OPAL95 opacity table file ended inside table ', nn
      call kap_terminal(lines)
      close(table_unit)
      ierr = 1
end subroutine read_opal95_tables

end module opal95_reader

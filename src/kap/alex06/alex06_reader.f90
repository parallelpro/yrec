!----------------------------------------------------------------------
! alex06_reader
!----------------------------------------------------------------------
! Reading the Alexander 2006 (Ferguson et al. 2005) file into an
! alex06_table_set (the family's only I/O; alex06_opacity builds the
! table at each composition from it and looks the opacity up).
!     read_alex06_tables   the file, checked against the set's grids
module alex06_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: alex06_table_set, n_alex06_d, n_alex06_t, n_alex06_x, &
           n_alex06_z
      use kap_messages, only: kap_report, message_length, unwritten_line
      use opacity_table_files, only: open_table, read_failed, iomsg_length
      implicit none
      private
      public :: read_alex06_tables

contains

!----------------------------------------------------------------------
! read_alex06_tables
!----------------------------------------------------------------------
! Reads the file into tbl%full_opacity, each table with read_alex06_table,
! which checks its header (X, Z and the log R columns) and each row's
! log T against the set's grids; a mismatch returns ierr = 1. The file
! holds, for each X = 0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8, 0.9 and then
! 1 - Z, the 16 Z = 0, 1e-5, 3e-5, 1e-4, 3e-4, 0.001, 0.002, 0.004,
! 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.08, 0.1 (15 for X = 1 - Z: the
! X = 0.9 table is Z = 0.1's), each in the 85 log T from 4.5 down to 2.7
! (0.05 dex above 3.5 and below 2.9, 0.01 dex between) and the 19
! log R = log rho - 3 log T6 from -8 to 1. The table of X index ix and Z
! index iz goes to slot n_alex06_z (ix - 1) + iz.
subroutine read_alex06_tables(tbl, alex06_table_path, ierr)
      type(alex06_table_set), intent(inout) :: tbl
      character(len=*), intent(in) :: alex06_table_path
      integer, intent(out) :: ierr

      integer :: table_unit, ix, iz

      call open_table(alex06_table_path, 'Alexander 2006', table_unit, ierr)
      if (ierr /= 0) return
      do ix = 1,n_alex06_x-1
         do iz = 1,n_alex06_z
            call read_alex06_table(tbl, table_unit, alex06_table_path, (ix-1)*n_alex06_z+iz, &
                 tbl%grid_x(ix), tbl%grid_z(iz), ierr)
            if (ierr /= 0) return
         end do
      end do
!     the X = 1 - Z tables, one per Z but the last (Z = 0.1's is the
!     X = 0.9 table)
      do iz = 1,n_alex06_z-1
         call read_alex06_table(tbl, table_unit, alex06_table_path, (n_alex06_x-1)*n_alex06_z+iz, &
              1.0d0 - tbl%grid_z(iz), tbl%grid_z(iz), ierr)
         if (ierr /= 0) return
      end do
      close(unit=table_unit)
!     Initialize fixed Z,X table, and the cache that names its
!     composition (X = 0, Z = 0: a second kap_init starts from it too)
      tbl%opacity = tbl%full_opacity(1,:,:)
      tbl%cached_x = 0.0d0
      tbl%cached_z = 0.0d0
end subroutine read_alex06_tables

! One table of the file, into tbl%full_opacity(slot,:,:): its header,
! which must state expected_x and expected_z; its log R columns, which
! must be the set's; and its rows from the highest log T down, each of
! which must be at the set's log T. A mismatch or a failed read is
! reported, with the unit closed and ierr = 1.
subroutine read_alex06_table(tbl, table_unit, path, slot, expected_x, expected_z, ierr)
      type(alex06_table_set), intent(inout) :: tbl
      integer, intent(in) :: table_unit, slot
      character(len=*), intent(in) :: path
      real(dp), intent(in) :: expected_x, expected_z
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success

      real(dp) :: header_x, header_z, header_logr(n_alex06_d), row_logt
      integer :: row, column, ios
      character(len=iomsg_length) :: iomsg
      character(len=message_length) :: lines(1)

      read(table_unit,'(30x,f9.6,7x,f9.6//)',iostat=ios,iomsg=iomsg) header_x, header_z
      if (read_failed(ios, iomsg, 'Alexander 2006', path, table_unit, ierr)) return
      if (header_x /= expected_x .or. header_z /= expected_z) then
         lines = unwritten_line
         write(lines,15) trim(path), expected_x, header_x, expected_z, header_z
   15    format('kap_init: the Alexander 2006 table ',a,' expects X =',f7.4, &
              ' and has',f7.4,', expects Z =',f7.4,' and has',f7.4)
         call refuse(lines)
         return
      end if
      read(table_unit,'(6x,19f7.3)',iostat=ios,iomsg=iomsg) header_logr
      if (read_failed(ios, iomsg, 'Alexander 2006', path, table_unit, ierr)) return
      do column = 1,n_alex06_d
         if (header_logr(column) /= tbl%grid_logr(column)) then
            lines = unwritten_line
            write(lines,25) trim(path), tbl%grid_logr(column), column, header_logr(column)
   25       format('kap_init: the Alexander 2006 table ',a,' expects log R',f7.3, &
                 ' in column ',i0,' and has',f7.3)
            call refuse(lines)
            return
         end if
      end do
      do row = n_alex06_t,1,-1
         read(table_unit,'(f5.3,1x,19f7.3)',iostat=ios,iomsg=iomsg) row_logt, &
              tbl%full_opacity(slot,row,:)
         if (read_failed(ios, iomsg, 'Alexander 2006', path, table_unit, ierr)) return
         if (row_logt /= tbl%grid_logt(row)) then
            lines = unwritten_line
            write(lines,35) trim(path), tbl%grid_logt(row), row, row_logt
   35       format('kap_init: the Alexander 2006 table ',a,' expects log T',f7.3, &
                 ' in row ',i0,' and has',f7.3)
            call refuse(lines)
            return
         end if
      end do

contains

subroutine refuse(lines)
      character(len=*), intent(in) :: lines(:)

      call kap_report(lines)
      close(table_unit)
      ierr = 1
end subroutine refuse

end subroutine read_alex06_table

end module alex06_reader

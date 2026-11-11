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
      use kap_messages, only: kap_terminal, message_length, unwritten_line
      use opacity_table_files, only: open_table, read_failed
      implicit none
      private
      public :: read_alex06_tables

contains

!----------------------------------------------------------------------
! read_alex06_tables (formerly readalex06)
!----------------------------------------------------------------------
! Reads the file into tbl%full_opacity, checking each table's header (X,
! Z and the log R columns) and each row's log T against the set's grids;
! a mismatch returns ierr = 1. The file holds, for each X = 0, 0.1, 0.2,
! 0.35, 0.5, 0.7, 0.8, 0.9 and then 1 - Z, the 16 Z = 0, 1e-5, 3e-5,
! 1e-4, 3e-4, 0.001, 0.002, 0.004, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06,
! 0.08, 0.1 (15 for X = 1 - Z: the X = 0.9 table is Z = 0.1's), each in
! the 85 log T from 4.5 down to 2.7 (0.05 dex above 3.5 and below 2.9,
! 0.01 dex between) and the 19 log R = log rho - 3 log T6 from -8 to 1.
subroutine read_alex06_tables(tbl, alex06_table_path, ierr)
      type(alex06_table_set), intent(inout) :: tbl
      integer, intent(out) :: ierr

! the table file's unit
      integer :: table_unit
      character(len=256), intent(in) :: alex06_table_path

! a table's header X, Z and log R columns, and a row's log T, as read
      real(dp) :: row_logr_check(n_alex06_d), header_x, header_z, row_temp
      integer :: i, ii, j, jj, k, kk
      character(len=message_length) :: lines(2)
      integer :: ios
      character(len=256) :: iomsg

      call open_table(alex06_table_path, 'Alexander 2006', table_unit, ierr)
      if (ierr /= 0) return
!     Read in initial X and Z; ensure that they have the expected values.
      do i = 1,n_alex06_x-1
         do ii = 1, n_alex06_z
!           the table of X index i and Z index ii goes to slot
!           n_alex06_z (i - 1) + ii; its header gives X and Z
            read(table_unit,10,iostat=ios,iomsg=iomsg) header_x, header_z
            if (read_failed(ios, iomsg, 'Alexander 2006', alex06_table_path, table_unit, ierr)) return
   10       format(30x,f9.6,7x,f9.6//)
            if (header_x /= tbl%grid_x(i) .or. header_z /= tbl%grid_z(ii)) then
               lines = unwritten_line
               write(lines,15) tbl%grid_x(i), header_x, tbl%grid_z(ii), header_z
               call kap_terminal(lines)
   15          format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                    1x,'EXPECTED AND ACTUAL X,Z',4f7.2,' RUN STOPPED')
               close(table_unit)
               ierr = 1
               return
            end if
            read(table_unit,20,iostat=ios,iomsg=iomsg) (row_logr_check(k),k=1,n_alex06_d)
            if (read_failed(ios, iomsg, 'Alexander 2006', alex06_table_path, table_unit, ierr)) return
   20       format(6x,19f7.3)
            do kk=1,n_alex06_d
               if (row_logr_check(kk) /= tbl%grid_logr(kk)) then
                  lines = unwritten_line
                  write(lines,25) kk, tbl%grid_logr(kk), row_logr_check(kk)
                  call kap_terminal(lines)
   25             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL R',i3,2f7.3,' RUN STOPPED')
                  close(table_unit)
                  ierr = 1
                  return
               end if
            end do
            jj = (i-1)*n_alex06_z+ii
!           its rows, from the highest log T down, each checked against
!           the grid's log T
            do j = n_alex06_t,1,-1
               read(table_unit,30,iostat=ios,iomsg=iomsg) row_temp, (tbl%full_opacity(jj,j,k),k=1,n_alex06_d)
               if (read_failed(ios, iomsg, 'Alexander 2006', alex06_table_path, table_unit, ierr)) return
   30          format(f5.3,1x,19f7.3)
               if (row_temp /= tbl%grid_logt(j)) then
                  lines = unwritten_line
                  write(lines,35) j, tbl%grid_logt(j), row_temp
                  call kap_terminal(lines)
   35             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL T',i3,2f7.3,' RUN STOPPED')
                  close(table_unit)
                  ierr = 1
                  return
               end if
            end do
         end do
      end do
!     the X = 1 - Z tables, one per Z but the last (Z = 0.1's is the
!     X = 0.9 table)
      do ii = 1, n_alex06_z-1
!        slot n_alex06_z (n_alex06_x - 1) + ii; its header gives X and Z
         read(table_unit,10,iostat=ios,iomsg=iomsg) header_x, header_z
         if (read_failed(ios, iomsg, 'Alexander 2006', alex06_table_path, table_unit, ierr)) return
         if (header_x /= (1.0d0 - tbl%grid_z(ii)) .or. &
              header_z /= tbl%grid_z(ii)) then
            lines = unwritten_line
            write(lines,15) 1.0d0 - tbl%grid_z(ii), header_x, tbl%grid_z(ii), header_z
            call kap_terminal(lines)
            close(table_unit)
            ierr = 1
            return
         end if
         read(table_unit,20,iostat=ios,iomsg=iomsg) (row_logr_check(k),k=1,n_alex06_d)
         if (read_failed(ios, iomsg, 'Alexander 2006', alex06_table_path, table_unit, ierr)) return
         do kk=1,n_alex06_d
            if (row_logr_check(kk) /= tbl%grid_logr(kk)) then
               lines = unwritten_line
               write(lines,25) kk, tbl%grid_logr(kk), row_logr_check(kk)
               call kap_terminal(lines)
               close(table_unit)
               ierr = 1
               return
            end if
         end do
         jj = (n_alex06_x-1)*n_alex06_z+ii
!        its rows, checked as above
         do j = n_alex06_t,1,-1
            read(table_unit,30,iostat=ios,iomsg=iomsg) row_temp, (tbl%full_opacity(jj,j,k),k=1,n_alex06_d)
            if (read_failed(ios, iomsg, 'Alexander 2006', alex06_table_path, table_unit, ierr)) return
            if (row_temp /= tbl%grid_logt(j)) then
               lines = unwritten_line
               write(lines,35) j, tbl%grid_logt(j), row_temp
               call kap_terminal(lines)
               close(table_unit)
               ierr = 1
               return
            end if
         end do
      end do
      close(unit=table_unit)
!     Initialize fixed Z,X table, and the cache that names its
!     composition (X = 0, Z = 0: a second kap_init starts from it too)
      do i = 1,n_alex06_t
         do j = 1,n_alex06_d
            tbl%opacity(i,j)=tbl%full_opacity(1,i,j)
         end do
      end do
      tbl%cached_x = 0.0d0
      tbl%cached_z = 0.0d0
end subroutine read_alex06_tables

end module alex06_reader

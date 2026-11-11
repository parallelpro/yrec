!----------------------------------------------------------------------
! alex94_reader
!----------------------------------------------------------------------
! Reading the seven Alexander & Ferguson (1994) files into an
! alex94_table_set (the family's only I/O; alex94_opacity builds the
! fixed-Z and surface tables from them and looks the opacity up).
!     read_alex94_tables   the seven files, checked against the set's
!                          grids
module alex94_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: alex94_table_set, n_alex94_d, n_alex94_t, n_alex94_x, &
           n_alex94_z
      use kap_messages, only: kap_log, message_length, unwritten_line
      use opacity_table_files, only: open_table, read_failed
      implicit none
      private
      public :: read_alex94_tables

contains

!----------------------------------------------------------------------
! read_alex94_tables
!----------------------------------------------------------------------
! Reads the seven files (one per X = 0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.8)
! into tbl%full_opacity, checking each table's X, Z, row count, log T
! and first log R against the set's grids (a mismatch returns ierr = 1).
! Each file holds the 15
! Z = 0.1 down to 0, each in the 23 log T from 4.10 down to 3.00 and the
! 17 log R = log rho - 3 log T6 from -7 to 1.
subroutine read_alex94_tables(tbl, alex95_table_paths, ierr)
      type(alex94_table_set), intent(inout) :: tbl
      integer, intent(out) :: ierr

! the table file's unit
      integer :: table_unit
      character(len=256), intent(in) :: alex95_table_paths(:)

! a row as read (its index goes to tbl%index_t, as in the F77 code, so
! the last row read seeds the T cursor)
      real(dp) :: row_opacity_temp(n_alex94_d)
      integer :: table_index, i, ii, j, k, row_density_count
      real(dp) :: header_x, header_z, row_temp, row_logr0
      character(len=message_length) :: lines(2)
      integer :: ios
      character(len=256) :: iomsg

      do table_index=1,n_alex94_x
!        the file of the table_index-th X, with every Z
         call open_table(alex95_table_paths(table_index), 'Alexander 1994', table_unit, ierr)
         if (ierr /= 0) return
         do i = n_alex94_z,1,-1
!           the table of X index table_index and Z index i goes to slot
!           table_index + 7 (i - 1)
            ii = table_index + 7*(i-1)
!           Header information: X and Z
            read(table_unit,10,iostat=ios,iomsg=iomsg) header_x, header_z
            if (read_failed(ios, iomsg, 'Alexander 1994', alex95_table_paths(table_index), table_unit, ierr)) return
   10       format(18x,f6.2,2x,f7.2)
            if (header_x /= tbl%grid_x(table_index) .or. &
                 header_z /= tbl%grid_z(i)) then
               lines = unwritten_line
               write(lines,20) header_x, tbl%grid_x(table_index), &
                    header_z, tbl%grid_z(i)
               call kap_log(lines)
   20          format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                    1x,'EXPECTED AND ACTUAL X,Z',4f7.2,' RUN STOPPED')
               close(table_unit)
               ierr = 1
               return
            end if
!           its rows, from the highest log T down, each checked against
!           the grid's log T and first log R
            do j = n_alex94_t,1,-1
               read(table_unit,30,iostat=ios,iomsg=iomsg) tbl%index_t, row_density_count, &
                    row_temp, row_logr0, (row_opacity_temp(k),k=1,n_alex94_d)
               if (read_failed(ios, iomsg, 'Alexander 1994', alex95_table_paths(table_index), table_unit, ierr)) return
               do k = 1, n_alex94_d
                  tbl%full_opacity(ii,j,k) = row_opacity_temp(k)
               end do
   30          format(i2,i3,f6.3,f5.1,8f8.3/9f8.3)
               if (row_density_count /= n_alex94_d .or. row_temp /= tbl%grid_logt(j) &
                    .or. row_logr0 /= tbl%grid_logr(1)) then
                  lines = unwritten_line
                  write(lines,40) row_density_count, row_temp, &
                       tbl%grid_logt(j), row_logr0, tbl%grid_logr(1)
                  call kap_log(lines)
   40             format(1x,'ERROR IN ALEXANDER OPACITY TABLES:'/ &
                       1x,'EXPECTED AND ACTUAL T,RHO',i3,4f7.2, &
                       ' RUN STOPPED')
                  close(table_unit)
                  ierr = 1
                  return
               end if
            end do
         end do
         close(unit=table_unit)
      end do
      ierr = 0
end subroutine read_alex94_tables

end module alex94_reader

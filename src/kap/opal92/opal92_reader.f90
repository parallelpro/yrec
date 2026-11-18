!----------------------------------------------------------------------
! opal92_reader
!----------------------------------------------------------------------
! Reading OPAL 1992 table files into opal92_table_sets (the family's
! only I/O; opal92_opacity splines them, builds the surface tables and
! looks the opacity up).
!     read_opal92_tables   the table at opal_table_z1 and, with
!                          two_z_tables, the one at opal_table_z2
module opal92_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use opacity_tables, only: n_opal92_d, n_opal92_t, n_opal92_x, opal92_table_set
      use kap_messages, only: kap_report, message_length, unwritten_line
      use opacity_table_files, only: open_table, read_failed, table_z_differs, iomsg_length
      use math_lib
      implicit none
      private
      public :: read_opal92_tables

contains

!----------------------------------------------------------------------
! read_opal92_tables
!----------------------------------------------------------------------
! Reads the OPAL92 (Lawrence Livermore) opacity table(s) with
! read_opal92_table: the table at opal_table_z1 and, with
! two_z_tables, the one at opal_table_z2. A table that cannot be opened,
! is short, or states another Z than its control's (table_z) returns
! ierr.
subroutine read_opal92_tables(tables, opal92_table_path, opal92_table2_path, &
     two_z_tables, table_z, ierr)
! (1) at opal_table_z1 and, with two_z_tables, (2) at opal_table_z2
      type(opal92_table_set), intent(inout) :: tables(:)
      character(len=*), intent(in) :: opal92_table_path, opal92_table2_path
      logical, intent(in) :: two_z_tables
! opal_table_z1 and opal_table_z2
      real(dp), intent(in) :: table_z(:)
      integer, intent(out) :: ierr

      ierr = 0
      call read_opal92_table(opal92_table_path, tables(1), table_z(1), 'opal_table_z1', ierr)
      if (ierr /= 0) return
      if (two_z_tables) then
         call read_opal92_table(opal92_table2_path, tables(2), table_z(2), 'opal_table_z2', ierr)
      end if
end subroutine read_opal92_tables

! One table file, at path, into tbl: n_opal92_x X tables of n_opal92_t
! temperature rows each (a shorter file is an error), each at the Z of
! its control (control_z, of control_name). The rows give T6 = T/1e6 K,
! stored as log10 T6.
subroutine read_opal92_table(path, tbl, control_z, control_name, ierr)
      character(len=*), intent(in) :: path
      type(opal92_table_set), intent(inout) :: tbl
      real(dp), intent(in) :: control_z
      character(len=*), intent(in) :: control_name
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
      integer :: unit, ix, row, ios
      real(dp) :: header_z, t6
      character(len=iomsg_length) :: iomsg
      character(len=message_length) :: lines(1)

      call open_table(path, 'OPAL 1992', unit, ierr)
      if (ierr /= 0) return
      do ix = 1,n_opal92_x
!        the X table's header (X and Z)
         read(unit,'(33x,f7.4,2x,f7.4)',iostat=ios,iomsg=iomsg) tbl%grid_x(ix), header_z
         if (read_failed(ios, iomsg, 'OPAL 1992', path, unit, ierr)) return
         if (table_z_differs(header_z, control_z, control_name, 'OPAL 1992', path, &
              unit, ierr)) return
         read(unit,'()',iostat=ios,iomsg=iomsg)
         if (read_failed(ios, iomsg, 'OPAL 1992', path, unit, ierr)) return
!        its log R grid (log R = log rho - 3 log T6)
         read(unit,'(6x,17f7.1)',iostat=ios,iomsg=iomsg) tbl%grid_logr
         if (read_failed(ios, iomsg, 'OPAL 1992', path, unit, ierr)) return
!        its rows: T6 and log kappa at each log R
         do row = 1,n_opal92_t
            read(unit,'(18f7.3)',iostat=ios,iomsg=iomsg) t6, &
                 tbl%log10_opacity(row+(ix-1)*n_opal92_t,:)
            if (read_failed(ios, iomsg, 'OPAL 1992', path, unit, ierr)) return
            tbl%grid_logt(row) = log10(t6)
         end do
      end do
      close(unit,iostat=ios)
      if (ios /= 0) then
         lines = unwritten_line
         write(lines,'(a)') 'kap_init: cannot close the OPAL 1992 table ' // trim(path)
         call kap_report(lines)
         ierr = 1
         return
      end if
      tbl%num_temps = n_opal92_t
      tbl%num_x = n_opal92_x
end subroutine read_opal92_table

end module opal92_reader

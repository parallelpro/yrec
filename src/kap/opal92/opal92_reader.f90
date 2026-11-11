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
      use opacity_table_files, only: open_table, read_failed
      use math_lib
      implicit none
      private
      public :: read_opal92_tables

contains

!----------------------------------------------------------------------
! read_opal92_tables
!----------------------------------------------------------------------
! Reads the OPAL92 (Lawrence Livermore) opacity table(s). Each table is read by the internal read_opal92_table, first the
! table at opal_table_z1 and then, with use_two_z_tables, the one at
! opal_table_z2; a table that cannot be opened or is short returns ierr.
subroutine read_opal92_tables(tables, opal92_table_path, opal92_table2_path, &
     two_z_tables, ierr)
! (1) at opal_table_z1 and, with two_z_tables, (2) at opal_table_z2
      type(opal92_table_set), intent(inout) :: tables(:)
      character(len=256), intent(in) :: opal92_table_path, opal92_table2_path
      logical, intent(in) :: two_z_tables
      integer, intent(out) :: ierr
!     Local_grid_z is a read-list target that is never used afterwards.
      real(dp) :: local_grid_z(n_opal92_x)
      integer :: i, k, density_index, ios
      real(dp) :: grid_temp_k
      character(len=256) :: iomsg
      character(len=message_length) :: lines(1)

      ierr = 0
      call read_opal92_table(opal92_table_path, tables(1), ierr)
      if (ierr /= 0) return

      if (two_z_tables) then
         call read_opal92_table(opal92_table2_path, tables(2), ierr)
         if (ierr /= 0) return
      end if
      return

contains

! read_opal92_table: read the table at path into tbl: n_opal92_x X
! tables of n_opal92_t temperature rows each (a shorter file is an error).
subroutine read_opal92_table(path, tbl, ierr)
      character(len=256), intent(in) :: path
      type(opal92_table_set), intent(inout) :: tbl
      integer, intent(inout) :: ierr   ! left at the caller's 0 on success
      integer :: unit

      call open_table(path, 'OPAL 1992', unit, ierr)
      if (ierr /= 0) return
      do i=1,n_opal92_x
!        the X table's header (X and Z)
        read(unit,190,iostat=ios,iomsg=iomsg) tbl%grid_x(i), local_grid_z(i)
        if (read_failed(ios, iomsg, 'OPAL 1992', path, unit, ierr)) return
  190   format(33x,f7.4,2x,f7.4)
         read(unit,'()',iostat=ios,iomsg=iomsg)
         if (read_failed(ios, iomsg, 'OPAL 1992', path, unit, ierr)) return
!        its log R grid (log R = log rho - 3 log T6)
            read(unit, 200, iostat=ios, iomsg=iomsg) (tbl%grid_logr(density_index), density_index=1, n_opal92_d)
            if (read_failed(ios, iomsg, 'OPAL 1992', path, unit, ierr)) return
  200   format (6x, 17f7.1)
!        its rows: T6 and log kappa at each log R
         do k=1, n_opal92_t
         read(unit,196,iostat=ios,iomsg=iomsg) grid_temp_k, &
              (tbl%log10_opacity(k+(i-1)*n_opal92_t,density_index),density_index=1,n_opal92_d)
         if (read_failed(ios, iomsg, 'OPAL 1992', path, unit, ierr)) return
         tbl%grid_logt(k)=log10(grid_temp_k)
         end do
  196    format(18f7.3)
!
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

end subroutine read_opal92_tables

end module opal92_reader

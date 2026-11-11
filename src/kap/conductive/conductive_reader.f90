!----------------------------------------------------------------------
! conductive_reader
!----------------------------------------------------------------------
! Reading the Potekhin conductive table (condall06.d) into a
! conductive_table_set (the family's only I/O; conductive_opacity looks
! the conductive opacity up).
!     read_conductive_table   the table's log conductivities on its
!                             log T, log rho and log Z grid
module conductive_reader
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      use opacity_tables, only: cond_n_rho, cond_n_temp, cond_n_z, conductive_table_set
      use kap_messages, only: kap_terminal, message_length, unwritten_line
      use opacity_table_files, only: read_failed
      implicit none
      private
      public :: read_conductive_table

contains

!----------------------------------------------------------------------
! read_conductive_table
!----------------------------------------------------------------------
! Read the Potekhin table (condall06.d) at path into tbl, once per run
! (kap_init). A table that cannot be opened or is short returns ierr.
subroutine read_conductive_table(tbl, path, ierr)
      type(conductive_table_set), intent(inout) :: tbl
      character(len=*), intent(in) :: path
      integer, intent(out) :: ierr
      integer :: file_unit, z_index, t_index, r_index
      real(dp) :: z_grid_value
      character(len=message_length) :: lines(1)
      integer :: ios
      character(len=256) :: iomsg

      open(newunit=file_unit, file=path, status='OLD', iostat=ierr)
      if (ierr /= 0) then
         lines = unwritten_line
         write(lines,*) 'read_conductive_table: cannot open ', trim(path)
         call kap_terminal(lines)
         ierr = 1
         return
      end if
      read(file_unit,'(A)',iostat=ios,iomsg=iomsg) ! skip the first line
      if (read_failed(ios, iomsg, 'conductive', path, file_unit, ierr)) return
      do z_index=1,cond_n_z
         read(file_unit,*,iostat=ios,iomsg=iomsg) z_grid_value,(tbl%temp_grid(t_index),t_index=1,cond_n_temp)
         if (read_failed(ios, iomsg, 'conductive', path, file_unit, ierr)) return
         tbl%z_grid(z_index)=log10(z_grid_value)
         do r_index=1,cond_n_rho
            read(file_unit,*,iostat=ios,iomsg=iomsg) tbl%rho_grid(r_index), &
                 (tbl%log10_kappa_table(t_index,r_index,z_index),t_index=1,cond_n_temp)
            if (read_failed(ios, iomsg, 'conductive', path, file_unit, ierr)) return
         end do
      end do
      close(file_unit)
      tbl%loaded = .true.
end subroutine read_conductive_table

end module conductive_reader

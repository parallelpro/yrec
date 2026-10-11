!----------------------------------------------------------------------
! conductive_opacity_table
!----------------------------------------------------------------------
! Added 2026 (save-migration campaign): the Potekhin conductive-
! opacity table, promoted from SAVEd locals of
! kap/conductive/condopacp.f90 -- the table was read once into static
! locals on first call, invisible state that crashed when the blanket
! save was removed. Now it is named, resettable domain state like
! every other opacity table. The grid dimensions must stay consistent
! with the condall06.d table file (see condopacp.f90).
module conductive_opacity_table
      use math_lib
      implicit none
      private
      public :: read_conductive_table
      integer, parameter, public :: cond_n_temp = 19, cond_n_rho = 64, &
           cond_n_z = 15
! table_loaded_flag holds this value once condopacp has read the table.
      integer, parameter, public :: cond_table_loaded_marker = 12345

      type, public :: conductive_table_state
            double precision :: temp_grid(cond_n_temp)
            double precision :: rho_grid(cond_n_rho)
            double precision :: z_grid(cond_n_z)
            double precision :: log10_kappa_table(cond_n_temp, cond_n_rho, &
                 cond_n_z)
            integer :: table_loaded_flag = -1
      end type conductive_table_state

      type(conductive_table_state), public, save :: cond_table

contains

! read_conductive_table: read the Potekhin table (condall06.d) at path,
! once per run from kap_init. (2026-10: read_controls opened it on a
! unit and condopacp read it on its first lookup.)
subroutine read_conductive_table(path, ierr)
      character(len=*), intent(in) :: path
      integer, intent(out) :: ierr
      integer :: file_unit, z_index, t_index, r_index
      double precision :: z_grid_value

      open(newunit=file_unit, file=path, status='OLD', iostat=ierr)
      if (ierr /= 0) then
         write(*,*) 'read_conductive_table: cannot open ', trim(path)
         ierr = 1
         return
      end if
      read(file_unit,'(A)') ! skip the first line
      do z_index=1,cond_n_z
         read(file_unit,*) z_grid_value,(cond_table%temp_grid(t_index),t_index=1,cond_n_temp)
         cond_table%z_grid(z_index)=log10(z_grid_value)
         do r_index=1,cond_n_rho
            read(file_unit,*) cond_table%rho_grid(r_index), &
                 (cond_table%log10_kappa_table(t_index,r_index,z_index),t_index=1,cond_n_temp)
         enddo
      enddo
      close(file_unit)
      cond_table%table_loaded_flag=cond_table_loaded_marker
end subroutine read_conductive_table

end module conductive_opacity_table

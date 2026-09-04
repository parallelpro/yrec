!----------------------------------------------------------------------
! conductive_table_lib
!----------------------------------------------------------------------
! Added 2026 (save-migration campaign): the Potekhin conductive-
! opacity table, promoted from SAVEd locals of
! kap/conductive/condopacp.f90 -- the table was read once into static
! locals on first call, invisible state that crashed when the blanket
! save was removed. Now it is named, resettable domain state like
! every other opacity table. The grid dimensions must stay consistent
! with the condall06.d table file (see condopacp.f90).
module conductive_table_lib
      implicit none
      private
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

end module conductive_table_lib

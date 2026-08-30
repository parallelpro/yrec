!----------------------------------------------------------------------
! alprint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original alprint.f; only variable names, source form, and comment
! style were updated.
!
! Prints the internal Allard-format atmosphere tables and auxiliary
! information provided by alfilein and verified by altabinit, to the
! .short log file.
subroutine print_allard_tables

      use atm_table_lib
      use luout_lib
      implicit none
      integer, parameter :: nta = 250, nga = 25
      integer :: i, j

      atm_table%allard_teffl_min = atm_table%allard_teffl_grid(1) - &
           (atm_table%allard_teffl_grid(2)-atm_table%allard_teffl_grid(1))
      atm_table%allard_teffl_max = atm_table%allard_teffl_grid(atm_table%allard_num_teff) + &
           (atm_table%allard_teffl_grid(atm_table%allard_num_teff)-atm_table%allard_teffl_grid(atm_table%allard_num_teff-1))

      write(run_log_unit,*)
      write(run_log_unit,*)' **************** Allard Table  ***************'
      write(run_log_unit,*)
      write(run_log_unit,*)'nTeff= ',atm_table%allard_num_teff,' ,nGL=',atm_table%allard_num_gl
      write(run_log_unit,*)'nFeH= ',atm_table%allard_num_feh,' ,nAlpha=',atm_table%allard_num_alpha
      write(run_log_unit,*)
      write(run_log_unit,*)'TEFFLmin= ',atm_table%allard_teffl_min,' ,TEFFLmax= ',atm_table%allard_teffl_max
      write(run_log_unit,*)'GLXmin=',atm_table%allard_gl_min,' GLXmax=',atm_table%allard_gl_max
      write(run_log_unit,900)' PL @ Teff','GL ->',(atm_table%allard_gl_grid(j),j=1,atm_table%allard_num_gl)
  900      format(A10,16X,A,20(5X,F6.2,5X))
      write(run_log_unit,*)
      write(run_log_unit,*) 'PL at Teff,GL'
      write(run_log_unit,*)
      do i = 1, atm_table%allard_num_teff
         write(run_log_unit,910) i,10D0**atm_table%allard_teffl_grid(i), &
              atm_table%allard_gl_index_min(i),atm_table%allard_gl_index_max(i), &
              atm_table%allard_gl_row_min(i),atm_table%allard_gl_row_max(i), &
              (atm_table%allard_log10_pressure(i,j),j=1,atm_table%allard_num_gl)
  910     format(I5,F6.0,2I4,2F6.2,1P20E16.8)
      enddo
      write(run_log_unit,*)
      write(run_log_unit,*) 'PL at Tau=100'
      write(run_log_unit,*)
      do i = 1, atm_table%allard_num_teff
         write(run_log_unit,910) i,10D0**atm_table%allard_teffl_grid(i), &
              atm_table%allard_gl_index_min(i),atm_table%allard_gl_index_max(i), &
              atm_table%allard_gl_row_min(i),atm_table%allard_gl_row_max(i), &
              (atm_table%allard_log10_pressure_tau100(i,j),j=1,atm_table%allard_num_gl)
      enddo
      write(run_log_unit,*)
      write(run_log_unit,*) 'TL at Tau=100'
      write(run_log_unit,*)
      do i = 1, atm_table%allard_num_teff
         write(run_log_unit,910) i,10D0**atm_table%allard_teffl_grid(i), &
              atm_table%allard_gl_index_min(i),atm_table%allard_gl_index_max(i), &
              atm_table%allard_gl_row_min(i),atm_table%allard_gl_row_max(i), &
              (atm_table%allard_log10_temp_tau100(i,j),j=1,atm_table%allard_num_gl)
      enddo
      write(run_log_unit,*)

end subroutine print_allard_tables

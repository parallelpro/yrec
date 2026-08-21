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
subroutine alprint

      use luout_lib
      implicit none
      integer, parameter :: nta = 250, nga = 25


! common/alatm01/: Allard grid arrays and sizes; shared with alfilein,
! altabinit, alsurfp. Naming matches alfilein.f90.
      double precision :: allard_teffl_grid(nta), allard_gl_grid(nga), &
           allard_feh_grid(nga), allard_alpha_grid(nga), &
           allard_log10_pressure(nta,nga), allard_log10_pressure_tau100(nta,nga), &
           allard_log10_temp_tau100(nta,nga)
      logical :: allard_is_old_nextgen
      integer :: allard_num_teff, allard_num_gl, allard_num_feh, allard_num_alpha
      common /alatm01/ allard_teffl_grid, allard_gl_grid, allard_feh_grid, &
           allard_alpha_grid, allard_log10_pressure, allard_log10_pressure_tau100, &
           allard_log10_temp_tau100, allard_is_old_nextgen, allard_num_teff, &
           allard_num_gl, allard_num_feh, allard_num_alpha

! common/alatm02/: shared with altabinit, alsurfp. Naming matches
! altabinit.f90.
      double precision :: allard_gl_row_min(nta), allard_gl_row_max(nta)
      integer :: allard_gl_index_min(nta), allard_gl_index_max(nta)
      double precision :: allard_teffl_min, allard_teffl_max, allard_gl_min, &
           allard_gl_max
      common /alatm02/ allard_gl_row_min, allard_gl_row_max, allard_gl_index_min, &
           allard_gl_index_max, allard_teffl_min, allard_teffl_max, &
           allard_gl_min, allard_gl_max

! common/alatm03/: not used in this file; declared only to preserve
! layout. Naming matches alfilein.f90.
      double precision :: allard_target_feh, allard_target_alpha
      logical :: allard_use_tau100
      integer :: allard_table_unit
      common /alatm03/ allard_target_feh, allard_target_alpha, allard_use_tau100, &
           allard_table_unit

! common/alatm04/: unused everywhere in this file family; placeholders
! preserving the layout. Naming matches alfilein.f90/altabinit.f90.
      double precision :: alatm04_placeholder1, alatm04_placeholder2, &
           alatm04_placeholder3, alatm04_placeholder4
      common /alatm04/ alatm04_placeholder1, alatm04_placeholder2, &
           alatm04_placeholder3, alatm04_placeholder4

      save

      integer :: i, j

      allard_teffl_min = allard_teffl_grid(1) - &
           (allard_teffl_grid(2)-allard_teffl_grid(1))
      allard_teffl_max = allard_teffl_grid(allard_num_teff) + &
           (allard_teffl_grid(allard_num_teff)-allard_teffl_grid(allard_num_teff-1))

      write(short_file_unit,*)
      write(short_file_unit,*)' **************** Allard Table  ***************'
      write(short_file_unit,*)
      write(short_file_unit,*)'nTeff= ',allard_num_teff,' ,nGL=',allard_num_gl
      write(short_file_unit,*)'nFeH= ',allard_num_feh,' ,nAlpha=',allard_num_alpha
      write(short_file_unit,*)
      write(short_file_unit,*)'TEFFLmin= ',allard_teffl_min,' ,TEFFLmax= ',allard_teffl_max
      write(short_file_unit,*)'GLXmin=',allard_gl_min,' GLXmax=',allard_gl_max
      write(short_file_unit,900)' PL @ Teff','GL ->',(allard_gl_grid(j),j=1,allard_num_gl)
  900      format(A10,16X,A,20(5X,F6.2,5X))
      write(short_file_unit,*)
      write(short_file_unit,*) 'PL at Teff,GL'
      write(short_file_unit,*)
      do i = 1, allard_num_teff
         write(short_file_unit,910) i,10D0**allard_teffl_grid(i), &
              allard_gl_index_min(i),allard_gl_index_max(i), &
              allard_gl_row_min(i),allard_gl_row_max(i), &
              (allard_log10_pressure(i,j),j=1,allard_num_gl)
  910     format(I5,F6.0,2I4,2F6.2,1P20E16.8)
      enddo
      write(short_file_unit,*)
      write(short_file_unit,*) 'PL at Tau=100'
      write(short_file_unit,*)
      do i = 1, allard_num_teff
         write(short_file_unit,910) i,10D0**allard_teffl_grid(i), &
              allard_gl_index_min(i),allard_gl_index_max(i), &
              allard_gl_row_min(i),allard_gl_row_max(i), &
              (allard_log10_pressure_tau100(i,j),j=1,allard_num_gl)
      enddo
      write(short_file_unit,*)
      write(short_file_unit,*) 'TL at Tau=100'
      write(short_file_unit,*)
      do i = 1, allard_num_teff
         write(short_file_unit,910) i,10D0**allard_teffl_grid(i), &
              allard_gl_index_min(i),allard_gl_index_max(i), &
              allard_gl_row_min(i),allard_gl_row_max(i), &
              (allard_log10_temp_tau100(i,j),j=1,allard_num_gl)
      enddo
      write(short_file_unit,*)

end subroutine alprint

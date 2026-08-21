!----------------------------------------------------------------------
! altabinit
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original altabinit.f; only variable names, source form, and comment
! style were updated.
!
! PURPOSE
! To analyze the internal Allard-format tables provided by alfilein.f90,
! verifying them and creating needed additional tables and elements.
!
! Several steps are taken:
! 1. Find the minimum (TEFFLmin) and maximum (TEFFLmax) permissable values
!    of TEFFL. These are one row's width below the bottom and above the
!    top of the table. Because the first level of interpolation is in GL,
!    only a single minimum and maximum value of TEFFL are needed.
! 2. For each row in GL, find the index of the smallest element (iGLmin)
!    and the index of the largest element (iGLmax). There is one pair of
!    these for each TEFFL.
! 3. For each row in GL, find the the minimum (GLmin) and maximum (GLmax)
!    permissable value of GL. These are one column's width less than the
!    row minimum in GL and one column's width greater than the column
!    maximum. There is one pair of these for each TEFFL.
! 4. Validate the table. (a) Ensure that there are no invalid elements
!    inside the table, i.e., no invalid element between the row's iGLmin
!    and iGLmax. This check is made for every TEFFL. (b) Ensure that
!    every row has at least 4 valid entries and that there are at least
!    4 columns.
subroutine altabinit

      use luout_lib
      implicit none
      integer, parameter :: nta = 250
      integer, parameter :: nga = 25


! Shared: ALFILEIN, ALTABINIT and ALSURFP
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

! Shared: ALTABINIT and ALSURFP
      double precision :: allard_gl_row_min(nta), allard_gl_row_max(nta)
      integer :: allard_gl_index_min(nta), allard_gl_index_max(nta)
      double precision :: allard_teffl_min, allard_teffl_max, allard_gl_min, &
           allard_gl_max
      common /alatm02/ allard_gl_row_min, allard_gl_row_max, allard_gl_index_min, &
           allard_gl_index_max, allard_teffl_min, allard_teffl_max, &
           allard_gl_min, allard_gl_max

! Shared: ALFILEIN, ALSURFP and PARMIN
      double precision :: allard_target_feh, allard_target_alpha
      logical :: allard_use_tau100
      integer :: allard_table_unit
      common /alatm03/ allard_target_feh, allard_target_alpha, allard_use_tau100, &
           allard_table_unit

! MHP 8/25 Removed character file names from common block
! common/alatm04/: unused everywhere in this file family; placeholders
! preserving the shared storage layout.
      double precision :: alatm04_placeholder1, alatm04_placeholder2, &
           alatm04_placeholder3, alatm04_placeholder4
      common /alatm04/ alatm04_placeholder1, alatm04_placeholder2, &
           alatm04_placeholder3, alatm04_placeholder4

      double precision :: allard_al_teffl_min, allard_al_teffl_max
      common /alatm05/ allard_al_teffl_min, allard_al_teffl_max

      integer :: i, j, j1, j2
      logical :: table_is_bad

!     1. Find the minimum (TEFFLmin) and maximum (TEFFLmax) permissable values of TEFFL. These are
!        one row's width below the bottom and above the top of the table. Because the first level
!        of interpolation is in GL, only a single minimum and maximum value of TEFFL are needed.
      allard_al_teffl_min = allard_teffl_grid(1)
      allard_al_teffl_max = allard_teffl_grid(allard_num_teff)
      allard_teffl_min = allard_teffl_grid(1) - (allard_teffl_grid(2)-allard_teffl_grid(1))
      allard_teffl_max = allard_teffl_grid(allard_num_teff) + &
           (allard_teffl_grid(allard_num_teff)-allard_teffl_grid(allard_num_teff-1))

!     2. For each row in GL, find the index of the smallest element (iGLmin) and the index of the
!        largest element (iGLmax). There is one pair of these for each TEFFL.
      do i = 1, allard_num_teff
         do j = 1, allard_num_gl
            if (allard_log10_pressure(i,j) .gt. -998d0) then
               allard_gl_index_min(i) = j
               goto 100
            endif
         enddo
  100    continue
         do j = allard_num_gl, 1, -1
            if (allard_log10_pressure(i,j) .gt. -998d0) then
               allard_gl_index_max(i) = j
               goto 110
            endif
         enddo
  110    continue

      enddo

!     3. For each row in GL, find the the minimum (GLmin) and maximum (GLmax) permissable value of
!        GL. These are one column's width less than the row minimum in GL and one column's width
!        greater than the column maximum. There is one pair of these for each TEFFL.

      allard_gl_min = 999d0
      allard_gl_max = -999d0
      do i = 1, allard_num_teff
         j1 = allard_gl_index_min(i)
         j2 = allard_gl_index_max(i)
         allard_gl_row_min(i) = allard_gl_grid(j1) - 4d0*(allard_gl_grid(j1+1) - allard_gl_grid(j1))
         if (allard_gl_row_min(i) .lt. allard_gl_min) allard_gl_min = allard_gl_row_min(i)
         allard_gl_row_max(i) = allard_gl_grid(j2) + (allard_gl_grid(j2) - allard_gl_grid(j2-1))
         if (allard_gl_row_max(i) .gt. allard_gl_max) allard_gl_max = allard_gl_row_max(i)
      enddo

!     4. Validate the table. (a) Ensure that there are no invalid elements inside the table, i.e.,
!        no invalid element between the row's iGLmin and iGLmax. This check is made for every TEFFL.

      table_is_bad = .false.
      do i = 1, allard_num_teff
         j1 = allard_gl_index_min(i)
         j2 = allard_gl_index_max(i)
         do j = j1, j2
            if (allard_log10_pressure(i,j) .lt. -998d0) then
              table_is_bad = .true.
              write(short_file_unit,900) 'ALTABINIT: Bad input Allard Table: ', &
                 'TEFF, GL: ', 10d0**allard_teffl_grid(i), allard_gl_grid(j)
  900              format(2a,f5.0,f7.2)
            endif
         enddo
      enddo

!        (b) Ensure that every row has at least 4 valid entries and that there are at least 4 columns.
      if (allard_num_teff .lt. 4) then
         table_is_bad = .true.
         write(short_file_unit,910) 'ALTABINIT: Bad input Allard Table: ', &
            'Less than 4 rows: nTeff = ',allard_num_teff
  910         format(a,i3)
      endif
      do i = 1, allard_num_teff
         j1 = allard_gl_index_min(i)
         j2 = allard_gl_index_max(i)
         if ((j2 - j1 + 1) .lt. 4) then
            table_is_bad = .true.
            write(short_file_unit,920) 'ALTABINIT: Bad input Allard Table: ', &
              'Row with less that 4 elements: i,#,Teff,GLMin,GLmax: ', &
               i,j2-j1+1,10d0**allard_teffl_grid(i),allard_gl_row_min(i),allard_gl_row_max(i)
  920            format(2a,2i4,2x,3f7.2)
         endif
      enddo

      if (table_is_bad) goto 9999            ! If bad table, go to error exit

      return                  ! If good table, return

 9999      continue
       write(*,*)
       write(*,*)'******** ALTABINIT: Program Terminated ********'
       write(*,*)
       write(short_file_unit,*)
       write(short_file_unit,*)'******** ALTABINIT: Program Terminated ********'
       write(short_file_unit,*)
       call alprint
       stop

end subroutine altabinit

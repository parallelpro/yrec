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
subroutine altabinit(ierr)

      use atm_table_lib
      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: nta = 250
      integer, parameter :: nga = 25







      integer :: i, j, j1, j2
      logical :: table_is_bad

!     1. Find the minimum (TEFFLmin) and maximum (TEFFLmax) permissable values of TEFFL. These are
!        one row's width below the bottom and above the top of the table. Because the first level
!        of interpolation is in GL, only a single minimum and maximum value of TEFFL are needed.
      integer, intent(out) :: ierr

      ierr = 0

      atm_table%allard_al_teffl_min = atm_table%allard_teffl_grid(1)
      atm_table%allard_al_teffl_max = atm_table%allard_teffl_grid(atm_table%allard_num_teff)
      atm_table%allard_teffl_min = atm_table%allard_teffl_grid(1) - (atm_table%allard_teffl_grid(2)-atm_table%allard_teffl_grid(1))
      atm_table%allard_teffl_max = atm_table%allard_teffl_grid(atm_table%allard_num_teff) + &
           (atm_table%allard_teffl_grid(atm_table%allard_num_teff)-atm_table%allard_teffl_grid(atm_table%allard_num_teff-1))

!     2. For each row in GL, find the index of the smallest element (iGLmin) and the index of the
!        largest element (iGLmax). There is one pair of these for each TEFFL.
      do i = 1, atm_table%allard_num_teff
         do j = 1, atm_table%allard_num_gl
            if (atm_table%allard_log10_pressure(i,j) .gt. -998d0) then
               atm_table%allard_gl_index_min(i) = j
               exit
            endif
         enddo
         do j = atm_table%allard_num_gl, 1, -1
            if (atm_table%allard_log10_pressure(i,j) .gt. -998d0) then
               atm_table%allard_gl_index_max(i) = j
               exit
            endif
         enddo

      enddo

!     3. For each row in GL, find the the minimum (GLmin) and maximum (GLmax) permissable value of
!        GL. These are one column's width less than the row minimum in GL and one column's width
!        greater than the column maximum. There is one pair of these for each TEFFL.

      atm_table%allard_gl_min = 999d0
      atm_table%allard_gl_max = -999d0
      do i = 1, atm_table%allard_num_teff
         j1 = atm_table%allard_gl_index_min(i)
         j2 = atm_table%allard_gl_index_max(i)
         atm_table%allard_gl_row_min(i) = atm_table%allard_gl_grid(j1) - 4d0*(atm_table%allard_gl_grid(j1+1) - atm_table%allard_gl_grid(j1))
         if (atm_table%allard_gl_row_min(i) .lt. atm_table%allard_gl_min) atm_table%allard_gl_min = atm_table%allard_gl_row_min(i)
         atm_table%allard_gl_row_max(i) = atm_table%allard_gl_grid(j2) + (atm_table%allard_gl_grid(j2) - atm_table%allard_gl_grid(j2-1))
         if (atm_table%allard_gl_row_max(i) .gt. atm_table%allard_gl_max) atm_table%allard_gl_max = atm_table%allard_gl_row_max(i)
      enddo

!     4. Validate the table. (a) Ensure that there are no invalid elements inside the table, i.e.,
!        no invalid element between the row's iGLmin and iGLmax. This check is made for every TEFFL.

      table_is_bad = .false.
      do i = 1, atm_table%allard_num_teff
         j1 = atm_table%allard_gl_index_min(i)
         j2 = atm_table%allard_gl_index_max(i)
         do j = j1, j2
            if (atm_table%allard_log10_pressure(i,j) .lt. -998d0) then
              table_is_bad = .true.
              write(short_file_unit,900) 'ALTABINIT: Bad input Allard Table: ', &
                 'TEFF, GL: ', 10d0**atm_table%allard_teffl_grid(i), atm_table%allard_gl_grid(j)
  900              format(2a,f5.0,f7.2)
            endif
         enddo
      enddo

!        (b) Ensure that every row has at least 4 valid entries and that there are at least 4 columns.
      if (atm_table%allard_num_teff .lt. 4) then
         table_is_bad = .true.
         write(short_file_unit,910) 'ALTABINIT: Bad input Allard Table: ', &
            'Less than 4 rows: nTeff = ',atm_table%allard_num_teff
  910         format(a,i3)
      endif
      do i = 1, atm_table%allard_num_teff
         j1 = atm_table%allard_gl_index_min(i)
         j2 = atm_table%allard_gl_index_max(i)
         if ((j2 - j1 + 1) .lt. 4) then
            table_is_bad = .true.
            write(short_file_unit,920) 'ALTABINIT: Bad input Allard Table: ', &
              'Row with less that 4 elements: i,#,Teff,GLMin,GLmax: ', &
               i,j2-j1+1,10d0**atm_table%allard_teffl_grid(i),atm_table%allard_gl_row_min(i),atm_table%allard_gl_row_max(i)
  920            format(2a,2i4,2x,3f7.2)
         endif
      enddo

      if (.not. table_is_bad) then

      return                  ! If good table, return

      end if
       write(*,*)
       write(*,*)'******** ALTABINIT: Program Terminated ********'
       write(*,*)
       write(short_file_unit,*)
       write(short_file_unit,*)'******** ALTABINIT: Program Terminated ********'
       write(short_file_unit,*)
       call alprint
       ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the atm_lib
       ! facades stop when their caller passes no ierr.
       ierr = 1
       return

end subroutine altabinit

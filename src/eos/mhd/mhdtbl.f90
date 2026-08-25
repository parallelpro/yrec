!
!
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! MHDTBL
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
!----------------------------------------------------------------------
! mhdtbl
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mhdtbl.f; only variable names, source form, and comment
! style were updated.
!
!     FOR SET TABLES
! MHP 8/25 Removed un-used file names and file names from common blocks
! Opens the 8 MHD equation-of-state table files, calls mhdst to read
! and spline them, and closes the files again.
subroutine mhdtbl(zams_a_table_path, zams_b_table_path, zams_c_table_path, &
     centre1_table_path, centre2_table_path, centre3_table_path, &
     centre4_table_path, centre5_table_path, ierr)
      use star_info_lib, only: star
      use const_lib
      implicit none

!     COMMON/LUFNM/ FLAST, FFIRST, FRUN, FSTAND, FFERMI,
!    1    FDEBUG, FTRACK, FSHORT, FMILNE, FMODPT,
!    2    FSTOR, FPMOD, FPENV, FPATM, FDYN,
!    3    FLLDAT, FSNU, FSCOMP, FKUR,
!    4    FMHD1, FMHD2, FMHD3, FMHD4, FMHD5, FMHD6, FMHD7, FMHD8
      character(len=256), intent(in) :: zams_a_table_path, zams_b_table_path, &
           zams_c_table_path, centre1_table_path, centre2_table_path, &
           centre3_table_path, centre4_table_path, centre5_table_path
!     SET MHD EQUATION OF STATE
!     USES 3 ZAMS-TYPE TABLES AND OPTIONALLY 5 CENTER-TYPE TABLES.
!
      integer, intent(out) :: ierr

      ierr = 0

      open(unit=star%ctrl%unit_zams_a, file=zams_a_table_path, status='OLD', &
          form='UNFORMATTED')
      open(unit=star%ctrl%unit_zams_b, file=zams_b_table_path, status='OLD', &
          form='UNFORMATTED')
      open(unit=star%ctrl%unit_zams_c, file=zams_c_table_path, status='OLD', &
          form='UNFORMATTED')
      open(unit=star%ctrl%unit_centre1, file=centre1_table_path, status='OLD', &
          form='UNFORMATTED')
      open(unit=star%ctrl%unit_centre2, file=centre2_table_path, status='OLD', &
          form='UNFORMATTED')
      open(unit=star%ctrl%unit_centre3, file=centre3_table_path, status='OLD', &
          form='UNFORMATTED')
      open(unit=star%ctrl%unit_centre4, file=centre4_table_path, status='OLD', &
          form='UNFORMATTED')
      open(unit=star%ctrl%unit_centre5, file=centre5_table_path, status='OLD', &
          form='UNFORMATTED')

      call mhdst(star%ctrl%unit_zams_a, star%ctrl%unit_zams_b, star%ctrl%unit_zams_c, star%ctrl%unit_centre1, &
                 star%ctrl%unit_centre2, star%ctrl%unit_centre3, star%ctrl%unit_centre4, star%ctrl%unit_centre5, ierr)
      if (ierr /= 0) return
!     END MHD TABLE SETTING
      close(star%ctrl%unit_zams_a)
      close(star%ctrl%unit_zams_b)
      close(star%ctrl%unit_zams_c)
      close(star%ctrl%unit_centre1)
      close(star%ctrl%unit_centre2)
      close(star%ctrl%unit_centre3)
      close(star%ctrl%unit_centre4)
      close(star%ctrl%unit_centre5)
      return
end subroutine mhdtbl

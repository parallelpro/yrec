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
! Opens the 8 MHD equation-of-state table files (3 ZAMS-type tables
! and optionally 5 centre-type tables) on the unit numbers held in
! star%ctrl, calls mhdst to read them, and closes the files again.
subroutine mhdtbl(zams_a_table_path, zams_b_table_path, zams_c_table_path, &
     centre1_table_path, centre2_table_path, centre3_table_path, &
     centre4_table_path, centre5_table_path, ierr)
      use star_info_lib, only: star
      implicit none

      character(len=256), intent(in) :: zams_a_table_path, zams_b_table_path, &
           zams_c_table_path, centre1_table_path, centre2_table_path, &
           centre3_table_path, centre4_table_path, centre5_table_path
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

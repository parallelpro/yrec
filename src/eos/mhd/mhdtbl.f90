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
      implicit none

      character(len=256), intent(in) :: zams_a_table_path, zams_b_table_path, &
           zams_c_table_path, centre1_table_path, centre2_table_path, &
           centre3_table_path, centre4_table_path, centre5_table_path
      integer, intent(out) :: ierr
      integer :: unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, unit_centre2, unit_centre3, unit_centre4, unit_centre5
      ierr = 0

      open(newunit=unit_zams_a, file=zams_a_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_zams_b, file=zams_b_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_zams_c, file=zams_c_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre1, file=centre1_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre2, file=centre2_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre3, file=centre3_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre4, file=centre4_table_path, status='OLD', &
          form='UNFORMATTED')
      open(newunit=unit_centre5, file=centre5_table_path, status='OLD', &
          form='UNFORMATTED')

      call mhdst(unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, &
                 unit_centre2, unit_centre3, unit_centre4, unit_centre5, ierr)
      if (ierr /= 0) return
      close(unit_zams_a)
      close(unit_zams_b)
      close(unit_zams_c)
      close(unit_centre1)
      close(unit_centre2)
      close(unit_centre3)
      close(unit_centre4)
      close(unit_centre5)
      return
end subroutine mhdtbl

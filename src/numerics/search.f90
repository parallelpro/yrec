!----------------------------------------------------------------------
! search
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original search.f; only variable names, source form, and comment
! style were updated.
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
! SEARCH CONDUCTS A BINARY SEARCH FOR eval_point. SEARCH IS CALLED ONLY
! IF eval_point IS BETWEEN table_x(1) AND table_x(num_table_points).
!
! ON INPUT--
!
!   table_x CONTAINS THE ABSCISSAS OF THE DATA POINTS OF INTERPOLATION.
!
!   num_table_points IS THE DIMENSION OF table_x
!
!   eval_point IS THE VALUE WHOSE RELATIVE POSITION IN table_x IS
!   LOCATED BY SEARCH.
!
!
! ON OUTPUT--
!
!   found_flag IS SET EQUAL TO 1 IF eval_point IS FOUND IN table_x AND
!   IS SET EQUAL TO 0 OTHERWISE.
!
!   table_idx IS THE INDEX OF THE LARGEST VALUE IN table_x FOR WHICH
!   table_x(I) .LT. eval_point.
!
! AND
!
!   SEARCH DOES NOT ALTER table_x,num_table_points,eval_point.
!
! NOTE: unlike most YREC files, the original search.f has no blanket
! IMPLICIT REAL*8(A-H,O-Z) statement -- it relies on default Fortran
! implicit typing (I-N integer, else real) except where XTAB/S and
! FND/FIRST are explicitly declared. Types below are chosen to match
! that original default typing exactly (NUM/LCN/MIDDLE/LAST all fall
! in the I-N default-integer range).
!
subroutine search(table_x, num_table_points, eval_point, table_idx, &
     found_flag)
      use star_info_lib, only: json

      implicit none

      double precision, intent(in) :: table_x(json), eval_point
      integer, intent(in) :: num_table_points
      integer, intent(out) :: table_idx, found_flag
      integer :: first_idx, last_idx, middle_idx

      first_idx=1
      last_idx=num_table_points
      found_flag=0
!
! (Restructured 2026 from the original goto binary search at labels
! 10-60; comparisons are unchanged.)
      if (table_x(1) .eq. eval_point) then
         table_idx=1
         found_flag=1
         return
      end if
      if (table_x(num_table_points) .eq. eval_point) then
         table_idx=num_table_points
         found_flag=1
         return
      end if
!
      do
!
! IF (LAST-FIRST) .EQ. 1, S IS NOT IN XTAB.  SET POSITION EQUAL TO
! FIRST.
      if ((last_idx-first_idx) .eq. 1) then
         table_idx=first_idx
         return
      end if
!
      middle_idx=(first_idx+last_idx)/2
!
! CHECK IF S .EQ. XTAB(MIDDLE). IF NOT, CONTINUE THE SEARCH IN THE
! APPROPRIATE HALF OF THE VECTOR XTAB.
      if (table_x(middle_idx) .lt. eval_point) then
         first_idx=middle_idx
      else if (table_x(middle_idx) .eq. eval_point) then
         table_idx=middle_idx
         found_flag =1
         return
      else
         last_idx=middle_idx
      end if
      end do
end subroutine search

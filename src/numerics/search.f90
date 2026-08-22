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

      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: table_x(json), eval_point
      integer, intent(in) :: num_table_points
      integer, intent(out) :: table_idx, found_flag
      integer :: first_idx, last_idx, middle_idx

      first_idx=1
      last_idx=num_table_points
      found_flag=0
!
      if (table_x(1) .ne. eval_point) go to 10
      table_idx=1
      found_flag=1
      return
!
  10  if (table_x(num_table_points) .ne. eval_point) go to 20
      table_idx=num_table_points
      found_flag=1
      return
!
! IF (LAST-FIRST) .EQ. 1, S IS NOT IN XTAB.  SET POSITION EQUAL TO
! FIRST.
  20  if ((last_idx-first_idx) .eq. 1) go to 30
!
      middle_idx=(first_idx+last_idx)/2
!
! CHECK IF S .EQ. XTAB(MIDDLE). IF NOT, CONTINUE THE SEARCH IN THE
! APPROPRIATE HALF OF THE VECTOR XTAB.
! KC 2025-05-30 fixed "Arithmetic IF statement"
!       IF (XTAB(MIDDLE) - S) 40,50,60
      if (table_x(middle_idx) .lt. eval_point) then
         goto 40
      else if (table_x(middle_idx) .eq. eval_point) then
         goto 50
      else
         goto 60
      end if
!
  30  table_idx=first_idx
      return
  40  first_idx=middle_idx
      go to 20
  50  table_idx=middle_idx
      found_flag =1
      return
  60  last_idx=middle_idx
      go to 20
end subroutine search

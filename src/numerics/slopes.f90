!----------------------------------------------------------------------
! slopes
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original slopes.f; only variable names, source form, and comment
! style were updated.
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C.STATE UNIVERSITY
!
! SLOPES CALCULATES THE DERIVATIVE AT EACH OF THE DATA POINTS. THE
! SLOPES PROVIDED WILL INSURE THAT AN OSCULATORY QUADRATIC SPLINE WILL
! HAVE ONE ADDITIONAL KNOT BETWEEN TWO ADJACENT POINTS OF INTERPOLATION.
! CONVEXITY AND MONOTONICITY ARE PRESERVED WHEREVER THESE CONDITIONS
! ARE COMPATIBLE WITH THE DATA.
!
! ON INPUT--
!
!   table_x CONTAINS THE ABSCISSAS OF THE DATA POINTS.
!
!   table_y CONTAINS THE ORDINATES OF THE DATA POINTS.
!
!   num_points IS THE NUMBER OF DATA POINTS (DIMENSION OF table_x,
!   table_y).
!
!
! ON OUTPUT--
!
!   first_derivs CONTAINS THE VALUE OF THE FIRST DERIVATIVE AT EACH
!   DATA POINT.
!
! AND
!
!   SLOPES DOES NOT ALTER table_x,table_y,num_points.
!
! NOTE: like search.f, the original slopes.f has no blanket IMPLICIT
! REAL*8(A-H,O-Z) statement -- it relies on default Fortran implicit
! typing (I-N integer, else real) except where table_x/table_y/
! first_derivs and the M1/M2/... scratch variables are explicitly
! declared REAL*8. Types below match that original default typing.
!
!----------------------------------------------------------------------
!
subroutine slopes(table_x, table_y, first_derivs, num_points)

      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: table_x(json), table_y(json)
      double precision, intent(out) :: first_derivs(json)
      integer, intent(in) :: num_points
      integer :: num_points_m1, prev_idx, idx, next_idx
      double precision :: slope1, slope2, x_bar, x_hat, y_diff1, y_diff2, &
           y_x_mid, x_mid, slope1_saved, slope2_saved

      num_points_m1= num_points - 1
      prev_idx=1
      idx=2
      next_idx=3
!
! CALCULATE THE SLOPES OF THE TWO LINES JOINING THE FIRST THREE DATA
! POINTS
      y_diff1=table_y(2) - table_y(1)
      y_diff2=table_y(3) - table_y(2)
      slope1=y_diff1/(table_x(2) - table_x(1))
      slope1_saved=slope1
      slope2=y_diff2/(table_x(3)-table_x(2))
      slope2_saved=slope2
!
! IF ONE OF THE PRECEDING SLOPES IS ZERO OR IF THEY HAVE OPPOSITE SIGN,
! ASSIGN THE VALUE ZERO TO THE DERIVATIVE AT THE MIDDLE POINT.
  10  if(slope1.eq.0.d0 .or. slope2.eq.0.d0 .or. (slope1*slope2).le.0.d0) go to 20
      if (abs(slope1) .gt. abs(slope2)) go to 30
      go to 40
  20  first_derivs(idx)= 0.d0
      go to 50
!
! CALCULATE THE SLOPE BY EXTENDING THE LINE WITH SLOPE M1.
  30  x_bar=(y_diff2/slope1) + table_x(idx)
      x_hat= (x_bar + table_x(next_idx))/2.d0
      first_derivs(idx)=y_diff2/(x_hat - table_x(idx))
      go to 50
!
! CALCULATE THE SLOPE BY EXTENDING THE LINE WITH SLOPE M2.
  40  x_bar=(-y_diff1/slope2) + table_x(idx)
      x_hat=(table_x(prev_idx) + x_bar)/2.d0
      first_derivs(idx)=y_diff1/(table_x(idx) - x_hat)
!
! INCREMENT COUNTERS
  50  prev_idx=idx
      idx=next_idx
      next_idx=next_idx+1
      if (idx .gt. num_points_m1) go to 60
!
! CALCULATE THE SLOPES OF THE TWO LINES JOINING THREE CONSECUTIVE DATA
! POINTS.
      y_diff1=table_y(idx) - table_y(prev_idx)
      y_diff2=table_y(next_idx) - table_y(idx)
      slope1=y_diff1/(table_x(idx) - table_x(prev_idx))
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION
!       M2=YDIF2/(XTAB(I1) - XTAB(I))
      call safedivide(y_diff2, (table_x(next_idx) - table_x(idx)), slope2)
      go to 10
!
! CALCULATE THE SLOPE AT THE LAST POINT, XTAB(NUM).
  60  if ((slope1*slope2) .lt. 0.d0) go to 80
      x_mid= (table_x(num_points_m1)+table_x(num_points))/2.d0
      y_x_mid=first_derivs(num_points_m1)*(x_mid - table_x(num_points_m1)) + &
           table_y(num_points_m1)
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION
!       MTAB(NUM)=(YTAB(NUM)-YXMID)/(XTAB(NUM)-XMID)
      call safedivide((table_y(num_points)-y_x_mid), &
           (table_x(num_points)-x_mid), first_derivs(num_points))
      if ((first_derivs(num_points)*slope2) .lt. 0.d0) go to 70
      go to 90
  70  first_derivs(num_points)=0.d0
      go to 90
  80  first_derivs(num_points)=2.d0*slope2
!
! CALCULATE THE SLOPE AT THE FIRST POINT, XTAB(1).
  90  if((slope1_saved*slope2_saved) .lt. 0.d0) go to 110
      x_mid=(table_x(1) + table_x(2))/2.d0
      y_x_mid=first_derivs(2)*(x_mid - table_x(2)) + table_y(2)
      first_derivs(1)=(y_x_mid - table_y(1))/(x_mid - table_x(1))
      if ((first_derivs(1) * slope1_saved) .lt. 0.d0) go to 100
      return
  100 first_derivs(1)=0.d0
      return
  110 first_derivs(1)=2.d0*slope1_saved
      return
!
end subroutine slopes

!----------------------------------------------------------------------
! KC 2025-05-31 SAFEDIVIDE
!----------------------------------------------------------------------
subroutine safedivide(numerator, denominator, quotient)
      implicit none
      double precision, intent(in) :: numerator, denominator
      double precision, intent(out) :: quotient

      quotient = 0.d0
      if (numerator .ne. 0.d0 .and. denominator .ne. 0.d0) then
         quotient = numerator / denominator
      end if

      return
end subroutine safedivide

!----------------------------------------------------------------------
! quadratic_spline
!----------------------------------------------------------------------
! The shape-preserving osculatory quadratic spline of ACM Algorithm 574
! (McAllister and Roulier, ACM TOMS, September 1981): convexity and
! monotonicity are kept wherever the data allow (YREC's numerics_lib,
! split 2026-10). osplin is the entry point; slopes chooses the slopes
! at the points, and meval evaluates the spline through choose (the knot
! case of an interval), cases (its knots) and spline (the value), with
! search for the first interval.
! With the slopes that slopes computes, every slope is below twice both
! adjacent secants and never opposes their sign, so choose only returns
! cases 1 and 2.
module quadratic_spline
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: osplin

contains

!----------------------------------------------------------------------
! osplin
!----------------------------------------------------------------------
! The spline through xtab/ytab(1:n) at the k nondecreasing points
! xval, into yval; outside the table it extrapolates the end
! intervals' quadratics.
subroutine osplin(xval, yval, xtab, ytab, n, k)
      integer, intent(in) :: n, k
      real(dp), intent(in) :: xval(k)
      real(dp), intent(out) :: yval(k)
      real(dp), intent(in) :: xtab(n), ytab(n)

      real(dp) :: first_derivs(n)
! error tolerance used in subroutine 'choose' (via meval)
      real(dp), parameter :: eps = 1.d-04
      integer :: err

! calculate the slopes at each data point.
      call slopes(xtab, ytab, first_derivs, n)

! call meval to evaluate the spline at the run of points xval.
      call meval(xval, yval, xtab, ytab, first_derivs, n, k, eps, err)

      return
end subroutine osplin

!----------------------------------------------------------------------
! slopes
!----------------------------------------------------------------------
! Shape-preserving quadratic splines, D.F. McAllister and J.A. Roulier,
! coded by S.L. Dodd and M. Roulier, N.C. State University.
! The slope at each point, such that the osculatory quadratic spline
! needs one knot per interval and keeps the data's convexity and
! monotonicity. Interior points get a weighted mean of the two
! adjacent secants (zero where they differ in sign or one is zero); the
! end points extrapolate from their neighbour. (KC 2025-05-31:
! safedivide guards two divisions by zero.)
subroutine slopes(table_x, table_y, first_derivs, num_points)
      integer, intent(in) :: num_points
      real(dp), intent(in) :: table_x(num_points), table_y(num_points)
      real(dp), intent(out) :: first_derivs(num_points)
      integer :: num_points_m1, prev_idx, idx, next_idx
      real(dp) :: slope1, slope2, x_bar, x_hat, y_diff1, y_diff2, &
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
! MAIN LOOP OVER THE INTERIOR POINTS. (Restructured 2026 from the
! original goto flow at labels 10-50; arithmetic is unchanged.)
!
! IF ONE OF THE PRECEDING SLOPES IS ZERO OR IF THEY HAVE OPPOSITE SIGN,
! ASSIGN THE VALUE ZERO TO THE DERIVATIVE AT THE MIDDLE POINT.
      do
      if (slope1.eq.0.d0 .or. slope2.eq.0.d0 .or. (slope1*slope2).le.0.d0) then
         first_derivs(idx)= 0.d0
      else if (abs(slope1) .gt. abs(slope2)) then
!
! CALCULATE THE SLOPE BY EXTENDING THE LINE WITH SLOPE M1.
         x_bar=(y_diff2/slope1) + table_x(idx)
         x_hat= (x_bar + table_x(next_idx))/2.d0
         first_derivs(idx)=y_diff2/(x_hat - table_x(idx))
      else
!
! CALCULATE THE SLOPE BY EXTENDING THE LINE WITH SLOPE M2.
         x_bar=(-y_diff1/slope2) + table_x(idx)
         x_hat=(table_x(prev_idx) + x_bar)/2.d0
         first_derivs(idx)=y_diff1/(table_x(idx) - x_hat)
      end if
!
! INCREMENT COUNTERS
      prev_idx=idx
      idx=next_idx
      next_idx=next_idx+1
      if (idx .gt. num_points_m1) exit
!
! CALCULATE THE SLOPES OF THE TWO LINES JOINING THREE CONSECUTIVE DATA
! POINTS.
      y_diff1=table_y(idx) - table_y(prev_idx)
      y_diff2=table_y(next_idx) - table_y(idx)
      slope1=y_diff1/(table_x(idx) - table_x(prev_idx))
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION
!       M2=YDIF2/(XTAB(I1) - XTAB(I))
      call safedivide(y_diff2, (table_x(next_idx) - table_x(idx)), slope2)
      end do
!
! CALCULATE THE SLOPE AT THE LAST POINT, XTAB(NUM).
      if ((slope1*slope2) .lt. 0.d0) then
         first_derivs(num_points)=2.d0*slope2
      else
         x_mid= (table_x(num_points_m1)+table_x(num_points))/2.d0
         y_x_mid=first_derivs(num_points_m1)*(x_mid - table_x(num_points_m1)) + &
              table_y(num_points_m1)
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION
!       MTAB(NUM)=(YTAB(NUM)-YXMID)/(XTAB(NUM)-XMID)
         call safedivide((table_y(num_points)-y_x_mid), &
              (table_x(num_points)-x_mid), first_derivs(num_points))
         if ((first_derivs(num_points)*slope2) .lt. 0.d0) then
            first_derivs(num_points)=0.d0
         end if
      end if
!
! CALCULATE THE SLOPE AT THE FIRST POINT, XTAB(1).
      if ((slope1_saved*slope2_saved) .lt. 0.d0) then
         first_derivs(1)=2.d0*slope1_saved
         return
      end if
      x_mid=(table_x(1) + table_x(2))/2.d0
      y_x_mid=first_derivs(2)*(x_mid - table_x(2)) + table_y(2)
      first_derivs(1)=(y_x_mid - table_y(1))/(x_mid - table_x(1))
      if ((first_derivs(1) * slope1_saved) .lt. 0.d0) then
         first_derivs(1)=0.d0
      end if
      return
!
end subroutine slopes

!----------------------------------------------------------------------
! safedivide
!----------------------------------------------------------------------
! numerator/denominator, or zero when either is zero (KC 2025-05-31).
subroutine safedivide(numerator, denominator, quotient)
      real(dp), intent(in) :: numerator, denominator
      real(dp), intent(out) :: quotient

      quotient = 0.d0
      if (numerator .ne. 0.d0 .and. denominator .ne. 0.d0) then
         quotient = numerator / denominator
      end if

      return
end subroutine safedivide

!----------------------------------------------------------------------
! cases
!----------------------------------------------------------------------
! Shape-preserving quadratic splines, D.F. McAllister and J.A. Roulier,
! coded by S.L. Dodd and M. Roulier, N.C. State University.
! The knots and spline parameters on (x_left, x_right) for knot case
! spline_case: (v, z, w) for cases 1-3; case 4 adds the second knot y
! and e.
subroutine cases(x_left, y_left, slope_left, slope_right, x_right, &
     y_right, knot_e_x, knot_e_y, knot_v_x, knot_v_y, knot_w_x, &
     knot_w_y, knot_z_x, knot_z_y, knot_y_x, knot_y_y, spline_case)
      real(dp), intent(in) :: x_left, y_left, slope_left, slope_right, &
           x_right, y_right
      integer, intent(in) :: spline_case
      real(dp), intent(out) :: knot_e_x, knot_e_y, knot_v_x, knot_v_y, &
           knot_w_x, knot_w_y, knot_z_x, knot_z_y, knot_y_x, knot_y_y
      real(dp) :: mbar1, mbar2, mbar3, c1, d1, h1, j1, k1, ztwo

      if (spline_case .ne. 3 .and. spline_case .ne. 4) then
      if (spline_case .ne. 2) then

! CALCULATE THE PARAMETERS FOR CASE 1.
      knot_z_x=(y_left-y_right+slope_right*x_right-slope_left*x_left)/ &
           (slope_right-slope_left)
      ztwo=y_left+slope_left*(knot_z_x-x_left)
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=(y_left+ztwo)/2.d0
      knot_w_x=(knot_z_x+x_right)/2.d0
      knot_w_y=(ztwo+y_right)/2.d0
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
           (knot_z_x-knot_v_x)
      return

! CALCULATE THE PARAMETERS FOR CASE 2.
      end if
      knot_z_x=(x_left+x_right)/2.d0
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(knot_z_x+x_right)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=(knot_v_y+knot_w_y)/2.d0
      return

! CALCULATE THE PARAMETERS USED IN BOTH CASES 3 AND 4.
      end if
      c1=x_left+(y_right-y_left)/slope_left
      d1=x_right+(y_left-y_right)/slope_right
      h1=2.d0*c1-x_left
      j1=2.d0*d1-x_right
      mbar1=(y_right-y_left)/(h1-x_left)
      mbar2=(y_left-y_right)/(j1-x_right)

      if (spline_case .ne. 4) then

! CALCULATE THE PARAMETERS FOR CASE 3.
      k1=(y_left-y_right+x_right*mbar2-x_left*mbar1)/(mbar2-mbar1)
      if (abs(slope_left) .le. abs(slope_right)) then
      knot_z_x=(k1+x_right)/2.d0
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
      (knot_z_x-knot_v_x)
      return
      end if
      knot_z_x=(k1+x_left)/2.d0
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
           (knot_z_x-knot_v_x)
      return

! CALCULATE THE PARAMETERS FOR CASE 4.
      end if
      knot_y_x=(x_left+c1)/2.d0
      knot_v_x=(x_left+knot_y_x)/2.d0
      knot_v_y=slope_left*(knot_v_x-x_left) + y_left
      knot_z_x=(d1+x_right)/2.d0
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=slope_right*(knot_w_x-x_right) + y_right
      mbar3=(knot_w_y-knot_v_y)/(knot_w_x-knot_v_x)
      knot_y_y=mbar3*(knot_y_x-knot_v_x) + knot_v_y
      knot_z_y=mbar3*(knot_z_x-knot_v_x) + knot_v_y
      knot_e_x=(knot_y_x+knot_z_x)/2.d0
      knot_e_y=mbar3*(knot_e_x-knot_v_x) + knot_v_y
      return

end subroutine cases

!----------------------------------------------------------------------
! choose
!----------------------------------------------------------------------
! Shape-preserving quadratic splines, D.F. McAllister and J.A. Roulier,
! coded by S.L. Dodd and M. Roulier, N.C. State University.
! The knot case of the interval from the end slopes and the secant
! spq: 1 when the end tangents meet inside the interval, 2 when the
! midpoint is the knot (MHP: also when an end slope is zero, avoiding a
! division by zero in cases), 3 or 4 when an end slope exceeds twice
! the secant (one or two knots). eps_tol separates slopes close to the
! secant or twice it.
subroutine choose(x_left, y_left, slope_left, slope_right, x_right, &
     y_right, eps_tol, spline_case)
      real(dp), intent(in) :: x_left, y_left, slope_left, slope_right, &
           x_right, y_right, eps_tol
      integer, intent(out) :: spline_case
      real(dp) :: mref, mref1, mref2, spq, prod, prod1, prod2

! CALCULATE THE SLOPE spq OF THE LINE JOINING (x_left,y_left),(x_right,y_right).
      spq=(y_right-y_left)/(x_right-x_left)

! CHECK WHETHER OR NOT spq IS 0.
! ******MODIFICATION BY MARC PINSONNEAULT TO AVOID DIVISION BY ZERO
! ******IN SR CASES
      if (slope_left.eq.0.0d0 .or. slope_right.eq.0.0d0) then
         spline_case=2
         return
      end if
! ******
      if (spq .eq. 0.d0) then
         if ((slope_left*slope_right) .ge. 0.d0) then
            spline_case=2
         else
            spline_case=1
         end if
         return
      end if

      prod1=spq*slope_left
      prod2=spq*slope_right

! FIND THE ABSOLUTE VALUES OF THE SLOPES spq,slope_left,AND slope_right.
      mref=abs(spq)
      mref1=abs(slope_left)
      mref2=abs(slope_right)

! THE SIGN OF AT LEAST ONE OF THE SLOPES slope_left,slope_right DOES NOT
! AGREE WITH THE SIGN OF THE SLOPE spq.
      if ((prod1 .lt. 0.d0) .or. (prod2 .lt. 0.d0)) then
         if ((prod1 .lt. 0.d0) .and. (prod2 .lt. 0.d0)) then
            spline_case=2
            return
         end if
         if (prod1 .ge. 0.d0) then
            if (mref1 .gt. ((1.d0+eps_tol)*mref)) then
               spline_case=1
            else
               spline_case=2
            end if
            return
         end if
         if (mref2 .gt. ((1.d0+eps_tol)*mref)) then
            spline_case=1
         else
            spline_case=2
         end if
         return
      end if

! IF THE RELATIVE DEVIATION OF slope_left OR slope_right FROM spq IS LESS THAN
! eps_tol, THEN CHOOSE CASE 2 OR CASE 3.
      if (abs(spq-slope_left).gt.eps_tol*mref .and. abs(spq-slope_right).gt.eps_tol*mref) then
         prod=(mref-mref1)*(mref-mref2)
         if (prod .lt. 0.d0) then

! L1, THE LINE THROUGH (x_left,y_left) WITH SLOPE slope_left, AND L2, THE LINE
! THROUGH (x_right,y_right) WITH SLOPE slope_right, INTERSECT AT A POINT WHOSE
! ABSCISSA IS BETWEEN x_left AND x_right. THE ABSCISSA BECOMES A KNOT OF THE
! SPLINE.
            spline_case=1
            return
         end if
      end if

! IN CASES 3 AND 4, SIGN(slope_left)=SIGN(slope_right)=SIGN(spq).
! CHOOSE CASE 4 IF THE OTHER SLOPE IS GREATER THAN (2.-eps_tol)*mref
! (NEITHER L1 NOR L2 CROSSES THE MIDLINE: TWO KNOTS); OTHERWISE CASE 3
! (EITHER L1 OR L2 CROSSES THE MIDLINE, BUT NOT BOTH).
      if (mref1 .gt. (2.d0*mref)) then
         if (mref2 .gt. (2.d0-eps_tol)*mref) then
            spline_case=4
         else
            spline_case=3
         end if
         return
      end if
      if (mref2 .gt. (2.d0*mref)) then
         if (mref1 .gt. (2.d0-eps_tol)*mref) then
            spline_case=4
         else
            spline_case=3
         end if
         return
      end if

! BOTH L1 AND L2 CROSS THE LINE THROUGH (x_left+x_right/2.,y_left) AND
! (x_left+x_right/2.,y_right), WHICH IS THE MIDLINE OF THE RECTANGLE FORMED
! BY (x_left,y_left),(x_right,y_left),(x_right,y_right), AND (x_left,y_right),
! OR BOTH slope_left AND slope_right HAVE SIGNS DIFFERENT THAN THE SIGN OF
! spq, OR ONE OF slope_left AND slope_right HAS OPPOSITE SIGN FROM spq AND L1
! AND L2 INTERSECT TO THE LEFT OF x_left OR TO THE RIGHT OF x_right. THE
! POINT (x_left+x_right)/2. IS A KNOT OF THE SPLINE.
! (Restructured 2026 from the original goto decision tree; two
! unreachable branches at old labels 60/110 were dropped.)
      spline_case=2
      return

end subroutine choose

!----------------------------------------------------------------------
! search
!----------------------------------------------------------------------
! Shape-preserving quadratic splines, D.F. McAllister and J.A. Roulier,
! coded by S.L. Dodd and M. Roulier, N.C. State University.
! Binary search for eval_point in table_x(1:num_table_points), which
! it lies within: table_idx is the last point below it, found_flag 1
! if it equals a point.
subroutine search(table_x, num_table_points, eval_point, table_idx, &
     found_flag)
      integer, intent(in) :: num_table_points
      real(dp), intent(in) :: table_x(num_table_points), eval_point
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

!----------------------------------------------------------------------
! spline
!----------------------------------------------------------------------
! Shape-preserving quadratic splines, D.F. McAllister and J.A. Roulier,
! coded by S.L. Dodd and M. Roulier, N.C. State University.
! The spline's value at eval_point in (x_left, x_right) from cases'
! parameters; one knot (z) for cases 1-3, two (y, z) for case 4. MHP
! 6/87: linear interpolation where a knot coincides with an end.
function spline(eval_point, z1, z2, x_left, y_left, x_right, y_right, &
     y1, y2, e2, w2, v2, spline_case)
      real(dp) :: spline
      real(dp), intent(in) :: eval_point, z1, z2, x_left, y_left, &
           x_right, y_right, y1, y2, e2, w2, v2
      integer, intent(in) :: spline_case
      real(dp) :: linear_interp_frac

! (Restructured 2026 from the original arithmetic-IF goto fans at
! labels 10-100; arithmetic is unchanged.)
      if (spline_case .ne. 4) then
!
! CASES 1,2, OR 3.
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE KNOT.
      if (z1 .lt. eval_point) then
         if(x_right.ne.z1)then
         spline=(z2*(x_right-eval_point)**2+w2*2.d0*(eval_point-z1)*(x_right-eval_point) &
                 +y_right*(eval_point-z1)**2)/(x_right-z1)**2
         else
          linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
          spline = y_left + linear_interp_frac*(y_right - y_left)
         end if
      else if (z1 .eq. eval_point) then
         spline=z2
      else
         if(z1.ne.x_left)then
         spline=(y_left*(z1-eval_point)**2+v2*2.d0*(eval_point-x_left)*(z1-eval_point)+ &
                 z2*(eval_point-x_left)**2)/(z1-x_left)**2
         else
          linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
          spline = y_left + linear_interp_frac*(y_right - y_left)
         end if
      end if
      return
      end if
!
! CASE 4.
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE FIRST KNOT.
      if (y1 .lt. eval_point) then
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE SECOND KNOT.
         if (z1 .lt. eval_point) then
            if(x_right.ne.z1)then
            spline=(z2*(x_right-eval_point)**2+w2*2.d0*(eval_point-z1)*(x_right-eval_point) &
                    +y_right*(eval_point-z1)**2)/(x_right-z1)**2
            else
             linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
             spline = y_left + linear_interp_frac*(y_right - y_left)
            end if
         else if (z1 .eq. eval_point) then
            spline=z2
         else
            if(z1.ne.y1)then
            spline=(y2*(z1-eval_point)**2+e2*2.d0*(eval_point-y1)*(z1-eval_point)+z2*(eval_point &
                    -y1)**2)/(z1-y1)**2
            else
             linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
             spline = y_left + linear_interp_frac*(y_right - y_left)
            end if
         end if
      else if (y1 .eq. eval_point) then
         spline=y2
      else
         if(y1.ne.x_left)then
         spline=(y_left*(y1-eval_point)**2+v2*2.d0*(eval_point-x_left)*(y1-eval_point)+ &
                 y2*(eval_point-x_left)**2)/(y1-x_left)**2
         else
          linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
          spline = y_left + linear_interp_frac*(y_right - y_left)
         end if
      end if
      return
end function spline

!----------------------------------------------------------------------
! meval
!----------------------------------------------------------------------
! Shape-preserving quadratic splines, D.F. McAllister and J.A. Roulier,
! coded by S.L. Dodd and M. Roulier, N.C. State University.
! Evaluates the spline through table_x/table_y with slopes table_slope
! at the nondecreasing points eval_x into eval_y. err_code: 0 normal;
! 1 some points lie outside the table (the end quadratics are
! extrapolated); 2 eval_x decreases somewhere (nothing evaluated).
! llp 8/19/08: the tail uses num_table_points - 1, never 0.
subroutine meval(eval_x, eval_y, table_x, table_y, table_slope, &
     num_table_points, num_eval_points, eps_tol, err_code)
      integer, intent(in) :: num_table_points, num_eval_points
      real(dp), intent(in) :: eval_x(num_eval_points), &
           table_x(num_table_points), table_y(num_table_points), &
           table_slope(num_table_points), eps_tol
      real(dp), intent(out) :: eval_y(num_eval_points)
      integer, intent(out) :: err_code

      real(dp) :: spline_v1, spline_v2, spline_w1, spline_w2, &
           spline_z1, spline_z2, spline_y1, spline_y2, spline_e1, spline_e2
      integer :: start_idx, start_idx1, end_idx, end_idx1, found_flag
      integer :: i, ind, loop_bound, table_idx, table_idx1, spline_case, &
           num_minus1
      logical :: recompute_tail_params
      start_idx = 1
      end_idx = num_eval_points
      err_code = 0
      if (num_eval_points .ne. 1) then

! Determine if eval_x is nondecreasing.
      loop_bound = num_eval_points - 1
      do i = 1, loop_bound
        if (eval_x(i+1) .ge. eval_x(i)) cycle
        err_code = 2
        return
      end do

! If eval_x(i) .lt. table_x(1), then eval_x(i)=table_y(1).
! If eval_x(i) .gt. table_x(num_table_points), then
! eval_x(i)=table_y(num_table_points).
!
! Determine if any of the points in eval_x are less than the abscissa
! of the first data point.
      end if
      do i = 1, num_eval_points
        if (eval_x(i) .ge. table_x(1)) exit
        start_idx = i + 1
      end do

      loop_bound = num_eval_points + 1

! Determine if any of the points in eval_x are greater than the
! abscissa of the last data point.
      do i = 1, num_eval_points
        ind = loop_bound - i
        if (eval_x(ind) .le. table_x(num_table_points)) exit
        end_idx = ind - 1
      end do

! Calculate the images of points of evaluation whose abscissas are
! less than the abscissa of the first data point.
      if (start_idx .ne. 1) then
! Set the error parameter to indicate that extrapolation has occurred.
      err_code = 1
      call choose(table_x(1), table_y(1), table_slope(1), table_slope(2), &
           table_x(2), table_y(2), eps_tol, spline_case)
      call cases(table_x(1), table_y(1), table_slope(1), table_slope(2), &
           table_x(2), table_y(2), spline_e1, spline_e2, spline_v1, &
           spline_v2, spline_w1, spline_w2, spline_z1, spline_z2, &
           spline_y1, spline_y2, spline_case)
      start_idx1 = start_idx - 1
      do i = 1, start_idx1
       eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, table_x(1), &
            table_y(1), table_x(2), table_y(2), spline_y1, spline_y2, &
            spline_e2, spline_w2, spline_v2, spline_case)
      end do
      if (num_eval_points .eq. 1) then
         return
      end if
      end if

! search locates the interval in which the first in-range point of
! evaluation lies. The single-out-of-range-high case skips straight
! to the high-extrapolation tail.
      recompute_tail_params = .true.
      if (num_eval_points .ne. 1 .or. end_idx .eq. num_eval_points) then
      call search(table_x, num_table_points, eval_x(start_idx), &
           table_idx, found_flag)

      table_idx1 = table_idx + 1

! If the first in-range point of evaluation is equal to one of the
! data points, assign the appropriate value from table_y. Continue
! until a point of evaluation is found which is not equal to a data
! point.
      if (found_flag .ne. 0) then
      do
         eval_y(start_idx) = table_y(table_idx)
      start_idx1 = start_idx
      start_idx = start_idx + 1
      if (start_idx .gt. num_eval_points) then
         return
      end if
      if (eval_x(start_idx1) .ne. eval_x(start_idx)) exit
      end do

! Advance the table pointers until the next point of evaluation lies
! inside the current table interval, assigning images to any points
! of evaluation that coincide with data points along the way.
! (Restructured 2026 from the original arithmetic-IF web at labels
! 100/110/120; comparisons and order of assignment are unchanged.)
      do
         if (eval_x(start_idx) .lt. table_x(table_idx1)) exit
         if (eval_x(start_idx) .eq. table_x(table_idx1)) then
            do
               eval_y(start_idx) = table_y(table_idx1)
               start_idx1 = start_idx
               start_idx = start_idx + 1
               if (start_idx .gt. num_eval_points) then
                  return
               end if
               if (eval_x(start_idx) .ne. eval_x(start_idx1)) exit
            end do
         end if
         table_idx = table_idx1
         table_idx1 = table_idx1 + 1
      end do

! Calculate the images of all the points which lie within range of
! the data.
      end if
      if (table_idx .ne. 1 .or. err_code .ne. 1) then
      call choose(table_x(table_idx), table_y(table_idx), &
           table_slope(table_idx), table_slope(table_idx1), &
           table_x(table_idx1), table_y(table_idx1), eps_tol, spline_case)
      call cases(table_x(table_idx), table_y(table_idx), &
           table_slope(table_idx), table_slope(table_idx1), &
           table_x(table_idx1), table_y(table_idx1), spline_e1, spline_e2, &
           spline_v1, spline_v2, spline_w1, spline_w2, spline_z1, &
           spline_z2, spline_y1, spline_y2, spline_case)
      end if

      do i = start_idx, end_idx

! If eval_x(i) is beyond the current table interval, advance the
! pointers (and recompute the spline parameters) first; if it is a
! data point, its image is known; otherwise evaluate the spline.
! (Restructured 2026 from the original arithmetic-IF web at labels
! 150/160/170/180; comparisons and evaluation order are unchanged.)
      if (eval_x(i) .gt. table_x(table_idx1)) then
         do
            table_idx = table_idx1
            table_idx1 = table_idx + 1
            if (eval_x(i) .le. table_x(table_idx1)) exit
         end do
         if (eval_x(i) .lt. table_x(table_idx1)) then
! Call choose to determine the appropriate case and then call cases
! to compute the parameters of the spline.
            call choose(table_x(table_idx), table_y(table_idx), &
                 table_slope(table_idx), table_slope(table_idx1), &
                 table_x(table_idx1), table_y(table_idx1), eps_tol, spline_case)
            call cases(table_x(table_idx), table_y(table_idx), &
                 table_slope(table_idx), table_slope(table_idx1), &
                 table_x(table_idx1), table_y(table_idx1), spline_e1, spline_e2, &
                 spline_v1, spline_v2, spline_w1, spline_w2, spline_z1, &
                 spline_z2, spline_y1, spline_y2, spline_case)
         end if
      end if
      if (eval_x(i) .eq. table_x(table_idx1)) then
! If eval_x(i) is a data point, its image is known.
         eval_y(i) = table_y(table_idx1)
         cycle
      end if
      eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, &
           table_x(table_idx), table_y(table_idx), table_x(table_idx1), &
           table_y(table_idx1), spline_y1, spline_y2, spline_e2, &
           spline_w2, spline_v2, spline_case)
      end do

! Calculate the images of the points of evaluation whose abscissas
! are greater than the abscissa of the last data point.
      if (end_idx .eq. num_eval_points) then
         return
      end if
      if ((table_idx1 .eq. num_table_points) .and. (eval_x(end_idx) .ne. table_x(num_table_points))) then
         recompute_tail_params = .false.
      end if
      end if

! Previously, when we arrived at 200 or 210, NUM1 could be improperly
! set. The NUM1= lines below protect from that. llp 8/19/08

! Set the error parameter to indicate that extrapolation has occurred.
      if (recompute_tail_params) then
      err_code = 1
      num_minus1 = max(num_table_points - 1, 1)
      call choose(table_x(num_minus1), table_y(num_minus1), &
           table_slope(num_minus1), table_slope(num_table_points), &
           table_x(num_table_points), table_y(num_table_points), &
           eps_tol, spline_case)
      call cases(table_x(num_minus1), table_y(num_minus1), &
           table_slope(num_minus1), table_slope(num_table_points), &
           table_x(num_table_points), table_y(num_table_points), &
           spline_e1, spline_e2, spline_v1, spline_v2, spline_w1, &
           spline_w2, spline_z1, spline_z2, spline_y1, spline_y2, &
           spline_case)
      end if
      end_idx1 = end_idx + 1
      num_minus1 = max(num_table_points - 1, 1)
      do i = end_idx1, num_eval_points
       eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, &
            table_x(num_minus1), table_y(num_minus1), &
            table_x(num_table_points), table_y(num_table_points), &
            spline_y1, spline_y2, spline_e2, spline_w2, spline_v2, &
            spline_case)
      end do

      return
end subroutine meval

end module quadratic_spline

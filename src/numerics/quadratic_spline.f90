!----------------------------------------------------------------------
! quadratic_spline
!----------------------------------------------------------------------
! The shape-preserving osculatory quadratic spline of ACM Algorithm 574
! (McAllister and Roulier, ACM TOMS 7, 1981; coded by S. L. Dodd and
! M. Roulier, N.C. State University): wherever the data are monotone or
! convex, so is the spline. osplin is the entry point: slopes chooses
! the slope at each point, and meval evaluates the spline, with search
! for the first interval, choose for the knot case of an interval,
! cases for its knots and spline for the value.
! On an interval from (x_left, y_left) to (x_right, y_right) with end
! slopes slope_left and slope_right, the spline is two quadratic arcs
! joined at a knot z (case 1: where the end tangents meet; case 2: the
! midpoint; case 3: between the midpoint and the far end), or three
! arcs with knots y and z (case 4). Each arc is the quadratic Bezier
! curve of its two ends and the control point where their tangents
! meet (v at the left, w at the right, e between case 4's knots).
! With the slopes that slopes computes, every slope is below twice both
! adjacent secants and never opposes their sign, so choose only returns
! cases 1 and 2; 3 and 4 serve other slopes.
module quadratic_spline
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use numerics_status, only: numerics_ok, numerics_bad_input, &
           numerics_unsorted_points
      implicit none
      private
      public :: osplin

! choose's tolerance: an end slope within this fraction of the secant,
! or of twice the secant, counts as equal to it
      real(dp), parameter :: slope_tolerance = 1.d-04

contains

!----------------------------------------------------------------------
! osplin
!----------------------------------------------------------------------
! The spline through xtab/ytab(1:n), n >= 3 points with xtab
! increasing, at the k nondecreasing points xval, into yval; outside the
! table it extrapolates the end intervals' quadratics.
! ierr (numerics_status): numerics_bad_input for n < 3,
! numerics_unsorted_points when xval decreases somewhere; on either
! yval is left as it was.
pure subroutine osplin(xval, yval, xtab, ytab, n, k, ierr)
      integer, intent(in) :: n, k
      real(dp), intent(in) :: xval(k)
      real(dp), intent(inout) :: yval(k)
      real(dp), intent(in) :: xtab(n), ytab(n)
      integer, intent(out), optional :: ierr

      real(dp) :: first_derivs(n)
      integer :: err

      if (present(ierr)) ierr = numerics_ok
      if (n < 3) then
         if (present(ierr)) ierr = numerics_bad_input
         return
      end if
      if (k < 1) return
      call slopes(xtab, ytab, first_derivs, n)
      call meval(xval, yval, xtab, ytab, first_derivs, n, k, slope_tolerance, err)
      if (err == 2 .and. present(ierr)) ierr = numerics_unsorted_points
end subroutine osplin

!----------------------------------------------------------------------
! slopes
!----------------------------------------------------------------------
! The slope at each point, such that the spline needs one knot per
! interval and keeps the data's monotonicity and convexity. At an
! interior point: zero where the two adjacent secants m1 and m2 differ
! in sign or one is zero, else their harmonic mean 2 m1 m2/(m1 + m2),
! formed by extending the steeper secant to the next point's ordinate
! (the two branches round differently). At an end point: twice the end
! secant where the two end secants differ in sign; else the slope from
! the end point to the neighbour's tangent at the end interval's
! midpoint, zero where that opposes the end secant. Every division by
! an interval's width goes through safedivide, which gives zero for a
! zero-width interval (repeated abscissae) instead of an infinity or a
! NaN (KC 2025-05-31 guarded the last interval; 2026-10 the others).
pure subroutine slopes(table_x, table_y, first_derivs, num_points)
      integer, intent(in) :: num_points
      real(dp), intent(in) :: table_x(num_points), table_y(num_points)
      real(dp), intent(out) :: first_derivs(num_points)
      integer :: num_points_m1, prev_idx, idx, next_idx
      real(dp) :: slope1, slope2, x_bar, x_hat, y_diff1, y_diff2, &
           y_x_mid, x_mid, slope1_saved, slope2_saved

      num_points_m1 = num_points - 1
      prev_idx = 1
      idx = 2
      next_idx = 3
! the secants of the first two intervals
      y_diff1 = table_y(2) - table_y(1)
      y_diff2 = table_y(3) - table_y(2)
      call safedivide(y_diff1, (table_x(2) - table_x(1)), slope1)
      slope1_saved = slope1
      call safedivide(y_diff2, (table_x(3) - table_x(2)), slope2)
      slope2_saved = slope2
! the interior points
      do
         if (slope1 == 0.d0 .or. slope2 == 0.d0 .or. (slope1*slope2) <= 0.d0) then
            first_derivs(idx) = 0.d0
         else if (abs(slope1) > abs(slope2)) then
! the line of slope m1 through point idx meets y(next_idx) at x_bar
            x_bar = (y_diff2/slope1) + table_x(idx)
            x_hat = (x_bar + table_x(next_idx))/2.d0
            first_derivs(idx) = y_diff2/(x_hat - table_x(idx))
         else
! the line of slope m2 through point idx meets y(prev_idx) at x_bar
            x_bar = (-y_diff1/slope2) + table_x(idx)
            x_hat = (table_x(prev_idx) + x_bar)/2.d0
            first_derivs(idx) = y_diff1/(table_x(idx) - x_hat)
         end if
         prev_idx = idx
         idx = next_idx
         next_idx = next_idx + 1
         if (idx > num_points_m1) exit
         y_diff1 = table_y(idx) - table_y(prev_idx)
         y_diff2 = table_y(next_idx) - table_y(idx)
         call safedivide(y_diff1, (table_x(idx) - table_x(prev_idx)), slope1)
         call safedivide(y_diff2, (table_x(next_idx) - table_x(idx)), slope2)
      end do
! the last point, from the last two secants
      if ((slope1*slope2) < 0.d0) then
         first_derivs(num_points) = 2.d0*slope2
      else
         x_mid = (table_x(num_points_m1) + table_x(num_points))/2.d0
         y_x_mid = first_derivs(num_points_m1)*(x_mid - table_x(num_points_m1)) + &
              table_y(num_points_m1)
         call safedivide((table_y(num_points) - y_x_mid), &
              (table_x(num_points) - x_mid), first_derivs(num_points))
         if ((first_derivs(num_points)*slope2) < 0.d0) first_derivs(num_points) = 0.d0
      end if
! the first point, from the first two secants
      if ((slope1_saved*slope2_saved) < 0.d0) then
         first_derivs(1) = 2.d0*slope1_saved
         return
      end if
      x_mid = (table_x(1) + table_x(2))/2.d0
      y_x_mid = first_derivs(2)*(x_mid - table_x(2)) + table_y(2)
      call safedivide((y_x_mid - table_y(1)), (x_mid - table_x(1)), first_derivs(1))
      if ((first_derivs(1)*slope1_saved) < 0.d0) first_derivs(1) = 0.d0
end subroutine slopes

!----------------------------------------------------------------------
! safedivide
!----------------------------------------------------------------------
! numerator/denominator, or zero when either is zero (KC 2025-05-31).
pure subroutine safedivide(numerator, denominator, quotient)
      real(dp), intent(in) :: numerator, denominator
      real(dp), intent(out) :: quotient

      quotient = 0.d0
      if (numerator /= 0.d0 .and. denominator /= 0.d0) then
         quotient = numerator / denominator
      end if
end subroutine safedivide

!----------------------------------------------------------------------
! cases
!----------------------------------------------------------------------
! The knots and control points of the spline on (x_left, x_right) for
! knot case spline_case (choose's): the knot (knot_z_x, knot_z_y) and
! the control points v and w of its two arcs for cases 1-3; case 4 adds
! the first knot y and the control point e of the middle arc (unset in
! cases 1-3). Any case other than 2, 3 and 4 is case 1.
pure subroutine cases(x_left, y_left, slope_left, slope_right, x_right, &
     y_right, knot_e_x, knot_e_y, knot_v_x, knot_v_y, knot_w_x, &
     knot_w_y, knot_z_x, knot_z_y, knot_y_x, knot_y_y, spline_case)
      real(dp), intent(in) :: x_left, y_left, slope_left, slope_right, &
           x_right, y_right
      integer, intent(in) :: spline_case
      real(dp), intent(out) :: knot_e_x, knot_e_y, knot_v_x, knot_v_y, &
           knot_w_x, knot_w_y, knot_z_x, knot_z_y, knot_y_x, knot_y_y
      real(dp) :: mbar1, mbar2, mbar3, c1, d1, h1, j1, k1, ztwo

      select case (spline_case)
      case (2)
! the knot at the midpoint, the control points on the end tangents
         knot_z_x = (x_left+x_right)/2.d0
         knot_v_x = (x_left+knot_z_x)/2.d0
         knot_v_y = y_left+slope_left*(knot_v_x-x_left)
         knot_w_x = (knot_z_x+x_right)/2.d0
         knot_w_y = y_right+slope_right*(knot_w_x-x_right)
         knot_z_y = (knot_v_y+knot_w_y)/2.d0
      case (3, 4)
! where the end tangents reach the other end's ordinate (c1, d1), and
! the slopes mbar1, mbar2 of the lines from the ends to their
! reflections h1, j1
         c1 = x_left+(y_right-y_left)/slope_left
         d1 = x_right+(y_left-y_right)/slope_right
         h1 = 2.d0*c1-x_left
         j1 = 2.d0*d1-x_right
         mbar1 = (y_right-y_left)/(h1-x_left)
         mbar2 = (y_left-y_right)/(j1-x_right)
         if (spline_case == 3) then
! one knot, midway between the far end and where those lines meet
            k1 = (y_left-y_right+x_right*mbar2-x_left*mbar1)/(mbar2-mbar1)
            if (abs(slope_left) <= abs(slope_right)) then
               knot_z_x = (k1+x_right)/2.d0
            else
               knot_z_x = (k1+x_left)/2.d0
            end if
            knot_v_x = (x_left+knot_z_x)/2.d0
            knot_v_y = y_left+slope_left*(knot_v_x-x_left)
            knot_w_x = (x_right+knot_z_x)/2.d0
            knot_w_y = y_right+slope_right*(knot_w_x-x_right)
            knot_z_y = knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
                 (knot_z_x-knot_v_x)
         else
! two knots, joined by a line of slope mbar3
            knot_y_x = (x_left+c1)/2.d0
            knot_v_x = (x_left+knot_y_x)/2.d0
            knot_v_y = slope_left*(knot_v_x-x_left) + y_left
            knot_z_x = (d1+x_right)/2.d0
            knot_w_x = (x_right+knot_z_x)/2.d0
            knot_w_y = slope_right*(knot_w_x-x_right) + y_right
            mbar3 = (knot_w_y-knot_v_y)/(knot_w_x-knot_v_x)
            knot_y_y = mbar3*(knot_y_x-knot_v_x) + knot_v_y
            knot_z_y = mbar3*(knot_z_x-knot_v_x) + knot_v_y
            knot_e_x = (knot_y_x+knot_z_x)/2.d0
            knot_e_y = mbar3*(knot_e_x-knot_v_x) + knot_v_y
         end if
      case default
! case 1: the knot where the end tangents meet
         knot_z_x = (y_left-y_right+slope_right*x_right-slope_left*x_left)/ &
              (slope_right-slope_left)
         ztwo = y_left+slope_left*(knot_z_x-x_left)
         knot_v_x = (x_left+knot_z_x)/2.d0
         knot_v_y = (y_left+ztwo)/2.d0
         knot_w_x = (knot_z_x+x_right)/2.d0
         knot_w_y = (ztwo+y_right)/2.d0
         knot_z_y = knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
              (knot_z_x-knot_v_x)
      end select
end subroutine cases

!----------------------------------------------------------------------
! choose
!----------------------------------------------------------------------
! The knot case of the interval from its end slopes and its secant spq:
! 1 when the end tangents meet inside the interval, 2 when the midpoint
! is the knot (MHP: also when an end slope is zero, which would divide
! by zero in cases), 3 or 4 when an end slope exceeds twice the secant
! (one knot or two). eps_tol separates slopes close to the secant or to
! twice it.
pure subroutine choose(x_left, y_left, slope_left, slope_right, x_right, &
     y_right, eps_tol, spline_case)
      real(dp), intent(in) :: x_left, y_left, slope_left, slope_right, &
           x_right, y_right, eps_tol
      integer, intent(out) :: spline_case
      real(dp) :: mref, mref1, mref2, spq, prod, prod1, prod2

      spq = (y_right-y_left)/(x_right-x_left)
      if (slope_left == 0.0d0 .or. slope_right == 0.0d0) then
         spline_case = 2
         return
      end if
      if (spq == 0.d0) then
         if ((slope_left*slope_right) >= 0.d0) then
            spline_case = 2
         else
            spline_case = 1
         end if
         return
      end if

      prod1 = spq*slope_left
      prod2 = spq*slope_right
      mref = abs(spq)
      mref1 = abs(slope_left)
      mref2 = abs(slope_right)

! an end slope opposes the secant: case 2 when both do, else case 1 when
! the other end slope exceeds the secant by more than eps_tol
      if ((prod1 < 0.d0) .or. (prod2 < 0.d0)) then
         if ((prod1 < 0.d0) .and. (prod2 < 0.d0)) then
            spline_case = 2
            return
         end if
         if (prod1 >= 0.d0) then
            if (mref1 > ((1.d0+eps_tol)*mref)) then
               spline_case = 1
            else
               spline_case = 2
            end if
            return
         end if
         if (mref2 > ((1.d0+eps_tol)*mref)) then
            spline_case = 1
         else
            spline_case = 2
         end if
         return
      end if

! both end slopes away from the secant by more than eps_tol, one above
! it and one below: the end tangents meet inside the interval, case 1
      if (abs(spq-slope_left) > eps_tol*mref .and. abs(spq-slope_right) > eps_tol*mref) then
         prod = (mref-mref1)*(mref-mref2)
         if (prod < 0.d0) then
            spline_case = 1
            return
         end if
      end if

! an end slope above twice the secant: case 4 when the other is above
! (2 - eps_tol) times it (neither tangent crosses the midline: two
! knots), else case 3 (one crosses it)
      if (mref1 > (2.d0*mref)) then
         if (mref2 > (2.d0-eps_tol)*mref) then
            spline_case = 4
         else
            spline_case = 3
         end if
         return
      end if
      if (mref2 > (2.d0*mref)) then
         if (mref1 > (2.d0-eps_tol)*mref) then
            spline_case = 4
         else
            spline_case = 3
         end if
         return
      end if

! both tangents cross the midline x = (x_left + x_right)/2 within the
! rectangle of the two points: the midpoint is the knot
      spline_case = 2
end subroutine choose

!----------------------------------------------------------------------
! search
!----------------------------------------------------------------------
! Bisection for eval_point in table_x(1:num_table_points), which it
! lies within: table_idx is the last point below it, or the point it
! equals (found_flag 1).
pure subroutine search(table_x, num_table_points, eval_point, table_idx, &
     found_flag)
      integer, intent(in) :: num_table_points
      real(dp), intent(in) :: table_x(num_table_points), eval_point
      integer, intent(out) :: table_idx, found_flag
      integer :: first_idx, last_idx, middle_idx

      first_idx = 1
      last_idx = num_table_points
      found_flag = 0
      if (table_x(1) == eval_point) then
         table_idx = 1
         found_flag = 1
         return
      end if
      if (table_x(num_table_points) == eval_point) then
         table_idx = num_table_points
         found_flag = 1
         return
      end if
      do
         if ((last_idx-first_idx) == 1) then
            table_idx = first_idx
            return
         end if
         middle_idx = (first_idx+last_idx)/2
         if (table_x(middle_idx) < eval_point) then
            first_idx = middle_idx
         else if (table_x(middle_idx) == eval_point) then
            table_idx = middle_idx
            found_flag = 1
            return
         else
            last_idx = middle_idx
         end if
      end do
end subroutine search

!----------------------------------------------------------------------
! spline
!----------------------------------------------------------------------
! The spline's value at eval_point on (x_left, x_right), from cases'
! knots (y1, y2) (case 4 only) and (z1, z2) and control points v2, w2
! and e2 (case 4 only). MHP 6/87: the chord of the interval where a
! knot coincides with an end.
pure function spline(eval_point, z1, z2, x_left, y_left, x_right, y_right, &
     y1, y2, e2, w2, v2, spline_case)
      real(dp) :: spline
      real(dp), intent(in) :: eval_point, z1, z2, x_left, y_left, &
           x_right, y_right, y1, y2, e2, w2, v2
      integer, intent(in) :: spline_case

      if (spline_case /= 4) then
! one knot: the arc to its left or to its right
         if (z1 < eval_point) then
            if (x_right /= z1) then
               spline = quadratic_arc(z1, z2, w2, x_right, y_right)
            else
               spline = chord()
            end if
         else if (z1 == eval_point) then
            spline = z2
         else
            if (z1 /= x_left) then
               spline = quadratic_arc(x_left, y_left, v2, z1, z2)
            else
               spline = chord()
            end if
         end if
      else
! two knots: the left, middle or right arc
         if (y1 < eval_point) then
            if (z1 < eval_point) then
               if (x_right /= z1) then
                  spline = quadratic_arc(z1, z2, w2, x_right, y_right)
               else
                  spline = chord()
               end if
            else if (z1 == eval_point) then
               spline = z2
            else
               if (z1 /= y1) then
                  spline = quadratic_arc(y1, y2, e2, z1, z2)
               else
                  spline = chord()
               end if
            end if
         else if (y1 == eval_point) then
            spline = y2
         else
            if (y1 /= x_left) then
               spline = quadratic_arc(x_left, y_left, v2, y1, y2)
            else
               spline = chord()
            end if
         end if
      end if

contains

! the quadratic Bezier arc from (a, fa) to (b, fb) with control value c
! at eval_point
      pure function quadratic_arc(a, fa, c, b, fb) result(f)
         real(dp), intent(in) :: a, fa, c, b, fb
         real(dp) :: f
         f = (fa*(b-eval_point)**2+c*2.d0*(eval_point-a)*(b-eval_point) &
              +fb*(eval_point-a)**2)/(b-a)**2
      end function quadratic_arc

! the interval's chord at eval_point
      pure function chord() result(f)
         real(dp) :: f
         real(dp) :: linear_interp_frac
         linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
         f = y_left + linear_interp_frac*(y_right - y_left)
      end function chord
end function spline

!----------------------------------------------------------------------
! meval
!----------------------------------------------------------------------
! The spline through table_x/table_y with slopes table_slope at the
! nondecreasing points eval_x, into eval_y: the points below the table
! on the first interval's arcs, extrapolated; the points in range
! interval by interval (a point on a node takes the node's value); the
! points above on the last interval's arcs, extrapolated. err_code: 0;
! 1 when some points lie outside the table; 2 when eval_x decreases
! somewhere (nothing is evaluated).
pure subroutine meval(eval_x, eval_y, table_x, table_y, table_slope, &
     num_table_points, num_eval_points, eps_tol, err_code)
      integer, intent(in) :: num_table_points, num_eval_points
      real(dp), intent(in) :: eval_x(num_eval_points), &
           table_x(num_table_points), table_y(num_table_points), &
           table_slope(num_table_points), eps_tol
      real(dp), intent(inout) :: eval_y(num_eval_points)
      integer, intent(out) :: err_code

      real(dp) :: spline_v1, spline_v2, spline_w1, spline_w2, &
           spline_z1, spline_z2, spline_y1, spline_y2, spline_e1, spline_e2
      integer :: start_idx, start_idx1, end_idx, found_flag
      integer :: i, table_idx, table_idx1, spline_case, num_minus1
      logical :: recompute_tail_params

      start_idx = 1
      end_idx = num_eval_points
      err_code = 0
      if (num_eval_points /= 1) then
         do i = 1, num_eval_points - 1
            if (eval_x(i+1) >= eval_x(i)) cycle
            err_code = 2
            return
         end do
      end if
! start_idx: the first point at or above table_x(1); end_idx: the last
! at or below the last node
      do i = 1, num_eval_points
         if (eval_x(i) >= table_x(1)) exit
         start_idx = i + 1
      end do
      do i = num_eval_points, 1, -1
         if (eval_x(i) <= table_x(num_table_points)) exit
         end_idx = i - 1
      end do

! the points below the table
      if (start_idx /= 1) then
         err_code = 1
         call choose(table_x(1), table_y(1), table_slope(1), table_slope(2), &
              table_x(2), table_y(2), eps_tol, spline_case)
         call cases(table_x(1), table_y(1), table_slope(1), table_slope(2), &
              table_x(2), table_y(2), spline_e1, spline_e2, spline_v1, &
              spline_v2, spline_w1, spline_w2, spline_z1, spline_z2, &
              spline_y1, spline_y2, spline_case)
         do i = 1, start_idx - 1
            eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, table_x(1), &
                 table_y(1), table_x(2), table_y(2), spline_y1, spline_y2, &
                 spline_e2, spline_w2, spline_v2, spline_case)
         end do
         if (start_idx > num_eval_points) return
      end if

! the points in range, if any (else all the rest lie above the table)
      recompute_tail_params = .true.
      if (start_idx <= end_idx) then
! the interval of the first point in range; points on a node take the
! node's value
         call search(table_x, num_table_points, eval_x(start_idx), &
              table_idx, found_flag)
         table_idx1 = table_idx + 1
         if (found_flag /= 0) then
            do
               eval_y(start_idx) = table_y(table_idx)
               start_idx1 = start_idx
               start_idx = start_idx + 1
               if (start_idx > num_eval_points) return
               if (eval_x(start_idx1) /= eval_x(start_idx)) exit
            end do
! advance to the interval of the next point, giving the points on the
! nodes passed their values; stop after the last point in range (those
! left lie above the table, after one on its last node)
            do
               if (start_idx > end_idx) exit
               if (eval_x(start_idx) < table_x(table_idx1)) exit
               if (eval_x(start_idx) == table_x(table_idx1)) then
                  do
                     eval_y(start_idx) = table_y(table_idx1)
                     start_idx1 = start_idx
                     start_idx = start_idx + 1
                     if (start_idx > num_eval_points) return
                     if (eval_x(start_idx) /= eval_x(start_idx1)) exit
                  end do
               end if
               table_idx = table_idx1
               table_idx1 = table_idx1 + 1
            end do
         end if

! interval by interval: the arcs of the current interval (the first
! interval's are set when points below the table used them), then each
! point's value
         if (start_idx <= end_idx) then
            if (table_idx /= 1 .or. err_code /= 1) then
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
               if (eval_x(i) > table_x(table_idx1)) then
                  do
                     table_idx = table_idx1
                     table_idx1 = table_idx + 1
                     if (eval_x(i) <= table_x(table_idx1)) exit
                  end do
                  if (eval_x(i) < table_x(table_idx1)) then
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
               if (eval_x(i) == table_x(table_idx1)) then
                  eval_y(i) = table_y(table_idx1)
                  cycle
               end if
               eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, &
                    table_x(table_idx), table_y(table_idx), table_x(table_idx1), &
                    table_y(table_idx1), spline_y1, spline_y2, spline_e2, &
                    spline_w2, spline_v2, spline_case)
            end do
         end if
         if (end_idx == num_eval_points) return
! the last interval's arcs are set already when the last point in
! range lay inside that interval
         if ((table_idx1 == num_table_points) .and. &
              (eval_x(end_idx) /= table_x(num_table_points))) recompute_tail_params = .false.
      end if

! the points above the table, on the last interval's arcs (its first
! point is max(n - 1, 1); llp 8/19/08)
      num_minus1 = max(num_table_points - 1, 1)
      if (recompute_tail_params) then
         err_code = 1
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
      do i = end_idx + 1, num_eval_points
         eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, &
              table_x(num_minus1), table_y(num_minus1), &
              table_x(num_table_points), table_y(num_table_points), &
              spline_y1, spline_y2, spline_e2, spline_w2, spline_v2, &
              spline_case)
      end do
end subroutine meval

end module quadratic_spline

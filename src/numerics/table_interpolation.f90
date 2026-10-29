!----------------------------------------------------------------------
! table_interpolation
!----------------------------------------------------------------------
! Interpolation in tables and the searches that bracket a point
! (YREC's numerics_lib, split 2026-10):
!   inter3, interp          Lagrange weights (and their derivatives) for
!                           3 and 4 points; lagrange4 applies 4 weights
!   quint                   quadratic interpolation at equidistant points
!   polint                  Neville's algorithm (Numerical Recipes)
!   lir                     cubic (linear below 4 points) interpolation of
!                           a set of tabulated functions
!   locate, hunt, findex    the interval holding a point: by bisection,
!                           and from a starting guess (Numerical Recipes'
!                           hunt; findex, hunt with YCK's codes outside)
module table_interpolation
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use numerics_status, only: numerics_ok, numerics_bad_input, &
           numerics_coincident_points
      implicit none
      private
      public :: inter3, interp, lagrange4, quint, polint, lir, findex, locate, hunt

contains

!----------------------------------------------------------------------
! inter3
!----------------------------------------------------------------------
! The 3-point Lagrange weights at x_eval for the nodes x_nodes, and
! their derivatives in x_eval: f(x_eval) = sum(weight*f(x_nodes)). For
! the Kurucz atmosphere tables (gravity) and the SCV EOS.
pure subroutine inter3(x_nodes, weight, dweight, x_eval)
      real(dp), intent(in) :: x_nodes(3)
      real(dp), intent(out) :: weight(3), dweight(3)
      real(dp), intent(in) :: x_eval

      real(dp) :: diff32, diff31, diff21
      real(dp) :: denom1, denom2, denom3
      real(dp) :: dx1, dx2, dx3

      diff32 = x_nodes(3) - x_nodes(2)
      diff31 = x_nodes(3) - x_nodes(1)
      diff21 = x_nodes(2) - x_nodes(1)
      denom1 = diff21*diff31
      denom2 = -diff21*diff32
      denom3 = diff31*diff32
      dx1 = x_eval - x_nodes(1)
      dx2 = x_eval - x_nodes(2)
      dx3 = x_eval - x_nodes(3)
      weight(1) = (dx2*dx3)/denom1
      weight(2) = (dx1*dx3)/denom2
      weight(3) = (dx1*dx2)/denom3
      dweight(1) = (dx2 + dx3)/denom1
      dweight(2) = (dx1 + dx3)/denom2
      dweight(3) = (dx1 + dx2)/denom3
end subroutine inter3

!----------------------------------------------------------------------
! interp
!----------------------------------------------------------------------
! The 4-point Lagrange weights at x_eval, and when asked their
! derivatives in x_eval: f(x_eval) = sum(weight*f(x_nodes)), applied by
! lagrange4. For the opacity tables (Alexander 1994 and 2006, OPAL95),
! the SCV EOS, the microscopic diffusion's grid and rotation's
! interfaces.
pure subroutine interp(x_nodes, weight, dweight, x_eval)
      real(dp), intent(in) :: x_nodes(4)
      real(dp), intent(out) :: weight(4)
      real(dp), intent(out), optional :: dweight(4)
      real(dp), intent(in) :: x_eval

      real(dp) :: diff43, diff42, diff41, diff32, diff31, diff21
      real(dp) :: denom1, denom2, denom3, denom4
      real(dp) :: dx1, dx2, dx3, dx4

      diff43 = x_nodes(4) - x_nodes(3)
      diff42 = x_nodes(4) - x_nodes(2)
      diff41 = x_nodes(4) - x_nodes(1)
      diff32 = x_nodes(3) - x_nodes(2)
      diff31 = x_nodes(3) - x_nodes(1)
      diff21 = x_nodes(2) - x_nodes(1)
      denom1 = -diff21*diff31*diff41
      denom2 = diff21*diff32*diff42
      denom3 = -diff31*diff32*diff43
      denom4 = diff41*diff42*diff43
      dx1 = x_eval - x_nodes(1)
      dx2 = x_eval - x_nodes(2)
      dx3 = x_eval - x_nodes(3)
      dx4 = x_eval - x_nodes(4)
      weight(1) = (dx2*dx3*dx4)/denom1
      weight(2) = (dx1*dx3*dx4)/denom2
      weight(3) = (dx1*dx2*dx4)/denom3
      weight(4) = (dx1*dx2*dx3)/denom4
      if (present(dweight)) then
         dweight(1) = (dx3*dx4 + dx2*dx4 + dx2*dx3)/denom1
         dweight(2) = (dx3*dx4 + dx1*dx4 + dx1*dx3)/denom2
         dweight(3) = (dx2*dx4 + dx1*dx4 + dx1*dx2)/denom3
         dweight(4) = (dx2*dx3 + dx1*dx3 + dx1*dx2)/denom4
      end if
end subroutine interp

!----------------------------------------------------------------------
! lagrange4
!----------------------------------------------------------------------
! w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4), summed left to right:
! the 4-point weights of interp applied to a stencil y(k0:k0+3).
pure function lagrange4(w, y) result(s)
      real(dp), intent(in) :: w(4), y(4)
      real(dp) :: s

      s = w(1)*y(1)+w(2)*y(2)+w(3)*y(3)+w(4)*y(4)
end function lagrange4

!----------------------------------------------------------------------
! quint
!----------------------------------------------------------------------
! Quadratic interpolation at equidistant points: y(x) from y0 = y(x0),
! y1 = y(x0+h), y2 = y(x0+2h).
pure subroutine quint(x, x0, h, y0, y1, y2, y)
      real(dp), intent(in) :: x, x0, h, y0, y1, y2
      real(dp), intent(out) :: y

      real(dp) :: d1, d2, t

      d1 = y1 - y0
      d2 = y2 - 2.d0*y1 + y0
      t = (x - x0)/h
      y = y0 + t*d1 + 0.5d0*t*(t-1.d0)*d2
end subroutine quint

!----------------------------------------------------------------------
! polint
!----------------------------------------------------------------------
! Neville's algorithm (Numerical Recipes' polint): the value y at x
! of the polynomial through xa/ya(1:n) and an error estimate dy.
! Coincident abscissae, equal to within coincidence_tolerance of their
! size, return numerics_coincident_points (numerics_status).
pure subroutine polint(xa, ya, n, x, y, dy, ierr)
      integer, intent(in) :: n
      real(dp), intent(in) :: xa(n), ya(n)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: y, dy
      integer, intent(out) :: ierr

      real(dp) :: c(n), d(n)
      integer :: ns, i, j
      real(dp) :: dif, dift, ho, hp, w, den
! two abscissae closer than this, relative to their distances from x,
! count as coincident (a few ulp)
      real(dp), parameter :: coincidence_tolerance = 4.0d0*epsilon(1.0d0)

      ierr = numerics_ok
! the nearest abscissa, where the tableau starts
      ns = 1
      dif = abs(x-xa(1))
      do i = 1, n
         dift = abs(x-xa(i))
         if (dift < dif) then
            ns = i
            dif = dift
         end if
         c(i) = ya(i)
         d(i) = ya(i)
      end do
      y = ya(ns)
      ns = ns - 1
! each column of the tableau: the corrections c and d, and the one that
! leads to y along the path closest to the middle of the table
      do j = 1, n-1
         do i = 1, n-j
            ho = xa(i)-x
            hp = xa(i+j)-x
            w = c(i+1) - d(i)
            den = ho - hp
            if (abs(den) <= coincidence_tolerance*max(abs(ho), abs(hp))) then
               ierr = numerics_coincident_points
               return
            end if
            den = w/den
            d(i) = hp*den
            c(i) = ho*den
         end do
         if (2*ns < n-j) then
            dy = c(ns+1)
         else
            dy = d(ns)
            ns = ns-1
         end if
         y = y + dy
      end do
end subroutine polint

!----------------------------------------------------------------------
! lir
!----------------------------------------------------------------------
! For target_z among the points table_z(1:num_points), increasing or
! decreasing, the num_y tabulated functions table_y(1:num_y, point)
! interpolated into result_y(1:num_y): cubic through four points around
! target_z, interp's Lagrange weights with the one of the first point
! beyond target_z (the last point when target_z lies beyond the table)
! replaced so that the four sum to one, applied by lagrange4; linear between the two
! points around it below four points; the end points' cubic or line
! outside the table. On a node the tabulated values are copied.
! interp_flag is 1 inside the table, 0 outside it. y_stride is
! table_y's leading dimension, size(table_y, 1). A result of zero is
! stored as +0 (a negative zero from the table or the sum becomes +0).
! ierr (numerics_status, optional): numerics_bad_input for fewer than
! two points or a y_stride other than table_y's, numerics_coincident_points
! for equal end points; either leaves result_y unset and interp_flag 0.
! (Formerly lir_impl behind lir, and once in net_lib; a SAVEd search
! index served a resume-the-search option no caller used, and went in
! 2026-10.)
pure subroutine lir(target_z, table_z, result_y, table_y, num_y, y_stride, &
     num_points, interp_flag, ierr)
      real(dp), intent(in) :: target_z
      real(dp), intent(in) :: table_z(:)
      real(dp), intent(out) :: result_y(:)
      real(dp), intent(in) :: table_y(:,:)
      integer, intent(in) :: num_y, y_stride, num_points
      integer, intent(out) :: interp_flag
      integer, intent(out), optional :: ierr

      real(dp) :: weight(4), span, weight_sum
      real(dp) :: frac_lower, frac_upper
      integer :: above, first, corrected, i

      if (present(ierr)) ierr = numerics_ok
      interp_flag = 0
      if (num_points < 2 .or. y_stride /= size(table_y, 1)) then
         if (present(ierr)) ierr = numerics_bad_input
         return
      end if
      span = table_z(num_points) - table_z(1)
      if (span == 0.0d0) then
         if (present(ierr)) ierr = numerics_coincident_points
         return
      end if
! the first point beyond target_z (num_points + 1 when there is none);
! a node gives its values
      interp_flag = 1
      above = 1
      do
         if (above > num_points) then
            interp_flag = 0
            exit
         end if
         if (span < 0.0d0) then
            if (table_z(above) < target_z) exit
         else
            if (table_z(above) > target_z) exit
         end if
         if (table_z(above) == target_z) then
            do i = 1, num_y
               result_y(i) = table_y(i, above)
               if (result_y(i) == 0.d0) result_y(i) = 0.d0
            end do
            return
         end if
         above = above + 1
      end do
      if (above <= 1) interp_flag = 0

      if (num_points >= 4) then
! the four points first..first+3 (first + 2 is the first point beyond
! target_z, kept inside the table), and corrected, the stencil's slot of
! that first point beyond (or of the last point), whose weight takes up
! the four's round-off
         first = above - 2
         corrected = 3
         if (above <= 2) then
            first = 1
            corrected = above
         end if
         if (above >= num_points) then
            first = num_points - 3
            corrected = 4
         end if
         call interp(table_z(first:first+3), weight, x_eval=target_z)
         weight_sum = weight(1) + weight(2) + weight(3) + weight(4)
         weight(corrected) = (1.d0 + weight(corrected)) - weight_sum
         do i = 1, num_y
            result_y(i) = lagrange4(weight, table_y(i, first:first+3))
            if (result_y(i) == 0.d0) result_y(i) = 0.d0
         end do
         return
      end if

! linear between the points above - 1 and above (the end points
! outside the table)
      if (above == 1) above = 2
      if (above > num_points) above = num_points
      frac_lower = (table_z(above) - target_z)/(table_z(above) - table_z(above-1))
      frac_upper = 1.0d0 - frac_lower
      do i = 1, num_y
         result_y(i) = frac_lower*table_y(i, above-1) + frac_upper*table_y(i, above)
         if (result_y(i) == 0.d0) result_y(i) = 0.d0
      end do
end subroutine lir

!----------------------------------------------------------------------
! findex
!----------------------------------------------------------------------
! YCK 3/91, for the OPAL92 and Kurucz 1990 opacities: the index such
! that grid_x(index) <= x_eval < grid_x(index+1) in an increasing grid,
! starting from the incoming index as a guess; below the grid it
! returns -1, at or above its last point -n_grid. It is hunt with these
! two codes for "outside" (hunt's 0 and n), so the two warm-started
! searches are one; findex's own linear walk from the guess went in
! 2026-10 (the index found is the same; a NaN x_eval, which the walk
! sent to -n_grid, now gets whatever interval the bisection ends in).
pure subroutine findex(grid_x, n_grid, x_eval, index)
      integer, intent(in) :: n_grid
      real(dp), intent(in) :: grid_x(n_grid)
      real(dp), intent(in) :: x_eval
      integer, intent(inout) :: index

      call hunt(grid_x, n_grid, x_eval, index)
      if (index == 0) index = -1
      if (index == n_grid) index = -n_grid
end subroutine findex

!----------------------------------------------------------------------
! locate
!----------------------------------------------------------------------
! The j with x between xx(j) and xx(j+1), by bisection (Numerical
! Recipes' locate): xx(j) < x <= xx(j+1) for an increasing table,
! xx(j) >= x > xx(j+1) for a decreasing one. Outside the table j is 0
! or n, except within 1% of the table's span beyond an end, which gives
! the end interval (1 or n-1); so does x on an end node.
pure subroutine locate(xx, n, x, j)
      integer, intent(in) :: n
      real(dp), intent(in) :: xx(n), x
      integer, intent(out) :: j

      integer :: jl, ju, jm
! how far beyond an end, as a fraction of the table's span, x still
! takes the end interval
      real(dp), parameter :: edge_slack = 0.01d0
      real(dp) :: slack

      jl = 0
      ju = n+1
      do while (ju-jl > 1)
         jm = (ju+jl)/2
         if ((xx(n) > xx(1)).eqv.(x > xx(jm))) then
            jl = jm
         else
            ju = jm
         end if
      end do
      j = jl
      slack = edge_slack*abs(xx(n) - xx(1))
      if (j == 0 .and. abs(x - xx(1)) <= slack) then
         j = 1
      else if (j == n .and. abs(x - xx(n)) <= slack) then
         j = n-1
      end if
end subroutine locate

!----------------------------------------------------------------------
! hunt
!----------------------------------------------------------------------
! Numerical Recipes' hunt (Press et al. 1986): the index such that x
! lies between table(index) and table(index+1), for a monotone table,
! increasing or decreasing, starting from the incoming index as a guess
! (outside 1..n it is ignored). From the guess the search steps up or
! down by 1, 2, 4, ... until x is bracketed, then bisects: for a
! sequence of nearby x it costs O(1) per call. For an increasing table,
! table(index) <= x < table(index+1); for a decreasing one,
! table(index) > x >= table(index+1). Outside the table it returns 0
! (before the first point) or n (at or beyond the last). The result does
! not depend on the guess, only the cost does.
pure subroutine hunt(table, n, x, index)
      integer, intent(in) :: n
      real(dp), intent(in) :: table(n)
      real(dp), intent(in) :: x
      integer, intent(inout) :: index

      logical :: ascending
      integer :: index_hi, increment, index_mid

      ascending = table(n) > table(1)
      if (index <= 0 .or. index > n) then
! no useful guess: bisect the whole table
         index = 0
         index_hi = n + 1
      else if ((x >= table(index)) .eqv. ascending) then
! hunt up
         increment = 1
         do
            index_hi = index + increment
            if (index_hi > n) then
               index_hi = n + 1
               exit
            else if ((x >= table(index_hi)) .eqv. ascending) then
               index = index_hi
               increment = increment + increment
            else
               exit
            end if
         end do
      else
! hunt down
         increment = 1
         index_hi = index
         do
            index = index_hi - increment
            if (index < 1) then
               index = 0
               exit
            else if ((x < table(index)) .eqv. ascending) then
               index_hi = index
               increment = increment + increment
            else
               exit
            end if
         end do
      end if
! bisect between index and index_hi
      do
         if (index_hi - index == 1) return
         index_mid = (index_hi + index)/2
         if ((x >= table(index_mid)) .eqv. ascending) then
            index = index_mid
         else
            index_hi = index_mid
         end if
      end do
end subroutine hunt

end module table_interpolation

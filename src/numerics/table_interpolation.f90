!----------------------------------------------------------------------
! table_interpolation
!----------------------------------------------------------------------
! Interpolation in tables and the searches that bracket a point
! (YREC's numerics_lib, split 2026-10):
!   inter3, interp, intrp2  Lagrange weights (and their derivatives) for
!                           3 and 4 points; lagrange4 applies 4 weights
!   quint                   quadratic interpolation at equidistant points
!   polint                  Neville's algorithm (Numerical Recipes)
!   lir                     cubic (linear below 4 points) interpolation of
!                           a set of tabulated functions
!   findex, locate          the interval holding a point
module table_interpolation
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: inter3, interp, intrp2, lagrange4, quint, polint, lir, findex, locate

contains

!----------------------------------------------------------------------
! inter3
!----------------------------------------------------------------------
! The 3-point Lagrange weights at x_eval for the nodes x_nodes, and
! their derivatives in x_eval: f(x_eval) = sum(weight*f(x_nodes)). For
! density in the Livermore opacity tables.
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
! The 4-point Lagrange weights at x_eval and their derivatives in
! x_eval, as inter3. For the VandenBerg opacity tables (CAPPA).
pure subroutine interp(x_nodes, weight, dweight, x_eval)
      real(dp), intent(in) :: x_nodes(4)
      real(dp), intent(out) :: weight(4), dweight(4)
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
      dweight(1) = (dx3*dx4 + dx2*dx4 + dx2*dx3)/denom1
      dweight(2) = (dx3*dx4 + dx1*dx4 + dx1*dx3)/denom2
      dweight(3) = (dx2*dx4 + dx1*dx4 + dx1*dx2)/denom3
      dweight(4) = (dx2*dx3 + dx1*dx3 + dx1*dx2)/denom4
end subroutine interp

!----------------------------------------------------------------------
! intrp2
!----------------------------------------------------------------------
! interp without the derivatives.
pure subroutine intrp2(x_nodes, weight, x_eval)
      real(dp), intent(in) :: x_nodes(4)
      real(dp), intent(out) :: weight(4)
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
end subroutine intrp2

!----------------------------------------------------------------------
! lagrange4
!----------------------------------------------------------------------
! w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4), summed left to right:
! the 4-point weights of interp/intrp2 applied to a stencil y(k0:k0+3).
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
! Coincident abscissae return ierr = 1.
subroutine polint(xa, ya, n, x, y, dy, ierr)
      integer, intent(in) :: n
      real(dp), intent(in) :: xa(n), ya(n)
      real(dp), intent(in) :: x
      real(dp), intent(out) :: y, dy
      integer, intent(out) :: ierr

      real(dp) :: c(n), d(n)
      integer :: ns, i, j
      real(dp) :: dif, dift, ho, hp, w, den

      ierr = 0
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
      do j = 1, n-1
       do i = 1, n-j
          ho = xa(i)-x
          hp = xa(i+j)-x
          w = c(i+1) - d(i)
          den = ho - hp
          if (abs(den) < 1.0d-20) then
             write(*,*) 'polint: repeated abscissa (den ~ 0)'
             ierr = 1
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
! For target_z at the points table_z(1:num_points), increasing or
! decreasing, the num_y tabulated functions table_y(i,n) (first
! dimension y_stride) interpolated into result_y: cubic through the
! four nearest points, with the weight of the closest corrected so the
! four sum to one; linear below four points. interp_flag is 1 inside
! the table, 0 when extrapolating. On a node the tabulated values are
! copied. Fewer than two points, or equal end points, leave result_y
! unset. (Formerly lir_impl behind lir, and once in net_lib; a SAVEd
! search index served a resume-the-search option no caller used, and
! went in 2026-10.)
pure subroutine lir(target_z,table_z,result_y,table_y,num_y,y_stride, &
     num_points,interp_flag)
      real(dp), intent(in) :: target_z
      real(dp), intent(in) :: table_z(*)
      real(dp), intent(out) :: result_y(*)
      real(dp), intent(in) :: table_y(*)
      integer, intent(in) :: num_y, y_stride, num_points
      integer, intent(out) :: interp_flag

      real(dp) :: weight(4)
      integer :: search_idx
      integer :: linear_mode
      integer :: stride, stride_m1, y_strided, num_y_strided, table_end
      real(dp) :: diff
      integer :: pivot, closest
      real(dp) :: y1, y2, y3, y4, z1, z2, z3, z4, z12, z34
      integer :: y_idx, j, k, base_idx
      real(dp) :: yy

      linear_mode=0
      stride=1
! CHECK NT AND RESET IL IF NECESSARY
      if (num_points < 2) then
         return
      end if
      if(num_points < 4) linear_mode=1
! ADDRESSING CONSTANTS
      interp_flag=1
      stride_m1=stride-1
      y_strided=stride*y_stride
      num_y_strided=(num_y-1)*stride+1
      table_end=(num_points-1)*stride+1
      diff=table_z(table_end)-table_z(1)
! THE SEARCH STARTS AT THE FIRST POINT
      search_idx=1
! DETERMINE POSITION OF target_z WITHIN table_z
      do
         if(search_idx > table_end) then
            interp_flag=0
            exit
         end if
         if (diff == 0.0) then
            return
         end if
         if (diff < 0.0) then
            if (table_z(search_idx) < target_z) exit
         else
            if (table_z(search_idx) > target_z) exit
         end if
         if (table_z(search_idx) == target_z) then
! SET Y WHEN Z LIES ON A MESH POINT
            base_idx=(search_idx-1)*y_stride
            do y_idx=1,num_y_strided
               result_y(y_idx)=table_y(y_idx+base_idx)
               if(result_y(y_idx) == 0.d0) result_y(y_idx+stride_m1)=0.d0
            end do
            return
         end if
         search_idx=search_idx+stride
      end do
! CONTROL WHEN Z DOES NOT LIE ON A MESH POINT
      if(search_idx <= 1) interp_flag=0
      if (linear_mode /= 1) then
! CUBIC INTERPOLATION/EXTRAPOLATION
! PIVOTAL POINT (M) AND POINT (K) CLOSEST TO Z
      pivot=search_idx
      closest=3
      if (search_idx <= 1+stride) then
      pivot=1+stride+stride
      closest=search_idx
      end if
      if (search_idx >= table_end) then
      pivot=table_end-stride
      closest=4
      end if
! WEIGHTING FACTORS
      y1=table_z(pivot-stride*2)
      y2=table_z(pivot-stride)
      y3=table_z(pivot)
      y4=table_z(pivot+stride)
      z1=target_z-y1
      z2=target_z-y2
      z3=target_z-y3
      z4=target_z-y4
      z12=z1*z2
      z34=z3*z4
      weight(1)=z2*z34/((y1-y2)*(y1-y3)*(y1-y4))
      weight(2)=z1*z34/((y2-y1)*(y2-y3)*(y2-y4))
      weight(3)=z12*z4/((y3-y1)*(y3-y2)*(y3-y4))
      weight(4)=z12*z3/((y4-y1)*(y4-y2)*(y4-y3))
! CORRECT A(K)
      diff=weight(1)+weight(2)+weight(3)+weight(4)
      weight(closest)=(1.d0+weight(closest))-diff
! COMPUTE Y
      pivot=(pivot-1)/stride-3
      pivot=pivot*y_strided
      do y_idx=1,num_y_strided
         k=y_idx+pivot
         yy=0.d0
         do j=1,4
            k=k+y_strided
            diff=table_y(k)
            yy=yy+weight(j)*diff
         end do
         result_y(y_idx)=yy
         if(result_y(y_idx) == 0.d0) result_y(y_idx+stride_m1)=0.d0
      end do
      return
! LINEAR INTERPOLATION/EXTRAPOLATION
      end if
      if(search_idx == 1) search_idx=1+stride
      if(search_idx > table_end) search_idx=table_end
      z1=table_z(search_idx)
      y1=(z1-target_z)/(z1-table_z(search_idx-stride))
      y2=1.0d0-y1
      base_idx=(search_idx-1)*y_stride
      pivot=base_idx-y_strided
      do y_idx=1,num_y_strided,stride
         result_y(y_idx)=y1*table_y(y_idx+pivot)+y2*table_y(y_idx+base_idx)
         if(result_y(y_idx) == 0.d0) result_y(y_idx+stride_m1)=0.d0
      end do
end subroutine lir

!----------------------------------------------------------------------
! findex
!----------------------------------------------------------------------
! YCK 3/91, for the Cox and Kurucz opacities: index such that
! grid_x(index) <= x_eval < grid_x(index+1), searching outward from the
! incoming index (reset to 1 when outside [1, n_grid]). Below the grid
! it returns -1, at or above its last point -n_grid.
pure subroutine findex(grid_x, n_grid, x_eval, index)
      integer, intent(in) :: n_grid
      real(dp), intent(in) :: grid_x(n_grid)
      real(dp), intent(in) :: x_eval
      integer, intent(inout) :: index

      integer :: found_index, j

! find the 'index'
      if(index < 1.or.index > n_grid)index=1
      found_index=index
      if(x_eval < grid_x(found_index))then
         do j=found_index-1,1,-1
            if(grid_x(j) <= x_eval)then
               found_index=j
               index=found_index
               
               return
            endif
         end do
         found_index=-1
      else
         do j=found_index,n_grid-1
            if(grid_x(j+1) > x_eval)then
               found_index=j
               index=found_index
               
               return
            endif
         end do
         found_index = -n_grid
      endif
      index=found_index
end subroutine findex

!----------------------------------------------------------------------
! locate
!----------------------------------------------------------------------
! The j with x between xx(j) and xx(j+1), by bisection, for xx
! increasing or decreasing; j = 0 or n when outside, except within 1%
! of an end, which clamps to the end interval.
pure subroutine locate(xx, n, x, j)
      integer, intent(in) :: n
      real(dp), intent(in) :: xx(n), x
      integer, intent(out) :: j

      integer :: jl, ju, jm
! slack allowed when x falls just outside the table before the
! bracket is clamped to the end intervals
      real(dp), parameter :: below_table_slack = 0.99d0, &
           above_table_slack = 1.01d0

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
      if ((j == 0) .and. (x > below_table_slack*xx(1))) then
         j = 1
      end if
      if ((j == n) .and. (x < above_table_slack*xx(n))) then
         j = n-1
      end if
end subroutine locate

end module table_interpolation

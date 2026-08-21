!----------------------------------------------------------------------
! cspline
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original cspline.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! double precision version for opacities. Taken from Numerical
! Recipes, Press et al, p88. Given arrays x and y of length n
! containing a tabulated function, i.e. y(i) = f(x(i)) with
! x(1) < x(2) < ... < x(n), and given values yp1 and ypn for the first
! derivative of the interpolating function at points 1 and n,
! respectively, this routine returns an array y2 of length n which
! contains the second derivatives of the interpolating function at the
! tabulated points x(i). If yp1 and/or ypn are equal to 1.0e30 or
! larger, the routine is signalled to set the corresponding boundary
! condition for a natural spline, with zero second derivative on that
! boundary.
module cspline_mod
contains
subroutine cspline(x, y, n, yp1, ypn, y2)
      implicit none
      integer, parameter :: json = 5000

      integer, intent(in) :: n
      double precision, intent(in) :: x(n), y(n), yp1, ypn
      double precision, intent(out) :: y2(n)

      double precision :: u(json)
      integer :: i, k
      double precision :: sig, p, qn, un
      save

      if (yp1 .gt. 0.99d30) then
         y2(1) = 0.0d0
         u(1) = 0.0d0
      else
         y2(1) = -0.5d0
         u(1) = (3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      end if
      do i = 2, n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.0d0
         y2(i) = (sig-1.0d0)/p
         u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn .gt. 0.99d30) then
         qn = 0.0d0
         un = 0.0d0
      else
         qn = 0.5d0
         un = (3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      end if
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do k = n-1, 1, -1
         y2(k) = y2(k)*y2(k+1)+u(k)
      end do
      return
end subroutine cspline
end module cspline_mod

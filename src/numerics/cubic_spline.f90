!----------------------------------------------------------------------
! cubic_spline
!----------------------------------------------------------------------
! Cubic splines through tabulated points (YREC's numerics_lib, split
! 2026-10):
!   cspline, splint   Numerical Recipes' spline/splint: the second
!                     derivatives for given end slopes (a slope above
!                     0.99d30 makes that end natural), and the value at x,
!                     with its first derivative when asked
!   natural_spline    cspline with both ends natural (formerly kspline,
!                     for four points, and splinc)
!   compute_spline_coefficients
!                     the spline's cubic on each interval in powers of
!                     x - x(i), for tables that store it (OPAL92 and
!                     Kurucz 1990 opacities)
! (ysplin and intpol, YCK's 3/91 natural spline in slope form, wrong for
! 2 and 3 points, became wrappers of these in 2026-10 and went once their
! callers had moved to them.)
module cubic_spline
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use numerics_status, only: numerics_ok, numerics_coincident_points, &
           numerics_bad_input
      implicit none
      private
      public :: cspline, natural_spline, splint, compute_spline_coefficients

contains

!----------------------------------------------------------------------
! cspline
!----------------------------------------------------------------------
! Numerical Recipes' spline (p. 88; double precision for the
! opacities): the second derivatives y2 at the tabulated points x(1:n),
! n >= 2, strictly monotone (increasing, as splint needs, or
! decreasing), of the cubic spline through y with first derivatives yp1
! and ypn at the ends. yp1 or ypn above 0.99d30 makes that end natural
! (zero second derivative).
! ierr (numerics_status, optional): numerics_bad_input for n < 2 or
! abscissae that are not monotone, numerics_coincident_points for two
! equal neighbours; y2 is then not set.
pure subroutine cspline(x, y, n, yp1, ypn, y2, ierr)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), y(n), yp1, ypn
      real(dp), intent(out) :: y2(n)
      integer, intent(out), optional :: ierr

      real(dp) :: u(n)
      integer :: i, k
      real(dp) :: sig, p, qn, un

      if (present(ierr)) ierr = numerics_ok
      if (n < 2) then
         if (present(ierr)) ierr = numerics_bad_input
         return
      end if
      do i = 1, n - 1
         if (x(i+1) == x(i)) then
            if (present(ierr)) ierr = numerics_coincident_points
            return
         end if
         if ((x(i+1) > x(i)) .neqv. (x(2) > x(1))) then
            if (present(ierr)) ierr = numerics_bad_input
            return
         end if
      end do
      if (yp1 > 0.99d30) then
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
      if (ypn > 0.99d30) then
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
end subroutine cspline

!----------------------------------------------------------------------
! natural_spline
!----------------------------------------------------------------------
! cspline with natural boundary conditions at both ends (formerly
! kspline, its four-point form, and splinc); ierr (optional) is
! cspline's.
pure subroutine natural_spline(x, y, n, y2, ierr)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), y(n)
      real(dp), intent(out) :: y2(n)
      integer, intent(out), optional :: ierr
! cspline reads a boundary slope above 0.99d30 as natural
      real(dp), parameter :: natural_bc = 1.0d30

      call cspline(x, y, n, natural_bc, natural_bc, y2, ierr)
end subroutine natural_spline

!----------------------------------------------------------------------
! splint
!----------------------------------------------------------------------
! Numerical Recipes' splint (p. 89; MHP 10/02): the cubic spline
! through xa/ya, with the second derivatives y2a from cspline, at x, and
! its first derivative dy when asked; klo/khi optionally return the
! bracketing points. Outside the table it extrapolates the end
! interval's cubic. ierr (numerics_status): numerics_coincident_points
! when the bracketing points coincide (formerly also ksplint, for four
! points, and splintd2).
subroutine splint(xa, ya, n, y2a, x, y, ierr, klo, khi, dy)
      integer, intent(in) :: n
      real(dp), intent(in) :: xa(n), ya(n), y2a(n), x
      real(dp), intent(out) :: y
      integer, intent(out) :: ierr
      integer, intent(out), optional :: klo, khi
      real(dp), intent(out), optional :: dy

      integer :: lo, hi, k
      real(dp) :: h, a, b

      ierr = numerics_ok
      lo = 1
      hi = n
      do while (hi-lo > 1)
         k = (hi+lo)/2
         if (xa(k) > x) then
            hi = k
         else
            lo = k
         end if
      end do
      if (present(klo)) klo = lo
      if (present(khi)) khi = hi
      h = xa(hi) - xa(lo)
      if (h == 0d0) then
         ierr = numerics_coincident_points
         return
      end if
      a = (xa(hi)-x)/h
      b = (x - xa(lo))/h
      y = a*ya(lo)+b*ya(hi)+ &
           ((a**3-a)*y2a(lo)+(b**3-b)*y2a(hi))*(h**2)/6.0d0
      if (present(dy)) dy = (ya(hi)-ya(lo))/h &
           - ((3.0d0*a**2-1.0d0)*y2a(lo) - (3.0d0*b**2-1.0d0)*y2a(hi))*h/6.0d0
end subroutine splint

!----------------------------------------------------------------------
! compute_spline_coefficients
!----------------------------------------------------------------------
! The cubic spline through x/y(1:n) with second derivatives y2 (from
! cspline or natural_spline) in powers of x - x(i): on [x(i), x(i+1)]
!   s(x) = c(1,i) + c(2,i) d + c(3,i) d**2 + c(4,i) d**3,  d = x - x(i),
! with c(1,i) = y(i), c(2,i) = s'(x(i)), c(3,i) = y2(i)/2 and c(4,i) =
! (y2(i+1) - y2(i))/(6 h), h = x(i+1) - x(i). Column n holds the value,
! slope and half the second derivative at the last point (c(4,n) = 0),
! for a linear extrapolation beyond it.
pure subroutine compute_spline_coefficients(x, y, y2, n, c)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), y(n), y2(n)
      real(dp), intent(out) :: c(4,n)

      integer :: i
      real(dp) :: h

      do i = 1, n-1
         h = x(i+1) - x(i)
         c(1,i) = y(i)
         c(2,i) = (y(i+1)-y(i))/h - h*(2.0d0*y2(i) + y2(i+1))/6.0d0
         c(3,i) = 0.5d0*y2(i)
         c(4,i) = (y2(i+1) - y2(i))/(6.0d0*h)
      end do
      c(1,n) = y(n)
      if (n > 1) then
         h = x(n) - x(n-1)
         c(2,n) = (y(n)-y(n-1))/h + h*(y2(n-1) + 2.0d0*y2(n))/6.0d0
      else
         c(2,n) = 0.0d0
      end if
      c(3,n) = 0.5d0*y2(n)
      c(4,n) = 0.0d0
end subroutine compute_spline_coefficients

end module cubic_spline

!----------------------------------------------------------------------
! cubic_spline
!----------------------------------------------------------------------
! Cubic splines through tabulated points (YREC's numerics_lib, split
! 2026-10):
!   cspline, splint   Numerical Recipes' spline/splint: the second
!                     derivatives for given end slopes (a slope of 1e30 or
!                     more means a natural end), and the value at x
!   natural_spline    cspline with both ends natural (formerly kspline,
!                     for four points, and splinc)
!   ysplin, intpol    YCK's (3/91) natural spline in coefficient form,
!                     for the Cox and Kurucz opacity tables: the cubic's
!                     coefficients on every interval, and the value and
!                     slope at x
module cubic_spline
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use io_units, only: run_log_unit, terminal_unit
      implicit none
      private
      public :: cspline, natural_spline, splint, ysplin, intpol

contains

!----------------------------------------------------------------------
! cspline
!----------------------------------------------------------------------
! Numerical Recipes' spline (p. 88; double precision for the
! opacities): the second derivatives y2 at the tabulated points
! x(1) < ... < x(n) of the cubic spline through y with first
! derivatives yp1 and ypn at the ends. yp1 or ypn of 0.99d30 or more
! makes that end natural (zero second derivative).
subroutine cspline(x, y, n, yp1, ypn, y2)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), y(n), yp1, ypn
      real(dp), intent(out) :: y2(n)

      real(dp) :: u(n)
      integer :: i, k
      real(dp) :: sig, p, qn, un

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

!----------------------------------------------------------------------
! natural_spline
!----------------------------------------------------------------------
! cspline with natural boundary conditions at both ends (formerly
! kspline, its four-point form, and splinc).
subroutine natural_spline(x, y, n, y2)
      integer, intent(in) :: n
      real(dp), intent(in) :: x(n), y(n)
      real(dp), intent(out) :: y2(n)
! cspline reads a boundary slope of 0.99d30 or more as natural
      real(dp), parameter :: natural_bc = 1.0d30

      call cspline(x, y, n, natural_bc, natural_bc, y2)
end subroutine natural_spline

!----------------------------------------------------------------------
! splint
!----------------------------------------------------------------------
! Numerical Recipes' splint (p. 89; MHP 10/02): the cubic spline
! through xa/ya, with the second derivatives y2a from cspline, at x;
! klo/khi optionally return the bracketing points. Coincident
! bracketing points return ierr = 1 (formerly also ksplint, for four
! points, and splintd2).
subroutine splint(xa, ya, n, y2a, x, y, ierr, klo, khi)
      integer, intent(in) :: n
      real(dp), intent(in) :: xa(n), ya(n), y2a(n), x
      real(dp), intent(out) :: y
      integer, intent(out) :: ierr
      integer, intent(out), optional :: klo, khi

      integer :: lo, hi, k
      real(dp) :: h, a, b

      ierr = 0
      lo = 1
      hi = n
      do while (hi-lo .gt. 1)
         k = (hi+lo)/2
         if (xa(k) .gt. x) then
            hi = k
         else
            lo = k
         end if
      end do
      if (present(klo)) klo = lo
      if (present(khi)) khi = hi
      h = xa(hi) - xa(lo)
      if (h .eq. 0d0) then
           write(run_log_unit,*) 'ERROR IN SPLINT ROUTINE.'
         ierr = 1
         return
      end if
      a = (xa(hi)-x)/h
      b = (x - xa(lo))/h
      y = a*ya(lo)+b*ya(hi)+ &
           ((a**3-a)*y2a(lo)+(b**3-b)*y2a(hi))*(h**2)/6.0d0
      return
end subroutine splint

!----------------------------------------------------------------------
! ysplin
!----------------------------------------------------------------------
! YCK 3/91: the natural cubic spline through xi(1:n) in coefficient
! form. On entry c(1,i) holds the values; on exit c(2,i), c(3,i) and
! c(4,i) hold the first-, second- and third-order coefficients of the
! cubic on [xi(i), xi(i+1)].
subroutine ysplin(xi, c, n)
      integer, parameter :: np = 100

      real(dp), intent(in) :: xi(n)
      real(dp), intent(inout) :: c(4,np)
      integer, intent(in) :: n

      real(dp) :: f(np), h(np), d(np), g, d3
      real(dp) :: const1, const2, const3, const4
      integer :: i

! set the divided difference at each subinterval.
      do i = 2, n
         h(i) = xi(i)-xi(i-1)
         f(i) = (c(1,i)-c(1,i-1))/h(i)
      end do
! set the tridiagonally dominant matrix and the matrix equation
! for the natural spline.
!!!!! find the coefficients for the second order terms !!!!!
      i = 2
      const1 = 3.0d0
      const2 = 1.5d0
      const3 = 2.0d0
      const4 = 1.5d0
      c(2,i) = const1*h(i)*f(i+1)+const2*h(i+1)*f(i)
      d(i) = const3*h(i)+const4*h(i+1)
      do i = 3, n-2
         const1 = 3.0d0
         const2 = 3.0d0
         const3 = 2.0d0
         const4 = 2.0d0
         c(2,i) = const1*h(i)*f(i+1)+const2*h(i+1)*f(i)
         d(i) = const3*h(i)+const4*h(i+1)
      end do
      i = n-1
      const1 = 1.5d0
      const2 = 3.0d0
      const3 = 1.5d0
      const4 = 2.0d0
      c(2,i) = const1*h(i)*f(i+1)+const2*h(i+1)*f(i)
      d(i) = const3*h(i)+const4*h(i+1)
! solve the matrix equation with gauss method
! elimination of the sub-diagonal
      do i = 2, n-2
         g = h(i+2)/d(i)
         d(i+1) = d(i+1)-g*h(i)
         c(2,i+1) = c(2,i+1)-g*c(2,i)
      end do
! elimination of the super-diagonal
      c(2,n-1) = c(2,n-1)/d(n-1)
      do i = n-2, 2, -1
         c(2,i) = (c(2,i)-h(i)*c(2,i+1))/d(i)
      end do
! treatment for the first and last row.
      c(2,1) = 1.5d0*f(2)-c(2,2)/2.0d0
      c(2,n) = 1.5d0*f(n)-c(2,n-1)/2.0d0
!!! now, we have the coefficients for the first order terms
!
! find the coefficients for the second, and the third order terms
      do i = 1, n-1
         d3 = c(2,i)+c(2,i+1)-2.0d0*f(i+1)
         c(3,i) = (f(i+1)-c(2,i)-d3)/h(i+1)
         c(4,i) = d3/(h(i+1)*h(i+1))
      end do
!!!   now, we have complete set of coefficients for each sub-interval.

      return
end subroutine ysplin

!----------------------------------------------------------------------
! intpol
!----------------------------------------------------------------------
! YCK 3/91: the value y_eval and slope dy_eval at x_eval of ysplin's
! spline through x_grid/y_grid (n_grid <= 100 points); outside the grid
! it extrapolates the end cubic.
subroutine intpol(x_grid, y_grid, n_grid, x_eval, y_eval, dy_eval, ierr)
      integer, parameter :: np=100

      real(dp), intent(in) :: x_grid(n_grid), y_grid(n_grid)
      integer, intent(in) :: n_grid
      real(dp), intent(in) :: x_eval
      real(dp), intent(out) :: y_eval, dy_eval

      real(dp) :: spline_coeff(4,np), dx_local, x_eval_copy
      real(dp) :: interp_value, interp_deriv
      integer :: i, k_lo, k_hi, k_mid

      integer, intent(out) :: ierr

      ierr = 0
! the coefficients for the zero-th order term
      do i=1,n_grid
         spline_coeff(1,i)=y_grid(i)
      end do
! find the coefficients for spline interploation
      call ysplin(x_grid, spline_coeff, n_grid)
      x_eval_copy=x_eval
! find the grid point,  k_lo, such that x_grid(k_lo)<=x_eval_copy, and
! abs(x_grid(k_lo)-x_eval_copy)<1.
      if(x_grid(1).gt.x_eval)then
         k_lo=1
         k_hi=2
      else if(x_grid(n_grid).lt.x_eval)then
         k_lo=n_grid-1
         k_hi=n_grid
      else
      k_lo=1
      k_hi=n_grid
    do while ((k_hi-k_lo).gt.1)
         k_mid = (k_hi+k_lo)/2
         if(x_grid(k_mid).gt.x_eval_copy)then
            k_hi=k_mid
         else
            k_lo=k_mid
         endif
    end do
      if((k_hi-k_lo).le.0)then
         write(terminal_unit, *) 'ERROR COX OP: INTERPOLATION'
         write(run_log_unit, *) 'ERROR COX OP: INTERPOLATION'
         ierr = 1
         return
      endif
      end if
! now, (k_lo,k_hi) is sub-range of x_grid which contains x_eval_copy.
      dx_local=x_eval_copy-x_grid(k_lo)
! go on to the spline interpolation routine.
! evaluates the interpolation value in the sub-range we determined.
      interp_value=((spline_coeff(4,k_lo)*dx_local+spline_coeff(3,k_lo))*dx_local &
           +spline_coeff(2,k_lo))*dx_local+spline_coeff(1,k_lo)
      interp_deriv=(3.0d0*spline_coeff(4,k_lo)*dx_local+2.0d0*spline_coeff(3,k_lo)) &
           *dx_local+spline_coeff(2,k_lo)
! return the results from the spline routine
      y_eval=interp_value
      dy_eval=interp_deriv

      return
end subroutine intpol

end module cubic_spline

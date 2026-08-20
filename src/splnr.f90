!----------------------------------------------------------------------
! splnr
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original splnr.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Single-precision (real, not double precision) natural/clamped cubic
! spline coefficient generator, the classic Numerical Recipes SPLINE
! routine kept in its original real precision (unlike cspline/splinj/
! splinc, which are real*8 ports of the same algorithm).
subroutine splnr(x, y, n, yp1, ypn, y2)
      implicit none
      integer, parameter :: nmax = 500

      integer, intent(in) :: n
      real, intent(in) :: x(n), y(n), yp1, ypn
      real, intent(out) :: y2(n)

      integer :: i, k
      real :: p, qn, sig, un, u(nmax)

      if (yp1.gt..99e30) then
        y2(1)=0.
        u(1)=0.
      else
        y2(1)=-0.5
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1))) &
             /(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      if (ypn.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      end do
      return
end subroutine splnr

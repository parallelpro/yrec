!----------------------------------------------------------------------
! splinj
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original splinj.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Natural cubic spline coefficient generator, sized exactly to n
! (unlike splinc, which is dimensioned to the json=5000 maximum).
subroutine splinj(x, y, y2, n)
      implicit none

      integer, intent(in) :: n
      double precision, intent(in) :: x(n), y(n)
      double precision, intent(out) :: y2(n)

      double precision :: u(n)
      integer :: i, k
      double precision :: sig, p, qn, un

! natural spline
      y2(1) = 0.0d0
      u(1) = 0.0d0
      do i = 2, n-1
         sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
         p = sig*y2(i-1)+2.0d0
         y2(i) = (sig-1.0d0)/p
         u(i) = (6.0d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) &
              /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      end do
      qn = 0.0d0
      un = 0.0d0
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)
      do k = n-1, 1, -1
         y2(k) = y2(k)*y2(k+1)+u(k)
      end do
      return
end subroutine splinj

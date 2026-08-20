!----------------------------------------------------------------------
! splinnr
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original splinnr.f; only variable names, source form, and comment
! style were updated. Not called from any already-converted file in
! this batch; a Numerical-Recipes-style cubic spline setup routine
! (single precision, as in the original -- preserved verbatim).
!
! Given the tabulated function y_values(1:num_points) on strictly
! increasing abscissas x_values(1:num_points), and first derivatives
! first_deriv_start/first_deriv_end at the endpoints (or a value
! greater than 0.99e30 to signal a natural spline boundary there),
! returns the second derivatives second_derivs(1:num_points) of the
! interpolating cubic spline.
subroutine splinnr(x_values, y_values, num_points, first_deriv_start, &
     first_deriv_end, second_derivs)

      implicit none
      integer, parameter :: nmax = 500

      integer, intent(in) :: num_points
      real, intent(in) :: first_deriv_start, first_deriv_end
      real, intent(in) :: x_values(num_points), y_values(num_points)
      real, intent(out) :: second_derivs(num_points)

      integer :: i, k
      real :: p, qn, sig, un, u(nmax)

      if (first_deriv_start.gt..99e30) then
        second_derivs(1)=0.
        u(1)=0.
      else
        second_derivs(1)=-0.5
        u(1)=(3./(x_values(2)-x_values(1)))*((y_values(2)-y_values(1))/ &
             (x_values(2)-x_values(1))-first_deriv_start)
      endif
      do 11 i=2,num_points-1
        sig=(x_values(i)-x_values(i-1))/(x_values(i+1)-x_values(i-1))
        p=sig*second_derivs(i-1)+2.
        second_derivs(i)=(sig-1.)/p
        u(i)=(6.*((y_values(i+1)-y_values(i))/(x_values(i+1)-x_values(i))- &
             (y_values(i)-y_values(i-1))/(x_values(i)-x_values(i-1)))/ &
             (x_values(i+1)-x_values(i-1))-sig*u(i-1))/p
11    continue
      if (first_deriv_end.gt..99e30) then
        qn=0.
        un=0.
      else
        qn=0.5
        un=(3./(x_values(num_points)-x_values(num_points-1)))* &
           (first_deriv_end-(y_values(num_points)-y_values(num_points-1))/ &
           (x_values(num_points)-x_values(num_points-1)))
      endif
      second_derivs(num_points)=(un-qn*u(num_points-1))/(qn*second_derivs(num_points-1)+1.)
      do 12 k=num_points-1,1,-1
        second_derivs(k)=second_derivs(k)*second_derivs(k+1)+u(k)
12    continue
      return
end subroutine splinnr

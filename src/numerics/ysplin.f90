!----------------------------------------------------------------------
! ysplin
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ysplin.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! YCK 3/91. Find the coefficients for the natural cubic spline.
! On entry c(1,i) holds the tabulated function values at xi(i); on
! exit c(2,i)/c(3,i)/c(4,i) hold the first/second/third-order
! coefficients of the cubic on each sub-interval.
subroutine ysplin(xi, c, n)
      implicit none
      integer, parameter :: np = 100

      double precision, intent(in) :: xi(n)
      double precision, intent(inout) :: c(4,np)
      integer, intent(in) :: n

      double precision :: f(np), h(np), d(np), g, d3
      double precision :: const1, const2, const3, const4
      integer :: i
      save

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

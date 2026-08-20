!----------------------------------------------------------------------
! boole
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original boole.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Performs 5-point Newton-Cotes (Boole's rule) integration of
! tabulated data.
!      x: vector of independent variable values (x)
!      y: vector of dependent variable values (f(x))
!      n: number of elements in x
!      n_grid: number of elements desired in the interpolated grid
!              (should be 1+4*integer)
! integral: returned, integrated value of the function (kept as a
!           length-1 array, exactly as in the original file)
!
! Data need to be evenly gridded for the integration formula to
! work. Data that are not already so are interpolated onto an even
! grid using spline interpolation (Numerical Recipes SPLINE, renamed
! splinj here, and SPLINT).
subroutine boole(x, y, n, n_grid, integral)
      implicit none
      double precision, parameter :: scalex = 1e-11
      double precision, parameter :: scaley = 1e7

      double precision, intent(in) :: x(n), y(n)
      integer, intent(in) :: n, n_grid
      double precision, intent(out) :: integral(1)

      double precision :: h_step(1), y2_deriv(n)
      double precision :: x_even(n_grid), y_even(n_grid)
      double precision :: x_scaled(n), y_scaled(n)
      integer :: i, num_quads, klo, khi

! rescale radius and cs vectors to have values ~1
      do i = 1, n
            x_scaled(i) = x(i)*scalex
            y_scaled(i) = y(i)*scaley
      end do
! deal with unevenly gridded datasets:
            call splinj(x_scaled, y_scaled, y2_deriv, n) ! get derivs of interp. fn.
            do i = 1, n_grid
                  x_even(i) = x_scaled(1)+(i-1)*(x_scaled(n)-x_scaled(1))/(n_grid-1)
                  call splint(x_scaled, y_scaled, n, y2_deriv, x_even(i), y_even(i), klo, khi)
            end do

! how many sets of four points do we have?
      num_quads = (n_grid-1)/4
      h_step = (x_scaled(n)-x_scaled(1))/(n_grid-1)
! for each set of 4 points before the last, apply formula and add up
      integral = 0.0d0
      do i = 0, num_quads-1
            integral = integral+2.0*h_step*(7*y_even(1+4*i)+ 32*y_even(2+4*i) + 12*y_even(3+4*i) + &
                  32*y_even(4+4*i) + 7*y_even(5+4*i))/45.0
      end do
! rescale result back to actual units:
      integral = integral/(scalex*scaley)
!      print*, 'Last point in vector:', n_grid
!      print*, 'Last point integrated:', 5+4*(num_quads-1)

!--------------------------------------------------------------
!                  open(unit=100,file='diagnostic3.out',status='old')
!                  do i=1,n
!                              write(100,1504) x(i), y(i), x_scaled(i), y_scaled(i),
!     *                         x_even(i), y_even(i), integral(1), h_step(1)
!                  end do
!1504                  format(1x,8e16.8)
!                  close(100)
!----------------------------------------------------------------

end subroutine boole

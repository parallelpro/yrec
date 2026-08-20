!----------------------------------------------------------------------
! polint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original polint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Neville's algorithm polynomial interpolation/extrapolation (as in
! Numerical Recipes POLINT), returning the interpolated value y at x
! and an error estimate dy.
! MHP 10/02: dimensions changed for consistency with fpft (fixed at
! 20 rather than sized to n).
subroutine polint(xa, ya, n, x, y, dy)
      implicit none

      double precision, intent(in) :: xa(20), ya(20)
      integer, intent(in) :: n
      double precision, intent(in) :: x
      double precision, intent(out) :: y, dy

      double precision :: c(20), d(20)
      integer :: ns, i, j
      double precision :: dif, dift, ho, hp, w, den
      save

      ns = 1
      dif = dabs(x-xa(1))
      do i = 1, n
       dift = dabs(x-xa(i))
       if (dift.lt.dif) then
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
          if (dabs(den).lt.1.0d-20) stop
          den = w/den
          d(i) = hp*den
          c(i) = ho*den
       end do
       if (2*ns.lt.n-j) then
          dy = c(ns+1)
       else
          dy = d(ns)
          ns = ns-1
       end if
       y = y + dy
      end do
      return
end subroutine polint

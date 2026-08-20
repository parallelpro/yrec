!----------------------------------------------------------------------
! osplin
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original osplin.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! ACM Algorithm 574: shape-preserving osculatory quadratic splines,
! by D.F. McAllister and J.A. Roulier, ACM Transactions on
! Mathematical Software, September 1981.
!
! xtab contains the abscissas of the points of interpolation.
! ytab contains the ordinates of the points of interpolation.
! n is the number of data points.
! k is the number of points at which the spline is to be evaluated
! (the points themselves are xval; the evaluated values come back in
! yval).
!
! Upon exit from subroutine 'slopes' -- first_derivs contains the
! computed first derivatives at each data point.
subroutine osplin(xval, yval, xtab, ytab, n, k)
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: xval(json)
      double precision, intent(out) :: yval(json)
      double precision, intent(in) :: xtab(json), ytab(json)
      integer, intent(in) :: n, k

      double precision :: first_derivs(json), eps
      integer :: err
      save

! calculate the slopes at each data point.
      call slopes(xtab, ytab, first_derivs, n)

! set the error tolerance eps, which is used in subroutine 'choose'.
      eps = 1.d-04
! call meval to evaluate the spline at the run of points xval.
      call meval(xval, yval, xtab, ytab, first_derivs, n, k, eps, err)

      return
end subroutine osplin

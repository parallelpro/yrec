!----------------------------------------------------------------------
! intpol
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original intpol.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! ***** the interpolation routine *****
! this routine contains two interpolation methods; hermite,
! and spline.  each of these has own merit and demerit.
! if the interpolant is smooth, both of these will give good
! results.  generally, spline gives more smooth interpolation.
! when the interpolant contains abrupt variation in gradient,
! however, spline get worse at that part, while hermit still
! gives reasonable result.  unfortunately, there is no criterion
! for selection between these two methods.  therefore, this
! routine gives the right for selection to user.
! YCK 3/91
!
! x_grid  ; array of abscissa points
! y_grid  ; array of ordinate points
! n_grid  ; size of the arrays
! x_eval  ; a x-value at which we want to find y-value
! k_lo    ; the grid point smaller than and closest to x_eval
! y_eval  ; the value we want
! dy_eval ; the derivative value at x_eval
subroutine intpol(x_grid, y_grid, n_grid, x_eval, y_eval, dy_eval)
      use luout_lib
      use numerics_lib
      implicit none
      integer, parameter :: np=100

      double precision, intent(in) :: x_grid(n_grid), y_grid(n_grid)
      integer, intent(in) :: n_grid
      double precision, intent(in) :: x_eval
      double precision, intent(out) :: y_eval, dy_eval

      double precision :: spline_coeff(4,np), dx_local, x_eval_copy
      double precision :: interp_value, interp_deriv
      integer :: i, k_lo, k_hi, k_mid
      data spline_coeff/400*0.0d0/
      save

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
         go to 522
      endif
      if(x_grid(n_grid).lt.x_eval)then
         k_lo=n_grid-1
         k_hi=n_grid
         go to 522
      endif
      k_lo=1
      k_hi=n_grid
    2 if((k_hi-k_lo).gt.1)then
         k_mid = (k_hi+k_lo)/2
         if(x_grid(k_mid).gt.x_eval_copy)then
            k_hi=k_mid
         else
            k_lo=k_mid
         endif
         go to 2
      endif
      if((k_hi-k_lo).le.0)then
         write(iowr, *) 'ERROR COX OP: INTERPOLATION'
         write(short_file_unit, *) 'ERROR COX OP: INTERPOLATION'
         stop
      endif
  522 continue
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

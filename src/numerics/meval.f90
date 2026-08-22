!----------------------------------------------------------------------
! meval
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original meval.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Controls the evaluation of an osculatory (shape-preserving)
! quadratic spline, by D.F. McAllister & J.A. Roulier, coded by
! S.L. Dodd & M. Roulier, N.C. State University. The caller may
! provide slopes at the points of interpolation directly, or use
! subroutine SLOPES to compute slopes consistent with the shape of
! the data.
!
! On input --
!   eval_x must be a nondecreasing vector of points at which the
!   spline will be evaluated.
!   table_x contains the abscissas of the data points to be
!   interpolated. table_x must be increasing.
!   table_y contains the ordinates of the data points to be
!   interpolated.
!   table_slope contains the slope of the spline at each point of
!   interpolation.
!   num_table_points is the number of data points (dimension of
!   table_x and table_y).
!   num_eval_points is the number of points of evaluation (dimension
!   of eval_x and eval_y).
!   eps_tol is a relative error tolerance used in subroutine CHOOSE to
!   distinguish the situation table_slope(i) or table_slope(i+1) is
!   relatively close to the slope or twice the slope of the linear
!   segment between table_x(i) and table_x(i+1). If this situation
!   occurs, roundoff may cause a change in convexity or monotonicity
!   of the resulting spline and a change in the case number provided
!   by CHOOSE. If eps_tol is not equal to zero, then eps_tol should be
!   greater than or equal to machine epsilon.
!
! On output --
!   eval_y contains the images of the points in eval_x.
!   err_code is an error code --
!   err_code=0 - meval ran normally.
!   err_code=1 - eval_x(i) is less than table_x(1) for at least one i,
!                or eval_x(i) is greater than table_x(num_table_points)
!                for at least one i. meval will extrapolate to provide
!                function values for these abscissas.
!   err_code=2 - eval_x(i+1) .lt. eval_x(i) for some i.
!
! meval does not alter eval_x, table_x, table_y, table_slope,
! num_table_points, num_eval_points.
!
! meval calls the following subroutines or functions:
!    search
!    cases
!    choose
!    spline
!----------------------------------------------------------------------
subroutine meval(eval_x, eval_y, table_x, table_y, table_slope, &
     num_table_points, num_eval_points, eps_tol, err_code)

      implicit none

      integer, parameter :: max_points = 5000

      double precision, intent(in) :: eval_x(max_points), &
           table_x(max_points), table_y(max_points), &
           table_slope(max_points), eps_tol
      double precision, intent(out) :: eval_y(max_points)
      integer, intent(in) :: num_table_points, num_eval_points
      integer, intent(out) :: err_code

      double precision, external :: spline
      double precision :: spline_v1, spline_v2, spline_w1, spline_w2, &
           spline_z1, spline_z2, spline_y1, spline_y2, spline_e1, spline_e2
      integer :: start_idx, start_idx1, end_idx, end_idx1, found_flag
      integer :: i, ind, loop_bound, table_idx, table_idx1, spline_case, &
           num_minus1
      start_idx = 1
      end_idx = num_eval_points
      err_code = 0
      if (.not. (num_eval_points .eq. 1)) then

! Determine if eval_x is nondecreasing.
      loop_bound = num_eval_points - 1
      do i = 1, loop_bound
        if (eval_x(i+1) .ge. eval_x(i)) cycle
        err_code = 2
        return
  10  continue
      end do

! If eval_x(i) .lt. table_x(1), then eval_x(i)=table_y(1).
! If eval_x(i) .gt. table_x(num_table_points), then
! eval_x(i)=table_y(num_table_points).
!
! Determine if any of the points in eval_x are less than the abscissa
! of the first data point.
      end if
  20 do i = 1, num_eval_points
        if (eval_x(i) .ge. table_x(1)) exit
        start_idx = i + 1
  30  continue
  end do

  40  loop_bound = num_eval_points + 1

! Determine if any of the points in eval_x are greater than the
! abscissa of the last data point.
      do i = 1, num_eval_points
        ind = loop_bound - i
        if (eval_x(ind) .le. table_x(num_table_points)) exit
        end_idx = ind - 1
  50  continue
      end do

! Calculate the images of points of evaluation whose abscissas are
! less than the abscissa of the first data point.
  60  if (start_idx .eq. 1) go to 80
! Set the error parameter to indicate that extrapolation has occurred.
      err_code = 1
      call choose(table_x(1), table_y(1), table_slope(1), table_slope(2), &
           table_x(2), table_y(2), eps_tol, spline_case)
      call cases(table_x(1), table_y(1), table_slope(1), table_slope(2), &
           table_x(2), table_y(2), spline_e1, spline_e2, spline_v1, &
           spline_v2, spline_w1, spline_w2, spline_z1, spline_z2, &
           spline_y1, spline_y2, spline_case)
      start_idx1 = start_idx - 1
      do i = 1, start_idx1
       eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, table_x(1), &
            table_y(1), table_x(2), table_y(2), spline_y1, spline_y2, &
            spline_e2, spline_w2, spline_v2, spline_case)
  70  continue
      end do
      if (num_eval_points .eq. 1) then
         return
      end if

! search locates the interval in which the first in-range point of
! evaluation lies.
  80  if ((num_eval_points .eq. 1) .and. (end_idx .ne. num_eval_points)) &
           go to 200
      call search(table_x, num_table_points, eval_x(start_idx), &
           table_idx, found_flag)

      table_idx1 = table_idx + 1

! If the first in-range point of evaluation is equal to one of the
! data points, assign the appropriate value from table_y. Continue
! until a point of evaluation is found which is not equal to a data
! point.
      if (found_flag .eq. 0) go to 130
      do
         eval_y(start_idx) = table_y(table_idx)
      start_idx1 = start_idx
      start_idx = start_idx + 1
      if (start_idx .gt. num_eval_points) then
         return
      end if
      if (.not. (eval_x(start_idx1) .eq. eval_x(start_idx))) exit
      end do

! KC 2025-05-30 fixed "Arithmetic IF statement"
! 100 IF (XVAL(START) - XTAB(LCN1)) 130,110,120
 100  if (eval_x(start_idx) .lt. table_x(table_idx1)) then
         goto 130
      else if (eval_x(start_idx) .eq. table_x(table_idx1)) then
         goto 110
      else
         goto 120
      end if

 110  eval_y(start_idx) = table_y(table_idx1)
      start_idx1 = start_idx
      start_idx = start_idx + 1
      if (start_idx .gt. num_eval_points) then
         return
      end if
      if (eval_x(start_idx) .ne. eval_x(start_idx1)) go to 120
      go to 110

 120  table_idx = table_idx1
      table_idx1 = table_idx1 + 1
      go to 100

! Calculate the images of all the points which lie within range of
! the data.
 130  if ((table_idx .eq. 1) .and. (err_code .eq. 1)) go to 140
      call choose(table_x(table_idx), table_y(table_idx), &
           table_slope(table_idx), table_slope(table_idx1), &
           table_x(table_idx1), table_y(table_idx1), eps_tol, spline_case)
      call cases(table_x(table_idx), table_y(table_idx), &
           table_slope(table_idx), table_slope(table_idx1), &
           table_x(table_idx1), table_y(table_idx1), spline_e1, spline_e2, &
           spline_v1, spline_v2, spline_w1, spline_w2, spline_z1, &
           spline_z2, spline_y1, spline_y2, spline_case)

 140 do i = start_idx, end_idx

! If eval_x(i) - table_x(table_idx1) is negative, do not recalculate
! the parameters for this section of the spline since they are
! already known.
! KC 2025-05-30 fixed "Arithmetic IF statement"
!       IF (XVAL(I) - XTAB(LCN1)) 150,160,170
      if (eval_x(i) .lt. table_x(table_idx1)) then
         goto 150
      else if (eval_x(i) .eq. table_x(table_idx1)) then
         goto 160
      else
         goto 170
      end if

 150  eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, &
           table_x(table_idx), table_y(table_idx), table_x(table_idx1), &
           table_y(table_idx1), spline_y1, spline_y2, spline_e2, &
           spline_w2, spline_v2, spline_case)
      cycle

! If eval_x(i) is a data point, its image is known.
 160  eval_y(i) = table_y(table_idx1)
      cycle

! Increment the pointers which give the location in the data vector.
 170  table_idx = table_idx1
      table_idx1 = table_idx + 1

! Determine that the routine is in the correct part of the spline.
! KC 2025-05-30 fixed "Arithmetic IF statement"
!       IF (XVAL(I) - XTAB(LCN1)) 180,160,170
      if (eval_x(i) .lt. table_x(table_idx1)) then
         goto 180
      else if (eval_x(i) .eq. table_x(table_idx1)) then
         goto 160
      else
         goto 170
      end if

! Call choose to determine the appropriate case and then call cases
! to compute the parameters of the spline.
 180  call choose(table_x(table_idx), table_y(table_idx), &
           table_slope(table_idx), table_slope(table_idx1), &
           table_x(table_idx1), table_y(table_idx1), eps_tol, spline_case)
      call cases(table_x(table_idx), table_y(table_idx), &
           table_slope(table_idx), table_slope(table_idx1), &
           table_x(table_idx1), table_y(table_idx1), spline_e1, spline_e2, &
           spline_v1, spline_v2, spline_w1, spline_w2, spline_z1, &
           spline_z2, spline_y1, spline_y2, spline_case)
      go to 150
 190  continue
 end do

! Calculate the images of the points of evaluation whose abscissas
! are greater than the abscissa of the last data point.
      if (end_idx .eq. num_eval_points) then
         return
      end if
      if ((table_idx1 .eq. num_table_points) .and. &
           (eval_x(end_idx) .ne. table_x(num_table_points))) go to 210

! Previously, when we arrived at 200 or 210, NUM1 could be improperly
! set. The NUM1= lines below protect from that. llp 8/19/08

! Set the error parameter to indicate that extrapolation has occurred.
 200  err_code = 1
      num_minus1 = max(num_table_points - 1, 1)
      call choose(table_x(num_minus1), table_y(num_minus1), &
           table_slope(num_minus1), table_slope(num_table_points), &
           table_x(num_table_points), table_y(num_table_points), &
           eps_tol, spline_case)
      call cases(table_x(num_minus1), table_y(num_minus1), &
           table_slope(num_minus1), table_slope(num_table_points), &
           table_x(num_table_points), table_y(num_table_points), &
           spline_e1, spline_e2, spline_v1, spline_v2, spline_w1, &
           spline_w2, spline_z1, spline_z2, spline_y1, spline_y2, &
           spline_case)
 210  end_idx1 = end_idx + 1
      num_minus1 = max(num_table_points - 1, 1)
      do i = end_idx1, num_eval_points
       eval_y(i) = spline(eval_x(i), spline_z1, spline_z2, &
            table_x(num_minus1), table_y(num_minus1), &
            table_x(num_table_points), table_y(num_table_points), &
            spline_y1, spline_y2, spline_e2, spline_w2, spline_v2, &
            spline_case)
 220  continue
      end do

 230  return
end subroutine meval

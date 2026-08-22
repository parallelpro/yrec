!----------------------------------------------------------------------
! cases
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original cases.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
! Computes the knots and other parameters of the spline on the
! interval (x_left,x_right).
!
! ON INPUT--
!
!   (x_left,y_left) AND (x_right,y_right) ARE THE COORDINATES OF THE
!   POINTS OF INTERPOLATION.
!
!   slope_left IS THE SLOPE AT (x_left,y_left).
!
!   slope_right IS THE SLOPE AT (x_right,y_right)
!
!   spline_case CONTROLS THE NUMBER AND LOCATION OF THE KNOTS.
!
!
! ON OUTPUT--
!
!   (knot_v_x,knot_v_y),(knot_w_x,knot_w_y),(knot_z_x,knot_z_y), AND
!   (knot_e_x,knot_e_y) ARE THE COORDINATES OF THE KNOTS AND OTHER
!   PARAMETERS OF THE SPLINE ON (x_left,x_right). (knot_e_x,knot_e_y)
!   AND (knot_y_x,knot_y_y) ARE USED ONLY IF spline_case=4.
!
! AND
!
!   CASES DOES NOT ALTER x_left,y_left,slope_left,slope_right,x_right,
!   y_right.
!
!----------------------------------------------------------------------
subroutine cases(x_left, y_left, slope_left, slope_right, x_right, &
     y_right, knot_e_x, knot_e_y, knot_v_x, knot_v_y, knot_w_x, &
     knot_w_y, knot_z_x, knot_z_y, knot_y_x, knot_y_y, spline_case)

      double precision :: x_left, y_left, slope_left, slope_right, &
           x_right, y_right, knot_v_x, knot_v_y, knot_z_x, knot_z_y, &
           knot_w_x, knot_w_y, knot_e_x, knot_e_y, &
           mbar1, mbar2, mbar3, c1, d1, h1, j1, knot_y_x, knot_y_y, &
           k1, ztwo
      integer :: spline_case
      if (.not. ((spline_case .eq. 3).or.(spline_case .eq. 4))) then
      if (.not. (spline_case .eq. 2)) then

! CALCULATE THE PARAMETERS FOR CASE 1.
      knot_z_x=(y_left-y_right+slope_right*x_right-slope_left*x_left)/ &
           (slope_right-slope_left)
      ztwo=y_left+slope_left*(knot_z_x-x_left)
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=(y_left+ztwo)/2.d0
      knot_w_x=(knot_z_x+x_right)/2.d0
      knot_w_y=(ztwo+y_right)/2.d0
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
           (knot_z_x-knot_v_x)
      return

! CALCULATE THE PARAMETERS FOR CASE 2.
      end if
  10  knot_z_x=(x_left+x_right)/2.d0
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(knot_z_x+x_right)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=(knot_v_y+knot_w_y)/2.d0
      return

! CALCULATE THE PARAMETERS USED IN BOTH CASES 3 AND 4.
      end if
  20  c1=x_left+(y_right-y_left)/slope_left
      d1=x_right+(y_left-y_right)/slope_right
      h1=2.d0*c1-x_left
      j1=2.d0*d1-x_right
      mbar1=(y_right-y_left)/(h1-x_left)
      mbar2=(y_left-y_right)/(j1-x_right)

      if (.not. (spline_case .eq. 4)) then

! CALCULATE THE PARAMETERS FOR CASE 3.
      k1=(y_left-y_right+x_right*mbar2-x_left*mbar1)/(mbar2-mbar1)
      if (.not. (abs(slope_left) .gt. abs(slope_right))) then
      knot_z_x=(k1+x_right)/2.d0
      knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
      (knot_z_x-knot_v_x)
      return
      end if
  30  knot_z_x=(k1+x_left)/2.d0
  40  knot_v_x=(x_left+knot_z_x)/2.d0
      knot_v_y=y_left+slope_left*(knot_v_x-x_left)
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=y_right+slope_right*(knot_w_x-x_right)
      knot_z_y=knot_v_y+((knot_w_y-knot_v_y)/(knot_w_x-knot_v_x))* &
           (knot_z_x-knot_v_x)
      return

! CALCULATE THE PARAMETERS FOR CASE 4.
      end if
  50  knot_y_x=(x_left+c1)/2.d0
      knot_v_x=(x_left+knot_y_x)/2.d0
      knot_v_y=slope_left*(knot_v_x-x_left) + y_left
      knot_z_x=(d1+x_right)/2.d0
      knot_w_x=(x_right+knot_z_x)/2.d0
      knot_w_y=slope_right*(knot_w_x-x_right) + y_right
      mbar3=(knot_w_y-knot_v_y)/(knot_w_x-knot_v_x)
      knot_y_y=mbar3*(knot_y_x-knot_v_x) + knot_v_y
      knot_z_y=mbar3*(knot_z_x-knot_v_x) + knot_v_y
      knot_e_x=(knot_y_x+knot_z_x)/2.d0
      knot_e_y=mbar3*(knot_e_x-knot_v_x) + knot_v_y
      return

end subroutine cases

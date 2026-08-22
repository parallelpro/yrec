!----------------------------------------------------------------------
! spline
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original spline.f; only variable names, source form, and comment
! style were updated.
!
!                                 SHAPE PRESERVING QUDRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A.ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!                                       N.C. STATE UNIVERSITY
!
! SPLINE FINDS THE IMAGE OF A POINT IN eval_point.
!
! ON INPUT--
!
!   eval_point CONTAINS THE VALUE AT WHICH THE SPLINE IS EVALUATED.
!
!   (x_left,y_left) ARE THE COORDINATES OF THE LEFT-HAND DATA POINT
!   USED IN THE EVALUATION OF eval_point.
!
!   (x_right,y_right) ARE THE COORDINATES OF THE RIGHT-HAND DATA POINT
!   USED IN THE EVALUATION OF eval_point.
!
!   z1,z2,y1,y2,e2,w2,v2 ARE THE PARAMETERS OF THE SPLINE.
!
!   spline_case CONTROLS THE EVALUATION OF THE SPLINE BY INDICATING
!   WHETHER ONE OR TWO KNOTS WERE PLACED IN THE INTERVAL
!   (x_left,x_right).
!
!
! ON OUTPUT--
!
!   SPLINE IS THE IMAGE OF eval_point.
!
! AND
!
!   SPLINE DOES NOT ALTER ANY OF THE INPUT PARAMETERS.
!
!----------------------------------------------------------------------
!
!  *****MODIFICATION DUE TO MARC PINSONNEAULT 6/87*****
!  IF DIVISION BY ZERO WOULD BE CAUSED,LINEAR INTERPOLATION IS USED
!  INSTEAD OF THE SPLINE.
! IF NCASE .EQ. 4, MORE THAN ONE KNOT WAS PLACED IN THE INTERVAL.
function spline(eval_point, z1, z2, x_left, y_left, x_right, y_right, &
     y1, y2, e2, w2, v2, spline_case)

      implicit none
      double precision :: spline
      double precision, intent(in) :: eval_point, z1, z2, x_left, y_left, &
           x_right, y_right, y1, y2, e2, w2, v2
      integer, intent(in) :: spline_case
      double precision :: linear_interp_frac

      if (.not. (spline_case .eq. 4)) then
!
! CASES 1,2, OR 3.
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE KNOT.
! KC 2025-05-30 fixed "Arithmetic IF statement"
!       IF (Z1 - XVALS) 10,20,30
      if (z1 .lt. eval_point) then
         goto 10
      else if (z1 .eq. eval_point) then
         spline=z2
         return
      else
         goto 30
      end if
!
  10  continue
      if(x_right.ne.z1)then
      spline=(z2*(x_right-eval_point)**2+w2*2.d0*(eval_point-z1)*(x_right-eval_point) &
              +y_right*(eval_point-z1)**2)/(x_right-z1)**2
      else
       linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
       spline = y_left + linear_interp_frac*(y_right - y_left)
      end if
      return
!
  20  spline=z2
      return
!
  30  continue
      if(z1.ne.x_left)then
      spline=(y_left*(z1-eval_point)**2+v2*2.d0*(eval_point-x_left)*(z1-eval_point)+ &
              z2*(eval_point-x_left)**2)/(z1-x_left)**2
      else
       linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
       spline = y_left + linear_interp_frac*(y_right - y_left)
      end if
      return
!
! CASE 4.
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE FIRST KNOT.
! KC 2025-05-30 fixed "Arithmetic IF statement"
!   40  IF (Y1 - XVALS) 70,60,50
      end if
  40  if (y1 .lt. eval_point) then
         goto 70
      else if (y1 .eq. eval_point) then
         spline=y2
         return
      else
         goto 50
      end if
!
  50  continue
      if(y1.ne.x_left)then
      spline=(y_left*(y1-eval_point)**2+v2*2.d0*(eval_point-x_left)*(y1-eval_point)+ &
              y2*(eval_point-x_left)**2)/(y1-x_left)**2
      else
       linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
       spline = y_left + linear_interp_frac*(y_right - y_left)
      end if
      return
!
  60  spline=y2
      return
!
! DETERMINE THE LOCATION OF XVALS RELATIVE TO THE SECOND KNOT.
! KC 2025-05-30 fixed "Arithmetic IF statement"
!   70  IF (Z1 - XVALS) 100,90,80
  70  if (z1 .lt. eval_point) then
         goto 100
      else if (z1 .eq. eval_point) then
         spline=z2
         return
      else
         goto 80
      end if
!
  80  continue
      if(z1.ne.y1)then
      spline=(y2*(z1-eval_point)**2+e2*2.d0*(eval_point-y1)*(z1-eval_point)+z2*(eval_point &
              -y1)**2)/(z1-y1)**2
      else
       linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
       spline = y_left + linear_interp_frac*(y_right - y_left)
      end if
      return
!
  90  spline=z2
      return
!
  100 continue
      if(x_right.ne.z1)then
      spline=(z2*(x_right-eval_point)**2+w2*2.d0*(eval_point-z1)*(x_right-eval_point) &
              +y_right*(eval_point-z1)**2)/(x_right-z1)**2
      else
       linear_interp_frac = (eval_point - x_left)/(x_right - x_left)
       spline = y_left + linear_interp_frac*(y_right - y_left)
      end if
      return
end function spline

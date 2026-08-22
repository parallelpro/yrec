!----------------------------------------------------------------------
! choose
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original choose.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!                                 SHAPE PRESERVING QUADRATIC SPLINES
!                                   BY D.F.MCALLISTER & J.A. ROULIER
!                                     CODED BY S.L.DODD & M.ROULIER
!
! CHOOSE DETERMINES THE CASE NEEDED FOR THE COMPUTATION OF THE PARAME-
! TERS OF THE QUADRATIC SPLINE AND RETURNS THE VALUE IN THE VARIABLE
! spline_case.
!
! ON INPUT--
!
!   (x_left,y_left) GIVES THE COORDINATES OF ONE OF THE POINTS OF
!   INTERPOLATION.
!
!   slope_left SPECIFIES THE DERIVATIVE CONDITION AT (x_left,y_left).
!
!   (x_right,y_right) GIVES THE COORDINATES OF ONE OF THE POINTS OF
!   INTERPOLATION.
!
!   slope_right SPECIFIES THE DERIVATIVE CONDITION AT (x_right,y_right).
!
!   eps_tol IS AN ERROR TOLERANCE USED TO DISTINGUISH CASES WHEN
!   slope_left OR slope_right IS RELATIVELY CLOSE TO THE SLOPE OR TWICE
!   THE SLOPE OF THE LINE SEGMENT JOINING (x_left,y_left) AND
!   (x_right,y_right). IF eps_tol IS NOT EQUAL TO ZERO, THEN eps_tol
!   SHOULD BE GREATER THAN OR EQUAL TO MACHINE EPSILON.
!
!
! ON OUTPUT--
!
!   spline_case CONTAINS THE VALUE WHICH CONTROLS HOW THE PARAMETERS OF
!   THE QUADRATIC SPLINE ARE EVALUATED.
!
! AND
!
!   CHOOSE DOES NOT ALTER x_left,y_left,x_right,y_right,slope_left,
!   slope_right,eps_tol.
!
!----------------------------------------------------------------------
subroutine choose(x_left, y_left, slope_left, slope_right, x_right, &
     y_right, eps_tol, spline_case)

      double precision :: x_left, y_left, slope_left, slope_right, &
           x_right, y_right, mref, mref1, mref2, spq, prod, &
           prod1, prod2, eps_tol
      integer :: spline_case

! CALCULATE THE SLOPE spq OF THE LINE JOINING (x_left,y_left),(x_right,y_right).
      spq=(y_right-y_left)/(x_right-x_left)

! CHECK WHETHER OR NOT spq IS 0.
! ******MODIFICATION BY MARC PINSONNEAULT TO AVOID DIVISION BY ZERO
! ******IN SR CASES
      if (.not. (slope_left.eq.0.0d0.or.slope_right.eq.0.0d0)) then
! ******
      if (spq .ne. 0.d0) go to 20
      if ((slope_left*slope_right) .ge. 0.d0) then
         spline_case=2
         return
      end if
      spline_case=1
      return
      end if
   9  continue
  10  spline_case=2
      return

  20  prod1=spq*slope_left
      prod2=spq*slope_right

! FIND THE ABSOLUTE VALUES OF THE SLOPES spq,slope_left,AND slope_right.
      mref=abs(spq)
      mref1=abs(slope_left)
      mref2=abs(slope_right)

! IF THE RELATIVE DEVIATION OF slope_left OR slope_right FROM spq IS LESS THAN
! eps_tol, THEN CHOOSE CASE 2 OR CASE 3.
      if ((abs(spq-slope_left).le.eps_tol*mref) .or. &
           (abs(spq-slope_right).le.eps_tol*mref)) &
          go to 30

! COMPARE THE SIGNS OF THE SLOPES spq,slope_left, AND slope_right.
      if ((prod1 .lt. 0.d0).or.(prod2 .lt. 0.d0)) go to 80
      prod=(mref-mref1)*(mref-mref2)
      if (prod .ge. 0.d0) go to 40

! L1, THE LINE THROUGH (x_left,y_left) WITH SLOPE slope_left, AND L2, THE LINE
! THROUGH (x_right,y_right) WITH SLOPE slope_right, INTERSECT AT A POINT WHOSE
! ABSCISSA IS BETWEEN x_left AND x_right. THE ABSCISSA BECOMES A KNOT OF THE
! SPLINE.
      spline_case=1
      return

  30  if ((prod1 .lt. 0.d0).or.(prod2 .lt. 0.d0)) go to 80
  40  if (mref1 .gt. (2.d0*mref)) go to 50
      if (mref2 .gt. (2.d0*mref)) then
         if (mref1 .gt. (2.d0-eps_tol)*mref) go to 70
         spline_case=3
         return
      end if

! BOTH L1 AND L2 CROSS THE LINE THROUGH (x_left+x_right/2.,y_left) AND
! (x_left+x_right/2.,y_right), WHICH IS THE MIDLINE OF THE RECTANGLE FORMED
! BY (x_left,y_left),(x_right,y_left),(x_right,y_right), AND (x_left,y_right),
! OR BOTH slope_left AND slope_right HAVE SIGNS DIFFERENT THAN THE SIGN OF
! spq, OR ONE OF slope_left AND slope_right HAS OPPOSITE SIGN FROM spq AND L1
! AND L2 INTERSECT TO THE LEFT OF x_left OR TO THE RIGHT OF x_right. THE
! POINT (x_left+x_right)/2. IS A KNOT OF THE SPLINE.
      spline_case=2
      return

! CHOOSE CASE 4 IF mref2 IS GREATER THAN (2.-eps_tol)*mref; OTHERWISE,
! CHOOSE CASE 3.
  50  if (mref2 .gt. (2.d0-eps_tol)*mref) go to 70
      spline_case=3
      return

! IN CASES 3 AND 4, SIGN(slope_left)=SIGN(slope_right)=SIGN(spq).
!
! EITHER L1 OR L2 CROSSES THE MIDLINE, BUT NOT BOTH.
! CHOOSE CASE 4 IF mref1 IS GREATER THAN (2.-eps_tol)*mref; OTHERWISE,
! CHOOSE CASE 3.
  60  if (mref1 .gt. (2.d0-eps_tol)*mref) go to 70
      spline_case=3
      return

! IF NEITHER L1 NOR L2 CROSSES THE MIDLINE, THE SPLINE REQUIRES TWO
! KNOTS BETWEEN x_left AND x_right.
  70  spline_case=4
      return

! THE SIGN OF AT LEAST ONE OF THE SLOPES slope_left,slope_right DOES NOT
! AGREE WITH THE SIGN OF THE SLOPE spq.
  80  if ((prod1 .lt. 0.d0).and.(prod2 .lt. 0.d0)) go to 130

      if (.not. (prod1 .lt. 0.d0)) then
      if (mref1 .gt. ((1.d0+eps_tol)*mref)) go to 120
      spline_case=2
      return

      end if
  90  if (mref2 .gt. ((1.d0+eps_tol)*mref)) go to 100
      spline_case=2
      return

  100 spline_case=1
      return

  110 if (mref1 .gt. ((1.d0+eps_tol)*mref)) go to 120
      spline_case=2
      return

  120 spline_case=1
      return

  130 spline_case=2
      return

end subroutine choose

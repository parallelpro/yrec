!----------------------------------------------------------------------
! intrp2
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original intrp2.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! interp without derivatives. Interp is the interpolation routine for
! the VandenBerg opacity tables (CAPPA), and it uses a 4-point
! Lagrangian interpolation scheme. Given the 4 abscissa points
! x_nodes and evaluation point x_eval, returns the Lagrangian weights
! (weight) such that a function's interpolated value at x_eval is
! sum(weight(:)*f(:)) for f sampled at x_nodes.
subroutine intrp2(x_nodes, weight, x_eval)

      implicit none
      double precision, intent(in) :: x_nodes(4)
      double precision, intent(out) :: weight(4)
      double precision, intent(in) :: x_eval

      double precision :: diff43, diff42, diff41, diff32, diff31, diff21
      double precision :: denom1, denom2, denom3, denom4
      double precision :: dx1, dx2, dx3, dx4
      save

! interp is the interpolation routine for the VandenBerg
! opacity tables(CAPPA), and it uses a 4-point Lagrangian
! interpolation scheme.
      diff43 = x_nodes(4) - x_nodes(3)
      diff42 = x_nodes(4) - x_nodes(2)
      diff41 = x_nodes(4) - x_nodes(1)
      diff32 = x_nodes(3) - x_nodes(2)
      diff31 = x_nodes(3) - x_nodes(1)
      diff21 = x_nodes(2) - x_nodes(1)
      denom1 = -diff21*diff31*diff41
      denom2 = diff21*diff32*diff42
      denom3 = -diff31*diff32*diff43
      denom4 = diff41*diff42*diff43
      dx1 = x_eval - x_nodes(1)
      dx2 = x_eval - x_nodes(2)
      dx3 = x_eval - x_nodes(3)
      dx4 = x_eval - x_nodes(4)
      weight(1) = (dx2*dx3*dx4)/denom1
      weight(2) = (dx1*dx3*dx4)/denom2
      weight(3) = (dx1*dx2*dx4)/denom3
      weight(4) = (dx1*dx2*dx3)/denom4
      return
end subroutine intrp2

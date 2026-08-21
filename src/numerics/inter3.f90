!----------------------------------------------------------------------
! inter3
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original inter3.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! inter3 is the interpolation routine for density in the Livermore
! opacity tables, and it uses a 3-point Lagrangian interpolation
! scheme. Given the 3 abscissa points x_nodes and evaluation point
! x_eval, returns the Lagrangian weights (weight) and their
! derivatives with respect to x_eval (dweight) such that a function's
! interpolated value/derivative at x_eval are sum(weight(:)*f(:)) and
! sum(dweight(:)*f(:)) for f sampled at x_nodes.
module inter3_mod
contains
subroutine inter3(x_nodes, weight, dweight, x_eval)

      implicit none
      double precision, intent(in) :: x_nodes(3)
      double precision, intent(out) :: weight(3), dweight(3)
      double precision, intent(in) :: x_eval

      double precision :: diff32, diff31, diff21
      double precision :: denom1, denom2, denom3
      double precision :: dx1, dx2, dx3
      save

! inter3 is the interpolation routine for density in the livermore
! opacity tables, and it uses a 3-point lagrangian interpolation scheme.
      diff32 = x_nodes(3) - x_nodes(2)
      diff31 = x_nodes(3) - x_nodes(1)
      diff21 = x_nodes(2) - x_nodes(1)
      denom1 = diff21*diff31
      denom2 = -diff21*diff32
      denom3 = diff31*diff32
      dx1 = x_eval - x_nodes(1)
      dx2 = x_eval - x_nodes(2)
      dx3 = x_eval - x_nodes(3)
      weight(1) = (dx2*dx3)/denom1
      weight(2) = (dx1*dx3)/denom2
      weight(3) = (dx1*dx2)/denom3
      dweight(1) = (dx2 + dx3)/denom1
      dweight(2) = (dx1 + dx3)/denom2
      dweight(3) = (dx1 + dx2)/denom3

      return
end subroutine inter3
end module inter3_mod

!----------------------------------------------------------------------
! numerics_lib
!----------------------------------------------------------------------
! The numerics domain's one public module (docs/GUIDELINES.md naming rule
! 2): it re-exports the generic numerical routines, which live in five
! modules by topic:
!   cubic_spline         cspline, natural_spline, splint, ysplin, intpol
!   quadratic_spline     osplin (shape-preserving, ACM 574)
!   table_interpolation  inter3, interp, intrp2, lagrange4, quint,
!                        polint, lir, findex, locate
!   linear_solvers       ludcmp, lubksb, simeqc, solve_tridiagonal
!   bulirsch_stoer       bsstep, for an ode_system
! None of them reads the star, keeps state between calls or prints.
! They report failure through ierr, as numerics_status's codes, and the
! caller reports it. numerics_termination is kept for now: the envelope
! integration's status is envelope_termination
! (core/envelope/envelope_integration_kernel).
module numerics_lib
      use cubic_spline, only: cspline, natural_spline, splint, &
           compute_spline_coefficients, ysplin, intpol
      use quadratic_spline, only: osplin
      use table_interpolation, only: inter3, interp, intrp2, lagrange4, quint, &
           polint, lir, findex, locate
      use linear_solvers, only: ludcmp, lubksb, simeqc, solve_tridiagonal
      use bulirsch_stoer, only: ode_system, bsstep, bsstep_max_equations, &
           bsstep_max_stages, bsstep_max_extrap_order
      use numerics_status, only: numerics_ok, numerics_singular, &
           numerics_coincident_points, numerics_bad_input, &
           numerics_step_underflow, numerics_unsorted_points, &
           numerics_derivatives_failed, numerics_status_text
      implicit none
      private
      public :: cspline, natural_spline, splint, compute_spline_coefficients
      public :: ysplin, intpol
      public :: osplin
      public :: inter3, interp, intrp2, lagrange4, quint, polint, lir, findex, locate
      public :: ludcmp, lubksb, simeqc, solve_tridiagonal
      public :: ode_system, bsstep, bsstep_max_equations, bsstep_max_stages, &
           bsstep_max_extrap_order
      public :: numerics_ok, numerics_singular, numerics_coincident_points, &
           numerics_bad_input, numerics_step_underflow, numerics_unsorted_points, &
           numerics_derivatives_failed, numerics_status_text
! kept for now: use envelope_termination (envelope_integration_kernel)
      integer, parameter, public :: numerics_termination = -2
end module numerics_lib

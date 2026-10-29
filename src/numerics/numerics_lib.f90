!----------------------------------------------------------------------
! numerics_lib
!----------------------------------------------------------------------
! The numerics domain's one public module (docs/GUIDELINES.md naming rule
! 2): it re-exports the generic numerical routines, which live in five
! modules by topic, and the status codes (numerics/README.md):
!   cubic_spline         cspline, natural_spline, splint,
!                        compute_spline_coefficients
!   quadratic_spline     osplin (shape-preserving, ACM 574)
!   table_interpolation  inter3, interp, lagrange4, quint, polint, lir,
!                        locate, hunt, findex
!   linear_solvers       ludcmp, lubksb, simeqc, solve_tridiagonal,
!                        tridiagonal_forward_sweep,
!                        tridiagonal_back_substitution
!   bulirsch_stoer       bsstep for an ode_system, and its limits
!   numerics_status      the codes returned in ierr, numerics_status_text
! The YREC-written routines with short legacy names have a descriptive
! alias, the same procedure under a second name (new code may use
! either; the legacy names stay for their callers):
!   lagrange3_weights => inter3     lagrange4_weights => interp
!   equidistant_quadratic => quint  interpolate_tabulated_functions => lir
!   bracket_index => findex
! None of them reads the star, keeps state between calls or prints.
! They report failure through ierr, as numerics_status's codes, and the
! caller reports it.
module numerics_lib
      use cubic_spline, only: cspline, natural_spline, splint, &
           compute_spline_coefficients
      use quadratic_spline, only: osplin
      use table_interpolation, only: inter3, interp, lagrange4, quint, &
           polint, lir, findex, locate, hunt
      use table_interpolation, only: lagrange3_weights => inter3, &
           lagrange4_weights => interp, equidistant_quadratic => quint, &
           interpolate_tabulated_functions => lir, bracket_index => findex
      use linear_solvers, only: ludcmp, lubksb, simeqc, solve_tridiagonal, &
           tridiagonal_forward_sweep, tridiagonal_back_substitution
      use bulirsch_stoer, only: ode_system, bsstep, bsstep_max_equations, &
           bsstep_max_stages, bsstep_max_extrap_order
      use numerics_status, only: numerics_ok, numerics_singular, &
           numerics_coincident_points, numerics_bad_input, &
           numerics_step_underflow, numerics_unsorted_points, &
           numerics_derivatives_failed, numerics_status_text
      implicit none
      private
      public :: cspline, natural_spline, splint, compute_spline_coefficients
      public :: osplin
      public :: inter3, interp, lagrange4, quint, polint, lir, findex, locate, hunt
      public :: lagrange3_weights, lagrange4_weights, equidistant_quadratic, &
           interpolate_tabulated_functions, bracket_index
      public :: ludcmp, lubksb, simeqc, solve_tridiagonal
      public :: tridiagonal_forward_sweep, tridiagonal_back_substitution
      public :: ode_system, bsstep, bsstep_max_equations, bsstep_max_stages, &
           bsstep_max_extrap_order
      public :: numerics_ok, numerics_singular, numerics_coincident_points, &
           numerics_bad_input, numerics_step_underflow, numerics_unsorted_points, &
           numerics_derivatives_failed, numerics_status_text
end module numerics_lib

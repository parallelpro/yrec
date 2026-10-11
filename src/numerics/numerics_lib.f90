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
!   bulirsch_stoer       bsstep
! None of them reads the star or keeps state between calls. They report
! failure through ierr; numerics_termination is the status a driver
! gives a failed envelope integration (the historical "solution
! diverged, stop", which exited 0), so that main and yrec_capi can tell
! it from a configuration error.
module numerics_lib
      use cubic_spline, only: cspline, natural_spline, splint, ysplin, intpol
      use quadratic_spline, only: osplin
      use table_interpolation, only: inter3, interp, intrp2, lagrange4, quint, &
           polint, lir, findex, locate
      use linear_solvers, only: ludcmp, lubksb, simeqc, solve_tridiagonal
      use bulirsch_stoer, only: bsstep
      implicit none
      private
      public :: cspline, natural_spline, splint, ysplin, intpol
      public :: osplin
      public :: inter3, interp, intrp2, lagrange4, quint, polint, lir, findex, locate
      public :: ludcmp, lubksb, simeqc, solve_tridiagonal
      public :: bsstep
      integer, parameter, public :: numerics_termination = -2
end module numerics_lib

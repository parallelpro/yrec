!----------------------------------------------------------------------
! quint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original quint.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Quadratic interpolation for equidistant points.
! y0=y(x0), y1=y(x1), y2=y(x2); h=x1-x0=x2-x1; computes y=y(x).
subroutine quint(x, x0, h, y0, y1, y2, y)
      implicit none
      double precision, intent(in) :: x, x0, h, y0, y1, y2
      double precision, intent(out) :: y

      double precision :: d1, d2, t
      save

      d1 = y1 - y0
      d2 = y2 - 2.d0*y1 + y0
      t = (x - x0)/h
      y = y0 + t*d1 + 0.5d0*t*(t-1.d0)*d2
      return
end subroutine quint

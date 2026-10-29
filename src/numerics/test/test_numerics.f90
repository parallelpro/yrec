!----------------------------------------------------------------------
! test_numerics
!----------------------------------------------------------------------
! numerics_lib's routines, in two parts.
! The pins: every public routine on small polynomial and rational data,
! error paths included, printed for a byte comparison with
! expected_test_numerics.out (they guard refactors):
!   cubic splines     cspline/splint, natural_spline (four points and
!                     seven), ysplin/intpol (eight, three and two points)
!   quadratic spline  osplin (McAllister-Roulier) on monotone, oscillating,
!                     stepped and concave data, extrapolation, single
!                     points, nodes, runs outside the table, a repeated
!                     abscissa and its error codes
!   interpolation     inter3, interp, intrp2, lagrange4, quint, polint,
!                     lir (cubic, linear, on a node, extrapolated)
!   search            findex, locate
!   linear algebra    ludcmp/lubksb, simeqc, solve_tridiagonal
!   ODE               bsstep (with mmid and ratext) on a harmonic
!                     ode_system: its acceptance branches, a retry, an
!                     integrand failure, step underflow, NaN errors and
!                     its limits
! The checks: each prints its error and ok or FAIL against its
! tolerance (a contract check prints 0 or 1), and the program stops
! with status 1 on any FAIL; numerics/README.md section 6 lists them.
! Each was confirmed to fail under a mutation of the code it covers
! (mutate_numerics.py).
module test_numerics_callbacks
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      use numerics_lib, only: ode_system
      implicit none

! dy1/dx = y2, dy2/dx = -y1, dy3/dx = -decay*y3: a rotation and a decay.
! Counts its calls; fails past x = 0.75 when fail is set; returns a NaN
! derivative for component nan_component (0: none).
      type, extends(ode_system) :: harmonic_system
         real(dp) :: decay = 0.3d0
         integer :: calls = 0
         logical :: fail = .false.
         integer :: nan_component = 0
      contains
         procedure :: derivatives => harmonic_derivatives
      end type harmonic_system

contains

subroutine harmonic_derivatives(self, x, y, dydx, ierr)
      class(harmonic_system), intent(inout) :: self
      real(dp), intent(in) :: x, y(:)
      real(dp), intent(out) :: dydx(:)
      integer, intent(out) :: ierr
      real(dp), parameter :: fail_x = 0.75d0
      ierr = 0
      self%calls = self%calls + 1
      if (self%fail .and. x > fail_x) then
         ierr = 3
         return
      end if
      dydx(1) = y(2)
      dydx(2) = -y(1)
      dydx(3) = -self%decay*y(3)
      if (self%nan_component > 0) dydx(self%nan_component) = ieee_value(1.0d0, ieee_quiet_nan)
end subroutine harmonic_derivatives

end module test_numerics_callbacks

program test_numerics
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      use math_lib
      use numerics_lib, only: bsstep, cspline, findex, inter3, interp, intpol, &
           intrp2, lagrange4, lir, locate, &
           lubksb, ludcmp, natural_spline, osplin, &
           polint, quint, simeqc, solve_tridiagonal, &
           splint, ysplin, compute_spline_coefficients, numerics_ok, &
           numerics_singular, numerics_coincident_points, numerics_bad_input, numerics_step_underflow, &
           numerics_unsorted_points, numerics_derivatives_failed, numerics_status_text
      use numerics_lib, only: tridiagonal_forward_sweep, tridiagonal_back_substitution
      use numerics_lib, only: hunt
      use test_numerics_callbacks, only: harmonic_system
      implicit none
! the pins' tables are larger than the data they hold
      integer, parameter :: table_size = 100
      real(dp) :: jx(table_size), jy(table_size), jy2(table_size), jval(table_size), &
           jout(table_size)
      real(dp) :: ja(table_size), jb(table_size), jc(table_size), jr(table_size), &
           jsol(table_size)
      real(dp) :: x6(6), y6(6), y26(6), v, dv, d
      real(dp) :: x4(4), y4(4), y24(4), w4(4), dw4(4), w3(3), dw3(3)
      real(dp) :: x20(20), y20(20)
      real(dp) :: a5(5,5), b5(5), sm(15)
      integer :: indx(5), i, idx, klo, khi, ierr, flag, guesses(7), failures
      real(dp) :: pts(10), pts2(10)

      failures = 0

! ---------------- cubic splines
      write(*,'(a)') "# cspline/splint: y = x**3 - 2x + 1 on six uneven points"
      x6 = [0.0d0, 0.5d0, 1.3d0, 2.0d0, 3.1d0, 4.0d0]
      y6 = x6**3 - 2.0d0*x6 + 1.0d0
      call cspline(x6, y6, 6, 1.0d30, 1.0d30, y26)
      call show("natural y2", y26)
      call cspline(x6, y6, 6, -2.0d0, 46.0d0, y26)
      call show("clamped y2", y26)
      pts(1:5) = [-0.25d0, 0.5d0, 1.7d0, 3.95d0, 4.5d0]
      do i = 1, 5
         v = pts(i)
         call splint(x6, y6, 6, y26, v, d, ierr, klo, khi)
         write(*,'(a,f6.2,a,1pe24.15,3i3)') "splint at", v, " =", d, klo, khi, ierr
      end do
! bisection gives a zero-width interval only between equal last points
      x6(6) = x6(5)
      call splint(x6, y6, 6, y26, 5.0d0, d, ierr, klo, khi)
      write(*,'(a,3i3)') "splint, repeated abscissa: klo khi ierr", klo, khi, ierr
      x6(6) = 4.0d0

      write(*,'(a)') "# natural_spline/splint: four points"
      x4 = [1.0d0, 2.0d0, 4.0d0, 7.0d0]
      y4 = 1.0d0/x4
      call natural_spline(x4, y4, 4, y24)
      call show("natural_spline y2", y24)
      pts(1:3) = [1.5d0, 4.0d0, 6.5d0]
      do i = 1, 3
         v = pts(i)
         call splint(x4, y4, 4, y24, v, d, ierr)
         write(*,'(a,f6.2,a,1pe24.15,i3)') "splint at", v, " =", d, ierr
      end do
      x4(4) = x4(3)
      call splint(x4, y4, 4, y24, 8.0d0, d, ierr)
      write(*,'(a,i3)') "splint (4 points), repeated abscissa: ierr", ierr

      write(*,'(a)') "# natural_spline/splint: seven points in arrays of 100"
      jx = 0.0d0
      jy = 0.0d0
      do i = 1, 7
         jx(i) = 0.3d0*i + 0.01d0*i*i
         jy(i) = jx(i)**2/(1.0d0 + jx(i))
      end do
      call natural_spline(jx, jy, 7, jy2)
      call show("natural_spline y2", jy2(1:7))
      call splint(jx, jy, 7, jy2, 1.234d0, d, ierr, klo, khi)
      write(*,'(a,1pe24.15,3i3)') "splint at 1.234 =", d, klo, khi, ierr

      write(*,'(a)') "# ysplin/intpol: eight points"
      do i = 1, 8
         jx(i) = 0.5d0*i
         jy(i) = (jx(i) - 2.0d0)**3 + jx(i)
      end do
      pts(1:5) = [0.2d0, 0.5d0, 1.9d0, 3.3d0, 4.6d0]
      do i = 1, 5
         v = pts(i)
         call intpol(jx(1:8), jy(1:8), 8, v, d, dv, ierr)
         write(*,'(a,f6.2,a,2(1pe24.15),i3)') "intpol at", v, " =", d, dv, ierr
      end do
      call ysplin_check()
      write(*,'(a)') "# ysplin/intpol: three and two points (the natural spline's slopes)"
      jx(1:3) = [0.0d0, 1.0d0, 3.0d0]
      jy(1:3) = [0.0d0, 1.0d0, 0.0d0]
      call intpol(jx(1:3), jy(1:3), 3, 1.0d0, d, dv, ierr)
      write(*,'(a,2(1pe24.15),i3)') "intpol (0,0),(1,1),(3,0) at 1: value, slope 1/2", d, dv, ierr
      call intpol(jx(1:3), jy(1:3), 3, 2.0d0, d, dv, ierr)
      write(*,'(a,2(1pe24.15),i3)') "intpol (0,0),(1,1),(3,0) at 2", d, dv, ierr
      jx(1:2) = [0.0d0, 2.0d0]
      jy(1:2) = [1.0d0, 5.0d0]
      call intpol(jx(1:2), jy(1:2), 2, 0.5d0, d, dv, ierr)
      write(*,'(a,2(1pe24.15),i3)') "intpol (0,1),(2,5) at 0.5: value 2, slope 2", d, dv, ierr

! ---------------- quadratic spline
      write(*,'(a)') "# osplin: monotone convex data, points below, on, between and above"
      do i = 1, 8
         jx(i) = dble(i - 1)
         jy(i) = jx(i)**2
      end do
      call osplin_case(8, [-0.5d0, 0.0d0, 0.0d0, 0.4d0, 2.0d0, 2.5d0, 3.7d0, 6.9d0, 7.0d0, 7.6d0])
      write(*,'(a)') "# osplin: oscillating data (slopes of opposite sign, zero slopes)"
      jy(1:8) = [0.0d0, 1.0d0, 0.5d0, 2.0d0, -1.0d0, -1.0d0, 3.0d0, 2.5d0]
      call osplin_case(8, [0.3d0, 0.9d0, 1.5d0, 2.2d0, 3.6d0, 4.5d0, 5.2d0, 6.4d0, 6.8d0, 6.99d0])
      write(*,'(a)') "# osplin: steps (one and two knots per interval)"
      jy(1:8) = [0.0d0, 0.01d0, 0.02d0, 5.0d0, 5.01d0, 5.02d0, 10.0d0, 10.1d0]
      call osplin_case(8, [0.5d0, 1.5d0, 2.1d0, 2.5d0, 2.9d0, 3.5d0, 4.5d0, 5.5d0, 6.2d0, 6.8d0])
      write(*,'(a)') "# osplin: concave data with a near-linear stretch"
      jy(1:8) = [0.0d0, 3.0d0, 5.0d0, 6.0d0, 6.5d0, 7.0d0, 7.5d0, 7.6d0]
      call osplin_case(8, [0.25d0, 0.75d0, 1.25d0, 2.5d0, 3.5d0, 4.25d0, 4.75d0, 5.5d0, 6.5d0, 6.9d0])
      write(*,'(a)') "# osplin: single points (below, inside, on a node, above)"
      jy(1:8) = [0.0d0, 1.0d0, 4.0d0, 9.0d0, 16.0d0, 25.0d0, 36.0d0, 49.0d0]
      call osplin_case(8, [-1.0d0])
      call osplin_case(8, [2.6d0])
      call osplin_case(8, [3.0d0])
      call osplin_case(8, [8.5d0])
      write(*,'(a)') "# osplin: nodes after a node, the last node, the last interval's tail"
      call osplin_case(8, [2.0d0, 3.0d0, 3.0d0, 4.5d0])
      call osplin_case(8, [2.0d0, 3.0d0])
      call osplin_case(8, [7.0d0, 7.6d0])
      call osplin_case(8, [6.5d0, 7.5d0])
      write(*,'(a)') "# osplin: decreasing points (meval refuses, output untouched)"
      call osplin_case(8, [2.0d0, 1.0d0])
      write(*,'(a)') "# osplin: runs outside the table, a node then points above"
      call osplin_case(8, [-2.0d0, -1.0d0, -0.5d0])
      call osplin_case(8, [7.5d0, 8.0d0, 9.0d0])
      call osplin_case(8, [-1.0d0, 7.5d0, 8.0d0])
      call osplin_case(8, [7.0d0, 7.5d0, 8.0d0])
      call osplin_case(8, [6.0d0, 7.0d0, 7.0d0, 8.0d0])
      write(*,'(a)') "# osplin: a repeated abscissa inside the table (finite slopes around it)"
      jx(1:6) = [0.0d0, 1.0d0, 2.0d0, 2.0d0, 3.0d0, 4.0d0]
      jy(1:6) = [0.0d0, 1.0d0, 4.0d0, 4.0d0, 9.0d0, 16.0d0]
      call osplin_case(6, [0.5d0, 1.5d0, 2.5d0, 3.5d0])
      do i = 1, 8
         jx(i) = dble(i - 1)
      end do
      jy(1:8) = [0.0d0, 1.0d0, 4.0d0, 9.0d0, 16.0d0, 25.0d0, 36.0d0, 49.0d0]
      write(*,'(a)') "# osplin: error codes (n < 3; decreasing points; none; a zero-width end interval)"
      jout = -1.0d0
      call osplin(jval, jout, jx, jy, 2, 2, ierr)
      write(*,'(a,i3,2f6.1)') "osplin n = 2: ierr, yval untouched", ierr, jout(1:2)
      jval(1:2) = [2.0d0, 1.0d0]
      call osplin(jval, jout, jx, jy, 8, 2, ierr)
      write(*,'(a,i3,2f6.1)') "osplin decreasing: ierr, yval untouched", ierr, jout(1:2)
      jval(1:2) = [1.0d0, 2.0d0]
      call osplin(jval, jout, jx, jy, 8, 2, ierr)
      write(*,'(a,i3,2f6.1)') "osplin on nodes: ierr", ierr, jout(1:2)
      jx(8) = jx(7)
      jval(1:2) = [6.5d0, 7.5d0]
      jout = -1.0d0
      call osplin(jval, jout, jx, jy, 8, 2, ierr)
      write(*,'(a,i3,2f6.1)') "osplin, repeated last abscissa, a point above: ierr, yval untouched", &
           ierr, jout(1:2)
      jx(8) = 7.0d0

! ---------------- interpolation
      write(*,'(a)') "# inter3/interp/intrp2/lagrange4/quint"
      call inter3([1.0d0, 1.5d0, 2.5d0], w3, dw3, 1.8d0)
      call show("inter3 w", w3)
      call show("inter3 dw", dw3)
      x4 = [0.0d0, 1.0d0, 3.0d0, 4.5d0]
      call interp(x4, w4, dw4, 2.2d0)
      call show("interp w", w4)
      call show("interp dw", dw4)
      call intrp2(x4, w4, 5.0d0)
      call show("intrp2 w (extrapolated)", w4)
      write(*,'(a,1pe24.15)') "lagrange4", lagrange4(w4, [1.0d0, -2.0d0, 3.5d0, 0.25d0])
      call quint(0.37d0, 0.2d0, 0.1d0, 1.0d0, 1.3d0, 1.9d0, v)
      write(*,'(a,1pe24.15)') "quint", v
      write(*,'(a)') "# polint: five and three points, repeated abscissa"
      x20 = 0.0d0
      y20 = 0.0d0
      x20(1:5) = [0.1d0, 0.3d0, 0.35d0, 0.8d0, 1.2d0]
      y20(1:5) = 1.0d0/(1.0d0 + x20(1:5)**2)
      call polint(x20, y20, 5, 0.5d0, v, dv, ierr)
      write(*,'(a,2(1pe24.15),i3)') "polint n=5", v, dv, ierr
      call polint(x20, y20, 3, 0.0d0, v, dv, ierr)
      write(*,'(a,2(1pe24.15),i3)') "polint n=3 extrapolated", v, dv, ierr
      call polint(1.0d-25*x20, y20, 5, 0.5d-25, v, dv, ierr)
      write(*,'(a,2(1pe24.15),i3)') "polint n=5, abscissae scaled by 1e-25", v, dv, ierr
      x20(2) = x20(1)
      call polint(x20, y20, 3, 0.2d0, v, dv, ierr)
      write(*,'(a,i3)') "polint, repeated abscissa: ierr", ierr
      call lir_check()

! ---------------- search
      write(*,'(a)') "# findex: guesses below, above and outside the grid"
      x6 = [0.0d0, 0.5d0, 1.3d0, 2.0d0, 3.1d0, 4.0d0]
      pts(1:7) = [1.0d0, 1.0d0, 3.5d0, -0.1d0, 4.0d0, 0.5d0, 2.5d0]
      guesses(1:7) = [5, 1, 9, 3, 2, 6, 0]
      do i = 1, 7
         v = pts(i)
         idx = guesses(i)
         call findex(x6, 6, v, idx)
         write(*,'(a,f6.2,a,i4)') "findex", v, " ->", idx
      end do
      write(*,'(a)') "# locate: increasing and decreasing tables, edge slack"
      pts(1:6) = [1.6d0, 0.0d0, -0.01d0, 4.03d0, 4.5d0, -1.0d0]
      do i = 1, 6
         v = pts(i)
         call locate(x6, 6, v, idx)
         write(*,'(a,f6.2,a,i4)') "locate up", v, " ->", idx
      end do
      pts(1:2) = [0.995d0, 5.02d0]
      do i = 1, 2
         v = pts(i)
         call locate(x6 + 1.0d0, 6, v, idx)
         write(*,'(a,f6.3,a,i4)') "locate within slack", v, " ->", idx
      end do
      y6 = x6(6:1:-1) + 1.0d0
      pts(1:3) = [4.0d0, 5.0d0, 0.5d0]
      do i = 1, 3
         v = pts(i)
         call locate(y6, 6, v, idx)
         write(*,'(a,f6.2,a,i4)') "locate down", v, " ->", idx
      end do

! ---------------- linear algebra
      write(*,'(a)') "# ludcmp/lubksb: 4x4 in 5x5 storage, then a singular matrix"
      a5 = 0.0d0
      a5(1:4,1:4) = reshape([2.0d0, 1.0d0, 0.0d0, 4.0d0, &
                             1.0d0, -3.0d0, 2.0d0, 1.0d0, &
                             3.0d0, 0.5d0, 1.0d0, -2.0d0, &
                             0.0d0, 2.0d0, 5.0d0, 1.0d0], [4, 4])
      b5 = [1.0d0, 2.0d0, 3.0d0, 4.0d0, 0.0d0]
      call ludcmp(a5, 4, 5, indx, d, ierr)
      write(*,'(a,f5.1,i3,4i3)') "ludcmp d ierr indx", d, ierr, indx(1:4)
      call show("ludcmp LU", reshape(a5(1:4,1:4), [16]))
      call lubksb(a5, 4, 5, indx, b5)
      call show("lubksb x", b5(1:4))
      a5 = 0.0d0
      a5(1,1) = 1.0d0
      call ludcmp(a5, 3, 5, indx, d, ierr)
      write(*,'(a,i3)') "ludcmp singular: ierr", ierr
      a5 = 0.0d0
      a5(1:3,1:3) = reshape([1.0d0, 4.0d0, 7.0d0, 2.0d0, 5.0d0, 8.0d0, 3.0d0, 6.0d0, 9.0d0], [3, 3])
      call ludcmp(a5, 3, 5, indx, d, ierr)
      write(*,'(a,i3)') "ludcmp, rows (1,2,3),(4,5,6),(7,8,9) (a round-off pivot): ierr", ierr
      a5 = 1.0d0
      a5(1:3,1) = ieee_value(1.0d0, ieee_quiet_nan)
      call ludcmp(a5, 3, 5, indx, d, ierr)
      write(*,'(a,i3)') "ludcmp, a NaN column: ierr", ierr
      write(*,'(a)') "# simeqc: 3 unknowns (rows interchanged), two right-hand sides; then singular"
      sm(1:15) = [1.0d0, 0.5d0, 4.0d0, 3.0d0, 1.0d0, 2.0d0, 2.0d0, 6.0d0, 1.0d0, &
                  1.0d0, 2.0d0, 3.0d0, 0.5d0, -1.0d0, 2.0d0]
      call simeqc(sm(1:15), 5, 3, ierr)
      call show("simeqc", sm(1:15))
      write(*,'(a,i3)') "simeqc ierr", ierr
      sm(1:6) = [1.0d0, 2.0d0, 2.0d0, 4.0d0, 1.0d0, 1.0d0]
      call simeqc(sm(1:6), 3, 2, ierr)
      write(*,'(a,i3)') "simeqc singular: ierr", ierr
      call tridiagonal_check()

! ---------------- ODE
      call bsstep_check()

      write(*,'(a)') "test_numerics: done"

! ---------------- checks
      write(*,'(a)') "# checks"
      call cubic_spline_checks()
      call quadratic_spline_checks()
      call lagrange_checks()
      call search_checks()
      call linear_algebra_checks()
      call bsstep_checks()
      call status_checks()
      if (failures > 0) then
         write(*,'(i0,a)') failures, " check(s) failed"
         stop 1
      end if

contains

subroutine show(label, values)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: values(:)
      write(*,'(a)') label
      write(*,'(3(1pe24.15))') values
end subroutine show

! ysplin on its own: the coefficients of the cubic on each interval
subroutine ysplin_check()
      real(dp) :: c(4,100)
      c = 0.0d0
      c(1,1:8) = jy(1:8)
      call ysplin(jx(1:8), c, 8)
      call show("ysplin c(2:4,1:7)", reshape(c(2:4,1:7), [21]))
end subroutine ysplin_check

! osplin on jx/jy(1:n) at the given nondecreasing points
subroutine osplin_case(n, points)
      integer, intent(in) :: n
      real(dp), intent(in) :: points(:)
      integer :: m
      m = size(points)
      jval = 0.0d0
      jval(1:m) = points
      jout = 0.0d0
      call osplin(jval, jout, jx, jy, n, m)
      call show("osplin", jout(1:m))
end subroutine osplin_case

! lir: cubic, on a node, extrapolated low and high, linear (3 points),
! a decreasing table, fewer than two points, a flat table
subroutine lir_check()
      real(dp) :: tz(6), ty(2,6), res(2)
      write(*,'(a)') "# lir"
      tz = [1.0d0, 2.0d0, 3.5d0, 4.0d0, 6.0d0, 7.5d0]
      ty(1,:) = tz**2
      ty(2,:) = 10.0d0 - tz
      pts(1:5) = [2.7d0, 4.0d0, 0.5d0, 8.0d0, 1.2d0]
      do i = 1, 5
         v = pts(i)
         res = -1.0d0
         call lir(v, tz, res, ty, 2, 2, 6, flag)
         write(*,'(a,f6.2,a,2(1pe24.15),i3)') "lir cubic at", v, " =", res, flag
      end do
      res = -1.0d0
      call lir(2.7d0, tz, res, ty, 2, 2, 3, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir linear (3 points)", res, flag
      res = -1.0d0
      call lir(2.7d0, tz(6:1:-1), res, ty(:,6:1:-1), 2, 2, 6, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir decreasing table", res, flag
      res = -1.0d0
      flag = -7
      call lir(2.7d0, tz, res, ty, 2, 2, 1, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir one point (no change)", res, flag
      res = -1.0d0
      call lir(2.7d0, [3.0d0, 3.0d0, 3.0d0, 3.0d0], res, ty, 2, 2, 4, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir flat table (no change)", res, flag
      res = -1.0d0
      ty(1,3) = 0.0d0
      call lir(3.5d0, tz, res, ty, 2, 2, 6, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir zero on a node", res, flag
end subroutine lir_check

! solve_tridiagonal on one system, on its first n-2 rows, and on a
! singular system
subroutine tridiagonal_check()
      integer, parameter :: n = 6
      write(*,'(a)') "# solve_tridiagonal"
      ja = 0.0d0
      jb = 0.0d0
      jc = 0.0d0
      jr = 0.0d0
      do i = 1, n
         ja(i) = -1.0d0 + 0.1d0*i
         jb(i) = 4.0d0 + 0.5d0*i
         jc(i) = -1.5d0 + 0.05d0*i*i
         jr(i) = 1.0d0 + i*(7.0d0 - i)
      end do
      call solve_tridiagonal(ja, jb, jc, jr, n, jsol, ierr)
      call show("solve_tridiagonal", jsol(1:n))
      call solve_tridiagonal(ja, jb, jc, jr, n - 2, jsol, ierr)
      call show("solve_tridiagonal, first n-2 rows", jsol(1:n-2))
      jb(2) = ja(2)*jc(1)/jb(1)
      call solve_tridiagonal(ja, jb, jc, jr, n, jsol, ierr)
      write(*,'(a,i3)') "solve_tridiagonal singular: ierr", ierr
      jb(1) = 0.0d0
      call solve_tridiagonal(ja, jb, jc, jr, n, jsol, ierr)
      write(*,'(a,i3)') "solve_tridiagonal, zero first pivot: ierr", ierr
end subroutine tridiagonal_check

! bsstep: four accepted steps of the harmonic system from x = 0; large
! trial steps accepted at the sixth stage (h grows), the seventh (h
! shrinks) and after a retry; then an integrand failure and a step
! that underflows (against 0.7, as the envelope's hydrogen fraction)
subroutine bsstep_check()
      real(dp), parameter :: underflow_reference = 0.7d0
      type(harmonic_system) :: sys
      real(dp) :: y(3), dydx(3), yscale(3), err(3), x, h, h_did, h_next
      integer :: step, jerr, limits(3,6)
      write(*,'(a)') "# bsstep (mmid, ratext)"
      y = [0.0d0, 1.0d0, 1.0d0]
      x = 0.0d0
      h = 0.5d0
      yscale = 1.0d0
      do step = 1, 4
         call sys%derivatives(x, y, dydx, ierr)
         call bsstep(sys, y, dydx, 3, x, h, 1.0d-10, yscale, h_did, h_next, err, &
              11, 7, underflow_reference, ierr)
         write(*,'(a,i2,a,i3,a,i6)') "step", step, " ierr", ierr, " deriv calls", sys%calls
         call show("  y, x, h_did, h_next", [y, x, h_did, h_next])
         call show("  step_err", err)
         h = h_next
      end do
      pts(1:6) = [1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 4.0d0]
      do step = 1, 6
         h = pts(step)
         call bsstep(sys, y, dydx, 3, x, h, 1.0d-13, yscale, h_did, h_next, err, &
              11, 7, underflow_reference, ierr)
         write(*,'(a,f5.2,a,i3)') "large step, h_try", h, " ierr", ierr
         call show("  y, x, h_did, h_next", [y, x, h_did, h_next])
      end do
      sys%fail = .true.
      call bsstep(sys, y, dydx, 3, x, h, 1.0d-10, yscale, h_did, h_next, err, &
           11, 7, underflow_reference, ierr, jerr)
      write(*,'(a,2i3)') "bsstep, integrand failure: ierr, derivatives_ierr", ierr, jerr
      sys%fail = .false.
! zero tolerance: no step is accepted, so h shrinks until it vanishes
! against the underflow reference
      call bsstep(sys, y, dydx, 3, x, 0.5d0, 0.0d0, yscale, h_did, h_next, err, &
           11, 7, underflow_reference, ierr)
      write(*,'(a,i3,a,i8)') "bsstep, step underflow: ierr", ierr, " deriv calls", sys%calls
! NaN errors: with every component NaN (one equation) no stage is
! accepted and the step ends on underflow; a NaN in one of three
! components is passed over
      sys%nan_component = 1
      y = [0.0d0, 1.0d0, 1.0d0]
      x = 0.0d0
      call bsstep(sys, y, dydx, 1, x, 0.5d0, 1.0d-10, yscale, h_did, h_next, err, &
           11, 7, underflow_reference, ierr)
      write(*,'(a,i3)') "bsstep, one equation, NaN: ierr", ierr
      sys%nan_component = 3
      y = [0.0d0, 1.0d0, 1.0d0]
      x = 0.0d0
      call bsstep(sys, y, dydx, 3, x, 0.5d0, 1.0d-10, yscale, h_did, h_next, err, &
           11, 7, underflow_reference, ierr)
      write(*,'(a,i3,a,2(1pe24.15),l2)') "bsstep, NaN in component 3 of 3: ierr", ierr, &
           " y(1:2), y(3) is NaN", y(1:2), y(3) /= y(3)
! limits out of range: no step, numerics_bad_input (num_eqs,
! max_stage_index, extrap_order)
      limits = reshape([0, 11, 7,  16, 11, 7,  3, 12, 7,  3, 11, 1,  3, 11, 8,  3, 5, 7], [3, 6])
      do step = 1, 6
         sys%calls = 0
         x = 0.0d0
         call bsstep(sys, y, dydx, limits(1,step), x, 0.5d0, 1.0d-10, yscale, h_did, &
              h_next, err, limits(2,step), limits(3,step), underflow_reference, ierr)
         write(*,'(a,3i4,a,i3,a,i3,f5.1)') "bsstep limits", limits(:,step), ": ierr", ierr, &
              " deriv calls", sys%calls, x
      end do
end subroutine bsstep_check

! ======================================================================
! The checks: each prints its relative error (or 0/1 for a contract)
! and ok or FAIL against its tolerance.
! ======================================================================

! one check: its error against its tolerance
subroutine report(label, rel_err, tolerance)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: rel_err, tolerance
      if (rel_err <= tolerance) then
         write(*,'(a,t74,1pe10.2,a)') label, rel_err, "  ok"
      else
         write(*,'(a,t74,1pe10.2,a)') label, rel_err, "  FAIL"
         failures = failures + 1
      end if
end subroutine report

! a contract check: holds or not
subroutine report_true(label, holds)
      character(len=*), intent(in) :: label
      logical, intent(in) :: holds
      call report(label, merge(0.0d0, 1.0d0, holds), 0.0d0)
end subroutine report_true

! ---------------------------------------------------------------- splines
! The cubic f = x^3 - 2x + 1, its slope and a smooth rational g = 1/(1+x^2)
elemental function cubic(x) result(f)
      real(dp), intent(in) :: x
      real(dp) :: f
      f = x**3 - 2.0d0*x + 1.0d0
end function cubic

elemental function cubic_slope(x) result(f)
      real(dp), intent(in) :: x
      real(dp) :: f
      f = 3.0d0*x**2 - 2.0d0
end function cubic_slope

elemental function rational(x) result(f)
      real(dp), intent(in) :: x
      real(dp) :: f
      f = 1.0d0/(1.0d0 + x**2)
end function rational

elemental function rational_slope(x) result(f)
      real(dp), intent(in) :: x
      real(dp) :: f
      f = -2.0d0*x/(1.0d0 + x**2)**2
end function rational_slope

subroutine cubic_spline_checks()
      real(dp) :: x(6), y(6), y2(6), c(4,6), xe, ye, dye, err, derr, cerr, d
      real(dp) :: line_err, e_coarse, e_fine, de_coarse, de_fine
      integer :: i, k, ierr, klo
      x = [0.0d0, 0.5d0, 1.3d0, 2.0d0, 3.1d0, 4.0d0]
      y = [(cubic(x(i)), i = 1, 6)]
! clamped with the cubic's own end slopes: the spline is the cubic
      call cspline(x, y, 6, cubic_slope(x(1)), cubic_slope(x(6)), y2)
      call compute_spline_coefficients(x, y, y2, 6, c)
      err = 0.0d0
      derr = 0.0d0
      cerr = 0.0d0
      do k = 0, 80
         xe = 4.0d0*k/80.0d0
         call splint(x, y, 6, y2, xe, ye, ierr, klo=klo, dy=dye)
         err = max(err, abs(ye - cubic(xe))/49.0d0)
         derr = max(derr, abs(dye - cubic_slope(xe))/46.0d0)
         d = xe - x(klo)
         cerr = max(cerr, abs(((c(4,klo)*d + c(3,klo))*d + c(2,klo))*d + c(1,klo) - ye)/49.0d0)
      end do
      call report("cspline/splint: clamped spline exact on a cubic", err, 1.0d-14)
      call report("splint: its dy exact on the clamped cubic", derr, 1.0d-14)
      call report("compute_spline_coefficients: the power form is splint's", cerr, 1.0d-14)
! natural ends: zero second derivatives there, and exact on a line
      y = 2.0d0*x + 1.0d0
      call natural_spline(x, y, 6, y2)
      call report("natural_spline: zero second derivatives at the ends", &
           max(abs(y2(1)), abs(y2(6))), 0.0d0)
      line_err = 0.0d0
      do k = -4, 84
         xe = 4.0d0*k/80.0d0
         call splint(x, y, 6, y2, xe, ye, ierr, dy=dye)
         line_err = max(line_err, abs(ye - (2.0d0*xe + 1.0d0))/9.0d0, abs(dye - 2.0d0)/2.0d0)
      end do
      call report("natural_spline/splint: exact on a line, extrapolated too", line_err, 1.0d-14)
! the interior error of the clamped spline on 1/(1+x^2): h^4 in the
! value, h^3 in the slope
      call rational_errors(10, e_coarse, de_coarse)
      call rational_errors(20, e_fine, de_fine)
      call report("cspline: the error falls as h^4 (order to 0.2)", &
           abs(log(e_coarse/e_fine)/log(2.0d0) - 4.0d0), 0.2d0)
      call report("splint: the slope's error falls as h^3 (order to 0.3)", &
           abs(log(de_coarse/de_fine)/log(2.0d0) - 3.0d0), 0.3d0)
      call natural_spline_small_n_checks()
end subroutine cubic_spline_checks

! the clamped spline's largest errors in value and slope on [0, 2] with
! m equal intervals, sampled at 10 points per interval
subroutine rational_errors(m, err, derr)
      integer, intent(in) :: m
      real(dp), intent(out) :: err, derr
      real(dp) :: x(m+1), y(m+1), y2(m+1), xe, ye, dye
      integer :: i, k, ierr
      do i = 1, m+1
         x(i) = 2.0d0*(i-1)/m
         y(i) = rational(x(i))
      end do
      call cspline(x, y, m+1, rational_slope(x(1)), rational_slope(x(m+1)), y2)
      err = 0.0d0
      derr = 0.0d0
      do k = 0, 10*m
         xe = 2.0d0*k/(10*m)
         call splint(x, y, m+1, y2, xe, ye, ierr, dy=dye)
         err = max(err, abs(ye - rational(xe)))
         derr = max(derr, abs(dye - rational_slope(xe)))
      end do
end subroutine rational_errors

! ysplin and intpol (kept for now) against natural_spline and splint
! for n = 2, 3, 4 and 8, and the natural spline's slope at the middle of
! (0,0), (1,1), (3,0), which is 1/2
subroutine natural_spline_small_n_checks()
      real(dp) :: x(8), y(8), y2(8), c(4,8), ye, dye, yi, dyi, err, xe
      integer :: sizes(4), m, n, i, k, ierr
      character(len=80) :: label
      sizes = [2, 3, 4, 8]
      do m = 1, 4
         n = sizes(m)
         do i = 1, n
            x(i) = 0.5d0*i + 0.1d0*i*i
            y(i) = rational(x(i))
         end do
         call natural_spline(x(1:n), y(1:n), n, y2(1:n))
         c = 0.0d0
         c(1,1:n) = y(1:n)
         call ysplin(x(1:n), c, n)
         err = 0.0d0
         do i = 1, n
            call splint(x(1:n), y(1:n), n, y2(1:n), x(i), ye, ierr, dy=dye)
            err = max(err, abs(c(2,i) - dye))
         end do
         do k = 0, 20
            xe = x(1) - 0.2d0 + (x(n) - x(1) + 0.4d0)*k/20.0d0
            call splint(x(1:n), y(1:n), n, y2(1:n), xe, ye, ierr, dy=dye)
            call intpol(x(1:n), y(1:n), n, xe, yi, dyi, ierr)
            err = max(err, abs(yi - ye), abs(dyi - dye))
         end do
         write(label,'(a,i0)') "ysplin/intpol: the natural spline, n = ", n
         call report(trim(label), err, 1.0d-15)
      end do
      x(1:3) = [0.0d0, 1.0d0, 3.0d0]
      y(1:3) = [0.0d0, 1.0d0, 0.0d0]
      c = 0.0d0
      c(1,1:3) = y(1:3)
      call ysplin(x(1:3), c, 3)
      call report("ysplin: slope 1/2 at the middle of (0,0),(1,1),(3,0)", &
           abs(c(2,2) - 0.5d0), 1.0d-15)
end subroutine natural_spline_small_n_checks

! ---------------------------------------------------------------- osplin
subroutine quadratic_spline_checks()
      integer, parameter :: m = 400
      real(dp) :: x(9), y(9), xe(m), ye(m), diff_min, curv_min, err, keep(2), &
           xr(5), yr(5), inside(3)
      integer :: i, k, ierr, ierr_below, ierr_inside
      x = [0.0d0, 0.4d0, 1.0d0, 1.5d0, 2.5d0, 3.0d0, 4.2d0, 5.0d0, 6.5d0]
! monotone and convex data: (x+1)^2/2 + 1/(x+1)^... take f = x^2/2 + x/4
      y = 0.5d0*x**2 + 0.25d0*x
      do k = 1, m
         xe(k) = x(1) + (x(9) - x(1))*(k-1)/(m-1)
      end do
      call osplin(xe, ye, x, y, 9, m, ierr)
      diff_min = minval(ye(2:m) - ye(1:m-1))
      curv_min = minval(ye(3:m) - 2.0d0*ye(2:m-1) + ye(1:m-2))
      call report("osplin: monotone data give a monotone spline", &
           max(0.0d0, -diff_min)/maxval(abs(y)), 0.0d0)
      call report("osplin: convex data give a convex spline", &
           max(0.0d0, -curv_min)/maxval(abs(y)), 1.0d-14)
! the nodes themselves
      call osplin(x, ye(1:9), x, y, 9, 9, ierr)
      call report("osplin: interpolates the nodes", maxval(abs(ye(1:9) - y)), 0.0d0)
! a straight line is reproduced, extrapolation included
      y = 3.0d0 - 0.5d0*x
      do k = 1, m
         xe(k) = -1.0d0 + 9.0d0*(k-1)/(m-1)
      end do
      call osplin(xe, ye, x, y, 9, m, ierr)
      call report("osplin: exact on a line, extrapolated too", &
           maxval(abs(ye - (3.0d0 - 0.5d0*xe)))/3.0d0, 1.0d-13)
! every sorted run of points equals its points one at a time
      y = [0.0d0, 1.0d0, 0.5d0, 2.0d0, -1.0d0, -1.0d0, 3.0d0, 2.5d0, 4.0d0]
      call report("osplin: every sorted subset of 16 points equals them one at a time", &
           subsets_against_points(x, y), 0.0d0)
! the error contract: too few points, decreasing points
      keep = [-7.0d0, -7.0d0]
      call osplin([1.0d0, 2.0d0], keep, x, y, 2, 2, ierr)
      call report_true("osplin: n < 3 is numerics_bad_input, yval untouched", &
           ierr == numerics_bad_input .and. all(keep == -7.0d0))
      call osplin([2.0d0, 1.0d0], keep, x, y, 9, 2, ierr)
      call report_true("osplin: decreasing points are numerics_unsorted_points", &
           ierr == numerics_unsorted_points .and. all(keep == -7.0d0))
! a repeated end abscissa: no extrapolation from the zero-width end
! interval (it gave NaN), but the points inside are evaluated
      xr = [0.0d0, 1.0d0, 2.0d0, 3.0d0, 3.0d0]
      yr = xr**2
      call osplin([2.5d0, 3.5d0], keep, xr, yr, 5, 2, ierr)
      xr = [0.0d0, 0.0d0, 1.0d0, 2.0d0, 3.0d0]
      yr = xr**2
      call osplin([-0.5d0, 1.5d0], keep, xr, yr, 5, 2, ierr_below)
      call osplin([0.0d0, 0.5d0, 3.0d0], inside, xr, yr, 5, 3, ierr_inside)
      call report_true("osplin: no extrapolation from a zero-width end interval", &
           ierr == numerics_coincident_points .and. ierr_below == numerics_coincident_points &
           .and. all(keep == -7.0d0) .and. ierr_inside == numerics_ok &
           .and. inside(1) == 0.0d0 .and. inside(3) == 9.0d0 .and. abs(inside(2)) < 1.0d0)
end subroutine quadratic_spline_checks

! The fraction of points whose value in a run differs from their value
! alone (bit for bit), over every nonempty sorted subset of 16 candidate
! points on the 9-point table x/y: two below the table, the first node
! and a point in its interval, a node twice and a point after it, the
! third-last node and a point after it, the last-but-one node twice and
! a point in the last interval, the last node twice and two above.
function subsets_against_points(x, y) result(mismatch_fraction)
      real(dp), intent(in) :: x(9), y(9)
      real(dp) :: mismatch_fraction
      integer, parameter :: n_cand = 16
      real(dp), parameter :: candidates(n_cand) = [-2.0d0, -0.5d0, 0.0d0, 0.2d0, &
           1.0d0, 1.0d0, 1.2d0, 4.2d0, 4.6d0, 5.0d0, 5.0d0, 5.7d0, 6.5d0, 6.5d0, &
           7.0d0, 8.0d0]
      real(dp) :: alone(n_cand), points(n_cand), run_y(n_cand), y1(1)
      integer :: c, mask, m, picked(n_cand), jerr, mismatches, compared
      do c = 1, n_cand
         call osplin(candidates(c:c), y1, x, y, 9, 1, jerr)
         alone(c) = y1(1)
      end do
      mismatches = 0
      compared = 0
      do mask = 1, 2**n_cand - 1
         m = 0
         do c = 1, n_cand
            if (.not. btest(mask, c - 1)) cycle
            m = m + 1
            picked(m) = c
            points(m) = candidates(c)
         end do
         call osplin(points(1:m), run_y(1:m), x, y, 9, m, jerr)
         do c = 1, m
            compared = compared + 1
            if (.not. (run_y(c) == alone(picked(c)))) mismatches = mismatches + 1
         end do
      end do
      mismatch_fraction = real(mismatches, dp)/real(compared, dp)
end function subsets_against_points

! ---------------------------------------------------------------- Lagrange
subroutine lagrange_checks()
      real(dp) :: x3(3), w3(3), dw3(3), x4(4), w4(4), dw4(4), w4b(4), xe, err, derr
      real(dp) :: xa(5), ya(5), ye, dy, y0, y1, y2, sum_err
      integer :: k, ierr
      x3 = [1.0d0, 1.5d0, 2.5d0]
      x4 = [0.0d0, 1.0d0, 3.0d0, 4.5d0]
      err = 0.0d0
      derr = 0.0d0
      sum_err = 0.0d0
      do k = 0, 12
         xe = -0.5d0 + 0.5d0*k
         call inter3(x3, w3, dw3, xe)
         err = max(err, abs(sum(w3*quadratic(x3)) - quadratic(xe))/10.0d0)
         derr = max(derr, abs(sum(dw3*quadratic(x3)) - quadratic_slope(xe))/10.0d0)
      end do
      call report("inter3: exact on a quadratic", err, 1.0d-13)
      call report("inter3: its derivative weights exact on a quadratic", derr, 1.0d-13)
      err = 0.0d0
      derr = 0.0d0
      do k = 0, 12
         xe = -0.5d0 + 0.5d0*k
         call interp(x4, w4, dw4, xe)
         err = max(err, abs(lagrange4(w4, cubic(x4)) - cubic(xe))/100.0d0)
         derr = max(derr, abs(lagrange4(dw4, cubic(x4)) - cubic_slope(xe))/100.0d0)
         sum_err = max(sum_err, abs(sum(w4) - 1.0d0), abs(sum(dw4)))
         call intrp2(x4, w4b, xe)
         sum_err = max(sum_err, maxval(abs(w4b - w4)))
      end do
      call report("interp/lagrange4: exact on a cubic", err, 1.0d-13)
      call report("interp: its derivative weights exact on a cubic", derr, 1.0d-13)
      call report("interp: weights sum to 1, derivatives to 0; intrp2's", sum_err, 1.0d-13)
      call interp(x4, w4, dweight=dw4, x_eval=2.2d0)
      call intrp2(x4, w4b, 2.2d0)
      call report_true("intrp2: interp's weights bit for bit", all(w4b == w4))
      err = 0.0d0
      do k = 0, 8
         xe = 0.15d0 + 0.05d0*k
         y0 = quadratic(0.2d0)
         y1 = quadratic(0.3d0)
         y2 = quadratic(0.4d0)
         call quint(xe, 0.2d0, 0.1d0, y0, y1, y2, ye)
         err = max(err, abs(ye - quadratic(xe)))
      end do
      call report("quint: exact on a quadratic", err, 1.0d-14)
      xa = [0.1d0, 0.3d0, 0.35d0, 0.8d0, 1.2d0]
      ya = xa**4 - xa
      err = 0.0d0
      do k = 0, 10
         xe = 0.13d0*k
         call polint(xa, ya, 5, xe, ye, dy, ierr)
         err = max(err, abs(ye - (xe**4 - xe)))
      end do
      call report("polint: exact on a quartic through five points", err, 1.0d-13)
! the same, on abscissae scaled by 1e-25: no false coincidence
      call polint(1.0d-25*xa, ya, 5, 1.0d-25*0.5d0, ye, dy, ierr)
      call report("polint: scale-free (abscissae of 1e-25)", &
           merge(abs(ye - (0.5d0**4 - 0.5d0)), 1.0d0, ierr == 0), 1.0d-13)
      xa(2) = xa(1)
      call polint(xa, ya, 3, 0.2d0, ye, dy, ierr)
      call report_true("polint: coincident abscissae are numerics_coincident_points", &
           ierr == numerics_coincident_points)
      call lir_checks()
end subroutine lagrange_checks

elemental function quadratic(x) result(f)
      real(dp), intent(in) :: x
      real(dp) :: f
      f = 2.0d0*x**2 - 3.0d0*x + 0.5d0
end function quadratic

elemental function quadratic_slope(x) result(f)
      real(dp), intent(in) :: x
      real(dp) :: f
      f = 4.0d0*x - 3.0d0
end function quadratic_slope

! lir: exact on cubics (4 points or more), on lines (fewer), on nodes;
! a decreasing table as the same table increasing; its error contract
subroutine lir_checks()
      real(dp) :: tz(6), ty(2,6), tr(2,6), res(2), res_down(2), err, err_lin, err_down
      integer :: k, flag, ierr
      tz = [1.0d0, 2.0d0, 3.5d0, 4.0d0, 6.0d0, 7.5d0]
      do k = 1, 6
         ty(1,k) = cubic(tz(k))
         ty(2,k) = 10.0d0 - tz(k)
         tr(1,k) = rational(tz(k))
         tr(2,k) = rational(0.5d0*tz(k))
      end do
      err = 0.0d0
      err_down = 0.0d0
      do k = 0, 18
         call lir(0.5d0*k, tz, res, ty, 2, 2, 6, flag)
         err = max(err, abs(res(1) - cubic(0.5d0*k))/cubic(9.0d0), &
              abs(res(2) - (10.0d0 - 0.5d0*k))/10.0d0)
! on 1/(1+z^2), where the stencil matters
         call lir(0.5d0*k, tz, res, tr, 2, 2, 6, flag)
         call lir(0.5d0*k, tz(6:1:-1), res_down, tr(:,6:1:-1), 2, 2, 6, flag)
         err_down = max(err_down, maxval(abs(res_down - res)))
      end do
      call report("lir: exact on a cubic and a line (6 points)", err, 1.0d-13)
      call report("lir: a decreasing table as the same table increasing", err_down, 1.0d-15)
      err_lin = 0.0d0
      do k = 0, 10
         call lir(0.5d0*k, tz, res, ty, 2, 2, 3, flag)
         err_lin = max(err_lin, abs(res(2) - (10.0d0 - 0.5d0*k))/10.0d0)
      end do
      call report("lir: linear below four points, exact on a line", err_lin, 1.0d-14)
      call lir(3.5d0, tz, res, ty, 2, 2, 6, flag)
      call report_true("lir: on a node, the tabulated values (flag 1)", &
           all(res == ty(:,3)) .and. flag == 1)
      call lir(9.0d0, tz, res, ty, 2, 2, 6, flag)
      call report_true("lir: outside the table, flag 0", flag == 0)
      res = -1.0d0
      call lir(2.7d0, tz, res, ty, 2, 2, 1, flag, ierr)
      call report_true("lir: one point is numerics_bad_input, result untouched", &
           ierr == numerics_bad_input .and. all(res == -1.0d0))
      call lir(2.7d0, [3.0d0, 3.0d0, 3.0d0, 3.0d0], res, ty, 2, 2, 4, flag, ierr)
      call report_true("lir: a flat table is numerics_coincident_points", &
           ierr == numerics_coincident_points .and. all(res == -1.0d0))
end subroutine lir_checks

! ---------------------------------------------------------------- searches
! the search checks' points: 61 across and beyond the table (k = 0..60),
! then its seven nodes (k = 61..67)
pure function sample(k) result(xe)
      integer, intent(in) :: k
      real(dp) :: xe
      real(dp), parameter :: nodes(7) = [-2.0d0, -0.5d0, 0.0d0, 0.5d0, 1.3d0, 2.0d0, 4.0d0]
      if (k <= 60) then
         xe = -2.5d0 + 7.0d0*k/60.0d0
      else
         xe = nodes(k - 60)
      end if
end function sample

subroutine search_checks()
      real(dp) :: up(7), down(7), xe
      integer :: k, j, guess
      logical :: holds
      up = [-2.0d0, -0.5d0, 0.0d0, 0.5d0, 1.3d0, 2.0d0, 4.0d0]
      down = up(7:1:-1)
! locate: xx(j) < x <= xx(j+1) inside (xx(j) >= x > xx(j+1) for a
! decreasing table), 1 and n-1 on the end nodes, the end interval within
! 1% of the span beyond an end, 0 and n beyond that
      holds = .true.
      do k = 1, 60
         xe = -2.0d0 + 6.0d0*k/60.0d0
         call locate(up, 7, xe, j)
         holds = holds .and. j >= 1 .and. j <= 6
         if (holds) holds = up(j) < xe .and. xe <= up(j+1)
         call locate(down, 7, xe, j)
         holds = holds .and. j >= 1 .and. j <= 6
         if (holds) holds = down(j) >= xe .and. xe > down(j+1)
      end do
      call report_true("locate: brackets x in increasing and decreasing tables", holds)
      call locate(up, 7, up(1), j)
      holds = j == 1
      call locate(up, 7, up(7), j)
      holds = holds .and. j == 6
      call locate(down, 7, down(1), j)
      holds = holds .and. j == 1
      call locate(down, 7, down(7), j)
      holds = holds .and. j == 6
      call report_true("locate: the end nodes give 1 and n-1, either direction", holds)
      call locate(up, 7, -2.05d0, j)
      holds = j == 1
      call locate(up, 7, 4.05d0, j)
      holds = holds .and. j == 6
      call locate(down, 7, 4.05d0, j)
      holds = holds .and. j == 1
      call locate(down, 7, -2.05d0, j)
      holds = holds .and. j == 6
      call report_true("locate: within 1% of the span beyond an end, the end interval", holds)
      call locate(up, 7, -2.1d0, j)
      holds = j == 0
      call locate(up, 7, 4.1d0, j)
      holds = holds .and. j == 7
      call locate(down, 7, 4.1d0, j)
      holds = holds .and. j == 0
      call locate(down, 7, -2.1d0, j)
      holds = holds .and. j == 7
      call report_true("locate: beyond the slack, 0 and n, either direction", holds)
! findex and hunt from every starting guess, at the nodes too
      holds = .true.
      do k = 0, 67
         xe = sample(k)
         do guess = -1, 8
            j = guess
            call findex(up, 7, xe, j)
            if (xe < up(1)) then
               holds = holds .and. j == -1
            else if (xe >= up(7)) then
               holds = holds .and. j == -7
            else
               holds = holds .and. up(j) <= xe .and. xe < up(j+1)
            end if
         end do
      end do
      call report_true("findex: brackets x from any guess; -1 below, -n at the top", holds)
      holds = .true.
      do k = 0, 67
         xe = sample(k)
         do guess = -1, 8
            j = guess
            call hunt(up, 7, xe, j)
            if (xe < up(1)) then
               holds = holds .and. j == 0
            else if (xe >= up(7)) then
               holds = holds .and. j == 7
            else
               holds = holds .and. up(j) <= xe .and. xe < up(j+1)
            end if
            j = guess
            call hunt(down, 7, xe, j)
            if (xe >= down(1)) then
               holds = holds .and. j == 0
            else if (xe < down(7)) then
               holds = holds .and. j == 7
            else
               holds = holds .and. down(j) > xe .and. xe >= down(j+1)
            end if
         end do
      end do
      call report_true("hunt: brackets x from any guess, either direction", holds)
end subroutine search_checks

! ---------------------------------------------------------------- linear algebra
subroutine linear_algebra_checks()
      integer, parameter :: n = 6
      real(dp) :: a(n,n), lu(n,n), b(n), x(n), d, res, sub(n), diag(n), sup(n), &
           xt(n), xs(n), gam(n), dense(n,n), sm(3*5), a3(3,3), ref3(3,3), xa(3), xb(3), &
           ref(3,2)
      real(dp) :: nan
      integer :: indx(n), i, j, ierr, indx3(3)
      do j = 1, n
         do i = 1, n
            a(i,j) = 1.0d0/(i + 2*j - 1) + merge(3.0d0, 0.0d0, i == j)
         end do
         b(j) = 1.0d0 - 0.3d0*j
      end do
! LU: the residual of a solve
      lu = a
      call ludcmp(lu, n, n, indx, d, ierr)
      x = b
      call lubksb(lu, n, n, indx, x)
      res = maxval(abs(matmul(a, x) - b))/(maxval(abs(a))*maxval(abs(x)))
      call report("ludcmp/lubksb: the residual of a 6x6 solve", res, 1.0d-15)
! rank-deficient matrices: rows (1,2,3), (4,5,6), (7,8,9), whose last
! pivot is round-off (1.9e-17 of its row), and row 3 = row 1 + row 2,
! whose last pivot is zero; and a zero row
      a3 = reshape([1.0d0, 4.0d0, 7.0d0,  2.0d0, 5.0d0, 8.0d0,  3.0d0, 6.0d0, 9.0d0], [3, 3])
      call ludcmp(a3, 3, 3, indx3, d, ierr)
      call report_true("ludcmp: a round-off pivot (rows 1..9) is numerics_singular", &
           ierr == numerics_singular)
      a3 = reshape([1.0d0, 4.0d0, 5.0d0,  2.0d0, 5.0d0, 7.0d0,  3.0d0, 6.0d0, 9.0d0], [3, 3])
      call ludcmp(a3, 3, 3, indx3, d, ierr)
      call report_true("ludcmp: a zero pivot (row 3 = row 1 + row 2) is numerics_singular", &
           ierr == numerics_singular)
      a3 = reshape([1.0d0, 0.0d0, 2.0d0,  2.0d0, 0.0d0, 1.0d0,  3.0d0, 0.0d0, 1.0d0], [3, 3])
      ref3 = a3
      call ludcmp(a3, 3, 3, indx3, d, ierr)
      call report_true("ludcmp: a zero row is numerics_singular, the matrix untouched", &
           ierr == numerics_singular .and. all(a3 == ref3))
      nan = ieee_value(1.0d0, ieee_quiet_nan)
      a3 = reshape([nan, nan, nan,  2.0d0, 5.0d0, 7.0d0,  3.0d0, 6.0d0, 1.0d0], [3, 3])
      call ludcmp(a3, 3, 3, indx3, d, ierr)
      call report_true("ludcmp: a NaN column is numerics_bad_input", ierr == numerics_bad_input)
! simeqc against LU on two right-hand sides, for a matrix whose first
! two columns need a row interchange (their largest elements are in
! rows 3 and 1 after it)
      a3 = reshape([1.0d0, 0.5d0, 4.0d0,  3.0d0, 1.0d0, 2.0d0,  2.0d0, 6.0d0, 1.0d0], [3, 3])
      sm(1:9) = reshape(a3, [9])
      sm(10:12) = [1.0d0, 2.0d0, 3.0d0]
      sm(13:15) = [0.5d0, -1.0d0, 2.0d0]
      ref(:,1) = sm(10:12)
      ref(:,2) = sm(13:15)
      call simeqc(sm, 5, 3, ierr)
      call ludcmp(a3, 3, 3, indx3, d, ierr)
      xa = ref(:,1)
      xb = ref(:,2)
      call lubksb(a3, 3, 3, indx3, xa)
      call lubksb(a3, 3, 3, indx3, xb)
      call report("simeqc: LU's solution for two right-hand sides", &
           max(maxval(abs(sm(10:12) - xa)), maxval(abs(sm(13:15) - xb)))/maxval(abs(xa)), 1.0d-15)
      sm(1:6) = [1.0d0, 2.0d0, 2.0d0, 4.0d0, 1.0d0, 1.0d0]
      call simeqc(sm(1:6), 3, 2, ierr)
      call report_true("simeqc: a singular matrix is numerics_singular", ierr == numerics_singular)
! rows (1,2,3), (4,5,6), (7,8,9): the last pivot is round-off
      sm(1:12) = [1.0d0, 4.0d0, 7.0d0,  2.0d0, 5.0d0, 8.0d0,  3.0d0, 6.0d0, 9.0d0, &
           1.0d0, 1.0d0, 1.0d0]
      call simeqc(sm(1:12), 4, 3, ierr)
      call report_true("simeqc: a round-off pivot (rows 1..9) is numerics_singular", &
           ierr == numerics_singular)
! a well-conditioned matrix with a row scaled by 1e-30 is solved
      sm(1:6) = [2.0d0, 1.0d-30,  1.0d0, 3.0d-30,  3.0d0, 4.0d-30]
      call simeqc(sm(1:6), 3, 2, ierr)
      call report_true("simeqc: a row scaled by 1e-30 is not singular, and solved", &
           ierr == numerics_ok .and. abs(sm(5) - 1.0d0) < 1.0d-15 .and. abs(sm(6) - 1.0d0) < 1.0d-15)
      sm(1:6) = [nan, nan,  1.0d0, 2.0d0,  1.0d0, 1.0d0]
      call simeqc(sm(1:6), 3, 2, ierr)
      call report_true("simeqc: a NaN column is numerics_bad_input", ierr == numerics_bad_input)
! the tridiagonal solve against a dense one, and its two halves
      do i = 1, n
         sub(i) = -1.0d0 + 0.1d0*i
         diag(i) = 4.0d0 + 0.5d0*i
         sup(i) = -1.5d0 + 0.05d0*i*i
         b(i) = 1.0d0 + i*(7.0d0 - i)
      end do
      dense = 0.0d0
      do i = 1, n
         dense(i,i) = diag(i)
         if (i > 1) dense(i,i-1) = sub(i)
         if (i < n) dense(i,i+1) = sup(i)
      end do
      call solve_tridiagonal(sub, diag, sup, b, n, xt, ierr)
      call ludcmp(dense, n, n, indx, d, ierr)
      x = b
      call lubksb(dense, n, n, indx, x)
      call report("solve_tridiagonal: a dense solve's solution", &
           maxval(abs(xt - x))/maxval(abs(x)), 1.0d-15)
      call tridiagonal_forward_sweep(sub, diag, sup, b, n, xs, gam, ierr)
      call tridiagonal_back_substitution(gam, n, xs)
      call report("tridiagonal sweep + back substitution: a dense solve's solution", &
           maxval(abs(xs - x))/maxval(abs(x)), 1.0d-15)
      diag(1) = 0.0d0
      call solve_tridiagonal(sub, diag, sup, b, n, xt, ierr)
      call report_true("solve_tridiagonal: a zero first pivot is numerics_singular", &
           ierr == numerics_singular)
      diag(1) = 4.5d0
      diag(2) = sub(2)*sup(1)/diag(1)
      call solve_tridiagonal(sub, diag, sup, b, n, xt, ierr)
      call report_true("solve_tridiagonal: a zero later pivot is numerics_singular", &
           ierr == numerics_singular)
end subroutine linear_algebra_checks

! ---------------------------------------------------------------- bsstep
! The harmonic system from y = (0, 1, 1) at x = 0 is (sin x, cos x,
! exp(-0.3 x)); bsstep's steps to x = 5 at tolerance 1e-10.
subroutine bsstep_checks()
      double precision, parameter :: underflow_reference = 0.7d0
      type(harmonic_system) :: sys
      real(dp) :: y(3), dydx(3), yscale(3), err(3), x, h, h_did, h_next, x_before, &
           big_y(16), big_dydx(16), big_scale(16), big_err(16), worst
      integer :: ierr, jerr, k, limits(3,6)
      logical :: advanced
      y = [0.0d0, 1.0d0, 1.0d0]
      x = 0.0d0
      h = 0.3d0
      yscale = 1.0d0
      worst = 0.0d0
      advanced = .true.
      do while (x < 5.0d0)
         call sys%derivatives(x, y, dydx, ierr)
         x_before = x
         call bsstep(sys, y, dydx, 3, x, min(h, 5.0d0 - x), 1.0d-10, yscale, h_did, &
              h_next, err, 11, 7, underflow_reference, ierr)
         advanced = advanced .and. ierr == 0 .and. x == x_before + h_did
         worst = max(worst, abs(y(1) - sin(x)), abs(y(2) - cos(x)), &
              abs(y(3) - exp(-0.3d0*x)))
         h = h_next
      end do
      call report("bsstep: sin x, cos x and exp(-0.3x) to x = 5", worst, 1.0d-9)
      call report_true("bsstep: x advances by h_did on every step", advanced)
! the limits: none is used, the derivatives are not called
      limits = reshape([0, 11, 7,  16, 11, 7,  3, 12, 7,  3, 11, 1,  3, 11, 8,  3, 5, 7], [3, 6])
      advanced = .true.
      big_y = 1.0d0
      big_dydx = 0.0d0
      big_scale = 1.0d0
      do k = 1, 6
         sys%calls = 0
         x = 0.0d0
         call bsstep(sys, big_y, big_dydx, limits(1,k), x, 0.5d0, 1.0d-10, big_scale, h_did, &
              h_next, big_err, limits(2,k), limits(3,k), underflow_reference, ierr)
         advanced = advanced .and. ierr == numerics_bad_input .and. sys%calls == 0 .and. x == 0.0d0
      end do
      call report_true("bsstep: limits out of range are numerics_bad_input, no step", advanced)
! the derivatives' failure and a vanishing step
      y = [0.0d0, 1.0d0, 1.0d0]
      x = 0.0d0
      call sys%derivatives(x, y, dydx, ierr)
      sys%fail = .true.
      call bsstep(sys, y, dydx, 3, x, 2.0d0, 1.0d-10, yscale, h_did, h_next, err, &
           11, 7, underflow_reference, ierr, jerr)
      call report_true("bsstep: failed derivatives are numerics_derivatives_failed", &
           ierr == numerics_derivatives_failed .and. jerr == 3)
      sys%fail = .false.
      call bsstep(sys, y, dydx, 3, x, 0.5d0, 0.0d0, yscale, h_did, h_next, err, &
           11, 7, underflow_reference, ierr)
      call report_true("bsstep: zero tolerance ends in numerics_step_underflow", &
           ierr == numerics_step_underflow)
end subroutine bsstep_checks

! every status code has its own description
subroutine status_checks()
      integer :: codes(7), i, j
      logical :: holds
      codes = [numerics_ok, numerics_singular, numerics_coincident_points, &
           numerics_bad_input, numerics_step_underflow, numerics_unsorted_points, &
           numerics_derivatives_failed]
      holds = .true.
      do i = 1, 7
         holds = holds .and. numerics_status_text(codes(i)) /= numerics_status_text(-99)
         do j = 1, i-1
            holds = holds .and. codes(i) /= codes(j) .and. &
                 numerics_status_text(codes(i)) /= numerics_status_text(codes(j))
         end do
      end do
      call report_true("numerics_status: distinct codes, each with its own text", holds)
end subroutine status_checks

end program test_numerics

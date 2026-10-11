!----------------------------------------------------------------------
! test_numerics
!----------------------------------------------------------------------
! Standalone test of numerics_lib: calls every public routine on small
! hand-made data and prints what it returns, error paths included, for
! byte-comparison against expected_test_numerics.out. The data are
! polynomials and rationals, so no elementary function enters the
! inputs. Families:
!   cubic splines     cspline/splint, kspline/ksplint, splinc/splintd2,
!                     ysplin/intpol
!   quadratic spline  osplin (McAllister-Roulier), through all four of
!                     its knot cases, extrapolation and single points
!   interpolation     inter3, interp, intrp2, lagrange4, quint, polint,
!                     lir (cubic, linear, on a node, extrapolated)
!   search            findex, locate
!   linear algebra    ludcmp/lubksb, simeqc, tridiag_gs, ctridi
!   ODE               bsstep (with mmid and ratext), its integrand
!                     failure and its step underflow
module test_numerics_callbacks
      implicit none
contains

! dy1/dx = y2, dy2/dx = -y1, dy3/dx = -L y3: a rotation and a decay.
! Counts its calls in saha_state; fails past x = fail_x when
! conductive_opacity_flag is set.
subroutine harmonic_deriv(x, y, dydx, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, log10_gravity, &
     in_atmosphere, want_derivatives, conductive_opacity_flag, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     saha_state, ierr)
      double precision, intent(inout) :: x, y(3), dydx(3)
      double precision, intent(inout) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(inout) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag
      double precision, intent(inout) :: log10_radius, log10_teff, &
           hydrogen_fraction, metal_fraction
      integer, intent(inout) :: saha_state
      integer, intent(out) :: ierr
      double precision, parameter :: fail_x = 0.75d0
      ierr = 0
      saha_state = saha_state + 1
      if (conductive_opacity_flag .and. x > fail_x) then
         ierr = 3
         return
      end if
      dydx(1) = y(2)
      dydx(2) = -y(1)
      dydx(3) = -luminosity_linear*y(3)
end subroutine harmonic_deriv

end module test_numerics_callbacks

program test_numerics
      use numerics_lib
      use integration_parameters, only: max_stage_index, extrap_order
      use io_units, only: run_log_unit, terminal_unit
      use star_info_lib, only: json
      use test_numerics_callbacks
      implicit none
      double precision :: jx(json), jy(json), jy2(json), jval(json), jout(json)
      double precision :: ja(json), jb(json), jc(json), jr(json), jsol(json)
      double precision :: x6(6), y6(6), y26(6), v, dv, d
      double precision :: x4(4), y4(4), y24(4), w4(4), dw4(4), w3(3), dw3(3)
      double precision :: x20(20), y20(20)
      double precision :: a5(5,5), b5(5), sm(15)
      integer :: indx(5), i, idx, klo, khi, ierr, flag, guesses(7)
      double precision :: pts(10), pts2(10)

      terminal_unit = 6
      run_log_unit = 6

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
         call splint(x6, y6, 6, y26, v, d, klo, khi, ierr)
         write(*,'(a,f6.2,a,1pe24.15,3i3)') "splint at", v, " =", d, klo, khi, ierr
      end do
! bisection gives a zero-width interval only between equal last points
      x6(6) = x6(5)
      call splint(x6, y6, 6, y26, 5.0d0, d, klo, khi, ierr)
      write(*,'(a,3i3)') "splint, repeated abscissa: klo khi ierr", klo, khi, ierr
      x6(6) = 4.0d0

      write(*,'(a)') "# kspline/ksplint: four points"
      x4 = [1.0d0, 2.0d0, 4.0d0, 7.0d0]
      y4 = 1.0d0/x4
      call kspline(x4, y4, y24)
      call show("kspline y2", y24)
      pts(1:3) = [1.5d0, 4.0d0, 6.5d0]
      do i = 1, 3
         v = pts(i)
         call ksplint(x4, y4, y24, v, d, ierr)
         write(*,'(a,f6.2,a,1pe24.15,i3)') "ksplint at", v, " =", d, ierr
      end do
      x4(4) = x4(3)
      call ksplint(x4, y4, y24, 8.0d0, d, ierr)
      write(*,'(a,i3)') "ksplint, repeated abscissa: ierr", ierr

      write(*,'(a)') "# splinc/splintd2: seven points in json-sized arrays"
      jx = 0.0d0
      jy = 0.0d0
      do i = 1, 7
         jx(i) = 0.3d0*i + 0.01d0*i*i
         jy(i) = jx(i)**2/(1.0d0 + jx(i))
      end do
      call splinc(jx, jy, jy2, 7)
      call show("splinc y2", jy2(1:7))
      call splintd2(jx, jy, 7, jy2, 1.234d0, d, klo, khi, ierr)
      write(*,'(a,1pe24.15,3i3)') "splintd2 at 1.234 =", d, klo, khi, ierr

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
      write(*,'(a)') "# simeqc: 3 unknowns, two right-hand sides; then singular"
      sm(1:15) = [4.0d0, 1.0d0, 2.0d0, 1.0d0, 5.0d0, 1.0d0, 2.0d0, 1.0d0, 6.0d0, &
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

contains

subroutine show(label, values)
      character(len=*), intent(in) :: label
      double precision, intent(in) :: values(:)
      write(*,'(a)') label
      write(*,'(3(1pe24.15))') values
end subroutine show

! ysplin on its own: the coefficients of the cubic on each interval
subroutine ysplin_check()
      double precision :: c(4,100)
      c = 0.0d0
      c(1,1:8) = jy(1:8)
      call ysplin(jx(1:8), c, 8)
      call show("ysplin c(2:4,1:7)", reshape(c(2:4,1:7), [21]))
end subroutine ysplin_check

! osplin on jx/jy(1:n) at the given nondecreasing points
subroutine osplin_case(n, points)
      integer, intent(in) :: n
      double precision, intent(in) :: points(:)
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
      double precision :: tz(6), ty(2,6), res(2)
      write(*,'(a)') "# lir"
      tz = [1.0d0, 2.0d0, 3.5d0, 4.0d0, 6.0d0, 7.5d0]
      ty(1,:) = tz**2
      ty(2,:) = 10.0d0 - tz
      pts(1:5) = [2.7d0, 4.0d0, 0.5d0, 8.0d0, 1.2d0]
      do i = 1, 5
         v = pts(i)
         res = -1.0d0
         call lir(v, tz, res, ty, 2, 2, 6, 1, flag)
         write(*,'(a,f6.2,a,2(1pe24.15),i3)') "lir cubic at", v, " =", res, flag
      end do
      res = -1.0d0
      call lir(2.7d0, tz, res, ty, 2, 2, 3, 1, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir linear (3 points)", res, flag
      res = -1.0d0
      call lir(2.7d0, tz(6:1:-1), res, ty(:,6:1:-1), 2, 2, 6, 1, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir decreasing table", res, flag
      res = -1.0d0
      flag = -7
      call lir(2.7d0, tz, res, ty, 2, 2, 1, 1, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir one point (no change)", res, flag
      res = -1.0d0
      call lir(2.7d0, [3.0d0, 3.0d0, 3.0d0, 3.0d0], res, ty, 2, 2, 4, 1, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir flat table (no change)", res, flag
      res = -1.0d0
      ty(1,3) = 0.0d0
      call lir(3.5d0, tz, res, ty, 2, 2, 6, 1, flag)
      write(*,'(a,2(1pe24.15),i3)') "lir zero on a node", res, flag
end subroutine lir_check

! tridiag_gs and ctridi on one system (ctridi twice: no carried
! state), and a singular system
subroutine tridiagonal_check()
      integer, parameter :: n = 6
      write(*,'(a)') "# tridiag_gs/ctridi"
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
      call tridiag_gs(ja, jb, jc, jr, n, jsol, ierr)
      call show("tridiag_gs", jsol(1:n))
      call ctridi(n, ja, jb, jc, jr, jsol, ierr)
      call show("ctridi", jsol(1:n))
      call ctridi(n - 2, ja, jb, jc, jr, jsol, ierr)
      call show("ctridi n-2", jsol(1:n-2))
      jb(2) = ja(2)*jc(1)/jb(1)
      call tridiag_gs(ja, jb, jc, jr, n, jsol, ierr)
      write(*,'(a,i3)') "tridiag_gs singular: ierr", ierr
      call ctridi(n, ja, jb, jc, jr, jsol, ierr)
      write(*,'(a,i3)') "ctridi singular: ierr", ierr
end subroutine tridiagonal_check

! bsstep: four accepted steps of the harmonic system from x = 0; large
! trial steps accepted at the sixth stage (h grows), the seventh (h
! shrinks) and after a retry; then an integrand failure and a step
! that underflows
subroutine bsstep_check()
      double precision :: y(3), dydx(3), yscale(3), err(3), x, h, h_did, h_next
      double precision :: lum, fp, ft, gl, rl, teffl, xh, z
      logical :: latmo, lderiv, locond
      integer :: ksaha, step
      write(*,'(a)') "# bsstep (mmid, ratext)"
      max_stage_index = 11
      extrap_order = 7
      lum = 0.3d0
      fp = 1.0d0
      ft = 1.0d0
      gl = 4.4d0
      rl = 10.0d0
      teffl = 3.76d0
      xh = 0.7d0
      z = 0.02d0
      latmo = .false.
      lderiv = .false.
      locond = .false.
      ksaha = 0
      y = [0.0d0, 1.0d0, 1.0d0]
      x = 0.0d0
      h = 0.5d0
      yscale = 1.0d0
      do step = 1, 4
         call harmonic_deriv(x, y, dydx, lum, fp, ft, gl, latmo, lderiv, locond, &
              rl, teffl, xh, z, ksaha, ierr)
         call bsstep(y, dydx, 3, x, h, 1.0d-10, yscale, h_did, h_next, harmonic_deriv, &
              lum, fp, ft, gl, latmo, lderiv, locond, rl, teffl, xh, z, ksaha, err, ierr)
         write(*,'(a,i2,a,i3,a,i6)') "step", step, " ierr", ierr, " deriv calls", ksaha
         call show("  y, x, h_did, h_next", [y, x, h_did, h_next])
         call show("  step_err", err)
         h = h_next
      end do
      pts(1:6) = [1.0d0, 1.5d0, 2.0d0, 2.5d0, 3.0d0, 4.0d0]
      do step = 1, 6
         h = pts(step)
         call bsstep(y, dydx, 3, x, h, 1.0d-13, yscale, h_did, h_next, harmonic_deriv, &
              lum, fp, ft, gl, latmo, lderiv, locond, rl, teffl, xh, z, ksaha, err, ierr)
         write(*,'(a,f5.2,a,i3)') "large step, h_try", h, " ierr", ierr
         call show("  y, x, h_did, h_next", [y, x, h_did, h_next])
      end do
      locond = .true.
      call bsstep(y, dydx, 3, x, h, 1.0d-10, yscale, h_did, h_next, harmonic_deriv, &
           lum, fp, ft, gl, latmo, lderiv, locond, rl, teffl, xh, z, ksaha, err, ierr)
      write(*,'(a,i3)') "bsstep, integrand failure: ierr", ierr
      locond = .false.
! zero tolerance: no step is accepted, so h shrinks until it vanishes
! against the hydrogen fraction (the inherited underflow test)
      call bsstep(y, dydx, 3, x, 0.5d0, 0.0d0, yscale, h_did, h_next, harmonic_deriv, &
           lum, fp, ft, gl, latmo, lderiv, locond, rl, teffl, xh, z, ksaha, err, ierr)
      write(*,'(a,i3,a,i8)') "bsstep, step underflow: ierr", ierr, " deriv calls", ksaha
end subroutine bsstep_check

end program test_numerics

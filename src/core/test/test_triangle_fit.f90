!----------------------------------------------------------------------
! test_triangle_fit
!----------------------------------------------------------------------
! The envelope triangle's kernels (core/structure/triangle_fit.f90) on
! hand-made triangles; no star, no tables. Sections:
!   fit         fit_triangle_coefficients, with the fitting-point log P,
!               log T and log R linear functions of the vertex's
!               (log Teff, log L): the three planes (log R and log L,
!               and the reverse relation log Teff) reproduce them at the
!               vertices and at interior points to round-off (the
!               residual over |c1 log P| + |c2 log T| + |c3|, tolerance
!               1e-12), and their slopes match the inverse of the linear
!               map (relative, 1e-10). A plane also passes through
!               arbitrary values at the vertices.
!   degenerate  the fit when the vertices' (log P, log T) are collinear,
!               and when T_2 = T_1 for a triangle that is not: nothing
!               checks the denominators, and all nine coefficients come
!               back infinite or NaN. A zero-area triangle in
!               (log Teff, log L) contains the points on its line, and
!               no reflection brings a point off the line inside; nor
!               any point inside a triangle whose orientation is the
!               opposite of its vertices' (a negative tri_delta_teffl):
!               there reflect_until_inside would not return.
!   inside      vertex_to_reflect for points inside, beyond each edge,
!               beyond two edges, on each edge, at each vertex, and
!               just beyond an edge (with and without a tolerance), for
!               a counterclockwise triangle and the same one clockwise.
!               The value printed is the smallest edge_side.
!   reflection  reflect_vertex is the half turn about the opposite
!               edge's midpoint: the other vertices stay, the old and
!               new vertex straddle the midpoint, the area is kept and
!               the orientation reverses, and a point just beyond the
!               edge is inside after the reflection. It is the mirror
!               image across the edge only when the vertex projects
!               onto the edge's midpoint (the start triangle's vertex 3).
!   walk        reflect_until_inside from the start triangle to near and
!               far points, and along a track of 40 models: it returns
!               with the point inside, the orientation matching the
!               vertices, the moved vertices marked, and the area kept
!               (relative 1e-12). The value printed is the number of
!               reflections. Each walk is first done step by step and
!               cut off (bounded_walk), so a kernel whose walk does not
!               end fails the check instead of hanging the test.
! The triangles of the inside and reflection sections have dyadic
! coordinates, so their edge tests are exact and points on an edge give
! exactly zero. Each check prints its value (1pe9.1) and ok or FAIL;
! any FAIL stops the program with status 1.
program test_triangle_fit
      use, intrinsic :: iso_fortran_env, only: dp => real64, output_unit
      use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
      use, intrinsic :: ieee_exceptions, only: ieee_set_flag, ieee_all
      use triangle_fit, only: start_triangle, edge_side, &
           vertex_to_reflect, reflect_vertex, reflect_until_inside, &
           fit_triangle_coefficients
      implicit none

! The default triangle size (tri_delta_teffl, tri_delta_logl) and the
! reference log Teff of the linear maps.
      real(dp), parameter :: d_teff = 0.01d0, d_logl = 0.08d0
      real(dp), parameter :: teff_ref = 3.76d0
      real(dp), parameter :: fit_tolerance = 1.0d-12, &
           slope_tolerance = 1.0d-10, area_tolerance = 1.0d-12
      integer :: nfail

      nfail = 0
      call fit_section()
      call degenerate_section()
      call inside_section()
      call reflection_section()
      call walk_section()
      if (nfail > 0) then
         write(*,'(a,i4,a)') 'test_triangle_fit: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_triangle_fit: done'

contains

! ---------------------------------------------------------------
! One check: its label, its value and ok or FAIL, flushed at once (a
! walk that does not end leaves the checks before it on record). value
! + 0 prints a zero as +0.
subroutine check(label, value, ok)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: value
      logical, intent(in) :: ok
      if (ok) then
         write(*,'(2x,a,t68,1pe9.1,2x,a)') label, value + 0.0d0, 'ok'
      else
         write(*,'(2x,a,t68,1pe9.1,2x,a)') label, value + 0.0d0, 'FAIL'
         nfail = nfail + 1
      end if
      flush(output_unit)
end subroutine check

! Twice the signed area of the triangle: positive when 1 -> 2 -> 3 runs
! counterclockwise.
pure function twice_area(vt, vl) result(a2)
      real(dp), intent(in) :: vt(3), vl(3)
      real(dp) :: a2
      a2 = (vt(2) - vt(1))*(vl(3) - vl(1)) - (vt(3) - vt(1))*(vl(2) - vl(1))
end function twice_area

! The smallest edge_side of a point.
pure function min_side(vt, vl, orientation, x, y) result(s)
      real(dp), intent(in) :: vt(3), vl(3), orientation, x, y
      real(dp) :: s
      s = min(edge_side(vt, vl, orientation, 1, x, y), &
           edge_side(vt, vl, orientation, 2, x, y), &
           edge_side(vt, vl, orientation, 3, x, y))
end function min_side

! A linear function of (log Teff, log L): a(1) + a(2) (log Teff -
! teff_ref) + a(3) log L.
pure function linear(a, x, y) result(f)
      real(dp), intent(in) :: a(3), x, y
      real(dp) :: f
      f = a(1) + a(2)*(x - teff_ref) + a(3)*y
end function linear

! ===============================================================
! fit
! ===============================================================
subroutine fit_section()
      real(dp) :: vt(3), vl(3), orientation
      real(dp), parameter :: &
           p_a(3) = [13.52d0, -3.1d0, 0.42d0], &
           t_a(3) = [6.08d0, 1.7d0, 0.31d0], &
           r_a(3) = [10.84d0, -2.0d0, 0.5d0], &
           p_b(3) = [11.9d0, -5.0d0, 0.8d0], &
           t_b(3) = [5.6d0, 2.2d0, 0.15d0], &
           r_b(3) = [11.95d0, -1.6d0, 0.47d0], &
           r_c(3) = [10.9d0, 0.7d0, -0.25d0]
      real(dp) :: p(3), t(3), r(3), c(9), resid
      integer :: k

      write(*,'(a)') '# fit: planes through linear functions of (log Teff, log L)'
! A: the start triangle about the Sun's (log Teff, log L).
      call start_triangle(3.7617d0, 0.0d0, d_teff, d_logl, vt, vl, orientation)
      call fit_case('A (start triangle, Sun)', vt, vl, p_a, t_a, r_a)
! B: about a red giant, then vertex 1 reflected, so clockwise.
      call start_triangle(3.65d0, 2.1d0, d_teff, d_logl, vt, vl, orientation)
      call reflect_vertex(vt, vl, orientation, 1)
      call fit_case('B (giant, clockwise)', vt, vl, p_b, t_b, r_b)
! C: a scalene triangle, with another linear log R.
      vt = [3.70d0, 3.715d0, 3.705d0]
      vl = [0.50d0, 0.47d0, 0.56d0]
      call fit_case('C (scalene)', vt, vl, p_a, t_a, r_c)
! The log R plane through arbitrary values at the vertices of C.
      do k = 1,3
         p(k) = linear(p_a, vt(k), vl(k))
         t(k) = linear(t_a, vt(k), vl(k))
      end do
      r = [10.91d0, 10.87d0, 10.95d0]
      call fit_triangle_coefficients(p, t, r, vl, vt, c)
      resid = 0.0d0
      do k = 1,3
         resid = max(resid, abs(c(1)*p(k) + c(2)*t(k) + c(3) - r(k)) &
              /(abs(c(1)*p(k)) + abs(c(2)*t(k)) + abs(c(3))))
      end do
      call check('C: log R plane through arbitrary values at the vertices', &
           resid, resid <= fit_tolerance)
end subroutine fit_section

! The fitting-point log P, log T and log R as the linear functions pa,
! ta, ra at the vertices; fit; then the residuals of the three planes
! at the vertices and at four interior points, and the slopes against
! the inverse of the map (log Teff, log L) -> (log P, log T).
subroutine fit_case(name, vt, vl, pa, ta, ra)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: vt(3), vl(3), pa(3), ta(3), ra(3)
      real(dp), parameter :: weights(3,4) = reshape([ &
           1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, &
           0.6d0, 0.3d0, 0.1d0, &
           0.1d0, 0.2d0, 0.7d0, &
           0.45d0, 0.45d0, 0.1d0], [3,4])
      real(dp) :: p(3), t(3), r(3), c(9), x, y, px, tx
      real(dp) :: res_vtx(3), res_int(3), exact(9), det, slope_err
      integer :: k, m, j

      do k = 1,3
         p(k) = linear(pa, vt(k), vl(k))
         t(k) = linear(ta, vt(k), vl(k))
         r(k) = linear(ra, vt(k), vl(k))
      end do
      call fit_triangle_coefficients(p, t, r, vl, vt, c)
      res_vtx = 0.0d0
      do k = 1,3
         call add_residuals(c, p(k), t(k), r(k), vl(k), vt(k), res_vtx)
      end do
      res_int = 0.0d0
      do m = 1,4
         x = sum(weights(:,m)*vt)
         y = sum(weights(:,m)*vl)
         px = linear(pa, x, y)
         tx = linear(ta, x, y)
         call add_residuals(c, px, tx, linear(ra, x, y), y, x, res_int)
      end do
      call check(name//': log R plane at the vertices', res_vtx(1), &
           res_vtx(1) <= fit_tolerance)
      call check(name//': log R plane at interior points', res_int(1), &
           res_int(1) <= fit_tolerance)
      call check(name//': log L plane at the vertices', res_vtx(2), &
           res_vtx(2) <= fit_tolerance)
      call check(name//': log L plane at interior points', res_int(2), &
           res_int(2) <= fit_tolerance)
      call check(name//': log Teff plane at the vertices', res_vtx(3), &
           res_vtx(3) <= fit_tolerance)
      call check(name//': log Teff plane at interior points', res_int(3), &
           res_int(3) <= fit_tolerance)
! d(log Teff, log L)/d(log P, log T) is the inverse of the map's
! Jacobian [[pa(2), pa(3)], [ta(2), ta(3)]]; log R's slopes follow by
! the chain rule.
      det = pa(2)*ta(3) - pa(3)*ta(2)
      exact = 0.0d0
      exact(7) = ta(3)/det
      exact(8) = -pa(3)/det
      exact(4) = -ta(2)/det
      exact(5) = pa(2)/det
      exact(1) = ra(2)*exact(7) + ra(3)*exact(4)
      exact(2) = ra(2)*exact(8) + ra(3)*exact(5)
      slope_err = 0.0d0
      do j = 1,9
         if (mod(j,3) == 0) cycle
         slope_err = max(slope_err, abs(c(j) - exact(j))/abs(exact(j)))
      end do
      call check(name//': slopes against the inverse map', slope_err, &
           slope_err <= slope_tolerance)
end subroutine fit_case

! The residuals of the three planes at (log P, log T) = (px, tx), where
! log R, log L and log Teff are rx, lx and ex, over the size of the
! plane's terms; res keeps the largest of each.
subroutine add_residuals(c, px, tx, rx, lx, ex, res)
      real(dp), intent(in) :: c(9), px, tx, rx, lx, ex
      real(dp), intent(inout) :: res(3)
      real(dp) :: f(3)
      integer :: i, i0
      f = [rx, lx, ex]
      do i = 1,3
         i0 = 3*(i-1)
         res(i) = max(res(i), abs(c(i0+1)*px + c(i0+2)*tx + c(i0+3) - f(i)) &
              /(abs(c(i0+1)*px) + abs(c(i0+2)*tx) + abs(c(i0+3))))
      end do
end subroutine add_residuals

! ===============================================================
! degenerate
! ===============================================================
subroutine degenerate_section()
      real(dp) :: vt(3), vl(3), orientation, c(9), p(3), t(3), r(3), s
      integer :: n_steps
      logical :: inside

      write(*,'(a)') '# degenerate: collinear fitting points, a zero-area triangle'
      call start_triangle(3.7617d0, 0.0d0, d_teff, d_logl, vt, vl, orientation)
      r = [10.80d0, 10.82d0, 10.81d0]
! Collinear: log T = 6 + (log P - 13.5)/2 at all three, exactly.
      p = [13.5d0, 13.75d0, 14.0d0]
      t = [6.0d0, 6.125d0, 6.25d0]
      call fit_triangle_coefficients(p, t, r, vl, vt, c)
      call check('collinear (log P, log T): non-finite coefficients', &
           real(count(.not.ieee_is_finite(c)), dp), &
           count(.not.ieee_is_finite(c)) == 9)
! Not collinear, but T_2 = T_1: the plane exists, the elimination's
! first pivot is zero.
      p = [13.5d0, 13.75d0, 13.5d0]
      t = [6.0d0, 6.0d0, 6.25d0]
      call fit_triangle_coefficients(p, t, r, vl, vt, c)
      call check('T_2 = T_1, plane exists: non-finite coefficients', &
           real(count(.not.ieee_is_finite(c)), dp), &
           count(.not.ieee_is_finite(c)) == 9)
      call ieee_set_flag(ieee_all, .false.)
! A zero-area triangle: tri_delta_logl = 0 puts the three vertices on
! the line log L = 0.
      call start_triangle(3.7617d0, 0.0d0, d_teff, 0.0d0, vt, vl, orientation)
      s = min_side(vt, vl, orientation, 3.7600d0, 0.0d0)
      call check('zero area: a point on its line is inside', s, &
           vertex_to_reflect(vt, vl, orientation, 3.7600d0, 0.0d0, 0.0d0) == 0)
      s = min_side(vt, vl, orientation, 3.9d0, 0.0d0)
      call check('zero area: so is one on the line far outside', s, &
           vertex_to_reflect(vt, vl, orientation, 3.9d0, 0.0d0, 0.0d0) == 0)
! The walk of reflect_until_inside, cut off after 1000 reflections.
      call bounded_walk(vt, vl, orientation, 3.7617d0, 0.01d0, 1000, &
           inside, n_steps)
      call check('zero area: a point off its line, reflections tried', &
           real(n_steps, dp), .not.inside)
! A negative tri_delta_teffl: the start triangle runs clockwise but its
! orientation is +1, so the three sides sum to minus twice the area and
! one is always negative: not even the triangle's own centre is inside.
      call start_triangle(3.7617d0, 0.0d0, -d_teff, d_logl, vt, vl, orientation)
      call bounded_walk(vt, vl, orientation, 3.7617d0, -0.01d0, 1000, &
           inside, n_steps)
      call check('orientation against the vertices: reflections tried', &
           real(n_steps, dp), .not.inside)
end subroutine degenerate_section

! ===============================================================
! inside
! ===============================================================
subroutine inside_section()
      real(dp), parameter :: tiny_step = 2.0d0**(-40)
      real(dp) :: vt(3), vl(3), wt(3), wl(3)
      real(dp) :: px(13), py(13)
      integer :: expect_ccw(13), expect_cw(13), j
      character(len=40) :: what(13)
      character(len=1) :: digit

! Counterclockwise, dyadic: 1 lower left, 2 lower right, 3 top.
      vt = [3.75d0, 3.765625d0, 3.7578125d0]
      vl = [-0.0625d0, -0.0625d0, 0.0625d0]
! The same triangle clockwise: (1, 3, 2).
      wt = [vt(1), vt(3), vt(2)]
      wl = [vl(1), vl(3), vl(2)]
      what(1) = 'inside';                 px(1) = 3.7578125d0;  py(1) = 0.0d0
      what(2) = 'beyond edge 2-3';        px(2) = 3.765625d0;   py(2) = 0.0d0
      what(3) = 'beyond edge 3-1';        px(3) = 3.75d0;       py(3) = 0.0d0
      what(4) = 'beyond edge 1-2';        px(4) = 3.7578125d0;  py(4) = -0.125d0
      what(5) = 'beyond edges 2-3 and 1-2'; px(5) = 3.78125d0;  py(5) = -0.125d0
      what(6) = 'on edge 1-2';            px(6) = 3.7578125d0;  py(6) = -0.0625d0
      what(7) = 'on edge 2-3';            px(7) = 3.76171875d0; py(7) = 0.0d0
      what(8) = 'on edge 3-1';            px(8) = 3.75390625d0; py(8) = 0.0d0
      what(9) = 'at vertex 1';            px(9) = vt(1);        py(9) = vl(1)
      what(10) = 'at vertex 2';           px(10) = vt(2);       py(10) = vl(2)
      what(11) = 'at vertex 3';           px(11) = vt(3);       py(11) = vl(3)
      what(12) = '2^-40 beyond edge 1-2'; px(12) = 3.7578125d0; py(12) = -0.0625d0 - tiny_step
      what(13) = '(same, tolerance 2^-45)'; px(13) = px(12);    py(13) = py(12)
! The vertex vertex_to_reflect returns (0: inside), in the clockwise
! triangle's labels for cw (vertex 2 of ccw is vertex 3 of cw).
      expect_ccw = [0, 1, 2, 3, 1, 0, 0, 0, 0, 0, 0, 3, 0]
      expect_cw =  [0, 1, 3, 2, 1, 0, 0, 0, 0, 0, 0, 2, 0]

      write(*,'(a)') '# inside: vertex_to_reflect (0 inside), value the smallest edge_side'
      do j = 1,13
         write(digit,'(i1)') expect_ccw(j)
         call one_point('ccw, '//trim(what(j))//' -> '//digit, vt, vl, &
              +1.0d0, px(j), py(j), merge(2.0d0**(-45), 0.0d0, j == 13), &
              expect_ccw(j))
      end do
      do j = 1,13
         write(digit,'(i1)') expect_cw(j)
         call one_point('cw, '//trim(what(j))//' -> '//digit, wt, wl, &
              -1.0d0, px(j), py(j), merge(2.0d0**(-45), 0.0d0, j == 13), &
              expect_cw(j))
      end do
end subroutine inside_section

subroutine one_point(label, vt, vl, orientation, x, y, tolerance, expected)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: vt(3), vl(3), orientation, x, y, tolerance
      integer, intent(in) :: expected
      call check(label, min_side(vt, vl, orientation, x, y), &
           vertex_to_reflect(vt, vl, orientation, x, y, tolerance) == expected)
end subroutine one_point

! ===============================================================
! reflection
! ===============================================================
subroutine reflection_section()
      real(dp) :: vt(3), vl(3), ut(3), ul(3), orientation, o
      real(dp) :: mt, ml, qt, ql, dev, a_old, a_new, s
      real(dp) :: et, el, len2, proj
      integer :: i1, i2, i3, before, after
      character(len=1) :: digit

      write(*,'(a)') '# reflection: reflect_vertex, the half turn about the opposite edge''s midpoint'
      vt = [3.75d0, 3.765625d0, 3.7578125d0]
      vl = [-0.0625d0, -0.0625d0, 0.0625d0]
      orientation = +1.0d0
      do i1 = 1,3
         write(digit,'(i1)') i1
         i2 = mod(i1,3) + 1
         i3 = mod(i2,3) + 1
         ut = vt
         ul = vl
         o = orientation
         call reflect_vertex(ut, ul, o, i1)
! The edge's midpoint is halfway between the old and new vertex; the
! other two vertices do not move.
         mt = 0.5d0*(vt(i2) + vt(i3))
         ml = 0.5d0*(vl(i2) + vl(i3))
         dev = max(abs(0.5d0*(vt(i1) + ut(i1)) - mt), abs(0.5d0*(vl(i1) + ul(i1)) - ml), &
              abs(ut(i2) - vt(i2)), abs(ut(i3) - vt(i3)), &
              abs(ul(i2) - vl(i2)), abs(ul(i3) - vl(i3)))
         call check('vertex '//digit//': half turn about the edge midpoint', &
              dev, dev == 0.0d0)
         a_old = twice_area(vt, vl)
         a_new = twice_area(ut, ul)
         call check('vertex '//digit//': area kept, orientation reversed', &
              abs(a_new + a_old)/abs(a_old), abs(a_new + a_old) <= area_tolerance*abs(a_old) &
              .and. o == -orientation .and. sign(1.0d0, a_new) == o)
! A point a quarter of the way from the edge's midpoint to the new
! vertex: beyond the edge before, inside after.
         qt = mt + 0.25d0*(ut(i1) - mt)
         ql = ml + 0.25d0*(ul(i1) - ml)
         before = vertex_to_reflect(vt, vl, orientation, qt, ql, 0.0d0)
         after = vertex_to_reflect(ut, ul, o, qt, ql, 0.0d0)
         s = min_side(ut, ul, o, qt, ql)
         call check('vertex '//digit//': a point beyond the edge is inside after', &
              s, before == i1 .and. after == 0)
      end do
! The half turn against the mirror image across the edge's line (in
! the (log Teff, log L) plane as drawn): the same for the start
! triangle's vertex 3, which lies above the midpoint of edge 1-2, not
! for its vertex 1.
      call start_triangle(3.7617d0, 0.0d0, d_teff, d_logl, vt, vl, orientation)
      do i1 = 3,1,-2
         write(digit,'(i1)') i1
         i2 = mod(i1,3) + 1
         i3 = mod(i2,3) + 1
         et = vt(i3) - vt(i2)
         el = vl(i3) - vl(i2)
         len2 = et*et + el*el
         proj = ((vt(i1) - vt(i2))*et + (vl(i1) - vl(i2))*el)/len2
! The mirror image of vertex i1: twice its foot on the edge's line,
! minus itself.
         mt = 2.0d0*(vt(i2) + proj*et) - vt(i1)
         ml = 2.0d0*(vl(i2) + proj*el) - vl(i1)
         ut = vt
         ul = vl
         o = orientation
         call reflect_vertex(ut, ul, o, i1)
         dev = sqrt((ut(i1) - mt)**2 + (ul(i1) - ml)**2)
         if (i1 == 3) then
            call check('start triangle, vertex 3: half turn = mirror image', &
                 dev, dev <= 1.0d-14)
         else
            call check('start triangle, vertex 1: half turn /= mirror image', &
                 dev, dev > 1.0d-3)
         end if
      end do
end subroutine reflection_section

! ===============================================================
! walk
! ===============================================================
subroutine walk_section()
      real(dp), parameter :: sun_teff = 3.7617d0
      real(dp) :: targets(2,7), vt(3), vl(3), orientation, a0
      real(dp) :: x, y, area_err
      logical :: reflected(3), all_inside, same
      integer :: n_reflections, j, total
      character(len=48) :: label

      write(*,'(a)') '# walk: reflect_until_inside from the start triangle, value the reflections'
      targets(:,1) = [sun_teff + 0.001d0, 0.0d0]
      targets(:,2) = [sun_teff + 0.008d0, 0.03d0]
      targets(:,3) = [sun_teff - 0.03d0, 0.3d0]
      targets(:,4) = [3.65d0, 2.0d0]
      targets(:,5) = [3.9d0, -1.0d0]
      targets(:,6) = [sun_teff, -3.0d0]
      targets(:,7) = [4.2d0, 3.0d0]
      do j = 1,7
         call start_triangle(sun_teff, 0.0d0, d_teff, d_logl, vt, vl, orientation)
         a0 = twice_area(vt, vl)
         call checked_walk(vt, vl, orientation, targets(1,j), targets(2,j), &
              reflected, n_reflections, same)
         write(label,'(a,f6.3,a,f6.2,a)') 'to (', targets(1,j), ',', targets(2,j), ')'
         call walk_checks(trim(label), vt, vl, orientation, targets(1,j), &
              targets(2,j), a0, reflected, n_reflections, same)
      end do
! A track of 40 models from the Sun's point to (3.66, 1.5), the triangle
! carried from model to model as compute_surface_boundary does.
      call start_triangle(sun_teff, 0.0d0, d_teff, d_logl, vt, vl, orientation)
      a0 = twice_area(vt, vl)
      total = 0
      all_inside = .true.
      do j = 1,40
         x = sun_teff + (3.66d0 - sun_teff)*j/40.0d0
         y = 1.5d0*j/40.0d0
         call checked_walk(vt, vl, orientation, x, y, reflected, &
              n_reflections, same)
         total = total + n_reflections
         all_inside = all_inside .and. same .and. &
              vertex_to_reflect(vt, vl, orientation, x, y, 0.0d0) == 0 .and. &
              sign(1.0d0, twice_area(vt, vl)) == orientation
      end do
      area_err = abs(abs(twice_area(vt, vl))/abs(a0) - 1.0d0)
      call check('track of 40 models: inside at every model', real(total, dp), all_inside)
      call check('track of 40 models: area kept', area_err, area_err <= area_tolerance)
end subroutine walk_section

! The walk of reflect_until_inside step by step, cut off after
! max_steps reflections: ended is true when the triangle contains the
! point, and n_steps counts the reflections.
subroutine bounded_walk(vt, vl, orientation, x, y, max_steps, ended, n_steps)
      real(dp), intent(inout) :: vt(3), vl(3), orientation
      real(dp), intent(in) :: x, y
      integer, intent(in) :: max_steps
      logical, intent(out) :: ended
      integer, intent(out) :: n_steps
      integer :: i1
      ended = .false.
      n_steps = 0
      do
         i1 = vertex_to_reflect(vt, vl, orientation, x, y, 0.0d0)
         if (i1 == 0) then
            ended = .true.
            return
         end if
         if (n_steps == max_steps) return
         call reflect_vertex(vt, vl, orientation, i1)
         n_steps = n_steps + 1
      end do
end subroutine bounded_walk

! The walk to (x, y) with reflect_until_inside, once bounded_walk has
! shown that it ends within 100000 reflections (a broken kernel then
! fails here instead of hanging). same is true when it ended and
! reflect_until_inside ends in the same triangle after as many
! reflections.
subroutine checked_walk(vt, vl, orientation, x, y, reflected, &
     n_reflections, same)
      real(dp), intent(inout) :: vt(3), vl(3), orientation
      real(dp), intent(in) :: x, y
      logical, intent(out) :: reflected(3), same
      integer, intent(out) :: n_reflections
      real(dp) :: bt(3), bl(3), bo
      logical :: ended
      integer :: n_steps
      bt = vt
      bl = vl
      bo = orientation
      call bounded_walk(bt, bl, bo, x, y, 100000, ended, n_steps)
      reflected = .false.
      n_reflections = n_steps
      same = .false.
      if (.not.ended) return
      call reflect_until_inside(vt, vl, orientation, x, y, 0.0d0, &
           reflected, n_reflections)
      same = all(vt == bt) .and. all(vl == bl) .and. orientation == bo &
           .and. n_reflections == n_steps
end subroutine checked_walk

! After a walk to (x, y): it ended as bounded_walk's did (same), the
! point is inside, the orientation matches the vertices, every moved
! vertex is marked reflected, and the area is that of the start
! triangle (a0, twice the signed area).
subroutine walk_checks(label, vt, vl, orientation, x, y, a0, reflected, &
     n_reflections, same)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: vt(3), vl(3), orientation, x, y, a0
      logical, intent(in) :: reflected(3), same
      integer, intent(in) :: n_reflections
      real(dp) :: st(3), sl(3), so, area_err
      logical :: marks_ok
      integer :: k

      call start_triangle(3.7617d0, 0.0d0, d_teff, d_logl, st, sl, so)
      marks_ok = count(reflected) <= n_reflections
      do k = 1,3
         if (vt(k) /= st(k) .or. vl(k) /= sl(k)) marks_ok = marks_ok .and. reflected(k)
      end do
      call check(label//': inside, orientation matches', &
           real(n_reflections, dp), &
           vertex_to_reflect(vt, vl, orientation, x, y, 0.0d0) == 0 .and. &
           sign(1.0d0, twice_area(vt, vl)) == orientation .and. marks_ok &
           .and. same)
      area_err = abs(abs(twice_area(vt, vl))/abs(a0) - 1.0d0)
      call check(label//': area kept', area_err, area_err <= area_tolerance)
end subroutine walk_checks

end program test_triangle_fit

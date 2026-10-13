!----------------------------------------------------------------------
! triangle_fit
!----------------------------------------------------------------------
! The envelope triangle of the surface boundary condition: star-blind
! kernels for compute_surface_boundary (surface_boundary.f90). The
! triangle has three vertices in the (log Teff, log L) plane; at each
! vertex an envelope integration gives log P, log T and log R at the
! fitting point.
!   start_triangle             a new triangle about a point
!   edge_side                  which side of the edge opposite a vertex
!                              a point lies on
!   vertex_to_reflect          the point-in-triangle test
!   reflect_vertex             a vertex through the midpoint of the
!                              opposite edge
!   reflect_until_inside       reflect vertices until the triangle
!                              contains the point
!   fit_triangle_coefficients  log R, log L and log Teff as planes in
!                              the fitting point's (log P, log T)
! They take what they read as arguments, return what they compute in
! arguments, and keep nothing between calls.
!
! The vertices are vtx_log_teff(k), vtx_log_l(k), k = 1, 2, 3, and
! orientation is +1 when 1 -> 2 -> 3 runs counterclockwise with
! log Teff increasing to the right and log L up, -1 when clockwise.
! core/test/test_triangle_fit checks every procedure.
module triangle_fit
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: start_triangle, edge_side, vertex_to_reflect, &
           reflect_vertex, reflect_until_inside, fit_triangle_coefficients

contains

!----------------------------------------------------------------------
! start_triangle
!----------------------------------------------------------------------
! A new triangle, counterclockwise (orientation +1), about the point
! (log_teff, log_l), with vertices (log Teff, log L) at
!   1: (Teff - dT/2, L - dL/2), 2: (Teff + dT/2, L - dL/2),
!   3: (Teff, L + dL/2),
! dT = delta_log_teff and dL = delta_log_l.
pure subroutine start_triangle(log_teff, log_l, delta_log_teff, &
     delta_log_l, vtx_log_teff, vtx_log_l, orientation)
      real(dp), intent(in) :: log_teff, log_l, delta_log_teff, delta_log_l
      real(dp), intent(out) :: vtx_log_teff(3), vtx_log_l(3)
      real(dp), intent(out) :: orientation

      orientation = +1.0d0
      vtx_log_teff(3) = log_teff
      vtx_log_teff(1) = vtx_log_teff(3) - 0.5d0*delta_log_teff
      vtx_log_teff(2) = vtx_log_teff(1) + delta_log_teff
      vtx_log_l(3) = log_l + 0.5d0*delta_log_l
      vtx_log_l(1) = vtx_log_l(3) - delta_log_l
      vtx_log_l(2) = vtx_log_l(1)
end subroutine start_triangle

!----------------------------------------------------------------------
! edge_side
!----------------------------------------------------------------------
! For vertex i1, with i2 and i3 the other two in cyclic order,
!   side = orientation [(Teff_i3 - Teff_i2)(L - L_i2)
!                       - (L_i3 - L_i2)(Teff - Teff_i2)],
! the cross product of the edge i2 -> i3 with the point
! (log_teff, log_l) minus vertex i2: positive when the point lies on
! i1's side of that edge, zero on the edge's line, negative beyond it.
! The three sides of a point sum to twice the triangle's area, so for a
! triangle of zero area one of them is negative unless all are zero.
pure function edge_side(vtx_log_teff, vtx_log_l, orientation, i1, &
     log_teff, log_l) result(side)
      real(dp), intent(in) :: vtx_log_teff(3), vtx_log_l(3), orientation
      integer, intent(in) :: i1
      real(dp), intent(in) :: log_teff, log_l
      real(dp) :: side
      integer :: i2, i3

      i2 = mod(i1,3) + 1
      i3 = mod(i2,3) + 1
      side = orientation*((vtx_log_l(i2)-vtx_log_l(i3))*(log_teff-vtx_log_teff(i2)) + &
             (vtx_log_teff(i3)-vtx_log_teff(i2))*(log_l-vtx_log_l(i2)) )
end function edge_side

!----------------------------------------------------------------------
! vertex_to_reflect
!----------------------------------------------------------------------
! The point-in-triangle test: the first vertex i1, scanning 1, 2, 3,
! whose edge_side is below -tolerance (the point lies beyond the edge
! opposite i1), or 0 when there is none: the point is inside, on an
! edge or at a vertex (or within tolerance of one).
pure function vertex_to_reflect(vtx_log_teff, vtx_log_l, orientation, &
     log_teff, log_l, tolerance) result(i1)
      real(dp), intent(in) :: vtx_log_teff(3), vtx_log_l(3), orientation
      real(dp), intent(in) :: log_teff, log_l, tolerance
      integer :: i1
      integer :: i

      i1 = 0
      do i = 1,3
         if (edge_side(vtx_log_teff, vtx_log_l, orientation, i, &
              log_teff, log_l).lt.-tolerance) then
            i1 = i
            return
         endif
      end do
end function vertex_to_reflect

!----------------------------------------------------------------------
! reflect_vertex
!----------------------------------------------------------------------
! Reflect vertex i1 through the midpoint of the opposite edge,
! x_i1 -> x_i2 + x_i3 - x_i1: the triangle turned half a turn about
! that midpoint. The new triangle shares the edge, lies on its other
! side, and has the same area; with the labels kept, the orientation
! reverses. (It is the mirror image across the edge only when x_i1
! projects onto the edge's midpoint.)
pure subroutine reflect_vertex(vtx_log_teff, vtx_log_l, orientation, i1)
      real(dp), intent(inout) :: vtx_log_teff(3), vtx_log_l(3)
      real(dp), intent(inout) :: orientation
      integer, intent(in) :: i1
      integer :: i2, i3

      i2 = mod(i1,3) + 1
      i3 = mod(i2,3) + 1
      orientation = -orientation
      vtx_log_teff(i1) = vtx_log_teff(i2) + vtx_log_teff(i3) - vtx_log_teff(i1)
      vtx_log_l(i1) = vtx_log_l(i2) + vtx_log_l(i3) - vtx_log_l(i1)
end subroutine reflect_vertex

!----------------------------------------------------------------------
! reflect_until_inside
!----------------------------------------------------------------------
! Reflect vertices until the triangle contains the point
! (log_teff, log_l): while vertex_to_reflect finds a vertex, reflect it
! and scan again from vertex 1. reflected(k) is true for every vertex
! that moved, and n_reflections counts the reflections. The reflections
! tile the plane with copies of the triangle, and for a triangle of
! nonzero area whose orientation matches its vertices the walk ends in
! the copy that contains the point. For a triangle of zero area it never
! ends unless the point lies on the triangle's line (see edge_side).
pure subroutine reflect_until_inside(vtx_log_teff, vtx_log_l, &
     orientation, log_teff, log_l, tolerance, reflected, n_reflections)
      real(dp), intent(inout) :: vtx_log_teff(3), vtx_log_l(3)
      real(dp), intent(inout) :: orientation
      real(dp), intent(in) :: log_teff, log_l, tolerance
      logical, intent(out) :: reflected(3)
      integer, intent(out) :: n_reflections
      integer :: i1

      reflected = .false.
      n_reflections = 0
      do
         i1 = vertex_to_reflect(vtx_log_teff, vtx_log_l, orientation, &
              log_teff, log_l, tolerance)
         if (i1 == 0) exit
         call reflect_vertex(vtx_log_teff, vtx_log_l, orientation, i1)
         reflected(i1) = .true.
         n_reflections = n_reflections + 1
      end do
end subroutine reflect_until_inside

!----------------------------------------------------------------------
! fit_triangle_coefficients
!----------------------------------------------------------------------
! The planes through the three vertices, c = coeffs:
!   log R    = c(1) log P + c(2) log T + c(3)
!   log L    = c(4) log P + c(5) log T + c(6)
!   log Teff = c(7) log P + c(8) log T + c(9)
! with log P, log T and log R at the fitting point (vtx_log_p,
! vtx_log_t, vtx_log_r) and log L, log Teff the vertex's own. For each
! f, the plane f = a log P + b log T + c through the three vertices
! (P_k, T_k: vertex k's fitting-point log P and log T):
!   temp1 = (T_3 - T_1)/(T_2 - T_1)
!   a = [(f_3 - f_1) - temp1 (f_2 - f_1)]/[(P_3 - P_1) - temp1 (P_2 - P_1)]
!   b = [(f_2 - f_1) - a (P_2 - P_1)]/(T_2 - T_1)
!   c = f_1 - a P_1 - b T_1
! Gaussian elimination with T_2 - T_1 as the first pivot, without
! pivoting, and no check of the two denominators: when the vertices'
! (P_k, T_k) are collinear, or T_2 = T_1, a denominator is zero and
! every coefficient comes back infinite or NaN (a build that traps
! floating-point exceptions, make devel, stops there instead).
pure subroutine fit_triangle_coefficients(vtx_log_p, vtx_log_t, &
     vtx_log_r, vtx_log_l, vtx_log_teff, coeffs)
      real(dp), intent(in) :: vtx_log_p(3), vtx_log_t(3), vtx_log_r(3)
      real(dp), intent(in) :: vtx_log_l(3), vtx_log_teff(3)
      real(dp), intent(out) :: coeffs(9)
      real(dp) :: temp, temp1, temp2

      temp = 1.0d0/(vtx_log_t(2)-vtx_log_t(1))
      temp1 = (vtx_log_t(3)-vtx_log_t(1))*temp
      temp2 = 1.0d0/(vtx_log_p(3)-vtx_log_p(1)-temp1*(vtx_log_p(2)-vtx_log_p(1)))
      coeffs(1) = (vtx_log_r(3)-vtx_log_r(1)-temp1*(vtx_log_r(2)-vtx_log_r(1)))*temp2
      coeffs(4) = (vtx_log_l(3)-vtx_log_l(1)-temp1*(vtx_log_l(2)-vtx_log_l(1)))*temp2
      coeffs(7) = (vtx_log_teff(3)-vtx_log_teff(1)-temp1*(vtx_log_teff(2)-vtx_log_teff(1)))*temp2
      coeffs(2) = (vtx_log_r(2)-vtx_log_r(1)-coeffs(1)*(vtx_log_p(2)-vtx_log_p(1)))*temp
      coeffs(5) = (vtx_log_l(2)-vtx_log_l(1)-coeffs(4)*(vtx_log_p(2)-vtx_log_p(1)))*temp
      coeffs(8) = (vtx_log_teff(2)-vtx_log_teff(1)-coeffs(7)*(vtx_log_p(2)-vtx_log_p(1)))*temp
      coeffs(3) = vtx_log_r(1) - coeffs(1)*vtx_log_p(1) - coeffs(2)*vtx_log_t(1)
      coeffs(6) = vtx_log_l(1) - coeffs(4)*vtx_log_p(1) - coeffs(5)*vtx_log_t(1)
      coeffs(9) = vtx_log_teff(1) - coeffs(7)*vtx_log_p(1) - coeffs(8)*vtx_log_t(1)
end subroutine fit_triangle_coefficients

end module triangle_fit

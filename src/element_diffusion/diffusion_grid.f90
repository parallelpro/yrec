!----------------------------------------------------------------------
! diffusion_grid
!----------------------------------------------------------------------
! The equally spaced grids in radius on which both pipelines solve, and
! the interpolation between them and the model's zones. Over the
! settling region r(zone_begin)..r(zone_end) (any radius unit; the
! pipelines use Bahcall & Loeb's):
!   the centres  r_1 = r(zone_begin), r_i = r_(i-1) + dr, i = 2..n,
!                so r_n is r(zone_end) to rounding;
!   the midpoints  r_1 + dr/2, then every dr: n - 1 of them (a second
!                is placed, unused, when n = 2).
! size_equal_grid picks n and dr; place_equal_grids places the points;
! build_equal_grid_map finds, for every point, the model zones it is
! interpolated from (map_field then interpolates any field):
!   midpoint 1   linearly between the two region zones around it;
!   centre 1, n  the model's values at zone_begin and zone_end;
!   the others   4-point Lagrange (numerics_lib's interp and lagrange4)
!                on the four zones around the point (locate_stencil).
! apply_grid_changes maps the species' changes on the centres back to
! the model's zones (both pipelines' map back), with locate_stencil on
! the equal grid as the table, and change_zone applies a zone's change.
module diffusion_grid
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json, i_h1, i_he4, i_metals, i_he3, i_c12, i_o18
      use numerics_lib, only: interp, lagrange4
      implicit none
      private
      public :: equal_grid_map, size_equal_grid, place_equal_grids, &
           build_equal_grid_map, map_field, locate_stencil, apply_grid_changes, &
           change_zone

! the three layers of the model whose zones change by different rules
! (change_zone): the convective core (and zone_begin), the region's
! interior, the surface convection zone (and zone_end)
      integer, parameter :: layer_core = 1, layer_interior = 2, layer_envelope = 3

! Where each point of the two grids takes its values from: midpoint 1
! from zones first_mid_zone-1 and first_mid_zone with the fraction
! first_mid_fraction of the way; centre or midpoint i (2..n-1) from the
! four zones from *_start(i) on, with the Lagrange weights *_weight(:,i).
      type :: equal_grid_map
         integer :: num_points = 0, zone_begin = 0, zone_end = 0
         integer :: first_mid_zone = 0
         real(dp) :: first_mid_fraction = 0.0d0
         integer :: centre_start(json), mid_start(json)
         real(dp) :: centre_weight(4,json), mid_weight(4,json)
      end type equal_grid_map

contains

! The equal grid's point count num_points and spacing grid_spacing over
! radius(zone_begin..zone_end): num_points is the span over the
! narrowest model zone among zones spacing_from..zone_end (each zone k
! spanning radius(k-1)..radius(k)), rounded up, at least 2 and at most
! json, and grid_spacing = span/(num_points - 1). (The F77 meant the
! intervals to be no wider than that zone, "DR <= DRMIN", but counts the
! intervals as points, so the spacing is wider by num_points/
! (num_points - 1): one interval short. Changing it moves results.)
pure subroutine size_equal_grid(radius, zone_begin, zone_end, spacing_from, &
     num_points, grid_spacing)
      real(dp), intent(in) :: radius(:)
      integer, intent(in) :: zone_begin, zone_end, spacing_from
      integer, intent(out) :: num_points
      real(dp), intent(out) :: grid_spacing
      real(dp) :: span, narrowest
      integer :: k

      span = radius(zone_end) - radius(zone_begin)
      narrowest = span
      do k = spacing_from, zone_end
         narrowest = min(narrowest, radius(k) - radius(k-1))
      end do
      num_points = int(span/narrowest)
      if (mod(span, narrowest) /= 0.0d0) num_points = num_points + 1
      if (num_points == 1) num_points = 2
      num_points = min(num_points, json)
      grid_spacing = span/real(num_points-1, dp)
end subroutine size_equal_grid

! The centres and midpoints of the grid of num_points points that starts
! at start_radius with spacing grid_spacing (the module's header).
pure subroutine place_equal_grids(start_radius, grid_spacing, num_points, &
     centre_radius, mid_radius)
      real(dp), intent(in) :: start_radius, grid_spacing
      integer, intent(in) :: num_points
      real(dp), intent(out) :: centre_radius(:), mid_radius(:)
      integer :: i

      mid_radius(1) = start_radius + 0.5d0*grid_spacing
      if (num_points == 2) then
         mid_radius(2) = mid_radius(1) + grid_spacing
      else
         do i = 2, num_points-1
            mid_radius(i) = mid_radius(i-1) + grid_spacing
         end do
      end if
      centre_radius(1) = start_radius
      do i = 2, num_points
         centre_radius(i) = centre_radius(i-1) + grid_spacing
      end do
end subroutine place_equal_grids

! The map from the model's zones (radius(1..num_zones), the region
! zone_begin..zone_end) to the grids' num_points centres and midpoints.
! The zones around midpoint 1 are searched over model zones
! 2..num_points, the grid's point count, as in the F77: when the region
! starts at or above zone num_points (a convective core's top beyond
! it), the search finds nothing and midpoint 1 is extrapolated from
! zones num_points-1 and num_points, inside the core
! (docs/bug_sweep.md M247; searching zone_begin+1..zone_end instead
! moves a pinned deck, so it waits for a measured commit).
pure subroutine build_equal_grid_map(radius, zone_begin, zone_end, num_zones, &
     num_points, centre_radius, mid_radius, map)
      real(dp), intent(in) :: radius(:)
      integer, intent(in) :: zone_begin, zone_end, num_zones, num_points
      real(dp), intent(in) :: centre_radius(:), mid_radius(:)
      type(equal_grid_map), intent(out) :: map
      integer :: iu

      map%num_points = num_points
      map%zone_begin = zone_begin
      map%zone_end = zone_end
      do iu = 2, num_points
         if (radius(iu) >= mid_radius(1)) exit
      end do
      if (iu > num_points) iu = num_points
      map%first_mid_zone = iu
      map%first_mid_fraction = (mid_radius(1) - radius(iu-1))/(radius(iu) - radius(iu-1))
      call find_interior_stencils(radius, zone_begin, zone_end, num_zones, num_points, &
           mid_radius, map%mid_start, map%mid_weight)
      call find_interior_stencils(radius, zone_begin, zone_end, num_zones, num_points, &
           centre_radius, map%centre_start, map%centre_weight)
end subroutine build_equal_grid_map

! The stencils and weights of points 2..num_points-1 of one grid, the
! searches moving out from zone_begin + 1 as the points do.
pure subroutine find_interior_stencils(radius, zone_begin, zone_end, num_zones, &
     num_points, grid_radius, stencil_start, weight)
      real(dp), intent(in) :: radius(:)
      integer, intent(in) :: zone_begin, zone_end, num_zones, num_points
      real(dp), intent(in) :: grid_radius(:)
      integer, intent(inout) :: stencil_start(:)
      real(dp), intent(inout) :: weight(:,:)
      integer :: i, k0, search_from
      logical :: found

      search_from = zone_begin+1
      do i = 2, num_points-1
         call locate_stencil(radius, num_zones, grid_radius(i), search_from, zone_end, &
              k0, found)
         stencil_start(i) = k0
         call interp(radius(k0:k0+3), weight(:,i), x_eval=grid_radius(i))
      end do
end subroutine find_interior_stencils

! The model field's values at the grid points (map): at the centres
! (1..num_points) and at the midpoints (1..num_points-1), either or both.
pure subroutine map_field(map, field, centre_value, mid_value)
      type(equal_grid_map), intent(in) :: map
      real(dp), intent(in) :: field(:)
      real(dp), intent(inout), optional :: centre_value(:), mid_value(:)
      integer :: i, iu, k0

      if (present(mid_value)) then
         iu = map%first_mid_zone
         mid_value(1) = field(iu-1) + map%first_mid_fraction*(field(iu) - field(iu-1))
         do i = 2, map%num_points-1
            k0 = map%mid_start(i)
            mid_value(i) = lagrange4(map%mid_weight(:,i), field(k0:k0+3))
         end do
      end if
      if (present(centre_value)) then
         centre_value(1) = field(map%zone_begin)
         do i = 2, map%num_points-1
            k0 = map%centre_start(i)
            centre_value(i) = lagrange4(map%centre_weight(:,i), field(k0:k0+3))
         end do
         centre_value(map%num_points) = field(map%zone_end)
      end if
end subroutine map_field

! The first of the four table points (table(1..num_table), increasing)
! around target, for 4-point Lagrange interpolation: the first point j
! of search_from..search_to at or above target, less 2, kept within
! 1..num_table-3 (found true); without such a point, num_table - 3
! (found false), which is below 1 for a table of fewer than four
! points (docs/bug_sweep.md O21). search_from moves to j (to num_table
! without one), where the search for the next, larger target starts.
pure subroutine locate_stencil(table, num_table, target, search_from, search_to, &
     k0, found)
      real(dp), intent(in) :: table(:)
      integer, intent(in) :: num_table, search_to
      real(dp), intent(in) :: target
      integer, intent(inout) :: search_from
      integer, intent(out) :: k0
      logical, intent(out) :: found
      integer :: j

      do j = search_from, search_to
         if (table(j) >= target) then
            k0 = min(max(j-2, 1), num_table-3)
            search_from = j
            found = .true.
            return
         end if
      end do
      k0 = num_table - 3
      search_from = num_table
      found = .false.
end subroutine locate_stencil

! ---------------------------------------------------------------
! The changes of X (hydrogen_change), Z (metal_change, with_metals) and
! the light elements (light_change, rows light_rows, with_light) at the
! num_points centres eq_radius, applied to the model's zones (radius):
! every zone up to zone_begin takes the first centre's change, every
! zone from zone_end on (to num_zones) the last centre's, and the zones
! between the change interpolated at their radius (4-point Lagrange on
! locate_stencil's points). A grid of fewer than four points has no
! stencil of its own: with clamp_short_grid (the original routines) the
! stencil starts at point 1; without it (microdiff, the F77's "JVS fix")
! only a three-point grid's does, and a two-point grid's reads outside
! the grid (docs/bug_sweep.md O21). change_zone applies each change.
pure subroutine apply_grid_changes(eq_radius, num_points, radius, zone_begin, zone_end, &
     num_zones, hydrogen_change, metal_change, light_change, light_rows, with_metals, &
     with_light, clamp_short_grid, composition)
      real(dp), intent(in) :: eq_radius(:), radius(:)
      integer, intent(in) :: num_points, zone_begin, zone_end, num_zones
      real(dp), intent(in) :: hydrogen_change(:), metal_change(:), light_change(:,:)
      integer, intent(in) :: light_rows(:)
      logical, intent(in) :: with_metals, with_light, clamp_short_grid
      real(dp), intent(inout) :: composition(:,:)
      real(dp) :: weight(4), dx, dz, dlight(size(light_rows))
      integer :: i, k, k0, search_from
      logical :: found

      do i = 1, zone_begin
         call change_zone(composition(:,i), hydrogen_change(1), metal_change(1), &
              light_change(:,1), light_rows, with_metals, with_light, layer_core)
      end do
      search_from = 2
      do i = zone_begin+1, zone_end-1
         call locate_stencil(eq_radius, num_points, radius(i), search_from, num_points, &
              k0, found)
         if (found .and. clamp_short_grid .and. k0 < 1) k0 = 1
         if (found .and. .not. clamp_short_grid .and. k0 == 0) k0 = 1
         call interp(eq_radius(k0:k0+3), weight, x_eval=radius(i))
         dx = lagrange4(weight, hydrogen_change(k0:k0+3))
         dz = 0.0d0
         if (with_metals) dz = lagrange4(weight, metal_change(k0:k0+3))
         do k = 1, size(light_rows)
            dlight(k) = 0.0d0
            if (with_light) dlight(k) = lagrange4(weight, light_change(k,k0:k0+3))
         end do
         call change_zone(composition(:,i), dx, dz, dlight, light_rows, with_metals, &
              with_light, layer_interior)
      end do
      do i = zone_end, num_zones
         call change_zone(composition(:,i), hydrogen_change(num_points), &
              metal_change(num_points), light_change(:,num_points), light_rows, &
              with_metals, with_light, layer_envelope)
      end do
end subroutine apply_grid_changes

! ---------------------------------------------------------------
! One zone's mass fractions xa(species) after settling changes X by dx,
! Z by dz (with_metals) and the light elements (rows light_rows) by
! dlight (with_light), in the layer given (layer_core, layer_interior,
! layer_envelope). Hydrogen rises and the metals sink, so the guards
! differ by layer (the F77's):
!   X  at least 0 in the core; elsewhere at most 1 - Z - He3 (the Z
!      before its change);
!   Z  at most 1 - X - He3 (the new X), but at least 0 in the envelope;
!      the CNO tracers (i_c12..i_o18) scaled with it;
!   Y  1 - X - Z - He3;
!   the light elements at least 0 in the core, unguarded elsewhere.
pure subroutine change_zone(xa, dx, dz, dlight, light_rows, with_metals, with_light, layer)
      real(dp), intent(inout) :: xa(:)
      real(dp), intent(in) :: dx, dz, dlight(:)
      integer, intent(in) :: light_rows(:), layer
      logical, intent(in) :: with_metals, with_light
      real(dp) :: x_max, z_max, z_new, z_ratio
      integer :: k, row

      if (layer == layer_core) then
         xa(i_h1) = max(xa(i_h1) + dx, 0.0d0)
      else
         x_max = 1.0d0 - xa(i_metals) - xa(i_he3)
         xa(i_h1) = min(xa(i_h1) + dx, x_max)
      end if
      if (with_metals) then
         if (layer == layer_envelope) then
            z_new = max(xa(i_metals) + dz, 0.0d0)
         else
            z_max = 1.0d0 - xa(i_h1) - xa(i_he3)
            z_new = min(xa(i_metals) + dz, z_max)
         end if
         z_ratio = z_new/xa(i_metals)
         xa(i_metals) = z_new
         xa(i_c12:i_o18) = z_ratio*xa(i_c12:i_o18)
      end if
      xa(i_he4) = 1.0d0 - xa(i_h1) - xa(i_metals) - xa(i_he3)
      if (with_light) then
         do k = 1, size(light_rows)
            row = light_rows(k)
            if (layer == layer_core) then
               xa(row) = max(xa(row) + dlight(k), 0.0d0)
            else
               xa(row) = xa(row) + dlight(k)
            end if
         end do
      end if
end subroutine change_zone

end module diffusion_grid

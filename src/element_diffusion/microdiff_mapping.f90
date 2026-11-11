!----------------------------------------------------------------------
! microdiff_mapping
!----------------------------------------------------------------------
! microdiff's grids and its maps between them and the model:
!   map_microdiff_to_equal_grid  the two equally spaced grids over the
!       settling region (centres and midpoints, Bahcall & Loeb radius
!       units; diffusion_grid) and the structure and the mass fractions
!       on them (the microdiff_grid records);
!   microdiff_etm  the species' changes on the centre grid back on the
!       model's zones (change_zone), with Y from X + Y + Z + He3 = 1.
module microdiff_mapping
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use element_diffusion_configuration, only: diffusion_settings
      use composition_layout, only: json, i_h1, i_he4, i_metals, i_he3, i_c12, i_o18
      use species_table, only: num_light_diffused
      implicit none
      private
      public :: microdiff_grid, map_microdiff_to_equal_grid, microdiff_etm

! One equally spaced grid of microdiff (Bahcall & Loeb units for the
! radius, the enclosed mass and d ln P/dr; g/cm^3 and K for the density
! and the temperature): the structure and the mass fractions at its
! points, light holding a row for each light element species_table
! lists (num_light_diffused: Li6, Li7, Be9).
      type :: microdiff_grid
         real(dp) :: mass(json), radius(json), density(json), &
              temperature(json), dlnp_dr(json), del_grad(json), &
              hydrogen(json), helium(json), metal(json), &
              light(num_light_diffused,json)
      end type microdiff_grid

! the three layers of the model that microdiff_etm changes by different
! rules: the convective core (and zone_begin), the region's interior,
! the surface convection zone (and zone_end)
      integer, parameter :: layer_core = 1, layer_interior = 2, layer_envelope = 3

contains

! ---------------------------------------------------------------
! The two grids over the settling region zone_begin..zone_end of a model
! of num_zones zones, and the structure and the mass fractions on them.
! In: the mass fractions composition(species, zone) (with the light
! elements' rows light_element_id(1..num_light)), and d ln P/dr, the
! radius and the enclosed mass in Bahcall & Loeb units, del_grad, the
! density (g/cm^3) and the temperature (K). Out: the grids' point count
! num_eq_points and spacing grid_spacing, from the narrowest of the
! region's outer num_spacing_zones zones (size_equal_grid); the centre
! and midpoint grids.
subroutine map_microdiff_to_equal_grid(num_light, light_element_id, composition, &
     dlnp_dr, radius_bl, enclosed_mass, del_grad, zone_begin, zone_end, &
     num_zones, grid_spacing, num_eq_points, density, temperature, &
     centre_grid, mid_grid)
      use diffusion_grid, only: equal_grid_map, size_equal_grid, place_equal_grids, &
           build_equal_grid_map, map_field
      integer, intent(in) :: num_light
      integer, intent(in) :: light_element_id(:)
      real(dp), intent(in) :: composition(:,:), dlnp_dr(:), radius_bl(:), &
           enclosed_mass(:), del_grad(:), density(:), temperature(:)
      integer, intent(in) :: zone_begin, zone_end, num_zones
      real(dp), intent(out) :: grid_spacing
      integer, intent(out) :: num_eq_points
      type(microdiff_grid), intent(out) :: centre_grid, mid_grid
! the zones each grid point is interpolated from (diffusion_grid)
      type(equal_grid_map) :: map
! the grid's spacing comes from the region's outer zones, below the
! surface convection zone
      integer, parameter :: num_spacing_zones = 20
      integer :: k

      call size_equal_grid(radius_bl, zone_begin, zone_end, &
           max(zone_end-num_spacing_zones, zone_begin+1), num_eq_points, grid_spacing)
      call place_equal_grids(radius_bl(zone_begin), grid_spacing, num_eq_points, &
           centre_grid%radius, mid_grid%radius)
      call build_equal_grid_map(radius_bl, zone_begin, zone_end, num_zones, &
           num_eq_points, centre_grid%radius, mid_grid%radius, map)
      call map_field(map, enclosed_mass, centre_grid%mass, mid_grid%mass)
      call map_field(map, density, centre_grid%density, mid_grid%density)
      call map_field(map, temperature, centre_grid%temperature, mid_grid%temperature)
      call map_field(map, dlnp_dr, centre_grid%dlnp_dr, mid_grid%dlnp_dr)
      call map_field(map, del_grad, centre_grid%del_grad, mid_grid%del_grad)
      call map_field(map, composition(i_h1,:), centre_grid%hydrogen, mid_grid%hydrogen)
      call map_field(map, composition(i_he4,:), centre_grid%helium, mid_grid%helium)
      call map_field(map, composition(i_metals,:), centre_grid%metal, mid_grid%metal)
      if (diffusion_settings%diffuse_lithium) then
         do k = 1, num_light
            call map_field(map, composition(light_element_id(k),:), &
                 centre_grid%light(k,:), mid_grid%light(k,:))
         end do
      end if
end subroutine map_microdiff_to_equal_grid

! ---------------------------------------------------------------
! The changes of X (hydrogen_change), Z (metal_change, with
! use_diffusion_z) and the light elements (light_change, rows
! light_element_id, with diffuse_lithium) at the num_eq_points centres
! eq_radius, applied to the model's zones (radius radius_bl): every zone
! up to zone_begin takes the first centre's change, every zone from
! zone_end on the last centre's, and the zones between take the change
! interpolated at their radius (4-point Lagrange, diffusion_grid's
! locate_stencil). change_zone applies a change.
subroutine microdiff_etm(eq_radius, hydrogen_change, metal_change, light_change, &
     zone_begin, zone_end, num_eq_points, composition, radius_bl, num_zones, &
     num_light, light_element_id)
      use numerics_lib, only: interp, lagrange4
      use diffusion_grid, only: locate_stencil
      real(dp), intent(in) :: eq_radius(:), hydrogen_change(:), metal_change(:)
      integer, intent(in) :: num_light
      real(dp), intent(in) :: light_change(:,:)
      integer, intent(in) :: zone_begin, zone_end, num_eq_points
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(in) :: radius_bl(:)
      integer, intent(in) :: num_zones
      integer, intent(in) :: light_element_id(:)
      real(dp) :: weight(4), dx, dz, dlight(num_light)
      integer :: i, k, k0, search_from
      logical :: found

      do i = 1, zone_begin
         call change_zone(composition(:,i), hydrogen_change(1), metal_change(1), &
              light_change(:,1), light_element_id, layer_core)
      end do
      search_from = 2
      do i = zone_begin+1, zone_end-1
! The four centres around the zone. A grid of three points takes the
! stencil from point 1, the one it has; one of two points does not, and
! reads outside the grid (docs/bug_sweep.md O21).
         call locate_stencil(eq_radius, num_eq_points, radius_bl(i), search_from, &
              num_eq_points, k0, found)
         if (found .and. k0 == 0) k0 = 1
         call interp(eq_radius(k0:k0+3), weight, x_eval=radius_bl(i))
         dx = lagrange4(weight, hydrogen_change(k0:k0+3))
         dz = 0.0d0
         if (diffusion_settings%use_diffusion_z) dz = lagrange4(weight, metal_change(k0:k0+3))
         do k = 1, num_light
            dlight(k) = 0.0d0
            if (diffusion_settings%diffuse_lithium) dlight(k) = lagrange4(weight, &
                 light_change(k,k0:k0+3))
         end do
         call change_zone(composition(:,i), dx, dz, dlight, light_element_id, layer_interior)
      end do
      do i = zone_end, num_zones
         call change_zone(composition(:,i), hydrogen_change(num_eq_points), &
              metal_change(num_eq_points), light_change(:,num_eq_points), &
              light_element_id, layer_envelope)
      end do
end subroutine microdiff_etm

! ---------------------------------------------------------------
! One zone's mass fractions xa(species) after settling changes X by dx,
! Z by dz (with use_diffusion_z) and the light elements (rows
! light_element_id) by dlight (with diffuse_lithium), in the given layer.
! Hydrogen rises and the metals sink, so the guards differ by layer (the
! F77's):
!   X  at least 0 in the core; elsewhere at most 1 - Z - He3 (the Z
!      before its change);
!   Z  at most 1 - X - He3 (the new X), but at least 0 in the envelope;
!      the CNO tracers (i_c12..i_o18) scaled with it;
!   Y  1 - X - Z - He3;
!   the light elements at least 0 in the core, unguarded elsewhere.
pure subroutine change_zone(xa, dx, dz, dlight, light_element_id, layer)
      real(dp), intent(inout) :: xa(:)
      real(dp), intent(in) :: dx, dz, dlight(:)
      integer, intent(in) :: light_element_id(:), layer
      real(dp) :: x_max, z_max, z_new, z_ratio
      integer :: k, row

      if (layer == layer_core) then
         xa(i_h1) = max(xa(i_h1) + dx, 0.0d0)
      else
         x_max = 1.0d0 - xa(i_metals) - xa(i_he3)
         xa(i_h1) = min(xa(i_h1) + dx, x_max)
      end if
      if (diffusion_settings%use_diffusion_z) then
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
      if (diffusion_settings%diffuse_lithium) then
         do k = 1, size(light_element_id)
            row = light_element_id(k)
            if (layer == layer_core) then
               xa(row) = max(xa(row) + dlight(k), 0.0d0)
            else
               xa(row) = xa(row) + dlight(k)
            end if
         end do
      end if
end subroutine change_zone

end module microdiff_mapping

! Sources
! - The equal grids and the maps between them and the model: MHP (3/94,
!   the metals), G. Somers (5/2015, the light elements), JvS (the guards
!   for grids of one to three points).
! - 2026: the grids built by diffusion_grid, shared with the original
!   routines; microdiff_grid's light rows sized by species_table
!   (docs/bug_sweep.md M252).

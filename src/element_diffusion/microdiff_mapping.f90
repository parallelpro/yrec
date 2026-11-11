!----------------------------------------------------------------------
! microdiff_mapping
!----------------------------------------------------------------------
! microdiff's grids and its maps between them and the model:
!   map_microdiff_to_equal_grid  the two equally spaced grids over the
!       settling region (centres and midpoints, Bahcall & Loeb radius
!       units; diffusion_grid) and the structure and the mass fractions
!       on them (the microdiff_grid records);
!   microdiff_etm  the species' changes on the centre grid back on the
!       model's zones (diffusion_grid's apply_grid_changes), with Y from
!       X + Y + Z + He3 = 1.
module microdiff_mapping
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use element_diffusion_configuration, only: diffusion_settings
      use composition_layout, only: json, i_h1, i_he4, i_metals
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
! eq_radius, applied to the model's zones (radius radius_bl) by
! diffusion_grid's apply_grid_changes: the core's zones take the first
! centre's change, the envelope's the last's, the zones between the
! change at their radius; a two-point grid keeps the F77's guard
! (docs/bug_sweep.md O21).
subroutine microdiff_etm(eq_radius, hydrogen_change, metal_change, light_change, &
     zone_begin, zone_end, num_eq_points, composition, radius_bl, num_zones, &
     num_light, light_element_id)
      use diffusion_grid, only: apply_grid_changes
      real(dp), intent(in) :: eq_radius(:), hydrogen_change(:), metal_change(:)
      integer, intent(in) :: num_light
      real(dp), intent(in) :: light_change(:,:)
      integer, intent(in) :: zone_begin, zone_end, num_eq_points
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(in) :: radius_bl(:)
      integer, intent(in) :: num_zones
      integer, intent(in) :: light_element_id(:)

      call apply_grid_changes(eq_radius, num_eq_points, radius_bl, zone_begin, zone_end, &
           num_zones, hydrogen_change, metal_change, light_change(1:num_light,:), &
           light_element_id(1:num_light), diffusion_settings%use_diffusion_z, &
           diffusion_settings%diffuse_lithium, .false., composition)
end subroutine microdiff_etm

end module microdiff_mapping

! Sources
! - The equal grids and the maps between them and the model: MHP (3/94,
!   the metals), G. Somers (5/2015, the light elements), JvS (the guards
!   for grids of one to three points).
! - 2026: the grids built by diffusion_grid, shared with the original
!   routines; microdiff_grid's light rows sized by species_table
!   (docs/bug_sweep.md M252).

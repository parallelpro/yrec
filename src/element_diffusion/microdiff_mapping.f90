!----------------------------------------------------------------------
! microdiff_mapping
!----------------------------------------------------------------------
! microdiff's grids. map_microdiff_to_equal_grid builds the equally
! spaced radial grids (zone centers, eq; midpoints, eq_mid; Bahcall &
! Loeb units) and interpolates the structure and the mass fractions
! onto them (4-point Lagrange, numerics_lib). microdiff_etm carries the
! abundance changes back to the model grid.
module microdiff_mapping
      use element_diffusion_configuration, only: diffusion_settings
      use composition_layout, only: json
      use species_table, only: num_light_diffused
      implicit none

! One equally spaced radial grid of the microdiff settling pipeline
! (Bahcall & Loeb units): structure plus interpolated composition.
! light has a row for each light element species_table lists
! (num_light_diffused: Li6, Li7, Be9).
      type :: microdiff_grid
         double precision :: mass(json), radius(json), density(json), &
              temperature(json), dlnp_dr(json), del_grad(json), &
              hydrogen(json), helium(json), metal(json), &
              light(num_light_diffused,json)
      end type microdiff_grid
contains

subroutine map_microdiff_to_equal_grid(num_light, light_element_id, composition, &
     dlnp_dr, radius_bl, enclosed_mass, del_grad, zone_begin, zone_end, &
     num_zones, grid_spacing, num_eq_points, density_orig, &
     temperature_orig, eq, eq_mid)

      use composition_layout, only: json, i_h1, i_he4, i_metals
      use diffusion_grid, only: equal_grid_map, size_equal_grid, place_equal_grids, &
           build_equal_grid_map, map_field
      implicit none

      integer, intent(in) :: num_light
      integer, intent(in) :: light_element_id(num_light)
      double precision, intent(in) :: composition(15,json), dlnp_dr(json), &
           radius_bl(json), enclosed_mass(json), del_grad(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(out) :: grid_spacing
      integer, intent(out) :: num_eq_points
      double precision, intent(in) :: density_orig(json), temperature_orig(json)
      type(microdiff_grid), intent(out) :: eq, eq_mid
! --- locals ---
! the zones each grid point is interpolated from (diffusion_grid)
      type(equal_grid_map) :: map
! the grid's spacing is the narrowest of the region's outer zones (the
! 20 below the surface convection zone)
      integer, parameter :: num_spacing_zones = 20
      integer :: kk

      call size_equal_grid(radius_bl, zone_begin, zone_end, &
           max(zone_end-num_spacing_zones, zone_begin+1), num_eq_points, grid_spacing)
      call place_equal_grids(radius_bl(zone_begin), grid_spacing, num_eq_points, &
           eq%radius, eq_mid%radius)
      call build_equal_grid_map(radius_bl, zone_begin, zone_end, num_zones, &
           num_eq_points, eq%radius, eq_mid%radius, map)
      call map_field(map, enclosed_mass, eq%mass, eq_mid%mass)
      call map_field(map, density_orig, eq%density, eq_mid%density)
      call map_field(map, temperature_orig, eq%temperature, eq_mid%temperature)
      call map_field(map, dlnp_dr, eq%dlnp_dr, eq_mid%dlnp_dr)
      call map_field(map, del_grad, eq%del_grad, eq_mid%del_grad)
      call map_field(map, composition(i_h1,:), eq%hydrogen, eq_mid%hydrogen)
      call map_field(map, composition(i_he4,:), eq%helium, eq_mid%helium)
      call map_field(map, composition(i_metals,:), eq%metal, eq_mid%metal)
      if (diffusion_settings%diffuse_lithium) then
         do kk = 1, num_light
            call map_field(map, composition(light_element_id(kk),:), eq%light(kk,:), &
                 eq_mid%light(kk,:))
         end do
      end if
      return
end subroutine map_microdiff_to_equal_grid

! ---------------------------------------------------------------
! The reverse of map_microdiff_to_equal_grid: interpolates the changes
! in hydrogen, metal and light-element abundance from the equally spaced
! grid back onto the model grid, applies the floors and ceilings,
! renormalizes helium by X+Y+Z=1, and converts DT/STOT/HS1/HQPR back out
! of Bahcall & Loeb units (radius_bl is not converted back; nothing
! reads it afterwards).
!  TRANSFORM BACK TO ORIGINAL GRID OF MODEL POINTS FROM EQUALLY
!  SPACED GRID.
subroutine microdiff_etm(eq_radius, eq_delta_hydrogen, &
     eq_delta_metal, eq_delta_light, zone_begin, zone_end, num_eq_points, &
     composition, radius_bl, num_zones, num_light, light_element_id)
      use composition_layout, only: json, i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_o18
      use numerics_lib, only: interp, lagrange4
      use diffusion_grid, only: locate_stencil
      implicit none

      double precision, intent(in) :: eq_radius(json), &
           eq_delta_hydrogen(json), eq_delta_metal(json)
      integer, intent(in) :: num_light
      double precision, intent(in) :: eq_delta_light(num_light,json)
      integer, intent(in) :: zone_begin, zone_end, num_eq_points
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: radius_bl(json)
      integer, intent(in) :: num_zones
      integer, intent(in) :: light_element_id(num_light)
! --- locals ---
      double precision :: facinterp(4)
      integer :: i, ii, j, jmin, k0, kk
      logical :: found
      double precision :: x_min_floor, z_max, zz, zz2, radmod, dxmod, &
           x_max, dzmod

!  TRANSFORM BACK TO ORIGINAL GRID OF MODEL POINTS FROM EQUALLY
!  SPACED GRID.
      x_min_floor = 0.0d0
      do i = zone_begin,1,-1
         composition(i_h1,i)=max(composition(i_h1,i) + eq_delta_hydrogen(1),x_min_floor)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
! NOTE THAT BECAUSE METALS SINK, AND HYDROGEN RISES, THE FAILSAFES
! ARE OPPOSITE (GUARDING AGAINST NEGATIVE X AND Z>1 RESPECTIVELY).
      if(diffusion_settings%use_diffusion_z)then
         do i = zone_begin,1,-1
            z_max = 1.0d0 - composition(i_h1,i) - composition(i_he3,i)
            zz=min(composition(i_metals,i)+eq_delta_metal(1),z_max)
            zz2 = zz/composition(i_metals,i)
            composition(i_metals,i) = zz
            do j = i_c12,i_o18
               composition(j,i) = zz2*composition(j,i)
            end do
            composition(i_he4,i)=1.0d0-composition(i_h1,i)-composition(i_metals,i)-composition(i_he3,i)
         end do
      else
         do i = zone_begin,1,-1
            composition(i_he4,i)=1.0d0-composition(i_h1,i)-composition(i_metals,i)-composition(i_he3,i)
         end do
      endif
! G SOMERS 5/15; ADD LIGHT ELEMENT DIFFUSION
      if(diffusion_settings%diffuse_lithium)then
         do kk = 1,num_light
            x_min_floor = 0.0d0
            do i = zone_begin,1,-1
               ii = light_element_id(kk)
               composition(ii,i)=max(composition(ii,i) + eq_delta_light(kk,1),x_min_floor)
            end do
         end do
      endif
!
      jmin=2
      do i=zone_begin+1,zone_end-1
!  THE FOUR EQUALLY SPACED GRID POINTS AROUND THE MODEL POINT
!  (diffusion_grid). A grid of three points gets the stencil from point
!  1 (the one it has); one of two points does not, and reads outside the
!  grid (docs/bug_sweep.md O21).
         call locate_stencil(eq_radius, num_eq_points, radius_bl(i), jmin, &
              num_eq_points, k0, found)
         if (found .and. k0 == 0) k0 = 1
         radmod=radius_bl(i)
!  FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
         call interp(eq_radius(k0:k0+3), facinterp, x_eval=radmod)
!  PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR CHANGE IN X.
         dxmod = lagrange4(facinterp, eq_delta_hydrogen(k0:k0+3))
         x_max = 1.0d0 - composition(i_metals,i) - composition(i_he3,i)
         composition(i_h1,i)=min(composition(i_h1,i) + dxmod,x_max)
! MHP 3/94 ADDED METAL DIFFUSION
         if(diffusion_settings%use_diffusion_z)then
            z_max = 1.0d0 - composition(i_h1,i) - composition(i_he3,i)
            dzmod = lagrange4(facinterp, eq_delta_metal(k0:k0+3))
            zz = min(composition(i_metals,i)+dzmod,z_max)
            zz2 = zz/composition(i_metals,i)
            composition(i_metals,i)=zz
            do j = i_c12,i_o18
               composition(j,i) = zz2*composition(j,i)
            enddo
            composition(i_he4,i)=1.0d0-composition(i_h1,i)-composition(i_metals,i)-composition(i_he3,i)
         else
            composition(i_he4,i)=1.0d0-composition(i_h1,i)-composition(i_metals,i)-composition(i_he3,i)
         endif
! GES 5/15 ADDED LIGHT ELEMENT DIFFUSION
         if(diffusion_settings%diffuse_lithium)then
            do kk = 1,num_light
               ii = light_element_id(kk)
               dxmod = lagrange4(facinterp, eq_delta_light(kk,k0:k0+3))
               composition(ii,i) = composition(ii,i) + dxmod
            enddo
         endif
      end do
!
      do i = zone_end,num_zones
         x_max = 1.0d0 - composition(i_metals,i) - composition(i_he3,i)
         composition(i_h1,i)=min(composition(i_h1,i) + eq_delta_hydrogen(num_eq_points),x_max)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
      if(diffusion_settings%use_diffusion_z)then
         do i = zone_end,num_zones
            zz = max(composition(i_metals,i)+eq_delta_metal(num_eq_points),0.0d0)
            zz2 = zz/composition(i_metals,i)
            composition(i_metals,i) = zz
            do j = i_c12,i_o18
               composition(j,i) = zz2*composition(j,i)
            end do
            composition(i_he4,i)=1.0d0-composition(i_h1,i)-composition(i_metals,i)-composition(i_he3,i)
         end do
      else
         do i = zone_end,num_zones
            composition(i_he4,i)=1.0d0-composition(i_h1,i)-composition(i_metals,i)-composition(i_he3,i)
         end do
      endif
! GES 5/15 LIGHT ELEMENT DIFFUSION
      if(diffusion_settings%diffuse_lithium)then
         do kk = 1,num_light
            do i = zone_end,num_zones
               ii = light_element_id(kk)
               composition(ii,i)=composition(ii,i) + eq_delta_light(kk,num_eq_points)
            end do
         enddo
      endif
! (The timestep, dlnP/dr, mass coordinate and total mass are microdiff's
! copies in Bahcall and Loeb units, which it discards: nothing to
! convert back; docs/bug_sweep.md O16.)
      return
end subroutine microdiff_etm

end module microdiff_mapping

!----------------------------------------------------------------------
! settling_grid
!----------------------------------------------------------------------
! The original settling routines' grid mapping: the model onto an equally
! spaced grid in radius and back.
module settling_grid
      use element_diffusion_configuration, only: diffusion_settings
      use element_diffusion_scratch, only: diffusion_scratch
      implicit none
      private
      public :: model_to_equal, equal_to_model

contains

! ---------------------------------------------------------------
! MHP 3/94 added metal diffusion.
subroutine model_to_equal(diff_scr, diffusion_coeff1, diffusion_coeff2, composition, &
     radius, enclosed_mass, &
     zone_begin, zone_end, num_zones, diffusion_coeff1_dx, &
     diffusion_coeff2_dx, grid_spacing, &
     equal_diffusion_coeff1, equal_diffusion_coeff1_mid, &
     equal_diffusion_coeff2_mid, equal_mass, equal_mass_mid, &
     equal_diffusion_coeff1_dx_mid, equal_diffusion_coeff2_dx_mid, &
     equal_radius, equal_hydrogen_fraction, equal_hydrogen_fraction_mid, &
     num_equal_points)
      use composition_layout, only: json, i_h1, i_metals
      use diffusion_grid, only: equal_grid_map, size_equal_grid, place_equal_grids, &
           build_equal_grid_map, map_field
      implicit none
      type(diffusion_scratch), intent(inout) :: diff_scr

      double precision, intent(in) :: diffusion_coeff1(json), &
           diffusion_coeff2(json)
      double precision, intent(in) :: composition(15,json)
      double precision, intent(in) :: radius(json), enclosed_mass(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(in) :: diffusion_coeff1_dx(json), &
           diffusion_coeff2_dx(json)
      double precision, intent(out) :: grid_spacing
      double precision, intent(out) :: equal_diffusion_coeff1(json), &
           equal_diffusion_coeff1_mid(json), equal_diffusion_coeff2_mid(json), &
           equal_mass(json), equal_mass_mid(json), &
           equal_diffusion_coeff1_dx_mid(json), equal_diffusion_coeff2_dx_mid(json), &
           equal_radius(json), equal_hydrogen_fraction(json), &
           equal_hydrogen_fraction_mid(json)
      integer, intent(out) :: num_equal_points




! the zones each grid point is interpolated from (diffusion_grid), and
! the midpoints' radii
      type(equal_grid_map) :: map
      double precision :: mid_radius(json)

! the grid's spacing is the narrowest of the region's zones
      call size_equal_grid(radius, zone_begin, zone_end, zone_begin+1, &
           num_equal_points, grid_spacing)
      call place_equal_grids(radius(zone_begin), grid_spacing, num_equal_points, &
           equal_radius, mid_radius)
      call build_equal_grid_map(radius, zone_begin, zone_end, num_zones, &
           num_equal_points, equal_radius, mid_radius, map)
      call map_field(map, enclosed_mass, equal_mass, equal_mass_mid)
      call map_field(map, diffusion_coeff1, equal_diffusion_coeff1, equal_diffusion_coeff1_mid)
      call map_field(map, diffusion_coeff2, mid_value=equal_diffusion_coeff2_mid)
      call map_field(map, diffusion_coeff1_dx, mid_value=equal_diffusion_coeff1_dx_mid)
      call map_field(map, diffusion_coeff2_dx, mid_value=equal_diffusion_coeff2_dx_mid)
      call map_field(map, composition(i_h1,:), equal_hydrogen_fraction, &
           equal_hydrogen_fraction_mid)
! The metals (as fully ionized iron): their coefficients, their
! derivatives in Z, and Z itself (slot i_metals; the inherited code read
! N15, slot 8, until the 2026 bug sweep).
      if(diffusion_settings%use_diffusion_z)then
         call map_field(map, diff_scr%src_grid_metal_diffusion_coeff1, &
              diff_scr%metal_diffusion_coeff1, diff_scr%metal_diffusion_coeff1_mid)
         call map_field(map, diff_scr%src_grid_metal_diffusion_coeff2, &
              mid_value=diff_scr%metal_diffusion_coeff2_mid)
         call map_field(map, diff_scr%src_grid_metal_diffusion_coeff1_dz, &
              mid_value=diff_scr%eq_metal_diffusion_coeff1_mid)
         call map_field(map, diff_scr%src_grid_metal_diffusion_coeff2_dz, &
              mid_value=diff_scr%eq_metal_diffusion_coeff2_mid)
         call map_field(map, composition(i_metals,:), diff_scr%metal_abundance_change, &
              diff_scr%metal_abundance_change_mid)
      endif
      return
end subroutine model_to_equal

! ---------------------------------------------------------------
subroutine equal_to_model(diff_scr, timestep, equal_radius, equal_hydrogen_fraction, &
     zone_begin, zone_end, num_equal_points, composition, &
     aux_radial_quantity, radius, enclosed_mass, num_zones, total_mass)
      use composition_layout, only: json, i_h1, i_he4, i_metals, i_he3, &
           i_c12, i_o18
      use numerics_lib, only: interp, lagrange4
      use bahcall_loeb_units, only: bl_radius_scale, bl_mass_scale, bl_time_scale
      use diffusion_grid, only: locate_stencil
      implicit none
      type(diffusion_scratch), intent(inout) :: diff_scr

      double precision, intent(inout) :: timestep
      double precision, intent(in) :: equal_radius(json), &
           equal_hydrogen_fraction(json)
      integer, intent(in) :: zone_begin, zone_end, num_equal_points
      double precision, intent(inout) :: composition(15,json)
! aux_radial_quantity (originally HQPR): only scaled by bl_radius_scale
! here (parallel to radius), its physical meaning is not otherwise
! exercised in this file.
      double precision, intent(inout) :: aux_radial_quantity(json)
! radius: the caller's Bahcall-Loeb radius copy, read only
      double precision, intent(in) :: radius(json)
      double precision, intent(inout) :: enclosed_mass(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass




      double precision :: interp_factors(4)
      double precision :: hydrogen_floor
      integer :: zone_index
      double precision :: metal_max, metal_new, metal_scale_ratio
      integer :: j, search_start_index, k0
      logical :: found
      double precision :: target_radius
      double precision :: delta_x, hydrogen_max
      double precision :: delta_z, metal_floor

! TRANSFORM BACK TO ORIGINAL GRID OF MODEL POINTS FROM EQUALLY
! SPACED GRID.
      hydrogen_floor = 0.0D0
      do zone_index = zone_begin,1,-1
         composition(i_h1,zone_index)=max(composition(i_h1,zone_index) + &
              equal_hydrogen_fraction(1),hydrogen_floor)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
! NOTE THAT BECAUSE METALS SINK, AND HYDROGEN RISES, THE FAILSAFES
! ARE OPPOSITE (GUARDING AGAINST NEGATIVE X AND Z>1 RESPECTIVELY).
      if(diffusion_settings%use_diffusion_z)then
         do zone_index = zone_begin,1,-1
            metal_max = 1.0D0 - composition(i_h1,zone_index) - composition(i_he3,zone_index)
            metal_new=min(composition(i_metals,zone_index)+diff_scr%metal_abundance_change(1),metal_max)
            metal_scale_ratio = metal_new/composition(i_metals,zone_index)
            composition(i_metals,zone_index) = metal_new
            do j = i_c12,i_o18
               composition(j,zone_index) = metal_scale_ratio*composition(j,zone_index)
            end do
            composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                 composition(i_metals,zone_index)-composition(i_he3,zone_index)
         end do
      else
         do zone_index = zone_begin,1,-1
            composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                 composition(i_metals,zone_index)-composition(i_he3,zone_index)
         end do
      endif
      search_start_index=2
      do zone_index=zone_begin+1,zone_end-1
! THE FOUR EQUALLY SPACED GRID POINTS AROUND THE MODEL POINT
! (diffusion_grid), kept inside the grid when it has fewer than four
! points (docs/bug_sweep.md O21).
         call locate_stencil(equal_radius, num_equal_points, radius(zone_index), &
              search_start_index, num_equal_points, k0, found)
         if (found .and. k0 < 1) k0 = 1
         target_radius=radius(zone_index)
! FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
         call interp(equal_radius(k0:k0+3), interp_factors, x_eval=target_radius)
! PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR CHANGE IN X.
         delta_x = lagrange4(interp_factors, equal_hydrogen_fraction(k0:k0+3))
         hydrogen_max = 1.0D0 - composition(i_metals,zone_index) - composition(i_he3,zone_index)
         composition(i_h1,zone_index)=min(composition(i_h1,zone_index) + delta_x,hydrogen_max)
! MHP 3/94 ADDED METAL DIFFUSION
         if(diffusion_settings%use_diffusion_z)then
            metal_max = 1.0D0 - composition(i_h1,zone_index) - composition(i_he3,zone_index)
            delta_z = lagrange4(interp_factors, diff_scr%metal_abundance_change(k0:k0+3))
            metal_new = min(composition(i_metals,zone_index)+delta_z,metal_max)
            metal_scale_ratio = metal_new/composition(i_metals,zone_index)
            composition(i_metals,zone_index)=metal_new
            do j = i_c12,i_o18
               composition(j,zone_index) = metal_scale_ratio*composition(j,zone_index)
            end do
            composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                 composition(i_metals,zone_index)-composition(i_he3,zone_index)
         else
            composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                 composition(i_metals,zone_index)-composition(i_he3,zone_index)
         endif
      end do
      do zone_index = zone_end,num_zones
         hydrogen_max = 1.0D0 - composition(i_metals,zone_index) - composition(i_he3,zone_index)
         composition(i_h1,zone_index)=min(composition(i_h1,zone_index) + &
              equal_hydrogen_fraction(num_equal_points),hydrogen_max)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
      if(diffusion_settings%use_diffusion_z)then
         metal_floor = 0.0D0
         do zone_index = zone_end,num_zones
            metal_new = max(composition(i_metals,zone_index)+ &
                 diff_scr%metal_abundance_change(num_equal_points),metal_floor)
            metal_scale_ratio = metal_new/composition(i_metals,zone_index)
            composition(i_metals,zone_index) = metal_new
            do j = i_c12,i_o18
               composition(j,zone_index) = metal_scale_ratio*composition(j,zone_index)
            end do
            composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                 composition(i_metals,zone_index)-composition(i_he3,zone_index)
         end do
      else
         do zone_index = zone_end,num_zones
            composition(i_he4,zone_index)=1.0D0-composition(i_h1,zone_index)- &
                 composition(i_metals,zone_index)-composition(i_he3,zone_index)
         end do
      endif
      do zone_index=1,num_zones

         enclosed_mass(zone_index)=enclosed_mass(zone_index)/bl_mass_scale
         aux_radial_quantity(zone_index)=aux_radial_quantity(zone_index)*bl_radius_scale
      end do
      timestep=timestep*bl_time_scale
      total_mass=total_mass/bl_mass_scale
      return
end subroutine equal_to_model

end module settling_grid

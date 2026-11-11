!----------------------------------------------------------------------
! settling_grid
!----------------------------------------------------------------------
! The original routines' maps between the model and the equal grids
! (diffusion_grid builds them, as for microdiff):
!   model_to_equal  the grids over the settling region (intervals from
!       the narrowest of all the region's zones) and the coefficients,
!       their derivatives in X (and Z) and the abundances on them;
!   equal_to_model  the changes of X and Z back on the model's zones
!       (diffusion_grid's apply_grid_changes), and the timestep, the
!       enclosed mass, d ln P/dr and the total mass back from Bahcall &
!       Loeb units, in place (the round trip can change their last bit).
module settling_grid
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: json
      use element_diffusion_configuration, only: diffusion_settings
      use element_diffusion_scratch, only: diffusion_scratch
      implicit none
      private
      public :: model_to_equal, equal_to_model

contains

! ---------------------------------------------------------------
! The equal grids over zone_begin..zone_end (radius in Bahcall & Loeb
! units): their point count num_equal_points, spacing grid_spacing and
! centres equal_radius, and on them hydrogen's D1 (centres and
! midpoints), D2 and the two derivatives in X (midpoints), X, and the
! enclosed mass; with use_diffusion_z the metals' coefficients and Z
! (diff_scr's grid fields, from its zone fields).
subroutine model_to_equal(diff_scr, diffusion_coeff1, diffusion_coeff2, composition, &
     radius, enclosed_mass, &
     zone_begin, zone_end, num_zones, diffusion_coeff1_dx, &
     diffusion_coeff2_dx, grid_spacing, &
     equal_diffusion_coeff1, equal_diffusion_coeff1_mid, &
     equal_diffusion_coeff2_mid, equal_mass, equal_mass_mid, &
     equal_diffusion_coeff1_dx_mid, equal_diffusion_coeff2_dx_mid, &
     equal_radius, equal_hydrogen_fraction, equal_hydrogen_fraction_mid, &
     num_equal_points)
      use composition_layout, only: i_h1, i_metals
      use diffusion_grid, only: equal_grid_map, size_equal_grid, place_equal_grids, &
           build_equal_grid_map, map_field
      type(diffusion_scratch), intent(inout) :: diff_scr
      real(dp), intent(in) :: diffusion_coeff1(:), diffusion_coeff2(:)
      real(dp), intent(in) :: composition(:,:)
      real(dp), intent(in) :: radius(:), enclosed_mass(:)
      integer, intent(in) :: zone_begin, zone_end, num_zones
      real(dp), intent(in) :: diffusion_coeff1_dx(:), diffusion_coeff2_dx(:)
      real(dp), intent(out) :: grid_spacing
      real(dp), intent(inout) :: equal_diffusion_coeff1(:), &
           equal_diffusion_coeff1_mid(:), equal_diffusion_coeff2_mid(:), &
           equal_mass(:), equal_mass_mid(:), &
           equal_diffusion_coeff1_dx_mid(:), equal_diffusion_coeff2_dx_mid(:), &
           equal_radius(:), equal_hydrogen_fraction(:), &
           equal_hydrogen_fraction_mid(:)
      integer, intent(out) :: num_equal_points
! the zones each grid point is interpolated from, and the midpoints
      type(equal_grid_map) :: map
      real(dp) :: mid_radius(json)

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
! the metals, as fully ionized iron
      if (diffusion_settings%use_diffusion_z) then
         call map_field(map, diff_scr%zone_metal_coeff1, diff_scr%metal_coeff1, &
              diff_scr%metal_coeff1_mid)
         call map_field(map, diff_scr%zone_metal_coeff2, mid_value=diff_scr%metal_coeff2_mid)
         call map_field(map, diff_scr%zone_metal_coeff1_dz, mid_value=diff_scr%metal_coeff1_dz_mid)
         call map_field(map, diff_scr%zone_metal_coeff2_dz, mid_value=diff_scr%metal_coeff2_dz_mid)
         call map_field(map, composition(i_metals,:), diff_scr%eq_metal_fraction, &
              diff_scr%eq_metal_fraction_mid)
      end if
end subroutine model_to_equal

! ---------------------------------------------------------------
! The changes of X (equal_hydrogen_fraction) and Z (diff_scr's
! eq_metal_fraction, with use_diffusion_z) at the num_equal_points
! centres equal_radius applied to the model's zones (radius, Bahcall &
! Loeb units; diffusion_grid's apply_grid_changes, a grid of fewer than
! four points taking its stencil from point 1); then the timestep, the
! enclosed mass, dlnp_dr and the total mass converted back from Bahcall
! & Loeb units in place.
subroutine equal_to_model(diff_scr, timestep, equal_radius, equal_hydrogen_fraction, &
     zone_begin, zone_end, num_equal_points, composition, &
     dlnp_dr, radius, enclosed_mass, num_zones, total_mass)
      use bahcall_loeb_units, only: bl_radius_scale, bl_mass_scale, bl_time_scale
      use diffusion_grid, only: apply_grid_changes
      type(diffusion_scratch), intent(inout) :: diff_scr
      real(dp), intent(inout) :: timestep
      real(dp), intent(in) :: equal_radius(:), equal_hydrogen_fraction(:)
      integer, intent(in) :: zone_begin, zone_end, num_equal_points
      real(dp), intent(inout) :: composition(:,:)
      real(dp), intent(inout) :: dlnp_dr(:)
      real(dp), intent(in) :: radius(:)
      real(dp), intent(inout) :: enclosed_mass(:)
      integer, intent(in) :: num_zones
      real(dp), intent(inout) :: total_mass
! (the original routines settle no light element)
      integer, parameter :: no_rows(0) = [integer ::]
      real(dp) :: no_change(0,json)
      integer :: k

      call apply_grid_changes(equal_radius, num_equal_points, radius, zone_begin, zone_end, &
           num_zones, equal_hydrogen_fraction, diff_scr%eq_metal_fraction, no_change, &
           no_rows, diffusion_settings%use_diffusion_z, .false., .true., composition)
      do k = 1, num_zones
         enclosed_mass(k) = enclosed_mass(k)/bl_mass_scale
         dlnp_dr(k) = dlnp_dr(k)*bl_radius_scale
      end do
      timestep = timestep*bl_time_scale
      total_mass = total_mass/bl_mass_scale
end subroutine equal_to_model

end module settling_grid

! Sources
! - MHP (3/94: the metals); JvS (the guards for grids of one to three
!   points). 2026: the grids built by diffusion_grid, shared with
!   microdiff; the literal composition slots named (docs/bug_sweep.md
!   M252's commit).

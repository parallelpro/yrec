!----------------------------------------------------------------------
! equal_to_model
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original equal_to_model.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
subroutine equal_to_model(timestep, equal_radius, equal_hydrogen_fraction, &
     zone_begin, zone_end, num_equal_points, composition, &
     aux_radial_quantity, radius, enclosed_mass, temperature, num_zones, &
     total_mass)
      use star_info_lib, only: star, json
      use numerics_lib
      implicit none

      double precision, intent(inout) :: timestep
      double precision, intent(in) :: equal_radius(json), &
           equal_hydrogen_fraction(json)
      integer, intent(in) :: zone_begin, zone_end, num_equal_points
      double precision, intent(inout) :: composition(15,json)
! aux_radial_quantity (originally HQPR): only scaled by star%rot%bl_radius_scale
! here (parallel to radius), its physical meaning is not otherwise
! exercised in this file.
      double precision, intent(inout) :: aux_radial_quantity(json)
      double precision, intent(inout) :: radius(json)
      double precision, intent(inout) :: enclosed_mass(json), temperature(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass




      double precision :: radius_table(4), interp_factors(4)
      double precision :: hydrogen_floor
      integer :: zone_index
      double precision :: metal_max, metal_new, metal_scale_ratio
      integer :: j, search_start_index, k0, k
      double precision :: target_radius
      double precision :: delta_x, hydrogen_max
      double precision :: delta_z, metal_floor

! TRANSFORM BACK TO ORIGINAL GRID OF MODEL POINTS FROM EQUALLY
! SPACED GRID.
      hydrogen_floor = 0.0D0
      do zone_index = zone_begin,1,-1
         composition(1,zone_index)=max(composition(1,zone_index) + &
              equal_hydrogen_fraction(1),hydrogen_floor)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
! NOTE THAT BECAUSE METALS SINK, AND HYDROGEN RISES, THE FAILSAFES
! ARE OPPOSITE (GUARDING AGAINST NEGATIVE X AND Z>1 RESPECTIVELY).
      if(star%job%use_diffusion_z)then
         do zone_index = zone_begin,1,-1
            metal_max = 1.0D0 - composition(1,zone_index) - composition(4,zone_index)
            metal_new=min(composition(3,zone_index)+star%rot%metal_abundance_change(1),metal_max)
            metal_scale_ratio = metal_new/composition(3,zone_index)
            composition(3,zone_index) = metal_new
            do j = 5,11
               composition(j,zone_index) = metal_scale_ratio*composition(j,zone_index)
            end do
            composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                 composition(3,zone_index)-composition(4,zone_index)
         end do
      else
         do zone_index = zone_begin,1,-1
            composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                 composition(3,zone_index)-composition(4,zone_index)
         end do
      endif
      search_start_index=2
      do zone_index=zone_begin+1,zone_end-1
         do j=search_start_index,num_equal_points

! FIND EQUALLY SPACED GRID POINTS CLOSEST TO THE MODEL POINT.
            if (j .eq. 0) print*, 'line 47 etm'
            if(equal_radius(j).ge.radius(zone_index))then

! ENSURE THAT FIRST INTERP. POINT NO LESS THAN FIRST EQUALLY SPACED POINT.
               k0 = max(j-2,1)
! ENSURE THAT LAST INTERP. POINT NO GREATER THAN LAST EQUALLY SPACED POINT.
               k0 = min(k0,num_equal_points-3)
! JVS fix for NPT = 3?
               if (k0 .eq. 0) k0=1
               search_start_index=j
               exit
            endif
         end do
         if (j > num_equal_points) then
         k0 = num_equal_points-3
         search_start_index=num_equal_points
         end if

         do k=1,4
            radius_table(k)=equal_radius(k0+k-1)
         end do
         target_radius=radius(zone_index)
! FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
         call intrp2(radius_table,interp_factors,target_radius)
! PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR CHANGE IN X.
         delta_x = interp_factors(1)*equal_hydrogen_fraction(k0)+ &
              interp_factors(2)*equal_hydrogen_fraction(k0+1)+ &
              interp_factors(3)*equal_hydrogen_fraction(k0+2)+ &
              interp_factors(4)*equal_hydrogen_fraction(k0+3)
         hydrogen_max = 1.0D0 - composition(3,zone_index) - composition(4,zone_index)
         composition(1,zone_index)=min(composition(1,zone_index) + delta_x,hydrogen_max)
! MHP 3/94 ADDED METAL DIFFUSION
         if(star%job%use_diffusion_z)then
            metal_max = 1.0D0 - composition(1,zone_index) - composition(4,zone_index)
            delta_z = interp_factors(1)*star%rot%metal_abundance_change(k0)+ &
                 interp_factors(2)*star%rot%metal_abundance_change(k0+1)+ &
                 interp_factors(3)*star%rot%metal_abundance_change(k0+2)+ &
                 interp_factors(4)*star%rot%metal_abundance_change(k0+3)
            metal_new = min(composition(3,zone_index)+delta_z,metal_max)
            metal_scale_ratio = metal_new/composition(3,zone_index)
            composition(3,zone_index)=metal_new
            do j = 5,11
               composition(j,zone_index) = metal_scale_ratio*composition(j,zone_index)
            end do
            composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                 composition(3,zone_index)-composition(4,zone_index)
         else
            composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                 composition(3,zone_index)-composition(4,zone_index)
         endif
      end do
      do zone_index = zone_end,num_zones
         hydrogen_max = 1.0D0 - composition(3,zone_index) - composition(4,zone_index)
         composition(1,zone_index)=min(composition(1,zone_index) + &
              equal_hydrogen_fraction(num_equal_points),hydrogen_max)
      end do
! MHP 3/94 ADDED METAL DIFFUSION
      if(star%job%use_diffusion_z)then
         metal_floor = 0.0D0
         do zone_index = zone_end,num_zones
            metal_new = max(composition(3,zone_index)+ &
                 star%rot%metal_abundance_change(num_equal_points),metal_floor)
            metal_scale_ratio = metal_new/composition(3,zone_index)
            composition(3,zone_index) = metal_new
            do j = 5,11
               composition(j,zone_index) = metal_scale_ratio*composition(j,zone_index)
            end do
            composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                 composition(3,zone_index)-composition(4,zone_index)
         end do
      else
         do zone_index = zone_end,num_zones
            composition(2,zone_index)=1.0D0-composition(1,zone_index)- &
                 composition(3,zone_index)-composition(4,zone_index)
         end do
      endif
      do zone_index=1,num_zones

         radius(zone_index)=radius(zone_index)/star%rot%bl_radius_scale
         temperature(zone_index)=temperature(zone_index)/star%rot%bl_temp_scale
         enclosed_mass(zone_index)=enclosed_mass(zone_index)/star%rot%bl_mass_scale
         aux_radial_quantity(zone_index)=aux_radial_quantity(zone_index)*star%rot%bl_radius_scale
      end do
      timestep=timestep*star%rot%bl_time_scale
      total_mass=total_mass/star%rot%bl_mass_scale
      return
end subroutine equal_to_model

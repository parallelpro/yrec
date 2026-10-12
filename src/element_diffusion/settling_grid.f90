!----------------------------------------------------------------------
! settling_grid
!----------------------------------------------------------------------
! The original settling routines' grid mapping: the model onto an equally
! spaced grid in radius and back.
module settling_grid
      use element_diffusion_scratch, only: diff_scr
      implicit none
      private
      public :: model_to_equal, equal_to_model

contains

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! model_to_equal
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original model_to_equal.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! MHP 3/94 added metal diffusion.
subroutine model_to_equal(diffusion_coeff1, diffusion_coeff2, composition, &
     radius, enclosed_mass, &
     zone_begin, zone_end, num_zones, diffusion_coeff1_dx, &
     diffusion_coeff2_dx, grid_spacing, &
     equal_diffusion_coeff1, equal_diffusion_coeff1_mid, &
     equal_diffusion_coeff2_mid, equal_mass, equal_mass_mid, &
     equal_diffusion_coeff1_dx_mid, equal_diffusion_coeff2_dx_mid, &
     equal_radius, equal_hydrogen_fraction, equal_hydrogen_fraction_mid, &
     num_equal_points)
      use star_info_lib, only: star, json, i_metals
      use numerics_lib
      implicit none

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




      double precision :: interp_factors(4), radius_table(4)
      double precision :: total_radius_span, min_radius_spacing
      integer :: zone_index
      integer :: interp_search_index
      double precision :: interp_fraction
      integer :: search_start_index, j, k0, k
      double precision :: target_radius

! FIND MINIMUM MODEL POINT SPACING IN RADIUS.

      total_radius_span = radius(zone_end)-radius(zone_begin)
      min_radius_spacing = total_radius_span
      do zone_index = zone_begin+1,zone_end
         min_radius_spacing=min(min_radius_spacing,radius(zone_index)-radius(zone_index-1))
      end do
! ASSIGN THE MINIMUM NUMBER OF EQUALLY SPACED GRID POINTS SUCH THAT
! DR <= DRMIN.
      num_equal_points=int(total_radius_span/min_radius_spacing)
      if(mod(total_radius_span,min_radius_spacing).ne.0.0D0)num_equal_points=num_equal_points+1
! JVS add additional trap to deal with models with NPT=1
      if (num_equal_points .eq. 1) num_equal_points=num_equal_points+1
! ENSURE THAT NUMBER OF POINTS DOES NOT EXCEED JSON.

      num_equal_points=min(num_equal_points,json)
      grid_spacing = total_radius_span/dfloat(num_equal_points-1)

! SET UP VECTOR OF EQUALLY SPACED RADII AT ZONE MIDPOINTS.
      equal_radius(1)=radius(zone_begin)+0.5D0*grid_spacing
! JVS added logic trap (IF statement)
      if (num_equal_points .eq. 2) then
            equal_radius(2)=equal_radius(1)+grid_spacing
      else
            do zone_index = 2,num_equal_points-1
               equal_radius(zone_index)=equal_radius(zone_index-1)+grid_spacing
            end do
      endif

! NOW USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND RUN OF VARIABLES
! AT EQUALLY SPACED ZONE MIDPOINTS.
!
! FIRST POINT : LINEAR INTERPOLATION BETWEEN STARTING POINT AND 2ND PT.
      do interp_search_index=2,num_equal_points
         if(radius(interp_search_index).ge.equal_radius(1))exit
      end do
      if (interp_search_index > num_equal_points) then
      interp_search_index=num_equal_points
      end if
      interp_fraction=(equal_radius(1)-radius(interp_search_index-1))/ &
           (radius(interp_search_index)-radius(interp_search_index-1))
      equal_mass_mid(1) = enclosed_mass(interp_search_index-1)+interp_fraction* &
           (enclosed_mass(interp_search_index)-enclosed_mass(interp_search_index-1))
      equal_diffusion_coeff1_mid(1)=diffusion_coeff1(interp_search_index-1)+ &
           interp_fraction*(diffusion_coeff1(interp_search_index)- &
           diffusion_coeff1(interp_search_index-1))
      equal_diffusion_coeff2_mid(1)=diffusion_coeff2(interp_search_index-1)+ &
           interp_fraction*(diffusion_coeff2(interp_search_index)- &
           diffusion_coeff2(interp_search_index-1))
      equal_diffusion_coeff1_dx_mid(1)=diffusion_coeff1_dx(interp_search_index-1)+ &
           interp_fraction*(diffusion_coeff1_dx(interp_search_index)- &
           diffusion_coeff1_dx(interp_search_index-1))
      equal_diffusion_coeff2_dx_mid(1)=diffusion_coeff2_dx(interp_search_index-1)+ &
           interp_fraction*(diffusion_coeff2_dx(interp_search_index)- &
           diffusion_coeff2_dx(interp_search_index-1))
      equal_hydrogen_fraction_mid(1) = composition(1,interp_search_index-1)+ &
           interp_fraction*(composition(1,interp_search_index)- &
           composition(1,interp_search_index-1))
      if(star%job%use_diffusion_z)then
         diff_scr%metal_diffusion_coeff1_mid(1)=diff_scr%src_grid_metal_diffusion_coeff1(interp_search_index-1)+ &
              interp_fraction*(diff_scr%src_grid_metal_diffusion_coeff1(interp_search_index)- &
              diff_scr%src_grid_metal_diffusion_coeff1(interp_search_index-1))
         diff_scr%metal_diffusion_coeff2_mid(1)=diff_scr%src_grid_metal_diffusion_coeff2(interp_search_index-1)+ &
              interp_fraction*(diff_scr%src_grid_metal_diffusion_coeff2(interp_search_index)- &
              diff_scr%src_grid_metal_diffusion_coeff2(interp_search_index-1))
         diff_scr%eq_metal_diffusion_coeff1_mid(1)=diff_scr%src_grid_metal_diffusion_coeff1_dz(interp_search_index-1)+ &
              interp_fraction*(diff_scr%src_grid_metal_diffusion_coeff1_dz(interp_search_index)- &
              diff_scr%src_grid_metal_diffusion_coeff1_dz(interp_search_index-1))
         diff_scr%eq_metal_diffusion_coeff2_mid(1)=diff_scr%src_grid_metal_diffusion_coeff2_dz(interp_search_index-1)+ &
              interp_fraction*(diff_scr%src_grid_metal_diffusion_coeff2_dz(interp_search_index)- &
              diff_scr%src_grid_metal_diffusion_coeff2_dz(interp_search_index-1))
! 2026 (bugsweep Batch 2): the metal abundance carried to the equal
! grid is Z (slot i_metals = 3); the inherited code read slot 8 (N15)
! here and at the four sibling sites below. equal_to_model applies
! the change back to slot 3, so this was an N15-for-Z mismatch on the
! old (use_new_diffusion_routines = F) settling path.
         diff_scr%metal_abundance_change_mid(1) = composition(i_metals,interp_search_index-1)+ &
              interp_fraction*(composition(i_metals,interp_search_index)- &
              composition(i_metals,interp_search_index-1))
      endif
! FOR OTHER POINTS: FIRST FIND 4 NEAREST (IN RADIUS) MODEL POINTS
! AND THEN FIND LAGRANGIAN INTERPOLATION FACTORS. APPLY THEM TO FIND
! MODEL QUANTITIES AT THE EQUALLY SPACED GRID POINTS.
! JMIN IS THE UPPERMOST MODEL POINT ABOVE THE PREVIOUS EQUALLY SPACED
! GRID POINT (IN RADIUS).
      search_start_index=zone_begin+1
      do zone_index=2,num_equal_points-1
         do j = search_start_index,zone_end
! FIND 4 MODEL POINTS CLOSEST TO THE EQUALLY SPACED GRID POINT.
            if(radius(j).ge.equal_radius(zone_index))then
! ENSURE THAT FIRST INTERPOLATION POINT NO LESS THAN FIRST MODEL POINT.
               k0 = max(j-2,1)
! ENSURE THAT LAST INTERPOLATION POINT NO GREATER THAN LAST MODEL POINT.
               k0 = min(k0,num_zones-3)
               search_start_index=j
               exit
            endif
         end do
         if (j > zone_end) then
         k0 = num_zones-3
         search_start_index=num_zones
         end if
         do k=1,4
            radius_table(k)=radius(k0+k-1)
         end do
         target_radius=equal_radius(zone_index)
! FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS FOR POINT target_radius
! GIVEN THE 4 TABLE RADII IN radius_table.
         call intrp2(radius_table,interp_factors,target_radius)
! PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR DESIRED QUANTITIES:
! MASS WITHIN THE RADIUS ER
         equal_mass_mid(zone_index) = lagrange4(interp_factors, enclosed_mass(k0:k0+3))
! D1
         equal_diffusion_coeff1_mid(zone_index) = lagrange4(interp_factors, &
              diffusion_coeff1(k0:k0+3))
! D2
         equal_diffusion_coeff2_mid(zone_index) = lagrange4(interp_factors, &
              diffusion_coeff2(k0:k0+3))
! DERIVATIVE OF D1 WRT X
         equal_diffusion_coeff1_dx_mid(zone_index) = lagrange4(interp_factors, &
              diffusion_coeff1_dx(k0:k0+3))
! DERIVATIVE OF D2 WRT X
         equal_diffusion_coeff2_dx_mid(zone_index) = lagrange4(interp_factors, &
              diffusion_coeff2_dx(k0:k0+3))
! MASS FRACTION OF HYDROGEN
         equal_hydrogen_fraction_mid(zone_index) = lagrange4(interp_factors, composition(1,k0:k0+3))
         if(star%job%use_diffusion_z)then
! METAL DIFFUSION-TREATED AS FULLY IONIZED IRON.
! D1
         diff_scr%metal_diffusion_coeff1_mid(zone_index) = lagrange4(interp_factors, &
              diff_scr%src_grid_metal_diffusion_coeff1(k0:k0+3))
! D2
         diff_scr%metal_diffusion_coeff2_mid(zone_index) = lagrange4(interp_factors, &
              diff_scr%src_grid_metal_diffusion_coeff2(k0:k0+3))
! DERIVATIVE OF D1 WRT Z
         diff_scr%eq_metal_diffusion_coeff1_mid(zone_index) = lagrange4(interp_factors, &
              diff_scr%src_grid_metal_diffusion_coeff1_dz(k0:k0+3))
! DERIVATIVE OF D2 WRT Z
         diff_scr%eq_metal_diffusion_coeff2_mid(zone_index) = lagrange4(interp_factors, &
              diff_scr%src_grid_metal_diffusion_coeff2_dz(k0:k0+3))
! MASS FRACTION OF METALS
         diff_scr%metal_abundance_change_mid(zone_index) = lagrange4(interp_factors, &
              composition(i_metals,k0:k0+3))
         endif
      end do
! SET UP VECTOR OF EQUALLY SPACED RADII AT ZONE CENTERS.

      equal_radius(1)=radius(zone_begin)
      do zone_index = 2,num_equal_points
         equal_radius(zone_index)=equal_radius(zone_index-1)+grid_spacing
      end do


! NOW USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND RUN OF VARIABLES
! AT EQUALLY SPACED ZONE CENTERS.
!
! FIRST POINT : BY DEFINITION, AT STARTING POINT.
      equal_mass(1) = enclosed_mass(zone_begin)
      equal_hydrogen_fraction(1) = composition(1,zone_begin)
      equal_diffusion_coeff1(1)=diffusion_coeff1(zone_begin)
      if(star%job%use_diffusion_z)then
         diff_scr%metal_abundance_change(1) = composition(i_metals,zone_begin)
         diff_scr%metal_diffusion_coeff1(1) = diff_scr%src_grid_metal_diffusion_coeff1(zone_begin)
      endif
! FOR OTHER POINTS: FIRST FIND 4 NEAREST (IN RADIUS) MODEL POINTS
! AND THEN FIND LAGRANGIAN INTERPOLATION FACTORS. APPLY THEM TO FIND
! MODEL QUANTITIES AT THE EQUALLY SPACED GRID POINTS.
! JMIN IS THE UPPERMOST MODEL POINT ABOVE THE PREVIOUS EQUALLY SPACED
! GRID POINT (IN RADIUS).
      search_start_index=zone_begin+1
      do zone_index=2,num_equal_points-1
         do j = search_start_index,zone_end
! FIND 4 MODEL POINTS CLOSEST TO THE EQUALLY SPACED GRID POINT.
            if(radius(j).ge.equal_radius(zone_index))then
! ENSURE THAT FIRST INTERPOLATION POINT NO LESS THAN FIRST MODEL POINT.
               k0 = max(j-2,1)
! ENSURE THAT LAST INTERPOLATION POINT NO GREATER THAN LAST MODEL POINT.
               k0 = min(k0,num_zones-3)
               search_start_index=j
               exit
            endif
         end do
         if (j > zone_end) then
         k0 = num_zones-3
         search_start_index=num_zones
         end if
         do k=1,4
            radius_table(k)=radius(k0+k-1)
         end do

         target_radius=equal_radius(zone_index)

! FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS FOR POINT target_radius
! GIVEN THE 4 TABLE RADII IN radius_table.
         call intrp2(radius_table,interp_factors,target_radius)
! PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR DESIRED QUANTITIES:
! MASS WITHIN THE RADIUS ER
         equal_mass(zone_index) = lagrange4(interp_factors, enclosed_mass(k0:k0+3))
! D1
         equal_diffusion_coeff1(zone_index) = lagrange4(interp_factors, diffusion_coeff1(k0:k0+3))
! MASS FRACTION OF HYDROGEN

         equal_hydrogen_fraction(zone_index) = lagrange4(interp_factors, composition(1,k0:k0+3))
! METAL DIFFUSION
         if(star%job%use_diffusion_z)then
! D1
         diff_scr%metal_diffusion_coeff1(zone_index) = lagrange4(interp_factors, &
              diff_scr%src_grid_metal_diffusion_coeff1(k0:k0+3))
! MASS FRACTION OF METALS
         diff_scr%metal_abundance_change(zone_index) = lagrange4(interp_factors, &
              composition(i_metals,k0:k0+3))
         endif
      end do
! LAST POINT : BY DEFINITION, AT ENDING POINT.

      equal_mass(num_equal_points) = enclosed_mass(zone_end)
      equal_hydrogen_fraction(num_equal_points) = composition(1,zone_end)
      equal_diffusion_coeff1(num_equal_points)=diffusion_coeff1(zone_end)

      if(star%job%use_diffusion_z)then
         diff_scr%metal_abundance_change(num_equal_points) = composition(i_metals,zone_end)
         diff_scr%metal_diffusion_coeff1(num_equal_points)=diff_scr%src_grid_metal_diffusion_coeff1(zone_end)
      endif
      return
end subroutine model_to_equal

! ---------------------------------------------------------------
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
     aux_radial_quantity, radius, enclosed_mass, num_zones, total_mass)
      use star_info_lib, only: star, json
      use numerics_lib
      implicit none

      double precision, intent(inout) :: timestep
      double precision, intent(in) :: equal_radius(json), &
           equal_hydrogen_fraction(json)
      integer, intent(in) :: zone_begin, zone_end, num_equal_points
      double precision, intent(inout) :: composition(15,json)
! aux_radial_quantity (originally HQPR): only scaled by diff_scr%bl_radius_scale
! here (parallel to radius), its physical meaning is not otherwise
! exercised in this file.
      double precision, intent(inout) :: aux_radial_quantity(json)
! radius: the caller's Bahcall-Loeb radius copy; read only (until 2026
! W2 it and a temperature copy were divided back out of BL units here,
! which nothing read afterwards).
      double precision, intent(in) :: radius(json)
      double precision, intent(inout) :: enclosed_mass(json)
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
            metal_new=min(composition(3,zone_index)+diff_scr%metal_abundance_change(1),metal_max)
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
            if(equal_radius(j).ge.radius(zone_index))then

! ENSURE THAT FIRST INTERP. POINT NO LESS THAN FIRST EQUALLY SPACED POINT.
               k0 = max(j-2,1)
! ENSURE THAT LAST INTERP. POINT NO GREATER THAN LAST EQUALLY SPACED POINT.
               k0 = min(k0,num_equal_points-3)
! keep the 4-point stencil inside the table when it has fewer than
! four points (only reachable for num_equal_points <= 3)
               if (k0 .lt. 1) k0 = 1
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
         delta_x = lagrange4(interp_factors, equal_hydrogen_fraction(k0:k0+3))
         hydrogen_max = 1.0D0 - composition(3,zone_index) - composition(4,zone_index)
         composition(1,zone_index)=min(composition(1,zone_index) + delta_x,hydrogen_max)
! MHP 3/94 ADDED METAL DIFFUSION
         if(star%job%use_diffusion_z)then
            metal_max = 1.0D0 - composition(1,zone_index) - composition(4,zone_index)
            delta_z = lagrange4(interp_factors, diff_scr%metal_abundance_change(k0:k0+3))
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
                 diff_scr%metal_abundance_change(num_equal_points),metal_floor)
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

         enclosed_mass(zone_index)=enclosed_mass(zone_index)/diff_scr%bl_mass_scale
         aux_radial_quantity(zone_index)=aux_radial_quantity(zone_index)*diff_scr%bl_radius_scale
      end do
      timestep=timestep*diff_scr%bl_time_scale
      total_mass=total_mass/diff_scr%bl_mass_scale
      return
end subroutine equal_to_model

end module settling_grid

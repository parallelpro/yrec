!----------------------------------------------------------------------
! microdiff_mte
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original microdiff_mte.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! "MTE" = model-to-equally-spaced-grid transform: builds the equally
! spaced radial grid (at both zone centers and zone midpoints, in
! Bahcall & Loeb units) used by microdiff_run.f90/microdiff_cod.f90,
! and interpolates density, temperature, dlnP/dr, the "del" temperature
! gradient, and the H/He/metal/light-element mass fractions onto it via
! 4-point Lagrangian interpolation (see interp.f, not part of this
! batch). Part of the microdiff.f90 pipeline (see also
! microdiff_setup.f90, microdiff_cod.f90, microdiff_run.f90,
! microdiff_etm.f90).
subroutine microdiff_mte(num_light, light_element_id, composition, &
     dlnp_dr, radius_bl, enclosed_mass, zone_begin, zone_end, num_zones, &
     grid_spacing, num_eq_points, density_orig, temperature_orig, eq_mass, &
     eq_radius, eq_density, eq_temperature, eq_dlnp_dr, eq_del_grad, &
     eq_hydrogen, eq_helium, eq_metal, eq_light, eq_mass_mid, eq_radius_mid, &
     eq_density_mid, eq_temperature_mid, eq_dlnp_dr_mid, eq_del_grad_mid, &
     eq_hydrogen_mid, eq_helium_mid, eq_metal_mid, eq_light_mid)

      use const_lib
      use star_info_lib, only: star, i_grad_actual, json
      use numerics_lib
      implicit none

      integer, intent(in) :: num_light
      integer, intent(in) :: light_element_id(num_light)
      double precision, intent(in) :: composition(15,json), dlnp_dr(json), &
           radius_bl(json), enclosed_mass(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(out) :: grid_spacing
      integer, intent(out) :: num_eq_points
      double precision, intent(in) :: density_orig(json), temperature_orig(json)
      double precision, intent(out) :: eq_mass(json), eq_radius(json), &
           eq_density(json), eq_temperature(json), eq_dlnp_dr(json), &
           eq_del_grad(json), eq_hydrogen(json), eq_helium(json), &
           eq_metal(json)
      double precision, intent(out) :: eq_light(num_light,json)
      double precision, intent(out) :: eq_mass_mid(json), eq_radius_mid(json), &
           eq_density_mid(json), eq_temperature_mid(json), &
           eq_dlnp_dr_mid(json), eq_del_grad_mid(json), eq_hydrogen_mid(json), &
           eq_helium_mid(json), eq_metal_mid(json)
      double precision, intent(out) :: eq_light_mid(num_light,json)




      integer :: half_json
      integer :: i, ii, iu, j, jmin, k, k0, kk
      double precision :: drtot, drmin, fx, tabler(4), gridrad, &
           facinterp(4), facderiv(4)

!  FIND MINIMUM MODEL POINT SPACING IN RADIUS.
      drtot = radius_bl(zone_end)-radius_bl(zone_begin)
      drmin = drtot
!  ID THE MINIMUM SPACING IN THE 20 LAYERS BELOW THE SURFACE CZ.
      ii = max(zone_end-20,zone_begin+1)
      do i = ii,zone_end
         drmin=min(drmin,radius_bl(i)-radius_bl(i-1))
      end do
!  ASSIGN THE MINIMUM NUMBER OF EQUALLY SPACED GRID POINTS SUCH THAT
!  DR <= DRMIN.
      num_eq_points=int(drtot/drmin)
      if(mod(drtot,drmin).ne.0.0d0)num_eq_points=num_eq_points+1
! JVS add additional trap to deal with models with NPT=1
      if (num_eq_points .eq. 1) num_eq_points=num_eq_points+1
!  ENSURE THAT NUMBER OF POINTS DOES NOT EXCEED JSON.

      half_json = 5000
      num_eq_points=min(num_eq_points,half_json)
      grid_spacing = drtot/dfloat(num_eq_points-1)
!  SET UP VECTOR OF EQUALLY SPACED RADII AT ZONE MIDPOINTS.
      eq_radius_mid(1)=radius_bl(zone_begin)+0.5d0*grid_spacing
! JVS added logic trap (IF statement)
      if (num_eq_points .eq. 2) then
         eq_radius_mid(2)=eq_radius_mid(1)+grid_spacing
      else
         do i = 2,num_eq_points-1! old piece
            if(i-1 .eq. 0) print*, 'mte line 47'
            eq_radius_mid(i)=eq_radius_mid(i-1)+grid_spacing  ! old piece
         end do
      endif

!  NOW USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND RUN OF VARIABLES
!  AT EQUALLY SPACED ZONE MIDPOINTS.
!
!  FIRST POINT : LINEAR INTERPOLATION BETWEEN STARTING POINT AND 2ND PT.
      do iu=2,num_eq_points
         if(radius_bl(iu).ge.eq_radius_mid(1))exit
      end do
      if (iu > num_eq_points) then
      iu=num_eq_points
      end if
      fx=(eq_radius_mid(1)-radius_bl(iu-1))/(radius_bl(iu)-radius_bl(iu-1))
      eq_mass_mid(1) = enclosed_mass(iu-1)+fx*(enclosed_mass(iu)-enclosed_mass(iu-1))
      eq_density_mid(1) = density_orig(iu-1)+fx*(density_orig(iu)-density_orig(iu-1))
      eq_temperature_mid(1) = temperature_orig(iu-1)+ &
           fx*(temperature_orig(iu)-temperature_orig(iu-1))
      eq_dlnp_dr_mid(1) = dlnp_dr(iu-1)+fx*(dlnp_dr(iu)-dlnp_dr(iu-1))
      eq_del_grad_mid(1) = star%diag%del_grad(i_grad_actual,iu-1)+fx*(star%diag%del_grad(i_grad_actual,iu)-star%diag%del_grad(i_grad_actual,iu-1))
      eq_hydrogen_mid(1) = composition(1,iu-1)+fx*(composition(1,iu)-composition(1,iu-1))
      eq_helium_mid(1) = composition(2,iu-1)+fx*(composition(2,iu)-composition(2,iu-1))
      eq_metal_mid(1) = composition(3,iu-1)+fx*(composition(3,iu)-composition(3,iu-1))
      if(ldifli)then
         do kk=1,num_light
            ii = light_element_id(kk)
            eq_light_mid(kk,1) = composition(ii,iu-1)+ &
                 fx*(composition(ii,iu)-composition(ii,iu-1))
         end do
      endif

!  FOR OTHER POINTS: FIRST FIND 4 NEAREST (IN RADIUS) MODEL POINTS
!  AND THEN FIND LAGRANGIAN INTERPOLATION FACTORS. APPLY THEM TO FIND
!  MODEL QUANTITIES AT THE EQUALLY SPACED GRID POINTS.
!  JMIN IS THE UPPERMOST MODEL POINT ABOVE THE PREVIOUS EQUALLY SPACED
!  GRID POINT (IN RADIUS).
      jmin=zone_begin+1
      do i=2,num_eq_points-1
         do j = jmin,zone_end
!  FIND 4 MODEL POINTS CLOSEST TO THE EQUALLY SPACED GRID POINT.
            if(radius_bl(j).ge.eq_radius_mid(i))then
!  ENSURE THAT FIRST INTERPOLATION POINT NO LESS THAN FIRST MODEL POINT.
               k0 = max(j-2,1)
!  ENSURE THAT LAST INTERPOLATION POINT NO GREATER THAN LAST MODEL POINT.
               k0 = min(k0,num_zones-3)
               jmin=j
               exit
            endif
         end do
         if (j > zone_end) then
         k0 = num_zones-3
         jmin=num_zones
         end if
         do k=1,4
            tabler(k)=radius_bl(k0+k-1)
         end do
         gridrad=eq_radius_mid(i)
!  FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
!  FACINTERP=INTERPOLATION FACTORS FOR POINT GRIDRAD GIVEN THE 4 TABLE
!  RADII IN TABLER; FACDERIV=SAME FOR DERIVATIVES AT POINT GRIDRAD.
         call interp(tabler,facinterp,facderiv,gridrad)
!  PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR DESIRED QUANTITIES:
!  MASS WITHIN THE RADIUS ER
         eq_mass_mid(i) = facinterp(1)*enclosed_mass(k0)+facinterp(2)*enclosed_mass(k0+1)+ &
                   facinterp(3)*enclosed_mass(k0+2)+facinterp(4)*enclosed_mass(k0+3)
!  RELAVENT PHYSICAL VARIABLES
         eq_density_mid(i) = facinterp(1)*density_orig(k0)+facinterp(2)*density_orig(k0+1)+ &
                   facinterp(3)*density_orig(k0+2)+facinterp(4)*density_orig(k0+3)
         eq_temperature_mid(i) = facinterp(1)*temperature_orig(k0)+facinterp(2)*temperature_orig(k0+1)+ &
                   facinterp(3)*temperature_orig(k0+2)+facinterp(4)*temperature_orig(k0+3)
         eq_dlnp_dr_mid(i) = facinterp(1)*dlnp_dr(k0)+facinterp(2)*dlnp_dr(k0+1)+ &
                   facinterp(3)*dlnp_dr(k0+2)+facinterp(4)*dlnp_dr(k0+3)
         eq_del_grad_mid(i) = facinterp(1)*star%diag%del_grad(i_grad_actual,k0)+facinterp(2)*star%diag%del_grad(i_grad_actual,k0+1)+ &
                   facinterp(3)*star%diag%del_grad(i_grad_actual,k0+2)+facinterp(4)*star%diag%del_grad(i_grad_actual,k0+3)
!  MASS FRACTION OF HYDROGEN
         eq_hydrogen_mid(i)=facinterp(1)*composition(1,k0) &
              +facinterp(2)*composition(1,k0+1) &
              +facinterp(3)*composition(1,k0+2) &
              +facinterp(4)*composition(1,k0+3)
!  MASS FRACTION OF HELIUM
         eq_helium_mid(i)=facinterp(1)*composition(2,k0) &
              +facinterp(2)*composition(2,k0+1) &
              +facinterp(3)*composition(2,k0+2) &
              +facinterp(4)*composition(2,k0+3)
!  MASS FRACTION OF METALS
         eq_metal_mid(i)=facinterp(1)*composition(3,k0) &
              +facinterp(2)*composition(3,k0+1) &
              +facinterp(3)*composition(3,k0+2) &
              +facinterp(4)*composition(3,k0+3)
!  MASS FRACTION OF LIGHT ELEMENTS
         if(ldifli)then
            do kk=1,num_light
               ii = light_element_id(kk)
               eq_light_mid(kk,i)=facinterp(1)*composition(ii,k0) &
                       +facinterp(2)*composition(ii,k0+1) &
                       +facinterp(3)*composition(ii,k0+2) &
                       +facinterp(4)*composition(ii,k0+3)
            end do
         endif
      end do

!  SET UP VECTOR OF EQUALLY SPACED RADII AT ZONE CENTERS.
      eq_radius(1)=radius_bl(zone_begin)
      do i = 2,num_eq_points
         eq_radius(i)=eq_radius(i-1)+grid_spacing
      end do

!  NOW USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND RUN OF VARIABLES
!  AT EQUALLY SPACED ZONE CENTERS.
!
!  FIRST POINT : BY DEFINITION, AT STARTING POINT.
! G Somers; added interpolation for Xfrac, Rho, T, and HQPR.
      eq_mass(1) = enclosed_mass(zone_begin)
      eq_density(1) = density_orig(zone_begin)
      eq_temperature(1) = temperature_orig(zone_begin)
      eq_dlnp_dr(1) = dlnp_dr(zone_begin)
      eq_del_grad(1) = star%diag%del_grad(i_grad_actual,zone_begin)
      eq_hydrogen(1) = composition(1,zone_begin)
      eq_helium(1) = composition(2,zone_begin)
      eq_metal(1) = composition(3,zone_begin)
      if(ldifli)then
         do kk=1,num_light
            ii = light_element_id(kk)
            eq_light(kk,1) = composition(ii,zone_begin)
         end do
      endif
!  FOR OTHER POINTS: FIRST FIND 4 NEAREST (IN RADIUS) MODEL POINTS
!  AND THEN FIND LAGRANGIAN INTERPOLATION FACTORS. APPLY THEM TO FIND
!  MODEL QUANTITIES AT THE EQUALLY SPACED GRID POINTS.
!  JMIN IS THE UPPERMOST MODEL POINT ABOVE THE PREVIOUS EQUALLY SPACED
!  GRID POINT (IN RADIUS).
      jmin=zone_begin+1
      do i=2,num_eq_points-1
         do j = jmin,zone_end
!  FIND 4 MODEL POINTS CLOSEST TO THE EQUALLY SPACED GRID POINT.
            if(radius_bl(j).ge.eq_radius(i))then
!  ENSURE THAT FIRST INTERPOLATION POINT NO LESS THAN FIRST MODEL POINT.
               k0 = max(j-2,1)
!  ENSURE THAT LAST INTERPOLATION POINT NO GREATER THAN LAST MODEL POINT.
               k0 = min(k0,num_zones-3)
               jmin=j
               exit
            endif
         end do
         if (j > zone_end) then
         k0 = num_zones-3
         jmin=num_zones
         end if
         do k=1,4
            tabler(k)=radius_bl(k0+k-1)
         end do
         gridrad=eq_radius(i)
!  FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
!  FACINTERP=INTERPOLATION FACTORS FOR POINT GRIDRAD GIVEN THE 4 TABLE
!  RADII IN TABLER; FACDERIV=SAME FOR DERIVATIVES AT POINT GRIDRAD.
         call interp(tabler,facinterp,facderiv,gridrad)
!  PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR DESIRED QUANTITIES:
!  MASS WITHIN THE RADIUS ER
         eq_mass(i) = facinterp(1)*enclosed_mass(k0)+facinterp(2)*enclosed_mass(k0+1)+ &
                   facinterp(3)*enclosed_mass(k0+2)+facinterp(4)*enclosed_mass(k0+3)
         eq_density(i) = facinterp(1)*density_orig(k0)+facinterp(2)*density_orig(k0+1)+ &
                   facinterp(3)*density_orig(k0+2)+facinterp(4)*density_orig(k0+3)
         eq_temperature(i) = facinterp(1)*temperature_orig(k0)+facinterp(2)*temperature_orig(k0+1)+ &
                   facinterp(3)*temperature_orig(k0+2)+facinterp(4)*temperature_orig(k0+3)
         eq_dlnp_dr(i) = facinterp(1)*dlnp_dr(k0)+facinterp(2)*dlnp_dr(k0+1)+ &
                   facinterp(3)*dlnp_dr(k0+2)+facinterp(4)*dlnp_dr(k0+3)
         eq_del_grad(i) = facinterp(1)*star%diag%del_grad(i_grad_actual,k0)+facinterp(2)*star%diag%del_grad(i_grad_actual,k0+1)+ &
                   facinterp(3)*star%diag%del_grad(i_grad_actual,k0+2)+facinterp(4)*star%diag%del_grad(i_grad_actual,k0+3)
         eq_hydrogen(i)=facinterp(1)*composition(1,k0) &
              +facinterp(2)*composition(1,k0+1) &
              +facinterp(3)*composition(1,k0+2) &
              +facinterp(4)*composition(1,k0+3)
         eq_helium(i)=facinterp(1)*composition(2,k0) &
              +facinterp(2)*composition(2,k0+1) &
              +facinterp(3)*composition(2,k0+2) &
              +facinterp(4)*composition(2,k0+3)
         eq_metal(i)=facinterp(1)*composition(3,k0) &
              +facinterp(2)*composition(3,k0+1) &
              +facinterp(3)*composition(3,k0+2) &
              +facinterp(4)*composition(3,k0+3)
         if(ldifli)then
            do kk=1,num_light
               ii = light_element_id(kk)
               eq_light(kk,i)=facinterp(1)*composition(ii,k0) &
                    +facinterp(2)*composition(ii,k0+1) &
                    +facinterp(3)*composition(ii,k0+2) &
                    +facinterp(4)*composition(ii,k0+3)
            end do
         endif
      end do
!  LAST POINT : BY DEFINITION, AT ENDING POINT.
      eq_mass(num_eq_points) = enclosed_mass(zone_end)
      eq_density(num_eq_points) = density_orig(zone_end)
      eq_temperature(num_eq_points) = temperature_orig(zone_end)
      eq_dlnp_dr(num_eq_points) = dlnp_dr(zone_end)
      eq_del_grad(num_eq_points) = star%diag%del_grad(i_grad_actual,zone_end)
      eq_hydrogen(num_eq_points) = composition(1,zone_end)
      eq_helium(num_eq_points) = composition(2,zone_end)
      eq_metal(num_eq_points) = composition(3,zone_end)
      if(ldifli)then
         do kk=1,num_light
            ii = light_element_id(kk)
            eq_light(kk,num_eq_points) = composition(ii,zone_end)
         end do
      endif
      return
end subroutine microdiff_mte

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
! Bahcall & Loeb units) used by microdiff_run.f90/microdiff_coefficients.f90,
! and interpolates density, temperature, dlnP/dr, the "del" temperature
! gradient, and the H/He/metal/light-element mass fractions onto it via
! 4-point Lagrangian interpolation (interp, numerics_lib). Part of
! the microdiff.f90 pipeline (see also
! microdiff_setup.f90, microdiff_coefficients.f90, microdiff_run.f90,
! microdiff_etm.f90).
! 2026 de-tramp (ROADMAP item 3): 33 arguments -> 15. The 22
! equally-spaced-grid output arrays are two microdiff_grid records
! (zone centers, eq; zone midpoints, eq_mid), defined here and
! threaded through microdiff_run/microdiff_coefficients, which are
! called once per grid instance.
module microdiff_mte_lib
      use star_info_lib, only: json
      implicit none

! One equally spaced radial grid of the microdiff settling pipeline
! (Bahcall & Loeb units): structure plus interpolated composition.
! light is dimensioned for the pipeline's three light elements
! (Li6/Li7/Be9; num_light in microdiff.f90).
      type :: microdiff_grid
         double precision :: mass(json), radius(json), density(json), &
              temperature(json), dlnp_dr(json), del_grad(json), &
              hydrogen(json), helium(json), metal(json), &
              light(3,json)
      end type microdiff_grid
contains

subroutine microdiff_mte(num_light, light_element_id, composition, &
     dlnp_dr, radius_bl, enclosed_mass, zone_begin, zone_end, num_zones, &
     grid_spacing, num_eq_points, density_orig, temperature_orig, &
     eq, eq_mid)

      use star_info_lib, only: star, json
      use numerics_lib, only: interp, intrp2, lagrange4
      implicit none

      integer, intent(in) :: num_light
      integer, intent(in) :: light_element_id(num_light)
      double precision, intent(in) :: composition(15,json), dlnp_dr(json), &
           radius_bl(json), enclosed_mass(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones
      double precision, intent(out) :: grid_spacing
      integer, intent(out) :: num_eq_points
      double precision, intent(in) :: density_orig(json), temperature_orig(json)
      type(microdiff_grid), intent(out) :: eq, eq_mid
! --- locals ---
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
      num_eq_points=min(num_eq_points,json)
      grid_spacing = drtot/dfloat(num_eq_points-1)
!  SET UP VECTOR OF EQUALLY SPACED RADII AT ZONE MIDPOINTS.
      eq_mid%radius(1)=radius_bl(zone_begin)+0.5d0*grid_spacing
! JVS added logic trap (IF statement)
      if (num_eq_points .eq. 2) then
         eq_mid%radius(2)=eq_mid%radius(1)+grid_spacing
      else
         do i = 2,num_eq_points-1
            eq_mid%radius(i)=eq_mid%radius(i-1)+grid_spacing
         end do
      endif

!  NOW USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND RUN OF VARIABLES
!  AT EQUALLY SPACED ZONE MIDPOINTS.
!
!  FIRST POINT : LINEAR INTERPOLATION BETWEEN STARTING POINT AND 2ND PT.
      do iu=2,num_eq_points
         if(radius_bl(iu).ge.eq_mid%radius(1))exit
      end do
      if (iu > num_eq_points) then
      iu=num_eq_points
      end if
      fx=(eq_mid%radius(1)-radius_bl(iu-1))/(radius_bl(iu)-radius_bl(iu-1))
      eq_mid%mass(1) = enclosed_mass(iu-1)+fx*(enclosed_mass(iu)-enclosed_mass(iu-1))
      eq_mid%density(1) = density_orig(iu-1)+fx*(density_orig(iu)-density_orig(iu-1))
      eq_mid%temperature(1) = temperature_orig(iu-1)+ &
           fx*(temperature_orig(iu)-temperature_orig(iu-1))
      eq_mid%dlnp_dr(1) = dlnp_dr(iu-1)+fx*(dlnp_dr(iu)-dlnp_dr(iu-1))
      eq_mid%del_grad(1) = star%gradT(iu-1)+fx*(star%gradT(iu)-star%gradT(iu-1))
      eq_mid%hydrogen(1) = composition(1,iu-1)+fx*(composition(1,iu)-composition(1,iu-1))
      eq_mid%helium(1) = composition(2,iu-1)+fx*(composition(2,iu)-composition(2,iu-1))
      eq_mid%metal(1) = composition(3,iu-1)+fx*(composition(3,iu)-composition(3,iu-1))
      if(star%ctrl%diffuse_lithium)then
         do kk=1,num_light
            ii = light_element_id(kk)
            eq_mid%light(kk,1) = composition(ii,iu-1)+ &
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
            if(radius_bl(j).ge.eq_mid%radius(i))then
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
         gridrad=eq_mid%radius(i)
!  FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
!  FACINTERP=INTERPOLATION FACTORS FOR POINT GRIDRAD GIVEN THE 4 TABLE
!  RADII IN TABLER; FACDERIV=SAME FOR DERIVATIVES AT POINT GRIDRAD.
         call interp(tabler,facinterp,facderiv,gridrad)
!  PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR DESIRED QUANTITIES:
!  MASS WITHIN THE RADIUS ER
         eq_mid%mass(i) = lagrange4(facinterp, enclosed_mass(k0:k0+3))
!  RELAVENT PHYSICAL VARIABLES
         eq_mid%density(i) = lagrange4(facinterp, density_orig(k0:k0+3))
         eq_mid%temperature(i) = lagrange4(facinterp, temperature_orig(k0:k0+3))
         eq_mid%dlnp_dr(i) = lagrange4(facinterp, dlnp_dr(k0:k0+3))
         eq_mid%del_grad(i) = lagrange4(facinterp, star%gradT(k0:k0+3))
!  MASS FRACTION OF HYDROGEN
         eq_mid%hydrogen(i)=lagrange4(facinterp, composition(1,k0:k0+3))
!  MASS FRACTION OF HELIUM
         eq_mid%helium(i)=lagrange4(facinterp, composition(2,k0:k0+3))
!  MASS FRACTION OF METALS
         eq_mid%metal(i)=lagrange4(facinterp, composition(3,k0:k0+3))
!  MASS FRACTION OF LIGHT ELEMENTS
         if(star%ctrl%diffuse_lithium)then
            do kk=1,num_light
               ii = light_element_id(kk)
               eq_mid%light(kk,i)=lagrange4(facinterp, composition(ii,k0:k0+3))
            end do
         endif
      end do

!  SET UP VECTOR OF EQUALLY SPACED RADII AT ZONE CENTERS.
      eq%radius(1)=radius_bl(zone_begin)
      do i = 2,num_eq_points
         eq%radius(i)=eq%radius(i-1)+grid_spacing
      end do

!  NOW USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND RUN OF VARIABLES
!  AT EQUALLY SPACED ZONE CENTERS.
!
!  FIRST POINT : BY DEFINITION, AT STARTING POINT.
! G Somers; added interpolation for Xfrac, Rho, T, and HQPR.
      eq%mass(1) = enclosed_mass(zone_begin)
      eq%density(1) = density_orig(zone_begin)
      eq%temperature(1) = temperature_orig(zone_begin)
      eq%dlnp_dr(1) = dlnp_dr(zone_begin)
      eq%del_grad(1) = star%gradT(zone_begin)
      eq%hydrogen(1) = composition(1,zone_begin)
      eq%helium(1) = composition(2,zone_begin)
      eq%metal(1) = composition(3,zone_begin)
      if(star%ctrl%diffuse_lithium)then
         do kk=1,num_light
            ii = light_element_id(kk)
            eq%light(kk,1) = composition(ii,zone_begin)
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
            if(radius_bl(j).ge.eq%radius(i))then
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
         gridrad=eq%radius(i)
!  FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
!  FACINTERP=INTERPOLATION FACTORS FOR POINT GRIDRAD GIVEN THE 4 TABLE
!  RADII IN TABLER; FACDERIV=SAME FOR DERIVATIVES AT POINT GRIDRAD.
         call interp(tabler,facinterp,facderiv,gridrad)
!  PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR DESIRED QUANTITIES:
!  MASS WITHIN THE RADIUS ER
         eq%mass(i) = lagrange4(facinterp, enclosed_mass(k0:k0+3))
         eq%density(i) = lagrange4(facinterp, density_orig(k0:k0+3))
         eq%temperature(i) = lagrange4(facinterp, temperature_orig(k0:k0+3))
         eq%dlnp_dr(i) = lagrange4(facinterp, dlnp_dr(k0:k0+3))
         eq%del_grad(i) = lagrange4(facinterp, star%gradT(k0:k0+3))
         eq%hydrogen(i)=lagrange4(facinterp, composition(1,k0:k0+3))
         eq%helium(i)=lagrange4(facinterp, composition(2,k0:k0+3))
         eq%metal(i)=lagrange4(facinterp, composition(3,k0:k0+3))
         if(star%ctrl%diffuse_lithium)then
            do kk=1,num_light
               ii = light_element_id(kk)
               eq%light(kk,i)=lagrange4(facinterp, composition(ii,k0:k0+3))
            end do
         endif
      end do
!  LAST POINT : BY DEFINITION, AT ENDING POINT.
      eq%mass(num_eq_points) = enclosed_mass(zone_end)
      eq%density(num_eq_points) = density_orig(zone_end)
      eq%temperature(num_eq_points) = temperature_orig(zone_end)
      eq%dlnp_dr(num_eq_points) = dlnp_dr(zone_end)
      eq%del_grad(num_eq_points) = star%gradT(zone_end)
      eq%hydrogen(num_eq_points) = composition(1,zone_end)
      eq%helium(num_eq_points) = composition(2,zone_end)
      eq%metal(num_eq_points) = composition(3,zone_end)
      if(star%ctrl%diffuse_lithium)then
         do kk=1,num_light
            ii = light_element_id(kk)
            eq%light(kk,num_eq_points) = composition(ii,zone_end)
         end do
      endif
      return
end subroutine microdiff_mte

end module microdiff_mte_lib

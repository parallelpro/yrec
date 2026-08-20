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
      implicit none
      integer, parameter :: json = 5000

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

! MHP 3/94 ADDED METAL DIFFUSION.
! common/gravs3/: only use_diffusion_z is used here. Naming matches
! setupopac.f90.
      double precision :: fgry, fgrz
      logical :: lthoul, use_diffusion_z
      common/gravs3/ fgry, fgrz, lthoul, use_diffusion_z

! common/gravsz/: the metal-diffusion analogs of the diffusion_coeff1/
! diffusion_coeff2/diffusion_coeff1_dx/diffusion_coeff2_dx dummy
! arguments above, on the original (non-equally-spaced) model grid.
! Not referenced in any already-converted file; names are chosen to
! parallel common/gravez/'s already-established names for the
! equally-spaced-grid analogs (metal_diffusion_coeff1 etc, in
! lax_wendrof1.f90/lax_wendrof2.f90), prefixed with "src_grid_" to
! avoid colliding with those names within this same file.
      double precision :: src_grid_metal_diffusion_coeff1(json), &
           src_grid_metal_diffusion_coeff2(json), &
           src_grid_metal_diffusion_coeff1_dz(json), &
           src_grid_metal_diffusion_coeff2_dz(json)
      common/gravsz/ src_grid_metal_diffusion_coeff1, &
           src_grid_metal_diffusion_coeff2, &
           src_grid_metal_diffusion_coeff1_dz, &
           src_grid_metal_diffusion_coeff2_dz

! common/gravez/: metal (Z) diffusion work arrays on the equally spaced
! grid; only metal_diffusion_coeff1/metal_diffusion_coeff1_mid/
! metal_diffusion_coeff2_mid/eq_metal_diffusion_coeff1_mid/
! eq_metal_diffusion_coeff2_mid/metal_abundance_change/
! metal_abundance_change_mid are used here. Naming matches
! lax_wendrof1.f90/lax_wendrof2.f90.
      double precision :: metal_diffusion_coeff1(json), &
           metal_diffusion_coeff1_mid(json), metal_diffusion_coeff2_mid(json), &
           eq_metal_diffusion_coeff1_mid(json), &
           eq_metal_diffusion_coeff2_mid(json), metal_abundance_change(json), &
           metal_abundance_change_mid(json)
      common/gravez/ metal_diffusion_coeff1, metal_diffusion_coeff1_mid, &
           metal_diffusion_coeff2_mid, eq_metal_diffusion_coeff1_mid, &
           eq_metal_diffusion_coeff2_mid, metal_abundance_change, &
           metal_abundance_change_mid

      double precision :: interp_factors(4), deriv_factors(4), radius_table(4)
      save

      double precision :: total_radius_span, min_radius_spacing
      integer :: zone_index
      integer :: interp_search_index
      double precision :: interp_fraction
      double precision :: dr1, dr2, fac1, fac2, delr
      integer :: search_start_index, j, k0, k
      double precision :: target_radius

! FIND MINIMUM MODEL POINT SPACING IN RADIUS.

      total_radius_span = radius(zone_end)-radius(zone_begin)
      min_radius_spacing = total_radius_span
      do 7 zone_index = zone_begin+1,zone_end
         min_radius_spacing=min(min_radius_spacing,radius(zone_index)-radius(zone_index-1))
    7 continue
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
            do 10 zone_index = 2,num_equal_points-1 ! old piece
            if(zone_index-1 .eq. 0) print*, 'mte line 47'
               equal_radius(zone_index)=equal_radius(zone_index-1)+grid_spacing  ! old piece
  10             continue          ! old piece
      endif

! NOW USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND RUN OF VARIABLES
! AT EQUALLY SPACED ZONE MIDPOINTS.
!
! FIRST POINT : LINEAR INTERPOLATION BETWEEN STARTING POINT AND 2ND PT.
      do 15 interp_search_index=2,num_equal_points
         if(radius(interp_search_index).ge.equal_radius(1))goto 17
   15 continue
      interp_search_index=num_equal_points
   17 continue
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
      if(use_diffusion_z)then
         metal_diffusion_coeff1_mid(1)=src_grid_metal_diffusion_coeff1(interp_search_index-1)+ &
              interp_fraction*(src_grid_metal_diffusion_coeff1(interp_search_index)- &
              src_grid_metal_diffusion_coeff1(interp_search_index-1))
         metal_diffusion_coeff2_mid(1)=src_grid_metal_diffusion_coeff2(interp_search_index-1)+ &
              interp_fraction*(src_grid_metal_diffusion_coeff2(interp_search_index)- &
              src_grid_metal_diffusion_coeff2(interp_search_index-1))
         eq_metal_diffusion_coeff1_mid(1)=src_grid_metal_diffusion_coeff1_dz(interp_search_index-1)+ &
              interp_fraction*(src_grid_metal_diffusion_coeff1_dz(interp_search_index)- &
              src_grid_metal_diffusion_coeff1_dz(interp_search_index-1))
         eq_metal_diffusion_coeff2_mid(1)=src_grid_metal_diffusion_coeff2_dz(interp_search_index-1)+ &
              interp_fraction*(src_grid_metal_diffusion_coeff2_dz(interp_search_index)- &
              src_grid_metal_diffusion_coeff2_dz(interp_search_index-1))
         metal_abundance_change_mid(1) = composition(8,interp_search_index-1)+ &
              interp_fraction*(composition(8,interp_search_index)- &
              composition(8,interp_search_index-1))
      endif
! CENTER DERIVATIVE.
      dr1=equal_radius(1)-radius(interp_search_index-1)
      dr2=radius(interp_search_index)-equal_radius(1)
      if(dr2.gt.dr1)then
         fac1=1.0D0
         fac2=dr1/dr2
         delr=2.0D0*dr1
      else if(dr1.gt.dr2)then
         fac1=dr2/dr1
         fac2=1.0D0
         delr=2.0D0*dr2
      else
         fac1=1.0D0
         fac2=1.0D0
         delr=dr1+dr2
      endif
! FOR OTHER POINTS: FIRST FIND 4 NEAREST (IN RADIUS) MODEL POINTS
! AND THEN FIND LAGRANGIAN INTERPOLATION FACTORS. APPLY THEM TO FIND
! MODEL QUANTITIES AT THE EQUALLY SPACED GRID POINTS.
! JMIN IS THE UPPERMOST MODEL POINT ABOVE THE PREVIOUS EQUALLY SPACED
! GRID POINT (IN RADIUS).
      search_start_index=zone_begin+1
      do 50 zone_index=2,num_equal_points-1
         do 20 j = search_start_index,zone_end
! FIND 4 MODEL POINTS CLOSEST TO THE EQUALLY SPACED GRID POINT.
            if(radius(j).ge.equal_radius(zone_index))then
! ENSURE THAT FIRST INTERPOLATION POINT NO LESS THAN FIRST MODEL POINT.
               k0 = max(j-2,1)
! ENSURE THAT LAST INTERPOLATION POINT NO GREATER THAN LAST MODEL POINT.
               k0 = min(k0,num_zones-3)
               search_start_index=j
               goto 30
            endif
   20    continue
         k0 = num_zones-3
         search_start_index=num_zones
   30    continue
         do 40 k=1,4
            radius_table(k)=radius(k0+k-1)
   40    continue
         target_radius=equal_radius(zone_index)
! FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
! FACINTERP=INTERPOLATION FACTORS FOR POINT GRIDRAD GIVEN THE 4 TABLE
! RADII IN TABLER; FACDERIV=SAME FOR DERIVATIVES AT POINT GRIDRAD.
         call interp(radius_table,interp_factors,deriv_factors,target_radius)
! PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR DESIRED QUANTITIES:
! MASS WITHIN THE RADIUS ER
         equal_mass_mid(zone_index) = interp_factors(1)*enclosed_mass(k0)+ &
              interp_factors(2)*enclosed_mass(k0+1)+ &
              interp_factors(3)*enclosed_mass(k0+2)+ &
              interp_factors(4)*enclosed_mass(k0+3)
! D1
         equal_diffusion_coeff1_mid(zone_index) = interp_factors(1)* &
              diffusion_coeff1(k0)+interp_factors(2)*diffusion_coeff1(k0+1)+ &
              interp_factors(3)*diffusion_coeff1(k0+2)+ &
              interp_factors(4)*diffusion_coeff1(k0+3)
! D2
         equal_diffusion_coeff2_mid(zone_index) = interp_factors(1)* &
              diffusion_coeff2(k0)+interp_factors(2)*diffusion_coeff2(k0+1)+ &
              interp_factors(3)*diffusion_coeff2(k0+2)+ &
              interp_factors(4)*diffusion_coeff2(k0+3)
! DERIVATIVE OF D1 WRT X
         equal_diffusion_coeff1_dx_mid(zone_index)=interp_factors(1)* &
              diffusion_coeff1_dx(k0)+interp_factors(2)*diffusion_coeff1_dx(k0+1) &
              +interp_factors(3)*diffusion_coeff1_dx(k0+2)+ &
              interp_factors(4)*diffusion_coeff1_dx(k0+3)
! DERIVATIVE OF D2 WRT X
         equal_diffusion_coeff2_dx_mid(zone_index)=interp_factors(1)* &
              diffusion_coeff2_dx(k0)+interp_factors(2)*diffusion_coeff2_dx(k0+1) &
              +interp_factors(3)*diffusion_coeff2_dx(k0+2)+ &
              interp_factors(4)*diffusion_coeff2_dx(k0+3)
! MASS FRACTION OF HYDROGEN
         equal_hydrogen_fraction_mid(zone_index)=interp_factors(1)* &
              composition(1,k0)+interp_factors(2)*composition(1,k0+1)+ &
              interp_factors(3)*composition(1,k0+2)+ &
              interp_factors(4)*composition(1,k0+3)
         if(use_diffusion_z)then
! METAL DIFFUSION-TREATED AS FULLY IONIZED IRON.
! D1
         metal_diffusion_coeff1_mid(zone_index)=interp_factors(1)* &
              src_grid_metal_diffusion_coeff1(k0)+interp_factors(2)* &
              src_grid_metal_diffusion_coeff1(k0+1)+interp_factors(3)* &
              src_grid_metal_diffusion_coeff1(k0+2)+interp_factors(4)* &
              src_grid_metal_diffusion_coeff1(k0+3)
! D2
         metal_diffusion_coeff2_mid(zone_index)=interp_factors(1)* &
              src_grid_metal_diffusion_coeff2(k0)+interp_factors(2)* &
              src_grid_metal_diffusion_coeff2(k0+1)+interp_factors(3)* &
              src_grid_metal_diffusion_coeff2(k0+2)+interp_factors(4)* &
              src_grid_metal_diffusion_coeff2(k0+3)
! DERIVATIVE OF D1 WRT Z
         eq_metal_diffusion_coeff1_mid(zone_index)=interp_factors(1)* &
              src_grid_metal_diffusion_coeff1_dz(k0)+interp_factors(2)* &
              src_grid_metal_diffusion_coeff1_dz(k0+1) &
              +interp_factors(3)*src_grid_metal_diffusion_coeff1_dz(k0+2)+ &
              interp_factors(4)*src_grid_metal_diffusion_coeff1_dz(k0+3)
! DERIVATIVE OF D2 WRT Z
         eq_metal_diffusion_coeff2_mid(zone_index)=interp_factors(1)* &
              src_grid_metal_diffusion_coeff2_dz(k0)+interp_factors(2)* &
              src_grid_metal_diffusion_coeff2_dz(k0+1) &
              +interp_factors(3)*src_grid_metal_diffusion_coeff2_dz(k0+2)+ &
              interp_factors(4)*src_grid_metal_diffusion_coeff2_dz(k0+3)
! MASS FRACTION OF METALS
         metal_abundance_change_mid(zone_index)=interp_factors(1)* &
              composition(8,k0)+interp_factors(2)*composition(8,k0+1)+ &
              interp_factors(3)*composition(8,k0+2)+ &
              interp_factors(4)*composition(8,k0+3)
         endif
   50 continue
! SET UP VECTOR OF EQUALLY SPACED RADII AT ZONE CENTERS.

      equal_radius(1)=radius(zone_begin)
      do 60 zone_index = 2,num_equal_points
         equal_radius(zone_index)=equal_radius(zone_index-1)+grid_spacing
   60 continue


! NOW USE 4-POINT LAGRANGIAN INTERPOLATION TO FIND RUN OF VARIABLES
! AT EQUALLY SPACED ZONE CENTERS.
!
! FIRST POINT : BY DEFINITION, AT STARTING POINT.
      equal_mass(1) = enclosed_mass(zone_begin)
      equal_hydrogen_fraction(1) = composition(1,zone_begin)
      equal_diffusion_coeff1(1)=diffusion_coeff1(zone_begin)
      if(use_diffusion_z)then
         metal_abundance_change(1) = composition(8,zone_begin)
         metal_diffusion_coeff1(1) = src_grid_metal_diffusion_coeff1(zone_begin)
      endif
! FOR OTHER POINTS: FIRST FIND 4 NEAREST (IN RADIUS) MODEL POINTS
! AND THEN FIND LAGRANGIAN INTERPOLATION FACTORS. APPLY THEM TO FIND
! MODEL QUANTITIES AT THE EQUALLY SPACED GRID POINTS.
! JMIN IS THE UPPERMOST MODEL POINT ABOVE THE PREVIOUS EQUALLY SPACED
! GRID POINT (IN RADIUS).
      search_start_index=zone_begin+1
      do 70 zone_index=2,num_equal_points-1
         do 80 j = search_start_index,zone_end
! FIND 4 MODEL POINTS CLOSEST TO THE EQUALLY SPACED GRID POINT.
            if(radius(j).ge.equal_radius(zone_index))then
! ENSURE THAT FIRST INTERPOLATION POINT NO LESS THAN FIRST MODEL POINT.
               k0 = max(j-2,1)
! ENSURE THAT LAST INTERPOLATION POINT NO GREATER THAN LAST MODEL POINT.
               k0 = min(k0,num_zones-3)
               search_start_index=j
               goto 90
            endif
   80    continue
         k0 = num_zones-3
         search_start_index=num_zones
   90    continue
         do 100 k=1,4
            radius_table(k)=radius(k0+k-1)
  100    continue

         target_radius=equal_radius(zone_index)

! FIND 4 POINT LAGRANGIAN INTERPOLATION FACTORS.
! FACINTERP=INTERPOLATION FACTORS FOR POINT GRIDRAD GIVEN THE 4 TABLE
! RADII IN TABLER; FACDERIV=SAME FOR DERIVATIVES AT POINT GRIDRAD.
         call interp(radius_table,interp_factors,deriv_factors,target_radius)
! PERFORM 4 POINT LAGRANGIAN INTERPOLATION FOR DESIRED QUANTITIES:
! MASS WITHIN THE RADIUS ER
         equal_mass(zone_index) = interp_factors(1)*enclosed_mass(k0)+ &
              interp_factors(2)*enclosed_mass(k0+1)+ &
              interp_factors(3)*enclosed_mass(k0+2)+ &
              interp_factors(4)*enclosed_mass(k0+3)
! D1
         equal_diffusion_coeff1(zone_index) = interp_factors(1)* &
              diffusion_coeff1(k0)+interp_factors(2)*diffusion_coeff1(k0+1)+ &
              interp_factors(3)*diffusion_coeff1(k0+2)+ &
              interp_factors(4)*diffusion_coeff1(k0+3)
! MASS FRACTION OF HYDROGEN

         equal_hydrogen_fraction(zone_index)=interp_factors(1)* &
              composition(1,k0)+interp_factors(2)*composition(1,k0+1)+ &
              interp_factors(3)*composition(1,k0+2)+ &
              interp_factors(4)*composition(1,k0+3)
! METAL DIFFUSION
         if(use_diffusion_z)then
! D1
         metal_diffusion_coeff1(zone_index)=interp_factors(1)* &
              src_grid_metal_diffusion_coeff1(k0)+interp_factors(2)* &
              src_grid_metal_diffusion_coeff1(k0+1)+interp_factors(3)* &
              src_grid_metal_diffusion_coeff1(k0+2)+interp_factors(4)* &
              src_grid_metal_diffusion_coeff1(k0+3)
! MASS FRACTION OF METALS
         metal_abundance_change(zone_index)=interp_factors(1)* &
              composition(8,k0)+interp_factors(2)*composition(8,k0+1)+ &
              interp_factors(3)*composition(8,k0+2)+ &
              interp_factors(4)*composition(8,k0+3)
         endif
   70 continue
! LAST POINT : BY DEFINITION, AT ENDING POINT.

      equal_mass(num_equal_points) = enclosed_mass(zone_end)
      equal_hydrogen_fraction(num_equal_points) = composition(1,zone_end)
      equal_diffusion_coeff1(num_equal_points)=diffusion_coeff1(zone_end)

      if(use_diffusion_z)then
         metal_abundance_change(num_equal_points) = composition(8,zone_end)
         metal_diffusion_coeff1(num_equal_points)=src_grid_metal_diffusion_coeff1(zone_end)
      endif
      return
end subroutine model_to_equal

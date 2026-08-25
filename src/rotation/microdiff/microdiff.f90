!----------------------------------------------------------------------
! microdiff
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original microdiff.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
!  GARRETT SOMERS - 05/2015
!  MICRODIFF IS AN UPDATED HELIUM, METAL, AND LIGHT ELEMENT DIFFUSION
!  ROUTINE. MUCH IS ADAPTED FROM GRSETT.F, WHICH IS RETAINED FOR BACKWARDS
!  COMPATABILITY. FOLLOWS THE FORMULISM OF BAHCALL AND LOEB 1989.
!
!  MICRODIFF PERFORMS THE FOLLOWING OPERATIONS:
!  1) TRANSFORMS PYHSICAL VARIABLES TO AN EQUALLY SPACED GRID IN RADIUS
!     (AND MANIPULATES THE MODEL VARIABLES INTO BAHCALL AND LOEB UNITS).
!  2) CALCULATES DIFFUSION COEFFICIENTS ALONG THE GRID FOR Y, Z, AND FOR
!     LIGHT ELEMENTS. THE DIFFUSION EQUATION HAS TWO TERMS : ONE THAT DEPENDS
!     ON DLNP/DR*DX/DR AND ONE THAT DEPENDS ON D^2X/DR^2. THE DIFFUSION
!     COEFFICIENTS ARE THEMSELVES FUNCTIONS OF BOTH THE THERMAL STRUCTURE,
!     ABUNDANCE, AND THE HYDROGEN GRADIENT.
!  3) SOLVES THIS DIFFUSION EQUATION IN TWO STEPS. THE FIRST TERM IS SOLVED
!     EXPLICITLY USING THE TWO-STEP LAX-WENDROFF SCHEME (NUMERICAL RECIPES,
!     PRESS ET AL. 1986,CAMBRIDGE UNIVERSITY PRESS,P.633). DIFFUSION
!     COEFFICIENTS ARE RECALCUALTED AT THE GRID MID-POINTS IN THE PROCESS.
!  4) THE SECOND TERM IS SOLVED IMPLICITLY, USING THE MODIFIED RUN OF
!     ABUNDANCE FOUND IN STEP 3) AS THE 'INITIAL' RUN OF ABUNDANCE.
!     THE IMPLICIT SCHEME ITERATES ON THE DIFFUSION COEFFICIENTS UNTIL
!     THE SOLUTION CONVERGES TO WITHIN THE USER-SPECIFIED TOLERANCE GRTOL.
!  5) UPDATES ABUNDANCE ARRAYS, TRANSFORMS BACK TO THE MODEL GRID, AND EXITS
!
!  timestep (DT) - TIMESTEP (SEC)
!  composition (HCOMP) - RUN OF MASS FRACTIONS OF SPECIES.HCOMP(1,..)=X,
!     HCOMP(2,...)=Y HCOMP(3,...)=Z
!  dlnp_dr (HQPR) - DLNP/DR
!  log_radius (HR) - LOG RADIUS (CM)
!  enclosed_mass (HS1) - MASS (GM), UNLOGGED
!  log_temperature (HT) - LOG TEMPERATURE (K)
!  convective_flag (LC) - FLAG T/F FOR CONVECTION
!  num_zones (M) - NUMBER OF MODEL POINTS
!  [COMMON] SDEL(2,...) - DEL (=DLNT/DLNP)
!
!  OUTPUT VARIABLES :
!
!  NEW RUN OF HCOMP(1-N,...).
!
!  CONVERT MODEL QUANTITIES TO BAHCALL AND LOEB UNITS; LOCATE
!  BOUNDARIES OF CONVECTIVE CORE AND ENVELOPE IF APPLICABLE;
!  COLLECT PHYSICAL VARIABLES FOR DIFFCO CALCULATION.
!
! Dispatcher for the microdiff.f90 element-settling pipeline: calls
! microdiff_setup.f90 (locate diffusion region, unit conversion),
! microdiff_mte.f90 (build the equally spaced grid), then
! microdiff_run.f90 (via microdiff_cod.f90) once each for hydrogen,
! heavy metals, and each light element in turn, and finally
! microdiff_etm.f90 (transform back to the model grid).
subroutine microdiff(timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, convective_flag, &
     num_zones, total_mass)
      use star_info_lib, only: star
      use star_info_lib, only: star, json
      use luout_lib
      implicit none
! SET NLIGHT TO THE NUMBER OF LIGHT ELEMENTS TO BE DIFFUSED.
      integer, parameter :: num_light = 3

      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass






      double precision :: radius_bl(json), temperature_bl(json)
      double precision :: density_orig(json), temperature_orig(json)
      double precision :: light_element_weight(num_light), &
           light_element_charge(num_light)
      integer :: light_element_id(num_light)
      double precision :: eq_mass(json), eq_radius(json), eq_density(json), &
           eq_temperature(json), eq_dlnp_dr(json), eq_del_grad(json), &
           eq_hydrogen(json), eq_helium(json), eq_metal(json), &
           eq_light(num_light,json), species_fraction(3,json), &
           hydrogen_dlnc_dr(json)
      double precision :: eq_mass_mid(json), eq_radius_mid(json), &
           eq_density_mid(json), eq_temperature_mid(json), &
           eq_dlnp_dr_mid(json), eq_del_grad_mid(json), &
           eq_hydrogen_mid(json), eq_helium_mid(json), eq_metal_mid(json), &
           eq_light_mid(num_light,json), species_fraction_mid(3,json), &
           hydrogen_dlnc_dr_mid(json), eq_delta_hydrogen(json), &
           eq_delta_metal(json), eq_delta_light(num_light,json)
!
!     SET UP VECTORS FOR LIGHT ELEMENT DIFFUSION. THE LENGTH
!     OF EACH SHOULD EQUAL NLIGHT.
      data light_element_id/13,14,15/
      data light_element_weight/6.015d0,7.016d0,9.012d0/
      data light_element_charge/3.0d0,3.0d0,4.0d0/
      integer :: i, ii, zone_begin, zone_end, num_eq_points, species_col
      logical :: fully_convective_flag
      double precision :: grid_spacing, atomic_weight_diffused, &
           atomic_charge_diffused

!  CONVERT MODEL QUANTITIES TO BAHCALL AND LOEB UNITS; LOCATE
!  BOUNDARIES OF CONVECTIVE CORE AND ENVELOPE IF APPLICABLE;
!  COLLECT PHYSICAL VARIABLES FOR DIFFCO CALCULATION.
!
      call microdiff_setup(timestep, dlnp_dr, log_radius, log_density, &
           enclosed_mass, log_temperature, convective_flag, num_zones, &
           total_mass, composition, radius_bl, temperature_bl, zone_begin, &
           zone_end, fully_convective_flag, density_orig, temperature_orig)
!
! SKIP SETTLING FOR FULLY CONVECTIVE MODELS.
      if(fully_convective_flag) return
!  TRANSFORM TO AN EQUALLY SPACED GRID IN RADIUS.
!  NOTE : PREFIX E ALONE=VARIABLE AT EQUALLY SPACED GRID POINTS.
!         PREFIX E + SUFFIX _H= VARIABLE AT MIDPOINT BETWEEN EQUALLY
!         SPACED GRID POINTS.(BOTH ARE NEEDED FOR THE DIFFUSION TECHNIQUE).
!
      call microdiff_mte(num_light, light_element_id, composition, &
           dlnp_dr, radius_bl, enclosed_mass, zone_begin, zone_end, &
           num_zones, grid_spacing, num_eq_points, density_orig, &
           temperature_orig, eq_mass, eq_radius, eq_density, &
           eq_temperature, eq_dlnp_dr, eq_del_grad, eq_hydrogen, eq_helium, &
           eq_metal, eq_light, eq_mass_mid, eq_radius_mid, eq_density_mid, &
           eq_temperature_mid, eq_dlnp_dr_mid, eq_del_grad_mid, &
           eq_hydrogen_mid, eq_helium_mid, eq_metal_mid, eq_light_mid)
!
!  PARTIALLY CONSTRUCT THE VECTORS ESPEC AND ESPEC_H, WHICH ARE THE
!  THREE SPECIES USED IN THE THOUL CALCULATION. 1=X, 2=Y, 3=METAL.
      do i=1,num_eq_points-1
         species_fraction(1,i) = eq_hydrogen(i)
         species_fraction(2,i) = eq_helium(i)
         species_fraction_mid(1,i) = eq_hydrogen_mid(i)
         species_fraction_mid(2,i) = eq_helium_mid(i)
      enddo
      species_fraction(1,num_eq_points) = eq_hydrogen(num_eq_points)
      species_fraction(2,num_eq_points) = eq_helium(num_eq_points)
!
!----------------------------------------------------------------------
!
!  WITH THE EQUAL GRID CALCULATED, RUN THE DIFFUSION ROUTINE FOR EACH
!  DESIRED ELEMENT. PASS IN EQUAL GRID AND THREE PARAMETERS:
!     J = THOUL COLUMN OF DIFFUSED SPECIES (1 FOR X, 3 OTHERWISE)
!     THEA = ATOMIC WEIGHT OF DIFFUSED SPECIES
!     THEZ = ATOMIC CHARGE OF DIFFUSED SPECIES
!  IN THE CASE OF HYDROGEN DIFFUSION, PASS IN THEA/THEZ FOR IRON SINCE
!  THESE VALUES ARE USED FOR THE METAL ABUNDANCE IN THOUL.
!
!     DIFFUSE HYDROGEN.
      if(star%job%diffuse_helium_active)then
         species_col = 1
         atomic_weight_diffused = 55.86d0
         atomic_charge_diffused = 26.0d0
!        SET ESPEC(3,*) TO THE HEAVY METAL ABUNDANCE.
         do i=1,num_eq_points-1
            species_fraction(3,i) = eq_metal(i)
            species_fraction_mid(3,i) = eq_metal_mid(i)
         enddo
         species_fraction(3,num_eq_points) = eq_metal(num_eq_points)
!        PASS THE EVEN GRID INTO MICRODIFF_RUN TO PERFORM THE DIFFUSION.
         call microdiff_run(grid_spacing, timestep, total_mass, &
              num_eq_points, eq_mass, eq_radius, eq_density, &
              eq_temperature, eq_dlnp_dr, eq_del_grad, species_fraction, &
              hydrogen_dlnc_dr, eq_mass_mid, eq_radius_mid, eq_density_mid, &
              eq_temperature_mid, eq_dlnp_dr_mid, eq_del_grad_mid, &
              species_fraction_mid, hydrogen_dlnc_dr_mid, &
              atomic_weight_diffused, atomic_charge_diffused, species_col)
!        STORE THE RUN OF CHANGES TO HYDROGEN
         do i = 1,num_eq_points
            eq_delta_hydrogen(i) = species_fraction(1,i)
         enddo
!        RESTORE THE ORIGINAL ESPEC HYDROGEN VECTORS.
         do i = 1,num_eq_points-1
            species_fraction(1,i) = eq_hydrogen(i)
            species_fraction_mid(1,i) = eq_hydrogen_mid(i)
         enddo
         species_fraction(1,num_eq_points) = eq_hydrogen(num_eq_points)
      endif
!
!----------------------------------------------------------------------
!
!     DIFFUSE HEAVY METALS.
      if(star%job%use_diffusion_z)then
         species_col = 3
         atomic_weight_diffused = 55.86d0
         atomic_charge_diffused = 26.0d0
!        SET ESPEC(3,*) TO THE HEAVY METAL ABUNDANCE.
         do i=1,num_eq_points-1
            species_fraction(3,i) = eq_metal(i)
            species_fraction_mid(3,i) = eq_metal_mid(i)
         enddo
         species_fraction(3,num_eq_points) = eq_metal(num_eq_points)
!        PASS THE EVEN GRID INTO MICRODIFF_RUN TO PERFORM THE DIFFUSION.
         call microdiff_run(grid_spacing, timestep, total_mass, &
              num_eq_points, eq_mass, eq_radius, eq_density, &
              eq_temperature, eq_dlnp_dr, eq_del_grad, species_fraction, &
              hydrogen_dlnc_dr, eq_mass_mid, eq_radius_mid, eq_density_mid, &
              eq_temperature_mid, eq_dlnp_dr_mid, eq_del_grad_mid, &
              species_fraction_mid, hydrogen_dlnc_dr_mid, &
              atomic_weight_diffused, atomic_charge_diffused, species_col)
!        STORE THE RUN OF CHANGES TO HEAVY METALS
         do i=1,num_eq_points
            eq_delta_metal(i) = species_fraction(3,i)
         enddo
!        THE ORIGINAL HEAVY VECTOR DOESN'T NEED RESTORATION, SINCE IT WILL
!        BE OVERWRITTEN BELOW.
      endif
!
!----------------------------------------------------------------------
!
!     DIFFUSE LIGHT ELEMENTS.
      if(star%ctrl%ldifli)then
!        ITERATE OVER THE DIFFUSION ROUTINES FOR EACH LIGHT ELEMENT.
         species_col = 3
         do ii = 1,num_light
            atomic_weight_diffused = light_element_weight(ii)
            atomic_charge_diffused = light_element_charge(ii)
!           SET ESPEC(3,*) TO THE LIGHT ELEMENT ABUNDANCE.
            do i=1,num_eq_points-1
               species_fraction(3,i) = eq_light(ii,i)
               species_fraction_mid(3,i) = eq_light_mid(ii,i)
            enddo
            species_fraction(3,i) = eq_light(ii,num_eq_points)
!           PASS THE EVEN GRID INTO MICRODIFF_RUN TO PERFORM THE DIFFUSION.
            call microdiff_run(grid_spacing, timestep, total_mass, &
                 num_eq_points, eq_mass, eq_radius, eq_density, &
                 eq_temperature, eq_dlnp_dr, eq_del_grad, species_fraction, &
                 hydrogen_dlnc_dr, eq_mass_mid, eq_radius_mid, &
                 eq_density_mid, eq_temperature_mid, eq_dlnp_dr_mid, &
                 eq_del_grad_mid, species_fraction_mid, hydrogen_dlnc_dr_mid, &
                 atomic_weight_diffused, atomic_charge_diffused, species_col)
!           STORE THE NEW RUN OF LIGHT METALS
            do i=1,num_eq_points
              eq_delta_light(ii,i) = species_fraction(3,i)
            enddo
         enddo
      endif
!
!----------------------------------------------------------------------

! TRANSFORM BACK TO ORIGINAL MODEL GRID; UPDATE HELIUM ARRAY USING
! X+Y+Z=1.  PRINT DIAGNOSTIC OUTPUT.
!
      call microdiff_etm(timestep, eq_radius, eq_delta_hydrogen, &
           eq_delta_metal, eq_delta_light, zone_begin, zone_end, &
           num_eq_points, composition, dlnp_dr, radius_bl, enclosed_mass, &
           temperature_bl, num_zones, total_mass, num_light, light_element_id)
      return
end subroutine microdiff

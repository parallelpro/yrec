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
!  ROUTINE. MUCH IS ADAPTED FROM gravitational_settling.f90, WHICH IS
!  RETAINED FOR BACKWARDS COMPATABILITY. FOLLOWS THE FORMULISM OF
!  BAHCALL AND LOEB 1989.
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
!  del_grad - DEL (=DLNT/DLNP) per zone (2026 W2: an explicit argument;
!     mix passes star%gradT, rotmix passes mix_scr%delm, which it used
!     to copy into star%gradT around this call)
!  convective_flag (LC) - FLAG T/F FOR CONVECTION
!  num_zones (M) - NUMBER OF MODEL POINTS
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
! microdiff_run.f90 (via microdiff_coefficients.f90) once each for hydrogen,
! heavy metals, and each light element in turn, and finally
! microdiff_etm.f90 (transform back to the model grid).
subroutine microdiff(timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, &
     convective_flag, num_zones, total_mass, ierr)
      use microdiff_mte_lib
      use microdiff_run_lib
      use star_info_lib, only: star, json
      use species_table_lib, only: thoul_fe, thoul_col_h, thoul_col_he, &
           thoul_col_metal, num_light_diffused, light_diffused, light_diffused_row
      implicit none
! SET NLIGHT TO THE NUMBER OF LIGHT ELEMENTS TO BE DIFFUSED.
      integer, parameter :: num_light = num_light_diffused

      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json), del_grad(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      integer, intent(out) :: ierr
! --- locals ---
      double precision :: radius_bl(json), temperature_bl(json)
      double precision :: density_orig(json), temperature_orig(json)
! the equally spaced grids at zone centers (eq) and zone midpoints
! (eq_mid), built by microdiff_mte -- see microdiff_mte_lib
      type(microdiff_grid) :: eq, eq_mid
      double precision :: species_fraction(3,json), &
           hydrogen_dlnc_dr(json)
      double precision :: species_fraction_mid(3,json), &
           hydrogen_dlnc_dr_mid(json), eq_delta_hydrogen(json), &
           eq_delta_metal(json), eq_delta_light(num_light,json)
!
!     THE LIGHT ELEMENTS DIFFUSED (ROW, WEIGHT, CHARGE) ARE
!     light_diffused_row / light_diffused OF species_table_lib.
      integer :: i, ii, zone_begin, zone_end, num_eq_points, species_col
      logical :: settling_skipped_flag
      double precision :: grid_spacing, atomic_weight_diffused, &
           atomic_charge_diffused

!  CONVERT MODEL QUANTITIES TO BAHCALL AND LOEB UNITS; LOCATE
!  BOUNDARIES OF CONVECTIVE CORE AND ENVELOPE IF APPLICABLE;
!  COLLECT PHYSICAL VARIABLES FOR DIFFCO CALCULATION.
!
      ierr = 0
      call microdiff_setup(timestep, dlnp_dr, log_radius, log_density, &
           enclosed_mass, log_temperature, convective_flag, num_zones, &
           total_mass, composition, radius_bl, temperature_bl, zone_begin, &
           zone_end, settling_skipped_flag, density_orig, temperature_orig)
!
! SKIP SETTLING FOR FULLY CONVECTIVE (OR H/HE-EXHAUSTED) MODELS.
      if(settling_skipped_flag) return
!  TRANSFORM TO AN EQUALLY SPACED GRID IN RADIUS.
!  NOTE : PREFIX E ALONE=VARIABLE AT EQUALLY SPACED GRID POINTS.
!         PREFIX E + SUFFIX _H= VARIABLE AT MIDPOINT BETWEEN EQUALLY
!         SPACED GRID POINTS.(BOTH ARE NEEDED FOR THE DIFFUSION TECHNIQUE).
!
      call microdiff_mte(num_light, light_diffused_row, composition, &
           dlnp_dr, radius_bl, enclosed_mass, del_grad, zone_begin, zone_end, &
           num_zones, grid_spacing, num_eq_points, density_orig, &
           temperature_orig, eq, eq_mid)
!
!  PARTIALLY CONSTRUCT THE VECTORS ESPEC AND ESPEC_H, WHICH ARE THE
!  THREE SPECIES USED IN THE THOUL CALCULATION. 1=X, 2=Y, 3=METAL.
      do i=1,num_eq_points-1
         species_fraction(thoul_col_h,i) = eq%hydrogen(i)
         species_fraction(thoul_col_he,i) = eq%helium(i)
         species_fraction_mid(thoul_col_h,i) = eq_mid%hydrogen(i)
         species_fraction_mid(thoul_col_he,i) = eq_mid%helium(i)
      enddo
      species_fraction(thoul_col_h,num_eq_points) = eq%hydrogen(num_eq_points)
      species_fraction(thoul_col_he,num_eq_points) = eq%helium(num_eq_points)
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
         species_col = thoul_col_h
         atomic_weight_diffused = thoul_fe%weight
         atomic_charge_diffused = thoul_fe%charge
!        SET ESPEC(3,*) TO THE HEAVY METAL ABUNDANCE.
         do i=1,num_eq_points-1
            species_fraction(thoul_col_metal,i) = eq%metal(i)
            species_fraction_mid(thoul_col_metal,i) = eq_mid%metal(i)
         enddo
         species_fraction(thoul_col_metal,num_eq_points) = eq%metal(num_eq_points)
!        PASS THE EVEN GRID INTO MICRODIFF_RUN TO PERFORM THE DIFFUSION.
         call microdiff_run(grid_spacing, timestep, total_mass, &
              num_eq_points, eq, species_fraction, &
              hydrogen_dlnc_dr, eq_mid, &
              species_fraction_mid, hydrogen_dlnc_dr_mid, &
              atomic_weight_diffused, atomic_charge_diffused, species_col, ierr)
         if (ierr /= 0) return
!        STORE THE RUN OF CHANGES TO HYDROGEN
         do i = 1,num_eq_points
            eq_delta_hydrogen(i) = species_fraction(thoul_col_h,i)
         enddo
!        RESTORE THE ORIGINAL ESPEC HYDROGEN VECTORS.
         do i = 1,num_eq_points-1
            species_fraction(thoul_col_h,i) = eq%hydrogen(i)
            species_fraction_mid(thoul_col_h,i) = eq_mid%hydrogen(i)
         enddo
         species_fraction(thoul_col_h,num_eq_points) = eq%hydrogen(num_eq_points)
      endif
!
!----------------------------------------------------------------------
!
!     DIFFUSE HEAVY METALS.
      if(star%job%use_diffusion_z)then
         species_col = thoul_col_metal
         atomic_weight_diffused = thoul_fe%weight
         atomic_charge_diffused = thoul_fe%charge
!        SET ESPEC(3,*) TO THE HEAVY METAL ABUNDANCE.
         do i=1,num_eq_points-1
            species_fraction(thoul_col_metal,i) = eq%metal(i)
            species_fraction_mid(thoul_col_metal,i) = eq_mid%metal(i)
         enddo
         species_fraction(thoul_col_metal,num_eq_points) = eq%metal(num_eq_points)
!        PASS THE EVEN GRID INTO MICRODIFF_RUN TO PERFORM THE DIFFUSION.
         call microdiff_run(grid_spacing, timestep, total_mass, &
              num_eq_points, eq, species_fraction, &
              hydrogen_dlnc_dr, eq_mid, &
              species_fraction_mid, hydrogen_dlnc_dr_mid, &
              atomic_weight_diffused, atomic_charge_diffused, species_col, ierr)
         if (ierr /= 0) return
!        STORE THE RUN OF CHANGES TO HEAVY METALS
         do i=1,num_eq_points
            eq_delta_metal(i) = species_fraction(thoul_col_metal,i)
         enddo
!        THE ORIGINAL HEAVY VECTOR DOESN'T NEED RESTORATION, SINCE IT WILL
!        BE OVERWRITTEN BELOW.
      endif
!
!----------------------------------------------------------------------
!
!     DIFFUSE LIGHT ELEMENTS.
      if(star%ctrl%diffuse_lithium)then
!        ITERATE OVER THE DIFFUSION ROUTINES FOR EACH LIGHT ELEMENT.
         species_col = thoul_col_metal
         do ii = 1,num_light
            atomic_weight_diffused = light_diffused(ii)%weight
            atomic_charge_diffused = light_diffused(ii)%charge
!           SET ESPEC(3,*) TO THE LIGHT ELEMENT ABUNDANCE.
            do i=1,num_eq_points-1
               species_fraction(thoul_col_metal,i) = eq%light(ii,i)
               species_fraction_mid(thoul_col_metal,i) = eq_mid%light(ii,i)
            enddo
            species_fraction(thoul_col_metal,num_eq_points) = eq%light(ii,num_eq_points)
!           PASS THE EVEN GRID INTO MICRODIFF_RUN TO PERFORM THE DIFFUSION.
            call microdiff_run(grid_spacing, timestep, total_mass, &
                 num_eq_points, eq, species_fraction, &
                 hydrogen_dlnc_dr, eq_mid, &
                 species_fraction_mid, hydrogen_dlnc_dr_mid, &
                 atomic_weight_diffused, atomic_charge_diffused, species_col, ierr)
            if (ierr /= 0) return
!           STORE THE NEW RUN OF LIGHT METALS
            do i=1,num_eq_points
              eq_delta_light(ii,i) = species_fraction(thoul_col_metal,i)
            enddo
         enddo
      endif
!
!----------------------------------------------------------------------

! TRANSFORM BACK TO ORIGINAL MODEL GRID; UPDATE HELIUM ARRAY USING
! X+Y+Z=1.  PRINT DIAGNOSTIC OUTPUT.
!
      call microdiff_etm(timestep, eq%radius, eq_delta_hydrogen, &
           eq_delta_metal, eq_delta_light, zone_begin, zone_end, &
           num_eq_points, composition, dlnp_dr, radius_bl, enclosed_mass, &
           num_zones, total_mass, num_light, light_diffused_row)
      return
end subroutine microdiff

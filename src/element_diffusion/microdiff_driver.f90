!----------------------------------------------------------------------
! microdiff_driver
!----------------------------------------------------------------------
! The newer settling routines (use_new_diffusion_routines true): microdiff
! and its setup.
module microdiff_driver
      use element_diffusion_configuration, only: diffusion_settings
      use element_diffusion_scratch, only: diffusion_scratch
      implicit none
      private
      public :: microdiff

contains

! ---------------------------------------------------------------
!  GARRETT SOMERS - 05/2015
!  MICRODIFF IS AN UPDATED HELIUM, METAL, AND LIGHT ELEMENT DIFFUSION
!  ROUTINE. MUCH IS ADAPTED FROM gravitational_settling, WHICH IS
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
!  del_grad - DEL (=DLNT/DLNP) per zone (core passes star%gradT, burn_and_settle
!     its mid-step gradient)
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
! The pipeline: microdiff_setup (the diffusion region, the unit
! conversion), map_microdiff_to_equal_grid (microdiff_mapping: the
! equally spaced grid), run_microdiff (microdiff_run, with the
! coefficients from microdiff_coefficients) once each for hydrogen, the
! metals and each light element, and microdiff_etm (microdiff_mapping:
! back to the model grid).
subroutine microdiff(diff_scr, timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, &
     convective_flag, num_zones, total_mass, diffuse_helium, &
     suspension_reported, ierr)
      use microdiff_mapping
      use microdiff_run
      use composition_layout, only: json
      use species_table, only: thoul_fe, thoul_col_h, thoul_col_he, &
           thoul_col_metal, num_light_diffused, light_diffused, light_diffused_row
      implicit none
      type(diffusion_scratch), intent(inout) :: diff_scr
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
      ! settle hydrogen and helium; the run log's suspension latch
      logical, intent(in) :: diffuse_helium
      logical, intent(inout) :: suspension_reported
! --- locals ---
      double precision :: radius_bl(json), temperature_bl(json)
      double precision :: density_orig(json), temperature_orig(json)
! the equally spaced grids at zone centers (eq) and zone midpoints
! (eq_mid), built by microdiff_mte -- see microdiff_mte
      type(microdiff_grid) :: eq, eq_mid
      double precision :: species_fraction(3,json), &
           hydrogen_dlnc_dr(json)
      double precision :: species_fraction_mid(3,json), &
           hydrogen_dlnc_dr_mid(json), eq_delta_hydrogen(json), &
           eq_delta_metal(json), eq_delta_light(num_light,json)
!
!     THE LIGHT ELEMENTS DIFFUSED (ROW, WEIGHT, CHARGE) ARE
!     light_diffused_row / light_diffused OF species_table.
      integer :: i, ii, zone_begin, zone_end, num_eq_points, species_col
      logical :: settling_skipped_flag
      double precision :: grid_spacing, atomic_weight_diffused, &
           atomic_charge_diffused

!  CONVERT MODEL QUANTITIES TO BAHCALL AND LOEB UNITS; LOCATE
!  BOUNDARIES OF CONVECTIVE CORE AND ENVELOPE IF APPLICABLE;
!  COLLECT PHYSICAL VARIABLES FOR DIFFCO CALCULATION.
!
      ierr = 0
      call microdiff_setup(diff_scr, timestep, dlnp_dr, log_radius, log_density, &
           enclosed_mass, log_temperature, convective_flag, num_zones, &
           total_mass, composition, radius_bl, temperature_bl, zone_begin, &
           zone_end, settling_skipped_flag, density_orig, temperature_orig, &
           suspension_reported)
!
! SKIP SETTLING FOR FULLY CONVECTIVE (OR H/HE-EXHAUSTED) MODELS.
      if(settling_skipped_flag) return
!  TRANSFORM TO AN EQUALLY SPACED GRID IN RADIUS.
!  NOTE : PREFIX E ALONE=VARIABLE AT EQUALLY SPACED GRID POINTS.
!         PREFIX E + SUFFIX _H= VARIABLE AT MIDPOINT BETWEEN EQUALLY
!         SPACED GRID POINTS.(BOTH ARE NEEDED FOR THE DIFFUSION TECHNIQUE).
!
      call map_microdiff_to_equal_grid(num_light, light_diffused_row, composition, &
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
      if(diffuse_helium)then
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
         call run_microdiff(diff_scr, grid_spacing, timestep, total_mass, &
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
      if(diffusion_settings%use_diffusion_z)then
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
         call run_microdiff(diff_scr, grid_spacing, timestep, total_mass, &
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
      if(diffusion_settings%diffuse_lithium)then
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
            call run_microdiff(diff_scr, grid_spacing, timestep, total_mass, &
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
      call microdiff_etm(diff_scr, timestep, eq%radius, eq_delta_hydrogen, &
           eq_delta_metal, eq_delta_light, zone_begin, zone_end, &
           num_eq_points, composition, dlnp_dr, radius_bl, enclosed_mass, &
           num_zones, total_mass, num_light, light_diffused_row)
      return
end subroutine microdiff

! ---------------------------------------------------------------
! The first stage of microdiff: locates the diffusion region (between any central
! convective core and any surface convection zone, and bounded by
! hydrogen/helium exhaustion), and converts the needed model quantities
! to Bahcall & Loeb (1990) natural units.
!
!     OUTPUT VARIABLES :
!
!     radius_bl - VECTOR OF UNLOGGED RADII IN BAHCALL AND LOEB UNITS.
!     temperature_bl - VECTOR OF UNLOGGED TEMPERATURES IN BAHCALL AND
!          LOEB UNITS.
!     zone_begin - FIRST ZONE FOR DIFFUSION PURPOSES (EITHER THE FIRST
!          MODEL POINT OR THE OUTERMOST POINT OF A CENTRAL CONVECTION
!          ZONE).
!     zone_end - LAST ZONE FOR DIFFUSION PURPOSES (EITHER THE LAST
!          MODEL POINT OR THE INNERMOST POINT OF A SURFACE CONVECTION
!          ZONE).
!     THE VECTORS dlnp_dr AND enclosed_mass, AND THE SCALARS timestep
!     AND total_mass, ARE ALSO CONVERTED TO BAHCALL AND LOEB UNITS.
!     CONSTANTS DEFINED :
!     ln10 = CONVERSION FACTOR FROM LN TO LOG10
!     CRSUN_BAH = SOLAR RADIUS (CM)
!     CSECYR_BAH = NUMBER OF SECONDS IN A YEAR.
!
! enclosed_mass (HS1) is the enclosed mass at the model points,
! unlogged, not the per-shell mass.
subroutine microdiff_setup(diff_scr, timestep, dlnp_dr, log_radius, log_density, &
     enclosed_mass, log_temperature, convective_flag, num_zones, &
     total_mass, composition, radius_bl, temperature_bl, zone_begin, &
     zone_end, settling_skipped_flag, density_orig, temperature_orig, &
     suspension_reported)

      use composition_layout, only: json, i_h1, i_he4
      use bahcall_loeb_units, only: set_bahcall_loeb_scales
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use math_lib
      implicit none
      type(diffusion_scratch), intent(inout) :: diff_scr

      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      double precision, intent(in) :: composition(15,json)
      double precision, intent(out) :: radius_bl(json), temperature_bl(json)
      integer, intent(out) :: zone_begin, zone_end
      logical, intent(out) :: settling_skipped_flag
      double precision, intent(out) :: density_orig(json), &
           temperature_orig(json)
! the run log's latch: a suspension of settling is reported once
      logical, intent(inout) :: suspension_reported
      integer :: i

!     settling_skipped_flag=T IF SETTLING IS SKIPPED THIS MODEL (FULLY
!     CONVECTIVE, HYDROGEN-EXHAUSTED OR HELIUM-EXHAUSTED).
      settling_skipped_flag=.false.
!     CHECK FOR CONVECTIVE CORE.
      if(convective_flag(1))then
         do i=2,num_zones
            if(.not.convective_flag(i))exit
         end do
         if (i > num_zones) then
!        DIFFUSION NOT COMPUTED FOR FULLY CONVECTIVE MODELS.
         settling_skipped_flag=.true.
! print once per suspension; every model only under
! report_solver_diagnostics
         if (solver_diagnostics() .or. &
              .not. suspension_reported) then
            write(run_log_unit,15)
   15       format(1x,' FULLY CONVECTIVE MODEL - NO SETTLING')
            suspension_reported = .true.
         end if
         return
         end if
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_begin = i-1
      else
         zone_begin = 1
      endif
! MHP 6/90 CHECK FOR HYDROGEN-EXHAUSTED CORE.
      do i = zone_begin,num_zones
         if(composition(i_h1,i).gt.diffusion_settings%hydrogen_diffusion_floor)exit
      end do
      if (i > num_zones) then
!     HYDROGEN-FREE MODEL - EXIT.
! print once per suspension; every model only under
! report_solver_diagnostics
      if (solver_diagnostics() .or. &
           .not. suspension_reported) then
         write(run_log_unit,16)diffusion_settings%hydrogen_diffusion_floor
   16    format(1x,'X BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
         suspension_reported = .true.
      end if
      settling_skipped_flag = .true.
      return
      end if
      zone_begin = i
!     CHECK FOR CONVECTIVE ENVELOPE.
      if(convective_flag(num_zones))then
         do i=num_zones-1,2,-1
            if(.not.convective_flag(i))exit
         end do
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_end = i+1
      else
         zone_end = num_zones
      endif
!     CHECK FOR HELIUM-EXHAUSTED SURFACE.
!     OUTER POINT IS SET WHEREVER Y>YMIN.
      do i=zone_end,1,-1
         if(composition(i_he4,i).gt.diffusion_settings%helium_diffusion_min) exit
      end do
      if (i < 1) then
!     HELIUM-EXHAUSTED MODEL - EXIT.
! print once per suspension; every model only under
! report_solver_diagnostics
      if (solver_diagnostics() .or. &
           .not. suspension_reported) then
         write(run_log_unit,17)diffusion_settings%helium_diffusion_min
   17    format(1x,'Y BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
         suspension_reported = .true.
      end if
      settling_skipped_flag = .true.
      return
      end if
      zone_end = i
! all suspension checks passed: settling proceeds this model
      if (suspension_reported) then
         write(run_log_unit,916)
  916    format(1x,' SETTLING RESUMED')
         suspension_reported = .false.
      end if
!     SET THE star%bl_*_scale CONVERSION FACTORS.
      call set_bahcall_loeb_scales(diff_scr)
!     CONVERT LOG(RADIUS) AND LOG(TEMPERATURE) TO NATURAL UNITS.
!     ALSO CONVERT NATURAL UNITS TO BAHCALL AND LOEB UNITS.
      do i=1,num_zones
         radius_bl(i)=exp(ln10*log_radius(i))*diff_scr%bl_radius_scale
         temperature_bl(i)=exp(ln10*log_temperature(i))*diff_scr%bl_temp_scale
         enclosed_mass(i)=enclosed_mass(i)*diff_scr%bl_mass_scale
         dlnp_dr(i)=dlnp_dr(i)/diff_scr%bl_radius_scale
      end do
      timestep=timestep/diff_scr%bl_time_scale
      total_mass=total_mass*diff_scr%bl_mass_scale
!
! COLLECT THE NECESSARY QUANTITIES (NAMELY RHO AND T) FOR LATER
! TRANSFORMATION TO THE EQUALLY SPACED GRID.
      do i = 1,num_zones
         density_orig(i) = exp(ln10*log_density(i))
         temperature_orig(i) = exp(ln10*log_temperature(i))
      enddo
      return
end subroutine microdiff_setup

end module microdiff_driver

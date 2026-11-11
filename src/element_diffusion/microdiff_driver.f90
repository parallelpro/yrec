!----------------------------------------------------------------------
! microdiff_driver
!----------------------------------------------------------------------
! The newer settling routines (use_new_diffusion_routines true): microdiff
! and its setup.
module microdiff_driver
      use element_diffusion_configuration, only: diffusion_settings
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
! The pipeline (diffuse_elements has found the settling region):
! microdiff_setup (the unit conversion), map_microdiff_to_equal_grid (microdiff_mapping: the
! equally spaced grid), run_microdiff (microdiff_run, with the
! coefficients from microdiff_coefficients) once each for hydrogen, the
! metals and each light element, and microdiff_etm (microdiff_mapping:
! back to the model grid).
subroutine microdiff(timestep, composition, dlnp_dr, log_radius, &
     log_density, enclosed_mass, log_temperature, del_grad, &
     num_zones, total_mass, zone_begin, zone_end, ierr)
      use microdiff_mapping
      use microdiff_run
      use composition_layout, only: json
      use species_table, only: thoul_fe, thoul_col_h, thoul_col_he, &
           thoul_col_metal, num_light_diffused, light_diffused, light_diffused_row
      implicit none
! SET NLIGHT TO THE NUMBER OF LIGHT ELEMENTS TO BE DIFFUSED.
      integer, parameter :: num_light = num_light_diffused

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json)
      double precision, intent(in) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json), del_grad(json)
      integer, intent(in) :: num_zones
      double precision, intent(in) :: total_mass
! the settling region (settling_region)
      integer, intent(in) :: zone_begin, zone_end
      integer, intent(out) :: ierr
! --- locals ---
! The timestep, dlnP/dr, mass coordinate and total mass in Bahcall and
! Loeb units: copies, so that the model's own are left exactly as they
! were (they were scaled and scaled back in place, which could move
! their last bits; docs/bug_sweep.md O16).
      double precision :: timestep_bl, dlnp_dr_bl(json), enclosed_mass_bl(json), &
           total_mass_bl
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
      integer :: i, ii, num_eq_points, species_col
      double precision :: grid_spacing, atomic_weight_diffused, &
           atomic_charge_diffused
! what each species' solve did (microdiff_run), for its messages
      type(species_solve_report) :: report

!  CONVERT MODEL QUANTITIES TO BAHCALL AND LOEB UNITS; LOCATE
!  BOUNDARIES OF CONVECTIVE CORE AND ENVELOPE IF APPLICABLE;
!  COLLECT PHYSICAL VARIABLES FOR DIFFCO CALCULATION.
!
      ierr = 0
      timestep_bl = timestep
      dlnp_dr_bl = dlnp_dr
      enclosed_mass_bl = enclosed_mass
      total_mass_bl = total_mass
      call microdiff_setup(timestep_bl, dlnp_dr_bl, log_radius, log_density, &
           enclosed_mass_bl, log_temperature, num_zones, &
           total_mass_bl, radius_bl, temperature_bl, density_orig, temperature_orig)
!  TRANSFORM TO AN EQUALLY SPACED GRID IN RADIUS.
!  NOTE : PREFIX E ALONE=VARIABLE AT EQUALLY SPACED GRID POINTS.
!         PREFIX E + SUFFIX _H= VARIABLE AT MIDPOINT BETWEEN EQUALLY
!         SPACED GRID POINTS.(BOTH ARE NEEDED FOR THE DIFFUSION TECHNIQUE).
!
      call map_microdiff_to_equal_grid(num_light, light_diffused_row, composition, &
           dlnp_dr_bl, radius_bl, enclosed_mass_bl, del_grad, zone_begin, zone_end, &
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
!     DIFFUSE HYDROGEN (always: diffuse_elements is called only while
!     hydrogen and helium settle, and the metal and light-element calls
!     read the hydrogen gradient this call computes).
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
         call run_microdiff(grid_spacing, timestep_bl, total_mass_bl, &
              num_eq_points, eq, species_fraction, &
              hydrogen_dlnc_dr, eq_mid, &
              species_fraction_mid, hydrogen_dlnc_dr_mid, &
              atomic_weight_diffused, atomic_charge_diffused, species_col, report, ierr)
         call report_species_solve(report, ierr)
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
         call run_microdiff(grid_spacing, timestep_bl, total_mass_bl, &
              num_eq_points, eq, species_fraction, &
              hydrogen_dlnc_dr, eq_mid, &
              species_fraction_mid, hydrogen_dlnc_dr_mid, &
              atomic_weight_diffused, atomic_charge_diffused, species_col, report, ierr)
         call report_species_solve(report, ierr)
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
            call run_microdiff(grid_spacing, timestep_bl, total_mass_bl, &
                 num_eq_points, eq, species_fraction, &
                 hydrogen_dlnc_dr, eq_mid, &
                 species_fraction_mid, hydrogen_dlnc_dr_mid, &
                 atomic_weight_diffused, atomic_charge_diffused, species_col, report, ierr)
            call report_species_solve(report, ierr)
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
      call microdiff_etm(eq%radius, eq_delta_hydrogen, &
           eq_delta_metal, eq_delta_light, zone_begin, zone_end, &
           num_eq_points, composition, radius_bl, &
           num_zones, num_light, light_diffused_row)
      return
end subroutine microdiff

! ---------------------------------------------------------------
! The messages of one species' solve (element_diffusion_messages): the
! numerics routine that failed (ierr nonzero); else each pass of the
! implicit solve under solver diagnostics (pass 1's change, then zero),
! and a solve that missed settling_tolerance in its
! settling_num_iterations passes, with the last pass's change.
subroutine report_species_solve(report, ierr)
      use microdiff_run, only: species_solve_report
      use element_diffusion_messages, only: note_solve_pass, report_unconverged, &
           report_numerics_failure
      type(species_solve_report), intent(in) :: report
      integer, intent(in) :: ierr
      integer :: pass

      if (ierr /= 0) then
         call report_numerics_failure(trim(report%failed_routine), ierr)
         return
      end if
      do pass = 1, report%num_passes
         if (pass == 1) then
            call note_solve_pass(pass, report%change, report%change_shell)
         else
            call note_solve_pass(pass, 0.0d0, 0)
         end if
      end do
      if (report%converged) return
      if (report%num_passes <= 1) then
         call report_unconverged('MICRODIFF', 'D', diffusion_settings%settling_tolerance, &
              diffusion_settings%settling_num_iterations, report%change, report%change_shell)
      else
         call report_unconverged('MICRODIFF', 'D', diffusion_settings%settling_tolerance, &
              diffusion_settings%settling_num_iterations, 0.0d0, 0)
      end if
end subroutine report_species_solve

! ---------------------------------------------------------------
! The first stage of microdiff: the model's quantities in Bahcall &
! Loeb (1990) units (bahcall_loeb_units), on microdiff's copies.
!   timestep, dlnp_dr, enclosed_mass, total_mass  converted in place
!       (s, 1/cm, g, g in; Bahcall & Loeb units out);
!   radius_bl, temperature_bl   the unlogged radius and temperature in
!       Bahcall & Loeb units;
!   density_orig, temperature_orig  the unlogged density (g/cm^3) and
!       temperature (K), for the equally spaced grids.
subroutine microdiff_setup(timestep, dlnp_dr, log_radius, log_density, &
     enclosed_mass, log_temperature, num_zones, &
     total_mass, radius_bl, temperature_bl, density_orig, temperature_orig)

      use composition_layout, only: json
      use bahcall_loeb_units, only: bl_radius_scale, bl_mass_scale, &
           bl_temp_scale, bl_time_scale
      use physical_constants, only: ln10
      use math_lib
      implicit none
      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      double precision, intent(out) :: radius_bl(json), temperature_bl(json)
      double precision, intent(out) :: density_orig(json), &
           temperature_orig(json)
      integer :: i

!     CONVERT LOG(RADIUS) AND LOG(TEMPERATURE) TO NATURAL UNITS.
!     ALSO CONVERT NATURAL UNITS TO BAHCALL AND LOEB UNITS.
      do i=1,num_zones
         radius_bl(i)=exp(ln10*log_radius(i))*bl_radius_scale
         temperature_bl(i)=exp(ln10*log_temperature(i))*bl_temp_scale
         enclosed_mass(i)=enclosed_mass(i)*bl_mass_scale
         dlnp_dr(i)=dlnp_dr(i)/bl_radius_scale
      end do
      timestep=timestep/bl_time_scale
      total_mass=total_mass*bl_mass_scale
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

!----------------------------------------------------------------------
! ndifcom
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ndifcom.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Dummy-argument names for the quantities passed straight through to
! mixgrid.f90 (log_density/log_luminosity/log_pressure/log_radius/
! log_mass/enclosed_mass/log_total_mass) match that already-converted
! callee's own parameter names, and the whole argument list matches the
! already-converted caller bursmix.f90's call to ndifcom.
!
!  DIFCOM CALCULATES THE DIFFUSION OF COMPOSITION DUE TO ANGULAR MOMENTUM
!  TRANSPORT.  THIS IS DONE BY TRANSFORMING TO AN EQUALLY SPACED GRID IN
!  RADIUS, SOLVING A DIFFUSION EQUATION, AND TRANSFORMING BACK.
!
!  INPUT VARIABLES:
!  timestep (DT) : DIFFUSION TIMESTEP(SEC).
!  diffusion_coeff (COD2) : DIFFUSION COEFFICIENTS FOR COMPOSITION
!     TRANSPORT AT THE ORIGINAL MODEL POINTS. Passed straight through
!     to mixgrid.
!  equally_spaced_diffusion_coeff (ECOD2) : DIFFUSION COEFFICIENTS FOR
!     COMPOSITION TRANSPORT AT THE EQUALLY SPACED GRID POINTS.
!  equally_spaced_mass (EM) : MASSES OF THE EQUALLY SPACED GRID
!     POINTS(GM).
!     NOTE: FOR CONVECTIVE BOUNDARIES THE MASS OF THE LAST GRID POINT IS
!           THE MASS OF THE ENTIRE CONVECTION ZONE.
!  composition (HCOMP) : ARRAY OF MASS FRACTION OF ALL OF THE SPECIES AT
!     THE ORIGINAL MODEL POINTS.
!  enclosed_mass (HS1) : LOCATION IN MASS(UNLOGGED) OF THE ORIGINAL MODEL
!     POINTS.
!  shell_mass (HS2) : MASSES OF THE ORIGINAL MODEL POINTS(UNLOGGED).
!  velocity (HV) : RUN OF DIFFUSION VELOCITIES.
!  zone_begin,zone_end (IBEG,IEND) :THE FIRST/LAST UNSTABLE POINTS IN
!     THE REGION.
!     NOTE: FOR CONVECTIVE BOUNDARIES THESE ARE ONLY THE FIRST CONVECTIVE
!           POINTS ADJACENT TO AN UNSTABLE RADIATIVE REGION.
!  zone_min (IMIN) : THE INNERMOST RADIATIVE ZONE OUTSIDE OF ZONE 1.
!  zone_max (IMAX) : the outermost radiative zone.
!  convective_flag (LCZ) : FLAG WHICH TELLS WHICH OF THE ORIGINAL MODEL
!     POINTS ARE CONVECTIVE FOR ANGULAR MOMENTUM TRANSPORT PURPOSES
!     (I.E. INCLUDES OVERSHOOT REGIONS). convective_flag=T IF
!     CONVECTIVE.
!  final_iteration_flag (LOK) : FLAG SET T WHEN FINAL ITERATION IS
!     BEING PERFORMED.
!  num_zones (M) : NUMBER OF MODEL POINTS.
!
!  OUTPUT VARIABLES:
!  composition IS UPDATED IN DIFCOM TO GIVE THE NEW RUN OF COMPOSITION
!  AFTER ANGULAR MOMENTUM TRANSPORT.
!
!  BEFORE THE LAST ITERATION(LOK=F),ONLY DIFFUSION OF H,HE4,HE3 CALCULATED
!  TO CALCULATE CHANGE IN MU GRADIENTS CAUSED BY DIFFUSION.
!
!  WHEN LOK IS T, ONLY DIFFUSION OF SPECIES HEAVIER THAN HE4 IS PERFORMED
!  (AS DIFFUSION OF LIGHTER SPECIES HAS ALREADY BEEN DONE).
!  THE EQUALLY SPACED GRID POINTS FROM THE LAST UNSTABLE REGION SOLVED
!  ARE USED INITIALLY, AND THEN THE PROGRAM CHECKS FOR ADDITIONAL UNSTABLE
!  REGIONS INTERIOR TO THAT.
! CALL MIXCOM FOR PREVIOUSLY IDENTIFIED UNSTABLE REGION IF THIS IS NOT
! THE CONVERGED CALL
subroutine ndifcom(timestep, diffusion_coeff, equally_spaced_diffusion_coeff, &
     equally_spaced_mass, log_density, log_luminosity, log_pressure, &
     log_radius, log_mass, enclosed_mass, shell_mass, log_total_mass, &
     velocity, zone_begin, zone_end, zone_max, zone_min, convective_flag, &
     final_iteration_flag, num_zones, composition, species_begin, species_end)
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: timestep
      double precision, intent(in) :: diffusion_coeff(json)
      double precision, intent(inout) :: equally_spaced_diffusion_coeff(json), &
           equally_spaced_mass(json)
      double precision, intent(in) :: log_density(json), log_luminosity(json), &
           log_pressure(json), log_radius(json), log_mass(json), &
           enclosed_mass(json), shell_mass(json), log_total_mass
      double precision, intent(in) :: velocity(json)
      integer, intent(inout) :: zone_begin, zone_end
      integer, intent(in) :: zone_max, zone_min
      logical, intent(in) :: convective_flag(json), final_iteration_flag
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: species_begin, species_end

      save

      logical :: unstable_region_active, unstable_zone_found, &
           two_zone_region
      integer :: search_start, zone_idx

      if (.not.final_iteration_flag) then
         call mixcom(timestep, equally_spaced_diffusion_coeff, &
              equally_spaced_mass, shell_mass, zone_begin, zone_end, &
              convective_flag, final_iteration_flag, num_zones, composition, &
              species_begin, species_end)
      else
! FIND UNSTABLE REGIONS IN ORDER, AND CALL MIXGRID TO SET UP THE
! EQUALLY SPACED GRID AND MIXCOM TO MIX THEM IN ORDER
!  EACH UNSTABLE REGION IS SOLVED SEPARATELY STARTING HERE.
!  unstable_zone_found IS SET T IF A NON-ZERO VELOCITY IS ENCOUNTERED.
!  zone_begin IS THE ZONE BELOW THE FIRST NON-ZERO V;zone_end IS THE
!  ZONE ABOVE THE LAST NON-ZERO V.
         unstable_region_active = .false.
         search_start = zone_min
   60    continue
         unstable_zone_found = .false.
         do zone_idx = search_start,zone_max
            if (velocity(zone_idx).gt.0.0d0) then
               unstable_zone_found = .true.
               if (.not.unstable_region_active) then
!  START OF UNSTABLE REGION
                  unstable_region_active = .true.
                  zone_begin = zone_idx - 1
               end if
            else if (unstable_region_active) then
!  END OF UNSTABLE REGION
               zone_end = zone_idx - 1
               unstable_region_active = .false.
               search_start = zone_idx + 1
               goto 80
            end if
         end do
!  IF THE LAST INTERFACE IS UNSTABLE (NON-ZERO V) ENSURE THAT zone_end
!  IS SET PROPERLY.
         if (unstable_region_active) zone_end = zone_max
         search_start = zone_max + 1
   80    continue
!  IF NO NON-ZERO V'S ENCOUNTERED, EXIT.
         if (.not.unstable_zone_found) goto 90
!  TRANSFORM TO EQUAL GRID SPACING IN R FOR THE REGION.
         call mixgrid(diffusion_coeff, log_density, log_luminosity, &
              log_pressure, log_radius, log_mass, enclosed_mass, shell_mass, &
              log_total_mass, zone_begin, zone_end, convective_flag, &
              num_zones, equally_spaced_diffusion_coeff, equally_spaced_mass, &
              two_zone_region)
!  two_zone_region=T IF TWO ZONES IN UNSTABLE REGION;
!  SKIP IF THIS OCCURS.
         if (two_zone_region) then
            if (search_start.le.zone_max) then
               goto 60
            else
               goto 90
            end if
         end if
!  PERFORM COMPOSITION DIFFUSION.
!  UNTIL THE FINAL ITERATION, ONLY COMPOSITION DIFFUSION OF SPECIES WHICH
!  AFFECT GRADIENTS IN MEAN MOLECULAR WEIGHT IS COMPUTED (H,HE3,HE4).
!  ON THE FINAL ITERATION, DIFFUSION OF ALL SPECIES IS PERFORMED.
         call mixcom(timestep, equally_spaced_diffusion_coeff, &
              equally_spaced_mass, shell_mass, zone_begin, zone_end, &
              convective_flag, final_iteration_flag, num_zones, composition, &
              species_begin, species_end)
!  RETURN FOR NEXT REGION IF APPLICABLE
         if (search_start.le.zone_max) goto 60
   90    continue
      end if
      return
end subroutine ndifcom

!----------------------------------------------------------------------
! rotmix
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original rotmix.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Burns each burn region (burn_step's burn_regions: shell by shell in
! the radiative regions, each convection zone as a unit, homogenized) by
! burn_network_region. Then (if diffuse_helium_active) performs gravitational settling of
! helium (Bahcall & Loeb method), and finally renormalizes the
! composition to guard against small negative/overflowing abundances.
! rate_by_zone holds the composition rate vectors of the model being
! burned (network_rates_of_model: mid_timestep_model passes the rotation
! sub-step's midpoint). (2026-10: they used to come from
! rot_scr%reaction_rate_by_zone, the rates of the last composition step,
! i.e. of the end-of-step structure.)
subroutine rotmix(timestep, composition, shell_mass, log_temperature, &
     rate_by_zone, num_zones, radiative_zone_bounds, convective_zone_bounds, &
     num_radiative_zones, num_convective_zones, log_total_mass, &
     log_density, log_mass, log_radius, log_pressure, convective_flag, &
     enclosed_mass, ierr)
      use rotation_scratch_lib

      use star_info_lib, only: star, json
      use io_units
      use run_log, only: solver_diagnostics
      use physical_constants
      use math_lib
      use burn_step, only: burn_network_region, burn_regions
      use element_diffusion_lib, only: diffuse_elements
      use net_lib, only: num_rate_rows
      implicit none

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: shell_mass(json), log_temperature(json)
      double precision, intent(in) :: rate_by_zone(num_rate_rows, json)
      integer, intent(in) :: num_zones
      integer, intent(in) :: radiative_zone_bounds(13,2), &
           convective_zone_bounds(12,2)
      integer, intent(in) :: num_radiative_zones, num_convective_zones
      double precision, intent(in) :: log_total_mass
      double precision, intent(in) :: log_density(json), log_mass(json), &
           log_radius(json), log_pressure(json)
      logical, intent(in) :: convective_flag(json)
! inout: the settling routines convert it to Bahcall & Loeb units in
! place and back, which can change its last bit
      double precision, intent(inout) :: enclosed_mass(json)

      integer :: num_species
      double precision :: dlnp_dr_settling(json)
      logical :: am_transport_convective_flag(json)
      double precision :: total_mass
      integer, allocatable :: regions(:,:)
      integer :: zone_idx, region_idx, species_idx, outer_boundary_zone, num_settling_substeps, &
           substep_idx
      double precision :: mass_fraction_above, settling_timescale, &
           max_settling_dt, settling_dt

! NSPEC IS THE NUMBER OF SPECIES BEING TRACKED.
      integer, intent(out) :: ierr

      ierr = 0

      if (star%job%use_extended_composition) then
         num_species = 15
      else
         num_species = 11
      end if
! THE NEW ABUNDANCES AT THE END OF THE (SUB-)STEP, timestep (s).
      regions = burn_regions(log_temperature, radiative_zone_bounds, &
           num_radiative_zones, convective_zone_bounds, num_convective_zones)
      do region_idx = 1, size(regions,2)
         call burn_network_region(log_temperature, regions(1,region_idx), &
              regions(2,region_idx), rate_by_zone, shell_mass, composition, &
              timestep, ierr)
         if (ierr /= 0) return
      end do
!
! MICROSCOPIC DIFFUSION OF HELIUM.
!
!***BC 5/92 ROTMIX MODIFIED TO INCLUDE CALL TO GRAVITATIONAL SETTLING
!      ROUTINE USING THE BAHCALL AND LOEB METHOD.
! FIRST DEFINE VARIABLES NEEDED FOR SETTLING -
! HQPR=VECTOR OF D LN P/DR.
! STOT=TOTAL STELLAR MASS(UNLOGGED).
      if (star%job%diffuse_helium_active) then
      settling: do
         if (composition(1,1).lt.star%ctrl%hydrogen_diffusion_floor) then
            star%job%diffuse_helium_active=.false.
            exit settling
         end if
! MHP 6/90 CHANGE ADDED : THE TIMESTEP FOR SETTLING IS RESTRICTED TO
!   A FRACTION OF THE TIMESCALE FOR SETTLING AT THE OUTER BOUNDARY.
!   THE OUTER BOUNDARY IS EITHER THE SURFACE CONVECTION ZONE OR THE
!   FIRST POINT WHERE THE HELIUM ABUNDANCE RISES ABOVE A USER-SPECIFIED
!   MINIMUM VALUE (YMIN).
!
!   LOCATE OUTER BOUNDARY.
         if (.not.convective_flag(num_zones)) then
            write(*,911)
            write(run_log_unit,911)
  911       format(1x,'NO SURFACE CZ - DIFFUSION NOT MEANINGFUL'/ &
                 'STOPPED IN SUBROUTINE MIX')
            ! 2026 (docs/ROADMAP.md stage 3): stop converted to ierr; the driver-side
            ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
            ! preserve the historical stop on a nonzero return.
            ierr = 1
            return
         end if
! 7/92 MHP STATEMENT ADDED FOR EXIT IF OVERSHOOT CAUSES A FULLY CONVECTIVE
! CASE.  AVOIDS OCCASIONAL PRE-MS CRASH IN FIRST RADIATIVE MODEL.
         if (convective_zone_bounds(num_convective_zones,1).gt.2 .or. convective_zone_bounds(num_convective_zones,2).ne.num_zones) then
         outer_boundary_zone = radiative_zone_bounds(num_radiative_zones,2)
         do zone_idx = outer_boundary_zone,1,-1
            if (composition(2,zone_idx).gt.star%ctrl%helium_diffusion_min) exit
         end do
         if (zone_idx < (1)) then
!   Y<YMIN FOR THE WHOLE STAR IF THE CODE GETS HERE.
         exit settling
         end if
         total_mass=exp(ln10*log_total_mass)
         do zone_idx = 1,num_zones
            dlnp_dr_settling(zone_idx)=-exp(ln10*(log_density(zone_idx)+ &
                 cgl+log_mass(zone_idx)-2.0d0*log_radius(zone_idx)- &
                 log_pressure(zone_idx)))
         end do
! ***BC 6/92 only check for timestep cutting if JMAX is large
!
         if (outer_boundary_zone.ge.2) then
!  FM IS THE MASS FRACTION ABOVE THE OUTER POINT.
            mass_fraction_above = (total_mass- &
                 enclosed_mass(outer_boundary_zone))/total_mass
!  TSCALE IS THE TIMESCALE FOR SETTLING OF HELIUM AT THE OUTER
!  BOUNDARY (MICHAUD ET AL 1984, APJ V.282,P.206)
            settling_timescale = 4.348d21*seconds_per_year* &
                 mass_fraction_above/ &
                 exp(ln10*1.5d0*log_temperature(outer_boundary_zone))
!  RESTRICT TIMESTEP TO THE MINIMUM OF THE MODEL TIMESTEP AND
!  A USER SPECIFIED FRACTION (DT_GS) OF THE SETTLING TIMESCALE.
! (2026: this settling diagnostic wrote to never-opened unit 69,
! leaving stray fort.69 files -- now solver forensics in the run
! log, like its peers)
            if (solver_diagnostics()) write(run_log_unit,*) 'JMAX=', &
                 outer_boundary_zone,' FM= ', &
                 mass_fraction_above,' TSCALE=',settling_timescale
            max_settling_dt = star%ctrl%settling_timestep_fraction*settling_timescale
            num_settling_substeps = int(timestep/max_settling_dt)
            if (mod(max_settling_dt,timestep).ne.0.0d0.or. &
                 num_settling_substeps.eq.0) &
                 num_settling_substeps=num_settling_substeps+1
            if (solver_diagnostics()) write(run_log_unit,*) &
                 'DTMAX=',max_settling_dt,' DELTS=',timestep, &
                 ' NSTEP=',num_settling_substeps
            settling_dt = timestep/dfloat(num_settling_substeps)
         else
            num_settling_substeps=1
            settling_dt = timestep
         end if
         do region_idx = 1,num_radiative_zones
            do zone_idx = radiative_zone_bounds(region_idx,1), &
                 radiative_zone_bounds(region_idx,2)
               am_transport_convective_flag(zone_idx) = .false.
            end do
         end do
         do region_idx = 1,num_convective_zones
            do zone_idx = convective_zone_bounds(region_idx,1), &
                 convective_zone_bounds(region_idx,2)
               am_transport_convective_flag(zone_idx) = .true.
            end do
         end do
         do substep_idx = 1,num_settling_substeps
! gravitational settling, with the mid-step temperature gradient
            call diffuse_elements(settling_dt,composition,dlnp_dr_settling, &
                 log_radius,log_density,enclosed_mass,log_temperature, &
                 mix_scr%delm,am_transport_convective_flag,num_zones, &
                 total_mass, ierr)
            if (ierr /= 0) return
         end do
         end if
      exit settling
      end do settling
      end if
! END OF GRAVITATIONAL SETTLING
!
! RENORMALIZE COMPOSITION TO GUARD AGAINST ANOMALIES (I.E. SMALL NEGATIVE
! ABUNDANCES...).
      do zone_idx = 1,num_zones
         do species_idx = 1,num_species
            composition(species_idx,zone_idx) = &
                 max(composition(species_idx,zone_idx),0.0d0)
            composition(species_idx,zone_idx) = &
                 min(composition(species_idx,zone_idx),1.0d0)
         end do
         composition(3,zone_idx) = min(composition(3,zone_idx), &
              1.0d0-composition(1,zone_idx)-composition(4,zone_idx))
         composition(2,zone_idx) = 1.0d0 - composition(1,zone_idx) - &
              composition(3,zone_idx) - composition(4,zone_idx)
      end do
      return
end subroutine rotmix

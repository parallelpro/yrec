!----------------------------------------------------------------------
! mdot
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mdot.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MASS LOSS ROUTINE
!
! PRESERVED CALL-SITE BUG (not fixed, per project policy): mdot.f's
! only caller, massloss.f90, calls this routine with 24 actual
! arguments (leading with log_luminosity_lsun/BL) while this
! subroutine declares only 23 dummy arguments (starting with
! timestep/DELTS, no luminosity argument at all). This is a pre-
! existing argument-count mismatch in the original F77 source
! (mdot.f/massloss.f) -- every actual argument after the first is
! therefore received one position off from what the caller intended.
! Reproduced exactly; NOT corrected here.
subroutine mdot(timestep, composition, log_density, specific_angular_momentum, &
     log_pressure, log_radius, log_mass, zone_mass_grams, shell_mass, &
     log_total_mass, log_temperature, envelope_boundary_zone, &
     new_surface_bc_needed, num_zones, omega, mean_molecular_weight_local, &
     total_radius_cm, total_mass_msun, mass_loss_rate_msun_yr, &
     accretion_specific_energy, mean_thermal_energy, &
     cz_total_mass_below_fitting, old_log_envelope_mass_fraction)
      use star_info_lib, only: star, json
      use phys_const_lib
      implicit none

      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: log_density(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      double precision, intent(inout) :: log_pressure(json)
      double precision, intent(inout) :: log_radius(json)
      double precision, intent(inout) :: log_mass(json)
      double precision, intent(inout) :: zone_mass_grams(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(inout) :: log_total_mass
      double precision, intent(inout) :: log_temperature(json)
      integer, intent(in) :: envelope_boundary_zone
      logical, intent(out) :: new_surface_bc_needed
      integer, intent(in) :: num_zones
      double precision, intent(in) :: omega(json)
      double precision, intent(out) :: mean_molecular_weight_local
      double precision, intent(in) :: total_radius_cm
      double precision, intent(inout) :: total_mass_msun
      double precision, intent(inout) :: mass_loss_rate_msun_yr
      double precision, intent(in) :: accretion_specific_energy, &
           mean_thermal_energy
      double precision, intent(in) :: cz_total_mass_below_fitting
      double precision, intent(out) :: old_log_envelope_mass_fraction
! --- locals ---
      double precision :: omega_ratio_sq, omega_max_ratio_sq
      double precision :: mass_loss_rate_cgs
      double precision :: cz_mass_below_fitting, cz_mass_grams
      double precision :: timestep_limit_total_mass, timestep_limit_cz_mass, &
           timestep_limit
      double precision :: timestep_limit_envelope
      double precision :: disk_age_test
      logical :: disk_exhausted_flag
      double precision :: delta_mass_cgs
      double precision :: new_thermal_energy, log_thermal_energy_ratio
      double precision :: local_temperature, local_pressure, local_density, &
           local_beta
      double precision :: delta_log_specific_entropy, delta_ln_mass, &
           delta_log_radius
      integer :: zone_idx, species_idx
      double precision :: boundary_radius_cm, radius_scale_factor, &
           radius_before_cm, radius_after_cm
      double precision :: surface_moment_of_inertia_per_mass, &
           delta_angular_momentum, surface_omega_local
      double precision :: cz_total_am, cz_new_am, am_ratio
      double precision :: mass_scale_factor
      double precision :: delta_mass_msun
      double precision :: total_mass_grams_old, total_mass_grams_new
      double precision :: mixed_abundance

      old_log_envelope_mass_fraction = log_mass(num_zones) - log_total_mass
! MHP 8/10- CHECK FOR SCALED SOLAR WIND MASS LOSS
      if(star%rot%use_rotation_scaled_solar_wind .and. star%job%rotation_active) then
         omega_ratio_sq = (omega(num_zones)/star%rot%wind_reference_omega)**2
         omega_max_ratio_sq = (star%rot%wind_max_omega/star%rot%wind_reference_omega)**2
         mass_loss_rate_msun_yr = star%rot%solar_wind_mass_loss_rate_msun_yr* &
              min(omega_ratio_sq,omega_max_ratio_sq)
         write(*,*)omega(num_zones),mass_loss_rate_msun_yr
      endif
!      LMDOT = .TRUE.
!      DMDT = 4.0D-10
! MAXIMUM FRACTION OF CZ MASS REMOVED PER TIMESTEP
!      FRAC = 0.1D0
!      IF(.NOT.LMDOT)GOTO 9999
! CONVERT FROM SOLAR MASSES/YEAR TO GM/SEC
      mass_loss_rate_cgs = abs(mass_loss_rate_msun_yr)*star%solar_mass_cgs/ &
           seconds_per_year
! THE SUM OF THE MASSES OF ALL SHELLS (E.G. TO JENV - 1)
! SHOULD BE USED RATHER THAN THE SUM OF THE MASSES DOWN
! TO THE MIDPOINT OF THE BOTTOM SHELL
!      DMCZT = 10.0D0**HSTOT - HS1(JENV)
      cz_mass_below_fitting = cz_total_mass_below_fitting
! COMPUTE MASS OF SURFACE CONVECTION ZONE BELOW FITTING POINT
      if(envelope_boundary_zone.eq.num_zones) stop 911
      if(envelope_boundary_zone.gt.1) then
         cz_mass_grams = zone_mass_grams(num_zones)- &
              zone_mass_grams(envelope_boundary_zone)
      else
         cz_mass_grams = zone_mass_grams(num_zones)
      endif
! COMPUTE MAXIMUM TIMESTEP BASED ON NOT REMOVING TOO MUCH MASS FROM THE
! SURFACE CZ (FCZDMDT) OR AS A FRACTION OF THE TOTAL MASS (FTOTDMDT)
      timestep_limit_total_mass = star%ctrl%ftotdmdt*total_mass_msun*star%solar_mass_cgs/ &
           mass_loss_rate_cgs
      timestep_limit_cz_mass = star%ctrl%fczdmdt*cz_mass_grams/mass_loss_rate_cgs
      timestep_limit = min(timestep_limit_total_mass,timestep_limit_cz_mass)
! RESTRICT TIMESTEP TO ADD NO MORE THAN 1/2 OF THE CURRENT MASS
! BEYOND THE FITTING POINT TO THE STAR.
      if(mass_loss_rate_cgs.gt.0.0d0)then
         timestep_limit_envelope = 0.5d0*(10.0d0**log_total_mass- &
              10.0d0**log_mass(num_zones))/mass_loss_rate_cgs
         timestep_limit = min(timestep_limit,timestep_limit_envelope)
      endif
      if(timestep_limit.lt.timestep)then
         write(*,10)timestep,timestep_limit
 10      format(' TIMESTEP REDUCED FOR MASS LOSS - OLD,NEW',1P2E12.3)
         timestep = timestep_limit
      endif
! mhp 8/10 turn mass loss off when disk exhausted only when dm /dt > 0, e.g. accretion
      if(star%job%disk_locking_active .and. mass_loss_rate_msun_yr.gt.0.0d0)then
         disk_age_test = star%disk_lifetime + 1.0d-9*timestep/seconds_per_year
         if(disk_age_test.gt.star%job%disk_temperature)then
            timestep = (star%job%disk_temperature-star%disk_lifetime)*1.0d9*seconds_per_year
            disk_exhausted_flag = .true.
         else
            disk_exhausted_flag = .false.
         endif
      else
         disk_exhausted_flag = .false.
      endif
! FINAL AMOUNT OF MASS LOSS INFERRED IN CGS UNITS.
      delta_mass_cgs = mass_loss_rate_msun_yr*star%solar_mass_cgs*timestep/ &
           seconds_per_year
! COMPUTE THE MEAN THERMAL ENERGY
! CONTENT OF THE CONVECTION ZONE FOR MODELS WITH ACCRETION.
      if(delta_mass_cgs.gt.0.0d0)then
         new_thermal_energy = (accretion_specific_energy*delta_mass_cgs+ &
              mean_thermal_energy*cz_total_mass_below_fitting)/ &
              (delta_mass_cgs+cz_total_mass_below_fitting)
         log_thermal_energy_ratio = log10(new_thermal_energy/mean_thermal_energy)
!         WRITE(*,912) JENV,EACC,ETHAV,ETHNEW,FAC,SUMDM,DELM
! 912     FORMAT(I5,1P6E12.3)
! OVERALL SCALE FACTOR IN R
         local_temperature = 10.0d0**log_temperature(envelope_boundary_zone)
         local_pressure = 10.0d0**log_pressure(envelope_boundary_zone)
         local_density = 10.0d0**log_density(envelope_boundary_zone)
         local_beta = 1.0d0-(radiation_constant_over_3*local_temperature**4/ &
              local_pressure)
         mean_molecular_weight_local = local_pressure*local_beta/ &
              (local_density*local_temperature)
!         DLOGEN = (DELM/DMCZ)*(SACC-SCEN)/RMU
         delta_log_specific_entropy = (delta_mass_cgs/cz_mass_grams)* &
              star%rot%envelope_specific_entropy/mean_molecular_weight_local
!         DLNM = LOG(DMCZ+DELM)-LOG(DMCZ)
         delta_ln_mass = 0.0d0
         delta_log_radius = (cc23*delta_log_specific_entropy- &
              cc13*delta_ln_mass)/ln10
         star%rot%delta_log_temperature = delta_ln_mass/ln10 - delta_log_radius
         star%rot%delta_log_pressure = 2.0d0*delta_ln_mass/ln10 - 4.0d0*delta_log_radius
!         WRITE(*,*)DLOGR,DLOGP,DLOGT
         if(envelope_boundary_zone.eq.1)then
            do zone_idx = envelope_boundary_zone,num_zones
               log_radius(zone_idx) = log_radius(zone_idx)+delta_log_radius
               log_pressure(zone_idx) = log_pressure(zone_idx)+star%rot%delta_log_pressure
               log_temperature(zone_idx) = &
                    log_temperature(zone_idx)+star%rot%delta_log_temperature
            end do
         else
            log_pressure(envelope_boundary_zone) = &
                 log_pressure(envelope_boundary_zone)+star%rot%delta_log_pressure
            log_temperature(envelope_boundary_zone) = &
                 log_temperature(envelope_boundary_zone)+star%rot%delta_log_temperature
            boundary_radius_cm = 10.0d0**log_radius(envelope_boundary_zone)
            radius_scale_factor = 10.0d0**delta_log_radius
            do zone_idx = envelope_boundary_zone+1,num_zones
               radius_before_cm = 10.0d0**log_radius(zone_idx)
               radius_after_cm = boundary_radius_cm+radius_scale_factor* &
                    (radius_before_cm-boundary_radius_cm)
               log_radius(zone_idx) = log10(radius_after_cm)
            end do
         endif
      endif
! FOR MODELS WITH ROTATION AND MASS LOSS REMOVE ANGULAR MOMENTUM.
! IN THIS CASE WE ASSUME A THERMAL WIND WHERE THE SURFACE MATTER
! CARRIES AWAY ONLY ITS LOCAL ANGULAR MOMENTUM PER UNIT MASS.
      if(star%job%rotation_active .and. delta_mass_cgs.lt.0.0d0)then
! MOMENT OF INERTIA PER UNIT MASS AT THE SURFACE.
         surface_moment_of_inertia_per_mass = 2.0d0/3.0d0*total_radius_cm**2
         if(star%ctrl%walpcz.ge.0.0d0)then
! SOLID BODY CZ ROTATION
            delta_angular_momentum = omega(num_zones)* &
                 surface_moment_of_inertia_per_mass*delta_mass_cgs
         else if(star%ctrl%walpcz.le.-2.0d0)then
! CONSTANT J/M
            delta_angular_momentum = specific_angular_momentum(num_zones)* &
                 delta_mass_cgs
         else
            surface_omega_local = omega(num_zones)*10.0d0** &
                 (log_radius(num_zones)*star%ctrl%walpcz)/total_radius_cm**star%ctrl%walpcz
            delta_angular_momentum = surface_omega_local* &
                 surface_moment_of_inertia_per_mass*delta_mass_cgs
         endif
!  FIND TOTAL STARTING ANGULAR MOMENTUM OF C.Z.
         cz_total_am = specific_angular_momentum(envelope_boundary_zone)* &
              shell_mass(envelope_boundary_zone)
         do zone_idx = envelope_boundary_zone+1,num_zones
            cz_total_am = cz_total_am + &
                 specific_angular_momentum(zone_idx)*shell_mass(zone_idx)
         end do
         cz_new_am = (cz_total_am+delta_angular_momentum)
         am_ratio = cz_new_am/cz_total_am
         write(*,22)delta_mass_cgs,delta_angular_momentum,omega(num_zones), &
              cz_new_am,cz_total_am,cz_mass_below_fitting
 22      format(1P6E13.4)
         do zone_idx = envelope_boundary_zone,num_zones
            specific_angular_momentum(zone_idx) = &
                 specific_angular_momentum(zone_idx)*am_ratio
         end do
      endif
! REMOVE OR ADD MASS IN THE CONVECTION ZONE.
! HS1 IS THE LOCATION OF THE ZONE CENTERS IN GM;
! HS IS THE BASE-10 LOG OF HS1
      if(delta_mass_cgs.lt.0.0d0)then
         mass_scale_factor = (cz_mass_grams+delta_mass_cgs)/cz_mass_grams
!      WRITE(*,*)JENV,M,DMCZ,DELM,FX,HS1(M),HS(M)
         do zone_idx = envelope_boundary_zone+1,num_zones
            zone_mass_grams(zone_idx) = zone_mass_grams(envelope_boundary_zone)+ &
                 mass_scale_factor*(zone_mass_grams(zone_idx)- &
                 zone_mass_grams(envelope_boundary_zone))
            log_mass(zone_idx) = log10(zone_mass_grams(zone_idx))
         end do
      endif
!      WRITE(*,*)HS1(M),HS(M)
! CORRECT TOTAL MASS IN SOLAR UNITS (SMASS) AND
! LOG OF TOTAL MASS IN GRAMS (HSTOT,STOTAL)
      total_mass_msun = total_mass_msun + delta_mass_cgs/star%solar_mass_cgs
      star%rot%updated_mass_msun = total_mass_msun
      delta_mass_msun = delta_mass_cgs/star%solar_mass_cgs
      write(*,20)total_mass_msun,delta_mass_msun
 20   format('MASS LOSS APPLIED - NEW M,DEL M',1P2E19.10)
      total_mass_grams_old = 10.0d0**log_total_mass
      total_mass_grams_new = total_mass_grams_old + delta_mass_cgs
      log_total_mass = log10(total_mass_grams_new)
      star%stotal = log_total_mass
! CORRECT MASS CONTENTS OF INDIVIDUAL SHELLS (HS2, IN GM)
      do zone_idx = envelope_boundary_zone,num_zones-1
       shell_mass(zone_idx) = 0.5d0*(zone_mass_grams(zone_idx+1)- &
            zone_mass_grams(zone_idx-1))
      end do
      shell_mass(num_zones) = total_mass_grams_new-0.5d0* &
           (zone_mass_grams(num_zones)+zone_mass_grams(num_zones-1))
! 07/02RESET SENV
      star%senv = log_mass(num_zones) - log_total_mass
! RECOMPUTE SURFACE BOUNDRY CONDITION
      new_surface_bc_needed = .true.
! REMIX THE SURFACE CONVECTION ZONE IF MDOT IS POSITIVE.
!      WRITE(*,911)DELM,DMCZT,COMPACC(12),HCOMP(12,M),
!     *            COMPACC(4),HCOMP(4,M)
! 911  FORMAT(1P6E12.3)
      if(delta_mass_cgs.gt.0.0d0)then
         do species_idx = 1,11
            mixed_abundance = (composition(species_idx,envelope_boundary_zone)* &
                 cz_mass_below_fitting+star%ctrl%accreted_composition(species_idx)* &
                 delta_mass_cgs)/(delta_mass_cgs+cz_mass_below_fitting)
            do zone_idx = envelope_boundary_zone,num_zones
               composition(species_idx,zone_idx) = mixed_abundance
            end do
         end do
         star%accreted_mass_fraction = delta_mass_cgs
      endif
!  9999 CONTINUE
      if(disk_exhausted_flag) star%job%use_mass_accretion = .false.
!      WRITE(*,*)M,HS(M),HSTOT
      return
end subroutine mdot

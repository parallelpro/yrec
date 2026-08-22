!----------------------------------------------------------------------
! massloss
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original massloss.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 10/01
! MASS LOSS DRIVER ROUTINE
!
! This is experimental code and valid for Allard atmospheres only.
!   llp  06/15/2009
!
! Computes the mass-loss (or accretion) rate for the current model
! (constant default rate, or a Reimers-law rate if use_reimers_law is
! set), the mean thermal-energy content of the surface convective
! envelope and of the accreted matter (needed by mdot.f90's radius/
! pressure/temperature adjustment), and the specific entropy of
! accreted matter (via an Allard-atmosphere surface-pressure/
! temperature solve, ALSURFP+EQSTAT), then calls mdot.f90 to actually
! apply the mass change.
!
! PRESERVED CALL-SITE BUG (not fixed, per project policy): the call to
! mdot below passes 24 actual arguments (leading with
! log_luminosity_lsun) to a subroutine, mdot.f90, that declares only
! 23 dummy arguments (no luminosity argument) -- see mdot.f90's header
! for the full description of this pre-existing mismatch, reproduced
! exactly here.
!
! FORMER CROSS-FILE COMMON LAYOUT MISMATCH (now resolved by the
! rotdiff_lib module conversion): this file's own former common/
! masschg2/ declared its 4th/5th members as (DLOGT, DLOGP) -- i.e.
! DLOGT before DLOGP -- while mdot.f90's declaration of the very same
! block had them in the opposite order (DLOGP, DLOGT). Since COMMON
! storage is positional, this meant the two files disagreed about
! which physical quantity occupied which storage slot -- harmless only
! because neither variable was actually read or set in this file's
! body (unused layout placeholders here). Name-based module access
! makes the position irrelevant, so the historical swap is moot now;
! documented here for the record rather than silently dropped.
subroutine massloss(log_luminosity_lsun, age_gyr, timestep, composition, &
     log_density, specific_angular_momentum, log_pressure, log_radius, &
     log_mass, zone_mass_grams, shell_mass, log_total_mass, log_temperature, &
     envelope_boundary_zone, new_surface_bc_needed, num_zones, omega, &
     total_mass_msun, log_teff, old_log_envelope_mass_fraction, &
     new_atmosphere_fit_needed)
      use atm_lib
      use atm_table_lib
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use const_lib
      use eos_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: log_luminosity_lsun, age_gyr
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
      double precision, intent(inout) :: total_mass_msun
      double precision, intent(in) :: log_teff
      double precision, intent(out) :: old_log_envelope_mass_fraction
      logical, intent(out) :: new_atmosphere_fit_needed











! MHP 5/02 EFFICIENCY FACTOR FOR THE THERMAL ENERGY CONTENT
! OF ACCRETED MATTER.
      double precision :: accretion_efficiency
      data accretion_efficiency/1.0d0/
! --- locals ---
      double precision :: mass_loss_rate_msun_yr
      logical :: apply_mass_change
      double precision :: log10_radius, total_radius_cm
      double precision :: total_mass_grams, age_seconds
      double precision :: surface_gravity_cgs
      double precision :: cz_total_mass_below_fitting
      double precision :: local_temperature, local_pressure, local_density, &
           local_beta, mean_molecular_weight_local
      double precision :: accretion_specific_energy, &
           accretion_specific_energy0
      double precision :: sum_thermal_energy, thermal_energy_accreted_bar
      double precision :: local_radius_cm, delta_accretion_energy, &
           thermal_energy_accreted, local_entropy, thermal_energy_per_gram
      double precision :: mean_thermal_energy
      logical :: print_flag
      double precision :: log10_gravity
      logical :: allard_surface_failed
      double precision :: temperature_local, pressure_local, log10_temperature_local, &
           log10_pressure_local, density_local, log10_density_local
      double precision :: hydrogen_fraction_local, metal_fraction_local
      logical :: eos_deriv_flag, eos_atmosphere_flag
      integer :: saha_flag
      double precision :: beta_local, beta_ion, beta14, ion_fraction(3), &
           mean_molecular_weight_eos, amu_eos, emu_eos, eta_eos, qdt_eos, &
           qdp_eos, qcp_eos, dela_eos, qdtt_eos, qdtp_eos, qat_eos, qap_eos, &
           qcpt_eos, qcpp_eos
      double precision :: mass_loss_rate_cgs, pressure_from_wind, &
           temperature_from_wind, accretion_specific_entropy2
      integer :: zone_idx

! INITIALIZE MASS LOSS AT DEFAULT RATE
      mass_loss_rate_msun_yr = mass_accretion_rate
      if(use_mass_accretion)then
         apply_mass_change = .true.
      else
         apply_mass_change = .false.
         new_atmosphere_fit_needed = .false.
         goto 9999
      endif
!      IF(.NOT.LDOMDOT)RETURN
! TEFFL IS THE BASE 10 LOG OF THE EFFECTIVE TEMPERATURE
! COMPUTE GLOBAL QUANTITIES (RADIUS,MASS,AGE) IN CGS UNITS.
! RADIUS
      log10_radius = 0.5d0*(log_luminosity_lsun+log10_solar_luminosity-c4pil- &
           csigl-4.0d0*log_teff)
      total_radius_cm = 10.0d0**log10_radius
! MASS
      total_mass_grams = total_mass_msun*solar_mass_cgs
! AGE
      age_seconds = age_gyr*1.0d9*seconds_per_year
! USE A REIMERS FORMULA TO COMPUTE MDOT IF DESIRED; OVERWRITES
! CONSTANT MDOT.  IN THIS EXPRESSION MDOT=K*L/G/R.
      if(apply_mass_change .and. lreimer)then
         surface_gravity_cgs = 10.0d0**(cgl)*total_mass_grams/total_radius_cm**2
         mass_loss_rate_msun_yr = creim*10.0d0**(log_luminosity_lsun+ &
              log10_solar_luminosity)/surface_gravity_cgs/total_radius_cm
      endif
! 02/12 MHP TAUCZ NOW COMPUTED PRIOR TO CALL IN MIXCZ
! CONVECTIVE OVERTURN TIMESCALE
      if(envelope_boundary_zone.lt.num_zones)then
!         TAUCZ = 0.0D0
!         DO I = JENV+1,M
!            V = 0.5D0*(SVEL(I-1)+SVEL(I))
!            DR = 10.0D0**(HR(I))-10.0D0**(HR(I-1))
!            IF(V.GT.0.0D0)THEN
!               TAUCZ = TAUCZ + DR/V
!            ENDIF
!         END DO
         write(*,*)star%turnover%convective_turnover_timescale/seconds_per_year, &
              total_radius_cm/solar_radius_cgs
         star%light_burn%jcz = envelope_boundary_zone
      else
         star%turnover%convective_turnover_timescale = 0.0d0
         star%light_burn%jcz = num_zones
      endif
! MHP 8/10
! THE RUNNING TOTAL OF THE MASS OF THE
! CONVECTIVE ENVELOPE, NEEDED ELSEWHERE AS WELL
      cz_total_mass_below_fitting = 0.0d0
      do zone_idx = envelope_boundary_zone, num_zones
         cz_total_mass_below_fitting = cz_total_mass_below_fitting+ &
              shell_mass(zone_idx)
      end do
! COMPUTE THE MEAN MOLECULAR WEIGHT FOR THE CZ
      if(mass_loss_rate_msun_yr.gt.0.0d0)then
         new_atmosphere_fit_needed = .true.
         local_temperature = 10.0d0**log_temperature(envelope_boundary_zone)
         local_pressure = 10.0d0**log_pressure(envelope_boundary_zone)
         local_density = 10.0d0**log_density(envelope_boundary_zone)
         local_beta = 1.0d0-(radiation_constant_over_3*local_temperature**4/ &
              local_pressure)
         mean_molecular_weight_local = local_pressure*local_beta/ &
              (local_density*local_temperature)
! FROM THE VIRIAL THEOREM, A MAXIMUM OF
! 1/2 OF THE GRAVITATIONAL POTENTIAL ENERGY OF
! ACCRETED MATTER WILL BE CONVERTED INTO KINETIC
! ENERGY.  THE PARAMETER FACC IS USED TO INFER
! THE KINETIC ENERGY PER GRAM WHICH IS COMPARED
! WITH THE MEAN THERMAL ENERGY CONTENT OF THE
! CONVECTION ZONE.
         accretion_specific_energy = accretion_efficiency*total_mass_grams* &
              (10.0d0**cgl)/total_radius_cm
         accretion_specific_energy0 = accretion_specific_energy
! DETERMINE THE MASS-WEIGHTED THERMAL ENERGY (PER GM)
! IN EACH SHELL OF THE CONVECTIVE ENV.
         sum_thermal_energy = 0.0d0
!         SUMDM = 0.0D0
         thermal_energy_accreted_bar = 0.0d0
         star%rot%envelope_specific_entropy = 0.0d0
         do zone_idx = envelope_boundary_zone, num_zones
            local_temperature = 10.0d0**log_temperature(zone_idx)
            local_pressure = 10.0d0**log_pressure(zone_idx)
            local_density = 10.0d0**log_density(zone_idx)
            local_radius_cm = 10.0d0**log_radius(zone_idx)
            delta_accretion_energy = 0.5d0*accretion_specific_energy* &
                 (total_radius_cm/local_radius_cm-1.0d0)
            thermal_energy_accreted = accretion_specific_energy+ &
                 delta_accretion_energy
            thermal_energy_accreted_bar = thermal_energy_accreted_bar+ &
                 shell_mass(zone_idx)*thermal_energy_accreted
            local_beta = 1.0d0-(radiation_constant_over_3*local_temperature**4/ &
                 local_pressure)
            mean_molecular_weight_local = local_pressure*local_beta/ &
                 (local_density*local_temperature)
            local_entropy = mean_molecular_weight_local* &
                 (1.5d0*log(local_temperature)-log(local_density))
            star%rot%envelope_specific_entropy = star%rot%envelope_specific_entropy+ &
                 local_entropy*shell_mass(zone_idx)
! THE THERMAL ENERGY PER GM IN THE JTH SHELL IS
            thermal_energy_per_gram = local_pressure*local_beta/local_density
! THE RUNNING TOTAL OF THERMAL ENERGY THROUGHOUT THE
! CONVECTIVE ENVELOPE IS
            sum_thermal_energy = sum_thermal_energy+thermal_energy_per_gram* &
                 shell_mass(zone_idx)
! THE RUNNING TOTAL OF THE MASS OF THE
! CONVECTIVE ENVELOPE IS
!            SUMDM = SUMDM + HS2(J)
         end do
! THE AVERAGE THERMAL ENERGY OF THE UNPERTURBED CE (W/O ACCN) IS
         mean_thermal_energy = sum_thermal_energy/cz_total_mass_below_fitting
         accretion_specific_energy = thermal_energy_accreted_bar/ &
              cz_total_mass_below_fitting
         star%rot%envelope_specific_entropy = star%rot%envelope_specific_entropy/ &
              cz_total_mass_below_fitting
         print_flag = .false.
         log10_gravity = cgl+log_total_mass-2.0d0*log10_radius
! This is experimental code and valid for Allard atmospheres only.
!   llp  06/15/2009
         call atm_get_surface_pt(log_teff,log10_gravity,print_flag, &
              allard_surface_failed)
         temperature_local = 10.0d0**atm_table%atm_log10_temperature
         pressure_local = 10.0d0**atm_table%atm_log10_pressure
         log10_temperature_local = atm_table%atm_log10_temperature
         log10_pressure_local = atm_table%atm_log10_pressure
         hydrogen_fraction_local = composition(1,num_zones)
         metal_fraction_local = composition(3,num_zones)
         eos_deriv_flag = .false.
         saha_flag = 1
         eos_atmosphere_flag = .true.
         call eos_get(log10_temperature_local,temperature_local, &
              log10_pressure_local,pressure_local,log10_density_local, &
              density_local,hydrogen_fraction_local,metal_fraction_local, &
              beta_local,beta_ion,beta14,ion_fraction,mean_molecular_weight_eos, &
              amu_eos,emu_eos,eta_eos,qdt_eos,qdp_eos,qcp_eos,dela_eos, &
              qdtt_eos,qdtp_eos,qat_eos,qap_eos,qcpt_eos,qcpp_eos, &
              eos_deriv_flag,eos_atmosphere_flag,saha_flag)
         beta_local = 1.0d0-(radiation_constant_over_3*temperature_local**4/ &
              pressure_local)
         mean_molecular_weight_local = pressure_local*beta_local/ &
              (density_local*temperature_local)
         star%rot%accretion_specific_entropy = mean_molecular_weight_local* &
              (1.5d0*log(temperature_local)-log(density_local))
!         WRITE(*,911)TL,PL,SACC,SCEN
!  911     FORMAT(' TSUR,PSUR ',2F8.5,' SACC ',1PE12.3,' SCORE ',E12.3)
! ALTERNATE EXPRESSION FOR SURFACE PRESSURE AND LUMINOSITY, FROM STAHLER 1988
         if(accretion_efficiency.gt.0.0d0)then
            mass_loss_rate_cgs = mass_loss_rate_msun_yr*solar_mass_cgs/ &
                 seconds_per_year
            pressure_from_wind = mass_loss_rate_cgs/c4pi* &
                 sqrt(2.0d0*accretion_efficiency*10.0d0**log10_gravity/ &
                 total_radius_cm**3)
            temperature_from_wind = (mass_loss_rate_cgs* &
                 accretion_specific_energy0*0.75d0/c4pi/total_radius_cm**2/ &
                 csig)**0.25d0
            log10_pressure_local = log10(pressure_from_wind)
            log10_temperature_local = log10(temperature_from_wind)
            call eos_get(log10_temperature_local,temperature_from_wind, &
                 log10_pressure_local,pressure_from_wind,log10_density_local, &
                 density_local,hydrogen_fraction_local,metal_fraction_local, &
                 beta_local,beta_ion,beta14,ion_fraction,mean_molecular_weight_eos, &
                 amu_eos,emu_eos,eta_eos,qdt_eos,qdp_eos,qcp_eos,dela_eos, &
                 qdtt_eos,qdtp_eos,qat_eos,qap_eos,qcpt_eos,qcpp_eos, &
                 eos_deriv_flag,eos_atmosphere_flag,saha_flag)
            beta_local = 1.0d0-(radiation_constant_over_3*temperature_from_wind**4/ &
                 pressure_from_wind)
            mean_molecular_weight_local = pressure_from_wind*beta_local/ &
                 (density_local*temperature_from_wind)
            accretion_specific_entropy2 = mean_molecular_weight_local* &
                 (1.5d0*log(temperature_from_wind)-log(density_local))
!            WRITE(*,911)TL,PL,SACC2,SCEN
!            SACC = MAX(SACC,SACC2)
            star%rot%envelope_specific_entropy = 0.0d0
!            SCEN = SACC2 - SACC
         else
            star%rot%envelope_specific_entropy = 0.0d0
         endif
      endif
! CALL MASS LOSS OR ACCRETION ROUTINE
      if(apply_mass_change)call mdot(log_luminosity_lsun,timestep,composition, &
           log_density,specific_angular_momentum,log_pressure,log_radius,log_mass, &
           zone_mass_grams,shell_mass,log_total_mass,log_temperature, &
           envelope_boundary_zone,new_surface_bc_needed,num_zones,omega, &
           mean_molecular_weight_local,total_radius_cm,total_mass_msun, &
           mass_loss_rate_msun_yr,accretion_specific_energy,mean_thermal_energy, &
           cz_total_mass_below_fitting,old_log_envelope_mass_fraction)
 9999 continue
      return
end subroutine massloss

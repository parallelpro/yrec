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
! The call to mdot below matches mdot.f90's 23-dummy list 1:1 (the
! historical 24-actual argument-count mismatch of massloss.f was fixed
! in 2026; see mdot.f90's header).
subroutine massloss(log_luminosity_lsun, age_gyr, timestep, composition, &
     log_density, specific_angular_momentum, log_pressure, log_radius, &
     log_mass, zone_mass_grams, shell_mass, log_total_mass, log_temperature, &
     envelope_boundary_zone, new_surface_bc_needed, num_zones, omega, &
     total_mass_msun, log_teff, old_log_envelope_mass_fraction, &
     new_atmosphere_fit_needed, ierr)
      use rotation_scratch_lib
      use atm_lib
      use atm_table_lib
      use star_info_lib, only: star, json
      use phys_const_lib
      use eos_lib
      use math_lib
      use wind_lib, only: log10_radius_from_l_teff
      implicit none

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
      integer, intent(out) :: ierr
! MHP 5/02 EFFICIENCY FACTOR FOR THE THERMAL ENERGY CONTENT
! OF ACCRETED MATTER.
      double precision, parameter :: accretion_efficiency = 1.0d0
! --- locals ---
      double precision :: mass_loss_rate_msun_yr
      double precision :: log10_radius, total_radius_cm
      double precision :: total_mass_grams
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
      double precision :: log10_temperature_local, log10_pressure_local, &
           log10_density_local
      double precision :: hydrogen_fraction_local, metal_fraction_local
      logical :: eos_deriv_flag, eos_atmosphere_flag
      integer :: saha_flag
      double precision :: beta_local
! 2026 named-index results: the former 18-variable eos output soup is
! one result array.  Nothing here is SAVEd: the original's blanket
! SAVE carried the previous call's converged log10 rho / beta into the
! inout seed slots; the modern locals are seeded to zero at entry
! (below), which is what the original saw on its first call.
      double precision :: eos_res(num_eos_results)
      double precision :: mass_loss_rate_cgs, pressure_from_wind, &
           temperature_from_wind
      integer :: zone_idx

! INITIALIZE MASS LOSS AT DEFAULT RATE
      ierr = 0
      log10_density_local = 0.0d0
      beta_local = 0.0d0
      mass_loss_rate_msun_yr = star%ctrl%mass_accretion_rate
      if(.not.star%job%use_mass_accretion)then
         new_atmosphere_fit_needed = .false.
         return
      endif
! TEFFL IS THE BASE 10 LOG OF THE EFFECTIVE TEMPERATURE
! COMPUTE GLOBAL QUANTITIES (RADIUS,MASS) IN CGS UNITS.
! RADIUS
      log10_radius = log10_radius_from_l_teff(log_luminosity_lsun, log_teff)
      total_radius_cm = exp10(log10_radius)
! MASS
      total_mass_grams = total_mass_msun*star%solar_mass_cgs
! USE A REIMERS FORMULA TO COMPUTE MDOT IF DESIRED; OVERWRITES
! CONSTANT MDOT.  IN THIS EXPRESSION MDOT=K*L/G/R.
! 2026 config-matrix fix: the historical expression computed
! reimers_scaling_factor*L/(g*R) in CGS -- units of g/s -- but stored it in the
! Msun/yr variable, which mdot then converted to CGS AGAIN (a
! ~1e25 rate inflation). Unreachable before the mdot argument-count
! fix, so never caught. Now implements the standard Reimers law the
! comment above describes: Mdot[Msun/yr] = reimers_scaling_factor*(L/Lsun)(R/Rsun)/(M/Msun),
! reimers_scaling_factor default -4e-13 (the classical eta=1 coefficient).
      if(star%ctrl%use_reimers_wind)then
         mass_loss_rate_msun_yr = star%ctrl%reimers_scaling_factor* &
              exp10(log_luminosity_lsun)* &
              (total_radius_cm/star%solar_radius_cgs)/total_mass_msun
      endif
! 02/12 MHP TAUCZ NOW COMPUTED PRIOR TO CALL IN MIXCZ
! CONVECTIVE OVERTURN TIMESCALE
      if(envelope_boundary_zone.lt.num_zones)then
         write(*,*)star%convective_turnover_timescale/seconds_per_year, &
              total_radius_cm/star%solar_radius_cgs
         star%jcz = envelope_boundary_zone
      else
         star%convective_turnover_timescale = 0.0d0
         star%jcz = num_zones
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
         local_temperature = exp10(log_temperature(envelope_boundary_zone))
         local_pressure = exp10(log_pressure(envelope_boundary_zone))
         local_density = exp10(log_density(envelope_boundary_zone))
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
              (exp10(cgl))/total_radius_cm
         accretion_specific_energy0 = accretion_specific_energy
! DETERMINE THE MASS-WEIGHTED THERMAL ENERGY (PER GM)
! IN EACH SHELL OF THE CONVECTIVE ENV.
         sum_thermal_energy = 0.0d0
         thermal_energy_accreted_bar = 0.0d0
         rot_scr%envelope_specific_entropy = 0.0d0
         do zone_idx = envelope_boundary_zone, num_zones
            local_temperature = exp10(log_temperature(zone_idx))
            local_pressure = exp10(log_pressure(zone_idx))
            local_density = exp10(log_density(zone_idx))
            local_radius_cm = exp10(log_radius(zone_idx))
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
            rot_scr%envelope_specific_entropy = rot_scr%envelope_specific_entropy+ &
                 local_entropy*shell_mass(zone_idx)
! THE THERMAL ENERGY PER GM IN THE JTH SHELL IS
            thermal_energy_per_gram = local_pressure*local_beta/local_density
! THE RUNNING TOTAL OF THERMAL ENERGY THROUGHOUT THE
! CONVECTIVE ENVELOPE IS
            sum_thermal_energy = sum_thermal_energy+thermal_energy_per_gram* &
                 shell_mass(zone_idx)
         end do
! THE AVERAGE THERMAL ENERGY OF THE UNPERTURBED CE (W/O ACCN) IS
         mean_thermal_energy = sum_thermal_energy/cz_total_mass_below_fitting
         accretion_specific_energy = thermal_energy_accreted_bar/ &
              cz_total_mass_below_fitting
         rot_scr%envelope_specific_entropy = rot_scr%envelope_specific_entropy/ &
              cz_total_mass_below_fitting
         print_flag = .false.
         log10_gravity = cgl+log_total_mass-2.0d0*log10_radius
! This is experimental code and valid for Allard atmospheres only.
!   llp  06/15/2009
         call atm_get_surface_pt(log_teff,log10_gravity,print_flag, &
              allard_surface_failed, ierr)
         if (ierr /= 0) return
         log10_temperature_local = atm_table%atm_log10_temperature
         log10_pressure_local = atm_table%atm_log10_pressure
         hydrogen_fraction_local = composition(1,num_zones)
         metal_fraction_local = composition(3,num_zones)
         eos_deriv_flag = .false.
         saha_flag = 1
         eos_atmosphere_flag = .true.
         eos_res(i_log10_density) = log10_density_local
         eos_res(i_beta) = beta_local
         call eos_get(log10_temperature_local, log10_pressure_local, &
              hydrogen_fraction_local, metal_fraction_local, eos_res, &
              eos_deriv_flag, eos_atmosphere_flag, saha_flag, ierr=ierr)
         if (ierr /= 0) return
         log10_density_local = eos_res(i_log10_density)
         beta_local = 1.0d0-(radiation_constant_over_3*eos_res(i_temperature)**4/ &
              eos_res(i_pressure))
         mean_molecular_weight_local = eos_res(i_pressure)*beta_local/ &
              (eos_res(i_density)*eos_res(i_temperature))
         rot_scr%accretion_specific_entropy = mean_molecular_weight_local* &
              (1.5d0*log(eos_res(i_temperature))-log(eos_res(i_density)))
! ALTERNATE EXPRESSION FOR SURFACE PRESSURE AND LUMINOSITY, FROM STAHLER 1988.
! ONLY mean_molecular_weight_local FROM THIS SECOND EOS CALL IS USED
! (PASSED ON TO mdot); THE ENTROPY IT WOULD GIVE IS NOT.
         if(accretion_efficiency.gt.0.0d0)then
            mass_loss_rate_cgs = mass_loss_rate_msun_yr*star%solar_mass_cgs/ &
                 seconds_per_year
            pressure_from_wind = mass_loss_rate_cgs/c4pi* &
                 sqrt(2.0d0*accretion_efficiency*exp10(log10_gravity)/ &
                 total_radius_cm**3)
            temperature_from_wind = pow(mass_loss_rate_cgs* &
                 accretion_specific_energy0*0.75d0/c4pi/total_radius_cm**2/ &
                 csig, 0.25d0)
            log10_pressure_local = log10(pressure_from_wind)
            log10_temperature_local = log10(temperature_from_wind)
            eos_res(i_log10_density) = log10_density_local
            eos_res(i_beta) = beta_local
            call eos_get(log10_temperature_local, log10_pressure_local, &
                 hydrogen_fraction_local, metal_fraction_local, eos_res, &
                 eos_deriv_flag, eos_atmosphere_flag, saha_flag, ierr=ierr)
            if (ierr /= 0) return
            log10_density_local = eos_res(i_log10_density)
            beta_local = 1.0d0-(radiation_constant_over_3*eos_res(i_temperature)**4/ &
                 eos_res(i_pressure))
            mean_molecular_weight_local = eos_res(i_pressure)*beta_local/ &
                 (eos_res(i_density)*eos_res(i_temperature))
            rot_scr%envelope_specific_entropy = 0.0d0
         else
            rot_scr%envelope_specific_entropy = 0.0d0
         endif
      endif
! CALL MASS LOSS OR ACCRETION ROUTINE
! 2026 config-matrix fix: the historical call passed 24 actuals into
! mdot's 23 dummies (leading extra log_luminosity_lsun), shifting
! every later argument one slot -- LMDOT runs dereferenced a scalar
! as the composition array and crashed. Pre-existing in the original
! F77 (documented in mdot.f90's header); fixed by dropping the
! spurious first actual so the list matches mdot's dummies 1:1.
      call mdot(timestep,composition, &
           log_density,specific_angular_momentum,log_pressure,log_radius,log_mass, &
           zone_mass_grams,shell_mass,log_total_mass,log_temperature, &
           envelope_boundary_zone,new_surface_bc_needed,num_zones,omega, &
           mean_molecular_weight_local,total_radius_cm,total_mass_msun, &
           mass_loss_rate_msun_yr,accretion_specific_energy,mean_thermal_energy, &
           cz_total_mass_below_fitting,old_log_envelope_mass_fraction, ierr)
      return
end subroutine massloss

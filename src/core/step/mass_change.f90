!----------------------------------------------------------------------
! mass_change
!----------------------------------------------------------------------
! The mass-change step of a timestep (formerly wind/massloss and
! wind/mdot, MHP 10/01): mass loss or accretion at the surface, when
! use_mass_accretion is set. The rate is mass_accretion_rate (Msun/yr;
! positive for accretion), or the Reimers law with use_reimers_wind, or
! the rotation-scaled solar wind (MHP 8/10). The step
!   - limits the timestep so that the step removes at most a fraction
!     of the star and of its surface convection zone (ftotdmdt,
!     fczdmdt) and adds at most half the mass outside the outer zone,
!     and ends it at the disk-locking age when accretion stops there;
!   - with accretion, adjusts R, P and T of the convection zone for the
!     entropy of the accreted matter, and mixes the accreted
!     composition into it;
!   - with mass loss from a rotating star, removes the angular momentum
!     the lost matter carries (a thermal wind);
!   - rescales the convection zone's masses and the total mass.
! The thermal-energy and entropy bookkeeping for accretion is
! experimental and valid for Allard atmospheres only (llp 06/15/2009).
! The model is read from and written to star; the timestep is star%dt.
module mass_change
      use atm_lib, only: atm_get_surface_pt
      use atmosphere_structure, only: atm_point
      use star_info_lib, only: star, json
      use physical_constants
      use eos_lib
      use math_lib
      implicit none
      private
      public :: change_mass

contains

!----------------------------------------------------------------------
! change_mass (formerly massloss)
!----------------------------------------------------------------------
! The rate, and for accretion the thermal energy and entropy of the
! convection zone and of the accreted matter; then apply_mass_change.
! old_log_envelope_mass_fraction: log10 of the mass inside the outer
! zone over the total mass, before the change (set only when the step
! runs). new_atmosphere_fit_needed: accretion was applied.
subroutine change_mass(old_log_envelope_mass_fraction, new_atmosphere_fit_needed, ierr)
      double precision, intent(out) :: old_log_envelope_mass_fraction
      logical, intent(out) :: new_atmosphere_fit_needed
      integer, intent(out) :: ierr
! MHP 5/02 EFFICIENCY FACTOR FOR THE THERMAL ENERGY CONTENT
! OF ACCRETED MATTER.
      double precision, parameter :: accretion_efficiency = 1.0d0
      logical :: disk_exhausted
      integer :: nz, jenv, zone_idx
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
      double precision :: log10_gravity
      logical :: allard_surface_failed
      double precision :: log10_temperature_local, log10_pressure_local, &
           log10_density_local
      double precision :: hydrogen_fraction_local, metal_fraction_local
      logical :: eos_deriv_flag, eos_atmosphere_flag
      integer :: saha_flag
      double precision :: beta_local
! Nothing here is SAVEd: the original's blanket SAVE carried the
! previous call's converged log10 rho / beta into the eos seeds; they
! start at zero (below), which is what the original saw on its first
! call.
      type(eos_result) :: eos_res
      double precision :: mass_loss_rate_cgs, pressure_from_wind, &
           temperature_from_wind

! INITIALIZE MASS LOSS AT DEFAULT RATE
      ierr = 0
      new_atmosphere_fit_needed = .false.
      log10_density_local = 0.0d0
      beta_local = 0.0d0
      mass_loss_rate_msun_yr = star%ctrl%mass_accretion_rate
      if(.not.star%job%use_mass_accretion) return
      nz = star%nz
      jenv = star%envelope_cz_bottom_index
      call set_mass_change_rate()
      call measure_convective_envelope()
      if(mass_loss_rate_msun_yr.gt.0.0d0)then
         new_atmosphere_fit_needed = .true.
         call average_envelope_energies()
         call find_accreted_surface_state()
         if (ierr /= 0) return
         call evaluate_stahler_surface()
         if (ierr /= 0) return
! the envelope entropy is reset after all, so the accretion step makes
! no entropy correction to R, P and T
         star%envelope_specific_entropy = 0.0d0
      endif
      call apply_mass_change(total_radius_cm, mass_loss_rate_msun_yr, &
           accretion_specific_energy, mean_thermal_energy, &
           cz_total_mass_below_fitting, old_log_envelope_mass_fraction, &
           disk_exhausted, ierr)
! the disk-locking age was reached this step: accretion stops for the
! rest of the run
      if(disk_exhausted) star%job%use_mass_accretion = .false.

contains

! ---------------------------------------------------------------
! The star's radius and mass in CGS, and the rate: the constant
! mass_accretion_rate, or the Reimers law with use_reimers_wind.
subroutine set_mass_change_rate()
! COMPUTE GLOBAL QUANTITIES (RADIUS,MASS) IN CGS UNITS.
      log10_radius = log10_radius_from_l_teff(star%log_L, star%log_Teff, &
           star%log10_solar_luminosity)
      total_radius_cm = exp10(log10_radius)
      total_mass_grams = star%star_mass*star%solar_mass_cgs
! USE A REIMERS FORMULA TO COMPUTE MDOT IF DESIRED; OVERWRITES
! CONSTANT MDOT: Mdot[Msun/yr] = reimers_scaling_factor*(L/Lsun)(R/Rsun)/(M/Msun),
! reimers_scaling_factor default -4e-13 (the classical eta=1 coefficient).
! (2026 config-matrix fix: the historical expression computed
! reimers_scaling_factor*L/(g*R) in CGS -- units of g/s -- but stored it
! in the Msun/yr variable, which mdot then converted to CGS AGAIN, a
! ~1e25 rate inflation. Unreachable before the mdot argument-count fix,
! so never caught.)
      if(star%ctrl%use_reimers_wind)then
         mass_loss_rate_msun_yr = star%ctrl%reimers_scaling_factor* &
              exp10(star%log_L)* &
              (total_radius_cm/star%solar_radius_cgs)/star%star_mass
      endif
end subroutine set_mass_change_rate

! ---------------------------------------------------------------
! With a surface convection zone, the turnover timescale is printed and
! jcz set; then the mass of the convective envelope.
subroutine measure_convective_envelope()
! 02/12 MHP TAUCZ NOW COMPUTED PRIOR TO CALL IN MIXCZ
      if(jenv.lt.nz)then
         write(*,*)star%convective_turnover_timescale/seconds_per_year, &
              total_radius_cm/star%solar_radius_cgs
         star%jcz = jenv
      else
         star%convective_turnover_timescale = 0.0d0
         star%jcz = nz
      endif
! MHP 8/10 THE MASS OF THE CONVECTIVE ENVELOPE
      cz_total_mass_below_fitting = 0.0d0
      do zone_idx = jenv, nz
         cz_total_mass_below_fitting = cz_total_mass_below_fitting+ &
              star%dm(zone_idx)
      end do
end subroutine measure_convective_envelope

! ---------------------------------------------------------------
! For accretion: the energy per gram the accreted matter brings, and the
! mass-weighted thermal energy and entropy of the convective envelope.
subroutine average_envelope_energies()
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
      star%envelope_specific_entropy = 0.0d0
      do zone_idx = jenv, nz
         local_temperature = exp10(star%logT(zone_idx))
         local_pressure = exp10(star%logP(zone_idx))
         local_density = exp10(star%logRho(zone_idx))
         local_radius_cm = exp10(star%logR(zone_idx))
         delta_accretion_energy = 0.5d0*accretion_specific_energy* &
              (total_radius_cm/local_radius_cm-1.0d0)
         thermal_energy_accreted = accretion_specific_energy+ &
              delta_accretion_energy
         thermal_energy_accreted_bar = thermal_energy_accreted_bar+ &
              star%dm(zone_idx)*thermal_energy_accreted
         local_beta = 1.0d0-(radiation_constant_over_3*local_temperature**4/ &
              local_pressure)
         mean_molecular_weight_local = local_pressure*local_beta/ &
              (local_density*local_temperature)
         local_entropy = mean_molecular_weight_local* &
              (1.5d0*log(local_temperature)-log(local_density))
         star%envelope_specific_entropy = star%envelope_specific_entropy+ &
              local_entropy*star%dm(zone_idx)
! THE THERMAL ENERGY PER GM IN THE JTH SHELL IS
         thermal_energy_per_gram = local_pressure*local_beta/local_density
! THE RUNNING TOTAL OF THERMAL ENERGY THROUGHOUT THE
! CONVECTIVE ENVELOPE IS
         sum_thermal_energy = sum_thermal_energy+thermal_energy_per_gram* &
              star%dm(zone_idx)
      end do
! THE AVERAGE THERMAL ENERGY OF THE UNPERTURBED CE (W/O ACCN) IS
      mean_thermal_energy = sum_thermal_energy/cz_total_mass_below_fitting
      accretion_specific_energy = thermal_energy_accreted_bar/ &
           cz_total_mass_below_fitting
      star%envelope_specific_entropy = star%envelope_specific_entropy/ &
           cz_total_mass_below_fitting
end subroutine average_envelope_energies

! ---------------------------------------------------------------
! The surface of an Allard atmosphere at the star's Teff and gravity, and
! the eos there: the state of the accreted matter.
subroutine find_accreted_surface_state()
! the entropy of the accreted matter, at the surface of an Allard
! atmosphere (when the lookup fails, atm_point keeps the last point)
      log10_gravity = cgl+star%log_total_mass-2.0d0*log10_radius
      call atm_get_surface_pt(4, star%log_Teff, log10_gravity, &
           atm_point%log10_pressure, atm_point%log10_temperature, &
           allard_surface_failed, ierr)
      if (ierr /= 0) return
      log10_temperature_local = atm_point%log10_temperature
      log10_pressure_local = atm_point%log10_pressure
      hydrogen_fraction_local = star%xa(1,nz)
      metal_fraction_local = star%xa(3,nz)
      eos_deriv_flag = .false.
      saha_flag = 1
      eos_atmosphere_flag = .true.
      eos_res%log10_density = log10_density_local
      eos_res%beta = beta_local
      call eos_get(log10_temperature_local, log10_pressure_local, &
           hydrogen_fraction_local, metal_fraction_local, eos_res, &
           eos_deriv_flag, eos_atmosphere_flag, saha_flag, ierr=ierr)
      if (ierr /= 0) return
      log10_density_local = eos_res%log10_density
      beta_local = 1.0d0-(radiation_constant_over_3*eos_res%temperature**4/ &
           eos_res%pressure)
      mean_molecular_weight_local = eos_res%pressure*beta_local/ &
           (eos_res%density*eos_res%temperature)
end subroutine find_accreted_surface_state

! ---------------------------------------------------------------
! The eos at the surface pressure and temperature of Stahler (1988) for
! the accretion rate (its results are unused; see inside).
subroutine evaluate_stahler_surface()
! ALTERNATE EXPRESSION FOR SURFACE PRESSURE AND LUMINOSITY, FROM STAHLER
! 1988. Nothing this second eos call returns is used (apply_mass_change
! recomputes mu before using it); the call is kept because the eos
! carries solver state from one call to the next.
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
         eos_res%log10_density = log10_density_local
         eos_res%beta = beta_local
         call eos_get(log10_temperature_local, log10_pressure_local, &
              hydrogen_fraction_local, metal_fraction_local, eos_res, &
              eos_deriv_flag, eos_atmosphere_flag, saha_flag, ierr=ierr)
         if (ierr /= 0) return
      endif
end subroutine evaluate_stahler_surface

end subroutine change_mass

!----------------------------------------------------------------------
! apply_mass_change (formerly mdot)
!----------------------------------------------------------------------
! Limits the timestep, then removes or adds the mass of the step at
! the surface: the R/P/T adjustment and mixing for accretion, the
! angular momentum carried off by mass loss, and the convection zone's
! masses. accretion_specific_energy and mean_thermal_energy are the
! thermal energies per gram of the accreted matter and of the
! convection zone, whose mass is cz_total_mass_below_fitting.
! disk_exhausted: the disk-locking age was reached (the timestep is
! clipped to land on it). ierr = 1: the envelope boundary is the
! surface zone.
subroutine apply_mass_change(total_radius_cm, mass_loss_rate_msun_yr, &
     accretion_specific_energy, mean_thermal_energy, &
     cz_total_mass_below_fitting, old_log_envelope_mass_fraction, &
     disk_exhausted, ierr)
      double precision, intent(in) :: total_radius_cm
      double precision, intent(inout) :: mass_loss_rate_msun_yr
      double precision, intent(in) :: accretion_specific_energy, &
           mean_thermal_energy
      double precision, intent(in) :: cz_total_mass_below_fitting
      double precision, intent(out) :: old_log_envelope_mass_fraction
      logical, intent(out) :: disk_exhausted
      integer, intent(out) :: ierr
      integer :: nz, jenv, zone_idx, species_idx
!      double precision :: omega_ratio_sq, omega_max_ratio_sq
      double precision :: mass_loss_rate_cgs
      double precision :: cz_mass_grams
      double precision :: timestep_limit_total_mass, timestep_limit_cz_mass, &
           timestep_limit
      double precision :: timestep_limit_envelope
      double precision :: disk_age_test
      double precision :: delta_mass_cgs
      double precision :: new_thermal_energy
      double precision :: local_temperature, local_pressure, local_density, &
           local_beta, mean_molecular_weight_local
      double precision :: delta_log_specific_entropy, delta_ln_mass, &
           delta_log_radius
      double precision :: boundary_radius_cm, radius_scale_factor, &
           radius_before_cm, radius_after_cm
      double precision :: surface_moment_of_inertia_per_mass, &
           delta_angular_momentum, surface_omega_local
      double precision :: cz_total_am, cz_new_am, am_ratio
      double precision :: mass_scale_factor
      double precision :: delta_mass_msun
      double precision :: total_mass_grams_old, total_mass_grams_new
      double precision :: mixed_abundance

      ierr = 0
      disk_exhausted = .false.
      nz = star%nz
      jenv = star%envelope_cz_bottom_index
      old_log_envelope_mass_fraction = star%log_mass(nz) - star%log_total_mass
! COMMENTED OUT (2026-10-12, the author's decision): the solar wind
! scaled by the surface rotation (MHP 8/10), Mdot = Mdot_sun
! min((omega/omega_ref)^2, (omega_max/omega_ref)^2). Nothing set its
! switch or its three parameters (the legacy LSOLWIND is read but never
! reached them), so it never ran. To revive it, uncomment this block,
! its locals and its star members (state/star_info_lib), and set them.
!      if(star%use_rotation_scaled_solar_wind .and. star%job%rotation_active) then
!         omega_ratio_sq = (star%omega(nz)/star%wind_reference_omega)**2
!         omega_max_ratio_sq = (star%wind_max_omega/star%wind_reference_omega)**2
!         mass_loss_rate_msun_yr = star%solar_wind_mass_loss_rate_msun_yr* &
!              min(omega_ratio_sq,omega_max_ratio_sq)
!         write(*,*)star%omega(nz),mass_loss_rate_msun_yr
!      endif
! CONVERT FROM SOLAR MASSES/YEAR TO GM/SEC
      mass_loss_rate_cgs = abs(mass_loss_rate_msun_yr)*star%solar_mass_cgs/ &
           seconds_per_year
      call weigh_surface_convection_zone()
      if (ierr /= 0) return
      call limit_mass_change_timestep()
! FINAL AMOUNT OF MASS LOSS INFERRED IN CGS UNITS.
      delta_mass_cgs = mass_loss_rate_msun_yr*star%solar_mass_cgs*star%dt/ &
           seconds_per_year
      call adjust_cz_for_accreted_entropy()
      call remove_wind_angular_momentum()
      call rescale_cz_masses()
      call update_total_mass()
      call mix_in_accreted_composition()

contains

! ---------------------------------------------------------------
! The mass of the surface convection zone; ierr = 1 when the envelope
! boundary is the surface zone.
subroutine weigh_surface_convection_zone()
! COMPUTE MASS OF SURFACE CONVECTION ZONE BELOW FITTING POINT
      if(jenv.eq.nz) then
         write(*,*) 'mdot: envelope boundary at the surface (911)'
         ierr = 1
         return
      end if
      if(jenv.gt.1) then
         cz_mass_grams = star%m(nz)-star%m(jenv)
      else
         cz_mass_grams = star%m(nz)
      endif
end subroutine weigh_surface_convection_zone

! ---------------------------------------------------------------
! Limits the timestep for the mass change (ftotdmdt, fczdmdt, half the
! outer mass) and, for accretion, ends it at the disk-locking age.
subroutine limit_mass_change_timestep()
! COMPUTE MAXIMUM TIMESTEP BASED ON NOT REMOVING TOO MUCH MASS FROM THE
! SURFACE CZ (FCZDMDT) OR AS A FRACTION OF THE TOTAL MASS (FTOTDMDT)
      timestep_limit_total_mass = star%ctrl%ftotdmdt*star%star_mass*star%solar_mass_cgs/ &
           mass_loss_rate_cgs
      timestep_limit_cz_mass = star%ctrl%fczdmdt*cz_mass_grams/mass_loss_rate_cgs
      timestep_limit = min(timestep_limit_total_mass,timestep_limit_cz_mass)
! RESTRICT TIMESTEP TO ADD NO MORE THAN 1/2 OF THE CURRENT MASS
! BEYOND THE FITTING POINT TO THE STAR.
      if(mass_loss_rate_cgs.gt.0.0d0)then
         timestep_limit_envelope = 0.5d0*(exp10(star%log_total_mass)- &
              exp10(star%log_mass(nz)))/mass_loss_rate_cgs
         timestep_limit = min(timestep_limit,timestep_limit_envelope)
      endif
      if(timestep_limit.lt.star%dt)then
         write(*,10)star%dt,timestep_limit
 10      format(' TIMESTEP REDUCED FOR MASS LOSS - OLD,NEW',1P2E12.3)
         star%dt = timestep_limit
      endif
! mhp 8/10 turn mass loss off when disk exhausted only when dm /dt > 0, e.g. accretion
      if(star%job%disk_locking_active .and. mass_loss_rate_msun_yr.gt.0.0d0)then
         disk_age_test = star%disk_gate_age_gyr + 1.0d-9*star%dt/seconds_per_year
         if(disk_age_test.gt.star%job%disk_locking_age_gyr)then
            star%dt = (star%job%disk_locking_age_gyr-star%disk_gate_age_gyr)*1.0d9*seconds_per_year
            disk_exhausted = .true.
         endif
      endif
end subroutine limit_mass_change_timestep

! ---------------------------------------------------------------
! With accretion (delta_mass_cgs > 0): the convection zone's R, P and T
! shifted for the entropy of the accreted matter.
subroutine adjust_cz_for_accreted_entropy()
! COMPUTE THE MEAN THERMAL ENERGY
! CONTENT OF THE CONVECTION ZONE FOR MODELS WITH ACCRETION.
      if(delta_mass_cgs.gt.0.0d0)then
         new_thermal_energy = (accretion_specific_energy*delta_mass_cgs+ &
              mean_thermal_energy*cz_total_mass_below_fitting)/ &
              (delta_mass_cgs+cz_total_mass_below_fitting)
! OVERALL SCALE FACTOR IN R
         local_temperature = exp10(star%logT(jenv))
         local_pressure = exp10(star%logP(jenv))
         local_density = exp10(star%logRho(jenv))
         local_beta = 1.0d0-(radiation_constant_over_3*local_temperature**4/ &
              local_pressure)
         mean_molecular_weight_local = local_pressure*local_beta/ &
              (local_density*local_temperature)
         delta_log_specific_entropy = (delta_mass_cgs/cz_mass_grams)* &
              star%envelope_specific_entropy/mean_molecular_weight_local
! THE MASS TERM ln((M_cz+dM)/M_cz) IS DELIBERATELY NOT APPLIED.
         delta_ln_mass = 0.0d0
         delta_log_radius = (cc23*delta_log_specific_entropy- &
              cc13*delta_ln_mass)/ln10
         star%delta_log_temperature = delta_ln_mass/ln10 - delta_log_radius
         star%delta_log_pressure = 2.0d0*delta_ln_mass/ln10 - 4.0d0*delta_log_radius
         if(jenv.eq.1)then
            do zone_idx = jenv,nz
               star%logR(zone_idx) = star%logR(zone_idx)+delta_log_radius
               star%logP(zone_idx) = star%logP(zone_idx)+star%delta_log_pressure
               star%logT(zone_idx) = star%logT(zone_idx)+star%delta_log_temperature
            end do
         else
            star%logP(jenv) = star%logP(jenv)+star%delta_log_pressure
            star%logT(jenv) = star%logT(jenv)+star%delta_log_temperature
            boundary_radius_cm = exp10(star%logR(jenv))
            radius_scale_factor = exp10(delta_log_radius)
            do zone_idx = jenv+1,nz
               radius_before_cm = exp10(star%logR(zone_idx))
               radius_after_cm = boundary_radius_cm+radius_scale_factor* &
                    (radius_before_cm-boundary_radius_cm)
               star%logR(zone_idx) = log10(radius_after_cm)
            end do
         endif
      endif
end subroutine adjust_cz_for_accreted_entropy

! ---------------------------------------------------------------
! With mass loss from a rotating star: the convection zone loses the
! angular momentum the lost matter carries (a thermal wind).
subroutine remove_wind_angular_momentum()
! FOR MODELS WITH ROTATION AND MASS LOSS REMOVE ANGULAR MOMENTUM.
! IN THIS CASE WE ASSUME A THERMAL WIND WHERE THE SURFACE MATTER
! CARRIES AWAY ONLY ITS LOCAL ANGULAR MOMENTUM PER UNIT MASS.
      if(star%job%rotation_active .and. delta_mass_cgs.lt.0.0d0)then
! MOMENT OF INERTIA PER UNIT MASS AT THE SURFACE.
         surface_moment_of_inertia_per_mass = 2.0d0/3.0d0*total_radius_cm**2
         if(star%ctrl%walpcz.ge.0.0d0)then
! SOLID BODY CZ ROTATION
            delta_angular_momentum = star%omega(nz)* &
                 surface_moment_of_inertia_per_mass*delta_mass_cgs
         else if(star%ctrl%walpcz.le.-2.0d0)then
! CONSTANT J/M
            delta_angular_momentum = star%j_rot(nz)*delta_mass_cgs
         else
            surface_omega_local = star%omega(nz)*exp10( &
                 star%logR(nz)*star%ctrl%walpcz)/pow(total_radius_cm, star%ctrl%walpcz)
            delta_angular_momentum = surface_omega_local* &
                 surface_moment_of_inertia_per_mass*delta_mass_cgs
         endif
!  FIND TOTAL STARTING ANGULAR MOMENTUM OF C.Z.
         cz_total_am = star%j_rot(jenv)*star%dm(jenv)
         do zone_idx = jenv+1,nz
            cz_total_am = cz_total_am + star%j_rot(zone_idx)*star%dm(zone_idx)
         end do
         cz_new_am = (cz_total_am+delta_angular_momentum)
         am_ratio = cz_new_am/cz_total_am
         write(*,22)delta_mass_cgs,delta_angular_momentum,star%omega(nz), &
              cz_new_am,cz_total_am,cz_total_mass_below_fitting
 22      format(1P6E13.4)
         do zone_idx = jenv,nz
            star%j_rot(zone_idx) = star%j_rot(zone_idx)*am_ratio
         end do
      endif
end subroutine remove_wind_angular_momentum

! ---------------------------------------------------------------
! With mass loss: the masses above the envelope boundary rescaled to the
! convection zone's new mass.
subroutine rescale_cz_masses()
! REMOVE OR ADD MASS IN THE CONVECTION ZONE: rescale the zone-center
! masses m (and log_mass) above the envelope boundary.
      if(delta_mass_cgs.lt.0.0d0)then
         mass_scale_factor = (cz_mass_grams+delta_mass_cgs)/cz_mass_grams
         do zone_idx = jenv+1,nz
            star%m(zone_idx) = star%m(jenv)+ &
                 mass_scale_factor*(star%m(zone_idx)-star%m(jenv))
            star%log_mass(zone_idx) = log10(star%m(zone_idx))
         end do
      endif
end subroutine rescale_cz_masses

! ---------------------------------------------------------------
! The new total mass, the shell masses dm of the convection zone and
! senv; the surface boundary condition is to be recomputed.
subroutine update_total_mass()
! CORRECT TOTAL MASS IN SOLAR UNITS AND ITS LOG IN GRAMS
      star%star_mass = star%star_mass + delta_mass_cgs/star%solar_mass_cgs
      delta_mass_msun = delta_mass_cgs/star%solar_mass_cgs
      write(*,20)star%star_mass,delta_mass_msun
 20   format('MASS LOSS APPLIED - NEW M,DEL M',1P2E19.10)
      total_mass_grams_old = exp10(star%log_total_mass)
      total_mass_grams_new = total_mass_grams_old + delta_mass_cgs
      star%log_total_mass = log10(total_mass_grams_new)
      star%stotal = star%log_total_mass
! CORRECT MASS CONTENTS OF INDIVIDUAL SHELLS (dm, IN GM)
      do zone_idx = jenv,nz-1
         star%dm(zone_idx) = 0.5d0*(star%m(zone_idx+1)-star%m(zone_idx-1))
      end do
      star%dm(nz) = total_mass_grams_new-0.5d0*(star%m(nz)+star%m(nz-1))
! 07/02RESET SENV
      star%senv = star%log_mass(nz) - star%log_total_mass
! RECOMPUTE SURFACE BOUNDRY CONDITION
      star%recompute_envelope_triangle = .true.
end subroutine update_total_mass

! ---------------------------------------------------------------
! With accretion: the accreted composition mixed into the convection
! zone, and the accreted mass recorded.
subroutine mix_in_accreted_composition()
! REMIX THE SURFACE CONVECTION ZONE IF MDOT IS POSITIVE.
      if(delta_mass_cgs.gt.0.0d0)then
         do species_idx = 1,11
            mixed_abundance = (star%xa(species_idx,jenv)* &
                 cz_total_mass_below_fitting+star%ctrl%accreted_composition(species_idx)* &
                 delta_mass_cgs)/(delta_mass_cgs+cz_total_mass_below_fitting)
            do zone_idx = jenv,nz
               star%xa(species_idx,zone_idx) = mixed_abundance
            end do
         end do
         star%accreted_mass = delta_mass_cgs
      endif
end subroutine mix_in_accreted_composition

end subroutine apply_mass_change

end module mass_change

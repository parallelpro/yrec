!----------------------------------------------------------------------
! mass_change
!----------------------------------------------------------------------
! The mass-change step of a timestep (formerly wind/massloss and
! wind/mdot, MHP 10/01): mass loss or accretion at the surface, when
! use_mass_accretion is set. The rate (Msun/yr, positive for accretion)
! is mass_accretion_rate, or the Reimers law with use_reimers_wind. (A
! solar wind scaled by the surface rotation, MHP 8/10, is commented out
! in apply_mass_change.) The step
!   - limits the timestep so that |dM| is at most the fraction ftotdmdt
!     of the star, the fraction fczdmdt of its surface convection zone
!     and half the mass outside the outer zone, and with accretion ends
!     the step at the disk-locking age, where the mass change stops;
!   - with accretion, shifts R, P and T of the convection zone for the
!     entropy of the accreted matter (a zero shift at present, see
!     change_mass) and mixes the accreted composition into the zone;
!   - with mass loss from a rotating star, removes the angular momentum
!     the lost matter carries (a thermal wind);
!   - with mass loss, compresses the zone masses of the convection
!     zone; then sets the new total mass and the shell masses dm.
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
! Sets the rate and, for accretion, the energy per gram of the accreted
! matter and the mean thermal energy and entropy of the convective
! envelope; then calls apply_mass_change. Outputs:
! old_log_envelope_mass_fraction = log10(m(nz)/M) before the change
! (set only when the step runs), and new_atmosphere_fit_needed, true
! when matter was accreted.
subroutine change_mass(old_log_envelope_mass_fraction, new_atmosphere_fit_needed, ierr)
      double precision, intent(out) :: old_log_envelope_mass_fraction
      logical, intent(out) :: new_atmosphere_fit_needed
      integer, intent(out) :: ierr
! The efficiency of the accreted matter's thermal energy: it brings
! accretion_efficiency G M/R per gram (MHP 5/02).
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
! Nothing here is saved between calls: the eos seeds log10 rho and beta
! start at zero on every call (below). The legacy blanket SAVE carried
! the previous call's values; zero is what it saw on its first call.
      type(eos_result) :: eos_res
      double precision :: mass_loss_rate_cgs, pressure_from_wind, &
           temperature_from_wind

! Start from the constant rate mass_accretion_rate (Msun/yr); the step
! runs only with use_mass_accretion.
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
! The envelope entropy is zeroed after all, so apply_mass_change shifts
! R, P and T by zero: the entropy correction is switched off.
         star%envelope_specific_entropy = 0.0d0
      endif
      call apply_mass_change(total_radius_cm, mass_loss_rate_msun_yr, &
           accretion_specific_energy, mean_thermal_energy, &
           cz_total_mass_below_fitting, old_log_envelope_mass_fraction, &
           disk_exhausted, ierr)
! The disk-locking age was reached this step: the mass change stops for
! the rest of the run.
      if(disk_exhausted) star%job%use_mass_accretion = .false.

contains

! ---------------------------------------------------------------
! The star's radius and mass in CGS, and the rate: the constant
! mass_accretion_rate, or the Reimers law with use_reimers_wind.
subroutine set_mass_change_rate()
! R from L = 4 pi R^2 sigma Teff^4, and M, in CGS.
      log10_radius = log10_radius_from_l_teff(star%log_L, star%log_Teff, &
           star%log10_solar_luminosity)
      total_radius_cm = exp10(log10_radius)
      total_mass_grams = star%star_mass*star%solar_mass_cgs
! The Reimers law replaces the constant rate:
!   Mdot [Msun/yr] = reimers_scaling_factor (L/Lsun) (R/Rsun)/(M/Msun),
! with reimers_scaling_factor = -4e-13 by default (eta = 1, mass loss).
! (2026 fix: the legacy code evaluated reimers_scaling_factor L/(g R) in
! CGS, a rate in g/s, but stored it as Msun/yr, and mdot converted it
! to CGS again, a rate about 1e25 too large.)
      if(star%ctrl%use_reimers_wind)then
         mass_loss_rate_msun_yr = star%ctrl%reimers_scaling_factor* &
              exp10(star%log_L)* &
              (total_radius_cm/star%solar_radius_cgs)/star%star_mass
      endif
end subroutine set_mass_change_rate

! ---------------------------------------------------------------
! With a surface convection zone, sets jcz to the zone's base (and, with
! the solver diagnostics on, prints the turnover timescale in yr and
! R/Rsun); without one, zeroes the turnover timescale and sets jcz = nz.
! Then sums the mass of the convective envelope.
subroutine measure_convective_envelope()
      use run_log, only: solver_diagnostics
! The turnover timescale was computed beforehand, for the previous model
! (compute_turnover_timescale; MHP 2/12). (The F77 code printed it
! unlabelled every step; it is now a diagnostic.)
      if(jenv.lt.nz)then
         if (solver_diagnostics()) write(*,'(a,1p,2e12.4)') &
              ' mass change: turnover time (yr), R/Rsun', &
              star%convective_turnover_timescale/seconds_per_year, &
              total_radius_cm/star%solar_radius_cgs
         star%jcz = jenv
      else
         star%convective_turnover_timescale = 0.0d0
         star%jcz = nz
      endif
! The mass of the convective envelope, sum of dm(jenv..nz) (MHP 8/10).
! dm(nz) reaches the total mass, so this includes the mass outside the
! fitting point, despite the name.
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
! By the virial theorem at most half of the gravitational potential
! energy of the accreted matter becomes kinetic energy. The energy per
! gram the accreted matter brings is taken as
!   e_acc = accretion_efficiency G M/R
! (twice that bound, G M/(2R), with accretion_efficiency = 1), and is
! compared with the mean thermal energy per gram of the convection zone.
      accretion_specific_energy = accretion_efficiency*total_mass_grams* &
           (exp10(cgl))/total_radius_cm
      accretion_specific_energy0 = accretion_specific_energy
! The mass-weighted means over the convective envelope (zones jenv..nz)
! of the accreted matter's energy per gram mixed down to radius r,
!   e(r) = e_acc [1 + (R/r - 1)/2]
! (it gains half the further potential energy released on the way
! down), of the thermal energy per gram and of the specific entropy.
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
! beta = P_gas/P = 1 - a T^4/(3 P). mean_molecular_weight_local holds
! P_gas/(rho T) = k/(mu m_u), the gas constant over mu, not mu itself.
! The ideal-gas specific entropy, up to a constant, is
!   s = (k/(mu m_u)) (3/2 ln T - ln rho).
         local_beta = 1.0d0-(radiation_constant_over_3*local_temperature**4/ &
              local_pressure)
         mean_molecular_weight_local = local_pressure*local_beta/ &
              (local_density*local_temperature)
         local_entropy = mean_molecular_weight_local* &
              (1.5d0*log(local_temperature)-log(local_density))
         star%envelope_specific_entropy = star%envelope_specific_entropy+ &
              local_entropy*star%dm(zone_idx)
! The thermal energy per gram of the shell, taken as P_gas/rho =
! kT/(mu m_u) (2/3 of the monatomic ideal-gas value 3/2 P_gas/rho).
         thermal_energy_per_gram = local_pressure*local_beta/local_density
         sum_thermal_energy = sum_thermal_energy+thermal_energy_per_gram* &
              star%dm(zone_idx)
      end do
! The means: the thermal energy per gram of the convective envelope
! before accretion, the accreted matter's energy per gram, and the
! specific entropy.
      mean_thermal_energy = sum_thermal_energy/cz_total_mass_below_fitting
      accretion_specific_energy = thermal_energy_accreted_bar/ &
           cz_total_mass_below_fitting
      star%envelope_specific_entropy = star%envelope_specific_entropy/ &
           cz_total_mass_below_fitting
end subroutine average_envelope_energies

! ---------------------------------------------------------------
! The surface of an Allard atmosphere at the star's Teff and gravity,
! and the eos there, meant as the state of the accreted matter. Only the
! eos seeds (log10 rho, beta) and the inputs evaluate_stahler_surface
! reuses are used afterwards: no entropy of the accreted matter is
! formed from this state.
subroutine find_accreted_surface_state()
! The Allard tables (atm_choice 4) at log g = log G + log M - 2 log R.
! The lookup writes atm_point (atmosphere_structure's current point);
! when it fails, atm_point keeps its last values (allard_surface_failed
! is not checked).
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
! The new seeds; beta and P_gas/(rho T) at the surface (the latter is
! not used afterwards).
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
! An alternative surface pressure and temperature, from Stahler (1988):
! the ram pressure and the radiating temperature of matter falling
! freely onto the star,
!   P_s = Mdot v_ff/(4 pi R^2),  v_ff = sqrt(2 accretion_efficiency g R),
!   T_s^4 = (3/4) Mdot e_acc/(4 pi R^2 sigma),  e_acc as above.
! Nothing this second eos call returns is used (apply_mass_change
! recomputes P_gas/(rho T) before using it); the call is kept because
! the eos carries solver state (the Saha solver's starting guesses)
! from one call to the next.
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
! Limits the timestep, then removes or adds the step's mass,
! dM = Mdot dt, at the surface: the R/P/T shift and the mixing for
! accretion, the angular momentum carried off by mass loss, and the
! convection zone's masses. Inputs: total_radius_cm (R), the rate
! mass_loss_rate_msun_yr (Msun/yr, positive for accretion),
! accretion_specific_energy and mean_thermal_energy (the energies per
! gram of the accreted matter and of the convection zone), and
! cz_total_mass_below_fitting (the convective envelope's mass).
! Outputs: old_log_envelope_mass_fraction = log10(m(nz)/M) before the
! change; disk_exhausted, true when the disk-locking age was reached
! (the timestep is clipped to land on it); ierr = 1 when the envelope
! boundary is the surface zone.
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
! Commented out (2026-10-12, by the author's decision): the solar wind
! scaled by the surface rotation (MHP 8/10),
!   Mdot = Mdot_sun min((omega/omega_ref)^2, (omega_max/omega_ref)^2).
! Nothing set its switch or its three parameters (the legacy LSOLWIND
! is read but never reached them), so it never ran. To revive it,
! uncomment this block, its locals and its star members
! (state/star_info_lib), and set them.
!      if(star%use_rotation_scaled_solar_wind .and. star%job%rotation_active) then
!         omega_ratio_sq = (star%omega(nz)/star%wind_reference_omega)**2
!         omega_max_ratio_sq = (star%wind_max_omega/star%wind_reference_omega)**2
!         mass_loss_rate_msun_yr = star%solar_wind_mass_loss_rate_msun_yr* &
!              min(omega_ratio_sq,omega_max_ratio_sq)
!         write(*,*)star%omega(nz),mass_loss_rate_msun_yr
!      endif
! The magnitude of the rate, in g/s.
      mass_loss_rate_cgs = abs(mass_loss_rate_msun_yr)*star%solar_mass_cgs/ &
           seconds_per_year
      call weigh_surface_convection_zone()
      if (ierr /= 0) return
      call limit_mass_change_timestep()
! The mass change of the step, dM = Mdot dt, in g (negative for loss).
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
! The zone's mass below the fitting point, cz_mass_grams =
! m(nz) - m(jenv), or m(nz) when the zone reaches the center (jenv = 1).
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
! The largest timestep that changes the mass by at most the fraction
! ftotdmdt of the star or fczdmdt of the surface convection zone,
!   dt_max = min(ftotdmdt M, fczdmdt M_cz)/|Mdot|.
! A zero rate sets no limit (it used to divide by zero, which traps in
! the devel builds).
      if(mass_loss_rate_cgs.gt.0.0d0)then
         timestep_limit_total_mass = star%ctrl%ftotdmdt*star%star_mass*star%solar_mass_cgs/ &
              mass_loss_rate_cgs
         timestep_limit_cz_mass = star%ctrl%fczdmdt*cz_mass_grams/mass_loss_rate_cgs
         timestep_limit = min(timestep_limit_total_mass,timestep_limit_cz_mass)
! Nor by more than half the mass outside the fitting point (the outer
! zone), (M - m(nz))/2, for loss and accretion alike.
         timestep_limit_envelope = 0.5d0*(exp10(star%log_total_mass)- &
              exp10(star%log_mass(nz)))/mass_loss_rate_cgs
         timestep_limit = min(timestep_limit,timestep_limit_envelope)
         if(timestep_limit.lt.star%dt)then
            write(*,10)star%dt,timestep_limit
 10         format(' TIMESTEP REDUCED FOR MASS LOSS - OLD,NEW',1P2E12.3)
            star%dt = timestep_limit
         endif
      endif
! With accretion only (Mdot > 0, MHP 8/10): a step that would pass the
! disk-locking age is clipped to end there, and the mass change stops.
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
! shifted for the entropy of the accreted matter. change_mass zeroes
! envelope_specific_entropy first, so every shift is zero at present.
subroutine adjust_cz_for_accreted_entropy()
! The mean thermal energy per gram of the convection zone once the
! accreted matter has mixed in (new_thermal_energy, not used further).
      if(delta_mass_cgs.gt.0.0d0)then
         new_thermal_energy = (accretion_specific_energy*delta_mass_cgs+ &
              mean_thermal_energy*cz_total_mass_below_fitting)/ &
              (delta_mass_cgs+cz_total_mass_below_fitting)
! The overall scale factor in R. P_gas/(rho T) = k/(mu m_u) at the base
! turns the envelope's mean entropy <s> into the dimensionless change
!   d s = (dM/M_cz) <s>/(k/(mu m_u)),
! and the homology of an n = 3/2 polytrope (R ~ K M^(-1/3), with
! ln K = (2/3) s + const for P = K rho^(5/3)) gives
!   d log R = (2/3 d s - 1/3 d ln M)/ln 10,
!   d log T = d ln M/ln 10 - d log R,  d log P = 2 d ln M/ln 10 - 4 d log R.
         local_temperature = exp10(star%logT(jenv))
         local_pressure = exp10(star%logP(jenv))
         local_density = exp10(star%logRho(jenv))
         local_beta = 1.0d0-(radiation_constant_over_3*local_temperature**4/ &
              local_pressure)
         mean_molecular_weight_local = local_pressure*local_beta/ &
              (local_density*local_temperature)
         delta_log_specific_entropy = (delta_mass_cgs/cz_mass_grams)* &
              star%envelope_specific_entropy/mean_molecular_weight_local
! The mass term d ln M = ln((M_cz + dM)/M_cz) is deliberately left out.
         delta_ln_mass = 0.0d0
         delta_log_radius = (cc23*delta_log_specific_entropy- &
              cc13*delta_ln_mass)/ln10
         star%delta_log_temperature = delta_ln_mass/ln10 - delta_log_radius
         star%delta_log_pressure = 2.0d0*delta_ln_mass/ln10 - 4.0d0*delta_log_radius
! A fully convective star (jenv = 1) shifts every zone. Otherwise only
! the base zone's P and T shift, and the radii above it stretch about
! the base, r' = r_b + 10^(d log R) (r - r_b).
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
! A thermal wind: the lost matter carries away only its own specific
! angular momentum, that of the surface under the convection zone's
! rotation law (walpcz), dJ = j_surface dM (dM < 0).
      if(star%job%rotation_active .and. delta_mass_cgs.lt.0.0d0)then
! The moment of inertia per unit mass of a thin shell at the surface,
! (2/3) R^2.
         surface_moment_of_inertia_per_mass = 2.0d0/3.0d0*total_radius_cm**2
         if(star%ctrl%walpcz.ge.0.0d0)then
! A solid-body convection zone: j_surface = (2/3) R^2 omega(nz).
            delta_angular_momentum = star%omega(nz)* &
                 surface_moment_of_inertia_per_mass*delta_mass_cgs
         else if(star%ctrl%walpcz.le.-2.0d0)then
! Constant specific angular momentum: j_surface = j_rot(nz).
            delta_angular_momentum = star%j_rot(nz)*delta_mass_cgs
         else
! A power law, omega extrapolated from the outer zone to R as
! omega_s = omega(nz) (r(nz)/R)^walpcz, i.e. omega r^walpcz constant
! (the rotation module's law is omega ~ r^walpcz, the inverse).
            surface_omega_local = star%omega(nz)*exp10( &
                 star%logR(nz)*star%ctrl%walpcz)/pow(total_radius_cm, star%ctrl%walpcz)
            delta_angular_momentum = surface_omega_local* &
                 surface_moment_of_inertia_per_mass*delta_mass_cgs
         endif
! The convection zone's angular momentum, J_cz = sum of j_rot dm over
! jenv..nz; every zone's j_rot is then scaled by (J_cz + dJ)/J_cz.
         cz_total_am = star%j_rot(jenv)*star%dm(jenv)
         do zone_idx = jenv+1,nz
            cz_total_am = cz_total_am + star%j_rot(zone_idx)*star%dm(zone_idx)
         end do
         cz_new_am = (cz_total_am+delta_angular_momentum)
         am_ratio = cz_new_am/cz_total_am
! The run log: dM, dJ, omega(nz), the new and old J_cz, and the
! convective envelope's mass.
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
! The zone-center masses above the base (and their logs) compress
! about m(jenv): m(k) = m(jenv) + f (m(k) - m(jenv)) for k > jenv, with
! f = (M_cz + dM)/M_cz and M_cz = cz_mass_grams. Accreted mass is not
! spread over the zones; it lands in dm(nz) (update_total_mass).
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
! The new total mass, in Msun (star_mass) and as log10 M in g
! (log_total_mass and its copy stotal).
      star%star_mass = star%star_mass + delta_mass_cgs/star%solar_mass_cgs
      delta_mass_msun = delta_mass_cgs/star%solar_mass_cgs
      write(*,20)star%star_mass,delta_mass_msun
 20   format('MASS LOSS APPLIED - NEW M,DEL M',1P2E19.10)
      total_mass_grams_old = exp10(star%log_total_mass)
      total_mass_grams_new = total_mass_grams_old + delta_mass_cgs
      star%log_total_mass = log10(total_mass_grams_new)
      star%stotal = star%log_total_mass
! The shell masses dm (g) of the convection zone: half the span between
! the neighbouring zone centers, the last shell reaching the total mass.
! A zone reaching the center (jenv = 1) starts from it: dm(1) =
! (m(1) + m(2))/2, i.e. m_0 = -m_1 as in rezoning. (It used to read
! m(0), outside the array.)
      do zone_idx = max(jenv,2),nz-1
         star%dm(zone_idx) = 0.5d0*(star%m(zone_idx+1)-star%m(zone_idx-1))
      end do
      if(jenv.eq.1) star%dm(1) = 0.5d0*(star%m(2)+star%m(1))
      star%dm(nz) = total_mass_grams_new-0.5d0*(star%m(nz)+star%m(nz-1))
! senv = log10(m(nz)/M), the mass fraction inside the fitting point.
      star%senv = star%log_mass(nz) - star%log_total_mass
! The envelope triangle (the surface boundary condition) is to be
! recomputed.
      star%recompute_envelope_triangle = .true.
end subroutine update_total_mass

! ---------------------------------------------------------------
! With accretion: the accreted composition mixed into the convection
! zone, and the accreted mass recorded.
subroutine mix_in_accreted_composition()
! For each of the 11 basic species, with M_cz = cz_total_mass_below_fitting
! and X_acc = accreted_composition,
!   X = (X(jenv) M_cz + X_acc dM)/(M_cz + dM)  in every zone jenv..nz.
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

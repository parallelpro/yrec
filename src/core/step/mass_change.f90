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
! The convection zone's new masses, the wind's angular momentum and the
! accreted composition are computed by the star-blind kernels of
! mass_change_kernel.f90.
module mass_change
      use atm_lib, only: atm_get_surface_pt
      use atmosphere_structure, only: atm_point
      use star_info_lib, only: json, star_info
      use mass_change_kernel, only: wind_angular_momentum, &
           convection_zone_mass, change_shell_masses, &
           wind_specific_angular_momentum, remove_wind_angular_momentum, &
           mix_accreted_composition
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
!   star reads: ctrl%mass_accretion_rate, ctrl%reimers_scaling_factor,
!     ctrl%use_reimers_wind, dm, envelope_cz_bottom_index,
!     log10_solar_luminosity, log_L, log_Teff, log_total_mass, logP,
!     logR, logRho, logT, nz, solar_mass_cgs, solar_radius_cgs,
!     star_mass, xa
!   star writes: convective_turnover_timescale,
!     envelope_specific_entropy, jcz, job%use_mass_accretion
subroutine change_mass(s, old_log_envelope_mass_fraction, new_atmosphere_fit_needed, ierr)
      type(star_info), intent(inout) :: s
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
      mass_loss_rate_msun_yr = s%ctrl%mass_accretion_rate
      if(.not.s%job%use_mass_accretion) return
      nz = s%nz
      jenv = s%envelope_cz_bottom_index
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
         s%envelope_specific_entropy = 0.0d0
      endif
      call apply_mass_change(s, total_radius_cm, mass_loss_rate_msun_yr, &
           accretion_specific_energy, mean_thermal_energy, &
           cz_total_mass_below_fitting, old_log_envelope_mass_fraction, &
           disk_exhausted, ierr)
! The disk-locking age was reached this step: the mass change stops for
! the rest of the run.
      if(disk_exhausted) s%job%use_mass_accretion = .false.

contains

! ---------------------------------------------------------------
! The star's radius and mass in CGS, and the rate: the constant
! mass_accretion_rate, or the Reimers law with use_reimers_wind.
subroutine set_mass_change_rate()
! R from L = 4 pi R^2 sigma Teff^4, and M, in CGS.
      log10_radius = log10_radius_from_l_teff(s%log_L, s%log_Teff, &
           s%log10_solar_luminosity)
      total_radius_cm = exp10(log10_radius)
      total_mass_grams = s%star_mass*s%solar_mass_cgs
! The Reimers law replaces the constant rate:
!   Mdot [Msun/yr] = reimers_scaling_factor (L/Lsun) (R/Rsun)/(M/Msun),
! with reimers_scaling_factor = -4e-13 by default (eta = 1, mass loss).
! (2026 fix: the legacy code evaluated reimers_scaling_factor L/(g R) in
! CGS, a rate in g/s, but stored it as Msun/yr, and mdot converted it
! to CGS again, a rate about 1e25 too large.)
      if(s%ctrl%use_reimers_wind)then
         mass_loss_rate_msun_yr = s%ctrl%reimers_scaling_factor* &
              exp10(s%log_L)* &
              (total_radius_cm/s%solar_radius_cgs)/s%star_mass
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
              s%convective_turnover_timescale/seconds_per_year, &
              total_radius_cm/s%solar_radius_cgs
         s%jcz = jenv
      else
         s%convective_turnover_timescale = 0.0d0
         s%jcz = nz
      endif
! The mass of the convective envelope, sum of dm(jenv..nz) (MHP 8/10).
! dm(nz) reaches the total mass, so this includes the mass outside the
! fitting point, despite the name.
      cz_total_mass_below_fitting = 0.0d0
      do zone_idx = jenv, nz
         cz_total_mass_below_fitting = cz_total_mass_below_fitting+ &
              s%dm(zone_idx)
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
      s%envelope_specific_entropy = 0.0d0
      do zone_idx = jenv, nz
         local_temperature = exp10(s%logT(zone_idx))
         local_pressure = exp10(s%logP(zone_idx))
         local_density = exp10(s%logRho(zone_idx))
         local_radius_cm = exp10(s%logR(zone_idx))
         delta_accretion_energy = 0.5d0*accretion_specific_energy* &
              (total_radius_cm/local_radius_cm-1.0d0)
         thermal_energy_accreted = accretion_specific_energy+ &
              delta_accretion_energy
         thermal_energy_accreted_bar = thermal_energy_accreted_bar+ &
              s%dm(zone_idx)*thermal_energy_accreted
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
         s%envelope_specific_entropy = s%envelope_specific_entropy+ &
              local_entropy*s%dm(zone_idx)
! The thermal energy per gram of the shell, taken as P_gas/rho =
! kT/(mu m_u) (2/3 of the monatomic ideal-gas value 3/2 P_gas/rho).
         thermal_energy_per_gram = local_pressure*local_beta/local_density
         sum_thermal_energy = sum_thermal_energy+thermal_energy_per_gram* &
              s%dm(zone_idx)
      end do
! The means: the thermal energy per gram of the convective envelope
! before accretion, the accreted matter's energy per gram, and the
! specific entropy.
      mean_thermal_energy = sum_thermal_energy/cz_total_mass_below_fitting
      accretion_specific_energy = thermal_energy_accreted_bar/ &
           cz_total_mass_below_fitting
      s%envelope_specific_entropy = s%envelope_specific_entropy/ &
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
      log10_gravity = cgl+s%log_total_mass-2.0d0*log10_radius
      call atm_get_surface_pt(4, s%log_Teff, log10_gravity, &
           atm_point%log10_pressure, atm_point%log10_temperature, &
           allard_surface_failed, ierr)
      if (ierr /= 0) return
      log10_temperature_local = atm_point%log10_temperature
      log10_pressure_local = atm_point%log10_pressure
      hydrogen_fraction_local = s%xa(1,nz)
      metal_fraction_local = s%xa(3,nz)
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
         mass_loss_rate_cgs = mass_loss_rate_msun_yr*s%solar_mass_cgs/ &
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
! convection zone's masses. The arithmetic is in the star-blind kernels
! of mass_change_kernel; this driver gathers the zone from star, calls
! them and stores the results. Inputs: total_radius_cm (R), the rate
! mass_loss_rate_msun_yr (Msun/yr, positive for accretion),
! accretion_specific_energy and mean_thermal_energy (the energies per
! gram of the accreted matter and of the convection zone), and
! cz_total_mass_below_fitting (the convective envelope's mass).
! Outputs: old_log_envelope_mass_fraction = log10(m(nz)/M) before the
! change; disk_exhausted, true when the disk-locking age was reached
! (the timestep is clipped to land on it); ierr = 1 when the envelope
! boundary is the surface zone.
!   star reads: ctrl%accreted_composition, ctrl%fczdmdt, ctrl%ftotdmdt,
!     ctrl%walpcz, disk_gate_age_gyr, envelope_cz_bottom_index,
!     envelope_specific_entropy, job%disk_locking_active,
!     job%disk_locking_age_gyr, job%rotation_active, logRho, nz, omega,
!     solar_mass_cgs
!   star writes: accreted_mass, delta_log_pressure,
!     delta_log_temperature, dm, dt, j_rot, log_mass, log_total_mass,
!     logP, logR, logT, m, recompute_envelope_triangle, senv, star_mass,
!     stotal, xa
subroutine apply_mass_change(s, total_radius_cm, mass_loss_rate_msun_yr, &
     accretion_specific_energy, mean_thermal_energy, &
     cz_total_mass_below_fitting, old_log_envelope_mass_fraction, &
     disk_exhausted, ierr)
      type(star_info), intent(inout) :: s
      double precision, intent(in) :: total_radius_cm
      double precision, intent(inout) :: mass_loss_rate_msun_yr
      double precision, intent(in) :: accretion_specific_energy, &
           mean_thermal_energy
      double precision, intent(in) :: cz_total_mass_below_fitting
      double precision, intent(out) :: old_log_envelope_mass_fraction
      logical, intent(out) :: disk_exhausted
      integer, intent(out) :: ierr
      integer :: nz, jenv, zone_idx
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
      logical :: wind_carries_angular_momentum
      double precision, allocatable :: cz_shell_masses_before(:)
      double precision :: wind_j_surface
      type(wind_angular_momentum) :: wind_am
      double precision :: delta_mass_msun
      double precision :: total_mass_grams_old, total_mass_grams_new

      ierr = 0
      disk_exhausted = .false.
      nz = s%nz
      jenv = s%envelope_cz_bottom_index
      old_log_envelope_mass_fraction = s%log_mass(nz) - s%log_total_mass
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
      mass_loss_rate_cgs = abs(mass_loss_rate_msun_yr)*s%solar_mass_cgs/ &
           seconds_per_year
      call weigh_surface_convection_zone()
      if (ierr /= 0) return
      call limit_mass_change_timestep()
! The mass change of the step, dM = Mdot dt, in g (negative for loss).
      delta_mass_cgs = mass_loss_rate_msun_yr*s%solar_mass_cgs*s%dt/ &
           seconds_per_year
! With mass loss from a rotating star, the lost matter carries off
! angular momentum (a thermal wind).
      wind_carries_angular_momentum = s%job%rotation_active .and. &
           delta_mass_cgs.lt.0.0d0
      call adjust_cz_for_accreted_entropy()
      call update_masses()
      call spin_down_cz_by_wind()
      call report_new_mass()
      call mix_in_accreted_composition()

contains

! ---------------------------------------------------------------
! The mass of the surface convection zone below the fitting point
! (convection_zone_mass); ierr = 1 when the envelope boundary is the
! surface zone.
subroutine weigh_surface_convection_zone()
      if(jenv.eq.nz) then
         write(*,*) 'mdot: envelope boundary at the surface (911)'
         ierr = 1
         return
      end if
      cz_mass_grams = convection_zone_mass(jenv, s%m(1:nz))
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
         timestep_limit_total_mass = s%ctrl%ftotdmdt*s%star_mass*s%solar_mass_cgs/ &
              mass_loss_rate_cgs
         timestep_limit_cz_mass = s%ctrl%fczdmdt*cz_mass_grams/mass_loss_rate_cgs
         timestep_limit = min(timestep_limit_total_mass,timestep_limit_cz_mass)
! Nor by more than half the mass outside the fitting point (the outer
! zone), (M - m(nz))/2, for loss and accretion alike.
         timestep_limit_envelope = 0.5d0*(exp10(s%log_total_mass)- &
              exp10(s%log_mass(nz)))/mass_loss_rate_cgs
         timestep_limit = min(timestep_limit,timestep_limit_envelope)
         if(timestep_limit.lt.s%dt)then
            write(*,10)s%dt,timestep_limit
 10         format(' TIMESTEP REDUCED FOR MASS LOSS - OLD,NEW',1P2E12.3)
            s%dt = timestep_limit
         endif
      endif
! With accretion only (Mdot > 0, MHP 8/10): a step that would pass the
! disk-locking age is clipped to end there, and the mass change stops.
      if(s%job%disk_locking_active .and. mass_loss_rate_msun_yr.gt.0.0d0)then
         disk_age_test = s%disk_gate_age_gyr + 1.0d-9*s%dt/seconds_per_year
         if(disk_age_test.gt.s%job%disk_locking_age_gyr)then
            s%dt = (s%job%disk_locking_age_gyr-s%disk_gate_age_gyr)*1.0d9*seconds_per_year
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
         local_temperature = exp10(s%logT(jenv))
         local_pressure = exp10(s%logP(jenv))
         local_density = exp10(s%logRho(jenv))
         local_beta = 1.0d0-(radiation_constant_over_3*local_temperature**4/ &
              local_pressure)
         mean_molecular_weight_local = local_pressure*local_beta/ &
              (local_density*local_temperature)
         delta_log_specific_entropy = (delta_mass_cgs/cz_mass_grams)* &
              s%envelope_specific_entropy/mean_molecular_weight_local
! The mass term d ln M = ln((M_cz + dM)/M_cz) is deliberately left out.
         delta_ln_mass = 0.0d0
         delta_log_radius = (cc23*delta_log_specific_entropy- &
              cc13*delta_ln_mass)/ln10
         s%delta_log_temperature = delta_ln_mass/ln10 - delta_log_radius
         s%delta_log_pressure = 2.0d0*delta_ln_mass/ln10 - 4.0d0*delta_log_radius
! A fully convective star (jenv = 1) shifts every zone. Otherwise only
! the base zone's P and T shift, and the radii above it stretch about
! the base, r' = r_b + 10^(d log R) (r - r_b).
         if(jenv.eq.1)then
            do zone_idx = jenv,nz
               s%logR(zone_idx) = s%logR(zone_idx)+delta_log_radius
               s%logP(zone_idx) = s%logP(zone_idx)+s%delta_log_pressure
               s%logT(zone_idx) = s%logT(zone_idx)+s%delta_log_temperature
            end do
         else
            s%logP(jenv) = s%logP(jenv)+s%delta_log_pressure
            s%logT(jenv) = s%logT(jenv)+s%delta_log_temperature
            boundary_radius_cm = exp10(s%logR(jenv))
            radius_scale_factor = exp10(delta_log_radius)
            do zone_idx = jenv+1,nz
               radius_before_cm = exp10(s%logR(zone_idx))
               radius_after_cm = boundary_radius_cm+radius_scale_factor* &
                    (radius_before_cm-boundary_radius_cm)
               s%logR(zone_idx) = log10(radius_after_cm)
            end do
         endif
      endif
end subroutine adjust_cz_for_accreted_entropy

! ---------------------------------------------------------------
! The new total mass, the masses of the convection zone
! (change_shell_masses: with mass loss the zone-center masses compress
! about the zone's base; the shell masses dm, the last reaching the new
! total) and senv; the surface boundary condition is to be recomputed.
! With a wind that carries angular momentum, the zone's shell masses
! before the change are kept for spin_down_cz_by_wind.
subroutine update_masses()
      if(wind_carries_angular_momentum) &
           cz_shell_masses_before = s%dm(jenv:nz)
! The new total mass, in Msun (star_mass) and as log10 M in g
! (log_total_mass and its copy stotal).
      s%star_mass = s%star_mass + delta_mass_cgs/s%solar_mass_cgs
      total_mass_grams_old = exp10(s%log_total_mass)
      total_mass_grams_new = total_mass_grams_old + delta_mass_cgs
      s%log_total_mass = log10(total_mass_grams_new)
      s%stotal = s%log_total_mass
      call change_shell_masses(jenv, delta_mass_cgs, total_mass_grams_new, &
           s%m(1:nz), s%log_mass(1:nz), s%dm(1:nz))
! senv = log10(m(nz)/M), the mass fraction inside the fitting point.
      s%senv = s%log_mass(nz) - s%log_total_mass
! The envelope triangle (the surface boundary condition) is to be
! recomputed.
      s%recompute_envelope_triangle = .true.
end subroutine update_masses

! ---------------------------------------------------------------
! With mass loss from a rotating star: the wind carries off the
! surface's specific angular momentum under the convection zone's
! rotation law (wind_specific_angular_momentum), and the zone keeps
! J_cz + dJ on its new masses (remove_wind_angular_momentum).
subroutine spin_down_cz_by_wind()
      if(wind_carries_angular_momentum)then
         wind_j_surface = wind_specific_angular_momentum(s%ctrl%walpcz, &
              total_radius_cm, s%omega(nz), s%logR(nz), s%j_rot(nz))
         call remove_wind_angular_momentum(wind_j_surface, delta_mass_cgs, &
              cz_shell_masses_before, s%dm(jenv:nz), s%j_rot(jenv:nz), &
              wind_am)
! The run log: dM, dJ, omega(nz), the new and old J_cz, and the
! convective envelope's mass.
         write(*,22)delta_mass_cgs,wind_am%lost,s%omega(nz), &
              wind_am%zone_after,wind_am%zone_before,cz_total_mass_below_fitting
 22      format(1P6E13.4)
      endif
end subroutine spin_down_cz_by_wind

! ---------------------------------------------------------------
! The run log: the new mass and the step's change, in Msun.
subroutine report_new_mass()
      delta_mass_msun = delta_mass_cgs/s%solar_mass_cgs
      write(*,20)s%star_mass,delta_mass_msun
 20   format('MASS LOSS APPLIED - NEW M,DEL M',1P2E19.10)
end subroutine report_new_mass

! ---------------------------------------------------------------
! With accretion: the accreted composition mixed into the convection
! zone (mix_accreted_composition, with M_cz = cz_total_mass_below_fitting
! and X_acc = accreted_composition), and the accreted mass recorded.
subroutine mix_in_accreted_composition()
      if(delta_mass_cgs.gt.0.0d0)then
         call mix_accreted_composition(cz_total_mass_below_fitting, &
              delta_mass_cgs, s%ctrl%accreted_composition, &
              s%xa(:,jenv:nz))
         s%accreted_mass = delta_mass_cgs
      endif
end subroutine mix_in_accreted_composition

end subroutine apply_mass_change

end module mass_change

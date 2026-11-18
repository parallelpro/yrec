!----------------------------------------------------------------------
! observables
!----------------------------------------------------------------------
! The per-model observables of the converged star, one subroutine per
! theme. compute_observables is called
! from evolve_step once per model, after the stitched model is built
! and before the output writer. Every quantity it computes is stored on
! star% (the observables members and the renormalized
! star%luminosity_breakdown_lsun), so the writer only reads. One of
! them is model state: renormalize_luminosity_breakdown rescales the
! luminosity channels, which the model file stores and the next card's
! first compute_timestep reads (its header). The turnover time, which
! the wind and the deuterium cap read, is refreshed by the step
! (step_driver's refresh_turnover_timescale, after compute_observables).
!
! The order of the themes matters: the central conditions precede the
! base of the surface convection zone (whose fully convective branch
! reads star%central_log10_*), and compute_surface_globals precedes the
! themes that read star%log_R_surface_rsun.
!
! The star carries the central eos call's Saha iteration counter
! (s%central_saha_state; eos_get adds to it and nothing reads it).
module observables
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: i_h1, i_metals, i_nu_f17, i_nu_pp, n_lum_channels, &
           n_nu_fluxes, star_info
      use energy_sources, only: num_energy_sources
      use structure_equations, only: luminosity_channel_sum
      use eos_lib, only: eos_get, eos_result
      use physical_constants, only: c4pi, c4pi3l, cc13, cc23, cgl, cpi, ln10, &
           log10_radius_from_l_teff
      implicit none
      private
      public :: compute_observables

! The eight solar-neutrino fluxes a detector captures (pp, pep, hep, 7Be,
! 8B, 13N, 15O, 17F: neutrino_flux_layout's slots i_nu_pp to i_nu_f17;
! the two slots after them are reaction rates, not neutrino sources).
      integer, parameter :: n_detected_fluxes = i_nu_f17 - i_nu_pp + 1
! The neutrino absorption cross sections of the gallium (71Ga) and
! chlorine (37Cl) detectors averaged over each flux's spectrum, in units
! of 1e-46 cm^2, in that order. With the fluxes in 1e10 cm^-2 s^-1 the
! detector's capture rate sum_i sigma_i phi_i is in SNU (1e-36 captures
! per target atom per second). The values are the F77 WRTOUT's (Sources).
      real(dp), parameter :: gallium_capture_cross_section(n_detected_fluxes) = &
           [1.18D1,2.15D2,7.14D4,7.17D1,2.40D4,6.04D1,1.137D2,1.139D2]
      real(dp), parameter :: chlorine_capture_cross_section(n_detected_fluxes) = &
           [0.0D0,1.6D1,4.26D4,2.4D0,1.14D4,1.7D0,6.8D0,6.9D0]

! Unit conversions and limits of the observables: Hz to uHz; seconds in
! a day; cm/s to km/s; the longest rotation period reported [days], which
! stands for a surface at rest.
      real(dp), parameter :: uhz_per_hz = 1.0d6, seconds_per_day = 8.64d4, &
           km_per_cm = 1.0d-5, max_rotation_period_days = 9999.0d0

contains

! ---------------------------------------------------------------
! Fills every per-model observable on star%, theme by theme, in the
! order the header explains.
!   star reads: none
!   star writes: none
subroutine compute_observables(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

      ierr = 0

      call renormalize_luminosity_breakdown(s)
      call locate_core_cz(s)
      call compute_central_conditions(s, ierr)
      if (ierr /= 0) return
      call locate_surface_cz_base(s)

! The per-model history columns (star% members read by write_history),
! by the formulas of the legacy .track output. The surface radius also
! feeds the run log's model line.
      call compute_surface_globals(s)
      call compute_seismic_observables(s)
      call compute_moment_of_inertia(s)
      call compute_snu_rates(s)
      call compute_rotation_observables(s)
      call compute_h_shell_boundaries(s)

end subroutine compute_observables

! ---------------------------------------------------------------
! Rescales the luminosity channels star%luminosity_breakdown_lsun so that
! they sum to the surface luminosity: each channel times
! L(nz)/sum of channels (structure_equations' luminosity_channel_sum,
! with the extra energy sources' channel when there are any). Skipped
! during the helium flash. This changes the model, after the step's
! compute_timestep has read the solver's channels: the rescaled values
! are what the model file stores and what the next card's first
! compute_timestep reads (the hydrogen and helium burning limits and the
! sign of L_grav). The factor is 1 to within the solver's last
! correction.
!   star reads: luminosity_lsun, nz, physics%helium_flash_active
!   star writes: luminosity_breakdown_lsun
subroutine renormalize_luminosity_breakdown(s)
      type(star_info), intent(inout) :: s
      real(dp) :: total_luminosity_sum, temp_value
      integer :: i

! The post-convergence twin of henyey_iterate's per-iteration
! renormalization, which also guards against a zero sum; kept separate
! on purpose (audit/readability-sweep-2026-09-03/reports/henyey-wave2.md,
! item 6).
      if(.not.s%physics%helium_flash_active) then
       total_luminosity_sum = luminosity_channel_sum(s%luminosity_breakdown_lsun, &
            num_energy_sources() > 0)
       temp_value = s%luminosity_lsun(s%nz)/total_luminosity_sum
       do i = 1,n_lum_channels
          s%luminosity_breakdown_lsun(i) = s%luminosity_breakdown_lsun(i)*temp_value
       end do
      end if
end subroutine renormalize_luminosity_breakdown

! ---------------------------------------------------------------
! The mass of the convective core in Msun, m at its top zone (zero
! without a convective core).
!   star reads: core_cz_top_index, m, solar_mass_cgs
!   star writes: core_cz_mass_msun
subroutine locate_core_cz(s)
      type(star_info), intent(inout) :: s
      if(s%core_cz_top_index > 1) then
       s%core_cz_mass_msun = s%m(s%core_cz_top_index)/s%solar_mass_cgs
      else
       s%core_cz_mass_msun = 0.0D0
      end if
end subroutine locate_core_cz

! ---------------------------------------------------------------
! The central P and T, extrapolated from the first zone, and the
! central density, beta and electron degeneracy parameter eta from the
! eos there. Stores star%central_*.
!   star reads: gradT, log_mass, logP, logRho, logT, xa
!   star writes: central_beta, central_degeneracy_eta,
!     central_log10_density, central_log10_pressure,
!     central_log10_temperature, central_saha_state
subroutine compute_central_conditions(s, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

      real(dp) :: temp_value
      real(dp) :: pressure_linear, log_pressure_center, &
           log_temperature_center, hydrogen_fraction_center, &
           metal_fraction_center
      logical :: is_atmosphere_point, compute_derivatives
! The eos results, as one eos_result record.
      type(eos_result) :: eos_res

      ierr = 0
      eos_res = eos_result()

! Hydrostatic equilibrium in a uniform-density center, with
! r_1 = (3 m_1/(4 pi rho_1))^(1/3), gives
!   P_c = P_1 + G m_1 rho_1/(2 r_1).
      temp_value =0.5D0*exp(ln10*(cc13*(c4pi3l+s%logRho(1)-s%log_mass(1))+s%logRho(1)+cgl+s%log_mass(1)))
      pressure_linear = exp(ln10*s%logP(1))
      log_pressure_center = log10(pressure_linear + temp_value)
! With the actual temperature gradient del = d ln T/d ln P of zone 1
! (star%gradT(1)), T_c = T_1 (1 + del (P_c - P_1)/P_1).
      log_temperature_center = s%logT(1) + log10(1.0D0+ temp_value*s%gradT(1)/pressure_linear)
      eos_res%log10_density = s%logRho(1)
      hydrogen_fraction_center = s%xa(i_h1,1)
      metal_fraction_center = s%xa(i_metals,1)
      is_atmosphere_point = .true.
      compute_derivatives = .false.
! The eos at (P_c, T_c) for the central density, beta and eta.
      call eos_get(log_temperature_center, log_pressure_center, &
           hydrogen_fraction_center, metal_fraction_center, eos_res, &
           compute_derivatives, is_atmosphere_point, s%central_saha_state, &
           composition_at_zone=s%xa(:,1), ierr=ierr)
      if (ierr /= 0) return
! Store the central values for later use.
      s%central_log10_pressure = log_pressure_center
      s%central_log10_temperature = log_temperature_center
      s%central_log10_density = eos_res%log10_density
      s%central_beta = eos_res%beta
      s%central_degeneracy_eta = eos_res%electron_degeneracy_parameter
end subroutine compute_central_conditions

! ---------------------------------------------------------------
! The base of the surface convection zone, interpolated between the
! zones around it: star%envelope_mass_msun (the convective envelope's
! mass), envelope_radius_rsun, and
! envelope_cz_temperature, _density, _pressure and _opacity. A fully
! convective star takes the central values (so compute_central_conditions
! runs first); without a surface convection zone all are zero. The
! step's turnover refresh, after compute_observables, replaces
! envelope_radius_rsun with the turnover walk's own base radius, the one
! the output writes.
!   star reads: central_log10_density, central_log10_pressure,
!     central_log10_temperature, envelope_cz_bottom_index, grada, gradr,
!     log10_solar_radius, log_total_mass, logP, logR, logRho, logT, m,
!     nz, opacity_zone, solar_mass_cgs, star_mass
!   star writes: envelope_cz_density, envelope_cz_opacity,
!     envelope_cz_pressure, envelope_cz_temperature, envelope_mass_msun,
!     envelope_radius_rsun
subroutine locate_surface_cz_base(s)
      use math_lib
      use convection_geometry, only: cz_base_fraction
      type(star_info), intent(inout) :: s
      real(dp) :: dd1, dd2, cz_base_mass, envelope_boundary_fx
      real(dp) :: envelope_cz_log_temperature, &
           envelope_cz_log_density, envelope_cz_log_pressure, &
           envelope_cz_log_radius

      if(s%envelope_cz_bottom_index < s%nz) then
       if(s%envelope_cz_bottom_index > 1) then
! The fraction f of the way from zone b-1 to the base zone b where
! gradr - grada crosses zero (cz_base_fraction), and the mass there.
! Linear in f: m, log r, the opacity, log T, log rho and log P.
            dd2 = s%gradr(s%envelope_cz_bottom_index-1)-s%grada(s%envelope_cz_bottom_index-1)
            dd1 = s%gradr(s%envelope_cz_bottom_index)-s%grada(s%envelope_cz_bottom_index)
            envelope_boundary_fx = cz_base_fraction(dd2, dd1)
            cz_base_mass = s%m(s%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (s%m(s%envelope_cz_bottom_index)-s%m(s%envelope_cz_bottom_index-1))
            s%envelope_mass_msun = (exp(ln10*s%log_total_mass) - cz_base_mass)/s%solar_mass_cgs
! The radius of the base (MHP 2/98).
            envelope_cz_log_radius = s%logR(s%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (s%logR(s%envelope_cz_bottom_index)-s%logR(s%envelope_cz_bottom_index-1))-s%log10_solar_radius
            s%envelope_radius_rsun = exp(ln10*envelope_cz_log_radius)
            s%envelope_cz_opacity = s%opacity_zone(s%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (s%opacity_zone(s%envelope_cz_bottom_index)-s%opacity_zone(s%envelope_cz_bottom_index-1))
            envelope_cz_log_temperature = s%logT(s%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (s%logT(s%envelope_cz_bottom_index)-s%logT(s%envelope_cz_bottom_index-1))
            envelope_cz_log_density = s%logRho(s%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (s%logRho(s%envelope_cz_bottom_index)-s%logRho(s%envelope_cz_bottom_index-1))
            envelope_cz_log_pressure = s%logP(s%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (s%logP(s%envelope_cz_bottom_index)-s%logP(s%envelope_cz_bottom_index-1))
            s%envelope_cz_temperature = exp(ln10*envelope_cz_log_temperature)
            s%envelope_cz_density = exp(ln10*envelope_cz_log_density)
            s%envelope_cz_pressure = exp(ln10*envelope_cz_log_pressure)
       else
          s%envelope_mass_msun = s%star_mass
          s%envelope_radius_rsun = 0.0D0
            s%envelope_cz_temperature = exp10(s%central_log10_temperature)
            s%envelope_cz_density = exp10(s%central_log10_density)
            s%envelope_cz_pressure = exp10(s%central_log10_pressure)
            s%envelope_cz_opacity = s%opacity_zone(1)
       end if
      else
       s%envelope_mass_msun = 0.0D0
       s%envelope_radius_rsun = 0.0D0
         s%envelope_cz_temperature = 0.0D0
         s%envelope_cz_density = 0.0D0
         s%envelope_cz_pressure = 0.0D0
         s%envelope_cz_opacity = 0.0D0
      end if
end subroutine locate_surface_cz_base

! ---------------------------------------------------------------
! The surface radius from L = 4 pi R^2 sigma Teff^4, log g = log G +
! log M - 2 log R (CGS), then log_R_surface_rsun = log10(R/Rsun).
!   star reads: log10_solar_luminosity, log10_solar_radius, log_L_lsun,
!     log_Teff, log_total_mass
!   star writes: log_g_surface, log_R_surface_rsun
subroutine compute_surface_globals(s)
      type(star_info), intent(inout) :: s
      s%log_R_surface_rsun = log10_radius_from_l_teff(s%log_L_lsun, s%log_Teff, &
           s%log10_solar_luminosity)
      s%log_g_surface = cgl + s%log_total_mass &
           - 2.0d0*s%log_R_surface_rsun
      s%log_R_surface_rsun = s%log_R_surface_rsun - s%log10_solar_radius
end subroutine compute_surface_globals

! ---------------------------------------------------------------
! The asteroseismic observables. The scaling relations, against the
! solar references nu_max_sun, delta_nu_sun and Teff_sun of the
! controls,
!   nu_max_uhz = nu_max_sun (g/g_sun) (Teff/Teff_sun)^(-1/2),
!   delta_nu_rho_uhz = delta_nu_sun (M/Msun)^(1/2) (R/Rsun)^(-3/2),
! with log g_sun from the run's own solar constants, so that
! Monte-Carlo-scaled solar values stay self-consistent; then the
! asymptotic delta_nu_uhz and delta_Pg from the stitched model.
!   star reads: log_g_surface, log_R_surface_rsun, log_Teff,
!     physics%delta_nu_sun, physics%nu_max_sun, physics%Teff_sun,
!     solar_mass_cgs, solar_radius_cgs, star_mass, stitched
!   star writes: delta_nu_rho_uhz, delta_nu_uhz, delta_Pg, nu_max_uhz
subroutine compute_seismic_observables(s)
      use math_lib
      use stitched_model, only: ip_logR, &
           ip_brunt_N2, ip_csound
      use stitched_kernels, only: acoustic_radius, buoyancy_radius
      type(star_info), intent(inout) :: s
      real(dp) :: log_g_solar
      real(dp) :: sound_travel_time, buoyancy_integral
      integer :: n_ie

      log_g_solar = cgl + log10(s%solar_mass_cgs) &
           - 2.0d0*log10(s%solar_radius_cgs)
      s%nu_max_uhz = s%physics%nu_max_sun &
           * exp10((s%log_g_surface - log_g_solar)) &
           * sqrt(s%physics%Teff_sun/exp10(s%log_Teff))
! The mean-density scaling.
      s%delta_nu_rho_uhz = s%physics%delta_nu_sun &
           * sqrt(s%star_mass) * exp10((-1.5d0*s%log_R_surface_rsun))

! The asymptotic p-mode large separation from the sound travel time,
!   delta_nu = [2 int dr/c]^(-1),
! over the stitched model's points 1..n_ie (the innermost point to the
! photosphere; the atmosphere above R is excluded, as in MESA's acoustic
! radius; stitched_kernels' acoustic_radius). In uHz (the integral is in
! seconds).
      n_ie = s%stitched%n_ie
      sound_travel_time = acoustic_radius(s%stitched%prof(ip_logR,1:n_ie), &
           s%stitched%prof(ip_csound,1:n_ie))
      if (sound_travel_time > 0.0d0) then
         s%delta_nu_uhz = uhz_per_hz/(2.0d0*sound_travel_time)
      else
         s%delta_nu_uhz = 0.0d0
      end if

! The asymptotic l = 1 g-mode period spacing from the Brunt-Vaisala
! frequency,
!   delta_Pg = 2 pi^2/sqrt(l(l+1)) / int N/r dr = sqrt(2) pi^2 / int N/r dr,
! over the innermost g-mode cavity only (stitched_kernels'
! buoyancy_radius). Zero when there is no cavity (a fully convective
! interior). In seconds.
      buoyancy_integral = buoyancy_radius(s%stitched%prof(ip_logR,1:n_ie), &
           s%stitched%prof(ip_brunt_N2,1:n_ie))
      if (buoyancy_integral > 0.0d0) then
         s%delta_Pg = cpi*cpi*sqrt(2.0d0)/buoyancy_integral
      else
         s%delta_Pg = 0.0d0
      end if
end subroutine compute_seismic_observables

! ---------------------------------------------------------------
! The total moment of inertia: without rotation the thin-shell sum
! I = sum (2/3) dm r^2, with rotation the sum of the shells' i_rot.
!   star reads: dm, i_rot, logR, nz, physics%rotation_active
!   star writes: total_moment_of_inertia
subroutine compute_moment_of_inertia(s)
      use math_lib
      type(star_info), intent(inout) :: s
      integer :: i

      s%total_moment_of_inertia = 0.0d0
      if (.not. s%physics%rotation_active) then
         do i = 1, s%nz
            s%total_moment_of_inertia = &
                 s%total_moment_of_inertia + &
                 cc23*s%dm(i)*exp(2.0d0*ln10*s%logR(i))
         end do
      else
         do i = 1, s%nz
            s%total_moment_of_inertia = &
                 s%total_moment_of_inertia + s%i_rot(i)
         end do
      end if
end subroutine compute_moment_of_inertia

! ---------------------------------------------------------------
! The chlorine and gallium capture rates (SNU) from the neutrino flux
! totals; without calc_neutrinos the rates and all ten fluxes are zero.
!   star reads: physics%calc_neutrinos
!   star writes: cl37_snu_rate, ga71_snu_rate,
!     neutrino_flux_total_e10_per_cm2_s
subroutine compute_snu_rates(s)
      type(star_info), intent(inout) :: s
      integer :: i

      s%cl37_snu_rate = 0.0d0
      s%ga71_snu_rate = 0.0d0
      if (s%physics%calc_neutrinos) then
         do i = 1, n_detected_fluxes
            s%cl37_snu_rate = s%cl37_snu_rate + &
                 chlorine_capture_cross_section(i)*s%neutrino_flux_total_e10_per_cm2_s(i)
            s%ga71_snu_rate = s%ga71_snu_rate + &
                 gallium_capture_cross_section(i)*s%neutrino_flux_total_e10_per_cm2_s(i)
         end do
      else
         do i = 1, n_nu_fluxes
            s%neutrino_flux_total_e10_per_cm2_s(i) = 0.0d0
         end do
      end if
end subroutine compute_snu_rates

! ---------------------------------------------------------------
! The surface rotation period, P = 2 pi/omega(nz) in days (at most
! max_rotation_period_days), the surface velocity omega(nz) R in km/s, and the moment of
! inertia of the surface convection zone (the i_rot sum with rotation,
! sum (2/3) dm r^2 without). Reads star%log_R_surface_rsun, so
! compute_surface_globals runs first.
!   star reads: dm, envelope_cz_bottom_index, i_rot, log10_solar_radius,
!     log_R_surface_rsun, logR, nz, omega, physics%rotation_active
!   star writes: cz_moment_of_inertia, rotation_period_days,
!     surf_velocity_kms
subroutine compute_rotation_observables(s)
      use math_lib
      type(star_info), intent(inout) :: s
      integer :: i

      s%cz_moment_of_inertia = 0.0d0
      if (s%physics%rotation_active) then
         s%rotation_period_days = min(max_rotation_period_days, &
              0.5d0*c4pi/s%omega(s%nz)/seconds_per_day)
         s%surf_velocity_kms = s%omega(s%nz)* &
              exp(ln10*(s%log_R_surface_rsun+s%log10_solar_radius))*km_per_cm
         if (s%envelope_cz_bottom_index < s%nz) then
            do i = s%envelope_cz_bottom_index, s%nz
               s%cz_moment_of_inertia = &
                    s%cz_moment_of_inertia + s%i_rot(i)
            end do
         end if
      else
         s%rotation_period_days = 0.0d0
         s%surf_velocity_kms = 0.0d0
         if (s%envelope_cz_bottom_index < s%nz) then
            do i = s%envelope_cz_bottom_index, s%nz
               s%cz_moment_of_inertia = &
                    s%cz_moment_of_inertia + &
                    cc23*s%dm(i)*exp(2.0d0*ln10*s%logR(i))
            end do
         end if
      end if
end subroutine compute_rotation_observables

! ---------------------------------------------------------------
! The bottom, middle and top of the hydrogen-burning shell: masses in
! Msun and radii as r/R (zero without a shell). Reads
! star%log_R_surface_rsun, so compute_surface_globals runs first.
!   star reads: h_shell_bot_index, h_shell_mid_index, h_shell_top_index,
!     has_h_shell, log10_solar_radius, log_R_surface_rsun, logR, m,
!     solar_mass_cgs
!   star writes: h_shell_bot_mass_msun, h_shell_bot_radius_fraction,
!     h_shell_mid_mass_msun, h_shell_mid_radius_fraction,
!     h_shell_top_mass_msun, h_shell_top_radius_fraction
subroutine compute_h_shell_boundaries(s)
      use math_lib
      type(star_info), intent(inout) :: s
      if (s%has_h_shell) then
         s%h_shell_bot_mass_msun = &
              s%m(s%h_shell_bot_index)/s%solar_mass_cgs
         s%h_shell_mid_mass_msun = &
              s%m(s%h_shell_mid_index)/s%solar_mass_cgs
         s%h_shell_top_mass_msun = &
              s%m(s%h_shell_top_index)/s%solar_mass_cgs
         s%h_shell_bot_radius_fraction = exp(ln10*( &
              s%logR(s%h_shell_bot_index) &
              - s%log_R_surface_rsun - s%log10_solar_radius))
         s%h_shell_mid_radius_fraction = exp(ln10*( &
              s%logR(s%h_shell_mid_index) &
              - s%log_R_surface_rsun - s%log10_solar_radius))
         s%h_shell_top_radius_fraction = exp(ln10*( &
              s%logR(s%h_shell_top_index) &
              - s%log_R_surface_rsun - s%log10_solar_radius))
      else
         s%h_shell_bot_mass_msun = 0.0d0
         s%h_shell_mid_mass_msun = 0.0d0
         s%h_shell_top_mass_msun = 0.0d0
         s%h_shell_bot_radius_fraction = 0.0d0
         s%h_shell_mid_radius_fraction = 0.0d0
         s%h_shell_top_radius_fraction = 0.0d0
      end if
end subroutine compute_h_shell_boundaries

! ---------------------------------------------------------------
! Sources
! 1. The themes were blocks of the F77 writer WRTOUT
!    (git show 6cd56739:src/wrtout.f), then
!    core/update_output_diagnostics.f90. Of the F77 routine's saved
!    locals only the central eos call's Saha iteration counter outlived a
!    call. WRTOUT's core-boundary radius, which read the previous call's
!    interpolation fraction through an FX/FX2 typo and was never used, is
!    not carried over.
! 2. The capture cross sections are WRTOUT's DATA statements (GASNUF and
!    CLSNUF), which name no source. Gallium: the hep, 7Be, 8B, 13N, 15O
!    and 17F values are those of Bahcall 1997, Phys. Rev. C 56, 3391
!    (7.14e4, 71.7, 2.40e4, 60.4, 113.7 and 113.9); that paper's pp and
!    pep values are 11.72 and 204, where the table has 11.8 and 215.
!    Chlorine: the 8B value is that of Bahcall et al. 1996, Phys. Rev. C
!    54, 411 (1.14e4); the others were not traced to a paper.
end module observables

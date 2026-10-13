!----------------------------------------------------------------------
! observables
!----------------------------------------------------------------------
! The per-model observables of the converged star, one subroutine per
! theme (formerly core/update_output_diagnostics.f90, and before that
! blocks of the legacy writer wrtout). compute_observables is called
! from evolve_step once per model, after the stitched model is built
! and before the output writer. Every quantity it computes is stored on
! star% (the observables members and the renormalized
! star%luminosity_breakdown), so the writer only reads; its one physics
! feedback is the turnover timescale (see compute_observables). The
! module also holds shell_masses_from_log_mass, the m and dm of a
! log-mass grid.
!
! The order of the themes matters: the central conditions precede the
! base of the surface convection zone (whose fully convective branch
! reads star%central_log10_*), and compute_surface_globals precedes the
! themes that read star%log_R_surface.
!
! One local of the legacy routine's blanket SAVE is real cross-call
! state: the Saha-table state of the central eos call, which the star
! carries (s%central_saha_state). The others are set before every use
! and are plain locals. (Legacy wrtout's core-boundary radius, which read the
! previous call's interpolation fraction through an FX/FX2 typo and was
! never used, is not carried over.)
module observables
      use turnover_timescale, only: compute_turnover_timescale
      use star_info_lib, only: i_h1, i_metals, i_lum_pp1, i_lum_pp2, i_lum_pp3, &
           i_lum_cno, i_lum_3alpha, i_lum_neu, i_lum_grav, i_lum_he_c, &
           n_lum_channels, star_info
      use eos_lib
      use physical_constants
      implicit none
      private
      public :: compute_observables
      public :: shell_masses_from_log_mass

! The capture rate per unit flux of the eight solar-neutrino fluxes (pp,
! pep, hep, 7Be, 8B, 13N, 15O, 17F, the order of neutrino_flux_layout)
! for the gallium (71Ga) and chlorine (37Cl) detectors: the detector's
! rate in SNU is sum_i c_i phi_i.
      double precision, parameter :: gasnuf_diag(8) = [1.18D1,2.15D2, &
           7.14D4,7.17D1,2.40D4,6.04D1,1.137D2,1.139D2]
      double precision, parameter :: clsnuf_diag(8) = [0.0D0,1.6D1, &
           4.26D4,2.4D0,1.14D4,1.7D0,6.8D0,6.9D0]

contains

! ---------------------------------------------------------------
! Fills every per-model observable on star%, theme by theme, in the
! order the header explains.
subroutine compute_observables(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

      ierr = 0

      call renormalize_luminosity_breakdown(s)
      call locate_core_cz(s)
      call compute_central_conditions(s, ierr)
      if (ierr /= 0) return
      call locate_surface_cz_base(s)
! The turnover timescale (with its _old and pphot0 copies) is physics
! state, read by the wind saturation and the deuterium limiter, and is
! refreshed here for every model.
      call refresh_turnover_timescale(s, ierr)
      if (ierr /= 0) return

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
! The shell-center masses m(1:nz) (g) and the shell masses dm(1:nz)
! of the log10 mass grid log_mass: dm(i) = (m(i+1) - m(i-1))/2, with a
! mirror point m(0) = -m(1) below the center, and the last shell reaching
! the total mass, dm(nz) = M - (m(nz-1) + m(nz))/2. Shared by
! read_starting_model (build_shell_masses: star%m, star%dm) and
! write_neutrino_flux_table (its local copies).
subroutine shell_masses_from_log_mass(log_mass, log_total_mass, nz, m, dm)
      use math_lib
      double precision, intent(in) :: log_mass(:), log_total_mass
      integer, intent(in) :: nz
      double precision, intent(inout) :: m(:), dm(:)
      double precision :: prev_mass, curr_mass, next_mass
      integer :: i

      next_mass = exp(ln10*log_mass(1))
      curr_mass = - next_mass
      do i = 2,nz
       prev_mass = curr_mass
       curr_mass = next_mass
       next_mass = exp(ln10*log_mass(i))
       m(i-1) = curr_mass
       dm(i-1) = 0.5d0*(next_mass-prev_mass)
      end do
      m(nz) = next_mass
      dm(nz) = exp(ln10*log_total_mass) - 0.5d0*(curr_mass+ &
           next_mass)
end subroutine shell_masses_from_log_mass

! ---------------------------------------------------------------
! Rescales the luminosity channels star%luminosity_breakdown so that
! they sum to the surface luminosity: each channel times
! L(nz)/sum of channels. Skipped during the helium flash. (This changes
! the model.)
subroutine renormalize_luminosity_breakdown(s)
      type(star_info), intent(inout) :: s
      double precision :: total_luminosity_sum, temp_value
      integer :: i

! The post-convergence twin of henyey_iterate's per-iteration
! renormalization, which also guards against a zero sum; kept separate
! on purpose (audit/readability-sweep-2026-09-03/reports/henyey-wave2.md,
! item 6).
      if(.not.s%ctrl%helium_flash_active) then
       total_luminosity_sum = s%luminosity_breakdown(i_lum_pp1)+s%luminosity_breakdown(i_lum_pp2)+ &
            s%luminosity_breakdown(i_lum_pp3)+s%luminosity_breakdown(i_lum_cno)+s%luminosity_breakdown(i_lum_3alpha)+ &
            s%luminosity_breakdown(i_lum_neu)+s%luminosity_breakdown(i_lum_grav)+s%luminosity_breakdown(i_lum_he_c)
       temp_value = s%luminosity_lsun(s%nz)/total_luminosity_sum
       do i = 1,n_lum_channels
          s%luminosity_breakdown(i) = s%luminosity_breakdown(i)*temp_value
       end do
      endif
end subroutine renormalize_luminosity_breakdown

! ---------------------------------------------------------------
! The mass of the convective core in Msun, m at its top zone (zero
! without a convective core).
subroutine locate_core_cz(s)
      type(star_info), intent(inout) :: s
      if(s%core_cz_top_index.gt.1) then
       s%core_cz_mass = s%m(s%core_cz_top_index)/s%solar_mass_cgs
      else
       s%core_cz_mass = 0.0D0
      endif
end subroutine locate_core_cz

! ---------------------------------------------------------------
! The central P and T, extrapolated from the first zone, and the
! central density, beta and electron degeneracy parameter eta from the
! eos there. Stores star%central_*.
subroutine compute_central_conditions(s, ierr)
      use math_lib
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

      double precision :: temp_value
      double precision :: pressure_linear, log_pressure_center, &
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
! zones around it: star%envelope_mass (the convective envelope's mass,
! Msun), envelope_radius (Rsun), and
! envelope_cz_temperature, _density, _pressure and _opacity. A fully
! convective star takes the central values (so compute_central_conditions
! runs first); without a surface convection zone all are zero.
subroutine locate_surface_cz_base(s)
      use math_lib
      use convection_geometry, only: cz_base_fraction
      type(star_info), intent(inout) :: s
      double precision :: dd1, dd2, cz_base_mass, envelope_boundary_fx
      double precision :: envelope_cz_log_temperature, &
           envelope_cz_log_density, envelope_cz_log_pressure, &
           envelope_cz_log_radius

      if(s%envelope_cz_bottom_index.lt.s%nz) then
       if(s%envelope_cz_bottom_index.gt.1) then
! The fraction f of the way from zone b-1 to the base zone b where
! gradr - grada crosses zero (cz_base_fraction), and the mass there.
! Linear in f: m, log r, the opacity, log T, log rho and log P.
            dd2 = s%gradr(s%envelope_cz_bottom_index-1)-s%grada(s%envelope_cz_bottom_index-1)
            dd1 = s%gradr(s%envelope_cz_bottom_index)-s%grada(s%envelope_cz_bottom_index)
            envelope_boundary_fx = cz_base_fraction(dd2, dd1)
            cz_base_mass = s%m(s%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (s%m(s%envelope_cz_bottom_index)-s%m(s%envelope_cz_bottom_index-1))
            s%envelope_mass = (exp(ln10*s%log_total_mass) - cz_base_mass)/s%solar_mass_cgs
! The radius of the base (MHP 2/98).
            envelope_cz_log_radius = s%logR(s%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (s%logR(s%envelope_cz_bottom_index)-s%logR(s%envelope_cz_bottom_index-1))-s%log10_solar_radius
            s%envelope_radius = exp(ln10*envelope_cz_log_radius)
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
          s%envelope_mass = s%star_mass
          s%envelope_radius = 0.0D0
            s%envelope_cz_temperature = exp10(s%central_log10_temperature)
            s%envelope_cz_density = exp10(s%central_log10_density)
            s%envelope_cz_pressure = exp10(s%central_log10_pressure)
            s%envelope_cz_opacity = s%opacity_zone(1)
       endif
      else
       s%envelope_mass = 0.0D0
       s%envelope_radius = 0.0D0
         s%envelope_cz_temperature = 0.0D0
         s%envelope_cz_density = 0.0D0
         s%envelope_cz_pressure = 0.0D0
         s%envelope_cz_opacity = 0.0D0
      endif
end subroutine locate_surface_cz_base

! ---------------------------------------------------------------
! Refreshes the turnover timescale and copies it and pphot to their
! _old / pphot0 versions, which the mid-step consumers (the wind, the
! deuterium limiter) read during the next step. compute_turnover_timescale
! also overwrites star%envelope_radius with its own base radius.
subroutine refresh_turnover_timescale(s, ierr)
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr
      call compute_turnover_timescale(s, s%envelope_radius, ierr)
      s%convective_turnover_timescale_old = &
           s%convective_turnover_timescale
      s%pphot0 = s%pphot
end subroutine refresh_turnover_timescale

! ---------------------------------------------------------------
! The surface radius from L = 4 pi R^2 sigma Teff^4, log g = log G +
! log M - 2 log R (CGS), then log_R_surface in Rsun.
subroutine compute_surface_globals(s)
      type(star_info), intent(inout) :: s
      s%log_R_surface = log10_radius_from_l_teff(s%log_L, s%log_Teff, &
           s%log10_solar_luminosity)
      s%log_g_surface = cgl + s%stotal &
           - 2.0d0*s%log_R_surface
      s%log_R_surface = s%log_R_surface - s%log10_solar_radius
end subroutine compute_surface_globals

! ---------------------------------------------------------------
! The asteroseismic observables. The scaling relations, against the
! solar references nu_max_sun, delta_nu_sun and Teff_sun of the
! controls,
!   nu_max = nu_max_sun (g/g_sun) (Teff/Teff_sun)^(-1/2),
!   delta_nu_rho = delta_nu_sun (M/Msun)^(1/2) (R/Rsun)^(-3/2),
! with log g_sun from the run's own solar constants, so that
! Monte-Carlo-scaled solar values stay self-consistent; then the
! asymptotic delta_nu and delta_Pg from the stitched model.
subroutine compute_seismic_observables(s)
      use math_lib
      use stitched_model, only: ip_logR, &
           ip_brunt_N2, ip_csound
      type(star_info), intent(inout) :: s
      double precision :: log_g_solar
      double precision :: acoustic_radius, buoyancy_radius, n2, dr, &
           r_inner, r_here, r_outer
      logical :: entered_g_mode_cavity
      integer :: k

      log_g_solar = cgl + log10(s%solar_mass_cgs) &
           - 2.0d0*log10(s%solar_radius_cgs)
      s%nu_max = s%ctrl%nu_max_sun &
           * exp10((s%log_g_surface - log_g_solar)) &
           * sqrt(s%ctrl%Teff_sun/exp10(s%log_Teff))
! The mean-density scaling.
      s%delta_nu_rho = s%ctrl%delta_nu_sun &
           * sqrt(s%star_mass) * exp10((-1.5d0*s%log_R_surface))

! The asymptotic p-mode large separation from the sound travel time,
!   delta_nu = [2 int_0^R dr/c]^(-1),
! by the trapezoidal rule over the stitched model's points 1..n_ie
! (center to photosphere; the atmosphere above R is excluded, as in
! MESA's acoustic radius). In uHz (the integral is in seconds).
      acoustic_radius = 0.0d0
      do k = 2, s%stitched%n_ie
         if (s%stitched%prof(ip_csound,k) <= 0.0d0 .or. &
             s%stitched%prof(ip_csound,k-1) <= 0.0d0) cycle
         dr = exp10(s%stitched%prof(ip_logR,k)) - exp10(s%stitched%prof(ip_logR,k-1))
         acoustic_radius = acoustic_radius + dr*0.5d0 &
              *(1.0d0/s%stitched%prof(ip_csound,k) + 1.0d0/s%stitched%prof(ip_csound,k-1))
      end do
      if (acoustic_radius > 0.0d0) then
         s%delta_nu = 1.0d6/(2.0d0*acoustic_radius)
      else
         s%delta_nu = 0.0d0
      end if

! The asymptotic l = 1 g-mode period spacing from the Brunt-Vaisala
! frequency,
!   delta_Pg = 2 pi^2/sqrt(l(l+1)) / int N/r dr = sqrt(2) pi^2 / int N/r dr,
! over the innermost g-mode cavity only: walking outward from the
! center, accumulate while N^2 > 0 and stop on first leaving the cavity
! (after WB's MESA treatment). Each point contributes N/r times half
! the span to its neighbours. Zero when there is no cavity (e.g. a fully
! convective interior). In seconds. The center point is skipped (N/r is
! 0/0 there).
      buoyancy_radius = 0.0d0
      entered_g_mode_cavity = .false.
      do k = 2, s%stitched%n_ie - 1
         n2 = s%stitched%prof(ip_brunt_N2,k)
         if (n2 > 0.0d0) then
            entered_g_mode_cavity = .true.
            r_inner = exp10(s%stitched%prof(ip_logR,k-1))
            r_here = exp10(s%stitched%prof(ip_logR,k))
            r_outer = exp10(s%stitched%prof(ip_logR,k+1))
            buoyancy_radius = buoyancy_radius &
                 + sqrt(n2)/r_here*0.5d0*(r_outer - r_inner)
         else if (entered_g_mode_cavity) then
            exit
         end if
      end do
      if (buoyancy_radius > 0.0d0) then
         s%delta_Pg = cpi*cpi*sqrt(2.0d0)/buoyancy_radius
      else
         s%delta_Pg = 0.0d0
      end if
end subroutine compute_seismic_observables

! ---------------------------------------------------------------
! The total moment of inertia: without rotation the thin-shell sum
! I = sum (2/3) dm r^2, with rotation the sum of the shells' i_rot.
subroutine compute_moment_of_inertia(s)
      use math_lib
      type(star_info), intent(inout) :: s
      integer :: i

      s%total_moment_of_inertia = 0.0d0
      if (.not. s%job%rotation_active) then
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
subroutine compute_snu_rates(s)
      type(star_info), intent(inout) :: s
      integer :: i

      s%cl37_snu_rate = 0.0d0
      s%ga71_snu_rate = 0.0d0
      if (s%ctrl%calc_neutrinos) then
         do i = 1, 8
            s%cl37_snu_rate = s%cl37_snu_rate + &
                 clsnuf_diag(i)*s%neutrino_flux_total(i)
            s%ga71_snu_rate = s%ga71_snu_rate + &
                 gasnuf_diag(i)*s%neutrino_flux_total(i)
         end do
      else
         do i = 1, 10
            s%neutrino_flux_total(i) = 0.0d0
         end do
      end if
end subroutine compute_snu_rates

! ---------------------------------------------------------------
! The surface rotation period, P = 2 pi/omega(nz) in days (at most
! 9999), the surface velocity omega(nz) R in km/s, and the moment of
! inertia of the surface convection zone (the i_rot sum with rotation,
! sum (2/3) dm r^2 without). Reads star%log_R_surface, so
! compute_surface_globals runs first.
subroutine compute_rotation_observables(s)
      use math_lib
      type(star_info), intent(inout) :: s
      integer :: i

      s%cz_moment_of_inertia = 0.0d0
      if (s%job%rotation_active) then
         s%rotation_period_days = min(9999.0d0, &
              0.5d0*c4pi/s%omega(s%nz)/8.64d4)
         s%surf_velocity_kms = s%omega(s%nz)* &
              exp(ln10*(s%log_R_surface+s%log10_solar_radius))*1.0d-5
         if (s%envelope_cz_bottom_index .lt. s%nz) then
            do i = s%envelope_cz_bottom_index, s%nz
               s%cz_moment_of_inertia = &
                    s%cz_moment_of_inertia + s%i_rot(i)
            end do
         end if
      else
         s%rotation_period_days = 0.0d0
         s%surf_velocity_kms = 0.0d0
         if (s%envelope_cz_bottom_index .lt. s%nz) then
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
! star%log_R_surface, so compute_surface_globals runs first.
subroutine compute_h_shell_boundaries(s)
      use math_lib
      type(star_info), intent(inout) :: s
      if (s%has_h_shell) then
         s%h_shell_bot_mass = &
              s%m(s%h_shell_zone_begin)/s%solar_mass_cgs
         s%h_shell_mid_mass = &
              s%m(s%h_shell_midpoint_zone)/s%solar_mass_cgs
         s%h_shell_top_mass = &
              s%m(s%h_shell_end_index)/s%solar_mass_cgs
         s%h_shell_bot_radius = exp(ln10*( &
              s%logR(s%h_shell_zone_begin) &
              - s%log_R_surface - s%log10_solar_radius))
         s%h_shell_mid_radius = exp(ln10*( &
              s%logR(s%h_shell_midpoint_zone) &
              - s%log_R_surface - s%log10_solar_radius))
         s%h_shell_top_radius = exp(ln10*( &
              s%logR(s%h_shell_end_index) &
              - s%log_R_surface - s%log10_solar_radius))
      else
         s%h_shell_bot_mass = 0.0d0
         s%h_shell_mid_mass = 0.0d0
         s%h_shell_top_mass = 0.0d0
         s%h_shell_bot_radius = 0.0d0
         s%h_shell_mid_radius = 0.0d0
         s%h_shell_top_radius = 0.0d0
      end if
end subroutine compute_h_shell_boundaries

end module observables

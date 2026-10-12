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
! state and lives at module level: ksaha_center, the Saha-table state
! of the central eos call. The others are set before every use and are
! plain locals. (Legacy wrtout's core-boundary radius, which read the
! previous call's interpolation fraction through an FX/FX2 typo and was
! never used, is not carried over.)
module observables
      use turnover_timescale, only: compute_turnover_timescale
      use star_info_lib
      use eos_lib
      use physical_constants
      implicit none
      private
      public :: compute_observables
      public :: shell_masses_from_log_mass

! Cross-call state (see the header): the central eos call's Saha-table
! state. Zero at process start, and zero again for a new run in the same
! process (observables_reset_pending).
      integer, save :: ksaha_center = 0

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
subroutine compute_observables(ierr)
      integer, intent(out) :: ierr

      ierr = 0

! A new run in the same process (yrec_reset) restarts the central eos
! call's Saha-table state.
      if (observables_reset_pending) then
         ksaha_center = 0
         observables_reset_pending = .false.
      end if

      call renormalize_luminosity_breakdown
      call locate_core_cz
      call compute_central_conditions(ierr)
      if (ierr /= 0) return
      call locate_surface_cz_base
! The turnover timescale (with its _old and pphot0 copies) is physics
! state, read by the wind saturation and the deuterium limiter, and is
! refreshed here for every model.
      call refresh_turnover_timescale(ierr)
      if (ierr /= 0) return

! The per-model history columns (star% members read by write_history),
! by the formulas of the legacy .track output. The surface radius also
! feeds the run log's model line.
      call compute_surface_globals
      call compute_seismic_observables
      call compute_moment_of_inertia
      call compute_snu_rates
      call compute_rotation_observables
      call compute_h_shell_boundaries

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
subroutine renormalize_luminosity_breakdown
      double precision :: total_luminosity_sum, temp_value
      integer :: i

! The post-convergence twin of henyey_iterate's per-iteration
! renormalization, which also guards against a zero sum; kept separate
! on purpose (audit/readability-sweep-2026-09-03/reports/henyey-wave2.md,
! item 6).
      if(.not.star%ctrl%helium_flash_active) then
       total_luminosity_sum = star%luminosity_breakdown(i_lum_pp1)+star%luminosity_breakdown(i_lum_pp2)+ &
            star%luminosity_breakdown(i_lum_pp3)+star%luminosity_breakdown(i_lum_cno)+star%luminosity_breakdown(i_lum_3alpha)+ &
            star%luminosity_breakdown(i_lum_neu)+star%luminosity_breakdown(i_lum_grav)+star%luminosity_breakdown(i_lum_he_c)
       temp_value = star%luminosity_lsun(star%nz)/total_luminosity_sum
       do i = 1,n_lum_channels
          star%luminosity_breakdown(i) = star%luminosity_breakdown(i)*temp_value
       end do
      endif
end subroutine renormalize_luminosity_breakdown

! ---------------------------------------------------------------
! The mass of the convective core in Msun, m at its top zone (zero
! without a convective core).
subroutine locate_core_cz
      if(star%core_cz_top_index.gt.1) then
       star%core_cz_mass = star%m(star%core_cz_top_index)/star%solar_mass_cgs
      else
       star%core_cz_mass = 0.0D0
      endif
end subroutine locate_core_cz

! ---------------------------------------------------------------
! The central P and T, extrapolated from the first zone, and the
! central density, beta and electron degeneracy parameter eta from the
! eos there. Stores star%central_*.
subroutine compute_central_conditions(ierr)
      use math_lib
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
      temp_value =0.5D0*exp(ln10*(cc13*(c4pi3l+star%logRho(1)-star%log_mass(1))+star%logRho(1)+cgl+star%log_mass(1)))
      pressure_linear = exp(ln10*star%logP(1))
      log_pressure_center = log10(pressure_linear + temp_value)
! With the actual temperature gradient del = d ln T/d ln P of zone 1
! (star%gradT(1)), T_c = T_1 (1 + del (P_c - P_1)/P_1).
      log_temperature_center = star%logT(1) + log10(1.0D0+ temp_value*star%gradT(1)/pressure_linear)
      eos_res%log10_density = star%logRho(1)
      hydrogen_fraction_center = star%xa(i_h1,1)
      metal_fraction_center = star%xa(i_metals,1)
      is_atmosphere_point = .true.
      compute_derivatives = .false.
! The eos at (P_c, T_c) for the central density, beta and eta.
      call eos_get(log_temperature_center, log_pressure_center, &
           hydrogen_fraction_center, metal_fraction_center, eos_res, &
           compute_derivatives, is_atmosphere_point, ksaha_center, &
           composition_at_zone=star%xa(:,1), ierr=ierr)
      if (ierr /= 0) return
! Store the central values for later use.
      star%central_log10_pressure = log_pressure_center
      star%central_log10_temperature = log_temperature_center
      star%central_log10_density = eos_res%log10_density
      star%central_beta = eos_res%beta
      star%central_degeneracy_eta = eos_res%electron_degeneracy_parameter
end subroutine compute_central_conditions

! ---------------------------------------------------------------
! The base of the surface convection zone, interpolated between the
! zones around it: star%envelope_mass (the convective envelope's mass,
! Msun), envelope_cz_log_radius and envelope_radius (Rsun), and
! envelope_cz_temperature, _density, _pressure and _opacity. A fully
! convective star takes the central values (so compute_central_conditions
! runs first); without a surface convection zone all are zero.
subroutine locate_surface_cz_base
      use math_lib
      use convection_geometry, only: cz_base_fraction
      double precision :: dd1, dd2, cz_base_mass, envelope_boundary_fx
      double precision :: envelope_cz_log_temperature, &
           envelope_cz_log_density, envelope_cz_log_pressure

      if(star%envelope_cz_bottom_index.lt.star%nz) then
       if(star%envelope_cz_bottom_index.gt.1) then
! The fraction f of the way from zone b-1 to the base zone b where
! gradr - grada crosses zero (cz_base_fraction), and the mass there.
! Linear in f: m, log r, the opacity, log T, log rho and log P.
            dd2 = star%gradr(star%envelope_cz_bottom_index-1)-star%grada(star%envelope_cz_bottom_index-1)
            dd1 = star%gradr(star%envelope_cz_bottom_index)-star%grada(star%envelope_cz_bottom_index)
            envelope_boundary_fx = cz_base_fraction(dd2, dd1)
            cz_base_mass = star%m(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%m(star%envelope_cz_bottom_index)-star%m(star%envelope_cz_bottom_index-1))
            star%envelope_mass = (exp(ln10*star%log_total_mass) - cz_base_mass)/star%solar_mass_cgs
! The radius of the base (MHP 2/98).
            star%envelope_cz_log_radius = star%logR(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%logR(star%envelope_cz_bottom_index)-star%logR(star%envelope_cz_bottom_index-1))-star%log10_solar_radius
            star%envelope_radius = exp(ln10*star%envelope_cz_log_radius)
            star%envelope_cz_opacity = star%opacity_zone(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%opacity_zone(star%envelope_cz_bottom_index)-star%opacity_zone(star%envelope_cz_bottom_index-1))
            envelope_cz_log_temperature = star%logT(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%logT(star%envelope_cz_bottom_index)-star%logT(star%envelope_cz_bottom_index-1))
            envelope_cz_log_density = star%logRho(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%logRho(star%envelope_cz_bottom_index)-star%logRho(star%envelope_cz_bottom_index-1))
            envelope_cz_log_pressure = star%logP(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%logP(star%envelope_cz_bottom_index)-star%logP(star%envelope_cz_bottom_index-1))
            star%envelope_cz_temperature = exp(ln10*envelope_cz_log_temperature)
            star%envelope_cz_density = exp(ln10*envelope_cz_log_density)
            star%envelope_cz_pressure = exp(ln10*envelope_cz_log_pressure)
       else
          star%envelope_mass = star%star_mass
          star%envelope_radius = 0.0D0
            star%envelope_cz_temperature = exp10(star%central_log10_temperature)
            star%envelope_cz_density = exp10(star%central_log10_density)
            star%envelope_cz_pressure = exp10(star%central_log10_pressure)
            star%envelope_cz_opacity = star%opacity_zone(1)
       endif
      else
       star%envelope_mass = 0.0D0
       star%envelope_radius = 0.0D0
         star%envelope_cz_temperature = 0.0D0
         star%envelope_cz_density = 0.0D0
         star%envelope_cz_pressure = 0.0D0
         star%envelope_cz_opacity = 0.0D0
      endif
end subroutine locate_surface_cz_base

! ---------------------------------------------------------------
! Refreshes the turnover timescale and copies it and pphot to their
! _old / pphot0 versions, which the mid-step consumers (the wind, the
! deuterium limiter) read during the next step. compute_turnover_timescale
! also overwrites star%envelope_radius with its own base radius.
subroutine refresh_turnover_timescale(ierr)
      integer, intent(out) :: ierr
      call compute_turnover_timescale(star%envelope_radius, ierr)
      star%convective_turnover_timescale_old = &
           star%convective_turnover_timescale
      star%pphot0 = star%pphot
end subroutine refresh_turnover_timescale

! ---------------------------------------------------------------
! The surface radius from L = 4 pi R^2 sigma Teff^4, log g = log G +
! log M - 2 log R (CGS), then log_R_surface in Rsun.
subroutine compute_surface_globals
      star%log_R_surface = log10_radius_from_l_teff(star%log_L, star%log_Teff, &
           star%log10_solar_luminosity)
      star%log_g_surface = cgl + star%stotal &
           - 2.0d0*star%log_R_surface
      star%log_R_surface = star%log_R_surface - star%log10_solar_radius
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
subroutine compute_seismic_observables
      use math_lib
      use stitched_model, only: n_ie, stx_prof, ip_logR, &
           ip_brunt_N2, ip_csound
      double precision :: log_g_solar
      double precision :: acoustic_radius, buoyancy_radius, n2, dr, &
           r_inner, r_here, r_outer
      logical :: entered_g_mode_cavity
      integer :: k

      log_g_solar = cgl + log10(star%solar_mass_cgs) &
           - 2.0d0*log10(star%solar_radius_cgs)
      star%nu_max = star%ctrl%nu_max_sun &
           * exp10((star%log_g_surface - log_g_solar)) &
           * sqrt(star%ctrl%Teff_sun/exp10(star%log_Teff))
! The mean-density scaling.
      star%delta_nu_rho = star%ctrl%delta_nu_sun &
           * sqrt(star%star_mass) * exp10((-1.5d0*star%log_R_surface))

! The asymptotic p-mode large separation from the sound travel time,
!   delta_nu = [2 int_0^R dr/c]^(-1),
! by the trapezoidal rule over the stitched model's points 1..n_ie
! (center to photosphere; the atmosphere above R is excluded, as in
! MESA's acoustic radius). In uHz (the integral is in seconds).
      acoustic_radius = 0.0d0
      do k = 2, n_ie
         if (stx_prof(ip_csound,k) <= 0.0d0 .or. &
             stx_prof(ip_csound,k-1) <= 0.0d0) cycle
         dr = exp10(stx_prof(ip_logR,k)) - exp10(stx_prof(ip_logR,k-1))
         acoustic_radius = acoustic_radius + dr*0.5d0 &
              *(1.0d0/stx_prof(ip_csound,k) + 1.0d0/stx_prof(ip_csound,k-1))
      end do
      if (acoustic_radius > 0.0d0) then
         star%delta_nu = 1.0d6/(2.0d0*acoustic_radius)
      else
         star%delta_nu = 0.0d0
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
      do k = 2, n_ie - 1
         n2 = stx_prof(ip_brunt_N2,k)
         if (n2 > 0.0d0) then
            entered_g_mode_cavity = .true.
            r_inner = exp10(stx_prof(ip_logR,k-1))
            r_here = exp10(stx_prof(ip_logR,k))
            r_outer = exp10(stx_prof(ip_logR,k+1))
            buoyancy_radius = buoyancy_radius &
                 + sqrt(n2)/r_here*0.5d0*(r_outer - r_inner)
         else if (entered_g_mode_cavity) then
            exit
         end if
      end do
      if (buoyancy_radius > 0.0d0) then
         star%delta_Pg = cpi*cpi*sqrt(2.0d0)/buoyancy_radius
      else
         star%delta_Pg = 0.0d0
      end if
end subroutine compute_seismic_observables

! ---------------------------------------------------------------
! The total moment of inertia: without rotation the thin-shell sum
! I = sum (2/3) dm r^2, with rotation the sum of the shells' i_rot.
subroutine compute_moment_of_inertia
      use math_lib
      integer :: i

      star%total_moment_of_inertia = 0.0d0
      if (.not. star%job%rotation_active) then
         do i = 1, star%nz
            star%total_moment_of_inertia = &
                 star%total_moment_of_inertia + &
                 cc23*star%dm(i)*exp(2.0d0*ln10*star%logR(i))
         end do
      else
         do i = 1, star%nz
            star%total_moment_of_inertia = &
                 star%total_moment_of_inertia + star%i_rot(i)
         end do
      end if
end subroutine compute_moment_of_inertia

! ---------------------------------------------------------------
! The chlorine and gallium capture rates (SNU) from the neutrino flux
! totals; without calc_neutrinos the rates and all ten fluxes are zero.
subroutine compute_snu_rates
      integer :: i

      star%cl37_snu_rate = 0.0d0
      star%ga71_snu_rate = 0.0d0
      if (star%ctrl%calc_neutrinos) then
         do i = 1, 8
            star%cl37_snu_rate = star%cl37_snu_rate + &
                 clsnuf_diag(i)*star%neutrino_flux_total(i)
            star%ga71_snu_rate = star%ga71_snu_rate + &
                 gasnuf_diag(i)*star%neutrino_flux_total(i)
         end do
      else
         do i = 1, 10
            star%neutrino_flux_total(i) = 0.0d0
         end do
      end if
end subroutine compute_snu_rates

! ---------------------------------------------------------------
! The surface rotation period, P = 2 pi/omega(nz) in days (at most
! 9999), the surface velocity omega(nz) R in km/s, and the moment of
! inertia of the surface convection zone (the i_rot sum with rotation,
! sum (2/3) dm r^2 without). Reads star%log_R_surface, so
! compute_surface_globals runs first.
subroutine compute_rotation_observables
      use math_lib
      integer :: i

      star%cz_moment_of_inertia = 0.0d0
      if (star%job%rotation_active) then
         star%rotation_period_days = min(9999.0d0, &
              0.5d0*c4pi/star%omega(star%nz)/8.64d4)
         star%surf_velocity_kms = star%omega(star%nz)* &
              exp(ln10*(star%log_R_surface+star%log10_solar_radius))*1.0d-5
         if (star%envelope_cz_bottom_index .lt. star%nz) then
            do i = star%envelope_cz_bottom_index, star%nz
               star%cz_moment_of_inertia = &
                    star%cz_moment_of_inertia + star%i_rot(i)
            end do
         end if
      else
         star%rotation_period_days = 0.0d0
         star%surf_velocity_kms = 0.0d0
         if (star%envelope_cz_bottom_index .lt. star%nz) then
            do i = star%envelope_cz_bottom_index, star%nz
               star%cz_moment_of_inertia = &
                    star%cz_moment_of_inertia + &
                    cc23*star%dm(i)*exp(2.0d0*ln10*star%logR(i))
            end do
         end if
      end if
end subroutine compute_rotation_observables

! ---------------------------------------------------------------
! The bottom, middle and top of the hydrogen-burning shell: masses in
! Msun and radii as r/R (zero without a shell). Reads
! star%log_R_surface, so compute_surface_globals runs first.
subroutine compute_h_shell_boundaries
      use math_lib
      if (star%has_h_shell) then
         star%h_shell_bot_mass = &
              star%m(star%h_shell_zone_begin)/star%solar_mass_cgs
         star%h_shell_mid_mass = &
              star%m(star%h_shell_midpoint_zone)/star%solar_mass_cgs
         star%h_shell_top_mass = &
              star%m(star%h_shell_end_index)/star%solar_mass_cgs
         star%h_shell_bot_radius = exp(ln10*( &
              star%logR(star%h_shell_zone_begin) &
              - star%log_R_surface - star%log10_solar_radius))
         star%h_shell_mid_radius = exp(ln10*( &
              star%logR(star%h_shell_midpoint_zone) &
              - star%log_R_surface - star%log10_solar_radius))
         star%h_shell_top_radius = exp(ln10*( &
              star%logR(star%h_shell_end_index) &
              - star%log_R_surface - star%log10_solar_radius))
      else
         star%h_shell_bot_mass = 0.0d0
         star%h_shell_mid_mass = 0.0d0
         star%h_shell_top_mass = 0.0d0
         star%h_shell_bot_radius = 0.0d0
         star%h_shell_mid_radius = 0.0d0
         star%h_shell_top_radius = 0.0d0
      end if
end subroutine compute_h_shell_boundaries

end module observables

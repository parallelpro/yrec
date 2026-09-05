!----------------------------------------------------------------------
! observables_lib
!----------------------------------------------------------------------
! 2026: the former core/update_output_diagnostics.f90 (phase four,
! step 5 -- the state-computing blocks that historically lived inside
! io/write_legacy_output.f90, moved so the writer only reads), restructured as a
! module with one subroutine per theme. The public entry point is
! compute_observables, called from evolve_step immediately before the
! output writer, once per output model -- exactly the cadence the
! blocks had inside wrtout. Every computed quantity is stored on
! star% (the observables members and star%luminosity_breakdown);
! the only physics feedback is the turnover-timescale refresh (see
! compute_observables).
!
! Theme order matches the original wrtout order and is load-bearing:
! central conditions must precede the surface-CZ base (its
! fully-convective branch reads star%central_log10_*), and the
! MESA-mode block reads star%log_R_surface set two themes
! earlier.
!
! Statement-level content is verbatim from the original; only the
! decomposition is new. One local of the old blanket-`save` version
! is genuinely cross-call state and lives at module level:
! ksaha_center, eqstat's saha-table state continuity across calls.
! All other former SAVEd locals are assigned before use on every path
! that reads them, so they are plain locals of their theme routine.
! (The original wrtout also carried a core-boundary radius computed
! from the previous call's surface-CZ interpolation fraction -- the
! FX/FX2 typo; that block was never consumed and is deleted.)
module observables_lib
      use star_info_lib
      use eos_lib
      use phys_const_lib
      implicit none
      private
      public :: compute_observables
      public :: shell_masses_from_log_mass

! cross-call state (see header): the saha table continuity. Static
! zero at process start; reset for repeated in-process runs via
! observables_reset_pending.
      integer, save :: ksaha_center = 0

! snu coefficients (as in wrtout) for the MESA-mode snu rates
      double precision, parameter :: gasnuf_diag(8) = [1.18D1,2.15D2, &
           7.14D4,7.17D1,2.40D4,6.04D1,1.137D2,1.139D2]
      double precision, parameter :: clsnuf_diag(8) = [0.0D0,1.6D1, &
           4.26D4,2.4D0,1.14D4,1.7D0,6.8D0,6.9D0]

contains

! ---------------------------------------------------------------
! Driver: one call fills every per-model observable on star%.
subroutine compute_observables(ierr)
      integer, intent(out) :: ierr

      ierr = 0

! 2026 (phase five, step C): see evolve_step's matching block.
      if (observables_reset_pending) then
         ksaha_center = 0
         observables_reset_pending = .false.
      end if

      call renormalize_luminosity_breakdown
      call locate_core_cz
      call compute_central_conditions(ierr)
      if (ierr /= 0) return
      call locate_surface_cz_base
! The turnover timescale (with its pphot/_old lag bookkeeping) is
! PHYSICS state -- the wind saturation and the deuterium limiter
! read it -- so it refreshes in BOTH output modes. (Pre-2026-
! restructure, legacy mode refreshed it via wrtout's own gettau
! call at write time; that call is gone.)
      call refresh_turnover_timescale(ierr)
      if (ierr /= 0) return

! ---- 2026 MESA-style output: fill the per-model history sources ----
! (star% members read by write_history). Formulas are the legacy
! .track v0 branch's, verbatim. Legacy mode skips this: wrtout
! computes the same values internally with pinned print ordering.
! surface radius/gravity run in BOTH modes since the run-log model
! line (2026 log redesign) reads log_R_surface.
      call compute_surface_globals
      call compute_seismic_observables
      call compute_moment_of_inertia
      call compute_snu_rates
      call compute_rotation_observables
      call compute_h_shell_boundaries

end subroutine compute_observables

! ---------------------------------------------------------------
! Unlogged shell masses m(1:nz) (grams, the shell centres) and the
! per-shell mass dm(1:nz) from the log10 mass grid: dm is half the
! span between the neighbouring shell centres, the last shell reaching
! to the total mass. Shared by read_starting_model (star%m, star%dm)
! and neutrino_flux_table (its local copies).
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
! RENORMALIZE LUMINOSITY TERMS TLUMX - SKIPPED FOR HE FLASH
! (MUTATES the model -- star%luminosity_breakdown.)
subroutine renormalize_luminosity_breakdown
      double precision :: total_luminosity_sum, temp_value
      integer :: i

! (post-convergence twin of henyey_iterate's per-iteration
! renormalization, which additionally guards total > 0 -- kept
! separate on purpose, see reports/henyey-wave2.md item 6)
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
!  CALCULATE MASS OF CENTRAL CONVECTION ZONE (SOLAR UNITS).
subroutine locate_core_cz
      if(star%core_cz_top_index.gt.1) then
       star%core_cz_mass = star%m(star%core_cz_top_index)/star%solar_mass_cgs
      else
       star%core_cz_mass = 0.0D0
      endif
end subroutine locate_core_cz

! ---------------------------------------------------------------
!  DETERMINE CENTRAL T,P, AND DENSITY USING THE FIRST SHELL VALUES.
!  CENTRAL ETA AND BETA ARE ALSO CALCULATED. Stores star%central_*.
subroutine compute_central_conditions(ierr)
      use math_lib
      integer, intent(out) :: ierr

      double precision :: temp_value
      double precision :: pressure_linear, log_pressure_center, &
           log_temperature_center, hydrogen_fraction_center, &
           metal_fraction_center
      logical :: is_atmosphere_point, compute_derivatives
! 2026 named-index results: the former 20-variable central-conditions
! soup is one eos result array (see eos_lib's index constants).
      double precision :: eos_res(num_eos_results)

      ierr = 0
      eos_res = 0.0d0

!  EXTRAPOLATE FROM INNER SHELL P AND T TO CENTRAL P AND T
      temp_value =0.5D0*exp(ln10*(cc13*(c4pi3l+star%logRho(1)-star%log_mass(1))+star%logRho(1)+cgl+star%log_mass(1)))
      pressure_linear = exp(ln10*star%logP(1))
      log_pressure_center = log10(pressure_linear + temp_value)
!  SDEL(2,1) IS THE ACTUAL T GRADIENT AT POINT 1( = DEL)
      log_temperature_center = star%logT(1) + log10(1.0D0+ temp_value*star%gradT(1)/pressure_linear)
      eos_res(i_log10_density) = star%logRho(1)
      hydrogen_fraction_center = star%xa(i_h1,1)
      metal_fraction_center = star%xa(i_metals,1)
      is_atmosphere_point = .true.
      compute_derivatives = .false.
!  CALL EQSTAT TO GET TRUE CENTRAL DENSITY, BETA, AND ETA.
      call eos_get(log_temperature_center, log_pressure_center, &
           hydrogen_fraction_center, metal_fraction_center, eos_res, &
           compute_derivatives, is_atmosphere_point, ksaha_center, &
           composition_at_zone=star%xa(:,1), ierr=ierr)
      if (ierr /= 0) return
! STORE CENTRAL RHO,P,T FOR LATER USE
      star%central_log10_pressure = log_pressure_center
      star%central_log10_temperature = log_temperature_center
      star%central_log10_density = eos_res(i_log10_density)
      star%central_beta = eos_res(i_beta)
      star%central_degeneracy_eta = eos_res(i_eta)
end subroutine compute_central_conditions

! ---------------------------------------------------------------
!  Surface-CZ base: interpolate the zone edge at the base of the
!  surface convection zone, filling star%envelope_mass/
!  envelope_radius/envelope_cz_*.
! (The fully-convective branch reads the central conditions, hence
! compute_central_conditions runs first.)
subroutine locate_surface_cz_base
      use math_lib
      double precision :: dd1, dd2, cz_base_mass, envelope_boundary_fx
      double precision :: envelope_cz_log_temperature, &
           envelope_cz_log_density, envelope_cz_log_pressure

      if(star%envelope_cz_bottom_index.lt.star%nz) then
       if(star%envelope_cz_bottom_index.gt.1) then
!  FIND MASS FRACTION OF THE ZONE EDGE AT BASE OF SURFACE C.Z.
            dd2 = star%gradr(star%envelope_cz_bottom_index-1)-star%grada(star%envelope_cz_bottom_index-1)
            dd1 = star%gradr(star%envelope_cz_bottom_index)-star%grada(star%envelope_cz_bottom_index)
            envelope_boundary_fx = dd2/(dd2-dd1)
            cz_base_mass = star%m(star%envelope_cz_bottom_index-1)+envelope_boundary_fx* &
                 (star%m(star%envelope_cz_bottom_index)-star%m(star%envelope_cz_bottom_index-1))
            star%envelope_mass = (exp(ln10*star%log_total_mass) - cz_base_mass)/star%solar_mass_cgs
! MHP 2/98 FIND RADIUS OF CZ BASE
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
! Refresh the turnover timescale for the history column, with the
! lag bookkeeping (the physics consumers -- getw, the wind, the
! deuterium limiter -- drive their own gettau calls; this one only
! freshens the reported value).
subroutine refresh_turnover_timescale(ierr)
      integer, intent(out) :: ierr
      call compute_turnover_timescale(star%envelope_radius, ierr)
      star%convective_turnover_timescale_old = &
           star%convective_turnover_timescale
      star%pphot0 = star%pphot
end subroutine refresh_turnover_timescale

! ---------------------------------------------------------------
! Surface radius and gravity from L and Teff.
subroutine compute_surface_globals
      star%log_R_surface = 0.5d0*(star%log_L + star%log10_solar_luminosity &
           - c4pil - csigl - 4.0d0*star%log_Teff)
      star%log_g_surface = cgl + star%stotal &
           - 2.0d0*star%log_R_surface
      star%log_R_surface = star%log_R_surface - star%log10_solar_radius
end subroutine compute_surface_globals

! ---------------------------------------------------------------
! Asteroseismic scaling-relation observables (2026): nu_max from
! (log g, Teff) and delta_nu_rho from the mean density, against the
! namelist solar references nu_max_sun / delta_nu_sun / Teff_sun.
! The solar log g is derived from the run's own solar constants, so
! Monte-Carlo-scaled solar values stay self-consistent.
subroutine compute_seismic_observables
      use math_lib
      use stitched_model_lib, only: n_ie, stx_prof, ip_logR, &
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
! mean-density scaling: sqrt( (M/Msun) / (R/Rsun)^3 )
      star%delta_nu_rho = star%ctrl%delta_nu_sun &
           * sqrt(star%star_mass) * exp10((-1.5d0*star%log_R_surface))

! Asymptotic p-mode large separation from the sound travel time:
! delta_nu = [2 int_0^R dr/c]^-1, trapezoidal over the stitched
! model's points 1..n_ie (center -> photosphere; the atmosphere
! points above R_star are excluded, as in MESA's acoustic radius).
! In uHz (the integral is in seconds).
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

! Asymptotic l=1 g-mode period spacing from the Brunt-Vaisala
! integral: delta_Pg = 2 pi^2 / sqrt(l(l+1)) / int(N/r dr), the
! integral taken over the INNERMOST g-mode cavity only -- walk
! outward from the center, accumulate while N^2 > 0, stop on first
! leaving the cavity after entering it (after WB's MESA treatment).
! Zero when no cavity exists (e.g. a fully convective interior).
! In seconds. The center point is skipped (N/r is 0/0 there).
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
! Total moment of inertia (thin-shell sum without rotation, i_rot
! sum with).
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
! Chlorine/gallium SNU capture rates from the neutrino flux totals.
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
! Surface rotation period/velocity and the CZ moment of inertia.
! (Reads star%log_R_surface -- compute_surface_globals first.)
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
! H-burning shell boundary masses and (surface-relative) radii.
! (Reads star%log_R_surface -- compute_surface_globals first.)
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

end module observables_lib

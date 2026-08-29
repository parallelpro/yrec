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
! star% (star%*, star%*, star%*,
! star%luminosity_breakdown); nothing here feeds back into the
! physics -- the physics consumers of the turnover timescale call
! gettau themselves at their own points (getw, starin).
!
! Theme order matches the original wrtout order and is load-bearing:
! central conditions must precede the surface-CZ base (its
! fully-convective branch reads star%central_log10_*), and the
! MESA-mode block reads star%log_R_surface set two themes
! earlier.
!
! Statement-level content is verbatim from the original; only the
! decomposition is new. Two locals of the old blanket-`save` version
! are genuinely cross-call state and live at module level:
!   * envelope_boundary_fx -- set by locate_surface_cz_base, but ALSO
!     read by locate_core_cz, which runs EARLIER in the call: it sees
!     the value from the PREVIOUS call. This is the documented FX/FX2
!     typo in the original wrtout.f, preserved exactly, not fixed
!     (its only consumer, core_boundary_radius, is itself dead).
!   * ksaha_center -- eqstat's saha-table state continuity across
!     calls.
! All other former SAVEd locals are assigned before use on every path
! that reads them, so they are plain locals of their theme routine;
! the observables_reset_pending block (repeated-run C API support)
! accordingly shrank to the two real carriers.
!
! NOT here: the legacy-mode gettau call stays in wrtout at its
! original spot. gettau's atm-side integration prints progress
! diagnostics into the .short stream, so hoisting it above wrtout's
! header writes reorders that stream (values identical, layout not)
! -- blocked by print interleaving, not by data flow. Documented
! residual.
module observables_lib
      use star_info_lib, only: star
      use star_info_lib
      use eos_lib
      use phys_const_lib
      implicit none
      private
      public :: compute_observables

! cross-call state (see header): the FX/FX2 stale carry and the saha
! table continuity. Static zero at process start; reset for repeated
! in-process runs via observables_reset_pending.
      double precision, save :: envelope_boundary_fx = 0.0d0
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
      use star_info_lib, only: star
      integer, intent(out) :: ierr

      ierr = 0

! 2026 (phase five, step C): see evolve_step's matching block. This
! includes envelope_boundary_fx, whose previous-call stale value is
! the documented FX/FX2 quirk -- a fresh process starts it at zero,
! so a repeated call must too.
      if (observables_reset_pending) then
         envelope_boundary_fx = 0.0d0
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
      call refresh_turnover_timescale

! ---- 2026 MESA-style output: fill the per-model history sources ----
! (star% members read by write_history). Formulas are the legacy
! .track v0 branch's, verbatim. Legacy mode skips this: wrtout
! computes the same values internally with pinned print ordering.
! surface radius/gravity run in BOTH modes since the run-log model
! line (2026 log redesign) reads log_R_surface.
      call compute_surface_globals
      call compute_moment_of_inertia
      call compute_snu_rates
      call compute_rotation_observables
      call compute_h_shell_boundaries

end subroutine compute_observables

! ---------------------------------------------------------------
! RENORMALIZE LUMINOSITY TERMS TLUMX - SKIPPED FOR HE FLASH
! (MUTATES the model -- star%luminosity_breakdown.)
subroutine renormalize_luminosity_breakdown
      use star_info_lib, only: star
      double precision :: total_luminosity_sum, temp_value
      integer :: i

      if(.not.star%ctrl%helium_flash_active) then
       total_luminosity_sum = star%luminosity_breakdown(i_lum_pp1)+star%luminosity_breakdown(i_lum_pp2)+ &
            star%luminosity_breakdown(i_lum_pp3)+star%luminosity_breakdown(i_lum_cno)+star%luminosity_breakdown(i_lum_3alpha)+ &
            star%luminosity_breakdown(i_lum_neu)+star%luminosity_breakdown(i_lum_grav)+star%luminosity_breakdown(i_lum_he_c)
       temp_value = star%luminosity_lsun(star%nz)/total_luminosity_sum
       do i = 1,8
          star%luminosity_breakdown(i) = star%luminosity_breakdown(i)*temp_value
       end do
      endif
end subroutine renormalize_luminosity_breakdown

! ---------------------------------------------------------------
!  CALCULATE MASS OF CENTRAL CONVECTION ZONE (SOLAR UNITS), plus the
!  (dead) core-boundary radius with the preserved FX/FX2 stale-carry
!  bug -- see the module header.
subroutine locate_core_cz
      use star_info_lib, only: star
      double precision :: core_boundary_fx2
! core_boundary_log_radius/core_boundary_radius (CORERL/CORER) are
! computed but never consumed (original behavior, preserved).
      double precision :: core_boundary_log_radius, core_boundary_radius

      if(star%core_cz_top_index.gt.1) then
       star%core_cz_mass = star%m(star%core_cz_top_index)/star%solar_mass_cgs
      else
       star%core_cz_mass = 0.0D0
      endif

! JVS 10/11 Be more care about the true boundary of the convective core
      if (star%core_cz_top_index.gt.1) then
! JVS 10/11 note: this formula reads envelope_boundary_fx (FX), which
! at this point has not yet been assigned in this call (it is set
! later, in locate_surface_cz_base) -- as module SAVE state it
! carries over whatever value it held at the end of the previous
! call. core_boundary_fx2 (FX2) is computed just above but is NOT
! what is used here -- this looks like a bug (FX2 vs FX typo) in the
! original wrtout.f, preserved exactly, not fixed.
       core_boundary_fx2 = (star%grada(star%core_cz_top_index+1)-star%gradr(star%core_cz_top_index))/ &
             (star%grada(star%core_cz_top_index+1)-star%gradr(star%core_cz_top_index))
       core_boundary_log_radius = star%logR(star%core_cz_top_index)+envelope_boundary_fx* &
            (star%logR(star%core_cz_top_index+1)-star%logR(star%core_cz_top_index))-star%log10_solar_radius
       core_boundary_radius = dexp(ln10*core_boundary_log_radius)
      else
       core_boundary_radius = 0.0D0
      endif
! JVS end
end subroutine locate_core_cz

! ---------------------------------------------------------------
!  DETERMINE CENTRAL T,P, AND DENSITY USING THE FIRST SHELL VALUES.
!  CENTRAL ETA AND BETA ARE ALSO CALCULATED. Stores star%central_*.
subroutine compute_central_conditions(ierr)
      use star_info_lib, only: star
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
      temp_value =0.5D0*dexp(ln10*(cc13*(c4pi3l+star%logRho(1)-star%log_mass(1))+star%logRho(1)+cgl+star%log_mass(1)))
      pressure_linear = dexp(ln10*star%logP(1))
      log_pressure_center = dlog10(pressure_linear + temp_value)
!  SDEL(2,1) IS THE ACTUAL T GRADIENT AT POINT 1( = DEL)
      log_temperature_center = star%logT(1) + dlog10(1.0D0+ temp_value*star%gradT(1)/pressure_linear)
      eos_res(i_log10_density) = star%logRho(1)
      hydrogen_fraction_center = star%xa(i_h1,1)
      metal_fraction_center = star%xa(i_metals,1)
      is_atmosphere_point = .true.
      compute_derivatives = .false.
!  CALL EQSTAT TO GET TRUE CENTRAL DENSITY, BETA, AND ETA.
      call eos_get_r(log_temperature_center, log_pressure_center, &
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
!  envelope_radius/envelope_cz_*. Sets the module-level
!  envelope_boundary_fx (the FX of the FX/FX2 story).
! MHP 02/12 FIXED MINOR GLITCH ON BASE OF THE CONVECTION ZONE
! PROPERTIES FOR FULLY CONVECTIVE STARS; TCENTER PCENTER RHOCENTER
! WERE BEING DEFINED AFTER THIS CODE SECTION (hence central
! conditions run first).
subroutine locate_surface_cz_base
      use star_info_lib, only: star
      double precision :: dd1, dd2, cz_base_mass
      double precision :: envelope_cz_log_temperature, &
           envelope_cz_log_density, envelope_cz_log_pressure

      if(star%envelope_cz_bottom_index.lt.star%nz) then
       if(star%envelope_cz_bottom_index.gt.1) then
!  FIND MASS FRACTION OF THE ZONE EDGE AT BASE OF SURFACE C.Z.
! JVS 10/11/13 SDEL(1,JENV) IN DENOMINATOR WAS A TYPO. CHANGED TO SDEL(3,JENV)
!            FX = (SDEL(3,JENV)-SDEL(1,JENV-1))/
!     *           (SDEL(3,JENV)-SDEL(1,JENV-1))
            dd2 = star%gradr(star%envelope_cz_bottom_index-1)-star%grada(star%envelope_cz_bottom_index-1)
            dd1 = star%gradr(star%envelope_cz_bottom_index)-star%grada(star%envelope_cz_bottom_index)
            envelope_boundary_fx = dd2/(dd2-dd1)
!            HSB = 0.5D0*(HS1(JENV)+HS1(JENV-1))
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
            star%envelope_cz_temperature = 10.0D0**star%central_log10_temperature
            star%envelope_cz_density = 10.0D0**star%central_log10_density
            star%envelope_cz_pressure = 10.0D0**star%central_log10_pressure
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
subroutine refresh_turnover_timescale
      use star_info_lib, only: star
      call compute_turnover_timescale(star%envelope_radius)
      star%convective_turnover_timescale_old = &
           star%convective_turnover_timescale
      star%pphot0 = star%pphot
end subroutine refresh_turnover_timescale

! ---------------------------------------------------------------
! Surface radius and gravity from L and Teff.
subroutine compute_surface_globals
      use star_info_lib, only: star
      star%log_R_surface = 0.5d0*(star%log_L + star%log10_solar_luminosity &
           - c4pil - csigl - 4.0d0*star%log_Teff)
      star%log_g_surface = cgl + star%stotal &
           - 2.0d0*star%log_R_surface
      star%log_R_surface = star%log_R_surface - star%log10_solar_radius
end subroutine compute_surface_globals

! ---------------------------------------------------------------
! Total moment of inertia (thin-shell sum without rotation, i_rot
! sum with).
subroutine compute_moment_of_inertia
      use star_info_lib, only: star
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
      use star_info_lib, only: star
      integer :: i

      star%cl37_snu_rate = 0.0d0
      star%ga71_snu_rate = 0.0d0
      if (star%ctrl%lsnu) then
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
      use star_info_lib, only: star
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
      use star_info_lib, only: star
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

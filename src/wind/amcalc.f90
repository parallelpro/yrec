!----------------------------------------------------------------------
! amcalc
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original amcalc.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! G Somers, 6/16
! Uses the delivered structural variables to determine the strength of
! the magnetic wind torque, according to the Matt et al. (2012)
! formulation. Sets structfactor (in common/cwind/), the structure-
! dependent factor combining mass, radius, luminosity, photospheric
! pressure, and (optionally) convective-turnover-timescale (Rossby)
! scaling, each raised to its own calibrated exponent.
subroutine amcalc(total_mass_msun, log_luminosity_lsun, log_teff)
      use turnover_lib
      use const_lib
      implicit none

      double precision, intent(in) :: total_mass_msun, log_luminosity_lsun, &
           log_teff

! common/const/: only solar_radius_cgs/solar_luminosity_cgs (via
! clsunl) are used here. Naming matches getw.f90.
      double precision :: solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude
      common/const/ solar_luminosity_cgs, log10_solar_luminosity, &
           ln_solar_luminosity, solar_mass_cgs, log10_solar_mass, &
           solar_radius_cgs, log10_solar_radius, solar_bolometric_magnitude



! common/cwind/: exr/exm/exl/expr/exm/extau and structfactor (the
! output) are used here. Naming matches getw.f90.
      double precision :: wind_saturation_omega, exmd, &
           wind_law_omega_exponent, extau, exr, exm, exl, expr, &
           constfactor, structfactor, excen, c_2
      logical :: ljdot0
      common/cwind/ wind_saturation_omega, exmd, wind_law_omega_exponent, &
           extau, exr, exm, exl, expr, constfactor, structfactor, excen, &
           c_2, ljdot0

! common/pmmwind/: Matt-type ("PMM") magnetic-wind calibration block.
! First converted here; also used (with these same names) in cowind.f90
! (no -- cowind.f90 doesn't reference it) and mcowind.f90/mwind.f90.
! pmm_exponent_a/b/c/d/m are dimensionless calibration exponents of the
! Matt+2012/2015-type wind-torque law (their exact individual physical
! roles are not confidently known from this file alone -- see parmin.f
! for how they combine into exw/exr/exm/exl/expr/excen). pmm_norm_jdot
! and pmm_norm_mdot are solar normalization constants (angular-
! momentum-loss and mass-loss rate scales, cgs); pmm_solar_pressure/
! pmm_solar_omega/pmm_solar_turnover_timescale are solar reference
! values used to non-dimensionalize the photospheric-pressure,
! rotation, and Rossby-scaling factors below. use_pmm_wind_law selects
! the Matt-type law (mcowind.f90/mwind.f90) over the legacy
! Kawaler-type law (cowind.f90/wind.f90); scale_by_rossby_number and
! scale_by_b_field are the associated Rossby/field-scaling toggles;
! wind_law_name is a short label (e.g. 'K97') for the adopted law.
      double precision :: pmm_exponent_a, pmm_exponent_b, pmm_exponent_c, &
           pmm_exponent_d, pmm_exponent_m, pmm_norm_jdot, pmm_norm_mdot, &
           pmm_solar_pressure, pmm_solar_omega, pmm_solar_turnover_timescale
      logical :: use_pmm_wind_law, scale_by_rossby_number, scale_by_b_field
      character*3 :: wind_law_name
      common/pmmwind/ pmm_exponent_a, pmm_exponent_b, pmm_exponent_c, &
           pmm_exponent_d, pmm_exponent_m, pmm_norm_jdot, pmm_norm_mdot, &
           pmm_solar_pressure, pmm_solar_omega, pmm_solar_turnover_timescale, &
           use_pmm_wind_law, scale_by_rossby_number, scale_by_b_field, &
           wind_law_name


      save

! --- locals ---
      double precision :: mass_factor, log10_radius, total_radius_cm, &
           radius_rsun, luminosity_lsun, photospheric_pressure_ratio, &
           turnover_ratio

! COLLECT THE RELEVANT STRUCTURE VARIABLES, AND
!     MASS
      mass_factor = total_mass_msun
!     RADIUS
      log10_radius = 0.5d0*(log_luminosity_lsun+log10_solar_luminosity-c4pil- &
           csigl-4.d0*log_teff)
      total_radius_cm = dexp(ln10*log10_radius)
      radius_rsun = total_radius_cm/solar_radius_cgs
!     LUMINOSITY
      luminosity_lsun = 10.**log_luminosity_lsun
!     PHOTOSPHERIC PRESSURE
      photospheric_pressure_ratio = 10.**(turnover%pphot0+turnover%fracstep*(turnover%pphot-turnover%pphot0))/ &
           (10.**pmm_solar_pressure)
!     CONVECTIVE OVERTURN TIMESCALE
      if(scale_by_rossby_number)then
         turnover_ratio = (turnover%convective_turnover_timescale_old+turnover%fracstep* &
              (turnover%convective_turnover_timescale-turnover%convective_turnover_timescale_old))/ &
              pmm_solar_turnover_timescale
      else
         turnover_ratio = 1.
      endif
!     COMBINE THEM ALL
      structfactor = mass_factor**exm * radius_rsun**exr * &
           luminosity_lsun**exl * photospheric_pressure_ratio**expr &
           * turnover_ratio**extau
      return
end subroutine amcalc

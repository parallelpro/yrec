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

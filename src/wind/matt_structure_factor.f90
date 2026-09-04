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
! formulation. Sets star%job%structfactor, the structure-
! dependent factor combining mass, radius, luminosity, photospheric
! pressure, and (optionally) convective-turnover-timescale (Rossby)
! scaling, each raised to its own calibrated exponent.
subroutine matt_structure_factor(total_mass_msun, log_luminosity_lsun, log_teff)
      use star_info_lib, only: star
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: total_mass_msun, log_luminosity_lsun, &
           log_teff
! --- locals ---
      double precision :: mass_factor, log10_radius, total_radius_cm, &
           radius_rsun, luminosity_lsun, photospheric_pressure_ratio, &
           turnover_ratio

! COLLECT THE RELEVANT STRUCTURE VARIABLES, AND
!     MASS
      mass_factor = total_mass_msun
!     RADIUS
      log10_radius = 0.5d0*(log_luminosity_lsun+star%log10_solar_luminosity-c4pil- &
           csigl-4.d0*log_teff)
      total_radius_cm = exp(ln10*log10_radius)
      radius_rsun = total_radius_cm/star%solar_radius_cgs
!     LUMINOSITY
      luminosity_lsun = exp10(log_luminosity_lsun)
!     PHOTOSPHERIC PRESSURE
      photospheric_pressure_ratio = exp10((star%pphot0+star%fracstep*(star%pphot-star%pphot0)))/ &
           (exp10(star%ctrl%pmm_solar_pressure))
!     CONVECTIVE OVERTURN TIMESCALE
      if(star%ctrl%scale_by_rossby_number)then
         turnover_ratio = (star%convective_turnover_timescale_old+star%fracstep* &
              (star%convective_turnover_timescale-star%convective_turnover_timescale_old))/ &
              star%ctrl%pmm_solar_turnover_timescale
      else
         turnover_ratio = 1.
      endif
!     COMBINE THEM ALL
      star%job%structfactor = pow(mass_factor, star%ctrl%exm) * pow(radius_rsun, star%ctrl%exr) * &
           pow(luminosity_lsun, star%ctrl%exl) * pow(photospheric_pressure_ratio, star%ctrl%expr) &
           * pow(turnover_ratio, star%ctrl%extau)
      return
end subroutine matt_structure_factor

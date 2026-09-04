!----------------------------------------------------------------------
! cowind
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original cowind.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! COWIND returns the change in the angular velocity of the surface
! convection zone of a star which is experiencing angular momentum
! loss (legacy Kawaler-type wind law; see wind_spindown_matt.f90 for the
! Matt-type replacement selected when use_pmm_wind_law is true).
!
! INPUT VARIABLES:
!  log_luminosity_lsun : total luminosity
!  full_timestep : diffusion timestep
!  cz_moment_of_inertia : total moment of inertia of the surface c.z.
!  iteration_number : iteration number (in the diffusion routines)
!  omega_surface : surface angular velocity at the end of the previous
!                  iteration
!  total_mass_msun : stellar mass (solar units)
!  log_teff : stellar effective temperature
!  omega_old : surface angular velocity at the beginning of the
!              timestep
!  star%ctrl%exmd, wind_law_omega_exponent, exr, exm, constfactor :
!  user parameters which vary the
!  strength of angular momentum loss and its dependence on surface
!  rotation rate and other stellar parameters.
!
! OUTPUT VARIABLES
!  domega_start,domega_end: change in angular velocity of the surface
!  c.z. for the timestep based on the angular velocity at the
!  beginning and end of the timestep respectively.
!
!  domega_start and domega_end are incorporated into the diffusion
!  equations; they are averaged to find the total angular momentum
!  loss for the timestep.
subroutine wind_spindown(log_luminosity_lsun, full_timestep, cz_moment_of_inertia, &
     iteration_number, omega_surface, total_mass_msun, log_teff, &
     omega_old, domega_start, domega_end, ierr)
      use star_info_lib, only: star
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log_luminosity_lsun, full_timestep, &
           cz_moment_of_inertia
      integer, intent(in) :: iteration_number
      double precision, intent(in) :: omega_surface, total_mass_msun, log_teff
      double precision, intent(in) :: omega_old
      double precision, intent(out) :: domega_start
      double precision, intent(inout) :: domega_end
      integer, intent(out) :: ierr
! --- locals ---
      double precision :: omega_saturation, log10_radius, total_radius_cm, &
           mass_loss_rate_msun_yr, wind_coefficient, omega_old_capped, &
           omega_new_capped, domega_end_this_iter

! MHP 3/09 IF WMAX > 1 THEN ASSUME THAT THE PARAMETER WMAX IS DEFINED BY
! WMAX = WMAX(SUN)*TAUCZ(SUN) AND THE SATURATION THRESHOLD WSAT = WMAX/TAUCZ(STAR)
      ierr = 0

      if(star%job%wind_saturation_omega.gt.1.0d0)then
         if(star%convective_turnover_timescale.gt.1.0d0)then
            omega_saturation = star%job%wind_saturation_omega/star%convective_turnover_timescale
         else
            write(*,911)star%job%wind_saturation_omega,star%convective_turnover_timescale
 911        format('ERROR IN WIND - TAUCZ NOT DEFINED ',1P2E12.3,'STOPPED')
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
            ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
            ! preserve the historical stop on a nonzero return.
            ierr = 1
            return
         endif
      else
         omega_saturation = star%job%wind_saturation_omega
      endif
!  FIND TOTAL RADIUS OF STAR.
      log10_radius=0.5d0*(log_luminosity_lsun+star%log10_solar_luminosity-c4pil- &
           csigl-4.d0*log_teff)
      total_radius_cm = exp(ln10*log10_radius)
! DMDOT IS THE MASS LOSS RATE IN SOLAR MASSES PER YEAR.
      mass_loss_rate_msun_yr = 2.0d-14
! DJ/DT = DT*CONSTFACTOR*(DMDOT/1.0D-14)**EXMD*OMEGA**EXW*(M/MSUN)**EXM
!         *(R/RSUN)**EXR
!  THE CONSTANT AND EXPONENTS ARE SET IN PARMIN BASED ON THE INPUT
!  INDEX ALFA;SEE PARMIN FOR DETAILS ON THE DEPENDENCE OF EACH ON ALFA.
      wind_coefficient = full_timestep/cz_moment_of_inertia*star%ctrl%constfactor* &
           pow((mass_loss_rate_msun_yr/1.0d-14), star%ctrl%exmd) &
           *pow((total_radius_cm/star%solar_radius_cgs), star%ctrl%exr)*pow(total_mass_msun, star%ctrl%exm)
      omega_old_capped = min(omega_old,omega_saturation)
      omega_new_capped = min(omega_surface,omega_saturation)
      domega_start = wind_coefficient*pow(omega_old_capped, &
           star%ctrl%wind_law_omega_exponent-1.0d0)*omega_old
      domega_end_this_iter = wind_coefficient*pow(omega_new_capped, &
           star%ctrl%wind_law_omega_exponent-1.0d0)*omega_surface
      if(iteration_number.eq.1) then
         domega_end = domega_end_this_iter
      else
         domega_end = 0.5d0*(domega_end+domega_end_this_iter)
      endif

      return
end subroutine wind_spindown

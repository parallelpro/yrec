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
! loss (legacy Kawaler-type wind law; see mcowind.f90 for the
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
!       IN COMMON BLOCKS
!  exmd,exw,exr,exm,constfactor : user parameters which vary the
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
subroutine cowind(log_luminosity_lsun, full_timestep, cz_moment_of_inertia, &
     iteration_number, omega_surface, total_mass_msun, log_teff, &
     omega_old, domega_start, domega_end, ierr)
      use star_info_lib, only: star
      use star_info_lib, only: star
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: log_luminosity_lsun, full_timestep, &
           cz_moment_of_inertia
      integer, intent(in) :: iteration_number
      double precision, intent(in) :: omega_surface, total_mass_msun, log_teff
      double precision, intent(in) :: omega_old
      double precision, intent(out) :: domega_start
      double precision, intent(inout) :: domega_end
! --- locals ---
      double precision :: omega_saturation, log10_radius, total_radius_cm, &
           mass_loss_rate_msun_yr, wind_coefficient, omega_old_capped, &
           omega_new_capped, domega_end_this_iter

!  FIND TOTAL RADIUS OF STAR.
! mhp 10/02 cgrav not used; omit
!      CGRAV = EXP(CLN*CGL)
! MHP 3/09 IF WMAX > 1 THEN ASSUME THAT THE PARAMETER WMAX IS DEFINED BY
! WMAX = WMAX(SUN)*TAUCZ(SUN) AND THE SATURATION THRESHOLD WSAT = WMAX/TAUCZ(STAR)
      integer, intent(out) :: ierr

      ierr = 0

      if(wind_saturation_omega.gt.1.0d0)then
         if(star%turnover%convective_turnover_timescale.gt.1.0d0)then
            omega_saturation = wind_saturation_omega/star%turnover%convective_turnover_timescale
!            WRITE(*,912)WSAT,TAUCZ
! 912        FORMAT('Omega sat, Tau',1p2e12.3)
         else
            write(*,911)wind_saturation_omega,star%turnover%convective_turnover_timescale
 911        format('ERROR IN WIND - TAUCZ NOT DEFINED ',1P2E12.3,'STOPPED')
            ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the driver-side
            ! call sites (core/main, core/crrect, core/starin, setup/hpoint)
            ! preserve the historical stop on a nonzero return.
            ierr = 1
            return
         endif
      else
         omega_saturation = wind_saturation_omega
      endif
      log10_radius=0.5d0*(log_luminosity_lsun+log10_solar_luminosity-c4pil- &
           csigl-4.d0*log_teff)
      total_radius_cm = exp(ln10*log10_radius)
! DMDOT IS THE MASS LOSS RATE IN SOLAR MASSES PER YEAR.
      mass_loss_rate_msun_yr = 2.0d-14
! DJ/DT = DT*CONSTFACTOR*(DMDOT/1.0D-14)**EXMD*OMEGA**EXW*(M/MSUN)**EXM
!         *(R/RSUN)**EXR
!  THE CONSTANT AND EXPONENTS ARE SET IN PARMIN BASED ON THE INPUT
!  INDEX ALFA;SEE PARMIN FOR DETAILS ON THE DEPENDENCE OF EACH ON ALFA.
      wind_coefficient = full_timestep/cz_moment_of_inertia*constfactor* &
           (mass_loss_rate_msun_yr/1.0d-14)**exmd &
           *(total_radius_cm/solar_radius_cgs)**exr*total_mass_msun**exm
      omega_old_capped = min(omega_old,omega_saturation)
      omega_new_capped = min(omega_surface,omega_saturation)
      domega_start = wind_coefficient*omega_old_capped** &
           (wind_law_omega_exponent-1.0d0)*omega_old
      domega_end_this_iter = wind_coefficient*omega_new_capped** &
           (wind_law_omega_exponent-1.0d0)*omega_surface
      if(iteration_number.eq.1) then
         domega_end = domega_end_this_iter
      else
         domega_end = 0.5d0*(domega_end+domega_end_this_iter)
      endif

      return
end subroutine cowind

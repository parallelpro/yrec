!----------------------------------------------------------------------
! mcowind
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mcowind.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! PRESERVED BUG (not fixed, per project policy of exact transliteration):
! the FSUN centrifugal-normalization line below uses the bare local
! variable "gl" where mwind.f90's equivalent line (and this file's own
! FCORR1/FCORR2 lines immediately after) use the common/const2/
! gravitational constant "cgl". "gl" is never assigned anywhere in this
! file and is not a member of any common block referenced here, so it
! is an uninitialized (SAVE'd) local -- this is a pre-existing defect
! in mcowind.f, reproduced exactly rather than silently corrected to
! cgl.
!
! COWIND RETURNS THE CHANGE IN THE ANGULAR VELOCITY OF THE SURFACE
! CONVECTION ZONE OF A STAR WHICH IS EXPERIENCING ANGULAR MOMENTUM LOSS.
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
!  exmd,exw,extau,exr,exm,factor : user parameters which vary the
!  strength of angular momentum loss and its dependence on surface
!  rotation rate and other stellar parameters.
!
! OUTPUT VARIABLES
!  domega_start,domega_end:change in angular velocity of the surface c.z. for the timestep
!  based on the angular velocity at the beginning and end of the timestep
!  respectively.
!
!  domega_start and domega_end are incorporated into the diffusion equations;
!  they are averaged to find the total angular momentum loss for the
!  timestep.
subroutine mcowind(log_luminosity_lsun, full_timestep, cz_moment_of_inertia, &
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








      save

! --- locals ---
      double precision :: omega_first, omega_now, current_turnover_timescale, &
           omega_saturation, wind_coefficient, gl, fsun, log10_radius, &
           fcorr1, fcorr2, fcen1, fcen2, omega_old_capped, omega_new_capped, &
           domega_end_this_iter

! RUN OLD WINDLAW IF LMWIND = FALSE
      integer, intent(out) :: ierr

      ierr = 0

      if(.not.use_pmm_wind_law)then
         call cowind(log_luminosity_lsun,full_timestep,cz_moment_of_inertia, &
              iteration_number,omega_surface,total_mass_msun,log_teff, &
              omega_old,domega_start,domega_end, ierr)
         if (ierr /= 0) return
         goto 9999
      endif
!
! G Somers 8/17 CREATE ROTATION DUMMY VARIABLES.
      omega_first = omega_old
      omega_now = omega_surface
!
! ADD ROSSBY SCALING TO WMAX = CRITICAL W FOR THE SUN.
      if(scale_by_rossby_number)then
! MHP 8/17 CORRECTED TAUCZ CALCULATION TO INTERPOLATE PROPERLY IN TIMESTEP
         current_turnover_timescale = star%turnover%convective_turnover_timescale_old+ &
              star%turnover%fracstep*(star%turnover%convective_turnover_timescale-star%turnover%convective_turnover_timescale_old)
         omega_saturation = wind_saturation_omega*pmm_solar_turnover_timescale/ &
              current_turnover_timescale
! G Somers 08/17 IF ADDING ADDITIONAL B SCALING, ADD ADDITIONAL TAUCZ TERM.
         if(scale_by_b_field)then
            omega_first = omega_first*current_turnover_timescale/ &
                 pmm_solar_turnover_timescale
            omega_now = omega_now*current_turnover_timescale/ &
                 pmm_solar_turnover_timescale
            omega_saturation = omega_saturation*current_turnover_timescale/ &
                 pmm_solar_turnover_timescale
         endif
      else
         omega_saturation = wind_saturation_omega
      endif
!
! G Somers 3/17, COMMENTING OUT OLD WIND LAW COEFFICIENT
! CALCULATION. ADDING NEW PMM-STYLE VERSION
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!      RL=0.5D0*(BL+CLSUNL-C4PIL-CSIGL-4.D0*TEFFL)
!      RTOT = EXP(CLN*RL)
! DMDOT IS THE MASS LOSS RATE IN SOLAR MASSES PER YEAR.
!      DMDOT = 2.0D-14
! DJ/DT = DT*FACTOR*(DMDOT/1.0D-14)**EXMD*OMEGA**EXW*(M/MSUN)**EXM
!         *(R/RSUN)**EXR
!  THE CONSTANT AND EXPONENTS ARE SET IN PARMIN BASED ON THE INPUT
!  INDEX ALFA;SEE PARMIN FOR DETAILS ON THE DEPENDENCE OF EACH ON ALFA.
!      C = DT/HICZ*FACTOR*(DMDOT/1.0D-14)**EXMD
!     *    *(RTOT/CRSUN)**EXR*SMASS**EXM
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! NEW CALCULATION. CALL AMCALC TO SET "STRUCTFACTOR"
!
      call amcalc(total_mass_msun,log_luminosity_lsun,log_teff)
!
! CALCULATE THE NEW WIND COEFFICIENT.
!
      wind_coefficient = (full_timestep/cz_moment_of_inertia)*constfactor* &
           structfactor
! MHP 8/17 ADDED CENTRIFUGAL REDUCTION TERM FROM MATT+2012 ApJ 754, L26
! NOTE THAT THIS IS IMPLEMENTED HERE RELATIVE TO THE SUN
!      C_2 = 0.0506
      fsun = 0.5*pmm_solar_omega**2*solar_radius_cgs**3/exp(ln10*gl)/solar_mass_cgs
!     RADIUS
      log10_radius = 0.5d0*(log_luminosity_lsun+log10_solar_luminosity-c4pil- &
           csigl-4.d0*log_teff)
      fcorr1 = 0.5*omega_old**2*exp(ln10*(3.0*log10_radius-cgl))/ &
           total_mass_msun/solar_mass_cgs
      fcorr2 = 0.5*omega_surface**2*exp(ln10*(3.0*log10_radius-cgl))/ &
           total_mass_msun/solar_mass_cgs
      fcen1 = ((c_2**2+fsun)/(c_2**2+fcorr1))**excen
      fcen2 = ((c_2**2+fsun)/(c_2**2+fcorr2))**excen
!
! G Somers, END
      omega_old_capped = min(omega_first,omega_saturation)
      omega_new_capped = min(omega_now,omega_saturation)
      domega_start = wind_coefficient*omega_old_capped** &
           (wind_law_omega_exponent-1.0d0)*omega_old*fcen1
      domega_end_this_iter = wind_coefficient*omega_new_capped** &
           (wind_law_omega_exponent-1.0d0)*omega_surface*fcen2
!      WIND1 = C*WP**(EXW-1.0D0)*WOLD
!      TEMP = C*WN**(EXW-1.0D0)*OMEGAS
      if(iteration_number.eq.1) then
         domega_end = domega_end_this_iter
      else
         domega_end = 0.5d0*(domega_end+domega_end_this_iter)
      endif
 9999 continue
      return
end subroutine mcowind

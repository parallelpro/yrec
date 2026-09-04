!----------------------------------------------------------------------
! mcowind
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mcowind.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! PRESERVED BUG (not fixed; see audit/readability-sweep-2026-09-03/
! SUMMARY.md): the FSUN centrifugal-normalization line below uses the
! bare local variable "gl" where matt_wind.f90's equivalent line (and
! wind_lib's matt_centrifugal_factor, called immediately after) use the
! phys_const_lib gravitational constant "cgl". "gl" is never assigned,
! so under -finit-local-zero it reads 0 and exp(ln10*gl) = 1 -- a
! pre-existing defect in mcowind.f, reproduced exactly rather than
! silently corrected to cgl (which would change numbers).
!
! MCOWIND RETURNS THE CHANGE IN THE ANGULAR VELOCITY OF THE SURFACE
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
!  star%ctrl%wind_law_omega_exponent, extau, exr, exm, constfactor, c_2,
!  excen : user parameters which vary the
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
subroutine wind_spindown_matt(log_luminosity_lsun, full_timestep, cz_moment_of_inertia, &
     iteration_number, omega_surface, total_mass_msun, log_teff, &
     omega_old, domega_start, domega_end, ierr)
      use star_info_lib, only: star
      use phys_const_lib
      use math_lib
      use wind_lib, only: log10_radius_from_l_teff, matt_centrifugal_factor
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
      double precision :: omega_first, omega_now, current_turnover_timescale, &
           omega_saturation, wind_coefficient, gl, fsun, log10_radius, &
           fcen1, fcen2, omega_old_capped, omega_new_capped, &
           domega_end_this_iter

      ierr = 0
! RUN THE LEGACY (KAWALER-TYPE) WIND LAW IF use_pmm_wind_law IS FALSE
      if(.not.star%ctrl%use_pmm_wind_law)then
         call wind_spindown(log_luminosity_lsun,full_timestep,cz_moment_of_inertia, &
              iteration_number,omega_surface,total_mass_msun,log_teff, &
              omega_old,domega_start,domega_end, ierr)
         return
      endif
!
! G Somers 8/17 CREATE ROTATION DUMMY VARIABLES.
      omega_first = omega_old
      omega_now = omega_surface
!
! ADD ROSSBY SCALING TO WMAX = CRITICAL W FOR THE SUN.
      if(star%ctrl%scale_by_rossby_number)then
! MHP 8/17 CORRECTED TAUCZ CALCULATION TO INTERPOLATE PROPERLY IN TIMESTEP
         current_turnover_timescale = star%convective_turnover_timescale_old+ &
              star%fracstep*(star%convective_turnover_timescale-star%convective_turnover_timescale_old)
         omega_saturation = star%job%wind_saturation_omega*star%ctrl%pmm_solar_turnover_timescale/ &
              current_turnover_timescale
! G Somers 08/17 IF ADDING ADDITIONAL B SCALING, ADD ADDITIONAL TAUCZ TERM.
         if(star%ctrl%scale_by_b_field)then
            omega_first = omega_first*current_turnover_timescale/ &
                 star%ctrl%pmm_solar_turnover_timescale
            omega_now = omega_now*current_turnover_timescale/ &
                 star%ctrl%pmm_solar_turnover_timescale
            omega_saturation = omega_saturation*current_turnover_timescale/ &
                 star%ctrl%pmm_solar_turnover_timescale
         endif
      else
         omega_saturation = star%job%wind_saturation_omega
      endif
!
! G Somers 3/17, PMM-STYLE WIND LAW: matt_structure_factor SETS
! star%job%structfactor (THE M, R, L, P_phot, tau_cz DEPENDENCE).
      call matt_structure_factor(total_mass_msun,log_luminosity_lsun,log_teff)
!
! CALCULATE THE NEW WIND COEFFICIENT.
!
      wind_coefficient = (full_timestep/cz_moment_of_inertia)*star%ctrl%constfactor* &
           star%job%structfactor
! MHP 8/17 ADDED CENTRIFUGAL REDUCTION TERM FROM MATT+2012 ApJ 754, L26
! NOTE THAT THIS IS IMPLEMENTED HERE RELATIVE TO THE SUN (star%ctrl%c_2).
      fsun = 0.5*star%ctrl%pmm_solar_omega**2*star%solar_radius_cgs**3/exp(ln10*gl)/star%solar_mass_cgs
!     RADIUS
      log10_radius = log10_radius_from_l_teff(log_luminosity_lsun, log_teff)
      fcen1 = matt_centrifugal_factor(omega_old, fsun, log10_radius, total_mass_msun)
      fcen2 = matt_centrifugal_factor(omega_surface, fsun, log10_radius, total_mass_msun)
!
! G Somers, END
      omega_old_capped = min(omega_first,omega_saturation)
      omega_new_capped = min(omega_now,omega_saturation)
      domega_start = wind_coefficient*pow(omega_old_capped, &
           star%ctrl%wind_law_omega_exponent-1.0d0)*omega_old*fcen1
      domega_end_this_iter = wind_coefficient*pow(omega_new_capped, &
           star%ctrl%wind_law_omega_exponent-1.0d0)*omega_surface*fcen2
      if(iteration_number.eq.1) then
         domega_end = domega_end_this_iter
      else
         domega_end = 0.5d0*(domega_end+domega_end_this_iter)
      endif
      return
end subroutine wind_spindown_matt

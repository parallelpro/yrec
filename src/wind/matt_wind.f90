!----------------------------------------------------------------------
! mwind
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mwind.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 10/02 unused LFIRST removed from call.
! If desired, removes angular momentum from the outer convection zone
! using a Weber-Davis magnetic wind model. Dispatches to the legacy
! Kawaler-type law (kawaler_wind.f90) when use_pmm_wind_law is false, otherwise
! uses the Matt-type ("PMM") law with optional Rossby/B-field scaling
! and the Matt+2012 centrifugal-reduction term (matt_structure_factor.f90's
! structfactor, fcen).
subroutine matt_wind(log_luminosity_lsun, full_timestep, cz_mass_bottom, &
     cz_mass_top, start_zone, end_zone, wind_loss_active, omega_surface, &
     total_mass_msun, log_teff, cz_moment_of_inertia, &
     specific_angular_momentum, ierr)
      use star_info_lib, only: star, json
      use phys_const_lib
      use math_lib
      use wind_lib, only: log10_radius_from_l_teff, matt_centrifugal_factor
      implicit none

      double precision, intent(in) :: log_luminosity_lsun, full_timestep, &
           cz_mass_bottom, cz_mass_top
      integer, intent(in) :: start_zone, end_zone
      logical, intent(in) :: wind_loss_active
      double precision, intent(in) :: omega_surface, total_mass_msun, log_teff, &
           cz_moment_of_inertia
      double precision, intent(inout) :: specific_angular_momentum(json)
! --- locals ---
      double precision :: current_turnover_timescale, omega_now, &
           omega_saturation
      double precision :: fsun, log10_radius, fcen
      double precision :: domega_test
      integer :: num_substeps
      double precision :: sub_timestep
      double precision :: omega_substep_start, omega_iter, omega_iter_new, &
           domega_relative_change, omega_iter_prev
      integer :: substep_idx, iter_count
      double precision :: cz_mass
      double precision :: delta_j_per_mass
      integer :: zone_idx

! IF DESIRED, REMOVE ANGULAR MOMENTUM FROM OUTER CONVECTION ZONE
! USING A WEBER-DAVIS MAGNETIC WIND MODEL
      integer, intent(out) :: ierr

      ierr = 0

      if(.not.star%ctrl%use_pmm_wind_law)then
         call kawaler_wind(log_luminosity_lsun,full_timestep,cz_mass_bottom, &
              cz_mass_top,start_zone,end_zone,wind_loss_active,omega_surface, &
              total_mass_msun,log_teff,cz_moment_of_inertia, &
              specific_angular_momentum, ierr)
         return
      endif
!
! ADD ROSSBY SCALING IF DESIRED.
      if(star%ctrl%scale_by_rossby_number)then
! MHP 8/17 CORRECTED TAUCZ CALCULATION TO INTERPOLATE PROPERLY IN TIMESTEP
         current_turnover_timescale = star%convective_turnover_timescale_old+ &
              star%fracstep*(star%convective_turnover_timescale-star%convective_turnover_timescale_old)
         if(star%ctrl%scale_by_b_field)then
! G Somers 8/17 CREATE ROTATION DUMMY VARIABLES.
            omega_now = omega_surface*current_turnover_timescale/ &
                 star%ctrl%pmm_solar_turnover_timescale
            omega_saturation = star%job%wind_saturation_omega
         else
            omega_now = omega_surface
            omega_saturation = star%job%wind_saturation_omega*star%ctrl%pmm_solar_turnover_timescale/ &
                 current_turnover_timescale
         endif
! If not scaling, just set dummies to the original values.
      else
         omega_now = omega_surface
         omega_saturation = star%job%wind_saturation_omega
      endif
!
! G Somers 3/17, MATT+2012 WIND LAW: matt_structure_factor SETS
! star%job%structfactor (THE M, R, L, P_phot, tau_cz DEPENDENCE).
      call matt_structure_factor(total_mass_msun,log_luminosity_lsun,log_teff)
!
! TEST : THE LOSS RATE DEPENDS ON OMEGA, AND FOR TIMESTESP THAT ARE
! TOO LARGE, ROTATION RATES THAT ARE TOO HIGH, AND THIN SURFACE C.Z.S
! NEGATIVE SURFACE ANGULAR VELOCITIES CAN BE PRODUCED.
! TO AVOID THIS, CHECK THAT THE TIMESTEP IS SMALL ENOUGH TO ALLOW
! A POSITIVE SOLUTION FOR OMEGA IN THE FIRST GUESS AT THE LOSS RATE.
! IF NOT, USE A SERIES OF SMALL STEPS.
! MHP 12/91 CAP LOSS RATE AT WSAT.
! MHP 8/17 ADDED CENTRIFUGAL REDUCTION TERM FROM MATT+2012 ApJ 754, L26
! NOTE THAT THIS IS IMPLEMENTED HERE RELATIVE TO THE SUN (star%ctrl%c_2).
      fsun = 0.5*star%ctrl%pmm_solar_omega**2*star%solar_radius_cgs**3/exp(ln10*cgl)/star%solar_mass_cgs
!     RADIUS
      log10_radius = log10_radius_from_l_teff(log_luminosity_lsun, log_teff)
      fcen = matt_centrifugal_factor(omega_surface, fsun, log10_radius, total_mass_msun)
      domega_test = wind_domega(full_timestep, omega_surface, omega_now, fcen)
      if(domega_test.gt.omega_surface)then
         num_substeps = int(domega_test/omega_surface)+1
         sub_timestep = full_timestep/dfloat(num_substeps)
      else
         num_substeps = 1
         sub_timestep = full_timestep
      endif
      omega_substep_start = omega_surface
      do substep_idx = 1,num_substeps
! THE CONSTANT AND EXPONENTS ARE SET IN PARMIN BASED ON THE INPUT
! INDEX ALFA;SEE PARMIN FOR DETAILS ON THE DEPENDENCE OF EACH ON ALFA.
! ITERATIVE SOLUTION : FOR FIRST GUESS, USE OMEGA=INITIAL OMEGA IN
! COMPUTING LOSS RATE. TO COMPUTE SUBSEQUENT RATES, USE THE AVERAGE OF THE
! STARTING OMEGA AND THE *PREVIOUS* ENDING OMEGA (I.E. AVERAGE EXPLICIT
! AND IMPLICIT SOLUTIONS).
         iter_count = 0
         omega_iter = omega_substep_start
         omega_iter_prev = omega_substep_start
         omega_fixed_point: do   ! (was label 5)
! G Somers 08/17 IF ADDING ADDITIONAL B SCALING, ADD ADDITIONAL TAUCZ TERM.
         if(star%ctrl%scale_by_b_field)then
            omega_now = omega_iter*current_turnover_timescale/ &
                 star%ctrl%pmm_solar_turnover_timescale
         else
            omega_now = omega_iter
         endif
         iter_count = iter_count + 1
! CENTRIFUGAL REDUCTION TERM (MATT+2012) RE-EVALUATED AT THE CURRENT OMEGA.
         fcen = matt_centrifugal_factor(omega_iter, fsun, log10_radius, total_mass_msun)
         omega_iter_new = omega_substep_start - wind_domega(sub_timestep, omega_iter, omega_now, fcen)
         domega_relative_change = 2.0d0*abs((omega_iter_prev-omega_iter_new)/ &
              (omega_iter_prev+omega_iter_new))
         if(domega_relative_change.gt.1.0d-6)then
            omega_iter = 0.5d0*(omega_substep_start+omega_iter_new)
            omega_iter_prev = omega_iter_new
            if(iter_count.le.20)cycle omega_fixed_point
         endif
         exit omega_fixed_point
         end do omega_fixed_point
         omega_substep_start = omega_iter_new
      end do
! DM IS THE TOTAL MASS IN THE CONVECTION ZONE.
      cz_mass = cz_mass_top - cz_mass_bottom
! FIND CHANGE IN ANGULAR MOMENTUM PER UNIT MASS AND SUBTRACT THIS
! NUMBER FROM THE J/M OF EACH SHELL IN THE SURFACE CONVECTION ZONE.
      delta_j_per_mass = (omega_surface-omega_iter_new)*cz_moment_of_inertia/ &
           cz_mass
      do zone_idx = start_zone,end_zone
         specific_angular_momentum(zone_idx) = &
              specific_angular_momentum(zone_idx) - delta_j_per_mass
      end do
      return

contains

!  Decrease in the surface angular velocity over a timestep dt from the
!  Matt-type wind torque at angular velocity omega: dt/I_cz *
!  CONSTFACTOR * STRUCTFACTOR * omega * min(omega_scaled,omega_sat)**
!  (EXW-1) * fcen, where omega_scaled is omega with the optional
!  Rossby/B-field scaling applied and fcen the centrifugal factor.
      double precision function wind_domega(dt, omega, omega_scaled, fcen)
      double precision, intent(in) :: dt, omega, omega_scaled, fcen
      wind_domega = (dt/cz_moment_of_inertia)*star%ctrl%constfactor* &
           star%job%structfactor*omega &
           *pow(min(omega_scaled,omega_saturation), (star%ctrl%wind_law_omega_exponent-1.0d0))*fcen
      end function wind_domega

end subroutine matt_wind

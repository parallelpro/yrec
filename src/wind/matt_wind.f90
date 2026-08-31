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
!      *                SJTOT,SMASS,TEFFL,HICZ,HJM,LFIRST)  ! KC 2025-05-31
      use star_info_lib, only: star, json
      use phys_const_lib
      use math_lib
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
! fcorr_local: a wind centrifugal-correction factor, unrelated to
! former common/ctol/'s fcorr -- renamed (2026) to avoid colliding
! with the const_lib fcorr added for that block.
      double precision :: fsun, log10_radius, fcorr_local, fcen
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
!      *                SJTOT,SMASS,TEFFL,HICZ,HJM)  ! KC 2025-05-31
              total_mass_msun,log_teff,cz_moment_of_inertia, &
              specific_angular_momentum, ierr)
         if (ierr /= 0) return
         continue
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
! G Somers 3/17, CHANGING WINDLAW TO NEW MATT+2012 METHOD.
!
! FIND TOTAL RADIUS OF STAR.
! G Somers - CGRAV, RTOT not used, so blacking out.
! DMDOT IS THE MASS LOSS RATE IN SOLAR MASSES PER YEAR.
! G Somers - Don't need DMDOT anymore.
!         DMDOT = 2.0D-14
! DJ/DT = DT*FACTOR*(DMDOT/1.0D-14)**EXMD*OMEGA**EXW*(M/MSUN)**EXM
!         *(R/RSUN)**EXR
!
! G Somers - New calculation. Call AMCALC to set "STRUCTFACTOR".
!
      call matt_structure_factor(total_mass_msun,log_luminosity_lsun,log_teff)
!
! TEST : THE LOSS RATE DEPENDS ON OMEGA, AND FOR TIMESTESP THAT ARE
! TOO LARGE, ROTATION RATES THAT ARE TOO HIGH, AND THIN SURFACE C.Z.S
! NEGATIVE SURFACE ANGULAR VELOCITIES CAN BE PRODUCED.
! TO AVOID THIS, CHECK THAT THE TIMESTEP IS SMALL ENOUGH TO ALLOW
! A POSITIVE SOLUTION FOR OMEGA IN THE FIRST GUESS AT THE LOSS RATE.
! IF NOT, USE A SERIES OF SMALL STEPS.
! MHP 12/91 CAP LOSS RATE AT WSAT.
!
! G Somers, Commenting out the old DWTEST, adding the PMM version.
!         DWTEST = (DELTS/HICZ)*FACTOR*(DMDOT/1.0D-14)**EXMD
!     *          *OMEGAS*(RTOT/CRSUN)**EXR*SMASS**EXM
!     *          *MIN(OMEGAS,WSAT)**(EXW-1.0D0)
! MHP 8/17 ADDED CENTRIFUGAL REDUCTION TERM FROM MATT+2012 ApJ 754, L26
! NOTE THAT THIS IS IMPLEMENTED HERE RELATIVE TO THE SUN
!      C_2 = 0.0506
      fsun = 0.5*star%ctrl%pmm_solar_omega**2*star%solar_radius_cgs**3/exp(ln10*cgl)/star%solar_mass_cgs
!     RADIUS
      log10_radius = 0.5d0*(log_luminosity_lsun+star%log10_solar_luminosity-c4pil- &
           csigl-4.d0*log_teff)
      fcorr_local = 0.5*omega_surface**2*exp(ln10*(3.0*log10_radius-cgl))/ &
           total_mass_msun/star%solar_mass_cgs
      fcen = pow(((star%ctrl%c_2**2+fsun)/(star%ctrl%c_2**2+fcorr_local)), star%ctrl%excen)
      domega_test = (full_timestep/cz_moment_of_inertia)*star%ctrl%constfactor* &
           star%job%structfactor*omega_surface &
           *pow(min(omega_now,omega_saturation), (star%ctrl%wind_law_omega_exponent-1.0d0))*fcen
!      DWTEST = (DELTS/HICZ)*CONSTFACTOR*STRUCTFACTOR*OMEGAS
!     *          *MIN(OMEGAS,WSAT)**(EXW-1.0D0)
! G Somers END
      if(domega_test.gt.omega_surface)then
         num_substeps = int(domega_test/omega_surface)+1
         sub_timestep = full_timestep/dfloat(num_substeps)
      else
         num_substeps = 1
         sub_timestep = full_timestep
      endif
!      WRITE(*,3)NSTEP
!    3 FORMAT(5X,I5)
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
! G Somers, Commenting out the old WNEW, adding the PMM version.
!         WNEW = WS - (DT/HICZ)*FACTOR*(DMDOT/1.0D-14)**EXMD
!     *          *W*(RTOT/CRSUN)**EXR*SMASS**EXM
!     *          *MIN(W,WSAT)**(EXW-1.0D0)
! MHP 8/17 ADDED CENTRIFUGAL REDUCTION TERM FROM MATT+2012 ApJ 754, L26
! NOTE THAT THIS IS IMPLEMENTED HERE RELATIVE TO THE SUN
         fcorr_local = 0.5*omega_iter**2*exp(ln10*(3.0*log10_radius-cgl))/ &
              total_mass_msun/star%solar_mass_cgs
         fcen = pow(((star%ctrl%c_2**2+fsun)/(star%ctrl%c_2**2+fcorr_local)), star%ctrl%excen)
         omega_iter_new = omega_substep_start - (sub_timestep/ &
              cz_moment_of_inertia)*star%ctrl%constfactor*star%job%structfactor*omega_iter &
              *pow(min(omega_now,omega_saturation), (star%ctrl%wind_law_omega_exponent-1.0d0))*fcen
!         WNEW = WS - (DT/HICZ)*CONSTFACTOR*STRUCTFACTOR*W
!     *          *MIN(W,WSAT)**(EXW-1.0D0)
! G Somers END
         domega_relative_change = 2.0d0*abs((omega_iter_prev-omega_iter_new)/ &
              (omega_iter_prev+omega_iter_new))
!       WRITE(*,4)WS,W,WNEW,DW,HICZ
!    4 FORMAT(1X,1P5E14.6)
         if(domega_relative_change.gt.1.0d-6)then
            omega_iter = 0.5d0*(omega_substep_start+omega_iter_new)
            omega_iter_prev = omega_iter_new
            if(iter_count.le.20)cycle omega_fixed_point
         endif
         exit omega_fixed_point
         end do omega_fixed_point
         omega_substep_start = omega_iter_new
      end do
!     CON = DELTS*FACTOR*(DMDOT/1.0D-14)**EXMD*OMEGAS**(EXW-1.0D0)
!    *        *(RTOT/CRSUN)**EXR*SMASS**EXM
!     FJDOT = CON*OMEGAS/(1.0D0+(EXW*CON/HICZ))
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
end subroutine matt_wind

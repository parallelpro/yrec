!----------------------------------------------------------------------
! wind
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original wind.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 10/02 unused LFIRST removed from call.
! Legacy Kawaler-type Weber-Davis magnetic wind angular-momentum-loss
! routine, called from mwind.f90 when use_pmm_wind_law is false. If
! desired, removes angular momentum from the outer convection zone,
! iterating (implicit-average) omega down over num_substeps sub-steps
! and then distributing the resulting specific-angular-momentum change
! uniformly (per unit mass) over the surface convection zone
! [start_zone,end_zone].
subroutine wind(log_luminosity_lsun, full_timestep, cz_mass_bottom, &
     cz_mass_top, start_zone, end_zone, wind_loss_active, omega_surface, &
     total_mass_msun, log_teff, cz_moment_of_inertia, &
     specific_angular_momentum, ierr)
!      *                SJTOT,SMASS,TEFFL,HICZ,HJM,LFIRST)  ! KC 2025-05-31
      use star_info_lib, only: star, json
      use const_lib
      implicit none

      double precision, intent(in) :: log_luminosity_lsun, full_timestep, &
           cz_mass_bottom, cz_mass_top
      integer, intent(in) :: start_zone, end_zone
      logical, intent(in) :: wind_loss_active
      double precision, intent(in) :: omega_surface, total_mass_msun, log_teff, &
           cz_moment_of_inertia
      double precision, intent(inout) :: specific_angular_momentum(json)
! --- locals ---
      double precision :: omega_saturation
      double precision :: gravity_cgs, log10_radius, total_radius_cm
      double precision :: mass_loss_rate_msun_yr
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
! MHP 3/09 IF WMAX > 1 THEN ASSUME THAT THE PARAMETER WMAX IS DEFINED BY
! WMAX = WMAX(SUN)*TAUCZ(SUN) AND THE SATURATION THRESHOLD WSAT = WMAX/TAUCZ(STAR)
! ONLY APPLY THIS IF ANGULAR MOMENTUM TRANSPORT ENABLED
      integer, intent(out) :: ierr

      ierr = 0

      if(.not.star%job%instability_transport_active)then
         omega_saturation = star%job%wind_saturation_omega
      else if(star%job%wind_saturation_omega.gt.1.0d0)then
         if(star%turnover%convective_turnover_timescale.gt.1.0d0)then
            omega_saturation = star%job%wind_saturation_omega/star%turnover%convective_turnover_timescale
!            WRITE(*,912)WSAT,TAUCZ
! 912        FORMAT('Omega sat, Tau',1p2e12.3)
         else
            write(*,911)star%job%wind_saturation_omega,star%turnover%convective_turnover_timescale
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
      if (wind_loss_active) then
! FIND TOTAL RADIUS OF STAR.
         gravity_cgs = dexp(ln10*cgl)
         log10_radius=0.5d0*(log_luminosity_lsun+star%log10_solar_luminosity-c4pil- &
              csigl-4.d0*log_teff)
         total_radius_cm = dexp(ln10*log10_radius)
! DMDOT IS THE MASS LOSS RATE IN SOLAR MASSES PER YEAR.
         mass_loss_rate_msun_yr = 2.0d-14
! DJ/DT = DT*CONSTFACTOR*(DMDOT/1.0D-14)**EXMD*OMEGA**EXW*(M/MSUN)**EXM
!         *(R/RSUN)**EXR
! TEST : THE LOSS RATE DEPENDS ON OMEGA, AND FOR TIMESTESP THAT ARE
! TOO LARGE, ROTATION RATES THAT ARE TOO HIGH, AND THIN SURFACE C.Z.S
! NEGATIVE SURFACE ANGULAR VELOCITIES CAN BE PRODUCED.
! TO AVOID THIS, CHECK THAT THE TIMESTEP IS SMALL ENOUGH TO ALLOW
! A POSITIVE SOLUTION FOR OMEGA IN THE FIRST GUESS AT THE LOSS RATE.
! IF NOT, USE A SERIES OF SMALL STEPS.
! MHP 12/91 CAP LOSS RATE AT WSAT.
         domega_test = (full_timestep/cz_moment_of_inertia)*star%ctrl%constfactor* &
              (mass_loss_rate_msun_yr/1.0d-14)**star%ctrl%exmd &
              *omega_surface*(total_radius_cm/star%solar_radius_cgs)**star%ctrl%exr* &
              total_mass_msun**star%ctrl%exm &
              *min(omega_surface,omega_saturation)**(star%ctrl%wind_law_omega_exponent-1.0d0)
         if(domega_test.gt.omega_surface)then
            num_substeps = int(domega_test/omega_surface)+1
            sub_timestep = full_timestep/dfloat(num_substeps)
         else
            num_substeps = 1
            sub_timestep = full_timestep
         endif
!         WRITE(*,3)NSTEP
!    3    FORMAT(5X,I5)
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
            iter_count = iter_count + 1
            omega_iter_new = omega_substep_start - (sub_timestep/ &
                 cz_moment_of_inertia)*star%ctrl%constfactor* &
                 (mass_loss_rate_msun_yr/1.0d-14)**star%ctrl%exmd &
                 *omega_iter*(total_radius_cm/star%solar_radius_cgs)**star%ctrl%exr* &
                 total_mass_msun**star%ctrl%exm &
                 *min(omega_iter,omega_saturation)**(star%ctrl%wind_law_omega_exponent-1.0d0)
            domega_relative_change = 2.0d0*abs((omega_iter_prev-omega_iter_new)/ &
                 (omega_iter_prev+omega_iter_new))
!         WRITE(*,4)WS,W,WNEW,DW,HICZ
!    4    FORMAT(1X,1P5E14.6)
            if(domega_relative_change.gt.1.0d-6)then
               omega_iter = 0.5d0*(omega_substep_start+omega_iter_new)
               omega_iter_prev = omega_iter_new
               if(iter_count.le.20)cycle omega_fixed_point
            endif
            exit omega_fixed_point
            end do omega_fixed_point
            omega_substep_start = omega_iter_new
         end do
!        CON = DELTS*CONSTFACTOR*(DMDOT/1.0D-14)**EXMD*OMEGAS**(EXW-1.0D0)
!    *           *(RTOT/CRSUN)**EXR*SMASS**EXM
!        FJDOT = CON*OMEGAS/(1.0D0+(EXW*CON/HICZ))
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
      endif
      return
end subroutine wind

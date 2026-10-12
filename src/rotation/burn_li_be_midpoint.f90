! ---------------------------------------------------------------
! COMMENTED OUT (2026-10-11, the author's decision): the Bulirsch-Stoer
! mixing-plus-burning path (burn_settle_mix -> burn_mix_extrapolated,
! burn_li_be_midpoint). A hard-coded switch kept it off and no control
! could enable it, so it never ran in the modernized code. The code is
! kept below for reference. To revive it, uncomment the three files, the
! two commented-out call sites in evolve_angular_momentum and its tables
! in rotation_scratch_lib; it has not been compiled since 2026-10 and may
! need updating.
! ---------------------------------------------------------------
!!----------------------------------------------------------------------
!! burn_li_be_midpoint
!!----------------------------------------------------------------------
!! Modernized (free-form, readable names) 2026 as part of the YREC
!! readability refactor. Logic and numerics are unchanged from the
!! original liburn2.f; only variable names, source form, and comment
!! style were updated. Validated against the Stage 0 regression suite
!! (examples/run_standard_solar_model).
!!
!! Determine lithium-6, lithium-7, and beryllium-9 burning. The burning
!! rates depend on the local T and rho, the abundance, and (in
!! convection zones) the mixing.
!!
!! In radiative regions burning is done implicitly, using a single
!! trapezoidal (midpoint-in-the-log) estimate of the burning rate over
!! the timestep -- unlike burn_li_be, this variant does not iterate
!! sub-step refinement with a rational-function extrapolation.
!!
!! In the convective region the program likewise uses a single
!! time-averaged (midpoint-in-the-log) rate over the timestep, again
!! without the iterative sub-step refinement used by burn_li_be.
!!
!! The degree of lithium burning in a surface CZ depends sensitively
!! on the temperature at its base -- so accurately locating it is
!! important. Determine the true location (cz_base_frac) of the base
!! of the CZ at the end of the timestep, and the location of the edge
!! of overshoot regions if applicable.
!!
!! 11/91 HR added to call.
!!
!! burn_li_be_midpoint (formerly liburn2)
!!   star reads: rate_be9, rate_be9_start, rate_li6, rate_li6_start,
!!     rate_li7, rate_li7_start
!!   star writes: none
!subroutine burn_li_be_midpoint(timestep, composition, radius, mass_coordinate, &
!     shell_mass, log_temperature, env_cz_zone, env_cz_zone_old, num_zones)
!      use rotation_scratch_lib
!      use star_info_lib, only: star, json, i_li6, i_li7, i_be9
!      use physical_constants
!      use math_lib
!      use convection_geometry, only: locate_surface_mixing_base
!      use burn_step, only: average_li_be_over_cz
!      implicit none
!
!      double precision, intent(in) :: timestep
!      double precision, intent(inout) :: composition(15,json)
!      double precision, intent(in) :: radius(json)
!      double precision, intent(in) :: mass_coordinate(json)
!      double precision, intent(inout) :: shell_mass(json)
!      double precision, intent(in) :: log_temperature(json)
!      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones
!
!      double precision :: li6_substep_depletion, li7_substep_depletion, &
!           be9_substep_depletion
!      integer :: zone_idx, min_zone, max_zone
!      double precision :: cz_base_frac
!      integer :: cz_base_zone, cz_base_zone_old
!      double precision :: log_rate_li6, log_rate_li7, log_rate_be9
!      double precision :: li6_cz_start, li7_cz_start, be9_cz_start
!      double precision :: log_rate_li6_cz_start, log_rate_li7_cz_start, &
!           log_rate_be9_cz_start
!      double precision :: li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end
!      double precision :: x_cz_start(3), x_cz_end(3), log_rate_cz_start(3), &
!           log_rate_cz_end(3)
!      double precision :: log_rate_li6_cz_end, log_rate_li7_cz_end, &
!           log_rate_be9_cz_end
!      double precision :: cz_log_rate_li6, cz_log_rate_li7, cz_log_rate_be9
!      double precision :: li6_depletion, li7_depletion, be9_depletion
!      double precision :: mass_coord_beg, mass_coord_end, delta_mass, &
!           radiative_frac
!
!! THE DEGREE OF LITHIUM BURNING IN A SURFACE CZ DEPENDS SENSITIVELY
!! ON THE TEMPERATURE AT ITS BASE - SO ACCURATELY LOCATING IS IMPORTANT.
!! DETERMINE THE TRUE LOCATION (FX) OF THE BASE OF THE CZ AT THE END
!! OF THE TIMESTEP, AND THE LOCATION OF THE EDGE OF OVERSHOOT REGIONS
!! IF APPLICABLE.
!      call locate_surface_mixing_base(radius, env_cz_zone, env_cz_zone_old, num_zones, &
!           cz_base_zone, cz_base_zone_old, cz_base_frac)
!! RADIATIVE INTERIOR.
!      min_zone = min(cz_base_zone,cz_base_zone_old)
!      max_zone = max(cz_base_zone,cz_base_zone_old)
!      do zone_idx = 1,min_zone-1
!         if(star%rate_be9(zone_idx).le.1.0d-32 .or. star%rate_be9_start(zone_idx).le.1.0d-32)cycle
!         if(composition(i_li6,zone_idx).lt.1.0d-24.and.composition(i_li7,zone_idx).lt.1.0d-24 &
!         .and.composition(i_be9,zone_idx).lt.1.0d-24)cycle
!         if(log_temperature(zone_idx).gt.7.0d0)then
!            composition(i_li6,zone_idx) = 0.0d0
!            composition(i_li7,zone_idx) = 0.0d0
!            composition(i_be9,zone_idx) = 0.0d0
!            cycle
!         endif
!         log_rate_li6 = 0.5d0*(log(star%rate_li6(zone_idx)) + log(star%rate_li6_start(zone_idx)))
!         log_rate_li7 = 0.5d0*(log(star%rate_li7(zone_idx)) + log(star%rate_li7_start(zone_idx)))
!         log_rate_be9 = 0.5d0*(log(star%rate_be9(zone_idx)) + log(star%rate_be9_start(zone_idx)))
!         li6_substep_depletion = timestep*exp(log_rate_li6)
!         li7_substep_depletion = timestep*exp(log_rate_li7)
!         be9_substep_depletion = timestep*exp(log_rate_be9)
!! THE REACTION RATES ARE OF THE FORM
!!    DLI/DT = F(RHO,T,X)*LI =>DLN LI = DT*F(RHO,T,X)
!! SOLVE FOR D LN (SPECIES)/DT AND ZERO OUT IF THE DEPLETION IS TOO HIGH.
!         if(li6_substep_depletion.lt.3.0d1)then
!            composition(i_li6,zone_idx) = composition(i_li6,zone_idx)/ &
!                 exp(li6_substep_depletion)
!         else
!            composition(i_li6,zone_idx) = 0.0d0
!         endif
!         if(li7_substep_depletion.lt.3.0d1)then
!            composition(i_li7,zone_idx) = composition(i_li7,zone_idx)/ &
!                 exp(li7_substep_depletion)
!         else
!            composition(i_li7,zone_idx) = 0.0d0
!         endif
!         if(be9_substep_depletion.lt.3.0d1)then
!            composition(i_be9,zone_idx) = composition(i_be9,zone_idx)/ &
!                 exp(be9_substep_depletion)
!         else
!            composition(i_be9,zone_idx) = 0.0d0
!         endif
!! WRITE NEW ABUNDANCES AND EXIT.
!         if(composition(i_li6,zone_idx).lt.1.0d-24)composition(i_li6,zone_idx)=0.0d0
!         if(composition(i_li7,zone_idx).lt.1.0d-24)composition(i_li7,zone_idx)=0.0d0
!         if(composition(i_be9,zone_idx).lt.1.0d-24)composition(i_be9,zone_idx)=0.0d0
!      end do
!! CONVECTION ZONE.
!!
!! SKIP IF WHOLE CZ IS BELOW THE BURNING THRESHOLD.
!      if (star%rate_be9_start(cz_base_zone_old).le.1.0d-32.or.star%rate_be9(cz_base_zone).le.1.0d-32) then
!         return
!      end if
!      call average_li_be_over_cz(composition, mass_coordinate, shell_mass, &
!           cz_base_zone, cz_base_zone_old, cz_base_frac, num_zones, &
!           x_cz_start, log_rate_cz_start, x_cz_end, cz_mass_end, log_rate_cz_end)
!      li6_cz_start = x_cz_start(1)
!      li7_cz_start = x_cz_start(2)
!      be9_cz_start = x_cz_start(3)
!      log_rate_li6_cz_start = log_rate_cz_start(1)
!      log_rate_li7_cz_start = log_rate_cz_start(2)
!      log_rate_be9_cz_start = log_rate_cz_start(3)
!      li6_cz_end = x_cz_end(1)
!      li7_cz_end = x_cz_end(2)
!      be9_cz_end = x_cz_end(3)
!      log_rate_li6_cz_end = log_rate_cz_end(1)
!      log_rate_li7_cz_end = log_rate_cz_end(2)
!      log_rate_be9_cz_end = log_rate_cz_end(3)
!! FOR THE STARTING ABUNDANCE, USE THE MIXED ABUNDANCE FOR THE
!! DEEPER CZ.
!      if(cz_base_zone.lt.cz_base_zone_old)then
!         li6_cz_start = li6_cz_end
!         li7_cz_start = li7_cz_end
!         be9_cz_start = be9_cz_end
!      endif
!! INITIALIZE ABUNDANCES.
!      li6_cz_end = li6_cz_start
!      li7_cz_end = li7_cz_start
!      be9_cz_end = be9_cz_start
!! TIME-AVERAGED RATES USING LINEAR INTERPOLATION IN THE LOG.
!      cz_log_rate_li6 = 0.5d0*(log_rate_li6_cz_end+log_rate_li6_cz_start)
!      cz_log_rate_li7 = 0.5d0*(log_rate_li7_cz_end+log_rate_li7_cz_start)
!      cz_log_rate_be9 = 0.5d0*(log_rate_be9_cz_end+log_rate_be9_cz_start)
!      li6_depletion = timestep*exp(cz_log_rate_li6)
!      li7_depletion = timestep*exp(cz_log_rate_li7)
!      be9_depletion = timestep*exp(cz_log_rate_be9)
!      if(li6_depletion.lt.3.0d1 .and. li6_cz_end.gt.1.0d-24)then
!         li6_cz_end = li6_cz_end/exp(li6_depletion)
!      else
!         li6_cz_end = 0.0d0
!      endif
!      if(li6_cz_end.lt.1.0d-24)li6_cz_end = 0.0d0
!      if(li7_depletion.lt.3.0d1 .and. li7_cz_end.gt.1.0d-24)then
!         li7_cz_end = li7_cz_end/exp(li7_depletion)
!      else
!         li7_cz_end = 0.0d0
!      endif
!      if(li7_cz_end.lt.1.0d-24)li7_cz_end = 0.0d0
!      if(be9_depletion.lt.3.0d1 .and. be9_cz_end.gt.1.0d-24)then
!         be9_cz_end = be9_cz_end/exp(be9_depletion)
!      else
!         be9_cz_end = 0.0d0
!      endif
!      if(be9_cz_end.lt.1.0d-24)be9_cz_end = 0.0d0
!      do zone_idx = max_zone,num_zones
!         composition(i_li6,zone_idx) = li6_cz_end
!         composition(i_li7,zone_idx) = li7_cz_end
!         composition(i_be9,zone_idx) = be9_cz_end
!      end do
!! NOW SOLVE FOR ABUNDANCES IN THE REGION WHICH BEGAN CONVECTIVE AND
!! ENDED RADIATIVE.
!      if (cz_base_zone.le.cz_base_zone_old) then
!         return
!      end if
!! FIND STARTING AND ENDING LOCATION IN MASS OF THE CZ BASE, AND
!! ASSUME THAT THE FRACTION OF TIME SPENT IN THE CZ IS PROPORTIONAL TO
!! THE LOCATION IN MASS (I.E. A POINT 1/3 FROM THE OLD TO THE NEW CZ
!! BASE IN MASS SPENT 1/3 OF THE TIME IN THE CZ AND 2/3 OUT OF IT).
!      if(cz_base_zone_old.gt.1)then
!         mass_coord_beg = 0.5d0*(mass_coordinate(cz_base_zone_old)+mass_coordinate(cz_base_zone_old-1))
!      else
!         mass_coord_beg = 0.0d0
!      endif
!      mass_coord_end = 0.5d0*(mass_coordinate(cz_base_zone)+mass_coordinate(cz_base_zone-1))
!! 2026 (bugsweep Batch 2): this block runs only for a retreating CZ
!! (cz_base_zone > cz_base_zone_old), so the original
!! delta_mass = beg - end was negative and FRAD = (m - beg)/delta_mass
!! ran from 0 at the old base to -1 at the new one. Per the comment
!! above FRAD must be 1 at the old base and 0 at the new base:
!! FRAD = (end - m)/(end - beg).
!      delta_mass = mass_coord_end - mass_coord_beg
!      do zone_idx = cz_base_zone_old,cz_base_zone-1
!! MHP 9/91 CHANGE TO AVOID DIVISION BY ZERO.
!! SKIP IF SHELL TEMPERATURE DROPS BELOW BURNING THRESHOLD.
!         if(star%rate_be9(zone_idx).le.1.0d-32)exit
!         radiative_frac = (mass_coord_end - mass_coordinate(zone_idx))/delta_mass
!! USE FRAD*RADIATIVE RATE AND (1-FRAD)*CONVECTIVE RATE.
!         li6_depletion = timestep*exp(radiative_frac*log(star%rate_li6(zone_idx))+ &
!              (1.0d0-radiative_frac)*log_rate_li6_cz_start)
!         li7_depletion = timestep*exp(radiative_frac*log(star%rate_li7(zone_idx))+ &
!              (1.0d0-radiative_frac)*log_rate_li7_cz_start)
!         be9_depletion = timestep*exp(radiative_frac*log(star%rate_be9(zone_idx))+ &
!              (1.0d0-radiative_frac)*log_rate_be9_cz_start)
!!***REMEMBER TO ADD FAILSAFES FOR LARGE DEPLETION***
!         composition(i_li6,zone_idx) = composition(i_li6,zone_idx)/exp(li6_depletion)
!         if(composition(i_li6,zone_idx).lt.1.0d-24)composition(i_li6,zone_idx)=0.0d0
!         composition(i_li7,zone_idx) = composition(i_li7,zone_idx)/exp(li7_depletion)
!         if(composition(i_li7,zone_idx).lt.1.0d-24)composition(i_li7,zone_idx)=0.0d0
!         composition(i_be9,zone_idx) = composition(i_be9,zone_idx)/exp(be9_depletion)
!         if(composition(i_be9,zone_idx).lt.1.0d-24)composition(i_be9,zone_idx)=0.0d0
!      end do
!      return
!end subroutine burn_li_be_midpoint

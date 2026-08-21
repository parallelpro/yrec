!----------------------------------------------------------------------
! liburn
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original liburn.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Determine lithium-6, lithium-7, and beryllium-9 burning in a
! non-rotating model (see liburn2.f90 for the simpler single-pass
! rotating-model variant). The burning rates depend on the local T and
! rho, the abundance, and (in convection zones) the mixing.
!
! In radiative regions burning is done implicitly, correcting the
! rates for the change in abundance as a function of time, by
! repeatedly halving the sub-step size and using a rational-function
! extrapolation (ratext) to converge on the zero-step-size limit.
!
! In the convective region the program iterates between burning and
! mixing in the same way, to correct the rates for the change in
! abundance as a function of time.
!
! 11/91 HR added to call.
subroutine liburn(timestep, composition, radius, mass_coordinate, &
     shell_mass, log_temperature, env_cz_zone, env_cz_zone_old, num_zones)
      use luout_lib
      use const_lib
      use nuclear_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: timestep
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: radius(json)
      double precision, intent(in) :: mass_coordinate(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(in) :: log_temperature(json)
      integer, intent(in) :: env_cz_zone, env_cz_zone_old, num_zones



! common/newrat/: lithium/beryllium burning rates at the end of the
! timestep, at the current (possibly overshoot-adjusted) depth.
      double precision :: rate_li6(json), rate_li7(json), rate_be9(json)
      common/newrat/ rate_li6, rate_li7, rate_be9

! common/oldmod/: previous-timestep model snapshot; only old_radius and
! old_composition are used here. Naming matches dburn.f90.
      double precision :: old_pressure(json), old_temperature(json), &
           old_radius(json), old_luminosity(json), old_density(json), &
           old_composition(15,json), old_shell_mass(json), old_teff
      logical :: old_convective_flag(json), old_cz_flag(json)
      integer :: old_num_zones
      common/oldmod/ old_pressure, old_temperature, old_radius, &
           old_luminosity, old_density, old_composition, old_shell_mass, &
           old_convective_flag, old_cz_flag, old_teff, old_num_zones

! common/oldrat/: lithium/beryllium burning rates at the start of the
! timestep.
      double precision :: rate_li6_start(json), rate_li7_start(json), &
           rate_be9_start(json)
      common/oldrat/ rate_li6_start, rate_li7_start, rate_be9_start

! common/scrtch/: only del_grad (originally SDEL) is used here, to
! evaluate del(ad)-del(rad) at the last convective point. The other
! members are declared only to preserve the storage layout.
! del_grad(1,*)/(3,*) are inferred (from the DELAM/DELRM analog in
! common/mdphy/ below) to be the radiative/adiabatic temperature
! gradients; the meaning of (2,*) is not used here and not confirmed.
      double precision :: sesum(json), seg(7,json), sbeta(json), seta(json)
      logical :: locons(json)
      double precision :: so(json), del_grad(3,json), sfxion(3,json), &
           svel(json), scp(json)
      common/scrtch/ sesum, seg, sbeta, seta, locons, so, del_grad, &
           sfxion, svel, scp

! common/dpmix/: only envelope_overshoot_active (originally LOVSTE) is
! used here -- uncertain interpretation, inferred as a flag gating
! whether the true (overshoot-adjusted) envelope convection-zone base
! location is used.
      double precision :: dpenv, alphac, alphae, alpham, betac
      integer :: iov1, iov2, iovim
      logical :: lovstc, envelope_overshoot_active, lovstm, lsemic, ladov, &
           lovmax
      common/dpmix/ dpenv, alphac, alphae, alpham, betac, iov1, iov2, &
           iovim, lovstc, envelope_overshoot_active, lovstm, lsemic, ladov, &
           lovmax

! common/liov/: pressure scale heights used to search downward from the
! CZ base for the true (overshoot-corrected) base location.
      double precision :: pressure_scale_height_start, &
           pressure_scale_height_end
      common/liov/ pressure_scale_height_start, pressure_scale_height_end

! common/rot/: only rotation_active (LROT) and instability_transport_
! active (LINSTB) are used here -- LINSTB's interpretation (rotational
! instability-driven mixing/transport enabled) is inferred from context
! and marked uncertain.
      double precision :: wnew, walpcz, acfpft
      integer :: itfp1, itfp2
      logical :: rotation_active, instability_transport_active, lwnew
      common/rot/ wnew, walpcz, acfpft, itfp1, itfp2, rotation_active, &
           instability_transport_active, lwnew

! common/mdphy/: only del_adiabatic_mix and del_radiative_mix
! (originally DELAM/DELRM) are used here -- the rotating-model
! counterpart of common/scrtch/'s del_grad.
      double precision :: amum(json), cpm(json), delm(json), &
           del_adiabatic_mix(json), del_radiative_mix(json), esumm(json), &
           om(json), qdtm(json), thdifm(json), velm(json), viscm(json), &
           epsm(json)
      common/mdphy/ amum, cpm, delm, del_adiabatic_mix, del_radiative_mix, &
           esumm, om, qdtm, thdifm, velm, viscm, epsm

! common/prevcz/: remembers the previous end-of-timestep values for use
! as the new beginning-of-timestep values.
      double precision :: cz_base_radius_prev, log_rate_li6_prev, &
           log_rate_li7_prev, log_rate_be9_prev
      integer :: envelope_cz_base_zone_prev
      common/prevcz/ cz_base_radius_prev, log_rate_li6_prev, &
           log_rate_li7_prev, log_rate_be9_prev, envelope_cz_base_zone_prev

! common/deuter/: not used in this file; declared only to preserve
! layout. Naming matches dburn.f90.
      double precision :: deuterium_burning_rate(json), &
           deuterium_burning_rate_start(json), accreted_mass_fraction
      integer :: jcz
      common/deuter/ deuterium_burning_rate, deuterium_burning_rate_start, &
           accreted_mass_fraction, jcz

! common/masschg/: use_mass_accretion, mass_accretion_rate, and
! accreted_composition are used here. Naming matches dburn.f90.
      double precision :: mass_accretion_rate, fczdmdt, ftotdmdt, &
           accreted_composition(15), creim
      logical :: lreimer, use_mass_accretion
      common/masschg/ mass_accretion_rate, fczdmdt, ftotdmdt, &
           accreted_composition, creim, lreimer, use_mass_accretion

      double precision :: li6_substep_depletion(json), &
           li7_substep_depletion(json), be9_substep_depletion(json)
      double precision :: light_element_save(3,json)
      integer :: substep_counts(11)
      double precision :: extrap_tol(3), extrap_y(3), extrap_err(3), &
           extrap_result(3)
      data substep_counts/2,4,6,8,12,16,24,32,48,64,96/
      data extrap_tol/1.0d-6,1.0d-6,1.0d-6/
      save

      integer :: zone_idx, refine_idx, substep_idx, species_idx
      integer :: cz_base_zone, cz_base_zone_old, min_zone, max_zone
      double precision :: del_diff, del_diff_below, cz_base_frac
      double precision :: search_radius, shell_radius, delta_radius, &
           cz_base_radius
      double precision :: substep_frac, substep_dt
      logical :: converged
      double precision :: log_rate_li6, log_rate_li7, log_rate_be9
      double precision :: extrap_x
      double precision :: li6_cz_start, li7_cz_start, be9_cz_start, &
           cz_mass_start
      double precision :: log_rate_li6_cz_start, log_rate_li7_cz_start, &
           log_rate_be9_cz_start
      double precision :: shell_mass_save
      double precision :: li6_cz_end, li7_cz_end, be9_cz_end, cz_mass_end
      double precision :: log_rate_li6_cz_end, log_rate_li7_cz_end, &
           log_rate_be9_cz_end
      logical :: accretion_active
      double precision :: li6_accreted, li7_accreted, be9_accreted
      double precision :: cz_log_rate_li6, cz_log_rate_li7, cz_log_rate_be9
      double precision :: li6_depletion, li7_depletion, be9_depletion
      double precision :: li6_added, li7_added, be9_added
      double precision :: mass_coord_beg, mass_coord_end, delta_mass, &
           radiative_frac

! SAVE ORIGINAL ABUNDANCES.
      do 10 zone_idx = 1,num_zones
         light_element_save(1,zone_idx) = composition(13,zone_idx)
         light_element_save(2,zone_idx) = composition(14,zone_idx)
         light_element_save(3,zone_idx) = composition(15,zone_idx)
   10 continue
! THE DEGREE OF LITHIUM BURNING IN A SURFACE CZ DEPENDS SENSITIVELY
! ON THE TEMPERATURE AT ITS BASE - SO ACCURATELY LOCATING IS IMPORTANT.
! DETERMINE THE TRUE LOCATION (FX) OF THE BASE OF THE CZ AT THE END
! OF THE TIMESTEP, AND THE LOCATION OF THE EDGE OF OVERSHOOT REGIONS
! IF APPLICABLE.
      if(env_cz_zone.gt.1.and.env_cz_zone.lt.num_zones)then
         if(rotation_active.and.instability_transport_active)then
            del_diff = del_adiabatic_mix(env_cz_zone) - &
                 del_radiative_mix(env_cz_zone)
            del_diff_below = del_adiabatic_mix(env_cz_zone-1) - &
                 del_radiative_mix(env_cz_zone-1)
         else
! EVALUATE DEL(AD) - DEL(RAD) AT THE LAST CONVECTIVE POINT AND THE ONE
! BELOW IT.
            del_diff = del_grad(3,env_cz_zone)-del_grad(1,env_cz_zone)
            del_diff_below = del_grad(3,env_cz_zone-1)-del_grad(1,env_cz_zone-1)
         endif
! USE LINEAR INTERPOLATION TO FIND THE DISTANCE OF THE TRUE LOCATION
! OF THE BASE FROM THE ZONE MIDPOINT. IF FX IS NEGATIVE,THEN THE TRUE
! BASE IS HIGHER; IF IT IS POSITIVE, THE TRUE BASE IS LOWER.
         cz_base_frac = max(-0.5d0,0.5d0-del_diff_below/ &
              (del_diff_below-del_diff))
         cz_base_frac = min(0.5d0,cz_base_frac)
         if(.not.envelope_overshoot_active)then
            cz_base_zone = env_cz_zone
            cz_base_zone_old = env_cz_zone_old
         else
! STARTING CZ DEPTH
            if(cz_base_radius_prev.eq.0.0d0)then
               cz_base_radius_prev = 0.5d0*(exp(ln10*old_radius(env_cz_zone_old)) &
                        +exp(ln10*old_radius(env_cz_zone_old-1)))
               search_radius = cz_base_radius_prev - pressure_scale_height_start
               do zone_idx = env_cz_zone_old-1,1,-1
                  shell_radius = exp(ln10*old_radius(zone_idx))
                  if(shell_radius.lt.search_radius)then
                     cz_base_zone_old = zone_idx + 1
                     goto 11
                  endif
               end do
               cz_base_zone_old = 1
   11          continue
            else
               cz_base_zone_old = envelope_cz_base_zone_prev
            endif
! ENDING CZ DEPTH : DETERMINE OVERSHOOT FROM TRUE CZ BASE.
            delta_radius = exp(ln10*radius(env_cz_zone))-exp(ln10*radius(env_cz_zone-1))
            cz_base_radius = 0.5d0*(exp(ln10*radius(env_cz_zone))+exp(ln10*radius(env_cz_zone-1))) &
                    -cz_base_frac*delta_radius
            cz_base_radius_prev = cz_base_radius
            search_radius = cz_base_radius - pressure_scale_height_end
            do zone_idx = env_cz_zone-1,1,-1
               shell_radius = exp(ln10*radius(zone_idx))
               if(shell_radius.lt.search_radius)then
                  cz_base_zone = zone_idx + 1
                  delta_radius = exp(ln10*old_radius(zone_idx+1))-shell_radius
                  cz_base_frac = 0.5d0-((search_radius-shell_radius)/delta_radius)
                  cz_base_frac = max(-0.5d0,cz_base_frac)
                  cz_base_frac = min(0.5d0,cz_base_frac)
                  goto 12
               endif
            end do
            cz_base_zone = 1
   12       continue
         endif
      else
         cz_base_zone = env_cz_zone
         cz_base_zone_old = env_cz_zone_old
      endif
      envelope_cz_base_zone_prev = cz_base_zone
! RADIATIVE INTERIOR.
      min_zone = min(cz_base_zone,cz_base_zone_old)
      max_zone = max(cz_base_zone,cz_base_zone_old)
      do 50 zone_idx = 1,min_zone-1
         if(rate_be9(zone_idx).le.1.0d-32 .or. rate_be9_start(zone_idx).le.1.0d-32)goto 60
         if(composition(13,zone_idx).lt.1.0d-24.and.composition(14,zone_idx).lt.1.0d-24 &
         .and.composition(15,zone_idx).lt.1.0d-24)goto 50
         if(log_temperature(zone_idx).gt.7.0d0)then
            composition(13,zone_idx) = 0.0d0
            composition(14,zone_idx) = 0.0d0
            composition(15,zone_idx) = 0.0d0
            goto 50
         endif
         do 40 refine_idx = 1,11
! PERFORM BURNING IN substep_counts SMALLER TIMESTEPS.
!  DEFINE FRAC - THE FRACTION OF THE TOTAL TIMESTEP IN EACH SMALL ONE -
!  AND DELT - THE FRACTIONAL TIMESTEP IN YEARS.
            substep_frac = 1.0d0/dfloat(substep_counts(refine_idx))
            substep_dt = substep_frac*timestep
            converged = .false.
! RESTORE INITIAL ABUNDANCES.
            composition(13,zone_idx) = light_element_save(1,zone_idx)
            composition(14,zone_idx) = light_element_save(2,zone_idx)
            composition(15,zone_idx) = light_element_save(3,zone_idx)
! STORE STARTING REACTION RATES.
            log_rate_li6 = log(rate_li6_start(zone_idx))-0.5d0*substep_frac* &
                   (log(rate_li6(zone_idx))-log(rate_li6_start(zone_idx)))
            log_rate_li7 = log(rate_li7_start(zone_idx))-0.5d0*substep_frac* &
                   (log(rate_li7(zone_idx))-log(rate_li7_start(zone_idx)))
            log_rate_be9 = log(rate_be9_start(zone_idx))-0.5d0*substep_frac* &
                   (log(rate_be9(zone_idx))-log(rate_be9_start(zone_idx)))
            do 20 substep_idx = 1,substep_counts(refine_idx)
! INCREMENT THE REACTION RATES.
               log_rate_li6 = log_rate_li6+substep_frac* &
                    (log(rate_li6(zone_idx))-log(rate_li6_start(zone_idx)))
               log_rate_li7 = log_rate_li7+substep_frac* &
                    (log(rate_li7(zone_idx))-log(rate_li7_start(zone_idx)))
               log_rate_be9 = log_rate_be9+substep_frac* &
                    (log(rate_be9(zone_idx))-log(rate_be9_start(zone_idx)))
               li6_substep_depletion(zone_idx) = substep_dt*exp(log_rate_li6)
               li7_substep_depletion(zone_idx) = substep_dt*exp(log_rate_li7)
               be9_substep_depletion(zone_idx) = substep_dt*exp(log_rate_be9)
! THE REACTION RATES ARE OF THE FORM
!    DLI/DT = F(RHO,T,X)*LI =>DLN LI = DT*F(RHO,T,X)
! SOLVE FOR D LN (SPECIES)/DT AND ZERO OUT IF THE DEPLETION IS TOO HIGH.
               if(li6_substep_depletion(zone_idx).lt.3.0d1)then
                  composition(13,zone_idx) = composition(13,zone_idx)/ &
                       exp(li6_substep_depletion(zone_idx))
               else
                  composition(13,zone_idx) = 0.0d0
               endif
               if(li7_substep_depletion(zone_idx).lt.3.0d1)then
                  composition(14,zone_idx) = composition(14,zone_idx)/ &
                       exp(li7_substep_depletion(zone_idx))
               else
                  composition(14,zone_idx) = 0.0d0
               endif
               if(be9_substep_depletion(zone_idx).lt.3.0d1)then
                  composition(15,zone_idx) = composition(15,zone_idx)/ &
                       exp(be9_substep_depletion(zone_idx))
               else
                  composition(15,zone_idx) = 0.0d0
               endif
   20       continue
! STORE ABUNDANCES AS A FUNCTION OF TIMESTEP IN VECTOR YEST.
            extrap_y(1) = composition(13,zone_idx)
            extrap_y(2) = composition(14,zone_idx)
            extrap_y(3) = composition(15,zone_idx)
            extrap_x = substep_frac**2
! USE A RATIONAL POLYNOMIAL EXTRAPOLATOR TO EXTRAPOLATE TO THE SOLUTION FOR
! ZERO TIMESTEP.
            call ratext(refine_idx,extrap_x,extrap_y,extrap_result,extrap_err,3,7)
! CHECK IF THE SOLUTION HAS CONVERGED.
            if(refine_idx.gt.1)then
               converged = .true.
               do 30 species_idx=1,3
                  if(extrap_y(species_idx).lt.1.0d-24)goto 30
                  extrap_err(species_idx) = abs(extrap_err(species_idx)/extrap_y(species_idx))
                  if(extrap_err(species_idx).gt.extrap_tol(species_idx))converged=.false.
   30          continue
               if(converged)then
!                 WRITE(ISHORT,912)I,J,(YEXT(K2)K2=1,3)
!   912             FORMAT(1X,'CONVERGED',I5,' LEVEL ',I2,' LI6 ',1P,
!      *            E12.5,' LI7 ',E12.5,' BE9 ',E12.5)
                  goto 45
               endif
            endif
   40    continue
! IF THE PROGRAM GETS HERE THEN IT FAILED TO CONVERGE TO WITHIN
! THE SPECIFIED TOLERANCE IN THE MAXIMUM NUMBER OF ITERATIONS.
         write(short_file_unit,911)zone_idx,(extrap_err(species_idx),species_idx=1,3)
  911    format(1x,'***LIBURN CONVERGENCE FAILURE IN SHELL ',i4, &
         'ERRORS '/1p3e10.3)
   45    continue
! WRITE NEW ABUNDANCES AND EXIT.
         composition(13,zone_idx)=extrap_result(1)
         if(composition(13,zone_idx).lt.1.0d-24)composition(13,zone_idx)=0.0d0
         composition(14,zone_idx)=extrap_result(2)
         if(composition(14,zone_idx).lt.1.0d-24)composition(14,zone_idx)=0.0d0
         composition(15,zone_idx)=extrap_result(3)
         if(composition(15,zone_idx).lt.1.0d-24)composition(15,zone_idx)=0.0d0
   50 continue
   60 continue
! CONVECTION ZONE.
!
! SKIP IF WHOLE CZ IS BELOW THE BURNING THRESHOLD.
      if(rate_be9_start(cz_base_zone_old).le.1.0d-32.or.rate_be9(cz_base_zone).le.1.0d-32)goto 200
! FIND RATES AT THE BEGINNING OF THE TIMESTEP (USING THE DEPTH AT THE START).
      li6_cz_start = 0.0d0
      li7_cz_start = 0.0d0
      be9_cz_start = 0.0d0
      cz_mass_start = 0.0d0
      do 65 zone_idx = cz_base_zone_old,num_zones
         li6_cz_start = li6_cz_start+composition(13,zone_idx)*shell_mass(zone_idx)
         li7_cz_start = li7_cz_start+composition(14,zone_idx)*shell_mass(zone_idx)
         be9_cz_start = be9_cz_start+composition(15,zone_idx)*shell_mass(zone_idx)
         cz_mass_start = cz_mass_start + shell_mass(zone_idx)
   65 continue
!    67 CONTINUE
      li6_cz_start = li6_cz_start/cz_mass_start
      li7_cz_start = li7_cz_start/cz_mass_start
      be9_cz_start = be9_cz_start/cz_mass_start
      if(log_rate_li6_prev.le.0.0d0)then
! COMPUTE MASS-WEIGHTED AVERAGE RATES AT THE START OF THE STEP.
         log_rate_li6_cz_start = 0.0d0
         log_rate_li7_cz_start = 0.0d0
         log_rate_be9_cz_start = 0.0d0
         do 68 zone_idx = cz_base_zone_old,num_zones
            log_rate_li6_cz_start = log_rate_li6_cz_start + rate_li6_start(zone_idx)*shell_mass(zone_idx)
            log_rate_li7_cz_start = log_rate_li7_cz_start + rate_li7_start(zone_idx)*shell_mass(zone_idx)
            log_rate_be9_cz_start = log_rate_be9_cz_start + rate_be9_start(zone_idx)*shell_mass(zone_idx)
   68    continue
         log_rate_li6_cz_start = log(log_rate_li6_cz_start/cz_mass_start)
         log_rate_li7_cz_start = log(log_rate_li7_cz_start/cz_mass_start)
         log_rate_be9_cz_start = log(log_rate_be9_cz_start/cz_mass_start)
      else
! USE THE RATE FROM THE END OF THE PREVIOUS TIMESTEP.
         log_rate_li6_cz_start = log_rate_li6_prev
         log_rate_li7_cz_start = log_rate_li7_prev
         log_rate_be9_cz_start = log_rate_be9_prev
      endif
! USE THE LOCATION OF THE TRUE EDGE OF THE CONVECTION ZONE (FX, FOUND AT
! BEGINNING OF SR) TO ADJUST THE BURNING RATE AND MASS OF THE BOTTOM POINT
! SUCH THAT IT INCLUDES THE ENTIRE C.Z.
      shell_mass_save = shell_mass(cz_base_zone)
      if(cz_base_zone.gt.1.and.cz_base_zone.lt.num_zones)then
         shell_mass(cz_base_zone) = shell_mass(cz_base_zone)+cz_base_frac* &
              (mass_coordinate(cz_base_zone)-mass_coordinate(cz_base_zone-1))
         rate_li6(cz_base_zone) = rate_li6(cz_base_zone)+0.5d0*cz_base_frac* &
              (rate_li6(cz_base_zone-1)-rate_li6(cz_base_zone))
         rate_li7(cz_base_zone) = rate_li7(cz_base_zone)+0.5d0*cz_base_frac* &
              (rate_li7(cz_base_zone-1)-rate_li7(cz_base_zone))
         rate_be9(cz_base_zone) = rate_be9(cz_base_zone)+0.5d0*cz_base_frac* &
              (rate_be9(cz_base_zone-1)-rate_be9(cz_base_zone))
      endif
! FIND RATES AT THE END OF THE TIMESTEP (USING THE DEPTH AT THE END).
! ALSO STORE INITIAL ABUNDANCES(FLI60,FLI70,FBE90).
! FM IS THE TOTAL MASS IN THE CZ.
      li6_cz_end = 0.0d0
      li7_cz_end = 0.0d0
      be9_cz_end = 0.0d0
      cz_mass_end = 0.0d0
      log_rate_li6_cz_end = 0.0d0
      log_rate_li7_cz_end = 0.0d0
      log_rate_be9_cz_end = 0.0d0
      do 70 zone_idx = cz_base_zone,num_zones
         log_rate_li6_cz_end = log_rate_li6_cz_end + rate_li6(zone_idx)*shell_mass(zone_idx)
         log_rate_li7_cz_end = log_rate_li7_cz_end + rate_li7(zone_idx)*shell_mass(zone_idx)
         log_rate_be9_cz_end = log_rate_be9_cz_end + rate_be9(zone_idx)*shell_mass(zone_idx)
         li6_cz_end = li6_cz_end+composition(13,zone_idx)*shell_mass(zone_idx)
         li7_cz_end = li7_cz_end+composition(14,zone_idx)*shell_mass(zone_idx)
         be9_cz_end = be9_cz_end+composition(15,zone_idx)*shell_mass(zone_idx)
         cz_mass_end = cz_mass_end + shell_mass(zone_idx)
   70 continue
!    75 CONTINUE
      li6_cz_end = li6_cz_end/cz_mass_end
      li7_cz_end = li7_cz_end/cz_mass_end
      be9_cz_end = be9_cz_end/cz_mass_end
! MASS-WEIGHTED AVERAGE RATES AT THE END OF THE STEP.
      log_rate_li6_cz_end = log(log_rate_li6_cz_end/cz_mass_end)
      log_rate_li7_cz_end = log(log_rate_li7_cz_end/cz_mass_end)
      log_rate_be9_cz_end = log(log_rate_be9_cz_end/cz_mass_end)
! RESTORE ORIGINAL MASS TO BASE OF C.Z.
      shell_mass(cz_base_zone) = shell_mass_save
! FOR THE STARTING ABUNDANCE, USE THE MIXED ABUNDANCE FOR THE
! DEEPER CZ.
      if(cz_base_zone.lt.cz_base_zone_old)then
         li6_cz_start = li6_cz_end
         li7_cz_start = li7_cz_end
         be9_cz_start = be9_cz_end
      endif
      do 120 refine_idx = 1,11
! PERFORM BURNING IN substep_counts SMALLER TIMESTEPS.
!  DEFINE FRAC - THE FRACTION OF THE TOTAL TIMESTEP IN EACH SMALL ONE -
!  AND DELT - THE FRACTIONAL TIMESTEP IN YEARS.
         substep_frac = 1.0d0/dfloat(substep_counts(refine_idx))
         substep_dt = substep_frac*timestep
         converged = .false.
! INITIALIZE ABUNDANCES.
         li6_cz_end = li6_cz_start
         li7_cz_end = li7_cz_start
         be9_cz_end = be9_cz_start
! INITIALIZE RATES.
         cz_log_rate_li6 = log_rate_li6_cz_start-0.5d0*substep_frac* &
              (log_rate_li6_cz_end-log_rate_li6_cz_start)
         cz_log_rate_li7 = log_rate_li7_cz_start-0.5d0*substep_frac* &
              (log_rate_li7_cz_end-log_rate_li7_cz_start)
         cz_log_rate_be9 = log_rate_be9_cz_start-0.5d0*substep_frac* &
              (log_rate_be9_cz_end-log_rate_be9_cz_start)
! BURN IN substep_counts SMALL TIMESTEPS.
! MHP 05/02 ADDED THE IMPACT OF MASS ACCRETION ON LIGHT
! ELEMENT ABUNDANCES.  INITIALIZE ARRAYS; SPECIES ARE
! ASSUMED TO EXPERIENCE BURNING FOR ON AVERAGE 1/2 OF
! THE STEP WHERE THEY ARE INITIALLY ACCRETED AND ARE
! THEN FULLY BURNED.
         if(use_mass_accretion.and.mass_accretion_rate.gt.0.0d0)then
            accretion_active = .true.
            li6_accreted = 0.0d0
            li7_accreted = 0.0d0
            be9_accreted = 0.0d0
         else
            accretion_active = .false.
         endif
         do 100 substep_idx = 1,substep_counts(refine_idx)
! INCREMENT RATES BY LINEAR INTERPOLATION IN THE LOG.
            cz_log_rate_li6 = cz_log_rate_li6 + substep_frac* &
                 (log_rate_li6_cz_end-log_rate_li6_cz_start)
            cz_log_rate_li7 = cz_log_rate_li7 + substep_frac* &
                 (log_rate_li7_cz_end-log_rate_li7_cz_start)
            cz_log_rate_be9 = cz_log_rate_be9 + substep_frac* &
                 (log_rate_be9_cz_end-log_rate_be9_cz_start)
            li6_depletion = substep_dt*exp(cz_log_rate_li6)
            li7_depletion = substep_dt*exp(cz_log_rate_li7)
            be9_depletion = substep_dt*exp(cz_log_rate_be9)
            if(li6_depletion.lt.3.0d1 .and. li6_cz_end.gt.1.0d-24)then
               li6_cz_end = li6_cz_end/exp(li6_depletion)
            else
               li6_cz_end = 0.0d0
            endif
            if(li6_cz_end.lt.1.0d-24)li6_cz_end = 0.0d0
            if(li7_depletion.lt.3.0d1 .and. li7_cz_end.gt.1.0d-24)then
               li7_cz_end = li7_cz_end/exp(li7_depletion)
            else
               li7_cz_end = 0.0d0
            endif
            if(li7_cz_end.lt.1.0d-24)li7_cz_end = 0.0d0
            if(be9_depletion.lt.3.0d1 .and. be9_cz_end.gt.1.0d-24)then
               be9_cz_end = be9_cz_end/exp(be9_depletion)
            else
               be9_cz_end = 0.0d0
            endif
            if(be9_cz_end.lt.1.0d-24)be9_cz_end = 0.0d0
! MHP 05/02 MASS ACCRETION.  THE NET LITHIUM ADDED IN
! A GIVEN STEP IS XXXADD WHILE THE LITHIUM ADDED IN PRIOR
! STEPS IS XXXA; THE TWO ARE ADDED AFTER BURNING.
! N.B. SINCE THIS IS A HALF-LIFE PROBLEM IT IS OK TO BURN
! THE COMPONENTS SEPARATELY.
            if(accretion_active)then
               li6_added = accreted_composition(13)*substep_frac/exp(0.5d0*li6_depletion)
               li6_accreted = li6_accreted/exp(li6_depletion) + li6_added
               li7_added = accreted_composition(14)*substep_frac/exp(0.5d0*li7_depletion)
               li7_accreted = li7_accreted/exp(li7_depletion) + li7_added
               be9_added = accreted_composition(15)*substep_frac/exp(0.5d0*be9_depletion)
               be9_accreted = be9_accreted/exp(be9_depletion) + be9_added
            endif
  100    continue
! DO A MASS-WEIGHTED AVERAGE OF THE ORIGINAL LIGHT
! ELEMENT CONTENT AND THE NET AMOUNT ADDED.
! FMASSACC = DMDT*DT/ORIGINAL CZ MASS
! NOTE: THIS FORMULATION ASSUMES THAT ALL ACCRETED MATTER
! LANDS IN A CZ WITH THE END-OF-TIMESTEP DEPTH.
! NOT STRICTLY TRUE, BUT NOT A BAD APPROXIMATION EITHER.
         if(accretion_active)then
            write(*,913)li6_cz_start,li6_cz_end,li6_accreted,li7_cz_start, &
                 li7_cz_end,li7_accreted,accreted_mass_fraction
 913        format(1p7e12.3)
            li6_cz_end = (li6_cz_end*cz_mass_end+li6_accreted*accreted_mass_fraction)/ &
                 (cz_mass_end+accreted_mass_fraction)
            li7_cz_end = (li7_cz_end*cz_mass_end+li7_accreted*accreted_mass_fraction)/ &
                 (cz_mass_end+accreted_mass_fraction)
            be9_cz_end = (be9_cz_end*cz_mass_end+be9_accreted*accreted_mass_fraction)/ &
                 (cz_mass_end+accreted_mass_fraction)
         endif
! STORE THE ABUNDANCE AS A FUNCTION OF THE TIME STEP.
! STORE ABUNDANCES AS A FUNCTION OF TIMESTEP IN VECTOR YEST.
         extrap_y(1) = li6_cz_end
         extrap_y(2) = li7_cz_end
         extrap_y(3) = be9_cz_end
         extrap_x = substep_frac**2
! USE A RATIONAL POLYNOMIAL EXTRAPOLATOR TO EXTRAPOLATE TO THE SOLUTION FOR
! ZERO TIMESTEP.
         call ratext(refine_idx,extrap_x,extrap_y,extrap_result,extrap_err,3,7)
! CHECK IF THE SOLUTION HAS CONVERGED.
         if(refine_idx.gt.1)then
            converged = .true.
            do 110 species_idx=1,3
               if(extrap_y(species_idx).lt.1.0d-24)goto 110
               extrap_err(species_idx) = abs(extrap_err(species_idx)/extrap_y(species_idx))
               if(extrap_err(species_idx).gt.extrap_tol(species_idx))converged=.false.
  110       continue
            if(converged)then
!              WRITE(ISHORT,912)JENV,J,(YEXT(K2),YERR(K2),K2=1,3)
               goto 125
            endif
         endif
  120 continue
! IF THE PROGRAM GETS HERE THEN IT FAILED TO CONVERGE TO WITHIN
! THE SPECIFIED TOLERANCE IN THE MAXIMUM NUMBER OF ITERATIONS.
      write(short_file_unit,911)cz_base_zone,(extrap_err(species_idx),species_idx=1,3)
  125 continue
! WRITE NEW ABUNDANCES AND EXIT.
      li6_cz_end = extrap_result(1)
      if(li6_cz_end.lt.1.0d-24)li6_cz_end=0.0d0
      li7_cz_end = extrap_result(2)
      if(li7_cz_end.lt.1.0d-24)li7_cz_end=0.0d0
      be9_cz_end = extrap_result(3)
      if(be9_cz_end.lt.1.0d-24)be9_cz_end=0.0d0
      do 130 zone_idx = max_zone,num_zones
         composition(13,zone_idx) = li6_cz_end
         composition(14,zone_idx) = li7_cz_end
         composition(15,zone_idx) = be9_cz_end
  130 continue
! STORE ENDING RATE FOR USE AT THE BEGINNING OF THE NEXT STEP.
      log_rate_li6_prev = log_rate_li6_cz_end
      log_rate_li7_prev = log_rate_li7_cz_end
      log_rate_be9_prev = log_rate_be9_cz_end
! NOW SOLVE FOR ABUNDANCES IN THE REGION WHICH BEGAN CONVECTIVE AND
! ENDED RADIATIVE.
      if(cz_base_zone.le.cz_base_zone_old)goto 200
! FIND STARTING AND ENDING LOCATION IN MASS OF THE CZ BASE, AND
! ASSUME THAT THE FRACTION OF TIME SPENT IN THE CZ IS PROPORTIONAL TO
! THE LOCATION IN MASS (I.E. A POINT 1/3 FROM THE OLD TO THE NEW CZ
! BASE IN MASS SPENT 1/3 OF THE TIME IN THE CZ AND 2/3 OUT OF IT).
      if(cz_base_zone_old.gt.1)then
         mass_coord_beg = 0.5d0*(mass_coordinate(cz_base_zone_old)+mass_coordinate(cz_base_zone_old-1))
      else
         mass_coord_beg = 0.0d0
      endif
      mass_coord_end = 0.5d0*(mass_coordinate(cz_base_zone)+mass_coordinate(cz_base_zone-1))
      delta_mass = mass_coord_beg - mass_coord_end
      do 140 zone_idx = cz_base_zone_old,cz_base_zone-1
! MHP 9/91 CHANGE TO AVOID DIVISION BY ZERO.
! SKIP IF SHELL TEMPERATURE DROPS BELOW BURNING THRESHOLD.
         if(rate_be9(zone_idx).le.1.0d-32)goto 145
         radiative_frac = (mass_coordinate(zone_idx)-mass_coord_beg)/delta_mass
! USE FRAD*RADIATIVE RATE AND (1-FRAD)*CONVECTIVE RATE.
         li6_depletion = timestep*exp(radiative_frac*log(rate_li6(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_li6_cz_start)
         li7_depletion = timestep*exp(radiative_frac*log(rate_li7(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_li7_cz_start)
         be9_depletion = timestep*exp(radiative_frac*log(rate_be9(zone_idx))+ &
              (1.0d0-radiative_frac)*log_rate_be9_cz_start)
!***REMEMBER TO ADD FAILSAFES FOR LARGE DEPLETION***
! KC 2025-05-31 PREVENT FLOATING POINT EXCEPTION
!          HCOMP(13,I) = HCOMP(13,I)/EXP(DLI6)
         call safedivexp(composition(13,zone_idx),li6_depletion)
         if(composition(13,zone_idx).lt.1.0d-24)composition(13,zone_idx)=0.0d0
!          HCOMP(14,I) = HCOMP(14,I)/EXP(DLI7)
         call safedivexp(composition(14,zone_idx),li7_depletion)
         if(composition(14,zone_idx).lt.1.0d-24)composition(14,zone_idx)=0.0d0
         composition(15,zone_idx) = composition(15,zone_idx)/exp(be9_depletion)
         if(composition(15,zone_idx).lt.1.0d-24)composition(15,zone_idx)=0.0d0
  140 continue
  145 continue
  200 continue
      return
end subroutine liburn

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
! KC 2025-05-31 SAFEDIVEXP
!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
subroutine safedivexp(numerator, exponent)
      double precision :: numerator, exponent

      if (exponent .lt. 709.7827d0) then
         numerator = numerator/exp(exponent)
      else
         numerator = 0.0d0
      endif

      return
end subroutine safedivexp

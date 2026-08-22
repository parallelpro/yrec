!----------------------------------------------------------------------
! getw
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getw.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 11/91 MHP revised subroutine (replaces old FINDW).
! GETW evolves the angular momentum distribution in time.
! It enforces local conservation of angular momentum in radiative regions,
! and enforces the user-specified rotation law in convective regions.
! Both meridional circulation and differential rotation with depth can
! cause the radial transport of angular momentum and material.
! GETW solves for the transport of angular momentum and the associated mixing
! if desired (LINSTB=T), and accounts for angular momentum loss caused by a
! magnetic stellar wind if desired (LJDOT0=T).
! If instabilities are being treated, the burning of the light elements
! lithium and beryllium are considered here if desired (LEXCOM=T).
! For a discussion of the overall technique used see Pinsonneault, Kawaler,
! Sofia, and Demarque (1989), Ap.J. vol. 338, p.424.
! 11/91 JENV0 added to call.
!
! NOTE: this file is only reached (from main.f) inside an IF(LROT)
! gate; within this file itself, whether SECULR (the rotation-
! diffusion instability solver) is called is gated by
! instability_transport_active (LINSTB), preserved exactly below.
subroutine getw(full_timestep, max_domega_step, wind_loss_active, &
     envelope_boundary_zone_prev, ierr)
      use star_info_lib, only: star
      use star_info_lib, only: star
      use net_lib
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use star_info_lib, only: star
      use luout_lib
      use const_lib
      use burn_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: full_timestep
      double precision, intent(inout) :: max_domega_step
      logical, intent(inout) :: wind_loss_active
      integer, intent(inout) :: envelope_boundary_zone_prev

! am_transport_convective_flag (originally LCZ): convective-for-AM-
! transport-purposes flag, set by OVROT below. Naming matches
! hpoint.f90/midmod.f90.
      logical :: am_transport_convective_flag(json)
      double precision :: specific_angular_momentum_saved(json)
! "_mid" arrays are the midpoint-in-time structure computed by MIDMOD
! each diffusion sub-step; names match midmod.f90's own dummy-argument
! names for these exact quantities.
      double precision :: eta_squared_mid(json), log_density_mid(json), &
           hg_mid(json), moment_of_inertia_mid(json), log_luminosity_mid(json), &
           log_pressure_mid(json), log_radius_mid(json), log_temperature_mid(json), &
           omega_mid(json), mean_radius_mid(json), qiw_mid(json)
      logical :: convective_flag_mid(json), am_transport_convective_flag_mid(json)
! MHP 6/00 added cod2, diffusion_velocity (HV) to allow BUR-ST mixing plus burning.
! Both are produced by SECULR (not yet converted) and consumed by
! BURSMIX; their exact roles are not otherwise exercised in this file.
      double precision :: cod2(json), diffusion_velocity(json)
      integer :: convective_zone_bounds(12,2), &
           convective_zone_bounds_burn(12,2), radiative_zone_bounds(13,2)
      double precision :: cz_mass_bottom, cz_mass_top
! LBURS is hardcoded .FALSE. here (preserved exactly, per project
! notes -- this permanently disables the BUR-ST extrapolation branch
! near the end of this routine regardless of LINSTB/LALLCZ).
      logical :: burs_extrapolation_active
      integer :: num_species_tracked
      logical :: disk_lock_engaged
      double precision :: omega_surface
      double precision :: moment_of_inertia_cz
! envelope_boundary_zone/core_boundary_zone (originally IMAX/IMIN) are
! reused across both the early (LINSTB=.false.) return branch, where
! envelope_boundary_zone is set from convective_zone_bounds directly,
! and the main diffusion loop, where MIDMOD sets both from the
! midpoint-in-time structure -- the two branches never execute in the
! same call (the early branch exits via GOTO 9999), so this matches
! midmod.f90's own naming for the same physical quantities throughout.
      integer :: envelope_boundary_zone, core_boundary_zone
      integer :: zone_index, species_index
      integer :: redo_count
      logical :: skip_diffusion_flag
      integer :: num_diffusion_steps, num_wind_diffusion_steps
      double precision :: sub_timestep
      double precision :: elapsed_substep_time
      logical :: first_call
      double precision :: delts_test
      double precision :: fx
      logical :: surface_cz_active, fully_convective_flag
      integer :: num_radiative_zones, num_convective_zones
      logical :: redo_needed_flag
      logical :: diffusion_solve_ok
      integer :: core_boundary_zone_cur, envelope_boundary_zone_cur
      integer :: num_convective_zones_burn
! mix_grads_flag (originally LIMIX) is set ("G Somers 6/14, SET LIMIX
! = .FALSE. SO THE CORRECT GRADS ARE USED") but the active MIXCZ call
! below does not actually take it as an argument -- vestigial, kept
! exactly as in the original.
      logical :: mix_grads_flag
      integer :: iend
      double precision :: surface_quad_term, surface_potential
      double precision :: log_radius_surface
      double precision :: omega_avg, domega_dr, delta_radius_step
      double precision :: radius_at_bcz

      integer, intent(out) :: ierr

      ierr = 0

      burs_extrapolation_active = .false.
! DETERMINE THE NUMBER OF DIFFERENT ELEMENTS AND ISOTOPES BEING TRACKED
! BY THE CODE (NSPEC).
      num_species_tracked = 11
      if(use_extended_composition) num_species_tracked = 15
!  CONVECTIVE AND RADIATIVE REGIONS ARE TREATED DIFFERENTLY FOR ANGULAR
!  MOMENTUM PURPOSES; IF CONVECTIVE OVERSHOOT IS BEING INCLUDED THE
!  OVERSHOOT REGION IS TREATED AS CONVECTIVE FOR ANGULAR MOMENTUM PURPOSES.
!  ***NOTE : IF LINSTB=T THERE MUST BE AT LEAST ONE SHELL IN THE OVERSHOOT
!  REGION, TO COUPLE THE RADIATIVE AND CONVECTIVE REGIONS AND TO AVOID
!  NUMERICAL PROBLEMS IN SOME OF THE TERMS WHICH ARISE IN THE SOLUTION OF
!  THE DIFFUSION EQUATIONS.
!
!  OVROT CALCULATES CONVECTIVE OVERSHOOT BY ARBITRARILY MIXING TO A
!  USER-SPECIFIED FRACTION OF A PRESSURE SCALE HEIGHT ABOVE OR BELOW A
!  CONVECTIVE REGION.  IT RETURNS A VECTOR LCZ, WHICH IS TRUE IF A SHELL IS
!  CONVECTIVE FOR ANGULAR MOMENTUM PURPOSES.
! 7/91 CHANGED CALL TO OVROT.
      call ovrot(star%composition,star%log_density,star%log_pressure,star%log_radius,star%log_mass, &
           star%log_temperature,star%convective_flag,star%num_zones, &
           am_transport_convective_flag,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones)
! MHP 9/94 ADDED DISK LOCKING OPTION.
      disk_lock_engaged = .false.
      if(disk_locking_active .and. disk_lifetime.le.disk_temperature) &
           disk_lock_engaged = .true.
      if(.not.instability_transport_active)then
!  STORE THE SURFACE ANGULAR VELOCITY FROM THE BEGINNING OF THE TIMESTEP.
         omega_surface = star%omega(star%num_zones)
!  ENFORCE SB ROTATION (OR UNIFORM ROTATION LAW IN ENTIRE STAR) IF DESIRED.
! JNT 2025/09 FOR 05/15 MAKE IMPJMOD=1 ACT LIKE LSOLID
         if(force_solid_body_rotation .or. (solid_body_mode_flag.eq.1))then
            do zone_index = 1,star%num_zones
               am_transport_convective_flag(zone_index) = .true.
            end do
         endif
!  GETROT TAKES THE ANGULAR MOMENTUM DISTRIBUTION AND FINDS THE
!  ROTATION CURVE THAT CORRESPONDS TO IT.  CONVECTIVE REGIONS HAVE
!  SOLID BODY ROTATION ENFORCED ON THEM.
         call getrot(star%log_density,star%specific_angular_momentum,star%log_radius, &
              star%log_mass,star%shell_mass,am_transport_convective_flag,star%num_zones, &
              star%eta_squared,star%moment_of_inertia,star%omega,star%qiw,star%mean_radius)
!  ANGULAR MOMENTUM LOSS WITHOUT INTERNAL ANGULAR MOMENTUM TRANSPORT.
!         DO I = 1,M
!            WOLD(I) = OMEGA(I)
!         END DO
         if(.not.disk_lock_engaged .and. ljdot0 .and. star%convective_flag(star%num_zones)) then
!  FIND MOMENT OF INERTIA OF THE SURFACE C.Z.
            moment_of_inertia_cz = 0.0D0
!  ENFORCE SB ROTATION (OR UNIFORM ROTATION LAW IN ENTIRE STAR) IF DESIRED.
            if(force_solid_body_rotation.or. (solid_body_mode_flag.eq.1))then
              envelope_boundary_zone = 1
            else
              envelope_boundary_zone = convective_zone_bounds(num_convective_zones,1)
            endif
            do zone_index = envelope_boundary_zone,star%num_zones
               moment_of_inertia_cz = moment_of_inertia_cz + star%moment_of_inertia(zone_index)
            end do
!  FIND LOWEST SHELL IN SURFACE CZ (IMAX)
!            IMAX = MXZONE(NZONE,1)
!  HSTOP IS THE MASS AT THE TOP OF THE C.Z.
!  HSBOT IS THE MASS AT THE BOTTOM OF THE C.Z.
            cz_mass_top = exp(ln10*star%log_total_mass)
            if(envelope_boundary_zone.gt.1)then
               cz_mass_bottom = 0.5D0*(star%enclosed_mass(envelope_boundary_zone)+ &
                    star%enclosed_mass(envelope_boundary_zone-1))
            else
               cz_mass_bottom = 0.0D0
            endif
            wind_loss_active = ljdot0
! MHP 10/02 UNUSED LFIRST REMOVED FROM CALL
! MHP 10/17 timestep average loss rate
!            FRACSTEP = 1.
            star%turnover%fracstep = 0.5
            call mwind(star%log10_luminosity,full_timestep,cz_mass_bottom, &
                 cz_mass_top,envelope_boundary_zone,star%num_zones,wind_loss_active, &
                 omega_surface, &
                 star%total_mass_msun,star%log_teff,moment_of_inertia_cz, &
                 star%specific_angular_momentum, ierr)
            if (ierr /= 0) return
!  FIND THE NEW RUN OF OMEGA IN THE CONVECTION ZONE AFTER THE WIND.
! MHP 10/02 REPLACED IEND WITH M IN CALL, DEFINED IEND=M
            iend = star%num_zones
! JNT 09/2025 FOR 05/15 REPLACE WCZ WITH WCZIMP
            call wczimp(star%log_density,star%specific_angular_momentum,star%log_radius, &
                 star%log_mass,star%shell_mass,envelope_boundary_zone,iend,star%eta_squared, &
                 star%moment_of_inertia,star%omega,star%qiw,star%mean_radius,star%num_zones)
         endif
      else
!  GETROT TAKES THE ANGULAR MOMENTUM DISTRIBUTION AND FINDS THE
!  ROTATION CURVE THAT CORRESPONDS TO IT.  CONVECTIVE REGIONS HAVE
!  SOLID BODY ROTATION ENFORCED ON THEM.
      call getrot(star%log_density,star%specific_angular_momentum,star%log_radius,star%log_mass, &
           star%shell_mass,am_transport_convective_flag,star%num_zones,star%eta_squared, &
           star%moment_of_inertia,star%omega,star%qiw,star%mean_radius)
      skip_diffusion_flag = .not.instability_transport_active .and. .not.wind_loss_active
      if (.not. (skip_diffusion_flag.or.full_timestep.le.0.0D0)) then
!  NOW LIMIT THE DIFFUSION TIMESTEP TO A MAXIMUM CHANGE IN OMEGA
!  FROM THE PREVIOUS MODEL.
      if(max_domega_step.eq.0.0D0) max_domega_step = max_domega_global
      num_diffusion_steps = int(max_domega_step/dtdif)
      if (mod(max_domega_step,dtdif).ne.0.0D0) num_diffusion_steps = num_diffusion_steps + 1
      num_wind_diffusion_steps = int(max_domega_global/dtdif)
      if (mod(max_domega_global,dtdif).ne.0.0D0) num_wind_diffusion_steps = num_wind_diffusion_steps + 1
!     NSTEP = MAX(NSTEP,NSTEP2/2)
      num_diffusion_steps = min(num_diffusion_steps,num_wind_diffusion_steps)
      sub_timestep = full_timestep/dfloat(num_diffusion_steps)
!  FIND BASIC PHYSICAL QUANTITIES NEEDED FOR BOTH SECULAR AND DYNAMICAL
!  INSTABILITES: ADIABATIC AND ACTUAL TEMPERATURE GRADIENTS,OPACITIES,
!  KINEMATIC VISCOSITIES,THERMOMETRIC DIFFUSIVITY, AND HEAT CAPACITY.
!       CALL PHYSIC(FP,FT,HCOMP,HD,HG,HL,HP,HR,HS,HT,LC,LCZ,M,TEFFL)  ! KC 2025-05-31
      call physic(star%pressure_rotation_factor,star%temperature_rotation_factor,star%composition,star%log_density,star%mean_gravity,star%luminosity_lsun,star%log_pressure, &
           star%log_radius,star%log_mass,star%log_temperature,star%convective_flag,star%num_zones,star%log_teff, ierr)
      if (ierr /= 0) return
! 8/17 DETERMINE TAUCZ AND PPHOT
!       CALL GETTAU(HCOMP,HR,HP,HD,HG,HS1,HT,FP,FT,TEFFL,  ! KC 2025-05-31
      call gettau(star%composition,star%log_radius,star%log_pressure,star%log_density,star%enclosed_mass, &
           star%log_temperature,star%pressure_rotation_factor,star%temperature_rotation_factor,star%log_teff, &
           star%log_total_mass,star%log10_luminosity,star%num_zones,star%convective_flag,radius_at_bcz)
! IF DT IS LESS THAN DELTS, THEN THE MODEL TIMESTEP IS TOO LONG FOR THE
! DIFFUSION CALCULATIONS.  IF THIS OCCURS,
! USE A SERIES OF SMALLER TIMESTEPS THAT DON'T VIOLATE THIS CONDITION.
! INTERPOLATE THE MODEL STRUCTURE BETWEEN THE OLD MODEL AND THE NEW ONE
! TO CALCULATE DIFFUSION AND WIND BETWEEN MODELS.
      elapsed_substep_time = 0.0D0
      first_call = .true.
!      DO 20 I = 1,M
!         WOLD(I) = OMEGA(I)
!   20 CONTINUE
!  ENTRY FOR SERIES OF DIFFUSION TIMESTEPS.
      do
         continue
      redo_count = 0
!  ENTRY FOR DIFFUSION TIMESTEP CUTTING.
      retry: do
! MHP 06/02 CHANGED TO ELIMINATE OCCASIONAL
! ALMOST ZERO TIMESTEP FROM ROUNDOFF ERROR
      delts_test = full_timestep - full_timestep*1.0D-6
      if(elapsed_substep_time+sub_timestep.ge.delts_test)then
         sub_timestep = full_timestep - elapsed_substep_time
         elapsed_substep_time = full_timestep
      else
         elapsed_substep_time = elapsed_substep_time + sub_timestep
      endif
! COPY OVER PRIOR THETA(TIME) TERM TO TEMPORARY SLOT
! FOR USE IN THE ADVECTION/DIFFUSION TREATMENT OF MAEDER&ZAHN 1998
      if(first_call)then
         star%rot%theta_prev(1) = star%rot%tho(1)
         star%rot%wmst(1) = star%run%old_omega(1)
         do zone_index = 2,star%num_zones
            star%rot%qwrmst(zone_index) = star%rot%qwrst(zone_index)
            star%rot%theta_prev(zone_index) = star%rot%tho(zone_index)
            star%rot%wmst(zone_index) = star%run%old_omega(zone_index)
         end do
! RECOMPUTE THETA
      else
         star%rot%wmst(1) = star%omega(1)
         do zone_index = 2,star%num_zones
            omega_avg = 0.5D0*(star%omega(zone_index)+star%omega(zone_index-1))
            delta_radius_step = 10.0D0**log_radius_mid(zone_index)- &
                 10.0D0**log_radius_mid(zone_index-1)
            domega_dr = (star%omega(zone_index)-star%omega(zone_index-1))/delta_radius_step
            star%rot%theta_prev(zone_index) = star%rot%theta_mean(zone_index)*omega_avg*domega_dr
            star%rot%qwrmst(zone_index) = domega_dr
            star%rot%wmst(zone_index) = star%omega(zone_index)
         end do
      endif
      fx = elapsed_substep_time/full_timestep
      star%turnover%fracstep = fx
! INTERPOLATE LINEARLY IN TIME FOR THE MODEL STRUCTURE BETWEEN THE
! START AND END OF THE TIMESTEP.
! JVS
      call midmod(full_timestep, sub_timestep, fx, first_call, &
           eta_squared_mid, log_density_mid, hg_mid, &
           moment_of_inertia_cz, moment_of_inertia_mid, &
           log_luminosity_mid, log_pressure_mid, log_radius_mid, &
           cz_mass_bottom, cz_mass_top, log_temperature_mid, &
           core_boundary_zone, envelope_boundary_zone, &
           fully_convective_flag, convective_flag_mid, &
           am_transport_convective_flag_mid, surface_cz_active, &
           omega_mid, mean_radius_mid, qiw_mid, radiative_zone_bounds, &
           convective_zone_bounds, num_radiative_zones, &
           num_convective_zones, ierr)
      if (ierr /= 0) return
! IF DESIRED, REMOVE ANGULAR MOMENTUM FROM OUTER CONVECTION ZONE
! USING A WEBER-DAVIS MAGNETIC WIND MODEL
      if(.not.disk_lock_engaged .and. wind_loss_active .and. surface_cz_active) then
         omega_surface = omega_mid(star%num_zones)
!  WIND CALCULATION WITHOUT INSTABILITIES
         if(fully_convective_flag.or..not.instability_transport_active) then
! MHP 10/02 UNUSED LFIRST REMOVED FROM CALL
            call mwind(star%log10_luminosity,sub_timestep,cz_mass_bottom,cz_mass_top, &
                 envelope_boundary_zone,star%num_zones,wind_loss_active,omega_surface, &
                 star%total_mass_msun,star%log_teff,moment_of_inertia_cz, &
                 star%specific_angular_momentum, ierr)
            if (ierr /= 0) return
!  FIND THE NEW OMEGA OF THE CONVECTION ZONE AFTER THE WIND.
! MHP 10/02 REPLACED IEND WITH M IN CALL, DEFINED IEND=M
            iend = star%num_zones
! JNT 09/2025 FOR 05/15 REPLACE WCZ WITH WCZIMP
            call wczimp(log_density_mid,star%specific_angular_momentum,log_radius_mid, &
                 star%log_mass,star%shell_mass,envelope_boundary_zone,iend,eta_squared_mid, &
                 moment_of_inertia_mid,omega_mid,qiw_mid,mean_radius_mid,star%num_zones)
         endif
      endif
!  NOW CHECK FOR INSTABILITIES IN RADIATIVE REGIONS
      if(instability_transport_active.and..not.fully_convective_flag) then
!  ENSURE THAT CONVECTIVE REGIONS ARE FULLY MIXED.
!  JVS 0212       CALL MIXCZ(HCOMP,HS2,LCZM,M)
! KC 2025-05-30 addressed warning messages from Makefile.legacy
! C G Somers 6/14, SET IMIX = .FALSE. SO THE CORRECT GRADS ARE USED.
!         IMIX = .FALSE.
!         CALL MIXCZ(HCOMP,HS2,HS1,LCZM,HRM,HPM,HDM,HGM,M,IMIX)
! G Somers 6/14, SET LIMIX = .FALSE. SO THE CORRECT GRADS ARE USED.
          mix_grads_flag = .false.
!         CALL MIXCZ(HCOMP,HS2,HS1,LCZM,HRM,HPM,HDM,HGM,M,LIMIX)  ! KC 2025-05-31
          call mixcz(star%composition,star%shell_mass,am_transport_convective_flag_mid,star%num_zones)
! G Somers END
!  NOW SOLVE FOR LONG-TIMESCALE(SECULAR) INSTABILITIES.
!  THESE ARE TREATED USING DIFFUSION EQUATIONS.
         do zone_index = 1,star%num_zones
            specific_angular_momentum_saved(zone_index) = star%specific_angular_momentum(zone_index)
   60    continue
         end do
! MHP 6/00 ADDED COD2,HV TO LIST RETURNED FROM SECULR
! FOR THE BUR-ST MIXING ROUTINES
!          CALL SECULR(DELTS,DT,HDM,HGM,HIM,HLM,HPM,HRM,HS,HS1,HS2,
!      *               HTM,HJMSAV,LCZM,LCM,M,MODEL,OMEGAM,
!      *               HJM,ETA2M,QIWM,R0M,HCOMP,LFIRST,IMIN,IMAX,BL,
!      *               HSTOT,SJTOT,SMASS,TEFFL,LREDO,IREDO,  ! KC 2025-05-31
         call seculr(sub_timestep,log_density_mid,hg_mid,moment_of_inertia_mid, &
              log_luminosity_mid,log_pressure_mid,log_radius_mid,star%log_mass, &
              star%enclosed_mass,star%shell_mass, &
              log_temperature_mid,specific_angular_momentum_saved, &
              am_transport_convective_flag_mid,star%num_zones,omega_mid, &
              star%specific_angular_momentum,eta_squared_mid,qiw_mid,mean_radius_mid, &
              star%composition,core_boundary_zone,envelope_boundary_zone, &
              star%log10_luminosity, &
              star%log_total_mass,star%total_mass_msun,star%log_teff,redo_needed_flag,redo_count, &
              moment_of_inertia_cz,cz_mass_bottom,cz_mass_top,omega_surface,surface_cz_active, &
!      *               MRZONE,MXZONE,NRZONE,NZONE,  ! KC 2025-05-31
              cod2,diffusion_velocity,diffusion_solve_ok, ierr)
         if (ierr /= 0) return
!  DIFFUSION TIMESTEP CUTTING REQUIRED IF LREDO IS TRUE.
!  RESET MEAN MOLECULAR WEIGHT,COMPOSITION,AND SPECIFIC ANGULAR MOMENTUM
!  TO THE VALUES THEY HAD PRIOR TO THE START OF THE STEP.
         if(redo_needed_flag)then
            fx = 2.0D0*sub_timestep/full_timestep
            elapsed_substep_time = elapsed_substep_time - 2.0D0*sub_timestep
            do zone_index = 1,star%num_zones
               star%specific_angular_momentum(zone_index) = specific_angular_momentum_saved(zone_index)
               star%mix_phys%amum(zone_index) = star%mix_phys%amum(zone_index) - fx*(star%thermo%mean_molecular_weight(zone_index)-star%rot%old_amu(zone_index))
               do species_index = 1,num_species_tracked
                  star%composition(species_index,zone_index) = star%prev%old_composition(species_index,zone_index)
   70          continue
               end do
   80       continue
            end do
            cycle retry
         endif
      endif
! MHP 05/02 TAKEN OUT OF LOOP SO THAT NUCLEAR BURNING
! IS ACCOUNTED FOR IN MODELS THAT ARE FULLY CONVECTIVE
! OR WITHOUT INSTABILITIES
! ADDED CHANGE FOR BURLICH-STORER TREATMENT OF MIXING PLUS
! BURNING - ONLY UPDATED IF NOT USED
      if(.not.burs_extrapolation_active)then
         do zone_index = 1,star%num_zones
            do species_index = 1,11
                  star%prev%old_composition(species_index,zone_index) = star%composition(species_index,zone_index)
            end do
         end do
      endif
      first_call = .false.
!  CLEAN UP EXTENDED COMP ARRAYS; ZERO VALUES THROUGHOUT MOST OF THE
!  INTERIOR CAN LEAD TO SMALL NEGATIVE VALUES DUE TO ROUNDOFF ERROR IN
!  THE DIFFUSION CALCULATIONS.
      if(use_extended_composition)then
!  PERFORM LIGHT ELEMENT BURNING.
! FIND SURFACE C.Z. DEPTH AT THE END OF THE TIME STEP.
         call convec(star%composition,log_density_mid,log_pressure_mid,log_radius_mid, &
              star%log_mass,log_temperature_mid,convective_flag_mid,star%num_zones, &
              radiative_zone_bounds,convective_zone_bounds, &
              convective_zone_bounds_burn,core_boundary_zone_cur, &
              envelope_boundary_zone_cur,num_radiative_zones,num_convective_zones, &
              num_convective_zones_burn)
! 11/91 CHANGED FOR LITHIUM BURNING WITH OVERSHOOT.
         if(lovstm .and. convective_flag_mid(star%num_zones))then
            star%light_burn%pressure_scale_height_end = alphae*exp(clndp*(log_pressure_mid(envelope_boundary_zone_cur)+ &
                 2.0D0*log_radius_mid(envelope_boundary_zone_cur) &
                 -log_density_mid(envelope_boundary_zone_cur)-cgl-star%log_mass(envelope_boundary_zone_cur)))
         else
            star%light_burn%pressure_scale_height_end = 0.0D0
         endif
! FIND LIGHT ELEMENT BURNING RATES AT THE END OF THE TIME STEP.
         call lirate88(star%composition,log_density_mid,log_temperature_mid,star%num_zones,2)
! STORE CURRENT "END OF STEP" QUANTITIES AS "BEGINNING" ONES FOR
! THE NEXT STEP.
! ADDED CHANGE FOR BURLICH-STORER TREATMENT OF MIXING PLUS
! BURNING - ONLY UPDATED IF NOT USED
         if(.not.instability_transport_active .or. .not.burs_extrapolation_active .or. fully_convective_flag)then
! COMPUTE BURNING.
!             CALL LIBURN(DT,HCOMP,HDM,HRM,HS1,HS2,HTM,JENV1,JENV0,M)  ! KC 2025-05-31
            call liburn(sub_timestep,star%composition,log_radius_mid,star%enclosed_mass, &
                 star%shell_mass,log_temperature_mid,envelope_boundary_zone_cur, &
                 envelope_boundary_zone_prev,star%num_zones)
            envelope_boundary_zone_prev = envelope_boundary_zone_cur
            star%light_burn%pressure_scale_height_start = star%light_burn%pressure_scale_height_end
            do zone_index = 1,star%num_zones
               star%light_burn%rate_li6_start(zone_index) = star%light_burn%rate_li6(zone_index)
               star%light_burn%rate_li7_start(zone_index) = star%light_burn%rate_li7(zone_index)
               star%light_burn%rate_be9_start(zone_index) = star%light_burn%rate_be9(zone_index)
  155       continue
            end do
         else
! COMPUTE BURNING.
!             CALL LIBURN2(DT,HCOMP,HDM,HRM,HS1,HS2,HTM,JENV1,JENV0,M)  ! KC 2025-05-31
            call liburn2(sub_timestep,star%composition,log_radius_mid,star%enclosed_mass, &
                 star%shell_mass,log_temperature_mid,envelope_boundary_zone_cur, &
                 envelope_boundary_zone_prev,star%num_zones)
         endif
!  ENSURE THAT CONVECTIVE REGIONS ARE FULLY MIXED.
!  JVS 0212       CALL MIXCZ(HCOMP,HS2,LCZM,M)
! KC 2025-05-30 addressed warning messages from Makefile.legacy
! C G Somers 6/14, SET IMIX = .FALSE. SO THE CORRECT GRADS ARE USED.
!         IMIX = .FALSE.
!         CALL MIXCZ(HCOMP,HS2,HS1,LCZM,HRM,HPM,HDM,HGM,M,IMIX)
! G Somers 6/14, SET LIMIX = .FALSE. SO THE CORRECT GRADS ARE USED.
          mix_grads_flag = .false.
!         CALL MIXCZ(HCOMP,HS2,HS1,LCZM,HRM,HPM,HDM,HGM,M,LIMIX)  ! KC 2025-05-31
          call mixcz(star%composition,star%shell_mass,am_transport_convective_flag_mid,star%num_zones)
! G Somers END
!  ZERO OUT LOW ABUNDANCES.
         do zone_index = 1,star%num_zones
            do species_index = 12,15
               if(star%composition(species_index,zone_index).lt.1.0D-24)star%composition(species_index,zone_index)=0.0D0
   90       continue
            end do
  100   continue
         end do
! MHP 6/00 ADDED OVERWRITE OF HCOMPP FOR LIGHT ELEMENTS
! ADDED CHANGE FOR BURLICH-STORER TREATMENT OF MIXING PLUS
! BURNING - ONLY UPDATED IF NOT USED
         if(.not.instability_transport_active.or.fully_convective_flag)then
            do zone_index = 1,star%num_zones
               do species_index = 12,15
                  star%prev%old_composition(species_index,zone_index) = star%composition(species_index,zone_index)
               end do
            end do
         else if(.not.burs_extrapolation_active)then
            do zone_index = 1,star%num_zones
               do species_index = 12,15
                  star%prev%old_composition(species_index,zone_index) = star%composition(species_index,zone_index)
               end do
            end do
         endif
      endif
! MHP 6/00 NOW ADDED THE OPTION OF PERFORMING A BUR-ST EXTRAPOLATION
! TO ZERO TIMESTEP FOR THE COMBINATION OF MIXING AND NUCLEAR BURNING.
      if(instability_transport_active .and. .not.fully_convective_flag)then
         if(burs_extrapolation_active)then
        call bursmix(cod2,sub_timestep,star%composition,log_density_mid,log_luminosity_mid, &
             log_pressure_mid,log_radius_mid,star%log_mass,star%enclosed_mass,star%shell_mass,star%log_total_mass, &
             log_temperature_mid,diffusion_velocity,envelope_boundary_zone,core_boundary_zone, &
             envelope_boundary_zone_prev,envelope_boundary_zone_cur,diffusion_solve_ok, &
             am_transport_convective_flag_mid,star%num_zones,radiative_zone_bounds, &
             convective_zone_bounds,num_radiative_zones,num_convective_zones, ierr)
        if (ierr /= 0) return
         endif
      endif
      exit retry
      end do retry
!  RETURN FOR NEXT SMALL DIFFUSION TIMESTEP IF NEEDED.
      if (.not. (elapsed_substep_time.lt.full_timestep)) exit
      end do
!  UPDATE OMEGA ARRAY TO REFLECT NEW ANGULAR MOMENTUM DISTRIBUTION.
      do zone_index = 1,star%num_zones
         star%omega(zone_index) = omega_mid(zone_index)
  110 continue
      end do
! MHP 3/96 ADDED CALL TO RECOMPUTE SELF-CONSISTENT SET OF OMEGAS
      call getrot(star%log_density,star%specific_angular_momentum,star%log_radius,star%log_mass, &
           star%shell_mass,am_transport_convective_flag,star%num_zones,star%eta_squared, &
           star%moment_of_inertia,star%omega,star%qiw,star%mean_radius)
      end if
      endif
      if(star%run%lprt0_placeholder)then
         star%log10_luminosity = log10(star%luminosity_lsun(star%num_zones))
         log_radius_surface = 0.5D0*(star%log10_luminosity + log10_solar_luminosity &
              - 4.0D0*star%log_teff - c4pil - csigl)
         fx = exp(ln10*3.0D0*(star%log_radius(star%num_zones)-log_radius_surface))
         surface_quad_term = fx*star%rot%quadrupole_moment(star%num_zones)
         surface_potential = exp(ln10*(cgl+star%log_total_mass-log_radius_surface))
         write(*,9911)surface_quad_term,surface_potential,-1.5D0*surface_quad_term/surface_potential
 9911    format(1X,'QUAD ',1PE12.3,' PHIS ',E12.3,' 3/2 QUAD/G ', &
     E12.3)
      endif
      return
end subroutine getw

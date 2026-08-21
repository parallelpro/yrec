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
subroutine getw(log_luminosity_lsun, full_timestep, max_domega_step, fp, ft, &
     composition, log_density, specific_angular_momentum, log_luminosity, &
     log_pressure, log_radius, log_mass, &
     enclosed_mass, shell_mass, log_total_mass, log_temperature, &
     convective_flag, wind_loss_active, num_zones, total_mass_msun, &
     log_teff, eta_squared, hg, moment_of_inertia, omega, qiw, mean_radius, &
     envelope_boundary_zone_prev)
      use temp2_lib
      use temp_lib
      use mdphy_lib
      use light_burn_lib
      use turnover_lib
      use oldmod_lib
      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(inout) :: log_luminosity_lsun
      double precision, intent(in) :: full_timestep
      double precision, intent(inout) :: max_domega_step
      double precision, intent(inout) :: fp(json), ft(json)
      double precision, intent(inout) :: composition(15,json)
      double precision, intent(inout) :: log_density(json), &
           specific_angular_momentum(json)
      double precision, intent(inout) :: log_luminosity(json), &
           log_pressure(json), log_radius(json)
      double precision, intent(in) :: log_mass(json)
      double precision, intent(in) :: enclosed_mass(json), shell_mass(json)
      double precision, intent(in) :: log_total_mass
      double precision, intent(inout) :: log_temperature(json)
      logical, intent(in) :: convective_flag(json)
      logical, intent(inout) :: wind_loss_active
      integer, intent(in) :: num_zones
      double precision, intent(in) :: total_mass_msun, log_teff
      double precision, intent(inout) :: eta_squared(json), hg(json), &
           moment_of_inertia(json), omega(json), qiw(json), mean_radius(json)
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









! MHP 8/17 added excen, c_2 to common block for Matt et al. 2012 cent. term
! common/cwind/: only ljdot0 is used here. Naming matches midmod.f90.
      double precision :: wind_saturation_omega, exmd, &
           wind_law_omega_exponent, extau, exr, exm, exl, expr, &
           constfactor, structfactor, excen, c_2
      logical :: ljdot0
      common/cwind/ wind_saturation_omega, exmd, wind_law_omega_exponent, &
           extau, exr, exm, exl, expr, constfactor, structfactor, excen, &
           c_2, ljdot0





! common/oldphy/: only old_amu is used here. Naming matches hpoint.f90.
      double precision :: old_delm(json), old_del_adiabatic_mix(json), &
           old_amu(json), old_om(json), old_cp(json), old_qdt(json), &
           old_vel(json), old_visc(json), old_thdif(json), old_esum(json), &
           old_del_radiative_mix(json), old_eps(json)
      common/oldphy/ old_delm, old_del_adiabatic_mix, old_amu, old_om, &
           old_cp, old_qdt, old_vel, old_visc, old_thdif, old_esum, &
           old_del_radiative_mix, old_eps


! common/oldrot/: only old_omega (WOLD) is used here. Naming matches
! hpoint.f90.
      double precision :: old_omega(json), old_specific_angular_momentum(json), &
           old_moment_of_inertia(json), old_hg(json), old_mean_radius(json), &
           old_eta_squared(json)
      common/oldrot/ old_omega, old_specific_angular_momentum, &
           old_moment_of_inertia, old_hg, old_mean_radius, old_eta_squared





! common/oldab/: not used in this file. Naming matches midmod.f90.
      double precision :: composition_snapshot(15,json)
      common/oldab/ composition_snapshot


! common/vmult/: not used in this file. Naming matches vcirc.f90.
      double precision :: difad_velocity_scale, mixing_velocity_scale, fo, &
           es_velocity_scale, gsf_velocity_scale, mu_gradient_scale, &
           secular_shear_velocity_scale, critical_reynolds
      common/vmult/ difad_velocity_scale, mixing_velocity_scale, fo, &
           es_velocity_scale, gsf_velocity_scale, mu_gradient_scale, &
           secular_shear_velocity_scale, critical_reynolds

! common/quadru/: only quadrupole_moment is used here. Naming matches
! vcirc.f90.
      double precision :: quadrupole_moment(json), local_gravity(json)
      common/quadru/ quadrupole_moment, local_gravity

! common/rotprt/: lprt0_placeholder is actually used here (gating the
! final diagnostic print), despite the "_placeholder" name inherited
! from wrtout.f90 (where this block is unused). Naming matches
! wrtout.f90.
      logical :: lprt0_placeholder
      common/rotprt/ lprt0_placeholder


! MHP 3/99 ADDED FLAG TO TREAT THE ENTIRE STAR AS 'CONVECTIVE'
! for angular momentum purposes.
! JNT 09/2025 for 05/15 add impjmod
! common/sbrot/: both used here. Naming matches midmod.f90.
      logical :: force_solid_body_rotation
      integer :: solid_body_mode_flag
      common/sbrot/ force_solid_body_rotation, solid_body_mode_flag

! MHP 06/02
! Time change of theta
! common/oldrot2/: tho/theta_mean/theta_prev/qwrst/wmst/qwrmst are
! used here. Naming matches hpoint.f90.
      double precision :: tho(json), theta_new(json), theta_mean(json), &
           del_grad_diff_interface(json), es_relaxation_factor(json), &
           theta_prev(json), qwrst(json), wmst(json), qwrmst(json)
      common/oldrot2/ tho, theta_new, theta_mean, del_grad_diff_interface, &
           es_relaxation_factor, theta_prev, qwrst, wmst, qwrmst


      save

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
      call ovrot(composition,log_density,log_pressure,log_radius,log_mass, &
           log_temperature,convective_flag,num_zones, &
           am_transport_convective_flag,radiative_zone_bounds, &
           convective_zone_bounds,num_radiative_zones,num_convective_zones)
! MHP 9/94 ADDED DISK LOCKING OPTION.
      disk_lock_engaged = .false.
      if(disk_locking_active .and. disk_lifetime.le.disk_temperature) &
           disk_lock_engaged = .true.
      if(.not.instability_transport_active)then
!  STORE THE SURFACE ANGULAR VELOCITY FROM THE BEGINNING OF THE TIMESTEP.
         omega_surface = omega(num_zones)
!  ENFORCE SB ROTATION (OR UNIFORM ROTATION LAW IN ENTIRE STAR) IF DESIRED.
! JNT 2025/09 FOR 05/15 MAKE IMPJMOD=1 ACT LIKE LSOLID
         if(force_solid_body_rotation .or. (solid_body_mode_flag.eq.1))then
            do zone_index = 1,num_zones
               am_transport_convective_flag(zone_index) = .true.
            end do
         endif
!  GETROT TAKES THE ANGULAR MOMENTUM DISTRIBUTION AND FINDS THE
!  ROTATION CURVE THAT CORRESPONDS TO IT.  CONVECTIVE REGIONS HAVE
!  SOLID BODY ROTATION ENFORCED ON THEM.
         call getrot(log_density,specific_angular_momentum,log_radius, &
              log_mass,shell_mass,am_transport_convective_flag,num_zones, &
              eta_squared,moment_of_inertia,omega,qiw,mean_radius)
!  ANGULAR MOMENTUM LOSS WITHOUT INTERNAL ANGULAR MOMENTUM TRANSPORT.
!         DO I = 1,M
!            WOLD(I) = OMEGA(I)
!         END DO
         if(.not.disk_lock_engaged .and. ljdot0 .and. convective_flag(num_zones)) then
!  FIND MOMENT OF INERTIA OF THE SURFACE C.Z.
            moment_of_inertia_cz = 0.0D0
!  ENFORCE SB ROTATION (OR UNIFORM ROTATION LAW IN ENTIRE STAR) IF DESIRED.
            if(force_solid_body_rotation.or. (solid_body_mode_flag.eq.1))then
              envelope_boundary_zone = 1
            else
              envelope_boundary_zone = convective_zone_bounds(num_convective_zones,1)
            endif
            do zone_index = envelope_boundary_zone,num_zones
               moment_of_inertia_cz = moment_of_inertia_cz + moment_of_inertia(zone_index)
            end do
!  FIND LOWEST SHELL IN SURFACE CZ (IMAX)
!            IMAX = MXZONE(NZONE,1)
!  HSTOP IS THE MASS AT THE TOP OF THE C.Z.
!  HSBOT IS THE MASS AT THE BOTTOM OF THE C.Z.
            cz_mass_top = exp(ln10*log_total_mass)
            if(envelope_boundary_zone.gt.1)then
               cz_mass_bottom = 0.5D0*(enclosed_mass(envelope_boundary_zone)+ &
                    enclosed_mass(envelope_boundary_zone-1))
            else
               cz_mass_bottom = 0.0D0
            endif
            wind_loss_active = ljdot0
! MHP 10/02 UNUSED LFIRST REMOVED FROM CALL
! MHP 10/17 timestep average loss rate
!            FRACSTEP = 1.
            turnover%fracstep = 0.5
            call mwind(log_luminosity_lsun,full_timestep,cz_mass_bottom, &
                 cz_mass_top,envelope_boundary_zone,num_zones,wind_loss_active, &
                 omega_surface, &
                 total_mass_msun,log_teff,moment_of_inertia_cz, &
                 specific_angular_momentum)
!  FIND THE NEW RUN OF OMEGA IN THE CONVECTION ZONE AFTER THE WIND.
! MHP 10/02 REPLACED IEND WITH M IN CALL, DEFINED IEND=M
            iend = num_zones
! JNT 09/2025 FOR 05/15 REPLACE WCZ WITH WCZIMP
            call wczimp(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,envelope_boundary_zone,iend,eta_squared, &
                 moment_of_inertia,omega,qiw,mean_radius,num_zones)
         endif
         goto 9999
      endif
!  GETROT TAKES THE ANGULAR MOMENTUM DISTRIBUTION AND FINDS THE
!  ROTATION CURVE THAT CORRESPONDS TO IT.  CONVECTIVE REGIONS HAVE
!  SOLID BODY ROTATION ENFORCED ON THEM.
      call getrot(log_density,specific_angular_momentum,log_radius,log_mass, &
           shell_mass,am_transport_convective_flag,num_zones,eta_squared, &
           moment_of_inertia,omega,qiw,mean_radius)
      skip_diffusion_flag = .not.instability_transport_active .and. .not.wind_loss_active
      if(skip_diffusion_flag.or.full_timestep.le.0.0D0) goto 9999
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
      call physic(fp,ft,composition,log_density,hg,log_luminosity,log_pressure, &
           log_radius,log_mass,log_temperature,convective_flag,num_zones,log_teff)
! 8/17 DETERMINE TAUCZ AND PPHOT
!       CALL GETTAU(HCOMP,HR,HP,HD,HG,HS1,HT,FP,FT,TEFFL,  ! KC 2025-05-31
      call gettau(composition,log_radius,log_pressure,log_density,enclosed_mass, &
           log_temperature,fp,ft,log_teff, &
           log_total_mass,log_luminosity_lsun,num_zones,convective_flag,radius_at_bcz)
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
   30 continue
      redo_count = 0
!  ENTRY FOR DIFFUSION TIMESTEP CUTTING.
   40 continue
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
         theta_prev(1) = tho(1)
         wmst(1) = old_omega(1)
         do zone_index = 2,num_zones
            qwrmst(zone_index) = qwrst(zone_index)
            theta_prev(zone_index) = tho(zone_index)
            wmst(zone_index) = old_omega(zone_index)
         end do
! RECOMPUTE THETA
      else
         wmst(1) = omega(1)
         do zone_index = 2,num_zones
            omega_avg = 0.5D0*(omega(zone_index)+omega(zone_index-1))
            delta_radius_step = 10.0D0**log_radius_mid(zone_index)- &
                 10.0D0**log_radius_mid(zone_index-1)
            domega_dr = (omega(zone_index)-omega(zone_index-1))/delta_radius_step
            theta_prev(zone_index) = theta_mean(zone_index)*omega_avg*domega_dr
            qwrmst(zone_index) = domega_dr
            wmst(zone_index) = omega(zone_index)
         end do
      endif
      fx = elapsed_substep_time/full_timestep
      turnover%fracstep = fx
! INTERPOLATE LINEARLY IN TIME FOR THE MODEL STRUCTURE BETWEEN THE
! START AND END OF THE TIMESTEP.
! JVS
      call midmod(full_timestep,sub_timestep,fx,composition,log_density,hg, &
           specific_angular_momentum,log_luminosity,log_pressure,log_radius, &
           log_mass,enclosed_mass,shell_mass, &
           log_total_mass,log_temperature,convective_flag,first_call,num_zones, &
           eta_squared_mid,log_density_mid,hg_mid,moment_of_inertia_cz, &
           moment_of_inertia_mid, &
           log_luminosity_mid,log_pressure_mid,log_radius_mid,cz_mass_bottom, &
           cz_mass_top,log_temperature_mid,core_boundary_zone, &
           envelope_boundary_zone,fully_convective_flag,convective_flag_mid, &
           am_transport_convective_flag_mid,surface_cz_active,omega_mid, &
           mean_radius_mid,qiw_mid,radiative_zone_bounds,convective_zone_bounds, &
           num_radiative_zones,num_convective_zones)
! IF DESIRED, REMOVE ANGULAR MOMENTUM FROM OUTER CONVECTION ZONE
! USING A WEBER-DAVIS MAGNETIC WIND MODEL
      if(.not.disk_lock_engaged .and. wind_loss_active .and. surface_cz_active) then
         omega_surface = omega_mid(num_zones)
!  WIND CALCULATION WITHOUT INSTABILITIES
         if(fully_convective_flag.or..not.instability_transport_active) then
! MHP 10/02 UNUSED LFIRST REMOVED FROM CALL
            call mwind(log_luminosity_lsun,sub_timestep,cz_mass_bottom,cz_mass_top, &
                 envelope_boundary_zone,num_zones,wind_loss_active,omega_surface, &
                 total_mass_msun,log_teff,moment_of_inertia_cz, &
                 specific_angular_momentum)
!  FIND THE NEW OMEGA OF THE CONVECTION ZONE AFTER THE WIND.
! MHP 10/02 REPLACED IEND WITH M IN CALL, DEFINED IEND=M
            iend = num_zones
! JNT 09/2025 FOR 05/15 REPLACE WCZ WITH WCZIMP
            call wczimp(log_density_mid,specific_angular_momentum,log_radius_mid, &
                 log_mass,shell_mass,envelope_boundary_zone,iend,eta_squared_mid, &
                 moment_of_inertia_mid,omega_mid,qiw_mid,mean_radius_mid,num_zones)
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
          call mixcz(composition,shell_mass,am_transport_convective_flag_mid,num_zones)
! G Somers END
!  NOW SOLVE FOR LONG-TIMESCALE(SECULAR) INSTABILITIES.
!  THESE ARE TREATED USING DIFFUSION EQUATIONS.
         do 60 zone_index = 1,num_zones
            specific_angular_momentum_saved(zone_index) = specific_angular_momentum(zone_index)
   60    continue
! MHP 6/00 ADDED COD2,HV TO LIST RETURNED FROM SECULR
! FOR THE BUR-ST MIXING ROUTINES
!          CALL SECULR(DELTS,DT,HDM,HGM,HIM,HLM,HPM,HRM,HS,HS1,HS2,
!      *               HTM,HJMSAV,LCZM,LCM,M,MODEL,OMEGAM,
!      *               HJM,ETA2M,QIWM,R0M,HCOMP,LFIRST,IMIN,IMAX,BL,
!      *               HSTOT,SJTOT,SMASS,TEFFL,LREDO,IREDO,  ! KC 2025-05-31
         call seculr(sub_timestep,log_density_mid,hg_mid,moment_of_inertia_mid, &
              log_luminosity_mid,log_pressure_mid,log_radius_mid,log_mass, &
              enclosed_mass,shell_mass, &
              log_temperature_mid,specific_angular_momentum_saved, &
              am_transport_convective_flag_mid,num_zones,omega_mid, &
              specific_angular_momentum,eta_squared_mid,qiw_mid,mean_radius_mid, &
              composition,core_boundary_zone,envelope_boundary_zone, &
              log_luminosity_lsun, &
              log_total_mass,total_mass_msun,log_teff,redo_needed_flag,redo_count, &
              moment_of_inertia_cz,cz_mass_bottom,cz_mass_top,omega_surface,surface_cz_active, &
!      *               MRZONE,MXZONE,NRZONE,NZONE,  ! KC 2025-05-31
              cod2,diffusion_velocity,diffusion_solve_ok)
!  DIFFUSION TIMESTEP CUTTING REQUIRED IF LREDO IS TRUE.
!  RESET MEAN MOLECULAR WEIGHT,COMPOSITION,AND SPECIFIC ANGULAR MOMENTUM
!  TO THE VALUES THEY HAD PRIOR TO THE START OF THE STEP.
         if(redo_needed_flag)then
            fx = 2.0D0*sub_timestep/full_timestep
            elapsed_substep_time = elapsed_substep_time - 2.0D0*sub_timestep
            do 80 zone_index = 1,num_zones
               specific_angular_momentum(zone_index) = specific_angular_momentum_saved(zone_index)
               mix_phys%amum(zone_index) = mix_phys%amum(zone_index) - fx*(shell_temp%mean_molecular_weight(zone_index)-old_amu(zone_index))
               do 70 species_index = 1,num_species_tracked
                  composition(species_index,zone_index) = prev_model%old_composition(species_index,zone_index)
   70          continue
   80       continue
            goto 40
         endif
      endif
! MHP 05/02 TAKEN OUT OF LOOP SO THAT NUCLEAR BURNING
! IS ACCOUNTED FOR IN MODELS THAT ARE FULLY CONVECTIVE
! OR WITHOUT INSTABILITIES
! ADDED CHANGE FOR BURLICH-STORER TREATMENT OF MIXING PLUS
! BURNING - ONLY UPDATED IF NOT USED
      if(.not.burs_extrapolation_active)then
         do zone_index = 1,num_zones
            do species_index = 1,11
                  prev_model%old_composition(species_index,zone_index) = composition(species_index,zone_index)
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
         call convec(composition,log_density_mid,log_pressure_mid,log_radius_mid, &
              log_mass,log_temperature_mid,convective_flag_mid,num_zones, &
              radiative_zone_bounds,convective_zone_bounds, &
              convective_zone_bounds_burn,core_boundary_zone_cur, &
              envelope_boundary_zone_cur,num_radiative_zones,num_convective_zones, &
              num_convective_zones_burn)
! 11/91 CHANGED FOR LITHIUM BURNING WITH OVERSHOOT.
         if(lovstm .and. convective_flag_mid(num_zones))then
            light_burn%pressure_scale_height_end = alphae*exp(clndp*(log_pressure_mid(envelope_boundary_zone_cur)+ &
                 2.0D0*log_radius_mid(envelope_boundary_zone_cur) &
                 -log_density_mid(envelope_boundary_zone_cur)-cgl-log_mass(envelope_boundary_zone_cur)))
         else
            light_burn%pressure_scale_height_end = 0.0D0
         endif
! FIND LIGHT ELEMENT BURNING RATES AT THE END OF THE TIME STEP.
         call lirate88(composition,log_density_mid,log_temperature_mid,num_zones,2)
! STORE CURRENT "END OF STEP" QUANTITIES AS "BEGINNING" ONES FOR
! THE NEXT STEP.
! ADDED CHANGE FOR BURLICH-STORER TREATMENT OF MIXING PLUS
! BURNING - ONLY UPDATED IF NOT USED
         if(.not.instability_transport_active .or. .not.burs_extrapolation_active .or. fully_convective_flag)then
! COMPUTE BURNING.
!             CALL LIBURN(DT,HCOMP,HDM,HRM,HS1,HS2,HTM,JENV1,JENV0,M)  ! KC 2025-05-31
            call liburn(sub_timestep,composition,log_radius_mid,enclosed_mass, &
                 shell_mass,log_temperature_mid,envelope_boundary_zone_cur, &
                 envelope_boundary_zone_prev,num_zones)
            envelope_boundary_zone_prev = envelope_boundary_zone_cur
            light_burn%pressure_scale_height_start = light_burn%pressure_scale_height_end
            do 155 zone_index = 1,num_zones
               light_burn%rate_li6_start(zone_index) = light_burn%rate_li6(zone_index)
               light_burn%rate_li7_start(zone_index) = light_burn%rate_li7(zone_index)
               light_burn%rate_be9_start(zone_index) = light_burn%rate_be9(zone_index)
  155       continue
         else
! COMPUTE BURNING.
!             CALL LIBURN2(DT,HCOMP,HDM,HRM,HS1,HS2,HTM,JENV1,JENV0,M)  ! KC 2025-05-31
            call liburn2(sub_timestep,composition,log_radius_mid,enclosed_mass, &
                 shell_mass,log_temperature_mid,envelope_boundary_zone_cur, &
                 envelope_boundary_zone_prev,num_zones)
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
          call mixcz(composition,shell_mass,am_transport_convective_flag_mid,num_zones)
! G Somers END
!  ZERO OUT LOW ABUNDANCES.
         do 100 zone_index = 1,num_zones
            do 90 species_index = 12,15
               if(composition(species_index,zone_index).lt.1.0D-24)composition(species_index,zone_index)=0.0D0
   90       continue
  100   continue
! MHP 6/00 ADDED OVERWRITE OF HCOMPP FOR LIGHT ELEMENTS
! ADDED CHANGE FOR BURLICH-STORER TREATMENT OF MIXING PLUS
! BURNING - ONLY UPDATED IF NOT USED
         if(.not.instability_transport_active.or.fully_convective_flag)then
            do zone_index = 1,num_zones
               do species_index = 12,15
                  prev_model%old_composition(species_index,zone_index) = composition(species_index,zone_index)
               end do
            end do
         else if(.not.burs_extrapolation_active)then
            do zone_index = 1,num_zones
               do species_index = 12,15
                  prev_model%old_composition(species_index,zone_index) = composition(species_index,zone_index)
               end do
            end do
         endif
      endif
! MHP 6/00 NOW ADDED THE OPTION OF PERFORMING A BUR-ST EXTRAPOLATION
! TO ZERO TIMESTEP FOR THE COMBINATION OF MIXING AND NUCLEAR BURNING.
      if(instability_transport_active .and. .not.fully_convective_flag)then
         if(burs_extrapolation_active)then
        call bursmix(cod2,sub_timestep,composition,log_density_mid,log_luminosity_mid, &
             log_pressure_mid,log_radius_mid,log_mass,enclosed_mass,shell_mass,log_total_mass, &
             log_temperature_mid,diffusion_velocity,envelope_boundary_zone,core_boundary_zone, &
             envelope_boundary_zone_prev,envelope_boundary_zone_cur,diffusion_solve_ok, &
             am_transport_convective_flag_mid,num_zones,radiative_zone_bounds, &
             convective_zone_bounds,num_radiative_zones,num_convective_zones)
         endif
      endif
!  RETURN FOR NEXT SMALL DIFFUSION TIMESTEP IF NEEDED.
      if(elapsed_substep_time.lt.full_timestep)goto 30
!  UPDATE OMEGA ARRAY TO REFLECT NEW ANGULAR MOMENTUM DISTRIBUTION.
      do 110 zone_index = 1,num_zones
         omega(zone_index) = omega_mid(zone_index)
  110 continue
! MHP 3/96 ADDED CALL TO RECOMPUTE SELF-CONSISTENT SET OF OMEGAS
      call getrot(log_density,specific_angular_momentum,log_radius,log_mass, &
           shell_mass,am_transport_convective_flag,num_zones,eta_squared, &
           moment_of_inertia,omega,qiw,mean_radius)
 9999 continue
      if(lprt0_placeholder)then
         log_luminosity_lsun = log10(log_luminosity(num_zones))
         log_radius_surface = 0.5D0*(log_luminosity_lsun + log10_solar_luminosity &
              - 4.0D0*log_teff - c4pil - csigl)
         fx = exp(ln10*3.0D0*(log_radius(num_zones)-log_radius_surface))
         surface_quad_term = fx*quadrupole_moment(num_zones)
         surface_potential = exp(ln10*(cgl+log_total_mass-log_radius_surface))
         write(*,9911)surface_quad_term,surface_potential,-1.5D0*surface_quad_term/surface_potential
 9911    format(1X,'QUAD ',1PE12.3,' PHIS ',E12.3,' 3/2 QUAD/G ', &
     E12.3)
      endif
      return
end subroutine getw

!----------------------------------------------------------------------
! seculr
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original seculr.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! SECULR CALCULATES THE TRANSFER OF ANGULAR MOMENTUM(AND THE ASSOCIATED
! TRANSPORT OF COMPOSITION) DUE TO SECULAR ROTATIONAL INSTABILITIES.
!
! MHP 6/00 ADDED MRZONE,MXZONE,NRZONE,NZONE FOR BS MIX PLUS BURN
! (long since removed from the active call list -- see the commented-
! out historical argument lines this file inherited from the .f
! version).
!
! Dummy-argument names/order/types for this routine are fixed by the
! existing call in getw.f90 (already converted, out of scope for this
! batch) and must not change; only the dummy names chosen here are
! free.
!
! A few local/dummy quantities are named differently by the various
! already-converted callees this file dispatches to. Where that
! happens the name used here is the one shared by the most callees (or
! the more descriptive of the two), and the alternate name(s) are
! noted alongside the declaration, following the convention set by
! crrect.f90:
!  - dynamical_shear_omega_limit (QWRMAX): vcirc.f90 calls the same
!    quantity dlnomega_dlnr_max.
!  - radius_unlogged (HRU): vcirc.f90/getfc.f90/getqua.f90/rotgrid.f90
!    call it radius; codiff.f90 calls it radius_mid_prev.
!  - mixing_diffusion_coeff (COD2, dummy arg): ndifcom.f90 calls it
!    diffusion_coeff.
!  - diffusion_velocity (HV, dummy arg): vcirc.f90 calls it
!    total_circulation_velocity.
!  - eq_mixing_diffusion_coeff (ECOD2): ndifcom.f90 calls it
!    equally_spaced_diffusion_coeff.
!  - eq_angular_momentum (EJ): eq2mod.f90 calls it angular_momentum.
!  - eq_am_diffusion_coeff (ECOD): dcoeft.f90 calls it diffusion_coeff.
!  - eq_delta_angular_momentum (DJ): eq2mod.f90 calls it
!    delta_angular_momentum; tridia.f90 calls it dj.
!  - eq_mass (EM): ndifcom.f90 calls it equally_spaced_mass.
!  - diffusion_converged (LOKAD): checkj.f90 calls it
!    already_converged_flag.
!  - wind_loss_explicit/wind_loss_implicit (WIND1/WIND2): mcowind.f90
!    calls the same pair domega_start/domega_end.
subroutine seculr(sub_timestep, log_density, local_gravity, &
     moment_of_inertia, luminosity, log_pressure, log_radius, log_mass, &
     enclosed_mass, shell_mass, log_temperature, &
     specific_angular_momentum_saved, am_transport_convective_flag, &
     num_zones, omega, specific_angular_momentum, eta_squared, qiw, &
     mean_radius, composition, zone_min, zone_max, log_luminosity_lsun, &
     log_total_mass, total_mass_msun, log_teff, redo_flag, cut_count, &
     cz_moment_of_inertia, cz_mass_bottom, cz_mass_top, omega_surface, &
     surface_cz_active, mixing_diffusion_coeff, diffusion_velocity, &
     diffusion_solve_ok)

      use oldmod_lib
      use luout_lib
      use const_lib
      use numerics_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(inout) :: sub_timestep
      double precision, intent(inout) :: log_density(json)
      double precision, intent(in) :: local_gravity(json)
      double precision, intent(inout) :: moment_of_inertia(json)
      double precision, intent(in) :: luminosity(json)
      double precision, intent(in) :: log_pressure(json)
      double precision, intent(inout) :: log_radius(json)
      double precision, intent(inout) :: log_mass(json)
      double precision, intent(in) :: enclosed_mass(json)
      double precision, intent(inout) :: shell_mass(json)
      double precision, intent(in) :: log_temperature(json)
      double precision, intent(in) :: specific_angular_momentum_saved(json)
      logical, intent(inout) :: am_transport_convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: omega(json)
      double precision, intent(inout) :: specific_angular_momentum(json)
      double precision, intent(inout) :: eta_squared(json)
      double precision, intent(inout) :: qiw(json)
      double precision, intent(inout) :: mean_radius(json)
      double precision, intent(inout) :: composition(15,json)
      integer, intent(in) :: zone_min, zone_max
      double precision, intent(in) :: log_luminosity_lsun
      double precision, intent(in) :: log_total_mass
      double precision, intent(in) :: total_mass_msun
      double precision, intent(in) :: log_teff
      logical, intent(out) :: redo_flag
      integer, intent(inout) :: cut_count
      double precision, intent(out) :: cz_moment_of_inertia
      double precision, intent(in) :: cz_mass_bottom
      double precision, intent(in) :: cz_mass_top
      double precision, intent(inout) :: omega_surface
      logical, intent(in) :: surface_cz_active
      double precision, intent(inout) :: mixing_diffusion_coeff(json)
      double precision, intent(inout) :: diffusion_velocity(json)
      logical, intent(inout) :: diffusion_solve_ok





! MHP 8/17 added excen, c_2 to common block for Matt et al. 2012 cent.
! term.
! common/cwind/: only ljdot0 is used here. Naming matches getw.f90/
! midmod.f90.
      double precision :: wind_saturation_omega, exmd, &
           wind_law_omega_exponent, extau, exr, exm, exl, expr, &
           constfactor, structfactor, excen, c_2
      logical :: ljdot0
      common/cwind/ wind_saturation_omega, exmd, wind_law_omega_exponent, &
           extau, exr, exm, exl, expr, constfactor, structfactor, excen, &
           c_2, ljdot0



! 3/92 common block added for centered dw/dlnr.
! common/intvar/: only interface_radius (RM) is used here. Naming
! matches setupv.f90/vcirc.f90/getfc.f90/codiff.f90.
      double precision :: interface_luminosity(json), delami(json), &
           delmi(json), dm(json), epsilm(json), interface_gravity_factor(json), &
           hs3(json), pm(json), qdtmi(json), interface_radius(json), tm(json)
      common/intvar/ interface_luminosity, delami, delmi, dm, epsilm, &
           interface_gravity_factor, hs3, pm, qdtmi, interface_radius, tm


      double precision :: es_circulation_velocity(json), &
           es_circulation_velocity_prev(json), secular_shear_velocity(json), &
           secular_shear_velocity_prev(json), hle(json), &
           gsf_circulation_velocity(json), gsf_circulation_velocity_prev(json), &
           mu_gradient_velocity(json)
      common/temp2/ es_circulation_velocity, es_circulation_velocity_prev, &
           secular_shear_velocity, secular_shear_velocity_prev, hle, &
           gsf_circulation_velocity, gsf_circulation_velocity_prev, &
           mu_gradient_velocity

! common/oldab/: composition_snapshot is used here to restore
! composition between diffusion iterations. Naming matches getw.f90/
! midmod.f90.
      double precision :: composition_snapshot(15,json)
      common/oldab/ composition_snapshot

! Add option for variable Fc.
! common/varfc/: use_diffusion_advection_transport (LDIFAD) gates the
! banded-matrix solve path; vfc/lvfc are set but otherwise unused here.
! Naming matches rotgrid.f90/checkj.f90/codiff.f90/dcoeft.f90/getfc.f90.
      double precision :: vfc(json)
      logical :: lvfc, use_diffusion_advection_transport
      common/varfc/ vfc, lvfc, use_diffusion_advection_transport


! common/rotprt/: not used in this file (LPRT is hardcoded .false.
! below rather than computed from lprt0_placeholder). Naming matches
! getw.f90/wrtout.f90.
      logical :: lprt0_placeholder
      common/rotprt/ lprt0_placeholder


! common/egrid/: only ntot is used here (as the equally-spaced point
! count passed to dcoeft/dadcoeft/tridia). Naming matches rotgrid.f90/
! checkj.f90.
      double precision :: chi(json), echi(json), es1(json), dchi
      integer :: ntot
      common/egrid/ chi, echi, es1, dchi, ntot

! MHP 9/93
! common/notran/: no_am_transport_in_core (originally LNOJ), "9/93
! MIXING WITHOUT TRANSPORT ADDED" -- forces convergence without
! angular-momentum transport in the core when set. Naming matches
! vcirc.f90.
      logical :: no_am_transport_in_core
      common/notran/ no_am_transport_in_core

! Time change of theta.
! common/oldrot2/: only wmst is used here (the surface omega from the
! start of the timestep). Naming matches getw.f90/hpoint.f90/
! setupv.f90/vcirc.f90.
      double precision :: tho(json), theta_new(json), theta_mean(json), &
           del_grad_diff_interface(json), es_relaxation_factor(json), &
           theta_prev(json), qwrst(json), wmst(json), qwrmst(json)
      common/oldrot2/ tho, theta_new, theta_mean, del_grad_diff_interface, &
           es_relaxation_factor, theta_prev, qwrst, wmst, qwrmst

! MHP 3/99 added flag to treat the entire star as 'convective' for
! angular momentum purposes.
! JNT 09/25 for 05/15 add impjmod
! common/sbrot/: both used here. Naming matches getw.f90/midmod.f90.
      logical :: force_solid_body_rotation
      integer :: solid_body_mode_flag
      common/sbrot/ force_solid_body_rotation, solid_body_mode_flag

      save
! MHP 6/00 added MRZONE,MXZONE for BS mixing plus burning (no longer
! part of the active call list).

! --- locals (JSON-dimensioned) ---
      double precision :: dynamical_shear_omega_limit(json), eq_mass(json), &
           am_diffusion_coeff(json), eq_mixing_diffusion_coeff(json), &
           eq_moment_of_inertia(json), eq_angular_momentum(json), &
           eq_am_diffusion_coeff(json), radius_unlogged(json), radius_mid(json), &
           dlnomega_dlnr(json), eq_delta_angular_momentum(json), &
           specific_angular_momentum_prev(json), omega_start(json), &
           eq_omega(json)
      integer :: print_zone_id(json)

! Tridiagonal-solve work arrays (Thomas algorithm) threaded from
! dcoeft (fills them) to tridia (solves): were originally shared via
! common/tridi/ (positional storage), converted (2026, GUIDELINES.md)
! to explicit arguments since this is real per-call data flow, not
! global configuration. tridia's solution(:) (the new omega
! distribution) was never actually read back here -- this call site
! only used dj/sumdj (eq_delta_angular_momentum/
! sum_delta_angular_momentum below) -- so it's captured in an unused
! local. surface_wind_loss_term is dcoeft's computed surface wind
! angular-momentum-loss term, previously smuggled to tridia via
! gamma_elim(n) (see tridia.f90's header note).
      double precision :: sub_diag(json), diag(json), super_diag(json), &
           rhs(json), unused_tridia_solution(json), surface_wind_loss_term

! --- other locals ---
      logical :: print_diffusion_flag
      logical :: disk_lock_active
      logical :: lcz_first_zone, lcz_last_zone
      double precision :: total_luminosity
      integer :: i, j
      double precision :: log_radius_center
! dlnr_weight (originally DR1) is the sum of the two one-sided inverse
! log-radius spacings used to build the centered dlnomega/dlnr
! estimate; not otherwise named elsewhere.
      double precision :: dlnr_weight
      integer :: iteration
      double precision :: omega_surface_start
      logical :: wind_loss_active
      integer :: solid_body_zone_start, solid_body_zone_end
      logical :: in_unstable_region
      integer :: scan_start_zone
      logical :: unstable_zone_found
      integer :: zone_begin, zone_end
      double precision :: grid_spacing
      logical :: single_interface_flag
      logical :: fix_omega_at_surface
      logical :: diffusion_converged
      double precision :: wind_loss_explicit, wind_loss_implicit
      double precision :: sum_delta_angular_momentum
      double precision :: solid_cz_mass_bottom, solid_cz_mass_top
      integer :: solid_start_zone
      integer :: species_begin, species_end
      integer :: print_zone_count
      double precision :: delta_h1, delta_he3, delta_c12, delta_c13, &
           delta_n14, delta_li7, delta_be9
! constant_diffusion_coeff_flag/constant_diffusion_coeff (originally
! LCODM/CODM, MHP 8/13 "treat entire domain as unstable if a constant
! diffusion coefficient is being added") are implicitly typed in the
! original and never assigned anywhere in seculr.f, nor declared in
! any common block reachable from it, nor referenced by any other file
! in the already-converted codebase -- they are permanently at their
! (undefined/zero-initialized) default, so this branch is dead code.
! Preserved exactly, with explicit types matching the original's
! implicit rules (L->logical, C->double precision).
      logical :: constant_diffusion_coeff_flag
      double precision :: constant_diffusion_coeff

!  PRINT OUT DETAILS OF THE DIFFUSION EVERY TIME SHORT OR LONG
!  OUTPUT OF THE MODEL IS GENERATED.
!
! G Somers 11/14, EXCLUDE THE IO TO THE FULL FILE. THE CODE
! WILL NO LONGER REPORT THE CHANGES TO J AT EACH POINT, BECAUSE
! THIS CAN BE EASILY INFERRED FROM THE EXTENDED SHORT FILE. IF
! DESIRED, THIS OUTPUT CAN BE RETURNED.
!      IF(.NOT.LPRT0)THEN
!         LPRT = MOD(MODEL+1,NPRTMOD).EQ.0
!      ELSE
!         LPRT = .TRUE.
!      ENDIF
!      IF(.NOT.LSCRIB.OR..NOT.LPRTIN) LPRT = .FALSE.
      print_diffusion_flag = .false.
! MHP 9/94
! DISK LOCKING CHECKED
      disk_lock_active = .false.
      if(disk_locking_active .and. disk_lifetime.le.disk_temperature) &
           disk_lock_active = .true.
!  FIRST AND LAST SHELL BOUNDARY CONDITIONS REQUIRE THAT THESE SHELLS BE
!  TREATED AS 'CONVECTIVE' FOR DIFFUSION PURPOSES. FIX THIS.
      lcz_first_zone = am_transport_convective_flag(1)
      lcz_last_zone = am_transport_convective_flag(num_zones)
!  SET UP FACTORS USED TO COMPUTE DIFFUSION VELOCITIES.
      call setupv(log_density,local_gravity,luminosity,log_pressure, &
           log_radius,enclosed_mass,log_temperature,zone_min,zone_max, &
           num_zones,radius_unlogged,dynamical_shear_omega_limit)
      total_luminosity = solar_luminosity_cgs*luminosity(num_zones)
!  COMPUTE ANGULAR VELOCITY GRADIENTS
      do 10 i = zone_min,zone_max
! CENTER LOGARITHMIC DERIVATIVE.
         log_radius_center = log10(interface_radius(i))
         dlnr_weight = 1.0D0/ln10/(log_radius(i)-log_radius_center)+ &
              1.0D0/ln10/(log_radius_center-log_radius(i-1))
         dlnomega_dlnr(i) = 0.25D0*(omega(i)-omega(i-1))*dlnr_weight
   10 continue
      call getqua(log_density,local_gravity,radius_unlogged,omega,num_zones)
      do i = 1,num_zones
         vfc(i) = 0.0D0
      end do
!  CHECK STABILITY OF THE MODEL AND COMPUTE CIRCULATION VELOCITIES.
!  MECHANISMS CONSIDERED:SECULAR AND DYNAMICAL SHEAR,EDDINGTON CIRCULATION,
!  AND THE GSF INSTABILITY.
      iteration = 1
! MHP 06/02 added timestep to call
      call vcirc(log_radius,radius_unlogged,zone_min,zone_max,iteration, &
           am_transport_convective_flag,num_zones,omega,unstable_zone_found, &
           dlnomega_dlnr,dynamical_shear_omega_limit,diffusion_velocity, &
           total_luminosity,sub_timestep,log_pressure)
      if(use_diffusion_advection_transport) &
           call getfc(log_density,radius_unlogged,diffusion_velocity, &
           zone_min,zone_max,omega)
!  STORE INITIAL ANGULAR MOMENTUM DISTRIBUTION.
      do 20 i =1,num_zones
         specific_angular_momentum_prev(i) = specific_angular_momentum_saved(i)
         omega_start(i) = omega(i)
   20 continue
!  STORE INITIAL SURFACE ANGULAR VELOCITY FOR USE IN ANGULAR MOMENTUM
!  LOSS CALCULATIONS.
!      WBEG = OMEGA(M)
      omega_surface_start = wmst(num_zones)
      diffusion_solve_ok = .false.
!  ON THE FIRST LEVEL OF ITERATION, THE UNPERTURBED MODEL IS USED TO
!  CALCULATE THE DIFFUSION VELOCITIES. ON THE SECOND AND SUBSEQUENT
!  ITERATIONS, THE ANGULAR MOMENTUM DISTRIBUTION FROM THE PREVIOUS
!  ITERATION IS USED TO GET THE VELOCITIES.  THIS 'NEW' VELOCITY IS THEN
!  AVERAGED WITH THE VELOCITY FOUND IN THE PREVIOUS ITERATION TO GET A
!  CORRECTED V AND THUS A MORE ACCURATE RUN OF DIFFUSION COEFFICIENTS.
      do 100 iteration = 1,itdif2
         omega_surface = omega(num_zones)
         if(iteration.gt.1)then
!  COMPUTE NEW RUN OF ANGULAR VELOCITIES (AVERAGE OF INITIAL AND
!  LATEST VALUES).
            do i = 1,num_zones
               omega(i) = 0.5D0*(omega(i)+omega_start(i))
            end do
!  COMPUTE NEW RUN OF ANGULAR VELOCITY GRADIENTS.
            do 30 i = zone_min,zone_max
! CENTER LOGARITHMIC DERIVATIVE.
               log_radius_center = log10(interface_radius(i))
               dlnr_weight = 1.0D0/ln10/(log_radius(i)-log_radius_center)+ &
                    1.0D0/ln10/(log_radius_center-log_radius(i-1))
               dlnomega_dlnr(i) = 0.25D0*(omega(i)-omega(i-1))*dlnr_weight
   30       continue
!  ON 2ND AND SUBSEQUENT ITERATIONS,COMPUTE CHARACTERISTIC VELOCITIES
!  FOR THE NEW RUN OF OMEGA AND COMPOSITION FOUND IN THE PREVIOUS
!  ITERATION.
! MHP 06/02 added pressure and timestep to call
            call vcirc(log_radius,radius_unlogged,zone_min,zone_max, &
                 iteration,am_transport_convective_flag,num_zones,omega, &
                 unstable_zone_found,dlnomega_dlnr,dynamical_shear_omega_limit, &
                 diffusion_velocity,total_luminosity,sub_timestep,log_pressure)
            if(use_diffusion_advection_transport) &
                 call getfc(log_density,radius_unlogged,diffusion_velocity, &
                 zone_min,zone_max,omega)
!  NOW THAT THE NEW DIFFUSION VELOCITIES HAVE BEEN COMPUTED, RESET THE
!  ANGULAR MOMENTUM AND COMPOSITION ARRAYS TO THEIR ORIGINAL VALUES.
            do 50 i = 1,num_zones
               specific_angular_momentum(i) = specific_angular_momentum_saved(i)
               do 40 j = 1,4
                  composition(j,i) = composition_snapshot(j,i)
   40          continue
   50       continue
! MHP 10/91 CHANGED TO REMIX CZ'S TO THEIR PROPER DEPTH!
! OTHERWISE, DRASTIC ERRORS OCCUR IN THE PRESENCE OF A DEEPENING CZ
! (THE LOCAL ABUNDANCE AT THE CZ BASE PRIOR TO THE DEEPENING OF THE
!  CZ IS ASSIGNED TO THE WHOLE CZ - X AND HE3 ONLY, AND ONLY FOR A
!  CZ WHICH IS BECOMING DEEPER.)
!  ENSURE THAT CONVECTIVE REGIONS ARE FULLY MIXED.
!  JVS 0212 CALL MIXCZ(HCOMP,HS2,LCZ,M)
! KC 2025-05-30 addressed warning messages from Makefile.legacy
! C G Somers 6/14, SET IMIX = .FALSE. SO THE CORRECT GRADS ARE USED.
!          IMIX = .FALSE.
!          CALL MIXCZ(HCOMP,HS2,HS1,LCZ,HR,HP,HD,HG,M,IMIX)
! G Somers 6/14, SET LIMIX = .FALSE. SO THE CORRECT GRADS ARE USED.
! mix_grads_flag (originally LIMIX) is set but the active MIXCZ call
! below does not take it as an argument -- vestigial, kept exactly as
! in the original.
         call mixcz(composition,shell_mass,am_transport_convective_flag, &
              num_zones)
! G Somers END
         endif
!  CALCULATE LOSS OF ANGULAR MOMENTUM DUE TO WIND FOR AN
!  ISOLATED SURFACE C.Z.(NO COUPLING WITH INTERIOR VIA DIFFUSION).
         if(surface_cz_active.and.diffusion_velocity(zone_max).eq.0.0D0 &
              .and. .not.disk_lock_active)then
!  FIND MOMENT OF INERTIA OF THE SURFACE C.Z.
            cz_moment_of_inertia = 0.0D0
            do 55 i = zone_max,num_zones
               cz_moment_of_inertia = cz_moment_of_inertia + moment_of_inertia(i)
   55       continue
            wind_loss_active = ljdot0
! MHP 10/02 UNUSED LFIRST REMOVED FROM CALL
            call mwind(log_luminosity_lsun,sub_timestep,cz_mass_bottom, &
                 cz_mass_top,zone_max,num_zones,wind_loss_active,omega_surface, &
                 total_mass_msun,log_teff,cz_moment_of_inertia, &
                 specific_angular_momentum)
! REMOVE TORQUE FROM ENTIRE STAR
! JNT 09/25 FOR 05/15 IMPJMOD=1 SAME AS LSOLID
         else if(surface_cz_active .and. (force_solid_body_rotation .or. &
              (solid_body_mode_flag.eq.1)))then
            wind_loss_active = ljdot0
            solid_cz_mass_bottom = 0.0D0
            solid_cz_mass_top = exp(ln10*log_total_mass)
            solid_start_zone = 1
            call mwind(log_luminosity_lsun,sub_timestep,solid_cz_mass_bottom, &
                 solid_cz_mass_top,solid_start_zone,num_zones,wind_loss_active, &
                 omega_surface,total_mass_msun,log_teff,cz_moment_of_inertia, &
                 specific_angular_momentum)
!            WRITE(*,*)HJM(1),HJM(M)
            solid_body_zone_start = 1
            solid_body_zone_end = num_zones
            call solid(log_density,specific_angular_momentum,log_radius, &
                 log_mass,shell_mass,solid_body_zone_start,solid_body_zone_end, &
                 eta_squared,moment_of_inertia,omega,qiw,mean_radius,num_zones)
!            WRITE(*,*)OMEGA(1),OMEGA(M)
         endif
!  IF LDO=F,NO INSTABILITIES OCCUR (STABLE AGAINST ALL MECHANISMS).
         if(.not.unstable_zone_found) goto 9999
!  TREAT CENTRAL AND SURFACE ZONES AS ALWAYS CONVECTIVE
!  (SHOULD BE FIXED TO GIVE BETTER CENTRAL/SURFACE B.C.)
         am_transport_convective_flag(1) = .true.
         am_transport_convective_flag(num_zones) = .true.
!  UNSTABLE REGION EXISTS.
!  FIND DIFFUSION COEFFICIENTS(COD) FOR ALL UNSTABLE REGIONS.
         call codiff(radius_unlogged,num_zones,radius_mid, &
              am_diffusion_coeff,mixing_diffusion_coeff)
!  EACH UNSTABLE REGION IS SOLVED SEPARATELY STARTING HERE.
!  LTEST IS SET T IF A NON-ZERO VELOCITY IS ENCOUNTERED.
!  IBEG IS THE ZONE BELOW THE FIRST NON-ZERO V;IEND IS THE ZONE ABOVE
!  THE LAST NON-ZERO V.
         in_unstable_region = .false.
         scan_start_zone = zone_min
   60    continue
         unstable_zone_found = .false.
! MHP 8/13 TREAT ENTIRE DOMAIN AS UNSTABLE IF A CONSTANT DIFFUSION
! COEFFICIENT IS BEING ADDED
         if(constant_diffusion_coeff_flag.and.constant_diffusion_coeff.gt.0.0D0)then
            unstable_zone_found = .true.
            zone_begin = zone_min - 1
            zone_end = zone_max
            scan_start_zone = zone_max + 1
            goto 80
         endif
         do 70 j = scan_start_zone,zone_max
            if(diffusion_velocity(j).gt.0.0D0) then
               unstable_zone_found = .true.
               if(.not.in_unstable_region) then
!  START OF UNSTABLE REGION
                  in_unstable_region = .true.
                  zone_begin = j - 1
               endif
            else if(in_unstable_region) then
!  END OF UNSTABLE REGION
               zone_end = j - 1
               in_unstable_region = .false.
               scan_start_zone = j + 1
               goto 80
            endif
   70    continue
!  IF THE LAST INTERFACE IS UNSTABLE (NON-ZERO V) ENSURE THAT IEND IS SET
!  PROPERLY.
         if(in_unstable_region) zone_end = zone_max
         scan_start_zone = zone_max + 1
   80    continue
!  IF NO NON-ZERO V'S ENCOUNTERED, EXIT.
         if(.not.unstable_zone_found) goto 90
! MHP 08/03 REMOVED OBSOLETE EQUAL ROUTINE
!         IF(M.GT.1)THEN
!  TRANSFORM TO EQUAL GRID SPACING IN CHI FOR THE REGION.
! CHI IS A NORMALIZED SUM OF THE VARIABLES USED TO PLACE POINTS
! IN THE HENYEY SCHEME, CHOSEN SUCH THAT THE GRID USED IN THE
! ANGULAR MOMENTUM EVOLUTION IS CLOSE TO THE GRID STORED FOR
! THE STRUCTURAL EVOLUTION.
         call rotgrid(am_diffusion_coeff,mixing_diffusion_coeff,log_density, &
              moment_of_inertia,specific_angular_momentum_saved,luminosity, &
              log_pressure,log_radius,radius_unlogged,log_mass,enclosed_mass, &
              shell_mass,log_total_mass,zone_begin,zone_end, &
              am_transport_convective_flag,num_zones,omega_start,grid_spacing, &
              eq_am_diffusion_coeff,eq_mixing_diffusion_coeff, &
              eq_moment_of_inertia,eq_angular_momentum,eq_mass,eq_omega, &
              single_interface_flag)
!         ELSE
!         CALL EQUAL(BL,COD,COD2,HD,HI,HJMSAV,HRU,HS,HS1,HS2,HSTOT,
!     *              IBEG,IEND,LCZ,M,TEFFL,DR,ECOD,ECOD2,EI,EJ,
!     *              EM,ES1,EW,WSAV,LDUM2,NTOT)
!         ENDIF
!  LDUM2=T IF TWO C.Z.'S ARE SEPARATED BY ONE RADIATIVE ZONE;
!  SKIP IF THIS OCCURS.
         if(single_interface_flag) then
            if(scan_start_zone.le.zone_max) then
               goto 60
            else
               goto 90
            endif
         endif
! MHP 3/09 SKIP ANGULAR MOMENTUM EVOLUTION FOR SOLID BODY MODELS
! JNT 09/25 FOR 05/15 IMPJMOD=1 IS THE SAME AS LSOLID
         if(.not.force_solid_body_rotation .and. &
              (solid_body_mode_flag.ne.1))then
!  CHECK IF SURFACE C.Z. IS PART OF THE UNSTABLE REGION.
!  IF SO,CALCULATE TERMS FOR DIFFUSION.
         if(surface_cz_active.and.zone_end.eq.zone_max)then
            if(disk_lock_active)then
               wind_loss_explicit = 0.0D0
               wind_loss_implicit = 0.0D0
               fix_omega_at_surface = .true.
            else
               fix_omega_at_surface = .false.
               if(ljdot0)then
                  cz_moment_of_inertia = eq_moment_of_inertia(ntot)
                  omega_surface = omega(num_zones)
                  call mcowind(log_luminosity_lsun,sub_timestep, &
                       cz_moment_of_inertia,iteration,omega_surface, &
                       total_mass_msun,log_teff,omega_surface_start, &
                       wind_loss_explicit,wind_loss_implicit)
               else
                  wind_loss_explicit = 0.0D0
                  wind_loss_implicit = 0.0D0
               endif
            endif
         else
            fix_omega_at_surface = .false.
            wind_loss_explicit = 0.0D0
            wind_loss_implicit = 0.0D0
         endif
         diffusion_converged = .false.
         if(.not.use_diffusion_advection_transport)then
!  SET UP DIFFUSION EQUATION ARRAYS TO SOLVE FOR OMEGA AT END OF TSTEP
            call dcoeft(eq_am_diffusion_coeff,grid_spacing,sub_timestep, &
                 eq_moment_of_inertia,eq_angular_momentum,eq_omega,ntot, &
                 wind_loss_explicit,wind_loss_implicit,fix_omega_at_surface, &
                 sub_diag,diag,super_diag,rhs,surface_wind_loss_term)
!  SOLVE MATRIX FOR THE RUN OF OMEGA AT THE END OF THE TIMESTEP AT THE
!  EQUALLY SPACED GRID POINTS.
            call tridia(ntot,eq_moment_of_inertia, &
                 eq_delta_angular_momentum,sum_delta_angular_momentum, &
                 sub_diag,diag,super_diag,rhs,unused_tridia_solution, &
                 surface_wind_loss_term)
!  TRANSFORM THE NEW ANGULAR MOMENTUM DISTRIBUTION BACK TO THE ORIGINAL GRID
!  POINTS IN THE UNSTABLE REGION.
         else
! SOLVE FOR OMEGA AND ITS DERIVATIVES IN A BAND MATRIX
            call dadcoeft(grid_spacing,sub_timestep,eq_moment_of_inertia, &
                 eq_omega,ntot,wind_loss_explicit,wind_loss_implicit, &
                 eq_delta_angular_momentum,eq_mixing_diffusion_coeff, &
                 sum_delta_angular_momentum,fix_omega_at_surface, &
                 diffusion_converged)
         endif
! MHP 08/03 REMOVED OBSOLETE EQUAL2 ROUTINE
!         IF(M.GT.1)THEN
! TRANSFORM BACK TO THE ORIGINAL GRID
            call eq2mod(eq_delta_angular_momentum,eq_angular_momentum, &
                 shell_mass,zone_begin,zone_end,am_transport_convective_flag, &
                 num_zones,sum_delta_angular_momentum,specific_angular_momentum)
!         ELSE
!         CALL EQUAL2(DJ,EJ,EM,ES1,HS1,IBEG,IEND,LCZ,M,NTOT,HJM)
!         ENDIF
!  PERFORM COMPOSITION DIFFUSION.
!  UNTIL THE FINAL ITERATION, ONLY COMPOSITION DIFFUSION OF SPECIES WHICH
!  AFFECT GRADIENTS IN MEAN MOLECULAR WEIGHT IS COMPUTED (H,HE3,HE4).
!  ON THE FINAL ITERATION, DIFFUSION OF ALL SPECIES IS PERFORMED.
! MHP 08/03 REMOVED OBSOLETE DIFCOM ROUTINE
!         IF(M.GT.1)THEN
         endif
         species_begin = 1
         species_end = 4
         call ndifcom(sub_timestep,mixing_diffusion_coeff, &
              eq_mixing_diffusion_coeff,eq_mass,log_density,luminosity, &
              log_pressure,log_radius,log_mass,enclosed_mass,shell_mass, &
              log_total_mass,diffusion_velocity,zone_begin,zone_end,zone_max, &
              zone_min,am_transport_convective_flag,diffusion_solve_ok, &
              num_zones,composition,species_begin,species_end)
!         ELSE
!         CALL DIFCOM(DR,DT,COD2,ECOD2,EM,ES1,HRU,HS,HS1,HS2,HV,
!     *               IBEG,IEND,IMIN,LCZ,LOK,M,NTOT,HCOMP)
!         ENDIF
!  RETURN FOR NEXT REGION IF APPLICABLE
         if(scan_start_zone.le.zone_max) goto 60
   90    continue
! CHECK SOLUTION,UPDATE OMEGA,AND SEE IF ANOTHER ITERATION IS NEEDED.
         am_transport_convective_flag(1) = lcz_first_zone
         am_transport_convective_flag(num_zones) = lcz_last_zone
! mhp 10/02 unused ecod, ecod2 no longer passed
!         WRITE(*,*)LSOLID,OMEGA(1),OMEGA(M)
         call checkj(log_density,specific_angular_momentum_prev, &
              specific_angular_momentum_saved,log_radius,log_mass,shell_mass, &
              diffusion_velocity,zone_min,zone_max,iteration, &
              am_transport_convective_flag,num_zones,sub_timestep,eta_squared, &
              moment_of_inertia,specific_angular_momentum,cut_count, &
              diffusion_solve_ok,redo_flag,omega,qiw,mean_radius,omega_start, &
              print_zone_id,print_zone_count,diffusion_converged)
!         WRITE(*,*)OMEGA(1),OMEGA(M)
! CHECK COMPOSITION DIFFUSION AND RECOMPUTE MEAN MOLECULAR WEIGHT.
         if(.not.redo_flag)call checkc(composition,iteration, &
              print_diffusion_flag,num_zones,sub_timestep,cut_count, &
              diffusion_solve_ok,redo_flag)
! MHP 9/93
         if(no_am_transport_in_core)diffusion_solve_ok = .true.
! IF LOK=T,CONVERGED.
         if(diffusion_solve_ok)goto 200
! IF LREDO=T, A PROBLEM REQUIRES TIMESTEP CUTTING.
         if(redo_flag)goto 9999
  100 continue
      diffusion_solve_ok = .true.
  200 continue
! PERFORM COMPOSITION DIFFUSION OF REMAINING SPECIES.
      am_transport_convective_flag(1) = .true.
      am_transport_convective_flag(num_zones) = .true.
! MHP 6/00 ADDED OPTION OF BS EXTRAPOLATION FOR HE3, CNO
!      LBSMIX=.TRUE.
!      IF(LBSMIX)THEN
!         CALL BSCOMP(DR,DT,COD2,ECOD2,EM,ES1,HRU,HS,HS1,HS2,
!     *               IBEG,IEND,IMIN,LCZ,LOK,M,NTOT,HCOMP,HV,
!     *           HD,HP,HR,HT,MRZONE,MXZONE,NRZONE,NZONE,HSTOT)
!      ELSE
! MHP 08/03 REMOVED OBSOLETE DIFCOM ROUTINE
!      IF(M.GT.1)THEN
      species_begin = 5
      if(use_extended_composition)then
         species_end = 15
      else
         species_end = 11
      endif
      call ndifcom(sub_timestep,mixing_diffusion_coeff, &
           eq_mixing_diffusion_coeff,eq_mass,log_density,luminosity, &
           log_pressure,log_radius,log_mass,enclosed_mass,shell_mass, &
           log_total_mass,diffusion_velocity,zone_begin,zone_end,zone_max, &
           zone_min,am_transport_convective_flag,diffusion_solve_ok, &
           num_zones,composition,species_begin,species_end)
!      ELSE
!      CALL DIFCOM(DR,DT,COD2,ECOD2,EM,ES1,HRU,HS,HS1,HS2,HV,
!     *            IBEG,IEND,IMIN,LCZ,LOK,M,NTOT,HCOMP)
!      ENDIF
! MHP 6/00
      am_transport_convective_flag(1) = lcz_first_zone
      am_transport_convective_flag(num_zones) = lcz_last_zone
! ADD I/O FOR MIXING.
      if(print_diffusion_flag)then
! MHP 8/03 ADDED HEADER LINE.
      write(imodpt,198)
 198  format('COMPOSITION CHANGE FROM ROTATIONAL MIXING'/ &
             5X,'H',8X,'HE3',7X,'C12',7X,'C13',7X,'N14', &
             7X,'LI7',7X,'BE9')
      do i = 1,print_zone_count
         delta_h1 = composition(1,print_zone_id(i))-composition_snapshot(1,print_zone_id(i))
         delta_he3 = composition(4,print_zone_id(i))-composition_snapshot(4,print_zone_id(i))
         delta_c12 = composition(5,print_zone_id(i))-composition_snapshot(5,print_zone_id(i))
         delta_c13 = composition(6,print_zone_id(i))-composition_snapshot(6,print_zone_id(i))
         delta_n14 = composition(7,print_zone_id(i))-composition_snapshot(7,print_zone_id(i))
         delta_li7 = composition(14,print_zone_id(i))-composition_snapshot(14,print_zone_id(i))
         delta_be9 = composition(15,print_zone_id(i))-composition_snapshot(15,print_zone_id(i))
         write(imodpt,199)print_zone_id(i),delta_h1,delta_he3,delta_c12, &
              delta_c13,delta_n14,delta_li7,delta_be9
 199     format(I5,1P7E10.3)
      end do
      endif
 9999 continue
! MHP 8/03 - OMITTED I/O, COULD REINTRODUCE IN ANOTHER FILE
!  DETERMINE COUPLING FACTOR (I.E. THE FRACTION OF THE TOTAL ANGULAR
!  MOMENTUM LOST FROM THE CORE RELATIVE TO ITS FRACTION OF THE TOTAL
!  MOMENT OF INERTIA).
!      DJCORE = 0.0D0
!      HICORE = 0.0D0
!      HJTOT = 0.0D0
!      DO I = 1,IMAX-1
!         DJCORE = DJCORE + (HJM(I)-HJMSAV(I))*HS2(I)
!         HICORE = HICORE + HI(I)
!         HJTOT = HJTOT + HJM(I)*HS2(I)
!      END DO
!      HIENV = 0.0D0
!      DJENV = 0.0D0
!      DO I = IMAX,M
!         DJENV = DJENV + (HJM(I)-HJMSAV(I))*HS2(I)
!         HIENV = HIENV + HI(I)
!         HJTOT = HJTOT + HJM(I)*HS2(I)
!      END DO
!      FI = HICORE/(HICORE+HIENV)
!      IF(ABS(DJCORE+DJENV).GT.1.0D-32)THEN
!         FJ = DJCORE/(DJCORE+DJENV)
!      ELSE
!         FJ = 1.0D0
!      ENDIF
!      WRITE(88,709)SAGE,HJTOT,OMEGA(M),FJ,FJ/FI,
!     *               HCOMP(14,M),HCOMP(15,M)
!  709 FORMAT(1P7E11.3)
      return
end subroutine seculr

!----------------------------------------------------------------------
! microdiff_setup
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original microdiff_setup.f; only variable names, source form, and
! comment style were updated. Validated against the Stage 0 regression
! suite (examples/run_standard_solar_model).
!
! First stage of the microdiff.f90 element-settling pipeline (see also
! microdiff_mte.f90, microdiff_cod.f90, microdiff_run.f90,
! microdiff_etm.f90): locates the diffusion region (between any central
! convective core and any surface convection zone, and bounded by
! hydrogen/helium exhaustion), and converts the needed model quantities
! to Bahcall & Loeb (1990) natural units.
!
!     OUTPUT VARIABLES :
!
!     radius_bl - VECTOR OF UNLOGGED RADII IN BAHCALL AND LOEB UNITS.
!     temperature_bl - VECTOR OF UNLOGGED TEMPERATURES IN BAHCALL AND
!          LOEB UNITS.
!     zone_begin - FIRST ZONE FOR DIFFUSION PURPOSES (EITHER THE FIRST
!          MODEL POINT OR THE OUTERMOST POINT OF A CENTRAL CONVECTION
!          ZONE).
!     zone_end - LAST ZONE FOR DIFFUSION PURPOSES (EITHER THE LAST
!          MODEL POINT OR THE INNERMOST POINT OF A SURFACE CONVECTION
!          ZONE).
!     THE VECTORS dlnp_dr AND enclosed_mass, AND THE SCALARS timestep
!     AND total_mass, ARE ALSO CONVERTED TO BAHCALL AND LOEB UNITS.
!     CONSTANTS DEFINED :
!     ln10 = CONVERSION FACTOR FROM LN TO LOG10
!     CRSUN_BAH = SOLAR RADIUS (CM)
!     CSECYR_BAH = NUMBER OF SECONDS IN A YEAR.
!
! Note: enclosed_mass (originally HS1) is the run of enclosed mass at
! the original model points, unlogged -- naming matches ndifcom.f90/
! bursmix.f90's HS1, not the per-shell mass (HS2 there).
subroutine microdiff_setup(timestep, dlnp_dr, log_radius, log_density, &
     enclosed_mass, log_temperature, convective_flag, num_zones, &
     total_mass, composition, radius_bl, temperature_bl, zone_begin, &
     zone_end, fully_convective_flag, density_orig, temperature_orig)

      use luout_lib
      use const_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(inout) :: timestep
      double precision, intent(inout) :: dlnp_dr(json)
      double precision, intent(in) :: log_radius(json), log_density(json)
      double precision, intent(inout) :: enclosed_mass(json)
      double precision, intent(in) :: log_temperature(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      double precision, intent(inout) :: total_mass
      double precision, intent(in) :: composition(15,json)
      double precision, intent(out) :: radius_bl(json), temperature_bl(json)
      integer, intent(out) :: zone_begin, zone_end
      logical, intent(out) :: fully_convective_flag
      double precision, intent(out) :: density_orig(json), &
           temperature_orig(json)


! common/const/: not used in this file. Naming matches rotgrid.f90/
! vcirc.f90.
      double precision :: clsun, clsunl, clnsun, cmsun, cmsunl, crsun, &
           crsunl, cmbol
      common/const/ clsun, clsunl, clnsun, cmsun, cmsunl, crsun, crsunl, cmbol



! common/confac/: CON_RAD/CON_MASS/CON_TEMP/CON_TIME, the Bahcall &
! Loeb unit-conversion factors set up here and shared with
! microdiff_etm.f90 (which converts back). Naming is local to this
! batch.
      double precision :: bl_radius_scale, bl_mass_scale, bl_temp_scale, &
           bl_time_scale
      common/confac/ bl_radius_scale, bl_mass_scale, bl_temp_scale, &
           bl_time_scale

! common/gravs2/: only hydrogen_diffusion_floor/helium_diffusion_min
! (XMIN/YMIN) are used here. Naming matches mix.f90/rotmix.f90.
      double precision :: settling_timestep_fraction, &
           hydrogen_diffusion_floor, helium_diffusion_min
      logical :: use_thoul_fit
      common/gravs2/ settling_timestep_fraction, hydrogen_diffusion_floor, &
           helium_diffusion_min, use_thoul_fit

! common/scrtch/: not used in this file. Naming matches liburn.f90/
! rotmix.f90.
      double precision :: sesum(json), seg(7,json), sbeta(json), seta(json)
      logical :: locons(json)
      double precision :: so(json), del_grad(3,json), sfxion(3,json), &
           svel(json), scp(json)
      common/scrtch/ sesum, seg, sbeta, seta, locons, so, del_grad, &
           sfxion, svel, scp

! common/gravst/: not used in this file. Naming matches mix.f90/
! rotmix.f90.
      double precision :: settling_tolerance
      integer :: coulomb_log_choice, settling_num_iterations
      logical :: diffuse_helium_active
      common/gravst/ settling_tolerance, coulomb_log_choice, &
           settling_num_iterations, diffuse_helium_active

! common/gravs3/: not used in this file. Naming matches eqstat.f90.
      double precision :: fgry, fgrz
      logical :: lthoul, use_diffusion_z
      common/gravs3/ fgry, fgrz, lthoul, use_diffusion_z

! MHP 8/94 ADDED I/O FOR DIFFUSION
!      COMMON/GSCOF2/TAPP(JSON),TATP(JSON),TCLP(JSON),TAPZP(JSON),
!     *              TATZP(JSON)
! common/gscof/: not used in this file; declared only to preserve
! layout. Not referenced in any already-converted file, so kept as
! lowercased originals pending a confirmed source.
      double precision :: app(json), atp(json), apzp(json), atzp(json)
      common/gscof/ app, atp, apzp, atzp

! FD 10/09 ADDED COMMON BLOCK FOR EXTRA MIXING. IT MIMIC SOME MIXING BY
! AFFECTING THE SETTLING COEFFICIENT DIRECTLY (in setup_grsett.f)
! common/cmixing/: not used in this file; declared only to preserve
! layout. Not referenced in any already-converted file, so kept as
! lowercased originals pending a confirmed source.
      double precision :: cstmixing, cstdiffmix
      common/cmixing/ cstmixing, cstdiffmix

      save

      integer :: i
      double precision :: crsun_bah, csecyr_bah

! DIMENSION dlnp_dr(json),log_radius(json),log_temperature(json), etc.
! all declared above with intent.

      crsun_bah=6.9598d10
      csecyr_bah=3.1558d7
!     fully_convective_flag=T FOR FULLY CONVECTIVE MODEL(AND IF TRUE, DIFFUSION IS SKIPPED).
      fully_convective_flag=.false.
!     CHECK FOR CONVECTIVE CORE.
      if(convective_flag(1))then
         do 10 i=2,num_zones
            if(.not.convective_flag(i))goto 20
   10    continue
!        DIFFUSION NOT COMPUTED FOR FULLY CONVECTIVE MODELS.
         fully_convective_flag=.true.
         write(short_file_unit,15)
   15    format(1x,' FULLY CONVECTIVE MODEL - NO SETTLING')
         goto 9999
   20    continue
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_begin = i-1
      else
         zone_begin = 1
      endif
! MHP 6/90 CHECK FOR HYDROGEN-EXHAUSTED CORE.
      do 23 i = zone_begin,num_zones
         if(composition(1,i).gt.hydrogen_diffusion_floor)goto 25
   23 continue
!     HYDROGEN-FREE MODEL - EXIT.
      write(short_file_unit,16)hydrogen_diffusion_floor
   16 format(1x,'X BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
      fully_convective_flag = .true.
      goto 9999
   25 continue
      zone_begin = i
!     CHECK FOR CONVECTIVE ENVELOPE.
      if(convective_flag(num_zones))then
         do 30 i=num_zones-1,2,-1
            if(.not.convective_flag(i))goto 40
   30    continue
   40    continue
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_end = i+1
      else
         zone_end = num_zones
      endif
!     CHECK FOR HELIUM-EXHAUSTED SURFACE.
!     OUTER POINT IS SET WHEREVER Y>YMIN.
      do 45 i=zone_end,1,-1
         if(composition(2,i).gt.helium_diffusion_min) goto 47
   45 continue
!     HYDROGEN-FREE MODEL - EXIT.
      write(short_file_unit,17)helium_diffusion_min
   17 format(1x,'Y BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
      fully_convective_flag = .true.
      goto 9999
   47 continue
      zone_end = i
!     bl_mass_scale=CONVERSION FACTOR FOR MASS.
!     bl_radius_scale=CONVERSION FACTOR FOR RADIUS.
!     bl_temp_scale=CONVERSION FACOTR FOR TEMPERATURE.
!     bl_time_scale=CONVERSION FACTOR FOR TIME.
      bl_radius_scale=1.0d0/crsun_bah
      bl_mass_scale=1.0d-2*bl_radius_scale**3
      bl_temp_scale=1.0d-7
!     INCLUDES FACTOR OF 2.2 FROM LN LAMBDA
      bl_time_scale=2.7d13*csecyr_bah
!     CONVERT LOG(RADIUS) AND LOG(TEMPERATURE) TO NATURAL UNITS.
!     ALSO CONVERT NATURAL UNITS TO BAHCALL AND LOEB UNITS.
      do 50 i=1,num_zones
         radius_bl(i)=exp(ln10*log_radius(i))*bl_radius_scale
         temperature_bl(i)=exp(ln10*log_temperature(i))*bl_temp_scale
         enclosed_mass(i)=enclosed_mass(i)*bl_mass_scale
         dlnp_dr(i)=dlnp_dr(i)/bl_radius_scale
!        SDEL(2,I)=0.4D0   !COMMENT OUT IN REAL CODE
   50 continue
      timestep=timestep/bl_time_scale
      total_mass=total_mass*bl_mass_scale
!
! COLLECT THE NECESSARY QUANTITIES (NAMELY RHO AND T) FOR LATER
! TRANSFORMATION TO THE EQUALLY SPACED GRID.
      do i = 1,num_zones
         density_orig(i) = exp(ln10*log_density(i))
         temperature_orig(i) = exp(ln10*log_temperature(i))
      enddo
 9999 continue
      return
end subroutine microdiff_setup

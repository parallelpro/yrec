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

      use star_info_lib, only: star, json
      use luout_lib
      use const_lib
      implicit none

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
         do i=2,num_zones
            if(.not.convective_flag(i))exit
         end do
         if (i > num_zones) then
!        DIFFUSION NOT COMPUTED FOR FULLY CONVECTIVE MODELS.
         fully_convective_flag=.true.
         write(short_file_unit,15)
   15    format(1x,' FULLY CONVECTIVE MODEL - NO SETTLING')
         continue
         return
         end if
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_begin = i-1
      else
         zone_begin = 1
      endif
! MHP 6/90 CHECK FOR HYDROGEN-EXHAUSTED CORE.
      do i = zone_begin,num_zones
         if(composition(1,i).gt.star%ctrl%hydrogen_diffusion_floor)exit
      end do
      if (i > num_zones) then
!     HYDROGEN-FREE MODEL - EXIT.
      write(short_file_unit,16)star%ctrl%hydrogen_diffusion_floor
   16 format(1x,'X BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
      fully_convective_flag = .true.
      continue
      return
      end if
      zone_begin = i
!     CHECK FOR CONVECTIVE ENVELOPE.
      if(convective_flag(num_zones))then
         do i=num_zones-1,2,-1
            if(.not.convective_flag(i))exit
         end do
!        COMPUTE OVERSHOOT (TO BE ADDED).
         zone_end = i+1
      else
         zone_end = num_zones
      endif
!     CHECK FOR HELIUM-EXHAUSTED SURFACE.
!     OUTER POINT IS SET WHEREVER Y>YMIN.
      do i=zone_end,1,-1
         if(composition(2,i).gt.star%ctrl%helium_diffusion_min) exit
      end do
      if (i < (1)) then
!     HYDROGEN-FREE MODEL - EXIT.
      write(short_file_unit,17)star%ctrl%helium_diffusion_min
   17 format(1x,'Y BELOW ',f9.6,' IN WHOLE MODEL-NO SETTLING')
      fully_convective_flag = .true.
      continue
      return
      end if
      zone_end = i
!     star%rot%bl_mass_scale=CONVERSION FACTOR FOR MASS.
!     star%rot%bl_radius_scale=CONVERSION FACTOR FOR RADIUS.
!     star%rot%bl_temp_scale=CONVERSION FACOTR FOR TEMPERATURE.
!     star%rot%bl_time_scale=CONVERSION FACTOR FOR TIME.
      star%rot%bl_radius_scale=1.0d0/crsun_bah
      star%rot%bl_mass_scale=1.0d-2*star%rot%bl_radius_scale**3
      star%rot%bl_temp_scale=1.0d-7
!     INCLUDES FACTOR OF 2.2 FROM LN LAMBDA
      star%rot%bl_time_scale=2.7d13*csecyr_bah
!     CONVERT LOG(RADIUS) AND LOG(TEMPERATURE) TO NATURAL UNITS.
!     ALSO CONVERT NATURAL UNITS TO BAHCALL AND LOEB UNITS.
      do i=1,num_zones
         radius_bl(i)=exp(ln10*log_radius(i))*star%rot%bl_radius_scale
         temperature_bl(i)=exp(ln10*log_temperature(i))*star%rot%bl_temp_scale
         enclosed_mass(i)=enclosed_mass(i)*star%rot%bl_mass_scale
         dlnp_dr(i)=dlnp_dr(i)/star%rot%bl_radius_scale
!        SDEL(2,I)=0.4D0   !COMMENT OUT IN REAL CODE
      end do
      timestep=timestep/star%rot%bl_time_scale
      total_mass=total_mass*star%rot%bl_mass_scale
!
! COLLECT THE NECESSARY QUANTITIES (NAMELY RHO AND T) FOR LATER
! TRANSFORMATION TO THE EQUALLY SPACED GRID.
      do i = 1,num_zones
         density_orig(i) = exp(ln10*log_density(i))
         temperature_orig(i) = exp(ln10*log_temperature(i))
      enddo
      return
end subroutine microdiff_setup

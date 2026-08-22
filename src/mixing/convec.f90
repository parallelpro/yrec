!----------------------------------------------------------------------
! convec
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original convec.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! CONVEC determines the boundaries of convective regions with and
! without overshoot.
! *IMPORTANT NOTE*
! In the evolution code, overshoot is assumed to be overmixing:
!    shells in the overshoot region are mixed with the adjacent
!    convection zone, but still use the radiative temperature
!    gradient.
!
! INPUT VARIABLES:
! composition,log_density,log_pressure,log_radius,log_mass,
!   log_temperature: mass fractions of species, and run of density,
!   pressure, radius, mass, and temperature respectively. Used to
!   determine the extent of overshoot regions.
! convective_flag: flag that is true if a shell is convective and
!   false if it is radiative.
! lovstc,envelope_overshoot_active,lovstm (in common/dpmix/): flags set
!   true if overshoot is to be computed for central, surface, and
!   intermediate convection zones respectively.
! num_zones: number of shells in the model.
!
! OUTPUT VARIABLES:
! core_cz_edge, envelope_cz_edge: ids of the edges of the surface and
!   central convection zones respectively. core_cz_edge=1 and
!   envelope_cz_edge=num_zones if the zones in question don't exist.
! mixed_zone_bounds_no_overshoot,mixed_zone_bounds: arrays giving the
!   starting and ending shells in mixed regions. The "_no_overshoot"
!   array refers to the boundaries without overshoot, and
!   mixed_zone_bounds refers to the boundaries with overshoot. The
!   zones are stored in pairs, i.e. elements (3,1) and (3,2) are the
!   first and last mixed shells in the third convection zone out from
!   the center.
! num_mixed_zones,num_mixed_zones_no_overshoot: number of distinct
!   mixed regions with and without overshoot. They can be different
!   because overshoot can cause two nearby convection zones to merge.
! radiative_zone_bounds,num_radiative_zones: the locations of
!   radiative regions and the number of contiguous radiative regions.
!   Needed because mix treats burning in radiative and convective
!   regions differently. These are defined in the same way as the
!   analogous convective variables.
subroutine convec(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, convective_flag, num_zones, &
     radiative_zone_bounds, mixed_zone_bounds, &
     mixed_zone_bounds_no_overshoot, core_cz_edge, envelope_cz_edge, &
     num_radiative_zones, num_mixed_zones, num_mixed_zones_no_overshoot)
      use const_lib
      use luout_lib
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: composition(15,json), &
           log_density(json), log_pressure(json), log_radius(json), &
           log_mass(json), log_temperature(json)
      logical, intent(inout) :: convective_flag(json)
      integer, intent(in) :: num_zones
      integer, intent(out) :: radiative_zone_bounds(13,2)
      integer, intent(inout) :: mixed_zone_bounds(12,2)
      integer, intent(out) :: mixed_zone_bounds_no_overshoot(12,2)
      integer, intent(out) :: core_cz_edge, envelope_cz_edge
      integer, intent(out) :: num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot
      integer :: j_idx, zone_idx, zone_start, k_idx, pair_idx
      logical :: in_convection_zone

! LOCATE BOUNDARIES OF STANDARD CONVECTION ZONES

      j_idx = 1
      in_convection_zone = .false.
      convective_flag(num_zones+1) = .false.
      do zone_idx = 1, num_zones + 1
         if (.not. (.not.convective_flag(zone_idx))) then
! CONVECTION ZONE
         if (in_convection_zone) cycle
! START OF CONVECTION ZONE
         in_convection_zone = .true.
         zone_start = zone_idx
         cycle
         end if
   10    if (.not.in_convection_zone) cycle
!   END OF CONVECTION ZONE
         in_convection_zone = .false.
         if (zone_start.ne.zone_idx-1) then
            mixed_zone_bounds(j_idx,1) = zone_start
            mixed_zone_bounds(j_idx,2) = zone_idx - 1
            j_idx = j_idx + 1
         end if
         if (j_idx.lt.12) cycle
         write(short_file_unit,661)
  661    format(' -----TOO MANY MIX ZONES')
         exit
   11 continue
      end do
   12 continue

! MHP 5/91 LOGIC CHANGE TO AVOID PROBLEMS IF NO CZ IN MODEL(NZONE=0)
! SKIP REST OF SR IF THERE ARE NO CONVECTION ZONES.
      if (j_idx.eq.1) then
         core_cz_edge = 1
         envelope_cz_edge = num_zones
         num_mixed_zones_no_overshoot = 0
         num_mixed_zones = 0
         num_radiative_zones = 1
         radiative_zone_bounds(1,1) = 1
         radiative_zone_bounds(1,2) = num_zones
         continue
         
         return
      end if
      num_mixed_zones = j_idx - 1
      do zone_idx = 1, num_mixed_zones
         mixed_zone_bounds_no_overshoot(zone_idx,1) = &
              mixed_zone_bounds(zone_idx,1)
         mixed_zone_bounds_no_overshoot(zone_idx,2) = &
              mixed_zone_bounds(zone_idx,2)
   20 continue
      end do
      num_mixed_zones_no_overshoot = num_mixed_zones
! FIND OUTER EDGE OF THE CONVECTIVE CORE (JCORE) AND INNER EDGE OF THE
! CONVECTIVE ENVELOPE (JENV) BEFORE OVERSHOOT.
      if (mixed_zone_bounds(1,1).eq.1) then
! CENTRAL CONVECTION ZONE EXISTS IF TRUE.
         if (mixed_zone_bounds(1,2).eq.num_zones) then
! FULLY CONVECTIVE STAR IF TRUE.
            core_cz_edge = 1
            envelope_cz_edge = 1
            num_radiative_zones = 0
            continue
            
            return
         else
            core_cz_edge = mixed_zone_bounds(1,2)
         end if
      else
         core_cz_edge = 1
      end if
      if (mixed_zone_bounds(num_mixed_zones,2).eq.num_zones) then
! SURFACE CONVECTION ZONE EXISTS IF TRUE.
         envelope_cz_edge = mixed_zone_bounds(num_mixed_zones,1)
      else
         envelope_cz_edge = num_zones
      end if
!  ADD CONVECTIVE OVERSHOOT IF NEEDED; THE SIZE OF THE OVERSHOOT REGION IS
!  COMPUTED AND THE EDGES IN MXZONE ARE MOVED TO THE EDGES OF THE
!  OVERSHOOT REGIONS.
      if (.not.lovstc .and. .not.envelope_overshoot_active .and. &
           .not.lovstm) goto 100
      call oversh(composition, log_density, log_pressure, log_radius, &
           log_mass, log_temperature, num_zones, mixed_zone_bounds, &
           mixed_zone_bounds_no_overshoot, num_mixed_zones)
!  CHECK FOR MERGERS OF NEARBY CONVECTION ZONES CAUSED BY OVERSHOOT.
      if (num_mixed_zones.eq.1) goto 100
      j_idx = 1
   85 continue
!  CHECK IF 'TOP' OF ONE REGION IS ABOVE 'BOTTOM' OF THE NEXT ONE.
      if (mixed_zone_bounds(j_idx,2).ge.mixed_zone_bounds(j_idx+1,1)) then
!  IF THIS OCCURS, TWO CONVECTION ZONES HAVE MERGED.
         write(short_file_unit,93) ((mixed_zone_bounds(k_idx,pair_idx), &
              pair_idx=1,2),k_idx=j_idx,j_idx+1), mixed_zone_bounds(j_idx,1), &
              mixed_zone_bounds(j_idx+1,2)
   93    format(2x,'CONVECTION ZONES MERGED DUE TO OVERSHOOT'/2x, &
              'OLD',2('[',i3,'-',i3,']'),' NEW','[',i3,'-',i3,']')
         mixed_zone_bounds(j_idx+1,1) = mixed_zone_bounds(j_idx,1)
         do k_idx = j_idx, num_mixed_zones-1
            do pair_idx = 1, 2
               mixed_zone_bounds(k_idx,pair_idx) = &
                    mixed_zone_bounds(k_idx+1,pair_idx)
   95       continue
            end do
   90    continue
         end do
         num_mixed_zones = num_mixed_zones - 1
         if (j_idx.le.num_mixed_zones-1) then
            goto 85
         else
            goto 100
         end if
      end if
      j_idx = j_idx + 1
      if (j_idx.le.num_mixed_zones-1) goto 85
  100 continue
! NOW DETERMINE THE NUMBER OF RADIATIVE REGIONS.
! CHECK FOR A RADIATIVE REGION BELOW THE FIRST CONVECTION ZONE.
      if (mixed_zone_bounds(1,1).gt.1) then
         num_radiative_zones = 1
         radiative_zone_bounds(1,1) = 1
         radiative_zone_bounds(1,2) = mixed_zone_bounds(1,1) - 1
      else
         num_radiative_zones = 0
      end if
! LOCATE ALL RADIATIVE REGIONS BETWEEN CONVECTION ZONES.
      if (num_mixed_zones.gt.1) then
         do zone_idx = 1, num_mixed_zones-1
            num_radiative_zones = num_radiative_zones + 1
            radiative_zone_bounds(num_radiative_zones,1) = &
                 mixed_zone_bounds(zone_idx,2) + 1
            radiative_zone_bounds(num_radiative_zones,2) = &
                 mixed_zone_bounds(zone_idx+1,1) - 1
  110    continue
         end do
      end if
! CHECK FOR A RADIATIVE REGION ABOVE THE LAST CONVECTION ZONE.
      if (mixed_zone_bounds(num_mixed_zones,2).lt.num_zones) then
         num_radiative_zones = num_radiative_zones + 1
         radiative_zone_bounds(num_radiative_zones,1) = &
              mixed_zone_bounds(num_mixed_zones,2) + 1
         radiative_zone_bounds(1,2) = num_zones
      end if
 9999 continue

      return
end subroutine convec

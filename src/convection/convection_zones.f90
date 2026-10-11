!----------------------------------------------------------------------
! convection_zones
!----------------------------------------------------------------------
! Where the model convects: the convection zones from the shells'
! convective flags, extended by overshoot, the radiative regions
! between them, and their homogenization (formerly convec, oversh,
! hsubp and mixcz). The overshoot switches and extents come from
! convection_settings.
module convection_zones
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use convection_configuration, only: convection_settings
      use eos_lib, only: eos_result, eos_get
      use io_units, only: run_log_unit
      use physical_constants, only: ln10, cgl
      use math_lib
      implicit none
      private
      public :: find_convection_zones, homogenize_convection_zones

contains

!----------------------------------------------------------------------
! find_convection_zones (formerly convec)
!----------------------------------------------------------------------
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
! overshoot_bottom, overshoot_top: the overshoot region below the
!   envelope, updated when overshoot is computed (the caller keeps
!   them; formerly star%iov1/iov2).
! radiative_zone_bounds,num_radiative_zones: the locations of
!   radiative regions and the number of contiguous radiative regions.
!   Needed because mix treats burning in radiative and convective
!   regions differently. These are defined in the same way as the
!   analogous convective variables.
subroutine find_convection_zones(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, convective_flag, num_zones, &
     radiative_zone_bounds, mixed_zone_bounds, &
     mixed_zone_bounds_no_overshoot, core_cz_edge, envelope_cz_edge, &
     num_radiative_zones, num_mixed_zones, num_mixed_zones_no_overshoot, &
     overshoot_bottom, overshoot_top, ierr)

      double precision, intent(in) :: composition(:,:), &
           log_density(:), log_pressure(:), log_radius(:), &
           log_mass(:), log_temperature(:)
      logical, intent(inout) :: convective_flag(:)
      integer, intent(in) :: num_zones
      integer, intent(out) :: radiative_zone_bounds(:,:)
      integer, intent(inout) :: mixed_zone_bounds(:,:)
      integer, intent(out) :: mixed_zone_bounds_no_overshoot(:,:)
      integer, intent(out) :: core_cz_edge, envelope_cz_edge
      integer, intent(out) :: num_radiative_zones, num_mixed_zones, &
           num_mixed_zones_no_overshoot
      integer, intent(inout) :: overshoot_bottom, overshoot_top
      integer, intent(out) :: ierr
      integer :: j_idx, zone_idx, zone_start, k_idx, pair_idx
      logical :: in_convection_zone

! LOCATE BOUNDARIES OF STANDARD CONVECTION ZONES

      ierr = 0
      j_idx = 1
      in_convection_zone = .false.
      convective_flag(num_zones+1) = .false.
      do zone_idx = 1, num_zones + 1
         if (convective_flag(zone_idx)) then
! CONVECTION ZONE
         if (in_convection_zone) cycle
! START OF CONVECTION ZONE
         in_convection_zone = .true.
         zone_start = zone_idx
         cycle
         end if
         if (.not.in_convection_zone) cycle
!   END OF CONVECTION ZONE
         in_convection_zone = .false.
         if (zone_start.ne.zone_idx-1) then
            mixed_zone_bounds(j_idx,1) = zone_start
            mixed_zone_bounds(j_idx,2) = zone_idx - 1
            j_idx = j_idx + 1
         end if
         if (j_idx.lt.size(mixed_zone_bounds, 1)) cycle
         write(run_log_unit,661)
  661    format(' -----TOO MANY MIX ZONES')
         exit
      end do

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
         return
      end if
      num_mixed_zones = j_idx - 1
      do zone_idx = 1, num_mixed_zones
         mixed_zone_bounds_no_overshoot(zone_idx,1) = &
              mixed_zone_bounds(zone_idx,1)
         mixed_zone_bounds_no_overshoot(zone_idx,2) = &
              mixed_zone_bounds(zone_idx,2)
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
      if (convection_settings%core_overshoot .or. convection_settings%envelope_overshoot .or. &
           convection_settings%intermediate_overshoot) then
      call extend_by_overshoot(composition, log_density, log_pressure, log_radius, &
           log_mass, log_temperature, num_zones, mixed_zone_bounds, &
           mixed_zone_bounds_no_overshoot, num_mixed_zones, overshoot_bottom, overshoot_top, ierr)
      if (ierr /= 0) return
!  CHECK FOR MERGERS OF NEARBY CONVECTION ZONES CAUSED BY OVERSHOOT.
      if (num_mixed_zones.ne.1) then
      j_idx = 1
      merge_scan: do
!  CHECK IF 'TOP' OF ONE REGION IS ABOVE 'BOTTOM' OF THE NEXT ONE.
      if (mixed_zone_bounds(j_idx,2).ge.mixed_zone_bounds(j_idx+1,1)) then
!  IF THIS OCCURS, TWO CONVECTION ZONES HAVE MERGED.
         write(run_log_unit,93) ((mixed_zone_bounds(k_idx,pair_idx), &
              pair_idx=1,2),k_idx=j_idx,j_idx+1), mixed_zone_bounds(j_idx,1), &
              mixed_zone_bounds(j_idx+1,2)
   93    format(2x,'CONVECTION ZONES MERGED DUE TO OVERSHOOT'/2x, &
              'OLD',2('[',i3,'-',i3,']'),' NEW','[',i3,'-',i3,']')
         mixed_zone_bounds(j_idx+1,1) = mixed_zone_bounds(j_idx,1)
         do k_idx = j_idx, num_mixed_zones-1
            do pair_idx = 1, 2
               mixed_zone_bounds(k_idx,pair_idx) = &
                    mixed_zone_bounds(k_idx+1,pair_idx)
            end do
         end do
         num_mixed_zones = num_mixed_zones - 1
         if (j_idx.le.num_mixed_zones-1) then
            cycle merge_scan
         else
            exit merge_scan
         end if
      end if
      j_idx = j_idx + 1
      if (j_idx.gt.num_mixed_zones-1) exit merge_scan
      end do merge_scan
      end if
      end if

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
         end do
      end if
! CHECK FOR A RADIATIVE REGION ABOVE THE LAST CONVECTION ZONE.
      if (mixed_zone_bounds(num_mixed_zones,2).lt.num_zones) then
         num_radiative_zones = num_radiative_zones + 1
         radiative_zone_bounds(num_radiative_zones,1) = &
              mixed_zone_bounds(num_mixed_zones,2) + 1
! 2026 audit fix: the original F77 (convec.f: MRZONE(1,2) = M) wrote
! the FIRST region's top here, clobbering it whenever earlier
! radiative regions exist and the outermost CZ is detached from the
! surface. Latent for decades because the top CZ almost always
! reaches the surface (then NRZONE = 1 and the two spellings agree).
         radiative_zone_bounds(num_radiative_zones,2) = num_zones
      end if

      return
end subroutine find_convection_zones

!----------------------------------------------------------------------
! extend_by_overshoot (formerly oversh, overshoot_boundaries)
!----------------------------------------------------------------------
! This routine computes the local pressure scale height at both edges
! of a given convective region and locates the boundaries of overshoot
! regions based on the user-specified extent.
subroutine extend_by_overshoot(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, num_zones, mixed_zone_bounds, &
     mixed_zone_bounds_no_overshoot, num_mixed_zones, overshoot_bottom, overshoot_top, ierr)

      double precision, intent(in) :: composition(:,:), &
           log_density(:), log_pressure(:), log_radius(:), &
           log_mass(:), log_temperature(:)
      integer, intent(in) :: num_zones
      integer, intent(inout) :: mixed_zone_bounds(:,:)
      integer, intent(in) :: mixed_zone_bounds_no_overshoot(:,:)
      integer, intent(in) :: num_mixed_zones
      integer, intent(out) :: overshoot_bottom, overshoot_top
      integer, intent(out) :: ierr
      logical :: up_overshoot_flag, down_overshoot_flag
      integer :: zone_idx, edge_idx, j_idx
      double precision :: pscale_up, pscale_down
      double precision :: cz_radius, overshoot_radius, radius

! overshoot_bottom/top: the overshoot region below the envelope, for
! the adiabatic gradient there (formerly star%iov1/iov2)
      overshoot_bottom = -1
      overshoot_top = -1
      ierr = 0
      do zone_idx = 1, num_mixed_zones
! DETERMINE IF THIS REGION IS A CORE CONVECTION ZONE, SURFACE CZ,
! OR INTERMEDIATE CZ. THERE ARE SEPARATE FLAGS GOVERNING WHETHER
! OVERSHOOT WILL BE PERFORMED IN EACH CASE, AND SEPARATE USER
! PARAMETERS GOVERNING THE DEGREE OF OVERSHOOT.
         if (mixed_zone_bounds(zone_idx,1).eq.1) then
! CONVECTIVE CORE
! CHECK FOR A FULLY CONVECTIVE STAR; SKIP THIS SR IF THERE IS ONE.
            if (mixed_zone_bounds(zone_idx,2).eq.num_zones) then
               write(run_log_unit,5)
    5          format(1x,'FULLY CONVECTIVE MODEL - NO OVERSHOOT')
               write(run_log_unit,200) (mixed_zone_bounds_no_overshoot( &
                    1,j_idx), j_idx=1,2)
               return
            end if
! SKIP IF NO CORE OVERSHOOT IS DESIRED.
            if (.not.convection_settings%core_overshoot) cycle
            up_overshoot_flag = .true.
            down_overshoot_flag = .false.
            edge_idx = mixed_zone_bounds(zone_idx,2)
            call compute_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_up, ierr)
            if (ierr /= 0) return
! PSCALU IS THE PRESSURE SCALE HEIGHT ABOVE THE CONVECTIVE REGION;
! ALPHAC IS THE DESIRED OVERSHOOT (IN SCALE HEIGHTS).
!            PSCALU = PSCALU*ALPHAC
! JVS 07/13 ALLOW FOR THE LIMITING OF OVERSHOOTING ABOVE THE CONVECTIVE
! CORES OF LOVE MASS STARS AS PER WOO & DEMARQUE 2001
            if (convection_settings%limit_core_overshoot) then
               pscale_up = min(pscale_up*convection_settings%overshoot_alpha_core, &
                    convection_settings%core_overshoot_limit*exp(ln10*log_radius(edge_idx)))
            else
               pscale_up = pscale_up*convection_settings%overshoot_alpha_core
            end if

         else if (mixed_zone_bounds(zone_idx,2).eq.num_zones) then
! CONVECTIVE ENVELOPE
! SKIP IF NO ENVELOPE OVERSHOOT IS DESIRED.
            if (.not.convection_settings%envelope_overshoot) cycle
            up_overshoot_flag = .false.
            down_overshoot_flag = .true.
            edge_idx = mixed_zone_bounds(zone_idx,1)
            call compute_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_down, ierr)
            if (ierr /= 0) return
! PSCALD IS THE PRESSURE SCALE HEIGHT BELOW THE CONVECTIVE REGION;
! ALPHAE IS THE DESIRED OVERSHOOT (IN SCALE HEIGHTS).
            pscale_down = pscale_down*convection_settings%overshoot_alpha_envelope
         else
! INTERMEDIATE CONVECTION ZONE (NOT INCLUDING CENTRAL OR SURFACE POINT).
! SKIP IF NO INTERMEDIATE CONVECTION.
            if (.not.convection_settings%intermediate_overshoot) cycle
            up_overshoot_flag = .true.
            down_overshoot_flag = .true.
! PSCALU AND PSCALD HAVE THE SAME MEANING AS ABOVE; OVERSHOOT BOTH BELOW
! AND ABOVE IS PERFORMED BY AN AMOUNT ALPHAM.
            edge_idx = mixed_zone_bounds(zone_idx,1)
            call compute_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_down, ierr)
            if (ierr /= 0) return
            pscale_down = pscale_down*convection_settings%overshoot_alpha_intermediate
            edge_idx = mixed_zone_bounds(zone_idx,2)
            call compute_scale_height(composition, log_density, log_pressure, log_radius, &
                 log_mass, log_temperature, edge_idx, pscale_up, ierr)
            if (ierr /= 0) return
            pscale_up = pscale_up*convection_settings%overshoot_alpha_intermediate
         end if
! COMPUTE EXTENSION OF CONVECTION ZONE BELOW SCHWARTZSCHILD BOUNDARY.
         if (down_overshoot_flag) then
            edge_idx = mixed_zone_bounds(zone_idx,1)
            cz_radius = exp(ln10*log_radius(edge_idx))
            overshoot_radius = cz_radius - pscale_down
! THE OVERSHOOT REGION IS EXTENDED THE RADIAL DISTANCE PSCALD DOWN; THE
! LAST POINT LESS THAN PSCALD FROM THE FORMAL EDGE OF THE CZ IS DEFINED
! AS THE NEW EDGE OF THE MIXED REGION.
            do j_idx = edge_idx-1, 1, -1
               radius = exp(ln10*log_radius(j_idx))
               if (radius.lt.overshoot_radius) exit
            end do
! IF THE LOOP COMPLETES, THE OVERSHOOT REGION EXTENDS BELOW THE FIRST
! POINT AND THE CZ WILL EXTEND TO THE CENTER: natural completion of the
! downward loop already leaves j_idx = 0 (the old explicit assignment
! before label 20 was redundant).
! FOR ROTATING MODELS, ENSURE THAT THERE IS AT LEAST ONE RADIATIVE POINT
! IN THE OVERSHOOT REGION.
            mixed_zone_bounds(zone_idx,1) = j_idx + 1
! 11/91 MHP CHANGED TO REQUIRE AN OVERSHOOT ZONE ONLY IF LINSTB=T.
            if (convection_settings%rotation_active .and. convection_settings%instability_transport_active .and. &
                 mixed_zone_bounds(zone_idx,1).eq. &
                 mixed_zone_bounds_no_overshoot(zone_idx,1)) &
                 mixed_zone_bounds(zone_idx,1) = &
                 mixed_zone_bounds(zone_idx,1) - 1
! DBG 8/94 STORE POSITION OF OVERSHOOT FOR ADIABATIC EXTENSION
            overshoot_top = edge_idx
            overshoot_bottom = mixed_zone_bounds(zone_idx,1)
         end if
! COMPUTE EXTENSION OF CONVECTION ZONE ABOVE SCHWARTZSCHILD BOUNDARY.
         if (up_overshoot_flag) then
            edge_idx = mixed_zone_bounds(zone_idx,2)
            cz_radius = exp(ln10*log_radius(edge_idx))
            overshoot_radius = cz_radius + pscale_up
! THE OVERSHOOT REGION IS EXTENDED THE RADIAL DISTANCE PSCALU UP; THE
! LAST POINT LESS THAN PSCALU FROM THE FORMAL EDGE OF THE CZ IS DEFINED
! AS THE NEW EDGE OF THE MIXED REGION.
            do j_idx = edge_idx+1, num_zones
               radius = exp(ln10*log_radius(j_idx))
               if (radius.gt.overshoot_radius) exit
            end do
! IF THE LOOP COMPLETES, THE OVERSHOOT REGION EXTENDS ABOVE THE LAST
! POINT AND THE CZ WILL EXTEND TO THE SURFACE: natural completion of
! the upward loop already leaves j_idx = num_zones + 1 (the old
! explicit assignment before label 40 was redundant).
            mixed_zone_bounds(zone_idx,2) = j_idx - 1
! 11/91 MHP CHANGED TO REQUIRE AN OVERSHOOT ZONE ONLY IF LINSTB=T.
            if (convection_settings%rotation_active .and. convection_settings%instability_transport_active .and. &
                 mixed_zone_bounds(zone_idx,2).eq. &
                 mixed_zone_bounds_no_overshoot(zone_idx,2)) &
                 mixed_zone_bounds(zone_idx,2) = &
                 mixed_zone_bounds(zone_idx,2) + 1
         end if
      end do
! OUTPUT : THE OLD AND NEW MIXED REGIONS ARE PRINTED OUT IN ISHORT.
      write(run_log_unit,200) ((mixed_zone_bounds_no_overshoot( &
           zone_idx,j_idx), j_idx=1,2), zone_idx=1,num_mixed_zones)
  200 format(1x,'MIXED REGIONS WITHOUT OVERSHOOT', &
           4('[',i4,'-',i4,' ]'))
      write(run_log_unit,210) ((mixed_zone_bounds(zone_idx,j_idx), &
           j_idx=1,2), zone_idx=1,num_mixed_zones)
  210 format(1x,'MIXED REGIONS WITH OVERSHOOT   ', &
           4('[',i4,'-',i4,' ]'))

      return
end subroutine extend_by_overshoot

!----------------------------------------------------------------------
! compute_scale_height (formerly hsubp)
!----------------------------------------------------------------------
! Calculates the local pressure scale height at the edge of a
! convective region (PD Dec. 1984; 2/91 MHP revised for YREC/Mark6).
! Calls the equation of state at the edge zone and returns the
! density-based pressure scale height (pscahe). The inputs are log10
! quantities; the result is in cm.
subroutine compute_scale_height(composition, log_density, log_pressure, &
     log_radius, log_mass, log_temperature, edge_zone, pscahe, ierr)

      double precision, intent(in) :: composition(:,:), log_density(:), &
           log_pressure(:), log_radius(:), log_mass(:), &
           log_temperature(:)
      integer, intent(in) :: edge_zone
      double precision, intent(out) :: pscahe
      integer, intent(out) :: ierr

      logical :: lderiv, latmo
      integer :: ksaha
      double precision :: hydrogen_fraction, metal_fraction, log10_pressure, &
           log10_temperature, log10_density
! 2026 named-index results: the former 20-variable eos output soup is
! one result array (see eos_lib's index constants); only the inout
! i_log10_density (the density consistent with P,T,X,Z) is consumed here.
      type(eos_result) :: eos_res
      double precision :: log10_mass, log10_radius

!  CALL EQUATION OF STATE TO GET THE DENSITY AT THE EDGE ZONE.
      ierr = 0
      lderiv = .false.
      latmo = .true.
      ksaha = 0
      hydrogen_fraction = composition(1,edge_zone)
      metal_fraction = composition(3,edge_zone)
      log10_pressure = log_pressure(edge_zone)
      log10_temperature = log_temperature(edge_zone)
      log10_density = log_density(edge_zone)
      eos_res = eos_result()
      eos_res%log10_density = log10_density
      call eos_get(log10_temperature, log10_pressure, &
           hydrogen_fraction, metal_fraction, eos_res, lderiv, latmo, &
           ksaha, composition_at_zone=composition(:,edge_zone), ierr=ierr)
      if (ierr /= 0) return
!  COMPUTE PRESSURE SCALE HEIGHT  H_P = P R^2/(G M rho).
      log10_mass = log_mass(edge_zone)
      log10_radius = log_radius(edge_zone)
      pscahe = exp(ln10*(log10_pressure + 2.0d0*log10_radius - &
           eos_res%log10_density - cgl - log10_mass))
      return
end subroutine compute_scale_height

!----------------------------------------------------------------------
! homogenize_convection_zones (formerly mixcz)
!----------------------------------------------------------------------
! Locates the boundaries of standard (non-overshoot) convection zones
! from convective_flag and homogenizes the composition within each
! one. Also renormalizes the composition to guard against small
! negative/overflowing abundances at the end.
!
! The convective-overturn-timescale (taucz) calculation that once lived
! here moved to the rotation code; its dummy arguments were dropped from
! the signature (KC 2025-05-31).
subroutine homogenize_convection_zones(composition, shell_mass, convective_flag, num_zones)


      double precision, intent(inout) :: composition(:,:)
      double precision, intent(in) :: shell_mass(:)
      logical, intent(inout) :: convective_flag(:)
      integer, intent(in) :: num_zones

      double precision :: species_sum(15)
      integer :: zone_bounds(24)
      integer :: j_idx, num_zones_plus1, zone_idx, zone_start, zone_end, &
           species_idx, inner_idx, num_species
      double precision :: weight_sum
      logical :: in_convection_zone

! COMPUTE STANDARD CONVECTION ZONES
      j_idx = 1
      in_convection_zone = .false.
      num_zones_plus1 = num_zones + 1
      convective_flag(num_zones_plus1) = .false.
      do zone_idx = 1, num_zones_plus1
         if (convective_flag(zone_idx)) then
! CONVECTION ZONE
         if (in_convection_zone) cycle
! START OF CONVECTION ZONE
         in_convection_zone = .true.
         zone_start = zone_idx
         cycle
         end if
         if (.not.in_convection_zone) cycle
!   END OF CONVECTION ZONE
         in_convection_zone = .false.
         zone_bounds(j_idx) = zone_start
         zone_bounds(j_idx+1) = zone_idx - 1
         j_idx = j_idx + 2
         if (j_idx.lt.24) cycle
         exit
      end do
! TERMINATE THE LIST (ONLY WHEN THE LOOP RAN TO COMPLETION; ON THE
! 12-ZONE OVERFLOW EXIT THE LIST IS LEFT FULL).
      if (zone_idx > num_zones_plus1) then
      zone_bounds(j_idx) = 0
      end if
      num_species = 11
      if (convection_settings%use_extended_composition) num_species = 15
! MIX ALL CONVECTIVE ZONES
      do j_idx = 1, 24, 2
         if (zone_bounds(j_idx).le.0) exit
         zone_start = zone_bounds(j_idx)
         zone_end = min0(num_zones, zone_bounds(j_idx+1))
         if (zone_start.ne.1 .and. zone_start.ge.zone_end) cycle
! INITIALIZE SUMS
         weight_sum = 0.0d0
         do species_idx = 1, num_species
            species_sum(species_idx) = 0.0d0
         end do
         do inner_idx = zone_start, zone_end
            weight_sum = weight_sum + shell_mass(inner_idx)
            do species_idx = 1, num_species
               species_sum(species_idx) = species_sum(species_idx) + &
                    composition(species_idx,inner_idx)*shell_mass(inner_idx)
            end do
         end do
         do species_idx = 1, num_species
            species_sum(species_idx) = species_sum(species_idx)/weight_sum
         end do
         do inner_idx = zone_start, zone_end
            do species_idx = 1, num_species
               composition(species_idx,inner_idx) = species_sum(species_idx)
            end do
         end do
      end do
! RENORMALIZE COMPOSITION IF NECESSARY
      do zone_idx = 1, num_zones
         composition(1,zone_idx) = dmax1(composition(1,zone_idx), 0.0d0)
         composition(3,zone_idx) = dmin1(composition(3,zone_idx), &
              1.0d0-composition(1,zone_idx))
         composition(9,zone_idx) = dmax1(composition(9,zone_idx), &
              0.99d-3*(composition(3,zone_idx)-convection_settings%envelope_metals_no_cno))
      end do

      return
end subroutine homogenize_convection_zones

end module convection_zones

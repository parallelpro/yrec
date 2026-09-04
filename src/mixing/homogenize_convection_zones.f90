!----------------------------------------------------------------------
! mixcz
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mixcz.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Locates the boundaries of standard (non-overshoot) convection zones
! from convective_flag and homogenizes the composition within each
! one. Also renormalizes the composition to guard against small
! negative/overflowing abundances at the end.
!
! The convective-overturn-timescale (taucz) calculation that once lived
! here moved to the rotation code; its dummy arguments were dropped from
! the signature (KC 2025-05-31).
subroutine homogenize_convection_zones(composition, shell_mass, convective_flag, num_zones)

      use star_info_lib, only: star, json
      implicit none

      double precision, intent(inout) :: composition(15,json)
      double precision, intent(in) :: shell_mass(json)
      logical, intent(inout) :: convective_flag(json)
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
      if (star%job%use_extended_composition) num_species = 15
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
              0.99d-3*(composition(3,zone_idx)-star%zenvm))
      end do

      return
end subroutine homogenize_convection_zones

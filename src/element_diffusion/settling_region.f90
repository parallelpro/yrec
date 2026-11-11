!----------------------------------------------------------------------
! settling_region
!----------------------------------------------------------------------
! Where a model settles: the zones zone_begin..zone_end (center first)
! between any convective core and the surface convection zone, bounded
! by hydrogen exhaustion below and helium exhaustion above. Both
! pipelines settle over it; diffuse_elements finds it once per step.
module settling_region
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: i_h1, i_he4
      implicit none
      private
      public :: find_settling_region

! the verdict: a region to settle, or why there is none
      integer, parameter, public :: region_found = 0, &
           region_fully_convective = 1, region_no_hydrogen = 2, &
           region_no_helium = 3

contains

! The settling region of a model of num_zones zones with convective
! flags convective_flag and mass fractions composition(species, zone):
!   zone_begin  the top zone of a convective core (zone 1 without one),
!               moved out to the first zone whose X exceeds
!               hydrogen_floor;
!   zone_end    the base zone of the surface convection zone (num_zones
!               without one), moved in to the first zone whose Y exceeds
!               helium_min.
! region is region_found, or the reason there is no region: the whole
! model convective, or X (from zone_begin out) or Y (from zone_end in)
! nowhere above its floor; zone_begin and zone_end are then not all set.
! (No overshoot widens the convective boundaries.)
pure subroutine find_settling_region(convective_flag, composition, num_zones, &
     hydrogen_floor, helium_min, zone_begin, zone_end, region)
      logical, intent(in) :: convective_flag(:)
      real(dp), intent(in) :: composition(:,:)
      integer, intent(in) :: num_zones
      real(dp), intent(in) :: hydrogen_floor, helium_min
      integer, intent(out) :: zone_begin, zone_end, region
      integer :: i

      zone_begin = 1
      zone_end = num_zones
! the top of a convective core
      if (convective_flag(1)) then
         do i = 2, num_zones
            if (.not. convective_flag(i)) exit
         end do
         if (i > num_zones) then
            region = region_fully_convective
            return
         end if
         zone_begin = i-1
      end if
! out of a hydrogen-exhausted core
      do i = zone_begin, num_zones
         if (composition(i_h1,i) > hydrogen_floor) exit
      end do
      if (i > num_zones) then
         region = region_no_hydrogen
         return
      end if
      zone_begin = i
! the base of the surface convection zone
      if (convective_flag(num_zones)) then
         do i = num_zones-1, 2, -1
            if (.not. convective_flag(i)) exit
         end do
         zone_end = i+1
      end if
! in from a helium-exhausted surface
      do i = zone_end, 1, -1
         if (composition(i_he4,i) > helium_min) exit
      end do
      if (i < 1) then
         region = region_no_helium
         return
      end if
      zone_end = i
      region = region_found
end subroutine find_settling_region

end module settling_region

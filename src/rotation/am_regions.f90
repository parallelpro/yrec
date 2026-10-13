!----------------------------------------------------------------------
! am_regions
!----------------------------------------------------------------------
! The convective and radiative regions for angular momentum transport
! (am_convective_regions): the convection zones with their overshoot.
module am_regions
      implicit none
      private
      public :: am_convective_regions

contains

! ---------------------------------------------------------------
! am_convective_regions: the regions for angular momentum transport and
! mixing. find_convection_zones gives the convection zones extended by
! overshoot (convective_zone_bounds) and the radiative regions between
! them (radiative_zone_bounds); am_transport_convective_flag is set
! .true. in the former and .false. in the latter. Unlike convective_flag
! without overshoot, which homogenize_convection_zones uses, it includes
! the overshoot regions.
subroutine am_convective_regions(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, convective_flag, num_zones, &
     am_transport_convective_flag, radiative_zone_bounds, &
     convective_zone_bounds, num_radiative_zones, num_convective_zones, &
     overshoot_bottom, overshoot_top, ierr)
      use composition_layout, only: json
      use convection_lib, only: find_convection_zones
      implicit none

      double precision, intent(in) :: composition(15,json), &
           log_density(json), log_pressure(json), log_radius(json), &
           log_mass(json), log_temperature(json)
! (find_convection_zones sets convective_flag(num_zones+1) to .false.)
      logical, intent(inout) :: convective_flag(json)
      integer, intent(in) :: num_zones
      logical, intent(out) :: am_transport_convective_flag(json)
      integer, intent(out) :: radiative_zone_bounds(13,2), &
           convective_zone_bounds(12,2)
      integer, intent(out) :: num_radiative_zones, num_convective_zones
! The overshoot region below the envelope (find_convection_zones).
      integer, intent(inout) :: overshoot_bottom, overshoot_top
      integer, intent(out) :: ierr

      integer :: convective_zone_bounds_raw(12,2)
      integer :: i, j, core_zone_boundary, envelope_zone_boundary, &
           num_convective_zones_raw

! The convection zones with and without overshoot, and the radiative
! regions.
      call find_convection_zones(composition,log_density,log_pressure,log_radius, &
           log_mass,log_temperature,convective_flag,num_zones, &
           radiative_zone_bounds,convective_zone_bounds, &
           convective_zone_bounds_raw, &
           core_zone_boundary,envelope_zone_boundary,num_radiative_zones, &
           num_convective_zones,num_convective_zones_raw, overshoot_bottom, overshoot_top, ierr)
      if (ierr /= 0) return
      do j = 1,num_radiative_zones
         do i = radiative_zone_bounds(j,1),radiative_zone_bounds(j,2)
! Skip a non-positive zone index, which would write out of bounds
! (KC 2025-05-30).
            if (i .gt. 0) am_transport_convective_flag(i) = .false.
         end do
      end do
      do j = 1,num_convective_zones
         do i = convective_zone_bounds(j,1),convective_zone_bounds(j,2)
            am_transport_convective_flag(i) = .true.
         end do
      end do

      return
end subroutine am_convective_regions

end module am_regions

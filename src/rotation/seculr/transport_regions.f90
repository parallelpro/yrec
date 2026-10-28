!----------------------------------------------------------------------
! transport_regions
!----------------------------------------------------------------------
! The transport regions of a rotation sub-step: the runs of interfaces
! where the total circulation and instability velocity U is positive.
! The angular momentum transport (secular_transport) and the final
! mixing of the composition (rotational_mixing) both walk them with
! next_transport_region. A star-blind kernel (tools/check_boundaries.py
! STAR_BLIND_FILES), checked by test_rotation_kernels.
module transport_regions
      implicit none
      private
      public :: next_transport_region

contains

! ---------------------------------------------------------------
! The next transport region at or above interface scan_start, up to
! zone_max. Interface k lies between zones k-1 and k, and velocity(k) is
! U there. On return, found says whether the scan met U > 0; zone_begin
! is the zone below the region's first interface with U > 0 and
! zone_end the zone above its last, or zone_max when the region reaches
! it; scan_start is where the next search starts (zone_max + 1 once the
! scan has reached zone_max, which ends the caller's walk). in_region
! says whether a region is open: false before the first call, and left
! true when a region runs up to zone_max.
pure subroutine next_transport_region(velocity, zone_max, scan_start, &
     in_region, zone_begin, zone_end, found)
      double precision, intent(in) :: velocity(:)
      integer, intent(in) :: zone_max
      integer, intent(inout) :: scan_start
      logical, intent(inout) :: in_region
      integer, intent(inout) :: zone_begin, zone_end
      logical, intent(out) :: found
      integer :: k

      found = .false.
      do k = scan_start, zone_max
         if (velocity(k).gt.0.0d0) then
            found = .true.
            if (.not.in_region) then
               in_region = .true.
               zone_begin = k - 1
            end if
         else if (in_region) then
            zone_end = k - 1
            in_region = .false.
            scan_start = k + 1
            exit
         end if
      end do
      if (k.gt.zone_max) then
! The scan reached zone_max: a region still open ends there.
         if (in_region) zone_end = zone_max
         scan_start = zone_max + 1
      end if
end subroutine next_transport_region

end module transport_regions

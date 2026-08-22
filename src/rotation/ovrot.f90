!----------------------------------------------------------------------
! ovrot
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original ovrot.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Isolates convective regions including overshoot (via CONVEC) and
! builds am_transport_convective_flag, the convective flag used for
! angular-momentum-transport/mixing purposes -- i.e. extended to
! include overshoot regions, as opposed to the plain convective_flag
! used elsewhere (see mixcz.f90). Set .false. over each radiative
! region and .true. over each (overshoot-extended) convective region.
subroutine ovrot(composition, log_density, log_pressure, log_radius, &
     log_mass, log_temperature, convective_flag, num_zones, &
     am_transport_convective_flag, radiative_zone_bounds, &
     convective_zone_bounds, num_radiative_zones, num_convective_zones)
      implicit none
      integer, parameter :: json = 5000

      double precision, intent(in) :: composition(15,json), &
           log_density(json), log_pressure(json), log_radius(json), &
           log_mass(json), log_temperature(json)
      logical, intent(in) :: convective_flag(json)
      integer, intent(in) :: num_zones
      logical, intent(out) :: am_transport_convective_flag(json)
      integer, intent(out) :: radiative_zone_bounds(13,2), &
           convective_zone_bounds(12,2)
      integer, intent(out) :: num_radiative_zones, num_convective_zones

      integer :: convective_zone_bounds_raw(12,2)
      save

      integer :: i, j, core_zone_boundary, envelope_zone_boundary, &
           num_convective_zones_raw

! ISOLATE CONVECTIVE REGIONS INCLUDING OVERSHOOT.
      call convec(composition,log_density,log_pressure,log_radius, &
           log_mass,log_temperature,convective_flag,num_zones, &
           radiative_zone_bounds,convective_zone_bounds, &
           convective_zone_bounds_raw, &
           core_zone_boundary,envelope_zone_boundary,num_radiative_zones, &
           num_convective_zones,num_convective_zones_raw)
      do j = 1,num_radiative_zones
         do i = radiative_zone_bounds(j,1),radiative_zone_bounds(j,2)
! KC 2025-05-30 ADDED IF CHECK TO AVOID RUNTIME ERROR.
            if (i .gt. 0) am_transport_convective_flag(i) = .false.
    5    continue
         end do
   10 continue
      end do
      do j = 1,num_convective_zones
         do i = convective_zone_bounds(j,1),convective_zone_bounds(j,2)
            am_transport_convective_flag(i) = .true.
   15    continue
         end do
   20 continue
      end do

      return
end subroutine ovrot

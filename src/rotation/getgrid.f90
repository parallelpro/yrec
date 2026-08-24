!----------------------------------------------------------------------
! getgrid
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getgrid.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This routine forms the sum of star%rot%chi, defined as
! star%rot%chi = log(M)/dm + L/(Ltot*dl) - log(P)/dp, where dm, dp, and dl are
! the point spacings for log mass, luminosity, and pressure
! respectively. This routine transforms these variables to an equally
! spaced grid from the first point to the last point. It returns a
! set of equally spaced star%rot%chi values and their location in mass.
subroutine getgrid(log_luminosity, log_pressure, log_mass, zone_begin, &
     zone_end, num_zones)
      use star_info_lib, only: star, json
      use const_lib
      use numerics_lib
      implicit none

! log_luminosity (originally HL): this file's own comment describes
! this term as "L/(LTOT*DL)" (not "LOG(L)/DL", unlike the LOG(M) and
! LOG(P) terms), but both already-converted callers (mixgrid.f90,
! rotgrid.f90) pass their log_luminosity array into this argument
! position, so it is named accordingly here.
      double precision, intent(in) :: log_luminosity(json), log_pressure(json), &
           log_mass(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones




      double precision :: log_mass_table(json)
      integer :: num_points_in_range, zone_index, range_index
      double precision :: mass_scale, luminosity_scale, pressure_scale

! USE THE MODEL CRITERIA FOR ASSIGNING THE SPACING BETWEEN GRID POINTS.
      mass_scale = chi_grid_scale(2)
      luminosity_scale = chi_grid_scale(9)*log_luminosity(num_zones)
      pressure_scale = chi_grid_scale(11)
      num_points_in_range = zone_end - zone_begin + 1
      do range_index = 1, num_points_in_range
         zone_index = zone_begin + range_index - 1
         star%rot%chi(range_index) = log_mass(zone_index)/mass_scale + &
              log_luminosity(zone_index)/luminosity_scale - &
              log_pressure(zone_index)/pressure_scale
      end do
! DEFINE NTOT EQUAL TO NTAB
      star%rot%ntot = num_points_in_range
! TOTAL NUMBER OF ZONES IS MODULUS OF FINAL CHI PLUS ONE.
!      NTOT = INT(CHI(NTAB)-CHI(1))+1
! FOR ROTATION PURPOSES, DEFINE THE MINIMUM NUMBER OF
! EQUALLY SPACED SHELLS AS 3.
      star%rot%ntot = max(star%rot%ntot,3)
! EQUALLY SPACED INCREMENT IN CHI
      star%rot%dchi = (star%rot%chi(num_points_in_range)-star%rot%chi(1))/dfloat(star%rot%ntot-1)
! ASSIGN VECTOR OF EQUALLY SPACED CHI
      star%rot%echi(1) = star%rot%chi(1)
      do zone_index = 2, star%rot%ntot
         star%rot%echi(zone_index) = star%rot%echi(zone_index-1)+star%rot%dchi
      end do
! NOW ASSIGN MASSES TO THE NEW EQUALLY SPACED GRID POINTS.
! PERFORM INTERPOLATION IN LOG MASS - CONSISTENT WITH CHI DEFINED
! IN LOG M, LOG P, ETC.
      do zone_index = 1, num_points_in_range
         log_mass_table(zone_index) = log_mass(zone_begin+zone_index-1)
      end do
      call osplin(star%rot%echi,star%rot%es1,star%rot%chi,log_mass_table,num_points_in_range,star%rot%ntot)
! TRANSFORM TO PHYSICAL MASS (GM)
      do zone_index = 1, star%rot%ntot
         star%rot%es1(zone_index) = exp(ln10*star%rot%es1(zone_index))
      end do
      return
end subroutine getgrid

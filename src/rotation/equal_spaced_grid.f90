!----------------------------------------------------------------------
! getgrid
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original getgrid.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! This routine forms the sum of rot_scr%chi, defined as
! rot_scr%chi = log(M)/dm + L/(Ltot*dl) - log(P)/dp, where dm, dp, and dl are
! the point spacings for log mass, luminosity, and pressure
! respectively. This routine transforms these variables to an equally
! spaced grid from the first point to the last point. It returns a
! set of equally spaced rot_scr%chi values and their location in mass.
subroutine equal_spaced_grid(log_luminosity, log_pressure, log_mass, zone_begin, &
     zone_end, num_zones)
      use rotation_scratch_lib
      use star_info_lib, only: star, json
      use phys_const_lib
      use numerics_lib
      use math_lib
      implicit none

! log_luminosity (originally HL): this file's own comment describes
! this term as "L/(LTOT*DL)" (not "LOG(L)/DL", unlike the LOG(M) and
! LOG(P) terms), but both already-converted callers (composition_grid.f90,
! am_transport_grid.f90) pass their log_luminosity array into this argument
! position, so it is named accordingly here.
      double precision, intent(in) :: log_luminosity(json), log_pressure(json), &
           log_mass(json)
      integer, intent(in) :: zone_begin, zone_end, num_zones




      double precision :: log_mass_table(json)
      integer :: num_points_in_range, zone_index, range_index
      double precision :: mass_scale, luminosity_scale, pressure_scale

! USE THE MODEL CRITERIA FOR ASSIGNING THE SPACING BETWEEN GRID POINTS.
      mass_scale = star%ctrl%chi_grid_scale(2)
      luminosity_scale = star%ctrl%chi_grid_scale(9)*log_luminosity(num_zones)
      pressure_scale = star%ctrl%chi_grid_scale(11)
      num_points_in_range = zone_end - zone_begin + 1
      do range_index = 1, num_points_in_range
         zone_index = zone_begin + range_index - 1
         rot_scr%chi(range_index) = log_mass(zone_index)/mass_scale + &
              log_luminosity(zone_index)/luminosity_scale - &
              log_pressure(zone_index)/pressure_scale
      end do
! DEFINE NTOT EQUAL TO NTAB
      rot_scr%ntot = num_points_in_range
! TOTAL NUMBER OF ZONES IS MODULUS OF FINAL CHI PLUS ONE.
!      NTOT = INT(CHI(NTAB)-CHI(1))+1
! FOR ROTATION PURPOSES, DEFINE THE MINIMUM NUMBER OF
! EQUALLY SPACED SHELLS AS 3.
      rot_scr%ntot = max(rot_scr%ntot,3)
! EQUALLY SPACED INCREMENT IN CHI
      rot_scr%dchi = (rot_scr%chi(num_points_in_range)-rot_scr%chi(1))/dfloat(rot_scr%ntot-1)
! ASSIGN VECTOR OF EQUALLY SPACED CHI
      rot_scr%echi(1) = rot_scr%chi(1)
      do zone_index = 2, rot_scr%ntot
         rot_scr%echi(zone_index) = rot_scr%echi(zone_index-1)+rot_scr%dchi
      end do
! NOW ASSIGN MASSES TO THE NEW EQUALLY SPACED GRID POINTS.
! PERFORM INTERPOLATION IN LOG MASS - CONSISTENT WITH CHI DEFINED
! IN LOG M, LOG P, ETC.
      do zone_index = 1, num_points_in_range
         log_mass_table(zone_index) = log_mass(zone_begin+zone_index-1)
      end do
      call osplin(rot_scr%echi,rot_scr%es1,rot_scr%chi,log_mass_table,num_points_in_range,rot_scr%ntot)
! TRANSFORM TO PHYSICAL MASS (GM)
      do zone_index = 1, rot_scr%ntot
         rot_scr%es1(zone_index) = exp(ln10*rot_scr%es1(zone_index))
      end do
      return
end subroutine equal_spaced_grid

!----------------------------------------------------------------------
! test_turnover
!----------------------------------------------------------------------
! ok/FAIL checks of the turnover walk (core/diagnostics/turnover_walk.f90,
! time_surface_convection_zone), each confirmed to fail under a mutation
! (test_turnover.mutations.toml). The walk is star-free, so the checks
! run on synthetic profiles of 40 points whose answers are known in
! closed form:
!   r_k = 1e10 + (k-1) 1e8 cm; log P falls by 0.3 dex per point;
!   the pressure scale height H(r) = h0 + h1 r and the convective
!   velocity v(r) = v0 + v1 r are linear in r (the density is P/H at
!   g = 1), so the walk's natural cubic splines reproduce them and
!     tau = H(r*)/v(r*)  at  r* - r_base = H(r*),  r* = (r_base + h0)/(1 - h1),
!   or, for a cell from r_base to r_top, tau = (r_top - r_base)/v at its
!   middle; gradr - grada is -0.3 below a zone and +0.1 at its first
!   point, so the base lies 3/4 of the way (in log r) from the last
!   radiative point to the first convective one, and +0.2 and -0.2 around
!   a cell's top, which lies halfway.
! The cases: a unique zone; the innermost-point branch; two cells on a
! walk with envelope points and on an interior-only walk; a fully
! convective and a fully radiative model; a zone of two points, which is
! skipped; a zone within a dex of a convective core's top, which is
! skipped on a walk with envelope points; a zone whose time exceeds 1e20
! s, which is skipped; a zone at the walk's last two points, with
! convective points beyond the walk that it must not read (bug_sweep
! M329); and a spline on coincident points, which returns the spline's
! status and names the spline.
! Prints one line per check and exits with status 1 when one fails.
program test_turnover
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use turnover_walk, only: turnover_result, time_surface_convection_zone, &
           timing_none, timing_scale_height, timing_bottom_cell, velocity_spline_failed
      use numerics_lib, only: numerics_coincident_points
      implicit none
      integer, parameter :: n = 40
      real(dp), parameter :: tol = 1.0e-10_dp
      real(dp), parameter :: v0 = 1.0e3_dp, v1 = 2.0e-6_dp
      real(dp) :: radius(n), mass(n), log_radius(n), log_pressure(n), log_density(n), &
           gravity(n), velocity(n), gradr(n), grada(n)
      logical :: convective(n)
      type(turnover_result) :: turnover
      real(dp) :: h0, h1, base_radius, top_radius, cell_time, scale_height_time
      integer :: ierr, num_failed = 0

! A unique zone, points 21 to 40, with H = 1.6e8 + 0.02 r (about 4e8 cm).
      h0 = 1.6e8_dp
      h1 = 0.02_dp
      call set_profile()
      call set_zone(21, n)
      call walk(n)
      base_radius = zone_base_radius(21)
      call turnover_check('unique zone: its base point and radius', &
           ierr == 0 .and. turnover%base_index == 21 .and. &
           abs(turnover%base_radius/base_radius - 1.0_dp) < tol)
      call turnover_check('unique zone: tau = H_P/v_c one scale height above the base', &
           turnover%timing == timing_scale_height .and. &
           abs(turnover%turnover_time/scale_height_tau(base_radius) - 1.0_dp) < tol)

! The innermost-point branch: H = 2e7 cm, less than the 2.5e7 cm of the
! zone's first point above its base; the zone is points 21 to 30. The
! walk times the zone at that point, tau = H_P/v_c there, and not by the
! bottom cell's crossing time, which differs from it by more than 10%
! (bug_sweep M324).
      h0 = 2.0e7_dp
      h1 = 0.0_dp
      call set_profile()
      call set_zone(21, 30)
      call walk(n)
      base_radius = zone_base_radius(21)
      top_radius = zone_top_radius(30)
      cell_time = cell_tau(base_radius, top_radius)
      scale_height_time = h0/velocity(21)
      call turnover_check('innermost point above H_P: H_P/v_c at that point (M324)', &
           ierr == 0 .and. turnover%timing == timing_scale_height .and. &
           abs(turnover%turnover_time/scale_height_time - 1.0_dp) < tol .and. &
           abs(cell_time/scale_height_time - 1.0_dp) > 0.1_dp)

! Two cells, points 11 to 15 and 26 to 40. With envelope points in the
! walk (30 interior points) the lower cell is timed by its crossing.
      h0 = 1.6e8_dp
      h1 = 0.02_dp
      call set_profile()
      call set_zone(11, 15)
      call set_zone(26, n)
      call walk(30)
      base_radius = zone_base_radius(11)
      top_radius = zone_top_radius(15)
      call turnover_check('two cells, envelope walk: the lower cell''s width over v_c at its middle', &
           ierr == 0 .and. turnover%base_index == 11 .and. &
           turnover%timing == timing_bottom_cell .and. &
           abs(turnover%turnover_time/cell_tau(base_radius, top_radius) - 1.0_dp) < tol)
! An interior-only walk takes the upper zone and times it as a unique one.
      call walk(n)
      base_radius = zone_base_radius(26)
      call turnover_check('two cells, interior walk: the upper zone, one scale height above its base', &
           ierr == 0 .and. turnover%base_index == 26 .and. &
           turnover%timing == timing_scale_height .and. &
           abs(turnover%turnover_time/scale_height_tau(base_radius) - 1.0_dp) < tol)

! A fully convective model, with H = 1.05e10 + 0.02 r: the base is the
! centre, at 10 cm.
      h0 = 1.05e10_dp
      h1 = 0.02_dp
      call set_profile()
      convective = .true.
      gradr = grada + 0.1_dp
      call walk(n)
      call turnover_check('fully convective: the base is the centre (10 cm), timed above it', &
           ierr == 0 .and. turnover%base_index == 1 .and. &
           abs(turnover%base_radius/10.0_dp - 1.0_dp) < tol .and. &
           turnover%timing == timing_scale_height .and. &
           abs(turnover%turnover_time/scale_height_tau(10.0_dp) - 1.0_dp) < tol)

! A fully radiative model.
      call set_profile()
      call walk(n)
      call turnover_check('fully radiative: no zone, tau = 0 and base radius 0', &
           ierr == 0 .and. turnover%timing == timing_none .and. &
           turnover%turnover_time == 0.0_dp .and. turnover%base_radius == 0.0_dp)

! A zone of two points, 11 and 12, below the zone 26 to 40: it is
! skipped, on a walk with envelope points too.
      h0 = 1.6e8_dp
      h1 = 0.02_dp
      call set_profile()
      call set_zone(11, 12)
      call set_zone(26, n)
      call walk(30)
      base_radius = zone_base_radius(26)
      call turnover_check('a zone of two points is skipped', &
           ierr == 0 .and. turnover%base_index == 26 .and. &
           turnover%timing == timing_scale_height .and. &
           abs(turnover%turnover_time/scale_height_tau(base_radius) - 1.0_dp) < tol)

! A convective core, points 1 to 5, and a zone 8 to 12 whose base is 0.9
! dex in pressure above the core's top: on a walk with envelope points it
! is skipped for the zone 26 to 40, 4.2 dex above the top the scan then
! leaves at point 12.
      call set_profile()
      convective(1:5) = .true.
      gradr(1:5) = grada(1:5) + 0.1_dp
      call set_zone(8, 12)
      call set_zone(26, n)
      call walk(30)
      call turnover_check('a zone within a dex of a convective core''s top is skipped', &
           ierr == 0 .and. turnover%base_index == 26 .and. &
           turnover%timing == timing_scale_height)

! A lower cell, 11 to 15, whose convective velocity is 1e-15 cm/s: its
! crossing time is above 1e20 s, and the walk goes on to the zone above.
      call set_profile()
      call set_zone(11, 15)
      call set_zone(26, n)
      velocity(9:17) = 1.0e-15_dp
      call walk(30)
      base_radius = zone_base_radius(26)
      call turnover_check('a zone slower than 1e20 s is skipped', &
           ierr == 0 .and. turnover%base_index == 26 .and. &
           abs(turnover%turnover_time/scale_height_tau(base_radius) - 1.0_dp) < tol)

! A walk of 38 of the 40 points whose last two are convective, as the
! two beyond it are: the zone has two points, not four, and is no zone.
      call set_profile()
      call set_zone(37, n)
      call time_surface_convection_zone(mass, convective, log_radius, log_pressure, &
           log_density, gravity, velocity, gradr, grada, 38, 38, turnover, ierr)
      call turnover_check('a zone at the walk''s last points: nothing beyond the walk is read', &
           ierr == 0 .and. turnover%timing == timing_none .and. &
           turnover%turnover_time == 0.0_dp)

! Four convective points, the first two at one radius, with H above the
! model's radius: no point lies one H_P above the base, and the bottom
! cell's spline in r meets the coincident points below its middle.
      call set_profile()
      radius(1:4) = [8.0e9_dp, 8.0e9_dp, 9.0e9_dp, 1.0e10_dp]
      log_radius(1:4) = log10(radius(1:4))
      log_density(1:4) = log_pressure(1:4) - 11.0_dp
      convective(1:4) = .true.
      gradr(1:4) = grada(1:4) + [0.2_dp, 0.2_dp, 0.2_dp, 0.1_dp]
      call time_surface_convection_zone(mass(1:4), convective(1:4), log_radius(1:4), &
           log_pressure(1:4), log_density(1:4), gravity(1:4), velocity(1:4), gradr(1:4), &
           grada(1:4), 4, 4, turnover, ierr)
      call turnover_check('coincident points: the spline''s status, and which spline', &
           ierr == numerics_coincident_points .and. &
           turnover%failed_spline == velocity_spline_failed)

      if (num_failed > 0) then
         write(*,'(a,i0,a)') 'test_turnover: ', num_failed, ' checks FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_turnover: done'

contains

! ---------------------------------------------------------------
! The radiative profile for the current h0 and h1: no convective point,
! gradr - grada = -0.3 everywhere.
subroutine set_profile()
      integer :: k
      do k = 1, n
         radius(k) = 1.0e10_dp + real(k - 1, dp)*1.0e8_dp
         mass(k) = real(k, dp)*1.0e30_dp
         log_pressure(k) = 17.0_dp - 0.3_dp*real(k - 1, dp)
      end do
      log_radius = log10(radius)
      log_density = log_pressure - log10(h0 + h1*radius)
      gravity = 1.0_dp
      velocity = v0 + v1*radius
      grada = 0.4_dp
      gradr = grada - 0.3_dp
      convective = .false.
end subroutine set_profile

! ---------------------------------------------------------------
! Points first to last convective: gradr - grada = +0.1 at the first
! point and +0.2 at the others (and at the last), over -0.3 below and
! -0.2 at the point above.
subroutine set_zone(first, last)
      integer, intent(in) :: first, last
      convective(first:last) = .true.
      gradr(first:last) = grada(first:last) + 0.2_dp
      gradr(first) = grada(first) + 0.1_dp
      if (last < n) gradr(last+1) = grada(last+1) - 0.2_dp
end subroutine set_zone

! ---------------------------------------------------------------
! The walk on the n points, the first num_interior of them the interior's.
subroutine walk(num_interior)
      integer, intent(in) :: num_interior
      call time_surface_convection_zone(mass, convective, log_radius, log_pressure, &
           log_density, gravity, velocity, gradr, grada, n, num_interior, turnover, ierr)
end subroutine walk

! ---------------------------------------------------------------
! The radius of the base of a zone whose first point is `first`: the
! zero of gradr - grada, 3/4 of the way in log r from the point below.
real(dp) function zone_base_radius(first)
      integer, intent(in) :: first
      zone_base_radius = 10.0_dp**(log_radius(first-1) + &
           0.75_dp*(log_radius(first) - log_radius(first-1)))
end function zone_base_radius

! The radius of the top of a cell whose last point is `last`: halfway in
! log r to the point above.
real(dp) function zone_top_radius(last)
      integer, intent(in) :: last
      zone_top_radius = 10.0_dp**(log_radius(last) + &
           0.5_dp*(log_radius(last+1) - log_radius(last)))
end function zone_top_radius

! tau = H(r*)/v(r*) where r* - r_base = H(r*).
real(dp) function scale_height_tau(r_base)
      real(dp), intent(in) :: r_base
      real(dp) :: r_star
      r_star = (r_base + h0)/(1.0_dp - h1)
      scale_height_tau = (h0 + h1*r_star)/(v0 + v1*r_star)
end function scale_height_tau

! tau = (r_top - r_base)/v at the cell's middle.
real(dp) function cell_tau(r_base, r_top)
      real(dp), intent(in) :: r_base, r_top
      cell_tau = (r_top - r_base)/(v0 + v1*0.5_dp*(r_base + r_top))
end function cell_tau

! ---------------------------------------------------------------
! One check: its label and ok, or FAIL.
subroutine turnover_check(label, condition)
      character(len=*), intent(in) :: label
      logical, intent(in) :: condition
      character(len=96) :: padded
      padded = label
      if (condition) then
         write(*,'(a,a)') padded, '  ok'
      else
         write(*,'(a,a)') padded, '  FAIL'
         num_failed = num_failed + 1
      end if
end subroutine turnover_check

end program test_turnover

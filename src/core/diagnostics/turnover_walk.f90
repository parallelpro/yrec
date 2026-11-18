!----------------------------------------------------------------------
! turnover_walk
!----------------------------------------------------------------------
! The convective turnover time near the base of the surface convection
! zone, on a run of points from the centre outward: a star-blind kernel
! (docs/GUIDELINES.md, star rule 6). compute_turnover_timescale
! (turnover_timescale.f90) gathers the points from the stitched model,
! calls time_surface_convection_zone, reports a failure and stores the
! result; nothing here reads the star, prints or keeps state between
! calls.
!
! The method (G. Somers 3/17). The search runs outward from the centre,
! past a convective core, to the base of the next convection zone, and
! inward from the surface to the base of the uppermost one.
! - Both find the same base: the zone is unique, and its turnover time is
!   tau = H_P/v_c one pressure scale height above the base, with
!   H_P = P/(rho g) and v_c the convective velocity, both splined to the
!   point where r - r_base = H_P (time_scale_height_above_base).
! - They differ (several convection cells lie in the envelope), or one
!   H_P above the base lies outside the points: tau is the crossing time
!   of the lowest cell, its width over v_c at its middle
!   (time_bottom_cell_crossing).
! A zone found is rejected, and the search restarts above its base, when
! it is thinner than three points, when on a walk with envelope points it
! lies within min_core_separation_dex in log P of a convective core's
! top, or when its turnover time exceeds max_turnover_time.
module turnover_walk
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      use physical_constants, only: ln10
      use numerics_lib, only: natural_spline, splint
      use convection_geometry, only: cz_base_fraction
      implicit none
      private
      public :: turnover_result, time_surface_convection_zone
      public :: timing_none, timing_scale_height, timing_bottom_cell
      public :: no_spline_failed, velocity_spline_failed, scale_height_spline_failed

! How the turnover time was found (turnover_turnover%timing): there is no
! surface convection zone; H_P/v_c one scale height above the base; the
! bottom cell's crossing time.
      integer, parameter :: timing_none = 0, timing_scale_height = 1, &
           timing_bottom_cell = 2
! The spline that failed (turnover_turnover%failed_spline), for the
! caller's message.
      integer, parameter :: no_spline_failed = 0, velocity_spline_failed = 1, &
           scale_height_spline_failed = 2

! What the walk returns:
!   turnover_time  the convective turnover time [s]; 0 without a surface
!                  convection zone;
!   base_radius    the radius of the zone's base [cm]: 10 cm for a fully
!                  convective model, 0 without a surface convection zone;
!   base_index     the innermost point of the zone timed (the last point
!                  for a radiative envelope);
!   timing         how the time was found (timing_*);
!   failed_spline  which spline failed when ierr is not 0 (*_spline_failed).
      type :: turnover_result
         real(dp) :: turnover_time = 0.0d0
         real(dp) :: base_radius = 0.0d0
         integer :: base_index = 0
         integer :: timing = timing_none
         integer :: failed_spline = no_spline_failed
      end type turnover_result

! The zone the search holds from one pass to the next:
!   core_top_index    the top of the convective region the scan left
!                     behind it (0 until a scan starts in one);
!   base_index        the base found from below (legacy IMAX);
!   upper_base_index  the base of the uppermost zone (legacy JMAX);
!   index_gap         upper_base_index - base_index: 0 for a unique zone;
!   fully_convective  the scan ran off the top inside a convective region;
!   has_mass          the surface convection zone holds mass (legacy
!                     LCZSUR).
      type :: zone_search
         integer :: core_top_index = 0
         integer :: base_index = 0
         integer :: upper_base_index = 0
         integer :: index_gap = 0
         logical :: fully_convective = .false.
         logical :: has_mass = .false.
      end type zone_search

! A zone whose base lies less than this in log P above the top of a
! convective core is a few convective shells on the core, not the surface
! zone [dex]; a turnover time above max_turnover_time means the zone
! found lies deep in the interior [s]; a zone is at least
! min_zone_points thick; the centre of a fully convective model is put at
! log r = central_log_radius, 10 cm [log10 cm]; the splines take
! stencil_size points.
      real(dp), parameter :: min_core_separation_dex = 1.0d0, &
           max_turnover_time = 1.0d20, central_log_radius = 1.0d0
      integer, parameter :: min_zone_points = 3, stencil_size = 4

contains

!----------------------------------------------------------------------
! time_surface_convection_zone (legacy tauintnew.f; G. Somers 3/17)
!----------------------------------------------------------------------
! The base of the surface convection zone and its turnover time, on the
! points 1 (centre) to num_points, the first num_interior_points of them
! the interior model's and the rest the envelope's. Per point: the mass
! inside it shell_mass (any unit), the convective flag, log10 of r (cm),
! P (dyn/cm^2) and rho (g/cm^3), the local gravity G m/r^2 (cm/s^2), the
! convective velocity (cm/s) and the radiative and adiabatic gradients.
! The result is described at turnover_result; ierr is the spline's
! status (numerics_status) when one fails, with turnover%failed_spline.
subroutine time_surface_convection_zone(shell_mass, convective_flag, log10_radius, &
     log10_pressure, log10_density, local_gravity, convective_velocity, &
     radiative_gradient, adiabatic_gradient, num_points, num_interior_points, &
     turnover, ierr)
      real(dp), intent(in) :: shell_mass(:), log10_radius(:), log10_pressure(:), &
           log10_density(:), local_gravity(:), convective_velocity(:), &
           radiative_gradient(:), adiabatic_gradient(:)
      logical, intent(in) :: convective_flag(:)
      integer, intent(in) :: num_points, num_interior_points
      type(turnover_result), intent(out) :: turnover
      integer, intent(out) :: ierr
      type(zone_search) :: zone
      integer :: search_start_index, scan_index
      logical :: scan_ended_radiative, timed

      ierr = 0
      search_start_index = 1
      search: do
         call scan_for_cz_base_from_below(convective_flag, num_points, &
              search_start_index, scan_index, scan_ended_radiative, zone%core_top_index)
         if (scan_index > num_points .and. scan_ended_radiative) then
!           A radiative envelope (perhaps over a convective core).
            zone%fully_convective = .false.
            zone%has_mass = .false.
            zone%base_index = num_points
            zone%upper_base_index = num_points
            zone%core_top_index = num_points
         else
            call find_cz_bases(convective_flag, num_points, scan_index, zone)
            if (zone_hugs_convective_core(log10_pressure, num_points, &
                 num_interior_points, zone)) then
               search_start_index = zone%base_index + 1
               cycle search
            end if
            if (.not. zone_is_thick_enough(convective_flag, zone%base_index)) then
!              Likely a single convective shell deep in the interior:
!              restart above it.
               if (.not. convective_flag(zone%base_index+1)) then
                  search_start_index = zone%base_index + 1
               else
                  search_start_index = zone%base_index + 2
               end if
               cycle search
            end if
            call assess_surface_cz_depth(shell_mass, num_points, &
                 num_interior_points, zone)
         end if

         if (zone%has_mass) then
            turnover%base_radius = cz_base_radius(log10_radius, radiative_gradient, &
                 adiabatic_gradient, zone)
!           A unique zone is timed one H_P above its base. Several cells
!           in a mostly radiative envelope, or no point one H_P above the
!           base, take the bottom cell's crossing time.
            timed = .false.
            if (zone%index_gap == 0) then
               call time_scale_height_above_base(log10_radius, log10_pressure, &
                    log10_density, local_gravity, convective_velocity, num_points, &
                    zone%base_index, turnover%base_radius, turnover%turnover_time, timed, &
                    turnover%failed_spline, ierr)
               if (ierr /= 0) return
            end if
            if (timed) then
               turnover%timing = timing_scale_height
            else
               call time_bottom_cell_crossing(convective_flag, log10_radius, &
                    convective_velocity, radiative_gradient, adiabatic_gradient, &
                    num_points, zone%base_index, turnover%base_radius, &
                    turnover%turnover_time, turnover%failed_spline, ierr)
               if (ierr /= 0) return
               turnover%timing = timing_bottom_cell
            end if
         else
            turnover%turnover_time = 0.0d0
            turnover%base_radius = 0.0d0
            turnover%timing = timing_none
         end if
         turnover%base_index = zone%base_index
!        A zone deep in the interior: search again above its base.
         if (turnover%turnover_time > max_turnover_time) then
            search_start_index = zone%base_index + 1
            cycle search
         end if
         exit search
      end do search
end subroutine time_surface_convection_zone

! ---------------------------------------------------------------
! Scans outward from search_start_index, past the convective region it
! starts in (whose top goes to core_top_index; it is left alone when the
! scan starts at a radiative point), to the base of the next convection
! zone, scan_index; num_points + 1 when there is none.
! scan_ended_radiative: the scan's last point was radiative.
pure subroutine scan_for_cz_base_from_below(convective_flag, num_points, &
     search_start_index, scan_index, scan_ended_radiative, core_top_index)
      logical, intent(in) :: convective_flag(:)
      integer, intent(in) :: num_points, search_start_index
      integer, intent(out) :: scan_index
      logical, intent(out) :: scan_ended_radiative
      integer, intent(inout) :: core_top_index
      integer :: i

      scan_ended_radiative = .not. convective_flag(search_start_index)
      do i = search_start_index, num_points
         if (.not. scan_ended_radiative .and. .not. convective_flag(i)) then
            core_top_index = i - 1
            scan_ended_radiative = .true.
         end if
         if (scan_ended_radiative .and. convective_flag(i)) exit
      end do
      scan_index = i
end subroutine scan_for_cz_base_from_below

! ---------------------------------------------------------------
! The base found from below (the scan's) and the base of the uppermost
! zone: the outermost point j + 1 that is convective over a radiative
! point j, so that outer points that are already radiative are skipped
! (JvS 02/2026). A scan that ran off the top inside a convective region
! is a fully convective model: both bases are the centre, and the zone
! holds mass.
pure subroutine find_cz_bases(convective_flag, num_points, scan_index, zone)
      logical, intent(in) :: convective_flag(:)
      integer, intent(in) :: num_points, scan_index
      type(zone_search), intent(inout) :: zone
      integer :: j

      if (scan_index > num_points) then
         zone%fully_convective = .true.
         zone%has_mass = .true.
         zone%core_top_index = 1
         zone%base_index = 1
         zone%upper_base_index = 1
      else
         do j = num_points-1, 1, -1
            if (convective_flag(j+1) .and. .not. convective_flag(j)) exit
         end do
         zone%base_index = scan_index
         zone%upper_base_index = j + 1
      end if
      zone%index_gap = zone%upper_base_index - zone%base_index
end subroutine find_cz_bases

! ---------------------------------------------------------------
! Whether the zone found is a few convective shells just above a
! convective core: on a walk with envelope points, its base less than
! min_core_separation_dex in log P above the core's top. (The test on
! positive indices keeps both reads inside the arrays; KC 2025-05-30.)
pure logical function zone_hugs_convective_core(log10_pressure, num_points, &
     num_interior_points, zone)
      real(dp), intent(in) :: log10_pressure(:)
      integer, intent(in) :: num_points, num_interior_points
      type(zone_search), intent(in) :: zone

      zone_hugs_convective_core = .false.
      if (zone%base_index > 0 .and. zone%core_top_index > 0) then
         zone_hugs_convective_core = num_points /= num_interior_points .and. &
              log10_pressure(zone%core_top_index) - log10_pressure(zone%base_index) &
              < min_core_separation_dex
      end if
end function zone_hugs_convective_core

! ---------------------------------------------------------------
! Whether the zone is at least min_zone_points thick: the base point and
! the points above it, min_zone_points in all, are convective.
pure logical function zone_is_thick_enough(convective_flag, base_index)
      logical, intent(in) :: convective_flag(:)
      integer, intent(in) :: base_index

      zone_is_thick_enough = all(convective_flag(base_index:base_index+min_zone_points-1))
end function zone_is_thick_enough

! ---------------------------------------------------------------
! An interior-only walk takes the uppermost zone's base and treats the
! zone as unique. Then whether the surface convection zone holds mass:
! the mass at its top (the outermost point) exceeds the mass at its base
! (midway between the uppermost zone's base point and the one below; 0 at
! the centre). The threshold is any mass.
pure subroutine assess_surface_cz_depth(shell_mass, num_points, num_interior_points, &
     zone)
      real(dp), intent(in) :: shell_mass(:)
      integer, intent(in) :: num_points, num_interior_points
      type(zone_search), intent(inout) :: zone
      real(dp) :: mass_at_cz_top, mass_at_cz_bottom

      if (zone%index_gap /= 0 .and. num_points == num_interior_points) then
         zone%base_index = zone%upper_base_index
         zone%index_gap = 0
      end if
      mass_at_cz_top = shell_mass(num_points)
      if (zone%upper_base_index > 1) then
         mass_at_cz_bottom = 0.5d0*(shell_mass(zone%upper_base_index) + &
              shell_mass(zone%upper_base_index-1))
      else
         mass_at_cz_bottom = 0.0d0
      end if
      zone%has_mass = mass_at_cz_top - mass_at_cz_bottom > 0.0d0
end subroutine assess_surface_cz_depth

! ---------------------------------------------------------------
! The radius of the zone's base [cm]: linear in log r to the zero of
! gradr - grada between the last radiative point and the base point,
! f = d_below/(d_below - d) (cz_base_fraction); the centre, 10 cm, for a
! fully convective model.
real(dp) function cz_base_radius(log10_radius, radiative_gradient, &
     adiabatic_gradient, zone)
      real(dp), intent(in) :: log10_radius(:), radiative_gradient(:), &
           adiabatic_gradient(:)
      type(zone_search), intent(in) :: zone
      real(dp) :: gap_below, gap_at_base, base_fraction, log_base_radius
      integer :: k

      if (zone%fully_convective) then
         log_base_radius = central_log_radius
      else
         k = zone%base_index
         gap_below = radiative_gradient(k-1) - adiabatic_gradient(k-1)
         gap_at_base = radiative_gradient(k) - adiabatic_gradient(k)
         base_fraction = cz_base_fraction(gap_below, gap_at_base)
         log_base_radius = log10_radius(k-1) + base_fraction* &
              (log10_radius(k) - log10_radius(k-1))
      end if
      cz_base_radius = exp(ln10*log_base_radius)
end function cz_base_radius

! ---------------------------------------------------------------
! The pressure scale height H_P = P/(rho g) at point k [cm].
pure real(dp) function scale_height_at(log10_pressure, log10_density, local_gravity, k)
      real(dp), intent(in) :: log10_pressure(:), log10_density(:), local_gravity(:)
      integer, intent(in) :: k

      scale_height_at = exp(ln10*(log10_pressure(k) - log10_density(k)))/local_gravity(k)
end function scale_height_at

! ---------------------------------------------------------------
! The height r - r_base of point k above the zone's base [cm].
pure real(dp) function height_above_base(log10_radius, base_radius, k)
      real(dp), intent(in) :: log10_radius(:), base_radius
      integer, intent(in) :: k

      height_above_base = exp(ln10*log10_radius(k)) - base_radius
end function height_above_base

! ---------------------------------------------------------------
! The first point of a four-point spline stencil around point k: the
! stencil k-2..k+1, moved one point outward when k is the second point of
! the zone, so that it stays above the base, and kept inside the points.
pure integer function stencil_start(k, base_index, num_points)
      integer, intent(in) :: k, base_index, num_points
      integer :: third_point

      third_point = k
      if (k == base_index + 1) third_point = k + 1
      third_point = max(3, min(third_point, num_points-1))
      stencil_start = third_point - 2
end function stencil_start

! ---------------------------------------------------------------
! The turnover time one pressure scale height above the base,
! tau = H_P/v_c, with H_P and v_c splined through four points to the
! place where r - r_base = H_P. timed is set when it is found: it is not
! when no point lies more than its own H_P above the base, and not when
! the base point itself does. That innermost-point case sets
! tau = H_P/v_c at the base point and leaves timed false, so that
! time_bottom_cell_crossing replaces the value (bug_sweep M324).
subroutine time_scale_height_above_base(log10_radius, log10_pressure, log10_density, &
     local_gravity, convective_velocity, num_points, base_index, base_radius, &
     turnover_time, timed, failed_spline, ierr)
      real(dp), intent(in) :: log10_radius(:), log10_pressure(:), log10_density(:), &
           local_gravity(:), convective_velocity(:)
      integer, intent(in) :: num_points, base_index
      real(dp), intent(in) :: base_radius
      real(dp), intent(inout) :: turnover_time
      logical, intent(inout) :: timed
      integer, intent(inout) :: failed_spline
      integer, intent(out) :: ierr
! The stencil's nodes x = r - r_base - H_P and the values on them.
      real(dp) :: excess_height(stencil_size), velocity(stencil_size), &
           scale_height(stencil_size)
      real(dp) :: velocity_at_scale_height, scale_height_there
      integer :: k, first, n

      ierr = 0
      if (height_above_base(log10_radius, base_radius, base_index) - &
           scale_height_at(log10_pressure, log10_density, local_gravity, base_index) > 0.0d0) then
         call time_innermost_point(log10_pressure, log10_density, local_gravity, &
              convective_velocity, base_index, turnover_time)
         return
      end if
      do k = base_index+1, num_points-1
!        The first point where r - r_base exceeds the local H_P.
         if (height_above_base(log10_radius, base_radius, k) - &
              scale_height_at(log10_pressure, log10_density, local_gravity, k) > 0.0d0) then
            first = stencil_start(k, base_index, num_points)
            do n = 1, stencil_size
               scale_height(n) = scale_height_at(log10_pressure, log10_density, &
                    local_gravity, first+n-1)
               excess_height(n) = height_above_base(log10_radius, base_radius, first+n-1) - &
                    scale_height(n)
               velocity(n) = convective_velocity(first+n-1)
            end do
            call spline_stencil_value(excess_height, velocity, 0.0d0, &
                 velocity_at_scale_height, ierr)
            if (ierr /= 0) then
               failed_spline = velocity_spline_failed
               return
            end if
            call spline_stencil_value(excess_height, scale_height, 0.0d0, &
                 scale_height_there, ierr)
            if (ierr /= 0) then
               failed_spline = scale_height_spline_failed
               return
            end if
            turnover_time = scale_height_there/velocity_at_scale_height
            timed = .true.
            exit
         end if
      end do
end subroutine time_scale_height_above_base

! ---------------------------------------------------------------
! The innermost point of the zone already lies more than its H_P above
! the base: tau = H_P/v_c there.
pure subroutine time_innermost_point(log10_pressure, log10_density, local_gravity, &
     convective_velocity, base_index, turnover_time)
      real(dp), intent(in) :: log10_pressure(:), log10_density(:), local_gravity(:), &
           convective_velocity(:)
      integer, intent(in) :: base_index
      real(dp), intent(out) :: turnover_time

      turnover_time = scale_height_at(log10_pressure, log10_density, local_gravity, &
           base_index)/convective_velocity(base_index)
end subroutine time_innermost_point

! ---------------------------------------------------------------
! The crossing time of the bottom convection cell, tau = width/v_c: the
! cell from the base to its upper edge (the zero of gradr - grada, linear
! in log r between its last convective point and the next), and v_c
! splined in r through four points to the cell's middle.
subroutine time_bottom_cell_crossing(convective_flag, log10_radius, &
     convective_velocity, radiative_gradient, adiabatic_gradient, num_points, &
     base_index, base_radius, turnover_time, failed_spline, ierr)
      logical, intent(in) :: convective_flag(:)
      real(dp), intent(in) :: log10_radius(:), convective_velocity(:), &
           radiative_gradient(:), adiabatic_gradient(:)
      integer, intent(in) :: num_points, base_index
      real(dp), intent(in) :: base_radius
      real(dp), intent(inout) :: turnover_time
      integer, intent(inout) :: failed_spline
      integer, intent(out) :: ierr
      real(dp) :: radius(stencil_size), velocity(stencil_size)
      real(dp) :: gap_at_top, gap_above, top_fraction, top_radius, cell_width, &
           log_middle_radius, middle_radius, velocity_at_middle
      integer :: k, top_index, first, n

      ierr = 0
!     The cell's last convective point.
      do k = base_index, num_points-1
         if (.not. convective_flag(k)) exit
      end do
      top_index = k - 1
      gap_at_top = radiative_gradient(top_index) - adiabatic_gradient(top_index)
      gap_above = radiative_gradient(top_index+1) - adiabatic_gradient(top_index+1)
      top_fraction = cz_base_fraction(gap_at_top, gap_above)
      top_radius = log10_radius(top_index) + top_fraction* &
           (log10_radius(top_index+1) - log10_radius(top_index))
      top_radius = exp(ln10*top_radius)
      cell_width = top_radius - base_radius
!     The first point above the cell's middle, r_base + width/2.
      log_middle_radius = log10(base_radius + 0.5d0*cell_width)
      do k = base_index, top_index+1
         if (log10_radius(k) > log_middle_radius) exit
      end do
      middle_radius = exp(ln10*log_middle_radius)
      first = stencil_start(k, base_index, num_points)
      do n = 1, stencil_size
         radius(n) = exp(ln10*log10_radius(first+n-1))
         velocity(n) = convective_velocity(first+n-1)
      end do
      call spline_stencil_value(radius, velocity, middle_radius, velocity_at_middle, ierr)
      if (ierr /= 0) then
         failed_spline = velocity_spline_failed
         return
      end if
      turnover_time = cell_width/velocity_at_middle
end subroutine time_bottom_cell_crossing

! ---------------------------------------------------------------
! The natural cubic spline through the stencil's stencil_size points
! (nodes, values) at x; ierr is splint's status.
pure subroutine spline_stencil_value(nodes, values, x, value_at_x, ierr)
      real(dp), intent(in) :: nodes(:), values(:), x
      real(dp), intent(out) :: value_at_x
      integer, intent(out) :: ierr
      real(dp) :: second_derivative(stencil_size)

      call natural_spline(nodes, values, stencil_size, second_derivative)
      call splint(nodes, values, stencil_size, second_derivative, x, value_at_x, ierr)
end subroutine spline_stencil_value

end module turnover_walk

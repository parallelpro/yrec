!----------------------------------------------------------------------
! mesh_spacing
!----------------------------------------------------------------------
! The arithmetic of the rezone's new mesh, on plain arrays: a star-blind
! kernel (docs/GUIDELINES.md, star rule 6). rezone (rezoning.f90) reads
! the model and the controls from the star, calls these in this order,
! and reports their failures; nothing here reads the star, prints or
! keeps state between calls.
!   accumulate_chi          the mesh coordinate chi of the old points;
!   sort_kept_points        the kept points in order, without duplicates;
!   space_points_along_chi  the new points, evenly spaced in chi between
!                           the kept points;
!   spline_onto_points      a profile at the new points;
!   insert_points_at_jumps  more points where a profile steps too far
!                           between two new points;
!   remove_crowded_points   new points too close in log m removed.
! A mesh is a run of points in log m (log10 of the mass inside the point,
! g), increasing outward. test_model_kernels checks each.
module mesh_spacing
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use numerics_lib, only: natural_spline, splint
      implicit none
      private
      public :: chi_scales
      public :: accumulate_chi, sort_kept_points, space_points_along_chi, &
           spline_onto_points, insert_points_at_jumps, remove_crowded_points

! The largest steps between two points of the new mesh, which scale chi:
!   log_mass               in log m [dex];
!   luminosity             in L, where L rises outward [the unit of the
!                          luminosity passed];
!   log_pressure_core      in log P below the fine zone [dex];
!   log_pressure_base      in log P inside the fine zone around the base
!                          of the surface convection zone [dex];
!   log_pressure_envelope  in log P above the fine zone [dex];
!   base_fine_width        the half width of the fine zone in log P [dex].
      type :: chi_scales
         real(dp) :: log_mass = 0.0d0
         real(dp) :: luminosity = 0.0d0
         real(dp) :: log_pressure_core = 0.0d0
         real(dp) :: log_pressure_base = 0.0d0
         real(dp) :: log_pressure_envelope = 0.0d0
         real(dp) :: base_fine_width = 0.0d0
      end type chi_scales

! A spacing limit at or below spacing_limit_off is off [the profile's unit].
      real(dp), parameter :: spacing_limit_off = 1.0D-15

contains

! ---------------------------------------------------------------
! chi, the running sum of the normalized changes from point to point
! (MHP 7/02):
!   chi(1) = 1,
!   chi(j) = chi(j-1) + Delta log m/scales%log_mass
!            + Delta L/scales%luminosity - Delta log P/dp_scale,
! the L term only where L increases outward. dp_scale depends on the
! point's log P above that of the base of the surface convection zone,
! base_log_pressure: log_pressure_base within base_fine_width of it,
! log_pressure_core deeper, and log_pressure_envelope above. A unit step
! in chi is the largest step allowed between two new points. The arrays
! hold the old points.
pure subroutine accumulate_chi(log_mass, luminosity, log_pressure, base_log_pressure, &
     scales, chi)
      real(dp), intent(in) :: log_mass(:), luminosity(:), log_pressure(:), &
           base_log_pressure
      type(chi_scales), intent(in) :: scales
      real(dp), intent(inout) :: chi(:)
      real(dp) :: pressure_test, dp_scale, dchi
      integer :: j

      chi(1) = 1.0D0
      do j = 2, size(log_mass)
         pressure_test = log_pressure(j) - base_log_pressure
         if (abs(pressure_test) < scales%base_fine_width) then
            dp_scale = scales%log_pressure_base
         else if (pressure_test > scales%base_fine_width) then
            dp_scale = scales%log_pressure_core
         else
            dp_scale = scales%log_pressure_envelope
         end if
         if (luminosity(j) > luminosity(j-1)) then
            dchi = (log_mass(j)-log_mass(j-1))/scales%log_mass + &
                 (luminosity(j)-luminosity(j-1))/scales%luminosity - &
                 (log_pressure(j)-log_pressure(j-1))/dp_scale
         else
            dchi = (log_mass(j)-log_mass(j-1))/scales%log_mass - &
                 (log_pressure(j)-log_pressure(j-1))/dp_scale
         end if
         chi(j) = chi(j-1)+dchi
      end do
end subroutine accumulate_chi

! ---------------------------------------------------------------
! The kept points flag_point(1:flag_count) sorted into ascending order
! (a bubble sort), then the duplicates removed: flag_count comes back as
! the number left. After a removal the scan moves on to the next slot,
! so of three equal points two are left.
pure subroutine sort_kept_points(flag_point, flag_count)
      integer, intent(inout) :: flag_point(:)
      integer, intent(inout) :: flag_count
      integer :: i, j, temp_swap
      logical :: sort_done

      if (flag_count == 1) return
      do
         sort_done = .true.
         do i = 1,flag_count-1
            if (flag_point(i+1) < flag_point(i)) then
               temp_swap = flag_point(i)
               flag_point(i) = flag_point(i+1)
               flag_point(i+1) = temp_swap
               sort_done = .false.
            end if
         end do
         if (sort_done) exit
      end do
      i = 2
      do
         if (flag_point(i) == flag_point(i-1)) then
            if (i < flag_count) then
               do j = i,flag_count-1
                  flag_point(j) = flag_point(j+1)
               end do
            end if
            flag_count = flag_count - 1
         end if
         i = i + 1
         if (i > flag_count) exit
      end do
end subroutine sort_kept_points

! ---------------------------------------------------------------
! The new mesh new_log_mass(1:new_num_zones): the first old point, then
! between successive kept points (flag_point, old point numbers in
! ascending order, the last one the outermost old point) int(Delta chi)
! + 1 intervals, evenly spaced in chi, so that every step in chi is below
! one; log m is read off a natural cubic spline of log m against chi, so
! the kept points themselves are points of the new mesh. needed_points is
! 0, or the number of points the mesh needs when that is more than
! new_log_mass holds; ierr is the spline's status (numerics_status). The
! mesh is not complete after either.
pure subroutine space_points_along_chi(chi, log_mass, flag_point, new_log_mass, &
     new_num_zones, needed_points, ierr)
      real(dp), intent(in) :: chi(:), log_mass(:)
      integer, intent(in) :: flag_point(:)
      real(dp), intent(inout) :: new_log_mass(:)
      integer, intent(out) :: new_num_zones, needed_points, ierr
      real(dp) :: spline_second_deriv(size(chi))
      real(dp) :: dchi, chi_prev, delta_chi, spline_eval_x, spline_eval_y
      integer :: i, j, point_prev_index, segment_point_count, num_old

      ierr = 0
      needed_points = 0
      num_old = size(chi)
      call natural_spline(chi, log_mass, num_old, spline_second_deriv)
      new_log_mass(1) = log_mass(1)
      chi_prev = chi(1)
      point_prev_index = 1
      new_num_zones = 1
      do i = 1, size(flag_point)
         delta_chi = chi(flag_point(i))-chi_prev
         segment_point_count = int(delta_chi)+1
         if (new_num_zones+segment_point_count > size(new_log_mass)) then
            needed_points = new_num_zones+segment_point_count
            return
         end if
         dchi = delta_chi/dble(segment_point_count)
         do j = point_prev_index+1,point_prev_index+segment_point_count
            spline_eval_x = chi_prev + dchi
            call splint(chi, log_mass, num_old, spline_second_deriv, spline_eval_x, &
                 spline_eval_y, ierr)
            if (ierr /= 0) return
            new_log_mass(j) = spline_eval_y
            chi_prev = spline_eval_x
         end do
         new_num_zones = new_num_zones + segment_point_count
         point_prev_index = new_num_zones
      end do
end subroutine space_points_along_chi

! ---------------------------------------------------------------
! A profile of the old mesh (log_mass, profile) at the points
! new_log_mass, by a natural cubic spline against log m; the first new
! point is the first old one and takes its value. ierr is the spline's
! status.
pure subroutine spline_onto_points(log_mass, profile, new_log_mass, values, ierr)
      real(dp), intent(in) :: log_mass(:), profile(:), new_log_mass(:)
      real(dp), intent(inout) :: values(:)
      integer, intent(out) :: ierr
      real(dp) :: spline_second_deriv(size(log_mass))
      integer :: i, num_old

      ierr = 0
      num_old = size(log_mass)
      values(1) = profile(1)
      call natural_spline(log_mass, profile, num_old, spline_second_deriv)
      do i = 2, size(new_log_mass)
         call splint(log_mass, profile, num_old, spline_second_deriv, new_log_mass(i), &
              values(i), ierr)
         if (ierr /= 0) return
      end do
end subroutine spline_onto_points

! ---------------------------------------------------------------
! Where a profile steps outward by more than max_step between two points
! of the new mesh, int(step/max_step) points are inserted between them,
! evenly in log m. values holds the profile at the new points
! new_log_mass(1:new_num_zones) on entry; the step from point j-1 to j is
! its rise, or with falling its fall (a profile that rises inward steps
! by nothing). Nothing is done when max_step is off (at most
! spacing_limit_off). The intervals are taken outermost first, so that an
! insertion moves only points outside the intervals still to be filled.
! needed_points is 0, or the number of points the mesh needs when that is
! more than new_log_mass holds (the mesh is then not complete).
pure subroutine insert_points_at_jumps(values, falling, max_step, new_log_mass, &
     new_num_zones, needed_points)
      real(dp), intent(in) :: values(:), max_step
      logical, intent(in) :: falling
      real(dp), intent(inout) :: new_log_mass(:)
      integer, intent(inout) :: new_num_zones
      integer, intent(out) :: needed_points
      integer :: jump_index(new_num_zones)
      integer :: i, j, k, jump_count, working_num_zones, num_new_points
      real(dp) :: step_over_max, point_insert_spacing

      needed_points = 0
      if (.not. (max_step > spacing_limit_off)) return
      jump_count = 0
      do j = new_num_zones,2,-1
         if (outward_step(values, j, falling) > max_step) then
            jump_count = jump_count + 1
            jump_index(jump_count) = j
         end if
      end do
      if (jump_count == 0) return
      working_num_zones = new_num_zones
      do i = 1,jump_count
         j = jump_index(i)
         step_over_max = outward_step(values, j, falling)/max_step
!        Shift the points from j outward up by num_new_points, then fill
!        the interval below j evenly in log m.
         num_new_points = int(step_over_max)
         point_insert_spacing = (new_log_mass(j)-new_log_mass(j-1))/ &
              dble(num_new_points+1)
         if (working_num_zones+num_new_points > size(new_log_mass)) then
            needed_points = working_num_zones+num_new_points
            return
         end if
         do k = working_num_zones+num_new_points,j+num_new_points,-1
            new_log_mass(k) = new_log_mass(k-num_new_points)
         end do
         do k = j + num_new_points -1, j -1, -1
            new_log_mass(k) = new_log_mass(k+1) - point_insert_spacing
         end do
         working_num_zones = working_num_zones + num_new_points
      end do
      new_num_zones = working_num_zones
end subroutine insert_points_at_jumps

! ---------------------------------------------------------------
! The step of a profile from point j-1 outward to point j: its rise, or
! with falling its fall.
pure real(dp) function outward_step(values, j, falling)
      real(dp), intent(in) :: values(:)
      integer, intent(in) :: j
      logical, intent(in) :: falling

      if (falling) then
         outward_step = values(j-1) - values(j)
      else
         outward_step = values(j) - values(j-1)
      end if
end function outward_step

! ---------------------------------------------------------------
! The new points that lie within min_spacing in log m of the last point
! kept are dropped; the first and the last point are always kept.
pure subroutine remove_crowded_points(new_log_mass, new_num_zones, min_spacing)
      real(dp), intent(inout) :: new_log_mass(:)
      integer, intent(inout) :: new_num_zones
      real(dp), intent(in) :: min_spacing
      integer :: j, k

      j = 1
      do k = 2,new_num_zones-1
         if (new_log_mass(k) - new_log_mass(j) > min_spacing) then
            j = j + 1
            new_log_mass(j) = new_log_mass(k)
         end if
      end do
      j = j + 1
      new_log_mass(j) = new_log_mass(new_num_zones)
      new_num_zones = j
end subroutine remove_crowded_points

end module mesh_spacing

!----------------------------------------------------------------------
! shell_boundaries
!----------------------------------------------------------------------
! Where the model's regions are: the outer edge of a convective core,
! the inner edge of the surface convection zone, and the beginning,
! middle and end of a hydrogen-burning shell.
!   find_shell_boundaries     the kernel: the profiles of X, L and the
!                             convective flag in, the zone indices out
!                             (no star);
!   locate_shell_boundaries   the driver: the kernel on the star's model,
!                             its results stored on the star.
! model_rescaling calls the kernel on the model it is rescaling.
module shell_boundaries
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use star_info_lib, only: star_info
      implicit none
      private
      public :: locate_shell_boundaries, find_shell_boundaries

contains

! ---------------------------------------------------------------
! The star's convective core, surface convection zone and hydrogen
! shell (find_shell_boundaries on the model's X, L/Lsun and convective
! flags; the shell's threshold is timestep_factors(itime_core_min)).
!   star reads: convective_flag, luminosity_lsun, nz,
!     physics%timestep_factors, xa
!   star writes: core_cz_top_index, envelope_cz_bottom_index,
!     h_shell_bot_index, h_shell_mid_index, h_shell_top_index,
!     has_h_shell
subroutine locate_shell_boundaries(s)
      use star_info_lib, only: i_h1
      use control_slots, only: itime_core_min
      type(star_info), intent(inout) :: s

      call find_shell_boundaries(s%xa(i_h1,:), s%luminosity_lsun, s%convective_flag, s%nz, &
           s%physics%timestep_factors(itime_core_min), s%core_cz_top_index, &
           s%envelope_cz_bottom_index, s%h_shell_bot_index, s%h_shell_top_index, &
           s%h_shell_mid_index, s%has_h_shell)
end subroutine locate_shell_boundaries

! ---------------------------------------------------------------
! The region boundaries of a model of num_points points (1 at the
! centre), from its hydrogen mass fraction h1, its luminosity (any
! unit) and its convective flags:
!   core_edge      the outermost point of the convective core (1 when
!                  there is none);
!   envelope_edge  the innermost point of the surface convection zone
!                  (num_points when there is none);
!   has_h_shell    whether the model burns hydrogen in a shell, that is,
!                  whether the central X is at most x_core_min;
! and, when it does (1 otherwise):
!   shell_begin    one past the last point, counted from the centre,
!                  with X <= x_core_min;
!   shell_mid      one past the last point with X <= x_core_min or X at
!                  most half the surface X;
!   shell_end      the first point beyond those where
!                  L(i) - L(i-1) < 1e-5 L(num_points), or where X is
!                  within 1e-5 of the surface X (num_points if none).
! A fully convective model (envelope_edge = 1) has core_edge = 1: the
! nuclear timestep reads core_edge, and a convective core reaching the
! surface would distort it.
pure subroutine find_shell_boundaries(h1, luminosity, is_convective, num_points, &
     x_core_min, core_edge, envelope_edge, shell_begin, shell_end, shell_mid, &
     has_h_shell)
      real(dp), intent(in) :: h1(:), luminosity(:)
      logical, intent(in) :: is_convective(:)
      integer, intent(in) :: num_points
      real(dp), intent(in) :: x_core_min
      integer, intent(out) :: core_edge, envelope_edge, shell_begin, &
           shell_end, shell_mid
      logical, intent(out) :: has_h_shell

! The end of the shell: the luminosity's step falls below this fraction
! of the surface luminosity, or X comes this close to the surface X.
      real(dp), parameter :: luminosity_change_tol = 1.0e-5_dp, &
           hydrogen_surface_tol = 1.0e-5_dp
      real(dp) :: half_surface_x, luminosity_end_threshold
      integer :: i

! The hydrogen shell.
      shell_begin = 1
      shell_mid = 1
      shell_end = 1
      has_h_shell = .false.
      if (h1(1) <= x_core_min) then
         has_h_shell = .true.
         half_surface_x = 0.5_dp*h1(num_points)
         luminosity_end_threshold = luminosity_change_tol*luminosity(num_points)
         do i = 1, num_points
            if (h1(i) <= x_core_min) then
               shell_begin = shell_begin + 1
               shell_mid = shell_mid + 1
            else if (h1(i) <= half_surface_x) then
               shell_mid = shell_mid + 1
            else if (luminosity(i) - luminosity(i-1) < luminosity_end_threshold) then
               exit
            else if (h1(num_points) - h1(i) < hydrogen_surface_tol) then
               exit
            end if
         end do
         if (i > num_points) then
            i = num_points
         end if
         shell_end = i
      end if
! The edge of the convective core.
      do i = 1, num_points
         if (.not. is_convective(i)) exit
      end do
      if (i > 1) then
         core_edge = i - 1
      else
         core_edge = 1
      end if
! The base of the surface convection zone.
      do i = num_points, 1, -1
         if (.not. is_convective(i)) exit
      end do
      if (i < num_points) then
         envelope_edge = i + 1
      else
         envelope_edge = num_points
      end if
      if (envelope_edge == 1) core_edge = 1
end subroutine find_shell_boundaries

! ---------------------------------------------------------------
! Sources
!   F77 FINDSH (git show 6cd56739:src/findsh.f). The shell's beginning
!   is tested against timestep_factors(itime_core_min), the threshold of
!   the core-burning criterion, rather than X <= 1e-10, for consistency
!   with that criterion (LLP 9/24/08).
end module shell_boundaries

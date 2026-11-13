!----------------------------------------------------------------------
! test_model_kernels
!----------------------------------------------------------------------
! ok/FAIL checks of model/'s star-free kernels, each confirmed to fail
! under a mutation (test_model_kernels.mutations.toml):
!   find_shell_boundaries (shell_boundaries.f90) on synthetic profiles of
!   40 points: the edge of a convective core, the base of a surface
!   convection zone, a fully convective model, and the hydrogen shell's
!   beginning, middle and end (by the luminosity's step and by the
!   surface X), with the central X on either side of x_core_min.
! (regridding's kernels have test_regrid.)
! Prints one line per check and exits with status 1 when one fails.
program test_model_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use shell_boundaries, only: find_shell_boundaries
      implicit none
      integer, parameter :: n = 40
      real(dp), parameter :: x_core_min = 1.0e-5_dp, surface_x = 0.7_dp
      real(dp) :: h1(n), luminosity(n)
      logical :: convective(n)
      integer :: core_edge, envelope_edge, shell_begin, shell_end, shell_mid, i
      logical :: has_h_shell
      integer :: num_failed = 0

! The zones: radiative throughout; a convective core of 10 points; a
! surface convection zone from point 31; fully convective.
      h1 = surface_x
      luminosity = [(real(i, dp), i = 1, n)]
      convective = .false.
      call boundaries()
      call model_kernel_check('find_shell_boundaries: radiative model: core_edge 1, envelope_edge n', &
           core_edge == 1 .and. envelope_edge == n)
      convective(1:10) = .true.
      call boundaries()
      call model_kernel_check('find_shell_boundaries: convective core of 10 points: core_edge 10', &
           core_edge == 10 .and. envelope_edge == n)
      convective = .false.
      convective(31:n) = .true.
      call boundaries()
      call model_kernel_check('find_shell_boundaries: surface convection zone from 31: envelope_edge 31', &
           envelope_edge == 31 .and. core_edge == 1)
      convective = .true.
      call boundaries()
      call model_kernel_check('find_shell_boundaries: fully convective: envelope_edge 1, core_edge 1', &
           envelope_edge == 1 .and. core_edge == 1)

! The hydrogen shell. Central X above x_core_min: none.
      convective = .false.
      h1 = surface_x
      h1(1) = 10.0_dp*x_core_min
      call boundaries()
      call model_kernel_check('find_shell_boundaries: central X above x_core_min: no H shell', &
           .not. has_h_shell .and. shell_begin == 1 .and. shell_mid == 1 .and. shell_end == 1)
! Points 1-10 depleted; 11-16 below half the surface X (16 just below);
! 17 on, above it; L rises by 1 per point to point 20, then is flat.
      h1(1:10) = 1.0e-12_dp
      h1(11:15) = [(0.05_dp*(i - 10), i = 11, 15)]
      h1(16) = 0.34_dp
      h1(17:n) = [(0.36_dp + 0.3_dp*(i - 17)/(n - 17), i = 17, n)]
      h1(n) = surface_x
      luminosity = [(real(min(i, 20), dp), i = 1, n)]
      call boundaries()
      call model_kernel_check('find_shell_boundaries: H shell begins past the last X <= x_core_min (11)', &
           has_h_shell .and. shell_begin == 11)
      call model_kernel_check('find_shell_boundaries: H shell middle past the last X <= X_surf/2 (17)', &
           shell_mid == 17)
      call model_kernel_check('find_shell_boundaries: H shell ends where the L step < 1e-5 L(n) (21)', &
           shell_end == 21)
! L rising to the surface; X within 1e-5 of the surface X from point 30
! (points 28 and 29 within 1e-3).
      luminosity = [(real(i, dp), i = 1, n)]
      h1(28) = surface_x - 5.0e-4_dp
      h1(29) = surface_x - 2.0e-4_dp
      h1(30:n-1) = surface_x - 5.0e-6_dp
      call boundaries()
      call model_kernel_check('find_shell_boundaries: H shell ends where X is within 1e-5 of X_surf (30)', &
           shell_end == 30)

      if (num_failed > 0) then
         write(*,'(a,i0,a)') 'test_model_kernels: ', num_failed, ' checks FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_model_kernels: done'

contains

subroutine boundaries()
      call find_shell_boundaries(h1, luminosity, convective, n, x_core_min, core_edge, &
           envelope_edge, shell_begin, shell_end, shell_mid, has_h_shell)
end subroutine boundaries

! One check: its label and ok, or FAIL.
subroutine model_kernel_check(label, condition)
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
end subroutine model_kernel_check

end program test_model_kernels

!----------------------------------------------------------------------
! test_model_kernels
!----------------------------------------------------------------------
! ok/FAIL checks of model/'s star-free kernels, each confirmed to fail
! under a mutation (test_model_kernels.mutations.toml):
!   find_shell_boundaries (shell_boundaries.f90) on synthetic profiles of
!   40 points: the edge of a convective core, the base of a surface
!   convection zone, a fully convective model, and the hydrogen shell's
!   beginning, middle and end (by the luminosity's step and by the
!   surface X), with the central X on either side of x_core_min;
!   compute_envelope_mixture (model_rescaling.f90) on a mixture of ten
!   metals with unequal weights: Y = 1 - X - Z - X(He3); the number
!   fractions sum to 1; the mass fractions they give back (n_i A_i/
!   sum n_j A_j) are X for H, 1 - X - Z for He and Z w_i/S for each
!   metal (relative 1e-13); 1/mu_ion is sum X_i/A_i of those mass
!   fractions; Z without C, N and O is Z (S - w_C - w_N - w_O)/S.
! (regridding's kernels have test_regrid.)
! Prints one line per check and exits with status 1 when one fails.
program test_model_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use shell_boundaries, only: find_shell_boundaries
      use model_rescaling, only: compute_envelope_mixture
      use composition_layout, only: n_mix_species, ix_h, ix_he, ix_c, ix_n, ix_o
      implicit none
      integer, parameter :: n = 40
      real(dp), parameter :: x_core_min = 1.0e-5_dp, surface_x = 0.7_dp
      real(dp) :: h1(n), luminosity(n)
      logical :: convective(n)
      integer :: core_edge, envelope_edge, shell_begin, shell_end, shell_mid, i
      logical :: has_h_shell
      integer :: num_failed = 0
! compute_envelope_mixture's case
      real(dp), parameter :: x_env = 0.7_dp, z_env = 0.02_dp, he3_env = 3.0e-5_dp, &
           mixture_tol = 1.0e-13_dp
      real(dp) :: seed(n_mix_species), number_fractions(n_mix_species), &
           mass_back(n_mix_species), y_env, metals_no_cno, inverse_atomic_weight, &
           metal_sum
      integer :: k

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

! compute_envelope_mixture: ten metals of weights 1..12 (H's and He's
! slots hold values the kernel must not read).
      seed = [(real(k, dp), k = 1, n_mix_species)]
      seed(ix_h) = 1.0e3_dp
      seed(ix_he) = -1.0e3_dp
      metal_sum = sum(seed) - seed(ix_h) - seed(ix_he)
      call compute_envelope_mixture(x_env, z_env, he3_env, seed, y_env, metals_no_cno, &
           inverse_atomic_weight, number_fractions)
      call model_kernel_check('compute_envelope_mixture: Y = 1 - X - Z - X(He3)', &
           y_env == 1.0_dp - x_env - z_env - he3_env)
      call model_kernel_check('compute_envelope_mixture: the number fractions sum to 1', &
           abs(sum(number_fractions) - 1.0_dp) < mixture_tol)
! The mass fractions back from the number fractions and the atomic
! weights, which sum_j n_j A_j = 1/(1/mu_ion) normalizes.
      mass_back = number_fractions*atomic_weights()
      mass_back = mass_back/sum(mass_back)
      call model_kernel_check('compute_envelope_mixture: X(H) and X(He) = 1 - X - Z back', &
           abs(mass_back(ix_h)/x_env - 1.0_dp) < mixture_tol .and. &
           abs(mass_back(ix_he)/(1.0_dp - x_env - z_env) - 1.0_dp) < mixture_tol)
      call model_kernel_check('compute_envelope_mixture: each metal Z w_i/S back', &
           all([(abs(mass_back(k)/(z_env*seed(k)/metal_sum) - 1.0_dp) < mixture_tol, &
           k = 1, n_mix_species)] .or. [(k == ix_h .or. k == ix_he, k = 1, n_mix_species)]))
! 1/mu_ion from the mass fractions the kernel should have used.
      mass_back = z_env*seed/metal_sum
      mass_back(ix_h) = x_env
      mass_back(ix_he) = 1.0_dp - x_env - z_env
      call model_kernel_check('compute_envelope_mixture: 1/mu_ion = sum X_i/A_i', &
           abs(inverse_atomic_weight/sum(mass_back/atomic_weights()) - 1.0_dp) < mixture_tol)
      call model_kernel_check('compute_envelope_mixture: Z without C, N, O = Z (S-wC-wN-wO)/S', &
           abs(metals_no_cno/(z_env*(metal_sum - seed(ix_c) - seed(ix_n) - seed(ix_o))/metal_sum) &
           - 1.0_dp) < mixture_tol)

      if (num_failed > 0) then
         write(*,'(a,i0,a)') 'test_model_kernels: ', num_failed, ' checks FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_model_kernels: done'

contains

! The atomic weights compute_envelope_mixture uses, in ix_* slot order
! (Na, Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He).
pure function atomic_weights() result(a)
      real(dp) :: a(n_mix_species)
      a = [23.0_dp, 26.99_dp, 24.32_dp, 55.86_dp, 28.1_dp, 12.015_dp, 1.008_dp, &
           16.0_dp, 14.01_dp, 39.96_dp, 20.19_dp, 4.004_dp]
end function atomic_weights

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

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
!   compute_envelope_mixture (rescaling_kernels.f90) on a mixture of ten
!   metals with unequal weights: Y = 1 - X - Z - X(He3); the number
!   fractions sum to 1; the mass fractions they give back (n_i A_i/
!   sum n_j A_j) are X for H, 1 - X - Z for He and Z w_i/S for each
!   metal (relative 1e-13); 1/mu_ion is sum X_i/A_i of those mass
!   fractions; Z without C, N and O is Z (S - w_C - w_N - w_O)/S;
!   step/convection_geometry's mixed_base_mass and mixed_base_fraction (the
!   Li/Be burn's mixed-region base, O9) on 40 zones with m(k) ~ k^2 and
!   rezoning's dm: the region above the base holds sum dm(k..n) + f (m(k) -
!   m(k-1)); the fraction comes back from the base's mass (1e-12); a base a
!   zone away clips to 0.5 or -0.5; no base, or the first or last zone,
!   gives 0.
!   rescaling_kernels' arithmetic of rescale_model on a 40-point model
!   of 0.8 Msun with m(k) = M (k/41)^1.5, whose hydrogen shell is points
!   21 to 25: hb_region_masses, whose core, shell and envelope add up to
!   M; a total mass of 0.9 Msun by stretching the mass above the shell
!   (mass_stretch_factor, stretch_masses_above: every point above the
!   base moves by the factor, the points below not at all); a core mass
!   raised by 0.02 Msun at fixed shell mass and fixed outer point
!   (move_core_mass); scale_hydrogen, shift_metals, ramp_metals and
!   set_helium_by_difference against their formulas.
!   mesh_spacing's arithmetic of the rezone's new mesh: accumulate_chi on
!   40 points with uniform steps (the mass, luminosity and pressure terms,
!   the three pressure regions); sort_kept_points (order, duplicates, and
!   the two it leaves of three equal points); space_points_along_chi on a
!   mesh whose log m is linear in chi (int(Delta chi) + 1 even intervals
!   per segment, the kept points on the mesh, the count needed when the
!   arrays are too small); spline_onto_points on a profile linear in
!   log m; insert_points_at_jumps (int(step/max) points, evenly in log m,
!   for a rise or with falling for a fall, none when the limit is off, the
!   count needed when they do not fit); remove_crowded_points.
! (regridding's kernels have test_regrid.)
! Prints one line per check and exits with status 1 when one fails.
program test_model_kernels
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use shell_boundaries, only: find_shell_boundaries
      use rescaling_kernels, only: compute_envelope_mixture
      use convection_geometry, only: mixed_base_mass, mixed_base_fraction
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
! the mixed-region base's case: zone-centre masses m(k) ~ k^2 (g) and the
! shell masses as rezoning sets them
      real(dp) :: zone_mass(n), shell_mass(n), base_mass, region_mass

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

! mixed_base_mass and mixed_base_fraction: 40 zones, m(k) = 1e30 k^2 g.
      zone_mass = [(1.0e30_dp*real(i*i, dp), i = 1, n)]
      shell_mass(1) = 0.5_dp*(zone_mass(1) + zone_mass(2))
      shell_mass(2:n-1) = 0.5_dp*(zone_mass(3:n) - zone_mass(1:n-2))
      shell_mass(n) = 1.0e30_dp*real(n*n + n, dp) - 0.5_dp*(zone_mass(n) + zone_mass(n-1))
      base_mass = mixed_base_mass(zone_mass, shell_mass, 25, 0.3_dp, n)
      region_mass = sum(shell_mass) - base_mass
      call model_kernel_check('mixed_base_mass: the region holds sum dm(k..n) + f (m(k) - m(k-1))', &
           abs(region_mass/(sum(shell_mass(25:n)) + 0.3_dp*(zone_mass(25) - zone_mass(24))) &
           - 1.0_dp) < 1.0e-13_dp)
      call model_kernel_check('mixed_base_fraction: f back from the base mass (0.3, -0.2)', &
           abs(mixed_base_fraction(zone_mass, shell_mass, 25, base_mass, n) - 0.3_dp) < 1.0e-12_dp &
           .and. abs(mixed_base_fraction(zone_mass, shell_mass, 25, &
           mixed_base_mass(zone_mass, shell_mass, 25, -0.2_dp, n), n) + 0.2_dp) < 1.0e-12_dp)
      call model_kernel_check('mixed_base_fraction: a base a zone away clips to 0.5, -0.5', &
           mixed_base_fraction(zone_mass, shell_mass, 25, &
           mixed_base_mass(zone_mass, shell_mass, 24, 0.0_dp, n), n) == 0.5_dp .and. &
           mixed_base_fraction(zone_mass, shell_mass, 25, &
           mixed_base_mass(zone_mass, shell_mass, 26, 0.0_dp, n), n) == -0.5_dp)
      call model_kernel_check('mixed_base_fraction: 0 with no base, at zone 1 or zone n', &
           mixed_base_fraction(zone_mass, shell_mass, 25, 0.0_dp, n) == 0.0_dp .and. &
           mixed_base_fraction(zone_mass, shell_mass, 1, base_mass, n) == 0.0_dp .and. &
           mixed_base_fraction(zone_mass, shell_mass, n, base_mass, n) == 0.0_dp)

      call rescaling_checks()
      call mesh_checks()

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

! ---------------------------------------------------------------
! rescaling_kernels on a 40-point model of 0.8 Msun.
subroutine rescaling_checks()
      use rescaling_kernels, only: hb_masses, scale_hydrogen, shift_metals, ramp_metals, &
           set_helium_by_difference, hb_region_masses, mass_stretch_factor, &
           stretch_masses_above, move_core_mass
      use composition_layout, only: i_h1, i_he4, i_metals, i_he3, i_c12, i_o18, i_h2, &
           n_species_extended
      real(dp), parameter :: solar_mass = 2.0e33_dp, total_msun = 0.8_dp, tol = 1.0e-10_dp
      integer, parameter :: shell_begin = 21, shell_end = 25
      real(dp) :: log_mass(n), old_mass(n), new_mass(n), log_total_mass, factor, &
           xa(n_species_extended,n), old_xa(n_species_extended,n), q(n), ramp_z
      type(hb_masses) :: masses
      logical :: ok
      integer :: j

      old_mass = [(total_msun*solar_mass*(real(j, dp)/real(n + 1, dp))**1.5_dp, j = 1, n)]
      log_mass = log10(old_mass)
      log_total_mass = log10(total_msun*solar_mass)

! The regions' masses.
      masses = hb_region_masses(log_mass, log_total_mass, shell_begin, shell_end, solar_mass)
      call model_kernel_check('hb_region_masses: core, shell and envelope add up to M', &
           abs((masses%core + masses%shell + masses%envelope)/total_msun - 1.0_dp) < tol .and. &
           abs(masses%core/(old_mass(shell_begin-1)/solar_mass) - 1.0_dp) < tol .and. &
           abs(masses%inner_envelope/((old_mass(n) - old_mass(shell_end))/solar_mass) - 1.0_dp) &
           < tol)

! A total mass of 0.9 Msun, the mass above the shell's end stretched.
      factor = mass_stretch_factor(log_mass(shell_end), 0.9_dp, total_msun, solar_mass)
      call stretch_masses_above(log_mass, log_total_mass, shell_end, factor)
      new_mass = 10.0_dp**log_mass
      call model_kernel_check('stretch_masses_above: the new total mass, the base and below unmoved', &
           abs(10.0_dp**log_total_mass/(0.9_dp*solar_mass) - 1.0_dp) < tol .and. &
           all(new_mass(1:shell_end) == 10.0_dp**log10(old_mass(1:shell_end))))
      call model_kernel_check('stretch_masses_above: every point above moves by (M_new - m_b)/(M_old - m_b)', &
           all(abs((new_mass(shell_end+1:n) - new_mass(shell_end))/ &
           (old_mass(shell_end+1:n) - old_mass(shell_end))/ &
           ((0.9_dp*solar_mass - old_mass(shell_end))/(total_msun*solar_mass - old_mass(shell_end))) &
           - 1.0_dp) < tol))

! The core mass raised by 0.02 Msun.
      log_mass = log10(old_mass)
      masses = hb_region_masses(log_mass, log10(total_msun*solar_mass), shell_begin, shell_end, &
           solar_mass)
      call move_core_mass(log_mass, shell_begin, shell_end, masses%core + 0.02_dp, masses%core, &
           solar_mass)
      new_mass = 10.0_dp**log_mass
      call model_kernel_check('move_core_mass: the new core mass, its points scaled alike', &
           abs(new_mass(shell_begin-1)/((masses%core + 0.02_dp)*solar_mass) - 1.0_dp) < tol .and. &
           all(abs(new_mass(1:shell_begin-1)/old_mass(1:shell_begin-1)/ &
           ((masses%core + 0.02_dp)/masses%core) - 1.0_dp) < tol))
      call model_kernel_check('move_core_mass: the shell keeps its mass', &
           all(abs((new_mass(shell_begin:shell_end) - new_mass(shell_begin-1))/ &
           (old_mass(shell_begin:shell_end) - old_mass(shell_begin-1)) - 1.0_dp) < tol))
      ok = abs(new_mass(n)/old_mass(n) - 1.0_dp) < tol
      do j = shell_end + 1, n - 1
         ok = ok .and. abs((new_mass(n) - new_mass(j))/(old_mass(n) - old_mass(j))/ &
              ((new_mass(n) - new_mass(shell_end))/(old_mass(n) - old_mass(shell_end))) - 1.0_dp) &
              < tol
      end do
      call model_kernel_check('move_core_mass: the outer point stays, the envelope shrinks toward it', ok)

! The composition: X rising from 0.5 at the centre to 0.7, Z from 0.03 to
! 0.02, with CNO isotopes and deuterium.
      do j = 1, n
         xa(:,j) = 0.0_dp
         xa(i_h1,j) = 0.5_dp + 0.2_dp*real(j - 1, dp)/real(n - 1, dp)
         xa(i_metals,j) = 0.03_dp - 0.01_dp*real(j - 1, dp)/real(n - 1, dp)
         xa(i_he3,j) = 1.0e-4_dp
         xa(i_c12:i_o18,j) = 1.0e-3_dp*real(j, dp)
         xa(i_h2,j) = 2.0e-5_dp
         xa(i_he4,j) = 0.123_dp
      end do
      old_xa = xa
      call scale_hydrogen(xa(i_h1,:), 0.6_dp, 0.7_dp)
      call model_kernel_check('scale_hydrogen: X_i X_new/X_surf, the surface at X_new', &
           abs(xa(i_h1,n)/0.6_dp - 1.0_dp) < tol .and. &
           all(abs(xa(i_h1,:)/(old_xa(i_h1,:)*0.6_dp/0.7_dp) - 1.0_dp) < tol))
      xa = old_xa
      call scale_hydrogen(xa(i_h1,:), 0.6_dp, 0.55_dp)
      call model_kernel_check('scale_hydrogen: no shell above X_new', &
           all(xa(i_h1,:) <= 0.6_dp) .and. xa(i_h1,n) == 0.6_dp .and. &
           abs(xa(i_h1,1)/(0.5_dp*0.6_dp/0.55_dp) - 1.0_dp) < tol)
      xa = old_xa
      call shift_metals(xa, -0.025_dp)
      call model_kernel_check('shift_metals: Z + dZ, not below 0; CNO by Z_new/Z_old; the rest kept', &
           all(abs(xa(i_metals,1:n/2) - (old_xa(i_metals,1:n/2) - 0.025_dp)) < 1.0e-15_dp) .and. &
           xa(i_metals,n) == 0.0_dp .and. xa(i_c12,n) == 0.0_dp .and. &
           all(abs(xa(i_o18,1:n/2)/old_xa(i_o18,1:n/2)/ &
           (xa(i_metals,1:n/2)/old_xa(i_metals,1:n/2)) - 1.0_dp) < tol) .and. &
           all(xa(i_h2,:) == old_xa(i_h2,:)) .and. all(xa(i_he3,:) == old_xa(i_he3,:)) .and. &
           all(xa(i_h1,:) == old_xa(i_h1,:)))
      xa = old_xa
      call set_helium_by_difference(xa)
      call model_kernel_check('set_helium_by_difference: X + Y + Z + X(He3) = 1', &
           all(abs(xa(i_h1,:) + xa(i_he4,:) + xa(i_metals,:) + xa(i_he3,:) - 1.0_dp) < 1.0e-15_dp))
! A ramp from Z = 0.04 inside q = 0.2 to the surface's 0.02 at q = 0.6.
      log_mass = log10(old_mass)
      q = old_mass/(total_msun*solar_mass)
      call ramp_metals(xa, log_mass, total_msun*solar_mass, 0.04_dp, 0.2_dp, 0.6_dp, 0.02_dp)
      ok = .true.
      do j = 1, n
         if (q(j) < 0.2_dp) then
            ramp_z = 0.04_dp
         else if (q(j) < 0.6_dp) then
            ramp_z = 0.04_dp + (q(j) - 0.2_dp)*(0.02_dp - 0.04_dp)/(0.6_dp - 0.2_dp)
         else
            ramp_z = old_xa(i_metals,j)
         end if
         ok = ok .and. abs(xa(i_metals,j) - ramp_z) < 1.0e-14_dp .and. &
              abs(xa(i_c12,j)/old_xa(i_c12,j)/(ramp_z/old_xa(i_metals,j)) - 1.0_dp) < tol
         if (q(j) < 0.6_dp) ok = ok .and. abs(xa(i_h1,j) + xa(i_he4,j) + xa(i_metals,j) + &
              xa(i_he3,j) - 1.0_dp) < 1.0e-15_dp
      end do
      call model_kernel_check('ramp_metals: the core value, the ramp, the surface kept; X = 1 - Y - Z - He3', &
           ok .and. count(q < 0.2_dp) > 3 .and. count(q >= 0.6_dp) > 3)
end subroutine rescaling_checks

! ---------------------------------------------------------------
! mesh_spacing: the rezone's new mesh.
subroutine mesh_checks()
      use mesh_spacing, only: chi_scales, accumulate_chi, sort_kept_points, &
           space_points_along_chi, spline_onto_points, insert_points_at_jumps, &
           remove_crowded_points
      real(dp), parameter :: tol = 1.0e-12_dp
      type(chi_scales) :: scales
      real(dp) :: log_mass(n), lum(n), log_pressure(n), chi(n), new_log_mass(200), &
           small_mesh(50), values(200), mesh(12), steps(6), crowded(6)
      integer :: kept(5), triple(4), count, new_count, needed, ierr, j
      logical :: ok

! 40 points: log m rises by 0.1 a point, log P falls by 0.05, L rises by
! 1 up to point 20 and is flat beyond. The fine zone is within 0.275 dex
! of log P(30): points 25 to 35.
      log_mass = [(30.0_dp + 0.1_dp*real(j - 1, dp), j = 1, n)]
      log_pressure = [(17.0_dp - 0.05_dp*real(j - 1, dp), j = 1, n)]
      lum = [(real(min(j, 20), dp), j = 1, n)]
      scales = chi_scales(log_mass=0.2_dp, luminosity=4.0_dp, log_pressure_core=0.1_dp, &
           log_pressure_base=0.02_dp, log_pressure_envelope=0.25_dp, base_fine_width=0.275_dp)
      call accumulate_chi(log_mass, lum, log_pressure, log_pressure(30), scales, chi)
      call model_kernel_check('accumulate_chi: chi(1) = 1; the log m and L terms where L rises', &
           chi(1) == 1.0_dp .and. abs(chi(2) - chi(1) - (0.5_dp + 0.25_dp + 0.5_dp)) < tol .and. &
           abs(chi(22) - chi(21) - (0.5_dp + 0.5_dp)) < tol)
      call model_kernel_check('accumulate_chi: the log P term below, inside and above the fine zone', &
           abs(chi(24) - chi(23) - (0.5_dp + 0.5_dp)) < tol .and. &
           abs(chi(25) - chi(24) - (0.5_dp + 2.5_dp)) < tol .and. &
           abs(chi(35) - chi(34) - (0.5_dp + 2.5_dp)) < tol .and. &
           abs(chi(36) - chi(35) - (0.5_dp + 0.2_dp)) < tol)

! The kept points.
      kept = [7, 3, 40, 3, 12]
      count = 5
      call sort_kept_points(kept, count)
      call model_kernel_check('sort_kept_points: ascending, a duplicate removed', &
           count == 4 .and. all(kept(1:4) == [3, 7, 12, 40]))
! Three equal points: the scan moves on after a removal, so two are left
! (reported with this round; the check holds today's result).
      triple = [5, 9, 5, 5]
      count = 4
      call sort_kept_points(triple, count)
      call model_kernel_check('sort_kept_points: of three equal points two are left', &
           count == 3 .and. all(triple(1:3) == [5, 5, 9]))

! The new mesh: chi rises by 2.5 a point, so log m is linear in chi;
! points 11 and 40 are kept. The segments are 25 and 72.5 long in chi:
! 26 and 73 intervals.
      chi = [(1.0_dp + 2.5_dp*real(j - 1, dp), j = 1, n)]
      call space_points_along_chi(chi, log_mass, [11, 40], new_log_mass, new_count, needed, ierr)
      ok = ierr == 0 .and. needed == 0 .and. new_count == 1 + 26 + 73
      if (ok) then
         ok = all(abs(new_log_mass(2:27) - new_log_mass(1:26) - 1.0_dp/26.0_dp) < 1.0e-11_dp) &
              .and. all(abs(new_log_mass(28:100) - new_log_mass(27:99) - 2.9_dp/73.0_dp) &
              < 1.0e-11_dp)
      end if
      call model_kernel_check('space_points_along_chi: int(Delta chi) + 1 even intervals a segment', ok)
      call model_kernel_check('space_points_along_chi: the first and the kept points are on the mesh', &
           new_log_mass(1) == log_mass(1) .and. abs(new_log_mass(27) - log_mass(11)) < 1.0e-11_dp &
           .and. abs(new_log_mass(100) - log_mass(40)) < 1.0e-11_dp)
      call space_points_along_chi(chi, log_mass, [11, 40], small_mesh, new_count, needed, ierr)
      call model_kernel_check('space_points_along_chi: the count needed when the arrays are too small', &
           needed == 100 .and. ierr == 0)

! A profile linear in log m at the new points.
      call space_points_along_chi(chi, log_mass, [11, 40], new_log_mass, new_count, needed, ierr)
      call spline_onto_points(log_mass, 0.7_dp - 0.01_dp*(log_mass - 30.0_dp), &
           new_log_mass(1:new_count), values(1:new_count), ierr)
      call model_kernel_check('spline_onto_points: a profile linear in log m is reproduced', &
           ierr == 0 .and. all(abs(values(1:new_count) - (0.7_dp - 0.01_dp* &
           (new_log_mass(1:new_count) - 30.0_dp))) < 1.0e-12_dp))

! Points at jumps: six points at log m = 1..6; the profile rises by 0.25
! into point 3 and by 0.34 into point 6; the limit is 0.1.
      steps = [0.1_dp, 0.1_dp, 0.35_dp, 0.36_dp, 0.36_dp, 0.7_dp]
      mesh = 0.0_dp
      mesh(1:6) = [(real(j, dp), j = 1, 6)]
      new_count = 6
      call insert_points_at_jumps(steps, .false., 0.1_dp, mesh, new_count, needed)
      call model_kernel_check('insert_points_at_jumps: int(step/max) points, evenly in log m', &
           needed == 0 .and. new_count == 11 .and. &
           all(abs(mesh(1:11) - [1.0_dp, 2.0_dp, 2.0_dp + 1.0_dp/3.0_dp, 2.0_dp + 2.0_dp/3.0_dp, &
           3.0_dp, 4.0_dp, 5.0_dp, 5.25_dp, 5.5_dp, 5.75_dp, 6.0_dp]) < tol))
! The same profile falling outward: nothing with the rise, the same
! points with falling.
      mesh(1:6) = [(real(j, dp), j = 1, 6)]
      new_count = 6
      call insert_points_at_jumps(-steps, .false., 0.1_dp, mesh, new_count, needed)
      ok = new_count == 6
      call insert_points_at_jumps(-steps, .true., 0.1_dp, mesh, new_count, needed)
      call model_kernel_check('insert_points_at_jumps: a fall counts with falling only', &
           ok .and. new_count == 11 .and. abs(mesh(8) - 5.25_dp) < tol)
      mesh(1:6) = [(real(j, dp), j = 1, 6)]
      new_count = 6
      call insert_points_at_jumps(steps, .false., 0.0_dp, mesh, new_count, needed)
      ok = new_count == 6 .and. needed == 0
      call insert_points_at_jumps(steps, .false., 0.1_dp, mesh(1:8), new_count, needed)
      call model_kernel_check('insert_points_at_jumps: none with the limit off; the count needed', &
           ok .and. needed == 9)

! Crowded points: within 0.01 of the last point kept; the last stays.
      crowded = [1.0_dp, 1.005_dp, 1.02_dp, 1.021_dp, 1.5_dp, 1.501_dp]
      new_count = 6
      call remove_crowded_points(crowded, new_count, 0.01_dp)
      call model_kernel_check('remove_crowded_points: within the limit of the last kept; the ends stay', &
           new_count == 4 .and. all(crowded(1:4) == [1.0_dp, 1.02_dp, 1.5_dp, 1.501_dp]))
end subroutine mesh_checks

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

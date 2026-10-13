!----------------------------------------------------------------------
! test_regrid
!----------------------------------------------------------------------
! The rezone's regridding kernels (core/model/regridding.f90) on a
! synthetic 1 Msun model.
!
! The old mesh has 60 points spaced evenly in log10(q/(1-q)), q = m/M,
! from q = 1e-4 to the outer point at q = 1 - 2e-4. Its composition
! (X, Y, Z and He3, with Y = 1 - X - Z) has a mixed convective core of
! X = 0.35 out to q = 0.08 and a step up to X = 0.55 at the core edge,
! between two old points, with a gradient above it; its rotation is
! solid in the core and slower outward, with j = (2/3) omega r^2. The
! new mesh is laid out as rezone lays out its points: it keeps the
! centre, the outer point and the point just above the step (a flagged
! point), has fewer points in the core and more above it, spaced evenly
! in the same coordinate at other intervals so that they fall between
! the old points, and three more points inside the step's interval,
! where add_points_for_x_jumps inserts them.
!
! The checks, each printed with its value and bound:
!   shell masses   on both meshes, sum dm = M to round-off (the weights
!                  telescope; bound n eps, n the number of points); on
!                  the new mesh, the centre, interior and outer rules of
!                  compute_shell_masses recomputed here (bound eps);
!   consistency    profiles linear in log m are reproduced at the new
!                  points (structure, composition and rotation; bound
!                  1e-13 of the profile's largest value, a few hundred
!                  round-offs: log m is near 33 and the steps between
!                  points near the surface near 1e-5); the points the
!                  two meshes share keep their values (bound eps);
!   composition    the relative change of each species' mass,
!                  sum_i X_i dm_i, from the old mesh to the new; the
!                  largest excursion of a mass fraction outside its old
!                  range, over its largest value (the spline is
!                  shape-preserving; bound 4 eps); the largest
!                  |X + Y + Z - 1| at the new points;
!   rotation       compute_angular_momentum_totals' J and rotational
!                  kinetic energy against the sums sum_i j_i dm_i and
!                  sum_i omega_i j_i dm_i/2 formed here (bound n eps);
!                  the relative change of each from the old mesh to the
!                  new.
! The regridding conserves none of the sums exactly, so their changes are
! measured and pinned (pin_*): a check fails when a change grows beyond
! its pin by more than 10%. They have two sources. A smooth profile's
! sum differs between the two meshes by their quadrature errors, which
! go as h^2: J changes by 1.3e-3 here and Z by 1.1e-5, and both fall
! about fourfold with twice the points on both meshes (the solar decks'
! models have about 1000). A step adds up to Delta X times half the mass of the old
! interval holding it, since the old mesh's sum puts the step at the
! interval's midpoint and the new points inside take interpolated
! values: here 0.2 x 0.022 M/2, 3e-3 of the hydrogen; X changes by
! 1.7e-4.
!
! Usage: test_regrid (no arguments). Prints, after the line
! 'test_regrid report', one line per check; exits with status 1 when a
! check fails.
program test_regrid
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use math_lib
      use physical_constants, only: ln10
      use regridding, only: structure_profiles, rotation_profiles, &
           compute_shell_masses, regrid_profile, regrid_composition, &
           regrid_structure, regrid_rotation, compute_angular_momentum_totals
      implicit none

      integer, parameter :: n_old = 60, n_species = 4
      integer, parameter :: i_x = 1, i_y = 2, i_z = 3, i_he3 = 4
! The model: total mass (g), surface radius (cm), central angular
! velocity (rad/s), the core edge in q, the mass fractions on either
! side of the step, and the q of the inner and outer points.
      real(dp), parameter :: total_mass = 1.989d33, surface_radius = 6.96d10, &
           omega_core = 2.9d-6, q_core_edge = 0.08d0, x_core = 0.35d0, &
           x_above_step = 0.55d0, q_inner = 1.0d-4, q_outer = 1.0d0 - 2.0d-4
! The new mesh: points from the centre to the point above the step, then
! from there to the outer point (both counts include the ends), and the
! points inserted in the step's interval.
      integer, parameter :: n_core_new = 13, n_envelope_new = 50, n_step_new = 3
! The pinned relative changes of the sums, |new - old|/old (the
! measured values, rounded up to three digits): of each species' mass,
! of X + Y + Z at the worst new point, of J and of the rotational kinetic
! energy. A check fails beyond 1.1 times its pin.
      real(dp), parameter :: pin_mass_change(n_species) = &
           [1.66d-4, 3.35d-4, 1.14d-5, 9.23d-4]
      real(dp), parameter :: pin_abundance_sum = 2.88d-6, &
           pin_angular_momentum = 1.35d-3, pin_rotational_ke = 1.12d-3
      character(len=3), parameter :: species_name(n_species) = ['X  ', 'Y  ', 'Z  ', 'He3']

! log_total_mass, and the total mass as compute_shell_masses forms it,
! 10^log_total_mass (which differs from total_mass in the last bits).
      real(dp) :: log_total_mass, mass_from_log, eps
      real(dp) :: old_log_mass(n_old), old_q(n_old), old_mass(n_old), old_dm(n_old)
      real(dp), allocatable :: new_log_mass(:), new_mass(:), new_dm(:)
      integer :: n_new, k_edge, nfail

      nfail = 0
      eps = epsilon(1.0d0)
      log_total_mass = log10(total_mass)
      mass_from_log = exp(ln10*log_total_mass)
      call build_old_mesh()
      call build_new_mesh()
      write(*,'(a)') 'test_regrid report'
      write(*,'(a,i4,a,i4,a,i3,a,i3)') 'old mesh', n_old, ' points, new mesh', n_new, &
           ' points; X steps between old points', k_edge, ' and', k_edge + 1
      write(*,'(2x,a5,43x,a9,a9,2x,a)') 'check', 'value', 'bound', 'status'
      call check_shell_masses()
      call check_linear_profiles()
      call check_species_masses()
      call check_rotation_sums()
      if (nfail > 0) then
         write(*,'(a,i4,a)') 'test_regrid: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_regrid: done'

contains

! ---------------------------------------------------------------
! The old mesh: n_old points evenly spaced in log10(q/(1-q)); k_edge is
! the last point inside the convective core (q <= q_core_edge).
subroutine build_old_mesh()
      real(dp) :: s_inner, s_outer, s
      integer :: i
      s_inner = log10(q_inner/(1.0d0 - q_inner))
      s_outer = log10(q_outer/(1.0d0 - q_outer))
      k_edge = 0
      do i = 1, n_old
         s = s_inner + (s_outer - s_inner)*dble(i - 1)/dble(n_old - 1)
         old_q(i) = 1.0d0/(1.0d0 + exp10(-s))
         old_log_mass(i) = log_total_mass + log10(old_q(i))
         if (old_q(i) <= q_core_edge) k_edge = i
      end do
      call compute_shell_masses(old_log_mass, log_total_mass, old_mass, old_dm)
end subroutine build_old_mesh

! ---------------------------------------------------------------
! The new mesh: n_core_new points evenly spaced in log10(q/(1-q)) from
! the centre to the old point k_edge + 1, n_step_new points inserted
! evenly in log m between the last two of them (the interval holding
! the step), and n_envelope_new points evenly spaced from the old point
! k_edge + 1 to the outer point.
subroutine build_new_mesh()
      real(dp) :: s_lo, s_hi, s, step_lo, step_hi
      integer :: i, j
      n_new = n_core_new + n_step_new + n_envelope_new - 1
      allocate(new_log_mass(n_new), new_mass(n_new), new_dm(n_new))
      j = 0
      s_lo = log10(old_q(1)/(1.0d0 - old_q(1)))
      s_hi = log10(old_q(k_edge+1)/(1.0d0 - old_q(k_edge+1)))
      do i = 1, n_core_new - 1
         s = s_lo + (s_hi - s_lo)*dble(i - 1)/dble(n_core_new - 1)
         j = j + 1
         new_log_mass(j) = log_total_mass + log10(1.0d0/(1.0d0 + exp10(-s)))
      end do
      new_log_mass(1) = old_log_mass(1)
      step_lo = new_log_mass(j)
      step_hi = old_log_mass(k_edge+1)
      do i = 1, n_step_new
         j = j + 1
         new_log_mass(j) = step_lo + (step_hi - step_lo)*dble(i)/dble(n_step_new + 1)
      end do
      s_lo = s_hi
      s_hi = log10(old_q(n_old)/(1.0d0 - old_q(n_old)))
      do i = 1, n_envelope_new
         s = s_lo + (s_hi - s_lo)*dble(i - 1)/dble(n_envelope_new - 1)
         j = j + 1
         new_log_mass(j) = log_total_mass + log10(1.0d0/(1.0d0 + exp10(-s)))
      end do
      new_log_mass(j - n_envelope_new + 1) = old_log_mass(k_edge+1)
      new_log_mass(n_new) = old_log_mass(n_old)
      call compute_shell_masses(new_log_mass, log_total_mass, new_mass, new_dm)
end subroutine build_new_mesh

! ---------------------------------------------------------------
! The composition at q = m/M: X, Y = 1 - X - Z, Z and He3.
function composition_at(q, in_core) result(xa)
      real(dp), intent(in) :: q
      logical, intent(in) :: in_core
      real(dp) :: xa(n_species)
      if (in_core) then
         xa(i_x) = x_core
      else
         xa(i_x) = x_above_step + 0.15d0*(1.0d0 - exp(-(q - q_core_edge)/0.1d0))
      end if
      xa(i_z) = 0.019d0 - 0.001d0*q*q
      xa(i_y) = 1.0d0 - xa(i_x) - xa(i_z)
      xa(i_he3) = 1.0d-5 + 3.0d-4*exp(-((q - 0.55d0)/0.12d0)**2)
end function composition_at

! ---------------------------------------------------------------
! The rotation at q = m/M: omega (solid in the core, slower
! outward), r, and j = (2/3) omega r^2 (a spherical shell's).
subroutine rotation_at(q, omega, radius, j)
      real(dp), intent(in) :: q
      real(dp), intent(out) :: omega, radius, j
      if (q <= q_core_edge) then
         omega = omega_core
      else
         omega = omega_core*(0.4d0 + 0.6d0*exp(-(q - q_core_edge)/0.3d0))
      end if
      radius = surface_radius*pow(q, 1.0d0/3.0d0)*(0.3d0 + 0.7d0*q**3)
      j = (2.0d0/3.0d0)*omega*radius**2
end subroutine rotation_at

! ---------------------------------------------------------------
! Print one check: ok when |value| <= bound.
subroutine report(label, value, bound)
      character(len=*), intent(in) :: label
      real(dp), intent(in) :: value, bound
      character(len=4) :: status
      character(len=48) :: padded
      padded = label
      if (abs(value) <= bound) then
         status = 'ok'
      else
         status = 'FAIL'
         nfail = nfail + 1
      end if
      write(*,'(2x,a48,1pe9.1,1pe9.1,2x,a)') padded, value, bound, trim(status)
end subroutine report

! ---------------------------------------------------------------
! The shell masses: the sums on both meshes, and the rules on the new.
subroutine check_shell_masses()
      real(dp) :: worst
      integer :: i
      call report('shell masses, old mesh: sum dm/M - 1', &
           sum(old_dm)/mass_from_log - 1.0d0, dble(n_old)*eps)
      call report('shell masses, new mesh: sum dm/M - 1', &
           sum(new_dm)/mass_from_log - 1.0d0, dble(n_new)*eps)
      call report('centre rule: dm_1 = (m_1 + m_2)/2', &
           (new_dm(1) - 0.5d0*(new_mass(1) + new_mass(2)))/new_dm(1), eps)
      worst = 0.0d0
      do i = 2, n_new - 1
         worst = max(worst, abs(new_dm(i) - 0.5d0*(new_mass(i+1) - new_mass(i-1)))/new_dm(i))
      end do
      call report('interior rule: dm_i = (m_i+1 - m_i-1)/2', worst, eps)
      call report('outer rule: dm_n = M - (m_n-1 + m_n)/2', &
           (new_dm(n_new) - (mass_from_log - 0.5d0*(new_mass(n_new-1) + new_mass(n_new)))) &
           /new_dm(n_new), eps)
end subroutine check_shell_masses

! ---------------------------------------------------------------
! Profiles linear in log m, a + b (log m - log m_1), with a different
! (a, b) per variable, through each kernel: the largest error at the new
! points over all its variables, relative to the variable's largest
! value. Then the points the meshes share (the centre, the point above
! the step and the outer point) keep the old values of the real model.
subroutine check_linear_profiles()
      type(structure_profiles) :: old_structure, new_structure
      type(rotation_profiles) :: old_rotation, new_rotation
      real(dp) :: old_xa(n_species, n_old), new_xa(n_species, n_new)
      real(dp) :: old_values(n_old), new_values(n_new), omega(n_old), &
           radius(n_old), worst
      integer :: k
      old_structure = structure_profiles(log_pressure=linear(old_log_mass, 17.3d0, -1.1d0), &
           log_radius=linear(old_log_mass, 9.0d0, 0.4d0), &
           luminosity_lsun=linear(old_log_mass, 0.2d0, 0.25d0), &
           log_temperature=linear(old_log_mass, 7.2d0, -0.3d0), &
           log_density=linear(old_log_mass, 2.1d0, -0.9d0))
      call regrid_structure(old_log_mass, old_structure, new_log_mass, new_structure)
      worst = max(linear_error(new_structure%log_pressure, 17.3d0, -1.1d0), &
           linear_error(new_structure%log_radius, 9.0d0, 0.4d0), &
           linear_error(new_structure%luminosity_lsun, 0.2d0, 0.25d0), &
           linear_error(new_structure%log_temperature, 7.2d0, -0.3d0), &
           linear_error(new_structure%log_density, 2.1d0, -0.9d0))
      call report('structure: linear in log m reproduced', worst, 1.0d-13)
      do k = 1, n_species
         old_xa(k,:) = linear(old_log_mass, 0.1d0*k, 0.02d0*(k - 2))
      end do
      call regrid_composition(old_log_mass, old_xa, new_log_mass, new_xa)
      worst = 0.0d0
      do k = 1, n_species
         worst = max(worst, linear_error(new_xa(k,:), 0.1d0*k, 0.02d0*(k - 2)))
      end do
      call report('composition: linear in log m reproduced', worst, 1.0d-13)
      old_rotation = rotation_profiles(omega=linear(old_log_mass, 3.0d-6, -1.0d-6), &
           specific_angular_momentum=linear(old_log_mass, 1.0d15, 4.0d14), &
           fp=linear(old_log_mass, 1.0d0, -1.0d-3), ft=linear(old_log_mass, 1.0d0, -5.0d-4), &
           mean_radius=linear(old_log_mass, 3.0d10, 2.0d10), &
           eta_squared=linear(old_log_mass, 0.05d0, 0.02d0))
      call regrid_rotation(old_log_mass, old_rotation, new_log_mass, new_rotation)
      worst = max(linear_error(new_rotation%omega, 3.0d-6, -1.0d-6), &
           linear_error(new_rotation%specific_angular_momentum, 1.0d15, 4.0d14), &
           linear_error(new_rotation%fp, 1.0d0, -1.0d-3), &
           linear_error(new_rotation%ft, 1.0d0, -5.0d-4), &
           linear_error(new_rotation%mean_radius, 3.0d10, 2.0d10), &
           linear_error(new_rotation%eta_squared, 0.05d0, 0.02d0))
      call report('rotation: linear in log m reproduced', worst, 1.0d-13)
! The shared points, with the model's X (which has the step) and j.
      call model_composition(old_xa)
      call regrid_profile(old_log_mass, old_xa(i_x,:), new_log_mass, new_values)
      worst = shared_point_error(old_xa(i_x,:), new_values)
      call model_rotation(omega, radius, old_values)
      call regrid_profile(old_log_mass, old_values, new_log_mass, new_values)
      worst = max(worst, shared_point_error(old_values, new_values))
      call report('shared points keep their values (X, j)', worst, eps)
end subroutine check_linear_profiles

! a + b (log m - log m_1) at the points log_mass.
function linear(log_mass, a, b) result(values)
      real(dp), intent(in) :: log_mass(:), a, b
      real(dp) :: values(size(log_mass))
      values = a + b*(log_mass - old_log_mass(1))
end function linear

! The largest |values - linear| at the new points over the largest
! |linear| there.
real(dp) function linear_error(values, a, b)
      real(dp), intent(in) :: values(:), a, b
      real(dp) :: exact(n_new)
      exact = linear(new_log_mass, a, b)
      linear_error = maxval(abs(values - exact))/maxval(abs(exact))
end function linear_error

! The largest relative change at the centre, at the point above the step
! and at the outer point, which both meshes have.
real(dp) function shared_point_error(old_values, new_values)
      real(dp), intent(in) :: old_values(:), new_values(:)
      shared_point_error = max(abs(new_values(1) - old_values(1))/abs(old_values(1)), &
           abs(new_values(n_core_new + n_step_new) - old_values(k_edge+1)) &
           /abs(old_values(k_edge+1)), &
           abs(new_values(n_new) - old_values(n_old))/abs(old_values(n_old)))
end function shared_point_error

! ---------------------------------------------------------------
! The model's composition, and its omega, r and j, at the old points.
subroutine model_composition(xa)
      real(dp), intent(out) :: xa(:,:)
      integer :: i
      do i = 1, n_old
         xa(:,i) = composition_at(old_q(i), i <= k_edge)
      end do
end subroutine model_composition

subroutine model_rotation(omega, radius, j_rot)
      real(dp), intent(out) :: omega(:), radius(:), j_rot(:)
      integer :: i
      do i = 1, n_old
         call rotation_at(old_q(i), omega(i), radius(i), j_rot(i))
      end do
end subroutine model_rotation

! ---------------------------------------------------------------
! The composition through regrid_composition: each species' mass before
! and after, the excursions outside the old ranges, and X + Y + Z.
subroutine check_species_masses()
      real(dp) :: old_xa(n_species, n_old), new_xa(n_species, n_new)
      real(dp) :: old_total, new_total, lo, hi, worst
      integer :: k
      call model_composition(old_xa)
      call regrid_composition(old_log_mass, old_xa, new_log_mass, new_xa)
      do k = 1, n_species
         old_total = sum(old_xa(k,:)*old_dm)
         new_total = sum(new_xa(k,:)*new_dm)
         call report('mass of '//trim(species_name(k))//': relative change', &
              (new_total - old_total)/old_total, 1.1d0*pin_mass_change(k))
      end do
      worst = 0.0d0
      do k = 1, n_species
         lo = minval(old_xa(k,:))
         hi = maxval(old_xa(k,:))
         worst = max(worst, maxval(max(new_xa(k,:) - hi, lo - new_xa(k,:), 0.0d0)) &
              /max(abs(lo), abs(hi)))
      end do
      call report('mass fractions outside their old range', worst, 4.0d0*eps)
      call report('max |X + Y + Z - 1| at the new points', &
           maxval(abs(new_xa(i_x,:) + new_xa(i_y,:) + new_xa(i_z,:) - 1.0d0)), &
           1.1d0*pin_abundance_sum)
end subroutine check_species_masses

! ---------------------------------------------------------------
! The rotation through regrid_rotation: J and the rotational kinetic
! energy before and after (compute_angular_momentum_totals).
subroutine check_rotation_sums()
      type(rotation_profiles) :: old_rotation, new_rotation
      real(dp) :: j_rot(n_old), omega(n_old), radius(n_old)
      real(dp) :: old_ke(n_old), new_ke(n_new)
      real(dp) :: old_j_total, old_ke_total, new_j_total, new_ke_total
      call model_rotation(omega, radius, j_rot)
      old_rotation = rotation_profiles(omega=omega, specific_angular_momentum=j_rot, &
           fp=1.0d0 - 1.0d-3*old_q, ft=1.0d0 - 5.0d-4*old_q, mean_radius=radius, &
           eta_squared=0.05d0*old_q)
      call regrid_rotation(old_log_mass, old_rotation, new_log_mass, new_rotation)
      call compute_angular_momentum_totals(old_rotation%specific_angular_momentum, &
           old_rotation%omega, old_dm, old_ke, old_j_total, old_ke_total)
      call compute_angular_momentum_totals(new_rotation%specific_angular_momentum, &
           new_rotation%omega, new_dm, new_ke, new_j_total, new_ke_total)
      call report('J and E_rot totals = sum j dm, sum omega j dm/2', &
           max(abs(new_j_total - sum(new_rotation%specific_angular_momentum*new_dm)) &
           /new_j_total, abs(new_ke_total - sum(0.5d0*new_rotation%omega &
           *new_rotation%specific_angular_momentum*new_dm))/new_ke_total), dble(n_new)*eps)
      call report('angular momentum J: relative change', &
           (new_j_total - old_j_total)/old_j_total, 1.1d0*pin_angular_momentum)
      call report('rotational kinetic energy: relative change', &
           (new_ke_total - old_ke_total)/old_ke_total, 1.1d0*pin_rotational_ke)
end subroutine check_rotation_sums

end program test_regrid

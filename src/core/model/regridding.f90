!----------------------------------------------------------------------
! regridding
!----------------------------------------------------------------------
! The arithmetic of moving a model onto a new mesh: a star-blind kernel
! (docs/GUIDELINES.md, star rule 6). rezone (rezoning.f90) decides where
! the new points go, gathers the model from the star, calls these, and
! stores the results; nothing here reads the star or keeps state between
! calls.
!
! A mesh is a run of points in log m (log10 of the mass inside the
! point, g), increasing outward; the model's variables are profiles on
! it, one value per point. The procedures:
!   compute_shell_masses   the points' masses m_i and shell masses dm_i
!   regrid_profile         one profile onto the new points
!   regrid_composition     every species of a composition matrix
!   regrid_structure       log P, log R, L, log T and log rho
!   regrid_rotation        omega, j, fp, ft, r0 and eta2
!   compute_angular_momentum_totals
!                          J = sum_i j_i dm_i and the rotational
!                          kinetic energy, shell by shell and in total
!
! Each profile is interpolated on its own, against log m, with the
! numerics domain's shape-preserving osculatory quadratic spline
! (osplin, ACM 574). It keeps the data's monotonicity and convexity, so
! a step in the composition gains no new extremum and every mass
! fraction stays inside the range of its old values (to round-off); it
! passes through the old points and reproduces a profile linear in
! log m. It does not conserve the mass of a species, sum_i X_i dm_i, or
! the angular momentum, sum_i j_i dm_i: those change by the difference
! between the two meshes' sums over the same profile, which goes as the
! square of the spacing for a smooth profile, plus, where a profile
! jumps between two old points (a step at a convective-core edge), up to
! the jump times half the mass between them, as the new points inside
! take interpolated values. Nor does it keep the mass fractions summing
! to one. core/test/test_regrid measures all three on a synthetic model.
module regridding
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: structure_profiles, rotation_profiles
      public :: compute_shell_masses, regrid_profile, regrid_composition, &
           regrid_structure, regrid_rotation, compute_angular_momentum_totals

! The structure variables at the points of a mesh: log P (P in
! dyn/cm^2), log R (cm), L (Lsun), log T (K) and log rho (g/cm^3).
      type :: structure_profiles
         real(dp), allocatable :: log_pressure(:), log_radius(:), &
              luminosity_lsun(:), log_temperature(:), log_density(:)
      end type structure_profiles

! The rotation variables at the points of a mesh: the angular velocity
! omega (rad/s), the specific angular momentum j (cm^2/s), the
! rotational factors fp and ft of the structure equations, and the
! mean radius r0 (cm) and distortion eta2 of the level surface.
      type :: rotation_profiles
         real(dp), allocatable :: omega(:), specific_angular_momentum(:), &
              fp(:), ft(:), mean_radius(:), eta_squared(:)
      end type rotation_profiles

contains

!----------------------------------------------------------------------
! compute_shell_masses
!----------------------------------------------------------------------
! The point masses m_i = 10^(log m_i) (g) of the mesh log_mass(1:nz) and
! their mass weights (g)
!   dm_i  = (m_(i+1) - m_(i-1))/2,  with m_0 = -m_1,
!   dm_nz = M - (m_(nz-1) + m_nz)/2,
! so dm_i is the mass between the midpoints on either side of point i
! (from the centre for i = 1), and dm_nz includes the envelope above the
! outer point (M = 10^log_total_mass). The weights telescope: their sum
! is M. mass and shell_mass have the size of log_mass.
pure subroutine compute_shell_masses(log_mass, log_total_mass, mass, shell_mass)
      use math_lib
      use physical_constants, only: ln10
      real(dp), intent(in) :: log_mass(:), log_total_mass
      real(dp), intent(out) :: mass(:), shell_mass(:)
      real(dp) :: mass_two_back, mass_prev, mass_curr
      integer :: i, nz

      nz = size(log_mass)
      mass_curr = exp(ln10*log_mass(1))
      mass_prev = - mass_curr
      do i = 2,nz
         mass_two_back = mass_prev
         mass_prev = mass_curr
         mass_curr = exp(ln10*log_mass(i))
         mass(i-1) = mass_prev
         shell_mass(i-1) = 0.5d0*(mass_curr-mass_two_back)
      end do
      mass(nz) = mass_curr
      shell_mass(nz) = exp(ln10*log_total_mass) - &
           0.5d0*(mass_prev+mass_curr)
end subroutine compute_shell_masses

!----------------------------------------------------------------------
! regrid_profile
!----------------------------------------------------------------------
! One profile onto new points: old_values at the old mesh old_log_mass
! (increasing, at least 3 points, the same size) evaluated at
! new_log_mass (nondecreasing) with osplin, into new_values (the size of
! new_log_mass). New points outside the old mesh get the end intervals'
! quadratics extrapolated; rezone's new mesh starts and ends on the old
! one's end points.
pure subroutine regrid_profile(old_log_mass, old_values, new_log_mass, new_values)
      use numerics_lib, only: osplin
      real(dp), intent(in) :: old_log_mass(:), old_values(:), new_log_mass(:)
      real(dp), intent(out) :: new_values(:)

      call osplin(new_log_mass, new_values, old_log_mass, old_values, &
           size(old_log_mass), size(new_log_mass))
end subroutine regrid_profile

!----------------------------------------------------------------------
! regrid_composition
!----------------------------------------------------------------------
! A composition matrix xa(species, point) onto new points: each species
! (row) of old_xa, on the old mesh old_log_mass, through regrid_profile
! at new_log_mass, into the same row of new_xa (as many rows as old_xa,
! a column per new point).
pure subroutine regrid_composition(old_log_mass, old_xa, new_log_mass, new_xa)
      real(dp), intent(in) :: old_log_mass(:), old_xa(:,:), new_log_mass(:)
      real(dp), intent(out) :: new_xa(:,:)
      integer :: k

      do k = 1,size(old_xa,1)
         call regrid_profile(old_log_mass, old_xa(k,:), new_log_mass, new_xa(k,:))
      end do
end subroutine regrid_composition

!----------------------------------------------------------------------
! regrid_structure
!----------------------------------------------------------------------
! The structure variables on the old mesh (old, at old_log_mass) at the
! new points new_log_mass, each through regrid_profile, into new.
pure subroutine regrid_structure(old_log_mass, old, new_log_mass, new)
      real(dp), intent(in) :: old_log_mass(:), new_log_mass(:)
      type(structure_profiles), intent(in) :: old
      type(structure_profiles), intent(out) :: new
      integer :: n

      n = size(new_log_mass)
      allocate(new%log_pressure(n), new%log_radius(n), new%luminosity_lsun(n), &
           new%log_temperature(n), new%log_density(n))
      call regrid_profile(old_log_mass, old%log_pressure, new_log_mass, new%log_pressure)
      call regrid_profile(old_log_mass, old%log_radius, new_log_mass, new%log_radius)
      call regrid_profile(old_log_mass, old%luminosity_lsun, new_log_mass, &
           new%luminosity_lsun)
      call regrid_profile(old_log_mass, old%log_temperature, new_log_mass, &
           new%log_temperature)
      call regrid_profile(old_log_mass, old%log_density, new_log_mass, new%log_density)
end subroutine regrid_structure

!----------------------------------------------------------------------
! regrid_rotation
!----------------------------------------------------------------------
! The rotation variables on the old mesh (old, at old_log_mass) at the
! new points new_log_mass, each through regrid_profile, into new. j is
! the quantity carried over: it is splined so that the angular momentum
! J = sum_i j_i dm_i changes only by the regridding of its profile. The
! others are starting values on the new mesh: rezone re-derives omega
! from j (omega_from_j, with the moments of inertia of the new shells)
! and then fp, ft, r0 and eta2 for that omega.
pure subroutine regrid_rotation(old_log_mass, old, new_log_mass, new)
      real(dp), intent(in) :: old_log_mass(:), new_log_mass(:)
      type(rotation_profiles), intent(in) :: old
      type(rotation_profiles), intent(out) :: new
      integer :: n

      n = size(new_log_mass)
      allocate(new%omega(n), new%specific_angular_momentum(n), new%fp(n), &
           new%ft(n), new%mean_radius(n), new%eta_squared(n))
      call regrid_profile(old_log_mass, old%omega, new_log_mass, new%omega)
      call regrid_profile(old_log_mass, old%specific_angular_momentum, new_log_mass, &
           new%specific_angular_momentum)
      call regrid_profile(old_log_mass, old%fp, new_log_mass, new%fp)
      call regrid_profile(old_log_mass, old%ft, new_log_mass, new%ft)
      call regrid_profile(old_log_mass, old%mean_radius, new_log_mass, new%mean_radius)
      call regrid_profile(old_log_mass, old%eta_squared, new_log_mass, new%eta_squared)
end subroutine regrid_rotation

!----------------------------------------------------------------------
! compute_angular_momentum_totals
!----------------------------------------------------------------------
! The angular momentum and rotational kinetic energy of a mesh's shells,
! from the specific angular momentum j (cm^2/s), the angular velocity
! omega (rad/s) and the shell masses dm (g):
!   J = sum_i j_i dm_i (g cm^2/s),
!   E_rot = sum_i E_i,  E_i = (1/2) omega_i j_i dm_i (erg),
! summed from the centre out. rotational_ke (the size of shell_mass)
! returns the E_i.
pure subroutine compute_angular_momentum_totals(specific_angular_momentum, omega, &
     shell_mass, rotational_ke, total_angular_momentum, total_rotational_ke)
      real(dp), intent(in) :: specific_angular_momentum(:), omega(:), shell_mass(:)
      real(dp), intent(out) :: rotational_ke(:)
      real(dp), intent(out) :: total_angular_momentum, total_rotational_ke
      real(dp) :: shell_angular_momentum
      integer :: i

      total_angular_momentum = 0.0d0
      total_rotational_ke = 0.0d0
      do i = 1,size(shell_mass)
         shell_angular_momentum = specific_angular_momentum(i)*shell_mass(i)
         rotational_ke(i) = 0.5d0*omega(i)*shell_angular_momentum
         total_angular_momentum = total_angular_momentum + shell_angular_momentum
         total_rotational_ke = total_rotational_ke + rotational_ke(i)
      end do
end subroutine compute_angular_momentum_totals

end module regridding

!----------------------------------------------------------------------
! envelope_refit_kernel
!----------------------------------------------------------------------
! The arithmetic of moving the fitting point outward, star-free: the
! steps that turn the points of a fresh envelope integration
! (envelope_points' envelope_structure_state, from the centre outward)
! into interior shells. The drivers (rebuild_envelope, and
! read_starting_model through append_envelope_points) read the star and
! call these with its arrays:
!   shift_envelope_onto_shell   make the envelope continuous with the
!                               last interior shell and drop its points
!                               at or below that shell's mass
!                               (rebuild_envelope);
!   append_envelope_shells      the points inside the new fitting mass
!                               become shells, and the first one outside
!                               is interpolated onto the fitting mass
!                               (both drivers);
!   refit_shell_masses          the new shells' masses m and dm
!                               (rebuild_envelope);
!   guess_new_shell_rotation    a first guess at the new shells'
!                               rotation (rebuild_envelope).
! Units: log is log10; masses in g, radii in cm, omega in rad/s, the
! specific angular momentum j in cm^2/s and the specific moment of
! inertia i in cm^2. The kernels print nothing; the one failure
! (a degenerate interpolation interval) returns ierr = 1, which the
! driver reports. test_envelope_refit checks them.
module envelope_refit_kernel
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use envelope_points, only: envelope_structure_state
      use composition_layout, only: i_h1, i_he3, i_he4, i_metals
      implicit none
      private
      public :: shift_envelope_onto_shell, append_envelope_shells, &
           refit_shell_masses, guess_new_shell_rotation, between_points

! The smallest interval in log m (dex) the interpolation onto the
! fitting mass divides by; a narrower one is an error.
      real(dp), parameter, public :: min_interpolation_interval = 1.0d-14
! An envelope point is kept by shift_envelope_onto_shell when its log m
! exceeds the last kept one's by more than this (dex).
      real(dp), parameter, public :: min_log_mass_step = 1.0d-10

contains

! ---------------------------------------------------------------
! The linear interpolation between two points, v0 at fraction 0 and v2
! at fraction 1: v = v0 + fraction (v2 - v0).
elemental real(dp) function between_points(v0, v2, fraction)
      real(dp), intent(in) :: v0, v2, fraction
      between_points = v0 + fraction*(v2 - v0)
end function between_points

! ---------------------------------------------------------------
! Matches the envelope to the last interior shell and keeps, shifted
! down by one, the points above its mass. Point 1 of the envelope lies
! at the mass of the last interior shell; its differences from that
! shell (shell_logP, shell_logRho, shell_logT, shell_logR) in log P,
! log rho, log T and log r become offsets that make the envelope
! continuous with the interior. Every later point whose log m exceeds
! the last kept one's (at first the shell's, shell_log_mass) by more
! than min_log_mass_step moves down one slot, with the offsets added.
! Log P, log rho, log r, log(m/M), log T, X, Z and the convective flag
! move; the other arrays are not read for the appended points.
! num_env_points becomes the number kept. (The slot is the point's own
! less one, whatever was dropped before it.)
pure subroutine shift_envelope_onto_shell(env, shell_logP, shell_logRho, shell_logT, &
     shell_logR, shell_log_mass, log_total_mass)
      type(envelope_structure_state), intent(inout) :: env
      real(dp), intent(in) :: shell_logP, shell_logRho, shell_logT, shell_logR
      real(dp), intent(in) :: shell_log_mass, log_total_mass
      real(dp) :: pressure_offset, density_offset, temperature_offset, radius_offset
      real(dp) :: last_kept_log_mass, log_mass
      integer :: k, num_kept

      pressure_offset = shell_logP - env%env_log10_pressure(1)
      density_offset = shell_logRho - env%env_log10_density(1)
      temperature_offset = shell_logT - env%env_log10_temperature(1)
      radius_offset = shell_logR - env%env_log10_radius(1)
      last_kept_log_mass = shell_log_mass
      num_kept = 0
      do k = 1, env%num_env_points - 1
         log_mass = log_total_mass+env%env_log10_mass(k+1)
         if (log_mass-last_kept_log_mass > min_log_mass_step) then
            env%env_log10_density(k) = env%env_log10_density(k+1)+density_offset
            env%env_log10_pressure(k) = env%env_log10_pressure(k+1)+pressure_offset
            env%env_log10_radius(k) = env%env_log10_radius(k+1)+radius_offset
            env%env_log10_mass(k) = env%env_log10_mass(k+1)
            env%env_log10_temperature(k) = env%env_log10_temperature(k+1)+temperature_offset
            env%env_hydrogen_fraction(k) = env%env_hydrogen_fraction(k+1)
            env%env_metal_fraction(k) = env%env_metal_fraction(k+1)
            env%env_convective_flag(k) = env%env_convective_flag(k+1)
            last_kept_log_mass = log_mass
            num_kept = num_kept + 1
         end if
      end do
      env%num_env_points = num_kept
end subroutine shift_envelope_onto_shell

! ---------------------------------------------------------------
! Appends the envelope points inside the new fitting mass
! log_fit_mass_fraction = log(m_fit/M) to the model as zones nz+1,
! nz+2, ..., and the first point outside it as a last zone interpolated
! linearly in log m onto the fitting mass (between_points), between the
! zone or point below it and that point; advances nz. The model is its
! arrays (log m in g, log P, log T, log r, log rho, L/Lsun, the
! convective flag, the composition xa(species, zone)), of which zones
! 1..nz are set on entry. Before the call the caller has shifted env so
! that its point 1 is the first point outside the last interior shell.
! - A point exactly at the fitting mass is taken as a zone with
!   accept_point_at_fit_mass, and interpolated (onto itself) without.
! - Each new zone keeps the luminosity of the last interior shell; it
!   takes X and Z from the envelope, He3 and the species up to
!   last_copied_species from the last interior shell, and
!   Y = 1 - X - Z - X(He3).
! - The interpolated zone is convective if either of its end points is.
! - When every envelope point, the photosphere included, lies inside the
!   fitting mass, all become zones and the model ends at the
!   photosphere.
! ierr = 1 when the interpolation's interval in log m is below
! min_interpolation_interval; nz is then unchanged (the zones above it
! are partly written).
pure subroutine append_envelope_shells(env, log_fit_mass_fraction, log_total_mass, &
     accept_point_at_fit_mass, last_copied_species, nz, log_mass, logP, logT, logR, &
     logRho, luminosity_lsun, convective_flag, xa, ierr)
      type(envelope_structure_state), intent(in) :: env
      real(dp), intent(in) :: log_fit_mass_fraction, log_total_mass
      logical, intent(in) :: accept_point_at_fit_mass
      integer, intent(in) :: last_copied_species
      integer, intent(inout) :: nz
      real(dp), intent(inout) :: log_mass(:), logP(:), logT(:), logR(:), logRho(:), &
           luminosity_lsun(:)
      logical, intent(inout) :: convective_flag(:)
      real(dp), intent(inout) :: xa(:,:)
      integer, intent(out) :: ierr

      integer :: zone, k, species
      logical :: inside_fit_mass
! the end points of the interpolation: 0 the zone or point below the
! fitting mass, 2 the first point outside it
      real(dp) :: log_mass0, log_mass2, fit_log_mass, fraction
      real(dp) :: logRho0, logP0, logR0, logT0, x0, z0
      logical :: convective0

      ierr = 0
      do k = 1, env%num_env_points
         zone = nz + k
! The luminosity is constant through the envelope.
         luminosity_lsun(zone) = luminosity_lsun(nz)
         do species = i_he3, last_copied_species
            xa(species,zone) = xa(species,nz)
         end do
         if (accept_point_at_fit_mass) then
            inside_fit_mass = env%env_log10_mass(k) <= log_fit_mass_fraction
         else
            inside_fit_mass = env%env_log10_mass(k) < log_fit_mass_fraction
         end if
         if (inside_fit_mass) then
! A point inside the new fitting mass becomes a zone as it is.
            logRho(zone) = env%env_log10_density(k)
            logP(zone) = env%env_log10_pressure(k)
            logR(zone) = env%env_log10_radius(k)
            log_mass(zone) = env%env_log10_mass(k) + log_total_mass
            logT(zone) = env%env_log10_temperature(k)
            xa(i_h1,zone) = env%env_hydrogen_fraction(k)
            xa(i_metals,zone) = env%env_metal_fraction(k)
            convective_flag(zone) = env%env_convective_flag(k)
         else
! The first point outside the fitting mass: the last zone, interpolated
! onto the fitting mass between the last interior shell (for the first
! envelope point) or the envelope point before, and this point.
            if (k == 1) then
               log_mass0 = log_mass(nz)
               logRho0 = logRho(nz)
               logP0 = logP(nz)
               logR0 = logR(nz)
               logT0 = logT(nz)
               x0 = xa(i_h1,nz)
               z0 = xa(i_metals,nz)
               convective0 = convective_flag(nz)
            else
               log_mass0 = env%env_log10_mass(k-1) + log_total_mass
               logRho0 = env%env_log10_density(k-1)
               logP0 = env%env_log10_pressure(k-1)
               logR0 = env%env_log10_radius(k-1)
               logT0 = env%env_log10_temperature(k-1)
               x0 = env%env_hydrogen_fraction(k-1)
               z0 = env%env_metal_fraction(k-1)
               convective0 = env%env_convective_flag(k-1)
            end if
            fit_log_mass = log_total_mass + log_fit_mass_fraction
            log_mass2 = env%env_log10_mass(k) + log_total_mass
            if (log_mass2-log_mass0 < min_interpolation_interval) then
               ierr = 1
               return
            end if
            fraction = (fit_log_mass-log_mass0)/(log_mass2-log_mass0)
            logRho(zone) = between_points(logRho0, env%env_log10_density(k), fraction)
            logP(zone) = between_points(logP0, env%env_log10_pressure(k), fraction)
            logR(zone) = between_points(logR0, env%env_log10_radius(k), fraction)
            log_mass(zone) = fit_log_mass
            logT(zone) = between_points(logT0, env%env_log10_temperature(k), fraction)
            xa(i_h1,zone) = between_points(x0, env%env_hydrogen_fraction(k), fraction)
            xa(i_metals,zone) = between_points(z0, env%env_metal_fraction(k), fraction)
            convective_flag(zone) = env%env_convective_flag(k) .or. convective0
         end if
         xa(i_he4,zone) = 1.0d0-xa(i_h1,zone)-xa(i_metals,zone)-xa(i_he3,zone)
         if (.not. inside_fit_mass) then
            nz = zone
            return
         end if
      end do
! Every envelope point lay inside the fitting mass.
      nz = nz + env%num_env_points
end subroutine append_envelope_shells

! ---------------------------------------------------------------
! The masses of the shells from first_new on, after zones were appended
! to log_mass (all of the size of the model): the enclosed mass
! m = 10^log m (g), and each zone's mass
!   dm(k) = (m(k+1) - m(k-1))/2,
! with the outermost zone holding the rest of the star,
!   dm(nz) = M - (m(nz) + m(nz-1))/2,  M = 10^log_total_mass.
! first_new is the last zone the model had before (its dm changes, as
! its outer neighbour is new), and at least 2.
pure subroutine refit_shell_masses(first_new, log_mass, log_total_mass, m, dm)
      use math_lib
      integer, intent(in) :: first_new
      real(dp), intent(in) :: log_mass(:), log_total_mass
      real(dp), intent(inout) :: m(:), dm(:)
      integer :: k, nz

      nz = size(log_mass)
      do k = first_new, nz
         m(k) = exp10(log_mass(k))
      end do
      do k = first_new, nz-1
         dm(k) = 0.5d0*(m(k+1)-m(k-1))
      end do
      dm(nz) = exp10(log_total_mass) - 0.5d0*(m(nz)+m(nz-1))
end subroutine refit_shell_masses

! ---------------------------------------------------------------
! The first guess at the rotation of the zones above last_old (the
! model's last zone before the refit), from that zone by the
! convection-zone law of exponent cz_rotation_law_exponent (legacy
! WALPCZ):
!   exponent <= -2   constant specific angular momentum j;
!   exponent >= 0    solid body, omega constant;
!   otherwise        omega = omega_old (r/r_old)^exponent, the law
!                    enforce_rotation_profile (omega_profile) imposes.
! Each new zone's specific moment of inertia is the thin shell's,
! i = (2/3) r^2, and j = omega i. The arrays have the size of the model.
pure subroutine guess_new_shell_rotation(last_old, cz_rotation_law_exponent, logR, &
     omega, i_rot, j_rot)
      use math_lib
      use physical_constants, only: cc23
      integer, intent(in) :: last_old
      real(dp), intent(in) :: cz_rotation_law_exponent
      real(dp), intent(in) :: logR(:)
      real(dp), intent(inout) :: omega(:), i_rot(:), j_rot(:)
      real(dp) :: omega_at_unit_radius
      integer :: k

      if (cz_rotation_law_exponent <= -2.0d0) then
         do k = last_old+1, size(logR)
            j_rot(k) = j_rot(last_old)
            i_rot(k) = cc23*exp10((2.0d0*logR(k)))
            omega(k) = j_rot(k)/i_rot(k)
         end do
      else if (cz_rotation_law_exponent >= 0.0d0) then
         do k = last_old+1, size(logR)
            omega(k) = omega(last_old)
            i_rot(k) = cc23*exp10((2.0d0*logR(k)))
            j_rot(k) = omega(last_old)*i_rot(k)
         end do
      else
         omega_at_unit_radius = omega(last_old)/exp10((logR(last_old)*cz_rotation_law_exponent))
         do k = last_old+1, size(logR)
            omega(k) = omega_at_unit_radius*exp10((logR(k)*cz_rotation_law_exponent))
            i_rot(k) = cc23*exp10((2.0d0*logR(k)))
            j_rot(k) = omega(k)*i_rot(k)
         end do
      end if
end subroutine guess_new_shell_rotation

end module envelope_refit_kernel

! Sources
! - The F77's GETNEWENV (getnewenv.f: the shift, the append loop, the
!   masses and the rotation guess) and STARIN's copy of the append loop;
!   the loop was shared as append_envelope_points in 2026 (readability
!   wave 2), and became these kernels in round 3 (2026-11).
! - GETNEWENV left the convective flag unshifted, one point off (fixed
!   2026-10), and interpolated X and Z of a zone between the last
!   interior shell and the first envelope point with the difference
!   reversed, X_shell + f (X_shell - X_env), away from the envelope's
!   value; the term vanishes when the envelope has the shell's X and Z,
!   as both drivers integrate it. between_points writes it as the other
!   variables' since round 3.

!----------------------------------------------------------------------
! test_envelope_refit
!----------------------------------------------------------------------
! The kernels that turn the points of a fresh envelope integration into
! interior shells (core/envelope/envelope_refit_kernel.f90), on a small
! synthetic model: six zones and an envelope of six points whose log P,
! log T, log r and log rho are linear in log m, so that an interpolation
! onto a mass has an exact answer. No deck, no tables. Sections:
!   append    append_envelope_shells: the points inside the fitting mass
!             become zones as they are (exactly); the last zone lies at
!             the fitting mass and on the linear profiles (continuity,
!             1e-11), in both branches (between two envelope points;
!             between the last interior shell and the first point); the zone
!             is convective if either end point is; a point exactly at
!             the fitting mass is a zone when accepted and the
!             interpolated zone when not; every point inside the fitting
!             mass becomes a zone; Y = 1 - X - Z - X(He3), the species
!             up to last_copied_species and the luminosity are the last
!             interior shell's; a degenerate interval returns ierr 1
!             with nz unchanged.
!   shift     shift_envelope_onto_shell: the kept points carry the
!             offsets that join point 1 to the shell; log(m/M), X, Z and
!             the convective flag move unshifted.
!   masses    refit_shell_masses: m = 10^log m; dm(k) = (m(k+1) -
!             m(k-1))/2 and the outermost zone holds the rest, so that
!             the masses sum to the star's (1e-14); the zones below
!             first_new are untouched.
!   rotation  guess_new_shell_rotation: constant j for an exponent <= -2,
!             solid body for >= 0, omega ~ r^exponent between; i = (2/3)
!             r^2 and j = omega i.
! Each check prints its value (1pe9.1) and ok or FAIL; any FAIL stops
! the program with status 1.
program test_envelope_refit
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use envelope_points, only: envelope_structure_state
      use composition_layout, only: i_h1, i_he3, i_he4, i_metals, n_species_basic, &
           n_species_extended
      use envelope_refit_kernel, only: append_envelope_shells, shift_envelope_onto_shell, &
           refit_shell_masses, guess_new_shell_rotation, between_points
      use physical_constants, only: cc23
      implicit none
      integer, parameter :: num_zones = 6, num_points = 6, max_zones = 20
! the model: log M (g), and the profiles' intercepts and slopes in log m
      real(dp), parameter :: log_total_mass = 33.0_dp
      real(dp), parameter :: slope(4) = [-40.0_dp, -9.0_dp, 0.3_dp, -30.0_dp]
      real(dp), parameter :: intercept(4) = [1330.0_dp, 301.0_dp, 0.9_dp, 985.0_dp]
! the zones' and the envelope points' log(m/M)
      real(dp), parameter :: zone_mass_fraction(num_zones) = &
           [-2.0_dp, -1.0_dp, -0.5_dp, -0.1_dp, -0.03_dp, -0.006_dp]
      real(dp), parameter :: point_mass_fraction(num_points) = &
           [-0.005_dp, -0.004_dp, -0.003_dp, -0.002_dp, -0.001_dp, 0.0_dp]
      real(dp), parameter :: shell_x = 0.70_dp, shell_z = 0.02_dp, shell_he3 = 1.0e-4_dp
      real(dp), parameter :: envelope_x = 0.72_dp, envelope_z = 0.015_dp
! the model's arrays, and the envelope
      real(dp) :: log_mass(max_zones), logP(max_zones), logT(max_zones), &
           logR(max_zones), logRho(max_zones), luminosity_lsun(max_zones), &
           xa(n_species_extended,max_zones)
      logical :: convective_flag(max_zones)
      type(envelope_structure_state), allocatable :: env
      integer :: nz, ierr, nfail

      nfail = 0
      allocate(env)
      write(*,'(a)') 'test_envelope_refit report'
      write(*,'(a)') '-- append'
      call check_append_between_points()
      call check_append_above_shell()
      call check_append_at_fit_mass()
      call check_append_all_inside()
      call check_append_degenerate()
      write(*,'(a)') '-- shift'
      call check_shift()
      write(*,'(a)') '-- masses'
      call check_masses()
      write(*,'(a)') '-- rotation'
      call check_rotation()
      if (nfail > 0) then
         write(*,'(a,i0,a)') 'test_envelope_refit: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_envelope_refit: done'

contains

! ---------------------------------------------------------------
! An ok/FAIL line.
subroutine report(name, value, pass)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      logical, intent(in) :: pass
      if (pass) then
         write(*,'(a,t68,1pe9.1,2x,a)') name, value, 'ok'
      else
         write(*,'(a,t68,1pe9.1,2x,a)') name, value, 'FAIL'
         nfail = nfail + 1
      end if
end subroutine report

! Profile j (1 log P, 2 log T, 3 log r, 4 log rho) at log(m/M).
pure real(dp) function profile(j, mass_fraction)
      integer, intent(in) :: j
      real(dp), intent(in) :: mass_fraction
      profile = intercept(j) + slope(j)*(mass_fraction + log_total_mass)
end function profile

! ---------------------------------------------------------------
! The model of num_zones zones, radiative, of the shell's composition,
! with tracers that tell the zones apart, and the envelope of
! num_points points from the centre outward, of the envelope's
! composition, convective from point 4 on.
subroutine set_model()
      integer :: k, species
      nz = num_zones
      log_mass = 0.0_dp
      logP = 0.0_dp
      logT = 0.0_dp
      logR = 0.0_dp
      logRho = 0.0_dp
      luminosity_lsun = 0.0_dp
      xa = 0.0_dp
      convective_flag = .false.
      do k = 1, num_zones
         log_mass(k) = zone_mass_fraction(k) + log_total_mass
         logP(k) = profile(1, zone_mass_fraction(k))
         logT(k) = profile(2, zone_mass_fraction(k))
         logR(k) = profile(3, zone_mass_fraction(k))
         logRho(k) = profile(4, zone_mass_fraction(k))
         luminosity_lsun(k) = 0.1_dp*k
         xa(i_h1,k) = shell_x
         xa(i_metals,k) = shell_z
         xa(i_he3,k) = shell_he3
         do species = i_he3+1, n_species_extended
            xa(species,k) = 1.0e-6_dp*species*k
         end do
         xa(i_he4,k) = 1.0_dp - shell_x - shell_z - shell_he3
      end do
      env%num_env_points = num_points
      do k = 1, num_points
         env%env_log10_mass(k) = point_mass_fraction(k)
         env%env_log10_pressure(k) = profile(1, point_mass_fraction(k))
         env%env_log10_temperature(k) = profile(2, point_mass_fraction(k))
         env%env_log10_radius(k) = profile(3, point_mass_fraction(k))
         env%env_log10_density(k) = profile(4, point_mass_fraction(k))
         env%env_hydrogen_fraction(k) = envelope_x
         env%env_metal_fraction(k) = envelope_z
         env%env_convective_flag(k) = k >= 4
      end do
end subroutine set_model

! The largest distance of zone's log P, log T, log r and log rho from
! the linear profiles at log(m/M).
real(dp) function off_profile(zone, mass_fraction)
      integer, intent(in) :: zone
      real(dp), intent(in) :: mass_fraction
      off_profile = max(abs(logP(zone) - profile(1, mass_fraction)), &
           abs(logT(zone) - profile(2, mass_fraction)), &
           abs(logR(zone) - profile(3, mass_fraction)), &
           abs(logRho(zone) - profile(4, mass_fraction)))
end function off_profile

! Whether zone is envelope point k as it is.
logical function is_point(zone, k)
      integer, intent(in) :: zone, k
      is_point = logP(zone) == env%env_log10_pressure(k) .and. &
           logT(zone) == env%env_log10_temperature(k) .and. &
           logR(zone) == env%env_log10_radius(k) .and. &
           logRho(zone) == env%env_log10_density(k) .and. &
           log_mass(zone) == env%env_log10_mass(k) + log_total_mass .and. &
           xa(i_h1,zone) == env%env_hydrogen_fraction(k) .and. &
           xa(i_metals,zone) == env%env_metal_fraction(k) .and. &
           (convective_flag(zone) .eqv. env%env_convective_flag(k))
end function is_point

! ---------------------------------------------------------------
! A fitting mass between envelope points 3 and 4.
subroutine check_append_between_points()
      real(dp), parameter :: fit = -0.0025_dp
      real(dp) :: x_expected, off
      integer :: species
      logical :: copied

      call set_model()
      call append_envelope_shells(env, fit, log_total_mass, .false., n_species_extended, &
           nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, convective_flag, xa, ierr)
      call report('append: three points inside, then the interpolated zone', &
           real(nz, dp), ierr == 0 .and. nz == num_zones + 4)
      call report('append: the points inside become zones as they are', 0.0_dp, &
           is_point(num_zones+1, 1) .and. is_point(num_zones+2, 2) .and. &
           is_point(num_zones+3, 3))
      call report('append: the last zone lies at the fitting mass', &
           abs(log_mass(nz) - (log_total_mass + fit)), &
           log_mass(nz) == log_total_mass + fit)
      off = off_profile(nz, fit)
      call report('append: continuity at the fitting mass (between two points)', off, &
           off <= 1.0e-11_dp)
! X between envelope points of one X is that X; the zone is convective
! because point 4 is.
      call report('append: convective if either end point is (points 3, 4)', 0.0_dp, &
           convective_flag(nz) .and. xa(i_h1,nz) == envelope_x .and. &
           xa(i_metals,nz) == envelope_z)
      x_expected = 1.0_dp - xa(i_h1,nz) - xa(i_metals,nz) - shell_he3
      copied = .true.
      do species = i_he3, n_species_extended
         copied = copied .and. all(xa(species,num_zones+1:nz) == xa(species,num_zones))
      end do
      call report('append: Y = 1 - X - Z - X(He3); the other species the shell''s', &
           abs(xa(i_he4,nz) - x_expected), xa(i_he4,nz) == x_expected .and. copied .and. &
           xa(i_he4,num_zones+1) == 1.0_dp - envelope_x - envelope_z - shell_he3)
      call report('append: the luminosity is the last shell''s', &
           luminosity_lsun(nz), all(luminosity_lsun(num_zones+1:nz) == luminosity_lsun(num_zones)))
! With the basic species set, the extended slots of the new zones are
! not written.
      call set_model()
      call append_envelope_shells(env, fit, log_total_mass, .false., n_species_basic, &
           nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, convective_flag, xa, ierr)
      call report('append: species above last_copied_species are not written', 0.0_dp, &
           all(xa(n_species_basic+1:,num_zones+1:nz) == 0.0_dp) .and. &
           xa(n_species_basic,nz) == xa(n_species_basic,num_zones))
end subroutine check_append_between_points

! ---------------------------------------------------------------
! A fitting mass between the last interior shell and the first envelope
! point.
subroutine check_append_above_shell()
      real(dp), parameter :: fit = -0.00525_dp
      real(dp) :: off

      call set_model()
      call append_envelope_shells(env, fit, log_total_mass, .false., n_species_extended, &
           nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, convective_flag, xa, ierr)
      off = off_profile(nz, fit)
      call report('append: continuity at the fitting mass (shell to first point)', off, &
           ierr == 0 .and. nz == num_zones + 1 .and. off <= 1.0e-11_dp .and. &
           log_mass(nz) == log_total_mass + fit)
      call report('append: radiative between a radiative shell and point', 0.0_dp, &
           .not. convective_flag(nz))
      call set_model()
      convective_flag(num_zones) = .true.
      call append_envelope_shells(env, fit, log_total_mass, .false., n_species_extended, &
           nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, convective_flag, xa, ierr)
      call report('append: convective if the shell is', 0.0_dp, convective_flag(nz))
end subroutine check_append_above_shell

! ---------------------------------------------------------------
! A fitting mass exactly at envelope point 3: interpolated (onto
! itself) unless the point is accepted as a zone.
subroutine check_append_at_fit_mass()
      real(dp) :: off

      call set_model()
      call append_envelope_shells(env, point_mass_fraction(3), log_total_mass, .false., &
           n_species_extended, nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, &
           convective_flag, xa, ierr)
      off = off_profile(nz, point_mass_fraction(3))
      call report('append: a point at the fitting mass, not accepted: interpolated', off, &
           ierr == 0 .and. nz == num_zones + 3 .and. off <= 1.0e-11_dp)
      call set_model()
      call append_envelope_shells(env, point_mass_fraction(3), log_total_mass, .true., &
           n_species_extended, nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, &
           convective_flag, xa, ierr)
      call report('append: a point at the fitting mass, accepted: a zone as it is', 0.0_dp, &
           ierr == 0 .and. is_point(num_zones+3, 3))
end subroutine check_append_at_fit_mass

! ---------------------------------------------------------------
! A fitting mass above every point: all become zones.
subroutine check_append_all_inside()
      call set_model()
      call append_envelope_shells(env, 1.0e-3_dp, log_total_mass, .false., &
           n_species_extended, nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, &
           convective_flag, xa, ierr)
      call report('append: every point inside the fitting mass becomes a zone', &
           real(nz, dp), ierr == 0 .and. nz == num_zones + num_points .and. &
           is_point(nz, num_points) .and. is_point(num_zones+1, 1))
end subroutine check_append_all_inside

! ---------------------------------------------------------------
! The first envelope point at the last shell's mass, and a fitting mass
! below it: an interval of zero width.
subroutine check_append_degenerate()
      call set_model()
      env%env_log10_mass(1) = zone_mass_fraction(num_zones)
      call append_envelope_shells(env, -0.0065_dp, log_total_mass, .false., &
           n_species_extended, nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, &
           convective_flag, xa, ierr)
      call report('append: a degenerate interval returns ierr 1, nz unchanged', &
           real(ierr, dp), ierr == 1 .and. nz == num_zones)
end subroutine check_append_degenerate

! ---------------------------------------------------------------
! The shift onto a shell that is offset from envelope point 1, with
! point 4 moved to within the mass step of point 3.
subroutine check_shift()
      real(dp), parameter :: offset(4) = [0.01_dp, -0.02_dp, 0.003_dp, 0.04_dp]
      type(envelope_structure_state), allocatable :: before
      real(dp) :: worst
      integer :: k
      logical :: moved

      call set_model()
      allocate(before)
      before = env
      call shift_envelope_onto_shell(env, before%env_log10_pressure(1) + offset(1), &
           before%env_log10_density(1) + offset(4), &
           before%env_log10_temperature(1) + offset(2), &
           before%env_log10_radius(1) + offset(3), &
           before%env_log10_mass(1) + log_total_mass, log_total_mass)
      worst = 0.0_dp
      moved = env%num_env_points == num_points - 1
      do k = 1, num_points - 1
         worst = max(worst, abs(env%env_log10_pressure(k) - before%env_log10_pressure(k+1) - offset(1)), &
              abs(env%env_log10_temperature(k) - before%env_log10_temperature(k+1) - offset(2)), &
              abs(env%env_log10_radius(k) - before%env_log10_radius(k+1) - offset(3)), &
              abs(env%env_log10_density(k) - before%env_log10_density(k+1) - offset(4)))
         moved = moved .and. env%env_log10_mass(k) == before%env_log10_mass(k+1) .and. &
              env%env_hydrogen_fraction(k) == before%env_hydrogen_fraction(k+1) .and. &
              env%env_metal_fraction(k) == before%env_metal_fraction(k+1) .and. &
              (env%env_convective_flag(k) .eqv. before%env_convective_flag(k+1))
      end do
      call report('shift: the kept points carry the offsets of point 1 to the shell', worst, &
           env%num_env_points == num_points - 1 .and. worst <= 1.0e-13_dp)
      call report('shift: log(m/M), X, Z and the convective flag move as they are', 0.0_dp, moved)
end subroutine check_shift

! ---------------------------------------------------------------
! The masses after an append: the model's own weights, with the centre's
! dm(1) = (m(1) + m(2))/2, telescope to M.
subroutine check_masses()
      real(dp) :: m(max_zones), dm(max_zones), total, worst
      integer :: k

      call set_model()
      m = -1.0_dp
      dm = -1.0_dp
      do k = 1, num_zones
         m(k) = 10.0_dp**log_mass(k)
      end do
      dm(1) = 0.5_dp*(m(1) + m(2))
      do k = 2, num_zones - 1
         dm(k) = 0.5_dp*(m(k+1) - m(k-1))
      end do
      dm(num_zones) = 12345.0_dp
      call append_envelope_shells(env, -0.0025_dp, log_total_mass, .false., &
           n_species_extended, nz, log_mass, logP, logT, logR, logRho, luminosity_lsun, &
           convective_flag, xa, ierr)
      call refit_shell_masses(num_zones, log_mass(1:nz), log_total_mass, m(1:nz), dm(1:nz))
      worst = 0.0_dp
      do k = num_zones, nz
         worst = max(worst, abs(m(k)/10.0_dp**log_mass(k) - 1.0_dp))
      end do
      call report('masses: m = 10^log m for the new zones', worst, worst <= 1.0e-14_dp)
      worst = 0.0_dp
      do k = num_zones, nz - 1
         worst = max(worst, abs(dm(k)/(0.5_dp*(m(k+1) - m(k-1))) - 1.0_dp))
      end do
      call report('masses: dm(k) = (m(k+1) - m(k-1))/2', worst, worst <= 1.0e-14_dp)
      total = sum(dm(1:nz))
      call report('masses: the zones'' masses sum to the star''s', &
           abs(total/10.0_dp**log_total_mass - 1.0_dp), &
           abs(total/10.0_dp**log_total_mass - 1.0_dp) <= 1.0e-14_dp)
      call report('masses: the zones below first_new are untouched', 0.0_dp, &
           dm(1) == 0.5_dp*(10.0_dp**log_mass(1) + 10.0_dp**log_mass(2)) .and. &
           all(m(nz+1:) == -1.0_dp) .and. all(dm(nz+1:) == -1.0_dp) .and. &
           dm(num_zones-1) == 0.5_dp*(m(num_zones) - m(num_zones-2)))
end subroutine check_masses

! ---------------------------------------------------------------
! The three laws of the new zones' rotation, on zones 5 and 6 above
! zone 4.
subroutine check_rotation()
      integer, parameter :: last_old = 4
      real(dp), parameter :: omega_old = 3.0e-6_dp, j_old = 2.0e16_dp
      real(dp) :: omega(num_zones), i_rot(num_zones), j_rot(num_zones), worst, expected
      integer :: k

      call set_model()
      omega = omega_old
      i_rot = 1.0_dp
      j_rot = j_old
      call guess_new_shell_rotation(last_old, -2.0_dp, logR(1:num_zones), omega, i_rot, j_rot)
      worst = 0.0_dp
      do k = last_old + 1, num_zones
         worst = max(worst, abs(omega(k)*cc23*10.0_dp**(2.0_dp*logR(k))/j_old - 1.0_dp))
      end do
      call report('rotation: exponent -2: constant j, omega = j/i', worst, &
           all(j_rot(last_old+1:) == j_old) .and. worst <= 1.0e-13_dp)

      omega = omega_old
      j_rot = j_old
      call guess_new_shell_rotation(last_old, 0.0_dp, logR(1:num_zones), omega, i_rot, j_rot)
      worst = 0.0_dp
      do k = last_old + 1, num_zones
         worst = max(worst, abs(j_rot(k)/(omega_old*cc23*10.0_dp**(2.0_dp*logR(k))) - 1.0_dp))
      end do
      call report('rotation: exponent 0: solid body, j = omega i', worst, &
           all(omega(last_old+1:) == omega_old) .and. worst <= 1.0e-13_dp)

      omega = omega_old
      j_rot = j_old
      call guess_new_shell_rotation(last_old, -1.0_dp, logR(1:num_zones), omega, i_rot, j_rot)
      worst = 0.0_dp
      do k = last_old + 1, num_zones
         expected = omega_old*10.0_dp**(-(logR(k) - logR(last_old)))
         worst = max(worst, abs(omega(k)/expected - 1.0_dp), &
              abs(j_rot(k)/(omega(k)*i_rot(k)) - 1.0_dp))
      end do
      call report('rotation: exponent -1: omega ~ 1/r, j = omega i', worst, &
           worst <= 1.0e-13_dp)
      worst = 0.0_dp
      do k = last_old + 1, num_zones
         worst = max(worst, abs(i_rot(k)/(cc23*10.0_dp**(2.0_dp*logR(k))) - 1.0_dp))
      end do
      call report('rotation: i = (2/3) r^2; the old zones untouched', worst, &
           worst <= 1.0e-13_dp .and. all(omega(:last_old) == omega_old) .and. &
           all(j_rot(:last_old) == j_old) .and. all(i_rot(:last_old) == 1.0_dp))
! between_points at its ends.
      call report('between_points: v0 at 0, v2 at 1, the mean at 1/2', &
           abs(between_points(2.0_dp, 6.0_dp, 0.5_dp) - 4.0_dp), &
           between_points(2.0_dp, 6.0_dp, 0.0_dp) == 2.0_dp .and. &
           between_points(2.0_dp, 6.0_dp, 1.0_dp) == 6.0_dp .and. &
           between_points(2.0_dp, 6.0_dp, 0.5_dp) == 4.0_dp)
end subroutine check_rotation

end program test_envelope_refit

!----------------------------------------------------------------------
! test_stitched
!----------------------------------------------------------------------
! ok/FAIL checks of the stitched model (core/diagnostics/), each
! confirmed to fail under a mutation (test_stitched.mutations.toml).
!   The kernels (stitched_kernels.f90) on synthetic arrays:
!     - accumulate_atmosphere_heights: each point's height is the sum of
!       the steps below it, the first step entering none;
!     - map_extended_grid: the envelope's repeat of the fitting point and
!       the atmosphere's repeat of the photosphere are skipped; an
!       envelope point just below the interior's top is skipped and one
!       above it kept; the map stops at the size of its arrays;
!     - thermal_brunt_n2 against g^2 (rho/P) delta (nabla_ad - nabla);
!     - difference_brunt_n2 on ln P and ln rho linear in r, for which the
!       centred difference is exact: N^2 = g (1/H_rho - 1/(Gamma_1 H_P));
!       a negative thermal value is kept; a thermal value of exactly zero
!       takes the difference; the end points copy their neighbours;
!     - acoustic_radius of a constant sound speed, (r_n - r_1)/c, and an
!       interval without a sound speed left out;
!     - buoyancy_radius of N = a r, a times the span the sum covers, and
!       only the innermost cavity.
!   The star's stitched model (stitched_model.f90), on the starting model
!   of the deck given as the argument, set up as run_yrec does:
!     - build_stitched_model: the three regions in order, the radius
!       increasing across both joins, and each region's points those of
!       the interior, the envelope and the atmosphere;
!     - the profile's brunt_N2 is the pulse points' N^2;
!     - a profile column without a value in a region is an error
!       (fill_stitched_arrays, ierr = 1);
!     - when the envelope integration fails (an atmosphere step size that
!       needs more steps than the integration allows), the grid is the
!       interior alone, and the failure is returned (bug_sweep M327).
! Four of the kernels' checks hold values that an open defect's fix will
! change, and say so: bug_sweep M359 (the envelope point above the
! interior's top), M361 (the zero thermal value), M363 (the innermost
! point's N^2) and M362 (the two integrals' inner ends).
!
! Usage: test_stitched DECK.inlist (test_stitched.py passes
! core/test/test_structure_tams.inlist), run from a directory with an
! empty output/, with YREC_OUTPUT=output and YREC_INPUT and YREC_START
! set. Prints one line per check and exits with status 1 when one fails.
program test_stitched
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use stitched_kernels, only: accumulate_atmosphere_heights, map_extended_grid, &
           thermal_brunt_n2, difference_brunt_n2, acoustic_radius, buoyancy_radius, &
           region_interior, region_envelope, region_atmosphere
      use stitched_model, only: build_stitched_model, fill_stitched_arrays, ip_logR, &
           ip_brunt_N2, ipul_N2
      use controls_reader, only: read_controls
      use run_files, only: open_run_files
      use run_setup, only: set_up_constants_and_tables
      use net_lib, only: set_net_settings
      use net_settings_builder, only: net_settings_from_star
      use convection_lib, only: set_convection_settings
      use convection_settings_builder, only: convection_settings_from_star
      use rotation_configuration, only: set_rotation_settings
      use rotation_settings_builder, only: rotation_settings_from_star
      use starting_model, only: read_starting_model
      use star_info_lib, only: star_info, control_nml_override, physics_nml_override
      implicit none
! the test's star; static, as it is too large for the stack
      type(star_info), save :: star
      character(len=512) :: deck
      integer :: num_failed = 0

      if (command_argument_count() /= 1) then
         write(*,'(a)') 'usage: test_stitched DECK.inlist'
         stop 2
      end if
      call get_command_argument(1, deck)

      call kernel_checks()
      call star_checks()

      if (num_failed > 0) then
         write(*,'(a,i0,a)') 'test_stitched: ', num_failed, ' checks FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_stitched: done'

contains

! ---------------------------------------------------------------
! The kernels on synthetic arrays.
subroutine kernel_checks()
      integer, parameter :: n = 12
      real(dp), parameter :: tol = 1.0e-12_dp
      real(dp), parameter :: pressure_height = 2.0e9_dp, density_height = 3.0e9_dp, &
           gamma1_value = 1.6_dp, newton_g = 6.67e-8_dp, sound_speed_value = 2.5e7_dp, &
           n_over_r = 3.0e-13_dp
      integer :: region(12), index(12), small_region(8), small_index(8)
      integer :: n_ie, n_ext, j
      real(dp) :: heights(4)
      real(dp) :: radius(n), mass(n), pressure(n), density(n), gamma1(n), brunt_n2(n), &
           expected_n2(n), log_radius(n), sound_speed(n)

! The atmosphere's heights: four points from the outside in, with steps
! 9 (which enters no height), 1, 2 and 3.
      heights = -1.0_dp
      call accumulate_atmosphere_heights([9.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], heights)
      call stitched_check('atmosphere heights: the sum of the steps below each point', &
           all(heights == [6.0_dp, 5.0_dp, 3.0_dp, 0.0_dp]))

! The index map: five interior points whose top is at log r = 10; three
! envelope points, the first the fitting point again; four atmosphere
! points, the deepest (height 0) the photosphere again.
      region = 0
      index = 0
      call map_extended_grid(5, 10.0_dp, [10.0_dp, 10.01_dp, 10.02_dp], heights, &
           region, index, n_ie, n_ext)
      call stitched_check('index map: the two repeated points at the joins are skipped', &
           n_ie == 7 .and. n_ext == 10 .and. &
           all(region(1:10) == [(region_interior, j = 1, 5), (region_envelope, j = 1, 2), &
           (region_atmosphere, j = 1, 3)]) .and. &
           all(index(1:10) == [1, 2, 3, 4, 5, 2, 3, 3, 2, 1]))
      call map_extended_grid(5, 10.0_dp, [10.0_dp - 1.0e-9_dp, 10.01_dp, 10.02_dp], &
           heights, region, index, n_ie, n_ext)
      call stitched_check('index map: an envelope point below the interior''s top is skipped', &
           n_ie == 7 .and. n_ext == 10 .and. index(6) == 2)
! bug_sweep M359: the envelope's first point is the fitting point again
! whatever its radius; today it is kept when its radius is above the
! interior's. That fix changes the expected n_ie to 7 and index(6) to 2.
      call map_extended_grid(5, 10.0_dp, [10.0_dp + 1.0e-6_dp, 10.01_dp, 10.02_dp], &
           heights, region, index, n_ie, n_ext)
      call stitched_check('index map: an envelope first point above the interior''s top is kept (M359)', &
           n_ie == 8 .and. n_ext == 11 .and. index(6) == 1 .and. &
           region(6) == region_envelope)
      small_region = 0
      small_index = 0
      call map_extended_grid(5, 10.0_dp, [10.0_dp, 10.01_dp, 10.02_dp], heights, &
           small_region, small_index, n_ie, n_ext)
      call stitched_check('index map: it stops at the size of its arrays', &
           n_ie == 7 .and. n_ext == 8 .and. small_region(8) == region_atmosphere)

! N^2, the thermal form.
      call stitched_check('thermal N^2 = g^2 (rho/P) delta (nabla_ad - nabla)', &
           abs(thermal_brunt_n2(2.0_dp, 3.0_dp, 6.0_dp, 1.5_dp, 0.4_dp, 0.25_dp) - 0.45_dp) &
           < tol)

! N^2 by differences: ln P and ln rho linear in r, with scale heights
! pressure_height and density_height, so that
!   N^2 = (G m/r^2) (1/density_height - 1/(Gamma_1 pressure_height))
! at every point. The thermal values on entry: positive, but for points
! 5 and 6 (negative, a convection zone) and point 8 (exactly zero).
      do j = 1, n
         radius(j) = 1.0e9_dp*real(j, dp)**1.3_dp
         mass(j) = 1.0e31_dp*real(j, dp)**2
      end do
      pressure = 1.0e17_dp*exp(-radius/pressure_height)
      density = 1.0e2_dp*exp(-radius/density_height)
      gamma1 = gamma1_value
      expected_n2 = newton_g*mass/radius**2* &
           (1.0_dp/density_height - 1.0_dp/(gamma1_value*pressure_height))
      brunt_n2 = 1.0_dp
      brunt_n2(5:6) = -1.0_dp
      brunt_n2(8) = 0.0_dp
      call difference_brunt_n2(radius, mass, pressure, density, gamma1, newton_g, brunt_n2)
      call stitched_check('N^2: the centred difference where the thermal form is positive', &
           all(abs(brunt_n2(2:4)/expected_n2(2:4) - 1.0_dp) < 1.0e-9_dp) .and. &
           all(abs(brunt_n2(9:n-1)/expected_n2(9:n-1) - 1.0_dp) < 1.0e-9_dp) .and. &
           abs(brunt_n2(7)/expected_n2(7) - 1.0_dp) < 1.0e-9_dp)
      call stitched_check('N^2: a negative thermal value is kept', &
           all(brunt_n2(5:6) == -1.0_dp))
! bug_sweep M361: a point whose thermal N^2 is exactly zero (nabla =
! nabla_ad in an efficient convection zone) takes the difference today;
! that fix changes the expected value to 0.
      call stitched_check('N^2: a thermal value of exactly zero takes the difference (M361)', &
           abs(brunt_n2(8)/expected_n2(8) - 1.0_dp) < 1.0e-9_dp)
! bug_sweep M363: the innermost point's N^2 is its neighbour's today;
! that fix changes the expected value to N^2(2) (r_1/r_2)^2.
      call stitched_check('N^2: the end points copy their neighbours (M363)', &
           brunt_n2(1) == brunt_n2(2) .and. brunt_n2(n) == brunt_n2(n-1))

! The acoustic radius of a constant sound speed: the trapezoid is exact,
! and starts at the first point. (bug_sweep M362: that fix adds
! r_1/c_1, and the expected value becomes r_n/c.)
      log_radius = log10(radius)
      sound_speed = sound_speed_value
      call stitched_check('acoustic radius of a constant sound speed: (r_n - r_1)/c (M362)', &
           abs(acoustic_radius(log_radius, sound_speed)/ &
           ((radius(n) - radius(1))/sound_speed_value) - 1.0_dp) < tol)
      sound_speed(4) = 0.0_dp
      call stitched_check('acoustic radius: an interval without a sound speed is left out', &
           abs(acoustic_radius(log_radius, sound_speed)/ &
           ((radius(n) - radius(5) + radius(3) - radius(1))/sound_speed_value) - 1.0_dp) < tol)

! The buoyancy radius of N = n_over_r r: every point of the cavity adds
! n_over_r times half the span to its neighbours, from point 2 to point
! n - 1. (bug_sweep M362: that fix adds the part inside (r_1 + r_2)/2,
! and the expected span starts at 0.)
      brunt_n2 = (n_over_r*radius)**2
      call stitched_check('buoyancy radius of N = a r: a times the span the sum covers (M362)', &
           abs(buoyancy_radius(log_radius, brunt_n2)/ &
           (n_over_r*0.5_dp*(radius(n) + radius(n-1) - radius(2) - radius(1))) - 1.0_dp) < tol)
! A convection zone at points 6 and 7 ends the cavity: points 2 to 5 count.
      brunt_n2(6:7) = -1.0e-8_dp
      call stitched_check('buoyancy radius: the innermost cavity only', &
           abs(buoyancy_radius(log_radius, brunt_n2)/ &
           (n_over_r*0.5_dp*(radius(6) + radius(5) - radius(2) - radius(1))) - 1.0_dp) < tol)
end subroutine kernel_checks

! ---------------------------------------------------------------
! The stitched model of the deck's starting model.
subroutine star_checks()
      integer :: ierr, j, n_ext, n_ie, nz
      logical :: regions_in_order, radius_increases, points_match

      call load_deck_model(deck)
! The fitting point's mass fraction, which the envelope is integrated
! down to: the step sets it before every solve (refresh_surface_bc), and
! a model just read has none.
      star%log_fit_point_mass_fraction = star%log_mass(star%nz) - star%log_total_mass
      call build_stitched_model(star, ierr)
      if (ierr /= 0) call die('build_stitched_model')
      n_ext = star%stitched%n_ext
      n_ie = star%stitched%n_ie
      nz = star%nz
      regions_in_order = n_ext > n_ie .and. n_ie > nz .and. &
           all(star%stitched%region(1:nz) == region_interior) .and. &
           all(star%stitched%region(nz+1:n_ie) == region_envelope) .and. &
           all(star%stitched%region(n_ie+1:n_ext) == region_atmosphere)
      call stitched_check('stitched model: interior, envelope and atmosphere, in order', &
           regions_in_order)
      radius_increases = all(star%stitched%prof(ip_logR,2:n_ext) > &
           star%stitched%prof(ip_logR,1:n_ext-1))
      call stitched_check('stitched model: the radius increases across both joins', &
           radius_increases)
      points_match = all(star%stitched%prof(ip_logR,1:nz) == star%logR(1:nz))
      do j = nz + 1, n_ie
         points_match = points_match .and. star%stitched%prof(ip_logR,j) == &
              star%env_points%env_log10_radius(star%stitched%index(j))
      end do
      call stitched_check('stitched model: each point is its region''s point', &
           points_match .and. star%stitched%index(nz) == nz .and. &
           star%stitched%index(n_ext) == 1)
      call stitched_check('stitched model: the profile''s brunt_N2 is the pulse points'' N^2', &
           all(star%stitched%prof(ip_brunt_N2,1:n_ext) == star%stitched%pulse(ipul_N2,1:n_ext)))

! A region without a case for a column: an unknown region code at point 1.
      star%stitched%region(1) = 0
      call fill_stitched_arrays(star, ierr)
      call stitched_check('a profile column without a value in a region: ierr = 1', ierr == 1)

! The envelope integration fails: an atmosphere step so small that the
! photosphere is not reached in the steps the integration allows.
      star%physics%atm_step_size = 1.0e-6_dp
      call build_stitched_model(star, ierr)
      call stitched_check('envelope integration fails: the grid is the interior alone', &
           star%stitched%n_ext == nz .and. star%stitched%n_ie == nz .and. &
           all(star%stitched%region(1:nz) == region_interior))
      call stitched_check('envelope integration fails: build_stitched_model returns the error', &
           ierr == 1)
end subroutine star_checks

! ---------------------------------------------------------------
! The deck's controls, tables, network and convection settings, and its
! first card's starting model, as run_yrec sets them up.
subroutine load_deck_model(path)
      character(len=*), intent(in) :: path
      integer :: ierr

      ierr = 0
      control_nml_override = path
      physics_nml_override = path
      call read_controls(star, ierr)
      if (ierr /= 0) call die('read_controls')
      call open_run_files(star, ierr)
      if (ierr /= 0) call die('open_run_files')
      call set_up_constants_and_tables(star, ierr)
      if (ierr /= 0) call die('set_up_constants_and_tables')
      call set_net_settings(net_settings_from_star(star), ierr)
      if (ierr /= 0) call die('set_net_settings')
      star%current_card = 1
      star%current_card_envelope_x = star%card_initial_x(1)
      star%current_card_envelope_z = star%card_initial_z(1)
      star%mixing_length_alpha = star%card_mixing_length_alpha(1)
      call set_convection_settings(convection_settings_from_star(star))
      call set_rotation_settings(rotation_settings_from_star(star))
      star%fit_point_move_requested = star%card_moves_fit_point(1)
      star%requested_log_fit_point_mass_fraction = star%card_log_fit_point_mass_fraction(1)
      star%total_angular_momentum = 0.0d0
      star%total_rotational_ke = 0.0d0
      call read_starting_model(star, ierr)
      if (ierr /= 0) call die('read_starting_model')
end subroutine load_deck_model

subroutine die(what)
      character(len=*), intent(in) :: what
      write(*,'(a)') 'test_stitched: FAIL ('//what//')'
      stop 1
end subroutine die

! ---------------------------------------------------------------
! One check: its label and ok, or FAIL.
subroutine stitched_check(label, condition)
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
end subroutine stitched_check

end program test_stitched

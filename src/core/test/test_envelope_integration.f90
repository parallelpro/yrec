!----------------------------------------------------------------------
! test_envelope_integration
!----------------------------------------------------------------------
! The envelope integration kernel (core/envelope/
! envelope_integration_kernel.f90) and its integrands
! (envelope_integrands.f90) on a solar surface: 1 Msun, 1 Lsun,
! Teff = 5777 K, the Eddington gray atmosphere, the fitting point at
! log(m/M) = -1e-4, with the microphysics and the integration controls
! of a deck (its controls and tables, as run_yrec sets them up; no
! model is read). Sections:
!   integrands  at a point of the atmosphere, d log P/d log tau equals
!               g tau/(kappa P) with kappa and P the point's (relative
!               1e-12), and T is the Eddington T(tau),
!               T^4 = (3/4) Teff^4 (tau + 2/3) (1e-12); at a point of the
!               envelope, d log m/d log P = -4 pi r^4 P/(G m^2) and
!               d log r/d log P = -r P/(G m rho), written here in linear
!               form (1e-12).
!   solve       the integration returns ierr = 0 and its events in order
!               (the atmosphere done, the vertex, the envelope done); it
!               ends at the photosphere's tau = 2/3 (log tau within the
!               tolerance fraction), where T = Teff (relative 1e-6), and
!               at the fitting mass (log(m/M) within the tolerance
!               fraction of it); the photospheric pressure is the
!               atmosphere's last point's.
!   depths      the atmosphere points' depth steps (atmo_delta_depth, the
!               stitched model's atmosphere radii): the first is measured
!               from tau = 0 with rho kappa = 0 there, dz = 2 tau/(rho
!               kappa) (M322: it read locals nothing had set); each later
!               one is the step in tau over the mean rho kappa of its two
!               points (relative 1e-6, with tau from the Eddington T(tau)).
!   reference   each Bulirsch-Stoer step of the atmosphere and of the
!               envelope against a fourth-order Runge-Kutta integration
!               of the same equations between the kernel's points (64
!               steps): the largest difference, relative to the step's
!               scale (as bsstep measures it), within 1000 times the
!               integration's error tolerance (the integrands interpolate
!               tables, so bsstep's estimate is no bound; the largest
!               seen is 55 times). (The atmosphere's points
!               carry log P and log T; tau follows from the Eddington
!               T(tau).) The vertex is the linear interpolation, in
!               log(m/M), of the last step's ends to the fitting mass.
!   thin        a fitting mass above max_log_fit_point_mass_fraction: no
!               envelope is integrated, and the vertex is the surface
!               point itself, exactly.
!   retry       a starting density so high that the starting optical
!               depth lies past the photosphere: an envint_start_retry
!               event (its log tau above the photosphere's), then a
!               normal integration. A surface of log g = -5, whose start
!               lies past the photosphere at any density: ierr = 1 and
!               an envint_start_exhausted event after max_start_tries
!               starting densities, the last max_start_tries - 1 decades
!               below the first.
!   errors      the error contracts: the atmosphere or the envelope
!               that cannot finish in max_steps steps (steps of 1e-6)
!               returns ierr = 1 and its maxstp event, last; the
!               photospheric pressure stays no_photosphere_log10_pressure
!               when the atmosphere failed. Bulirsch-Stoer limits out of range (nuse = 1)
!               return numerics_bad_input and an envint_bsstep_limits
!               event with imax and nuse. A step that underflows (an
!               error tolerance of zero, with the IEEE divide-by-zero and
!               invalid traps of a devel build held off for that one
!               integration) returns envelope_termination
!               and an envint_bsstep_failed event with bsstep's status.
!               The integrands failing inside a Bulirsch-Stoer step (a
!               first envelope step of 40 in log P toward a fitting mass
!               deep in the star: its midpoint substeps leave the
!               microphysics' tables) return their own positive code, an
!               error, not envelope_termination, and an
!               envint_bsstep_failed event with bsstep's
!               numerics_derivatives_failed and that code.
!   gray photosphere  envelope_integration's gray_photosphere (the
!               accretion surface after a failed Allard lookup, O19) at
!               the test's Teff and g: ierr = 0; its log P is the
!               kernel's gray photospheric pressure with the star's
!               settings (within 1e-10 dex), and T = Teff (relative
!               1e-6); the star's atmosphere point, atmosphere points,
!               photospheric pressure and atmosphere choice are left as
!               they were.
!   table       a tabulated surface, on a synthetic file in the Allard
!               layout that the program writes and atm_init reads (log P
!               a polynomial in log Teff and log g, 4000 to 7000 K, log g
!               3 to 5.5): no atmosphere is integrated (no points, no
!               event); the surface point is the table's (log P within
!               1e-7 of the polynomial, the file's digits; log T = log
!               Teff), and the envelope starts there. A surface above the
!               table (7300 K): the kernel switches to gray, and the
!               integration is the gray atmosphere's, bit for bit.
!               Through atm_get: on the table the star's photospheric
!               pressure is the table's and its choice in use stays;
!               above it the choice in use becomes gray. (atm's and
!               core's lines for the point above the table are part of
!               the report.)
!
! Usage: test_envelope_integration [DECK.inlist]   (run with YREC_INPUT
! and YREC_START set and YREC_OUTPUT pointing at a directory the deck
! may write to; the deck defaults to core/test/test_structure_tams.inlist
! from src/). Each check prints its value (1pe9.1) and ok or FAIL;
! any FAIL stops the program with status 1.
program test_envelope_integration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use controls_reader, only: read_controls
      use run_files, only: open_run_files
      use run_setup, only: set_up_constants_and_tables
      use convection_lib, only: set_convection_settings
      use convection_settings_builder, only: convection_settings_from_star
      use star_info_lib, only: star_info, control_nml_override, physics_nml_override
      use composition_layout, only: n_species_extended, i_h1, i_he4, i_metals
      use envelope_points, only: envelope_structure_state, &
           atmosphere_structure_state, atmosphere_point_state
      use envelope_integrands, only: point_physics, compute_atmosphere_derivs, &
           compute_envelope_derivs
      use envelope_integration_kernel, only: envint_config, envint_events, &
           fitting_vertex, integrate_envelope_atmosphere, envelope_termination, &
           max_log_fit_point_mass_fraction, no_photosphere_log10_pressure, &
           max_start_tries, envint_start_retry, envint_start_exhausted, &
           envint_atmosphere_done, envint_atmosphere_maxstp, &
           envint_bsstep_limits, envint_bsstep_failed, envint_vertex, &
           envint_envelope_done, envint_envelope_maxstp
      use atm_lib, only: atm_eddington_gray, atm_allard, hsra_t_tau_offset, &
           ttau_photosphere_x_limit, atm_init, atm_settings_t, atm_table_paths_t
      use numerics_lib, only: numerics_bad_input, numerics_step_underflow, &
           numerics_derivatives_failed
      use physical_constants, only: ln10, cgl, c4pil, csigl
      use rescaling_kernels, only: compute_envelope_mixture
      use settings_install, only: install_surface_mixture_settings
      use envelope_integration, only: gray_photosphere, atm_get
      implicit none
! the test's star (star_info_lib declares no instance); static, as it
! is too large for the stack
      type(star_info), save :: star
! the results of one integration, static for their size
      type(envelope_structure_state), save :: env_points
      type(atmosphere_structure_state), save :: atm_points
      type(atmosphere_point_state) :: atm_point
      type(envint_events) :: events
      type(envint_config) :: base_cfg

      real(dp), parameter :: x_surface = 0.7196_dp, z_surface = 0.016492_dp
      real(dp), parameter :: log10_fit_mass = -1.0e-4_dp
! How far a Bulirsch-Stoer step may stray from the reference, in units of
! the integration's error tolerance: the integrands interpolate tables,
! so they are not smooth and bsstep's error estimate is no bound (the
! largest seen is 55 tolerances, in the atmosphere).
      real(dp), parameter :: step_tolerance_factor = 1.0e3_dp
      real(dp), parameter :: log_mass_floor = 1.0e-6_dp
      real(dp) :: log10_teff, log10_star_mass, log10_radius_surface, &
           log10_gravity, composition(n_species_extended)
      real(dp) :: log10_photo_pressure
! the vertex the integrations fit (vertex 1)
      type(fitting_vertex) :: vertex
      logical :: switched_to_gray
      character(len=512) :: deck
      integer :: ierr, nfail

      if (command_argument_count() > 1) then
         write(*,'(a)') 'usage: test_envelope_integration [DECK.inlist]'
         stop 2
      end if
! Without an argument, test_structure's TAMS deck, as from src/ (where
! tools/coverage.py runs the test programs).
      deck = 'core/test/test_structure_tams.inlist'
      if (command_argument_count() == 1) call get_command_argument(1, deck)
      nfail = 0
      call load_deck(deck)
      call set_surface()
      write(*,'(a)') 'test_envelope_integration report'
      write(*,'(a)') '-- integrands'
      call check_integrands()
      write(*,'(a)') '-- solve'
      call check_solve()
      write(*,'(a)') '-- depths'
      call check_atmosphere_depths()
      write(*,'(a)') '-- reference'
      call check_atmosphere_steps()
      call check_envelope_steps()
      write(*,'(a)') '-- thin'
      call check_thin_envelope()
      write(*,'(a)') '-- retry'
      call check_start_retry()
      write(*,'(a)') '-- errors'
      call check_errors()
      write(*,'(a)') '-- gray photosphere'
      call check_gray_photosphere()
      write(*,'(a)') '-- table'
      call check_table_surface()
      if (nfail > 0) then
         write(*,'(a,i0,a)') 'test_envelope_integration: ', nfail, ' FAIL'
         stop 1
      end if
      write(*,'(a)') 'test_envelope_integration: done'

contains

! ---------------------------------------------------------------
! The deck's controls, tables and convection settings, as run_yrec sets
! them up before reading its first model.
subroutine load_deck(path)
      character(len=*), intent(in) :: path
      ierr = 0
      control_nml_override = path
      physics_nml_override = path
      call read_controls(star, ierr)
      if (ierr /= 0) call die('read_controls')
      call open_run_files(star, ierr)
      if (ierr /= 0) call die('open_run_files')
      call set_up_constants_and_tables(star, ierr)
      if (ierr /= 0) call die('set_up_constants_and_tables')
      star%current_card = 1
      star%current_card_envelope_x = star%card_initial_x(1)
      star%current_card_envelope_z = star%card_initial_z(1)
      star%mixing_length_alpha = star%card_mixing_length_alpha(1)
      call set_convection_settings(convection_settings_from_star(star))
end subroutine load_deck

subroutine die(what)
      character(len=*), intent(in) :: what
      write(*,'(a)') 'test_envelope_integration: FAIL ('//what//')'
      stop 1
end subroutine die

! ---------------------------------------------------------------
! The solar surface, and the kernel's settings as atm_get fills them
! from the star, with the gray atmosphere and the test's fitting mass.
subroutine set_surface()
      log10_teff = log10(5777.0_dp)
      log10_star_mass = log10(star%solar_mass_cgs)
! L = 4 pi R^2 sigma Teff^4, with L = Lsun.
      log10_radius_surface = 0.5_dp*(star%log10_solar_luminosity - &
           4.0_dp*log10_teff - c4pil - csigl)
      log10_gravity = cgl + log10_star_mass - 2.0_dp*log10_radius_surface
      composition = 0.0_dp
      composition(i_h1) = x_surface
      composition(i_metals) = z_surface
      composition(i_he4) = 1.0_dp - x_surface - z_surface
! The envelope's mixture for the EOS, as a run's first model read installs
! it (update_surface_mixture): without it the EOS works from the run
! plan's unset mixture, where the Saha solver does not converge (an
! error since M293).
      star%envelope_hydrogen_fraction = x_surface
      star%envelope_metal_fraction = z_surface
      star%envelope_he3_fraction = 0.0_dp
      call compute_envelope_mixture(x_surface, z_surface, 0.0_dp, &
           star%physics%mixture_weights_seed, star%envelope_helium_fraction, &
           star%envelope_metals_no_cno, star%envelope_inverse_atomic_weight, &
           star%envelope_number_fractions)
      call install_surface_mixture_settings(star)
      base_cfg%atm_choice = atm_eddington_gray
      base_cfg%atm_step_begin = star%physics%atm_step_begin
      base_cfg%atm_step_min = star%physics%atm_step_min
      base_cfg%atm_step_max = star%physics%atm_step_max
      base_cfg%env_step_begin = star%physics%env_step_begin
      base_cfg%env_step_min = star%physics%env_step_min
      base_cfg%env_step_max = star%physics%env_step_max
      base_cfg%atm_error_tol = star%physics%atm_error_tol
      base_cfg%env_error_tol = star%physics%env_error_tol
      base_cfg%atm_step_initial = star%physics%atm_step_initial
      base_cfg%hsra_log_t_offset = hsra_t_tau_offset()
      base_cfg%log10_total_mass = log10_star_mass
      base_cfg%surface_composition = composition
      base_cfg%log_fit_point_mass_fraction = log10_fit_mass
      base_cfg%envelope_store_log_temperature = star%physics%envelope_store_log_temperature
      base_cfg%log10_solar_luminosity = star%log10_solar_luminosity
      base_cfg%log10_gravitational_constant = cgl
      base_cfg%calc_envelope = .true.
      base_cfg%tolerance_fraction = star%physics%tolerance_fraction
      base_cfg%max_stage_index = star%physics%max_stage_index
      base_cfg%extrap_order = star%physics%extrap_order
end subroutine set_surface

! ---------------------------------------------------------------
! One integration with the settings cfg at the solar surface (or at
! another log Teff or log g), fitting vertex 1; the results go to the
! program's variables.
subroutine integrate(cfg, teff, gravity)
      type(envint_config), intent(in) :: cfg
! log Teff and log g of the surface, when not the solar ones
      real(dp), intent(in), optional :: teff, gravity
      type(envint_config) :: work
      real(dp) :: surface_teff, surface_gravity
      integer :: saha_iterations
      work = cfg
      surface_teff = log10_teff
      if (present(teff)) surface_teff = teff
      surface_gravity = log10_gravity
      if (present(gravity)) surface_gravity = gravity
      vertex = fitting_vertex(index=1)
      saha_iterations = 0
      env_points%num_env_points = 0
      atm_points%num_atm_points = 0
      call integrate_envelope_atmosphere(work, 1.0_dp, 1.0_dp, 1.0_dp, &
           surface_gravity, log10_radius_surface, surface_teff, x_surface, &
           z_surface, 0.0_dp, saha_iterations, env_points, atm_points, &
           atm_point, switched_to_gray, log10_photo_pressure, events, ierr, &
           vertex)
end subroutine integrate

! ---------------------------------------------------------------
! An ok/FAIL line.
subroutine report(name, value, pass)
      character(len=*), intent(in) :: name
      real(dp), intent(in) :: value
      logical, intent(in) :: pass
      if (pass) then
         write(*,'(a,t60,1pe9.1,2x,a)') name, value, 'ok'
      else
         write(*,'(a,t60,1pe9.1,2x,a)') name, value, 'FAIL'
         nfail = nfail + 1
      end if
end subroutine report

! The kind of event k of the last integration (0 when there is none).
integer function event_kind(k)
      integer, intent(in) :: k
      event_kind = 0
      if (k >= 1 .and. k <= events%num_events) event_kind = events%event(k)%kind
end function event_kind

! ---------------------------------------------------------------
! The atmosphere's right-hand side at a point (dP/dtau, T(tau)), and
! the envelope's (d log m/d log P, d log r/d log P), each against its
! equation in linear form.
subroutine check_integrands()
      type(point_physics) :: pt
      type(atmosphere_point_state) :: point
      real(dp) :: y1(1), dydx1(1), y3(3), dydx3(3), radius, tau, &
           g, p, expected, mass
      integer :: saha_iterations

      saha_iterations = 0
      y1(1) = 4.5_dp
      call compute_atmosphere_derivs(log10(0.1_dp), y1, dydx1, 1.0_dp, &
           log10_gravity, .true., .true., log10_teff, x_surface, z_surface, &
           saha_iterations, atm_eddington_gray, base_cfg%hsra_log_t_offset, &
           composition, point, pt, ierr)
      tau = 10.0_dp**point%log10_tau
      g = 10.0_dp**log10_gravity
      p = 10.0_dp**y1(1)
      expected = g*tau/(pt%opacity*p)
      call report('atmosphere: dlog P/dlog tau = g tau/(kappa P)', &
           abs(dydx1(1)/expected - 1.0_dp), &
           ierr == 0 .and. abs(dydx1(1)/expected - 1.0_dp) <= 1.0e-12_dp)
      expected = 0.75_dp*(10.0_dp**log10_teff)**4*(tau + 2.0_dp/3.0_dp)
      call report('atmosphere: T^4 = (3/4) Teff^4 (tau + 2/3)', &
           abs((10.0_dp**point%log10_temperature)**4/expected - 1.0_dp), &
           abs((10.0_dp**point%log10_temperature)**4/expected - 1.0_dp) <= 1.0e-12_dp)

! A point 0.4 dex below the photosphere in P, 2e-5 of the mass in, at
! 0.9 Rsun-ish radius and 8000 K.
      y3 = [-2.0e-5_dp, log10(8000.0_dp), log10_radius_surface - 1.0e-3_dp]
      call compute_envelope_derivs(5.6_dp, y3, dydx3, 1.0_dp, 1.0_dp, 1.0_dp, &
           .false., .true., log10_teff, x_surface, z_surface, &
           saha_iterations, log10_star_mass, cgl, composition, pt, ierr)
      radius = 10.0_dp**y3(3)
      mass = 10.0_dp**(y3(1) + log10_star_mass)
      p = 10.0_dp**5.6_dp
      expected = -4.0_dp*acos(-1.0_dp)*radius**4*p/(10.0_dp**cgl*mass**2)
      call report('envelope: dlog m/dlog P = -4 pi r^4 P/(G m^2)', &
           abs(dydx3(1)/expected - 1.0_dp), &
           ierr == 0 .and. abs(dydx3(1)/expected - 1.0_dp) <= 1.0e-12_dp)
      expected = -radius*p/(10.0_dp**cgl*mass*pt%density)
      call report('envelope: dlog r/dlog P = -r P/(G m rho)', &
           abs(dydx3(3)/expected - 1.0_dp), &
           abs(dydx3(3)/expected - 1.0_dp) <= 1.0e-12_dp)
      call report('envelope: dlog T/dlog P = grad', &
           abs(dydx3(2) - pt%temperature_gradient), &
           dydx3(2) == pt%temperature_gradient .and. dydx3(2) > 0.0_dp)
end subroutine check_integrands

! ---------------------------------------------------------------
! The integration with the deck's settings: its status, its events, and
! where it ends.
subroutine check_solve()
      real(dp) :: dtau, dteff, dmass

      call integrate(base_cfg)
      call report('ierr = 0', real(ierr, dp), ierr == 0)
      call report('events: atmosphere done, vertex, envelope done', &
           real(events%num_events, dp), events%num_events == 3 .and. &
           event_kind(1) == envint_atmosphere_done .and. &
           event_kind(2) == envint_vertex .and. &
           event_kind(3) == envint_envelope_done)
      call report('not switched to gray', 0.0_dp, .not. switched_to_gray)
      dtau = abs(atm_point%log10_tau - log10(2.0_dp/3.0_dp))
      call report('photosphere at tau = 2/3', dtau, dtau <= base_cfg%tolerance_fraction)
      dteff = abs(10.0_dp**(atm_point%log10_temperature - log10_teff) - 1.0_dp)
      call report('T = Teff at the photosphere', dteff, dteff <= 1.0e-6_dp)
      call report('photospheric P = the last atmosphere point''s', &
           abs(log10_photo_pressure - atm_point%log10_pressure), &
           log10_photo_pressure == atm_point%log10_pressure .and. &
           atm_points%num_atm_points > 2)
      dmass = abs(env_points%env_log10_mass(1) - log10_fit_mass)
      call report('envelope ends at the fitting mass', dmass, &
           env_points%num_env_points > 2 .and. &
           dmass <= base_cfg%tolerance_fraction*abs(log10_fit_mass))
! 1/mu_e of the recorded points: between 1/2 (ionized helium and metals)
! and 1 (ionized hydrogen) times the ionized fraction, and largest at
! the innermost, most ionized point.
      call report('the points record 1/mu_e (largest at the fitting point)', &
           env_points%env_electron_mean_weight_inverse(1), &
           env_points%env_electron_mean_weight_inverse(1) > 0.0_dp .and. &
           env_points%env_electron_mean_weight_inverse(1) < 1.0_dp .and. &
           env_points%env_electron_mean_weight_inverse(1) > &
           env_points%env_electron_mean_weight_inverse(env_points%num_env_points))
      call report('the vertex event holds the vertex', &
           abs(events%event(2)%values(1) - vertex%log10_pressure), &
           all(events%event(2)%values(1:4) == [vertex%log10_pressure, &
           vertex%log10_temperature, vertex%log10_radius, log10_fit_mass]))
end subroutine check_solve

! ---------------------------------------------------------------
! The depth steps of the last solve's atmosphere points against the
! quadrature dz = dtau/<rho kappa>, the first from tau = 0 with
! rho kappa = 0 there.
subroutine check_atmosphere_depths()
      real(dp) :: tau, tau_before, rho_kappa, rho_kappa_before, first, worst
      integer :: k

      call integrate(base_cfg)
      tau_before = 0.0_dp
      rho_kappa_before = 0.0_dp
      first = 0.0_dp
      worst = 0.0_dp
      do k = 1, atm_points%num_atm_points
         tau = 10.0_dp**log10_tau_of(atm_points%atmo_log10_temperature(k))
         rho_kappa = 10.0_dp**atm_points%atmo_log10_density(k)*atm_points%atmo_opacity(k)
         if (k == 1) then
            first = abs(atm_points%atmo_delta_depth(1)*rho_kappa/(2.0_dp*tau) - 1.0_dp)
         else
            worst = max(worst, abs(atm_points%atmo_delta_depth(k)* &
                 0.5_dp*(rho_kappa + rho_kappa_before)/(tau - tau_before) - 1.0_dp))
         end if
         tau_before = tau
         rho_kappa_before = rho_kappa
      end do
      call report('depth: the first step is measured from tau = 0', first, &
           atm_points%num_atm_points > 2 .and. first <= 1.0e-6_dp)
      call report('depth: each step is dtau over the mean rho kappa', worst, &
           atm_points%num_atm_points > 2 .and. worst <= 1.0e-6_dp)
end subroutine check_atmosphere_depths

! ---------------------------------------------------------------
! Each atmosphere step of the last solve against RK4 between its ends,
! with log tau from the Eddington T(tau).
subroutine check_atmosphere_steps()
      type(envint_config) :: cfg
      real(dp) :: x0, x1, y(1), worst, scale
      integer :: k

! Each difference over the step's scale, |log P| + |its change|.
      cfg = base_cfg
      call integrate(cfg)
      worst = 0.0_dp
      do k = 2, atm_points%num_atm_points
         x0 = log10_tau_of(atm_points%atmo_log10_temperature(k-1))
         x1 = log10_tau_of(atm_points%atmo_log10_temperature(k))
         y(1) = atm_points%atmo_log10_pressure(k-1)
         call rk4_atmosphere(x0, x1, y)
         scale = abs(atm_points%atmo_log10_pressure(k-1)) + &
              abs(atm_points%atmo_log10_pressure(k) - atm_points%atmo_log10_pressure(k-1))
         worst = max(worst, abs(y(1) - atm_points%atmo_log10_pressure(k))/scale)
      end do
      call report('atmosphere steps = RK4 (log P)', worst, &
           atm_points%num_atm_points > 2 .and. &
           worst <= step_tolerance_factor*cfg%atm_error_tol)
end subroutine check_atmosphere_steps

! log tau at log T on the Eddington T(tau): tau = (4/3)(T/Teff)^4 - 2/3.
real(dp) function log10_tau_of(log10_temperature)
      real(dp), intent(in) :: log10_temperature
      log10_tau_of = log10((4.0_dp/3.0_dp)*10.0_dp**(4.0_dp*(log10_temperature - &
           log10_teff)) - 2.0_dp/3.0_dp)
end function log10_tau_of

! RK4 of d log P/d log tau from x0 to x1 (64 steps).
subroutine rk4_atmosphere(x0, x1, y)
      real(dp), intent(in) :: x0, x1
      real(dp), intent(inout) :: y(1)
      real(dp) :: h, x, k1(1), k2(1), k3(1), k4(1)
      integer :: i
      h = (x1 - x0)/64.0_dp
      x = x0
      do i = 1, 64
         call atmosphere_rhs(x, y, k1)
         call atmosphere_rhs(x + 0.5_dp*h, y + 0.5_dp*h*k1, k2)
         call atmosphere_rhs(x + 0.5_dp*h, y + 0.5_dp*h*k2, k3)
         call atmosphere_rhs(x + h, y + h*k3, k4)
         y = y + h*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
         x = x0 + i*h
      end do
end subroutine rk4_atmosphere

subroutine atmosphere_rhs(x, y, dydx)
      real(dp), intent(in) :: x, y(1)
      real(dp), intent(out) :: dydx(1)
      type(point_physics) :: pt
      type(atmosphere_point_state) :: point
      integer :: saha_iterations, jerr
      saha_iterations = 0
      call compute_atmosphere_derivs(x, y, dydx, 1.0_dp, log10_gravity, &
           .true., .true., log10_teff, x_surface, z_surface, saha_iterations, &
           atm_eddington_gray, base_cfg%hsra_log_t_offset, composition, &
           point, pt, jerr)
      if (jerr /= 0) call die('atmosphere_rhs')
end subroutine atmosphere_rhs

! ---------------------------------------------------------------
! Each envelope step of the last solve against RK4 between its ends
! (env_points runs from the fitting point out), and the vertex as the
! interpolation of the last step's ends.
subroutine check_envelope_steps()
      real(dp) :: y(3), worst, scale, weight, expected(3), worst_vertex
      integer :: k, n, j

! Each difference over the step's scale, as bsstep measures its error:
! |y| + |the change over the step|. The points' log(m/M) is stored as
! log m - log M, which rounds at about 1e-14; log_mass_floor, added to
! its scale, keeps that round-off below the bound.
      n = env_points%num_env_points
      worst = 0.0_dp
      do k = n, 2, -1
         do j = 1, 3
            y(j) = point_value(k, j)
         end do
         call rk4_envelope(env_points%env_log10_pressure(k), &
              env_points%env_log10_pressure(k-1), y)
         do j = 1, 3
            scale = abs(point_value(k, j)) + abs(point_value(k-1, j) - point_value(k, j))
            if (j == 1) scale = scale + log_mass_floor
            worst = max(worst, abs(y(j) - point_value(k-1, j))/scale)
         end do
      end do
      call report('envelope steps = RK4 (log m, log T, log r)', worst, &
           n > 2 .and. worst <= step_tolerance_factor*base_cfg%env_error_tol)
! The vertex: v = v1 + w (v2 - v1), w = (log m_fit - y1)/(y2 - y1), with
! 1 the innermost point and 2 the next one out.
      weight = (log10_fit_mass - env_points%env_log10_mass(1))/ &
           (env_points%env_log10_mass(2) - env_points%env_log10_mass(1))
      expected(1) = env_points%env_log10_pressure(1) + weight* &
           (env_points%env_log10_pressure(2) - env_points%env_log10_pressure(1))
      expected(2) = env_points%env_log10_temperature(1) + weight* &
           (env_points%env_log10_temperature(2) - env_points%env_log10_temperature(1))
      expected(3) = env_points%env_log10_radius(1) + weight* &
           (env_points%env_log10_radius(2) - env_points%env_log10_radius(1))
      worst_vertex = maxval(abs([vertex%log10_pressure, vertex%log10_temperature, &
           vertex%log10_radius] - expected))
      call report('vertex = interpolation to the fitting mass', worst_vertex, &
           worst_vertex <= 1.0e-10_dp)
end subroutine check_envelope_steps

! Variable j (1 log m, 2 log T, 3 log r) of envelope point k.
real(dp) function point_value(k, j)
      integer, intent(in) :: k, j
      select case (j)
      case (1); point_value = env_points%env_log10_mass(k)
      case (2); point_value = env_points%env_log10_temperature(k)
      case default; point_value = env_points%env_log10_radius(k)
      end select
end function point_value

! RK4 of the envelope's equations in log P from x0 to x1 (64 steps).
subroutine rk4_envelope(x0, x1, y)
      real(dp), intent(in) :: x0, x1
      real(dp), intent(inout) :: y(3)
      real(dp) :: h, x, k1(3), k2(3), k3(3), k4(3)
      integer :: i
      h = (x1 - x0)/64.0_dp
      x = x0
      do i = 1, 64
         call envelope_rhs(x, y, k1)
         call envelope_rhs(x + 0.5_dp*h, y + 0.5_dp*h*k1, k2)
         call envelope_rhs(x + 0.5_dp*h, y + 0.5_dp*h*k2, k3)
         call envelope_rhs(x + h, y + h*k3, k4)
         y = y + h*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
         x = x0 + i*h
      end do
end subroutine rk4_envelope

subroutine envelope_rhs(x, y, dydx)
      real(dp), intent(in) :: x, y(3)
      real(dp), intent(out) :: dydx(3)
      type(point_physics) :: pt
      integer :: saha_iterations, jerr
      saha_iterations = 0
      call compute_envelope_derivs(x, y, dydx, 1.0_dp, 1.0_dp, 1.0_dp, &
           .false., .true., log10_teff, x_surface, z_surface, &
           saha_iterations, log10_star_mass, cgl, composition, pt, jerr)
      if (jerr /= 0) call die('envelope_rhs')
end subroutine envelope_rhs

! ---------------------------------------------------------------
! A fitting mass above the largest the kernel integrates to: the vertex
! is the surface point.
subroutine check_thin_envelope()
      type(envint_config) :: cfg
      cfg = base_cfg
      cfg%log_fit_point_mass_fraction = 0.5_dp*max_log_fit_point_mass_fraction
      call integrate(cfg)
      call report('thin: ierr = 0, no envelope points', real(ierr, dp), &
           ierr == 0 .and. env_points%num_env_points == 0)
      call report('thin: no envelope event (the atmosphere''s is the last)', &
           real(events%num_events, dp), events%num_events == 1 .and. &
           event_kind(1) == envint_atmosphere_done)
      call report('thin: vertex = the surface point', &
           abs(vertex%log10_pressure - atm_point%log10_pressure), &
           vertex%log10_pressure == atm_point%log10_pressure .and. &
           vertex%log10_temperature == atm_point%log10_temperature .and. &
           vertex%log10_radius == log10_radius_surface)
end subroutine check_thin_envelope

! ---------------------------------------------------------------
! A starting density of 1e-2 g/cm^3 puts the start past the photosphere
! (a density within the tables, so that the devel build's traps stay
! quiet).
subroutine check_start_retry()
      type(envint_config) :: cfg
      real(dp) :: decades
      integer :: last
      cfg = base_cfg
      cfg%atm_step_initial = 1.0e-2_dp
      call integrate(cfg)
      call report('retry: an envint_start_retry event first', &
           real(event_kind(1), dp), event_kind(1) == envint_start_retry)
      if (event_kind(1) == envint_start_retry) then
         call report('retry: its start lies past the photosphere', &
              events%event(1)%values(1) - events%event(1)%values(2), &
              events%event(1)%values(1) > events%event(1)%values(2) .and. &
              events%event(1)%values(2) == ttau_photosphere_x_limit(atm_eddington_gray))
      end if
      call report('retry: then the integration succeeds', real(ierr, dp), &
           ierr == 0 .and. event_kind(events%num_events) == envint_envelope_done)
! At log g = -5 the start lies past the photosphere whatever the
! density: tau = kappa P/g, and P does not fall below the radiation's.
      call integrate(base_cfg, gravity=-5.0_dp)
      last = events%num_events
      call report('retry: no start before the photosphere: ierr = 1', &
           real(ierr, dp), ierr == 1 .and. event_kind(last) == envint_start_exhausted &
           .and. last == max_start_tries)
      if (event_kind(last) == envint_start_exhausted) then
         decades = log10(base_cfg%atm_step_initial/events%event(last)%values(3))
         call report('retry: it stops after max_start_tries densities', &
              decades, events%event(last)%counts(1) == max_start_tries .and. &
              abs(decades - real(max_start_tries - 1, dp)) <= 1.0e-9_dp .and. &
              events%event(last)%values(1) > events%event(last)%values(2))
      end if
end subroutine check_start_retry

! ---------------------------------------------------------------
! The error contracts.
subroutine check_errors()
      use, intrinsic :: ieee_exceptions, only: ieee_get_halting_mode, &
           ieee_set_halting_mode, ieee_support_halting, ieee_set_flag, &
           ieee_divide_by_zero, ieee_invalid, ieee_overflow
      type(envint_config) :: cfg
      integer :: last
      logical :: halt_zero, halt_invalid, halt_overflow

      cfg = base_cfg
      cfg%atm_step_begin = 1.0e-6_dp
      cfg%atm_step_min = 1.0e-6_dp
      cfg%atm_step_max = 1.0e-6_dp
      call integrate(cfg)
      last = event_kind(events%num_events)
      call report('atmosphere maxstp: ierr = 1, its event last', real(ierr, dp), &
           ierr == 1 .and. last == envint_atmosphere_maxstp)
      call report('atmosphere maxstp: no photospheric pressure', &
           log10_photo_pressure, &
           log10_photo_pressure == no_photosphere_log10_pressure)

      cfg = base_cfg
      cfg%env_step_begin = 1.0e-6_dp
      cfg%env_step_min = 1.0e-6_dp
      cfg%env_step_max = 1.0e-6_dp
      call integrate(cfg)
      last = event_kind(events%num_events)
      call report('envelope maxstp: ierr = 1, its event last', real(ierr, dp), &
           ierr == 1 .and. last == envint_envelope_maxstp)
      call report('envelope maxstp: the photosphere was reached', &
           log10_photo_pressure, log10_photo_pressure == atm_point%log10_pressure &
           .and. log10_photo_pressure > no_photosphere_log10_pressure)

      cfg = base_cfg
      cfg%extrap_order = 1
      call integrate(cfg)
      last = event_kind(events%num_events)
      call report('bsstep limits: ierr = numerics_bad_input', real(ierr, dp), &
           ierr == numerics_bad_input .and. last == envint_bsstep_limits)
      if (last == envint_bsstep_limits) call report('bsstep limits: the event holds imax, nuse', &
           real(events%event(events%num_events)%counts(2), dp), &
           all(events%event(events%num_events)%counts == [cfg%max_stage_index, 1]))

! A zero error tolerance: bsstep divides the error estimate by it, so
! no stage is accepted and the step underflows. The division raises IEEE
! divide-by-zero (invalid for 0/0); a build that traps them (make
! devel_tests, -ffpe-trap) would stop there, so they do not halt during
! this one integration, and their flags are cleared after it.
      cfg = base_cfg
      cfg%atm_error_tol = 0.0_dp
      call ieee_get_halting_mode(ieee_divide_by_zero, halt_zero)
      call ieee_get_halting_mode(ieee_invalid, halt_invalid)
      if (ieee_support_halting(ieee_divide_by_zero)) &
           call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      if (ieee_support_halting(ieee_invalid)) &
           call ieee_set_halting_mode(ieee_invalid, .false.)
      call integrate(cfg)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      call ieee_set_flag(ieee_invalid, .false.)
      if (ieee_support_halting(ieee_divide_by_zero)) &
           call ieee_set_halting_mode(ieee_divide_by_zero, halt_zero)
      if (ieee_support_halting(ieee_invalid)) &
           call ieee_set_halting_mode(ieee_invalid, halt_invalid)
      last = event_kind(events%num_events)
      call report('step underflow: ierr = envelope_termination', real(ierr, dp), &
           ierr == envelope_termination .and. last == envint_bsstep_failed)
      if (last == envint_bsstep_failed) call report('step underflow: the event holds bsstep''s status', &
           real(events%event(events%num_events)%counts(1), dp), &
           events%event(events%num_events)%counts(1) == numerics_step_underflow)

! The integrands failing inside a step: a first envelope step of 40 in
! log P toward log(m/M) = -0.5 puts the step's midpoint substeps far past
! the microphysics' tables. The traps are held off as above (the failing
! points may overflow on their way to the error).
      cfg = base_cfg
      cfg%env_step_begin = 40.0_dp
      cfg%env_step_min = 40.0_dp
      cfg%env_step_max = 40.0_dp
      cfg%log_fit_point_mass_fraction = -0.5_dp
      call ieee_get_halting_mode(ieee_divide_by_zero, halt_zero)
      call ieee_get_halting_mode(ieee_invalid, halt_invalid)
      call ieee_get_halting_mode(ieee_overflow, halt_overflow)
      if (ieee_support_halting(ieee_divide_by_zero)) &
           call ieee_set_halting_mode(ieee_divide_by_zero, .false.)
      if (ieee_support_halting(ieee_invalid)) &
           call ieee_set_halting_mode(ieee_invalid, .false.)
      if (ieee_support_halting(ieee_overflow)) &
           call ieee_set_halting_mode(ieee_overflow, .false.)
      call integrate(cfg)
      call ieee_set_flag(ieee_divide_by_zero, .false.)
      call ieee_set_flag(ieee_invalid, .false.)
      call ieee_set_flag(ieee_overflow, .false.)
      if (ieee_support_halting(ieee_divide_by_zero)) &
           call ieee_set_halting_mode(ieee_divide_by_zero, halt_zero)
      if (ieee_support_halting(ieee_invalid)) &
           call ieee_set_halting_mode(ieee_invalid, halt_invalid)
      if (ieee_support_halting(ieee_overflow)) &
           call ieee_set_halting_mode(ieee_overflow, halt_overflow)
      last = event_kind(events%num_events)
      call report('derivatives failure in a step: ierr > 0 (an error)', real(ierr, dp), &
           ierr > 0 .and. last == envint_bsstep_failed)
      if (last == envint_bsstep_failed) call report('derivatives failure in a step: the event holds the code', &
           real(events%event(events%num_events)%counts(2), dp), &
           events%event(events%num_events)%counts(1) == numerics_derivatives_failed .and. &
           events%event(events%num_events)%counts(2) == ierr)
end subroutine check_errors

! ---------------------------------------------------------------
! gray_photosphere against the kernel's gray integration with the
! star's settings, and the star's integration state it must not touch.
subroutine check_gray_photosphere()
      real(dp) :: log10_pressure, log10_temperature, dteff
      real(dp), parameter :: unset = -7.0_dp
      logical :: untouched

      call integrate(base_cfg)
      star%log_total_mass = log10_star_mass
      star%log_L_lsun = 0.0_dp
      star%atm_choice_in_use = atm_eddington_gray
      star%atm_point%log10_pressure = unset
      star%atm_point%log10_temperature = unset
      star%atm_points%num_atm_points = -1
      star%photosphere_logP = unset
      call gray_photosphere(star, log10_teff, log10_gravity, x_surface, z_surface, &
           composition, log10_pressure, log10_temperature, ierr)
      call report('gray photosphere: ierr = 0', real(ierr, dp), ierr == 0)
      call report('gray photosphere: log P = the kernel''s', &
           abs(log10_pressure - log10_photo_pressure), &
           abs(log10_pressure - log10_photo_pressure) <= 1.0e-10_dp)
      dteff = abs(10.0_dp**(log10_temperature - log10_teff) - 1.0_dp)
      call report('gray photosphere: T = Teff', dteff, dteff <= 1.0e-6_dp)
      untouched = star%atm_point%log10_pressure == unset .and. &
           star%atm_point%log10_temperature == unset .and. &
           star%atm_points%num_atm_points == -1 .and. &
           star%photosphere_logP == unset .and. &
           star%atm_choice_in_use == atm_eddington_gray
      call report('gray photosphere: the star''s atmosphere untouched', 0.0_dp, untouched)
end subroutine check_gray_photosphere

! ---------------------------------------------------------------
! log P of the synthetic table at (log Teff, log g): cubic in log g and
! quadratic in log Teff, which the Allard lookup's 4-point Lagrange
! interpolations reproduce; about the solar photosphere's 5.0.
pure real(dp) function synthetic_log10_pressure(log_teff, log_g) result(log_p)
      real(dp), intent(in) :: log_teff, log_g
      real(dp) :: t, g
      t = log_teff - 3.76_dp
      g = log_g - 4.4_dp
      log_p = 5.0_dp + 0.62_dp*g - 2.1_dp*t + 0.05_dp*g**2 - 0.01_dp*g**3 + &
           0.3_dp*t**2 + 0.1_dp*g*t
end function synthetic_log10_pressure

! A file in the Allard layout (atm's allard_table_reader): 13 Teff from
! 4000 K by 250 K, log g 3.0 to 5.5 by 0.5, at [Fe/H] = [alpha/Fe] = 0;
! log P, and log P + 1 and log Teff + 0.1 at tau = 100 (not used).
subroutine write_synthetic_table(path)
      character(len=*), intent(in) :: path
      integer :: unit, row, column
      real(dp) :: teff, log_g
      open(newunit=unit, file=path, status='replace')
      write(unit,'(a)') ' 2 ALLARD # synthetic test table'
      write(unit,'(a)') 'Teff  logg  [M/H] [a/M] log P            log P(tau=100)   log T(tau=100)'
      do row = 1, 13
         teff = 4000.0_dp + 250.0_dp*(row - 1)
         do column = 1, 6
            log_g = 2.5_dp + 0.5_dp*column
            write(unit,'(f6.0,3f6.2,1p3d16.8)') teff, log_g, 0.0_dp, 0.0_dp, &
                 synthetic_log10_pressure(log10(teff), log_g), &
                 synthetic_log10_pressure(log10(teff), log_g) + 1.0_dp, log10(teff) + 0.1_dp
         end do
      end do
      close(unit)
end subroutine write_synthetic_table

! ---------------------------------------------------------------
! A tabulated surface and the fallback to gray above the table, in the
! kernel and through atm_get.
subroutine check_table_surface()
      character(len=*), parameter :: table_path = 'synthetic_allard_surface.all'
      real(dp), parameter :: unset = -7.0_dp
      type(atm_settings_t) :: settings
      type(envint_config) :: cfg
      type(fitting_vertex) :: fallback_vertex
      real(dp) :: expected, off, hot_log10_teff, fallback_pressure, fit_pressure
      integer :: n, fallback_points, saha_iterations

      call write_synthetic_table(table_path)
      settings%atm_choice = atm_allard
      call atm_init(settings, atm_table_paths_t('', table_path), ierr)
      if (ierr /= 0) call die('atm_init on the synthetic table')

! On the table.
      cfg = base_cfg
      cfg%atm_choice = atm_allard
      call integrate(cfg)
      call report('table: ierr = 0; no atmosphere points, no atmosphere event', &
           real(ierr, dp), ierr == 0 .and. atm_points%num_atm_points == 0 .and. &
           events%num_events == 2 .and. event_kind(1) == envint_vertex .and. &
           event_kind(2) == envint_envelope_done .and. .not. switched_to_gray)
      expected = synthetic_log10_pressure(log10_teff, log10_gravity)
      off = abs(atm_point%log10_pressure - expected)
      call report('table: the surface point is the table''s', &
           off, off <= 1.0e-7_dp .and. atm_point%log10_temperature == log10_teff .and. &
           log10_photo_pressure == atm_point%log10_pressure)
      n = env_points%num_env_points
      call report('table: the envelope starts at the surface point', &
           abs(env_points%env_log10_pressure(max(n, 1)) - atm_point%log10_pressure), &
           n > 2 .and. env_points%env_log10_pressure(max(n, 1)) == atm_point%log10_pressure &
           .and. env_points%env_log10_temperature(max(n, 1)) == log10_teff .and. &
           env_points%env_log10_radius(max(n, 1)) == log10_radius_surface .and. &
           env_points%env_log10_mass(max(n, 1)) == 0.0_dp)
      fit_pressure = vertex%log10_pressure

! Above the table: the kernel's fallback, against the gray atmosphere at
! the same surface.
      hot_log10_teff = log10(7300.0_dp)
      call integrate(cfg, teff=hot_log10_teff)
      call report('fallback: above the table: switched to gray, ierr = 0', real(ierr, dp), &
           ierr == 0 .and. switched_to_gray .and. event_kind(1) == envint_atmosphere_done)
      fallback_pressure = log10_photo_pressure
      fallback_points = atm_points%num_atm_points
      fallback_vertex = vertex
      cfg%atm_choice = atm_eddington_gray
      call integrate(cfg, teff=hot_log10_teff)
      call report('fallback: the integration is the gray atmosphere''s', &
           abs(fallback_pressure - log10_photo_pressure), ierr == 0 .and. &
           .not. switched_to_gray .and. fallback_pressure == log10_photo_pressure .and. &
           fallback_points == atm_points%num_atm_points .and. &
           fallback_vertex%log10_pressure == vertex%log10_pressure .and. &
           fallback_vertex%log10_temperature == vertex%log10_temperature .and. &
           fallback_vertex%log10_radius == vertex%log10_radius)

! Through atm_get, on a star whose choice in use is the table; the
! guess at the fitting point's pressure is the kernel's fitting point.
      star%log_total_mass = log10_star_mass
      star%log_fit_point_mass_fraction = log10_fit_mass
      star%atm_choice_in_use = atm_allard
      star%photosphere_logP = unset
      saha_iterations = 0
      call atm_get(star, 1.0_dp, 1.0_dp, 1.0_dp, log10_gravity, log10_radius_surface, &
           log10_teff, x_surface, z_surface, composition, fit_pressure, saha_iterations, &
           ierr)
      off = abs(star%photosphere_logP - expected)
      call report('atm_get: on the table: photosphere_logP the table''s', &
           off, ierr == 0 .and. off <= 1.0e-7_dp .and. star%atm_choice_in_use == atm_allard &
           .and. star%atm_points%num_atm_points == 0)
      call atm_get(star, 1.0_dp, 1.0_dp, 1.0_dp, log10_gravity, log10_radius_surface, &
           hot_log10_teff, x_surface, z_surface, composition, fit_pressure, &
           saha_iterations, ierr)
      call report('atm_get: above the table: the choice in use is gray', &
           real(star%atm_choice_in_use, dp), ierr == 0 .and. &
           star%atm_choice_in_use == atm_eddington_gray .and. &
           star%photosphere_logP == fallback_pressure)

! The run's own settings again (no table).
      settings%atm_choice = atm_eddington_gray
      call atm_init(settings, atm_table_paths_t('', ''), ierr)
end subroutine check_table_surface

end program test_envelope_integration

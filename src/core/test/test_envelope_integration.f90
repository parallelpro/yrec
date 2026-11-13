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
!               normal integration.
!   errors      the error contracts: the atmosphere or the envelope
!               that cannot finish in maxstp steps (steps of 1e-6)
!               returns ierr = 1 and its maxstp event, last; the
!               photospheric pressure stays -1e99 when the atmosphere
!               failed. Bulirsch-Stoer limits out of range (nuse = 1)
!               return numerics_bad_input and an envint_bsstep_limits
!               event with imax and nuse. A step that underflows (an
!               error tolerance of zero) returns envelope_termination
!               and an envint_bsstep_failed event with bsstep's status.
!
! Usage: test_envelope_integration DECK.inlist   (run with YREC_INPUT
! and YREC_START set and YREC_OUTPUT pointing at a directory the deck
! may write to). Each check prints its value (1pe9.1) and ok or FAIL;
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
           integrate_envelope_atmosphere, envelope_termination, &
           max_log_fit_point_mass_fraction, envint_start_retry, &
           envint_atmosphere_done, envint_atmosphere_maxstp, &
           envint_bsstep_limits, envint_bsstep_failed, envint_vertex, &
           envint_envelope_done, envint_envelope_maxstp
      use atm_lib, only: atm_eddington_gray, hsra_t_tau_offset, &
           ttau_photosphere_x_limit
      use numerics_lib, only: numerics_bad_input, numerics_step_underflow
      use physical_constants, only: ln10, cgl, c4pil, csigl
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
      real(dp) :: log10_photo_pressure, vtx_logp(3), vtx_logr(3), vtx_logt(3), &
           stored_state(4)
      logical :: switched_to_gray
      character(len=512) :: deck
      integer :: ierr, nfail

      if (command_argument_count() /= 1) then
         write(*,'(a)') 'usage: test_envelope_integration DECK.inlist'
         stop 2
      end if
      call get_command_argument(1, deck)
      nfail = 0
      call load_deck(deck)
      call set_surface()
      write(*,'(a)') 'test_envelope_integration report'
      write(*,'(a)') '-- integrands'
      call check_integrands()
      write(*,'(a)') '-- solve'
      call check_solve()
      write(*,'(a)') '-- reference'
      call check_atmosphere_steps()
      call check_envelope_steps()
      write(*,'(a)') '-- thin'
      call check_thin_envelope()
      write(*,'(a)') '-- retry'
      call check_start_retry()
      write(*,'(a)') '-- errors'
      call check_errors()
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
! One integration with the settings cfg, saving the boundary at vertex
! 1; the results go to the program's variables.
subroutine integrate(cfg)
      type(envint_config), intent(in) :: cfg
      type(envint_config) :: work
      real(dp) :: luminosity, fp, ft, log10_radius, teff, x, z
      integer :: stored_vertex, saha_iterations
      work = cfg
      luminosity = 1.0_dp
      fp = 1.0_dp
      ft = 1.0_dp
      log10_radius = log10_radius_surface
      teff = log10_teff
      x = x_surface
      z = z_surface
      stored_state = 0.0_dp
      stored_vertex = 0
      saha_iterations = 0
      vtx_logp = 0.0_dp
      vtx_logr = 0.0_dp
      vtx_logt = 0.0_dp
      env_points%num_env_points = 0
      atm_points%num_atm_points = 0
      call integrate_envelope_atmosphere(work, switched_to_gray, &
           log10_photo_pressure, luminosity, fp, ft, log10_gravity, &
           log10_star_mass, 1, .true., 0.0_dp, log10_radius, teff, x, z, &
           stored_state, stored_vertex, saha_iterations, vtx_logp, vtx_logr, &
           vtx_logt, env_points, atm_points, atm_point, events, ierr)
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
      real(dp) :: y1(1), dydx1(1), y3(3), dydx3(3), teff, radius, tau, &
           g, p, expected, log10_radius, mass
      integer :: saha_iterations

      teff = log10_teff
      saha_iterations = 0
      y1(1) = 4.5_dp
      call compute_atmosphere_derivs(log10(0.1_dp), y1, dydx1, 1.0_dp, &
           log10_gravity, .true., .true., teff, x_surface, z_surface, &
           saha_iterations, atm_eddington_gray, base_cfg%hsra_log_t_offset, &
           composition, point, pt, ierr)
      tau = 10.0_dp**point%log10_tau
      g = 10.0_dp**log10_gravity
      p = 10.0_dp**y1(1)
      expected = g*tau/(pt%qo*p)
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
      log10_radius = y3(3)
      call compute_envelope_derivs(5.6_dp, y3, dydx3, 1.0_dp, 1.0_dp, 1.0_dp, &
           .false., .true., log10_radius, teff, x_surface, z_surface, &
           saha_iterations, log10_star_mass, cgl, composition, pt, ierr)
      radius = 10.0_dp**y3(3)
      mass = 10.0_dp**(y3(1) + log10_star_mass)
      p = 10.0_dp**5.6_dp
      expected = -4.0_dp*acos(-1.0_dp)*radius**4*p/(10.0_dp**cgl*mass**2)
      call report('envelope: dlog m/dlog P = -4 pi r^4 P/(G m^2)', &
           abs(dydx3(1)/expected - 1.0_dp), &
           ierr == 0 .and. abs(dydx3(1)/expected - 1.0_dp) <= 1.0e-12_dp)
      expected = -radius*p/(10.0_dp**cgl*mass*pt%qd)
      call report('envelope: dlog r/dlog P = -r P/(G m rho)', &
           abs(dydx3(3)/expected - 1.0_dp), &
           abs(dydx3(3)/expected - 1.0_dp) <= 1.0e-12_dp)
      call report('envelope: dlog T/dlog P = grad', &
           abs(dydx3(2) - pt%current_gradients(3)), &
           dydx3(2) == pt%current_gradients(3) .and. dydx3(2) > 0.0_dp)
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
      call report('the vertex event holds the vertex', &
           abs(events%event(2)%values(1) - vtx_logp(1)), &
           all(events%event(2)%values(1:4) == [vtx_logp(1), vtx_logt(1), &
           vtx_logr(1), log10_fit_mass]))
end subroutine check_solve

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
      real(dp) :: teff
      integer :: saha_iterations, jerr
      teff = log10_teff
      saha_iterations = 0
      call compute_atmosphere_derivs(x, y, dydx, 1.0_dp, log10_gravity, &
           .true., .true., teff, x_surface, z_surface, saha_iterations, &
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
      worst_vertex = maxval(abs([vtx_logp(1), vtx_logt(1), vtx_logr(1)] - expected))
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
      real(dp) :: teff, log10_radius
      integer :: saha_iterations, jerr
      teff = log10_teff
      log10_radius = y(3)
      saha_iterations = 0
      call compute_envelope_derivs(x, y, dydx, 1.0_dp, 1.0_dp, 1.0_dp, &
           .false., .true., log10_radius, teff, x_surface, z_surface, &
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
      call report('thin: vertex = the surface point', &
           abs(vtx_logp(1) - atm_point%log10_pressure), &
           vtx_logp(1) == atm_point%log10_pressure .and. &
           vtx_logt(1) == atm_point%log10_temperature .and. &
           vtx_logr(1) == log10_radius_surface)
end subroutine check_thin_envelope

! ---------------------------------------------------------------
! A starting density of 1e3 g/cm^3 puts the start past the photosphere.
subroutine check_start_retry()
      type(envint_config) :: cfg
      cfg = base_cfg
      cfg%atm_step_initial = 1.0e3_dp
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
end subroutine check_start_retry

! ---------------------------------------------------------------
! The error contracts.
subroutine check_errors()
      type(envint_config) :: cfg
      integer :: last

      cfg = base_cfg
      cfg%atm_step_begin = 1.0e-6_dp
      cfg%atm_step_min = 1.0e-6_dp
      cfg%atm_step_max = 1.0e-6_dp
      call integrate(cfg)
      last = event_kind(events%num_events)
      call report('atmosphere maxstp: ierr = 1, its event last', real(ierr, dp), &
           ierr == 1 .and. last == envint_atmosphere_maxstp)
      call report('atmosphere maxstp: no photospheric pressure', &
           log10_photo_pressure, log10_photo_pressure == -1.0e99_dp)

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
           .and. log10_photo_pressure > -1.0e98_dp)

      cfg = base_cfg
      cfg%extrap_order = 1
      call integrate(cfg)
      last = event_kind(events%num_events)
      call report('bsstep limits: ierr = numerics_bad_input', real(ierr, dp), &
           ierr == numerics_bad_input .and. last == envint_bsstep_limits)
      if (last == envint_bsstep_limits) call report('bsstep limits: the event holds imax, nuse', &
           real(events%event(events%num_events)%counts(2), dp), &
           all(events%event(events%num_events)%counts == [cfg%max_stage_index, 1]))

      cfg = base_cfg
      cfg%atm_error_tol = 0.0_dp
      call integrate(cfg)
      last = event_kind(events%num_events)
      call report('step underflow: ierr = envelope_termination', real(ierr, dp), &
           ierr == envelope_termination .and. last == envint_bsstep_failed)
      if (last == envint_bsstep_failed) call report('step underflow: the event holds bsstep''s status', &
           real(events%event(events%num_events)%counts(1), dp), &
           events%event(events%num_events)%counts(1) == numerics_step_underflow)
end subroutine check_errors

end program test_envelope_integration

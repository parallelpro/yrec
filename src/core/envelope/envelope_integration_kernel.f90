!----------------------------------------------------------------------
! envelope_integration_kernel
!----------------------------------------------------------------------
! The integration of a model's outer layers, from the surface down to
! the fitting point. Its one entry, integrate_envelope_atmosphere, does
! three things in turn:
!   1. the surface: for a tabulated atmosphere (atm_choice 3, 4, 5) the
!      surface log P and log T from the table;
!   2. otherwise the atmosphere: log P against log tau, with T from the
!      T-tau relation, from tau near 0 down to the photosphere;
!   3. the envelope: log(m/M), log T and log r against log P, from the
!      surface down to the fitting mass, log(m/M) = log_fit_point_mass_fraction.
! On the way it records each atmosphere point in atm_points and each
! envelope point in env_points (the records it is given are its product).
! The physics of each point comes from point_physics, where the
! integrands of envelope_integrands leave it.
!
! The three parts are module procedures (find_table_surface,
! integrate_atmosphere, integrate_envelope, and the steps under each)
! that work on one record, integration_state: the equations with their
! surface, and the Bulirsch-Stoer integration in progress. Each takes
! the state, the settings and the records it writes as arguments.
!
! The kernel is star-blind: its code uses neither star_info_lib nor
! star%, which src/tools/check_boundaries.py checks (STAR_BLIND_FILES).
! Its settings arrive in an envint_config that the driver, atm_get in
! envelope_integration, fills from the star's controls and
! atm_choice_in_use; it does not change them. Its two effects on the
! star, the photospheric pressure (star%photosphere_logP) and the switch
! to a gray atmosphere when the Allard table lookup fails, come back as
! outputs that the driver applies.
!
! The kernel prints nothing. What it has to say (its notes, its solver
! diagnostics and the reason for a failure) comes back as events
! (envint_events), in the order they happened, and atm_get reports
! them (envelope_integration's report_envint_events); a failure also
! returns ierr.
module envelope_integration_kernel
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: n_species_extended
      use envelope_points, only: envelope_structure_state, &
           atmosphere_structure_state, atmosphere_point_state
      use envelope_integrands, only: compute_atmosphere_derivs, compute_envelope_derivs, &
           gamma1_from_eos_derivatives, point_physics
      use numerics_lib, only: bsstep, numerics_bad_input, numerics_derivatives_failed, &
           ode_system
      use eos_lib, only: eos_get, eos_result
      use kap_lib, only: i_kap, i_log10_kap, kap_get, num_kap_results
      use atm_lib, only: atm_choice_is_tabulated, atm_eddington_gray, atm_get_surface_pt, &
           ttau_photosphere_x_limit, ttau_start_log10_temperature
      use physical_constants, only: cdelrl, gas_constant, ln10, radiation_constant_over_3
      use math_lib
      implicit none
      private
      public :: integrate_envelope_atmosphere

! The kernel's settings (read only): the atmosphere choice; the Bulirsch-Stoer step
! sizes (begin, min, max, in log tau and in log P) and the error
! tolerances of the two integrations; the starting density of the
! atmosphere integration (atm_step_initial, legacy ATMD0, g/cm^3); the
! HSRA offset hsra_log_t_offset; the fitting mass
! log_fit_point_mass_fraction = log10(m_fit/M);
! envelope_store_log_temperature, the log T at which the envelope state
! is stored for the surface boundary; log10 Lsun; log10 G (cgs); log10
! of the total mass in g, the origin of the envelope's mass coordinate
! log(m/M) (star%log_total_mass); and calc_envelope (false: the
! photosphere only); and the Bulirsch-Stoer controls: the tolerance
! fraction (the atmosphere's step tolerance, and the envelope's stop at
! tolerance_fraction |log_fit_point_mass_fraction| of
! log_fit_point_mass_fraction), the number of stages and the
! extrapolation order (physics%tolerance_fraction, physics%max_stage_index,
! physics%extrap_order).
      type, public :: envint_config
            integer :: atm_choice
            real(dp) :: atm_step_begin, atm_step_min, atm_step_max
            real(dp) :: env_step_begin, env_step_min, env_step_max
            real(dp) :: atm_error_tol, env_error_tol, atm_step_initial
            real(dp) :: hsra_log_t_offset, log_fit_point_mass_fraction
            real(dp) :: envelope_store_log_temperature, log10_solar_luminosity
            real(dp) :: log10_gravitational_constant
            real(dp) :: log10_total_mass
! the outer shell's composition (star%xa layout), for the EOS's
! Debye-Hueckel correction
            real(dp) :: surface_composition(n_species_extended)
            logical :: calc_envelope
            real(dp) :: tolerance_fraction
            integer :: max_stage_index, extrap_order
      end type envint_config

! The largest fitting mass log_fit_point_mass_fraction = log10(m_fit/M)
! the code allows (rebuild_envelope clamps log_fit_point_mass_fraction
! to it): an envelope of at least 1 - 10^log_fit_point_mass_fraction =
! ln(10) 1e-12, about 2.3e-12 of the star. integrate_envelope skips the
! envelope integration when log_fit_point_mass_fraction lies above it,
! and integrates at it.
      real(dp), parameter, public :: max_log_fit_point_mass_fraction = -1.0d-12

! The photospheric log P the kernel returns when the atmosphere did not
! finish (an error before the photosphere): no pressure is this low, and
! a driver that tests log10_photo_pressure > no_photosphere_log10_pressure
! keeps the star's photosphere_logP after such an error.
      real(dp), parameter, public :: no_photosphere_log10_pressure = -1.0d99

! The status of an envelope integration whose Bulirsch-Stoer step
! vanished (bsstep's numerics_step_underflow: the solution diverged; the
! historical "solution diverged, stop", which exited 0): negative, so
! that main and yrec_capi tell it from an error. A failure of the integrands inside bsstep (the
! EOS or the opacity) is not this: it returns their own positive code, an
! error, as it does at a step's end (M132).
      integer, parameter, public :: envelope_termination = -2

! The events of an integration, for the driver to report. Each is one of
! the kinds below with the numbers its message shows (counts and
! values); unlisted slots are zero.
!   envint_start_retry          the atmosphere's starting log tau lay
!                               past the photosphere's, so the starting
!                               density was divided by 10 (values: the
!                               starting log tau, the photosphere's)
!   envint_atmosphere_done      the atmosphere reached the photosphere
!                               (counts: steps accepted, rejected;
!                               values: the summed relative error in P)
!   envint_atmosphere_maxstp    the atmosphere did not reach the
!                               photosphere in max_steps steps (counts:
!                               max_steps; ierr 1)
!   envint_bsstep_limits        the Bulirsch-Stoer limits are out of
!                               range (counts: max_stage_index,
!                               extrap_order; ierr numerics_bad_input)
!   envint_bsstep_failed        a Bulirsch-Stoer step failed (counts:
!                               bsstep's status, and the integrands'
!                               code when they failed; values: where
!                               the step started, log tau or log P, and
!                               the number of equations, 1 the
!                               atmosphere's and 3 the envelope's; ierr:
!                               that code, or envelope_termination when
!                               the step underflowed)
!   envint_vertex               the fitting point of the vertex being
!                               saved (values: log P, log T, log R,
!                               log(m/M))
!   envint_envelope_done        the envelope reached the fitting mass
!                               (counts: steps accepted, rejected;
!                               values: log(m/M) left to go, the last
!                               log(m/M), and the summed relative errors
!                               in log m, log T and log r)
!   envint_envelope_maxstp      the envelope did not reach the fitting
!                               mass in max_steps steps (counts:
!                               max_steps; ierr 1)
!   envint_start_exhausted      the atmosphere's starting log tau lay
!                               past the photosphere's at every one of
!                               max_start_tries starting densities
!                               (counts: max_start_tries; values: the
!                               last starting log tau, the
!                               photosphere's, the last starting
!                               density; ierr 1)
      integer, parameter, public :: envint_start_retry = 1, &
           envint_atmosphere_done = 2, envint_atmosphere_maxstp = 3, &
           envint_bsstep_limits = 4, envint_bsstep_failed = 5, &
           envint_vertex = 6, envint_envelope_done = 7, &
           envint_envelope_maxstp = 8, envint_start_exhausted = 9
      type, public :: envint_event
            integer :: kind = 0
            integer :: counts(2) = 0
            real(dp) :: values(5) = 0.0_dp
      end type envint_event
      type, public :: envint_events
            integer :: num_events = 0
            type(envint_event), allocatable :: event(:)
      end type envint_events

! The vertex of the surface boundary's triangle that an integration
! fits (compute_surface_boundary's three envelopes), when one is being
! fitted:
!   index               which vertex, 1 to 3 (in);
!   log10_pressure, log10_radius, log10_temperature
!                       the fitting point, at the fitting mass: written
!                       when the envelope reaches the fitting mass (for an
!                       envelope too thin to integrate, the surface
!                       point), and left as given otherwise;
!   stored_state, stored_vertex_index
!                       the envelope state (log P, log T, log r, log(m/M))
!                       where log T first reaches
!                       cfg%envelope_store_log_temperature, which the
!                       rezone moves the fitting point to, and the vertex
!                       it was stored for (0: none). They come in as the
!                       star's: an integration of the envelope forgets a
!                       state stored for its own vertex, and stores a new
!                       one when its T gets there.
      type, public :: fitting_vertex
            integer :: index = 0
            real(dp) :: log10_pressure = 0.0_dp, log10_radius = 0.0_dp, &
                 log10_temperature = 0.0_dp
            real(dp) :: stored_state(4) = 0.0_dp
            integer :: stored_vertex_index = 0
      end type fitting_vertex

! The most steps either integration may take (raised from 200 to give
! the atmosphere integration room; DBG).
      integer, parameter :: max_steps = 2000
! The floor of a Bulirsch-Stoer step's error scale |y| + |h dy/dx|
! (legacy TINY), so that a variable at zero has a scale.
      real(dp), parameter :: step_scale_floor = 1.0d-30
! A retried start of the atmosphere divides its starting density by
! this, and an atmosphere tries at most max_start_tries starting
! densities: 12 decades below the control's value the pressure of the
! start is the radiation's, which the density no longer lowers, so a
! start still past the photosphere there will stay past it.
      real(dp), parameter :: start_density_divisor = 10.0d0
      integer, parameter, public :: max_start_tries = 12
! An atmosphere point's gradient is the radiative one up to this margin
! above the adiabatic gradient, and the adiabatic one beyond it.
      real(dp), parameter :: radiative_gradient_margin = 1.0d-6
! The EOS derivatives are always requested: the integrands store them
! in the point's physics for atm_points and env_points.
      logical, parameter :: want_eos_derivatives = .true.

! The atmosphere (compute_atmosphere_derivs: log P against log tau) or
! the envelope (compute_envelope_derivs: log(m/M), log T and log r
! against log P) as an ode_system for bsstep, with what the integrands
! read and write:
!   in_atmosphere       which of the two (also the EOS's atmosphere flag);
!   the surface         L (Lsun), f_P, f_T, log g, log Teff, X and Z;
!   from the settings   the T(tau) relation in use (cfg%atm_choice, or
!                       gray after a failed Allard lookup), its HSRA
!                       offset, log10 of the total mass (g), log10 G
!                       (cgs) and the outer shell's composition;
!   saha_iterations     the EOS's iteration counter, which they add to;
!   atm_point           the atmosphere's current point, which the
!                       atmosphere integrand writes;
!   pt                  the physics of the point last evaluated, which
!                       both write (point_physics).
      type, extends(ode_system) :: envelope_equations
            logical :: in_atmosphere = .true.
            real(dp) :: luminosity_linear = 0.0_dp, &
                 pressure_rotation_factor = 0.0_dp, &
                 temperature_rotation_factor = 0.0_dp, log10_gravity = 0.0_dp, &
                 log10_teff = 0.0_dp, hydrogen_fraction = 0.0_dp, &
                 metal_fraction = 0.0_dp
            integer :: atm_choice = 0
            real(dp) :: hsra_log_t_offset = 0.0_dp, log10_total_mass = 0.0_dp, &
                 log10_gravitational_constant = 0.0_dp
            real(dp) :: surface_composition(n_species_extended) = 0.0_dp
            integer :: saha_iterations = 0
            type(atmosphere_point_state) :: atm_point
            type(point_physics) :: pt
      contains
            procedure :: derivatives => equation_derivatives
      end type envelope_equations

! The state of one integration, which the kernel's procedures work on.
! Everything starts at zero, so that what an integration did not reach
! reads as zero (the step counts of a table's surface, the mass left of
! an envelope too thin to integrate).
!   eqs                 the equations, with the surface and the point
!                       last evaluated (above);
!   tabulated_surface   a table gave the surface point: no atmosphere
!                       is integrated.
! The Bulirsch-Stoer integration in progress, the atmosphere's and then
! the envelope's:
!   num_eqs             1 (the atmosphere: y(1) = log P against x = log
!                       tau) or 3 (the envelope: y = log(m/M), log T,
!                       log r against x = log P);
!   x, y, dydx          the independent variable, the solution and its
!                       derivatives where the next step starts;
!   x_start, y_start    where the last step started;
!   h_step, h_did, h_next   the step tried, taken, and suggested next;
!   h_min, h_max        the smallest next step and the largest step;
!   tolerance, y_scale  the step's error tolerance and error scale;
!   step_err, err_sum   the last step's scaled error estimate and the
!                       sum over the steps;
!   num_ok, num_bad     the steps taken at the size tried, and smaller;
!   end_tolerance       how near its end an integration stops: in log
!                       tau of the photosphere, in log(m/M) of the
!                       fitting mass;
!   underflow_reference bsstep's scale for a vanishing step, the surface
!                       hydrogen fraction (numerics' bsstep).
! The atmosphere:
!   log10_tau_photosphere   where it ends;
!   prev_tau, prev_density, prev_opacity   the previous point of the
!                       depth quadrature (store_atmosphere_point).
! The envelope:
!   mass_left           log(m_fit/M) - log(m/M) after the last step;
!   state_stored        the vertex's stored state has been written;
!   pressure_limit_pending   the step that would cross
!                       log10_pressure_limit has not come yet.
      type :: integration_state
            type(envelope_equations) :: eqs
            logical :: tabulated_surface = .false.
            integer :: num_eqs = 0
            real(dp) :: x = 0.0_dp, y(3) = 0.0_dp, dydx(3) = 0.0_dp
            real(dp) :: x_start = 0.0_dp, y_start(3) = 0.0_dp
            real(dp) :: h_step = 0.0_dp, h_did = 0.0_dp, h_next = 0.0_dp, &
                 h_min = 0.0_dp, h_max = 0.0_dp
            real(dp) :: tolerance = 0.0_dp, y_scale(3) = 0.0_dp
            real(dp) :: step_err(3) = 0.0_dp, err_sum(3) = 0.0_dp
            integer :: num_ok = 0, num_bad = 0
            real(dp) :: end_tolerance = 0.0_dp, underflow_reference = 0.0_dp
            real(dp) :: log10_tau_photosphere = 0.0_dp
            real(dp) :: prev_tau = 0.0_dp, prev_density = 0.0_dp, &
                 prev_opacity = 0.0_dp
            real(dp) :: mass_left = 0.0_dp
            logical :: state_stored = .false., pressure_limit_pending = .false.
      end type integration_state

contains

! dy/dx of the atmosphere or of the envelope (as self%in_atmosphere
! says) at (x, y).
subroutine equation_derivatives(self, x, y, dydx, ierr)
      class(envelope_equations), intent(inout) :: self
      real(dp), intent(in) :: x, y(:)
      real(dp), intent(out) :: dydx(:)
      integer, intent(out) :: ierr
      if (self%in_atmosphere) then
         call compute_atmosphere_derivs(x, y, dydx, &
              self%pressure_rotation_factor, &
              self%log10_gravity, self%in_atmosphere, want_eos_derivatives, &
              self%log10_teff, &
              self%hydrogen_fraction, self%metal_fraction, self%saha_iterations, &
              self%atm_choice, self%hsra_log_t_offset, self%surface_composition, &
              self%atm_point, self%pt, ierr)
      else
         call compute_envelope_derivs(x, y, dydx, self%luminosity_linear, &
              self%pressure_rotation_factor, self%temperature_rotation_factor, &
              self%in_atmosphere, want_eos_derivatives, self%log10_teff, &
              self%hydrogen_fraction, self%metal_fraction, self%saha_iterations, &
              self%log10_total_mass, self%log10_gravitational_constant, &
              self%surface_composition, self%pt, ierr)
      end if
end subroutine equation_derivatives

! ---------------------------------------------------------------
! integrate_envelope_atmosphere (the body of legacy envint)
! ---------------------------------------------------------------
! Integrates the atmosphere and the envelope from the surface down to
! the fitting point: the three parts of the module header, in turn, on
! one integration_state; each integration by Bulirsch-Stoer steps
! (bsstep), at most max_steps of them.
! Inputs: the settings cfg; the surface: L in Lsun (luminosity_linear),
! the rotational distortion factors f_P and f_T (1 without rotation),
! log g, log R, log Teff, X and Z; log10_pressure_limit, a guess at the
! pressure of the fitting point that the first envelope step to cross it
! stops at (when no vertex is being fitted); saha_iterations, the EOS's
! iteration counter, which the integration adds to.
! Outputs: atm_points; env_points, from the centre outward; atm_point,
! the last atmosphere point (the photosphere, or a table's surface
! point); switched_to_gray; the photospheric log P
! (log10_photo_pressure, no_photosphere_log10_pressure unless the
! atmosphere finished); the events to report (envint_events); ierr.
! vertex, when present: the triangle's vertex this integration fits
! (fitting_vertex), whose fitting point and stored envelope state it
! writes.
subroutine integrate_envelope_atmosphere(cfg, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, log10_gravity, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     log10_pressure_limit, saha_iterations, env_points, atm_points, &
     atm_point, switched_to_gray, log10_photo_pressure, events, ierr, vertex)
      type(envint_config), intent(in) :: cfg
      real(dp), intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, &
           log10_gravity
      real(dp), intent(in) :: log10_radius, log10_teff
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(in) :: log10_pressure_limit
      integer, intent(inout) :: saha_iterations
! The results (envelope_points): the envelope and atmosphere points, and
! the atmosphere integration's current point (whose last values the
! caller keeps).
      type(envelope_structure_state), intent(inout) :: env_points
      type(atmosphere_structure_state), intent(inout) :: atm_points
      type(atmosphere_point_state), intent(inout) :: atm_point
      logical, intent(out) :: switched_to_gray
      real(dp), intent(out) :: log10_photo_pressure
! events: what the integration has to report, in order (envint_events).
      type(envint_events), intent(out) :: events
! ierr: nonzero when the table lookup, the EOS, the opacity or either
! integration fails; an integration's own failure is also an event.
      integer, intent(out) :: ierr
      type(fitting_vertex), intent(inout), optional :: vertex

      type(integration_state) :: st

! No photosphere until the atmosphere finishes, so that the driver does
! not overwrite star%photosphere_logP after an early error.
      log10_photo_pressure = no_photosphere_log10_pressure

! The equations' surface and settings. The atmosphere's current point
! starts as the caller's: a table's surface writes only its log P and
! log T.
      st%eqs%in_atmosphere = .true.
      st%eqs%luminosity_linear = luminosity_linear
      st%eqs%pressure_rotation_factor = pressure_rotation_factor
      st%eqs%temperature_rotation_factor = temperature_rotation_factor
      st%eqs%log10_gravity = log10_gravity
      st%eqs%log10_teff = log10_teff
      st%eqs%hydrogen_fraction = hydrogen_fraction
      st%eqs%metal_fraction = metal_fraction
      st%eqs%atm_choice = cfg%atm_choice
      st%eqs%hsra_log_t_offset = cfg%hsra_log_t_offset
      st%eqs%log10_total_mass = cfg%log10_total_mass
      st%eqs%log10_gravitational_constant = cfg%log10_gravitational_constant
      st%eqs%surface_composition = cfg%surface_composition
      st%eqs%saha_iterations = saha_iterations
      st%eqs%atm_point = atm_point
      st%underflow_reference = hydrogen_fraction

      call find_table_surface(st, switched_to_gray, ierr)
      if (ierr == 0) call integrate_atmosphere(st, cfg, atm_points, events, ierr)
      if (ierr == 0) then
! The photosphere (or the table's surface point) is reached, and its
! pressure is the star's photosphere_logP. Stop here when only that is
! wanted; otherwise integrate the envelope in pressure from there down
! to the fitting mass.
         log10_photo_pressure = st%eqs%atm_point%log10_pressure
         if (cfg%calc_envelope) call integrate_envelope(st, cfg, log10_radius, &
              log10_pressure_limit, env_points, events, ierr, vertex)
      end if
! What the caller keeps, after a failure too.
      saha_iterations = st%eqs%saha_iterations
      atm_point = st%eqs%atm_point
end subroutine integrate_envelope_atmosphere

! ===============================================================
! Part 1: the surface from a table
! ===============================================================

! ---------------------------------------------------------------
! For a tabulated atmosphere (atm_choice_is_tabulated: 3 Kurucz, 4
! Allard, 5 Castelli/Kurucz), the surface log P and log T from the
! table, into the atmosphere's current point (legacy ATMPL). When the
! point lies outside the Allard tables (Teff above their maximum, or
! log g outside their range), the choice in use becomes the Eddington
! gray atmosphere (atm_choice 0, legacy KTTAU = 0), which
! integrate_atmosphere then integrates; the driver reports the switch.
! Sets st%tabulated_surface; ierr when a table interpolation fails.
subroutine find_table_surface(st, switched_to_gray, ierr)
      type(integration_state), intent(inout) :: st
      logical, intent(out) :: switched_to_gray
      integer, intent(out) :: ierr
      logical :: allard_lookup_failed

      ierr = 0
      switched_to_gray = .false.
      st%tabulated_surface = .false.
      if (.not. atm_choice_is_tabulated(st%eqs%atm_choice)) return
      call atm_get_surface_pt(st%eqs%atm_choice, st%eqs%log10_teff, &
           st%eqs%log10_gravity, st%eqs%atm_point%log10_pressure, &
           st%eqs%atm_point%log10_temperature, allard_lookup_failed, ierr)
      if (ierr /= 0) return
      if (allard_lookup_failed) then
         st%eqs%atm_choice = atm_eddington_gray
         switched_to_gray = .true.
      else
         st%tabulated_surface = .true.
      end if
end subroutine find_table_surface

! ===============================================================
! Part 2: the atmosphere
! ===============================================================

! ---------------------------------------------------------------
! The atmosphere from a T-tau relation (atm_choice 0 Eddington gray,
! 1 Krishna Swamy 1966, 2 HSRA). It finds the starting point near
! tau = 0 (T from the relation, P from a guessed density, tau from P),
! then integrates log P in log tau down to the photosphere (tau = 2/3;
! 0.312 for Krishna Swamy), storing atm_points on the way. Very
! occasionally the starting tau lies past the photosphere; then the
! guessed density is divided by 10 and the start found again, at most
! max_start_tries times in all. Skipped when a table gave the surface;
! ierr 1 when no start lay before the photosphere, or after max_steps
! steps.
subroutine integrate_atmosphere(st, cfg, atm_points, events, ierr)
      type(integration_state), intent(inout) :: st
      type(envint_config), intent(in) :: cfg
      type(atmosphere_structure_state), intent(inout) :: atm_points
      type(envint_events), intent(inout) :: events
      integer, intent(out) :: ierr
      real(dp) :: start_density
      integer :: try
! The EOS at the starting point. It starts each integration as a new
! record (eos_results: the defaults that start the EOS's solvers afresh)
! and keeps its values from one try to the next, so that a retry's EOS
! call starts from the previous try's solution.
      type(eos_result) :: eos_at_start

      ierr = 0
! No atmosphere points until the integration stores them, so that a
! table's surface leaves atm_points empty rather than holding the
! previous integration's points.
      atm_points%num_atm_points = 0
      if (st%tabulated_surface) return

! One equation, the atmosphere's error tolerance, and the end point: log
! tau of the photosphere for the relation in use.
      st%eqs%in_atmosphere = .true.
      st%num_eqs = 1
      st%tolerance = cfg%atm_error_tol
      st%log10_tau_photosphere = ttau_photosphere_x_limit(st%eqs%atm_choice)
      st%err_sum(1) = 0.0d0
      start_density = cfg%atm_step_initial
      do try = 1, max_start_tries
         call find_atmosphere_start(st, start_density, eos_at_start, ierr)
         if (ierr /= 0) return
         if (st%x > st%log10_tau_photosphere) then
! The start lies past the end point: start again from a lower density,
! unless this was the last try.
            if (try == max_start_tries) exit
            call add_envint_event(events, envint_start_retry, &
                 values=[st%x, st%log10_tau_photosphere])
            start_density = start_density/start_density_divisor
         else
            call integrate_to_photosphere(st, cfg, atm_points, events, ierr)
            return
         end if
      end do
      call add_envint_event(events, envint_start_exhausted, counts=[max_start_tries, 0], &
           values=[st%x, st%log10_tau_photosphere, start_density])
      ierr = 1
end subroutine integrate_atmosphere

! ---------------------------------------------------------------
! The starting point of the atmosphere integration. With the starting
! T (the T-tau relation near tau = 0; for HSRA, the fit at its tau = 2/3
! anchor) and the guessed density rho0 (start_density), an ideal gas of
! mean molecular weight 1 plus radiation gives
!   P = (R rho0 + (a/3) T^3) T.
! The EOS and opacity at (P, T) give kappa, and with kappa and g taken
! constant above this point, tau = kappa P / g:
!   x = log tau = log P - log g + log kappa,
!   dydx(1) = d log P / d log tau = g_eff tau / (kappa P),
! with g_eff = g f_P as in compute_atmosphere_derivs, i.e. f_P up to
! rounding.
subroutine find_atmosphere_start(st, start_density, eos_at_start, ierr)
      type(integration_state), intent(inout) :: st
      real(dp), intent(in) :: start_density
      type(eos_result), intent(inout) :: eos_at_start
      integer, intent(out) :: ierr
      real(dp) :: log10_temperature, temperature, log10_pressure
      real(dp) :: opacity, log10_opacity
      real(dp) :: kap_res(num_kap_results)

      log10_temperature = ttau_start_log10_temperature(st%eqs%log10_teff, &
           st%eqs%atm_choice, st%eqs%hsra_log_t_offset)
      temperature = exp(ln10*log10_temperature)
! P = (R rho0 + (a/3) T^3) T.
      log10_pressure = log10((gas_constant*start_density + &
           radiation_constant_over_3*temperature**3)*temperature)

! The EOS at (P, T). The record is inout in eos_get (its density and its
! Saha unknowns are the solvers' starting points).
      eos_at_start%temperature = temperature
      call eos_get(log10_temperature, log10_pressure, st%eqs%hydrogen_fraction, &
           st%eqs%metal_fraction, eos_at_start, want_eos_derivatives, &
           st%eqs%in_atmosphere, st%eqs%saha_iterations, &
           composition_at_zone=st%eqs%surface_composition, ierr=ierr)
      if (ierr /= 0) return
! The opacity at the density the EOS returned (DBG 12/95).
      call kap_get(eos_at_start%log10_density, log10_temperature, &
           st%eqs%hydrogen_fraction, st%eqs%metal_fraction, kap_res, &
           eos_at_start%ion_fraction, ierr=ierr)
      if (ierr /= 0) return
      opacity = kap_res(i_kap)
      log10_opacity = kap_res(i_log10_kap)
! The starting optical depth, tau = kappa P / g, and d log P / d log tau.
      st%x = log10_pressure - st%eqs%log10_gravity + log10(opacity)
      st%y(1) = log10_pressure
      st%dydx(1) = exp(ln10*(st%eqs%log10_gravity+st%x-log10_opacity-log10_pressure))* &
           st%eqs%pressure_rotation_factor
end subroutine find_atmosphere_start

! ---------------------------------------------------------------
! Integrate d log P/d log tau from the starting optical depth to the
! photosphere, storing each point in atm_points; done when log tau is
! within end_tolerance of the photosphere's. ierr 1, and an event, after
! max_steps steps.
subroutine integrate_to_photosphere(st, cfg, atm_points, events, ierr)
      type(integration_state), intent(inout) :: st
      type(envint_config), intent(in) :: cfg
      type(atmosphere_structure_state), intent(inout) :: atm_points
      type(envint_events), intent(inout) :: events
      integer, intent(out) :: ierr
      integer :: step

      ierr = 0
! The depth quadrature of store_atmosphere_point starts above the first
! point, at tau = 0 with rho kappa = 0 there (M322).
      st%prev_tau = 0.0d0
      st%prev_density = 0.0d0
      st%prev_opacity = 0.0d0
      st%h_max = cfg%atm_step_max
      st%h_min = cfg%atm_step_min
      st%h_step = cfg%atm_step_begin
      st%end_tolerance = cfg%tolerance_fraction
      st%num_ok = 0
      st%num_bad = 0

      do step = 1, max_steps
         call limit_atmosphere_step(st)
         call take_bulirsch_stoer_step(st, cfg, events, ierr)
         if (ierr /= 0) return
         call store_atmosphere_point(st, cfg, step, atm_points)
         if (abs(st%x - st%log10_tau_photosphere) <= st%end_tolerance) then
            call add_envint_event(events, envint_atmosphere_done, &
                 counts=[st%num_ok, st%num_bad], values=[st%err_sum(1)])
            return
         end if
         if (st%h_next < st%h_min) st%h_next = st%h_min
         st%h_step = st%h_next
      end do
! The integration did not reach the photosphere in max_steps steps.
      call add_envint_event(events, envint_atmosphere_maxstp, counts=[max_steps, 0])
      ierr = 1
end subroutine integrate_to_photosphere

! ---------------------------------------------------------------
! The error scale of the atmosphere's next step, and the step kept from
! going past the photosphere and within the largest step.
subroutine limit_atmosphere_step(st)
      type(integration_state), intent(inout) :: st
      call set_step_error_scale(st)
      if ((st%x-st%log10_tau_photosphere)*(st%x+st%h_step-st%log10_tau_photosphere) < 0.0d0) &
           st%h_step = st%log10_tau_photosphere - st%x
      if (st%h_step > st%h_max) st%h_step = st%h_max
end subroutine limit_atmosphere_step

! ---------------------------------------------------------------
! The new atmosphere point, as the last call of
! compute_atmosphere_derivs left it in atm_point and pt, into slot k of
! atm_points:
! - the thermodynamic state, with the gas-pressure fraction
!   beta = 1 - (a/3) T^4 / P and Gamma1 from the EOS derivatives;
! - the geometric depth between this point and the previous one,
!     dz = (tau - tau_prev) / [(rho kappa + rho_prev kappa_prev)/2]
!   (Cox & Giuli 1968, Principles of Stellar Structure, p. 590), the
!   first point's from tau = 0 with rho kappa = 0 there. The heights the
!   stitched model sums from these steps (its atmosphere radii; the
!   first step enters none) are a plane-parallel layer's: the gravity of
!   the integration is the photosphere's and does not fall with height
!   (2 (r - R)/R = 2.5e-3 at the top of the standard solar model's
!   atmosphere), and the step is the trapezoid of 1/(rho kappa) in tau
!   (about h^2/12 = 1e-3 of a step). The validation tier measures the
!   two together: the heights depart from hydrostatic equilibrium by up
!   to 0.6% of a step, 6e-6 R over the atmosphere, which its contract
!   test's atmosphere tolerance holds (docs/PHYSICS_VALIDATION.md, V11;
!   bug_sweep M367);
! - the gradients (grad_rad, grad, grad_ad). The radiative gradient is
!   the one of the temperature-gradient routine (legacy TPGRAD),
!     grad_rad = 3 kappa L P / (64 pi sigma G M T^4) f_T/f_P,
!   with 10^cdelrl = 3/(64 pi sigma), L = luminosity_linear Lsun and M
!   the total mass; grad is grad_rad where grad_rad <= grad_ad + 1e-6
!   (radiative_gradient_margin), and grad_ad above.
subroutine store_atmosphere_point(st, cfg, k, atm_points)
      type(integration_state), intent(inout) :: st
      type(envint_config), intent(in) :: cfg
      integer, intent(in) :: k
      type(atmosphere_structure_state), intent(inout) :: atm_points
      real(dp) :: beta, tau, depth_step, radiative_gradient, gradient

      associate (point => st%eqs%atm_point, pt => st%eqs%pt)
      beta = 1.0d0 - radiation_constant_over_3* &
           exp(ln10*(4.0d0*point%log10_temperature-point%log10_pressure))
      atm_points%atmo_log10_pressure(k) = point%log10_pressure
      atm_points%atmo_log10_temperature(k) = point%log10_temperature
      atm_points%atmo_log10_density(k) = point%log10_density
      atm_points%atmo_beta(k) = beta
      atm_points%atmo_gamma1(k) = gamma1_from_eos_derivatives(pt%dlnrho_dlnp, &
           pt%dlnrho_dlnt, pt%specific_heat_cp, point%log10_pressure, &
           point%log10_density, point%log10_temperature)
      atm_points%atmo_dlnrho_dlnt(k) = pt%dlnrho_dlnt
      atm_points%atmo_ion_fraction(1,k) = point%ion_fraction(1)
      atm_points%atmo_ion_fraction(2,k) = point%ion_fraction(2)
      atm_points%atmo_ion_fraction(3,k) = point%ion_fraction(3)
      atm_points%atmo_opacity(k) = point%opacity
      atm_points%atmo_dlnkap_dlnrho(k) = pt%dlnkap_dlnrho
      atm_points%atmo_dlnkap_dlnt(k) = pt%dlnkap_dlnt
      atm_points%atmo_specific_heat_cp(k) = pt%specific_heat_cp

      tau = exp(ln10*point%log10_tau)
      depth_step = (tau - st%prev_tau)/ &
           (((pt%density*pt%opacity)+(st%prev_density*st%prev_opacity))/2.0d0)
      st%prev_tau = tau
      st%prev_density = pt%density
      st%prev_opacity = pt%opacity
      atm_points%atmo_delta_depth(k) = depth_step

      radiative_gradient = pt%opacity*st%eqs%luminosity_linear* &
           exp(ln10*(pt%log10_pressure-cfg%log10_total_mass-4.0d0*pt%log10_temperature &
           +cfg%log10_solar_luminosity-cfg%log10_gravitational_constant+ &
           cdelrl))*st%eqs%temperature_rotation_factor/st%eqs%pressure_rotation_factor
      if (radiative_gradient-pt%adiabatic_gradient <= radiative_gradient_margin) then
         gradient = radiative_gradient
      else
         gradient = pt%adiabatic_gradient
      end if
      atm_points%atmo_gradients(1,k) = radiative_gradient
      atm_points%atmo_gradients(2,k) = gradient
      atm_points%atmo_gradients(3,k) = pt%adiabatic_gradient
      atm_points%num_atm_points = k
      end associate
end subroutine store_atmosphere_point

! ===============================================================
! Part 3: the envelope
! ===============================================================

! ---------------------------------------------------------------
! The envelope: log(m/M), log T and log r against log P, from the
! surface point (the photosphere, or a table's) at log10_radius down to
! the fitting mass cfg%log_fit_point_mass_fraction, storing env_points,
! which is then reversed to run from the centre outward. When the
! fitting mass lies above max_log_fit_point_mass_fraction there is no
! envelope to integrate, and the vertex is the surface point itself.
! Either way an envint_envelope_done event closes it. ierr after
! max_steps steps.
subroutine integrate_envelope(st, cfg, log10_radius, log10_pressure_limit, &
     env_points, events, ierr, vertex)
      type(integration_state), intent(inout) :: st
      type(envint_config), intent(in) :: cfg
      real(dp), intent(in) :: log10_radius, log10_pressure_limit
      type(envelope_structure_state), intent(inout) :: env_points
      type(envint_events), intent(inout) :: events
      integer, intent(out) :: ierr
      type(fitting_vertex), intent(inout), optional :: vertex

      ierr = 0
      if (cfg%log_fit_point_mass_fraction > max_log_fit_point_mass_fraction) then
! An envelope too thin to integrate: the vertex is the surface point.
         if (present(vertex)) then
            vertex%log10_pressure = st%eqs%atm_point%log10_pressure
            vertex%log10_radius = log10_radius
            vertex%log10_temperature = st%eqs%atm_point%log10_temperature
         end if
      else
         call start_envelope_integration(st, cfg, log10_radius, ierr, vertex)
         if (ierr /= 0) return
! The surface point is point 1.
         env_points%num_env_points = 0
         call append_envelope_point(st%eqs, env_points)
         call integrate_to_fitting_mass(st, cfg, log10_pressure_limit, env_points, &
              events, ierr, vertex)
         if (ierr /= 0) return
         call reverse_envelope_points(env_points)
      end if
      call add_envint_event(events, envint_envelope_done, &
           counts=[st%num_ok, st%num_bad], &
           values=[st%mass_left, st%y(1), st%err_sum(1), st%err_sum(2), st%err_sum(3)])
end subroutine integrate_envelope

! ---------------------------------------------------------------
! The envelope integration's starting values and numerical parameters,
! and the derivatives at the surface point: x = log P, y = (log(m/M) =
! 0, log T, log r); three equations, the envelope's tolerance and steps,
! and done when log(m/M) is within tolerance_fraction
! |log_fit_point_mass_fraction| of the fitting mass. A stored envelope
! state that belongs to this vertex is forgotten; the pressure limit
! applies only when no vertex is being fitted.
subroutine start_envelope_integration(st, cfg, log10_radius, ierr, vertex)
      type(integration_state), intent(inout) :: st
      type(envint_config), intent(in) :: cfg
      real(dp), intent(in) :: log10_radius
      integer, intent(out) :: ierr
      type(fitting_vertex), intent(inout), optional :: vertex

      st%eqs%in_atmosphere = .false.
      st%err_sum = 0.0d0
      st%num_ok = 0
      st%num_bad = 0
      st%x = st%eqs%atm_point%log10_pressure
      st%y(1) = 0.0d0
      st%y(2) = st%eqs%atm_point%log10_temperature
      st%y(3) = log10_radius
      st%num_eqs = 3
      st%tolerance = cfg%env_error_tol
      st%h_max = cfg%env_step_max
      st%h_min = cfg%env_step_min
      st%h_step = cfg%env_step_begin
      st%end_tolerance = abs(cfg%tolerance_fraction*cfg%log_fit_point_mass_fraction)
      st%state_stored = .false.
      if (present(vertex)) then
         if (vertex%stored_vertex_index == vertex%index) vertex%stored_vertex_index = 0
         st%pressure_limit_pending = .false.
      else
         st%pressure_limit_pending = .true.
      end if
! dy/dx at the start of the first step.
      call st%eqs%derivatives(st%x, st%y, st%dydx, ierr)
end subroutine start_envelope_integration

! ---------------------------------------------------------------
! Integrate the envelope in pressure down to the fitting mass, storing
! each point in env_points. Done when log(m/M) is within end_tolerance
! of the fitting mass: then the vertex is interpolated there. Until
! then, for a vertex, the state where log T first reaches
! cfg%envelope_store_log_temperature is stored (log P, log T, log r,
! log(m/M)). ierr 1, and an event, after max_steps steps.
subroutine integrate_to_fitting_mass(st, cfg, log10_pressure_limit, env_points, &
     events, ierr, vertex)
      type(integration_state), intent(inout) :: st
      type(envint_config), intent(in) :: cfg
      real(dp), intent(in) :: log10_pressure_limit
      type(envelope_structure_state), intent(inout) :: env_points
      type(envint_events), intent(inout) :: events
      integer, intent(out) :: ierr
      type(fitting_vertex), intent(inout), optional :: vertex
      integer :: step

      ierr = 0
      do step = 1, max_steps
         call limit_envelope_step(st, cfg, log10_pressure_limit)
         call take_bulirsch_stoer_step(st, cfg, events, ierr)
         if (ierr /= 0) return
         st%mass_left = cfg%log_fit_point_mass_fraction - st%y(1)
! The new point, unless the step went past the fitting mass by more than
! end_tolerance.
         if (st%mass_left <= st%end_tolerance) call append_envelope_point(st%eqs, env_points)
         if (abs(st%mass_left) <= st%end_tolerance) then
            if (present(vertex)) call interpolate_envelope_vertex(st, cfg, events, vertex)
            return
         else if (.not.st%state_stored .and. present(vertex)) then
            if (st%y(2) >= cfg%envelope_store_log_temperature) then
               st%state_stored = .true.
               vertex%stored_vertex_index = vertex%index
               vertex%stored_state(1) = st%x
               vertex%stored_state(2) = st%y(2)
               vertex%stored_state(3) = st%y(3)
               vertex%stored_state(4) = st%y(1)
            end if
         end if
         if (st%h_next < st%h_min) st%h_next = st%h_min
         st%h_step = st%h_next
      end do
! The integration did not reach the fitting mass in max_steps steps.
      call add_envint_event(events, envint_envelope_maxstp, counts=[max_steps, 0])
      ierr = 1
end subroutine integrate_to_fitting_mass

! ---------------------------------------------------------------
! Where the envelope's next step starts, its error scale, and the step
! limited so that it neither overshoots the fitting mass nor exceeds the
! largest step or the pressure limit.
subroutine limit_envelope_step(st, cfg, log10_pressure_limit)
      type(integration_state), intent(inout) :: st
      type(envint_config), intent(in) :: cfg
      real(dp), intent(in) :: log10_pressure_limit
      real(dp) :: mass_after_step

      st%x_start = st%x
      st%y_start = st%y
      call set_step_error_scale(st)
      mass_after_step = st%y(1) + st%h_step*st%dydx(1)
      if (cfg%log_fit_point_mass_fraction - st%y(1) > 0.0d0 .or. &
           cfg%log_fit_point_mass_fraction - mass_after_step > 0.0d0) then
! The integration has passed the fitting mass, or the next step would:
! aim the step at it with the linear estimate h = (log(m_fit/M) - y(1))
! / dydx(1), a step back when it has passed.
         st%h_step = (cfg%log_fit_point_mass_fraction - st%y(1))/st%dydx(1)
      end if
! Keep |h_step| within the largest step.
      if (st%h_step < 0.0d0) then
         if (st%h_step < -st%h_max) st%h_step = -st%h_max
      else
         if (st%h_step > st%h_max) st%h_step = st%h_max
      end if
! log10_pressure_limit (legacy PLIM) estimates the pressure at the end
! of the integration, from the pressure of the last model point. The
! first step that would cross it stops at it.
      if (st%pressure_limit_pending) then
         if (st%x + st%h_step > log10_pressure_limit) then
            st%h_step = log10_pressure_limit - st%x
            st%pressure_limit_pending = .false.
         end if
      end if
end subroutine limit_envelope_step

! ---------------------------------------------------------------
! The point compute_envelope_derivs last evaluated (eqs%pt), appended to
! env_points with the surface X and Z and L (Lsun), which are constant
! through the envelope. To record another quantity at each point: a
! field of point_physics that the integrand fills, a member of
! envelope_structure_state, a line here and one in
! reverse_envelope_points (core/README.md, "Record a quantity at each
! envelope point").
subroutine append_envelope_point(eqs, env_points)
      type(envelope_equations), intent(in) :: eqs
      type(envelope_structure_state), intent(inout) :: env_points
      integer :: k

      k = env_points%num_env_points + 1
      env_points%num_env_points = k
      associate (pt => eqs%pt)
      env_points%env_log10_density(k) = pt%log10_density
      env_points%env_log10_pressure(k) = pt%log10_pressure
      env_points%env_log10_radius(k) = pt%log10_radius
      env_points%env_log10_mass(k) = pt%log10_mass_fraction
      env_points%env_log10_temperature(k) = pt%log10_temperature
      env_points%env_hydrogen_fraction(k) = eqs%hydrogen_fraction
      env_points%env_metal_fraction(k) = eqs%metal_fraction
      env_points%env_convective_flag(k) = pt%convective_velocity > 0.0d0
! The gradients (grad_rad, grad_ad, grad), the convective velocity and
! beta.
      env_points%env_gradients(1,k) = pt%radiative_gradient
      env_points%env_gradients(2,k) = pt%adiabatic_gradient
      env_points%env_gradients(3,k) = pt%temperature_gradient
      env_points%env_convective_velocity(k) = pt%convective_velocity
      env_points%env_beta(k) = pt%beta
! Gamma1 and c_P, the ionization fractions, the opacity and its
! derivatives, L (in Lsun) and (d ln rho/d ln T)_P.
      env_points%env_gamma1(k) = gamma1_from_eos_derivatives(pt%dlnrho_dlnp, &
           pt%dlnrho_dlnt, pt%specific_heat_cp, pt%log10_pressure, &
           pt%log10_density, pt%log10_temperature)
      env_points%env_specific_heat_cp(k) = pt%specific_heat_cp
      env_points%env_ion_fraction(1,k) = pt%ion_fraction(1)
      env_points%env_ion_fraction(2,k) = pt%ion_fraction(2)
      env_points%env_ion_fraction(3,k) = pt%ion_fraction(3)
      env_points%env_opacity(k) = pt%opacity
      env_points%env_dlnkap_dlnrho(k) = pt%dlnkap_dlnrho
      env_points%env_dlnkap_dlnt(k) = pt%dlnkap_dlnt
      env_points%env_luminosity(k) = eqs%luminosity_linear
      env_points%env_dlnrho_dlnt(k) = pt%dlnrho_dlnt
      end associate
end subroutine append_envelope_point

! ---------------------------------------------------------------
! The fitting point of the vertex (log P, log r, log T), interpolated
! linearly in log(m/M) to the fitting mass between the start and the end
! of the last step:
!   v = v_end + w (v_start - v_end),
!   w = (log(m_fit/M) - y(1)) / (y_start(1) - y(1)),
! and its event.
subroutine interpolate_envelope_vertex(st, cfg, events, vertex)
      type(integration_state), intent(in) :: st
      type(envint_config), intent(in) :: cfg
      type(envint_events), intent(inout) :: events
      type(fitting_vertex), intent(inout) :: vertex
      real(dp) :: weight

      weight = st%mass_left/(st%y_start(1)-st%y(1))
      vertex%log10_pressure = st%x + weight*(st%x_start - st%x)
      vertex%log10_radius = st%y(3) + weight*(st%y_start(3) - st%y(3))
      vertex%log10_temperature = st%y(2) + weight*(st%y_start(2) - st%y(2))
      call add_envint_event(events, envint_vertex, &
           values=[vertex%log10_pressure, vertex%log10_temperature, &
           vertex%log10_radius, cfg%log_fit_point_mass_fraction])
end subroutine interpolate_envelope_vertex

! ---------------------------------------------------------------
! Reverse env_points from surface-inward to centre-outward.
subroutine reverse_envelope_points(env_points)
      type(envelope_structure_state), intent(inout) :: env_points
      associate (n => env_points%num_env_points)
         env_points%env_log10_density(1:n) = env_points%env_log10_density(n:1:-1)
         env_points%env_log10_pressure(1:n) = env_points%env_log10_pressure(n:1:-1)
         env_points%env_log10_radius(1:n) = env_points%env_log10_radius(n:1:-1)
         env_points%env_log10_mass(1:n) = env_points%env_log10_mass(n:1:-1)
         env_points%env_log10_temperature(1:n) = env_points%env_log10_temperature(n:1:-1)
         env_points%env_hydrogen_fraction(1:n) = env_points%env_hydrogen_fraction(n:1:-1)
         env_points%env_metal_fraction(1:n) = env_points%env_metal_fraction(n:1:-1)
         env_points%env_convective_flag(1:n) = env_points%env_convective_flag(n:1:-1)
         env_points%env_opacity(1:n) = env_points%env_opacity(n:1:-1)
         env_points%env_dlnkap_dlnrho(1:n) = env_points%env_dlnkap_dlnrho(n:1:-1)
         env_points%env_dlnkap_dlnt(1:n) = env_points%env_dlnkap_dlnt(n:1:-1)
         env_points%env_luminosity(1:n) = env_points%env_luminosity(n:1:-1)
         env_points%env_dlnrho_dlnt(1:n) = env_points%env_dlnrho_dlnt(n:1:-1)
         env_points%env_gradients(:,1:n) = env_points%env_gradients(:,n:1:-1)
         env_points%env_convective_velocity(1:n) = env_points%env_convective_velocity(n:1:-1)
         env_points%env_beta(1:n) = env_points%env_beta(n:1:-1)
         env_points%env_gamma1(1:n) = env_points%env_gamma1(n:1:-1)
         env_points%env_specific_heat_cp(1:n) = env_points%env_specific_heat_cp(n:1:-1)
         env_points%env_ion_fraction(:,1:n) = env_points%env_ion_fraction(:,n:1:-1)
      end associate
end subroutine reverse_envelope_points

! ===============================================================
! The Bulirsch-Stoer step both integrations take, and the events
! ===============================================================

! ---------------------------------------------------------------
! y_scale (legacy YSCAL), the scale each variable's error is measured
! against: |y| + |h dy/dx| + step_scale_floor.
subroutine set_step_error_scale(st)
      type(integration_state), intent(inout) :: st
      integer :: i
      do i = 1, st%num_eqs
         st%y_scale(i) = abs(st%y(i)) + abs(st%h_step*st%dydx(i))+step_scale_floor
      end do
end subroutine set_step_error_scale

! ---------------------------------------------------------------
! One Bulirsch-Stoer step of the integration in progress, from x to
! x + h_step (h_step is the step tried, h_did the step taken and h_next
! the suggested next step), then the count of the step (taken at the
! size tried, or smaller), the summed error estimates, and dy/dx where
! the next step starts, which also leaves the new point in eqs%pt (and
! eqs%atm_point) for the stores. ierr as check_bsstep_status sets it,
! or the integrands' code at the new point.
subroutine take_bulirsch_stoer_step(st, cfg, events, ierr)
      type(integration_state), intent(inout) :: st
      type(envint_config), intent(in) :: cfg
      type(envint_events), intent(inout) :: events
      integer, intent(out) :: ierr
! bsstep's status, and the integrands' code when they failed inside it
      integer :: bs_status, derivatives_status
      integer :: i

      call bsstep(st%eqs, st%y, st%dydx, st%num_eqs, st%x, st%h_step, st%tolerance, &
           st%y_scale, st%h_did, st%h_next, st%step_err, cfg%max_stage_index, &
           cfg%extrap_order, st%underflow_reference, bs_status, &
           derivatives_ierr=derivatives_status)
      call check_bsstep_status(bs_status, derivatives_status, st, cfg, events, ierr)
      if (ierr /= 0) return
      if (st%h_did == st%h_step) then
         st%num_ok = st%num_ok + 1
      else
         st%num_bad = st%num_bad + 1
      end if
      do i = 1, st%num_eqs
         st%err_sum(i) = st%err_sum(i) + st%step_err(i)
      end do
      call st%eqs%derivatives(st%x, st%y, st%dydx, ierr)
end subroutine take_bulirsch_stoer_step

! ---------------------------------------------------------------
! bsstep's status after a step of the atmosphere or the envelope, into
! ierr and an event for the driver's message: a limit out of range is a
! configuration error (numerics_bad_input); a failure of the integrands
! (the EOS or the opacity at one of the step's midpoint substeps) returns
! their code, an error, as the same failure at the step's end does (M132);
! a step that underflowed (the solution diverged) returns
! envelope_termination. ierr is 0 when the step succeeded.
subroutine check_bsstep_status(bs_status, derivatives_status, st, cfg, events, ierr)
      integer, intent(in) :: bs_status, derivatives_status
      type(integration_state), intent(in) :: st
      type(envint_config), intent(in) :: cfg
      type(envint_events), intent(inout) :: events
      integer, intent(out) :: ierr

      ierr = 0
      if (bs_status == numerics_bad_input) then
         call add_envint_event(events, envint_bsstep_limits, &
              counts=[cfg%max_stage_index, cfg%extrap_order])
         ierr = bs_status
      else if (bs_status == numerics_derivatives_failed) then
         call add_envint_event(events, envint_bsstep_failed, &
              counts=[bs_status, derivatives_status], &
              values=[st%x, real(st%num_eqs, dp)])
         ierr = derivatives_status
      else if (bs_status /= 0) then
         call add_envint_event(events, envint_bsstep_failed, counts=[bs_status, 0], &
              values=[st%x, real(st%num_eqs, dp)])
         ierr = envelope_termination
      end if
end subroutine check_bsstep_status

! ---------------------------------------------------------------
! Append an event of the given kind, with its counts and values (the
! slots not given stay zero), to events.
subroutine add_envint_event(events, kind, counts, values)
      type(envint_events), intent(inout) :: events
      integer, intent(in) :: kind
      integer, intent(in), optional :: counts(:)
      real(dp), intent(in), optional :: values(:)
      type(envint_event), allocatable :: grown(:)
      integer :: n

      if (.not. allocated(events%event)) allocate(events%event(4))
      n = events%num_events + 1
      if (n > size(events%event)) then
         allocate(grown(2*size(events%event)))
         grown(1:n-1) = events%event(1:n-1)
         call move_alloc(grown, events%event)
      end if
      events%event(n) = envint_event()
      events%event(n)%kind = kind
      if (present(counts)) events%event(n)%counts(1:size(counts)) = counts
      if (present(values)) events%event(n)%values(1:size(values)) = values
      events%num_events = n
end subroutine add_envint_event

end module envelope_integration_kernel

! Sources
! - The F77's ENVINT (envint.f), with its integrands QATM and QENV
!   (envelope_integrands) and its table lookup ATMPL (atm_lib). The
!   fitting-mass tests: ENVINT skipped the envelope above -1e-12 and
!   reversed the points below it (.LT.), so at exactly -1e-12 it
!   integrated without reversing; fixed 2026-10. The f_P of the
!   atmosphere's first dydx was missing until 2026-10, as in envint.f, so
!   a rotating star's first step started from dydx = 1.
! - The integration was atm_lib's until the physics-purity pass of 2026
!   (docs/ROADMAP.md), which made it this star-blind kernel under the
!   driver atm_get; it prints nothing since round 2 (the events).
! - Round 3 (2026-11): the host's contained steps, which shared about 70
!   locals, became these module procedures on integration_state; the two
!   point stores became append_envelope_point, and the Gamma1 formula
!   envelope_integrands' gamma1_from_eos_derivatives. ENVINT's sum of
!   dr/v_conv over the envelope's convective points, which nothing read,
!   went (core/diagnostics/turnover_timescale.f90 computes the turnover
!   time from the whole model). ENVINT had a blanket SAVE, which went in
!   2026-08; since M322 the state starts at zero.

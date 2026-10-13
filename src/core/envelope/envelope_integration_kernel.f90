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
!      surface down to the fitting mass, log(m/M) = senv.
! On the way it records each atmosphere point in atmo_struct and each
! envelope point in env_struct (those module records are its product).
! The physics of each point comes from pt_scr (point_scratch_lib), where
! the integrands of envelope_integrands leave it.
!
! The kernel is star-blind: its code uses neither star_info_lib nor
! star%, which src/tools/check_boundaries.py checks (STAR_BLIND_FILES).
! Its settings arrive in an envint_config that the driver, atm_get in
! envelope_integration, fills from star%job and star%ctrl. Its two
! effects on the star, the photospheric pressure (star%pphot) and the
! switch to a gray atmosphere when the Allard table lookup fails, come
! back as outputs that the driver applies.
module envelope_integration_kernel
      use envelope_integrands, only: compute_atmosphere_derivs, compute_envelope_derivs
      use numerics_lib, only: ode_system
      implicit none

! The kernel's settings: the atmosphere choice; the Bulirsch-Stoer step
! sizes (begin, min, max, in log tau and in log P) and the error
! tolerances of the two integrations; the starting density of the
! atmosphere integration (atm_step_initial, legacy ATMD0, g/cm^3); the
! HSRA offset atm_hras; the fitting mass senv = log10(m_fit/M); tenv,
! the log T at which the envelope state is stored for the surface
! boundary; log10 Lsun; log10 of the total mass in g, the origin of the
! envelope's mass coordinate log(m/M) (star%stotal); and calc_envelope
! (false: the photosphere only).
      type, public :: envint_config
            integer :: atm_choice
            double precision :: atm_step_begin, atm_step_min, atm_step_max
            double precision :: env_step_begin, env_step_min, env_step_max
            double precision :: atm_error_tol, env_error_tol, atm_step_initial
            double precision :: atm_hras, senv, tenv, log10_solar_luminosity
            double precision :: log10_total_mass
            logical :: calc_envelope
      end type envint_config

! The largest fitting mass senv = log10(m_fit/M) the code allows
! (rebuild_envelope clamps senv to it): an envelope of at least
! 1 - 10^senv = ln(10) 1e-12, about 2.3e-12 of the star.
! integrate_envelope skips the envelope integration when senv lies
! above it.
      double precision, parameter, public :: senv_thin_envelope = -1.0d-12

! The atmosphere (compute_atmosphere_derivs: log P against log tau) or
! the envelope (compute_envelope_derivs: log(m/M), log T and log r
! against log P) as an ode_system for bsstep. It carries the values the
! integrands read, and the three they update (log10_radius, log10_teff,
! saha_state).
      type, extends(ode_system) :: envelope_equations
            logical :: atmosphere = .true.
            double precision :: luminosity_linear = 0.0d0, &
                 pressure_rotation_factor = 0.0d0, &
                 temperature_rotation_factor = 0.0d0, log10_gravity = 0.0d0
            logical :: in_atmosphere = .false., want_derivatives = .false., &
                 conductive_opacity_flag = .false.
            double precision :: log10_radius = 0.0d0, log10_teff = 0.0d0, &
                 hydrogen_fraction = 0.0d0, metal_fraction = 0.0d0
            integer :: saha_state = 0
! the T(tau) relation of the atmosphere integration (cfg%atm_choice,
! gray after a failed Allard lookup), its HSRA offset, and log10 of the
! total mass in g
            integer :: atm_choice = 0
            double precision :: atm_hras = 0.0d0, log10_total_mass = 0.0d0
      contains
            procedure :: derivatives => equation_derivatives
      end type envelope_equations

contains

! dy/dx of the atmosphere or of the envelope (as self%atmosphere says)
! at (x, y).
subroutine equation_derivatives(self, x, y, dydx, ierr)
      class(envelope_equations), intent(inout) :: self
      double precision, intent(in) :: x, y(:)
      double precision, intent(out) :: dydx(:)
      integer, intent(out) :: ierr
      if (self%atmosphere) then
         call compute_atmosphere_derivs(x, y, dydx, self%luminosity_linear, &
              self%pressure_rotation_factor, self%temperature_rotation_factor, &
              self%log10_gravity, self%in_atmosphere, self%want_derivatives, &
              self%conductive_opacity_flag, self%log10_radius, self%log10_teff, &
              self%hydrogen_fraction, self%metal_fraction, self%saha_state, &
              self%atm_choice, self%atm_hras, ierr)
      else
         call compute_envelope_derivs(x, y, dydx, self%luminosity_linear, &
              self%pressure_rotation_factor, self%temperature_rotation_factor, &
              self%log10_gravity, self%in_atmosphere, self%want_derivatives, &
              self%conductive_opacity_flag, self%log10_radius, self%log10_teff, &
              self%hydrogen_fraction, self%metal_fraction, self%saha_state, &
              self%log10_total_mass, ierr)
      end if
end subroutine equation_derivatives

! ---------------------------------------------------------------
! integrate_envelope_atmosphere (the body of legacy envint)
! ---------------------------------------------------------------
! Integrates the atmosphere and the envelope from the surface down to
! the fitting point (the three parts of the module header), each by
! Bulirsch-Stoer steps (bsstep), at most maxstp of them.
! Inputs: the settings cfg; L in Lsun (luminosity_linear); the
! rotational distortion factors f_P and f_T (1 without rotation);
! log g; log M; the surface X and Z; log Teff and the surface log R;
! log10_pressure_limit, a guess at the pressure of the fitting point
! that the first envelope step to cross it stops at (when the boundary
! is not being saved); and, with save_boundary_flag, the vertex of the
! surface-boundary triangle being computed (vertex_index).
! Outputs: atmo_struct; env_struct, from the centre outward; the
! photospheric log P (log10_photo_pressure, -1d99 unless the
! atmosphere finished); switched_to_gray; with save_boundary_flag, the
! vertex (vtx_logp, vtx_logr, vtx_logt) at the fitting mass and the
! state where T first reaches tenv (stored_envelope_state: log P,
! log T, log R, log(m/M); stored_vertex_index); ierr.
subroutine integrate_envelope_atmosphere(cfg, switched_to_gray, &
     log10_photo_pressure, &
     luminosity_linear, pressure_rotation_factor, &
     temperature_rotation_factor, log10_gravity, log10_star_mass, &
     vertex_index, print_flag, save_boundary_flag, log10_pressure_limit, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     stored_envelope_state, stored_vertex_index, &
     saha_state, vtx_logp, vtx_logr, vtx_logt, &
     ierr)

      use eos_lib
      use kap_lib
      use point_scratch_lib
      use atmosphere_structure
      use envelope_structure
      use atm_lib, only: atm_get_surface_pt, ttau_start_log10_temperature, &
           ttau_photosphere_x_limit
      use io_units
      use physical_constants
      use integration_parameters
      use numerics_lib
      use run_log, only: solver_diagnostics
      implicit none

      type(envint_config), intent(inout) :: cfg
      logical, intent(out) :: switched_to_gray
      double precision, intent(out) :: log10_photo_pressure
      double precision, intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, &
           log10_gravity
      double precision, intent(in) :: log10_star_mass
      integer, intent(in) :: vertex_index
! print_flag: also write the fitted vertex to the run log (with
! save_boundary_flag and the solver diagnostics on).
      logical, intent(in) :: print_flag
      logical, intent(in) :: save_boundary_flag
      double precision, intent(in) :: log10_pressure_limit
      double precision, intent(inout) :: log10_radius
      double precision, intent(inout) :: log10_teff
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      double precision, intent(inout) :: stored_envelope_state(4)
      integer, intent(inout) :: stored_vertex_index
      integer, intent(inout) :: saha_state
      double precision, intent(inout) :: vtx_logp(3), vtx_logr(3), vtx_logt(3)
! ierr: nonzero when the table lookup, the EOS, the opacity or either
! integration fails; the diagnostic is printed where the failure is.
      integer, intent(out) :: ierr

! The most steps either integration may take (raised from 200 to give
! the atmosphere integration room; DBG).
      integer, parameter :: maxstp = 2000
      logical :: tabulated_bc
      double precision, parameter :: tiny = 1.0d-30

      double precision :: ion_fraction(3)
! --- locals ---
      logical :: want_derivatives, in_atmosphere, conductive_opacity_flag
      logical :: allard_lookup_failed
      integer :: jj, jerr
      double precision :: log10_temperature, temperature, log10_pressure, &
           pressure, log10_density, density
      double precision :: atm_density_guess
      double precision :: beta, beta_inverse, beta14, specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision :: dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
           adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnt_dp, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      double precision :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
! The results of eos_get and kap_get. The EOS scalars above mirror the
! slots of eos_res: they are shared by the contained routines and
! survive the atm_retry loop, so find_atmosphere_start copies them into
! eos_res before the call and back after it.
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results)
! The independent variable: log tau in the atmosphere, log P in the
! envelope.
      double precision :: indep_var
! The equations bsstep integrates (see load_equations, unload_equations).
      type(envelope_equations) :: equations
      double precision :: chi_rho, chi_t, specific_heat_cv, gamma1
      double precision :: prev_tau, prev_opacity, &
           prev_density, delta_tau_step, &
           pulse_radiative_gradient, pulse_gradient, opacity_now, &
           density_now, tau_now
      integer :: num_eqs, num_ok, num_bad, step_index, i
      double precision :: tolerance, x_limit, h_max, h_min, h_step, &
           step_tolerance, h_did, h_next
      double precision :: y(3), dydx(3), y_scale(3), y_start(3), &
           err_sum(3), step_err(3)
      double precision :: initial_mass_coord
      logical :: store_flag_set, pressure_limit_test_flag
      logical :: cz_in_envelope
      double precision :: mass_diff_remaining, interp_weight
      double precision :: mass_next_step
      double precision :: x_start, taucz_env_accum, delta_radius_cz

      ierr = 0
      jerr = 0
      switched_to_gray = .false.
! -1d99 means "no photosphere": it stays unless the atmosphere
! finishes, so that the driver does not overwrite star%pphot after an
! early error.
      log10_photo_pressure = -1.0d99

      call prepare_surface_boundary
      if (ierr /= 0) return
      call integrate_atmosphere
      if (ierr /= 0) return
! The photosphere (or the table's surface point) is reached, and its
! pressure is pphot. Stop here when only pphot is wanted; otherwise
! integrate the envelope in pressure from there down to the fitting
! mass.
      log10_photo_pressure = atm_point%log10_pressure
      if (.not.cfg%calc_envelope) return

      call integrate_envelope
      if (ierr /= 0) return
! The convective turnover time is not computed here:
! core/diagnostics/turnover_timescale.f90 computes it from the whole
! model.
      return

contains

! ---------------------------------------------------------------
! The surface boundary: the EOS flags of the atmosphere, then, for a
! tabulated atmosphere (atm_choice 3 Kurucz, 4 Allard, 5
! Castelli/Kurucz), the surface log P and log T from the table. When
! the point lies outside the Allard tables, falls back to the Eddington
! gray atmosphere (atm_choice 0), which integrate_atmosphere then
! integrates. Sets tabulated_bc; ierr when a table interpolation fails.
subroutine prepare_surface_boundary
! The EOS derivatives are always requested: the integrands store them
! in pt_scr for atmo_struct and env_struct.
      want_derivatives = .true.
      in_atmosphere = .true.
      conductive_opacity_flag = .false.

! Interpolate the surface log P and log T in the table (legacy ATMPL).
      tabulated_bc = .false.
      if (cfg%atm_choice == 3 .or. cfg%atm_choice == 4 .or. cfg%atm_choice == 5) then
         call atm_get_surface_pt(cfg%atm_choice, log10_teff, log10_gravity, &
              atm_point%log10_pressure, atm_point%log10_temperature, &
              allard_lookup_failed, jerr)
         if (jerr /= 0) then
            ierr = jerr
            return
         end if
         if(allard_lookup_failed) then
            cfg%atm_choice=0
            switched_to_gray = .true.
! Teff is above the Allard tables' maximum, or log g outside their
! range: switch to the Eddington gray atmosphere (atm_choice 0,
! legacy KTTAU = 0).
            write(*,'(1x,a)') 'note: Allard table lookup out of range;' &
                 // ' switching to a gray atmosphere'
            write(run_log_unit,'(1x,a)') 'note: Allard table lookup' &
                 // ' out of range; switching to a gray atmosphere'
         else
            tabulated_bc = .true.
         endif
      endif
end subroutine prepare_surface_boundary

! ---------------------------------------------------------------
! The atmosphere from a T-tau relation (atm_choice 0 Eddington gray,
! 1 Krishna Swamy 1966, 2 HSRA). It finds the starting point near
! tau = 0 (T from the relation, P from a guessed density, tau from P),
! then integrates log P in log tau down to the photosphere (tau = 2/3;
! 0.312 for Krishna Swamy), storing atmo_struct on the way. When the
! starting tau lies past the photosphere, the guessed density is
! divided by 10 and the start found again (atm_retry). Skipped when a
! table gave the surface; ierr after maxstp steps.
subroutine integrate_atmosphere
      use math_lib
! No atmosphere points until the integration stores them, so that a
! table boundary leaves atmo_struct empty rather than holding the
! previous integration's points.
      atmo_struct%num_atm_points = 0
      if (.not. tabulated_bc) then
! The starting temperature: the T-tau relation near tau = 0 (for HSRA,
! the fit at its tau = 2/3 anchor).
      err_sum(1) = 0.0d0
      log10_temperature = ttau_start_log10_temperature(log10_teff, &
           cfg%atm_choice, cfg%atm_hras)
! Very occasionally the starting optical depth lies past the
! photosphere (indep_var > x_limit); then the starting density is
! divided by 10 and the start found again.
      atm_density_guess = cfg%atm_step_initial
      atm_retry: do
      call find_atmosphere_start()
      if (ierr /= 0) return

! Integrate log P in log tau from this starting tau to the
! photosphere: one equation, the atmosphere's error tolerance.
      num_eqs = 1
      tolerance = cfg%atm_error_tol
! The end point, log tau of the photosphere for the chosen relation.
      x_limit = ttau_photosphere_x_limit(cfg%atm_choice)

      if  (indep_var .gt. x_limit) then ! The start lies past the end point:
         write(run_log_unit,*)"ENVINT: X0>XLIM, X0,XLIM:",indep_var,x_limit
       write(run_log_unit,*)"ENVINT: get new X0 by dividing ATMD0 by 10"
       atm_density_guess = atm_density_guess / 10d0   ! divide the starting density by 10
       cycle atm_retry        ! and find the start again.
      endif

      call integrate_to_photosphere()
      if (ierr /= 0) return
      exit atm_retry
      end do atm_retry
      end if
end subroutine integrate_atmosphere

! ---------------------------------------------------------------
! The starting point of the atmosphere integration. With the starting
! T and the guessed density rho0 (atm_density_guess), an ideal gas of
! mean molecular weight 1 plus radiation gives
!   P = (R rho0 + (a/3) T^3) T.
! The EOS and opacity at (P, T) give kappa, and with kappa and g taken
! constant above this point, tau = kappa P / g:
!   indep_var = log tau = log P - log g + log kappa,
!   dydx(1) = d log P / d log tau = g_eff tau / (kappa P),
! with g_eff = g f_P as in compute_atmosphere_derivs, i.e. f_P up to
! rounding. (The f_P was missing until 2026-10, as in envint.f, so a
! rotating star's first step started from dydx = 1.)
subroutine find_atmosphere_start()
      use math_lib
      temperature = exp(ln10*log10_temperature)
! P = (R rho0 + (a/3) T^3) T.
      log10_pressure = log10((gas_constant*atm_density_guess + &
           radiation_constant_over_3*temperature**3)*temperature)

! The EOS at (P, T). eos_res is inout in eos_get, so every slot is
! copied in from its scalar before the call and back after it: a slot
! eos_get leaves alone keeps the scalar's value, and the scalars carry
! the values of the previous try into the next.
      eos_res%temperature = temperature
      eos_res%pressure = pressure
      eos_res%log10_density = log10_density
      eos_res%density = density
      eos_res%beta = beta
      eos_res%beta_inverse = beta_inverse
      eos_res%beta14 = beta14
      eos_res%ion_fraction = ion_fraction
      eos_res%specific_gas_constant = specific_gas_constant
      eos_res%ion_mean_weight_inverse = ion_mean_weight_inverse
      eos_res%electron_mean_weight_inverse = electron_mean_weight_inverse
      eos_res%electron_degeneracy_parameter = electron_degeneracy_parameter
      eos_res%dlnrho_dlnt = dlnrho_dlnt
      eos_res%dlnrho_dlnp = dlnrho_dlnp
      eos_res%specific_heat_cp = specific_heat_cp
      eos_res%adiabatic_gradient = adiabatic_gradient
      eos_res%dlnrho_dlnt_dt = dlnrho_dlnt_dt
      eos_res%dlnrho_dlnt_dp = dlnrho_dlnt_dp
      eos_res%adiabatic_gradient_dt = adiabatic_gradient_dt
      eos_res%adiabatic_gradient_dp = adiabatic_gradient_dp
      eos_res%specific_heat_cp_dt = specific_heat_cp_dt
      eos_res%specific_heat_cp_dp = specific_heat_cp_dp
      call eos_get(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, eos_res, want_derivatives, in_atmosphere, &
           saha_state, ierr=jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if
      temperature = eos_res%temperature
      pressure = eos_res%pressure
      log10_density = eos_res%log10_density
      density = eos_res%density
      beta = eos_res%beta
      beta_inverse = eos_res%beta_inverse
      beta14 = eos_res%beta14
      ion_fraction = eos_res%ion_fraction
      specific_gas_constant = eos_res%specific_gas_constant
      ion_mean_weight_inverse = eos_res%ion_mean_weight_inverse
      electron_mean_weight_inverse = eos_res%electron_mean_weight_inverse
      electron_degeneracy_parameter = eos_res%electron_degeneracy_parameter
      dlnrho_dlnt = eos_res%dlnrho_dlnt
      dlnrho_dlnp = eos_res%dlnrho_dlnp
      specific_heat_cp = eos_res%specific_heat_cp
      adiabatic_gradient = eos_res%adiabatic_gradient
      dlnrho_dlnt_dt = eos_res%dlnrho_dlnt_dt
      dlnrho_dlnt_dp = eos_res%dlnrho_dlnt_dp
      adiabatic_gradient_dt = eos_res%adiabatic_gradient_dt
      adiabatic_gradient_dp = eos_res%adiabatic_gradient_dp
      specific_heat_cp_dt = eos_res%specific_heat_cp_dt
      specific_heat_cp_dp = eos_res%specific_heat_cp_dp
! The opacity at the density the EOS returned (DBG 12/95).
      call kap_get(log10_density, log10_temperature, hydrogen_fraction, &
           metal_fraction, kap_res, ion_fraction, ierr=jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if
      opacity = kap_res(i_kap)
      log10_opacity = kap_res(i_log10_kap)
      dlnkap_dlnrho = kap_res(i_dlnkap_dlnrho)
      dlnkap_dlnt = kap_res(i_dlnkap_dlnt)
! The starting optical depth, tau = kappa P / g, and d log P / d log tau.
      indep_var = log10_pressure - log10_gravity + log10(opacity)
      y(1) = log10_pressure
      dydx(1) = exp(ln10*(log10_gravity+indep_var-log10_opacity-log10_pressure))* &
           pressure_rotation_factor
end subroutine find_atmosphere_start

! ---------------------------------------------------------------
! Integrate dP/dtau from the starting optical depth to the photosphere,
! storing each point in atmo_struct; ierr after maxstp steps.
subroutine integrate_to_photosphere()
      h_max = cfg%atm_step_max
      h_min = cfg%atm_step_min
      h_step = cfg%atm_step_begin
      step_tolerance = tolerance_fraction
      num_ok = 0
      num_bad = 0

! The integration loop. A step taken at the size tried counts as
! accepted (num_ok), a smaller one as rejected (num_bad).
      do step_index = 1,maxstp
       call limit_atmosphere_step()
       call take_atmosphere_step()
       if (ierr /= 0) return
       call store_atmosphere_point()
       if(h_did.eq.h_step) then
          num_ok = num_ok + 1
       else
          num_bad = num_bad + 1
       endif
       call store_atmosphere_gradients()


! Done when log tau is within step_tolerance of the photosphere.
       if(dabs(indep_var - x_limit).le.step_tolerance) then
          if (solver_diagnostics()) write(run_log_unit,35)num_ok,num_bad,err_sum(1)
   35       format(1X,'ATMOSPHERE INTEGRATION COMPLETE',1X, &
                 'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
                 I5/5X,'MAXIMUM RELATIVE ERROR IN P ',1PE22.13)
          exit
       endif
       if(h_next.lt.h_min) h_next = h_min
       h_step = h_next
      end do
      if (step_index .gt. maxstp) then
! The integration did not reach the photosphere in maxstp steps:
! report it and return an error.
      write(terminal_unit,50)
      write(run_log_unit,50)
   50 format(5X,'ATMOSPHERE INTEGRATION FAILED AFTER MAXSTP',1X, &
             'INTEGRATIONS.I QUIT.')
      jerr = 1
      ierr = jerr
      return
      end if
end subroutine integrate_to_photosphere

! ---------------------------------------------------------------
! The error scale of the step, and the step limited to the maximum size
! and to the end of the integration.
subroutine limit_atmosphere_step()
! y_scale (legacy YSCAL), the scale each variable's error is measured
! against: |y| + |h dy/dx| + tiny.
       do i = 1,num_eqs
          y_scale(i) = dabs(y(i)) + dabs(h_step*dydx(i))+tiny
       end do
! Keep the step from going past the end of the integration, and within
! the largest step.
       if((indep_var-x_limit)*(indep_var+h_step-x_limit).lt.0.0d0) h_step = x_limit - indep_var
       if(h_step.gt.h_max) h_step = h_max
end subroutine limit_atmosphere_step

! ---------------------------------------------------------------
! One Bulirsch-Stoer step of the atmosphere, then dP/dtau at its end.
subroutine take_atmosphere_step()
! One step from log tau to log tau + h_step: h_step is the step tried,
! h_did the step taken and h_next the predicted next step.
       call load_equations(.true.)
       call bsstep(equations,y,dydx,num_eqs,indep_var,h_step,tolerance,y_scale,h_did, &
            h_next,step_err,max_stage_index,extrap_order,hydrogen_fraction,jerr)
       call unload_equations()
! A bsstep failure (the solution diverged; bsstep has printed the
! diagnostic) returns the numerics-termination code.
       if (jerr /= 0) then
          ierr = numerics_termination
          return
       end if
! Add up the error estimate, and find d log P / d log tau at the start
! of the next step.
       err_sum(1) = err_sum(1) + step_err(1)
       call compute_atmosphere_derivs(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
            temperature_rotation_factor,log10_gravity,in_atmosphere, &
            want_derivatives,conductive_opacity_flag,log10_radius, &
            log10_teff,hydrogen_fraction,metal_fraction,saha_state, &
            cfg%atm_choice, cfg%atm_hras, jerr)
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
end subroutine take_atmosphere_step

! ---------------------------------------------------------------
! The thermodynamic state of the new atmosphere point, as the last call
! of compute_atmosphere_derivs left it in atm_point and pt_scr, into
! atmo_struct.
subroutine store_atmosphere_point()
      use math_lib
! The gas-pressure fraction beta = 1 - (a/3) T^4 / P; from the EOS
! derivatives, chi_rho = 1/(d ln rho/d ln P)_T,
! chi_T = -chi_rho (d ln rho/d ln T)_P,
! c_V = c_P - (P/(rho T)) chi_T^2/chi_rho and Gamma1 = chi_rho c_P/c_V.
       beta = 1.0d0 - radiation_constant_over_3*exp(ln10*(4.0d0*atm_point%log10_temperature-atm_point%log10_pressure))
       chi_rho = 1.0d0/pt_scr%qqdp
       chi_t = -chi_rho*pt_scr%qqdt
       specific_heat_cv = pt_scr%qqcp - exp(ln10*(atm_point%log10_pressure-atm_point%log10_density-atm_point%log10_temperature))*chi_t**2/chi_rho
       gamma1 = chi_rho*pt_scr%qqcp/specific_heat_cv
! Store the point in atmo_struct.
       atmo_struct%atmo_log10_pressure(step_index) = atm_point%log10_pressure
       atmo_struct%atmo_log10_temperature(step_index) = atm_point%log10_temperature
       atmo_struct%atmo_log10_density(step_index) = atm_point%log10_density
       atmo_struct%atmo_beta(step_index) = beta
       atmo_struct%atmo_gamma1(step_index) = gamma1
       atmo_struct%atmo_dlnrho_dlnt(step_index) = pt_scr%qqdt
       atmo_struct%atmo_ion_fraction(1,step_index) = atm_point%ion_fraction(1)
       atmo_struct%atmo_ion_fraction(2,step_index) = atm_point%ion_fraction(2)
       atmo_struct%atmo_ion_fraction(3,step_index) = atm_point%ion_fraction(3)
       atmo_struct%atmo_opacity(step_index) = atm_point%opacity
       atmo_struct%atmo_dlnkap_dlnrho(step_index) = pt_scr%qqod
       atmo_struct%atmo_dlnkap_dlnt(step_index) = pt_scr%qqot
       atmo_struct%atmo_specific_heat_cp(step_index) = pt_scr%qqcp
end subroutine store_atmosphere_point

! ---------------------------------------------------------------
! The depth step from the previous point and the temperature gradients
! of the new atmosphere point, into atmo_struct.
subroutine store_atmosphere_gradients()
      use math_lib
! The geometric depth between this tau point and the previous one,
!   dz = (tau - tau_prev) / [(rho kappa + rho_prev kappa_prev)/2]
! (Cox & Giuli 1968, Principles of Stellar Structure, p. 590).
! prev_tau, prev_density and prev_opacity carry over from step to step
! and are not set before the first step. They are not SAVEd, so the
! build's -finit-local-zero makes them zero on entry: the first point's
! depth is measured from tau = 0 with rho kappa = 0 there.
       opacity_now = pt_scr%qo
       density_now = exp(ln10*pt_scr%qdl)
       tau_now = exp(ln10*atm_point%log10_tau)
       delta_tau_step =  (tau_now - prev_tau)/(((density_now*opacity_now)+(prev_density*prev_opacity))/2)
       prev_opacity = opacity_now
! The radiative gradient, as at the start of the temperature-gradient
! routine (legacy TPGRAD),
!   grad_rad = 3 kappa L P / (64 pi sigma G M T^4) f_T/f_P,
! with 10^cdelrl = 3/(64 pi sigma), L = luminosity_linear Lsun and M the
! total mass. The gradient is grad_rad where grad_rad <= grad_ad + 1e-6,
! and grad_ad above.
       pulse_radiative_gradient = pt_scr%qo*luminosity_linear*exp(ln10*(pt_scr%qpl-log10_star_mass-4.0d0*pt_scr%qtl+cfg%log10_solar_luminosity-cgl+ &
              cdelrl))*temperature_rotation_factor/pressure_rotation_factor
       if (pulse_radiative_gradient-pt_scr%qdela .le. 1.0d-6) then
         pulse_gradient = pulse_radiative_gradient
       else
         pulse_gradient = pt_scr%qdela
       end if
! Store the depth step and the gradients (grad_rad, grad, grad_ad).
       atmo_struct%atmo_delta_depth(step_index) = delta_tau_step
       atmo_struct%atmo_gradients(1,step_index) = pulse_radiative_gradient
       atmo_struct%atmo_gradients(2,step_index) = pulse_gradient
       atmo_struct%atmo_gradients(3,step_index) = pt_scr%qdela
       atmo_struct%num_atm_points = step_index
       prev_tau = tau_now
       prev_density = density_now
end subroutine store_atmosphere_gradients

! ---------------------------------------------------------------
! The envelope: log(m/M), log T and log r against log P, from the
! surface point down to the fitting mass senv, storing env_struct, which
! is then reversed to run from the centre outward. When senv lies above
! senv_thin_envelope there is no envelope to integrate, and the vertex
! is the surface point itself. ierr after maxstp steps.
subroutine integrate_envelope
      use math_lib
! An envelope too thin to integrate: the vertex is the surface point.
      if(cfg%senv.gt.senv_thin_envelope) then
       if(save_boundary_flag) then
          vtx_logp(vertex_index) = atm_point%log10_pressure
          vtx_logr(vertex_index) = log10_radius
          vtx_logt(vertex_index) = atm_point%log10_temperature
       endif
      else
      call start_envelope_integration()
      if (ierr /= 0) return
      call store_first_envelope_point()
      call integrate_to_fitting_mass()
      if (ierr /= 0) return
      end if
      call invert_envelope_points()

      if (solver_diagnostics()) write(run_log_unit,215)num_ok,num_bad,mass_diff_remaining,y(1),(err_sum(jj),jj=1,3)
 215  format(1X,'ENVELOPE INTEGRATION COMPLETE',1X, &
           'NUMBER OF STEPS ACCEPTED',I5,' REJECTED', &
           I5/5X,'SENV-LAST M=',1PE22.13,'  LAST M=',E22.13/ &
           5X,'MAX RELATIVE ERRORS:M ',1PE14.5,'  T ',E14.5, &
           '  R ',E14.5)

end subroutine integrate_envelope

! ---------------------------------------------------------------
! The envelope integration's starting values and numerical parameters,
! and the derivatives at the photosphere.
subroutine start_envelope_integration()
! Start at atm_point (the photosphere, or the table's surface point):
! x = log P, y = (log(m/M) = 0, log T, log r). Three equations, the
! envelope's tolerance and steps, and done when log(m/M) is within
! tolerance_fraction |senv| of senv.
      in_atmosphere = .false.
      do i = 1,3
       err_sum(i) = 0.0d0
      end do
      num_ok = 0
      num_bad = 0
      indep_var = atm_point%log10_pressure
      initial_mass_coord = 0.0d0
      y(1) = initial_mass_coord
      y(2) = atm_point%log10_temperature
      y(3) = log10_radius
      num_eqs = 3
      tolerance = cfg%env_error_tol
      h_max = cfg%env_step_max
      h_min = cfg%env_step_min
      h_step = cfg%env_step_begin
      step_tolerance = dabs(tolerance_fraction*cfg%senv)
! Forget a stored envelope state that belongs to this vertex. The
! pressure limit applies only when the boundary is not being saved.
      if(stored_vertex_index.eq.vertex_index) stored_vertex_index = 0
      store_flag_set = .false.
      if (save_boundary_flag) then
         pressure_limit_test_flag = .false.
      else
         pressure_limit_test_flag = .true.
      end if
! dy/dx at the start of the first step.
      call compute_envelope_derivs(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
           temperature_rotation_factor,log10_gravity,in_atmosphere, &
           want_derivatives,conductive_opacity_flag,log10_radius, &
           log10_teff,hydrogen_fraction,metal_fraction,saha_state, &
           cfg%log10_total_mass, jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if
end subroutine start_envelope_integration

! ---------------------------------------------------------------
! The surface point as point 1 of env_struct, as compute_envelope_derivs
! left it in pt_scr; also starts the convection-zone bookkeeping.
subroutine store_first_envelope_point()
      use math_lib
      cz_in_envelope = .false.
      taucz_env_accum = 0.0d0
      env_struct%env_log10_density(1) = pt_scr%current_log10_density
      env_struct%env_log10_pressure(1) = pt_scr%current_log10_pressure
      env_struct%env_log10_radius(1) = pt_scr%current_log10_radius
      env_struct%env_log10_mass(1) = pt_scr%current_log10_mass
      env_struct%env_log10_temperature(1) = pt_scr%current_log10_temperature
      env_struct%env_hydrogen_fraction(1) = hydrogen_fraction
      env_struct%env_metal_fraction(1) = metal_fraction
      env_struct%env_convective_flag(1) = pt_scr%current_velocity.gt.0.0d0
! The gradients (grad_rad, grad_ad, grad), the convective velocity and
! beta.
      env_struct%env_gradients(1,1) = pt_scr%current_gradients(1)
      env_struct%env_gradients(2,1) = pt_scr%current_gradients(2)
      env_struct%env_gradients(3,1) = pt_scr%current_gradients(3)
      env_struct%env_convective_velocity(1) = pt_scr%current_velocity
      env_struct%env_beta(1) = pt_scr%current_beta
! Gamma1 and c_P (as in store_atmosphere_point), the ionization
! fractions, the opacity and its derivatives, L (in Lsun) and
! (d ln rho/d ln T)_P.
      chi_rho = 1.0d0/pt_scr%qqdp
      chi_t = -chi_rho*pt_scr%qqdt
      specific_heat_cv = pt_scr%qqcp - exp(ln10*(pt_scr%current_log10_pressure-pt_scr%current_log10_density- &
           pt_scr%current_log10_temperature))*chi_t**2/chi_rho
      env_struct%env_gamma1(1) = chi_rho*pt_scr%qqcp/specific_heat_cv
      env_struct%env_specific_heat_cp(1) = pt_scr%qqcp
      env_struct%env_ion_fraction(1,1) = pt_scr%current_ion_fraction(1)
      env_struct%env_ion_fraction(2,1) = pt_scr%current_ion_fraction(2)
      env_struct%env_ion_fraction(3,1) = pt_scr%current_ion_fraction(3)
      env_struct%env_opacity(1) = pt_scr%current_opacity
      env_struct%env_dlnkap_dlnrho(1) = pt_scr%qqod
      env_struct%env_dlnkap_dlnt(1) = pt_scr%qqot
      env_struct%env_luminosity(1) = luminosity_linear
      env_struct%env_dlnrho_dlnt(1) = pt_scr%qqdt
      if(env_struct%env_convective_flag(1))cz_in_envelope = .true.
      env_struct%num_env_points = 1
end subroutine store_first_envelope_point

! ---------------------------------------------------------------
! Integrate the envelope in pressure down to the fitting mass, storing
! each point in env_struct; ierr after maxstp steps.
subroutine integrate_to_fitting_mass()
      do step_index = 1,maxstp
       call limit_envelope_step()
       call take_envelope_step()
       if (ierr /= 0) return
       if(h_did.eq.h_step) then
          num_ok = num_ok + 1
       else
          num_bad = num_bad + 1
       endif
! Done when log(m/M) is within step_tolerance of senv: then interpolate
! the vertex at senv.
       mass_diff_remaining = cfg%senv - y(1)
       call store_envelope_point()
       if(dabs(mass_diff_remaining).le.step_tolerance)then
          call interpolate_envelope_vertex()
          exit
! Otherwise, with save_boundary_flag, store the state where log T
! first reaches tenv (log P, log T, log R, log(m/M)) for the surface
! boundary.
       else if(.not.store_flag_set) then
          if(y(2).ge.cfg%tenv .and. save_boundary_flag) then
             store_flag_set = .true.
             stored_vertex_index = vertex_index
             stored_envelope_state(1) = indep_var
             stored_envelope_state(2) = y(2)
             stored_envelope_state(3) = y(3)
             stored_envelope_state(4) = y(1)
          endif
       endif
       if(h_next.lt.h_min) h_next = h_min
       h_step = h_next
      end do
      if (step_index .gt. maxstp) then
! The integration did not reach the fitting mass in maxstp steps:
! report it and return an error.
      write(terminal_unit,911)
      write(run_log_unit,911)
 911  format(5X,'ENVELOPE INTEGRATION FAILED AFTER MAXSTP TRIES.',1X, &
           'I QUIT')
      jerr = 1
      ierr = jerr
      return
      end if
end subroutine integrate_to_fitting_mass

! ---------------------------------------------------------------
! The error scale of the step, and the step limited so it neither
! overshoots the fitting mass nor exceeds the maximum size or plim.
subroutine limit_envelope_step()
       x_start = indep_var
       do i = 1,num_eqs
          y_start(i) = y(i)
          y_scale(i) = dabs(y(i)) + dabs(h_step*dydx(i))+tiny
       end do
       mass_next_step = y(1) + h_step*dydx(1)
       if(cfg%senv - y(1).gt.0.0d0 .or. cfg%senv - mass_next_step.gt.0.0d0) then
! The integration has passed the fitting mass, or the next step would:
! aim the step at it with the linear estimate
! h = (senv - y(1)) / dydx(1), a step back when it has passed.
          h_step = (cfg%senv - y(1))/dydx(1)
       endif
! Keep |h_step| within the largest step.
       if(h_step.lt.0.0d0) then
          if(h_step.lt.-h_max) h_step = -h_max
       else
          if(h_step.gt.h_max)  h_step = h_max
       endif
! log10_pressure_limit (legacy PLIM) estimates the pressure at the end
! of the integration, from the pressure of the last model point. The
! first step that would cross it stops at it.
       if(pressure_limit_test_flag) then
          if(indep_var + h_step.gt.log10_pressure_limit) then
             h_step = log10_pressure_limit - indep_var
             pressure_limit_test_flag = .false.
          endif
       endif
end subroutine limit_envelope_step

! ---------------------------------------------------------------
! One Bulirsch-Stoer step of the envelope, then the derivatives at its end.
subroutine take_envelope_step()
! One step from log P to log P + h_step: h_step is the step tried,
! h_did the step taken and h_next the predicted next step.
       call load_equations(.false.)
       call bsstep(equations,y,dydx,num_eqs,indep_var,h_step,tolerance,y_scale,h_did, &
              h_next,step_err,max_stage_index,extrap_order,hydrogen_fraction,jerr)
       call unload_equations()
! A bsstep failure returns the numerics-termination code, as in the
! atmosphere.
       if (jerr /= 0) then
          ierr = numerics_termination
          return
       end if
! Add up the error estimates.
       do i = 1,3
          err_sum(i) = err_sum(i) + step_err(i)
       end do
! dy/dx at the start of the next step.
       call compute_envelope_derivs(indep_var,y,dydx,luminosity_linear,pressure_rotation_factor, &
              temperature_rotation_factor,log10_gravity,in_atmosphere, &
              want_derivatives,conductive_opacity_flag,log10_radius, &
              log10_teff,hydrogen_fraction,metal_fraction,saha_state, &
              cfg%log10_total_mass, jerr)
       if (jerr /= 0) then
          ierr = jerr
          return
       end if
end subroutine take_envelope_step

! ---------------------------------------------------------------
! The new envelope point into env_struct, unless the step went past the
! fitting mass by more than step_tolerance. Below the first convective
! point it also sums dr/v_conv over the convective points into
! taucz_env_accum, which nothing reads.
subroutine store_envelope_point()
      use math_lib
         if(mass_diff_remaining.le.step_tolerance)then
            env_struct%num_env_points = env_struct%num_env_points + 1
            env_struct%env_log10_density(env_struct%num_env_points) = pt_scr%current_log10_density
            env_struct%env_log10_pressure(env_struct%num_env_points) = pt_scr%current_log10_pressure
            env_struct%env_log10_radius(env_struct%num_env_points) = pt_scr%current_log10_radius
            env_struct%env_log10_mass(env_struct%num_env_points) = pt_scr%current_log10_mass
            env_struct%env_log10_temperature(env_struct%num_env_points) = pt_scr%current_log10_temperature
            env_struct%env_hydrogen_fraction(env_struct%num_env_points) = hydrogen_fraction
            env_struct%env_metal_fraction(env_struct%num_env_points) = metal_fraction
            env_struct%env_convective_flag(env_struct%num_env_points) = pt_scr%current_velocity.gt.0.0d0
! The gradients (grad_rad, grad_ad, grad), the convective velocity and
! beta.
            env_struct%env_gradients(1,env_struct%num_env_points) = pt_scr%current_gradients(1)
            env_struct%env_gradients(2,env_struct%num_env_points) = pt_scr%current_gradients(2)
            env_struct%env_gradients(3,env_struct%num_env_points) = pt_scr%current_gradients(3)
            env_struct%env_convective_velocity(env_struct%num_env_points) = pt_scr%current_velocity
            env_struct%env_beta(env_struct%num_env_points) = pt_scr%current_beta
! Gamma1, c_P, the ionization fractions, the opacity and its
! derivatives, L and (d ln rho/d ln T)_P, as for point 1.
            chi_rho = 1.0d0/pt_scr%qqdp
            chi_t = -chi_rho*pt_scr%qqdt
            specific_heat_cv = pt_scr%qqcp - exp(ln10*(pt_scr%current_log10_pressure-pt_scr%current_log10_density- &
                 pt_scr%current_log10_temperature))*chi_t**2/chi_rho
            env_struct%env_gamma1(env_struct%num_env_points) = chi_rho*pt_scr%qqcp/specific_heat_cv
            env_struct%env_specific_heat_cp(env_struct%num_env_points) = pt_scr%qqcp
            env_struct%env_ion_fraction(1,env_struct%num_env_points) = pt_scr%current_ion_fraction(1)
            env_struct%env_ion_fraction(2,env_struct%num_env_points) = pt_scr%current_ion_fraction(2)
            env_struct%env_ion_fraction(3,env_struct%num_env_points) = pt_scr%current_ion_fraction(3)
            env_struct%env_opacity(env_struct%num_env_points) = pt_scr%current_opacity
            env_struct%env_dlnkap_dlnrho(env_struct%num_env_points) = pt_scr%qqod
            env_struct%env_dlnkap_dlnt(env_struct%num_env_points) = pt_scr%qqot
            env_struct%env_luminosity(env_struct%num_env_points) = luminosity_linear
            env_struct%env_dlnrho_dlnt(env_struct%num_env_points) = pt_scr%qqdt

! From the first convective point on, add dr/v_conv of each convective
! point (dr from the previous point, which lies further out).
            if(.not.cz_in_envelope)then
               if(env_struct%env_convective_flag(env_struct%num_env_points))then
                  cz_in_envelope = .true.
               endif
            else if(pt_scr%current_velocity.gt.0.0d0)then
               delta_radius_cz = exp10(env_struct%env_log10_radius(env_struct%num_env_points-1)) - &
                    exp10(env_struct%env_log10_radius(env_struct%num_env_points))
               taucz_env_accum = taucz_env_accum + delta_radius_cz/pt_scr%current_velocity
            endif
         endif
end subroutine store_envelope_point

! ---------------------------------------------------------------
! The fitting-point vertex (log P, log R, log T), interpolated linearly
! in log(m/M) to senv between the start and the end of the last step:
!   v = v_end + w (v_start - v_end),  w = (senv - y(1)) / (y_start(1) - y(1)).
subroutine interpolate_envelope_vertex()
          if(save_boundary_flag) then
             interp_weight = mass_diff_remaining/(y_start(1)-y(1))
             vtx_logp(vertex_index) = indep_var + interp_weight*(x_start - indep_var)
             vtx_logr(vertex_index) = y(3) + interp_weight*(y_start(3) - y(3))
             vtx_logt(vertex_index) = y(2) + interp_weight*(y_start(2) - y(2))
             if(print_flag .and. solver_diagnostics())write(run_log_unit,230)vtx_logp(vertex_index),vtx_logt(vertex_index),vtx_logr(vertex_index),cfg%senv
          endif
 230     format(4X,3F16.12,8X,F16.12)
end subroutine interpolate_envelope_vertex

! ---------------------------------------------------------------
! Reverse env_struct from surface-inward to centre-outward, when the
! envelope was integrated (senv below senv_thin_envelope).
subroutine invert_envelope_points()
! (.le., the complement of the skip test above: at senv exactly
! senv_thin_envelope the envelope is integrated and must be reversed
! too; the F77 envint.f tested .LT., fixed 2026-10)
      if(cfg%senv.le.senv_thin_envelope)then
         associate (n => env_struct%num_env_points)
            env_struct%env_log10_density(1:n) = env_struct%env_log10_density(n:1:-1)
            env_struct%env_log10_pressure(1:n) = env_struct%env_log10_pressure(n:1:-1)
            env_struct%env_log10_radius(1:n) = env_struct%env_log10_radius(n:1:-1)
            env_struct%env_log10_mass(1:n) = env_struct%env_log10_mass(n:1:-1)
            env_struct%env_log10_temperature(1:n) = env_struct%env_log10_temperature(n:1:-1)
            env_struct%env_hydrogen_fraction(1:n) = env_struct%env_hydrogen_fraction(n:1:-1)
            env_struct%env_metal_fraction(1:n) = env_struct%env_metal_fraction(n:1:-1)
            env_struct%env_convective_flag(1:n) = env_struct%env_convective_flag(n:1:-1)
            env_struct%env_opacity(1:n) = env_struct%env_opacity(n:1:-1)
            env_struct%env_dlnkap_dlnrho(1:n) = env_struct%env_dlnkap_dlnrho(n:1:-1)
            env_struct%env_dlnkap_dlnt(1:n) = env_struct%env_dlnkap_dlnt(n:1:-1)
            env_struct%env_luminosity(1:n) = env_struct%env_luminosity(n:1:-1)
            env_struct%env_dlnrho_dlnt(1:n) = env_struct%env_dlnrho_dlnt(n:1:-1)
            env_struct%env_gradients(:,1:n) = env_struct%env_gradients(:,n:1:-1)
            env_struct%env_convective_velocity(1:n) = env_struct%env_convective_velocity(n:1:-1)
            env_struct%env_beta(1:n) = env_struct%env_beta(n:1:-1)
            env_struct%env_gamma1(1:n) = env_struct%env_gamma1(n:1:-1)
            env_struct%env_specific_heat_cp(1:n) = env_struct%env_specific_heat_cp(n:1:-1)
            env_struct%env_ion_fraction(:,1:n) = env_struct%env_ion_fraction(:,n:1:-1)
         end associate
      endif
end subroutine invert_envelope_points



! Copy the integration's current values into bsstep's equations, for
! the atmosphere (atmosphere = .true.) or the envelope.
subroutine load_equations(atmosphere)
      logical, intent(in) :: atmosphere
      equations%atmosphere = atmosphere
      equations%luminosity_linear = luminosity_linear
      equations%pressure_rotation_factor = pressure_rotation_factor
      equations%temperature_rotation_factor = temperature_rotation_factor
      equations%log10_gravity = log10_gravity
      equations%in_atmosphere = in_atmosphere
      equations%want_derivatives = want_derivatives
      equations%conductive_opacity_flag = conductive_opacity_flag
      equations%log10_radius = log10_radius
      equations%log10_teff = log10_teff
      equations%hydrogen_fraction = hydrogen_fraction
      equations%metal_fraction = metal_fraction
      equations%saha_state = saha_state
      equations%atm_choice = cfg%atm_choice
      equations%atm_hras = cfg%atm_hras
      equations%log10_total_mass = cfg%log10_total_mass
end subroutine load_equations

! Copy back the three values the integrands update (nothing else reads
! them during a step).
subroutine unload_equations()
      log10_radius = equations%log10_radius
      log10_teff = equations%log10_teff
      saha_state = equations%saha_state
end subroutine unload_equations

end subroutine integrate_envelope_atmosphere

end module envelope_integration_kernel

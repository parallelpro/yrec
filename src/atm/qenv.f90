!----------------------------------------------------------------------
! qenv
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original qenv.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Computes the derivatives (dS/dP, dT/dP, dR/dP) of the envelope
! structure equations with pressure as the independent variable:
! mass continuity, the temperature gradient (radiative or convective,
! via TPGRAD), and hydrostatic equilibrium. Also returns T,R,D,O and
! the ionization fractions for output purposes.
!
! Dummy-argument names log10_pressure_indep/luminosity_linear/
! pressure_rotation_factor/temperature_rotation_factor/log10_gravity/
! in_atmosphere/want_derivatives/conductive_opacity_flag/print_flag/
! log10_radius/log10_teff/hydrogen_fraction/metal_fraction match the
! sibling routine qatm.f90's naming at the equivalent call-site
! positions (both are called via EXTERNAL from ENVINT/BSSTEP).
! pressure_rotation_factor/temperature_rotation_factor (FPL/FTL) match
! tpgrad.f90's naming, discovered from that file's own analysis of the
! DELR formula.
subroutine qenv(log10_pressure_indep, y, dydx, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, log10_gravity, &
     in_atmosphere, want_derivatives, conductive_opacity_flag, print_flag, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     env_call_count, saha_state)

      use run_diag_lib
      use pulse_diag_lib
      use envelope_comp_lib
      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(in) :: log10_pressure_indep
      double precision, intent(inout) :: y(3), dydx(3)
      double precision, intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(in) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag, print_flag
      double precision, intent(inout) :: log10_radius
      double precision, intent(in) :: log10_teff
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: env_call_count, saha_state

      double precision :: ion_fraction(3)

      save

! --- locals ---
      double precision :: log10_pressure, log10_mass, log10_temperature
      double precision :: temperature, pressure, log10_density, density
      double precision :: beta, beta_inverse, beta14, specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision :: dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
           adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      double precision :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component
      double precision :: convective_velocity
      logical :: is_convective

      log10_pressure = log10_pressure_indep
      log10_mass = y(1) + env_comp%stotal
      log10_temperature = y(2)
      log10_radius = y(3)
      if(use_mhd_eos)then
         call meqos(log10_temperature,temperature,log10_pressure,pressure, &
              log10_density,density,hydrogen_fraction,metal_fraction,beta, &
              beta_inverse,beta14,ion_fraction,specific_gas_constant, &
              ion_mean_weight_inverse,electron_mean_weight_inverse, &
              electron_degeneracy_parameter,dlnrho_dlnt,dlnrho_dlnp, &
              specific_heat_cp,adiabatic_gradient,dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt,adiabatic_gradient_dt,adiabatic_gradient_dp, &
              specific_heat_cp_dt,specific_heat_cp_dp)
      else
         call eqstat(log10_temperature,temperature,log10_pressure,pressure, &
              log10_density,density,hydrogen_fraction,metal_fraction,beta, &
              beta_inverse,beta14,ion_fraction,specific_gas_constant, &
              ion_mean_weight_inverse,electron_mean_weight_inverse, &
              electron_degeneracy_parameter,dlnrho_dlnt,dlnrho_dlnp, &
              specific_heat_cp,adiabatic_gradient,dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt,adiabatic_gradient_dt,adiabatic_gradient_dp, &
              specific_heat_cp_dt,specific_heat_cp_dp,want_derivatives, &
              in_atmosphere,saha_state)
      endif
      call getopac(log10_density, log10_temperature, hydrogen_fraction, &
           metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, &
           dlnkap_dlnt, ion_fraction)
      iovim = -1
      call tpgrad(log10_temperature,temperature,log10_pressure,pressure, &
           density,log10_radius,log10_mass,luminosity_linear,opacity, &
           dlnrho_dlnt,dlnrho_dlnp,dlnkap_dlnt,dlnkap_dlnrho, &
           specific_heat_cp,actual_gradient,radiative_gradient, &
           adiabatic_gradient,dlnrho_dlnt_dt,dlnrho_dlnp_dt, &
           adiabatic_gradient_dt,adiabatic_gradient_dp,dgrad_dt_component, &
           dgrad_dp_component,dgrad_dr_component,specific_heat_cp_dt, &
           specific_heat_cp_dp,convective_velocity,want_derivatives, &
           is_convective,pressure_rotation_factor,temperature_rotation_factor, &
           log10_teff)
      dydx(1) = -dexp(ln10*(c4pil+4.0d0*log10_radius+log10_pressure-cgl- &
           log10_mass-log10_mass))/pressure_rotation_factor
      dydx(2) = actual_gradient
      dydx(3) = -dexp(ln10*(log10_pressure+log10_radius-cgl-log10_mass- &
           log10_density))*pressure_rotation_factor
      env_call_count = env_call_count + 1
! 07/02 ALWAYS STORE THE BASIC STRUCTURE VARIABLES.
      run_diag%current_log10_pressure = log10_pressure
      run_diag%current_log10_temperature = log10_temperature
      run_diag%current_log10_mass = log10_mass - env_comp%stotal
      run_diag%current_log10_radius = log10_radius
      run_diag%current_log10_density = log10_density
      run_diag%current_velocity = convective_velocity
! JVS 08/13 ALWAYS STORE GRADIENTS (FOR TRACKING CZ)
       run_diag%current_gradients(1) = radiative_gradient
       run_diag%current_gradients(2) = adiabatic_gradient
       run_diag%current_gradients(3) = actual_gradient
       run_diag%current_beta = beta ! added 03/14
! JVS 08/25 ALSO ALWAYS SAVE ADDITIONAL INFO FOR PROFILE
      run_diag%current_ion_fraction(1) = ion_fraction(1)
      run_diag%current_ion_fraction(2) = ion_fraction(2)
      run_diag%current_ion_fraction(3) = ion_fraction(3)
      pulse_diag%qqdp = dlnrho_dlnp
      pulse_diag%qqdt = dlnrho_dlnt
      pulse_diag%qqcp = specific_heat_cp

      if(print_flag .or. pulse_diag%lpumod) then
       run_diag%current_opacity = opacity
       run_diag%current_ion_fraction(1) = ion_fraction(1)
       run_diag%current_ion_fraction(2) = ion_fraction(2)
       run_diag%current_ion_fraction(3) = ion_fraction(3)
       pulse_diag%qtl = log10_temperature
       pulse_diag%qt = dexp(ln10*log10_temperature)
       pulse_diag%qpl = log10_pressure
       pulse_diag%qp = dexp(ln10*log10_pressure)
       pulse_diag%qdl = log10_density
       pulse_diag%qd = dexp(ln10*log10_density)
       pulse_diag%qo = opacity
       pulse_diag%qol = log10_opacity
       pulse_diag%qfs = dexp(ln10*(log10_mass-env_comp%stotal))
       pulse_diag%qqdp = dlnrho_dlnp
       pulse_diag%qqed = 0.0d0
       pulse_diag%qqod = dlnkap_dlnrho
       pulse_diag%qqot = dlnkap_dlnt
       pulse_diag%qdel = actual_gradient
       pulse_diag%qqdt = dlnrho_dlnt
       pulse_diag%qdela = adiabatic_gradient
       pulse_diag%qqcp = specific_heat_cp
       pulse_diag%qrmu = specific_gas_constant
       pulse_diag%qemu = electron_mean_weight_inverse
      endif

! KC 2025-05-31 THESE MUST BE RETAINED FOR EXTERNAL PROCEDURE COMPATIBILITY.
      if (.false.) print *, log10_gravity, conductive_opacity_flag

      return
end subroutine qenv

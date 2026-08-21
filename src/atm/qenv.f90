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
! common/pulse/: only pulsation_file_version is used here. Naming
! matches wrtmod.f90.
      double precision :: pulsation_mass_msun
      logical :: pulsation_output_active
      integer :: pulsation_file_version
      common/pulse/pulsation_mass_msun, pulsation_output_active, &
           pulsation_file_version
! common/pulse1/: only lpumod is used here; remaining members are
! unused placeholders. Naming matches wrtmod.f90.
      double precision :: pulse_dlnrho_dlnp(json), pulse_dlneps_dlnrho(json), &
           pulse_dlneps_dlnt(json), pulse_dlnkap_dlnrho(json), &
           pulse_dlnkap_dlnt(json), pulse_specific_heat(json), &
           pulse_mean_molecular_weight(json), pulse_dlnrho_dlnt(json), &
           pulse_electron_mean_molecular_weight(json)
      logical :: lpumod
      common/pulse1/pulse_dlnrho_dlnp, pulse_dlneps_dlnrho, &
           pulse_dlneps_dlnt, pulse_dlnkap_dlnrho, pulse_dlnkap_dlnt, &
           pulse_specific_heat, pulse_mean_molecular_weight, &
           pulse_dlnrho_dlnt, pulse_electron_mean_molecular_weight, lpumod
! common/pulse2/: values saved here for pulsation/diagnostic output
! when print_flag/lpumod is set. Naming matches wrtmod.f90 (kept close
! to the original cryptic names).
      double precision :: qqdp, qqed, qqet, qqod, qqot, qdel, qdela, qqcp, &
           qrmu, qtl, qpl, qdl, qo, qol, qt, qp, qqdt, qemu, qd, qfs
      common/pulse2/qqdp, qqed, qqet, qqod, qqot, qdel, qdela, qqcp, qrmu, &
           qtl, qpl, qdl, qo, qol, qt, qp, qqdt, qemu, qd, qfs
! common/comp/: only stotal is used here. Naming matches getopac.f90.
      double precision :: envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv(12), xnew, znew, stotal, senv
      common/comp/envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, amuenv, fxenv, xnew, znew, stotal, senv
! common/envprt/: all used/set here. Naming is local to this batch
! (shared with envint.f90's usage of this block).
      double precision :: current_log10_pressure, current_log10_temperature, &
           current_log10_radius, current_log10_mass, current_log10_density, &
           current_opacity, current_beta, current_gradients(3), &
           current_ion_fraction(3), current_velocity
      common/envprt/current_log10_pressure, current_log10_temperature, &
           current_log10_radius, current_log10_mass, current_log10_density, &
           current_opacity, current_beta, current_gradients, &
           current_ion_fraction, current_velocity
! common/mhd/: only use_mhd_eos is used here. Naming matches mhdtbl.f90.
      logical :: use_mhd_eos
      integer :: unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, &
           unit_centre2, unit_centre3, unit_centre4, unit_centre5
      common/mhd/use_mhd_eos, unit_zams_a, unit_zams_b, unit_zams_c, &
           unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

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
      log10_mass = y(1) + stotal
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
      current_log10_pressure = log10_pressure
      current_log10_temperature = log10_temperature
      current_log10_mass = log10_mass - stotal
      current_log10_radius = log10_radius
      current_log10_density = log10_density
      current_velocity = convective_velocity
! JVS 08/13 ALWAYS STORE GRADIENTS (FOR TRACKING CZ)
       current_gradients(1) = radiative_gradient
       current_gradients(2) = adiabatic_gradient
       current_gradients(3) = actual_gradient
       current_beta = beta ! added 03/14
! JVS 08/25 ALSO ALWAYS SAVE ADDITIONAL INFO FOR PROFILE
      current_ion_fraction(1) = ion_fraction(1)
      current_ion_fraction(2) = ion_fraction(2)
      current_ion_fraction(3) = ion_fraction(3)
      qqdp = dlnrho_dlnp
      qqdt = dlnrho_dlnt
      qqcp = specific_heat_cp

      if(print_flag .or. lpumod) then
       current_opacity = opacity
       current_ion_fraction(1) = ion_fraction(1)
       current_ion_fraction(2) = ion_fraction(2)
       current_ion_fraction(3) = ion_fraction(3)
       qtl = log10_temperature
       qt = dexp(ln10*log10_temperature)
       qpl = log10_pressure
       qp = dexp(ln10*log10_pressure)
       qdl = log10_density
       qd = dexp(ln10*log10_density)
       qo = opacity
       qol = log10_opacity
       qfs = dexp(ln10*(log10_mass-stotal))
       qqdp = dlnrho_dlnp
       qqed = 0.0d0
       qqod = dlnkap_dlnrho
       qqot = dlnkap_dlnt
       qdel = actual_gradient
       qqdt = dlnrho_dlnt
       qdela = adiabatic_gradient
       qqcp = specific_heat_cp
       qrmu = specific_gas_constant
       qemu = electron_mean_weight_inverse
      endif

! KC 2025-05-31 THESE MUST BE RETAINED FOR EXTERNAL PROCEDURE COMPATIBILITY.
      if (.false.) print *, log10_gravity, conductive_opacity_flag

      return
end subroutine qenv

!----------------------------------------------------------------------
! envelope_integrands
!----------------------------------------------------------------------
! The right-hand sides of the outer-layer integrations of
! envelope_integration_kernel: the atmosphere in optical depth
! (compute_atmosphere_derivs) and the envelope in pressure
! (compute_envelope_derivs). The kernel's envelope_equations calls
! either, with the arguments each uses; both leave the physics of the
! point in pt (point_physics, an argument) for the kernel to record.
! gamma1_from_eos_derivatives is the Gamma1 of a point the kernel records.
module envelope_integrands
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: point_physics, compute_atmosphere_derivs, compute_envelope_derivs, &
           gamma1_from_eos_derivatives

! The physics of the point an integrand just evaluated. Each integrand
! call fills it, and the integration kernel reads it to store the point
! in atm_points or env_points. Both integrands set the first three
! groups; the envelope's also sets the last (an atmosphere point's
! optical depth and ionization fractions are in atm_point, and the
! kernel computes its beta and gradients when it stores the point).
      type, public :: point_physics
! the point: log P (dyn/cm^2), log T (K), log rho and rho (g/cm^3)
            real(dp) :: log10_pressure = 0.0_dp, log10_temperature = 0.0_dp, &
                 log10_density = 0.0_dp, density = 0.0_dp
! from the EOS: (d ln rho/d ln P)_T, (d ln rho/d ln T)_P, c_P (erg/(g K))
! and the adiabatic gradient
            real(dp) :: dlnrho_dlnp = 0.0_dp, dlnrho_dlnt = 0.0_dp, &
                 specific_heat_cp = 0.0_dp, adiabatic_gradient = 0.0_dp
! the opacity (cm^2/g), (d ln kappa/d ln rho)_T and (d ln kappa/d ln T)_rho
            real(dp) :: opacity = 0.0_dp, dlnkap_dlnrho = 0.0_dp, &
                 dlnkap_dlnt = 0.0_dp
! an envelope point only: log r (cm) and log(m/M); the radiative
! gradient and the temperature gradient d log T/d log P; the convective
! velocity (cm/s, 0 where radiative); beta = P_gas/P and the ionization
! fractions of H, He+ and He++ (the EOS's)
            real(dp) :: log10_radius = 0.0_dp, log10_mass_fraction = 0.0_dp
            real(dp) :: radiative_gradient = 0.0_dp, temperature_gradient = 0.0_dp, &
                 convective_velocity = 0.0_dp
            real(dp) :: beta = 0.0_dp, ion_fraction(3) = 0.0_dp
      end type point_physics

contains

! ---------------------------------------------------------------
! compute_atmosphere_derivs (legacy qatm.f)
! ---------------------------------------------------------------
! The right-hand side of the atmosphere integration, hydrostatic
! equilibrium in optical depth, dP/dtau = g_eff/kappa:
!   d log P / d log tau = g_eff tau / (kappa P),   g_eff = g f_P,
! with x = log tau and y(1) = log P. T(tau) comes from the T-tau
! relation of atm_choice (atm/t_tau_relations.f90; for the
! Eddington gray atmosphere T^4 = (3/4) Teff^4 (tau + 2/3)), and rho
! and kappa from the EOS and opacity at (P, T). Also stores the point
! (tau, P, T, rho, kappa, the ionization fractions) in atm_point and
! its physics in pt.
subroutine compute_atmosphere_derivs(log10_optical_depth, y, dydx, &
     pressure_rotation_factor, log10_gravity, &
     in_atmosphere, want_derivatives, &
     log10_teff, hydrogen_fraction, metal_fraction, &
     saha_iterations, atm_choice, hsra_log_t_offset, composition, atm_point, pt, ierr)

      use eos_lib, only: eos_get, eos_result
      use kap_lib, only: i_dlnkap_dlnrho, i_dlnkap_dlnt, i_kap, kap_get, &
           num_kap_results
      use envelope_points, only: atmosphere_point_state
      use atm_lib, only: ttau_log10_temperature
      use physical_constants, only: ln10
      use math_lib
      implicit none

! One equation: x = log tau, y(1) = log P.
      real(dp), intent(in) :: log10_optical_depth
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: dydx(:)
      real(dp), intent(in) :: pressure_rotation_factor, log10_gravity
      logical, intent(in) :: in_atmosphere, want_derivatives
      real(dp), intent(in) :: log10_teff
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: saha_iterations
! atm_choice: the integration's T(tau) relation, the kernel's choice in
! use (gray after a failed Allard lookup), not the star's, which the
! driver resets only after the integration.
      integer, intent(in) :: atm_choice
! hsra_log_t_offset: the HSRA T(tau) relation's offset.
      real(dp), intent(in) :: hsra_log_t_offset
! composition: the outer shell's composition, for the EOS's
! Debye-Hueckel correction.
      real(dp), intent(in) :: composition(:)
! atm_point: the integration's current point, written here.
      type(atmosphere_point_state), intent(inout) :: atm_point
! pt: the point's physics, written here for the kernel.
      type(point_physics), intent(inout) :: pt
! ierr: nonzero when the EOS or the opacity fails; bsstep passes it up
! to the kernel.
      integer, intent(out) :: ierr

! --- locals ---
      real(dp) :: effective_gravity, optical_depth
      real(dp) :: log10_temperature, log10_pressure
! The EOS and opacity results at this point.
      type(eos_result) :: eos_res
      real(dp) :: kap_res(num_kap_results)
      ierr = 0
      effective_gravity = exp(ln10*log10_gravity)*pressure_rotation_factor
      atm_point%log10_tau = log10_optical_depth
      optical_depth = exp(ln10*atm_point%log10_tau)
! T(tau) from the T-tau relation.
      log10_temperature = ttau_log10_temperature(optical_depth, &
           log10_teff, atm_choice, hsra_log_t_offset)
      log10_pressure = y(1)
      call eos_get(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, eos_res, want_derivatives, in_atmosphere, &
           saha_iterations, composition_at_zone=composition, ierr=ierr)
      if (ierr /= 0) return
! The opacity at the density the EOS returned.
      call kap_get(eos_res%log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction, kap_res, &
           eos_res%ion_fraction, ierr=ierr)
      if (ierr /= 0) return
! d log P / d log tau = g_eff tau / (kappa P).
      dydx(1) = effective_gravity*optical_depth/ &
           (eos_res%pressure*kap_res(i_kap))
! Store the point for the kernel: log P and log T of the last point
! start the envelope integration, and the rest is output for
! atm_points.
      atm_point%log10_pressure = log10_pressure
      atm_point%log10_temperature = log10_temperature
      atm_point%log10_density = eos_res%log10_density
      atm_point%opacity = kap_res(i_kap)
      atm_point%ion_fraction(1) = eos_res%ion_fraction(1)
      atm_point%ion_fraction(2) = eos_res%ion_fraction(2)
      atm_point%ion_fraction(3) = eos_res%ion_fraction(3)
      pt%log10_pressure = log10_pressure
      pt%log10_temperature = log10_temperature
      pt%log10_density = eos_res%log10_density
      pt%density = exp(ln10*eos_res%log10_density)
      pt%dlnrho_dlnp = eos_res%dlnrho_dlnp
      pt%dlnrho_dlnt = eos_res%dlnrho_dlnt
      pt%specific_heat_cp = eos_res%specific_heat_cp
      pt%adiabatic_gradient = eos_res%adiabatic_gradient
      pt%opacity = kap_res(i_kap)
      pt%dlnkap_dlnrho = kap_res(i_dlnkap_dlnrho)
      pt%dlnkap_dlnt = kap_res(i_dlnkap_dlnt)

      return
end subroutine compute_atmosphere_derivs

! ---------------------------------------------------------------
! compute_envelope_derivs (legacy qenv.f)
! ---------------------------------------------------------------
! The right-hand side of the envelope integration, with x = log P and
! y = (log(m/M), log T, log r), M the total mass (log10_total_mass = log M):
!   d log m / d log P = -4 pi r^4 P / (G m^2 f_P)  (hydrostatic equilibrium)
!   d log T / d log P = grad                       (the temperature gradient)
!   d log r / d log P = -r P / (G m rho f_P)       (mass continuity with
!                                                   hydrostatic equilibrium)
! rho comes from the EOS at (P, T), and grad, radiative or convective,
! from compute_temperature_gradients (legacy TPGRAD). f_P and f_T
! (pressure_rotation_factor, temperature_rotation_factor; legacy FPL,
! FTL) are the rotational distortion factors, 1 without rotation. Also
! stores the point and its physics in pt for env_points.
subroutine compute_envelope_derivs(log10_pressure_indep, y, dydx, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, &
     in_atmosphere, want_derivatives, &
     log10_teff, hydrogen_fraction, metal_fraction, &
     saha_iterations, log10_total_mass, log10_gravitational_constant, composition, pt, ierr)

      use convection_lib, only: compute_temperature_gradients
      use eos_lib, only: eos_get, eos_result
      use kap_lib, only: i_dlnkap_dlnrho, i_dlnkap_dlnt, i_kap, kap_get, &
           num_kap_results
      use physical_constants, only: c4pil, ln10
      use math_lib
      implicit none

      real(dp), intent(in) :: log10_pressure_indep
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: dydx(:)
      real(dp), intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor
      logical, intent(in) :: in_atmosphere, want_derivatives
      real(dp), intent(in) :: log10_teff
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: saha_iterations
! log10_total_mass: log10 M (g), the origin of y(1) = log(m/M).
      real(dp), intent(in) :: log10_total_mass
! log10_gravitational_constant: log10 G (cgs).
      real(dp), intent(in) :: log10_gravitational_constant
! composition: the outer shell's composition, for the EOS's
! Debye-Hueckel correction.
      real(dp), intent(in) :: composition(:)
! pt: the point's physics, written here for the kernel.
      type(point_physics), intent(inout) :: pt
      integer, intent(out) :: ierr

! --- locals ---
      real(dp) :: log10_pressure, log10_mass, log10_temperature, log10_radius
! The EOS and opacity results at this point.
      type(eos_result) :: eos_res
      real(dp) :: kap_res(num_kap_results)
      real(dp) :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component
      real(dp) :: convective_velocity
      logical :: is_convective

! The error code of the gradient call, passed on in ierr.
      integer :: jerr

      ierr = 0
      log10_pressure = log10_pressure_indep
! log m = log(m/M) + log M.
      log10_mass = y(1) + log10_total_mass
      log10_temperature = y(2)
      log10_radius = y(3)
      call eos_get(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, eos_res, want_derivatives, in_atmosphere, &
           saha_iterations, composition_at_zone=composition, ierr=ierr)
      if (ierr /= 0) return
! The opacity at the density the EOS returned.
      call kap_get(eos_res%log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction, kap_res, &
           eos_res%ion_fraction, ierr=ierr)
      if (ierr /= 0) return
! The temperature gradient grad (actual_gradient), grad_rad and the
! convective velocity.
      call compute_temperature_gradients(log10_temperature, log10_pressure, &
           eos_res, kap_res, log10_radius, log10_mass, luminosity_linear, &
           actual_gradient, radiative_gradient, dgrad_dt_component, &
           dgrad_dp_component, dgrad_dr_component, convective_velocity, &
           want_derivatives, is_convective, pressure_rotation_factor, &
           temperature_rotation_factor, log10_teff, .false., jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if
! Hydrostatic equilibrium: d log m / d log P = -4 pi r^4 P / (G m^2 f_P).
      dydx(1) = -exp(ln10*(c4pil+4.0d0*log10_radius+log10_pressure-log10_gravitational_constant- &
           log10_mass-log10_mass))/pressure_rotation_factor
! d log T / d log P = grad.
      dydx(2) = actual_gradient
! Mass continuity: d log r / d log P = (d log r/d log m)(d log m/d log P)
! = -r P / (G m rho f_P), with the same 1/f_P as dydx(1) (Kippenhahn &
! Thomas 1970). (The legacy qenv multiplied by f_P instead, which broke
! mass continuity by a factor f_P^2 in rotating envelopes; bugsweep
! 2026.)
      dydx(3) = -exp(ln10*(log10_pressure+log10_radius-log10_gravitational_constant-log10_mass- &
           eos_res%log10_density))/pressure_rotation_factor
! Store the point for the kernel: the structure variables (log m as
! log(m/M)),
      pt%log10_pressure = log10_pressure
      pt%log10_temperature = log10_temperature
      pt%log10_density = eos_res%log10_density
      pt%density = exp(ln10*eos_res%log10_density)
      pt%log10_radius = log10_radius
      pt%log10_mass_fraction = log10_mass - log10_total_mass
! the EOS derivatives and the opacity with its derivatives,
      pt%dlnrho_dlnp = eos_res%dlnrho_dlnp
      pt%dlnrho_dlnt = eos_res%dlnrho_dlnt
      pt%specific_heat_cp = eos_res%specific_heat_cp
      pt%adiabatic_gradient = eos_res%adiabatic_gradient
      pt%opacity = kap_res(i_kap)
      pt%dlnkap_dlnrho = kap_res(i_dlnkap_dlnrho)
      pt%dlnkap_dlnt = kap_res(i_dlnkap_dlnt)
! the gradients and the convective velocity, which locate the convection
! zone,
      pt%radiative_gradient = radiative_gradient
      pt%temperature_gradient = actual_gradient
      pt%convective_velocity = convective_velocity
! and beta, the gas-pressure fraction, with the ionization fractions.
      pt%beta = eos_res%beta
      pt%ion_fraction(1) = eos_res%ion_fraction(1)
      pt%ion_fraction(2) = eos_res%ion_fraction(2)
      pt%ion_fraction(3) = eos_res%ion_fraction(3)

      return
end subroutine compute_envelope_derivs

! ---------------------------------------------------------------
! Gamma1 = (d ln P/d ln rho)_S of a point from its EOS derivatives:
!   chi_rho = 1/(d ln rho/d ln P)_T,  chi_T = -chi_rho (d ln rho/d ln T)_P,
!   c_V = c_P - (P/(rho T)) chi_T^2/chi_rho,  Gamma1 = chi_rho c_P/c_V
! (Guenther 1995, private communication, as shell_equations has it).
real(dp) function gamma1_from_eos_derivatives(dlnrho_dlnp, dlnrho_dlnt, &
     specific_heat_cp, log10_pressure, log10_density, log10_temperature) result(gamma1)
      use physical_constants, only: ln10
      use math_lib
      real(dp), intent(in) :: dlnrho_dlnp, dlnrho_dlnt, specific_heat_cp
      real(dp), intent(in) :: log10_pressure, log10_density, log10_temperature
      real(dp) :: chi_rho, chi_t, specific_heat_cv
      chi_rho = 1.0d0/dlnrho_dlnp
      chi_t = -chi_rho*dlnrho_dlnt
      specific_heat_cv = specific_heat_cp - exp(ln10*(log10_pressure-log10_density- &
           log10_temperature))*chi_t**2/chi_rho
      gamma1 = chi_rho*specific_heat_cp/specific_heat_cv
end function gamma1_from_eos_derivatives

end module envelope_integrands

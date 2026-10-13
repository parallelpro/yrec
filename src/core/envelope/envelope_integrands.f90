!----------------------------------------------------------------------
! envelope_integrands
!----------------------------------------------------------------------
! The right-hand sides of the outer-layer integrations of
! envelope_integration_kernel: the atmosphere in optical depth
! (compute_atmosphere_derivs) and the envelope in pressure
! (compute_envelope_derivs). Both take the same argument list, so that
! the kernel's envelope_equations can call either, and both leave the
! physics of the point in pt_scr (point_scratch_lib) for the kernel to
! record.
module envelope_integrands
      implicit none
      private
      public :: compute_atmosphere_derivs, compute_envelope_derivs

contains

! ---------------------------------------------------------------
! compute_atmosphere_derivs (legacy qatm.f)
! ---------------------------------------------------------------
! The right-hand side of the atmosphere integration, hydrostatic
! equilibrium in optical depth, dP/dtau = g_eff/kappa:
!   d log P / d log tau = g_eff tau / (kappa P),   g_eff = g f_P,
! with x = log tau and y(1) = log P. T(tau) comes from the T-tau
! relation of star%job%atm_choice (atm/t_tau_relations.f90; for the
! Eddington gray atmosphere T^4 = (3/4) Teff^4 (tau + 2/3)), and rho
! and kappa from the EOS and opacity at (P, T). Also stores the point
! (tau, P, T, rho, kappa, the ionization fractions) in atm_point and
! its physics in pt_scr.
subroutine compute_atmosphere_derivs(log10_optical_depth, y, dydx, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, log10_gravity, &
     in_atmosphere, want_derivatives, conductive_opacity_flag, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     saha_state, atm_choice, ierr)

      use eos_lib
      use kap_lib
      use atmosphere_structure, only: atm_point
      use atm_lib, only: ttau_log10_temperature
      use star_info_lib, only: star, json
      use point_scratch_lib
      use physical_constants
      use math_lib
      implicit none

! One equation: x = log tau, y(1) = log P.
      double precision, intent(in) :: log10_optical_depth
      double precision, intent(in) :: y(1)
      double precision, intent(out) :: dydx(1)
      double precision, intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(in) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag
      double precision, intent(in) :: log10_radius
      double precision, intent(inout) :: log10_teff
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: saha_state
! atm_choice: the integration's T(tau) relation, the kernel's (which
! falls back to gray, 0, after a failed Allard lookup, while
! star%job%atm_choice is reset only after the integration; reading the
! star's here gave HSRA for that one integration, a 2026 regression,
! fixed 2026-10; the F77 envint.f switched the common KTTAU qatm read)
      integer, intent(in) :: atm_choice
! ierr: nonzero when the EOS or the opacity fails; bsstep passes it up
! to the kernel.
      integer, intent(out) :: ierr

! --- locals ---
      double precision :: effective_gravity, optical_depth
      double precision :: log10_temperature, log10_pressure
! The EOS and opacity results at this point.
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results)
      ierr = 0
      effective_gravity = exp(ln10*log10_gravity)*pressure_rotation_factor
      atm_point%log10_tau = log10_optical_depth
      optical_depth = exp(ln10*atm_point%log10_tau)
! T(tau) from the T-tau relation.
      log10_temperature = ttau_log10_temperature(optical_depth, &
           log10_teff, atm_choice, star%atm_hras)
      log10_pressure = y(1)
      call eos_get(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, eos_res, want_derivatives, in_atmosphere, &
           saha_state, ierr=ierr)
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
! atmo_struct.
      atm_point%log10_pressure = log10_pressure
      atm_point%log10_temperature = log10_temperature
      atm_point%log10_density = eos_res%log10_density
      atm_point%opacity = kap_res(i_kap)
      atm_point%ion_fraction(1) = eos_res%ion_fraction(1)
      atm_point%ion_fraction(2) = eos_res%ion_fraction(2)
      atm_point%ion_fraction(3) = eos_res%ion_fraction(3)
      pt_scr%qtl = log10_temperature
      pt_scr%qt = exp(ln10*log10_temperature)
      pt_scr%qpl = log10_pressure
      pt_scr%qp = exp(ln10*log10_pressure)
      pt_scr%qdl = eos_res%log10_density
      pt_scr%qd = exp(ln10*eos_res%log10_density)
      pt_scr%qo = kap_res(i_kap)
      pt_scr%qol = kap_res(i_log10_kap)
      pt_scr%qqdp = eos_res%dlnrho_dlnp
      pt_scr%qqod = kap_res(i_dlnkap_dlnrho)
      pt_scr%qqot = kap_res(i_dlnkap_dlnt)
      pt_scr%qdel = 0.0d0
      pt_scr%qqdt = eos_res%dlnrho_dlnt
      pt_scr%qdela = eos_res%adiabatic_gradient
      pt_scr%qqcp = eos_res%specific_heat_cp
      pt_scr%qrmu = eos_res%specific_gas_constant
      pt_scr%qemu = eos_res%electron_mean_weight_inverse

! The arguments this integrand does not use stay in the shared argument
! list; this dead statement keeps the compiler from warning about them.
      if (.false.) print *, luminosity_linear, temperature_rotation_factor, conductive_opacity_flag, log10_radius

      return
end subroutine compute_atmosphere_derivs

! ---------------------------------------------------------------
! compute_envelope_derivs (legacy qenv.f)
! ---------------------------------------------------------------
! The right-hand side of the envelope integration, with x = log P and
! y = (log(m/M), log T, log r), M the total mass (star%stotal = log M):
!   d log m / d log P = -4 pi r^4 P / (G m^2 f_P)  (hydrostatic equilibrium)
!   d log T / d log P = grad                       (the temperature gradient)
!   d log r / d log P = -r P / (G m rho f_P)       (mass continuity with
!                                                   hydrostatic equilibrium)
! rho comes from the EOS at (P, T), and grad, radiative or convective,
! from compute_temperature_gradients (legacy TPGRAD). f_P and f_T
! (pressure_rotation_factor, temperature_rotation_factor; legacy FPL,
! FTL) are the rotational distortion factors, 1 without rotation. Also
! stores the point and its physics in pt_scr for env_struct.
subroutine compute_envelope_derivs(log10_pressure_indep, y, dydx, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, log10_gravity, &
     in_atmosphere, want_derivatives, conductive_opacity_flag, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     saha_state, ierr)

      use convection_lib, only: compute_temperature_gradients
      use eos_lib
      use kap_lib
      use star_info_lib, only: star, json
      use point_scratch_lib
      use physical_constants
      use math_lib
      implicit none

      double precision, intent(in) :: log10_pressure_indep
      double precision, intent(in) :: y(3)
      double precision, intent(out) :: dydx(3)
      double precision, intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(in) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag
      double precision, intent(inout) :: log10_radius
      double precision, intent(in) :: log10_teff
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: saha_state
      integer, intent(out) :: ierr

! --- locals ---
      double precision :: log10_pressure, log10_mass, log10_temperature
! The EOS and opacity results at this point.
      type(eos_result) :: eos_res
      double precision :: kap_res(num_kap_results)
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component
      double precision :: convective_velocity
      logical :: is_convective

! The error code of the gradient call, passed on in ierr.
      integer :: jerr

      ierr = 0
      log10_pressure = log10_pressure_indep
! log m = log(m/M) + log M.
      log10_mass = y(1) + star%stotal
      log10_temperature = y(2)
      log10_radius = y(3)
      call eos_get(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, eos_res, want_derivatives, in_atmosphere, &
           saha_state, ierr=ierr)
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
      dydx(1) = -exp(ln10*(c4pil+4.0d0*log10_radius+log10_pressure-cgl- &
           log10_mass-log10_mass))/pressure_rotation_factor
! d log T / d log P = grad.
      dydx(2) = actual_gradient
! Mass continuity: d log r / d log P = (d log r/d log m)(d log m/d log P)
! = -r P / (G m rho f_P), with the same 1/f_P as dydx(1) (Kippenhahn &
! Thomas 1970). (The legacy qenv multiplied by f_P instead, which broke
! mass continuity by a factor f_P^2 in rotating envelopes; bugsweep
! 2026.)
      dydx(3) = -exp(ln10*(log10_pressure+log10_radius-cgl-log10_mass- &
           eos_res%log10_density))/pressure_rotation_factor
! Store the point for the kernel: the structure variables (log m as
! log(m/M)),
      pt_scr%current_log10_pressure = log10_pressure
      pt_scr%current_log10_temperature = log10_temperature
      pt_scr%current_log10_mass = log10_mass - star%stotal
      pt_scr%current_log10_radius = log10_radius
      pt_scr%current_log10_density = eos_res%log10_density
      pt_scr%current_velocity = convective_velocity
! the gradients (grad_rad, grad_ad, grad), which locate the convection
! zone, and beta,
       pt_scr%current_gradients(1) = radiative_gradient
       pt_scr%current_gradients(2) = eos_res%adiabatic_gradient
       pt_scr%current_gradients(3) = actual_gradient
       pt_scr%current_beta = eos_res%beta ! the gas-pressure fraction
! the ionization fractions and EOS derivatives for the profile,
      pt_scr%current_ion_fraction(1) = eos_res%ion_fraction(1)
      pt_scr%current_ion_fraction(2) = eos_res%ion_fraction(2)
      pt_scr%current_ion_fraction(3) = eos_res%ion_fraction(3)
      pt_scr%qqdp = eos_res%dlnrho_dlnp
      pt_scr%qqdt = eos_res%dlnrho_dlnt
      pt_scr%qqcp = eos_res%specific_heat_cp

! and the opacity and the rest of the q* scratch (output only: no
! physics reads them).
      pt_scr%current_opacity = kap_res(i_kap)
      pt_scr%qtl = log10_temperature
      pt_scr%qt = exp(ln10*log10_temperature)
      pt_scr%qpl = log10_pressure
      pt_scr%qp = exp(ln10*log10_pressure)
      pt_scr%qdl = eos_res%log10_density
      pt_scr%qd = exp(ln10*eos_res%log10_density)
      pt_scr%qo = kap_res(i_kap)
      pt_scr%qol = kap_res(i_log10_kap)
      pt_scr%qqod = kap_res(i_dlnkap_dlnrho)
      pt_scr%qqot = kap_res(i_dlnkap_dlnt)
      pt_scr%qdel = actual_gradient
      pt_scr%qdela = eos_res%adiabatic_gradient
      pt_scr%qrmu = eos_res%specific_gas_constant
      pt_scr%qemu = eos_res%electron_mean_weight_inverse

! The arguments this integrand does not use stay in the shared argument
! list; this dead statement keeps the compiler from warning about them.
      if (.false.) print *, log10_gravity, conductive_opacity_flag

      return
end subroutine compute_envelope_derivs

end module envelope_integrands

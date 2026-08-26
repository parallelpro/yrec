!----------------------------------------------------------------------
! hra
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original hra.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Evaluates a fit to the Harvard Reference Atmosphere T-tau relation:
! given the optical depth tau, returns log10(T). Called (as an
! external double precision function) from atmosphere_derivs.f90 and atm_lib.f90.
double precision function harvard_t_tau(optical_depth)
      implicit none

      double precision, intent(in) :: optical_depth

      double precision :: log10_tau, log10_tau2, log10_tau3, log10_tau4, &
           log10_tau5, log10_tau6, log10_tau7, log10_tau8, log10_tau9
      double precision :: fit_value
      log10_tau = log10(optical_depth)
      fit_value = 3.81152046471d0
      fit_value = fit_value + 0.146133736471d0*log10_tau
      log10_tau2 = log10_tau*log10_tau
      fit_value = fit_value + 0.0267719174279d0*log10_tau2
      log10_tau3 = log10_tau2*log10_tau
      fit_value = fit_value - 0.029280655317d0*log10_tau3
      log10_tau4 = log10_tau3*log10_tau
      fit_value = fit_value - 0.0123814456666d0*log10_tau4
      log10_tau5 = log10_tau4*log10_tau
      fit_value = fit_value + 0.00285734990893d0*log10_tau5
      log10_tau6 = log10_tau5*log10_tau
      fit_value = fit_value + 0.0024575213331d0*log10_tau6
      log10_tau7 = log10_tau6*log10_tau
      fit_value = fit_value + 0.000521560431455d0*log10_tau7
      log10_tau8 = log10_tau7*log10_tau
      fit_value = fit_value + 4.71770176883d-5*log10_tau8
      log10_tau9 = log10_tau8*log10_tau
      fit_value = fit_value + 1.58685112637d-6*log10_tau9
      harvard_t_tau = fit_value
      return
end function harvard_t_tau

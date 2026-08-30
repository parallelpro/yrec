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
! in_atmosphere/want_derivatives/conductive_opacity_flag/
! log10_radius/log10_teff/hydrogen_fraction/metal_fraction match the
! sibling routine atmosphere_derivs.f90's naming at the equivalent call-site
! positions (both are called via EXTERNAL from atm_lib.f90's atm_get/BSSTEP).
! pressure_rotation_factor/temperature_rotation_factor (FPL/FTL) match
! temperature_gradients.f90's naming, discovered from that file's own analysis of the
! DELR formula.
subroutine envelope_derivs(log10_pressure_indep, y, dydx, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, log10_gravity, &
     in_atmosphere, want_derivatives, conductive_opacity_flag, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     env_call_count, saha_state, ierr)

      use temperature_gradients_lib
      use eos_lib
      use kap_lib
      use star_info_lib, only: star, json
      use point_scratch_lib
      use phys_const_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log10_pressure_indep
      double precision, intent(inout) :: y(3), dydx(3)
      double precision, intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(in) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag
      double precision, intent(inout) :: log10_radius
      double precision, intent(in) :: log10_teff
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: env_call_count, saha_state
      integer, intent(out) :: ierr

! --- locals ---
      double precision :: log10_pressure, log10_mass, log10_temperature
! 2026 named-index results: the eos/kap relay soup is two arrays
! (fresh each call -- this is an ODE integrand).
      double precision :: eos_res(num_eos_results), kap_res(num_kap_results)
      double precision :: actual_gradient, radiative_gradient, &
           dgrad_dt_component, dgrad_dp_component, dgrad_dr_component
      double precision :: convective_velocity
      logical :: is_convective

      ! 2026 ierr campaign: the callback protocol now carries ierr, so
      ! the gradient failure propagates (the former documented residual).
      integer :: jerr

      ierr = 0
      log10_pressure = log10_pressure_indep
      log10_mass = y(1) + star%stotal
      log10_temperature = y(2)
      log10_radius = y(3)
      call eos_get(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, eos_res, want_derivatives, in_atmosphere, &
           saha_state, ierr=ierr)
      if (ierr /= 0) return
! kap at eqstat's returned density -- the historical inout dataflow
      call kap_get(eos_res(i_log10_density), log10_temperature, &
           hydrogen_fraction, metal_fraction, kap_res, &
           eos_res(i_fxion:i_fxion+2), ierr=ierr)
      if (ierr /= 0) return
      star%iovim = -1
      call temperature_gradients(log10_temperature, log10_pressure, &
           eos_res, kap_res, log10_radius, log10_mass, luminosity_linear, &
           actual_gradient, radiative_gradient, dgrad_dt_component, &
           dgrad_dp_component, dgrad_dr_component, convective_velocity, &
           want_derivatives, is_convective, pressure_rotation_factor, &
           temperature_rotation_factor, log10_teff, jerr)
      if (jerr /= 0) then
         ierr = jerr
         return
      end if
      dydx(1) = -exp(ln10*(c4pil+4.0d0*log10_radius+log10_pressure-cgl- &
           log10_mass-log10_mass))/pressure_rotation_factor
      dydx(2) = actual_gradient
      dydx(3) = -exp(ln10*(log10_pressure+log10_radius-cgl-log10_mass- &
           eos_res(i_log10_density)))*pressure_rotation_factor
      env_call_count = env_call_count + 1
! 07/02 ALWAYS STORE THE BASIC STRUCTURE VARIABLES.
      pt_scr%current_log10_pressure = log10_pressure
      pt_scr%current_log10_temperature = log10_temperature
      pt_scr%current_log10_mass = log10_mass - star%stotal
      pt_scr%current_log10_radius = log10_radius
      pt_scr%current_log10_density = eos_res(i_log10_density)
      pt_scr%current_velocity = convective_velocity
! JVS 08/13 ALWAYS STORE GRADIENTS (FOR TRACKING CZ)
       pt_scr%current_gradients(1) = radiative_gradient
       pt_scr%current_gradients(2) = eos_res(i_grada)
       pt_scr%current_gradients(3) = actual_gradient
       pt_scr%current_beta = eos_res(i_beta) ! added 03/14
! JVS 08/25 ALSO ALWAYS SAVE ADDITIONAL INFO FOR PROFILE
      pt_scr%current_ion_fraction(1) = eos_res(i_fxion)
      pt_scr%current_ion_fraction(2) = eos_res(i_fxion+1)
      pt_scr%current_ion_fraction(3) = eos_res(i_fxion+2)
      pt_scr%qqdp = eos_res(i_dlnrho_dlnp)
      pt_scr%qqdt = eos_res(i_dlnrho_dlnt)
      pt_scr%qqcp = eos_res(i_cp)

! 2026 (.store convergence): these saves were gated on the print
! flag (or the retired pulse-derivative mode), which left
! current_opacity -- and hence the envelope opacity the stitched
! model materializes -- stale (typically zero) whenever the caller
! did not ask for printing. They are output-only scratch (no
! physics reads them), so save unconditionally.
      pt_scr%current_opacity = kap_res(i_kap)
      pt_scr%qtl = log10_temperature
      pt_scr%qt = exp(ln10*log10_temperature)
      pt_scr%qpl = log10_pressure
      pt_scr%qp = exp(ln10*log10_pressure)
      pt_scr%qdl = eos_res(i_log10_density)
      pt_scr%qd = exp(ln10*eos_res(i_log10_density))
      pt_scr%qo = kap_res(i_kap)
      pt_scr%qol = kap_res(i_log10_kap)
      pt_scr%qqod = kap_res(i_dlnkap_dlnrho)
      pt_scr%qqot = kap_res(i_dlnkap_dlnt)
      pt_scr%qdel = actual_gradient
      pt_scr%qdela = eos_res(i_grada)
      pt_scr%qrmu = eos_res(i_gas_constant)
      pt_scr%qemu = eos_res(i_mu_e_inv)

! KC 2025-05-31 THESE MUST BE RETAINED FOR EXTERNAL PROCEDURE COMPATIBILITY.
      if (.false.) print *, log10_gravity, conductive_opacity_flag

      return
end subroutine envelope_derivs

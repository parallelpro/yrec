!----------------------------------------------------------------------
! qatm
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original qatm.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! QATM CALCULATES THE DERIVATIVE D(P)/D(TAU), USING THE EDDINGTON
! APPROXIMATION FOR A T-TAU RELATION.
! IT ALSO RETURNS TL,O, AND FXION FOR OUTPUT PURPOSES
!   Q(TAU) = 0.6666667
subroutine qatm(log10_optical_depth, y, dydx, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, log10_gravity, &
     in_atmosphere, want_derivatives, conductive_opacity_flag, print_flag, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     atm_call_count, saha_state)

      use eos_lib
      use kap_lib
      use atm_table_lib
      use star_info_lib, only: star
      use const_lib
      implicit none
      integer, parameter :: json=5000

      double precision, intent(inout) :: log10_optical_depth
      double precision, intent(inout) :: y(3), dydx(3)
      double precision, intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(in) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag, print_flag
      double precision, intent(in) :: log10_radius
      double precision, intent(inout) :: log10_teff
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: atm_call_count, saha_state

! former common/nwlaol/: not used in this file; declared only to
! preserve layout.
! MHP 8/25 Removed unused variables
!      CHARACTER*256 FLAOL, FPUREZ
! MHP 8/25 Removed character file names from common block
      double precision :: fxion(3)

      save

! --- locals ---
      double precision :: effective_gravity, optical_depth
      double precision :: log10_temperature, temperature, log10_pressure, &
           pressure, log10_density, density
      double precision :: beta, beta_inverse, beta14, specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision :: dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
           adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      double precision :: opacity, log10_opacity, dlnkap_dlnrho, dlnkap_dlnt
      double precision :: ttaul0, ttaul1, yy
      double precision :: hra
      external hra

! EDDINGTON APPROXIMATION
      ttaul0(yy) = log10_teff - 0.031235d0 + 0.25d0*dlog10(yy + cc23)

! KRISHNA-SWAMY APPROXIMATION (BASED ON FIT TO SOLAR ATMOSPHER)
! SEE KRISHNA-SWAMY, AP.J. 1966, 145, 176.
      ttaul1(yy) = log10_teff - 0.031235d0 + 0.25d0*dlog10(yy + &
          1.39d0 - 0.815d0*exp(-2.54d0*yy) - 0.025d0*exp(-30.0d0*yy))

      effective_gravity = dexp(ln10*log10_gravity)*pressure_rotation_factor
      atm_table%atm_tau = log10_optical_depth
      optical_depth = dexp(ln10*atm_table%atm_tau)
! USE KTTAU TO IMPLIMENT FUTURE T TAU RELATIONS
      if (atm_choice .eq. 0) then
            log10_temperature = ttaul0(optical_depth)
      else if (atm_choice .eq. 1) then
            log10_temperature = ttaul1(optical_depth)
      else if (atm_choice .eq. 2) then
            log10_temperature = log10_teff + hra(optical_depth) - atm_hras
      end if
      log10_pressure = y(1)
      call eos_get(log10_temperature,temperature,log10_pressure,pressure, &
           log10_density,density,hydrogen_fraction,metal_fraction,beta, &
           beta_inverse,beta14,fxion,specific_gas_constant, &
           ion_mean_weight_inverse,electron_mean_weight_inverse, &
           electron_degeneracy_parameter,dlnrho_dlnt,dlnrho_dlnp, &
           specific_heat_cp,adiabatic_gradient,dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt,adiabatic_gradient_dt,adiabatic_gradient_dp, &
           specific_heat_cp_dt,specific_heat_cp_dp,want_derivatives, &
           in_atmosphere,saha_state)
      call kap_get(log10_density, log10_temperature, hydrogen_fraction, &
           metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, &
           dlnkap_dlnt, fxion)
      dydx(1) = effective_gravity*optical_depth/(pressure*opacity)
      atm_call_count = atm_call_count + 1
      atm_table%atm_log10_pressure = log10_pressure
      atm_table%atm_log10_temperature = log10_temperature
      if (print_flag .or. star%pulse%lpumod) then
       atm_table%atm_log10_density = log10_density
       atm_table%atm_opacity = opacity
       atm_table%atm_ion_fraction(1) = fxion(1)
       atm_table%atm_ion_fraction(2) = fxion(2)
       atm_table%atm_ion_fraction(3) = fxion(3)
       star%pulse%qtl = log10_temperature
       star%pulse%qt = dexp(ln10*log10_temperature)
       star%pulse%qpl = log10_pressure
       star%pulse%qp = dexp(ln10*log10_pressure)
       star%pulse%qdl = log10_density
       star%pulse%qd = dexp(ln10*log10_density)
       star%pulse%qo = opacity
       star%pulse%qol = log10_opacity
       star%pulse%qqdp = dlnrho_dlnp
       star%pulse%qqed = 0.0d0
       star%pulse%qqod = dlnkap_dlnrho
       star%pulse%qqot = dlnkap_dlnt
       star%pulse%qdel = 0.0d0
       star%pulse%qqdt = dlnrho_dlnt
       star%pulse%qdela = adiabatic_gradient
       star%pulse%qqcp = specific_heat_cp
       star%pulse%qrmu = specific_gas_constant
       star%pulse%qemu = electron_mean_weight_inverse
      endif

! KC 2025-05-31 THESE MUST BE RETAINED FOR EXTERNAL PROCEDURE COMPATIBILITY.
      if (.false.) print *, luminosity_linear, temperature_rotation_factor, conductive_opacity_flag, log10_radius

      return
end subroutine qatm

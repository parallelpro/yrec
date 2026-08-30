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
subroutine atmosphere_derivs(log10_optical_depth, y, dydx, luminosity_linear, &
     pressure_rotation_factor, temperature_rotation_factor, log10_gravity, &
     in_atmosphere, want_derivatives, conductive_opacity_flag, &
     log10_radius, log10_teff, hydrogen_fraction, metal_fraction, &
     atm_call_count, saha_state)

      use eos_lib
      use kap_lib
      use atm_table_lib
      use star_info_lib, only: star, json
      use phys_const_lib
      implicit none

      double precision, intent(inout) :: log10_optical_depth
      double precision, intent(inout) :: y(3), dydx(3)
      double precision, intent(in) :: luminosity_linear, &
           pressure_rotation_factor, temperature_rotation_factor, log10_gravity
      logical, intent(in) :: in_atmosphere, want_derivatives, &
           conductive_opacity_flag
      double precision, intent(in) :: log10_radius
      double precision, intent(inout) :: log10_teff
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      integer, intent(inout) :: atm_call_count, saha_state

! former common/nwlaol/: not used in this file; declared only to
! preserve layout.
! MHP 8/25 Removed unused variables
!      CHARACTER*256 FLAOL, FPUREZ
! MHP 8/25 Removed character file names from common block
! --- locals ---
      double precision :: effective_gravity, optical_depth
      double precision :: log10_temperature, log10_pressure
! 2026 named-index results: the eos/kap relay soup is two arrays
! (fresh each call -- this is an ODE integrand).
      double precision :: eos_res(num_eos_results), kap_res(num_kap_results)
      double precision :: ttaul0, ttaul1, yy
      double precision :: harvard_t_tau
      external harvard_t_tau

! EDDINGTON APPROXIMATION
      ttaul0(yy) = log10_teff - 0.031235d0 + 0.25d0*log10(yy + cc23)

! KRISHNA-SWAMY APPROXIMATION (BASED ON FIT TO SOLAR ATMOSPHER)
! SEE KRISHNA-SWAMY, AP.J. 1966, 145, 176.
      ttaul1(yy) = log10_teff - 0.031235d0 + 0.25d0*log10(yy + &
          1.39d0 - 0.815d0*exp(-2.54d0*yy) - 0.025d0*exp(-30.0d0*yy))

      effective_gravity = exp(ln10*log10_gravity)*pressure_rotation_factor
      atm_table%atm_tau = log10_optical_depth
      optical_depth = exp(ln10*atm_table%atm_tau)
! USE KTTAU TO IMPLIMENT FUTURE T TAU RELATIONS
      if (star%job%atm_choice .eq. 0) then
            log10_temperature = ttaul0(optical_depth)
      else if (star%job%atm_choice .eq. 1) then
            log10_temperature = ttaul1(optical_depth)
      else if (star%job%atm_choice .eq. 2) then
            log10_temperature = log10_teff + harvard_t_tau(optical_depth) - star%atm_hras
      end if
      log10_pressure = y(1)
      call eos_get(log10_temperature, log10_pressure, hydrogen_fraction, &
           metal_fraction, eos_res, want_derivatives, in_atmosphere, &
           saha_state)
! kap at eqstat's returned density -- the historical inout dataflow
      call kap_get(eos_res(i_log10_density), log10_temperature, &
           hydrogen_fraction, metal_fraction, kap_res, &
           eos_res(i_fxion:i_fxion+2))
      dydx(1) = effective_gravity*optical_depth/ &
           (eos_res(i_pressure)*kap_res(i_kap))
      atm_call_count = atm_call_count + 1
      atm_table%atm_log10_pressure = log10_pressure
      atm_table%atm_log10_temperature = log10_temperature
! 2026 (.store convergence): these saves were gated on the print
! flag (or the retired pulse-derivative mode), which meant the
! atmosphere structure the stitched model materializes carried
! stale density/opacity/ionization whenever the caller did not ask
! for printing. They are output-only scratch (no physics reads
! them), so save unconditionally: every integration leaves a fully
! populated atm_table behind.
      atm_table%atm_log10_density = eos_res(i_log10_density)
      atm_table%atm_opacity = kap_res(i_kap)
      atm_table%atm_ion_fraction(1) = eos_res(i_fxion)
      atm_table%atm_ion_fraction(2) = eos_res(i_fxion+1)
      atm_table%atm_ion_fraction(3) = eos_res(i_fxion+2)
      star%pulse%qtl = log10_temperature
      star%pulse%qt = exp(ln10*log10_temperature)
      star%pulse%qpl = log10_pressure
      star%pulse%qp = exp(ln10*log10_pressure)
      star%pulse%qdl = eos_res(i_log10_density)
      star%pulse%qd = exp(ln10*eos_res(i_log10_density))
      star%pulse%qo = kap_res(i_kap)
      star%pulse%qol = kap_res(i_log10_kap)
      star%pulse%qqdp = eos_res(i_dlnrho_dlnp)
      star%pulse%qqod = kap_res(i_dlnkap_dlnrho)
      star%pulse%qqot = kap_res(i_dlnkap_dlnt)
      star%pulse%qdel = 0.0d0
      star%pulse%qqdt = eos_res(i_dlnrho_dlnt)
      star%pulse%qdela = eos_res(i_grada)
      star%pulse%qqcp = eos_res(i_cp)
      star%pulse%qrmu = eos_res(i_gas_constant)
      star%pulse%qemu = eos_res(i_mu_e_inv)

! KC 2025-05-31 THESE MUST BE RETAINED FOR EXTERNAL PROCEDURE COMPATIBILITY.
      if (.false.) print *, luminosity_linear, temperature_rotation_factor, conductive_opacity_flag, log10_radius

      return
end subroutine atmosphere_derivs

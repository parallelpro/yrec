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

! common/nwlaol/: not used in this file; declared only to preserve
! layout. Naming matches getopac.f90.
      double precision :: olaol(12,104,52), oxa(12), ot(52), orho(104), &
           tollaol
      integer :: iolaol, numofxyz, numrho, numt, iopurez
      logical :: llaol, use_pure_z_table
! MHP 8/25 Removed unused variables
!      CHARACTER*256 FLAOL, FPUREZ
! MHP 8/25 Removed character file names from common block
      common/nwlaol/olaol, oxa, ot, orho, tollaol, &
           iolaol, numofxyz, numrho, numt, llaol, use_pure_z_table, iopurez
      double precision :: fxion(3)
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
! common/atmprt/: all used/set here. Naming matches alsurfp.f90.
      double precision :: atm_tau, atm_log10_pressure, &
           atm_log10_temperature, atm_log10_density, atm_opacity, &
           atm_ion_fraction(3)
      common/atmprt/atm_tau, atm_log10_pressure, atm_log10_temperature, &
           atm_log10_density, atm_opacity, atm_ion_fraction
! common/const1/: only ln10/cc23 are used here. Naming matches
! eqburn.f90.
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
      common/const1/ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
! common/atmos/: atm_choice/atm_hras are used here. Naming matches
! putstore.f90.
      double precision :: atm_hras
      integer :: atm_choice, atm_choice_initial
      logical :: use_ttau_relation
      common/atmos/atm_hras, atm_choice, atm_choice_initial, use_ttau_relation
! common/mhd/: only use_mhd_eos is used here. Naming matches mhdtbl.f90.
      logical :: use_mhd_eos
      integer :: unit_zams_a, unit_zams_b, unit_zams_c, unit_centre1, &
           unit_centre2, unit_centre3, unit_centre4, unit_centre5
      common/mhd/use_mhd_eos, unit_zams_a, unit_zams_b, unit_zams_c, &
           unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

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
      atm_tau = log10_optical_depth
      optical_depth = dexp(ln10*atm_tau)
! USE KTTAU TO IMPLIMENT FUTURE T TAU RELATIONS
      if (atm_choice .eq. 0) then
            log10_temperature = ttaul0(optical_depth)
      else if (atm_choice .eq. 1) then
            log10_temperature = ttaul1(optical_depth)
      else if (atm_choice .eq. 2) then
            log10_temperature = log10_teff + hra(optical_depth) - atm_hras
      end if
      log10_pressure = y(1)
! IF LMHD USE MHD EQUATION OF STATE.
      if (use_mhd_eos)then
       call meqos(log10_temperature,temperature,log10_pressure,pressure, &
            log10_density,density,hydrogen_fraction,metal_fraction,beta, &
            beta_inverse,beta14,fxion,specific_gas_constant, &
            ion_mean_weight_inverse,electron_mean_weight_inverse, &
            electron_degeneracy_parameter,dlnrho_dlnt,dlnrho_dlnp, &
            specific_heat_cp,adiabatic_gradient,dlnrho_dlnt_dt, &
            dlnrho_dlnp_dt,adiabatic_gradient_dt,adiabatic_gradient_dp, &
            specific_heat_cp_dt,specific_heat_cp_dp)
      else
       call eqstat(log10_temperature,temperature,log10_pressure,pressure, &
            log10_density,density,hydrogen_fraction,metal_fraction,beta, &
            beta_inverse,beta14,fxion,specific_gas_constant, &
            ion_mean_weight_inverse,electron_mean_weight_inverse, &
            electron_degeneracy_parameter,dlnrho_dlnt,dlnrho_dlnp, &
            specific_heat_cp,adiabatic_gradient,dlnrho_dlnt_dt, &
            dlnrho_dlnp_dt,adiabatic_gradient_dt,adiabatic_gradient_dp, &
            specific_heat_cp_dt,specific_heat_cp_dp,want_derivatives, &
            in_atmosphere,saha_state)
      endif
      call getopac(log10_density, log10_temperature, hydrogen_fraction, &
           metal_fraction, opacity, log10_opacity, dlnkap_dlnrho, &
           dlnkap_dlnt, fxion)
      dydx(1) = effective_gravity*optical_depth/(pressure*opacity)
      atm_call_count = atm_call_count + 1
      atm_log10_pressure = log10_pressure
      atm_log10_temperature = log10_temperature
      if (print_flag .or. lpumod) then
       atm_log10_density = log10_density
       atm_opacity = opacity
       atm_ion_fraction(1) = fxion(1)
       atm_ion_fraction(2) = fxion(2)
       atm_ion_fraction(3) = fxion(3)
       qtl = log10_temperature
       qt = dexp(ln10*log10_temperature)
       qpl = log10_pressure
       qp = dexp(ln10*log10_pressure)
       qdl = log10_density
       qd = dexp(ln10*log10_density)
       qo = opacity
       qol = log10_opacity
       qqdp = dlnrho_dlnp
       qqed = 0.0d0
       qqod = dlnkap_dlnrho
       qqot = dlnkap_dlnt
       qdel = 0.0d0
       qqdt = dlnrho_dlnt
       qdela = adiabatic_gradient
       qqcp = specific_heat_cp
       qrmu = specific_gas_constant
       qemu = electron_mean_weight_inverse
      endif

! KC 2025-05-31 THESE MUST BE RETAINED FOR EXTERNAL PROCEDURE COMPATIBILITY.
      if (.false.) print *, luminosity_linear, temperature_rotation_factor, conductive_opacity_flag, log10_radius

      return
end subroutine qatm

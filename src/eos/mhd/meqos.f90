!----------------------------------------------------------------------
! meqos
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original meqos.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Mihalas, Hummer, Dappen equation of state. Implemented by YC,
! March 7, 1990. Given log10(T), log10(P), X, and Z, calls the MHD
! table interpolator (mhdpx, converted separately; out of scope for
! this pass) and converts its output into the same thermodynamic
! derivative set used by the rest of the EOS machinery.
subroutine meqos(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, ierr)

! LATMO,KSAHA NEEDED FOR EQSAHA
      use mhd_eos_lib
      use star_info_lib
      use luout_lib
      implicit none

      double precision, intent(in) :: log10_temperature, log10_pressure
      double precision, intent(out) :: temperature, pressure, &
           log10_density, density
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      double precision, intent(out) :: beta, beta_inverse, beta14
      double precision, intent(out) :: ion_fraction(3)
      double precision, intent(out) :: specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp


      integer, parameter :: ivarx = 25
      double precision, parameter :: cnvs = 0.434294481d0
      double precision, parameter :: zero = 0.0d0
      integer :: ier_flag
      double precision :: mhdpx_r10
      double precision :: chi_rho, chi_t, log10_specific_heat_cp
      double precision :: specific_gas_constant_check
      integer :: ion_idx

      integer, intent(out) :: ierr

      ierr = 0

      ier_flag = 0
      temperature = 10.0d0**log10_temperature
      pressure = 10.0d0**log10_pressure
      call mhdpx(log10_pressure, log10_temperature, hydrogen_fraction, &
           mhdpx_r10, ierr)
      if (ierr /= 0) return
      log10_density = mhd_eos%mhd_output(1)
      density = 10.0d0**log10_density
      dlnrho_dlnp = 1.0d0/mhd_eos%mhd_output(4)
      chi_rho = mhd_eos%mhd_output(4)
      chi_t = mhd_eos%mhd_output(5)
      dlnrho_dlnt = -chi_t/chi_rho
      log10_specific_heat_cp = mhd_eos%mhd_output(9)
      specific_heat_cp = 10.0d0**log10_specific_heat_cp
      specific_heat_cp_dp = dlnrho_dlnp*mhd_eos%mhd_output(12)
      specific_heat_cp_dt = mhd_eos%mhd_output(13) + mhd_eos%mhd_output(12)*dlnrho_dlnt
      adiabatic_gradient = mhd_eos%mhd_output(8)
      adiabatic_gradient_dp = dlnrho_dlnp*mhd_eos%mhd_output(10)*cnvs/adiabatic_gradient
      adiabatic_gradient_dt = (mhd_eos%mhd_output(11) + mhd_eos%mhd_output(10)*dlnrho_dlnt)* &
           cnvs/adiabatic_gradient
      dlnrho_dlnp_dt = dlnrho_dlnt*(adiabatic_gradient_dp - 1.0d0 + &
           dlnrho_dlnp + specific_heat_cp_dp)
      dlnrho_dlnt_dt = dlnrho_dlnt*(adiabatic_gradient_dt + 1.0d0 + &
           dlnrho_dlnt + specific_heat_cp_dt)
      beta = 10.0d0**(mhd_eos%mhd_output(20) - mhd_eos%mhd_output(2))
      beta_inverse = 1.0d0/beta
      beta14 = 1.0d0 - beta
      do ion_idx = 1, 3
      ion_fraction(ion_idx) = mhd_eos%mhd_output(ion_idx+13)
      end do
      electron_degeneracy_parameter = mhd_eos%mhd_output(18)
      specific_gas_constant = 10.0d0**mhd_eos%mhd_output(20)/(density*temperature)
      call mu(temperature, pressure, density, hydrogen_fraction, &
           metal_fraction, specific_gas_constant_check, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, beta)
      if (abs((specific_gas_constant-specific_gas_constant_check)/ &
           specific_gas_constant).gt.5.0d-7) then
          write(terminal_unit,*)' ERROR(MHD) IN MEAN WEIGHTS ... '
          write(terminal_unit,*) specific_gas_constant, &
               specific_gas_constant_check
        write(terminal_unit,*) 'ERROR (MHD): CHECK MU'
          write(run_log_unit,*)' ERROR(MHD): IN MEAN WEIGHTS ... '
          write(run_log_unit,*) specific_gas_constant, &
               specific_gas_constant_check
        write(run_log_unit,*) 'ERROR (MHD): CHECK MU'
          ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
          ! facades stop when their caller passes no ierr.
          ierr = 1
          return
      end if
      return
!   999 CONTINUE
      write(run_log_unit,*) 'ERROR(MHD):... MHD TABLE FAIL'
      ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
      ! facades stop when their caller passes no ierr.
      ierr = 1
      return
end subroutine meqos

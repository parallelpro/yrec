!----------------------------------------------------------------------
! oeqos
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original oeqos.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! OPAL (1995) equation of state. Implemented by YC (Feb. 21, 1995).
! Given log10(T), log10(P), X, and Z, looks up density and the
! thermodynamic derivatives needed by the rest of the EOS machinery
! from the OPAL 1995 tables (via rhoofp/esac, converted separately;
! out of scope for this pass).
!
! Other routines in the original file (not converted here, still in
! their own .f files): MU, ESAC, T6RINTERP, READCO, QUAD, GMASS,
! RADSUB, RHOOFP.
subroutine oeqos(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, dlnrho_dlnt, &
     dlnrho_dlnp, specific_heat_cp, adiabatic_gradient, ierr, *)

      use opal_eos_lib
      use phys_const_lib
      use star_info_lib
      use luout_lib
      use math_lib
      implicit none

      double precision, intent(in) :: log10_temperature, log10_pressure
      double precision, intent(out) :: temperature, pressure, &
           log10_density, density
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      double precision, intent(out) :: beta, beta_inverse, beta14, &
           specific_gas_constant, ion_mean_weight_inverse, &
           electron_mean_weight_inverse, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient







      integer, parameter :: ivarx = 25
      double precision, parameter :: cnvs = 0.434294481d0
      double precision, parameter :: zero = 0.0d0
      double precision :: t_million_k, p_e12
      double precision :: hydrogen_fraction_work, metal_fraction_table
      double precision :: density_cgs
      double precision :: specific_gas_constant_check
      integer :: rad_flag, deriv_order
      double precision, external :: rhoofp

! CA3=2.5214D-15

      integer, intent(out) :: ierr

      ierr = 0

      temperature = exp10(log10_temperature)
      pressure = exp10(log10_pressure)
      t_million_k = temperature/1.0d6
      p_e12 = pressure/1.0d12
      if (t_million_k.lt.0.0050d0 .or. t_million_k.gt.100.0d0) return 1
      hydrogen_fraction_work = hydrogen_fraction
      metal_fraction_table = metal_fraction

      rad_flag = 1
      deriv_order = 10

      density_cgs = rhoofp(hydrogen_fraction_work, t_million_k, p_e12, rad_flag, ierr)
      if (ierr /= 0) return
      if (density_cgs.le.-998.0d0) then
         return 1
      end if
      density = density_cgs
      log10_density = log10(density)

      call esac(hydrogen_fraction_work, t_million_k, density_cgs, &
           deriv_order, rad_flag, ierr, *999)
      if (ierr /= 0) return

      if (abs((p_e12-opal_eos%eos_output(1))/p_e12).gt.0.5d-6) then
         write(run_log_unit,*) p_e12, opal_eos%eos_output(1)
         ! 2026 (ROADMAP.md stage 3): stop converted to ierr; the eos_lib
         ! facades stop when their caller passes no ierr.
         ierr = 1
         return
      end if
      dlnrho_dlnp = 1.0d0/opal_eos%eos_output(6)
      dlnrho_dlnt = -opal_eos%eos_output(7)/opal_eos%eos_output(6)

      specific_heat_cp = 1.0d6*opal_eos%eos_output(5)*opal_eos%eos_output(8)/opal_eos%eos_output(6)
      adiabatic_gradient = 1.0d0/opal_eos%eos_output(9)

      beta14 = (2.521971383d-3*t_million_k*t_million_k)* &
           (t_million_k*t_million_k/p_e12)
      beta = 1.0d0 - beta14
      beta_inverse = 1.0d0/beta
      specific_gas_constant = pressure*beta/(density*temperature)
      call mu(temperature, pressure, density, hydrogen_fraction, &
           metal_fraction, specific_gas_constant_check, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, beta)
      if (electron_mean_weight_inverse.le.0.0d0) then
      electron_mean_weight_inverse = 0.0d0
      ion_mean_weight_inverse = specific_gas_constant/gas_constant - &
           electron_mean_weight_inverse
      end if

      return
  999 write(run_log_unit, *) 'WARNNING... OPAL TBL FAIL'

      return 1
end subroutine oeqos

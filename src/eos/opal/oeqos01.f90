!----------------------------------------------------------------------
! oeqos01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original oeqos01.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 2001 OPAL equation of state. LLP, August 10, 2003. Given log10(T),
! log10(P), X, and Z, looks up density and the thermodynamic
! derivatives needed by the rest of the EOS machinery from the OPAL
! 2001 tables (via rhoofp01/esac01, converted separately; out of
! scope for this pass).
!
! Other routines associated with the 2001 OPAL EOS (not converted
! here, still in their own .f files): ESAC01, T6RINTEOS01,
! READCOEOS01, QUADEOS01, GMASS01, RADSUB01, RHOOFP01, EQBOUND01.
! Subroutine MU is also used (shared with the 1995/2006 OPAL EOS).
subroutine oeqos01(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, dlnrho_dlnt, &
     dlnrho_dlnp, specific_heat_cp, adiabatic_gradient, ierr, *)

      use opal_eos_lib
      use const_lib
      use star_info_lib
      use luout_lib
      implicit none

      double precision, intent(in) :: log10_temperature, log10_pressure
      double precision, intent(out) :: temperature, pressure, &
           log10_density, density
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      double precision, intent(out) :: beta, beta_inverse, beta14, &
           specific_gas_constant, ion_mean_weight_inverse, &
           electron_mean_weight_inverse, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 191
      integer, parameter :: ivarx = 25
      double precision, parameter :: cnvs = 0.434294481d0
      double precision, parameter :: zero = 0.0d0
      double precision :: t_million_k, p_e12
      double precision :: hydrogen_fraction_work, metal_fraction_table
      double precision :: density_cgs
      double precision :: specific_gas_constant_check
      integer :: rad_flag, deriv_order
      double precision, external :: rhoofp01

      integer, intent(out) :: ierr

      ierr = 0

      deriv_order = 9  ! gives all 1st and 2nd order data. See instructions
!                  in esac01.
!     NOTE: rad_flag=0 does not add radiation; rad_flag=1 adds radiation
      rad_flag = 1     ! does add radiation  corrections

      temperature = 10.0d0**log10_temperature
      pressure = 10.0d0**log10_pressure
      t_million_k = temperature/1.0d6
      p_e12 = pressure/1.0d12
      if (t_million_k.lt.0.0020d0 .or. t_million_k.gt.100.0d0) return 1
      hydrogen_fraction_work = hydrogen_fraction
      metal_fraction_table = metal_fraction

      density_cgs = rhoofp01(hydrogen_fraction_work, t_million_k, p_e12, &
           rad_flag, ierr)
      if (ierr /= 0) return
      if (density_cgs.le.-998.0d0) then
         return 1
      end if
      density = density_cgs
      log10_density = dlog10(density)

      call esac01(hydrogen_fraction_work, t_million_k, density_cgs, &
           deriv_order, rad_flag, ierr, *999)
      if (ierr /= 0) return

!     *         P12,EOS(1),T6,R,X,ZTAB
!         STOP ' ERROR IN OEQOS01 PTOT'
!      ENDIF
!      QDP=1.0D0/EOS(6)
      dlnrho_dlnp = 1.0d0/opal_eos%eos_output_01(5)
!      QDT= -EOS(7)/ EOS(6)
      dlnrho_dlnt = -opal_eos%eos_output_01(6)/opal_eos%eos_output_01(5)

!      QCP=1.0D6*EOS(5)*EOS(8)/EOS(6)
      specific_heat_cp = 1.0d6*opal_eos%eos_output_01(4)*opal_eos%eos_output_01(7)/opal_eos%eos_output_01(5)
!      DELA=1.0D0/EOS(9)
      adiabatic_gradient = 1.0d0/opal_eos%eos_output_01(8)

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
  999 write(short_file_unit, *) 'WARNNING... OPAL TBL FAIL'

      return 1
end subroutine oeqos01

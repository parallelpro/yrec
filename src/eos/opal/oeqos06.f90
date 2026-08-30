!----------------------------------------------------------------------
! oeqos06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original oeqos06.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! 2006 OPAL equation of state. LLP, October 17, 2006. Given log10(T),
! log10(P), X, and Z, looks up density and the thermodynamic
! derivatives needed by the rest of the EOS machinery from the OPAL
! 2006 tables (via rhoofp06/esac06, converted separately; out of
! scope for this pass).
!
! Other routines associated with the 2006 OPAL EOS (not converted
! here, still in their own .f files): ESAC06, T6RINTEOS06,
! READCOEOS06, QUADEOS06, GMASS06, RADSUB06, RHOOFP06, EQBOUND06.
! Subroutine MU is also used (shared with the 1995/2001 OPAL EOS).
subroutine oeqos06(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, dlnrho_dlnt, &
     dlnrho_dlnp, specific_heat_cp, adiabatic_gradient, ierr, *)

      use opal_eos_lib
      use phys_const_lib
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

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 197
      integer, parameter :: ivarx = 25
      double precision, parameter :: cnvs = 0.434294481d0
      double precision, parameter :: zero = 0.0d0
      double precision :: t_million_k, p_e12
      double precision :: hydrogen_fraction_work, metal_fraction_table
      double precision :: density_cgs
      double precision :: specific_gas_constant_check
      integer :: rad_flag, deriv_order
      double precision, external :: rhoofp06

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
      if (t_million_k.lt.0.001870d0 .or. t_million_k.gt.200.0d0) return 1
      hydrogen_fraction_work = hydrogen_fraction
      metal_fraction_table = metal_fraction

      density_cgs = rhoofp06(hydrogen_fraction_work, t_million_k, p_e12, &
           rad_flag, ierr)
      if (ierr /= 0) return
      if (density_cgs.le.-998.0d0) then
         return 1
      end if
      density = density_cgs
      log10_density = log10(density)

      call esac06(hydrogen_fraction_work, t_million_k, density_cgs, &
           deriv_order, rad_flag, ierr, *999)
      if (ierr /= 0) return

!     *         P12,EOS(1),T6,R,X,ZTAB
!         STOP ' ERROR IN OEQOS06 PTOT'
!      ENDIF

      dlnrho_dlnp = 1.0d0/opal_eos%eos_output_06(6)               ! O2006 EOS(6) is dlogP/dlogRho at const T6
      dlnrho_dlnt = -opal_eos%eos_output_06(7)/opal_eos%eos_output_06(6)       ! O2006 EOS(7) is dlogp/dlogT6 at const Rho
      specific_heat_cp = 1.0d6*opal_eos%eos_output_06(5)*opal_eos%eos_output_06(8)/opal_eos%eos_output_06(6)
                                      ! O2006 EOS(5) is the specific heat. dE/dT6
                                      !              at const Vol
                                      ! O2006 EOS(8) is gamma1
      adiabatic_gradient = 1.0d0/opal_eos%eos_output_06(9)         ! O2006 EOS(9) is gamma2/(gamma2-1)

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
end subroutine oeqos06

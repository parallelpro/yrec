!----------------------------------------------------------------------
! radsub
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original radsub.f; only variable names, source form, and comment
! style were updated.
!
! Adds the radiation-pressure/energy/entropy correction to the OPAL
! 1995 EOS values already in v%eos_output, and recomputes
! the derived thermodynamic derivatives (chi_rho, chi_t6, gamma1,
! gamma2/(gamma2-1), gamma3-1) including that correction. Called from
! esac.f90 when rad_flag=1.
!
! Physical constants (see comments below): na*k = 6.0221367e23 *
! 1.380658e-16 erg/degree K = 8.314511e7 erg/degree K = 8.314511e7 *
! 11604.5 erg/eV = 0.9648575e12 erg/eV. unit_factor = 0.9648575
! (latest physical constants); unit_factor_legacy = 0.9652 (older
! units, still used elsewhere in this EOS code). In these units,
! energy/density is in Mb-cc/gm and pressure is in units of 1e12
! bars = Mb. sigma (Stefan-Boltzmann) = 5.67051e-5 erg/(s cm**2 K**4),
! c = 2.99792458e10 cm/sec; rad_const_over_c = sigma/c, converted to
! units of 1e6 K (T6) via *1e24.
subroutine radsub(v, t6_temperature, density, total_moles, &
     mean_molecular_weight)

      use opal_eos_lib
      implicit none

      type(opal_eos_vintage95), intent(inout) :: v

      double precision, intent(in) :: t6_temperature, density
      double precision, intent(in) :: total_moles
      double precision, intent(in) :: mean_molecular_weight

      double precision :: rad_const_over_c, unit_factor, unit_factor_legacy, &
           molar_gas_constant_mbcc
      data unit_factor/0.9648575d0/, unit_factor_legacy/0.9652d0/, &
           rad_const_over_c/1.8914785d-3/, molar_gas_constant_mbcc/83.1446304d0/
! --- locals ---
      double precision :: moles_per_ev, radiation_pressure, &
           radiation_energy, radiation_entropy
      double precision :: total_pressure, total_energy, total_entropy
      double precision :: chi_rho, chi_t6
      double precision :: molar_specific_heat, gamma3_minus1, gamma1, &
           gamma2_over_gamma2_minus1
      double precision :: de_drho_at_t
      double precision :: unit_ratio

      moles_per_ev = total_moles*molar_gas_constant_mbcc  ! Mb-cc/unit T6

      radiation_pressure = 4.0d0/3.0d0*rad_const_over_c*t6_temperature**4   ! Mb
      radiation_energy = 3.0d0*radiation_pressure/density   ! Mb-cc/gm
      radiation_entropy = 4.0d0/3.0d0*radiation_energy/t6_temperature   ! Mb-cc/(gm-unit T6)
      total_pressure = v%eos_output(v%eos_index_inverse(i_opal_p)) + radiation_pressure
      total_energy = v%eos_output(v%eos_index_inverse(i_opal_e)) + radiation_energy
      total_entropy = v%eos_output(v%eos_index_inverse(i_opal_s)) + radiation_entropy
      chi_rho = v%eos_output(v%eos_index_inverse(i_opal_chi_rho))*v%eos_output(v%eos_index_inverse(i_opal_p))/ &
           total_pressure
      chi_t6 = (v%eos_output(v%eos_index_inverse(i_opal_p))*v%eos_output(v%eos_index_inverse(i_opal_chi_t)) &
           + 4.0d0*radiation_pressure)/total_pressure
      molar_specific_heat = (v%eos_output(v%eos_index_inverse(i_opal_cv))*moles_per_ev/ &
           mean_molecular_weight + 4.0d0*radiation_energy/t6_temperature)
      gamma3_minus1 = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1 = chi_rho + chi_t6*gamma3_minus1
      gamma2_over_gamma2_minus1 = gamma1/gamma3_minus1

!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)
!    x  /molenak
      de_drho_at_t = v%eos_output(v%eos_index_inverse(i_opal_dedrho)) - radiation_energy/density
      unit_ratio = unit_factor/unit_factor_legacy
      total_pressure = total_pressure*unit_ratio
! MHP 10/02 EN is never used; should this be ET=ET*REVISE???
!      EN=EN*REVISE
      total_entropy = total_entropy*unit_ratio
! DEDRHOA is never used; should this be DEDRHOAT=DEDRHOAT*REVISE????
!      DEDRHOA=DEDRHOA*REVISE
      v%eos_output(v%eos_index_inverse(i_opal_p)) = total_pressure
      v%eos_output(v%eos_index_inverse(i_opal_e)) = total_energy
      v%eos_output(v%eos_index_inverse(i_opal_s)) = total_entropy
      v%eos_output(v%eos_index_inverse(i_opal_dedrho)) = de_drho_at_t
      v%eos_output(v%eos_index_inverse(i_opal_cv)) = molar_specific_heat
      v%eos_output(v%eos_index_inverse(i_opal_chi_rho)) = chi_rho
      v%eos_output(v%eos_index_inverse(i_opal_chi_t)) = chi_t6
      v%eos_output(v%eos_index_inverse(i_opal_gamma1)) = gamma1
      v%eos_output(v%eos_index_inverse(i_opal_gamma2_ratio)) = gamma2_over_gamma2_minus1
      v%eos_output(v%eos_index_inverse(i_opal_gamma3m1)) = gamma3_minus1
      return
end subroutine radsub

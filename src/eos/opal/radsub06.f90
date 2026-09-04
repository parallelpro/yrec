!----------------------------------------------------------------------
! radsub06
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original radsub06.f; only variable names, source form, and comment
! style were updated.
!
! OPAL 2006 EOS analogue of radsub01.f90 (see there for the general
! "norad, then diff in the gammas" strategy), but restructured so the
! rad_flag/irad test is a real IF around the radiation-correction
! block (radsub.f90/radsub01.f90 always compute the radiation-added
! values; the caller, esac06.f90, currently only calls this routine
! when rad_flag=1, so rad_flag is effectively always nonzero here in
! practice, but the guard is preserved as written). total_pressure,
! total_energy, total_entropy, de_drho_at_t, chi_rho, chi_t6, and
! molar_specific_heat are computed once "without radiation" and then,
! if rad_flag/=0, OVERWRITTEN in place with the "with radiation"
! values -- unlike radsub01.f90 there is no separate set of "_norad"
! variables for these seven; only the three gammas get separate
! "_norad" values (needed for the later difference). Also note only
! gamma1 and gamma2_over_gamma2_minus1 get the "add the difference"
! treatment here -- there is no gamma3_minus1 slot in the OPAL 2006
! table (9 variables, not 10; see esac06.f90's header), and the
! corresponding line is commented out in the original.
subroutine radsub06(rad_flag, t6_temperature, density, total_moles, &
     mean_molecular_weight)

      use opal_eos_lib
      implicit none

      integer, intent(in) :: rad_flag
      double precision, intent(in) :: t6_temperature, density
      double precision, intent(in) :: total_moles
      double precision, intent(in) :: mean_molecular_weight

      double precision :: rad_const_over_c, molar_gas_constant_mbcc
      data rad_const_over_c/1.8914785d-3/, molar_gas_constant_mbcc/83.14510d0/
! --- locals ---
      double precision :: rat, moles_per_ev
      double precision :: radiation_pressure, radiation_energy, radiation_entropy
      double precision :: total_pressure, total_energy, total_entropy
      double precision :: de_drho_at_t, chi_rho, chi_t6, molar_specific_heat
      double precision :: gamma3_minus1_norad, gamma1_norad, &
           gamma2_over_gamma2_minus1_norad
      double precision :: gamma3_minus1, gamma1, gamma2_over_gamma2_minus1

      rat = rad_const_over_c

      moles_per_ev = total_moles*molar_gas_constant_mbcc  ! Mb-cc/unit T6

!-----Calculate EOS without radiation correction

      total_pressure = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1))
      total_energy = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(2))
      total_entropy = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(3))
      de_drho_at_t = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(4))
      chi_rho = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(6))*opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1))/ &
           total_pressure
      chi_t6 = (opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1))*opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(7)))/ &
           total_pressure
      molar_specific_heat = (opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(5))*moles_per_ev/ &
           mean_molecular_weight)
      gamma3_minus1_norad = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1_norad = chi_rho + chi_t6*gamma3_minus1_norad
      gamma2_over_gamma2_minus1_norad = gamma1_norad/gamma3_minus1_norad
!---- End  no radiation calculation

      if (rad_flag.ne.0) then
!---- Calculate EOS with radiation calculation
         radiation_pressure = 4.0d0/3.0d0*rat*t6_temperature**4   ! Mb
         radiation_energy = 3.0d0*radiation_pressure/density   ! Mb-cc/gm
         radiation_entropy = 4.0d0/3.0d0*radiation_energy/t6_temperature   ! Mb-cc/(gm-unit T6)
         total_pressure = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1)) + radiation_pressure
         total_energy = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(2)) + radiation_energy
         total_entropy = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(3)) + radiation_entropy
         de_drho_at_t = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(4)) - radiation_energy/density
         chi_rho = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(6))*opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1))/ &
              total_pressure
         chi_t6 = (opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1))*opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(7)) &
              + 4.0d0*radiation_pressure)/total_pressure
         molar_specific_heat = (opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(5))*moles_per_ev/ &
              mean_molecular_weight + 4.0d0*radiation_energy/t6_temperature)
         gamma3_minus1 = total_pressure*chi_t6/(molar_specific_heat*density* &
              t6_temperature)                                        ! DIRECT
         gamma1 = chi_rho + chi_t6*gamma3_minus1 !eq 16.16 Landau_Lifshitz (Stat. Mech) ! DIRECT
         gamma2_over_gamma2_minus1 = gamma1/gamma3_minus1             ! DIRECT

!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)
!    x  /molenak
!-----Add difference between EOS with and without radiation.  cvtt
!       calculation is not accurate enough to give accurate results using
!       eq. 16.16 Landau&Lifshitz (SEE line labeled DIRECT)
         opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(8)) = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(8)) + &
              gamma1 - gamma1_norad
         opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(9)) = opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(9)) + &
              gamma2_over_gamma2_minus1 - gamma2_over_gamma2_minus1_norad
!     eos(iri(10))=eos(iri(10))+gam3pt-gam3pt_norad
      end if
!-----End EOS calculations with radiation
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(1)) = total_pressure
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(2)) = total_energy
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(3)) = total_entropy
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(4)) = de_drho_at_t
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(5)) = molar_specific_heat
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(6)) = chi_rho
      opal_eos%eos_output_06(opal_eos%eos_index_inverse_06(7)) = chi_t6
      return
end subroutine radsub06

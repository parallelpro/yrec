!----------------------------------------------------------------------
! radsub01
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original radsub01.f; only variable names, source form, and comment
! style were updated.
!
! OPAL 2001 EOS analogue of radsub.f90, but with a different
! algorithm (preserved verbatim): rather than folding the radiation
! term directly into a single set of derivatives, this computes the
! thermodynamic derivatives BOTH without and with the radiation
! correction, then stores the "with radiation" pressure/energy/
! entropy/specific-heat/chi_rho/chi_t6 directly but only ADDS THE
! DIFFERENCE (with minus without) into gamma1/gamma2_over_gamma2_
! minus1/gamma3_minus1 -- because, per the original's comment, the
! molar_specific_heat calculation isn't accurate enough to give
! accurate gammas directly via the eq. 16.16 Landau & Lifshitz route
! (marked "DIRECT" below) when radiation dominates.
subroutine radsub01(t6_temperature, density, total_moles, &
     mean_molecular_weight)

      implicit none

      integer, parameter :: mx = 5, mv = 10, nr = 169, nt = 191

      double precision, intent(in) :: t6_temperature, density
      double precision, intent(in) :: total_moles
      double precision, intent(in) :: mean_molecular_weight

! common/eeos/: eos_output is the caller-facing interpolated result.
      double precision :: esact, eos_output(mv)
      common/eeos/ esact, eos_output

! common/beos/: eos_index_inverse is used here; the rest are placeholders.
      double precision :: z_table(mx)
      integer :: eos_index_inverse(10), eos_var_order(10), t6_index_lo(nr)
      common/beos/ z_table, eos_index_inverse, eos_var_order, t6_index_lo

      double precision :: rad_const_over_c, molar_gas_constant_mbcc
! NOTE: neither literal has a D-suffix in the original (radsub01.f) --
! preserved verbatim: both are parsed as single-precision constants
! and then widened to double precision, which is NOT bit-identical to
! the correctly-rounded double values.
      data rad_const_over_c/1.8914785e-3/, molar_gas_constant_mbcc/83.14510/

      save

! --- locals ---
      double precision :: rat
      double precision :: moles_per_ev, radiation_pressure, &
           radiation_energy, radiation_entropy
      double precision :: total_pressure, total_energy, total_entropy
      double precision :: chi_rho, chi_t6, molar_specific_heat
      double precision :: gamma3_minus1_norad, gamma1_norad, &
           gamma2_over_gamma2_minus1_norad
      double precision :: gamma3_minus1, gamma1, gamma2_over_gamma2_minus1

      rat = rad_const_over_c

      moles_per_ev = total_moles*molar_gas_constant_mbcc  ! Mb-cc/unit T6
      radiation_pressure = 4.0d0/3.0d0*rat*t6_temperature**4   ! Mb
      radiation_energy = 3.0d0*radiation_pressure/density   ! Mb-cc/gm
      radiation_entropy = 4.0d0/3.0d0*radiation_energy/t6_temperature   ! Mb-cc/(gm-unit T6)

!-----Calculate EOS without radiation correction
      total_pressure = eos_output(eos_index_inverse(1))
      total_energy = eos_output(eos_index_inverse(2))
      total_entropy = eos_output(eos_index_inverse(3))
      chi_rho = eos_output(eos_index_inverse(5))*eos_output(eos_index_inverse(1))/ &
           total_pressure
      chi_t6 = (eos_output(eos_index_inverse(1))*eos_output(eos_index_inverse(6)))/ &
           total_pressure
      molar_specific_heat = (eos_output(eos_index_inverse(4))*moles_per_ev/ &
           mean_molecular_weight)
      gamma3_minus1_norad = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)
      gamma1_norad = chi_rho + chi_t6*gamma3_minus1_norad
      gamma2_over_gamma2_minus1_norad = gamma1_norad/gamma3_minus1_norad
!---- End  no radiation calculation

!---- Calculate EOS with radiation calculation
      total_pressure = eos_output(eos_index_inverse(1)) + radiation_pressure
      total_energy = eos_output(eos_index_inverse(2)) + radiation_energy
      total_entropy = eos_output(eos_index_inverse(3)) + radiation_entropy
      chi_rho = eos_output(eos_index_inverse(5))*eos_output(eos_index_inverse(1))/ &
           total_pressure
      chi_t6 = (eos_output(eos_index_inverse(1))*eos_output(eos_index_inverse(6)) &
           + 4.0d0*radiation_pressure)/total_pressure
!     gam1t(jcs,i)=(p(jcs,i)*gam1(jcs,i)+4D0/3D0*pr)/pt(jcs,i)
!     gam2pt(jcs,i)=(gam2p(jcs,i)*p(jcs,i)+4D0*pr)/pt(jcs,i)
!     gam3pt(jcs,i)=gam1t(jcs,i)/gam2pt(jcs,i)
      molar_specific_heat = (eos_output(eos_index_inverse(4))*moles_per_ev/ &
           mean_molecular_weight + 4.0d0*radiation_energy/t6_temperature)
      gamma3_minus1 = total_pressure*chi_t6/(molar_specific_heat*density* &
           t6_temperature)                                        ! DIRECT
      gamma1 = chi_rho + chi_t6*gamma3_minus1 !eq 16.16 Landau_Lifshitz (Stat. Mech) ! DIRECT
      gamma2_over_gamma2_minus1 = gamma1/gamma3_minus1             ! DIRECT
!-----End Eos calculations with radiation

!     normalize cvt to 3/2 when gas is ideal,non-degenerate,
!     fully-ionized, and has no radiation correction
!     cvt=(eos(5)*molenak/tmass+4.*er/t6)
!    x  /molenak
      eos_output(eos_index_inverse(1)) = total_pressure
      eos_output(eos_index_inverse(2)) = total_energy
      eos_output(eos_index_inverse(3)) = total_entropy
      eos_output(eos_index_inverse(4)) = molar_specific_heat
      eos_output(eos_index_inverse(5)) = chi_rho
      eos_output(eos_index_inverse(6)) = chi_t6
!-----Add difference between EOS with and without radiation.  cvtt
!       calculation is not accurate enough to give accurate results using
!       eq. 16.16 Landau&Lifshitz (SEE line labeled DIRECT)
      eos_output(eos_index_inverse(7)) = eos_output(eos_index_inverse(7)) + &
           gamma1 - gamma1_norad
      eos_output(eos_index_inverse(8)) = eos_output(eos_index_inverse(8)) + &
           gamma2_over_gamma2_minus1 - gamma2_over_gamma2_minus1_norad
      eos_output(eos_index_inverse(9)) = eos_output(eos_index_inverse(9)) + &
           gamma3_minus1 - gamma3_minus1_norad
      return
end subroutine radsub01

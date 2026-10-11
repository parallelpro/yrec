!----------------------------------------------------------------------
! yale_eos
!----------------------------------------------------------------------
! The Yale equation of state of Prather (1976): analytic gas physics,
! no tables but the Fermi-Dirac integrals.
!     compute_saha_eos           partial ionization: the Saha equations
!                                for H, He and 11 metals, solved by
!                                Newton iteration on the free-electron
!                                number
!     compute_fully_ionized_eos  full ionization, with a degenerate,
!                                relativistic electron gas and the
!                                optional Debye-Hueckel correction
! The Fermi-Dirac table and the Debye-Hueckel state live in
! yale_eos_tables (yale).
module yale_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use yale_eos_tables, only: debye_huckel_coefficient, debye_huckel_nu, &
           debye_huckel_x, debye_huckel_y, debye_huckel_z, &
           debye_huckel_z_total, yale
      use physical_constants, only: clni, cmkh, gas_constant, ln10
      use io_units, only: run_log_unit
      use eos_mixture, only: ix_he, n_mix_species
      use eos_configuration, only: eos_settings
      use math_lib
      implicit none
      private
      public :: compute_saha_eos, compute_fully_ionized_eos, read_fermi_table

contains

!----------------------------------------------------------------------
! compute_saha_eos (formerly eqsaha)
!----------------------------------------------------------------------
! Prather partial-ionization equation of state: solves the coupled
! Saha ionization equations (11 metal/light-element species plus a
! 2-stage helium treatment) by Newton iteration on the mean number of
! free electrons per ion, then (if want_derivatives) differentiates
! the result with respect to log T and log P. Called from blend_eos.f90
! whenever the temperature is below the ionization cutoff, either
! directly (Prather EOS) or ahead of the SCV table lookup (to get
! derivatives SCV itself does not supply).
subroutine compute_saha_eos(saha_mass_fractions, log10_temperature, &
     log10_pressure, r, want_derivatives, in_atmosphere, saha_state)


      real(dp), intent(in) :: saha_mass_fractions(:)
      real(dp), intent(in) :: log10_temperature
      type(eos_result), intent(inout) :: r
      real(dp), intent(in) :: log10_pressure
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state

! --- locals ---
! ionization_temp/helium_ionization_temp_1/helium_ionization_temp_2:
! ionization potentials of the 11 tabulated species (plus He I/He II)
! expressed as temperatures (potential/k), Kelvin.
! saha_weight_term: ln of the statistical-weight ratio (ion/neutral)
! term entering each species' Saha constant.
! saha_ratio: the Saha equilibrium constant K_i = n_i+ n_e / n_i0 for
! each of the 11 species (beta_inverse-scaled, exponentiated below).
! species_ion_fraction: fractional ionization of each of the 11
! species. species_weighted_term/species_temp_deriv_term: mass- and
! electron-weighted terms entering the R3/derivative sums.
! mean_electrons_per_ion: E, the quantity being Newton-iterated on
! (mean number of free electrons per original particle, referenced to
! ion_mean_weight_inverse); max_electrons_per_ion caps it.
! helium_saha_ratio_1/2, helium_ion_fraction_1/2,
! helium_neutral_fraction: the analogous quantities for the 2-stage
! (He I / He II) helium ionization treatment, handled separately from
! the 11-species loop because helium's 2nd ionization couples directly
! into the electron balance.
      real(dp) :: ionization_temp_over_t(11), saha_ratio(11)
      real(dp) :: helium_mass_fraction, max_electrons_per_ion, &
           temperature_inverse
      real(dp) :: stemp
      integer :: nz1, i, nz0
      logical :: skip_helium_i, skip_helium_ii
      real(dp) :: helium_saha_ratio_1
      real(dp) :: helium_saha_ratio_2
      integer :: isaha
      real(dp) :: ep1, ep1e, ep12, eep1, c33, r3
      real(dp) :: div, temp8
      real(dp) :: species_ion_fraction(11)
      real(dp) :: species_weighted_term(11)
      logical :: converged
      real(dp) :: c22, c23, c32, r2
      real(dp) :: c11, c12, c13, c21, c31, r1
      real(dp) :: cr12, cr13, r1l, cr12l, fact, cr13l
      real(dp) :: cr23
      real(dp) :: delta_electrons_per_ion, deltx1
      real(dp) :: c12l, delx1l, fact1
      real(dp) :: c13l, deltel, fact2
      real(dp) :: deltx2
      real(dp) :: r3t, r3p
      real(dp) :: helium_neutral_fraction, sk0qt, temp, r2t, r2p
      real(dp) :: sk1qt, r1t, r1p
      real(dp) :: sqet, sqep
      real(dp) :: usum
      real(dp) :: sx1qt, sx1qp, sx2qt, sx2qp
      real(dp) :: sx0qt, sx0qp
      real(dp) :: species_temp_deriv_term(11)
      real(dp) :: betaut, ramu, qut
      real(dp) :: pdt, pdtq, qcpi
      real(dp) :: beta16
      real(dp) :: utsum, upsum
      real(dp) :: etemp
      real(dp) :: stemp1, stemp2, stemp3, stemp4, stemp5, &
           stemp6, stemp7
      real(dp) :: sqett, sqetp
      real(dp) :: sx1qtt, sx1qtp, sx2qtt, sx2qtp
      real(dp) :: btemp
      real(dp) :: qqutt, qqutp

! Single-precision literals, as in eqsaha.f: widened to double, they
! differ from the correctly rounded doubles in the last bits. The
! three saved values are the Saha solver's starting guesses, updated
! by each call.
      real(dp), parameter :: ionization_temp(11) = [real(dp) :: 59630., &
           69450., 88710., 91330., 94570., 130630., 157800., 158000., 168630., &
           183080., 250200.]
      real(dp), parameter :: helium_ionization_temp_1 = 285270.
      real(dp), parameter :: helium_ionization_temp_2 = 631370.
      real(dp), parameter :: saha_weight_term(13) = [real(dp) :: -.01, -.47, &
           .60, .49, .08, .11, .0, -.04, .64, .99, 1.03, .6, .0]
      integer, parameter :: nz = 11
      integer, parameter :: max_saha_iterations = 20
      real(dp), parameter :: saha_exponent_tol = 4.0D1
      real(dp), save :: mean_electrons_per_ion = 0.1D0
      real(dp), save :: helium_ion_fraction_1 = 0.0D0
      real(dp), save :: helium_ion_fraction_2 = 0.0D0
      real(dp), parameter :: saha_convergence_tol = 1.0D-09
      helium_mass_fraction = saha_mass_fractions(ix_he)
      max_electrons_per_ion = 1.0d0 + helium_mass_fraction
      temperature_inverse = 1.0d0/r%temperature
! Begin Saha routine
! compute Saha K's
      stemp = 2.50d0*log10_temperature - log10_pressure + cmkh
      nz1 = 1
      do i=1,nz
         ionization_temp_over_t(i) = ionization_temp(i)*temperature_inverse
         saha_ratio(i) = ln10*(saha_weight_term(i) + stemp) - &
              ionization_temp_over_t(i)
         if(saha_ratio(i) < -saha_exponent_tol) exit
         if (saha_ratio(i) <= +saha_exponent_tol) then
         saha_ratio(i) = r%beta_inverse*exp(saha_ratio(i))
         cycle
         end if
         saha_ratio(i) = 1.0d16
         nz1 = i + 1
      end do
      if (i > nz) then
      i = nz + 1
      end if
      nz0 = i - 1
      nz1 = min0(nz,nz1)
      skip_helium_i = .true.
      skip_helium_ii= .true.
      if (helium_mass_fraction >= 1.0d-10) then
      helium_saha_ratio_1 = ln10*(saha_weight_term(12) + stemp) - &
           helium_ionization_temp_1/r%temperature
      if (helium_saha_ratio_1 >= -saha_exponent_tol) then
      skip_helium_i = .false.
      helium_saha_ratio_1 = r%beta_inverse*exp(helium_saha_ratio_1)
      helium_saha_ratio_2 = ln10*(saha_weight_term(13) + stemp) - &
           helium_ionization_temp_2/r%temperature
      if(helium_saha_ratio_2 < -saha_exponent_tol) then
      helium_ion_fraction_2 = 0.0d0
      else
      skip_helium_ii= .false.
      helium_saha_ratio_2 = r%beta_inverse*exp(helium_saha_ratio_2)
      end if
      else
      helium_ion_fraction_1 = 0.0d0
      helium_ion_fraction_2 = 0.0d0
      end if
      else
      helium_ion_fraction_1 = 0.0d0
      helium_ion_fraction_2 = 0.0d0
      end if
! Begin iterations for solution of e
      do isaha=1,max_saha_iterations
         saha_state = saha_state + 1
         ep1 = 1.0d0/(1.0d0 + mean_electrons_per_ion)
         ep1e = ep1/mean_electrons_per_ion
         ep12 = ep1*ep1
         eep1 = mean_electrons_per_ion*ep1
         c33 = -mean_electrons_per_ion*(1.0d0 + mean_electrons_per_ion)
         r3 = mean_electrons_per_ion
         if (nz0 > 0) then
         do i=1,nz0
            div = 1.0d0/(saha_ratio(i) + eep1)
            temp8 = saha_ratio(i)*div
            species_ion_fraction(i) = temp8
            temp8 = temp8*saha_mass_fractions(i)
            r3 = r3 - temp8
            temp8 = temp8*eep1*div
            c33 = c33 - temp8
            species_weighted_term(i) = temp8
         end do
         end if
         c33 = c33*ep1e
         r3 = r3 - helium_mass_fraction*(helium_ion_fraction_1 + &
              helium_ion_fraction_2+helium_ion_fraction_2)
         converged = abs(r3) < saha_convergence_tol
         if (.not. skip_helium_i) then
         c22 = helium_saha_ratio_1 + eep1
         c23 = helium_ion_fraction_1*ep12
         c32 = helium_mass_fraction
         r2 = helium_saha_ratio_1*(1.0d0-helium_ion_fraction_1- &
              helium_ion_fraction_2) - helium_ion_fraction_1*eep1
         if (.not. skip_helium_ii) then
         c11 = eep1
         c12 = -helium_saha_ratio_2
         c13 = helium_ion_fraction_2*ep12
         c21 = helium_saha_ratio_1
         c31 = 2.0d0*helium_mass_fraction
         r1 = helium_ion_fraction_1*helium_saha_ratio_2 - &
              helium_ion_fraction_2*eep1
! Reduction of Saha matrix to diagonal form
         cr12 = - c21/c11
         cr13 = - c31/c11
         c22 = c22 + cr12*c12
         c32 = c32 + cr13*c12
         c23 = c23 + cr12*c13
         c33 = c33 + cr13*c13
         if (r1 /= 0.0d0) then
         r1l=log10(abs(r1))
         if (cr12 /= 0.0d0) then
         cr12l=log10(abs(cr12))
         fact=cr12l+r1l
         if (fact >= -38.0d0) then
         r2 = r2 + cr12*r1
         end if
         end if
         if (cr13 /= 0.0d0) then
         cr13l=log10(abs(cr13))
         fact=cr13l+r1l
         if (fact >= -38.0d0) then
         r3 = r3 + cr13*r1
         end if
         end if
         end if
! Entry for no fully ionized helium (SAHEX2 = 0.0)
         end if
         cr23 = - c32/c22
         c33 = c33 + cr23*c23
         r3 = r3 + cr23*r2
! Entry for neutral helium
         end if
         delta_electrons_per_ion = r3/c33
         if (.not. skip_helium_i) then
         deltx1 = (r2 - c23*delta_electrons_per_ion)/c22
         helium_ion_fraction_1 = helium_ion_fraction_1 + deltx1
         helium_ion_fraction_1 = max(0.0d0,min(1.0d0,helium_ion_fraction_1))
         if (.not. skip_helium_ii) then
! c12*deltx1, floored at 1e-38 in magnitude (written for the DEC-20's
! single-precision range)
         if (c12 /= 0.d0 .and. deltx1 /= 0.d0) then
         c12l=log10(abs(c12))
         delx1l=log10(abs(deltx1))
         fact1=c12l+delx1l
         if (fact1 < -38.0d0) then
         fact1=-38.0d0
         fact1=exp(ln10*fact1)*sign(1.0d0,c12)*sign(1.0d0,deltx1)
         else
         fact1=c12*deltx1
         end if
         else
         fact1=c12*deltx1
         end if
         if (c13 /= 0.d0 .and. delta_electrons_per_ion /= 0.d0) then
         c13l=log10(abs(c13))
         deltel=log10(abs(delta_electrons_per_ion))
         fact2=c13l+deltel
         if (fact2 < -38.0d0) then
         fact2=-38.0d0
         fact2=exp(ln10*fact2)*sign(1.0d0,c13)*sign(1.0d0,delta_electrons_per_ion)
         else
         fact2=c13*delta_electrons_per_ion
         end if
         else
         fact2=c13*delta_electrons_per_ion
         end if
         deltx2=(r1-fact1-fact2)/c11
         helium_ion_fraction_2 = helium_ion_fraction_2 + deltx2
         helium_ion_fraction_2 = max(0.0d0,min(1.0d0,helium_ion_fraction_2))
         end if
         end if
         mean_electrons_per_ion = mean_electrons_per_ion + delta_electrons_per_ion
         mean_electrons_per_ion = max(1.0d-11,min(max_electrons_per_ion, &
              mean_electrons_per_ion))
         converged = converged .and. abs(delta_electrons_per_ion) < saha_convergence_tol
         if(converged) exit
      end do
      if (isaha > max_saha_iterations) then
      write(run_log_unit,99) log10_temperature,log10_pressure, &
           mean_electrons_per_ion,delta_electrons_per_ion
 99   format(' -----SAHA FAILURE (TL,PL)=',2F10.6,'  (E,DE)=',2F20.12)
      end if
! System has been solved for e
      r%electron_mean_weight_inverse = mean_electrons_per_ion*r%ion_mean_weight_inverse
      r%specific_gas_constant = gas_constant*(r%electron_mean_weight_inverse + &
           r%ion_mean_weight_inverse)
      r%density = r%beta*r%pressure/(r%specific_gas_constant*r%temperature)
      r%log10_density = log10(r%density)
      r%ion_fraction(1) = 0.0d0
      if(nz0 >= 7) r%ion_fraction(1) = species_ion_fraction(7)
      r%ion_fraction(2) = helium_ion_fraction_1
      r%ion_fraction(3) = helium_ion_fraction_2
      if (in_atmosphere) then
         return
      end if
! Compute first derivatives
      r3t = 0.0d0
      r3p = 0.0d0
      if (nz0 > 0) then
      do i=nz1,nz0
         r3p = r3p - species_weighted_term(i)
         r3t = r3t -species_weighted_term(i)*(2.5d0 + r%beta14 + &
              ionization_temp_over_t(i))
      end do
      r3p = -r%beta_inverse*r3p
      end if
      if (.not. skip_helium_i) then
      helium_neutral_fraction = 1.0d0 - helium_ion_fraction_1 - helium_ion_fraction_2
      sk0qt = 2.5d0 + r%beta14 + helium_ionization_temp_1*temperature_inverse
      temp = helium_saha_ratio_1*helium_neutral_fraction
      r2t = temp*sk0qt
      r2p = -temp*r%beta_inverse
      if (.not. skip_helium_ii) then
      sk1qt = 2.5d0 + r%beta14 + helium_ionization_temp_2*temperature_inverse
      temp = helium_saha_ratio_2*helium_ion_fraction_1
      r1t = temp*sk1qt
      r1p = -temp*r%beta_inverse
      r2t = r2t + cr12*r1t
      r2p = r2p + cr12*r1p
      r3t = r3t + cr13*r1t
      r3p = r3p + cr13*r1p
      end if
      r3t = r3t + cr23*r2t
      r3p = r3p + cr23*r2p
      end if
      sqet = r3t/c33
      sqep = r3p/c33
      r%dlnrho_dlnt = -1.0d0 - r%beta14 - ep1*sqet
      r%dlnrho_dlnp = r%beta_inverse - ep1*sqep
! Compute internal energy temperature derivative(qut)
      usum = 0.0d0
      if (.not. skip_helium_i) then
      sx1qt = (r2t - c23*sqet)/c22
      sx1qp = (r2p - c23*sqep)/c22
      sx2qt = 0.0d0
      sx2qp = 0.0d0
      if (.not. skip_helium_ii) then
      sx2qt = (r1t - c12*sx1qt - c13*sqet)/c11
      sx2qp = (r1p - c12*sx1qp - c13*sqep)/c11
      end if
      sx0qt = -(sx1qt + sx2qt)
      sx0qp = -(sx1qp + sx2qp)
      usum = helium_mass_fraction*(sx1qt*helium_ionization_temp_1 + &
           sx2qt*(helium_ionization_temp_1 + helium_ionization_temp_2))
      end if
      if (nz0 > 0) then
      stemp = 2.50d0 + r%beta14 - ep1e*sqet
      do i=nz1,nz0
         species_temp_deriv_term(i) = species_weighted_term(i)*(stemp + &
              ionization_temp_over_t(i))
         usum = usum + species_temp_deriv_term(i)*ionization_temp(i)
      end do
      end if
      betaut = 0.75d0*(2.0d0 + r%beta14)
      ramu = gas_constant*r%ion_mean_weight_inverse
      qut =ramu*((1d0+mean_electrons_per_ion)*(betaut+3d0*r%beta_inverse*r%beta14) &
           +betaut*sqet+usum*temperature_inverse)
! Compute specific heat(qcp) and adiabatic gradient(dela)
      pdt = r%pressure/(r%density*r%temperature)
      pdtq = pdt*r%dlnrho_dlnt
      r%specific_heat_cp = qut - pdtq
      qcpi = 1.0d0/r%specific_heat_cp
      r%adiabatic_gradient = -pdtq*qcpi
      if (.not.want_derivatives) then
         return
      end if
! Compute derivatives of qdt,qcp,dela
      beta16 = 4.0d0*r%beta_inverse*r%beta14
      r3t = 0.0d0
      r3p = 0.0d0
      utsum = 0.0d0
      upsum = 0.0d0
      if (nz0 > 0) then
      etemp = (1.0d0+mean_electrons_per_ion+mean_electrons_per_ion)/ &
           (mean_electrons_per_ion*(1.0d0+mean_electrons_per_ion))**2
      stemp1 = 2.5d0 + r%beta14 - sqet*ep1e
      stemp2 = -r%beta_inverse - sqep*ep1e
      stemp3 = etemp*sqep*sqet - r%beta_inverse*r%beta14
      stemp4 = etemp*sqet**2 + beta16
      do i=nz1,nz0
         stemp5 = (stemp1 + ionization_temp_over_t(i))* &
              (1.0d0 - 2.0d0*species_ion_fraction(i))
         stemp6 = species_weighted_term(i)*(stemp4 - ionization_temp_over_t(i)) &
              + species_temp_deriv_term(i)*stemp5
         r3t = r3t - stemp6
         utsum = utsum + ionization_temp(i)*(stemp6 - species_temp_deriv_term(i))
         stemp7 = species_weighted_term(i)*(stemp3 + stemp2*stemp5)
         r3p = r3p - stemp7
         upsum = upsum + ionization_temp(i)*stemp7
      end do
      end if
      if (.not. skip_helium_i) then
      stemp = 2.0d0*helium_ion_fraction_1*ep1*sqet
      r2t = helium_saha_ratio_1*(helium_neutral_fraction*(beta16- &
           helium_ionization_temp_1*temperature_inverse+sk0qt**2)+2d0*sx0qt*sk0qt) &
           +ep12*(stemp*sqet - 2.0d0*sx1qt*sqet)
      r2p = helium_saha_ratio_1*(helium_neutral_fraction*r%beta_inverse* &
           (-r%beta14-sk0qt)+sx0qp*sk0qt-r%beta_inverse* &
           sx0qt) + ep12*(stemp*sqep - sx1qp*sqet - sx1qt*sqep)
      if (.not. skip_helium_ii) then
      stemp = 2.0d0*helium_ion_fraction_2*ep1*sqet
      r1t = helium_saha_ratio_2*(helium_ion_fraction_1*(beta16- &
           helium_ionization_temp_2*temperature_inverse+sk1qt**2)+2d0*sx1qt*sk1qt) &
           +ep12*(stemp*sqet - 2.0d0*sx2qt*sqet)
      r1p = helium_saha_ratio_2*(helium_ion_fraction_1*r%beta_inverse* &
           (-r%beta14-sk1qt)+sx1qp*sk1qt-r%beta_inverse* &
           sx1qt) + ep12*(stemp*sqep - sx2qp*sqet - sx2qt*sqep)
      r2t = r2t + cr12*r1t
      r2p = r2p + cr12*r1p
      r3p = r3p + cr13*r1p
      r3t = r3t + cr13*r1t
      end if
      r3p = r3p + cr23*r2p
      r3t = r3t + cr23*r2t
      end if
      sqett = r3t/c33
      sqetp = r3p/c33
      r%dlnrho_dlnt_dt = -beta16 + (ep1*sqet)**2 - ep1*sqett
      r%dlnrho_dlnt_dp = r%beta_inverse*r%beta14 + ep12*sqet*sqep - ep1*sqetp
      if (nz0 > 0) then
      stemp = ep1e*sqett
      stemp1 = ep1e*sqetp
      do i=nz1,nz0
         utsum = utsum - species_weighted_term(i)*ionization_temp(i)*stemp
         upsum = upsum - species_weighted_term(i)*ionization_temp(i)*stemp1
      end do
      end if
      if (.not. skip_helium_i) then
      sx1qtt = (r2t - c23*sqett)/c22
      sx1qtp = (r2p - c23*sqetp)/c22
      if (.not. skip_helium_ii) then
      sx2qtt = (r1t - c12*sx1qtt - c13*sqett)/c11
      sx2qtp = (r1p - c12*sx1qtp - c13*sqetp)/c11
      utsum = utsum + helium_mass_fraction*(helium_ionization_temp_1+ &
           helium_ionization_temp_2)*(sx2qtt - sx2qt)
      upsum = upsum + helium_mass_fraction*(helium_ionization_temp_1+ &
           helium_ionization_temp_2)*sx2qtp
      end if
      utsum = utsum + helium_mass_fraction*helium_ionization_temp_1*(sx1qtt - sx1qt)
      upsum = upsum + helium_mass_fraction*helium_ionization_temp_1*sx1qtp
      end if
      btemp = (1.0d0+mean_electrons_per_ion)*r%beta14*r%beta_inverse* &
           (24.0d0*r%beta_inverse - 9.0d0)
      qqutt = ramu*(btemp+0.75d0*(2.0d0+r%beta14*(8.0d0*r%beta_inverse+1d0))*sqet+ &
           betaut*sqett + utsum*temperature_inverse)
      qqutp = ramu*(-0.25d0*btemp+0.75d0*((2.0d0+r%beta14*(4.0d0*r%beta_inverse &
           +1.0d0))*sqep-r%beta14*r%beta_inverse*sqet) +betaut*sqetp+upsum*temperature_inverse)
      r%specific_heat_cp_dt = qcpi*(qqutt + pdtq*(1.0d0+r%dlnrho_dlnt) - pdt*r%dlnrho_dlnt_dt)
      r%specific_heat_cp_dp = qcpi*(qqutp - pdtq*(1.0d0-r%dlnrho_dlnp) - pdt*r%dlnrho_dlnt_dp)
      r%adiabatic_gradient_dt = -1.0d0 - r%dlnrho_dlnt - r%specific_heat_cp_dt + &
           r%dlnrho_dlnt_dt/r%dlnrho_dlnt
      r%adiabatic_gradient_dp = 1.0d0 - r%dlnrho_dlnp - r%specific_heat_cp_dp + &
           r%dlnrho_dlnt_dp/r%dlnrho_dlnt

end subroutine compute_saha_eos

!----------------------------------------------------------------------
! compute_fully_ionized_eos (formerly eqrelv)
!----------------------------------------------------------------------
! Fully-ionized-gas equation of state ("relativistic" electron gas via
! the tabulated Fermi-Dirac integrals in yale_eos_tables's
! yale%fermi_table_data, following Prather's thesis). Given log10(T), log10(P), and the mean molecular
! weights, iterates on density until the electron+ion+radiation
! pressure sum matches the given P, then returns density and (if
! in_atmosphere is false, via blend_eos.f90's call) thermodynamic
! derivatives. Called from blend_eos.f90 both for the T >= saha cutoff
! branch (fully ionized gas assumed) and to interpolate against the
! compute_saha_eos/compute_scv_envelope_eos result near the ionization cutoff.
!
!          Qdtt,qdtp = nat-log derivatives of qdt=d(log d )/d(log T )
!          qat,qap = nat-log derivatives of the adiabatic derivative
! DBG 7/92 modified to include Debye-Hueckel correction.
!     Correction is applied to pressure and all thermodynamic variables.
!     Treatment is also valid for electron degeneracy.  A ramp function
!     between the nondegenerate and degernate electron gas cases of the
!     D-H correction is used (this is not correct but, hopefully is
!     close enough.)  The constant term is debye_huckel_coefficient.
subroutine compute_fully_ionized_eos(log10_temperature, log10_pressure, r)


      real(dp), intent(in) :: log10_temperature
      type(eos_result), intent(inout) :: r
      real(dp), intent(in) :: log10_pressure
! --- locals ---
! Scratch algebra kept close to the original short names (rather than
! invented descriptive names) given the density of unlabeled physics
! here; see the inline comments (carried over from the original) for
! what each quantity represents. qpd (d ln P/d ln rho, a byproduct of
! the Newton iteration) is intentionally a separate local from the
! dlnrho_dlnp output argument (d ln rho/d ln P, its reciprocal) --
! they are related but distinct quantities in the original code (QPD
! vs QDP) and must not be aliased to the same name.
      real(dp) :: pinv, beta1, tl8, pl8, pr, ramu, ramut, emul, &
           emul8, ttud, ttcu, ttdcu, dml, almix, cmfdh0, cmfdh1, cmfdh, &
           xx, pr6, pr9, pr8, pr7, cden1, cden2, cden3, dpel2, dx1, dx2, &
           dx3, cl1, cl2, cl3, dcl1, dcl2, dcl3, t0, t1, slope, pel, &
           dpel1, dl8, d8, pa, pe, pdh, pt, pdhp, ptl, factor, pr8l, &
           pr9l, factl, fact, corr, peq, qped, qqpedd, pr1, pr2, pr3, &
           pr4, pr5, qpeqd, qpeqt, qpt, qpd, qqptt, qqpdt, qqpdd, qdpp, &
           ue, ueq, queqt, pd, pdt, ua, ur, udh, u, udhu, ui, pdu, qutd, &
           qudt, qutp, qqutt, qqudt, acp, bcp, qcpi
      integer :: jt1, jt2, jt3, nden, id1p, ii, kk, id1, id2, id3, nn, id
      real(dp) :: df1(5,3), df2(5,3), ff(5,3), ffe(5)

!  The electron pressure (tabulated in Fermi table) for a given (rho,T)
!  is found here; the density is adjusted so that the resulting total
!  pressure (pe + pa + pr) equals the initial P.  For details
!  consult Prather's thesis.
!  Suffix e=electron;R=radiation;a=nucleon;i.e. pa=pressure due to nuclei
      pinv = 1.0d0/r%pressure
      beta1 = 1.0d0 - r%beta
      tl8 = log10_temperature
      pl8 = log10_pressure
      pr = beta1*r%pressure
      ramu = gas_constant*r%ion_mean_weight_inverse
      ramut = ramu*r%temperature
      emul = log10(r%electron_mean_weight_inverse)
      emul8 = emul
!  Pe is tabulated as a function of X(=log d - log mue - 3/2log T) and
!  Y(= log T). Find index(JT1,JT2,JT3) for 3 pt interpolation in Y.
! MHP 10/02 corrected mixed type
      jt1 = int(5.0d0*(log10_temperature - 5.0d0))
      jt1 = max0(1,min0(jt1,18))  ! Force jt1 to be in the table
      jt2 = jt1 + 1
      jt3 = jt1 + 2
      ttud = 5.0d0*(tl8 - (5.0d0 + 0.20d0*real(jt1, dp)))
      ttcu = 0.50d0*ttud*(ttud - 1.0d0)
      ttdcu = ttud - 0.50d0
      dml = r%log10_density + emul
      nden = 1
      id1p = -100
! DBG 7/92 calculate D-H terms: Alpha-mixture, mass-FRAC0, mass-FRAC1
      if (eos_settings%use_debye_huckel_correction) then
         almix = 0.0d0
         do ii=1,18
            almix = almix + debye_huckel_nu(ii)*debye_huckel_z(ii)
         end do
         cmfdh0 = 2.0d0*debye_huckel_x+1.5d0*debye_huckel_y+0.5d0*(debye_huckel_z_total+almix)
         cmfdh1 = debye_huckel_x+debye_huckel_y+0.5d0*almix
      end if
      density_iter: do
! Begin iteration loop for correct density (was label 10)
!  find index (ID1,ID2,ID3) for 3-pt interpolation in X
      xx = dml - 1.50d0*tl8
      pr6= 20.0d0*(xx - yale%fermi_table_x_grid(1)) + 1.0d0
      kk = int(pr6)
      kk = min0(261,max0(1,kk))
      id1 = yale%fermi_table_x_lookup(kk)
!  If index unchanged from previous loop,skip this section
      if (id1 /= id1p) then
      id1p = id1
      id2 = id1 + 1
      id3 = id1 + 2
      pr9 = 1.0d0/( yale%fermi_table_x_grid(id3) - yale%fermi_table_x_grid(id2) )
      pr8 = 1.0d0/( yale%fermi_table_x_grid(id2) - yale%fermi_table_x_grid(id1) )
      pr7 = 1.0d0/( yale%fermi_table_x_grid(id1) - yale%fermi_table_x_grid(id3) )
      cden1 = -pr7*pr8
      cden2 = -pr8*pr9
      cden3 = -pr9*pr7
!  Interpolation in Y is Newtonian(equal spacing in Y)
      do nn = 1,3
       id = id1 + nn - 1
       df1(1,nn) = yale%fermi_table_data(1,id,jt2) - yale%fermi_table_data(1,id,jt1)
       df2(1,nn) = yale%fermi_table_data(1,id,jt3)-yale%fermi_table_data(1,id,jt2)-df1(1,nn)
       ff(1,nn) = yale%fermi_table_data(1,id,jt1) + ttud*df1(1,nn) + ttcu*df2(1,nn)
      end do
      dpel2 = 2.0d0*(ff(1,1)*cden1 + ff(1,2)*cden2 + ff(1,3)*cden3)
      end if
      dx1 = xx - yale%fermi_table_x_grid(id1)
      dx2 = xx - yale%fermi_table_x_grid(id2)
      dx3 = xx - yale%fermi_table_x_grid(id3)
!  Interpolation in X is Lagrangian(unequal table spacing in X)
      cl1 = dx2*dx3*cden1
      cl2 = dx3*dx1*cden2
      cl3 = dx1*dx2*cden3
      dcl1 = (dx2 + dx3)*cden1
      dcl2 = (dx3 + dx1)*cden2
      dcl3 = (dx1 + dx2)*cden3
      if (eos_settings%use_debye_huckel_correction) then
         r%electron_degeneracy_parameter = cl1*yale%fermi_table_eta(id1) + &
              cl2*yale%fermi_table_eta(id2) + cl3*yale%fermi_table_eta(id3)
         if (r%electron_degeneracy_parameter < eos_settings%debye_huckel_eta_min) then
            cmfdh = cmfdh0*sqrt(cmfdh0)
         else if (r%electron_degeneracy_parameter > eos_settings%debye_huckel_eta_max) then
            cmfdh = cmfdh0*sqrt(cmfdh1)
         else
            t0 = cmfdh0*sqrt(cmfdh0)
            t1 = cmfdh0*sqrt(cmfdh1)
            slope = (t1-t0)/(eos_settings%debye_huckel_eta_max-eos_settings%debye_huckel_eta_min)
            cmfdh = t0+(r%electron_degeneracy_parameter-eos_settings%debye_huckel_eta_min)*slope
         end if
      end if
!  Now find pe and its 1st and 2nd derivatives w.r.t. density
      pel =  cl1*ff(1,1) +  cl2*ff(1,2) +  cl3*ff(1,3)
      dpel1 = dcl1*ff(1,1) + dcl2*ff(1,2) + dcl3*ff(1,3)
      dl8 = dml - emul8
      d8 = exp(ln10*dl8)
      pa = ramut*d8
      pe = exp(ln10*pel)
      if (eos_settings%use_debye_huckel_correction) then
         pdh = debye_huckel_coefficient*cmfdh*sqrt(d8/r%temperature)*d8
      end if
      pt = pa + pe + pr
      if (eos_settings%use_debye_huckel_correction) then
         pt = pt + pdh
         pdhp = pdh/pt
      end if
      ptl = log10(pt)
      pr7 = (pa + pe*dpel1)/pt
      if (eos_settings%use_debye_huckel_correction) then
         pr7 = pr7 + 1.5d0*pdhp
      end if
      pr8 = (ptl - pl8)/pr7
      pr9 = ln10*( (pa + pe*(dpel1**2+dpel2/ln10))/(pt*pr7) - pr7 )
      if (eos_settings%use_debye_huckel_correction) then
         pr9 = pr9 + 2.25d0*pdhp/pr7
      end if
! Factor is the difference between the P predicted from the given rho
! and the true P; if too large increment rho and try again.
      factor=abs(pr8)
      if(factor < 1d-15) factor = 1d-15
      pr8l=log10(factor)
      pr9l=log10(abs(pr9))
      factl=log10(0.5d0)+pr9l+2.d0*pr8l
      fact=exp(ln10*factl)
      fact=sign(fact,pr9)
      corr=-pr8-fact
!  limit range of changes in density to a factor of 10**4   llp  1/31/07
      if(corr > 4d0) corr =  4d0
      if(corr < -4d0) corr = -4d0
      if(abs(corr) >= 1.0d-08) then
       dml = dml + corr
       nden = nden + 1
       if(nden <= 20) cycle density_iter
       write(run_log_unit,40) log10_temperature,log10_pressure,ptl,dml,corr
   40    format('EQRELV: Did not Converge: T,P,Pcalc,Dcalc,CORR', &
                4F10.6,F20.12)
       return
      end if
      exit density_iter
      end do density_iter
      r%log10_density = dl8
      r%density = d8
      do kk = 2,5
       do nn = 1,3
          id = id1 + nn - 1
          df1(kk,nn) = yale%fermi_table_data(kk,id,jt2) - yale%fermi_table_data(kk,id,jt1)
          df2(kk,nn)=(yale%fermi_table_data(kk,id,jt3)-yale%fermi_table_data(kk,id,jt2))-df1(kk,nn)
          ff(kk,nn)=yale%fermi_table_data(kk,id,jt1)+ttud*df1(kk,nn)+ttcu*df2(kk,nn)
       end do
       ffe(kk) = cl1*ff(kk,1) + cl2*ff(kk,2) + cl3*ff(kk,3)
      end do
! Derivatives of P
      peq = exp(ln10*ffe(3))
      qped = ffe(4)
      qqpedd =(dcl1*ff(4,1) + dcl2*ff(4,2) + dcl3*ff(4,3))*clni
      pr1 = df1(3,1)+ttdcu*df2(3,1)
      pr2 = df1(3,2)+ttdcu*df2(3,2)
      pr3 = df1(3,3)+ttdcu*df2(3,3)
      pr4 = dcl1*ff(3,1) + dcl2*ff(3,2) + dcl3*ff(3,3)
      pr5 =   5.0d0*(cl1*pr1 + cl2*pr2 + cl3*pr3)
      qpeqd = peq*pr4
      qpeqt = peq*(pr5 - 1.50d0*pr4)
      qpt = (pa + peq)*pinv + 4.0d0*beta1
      if (eos_settings%use_debye_huckel_correction) then
         qpt = qpt - 0.5d0*pdhp
      end if
      qpd = (pa + pe*qped)*pinv
      if (eos_settings%use_debye_huckel_correction) then
         qpd = qpd + 1.5d0*pdhp
      end if
      qqptt = (pa + qpeqt)*pinv + 16.0d0*beta1 - qpt**2
      if (eos_settings%use_debye_huckel_correction) then
         qqptt = qqptt+(0.5d0*qpt+0.25d0)*pdhp
      end if
      qqpdt = (pa + qpeqd)*pinv - qpt*qpd
      if (eos_settings%use_debye_huckel_correction) then
         qqpdt = qqpdt+(0.5d0*qpt-0.75d0)*pdhp
      end if
      qqpdd = (pa + pe*(qped**2 + qqpedd))*pinv - qpd**2
      if (eos_settings%use_debye_huckel_correction) then
        qqpdd = qqpdd + (-1.5d0*qpd+2.25d0)*pdhp
      end if
! Derivatives of density
      r%dlnrho_dlnp = 1.0d0/qpd
      r%dlnrho_dlnt = -qpt*r%dlnrho_dlnp
      pr1 = r%dlnrho_dlnp**2
      qdpp = -qqpdd*r%dlnrho_dlnp*pr1
      r%dlnrho_dlnt_dp = -qpt*qdpp - qqpdt*pr1
      r%dlnrho_dlnt_dt = -qpt*r%dlnrho_dlnt_dp - r%dlnrho_dlnp*(qqptt + r%dlnrho_dlnt*qqpdt)
! Internal energy term(u) and its derivatives
      ue  = exp(ln10*(ffe(2) + emul))
      ueq = exp(ln10*(ffe(5) + emul))
      pr1 = df1(5,1)+ttdcu*df2(5,1)
      pr2 = df1(5,2)+ttdcu*df2(5,2)
      pr3 = df1(5,3)+ttdcu*df2(5,3)
      pr4 = dcl1*ff(5,1) + dcl2*ff(5,2) + dcl3*ff(5,3)
      pr5 =   5.0d0*(cl1*pr1 + cl2*pr2 + cl3*pr3)
      queqt = ueq*(pr5 - 1.50d0*pr4)
      pd = r%pressure/r%density
      pdt = pd/r%temperature
      ua = 1.50d0*ramut
      ur = 3.0d0*beta1*pd
      if (eos_settings%use_debye_huckel_correction) then
         udh = 3.0d0*pdh/d8
      end if
      u = ua + ue + ur
      if (eos_settings%use_debye_huckel_correction) then
         u = u + udh
         udhu = udh/u
      end if
      ui = 1.0d0/u
      pdu = pd*ui
      qutd = (ua + ueq + 4.0d0*ur)*ui
      if(eos_settings%use_debye_huckel_correction) then
         qutd = qutd -0.5d0*udhu
      end if
      qudt = pdu*(1.0d0 - qpt)
      if (eos_settings%use_debye_huckel_correction) then
         qudt = qudt + 0.5d0*udhu
      end if
      qutp = qutd + r%dlnrho_dlnt*qudt
      qqutt = (ua + queqt + 16.0d0*ur)*ui - qutd**2
      if (eos_settings%use_debye_huckel_correction) then
         qqutt = qqutt + (0.5d0*qutd+0.25d0)*udhu
      end if
      qqudt = qudt*(qpt - qutd) - pdu*qqptt
      if (eos_settings%use_debye_huckel_correction) then
         qqudt = qqudt + (0.5d0*qudt-0.25d0)*udhu
      end if
! Specific heat(qcp) and its derivatives
      pr1 =-pdt*r%dlnrho_dlnt
      pr2 = u/r%temperature
      pr3 = -(1.0d0 + r%dlnrho_dlnt) +  r%dlnrho_dlnt_dt/r%dlnrho_dlnt
      pr4 = (1.0d0 - r%dlnrho_dlnp) +  r%dlnrho_dlnt_dp/r%dlnrho_dlnt
      acp = pr1*qpt
      bcp = pr2*qutd
      r%specific_heat_cp = acp + bcp
      qcpi = 1.0d0/r%specific_heat_cp
      r%specific_heat_cp_dt = acp*pr3 + pr1*(qqptt + r%dlnrho_dlnt*qqpdt) + &
              bcp*(qutp - 1.0d0) + pr2*(qqutt + r%dlnrho_dlnt*qqudt)
      r%specific_heat_cp_dp = acp*pr4 + r%dlnrho_dlnp*(pr1*qqpdt + bcp*qudt + pr2*qqudt)
      r%specific_heat_cp_dt = r%specific_heat_cp_dt*qcpi
      r%specific_heat_cp_dp = r%specific_heat_cp_dp*qcpi
! Adiabatic gradient(dela) and derivatives ( all derivs are ln(dela)/ln(
      r%adiabatic_gradient = pr1*qcpi
      r%adiabatic_gradient_dt = pr3 - r%specific_heat_cp_dt
      r%adiabatic_gradient_dp = pr4 - r%specific_heat_cp_dp
      r%electron_degeneracy_parameter = cl1*yale%fermi_table_eta(id1) + &
           cl2*yale%fermi_table_eta(id2) + cl3*yale%fermi_table_eta(id3)

end subroutine compute_fully_ionized_eos

!----------------------------------------------------------------------
! read_fermi_table
!----------------------------------------------------------------------
! Reads the Fermi-Dirac integral table of the degenerate electron gas
! (FERMI.TAB: 43 abscissas, 43 degeneracies, 860 cards of 5 values)
! into yale, and bins the abscissas into fermi_table_x_lookup, the
! interval index compute_fully_ionized_eos looks up. A non-increasing
! abscissa is a corrupt table: ierr = 1.
subroutine read_fermi_table(fermi_table_path, ierr)
      character(len=*), intent(in) :: fermi_table_path
      integer, intent(out) :: ierr
      integer :: fermi_unit, grid_idx, iden_idx, col_idx, card_idx
      integer :: bin_start, bin_width, bin_end

      ierr = 0
      open(newunit=fermi_unit, file=fermi_table_path, form='FORMATTED', status='OLD', &
           iostat=ierr)
      if (ierr /= 0) then
         write(run_log_unit,*) 'eos_init: cannot open the Fermi-Dirac table ', &
              trim(fermi_table_path)
         ierr = 1
         return
      end if
      read(fermi_unit,150)  (yale%fermi_table_x_grid(grid_idx), grid_idx = 1,43)
      read(fermi_unit,150)  (yale%fermi_table_eta(grid_idx), grid_idx = 1,43)
      read(fermi_unit,160) (grid_idx,iden_idx,(yale%fermi_table_data(col_idx,grid_idx,iden_idx), &
           col_idx=1,5),card_idx=1,860)
  150 format(8x,9f8.4)
  160 format(2i5,5f12.7)
      bin_start = 1
      do grid_idx=1,41
! a non-increasing abscissa is a corrupt table
       if (yale%fermi_table_x_grid(grid_idx+1) <= yale%fermi_table_x_grid(grid_idx)) then
          write(run_log_unit,1000) grid_idx
 1000       format(1x,39('>'),40('<')/1x,'read_fermi_table: abscissa',i4, &
           ' does not increase')
          close(fermi_unit)
          ierr = 1
          return
       endif
       bin_width = int((yale%fermi_table_x_grid(grid_idx+1) - yale%fermi_table_x_grid(grid_idx))*20.0d0 + 0.10d0)
       bin_end = bin_start + bin_width - 1
       if (grid_idx == 41)  bin_end = 261
       do iden_idx=bin_start,bin_end
          yale%fermi_table_x_lookup(iden_idx) = grid_idx
       end do
       bin_start = bin_end + 1
      end do

      close(fermi_unit)
end subroutine read_fermi_table

end module yale_eos

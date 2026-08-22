!----------------------------------------------------------------------
! eqsaha
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqsaha.f; only variable names, source form, and comment
! style were updated.
!
! Prather partial-ionization equation of state: solves the coupled
! Saha ionization equations (11 metal/light-element species plus a
! 2-stage helium treatment) by Newton iteration on the mean number of
! free electrons per ion, then (if want_derivatives) differentiates
! the result with respect to log T and log P. Called from eqstat2.f90
! whenever the temperature is below the ionization cutoff, either
! directly (Prather EOS) or ahead of the SCV table lookup (to get
! derivatives SCV itself does not supply).
subroutine eqsaha(saha_mass_fractions, log10_temperature, temperature, &
     log10_pressure, pressure, log10_density, density, beta, &
     beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     want_derivatives, in_atmosphere, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, saha_state)

      use luout_lib
      use const_lib
      implicit none

      double precision, intent(in) :: saha_mass_fractions(12)
      double precision, intent(inout) :: log10_temperature
      double precision, intent(in) :: temperature
      double precision, intent(inout) :: log10_pressure
      double precision, intent(in) :: pressure
      double precision, intent(out) :: log10_density, density
      double precision, intent(in) :: beta, beta_inverse, beta14
      double precision, intent(out) :: ion_fraction(3), &
           specific_gas_constant
      double precision, intent(in) :: ion_mean_weight_inverse
      double precision, intent(out) :: electron_mean_weight_inverse
      logical, intent(in) :: want_derivatives, in_atmosphere
      double precision, intent(out) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
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
      double precision :: ionization_temp(11), helium_ionization_temp_1, &
           helium_ionization_temp_2
      double precision :: saha_weight_term(13)
      double precision :: ionization_temp_over_t(11), saha_ratio(11)
      integer :: nz
      integer :: max_saha_iterations
      double precision :: saha_exponent_tol
      double precision :: mean_electrons_per_ion
      double precision :: saha_convergence_tol
      double precision :: helium_mass_fraction, max_electrons_per_ion, &
           temperature_inverse
      double precision :: stemp
      integer :: nz1, i, nz0
      logical :: skip_helium_i, skip_helium_ii
      double precision :: helium_saha_ratio_1
      double precision :: helium_saha_ratio_2
      integer :: isaha
      double precision :: ep1, ep1e, ep12, eep1, c33, r3
      double precision :: div, temp8
      double precision :: species_ion_fraction(11)
      double precision :: species_weighted_term(11)
      logical :: converged
      double precision :: c22, c23, c32, r2
      double precision :: c11, c12, c13, c21, c31, r1
      double precision :: cr12, cr13, r1l, cr12l, fact, cr13l
      double precision :: cr23
      double precision :: delta_electrons_per_ion, deltx1
      double precision :: c12l, delx1l, fact1
      double precision :: c13l, deltel, fact2
      double precision :: deltx2
      double precision :: helium_ion_fraction_1, helium_ion_fraction_2
      double precision :: rmub, r3t, r3p
      double precision :: helium_neutral_fraction, sk0qt, temp, r2t, r2p
      double precision :: sk1qt, r1t, r1p
      double precision :: sqet, sqep
      double precision :: usum
      double precision :: sx1qt, sx1qp, sx2qt, sx2qp
      double precision :: sx0qt, sx0qp
      double precision :: species_temp_deriv_term(11)
      double precision :: betaut, ramu, qut
      double precision :: pdt, pdtq, qcpi
      double precision :: beta16
      double precision :: utsum, upsum
      double precision :: etemp
      double precision :: stemp1, stemp2, stemp3, stemp4, stemp5, &
           stemp6, stemp7
      double precision :: sqett, sqetp
      double precision :: sx1qtt, sx1qtp, sx2qtt, sx2qtp
      double precision :: btemp
      double precision :: qqutt, qqutp

      data ionization_temp, helium_ionization_temp_1, &
           helium_ionization_temp_2/59630.,69450.,88710.,91330.,94570., &
           130630.,157800.,158000.,168630.,183080.,250200.,285270.,631370./
      data saha_weight_term/-.01,-.47,.60,.49,.08,.11,.0,-.04,.64,.99, &
           1.03,.6,.0/
      data nz/11/
      data max_saha_iterations/20/
      data saha_exponent_tol/4.0D1/
      data mean_electrons_per_ion, helium_ion_fraction_1, &
           helium_ion_fraction_2/0.1D0,0.0D0,0.0D0/
      data saha_convergence_tol/1.0D-09/
!  1    FXHE = FX(12)
      helium_mass_fraction = saha_mass_fractions(12)
      max_electrons_per_ion = 1.0d0 + helium_mass_fraction
      temperature_inverse = 1.0d0/temperature
! BEGIN SAHA ROUTINE
! COMPUTE SAHA K'S
      stemp = 2.50d0*log10_temperature - log10_pressure + cmkh
      nz1 = 1
      do i=1,nz
         ionization_temp_over_t(i) = ionization_temp(i)*temperature_inverse
         saha_ratio(i) = ln10*(saha_weight_term(i) + stemp) - &
              ionization_temp_over_t(i)
         if(saha_ratio(i).lt.-saha_exponent_tol) go to 12
         if(saha_ratio(i).gt.+saha_exponent_tol) go to 10
         saha_ratio(i) = beta_inverse*dexp(saha_ratio(i))
         go to 11
 10      saha_ratio(i) = 1.0d16
         nz1 = i + 1
 11   continue
      end do
      i = nz + 1
 12   nz0 = i - 1
      nz1 = min0(nz,nz1)
      skip_helium_i = .true.
      skip_helium_ii= .true.
      if(helium_mass_fraction.lt.1.0d-10) go to 13
      helium_saha_ratio_1 = ln10*(saha_weight_term(12) + stemp) - &
           helium_ionization_temp_1/temperature
      if(helium_saha_ratio_1.lt.-saha_exponent_tol) go to 13
      skip_helium_i = .false.
      helium_saha_ratio_1 = beta_inverse*dexp(helium_saha_ratio_1)
      helium_saha_ratio_2 = ln10*(saha_weight_term(13) + stemp) - &
           helium_ionization_temp_2/temperature
      if(helium_saha_ratio_2.lt.-saha_exponent_tol) go to 14
      skip_helium_ii= .false.
      helium_saha_ratio_2 = beta_inverse*dexp(helium_saha_ratio_2)
      go to 15
 13   helium_ion_fraction_1 = 0.0d0
 14   helium_ion_fraction_2 = 0.0d0
 15   continue
! BEGIN ITERATIONS FOR SOLUTION OF E
      do isaha=1,max_saha_iterations
         saha_state = saha_state + 1
         ep1 = 1.0d0/(1.0d0 + mean_electrons_per_ion)
         ep1e = ep1/mean_electrons_per_ion
         ep12 = ep1*ep1
         eep1 = mean_electrons_per_ion*ep1
         c33 = -mean_electrons_per_ion*(1.0d0 + mean_electrons_per_ion)
         r3 = mean_electrons_per_ion
         if(nz0.le.0) go to 22
         do i=1,nz0
            div = 1.0d0/(saha_ratio(i) + eep1)
            temp8 = saha_ratio(i)*div
            species_ion_fraction(i) = temp8
            temp8 = temp8*saha_mass_fractions(i)
            r3 = r3 - temp8
            temp8 = temp8*eep1*div
            c33 = c33 - temp8
            species_weighted_term(i) = temp8
 21      continue
         end do
 22      c33 = c33*ep1e
         r3 = r3 - helium_mass_fraction*(helium_ion_fraction_1 + &
              helium_ion_fraction_2+helium_ion_fraction_2)
         converged = dabs(r3).lt.saha_convergence_tol
         if(skip_helium_i) go to 25
         c22 = helium_saha_ratio_1 + eep1
         c23 = helium_ion_fraction_1*ep12
         c32 = helium_mass_fraction
         r2 = helium_saha_ratio_1*(1.0d0-helium_ion_fraction_1- &
              helium_ion_fraction_2) - helium_ion_fraction_1*eep1
         if(skip_helium_ii) go to 23
         c11 = eep1
         c12 = -helium_saha_ratio_2
         c13 = helium_ion_fraction_2*ep12
         c21 = helium_saha_ratio_1
         c31 = 2.0d0*helium_mass_fraction
         r1 = helium_ion_fraction_1*helium_saha_ratio_2 - &
              helium_ion_fraction_2*eep1
! REDUCTION OF SAHA MATRIX TO DIAGONAL FORM
         cr12 = - c21/c11
         cr13 = - c31/c11
         c22 = c22 + cr12*c12
         c32 = c32 + cr13*c12
         c23 = c23 + cr12*c13
         c33 = c33 + cr13*c13
         if (r1.eq.0.0d0) go to 210
         r1l=dlog10(dabs(r1))
         if (cr12.eq.0.0d0) go to 200
         cr12l=dlog10(dabs(cr12))
         fact=cr12l+r1l
         if (fact.lt.-38.0d0) go to 200
         r2 = r2 + cr12*r1
 200     continue
         if (cr13.eq.0.0d0) go to 210
         cr13l=dlog10(dabs(cr13))
         fact=cr13l+r1l
         if (fact.lt.-38.0d0) go to 210
         r3 = r3 + cr13*r1
 210     continue
! ENTRY FOR NO FULLY IONIZED HELIUM (SAHEX2 = 0.0)
 23      cr23 = - c32/c22
         c33 = c33 + cr23*c23
         r3 = r3 + cr23*r2
! ENTRY FOR NEUTRAL HELIUM
 25      delta_electrons_per_ion = r3/c33
         if(skip_helium_i) go to 26
         deltx1 = (r2 - c23*delta_electrons_per_ion)/c22
         helium_ion_fraction_1 = helium_ion_fraction_1 + deltx1
         helium_ion_fraction_1 = dmax1(0.0d0,dmin1(1.0d0,helium_ion_fraction_1))
         if(skip_helium_ii) go to 26
!CC   STATEMENT RECALCULATED FOR DEC-20 SYSTEM
         if (c12.eq.0.d0 .or. deltx1.eq.0.d0) go to 100
         c12l=dlog10(dabs(c12))
         delx1l=dlog10(dabs(deltx1))
         fact1=c12l+delx1l
         if (fact1.ge.-38.0d0) go to 100
         fact1=-38.0d0
         fact1=dexp(ln10*fact1)*dsign(1.0d0,c12)*dsign(1.0d0,deltx1)
         go to 105
 100     continue
         fact1=c12*deltx1
 105     continue
         if (c13.eq.0.d0 .or. delta_electrons_per_ion.eq.0.d0) go to 110
         c13l=dlog10(dabs(c13))
         deltel=dlog10(dabs(delta_electrons_per_ion))
         fact2=c13l+deltel
         if (fact2.ge.-38.0d0) go to 110
         fact2=-38.0d0
         fact2=dexp(ln10*fact2)*dsign(1.0d0,c13)*dsign(1.0d0,delta_electrons_per_ion)
         go to 115
 110     continue
         fact2=c13*delta_electrons_per_ion
 115     continue
         deltx2=(r1-fact1-fact2)/c11
         helium_ion_fraction_2 = helium_ion_fraction_2 + deltx2
         helium_ion_fraction_2 = dmax1(0.0d0,dmin1(1.0d0,helium_ion_fraction_2))
 26      mean_electrons_per_ion = mean_electrons_per_ion + delta_electrons_per_ion
         mean_electrons_per_ion = dmax1(1.0d-11,dmin1(max_electrons_per_ion, &
              mean_electrons_per_ion))
         converged = converged .and. dabs(delta_electrons_per_ion).lt.saha_convergence_tol
         if(converged) go to 29
 28   continue
      end do
      write(short_file_unit,99) log10_temperature,log10_pressure, &
           mean_electrons_per_ion,delta_electrons_per_ion
 99   format(' -----SAHA FAILURE (TL,PL)=',2F10.6,'  (E,DE)=',2F20.12)
! SYSTEM HAS BEEN SOLVED FOR E
 29   electron_mean_weight_inverse = mean_electrons_per_ion*ion_mean_weight_inverse
      specific_gas_constant = gas_constant*(electron_mean_weight_inverse + &
           ion_mean_weight_inverse)
      density = beta*pressure/(specific_gas_constant*temperature)
      log10_density = dlog10(density)
      ion_fraction(1) = 0.0d0
      if(nz0.ge.7) ion_fraction(1) = species_ion_fraction(7)
      ion_fraction(2) = helium_ion_fraction_1
      ion_fraction(3) = helium_ion_fraction_2
      if(in_atmosphere) go to 60
! COMPUTE FIRST DERIVATIVES
      rmub = specific_gas_constant*beta_inverse
      r3t = 0.0d0
      r3p = 0.0d0
      if(nz0.le.0) go to 31
      do i=nz1,nz0
         r3p = r3p - species_weighted_term(i)
!  30   R3T = R3T -FXS(I)*(2.5D0 + BETA14 + SAHATT(I))
         r3t = r3t -species_weighted_term(i)*(2.5d0 + beta14 + &
              ionization_temp_over_t(i))
 30   continue
      end do
      r3p = -beta_inverse*r3p
 31   if(skip_helium_i) go to 33
      helium_neutral_fraction = 1.0d0 - helium_ion_fraction_1 - helium_ion_fraction_2
      sk0qt = 2.5d0 + beta14 + helium_ionization_temp_1*temperature_inverse
      temp = helium_saha_ratio_1*helium_neutral_fraction
      r2t = temp*sk0qt
      r2p = -temp*beta_inverse
      if(skip_helium_ii) go to 32
      sk1qt = 2.5d0 + beta14 + helium_ionization_temp_2*temperature_inverse
      temp = helium_saha_ratio_2*helium_ion_fraction_1
      r1t = temp*sk1qt
      r1p = -temp*beta_inverse
      r2t = r2t + cr12*r1t
      r2p = r2p + cr12*r1p
      r3t = r3t + cr13*r1t
      r3p = r3p + cr13*r1p
 32   r3t = r3t + cr23*r2t
      r3p = r3p + cr23*r2p
 33   sqet = r3t/c33
      sqep = r3p/c33
      dlnrho_dlnt = -1.0d0 - beta14 - ep1*sqet
      dlnrho_dlnp = beta_inverse - ep1*sqep
! COMPUTE INTERNAL ENERGY TEMPERATURE DERIVATIVE(QUT)
      usum = 0.0d0
      if(skip_helium_i) go to 35
      sx1qt = (r2t - c23*sqet)/c22
      sx1qp = (r2p - c23*sqep)/c22
      sx2qt = 0.0d0
      sx2qp = 0.0d0
      if(skip_helium_ii) go to 34
      sx2qt = (r1t - c12*sx1qt - c13*sqet)/c11
      sx2qp = (r1p - c12*sx1qp - c13*sqep)/c11
 34   sx0qt = -(sx1qt + sx2qt)
      sx0qp = -(sx1qp + sx2qp)
      usum = helium_mass_fraction*(sx1qt*helium_ionization_temp_1 + &
           sx2qt*(helium_ionization_temp_1 + helium_ionization_temp_2))
 35   if(nz0.le.0) go to 37
      stemp = 2.50d0 + beta14 - ep1e*sqet
      do i=nz1,nz0
         species_temp_deriv_term(i) = species_weighted_term(i)*(stemp + &
              ionization_temp_over_t(i))
!  36   USUM = USUM + FXT(I)*SAHAT(I)
         usum = usum + species_temp_deriv_term(i)*ionization_temp(i)
 36   continue
      end do
 37   continue
      betaut = 0.75d0*(2.0d0 + beta14)
      ramu = gas_constant*ion_mean_weight_inverse
      qut =ramu*((1d0+mean_electrons_per_ion)*(betaut+3d0*beta_inverse*beta14) &
           +betaut*sqet+usum*temperature_inverse)
! COMPUTE SPECIFIC HEAT(QCP) AND ADIABATIC GRADIENT(DELA)
      pdt = pressure/(density*temperature)
      pdtq = pdt*dlnrho_dlnt
      specific_heat_cp = qut - pdtq
      qcpi = 1.0d0/specific_heat_cp
      adiabatic_gradient = -pdtq*qcpi
      if(.not.want_derivatives) go to 60
! COMPUTE DERIVATIVES OF QDT,QCP,DELA
      beta16 = 4.0d0*beta_inverse*beta14
      r3t = 0.0d0
      r3p = 0.0d0
      utsum = 0.0d0
      upsum = 0.0d0
      if(nz0.le.0) go to 42
      etemp = (1.0d0+mean_electrons_per_ion+mean_electrons_per_ion)/ &
           (mean_electrons_per_ion*(1.0d0+mean_electrons_per_ion))**2
      stemp1 = 2.5d0 + beta14 - sqet*ep1e
      stemp2 = -beta_inverse - sqep*ep1e
      stemp3 = etemp*sqep*sqet - beta_inverse*beta14
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
!  41   UPSUM = UPSUM + SAHAT(I)*STEMP7
         upsum = upsum + ionization_temp(i)*stemp7
 41   continue
      end do
 42   continue
      if(skip_helium_i) go to 44
      stemp = 2.0d0*helium_ion_fraction_1*ep1*sqet
      r2t = helium_saha_ratio_1*(helium_neutral_fraction*(beta16- &
           helium_ionization_temp_1*temperature_inverse+sk0qt**2)+2d0*sx0qt*sk0qt) &
           +ep12*(stemp*sqet - 2.0d0*sx1qt*sqet)
      r2p = helium_saha_ratio_1*(helium_neutral_fraction*beta_inverse* &
           (-beta14-sk0qt)+sx0qp*sk0qt-beta_inverse* &
           sx0qt) + ep12*(stemp*sqep - sx1qp*sqet - sx1qt*sqep)
      if(skip_helium_ii) go to 43
      stemp = 2.0d0*helium_ion_fraction_2*ep1*sqet
      r1t = helium_saha_ratio_2*(helium_ion_fraction_1*(beta16- &
           helium_ionization_temp_2*temperature_inverse+sk1qt**2)+2d0*sx1qt*sk1qt) &
           +ep12*(stemp*sqet - 2.0d0*sx2qt*sqet)
      r1p = helium_saha_ratio_2*(helium_ion_fraction_1*beta_inverse* &
           (-beta14-sk1qt)+sx1qp*sk1qt-beta_inverse* &
           sx1qt) + ep12*(stemp*sqep - sx2qp*sqet - sx2qt*sqep)
      r2t = r2t + cr12*r1t
      r2p = r2p + cr12*r1p
      r3p = r3p + cr13*r1p
      r3t = r3t + cr13*r1t
 43   r3p = r3p + cr23*r2p
      r3t = r3t + cr23*r2t
 44   sqett = r3t/c33
      sqetp = r3p/c33
      dlnrho_dlnt_dt = -beta16 + (ep1*sqet)**2 - ep1*sqett
      dlnrho_dlnp_dt = beta_inverse*beta14 + ep12*sqet*sqep - ep1*sqetp
      if(nz0.le.0) go to 46
      stemp = ep1e*sqett
      stemp1 = ep1e*sqetp
      do i=nz1,nz0
         utsum = utsum - species_weighted_term(i)*ionization_temp(i)*stemp
!  45   UPSUM = UPSUM - FXS(I)*SAHAT(I)*STEMP1
         upsum = upsum - species_weighted_term(i)*ionization_temp(i)*stemp1
 45   continue
      end do
 46   continue
      if(skip_helium_i) go to 48
      sx1qtt = (r2t - c23*sqett)/c22
      sx1qtp = (r2p - c23*sqetp)/c22
      if(skip_helium_ii) go to 47
      sx2qtt = (r1t - c12*sx1qtt - c13*sqett)/c11
      sx2qtp = (r1p - c12*sx1qtp - c13*sqetp)/c11
      utsum = utsum + helium_mass_fraction*(helium_ionization_temp_1+ &
           helium_ionization_temp_2)*(sx2qtt - sx2qt)
      upsum = upsum + helium_mass_fraction*(helium_ionization_temp_1+ &
           helium_ionization_temp_2)*sx2qtp
 47   utsum = utsum + helium_mass_fraction*helium_ionization_temp_1*(sx1qtt - sx1qt)
      upsum = upsum + helium_mass_fraction*helium_ionization_temp_1*sx1qtp
 48   continue
      btemp = (1.0d0+mean_electrons_per_ion)*beta14*beta_inverse* &
           (24.0d0*beta_inverse - 9.0d0)
      qqutt = ramu*(btemp+0.75d0*(2.0d0+beta14*(8.0d0*beta_inverse+1d0))*sqet+ &
           betaut*sqett + utsum*temperature_inverse)
      qqutp = ramu*(-0.25d0*btemp+0.75d0*((2.0d0+beta14*(4.0d0*beta_inverse &
           +1.0d0))*sqep-beta14*beta_inverse*sqet) +betaut*sqetp+upsum*temperature_inverse)
      specific_heat_cp_dt = qcpi*(qqutt + pdtq*(1.0d0+dlnrho_dlnt) - pdt*dlnrho_dlnt_dt)
      specific_heat_cp_dp = qcpi*(qqutp - pdtq*(1.0d0-dlnrho_dlnp) - pdt*dlnrho_dlnp_dt)
      adiabatic_gradient_dt = -1.0d0 - dlnrho_dlnt - specific_heat_cp_dt + &
           dlnrho_dlnt_dt/dlnrho_dlnt
      adiabatic_gradient_dp = 1.0d0 - dlnrho_dlnp - specific_heat_cp_dp + &
           dlnrho_dlnp_dt/dlnrho_dlnt
 60   continue

      return
end subroutine eqsaha

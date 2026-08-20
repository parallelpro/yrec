!----------------------------------------------------------------------
! condopacpint
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original condopacpint.f; only variable names, source form, and
! comment style were updated.
!
!     condOpacPInt -  Calculate  Conductive Opacities based on Chirs
!                     Burke's conductive opacity code and callin
!                     on Potehkin code and table.
!
!                      LLP    5/31/04
!
! Combines the fully-ionized-species Potekhin conductive opacity
! (condopacp.f90, called for pure H, He, and a Z=8 metal proxy) into a
! single mixture value, weighted by each species' actual ionization
! fraction (ion_fraction, from the eqsaha/eqrelv solve upstream).
! Called from getopac.f90 whenever conductive opacity is enabled.
subroutine condopacpint(log10_density, log10_temperature, &
     hydrogen_fraction, metal_fraction, conductive_opacity, &
     conductive_log10_opacity, conductive_dlnkap_dlnrho, &
     conductive_dlnkap_dlnt, ion_fraction, got_conductive_opacity)

      implicit none

      double precision, intent(in) :: log10_density, log10_temperature, &
           hydrogen_fraction, metal_fraction
      double precision, intent(out) :: conductive_opacity, &
           conductive_log10_opacity, conductive_dlnkap_dlnrho, &
           conductive_dlnkap_dlnt
      double precision, intent(in) :: ion_fraction(3)
      logical, intent(out) :: got_conductive_opacity

      save
      double precision :: atomic_weight_h1, atomic_weight_he4, &
           atomic_weight_ox
      data atomic_weight_h1,atomic_weight_he4,atomic_weight_ox &
           /1.008D0,4.004D0,16D0/

! --- locals ---
      double precision :: helium_fraction
      double precision :: number_fraction_h1, number_fraction_he4, &
           number_fraction_ox
      double precision :: weight_h1, weight_he4, weight_ox
      double precision :: extrapolation_indicator
      double precision :: log10_cond_h1, dlnkap_dlnrho_h1, dlnkap_dlnt_h1
      double precision :: log10_cond_he4, dlnkap_dlnrho_he4, dlnkap_dlnt_he4
      double precision :: log10_cond_ox, dlnkap_dlnrho_ox, dlnkap_dlnt_ox
      double precision :: cond_h1, cond_he4, cond_ox
      double precision :: mix_log10_cond
      double precision :: unused_deriv1, unused_deriv2

! Get fractional abundances in nubers of atoms of H1, He4 and metals (Ox)
      helium_fraction=1.0D0-hydrogen_fraction-metal_fraction
      number_fraction_h1  = (hydrogen_fraction/atomic_weight_h1)  / &
           (hydrogen_fraction/atomic_weight_h1 + &
           helium_fraction/atomic_weight_he4 + metal_fraction/atomic_weight_ox)
      number_fraction_he4 = (helium_fraction/atomic_weight_he4) / &
           (hydrogen_fraction/atomic_weight_h1 + &
           helium_fraction/atomic_weight_he4 + metal_fraction/atomic_weight_ox)
      number_fraction_ox  = (metal_fraction/atomic_weight_ox)  / &
           (hydrogen_fraction/atomic_weight_h1 + &
           helium_fraction/atomic_weight_he4 + metal_fraction/atomic_weight_ox)

! In addition we now weight the H1, He and Ox by oncly
! considering those atoms H1 and He4 atomos which are
! charged. FXION(1) is the fraction of H1 atoms that are
! (singly) charged. FXION(2) is the fraction of He4 atoms
! that are singly charged and FXION(3) the fractoin that
! are doubly charged.  It is presumed that the metals (Ox)
! are fully ionized in first approximation. Because of
! their low abundance (tyically less than 2% divided by
! the atomic mass
       weight_h1 = number_fraction_h1 * ion_fraction(1) ! (FXION(1)/1D0 is the fraction
       weight_he4 = number_fraction_he4 * (ion_fraction(2) + 2D0*ion_fraction(3))/2D0
       weight_ox = number_fraction_ox

! We initially presume no significant conductive opacity
       got_conductive_opacity=.FALSE.
!       OC = 1D50
!       OCL = LOG10(OC)
       conductive_log10_opacity = 9.9998D0
       conductive_opacity = 10D0**conductive_log10_opacity
       conductive_dlnkap_dlnrho = conductive_opacity
       conductive_dlnkap_dlnt = conductive_opacity

       extrapolation_indicator = log10_density - 3.0D0*(log10_temperature-6.0D0)
!     DO CONDUCTIVE OPACITY CORRECTION in fully ionized approximation
      if(log10_temperature.ge.3.0D0.and.log10_temperature.le.9.0D0 &
           .and.log10_density.le.9.75D0) then
         if (log10_density.ge.-6.0D0) then
           call condopacp(1.0D0,log10_temperature,log10_density, &
                log10_cond_h1,dlnkap_dlnrho_h1,dlnkap_dlnt_h1)
           call condopacp(2.0D0,log10_temperature,log10_density, &
                log10_cond_he4,dlnkap_dlnrho_he4,dlnkap_dlnt_he4)
           call condopacp(8.0D0,log10_temperature,log10_density, &
                log10_cond_ox,dlnkap_dlnrho_ox,dlnkap_dlnt_ox)

           cond_h1 = 10.0D0**(-log10_cond_h1)
           cond_he4 = 10.0D0**(-log10_cond_he4)
           cond_ox = 10.0D0**(-log10_cond_ox)

           mix_log10_cond = -log10(weight_h1*cond_h1 + weight_he4*cond_he4 + &
                weight_ox*cond_ox)
!           CONDL = LOG10( (X + 1D0) / 2D0 ) -
!     *             LOG10( X*CONDX + .5D0*Y*CONDY + .5D0*z*CONDZ)

           conductive_dlnkap_dlnrho = mix_log10_cond * (weight_h1*dlnkap_dlnrho_h1 &
                + weight_he4*dlnkap_dlnrho_he4 + weight_ox*dlnkap_dlnrho_ox)
!           QODC = 2D0 * CONDL / ( 1D0 + X ) *
!     *          (X*CONDX*QODCX +.5D0*Y*CONDY*QODCY +.5D0*Z*COBDZ*QODCZ)

            conductive_dlnkap_dlnt = mix_log10_cond * (weight_h1*dlnkap_dlnrho_h1 &
                 + weight_he4*dlnkap_dlnrho_he4 + weight_ox*dlnkap_dlnrho_ox)
!           QOTC = 2D0 * CONDL / ( 1D0 + X ) *
!     *          (X*CONDX*QOTCX +.5D0*Y*CONDY*QOTCY +.5D0*Z*CONDZ*QOTCZ)

         conductive_log10_opacity = -3.5194D0+3.0D0*log10_temperature-log10_density-mix_log10_cond
           conductive_dlnkap_dlnrho = -1.0D0-conductive_dlnkap_dlnrho
           conductive_dlnkap_dlnt = 3.0D0-conductive_dlnkap_dlnt
           got_conductive_opacity=.TRUE.
           conductive_opacity = 10.0D0**conductive_log10_opacity
         end if
         if (log10_density.lt.-6.0D0.and.extrapolation_indicator.ge.0.0D0) then
! Extrapolate Conductive opacity
           call condopacp(1.0D0,log10_temperature,-6.0D0,log10_cond_h1, &
                conductive_dlnkap_dlnrho,conductive_dlnkap_dlnt)
           call condopacp(2.0D0,log10_temperature,-6.0D0,log10_cond_he4, &
                unused_deriv1,unused_deriv2)
           call condopacp(8.0D0,log10_temperature,-6.0D0,log10_cond_ox, &
                unused_deriv1,unused_deriv2)

           cond_h1 = 10.0D0**(-log10_cond_h1)
           cond_he4 = 10.0D0**(-log10_cond_he4)
           cond_ox = 10.0D0**(-log10_cond_ox)

           mix_log10_cond = log10(weight_h1*cond_h1 + weight_he4*cond_he4 + &
                weight_ox*cond_ox)
!           CONDL = LOG10( (X + 1D0) / 2D0 ) -
!     *            LOG10( X*CONDX + .5D0*Y*CONDY + .5D0*Z*CONDZ)

           conductive_log10_opacity = -3.5194D0+3.0D0*log10_temperature+6.0D0-mix_log10_cond
           conductive_dlnkap_dlnrho = -1.0D0-conductive_dlnkap_dlnrho
           conductive_dlnkap_dlnt = 3.0D0-conductive_dlnkap_dlnt
           conductive_log10_opacity = conductive_log10_opacity + &
                conductive_dlnkap_dlnrho*(-6.0D0-log10_density)
           conductive_opacity = 10.0D0**conductive_log10_opacity
           got_conductive_opacity=.TRUE.
         end if
      end if

      return
end subroutine condopacpint

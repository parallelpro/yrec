!----------------------------------------------------------------------
! eqstat
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original eqstat.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! New front end to equation of state routines. Provided to allow
! numerical differentiation of the EOS (via eqstat2) by calling it at
! adjacent P and T points, when use_numerical_derivatives is set,
! instead of relying on eqstat2's own (analytic/table) derivatives.
!                                                                 LLP  10-22-06
subroutine eqstat(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
     in_atmosphere, saha_state)
!
!  Input Arguments: log10_temperature, log10_pressure, hydrogen_fraction,
!          metal_fraction, want_derivatives, in_atmosphere
!
!  Output Arguments: temperature, pressure, density, log10_density, beta,
!          beta_inverse, beta14, ion_fraction, specific_gas_constant,
!          ion_mean_weight_inverse, electron_mean_weight_inverse,
!          dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, adiabatic_gradient,
!          dlnrho_dlnt_dt, dlnrho_dlnp_dt, adiabatic_gradient_dt,
!          adiabatic_gradient_dp, specific_heat_cp_dt, specific_heat_cp_dp
!
!  Update (Input and Output) Arguments: saha_state
!

      use luout_lib
      use const_lib
      implicit none

      double precision, intent(inout) :: log10_temperature
      double precision, intent(out) :: temperature
      double precision, intent(inout) :: log10_pressure
      double precision, intent(out) :: pressure
      double precision, intent(inout) :: log10_density
      double precision, intent(out) :: density
      double precision, intent(in) :: hydrogen_fraction, metal_fraction
      double precision, intent(inout) :: beta
      double precision, intent(out) :: beta_inverse, beta14
      double precision, intent(inout) :: ion_fraction(3)
      double precision, intent(out) :: specific_gas_constant
      double precision, intent(inout) :: ion_mean_weight_inverse
      double precision, intent(out) :: electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      double precision, intent(inout) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient
      double precision, intent(out) :: dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state

      integer, parameter :: nts = 63, nps = 76


! common/comp/: not used in this file; declared only to preserve
! layout. Naming matches getopac.f90/meqos.f90/eqstat2.f90.
      double precision :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, zenvm, envelope_amu, &
           envelope_species_fractions(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, envelope_amu, envelope_species_fractions, xnew, znew, &
           stotal, senv




! common/debhu/: Debye-Huckel correction data; not used in this file.
! DBG 7/92 common block added to compute Debye-Huckel correction.
      double precision :: cdh, etadh0, etadh1, zdh(18), xxdy, yydh, zzdh, &
           dhnue(18)
      logical :: ldh
      common/debhu/ cdh, etadh0, etadh1, zdh, xxdy, yydh, zzdh, dhnue, ldh

! common/opaleos/: only use_numerical_derivatives is used here (the
! rest select which OPAL EOS table eqstat2 blends in; not read here).
! YCK 2/95, LLP 2001, LLP 2006 OPAL eos; LLP added the
! use_numerical_derivatives flag 7/07. It is part of the PARMIN
! PHYSICS namelist.
      logical :: use_opal95_eos, use_opal2001_eos, use_opal2006_eos, &
           use_numerical_derivatives
      integer :: iopale
      common/opaleos/ use_opal95_eos, iopale, use_opal2001_eos, &
           use_opal2006_eos, use_numerical_derivatives

! common/gravs3/: not used in this file; declared only to preserve
! layout. MHP 3/94 added metal diffusion. Naming matches setupopac.f90.
      double precision :: fgry, fgrz
      logical :: lthoul, use_diffusion_z
      common/gravs3/ fgry, fgrz, lthoul, use_diffusion_z

! common/scveos/: not used in this file; declared only to preserve
! layout. MHP 5/97 added common block for SCV EOS tables. Naming
! matches eqstat2.f90.
      double precision :: tlogx(nts), tablex(nts,nps,12), &
           tabley(nts,nps,12), smix(nts,nps), tablez(nts,nps,13), &
           tablenv(nts,nps,12)
      integer :: nptsx(nts), idt, idp
      logical :: use_scv_eos
      common/scveos/ tlogx, tablex, tabley, smix, tablez, tablenv, nptsx, &
           use_scv_eos, idt, idp

      save

!  want_derivatives: if true, provide derivatives needed for
!  relaxation, else don't
!
!  use_numerical_derivatives - if true derivatives are calculated
!  numerically, else get them the old way using Yale EOS.

      double precision :: log10_pressure_orig, log10_temperature_orig
      logical :: in_atmosphere_local
      integer :: saha_state_local
      logical :: want_derivatives_2
      double precision :: dpl, dtl, dpl2, dtl2
      double precision :: log10_density_1, density_1
      double precision :: ttl, ppl
      double precision :: dlnrho_dlnt_1, dlnrho_dlnp_1, specific_heat_cp_1, &
           adiabatic_gradient_1
      double precision :: dlnrho_dlnt_2, dlnrho_dlnp_2, specific_heat_cp_2, &
           adiabatic_gradient_2
      double precision :: dlnrho_dlnt_dt_1, dlnrho_dlnp_dt_1, &
           adiabatic_gradient_dt_1, adiabatic_gradient_dp_1, &
           specific_heat_cp_dt_1, specific_heat_cp_dp_1
      double precision :: dlnrho_dlnt_dt_throwaway, dlnrho_dlnp_dt_throwaway, &
           adiabatic_gradient_dt_throwaway, adiabatic_gradient_dp_throwaway, &
           specific_heat_cp_dt_throwaway, specific_heat_cp_dp_throwaway

      log10_pressure_orig = log10_pressure
      log10_temperature_orig = log10_temperature
      in_atmosphere_local = in_atmosphere
      saha_state_local = saha_state

      if (want_derivatives .and. use_numerical_derivatives) then
!        Get Numerical Derivatives of Current EOS    LLP  8/5/07
!        If both derivatives and numerical derivatives are requested.
!
!        A central difference approximation using values on both sides
!        of TL and then PL is used. The error term in this
!        approximation of the derivatives is of order h**2, where h is
!        the interval. The error term when one sided derivatives are
!        used is of order h.

         dpl = .150d0    ! the approximate table intervals
         dtl = .030d0

         dpl = dpl * .01d0   ! Empirically, scaling the intervals down
         dtl = dtl * .01d0   ! by 100 seemed to work best.

         dpl2 = 2d0 * dpl
         dtl2 = 2d0 * dtl

         log10_density_1 = log10_density   ! Initialize to called values.
         density_1 = density     ! This can avoid problems in initial table lookups

         ttl = log10_temperature + dtl   ! Get derivatives wrt T.
         temperature = 10.0d0**ttl
         want_derivatives_2 = .false.
         call eqstat2(ttl, temperature, log10_pressure, pressure, &
              log10_density_1, density_1, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt_1, dlnrho_dlnp_1, specific_heat_cp_1, &
              adiabatic_gradient_1, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere_local, saha_state_local)

         ttl = log10_temperature - dtl
         temperature = 10.0d0**ttl
         want_derivatives_2 = .false.
         call eqstat2(ttl, temperature, log10_pressure, pressure, &
              log10_density_1, density_1, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt_2, dlnrho_dlnp_2, specific_heat_cp_2, &
              adiabatic_gradient_2, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere_local, saha_state_local)
         dlnrho_dlnt_dt_1 = (dlnrho_dlnt_1 - dlnrho_dlnt_2)/dtl2/ln10
         specific_heat_cp_dt_1 = (dlog10(specific_heat_cp_1) - &
              dlog10(specific_heat_cp_2))/dtl2
         adiabatic_gradient_dt_1 = (dlog10(adiabatic_gradient_1) - &
              dlog10(adiabatic_gradient_2))/dtl2
         log10_temperature = log10_temperature_orig
         temperature = 10.0d0**log10_temperature_orig   ! Restore original T

         ppl = log10_pressure + dpl
         pressure = 10.0d0**ppl
         want_derivatives_2 = .false.
         call eqstat2(log10_temperature, temperature, ppl, pressure, &
              log10_density_1, density_1, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt_1, dlnrho_dlnp_1, specific_heat_cp_1, &
              adiabatic_gradient_1, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere_local, saha_state_local)
         ppl = log10_pressure - dpl
         pressure = 10.0d0**ppl
         want_derivatives_2 = .false.
         call eqstat2(log10_temperature, temperature, ppl, pressure, &
              log10_density_1, density_1, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt_2, dlnrho_dlnp_2, specific_heat_cp_2, &
              adiabatic_gradient_2, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere_local, saha_state_local)
         log10_pressure = log10_pressure_orig
         pressure = 10.0d0**log10_pressure_orig   ! Restore original P
         dlnrho_dlnp_dt_1 = (dlnrho_dlnt_1 - dlnrho_dlnt_2)/dpl2/ln10
         specific_heat_cp_dp_1 = (dlog10(specific_heat_cp_1) - &
              dlog10(specific_heat_cp_2))/dpl2
         adiabatic_gradient_dp_1 = (dlog10(adiabatic_gradient_1) - &
              dlog10(adiabatic_gradient_2))/dpl2

         dlnrho_dlnt_dt = dlnrho_dlnt_dt_1
         dlnrho_dlnp_dt = dlnrho_dlnp_dt_1
         adiabatic_gradient_dt = adiabatic_gradient_dt_1
         adiabatic_gradient_dp = adiabatic_gradient_dp_1
         specific_heat_cp_dp = specific_heat_cp_dp_1
         specific_heat_cp_dt = specific_heat_cp_dt_1
      end if

      log10_temperature = log10_temperature_orig   ! Restore original TL and PL.
      log10_pressure = log10_pressure_orig

      if (want_derivatives .and. .not.use_numerical_derivatives) then
         want_derivatives_2 = .true.   ! Need derivatives and have no numerical ones.
                          ! Call eqstat2 and request derivatives
         call eqstat2(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives_2, &
              in_atmosphere, saha_state)
      else
         want_derivatives_2 = .false.  ! We either already have numerical derivatives
                        ! or do not need any derivatives.
                          ! Call eqstat2 and request no derivatives
         call eqstat2(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere, saha_state)
!            Note that the dlnrho_dlnt_dt, dlnrho_dlnp_dt, etc OUTPUTS
!            are to dummy variables so they can not affect the
!            previously calculated second derivatives.
      end if

      return
end subroutine eqstat

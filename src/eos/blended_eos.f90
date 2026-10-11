!----------------------------------------------------------------------
! blended_eos
!----------------------------------------------------------------------
! The equation of state of every run that does not use MHD: the
! Yale EOS, with the SCV table at low temperature and an OPAL table
! blended in where it covers the point.
!     blend_eos            one evaluation at (log T, log P, X, Z):
!                          Saha (or SCV) below the ionization cutoff,
!                          fully ionized above it, interpolated across
!                          it, then the OPAL result adopted inside the
!                          table or blended in near its edge
!     compute_blended_eos  blend_eos, with the T and P derivatives of
!                          the gradients taken numerically (four more
!                          evaluations) when use_numerical_derivatives
!                          is set
module blended_eos
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use yale_eos, only: compute_saha_eos, compute_fully_ionized_eos
      use scv_eos, only: compute_scv_envelope_eos
      use opal95_eos, only: compute_opal95_eos, find_opal95_ramp
      use opal2001_eos, only: compute_opal2001_eos, find_opal2001_ramp
      use opal2006_eos, only: compute_opal2006_eos, find_opal2006_ramp
      use physical_constants, only: gas_constant, ln10, &
           radiation_constant_over_3
      use io_units, only: run_log_unit
      use eos_mixture, only: eos_mix, ion_mean_weight_excess, ix_c, ix_h, &
           ix_he, ix_na, ix_ne, ix_o, n_mix_species
      use eos_configuration, only: eos_settings
      use math_lib

      implicit none
      private
      public :: compute_blended_eos, blend_eos

contains

!----------------------------------------------------------------------
! compute_blended_eos (formerly eqstat)
!----------------------------------------------------------------------
! New front end to equation of state routines. Provided to allow
! numerical differentiation of the EOS (via blend_eos, below) by calling
! it at adjacent P and T points, when use_numerical_derivatives is
! set, instead of relying on blend_eos's own (analytic/table)
! derivatives.
!                                                                 LLP  10-22-06
subroutine compute_blended_eos(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
     in_atmosphere, saha_state, ierr)
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


      real(dp), intent(inout) :: log10_temperature
      real(dp), intent(out) :: temperature
      real(dp), intent(inout) :: log10_pressure
      real(dp), intent(out) :: pressure
      real(dp), intent(inout) :: log10_density
      real(dp), intent(out) :: density
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      real(dp), intent(inout) :: beta
      real(dp), intent(out) :: beta_inverse, beta14
      real(dp), intent(inout) :: ion_fraction(:)
      real(dp), intent(out) :: specific_gas_constant
      real(dp), intent(inout) :: ion_mean_weight_inverse
      real(dp), intent(out) :: electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      real(dp), intent(inout) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient
      real(dp), intent(out) :: dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state

!  want_derivatives: if true, provide derivatives needed for
!  relaxation, else don't
!
!  use_numerical_derivatives - if true derivatives are calculated
!  numerically, else get them the old way using Yale EOS.

      real(dp) :: log10_pressure_orig, log10_temperature_orig
      logical :: in_atmosphere_local
      integer :: saha_state_local
      logical :: want_derivatives_2
      real(dp) :: dpl, dtl, dpl2, dtl2
      real(dp) :: log10_density_1, density_1
      real(dp) :: ttl, ppl
      real(dp) :: dlnrho_dlnt_1, dlnrho_dlnp_1, specific_heat_cp_1, &
           adiabatic_gradient_1
      real(dp) :: dlnrho_dlnt_2, dlnrho_dlnp_2, specific_heat_cp_2, &
           adiabatic_gradient_2
      real(dp) :: dlnrho_dlnt_dt_1, dlnrho_dlnp_dt_1, &
           adiabatic_gradient_dt_1, adiabatic_gradient_dp_1, &
           specific_heat_cp_dt_1, specific_heat_cp_dp_1
      real(dp) :: dlnrho_dlnt_dt_throwaway, dlnrho_dlnp_dt_throwaway, &
           adiabatic_gradient_dt_throwaway, adiabatic_gradient_dp_throwaway, &
           specific_heat_cp_dt_throwaway, specific_heat_cp_dp_throwaway

      integer, intent(out) :: ierr

      ierr = 0

      log10_pressure_orig = log10_pressure
      log10_temperature_orig = log10_temperature
      in_atmosphere_local = in_atmosphere
      saha_state_local = saha_state

      if (want_derivatives .and. eos_settings%use_numerical_derivatives) then
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
         temperature = exp10(ttl)
         want_derivatives_2 = .false.
         call blend_eos(ttl, temperature, log10_pressure, pressure, &
              log10_density_1, density_1, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt_1, dlnrho_dlnp_1, specific_heat_cp_1, &
              adiabatic_gradient_1, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere_local, saha_state_local, ierr)
         if (ierr /= 0) return

         ttl = log10_temperature - dtl
         temperature = exp10(ttl)
         want_derivatives_2 = .false.
         call blend_eos(ttl, temperature, log10_pressure, pressure, &
              log10_density_1, density_1, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt_2, dlnrho_dlnp_2, specific_heat_cp_2, &
              adiabatic_gradient_2, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere_local, saha_state_local, ierr)
         if (ierr /= 0) return
         dlnrho_dlnt_dt_1 = (dlnrho_dlnt_1 - dlnrho_dlnt_2)/dtl2/ln10
         specific_heat_cp_dt_1 = (log10(specific_heat_cp_1) - &
              log10(specific_heat_cp_2))/dtl2
         adiabatic_gradient_dt_1 = (log10(adiabatic_gradient_1) - &
              log10(adiabatic_gradient_2))/dtl2
         log10_temperature = log10_temperature_orig
         temperature = exp10(log10_temperature_orig)   ! Restore original T

         ppl = log10_pressure + dpl
         pressure = exp10(ppl)
         want_derivatives_2 = .false.
         call blend_eos(log10_temperature, temperature, ppl, pressure, &
              log10_density_1, density_1, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt_1, dlnrho_dlnp_1, specific_heat_cp_1, &
              adiabatic_gradient_1, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere_local, saha_state_local, ierr)
         if (ierr /= 0) return
         ppl = log10_pressure - dpl
         pressure = exp10(ppl)
         want_derivatives_2 = .false.
         call blend_eos(log10_temperature, temperature, ppl, pressure, &
              log10_density_1, density_1, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt_2, dlnrho_dlnp_2, specific_heat_cp_2, &
              adiabatic_gradient_2, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere_local, saha_state_local, ierr)
         if (ierr /= 0) return
         log10_pressure = log10_pressure_orig
         pressure = exp10(log10_pressure_orig)   ! Restore original P
         dlnrho_dlnp_dt_1 = (dlnrho_dlnt_1 - dlnrho_dlnt_2)/dpl2/ln10
         specific_heat_cp_dp_1 = (log10(specific_heat_cp_1) - &
              log10(specific_heat_cp_2))/dpl2
         adiabatic_gradient_dp_1 = (log10(adiabatic_gradient_1) - &
              log10(adiabatic_gradient_2))/dpl2

         dlnrho_dlnt_dt = dlnrho_dlnt_dt_1
         dlnrho_dlnp_dt = dlnrho_dlnp_dt_1
         adiabatic_gradient_dt = adiabatic_gradient_dt_1
         adiabatic_gradient_dp = adiabatic_gradient_dp_1
         specific_heat_cp_dp = specific_heat_cp_dp_1
         specific_heat_cp_dt = specific_heat_cp_dt_1
      end if

      log10_temperature = log10_temperature_orig   ! Restore original TL and PL.
      log10_pressure = log10_pressure_orig

      if (want_derivatives .and. .not.eos_settings%use_numerical_derivatives) then
         want_derivatives_2 = .true.   ! Need derivatives and have no numerical ones.
                          ! Call blend_eos and request derivatives
         call blend_eos(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives_2, &
              in_atmosphere, saha_state, ierr)
         if (ierr /= 0) return
      else
         want_derivatives_2 = .false.  ! We either already have numerical derivatives
                        ! or do not need any derivatives.
                          ! Call blend_eos and request no derivatives
         call blend_eos(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt_throwaway, &
              dlnrho_dlnp_dt_throwaway, adiabatic_gradient_dt_throwaway, &
              adiabatic_gradient_dp_throwaway, specific_heat_cp_dt_throwaway, &
              specific_heat_cp_dp_throwaway, want_derivatives_2, &
              in_atmosphere, saha_state, ierr)
         if (ierr /= 0) return
!            Note that the dlnrho_dlnt_dt, dlnrho_dlnp_dt, etc OUTPUTS
!            are to dummy variables so they can not affect the
!            previously calculated second derivatives.
      end if

end subroutine compute_blended_eos

!----------------------------------------------------------------------
! blend_eos (formerly eqstat2)
!----------------------------------------------------------------------
! One evaluation of the EOS at (log T, log P, X, Z). Below the
! ionization cutoff the gas is solved with the Saha equations
! (compute_saha_eos), or looked up in the SCV table
! (compute_scv_envelope_eos) with Saha supplying what SCV lacks; above
! it, it is fully ionized (compute_fully_ionized_eos); near the cutoff
! the two are interpolated. Then, if an OPAL table (1995, 2001 or
! 2006) covers the point, its result is adopted, or blended in with
! the weight find_opal*_ramp returns near the table's edge.
! metal_fraction is intent(in): no callee writes it.
subroutine blend_eos(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
     in_atmosphere, saha_state, ierr)

      real(dp), intent(inout) :: log10_temperature
      real(dp), intent(out) :: temperature
      real(dp), intent(inout) :: log10_pressure
      real(dp), intent(out) :: pressure
      real(dp), intent(inout) :: log10_density
      real(dp), intent(out) :: density
      real(dp), intent(in) :: hydrogen_fraction
      real(dp), intent(in) :: metal_fraction
      real(dp), intent(inout) :: beta
      real(dp), intent(out) :: beta_inverse, beta14
      real(dp), intent(inout) :: ion_fraction(:)
      real(dp), intent(out) :: specific_gas_constant
      real(dp), intent(inout) :: ion_mean_weight_inverse
      real(dp), intent(out) :: electron_mean_weight_inverse, &
           electron_degeneracy_parameter
      real(dp), intent(inout) :: dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient
      real(dp), intent(out) :: dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
           adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state
! --- locals ---
      integer :: num_species, species_idx
      real(dp) :: inverse_atomic_weights(4), atomic_weights_full(n_mix_species)
      real(dp) :: saha_ramp_width, saha_ramp_scale
      data num_species/12/
      data inverse_atomic_weights/0.9921d0, 0.24975d0, 0.08322d0, 0.4995d0/
      data saha_ramp_width, saha_ramp_scale/0.500d0, 2.000d0/
! atomic weights in the ix_* species order (Na, Al, Mg, Fe, Si, C, H,
! O, N, Ar, Ne, He)
      data atomic_weights_full/23.0d0, 26.99d0, 24.32d0, 55.86d0, 28.1d0, &
           12.015d0, 1.008d0, 16.0d0, 14.01d0, 39.96d0, 20.19d0, 4.004d0/

      logical :: need_saha_solution, skip_relativistic_eos
      real(dp) :: saha_mass_fractions(n_mix_species)
      real(dp) :: metal_ratio, amu_correction, h_excess, y_excess
      real(dp) :: envelope_amu_over_amu
      real(dp) :: dfx1, dfx12, dfx4, amu_inverse, envelope_amu_frac
      logical :: use_envelope_mixture
      logical :: valid_table_point
      real(dp) :: opal_log10_density, opal_density, opal_beta, &
           opal_beta_inverse, opal_beta14, opal_specific_gas_constant, &
           opal_ion_mean_weight_inverse, opal_electron_mean_weight_inverse, &
           opal_dlnrho_dlnt, opal_dlnrho_dlnp, opal_specific_heat_cp, &
           opal_adiabatic_gradient
      real(dp) :: ramp_factor
      logical :: in_opal_table, needs_ramp
      logical :: opal_covers_point
! (T, rho) window inside which each OPAL vintage is consulted at all
! (the table edge itself is then handled by find_opal95_ramp/find_opal2001_ramp/
! find_opal2006_ramp). The three gates are not the same test: 1995 bounds
! log10 T by 8.0 where 2001/2006 bound T by 100e6 / 200e6 K, and 1995
! caps log10 rho at 5.0 against 7.0 for the other two -- the values
! are the original literals, kept as written.
      real(dp), parameter :: opal95_t_min = 5.0d3, &
           opal95_log10_t_max = 8.0d0, opal95_log10_rho_max = 5.0d0
      real(dp), parameter :: opal01_t_min = 2.0d3, &
           opal01_t_max = 100d6, opal01_log10_rho_max = 7.0d0
      real(dp), parameter :: opal06_t_min = 1.870d3, &
           opal06_t_max = 200d6, opal06_log10_rho_max = 7.0d0

! values saved across the compute_saha_eos/compute_fully_ionized_eos interpolation
! near the ionization cutoff
      real(dp) :: saha_log10_density, saha_dlnrho_dlnt, &
           saha_dlnrho_dlnp, saha_specific_heat_cp, saha_adiabatic_gradient
      real(dp) :: saha_dlnrho_dlnt_dt, saha_dlnrho_dlnp_dt, &
           saha_specific_heat_cp_dt, saha_specific_heat_cp_dp, &
           saha_adiabatic_gradient_dt, saha_adiabatic_gradient_dp
      real(dp) :: ramp_weight, ramp_weight_sq

!     need_saha_solution = T if T low enough to require a call to
!     compute_saha_eos to solve the Saha ionization equation; if F, a fully
!     ionized gas is assumed. skip_relativistic_eos = T if the Saha
!     equation must be solved but the fully-ionized (relativistic)
!     equation of state is not needed.
      integer, intent(out) :: ierr

      ierr = 0

      need_saha_solution = log10_temperature < eos_settings%saha_log10t_cutoff
      skip_relativistic_eos = log10_temperature <= &
           (eos_settings%saha_log10t_cutoff - saha_ramp_width)

!     MHP 3/94 metal diffusion added.  Assume all metals scale equally.
      if (eos_settings%use_diffusion_z) then
!        THE VECTOR eos_mix%fxenv IS DEFINED IN STARIN AS
!        (mass fraction of species/atomic wt/amuenv) for the composition
!        XENV,ZENV.
!        And amu is defined as the sum of the mass fractions divided by
!        their atomic weights.  For metal diffusion, all the metals are
!        assumed to change equally.
         metal_ratio = metal_fraction/eos_mix%envelope_metal_fraction
         amu_correction = (metal_ratio - 1.0d0)*eos_mix%amuenv
         ion_mean_weight_inverse = eos_mix%amuenv
         do species_idx = ix_na, ix_c
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 amu_correction*eos_mix%fxenv(species_idx)
         end do
         h_excess = (hydrogen_fraction - eos_mix%envelope_hydrogen_fraction)/ &
              atomic_weights_full(ix_h)
         ion_mean_weight_inverse = ion_mean_weight_inverse + h_excess
         do species_idx = ix_o, ix_ne
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 amu_correction*eos_mix%fxenv(species_idx)
         end do
         y_excess = (eos_mix%envelope_hydrogen_fraction + eos_mix%envelope_metal_fraction - &
              hydrogen_fraction - metal_fraction)/atomic_weights_full(ix_he)
         ion_mean_weight_inverse = ion_mean_weight_inverse + y_excess
         envelope_amu_over_amu = metal_ratio*eos_mix%amuenv/ &
              ion_mean_weight_inverse
         if (need_saha_solution) then
            do species_idx = ix_na, ix_c
               saha_mass_fractions(species_idx) = envelope_amu_over_amu* &
                    eos_mix%fxenv(species_idx)
            end do
            saha_mass_fractions(ix_h) = (eos_mix%fxenv(ix_h)* &
                 eos_mix%amuenv + h_excess)/ion_mean_weight_inverse
            do species_idx = ix_o, ix_ne
               saha_mass_fractions(species_idx) = envelope_amu_over_amu* &
                    eos_mix%fxenv(species_idx)
            end do
            saha_mass_fractions(ix_he) = (eos_mix%fxenv(ix_he)* &
                 eos_mix%amuenv + y_excess)/ion_mean_weight_inverse
         end if
      else
!        Ion mean weight from the envelope mixture plus the H/he/metal
!        EXCESS (SHARED WITH compute_mean_weights), THEN THE SAHA NUMBER FRACTIONS
         call ion_mean_weight_excess(hydrogen_fraction, metal_fraction, &
              ion_mean_weight_inverse, dfx1, dfx12, dfx4, use_envelope_mixture)
         if (use_envelope_mixture) then
!           Use envelope abundances
            if (need_saha_solution) then
               do species_idx = 1, num_species
                  saha_mass_fractions(species_idx) = &
                       eos_mix%fxenv(species_idx)
               end do
            end if
         else
            amu_inverse = 1.0d0/ion_mean_weight_inverse
            if (need_saha_solution) then
               envelope_amu_frac = eos_mix%amuenv*amu_inverse
               do species_idx = 1, num_species
                  saha_mass_fractions(species_idx) = envelope_amu_frac* &
                       eos_mix%fxenv(species_idx)
               end do
               saha_mass_fractions(ix_c) = saha_mass_fractions(ix_c) + &
                    dfx12*amu_inverse
               saha_mass_fractions(ix_h) = saha_mass_fractions(ix_h) + &
                    dfx1*amu_inverse
               saha_mass_fractions(ix_he) = saha_mass_fractions(ix_he) + &
                    dfx4*amu_inverse
            end if
         end if
      end if
!     Compute radiation pressure
      temperature = exp(ln10*log10_temperature)
      pressure = exp(ln10*log10_pressure)
      beta14 = radiation_constant_over_3*(temperature**2)**2/pressure
      beta = 1 - beta14

      if (beta < 0d0) then
         write(run_log_unit,*) &
              'blend_eos: BETA is negative: TL,PL,BETA,1-BETA', &
              log10_temperature, log10_pressure, beta, 1d0 - beta
      end if

      beta_inverse = 1.0d0/beta
      beta14 = 4.0d0*beta_inverse*min(1.0d0, beta14)
      electron_degeneracy_parameter = 0.0d0
      if (.not.need_saha_solution) then
!        Compute values for fully ionized gas and return
         electron_mean_weight_inverse = hydrogen_fraction*inverse_atomic_weights(1) &
              + (1.0d0 - hydrogen_fraction)*inverse_atomic_weights(4)
         specific_gas_constant = gas_constant*(ion_mean_weight_inverse + &
              electron_mean_weight_inverse)
         ion_fraction(1) = 1.0d0
         ion_fraction(2) = 0.0d0
         ion_fraction(3) = 1.0d0
         call compute_fully_ionized_eos(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, beta, &
              ion_mean_weight_inverse, electron_mean_weight_inverse, &
              electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp)
      else
!     Check if Saumon, Chabrier, and Van Horn equation of state needed.
!     THIS EOS REPLACES THE CALL TO compute_saha_eos, EXCEPT FOR DERIVATIVE PURPOSES.
      if (eos_settings%use_scv_eos) then
         call compute_saha_eos(saha_mass_fractions, log10_temperature, temperature, &
              log10_pressure, pressure, log10_density, density, beta, &
              beta_inverse, beta14, ion_fraction, specific_gas_constant, &
              ion_mean_weight_inverse, electron_mean_weight_inverse, &
              want_derivatives, in_atmosphere, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, saha_state)
         call compute_scv_envelope_eos(log10_temperature, temperature, pressure, &
              log10_density, density, hydrogen_fraction, metal_fraction, &
              beta, ion_fraction, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, valid_table_point, ierr)
         if (ierr /= 0) return

! The SCV path takes its derivatives from compute_saha_eos (analytic).

      else
!        CALL TO PRATHER EOS - Either because SCV was not requested or
!        it has failed.
         call compute_saha_eos(saha_mass_fractions, log10_temperature, temperature, &
              log10_pressure, pressure, log10_density, density, beta, &
              beta_inverse, beta14, ion_fraction, specific_gas_constant, &
              ion_mean_weight_inverse, electron_mean_weight_inverse, &
              want_derivatives, in_atmosphere, dlnrho_dlnt, dlnrho_dlnp, &
              specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
              dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, saha_state)
      end if
!     STORE compute_saha_eos VALUES FOR INTERPOLATION WITH compute_fully_ionized_eos VALUES
      saha_log10_density = log10_density
      saha_dlnrho_dlnt = dlnrho_dlnt
      saha_dlnrho_dlnp = dlnrho_dlnp
      saha_specific_heat_cp = specific_heat_cp
      saha_adiabatic_gradient = adiabatic_gradient
      if (want_derivatives) then
         saha_dlnrho_dlnt_dt = dlnrho_dlnt_dt
         saha_dlnrho_dlnp_dt = dlnrho_dlnp_dt
         saha_specific_heat_cp_dt = specific_heat_cp_dt
         saha_specific_heat_cp_dp = specific_heat_cp_dp
         saha_adiabatic_gradient_dt = adiabatic_gradient_dt
         saha_adiabatic_gradient_dp = adiabatic_gradient_dp
      end if
!     Bypass relativistic EOS (compute_fully_ionized_eos) if in low temperature region
      if (.not. skip_relativistic_eos) then
!     Compute values for fully ionized gas
      electron_mean_weight_inverse = hydrogen_fraction*inverse_atomic_weights(1) + &
           (1.0d0 - hydrogen_fraction)*inverse_atomic_weights(4)
      specific_gas_constant = gas_constant*(ion_mean_weight_inverse + &
           electron_mean_weight_inverse)
      ion_fraction(1) = 1.0d0
      ion_fraction(2) = 0.0d0
      ion_fraction(3) = 1.0d0
      call compute_fully_ionized_eos(log10_temperature, temperature, log10_pressure, &
           pressure, log10_density, density, beta, &
           ion_mean_weight_inverse, electron_mean_weight_inverse, &
           electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
           specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
           dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
           specific_heat_cp_dt, specific_heat_cp_dp)
!     Interpolate values
      ramp_weight = saha_ramp_scale*(eos_settings%saha_log10t_cutoff - log10_temperature)
      ramp_weight_sq = ramp_weight*ramp_weight
      ramp_factor = ramp_weight_sq*ramp_weight* &
           (6.0d0*ramp_weight_sq - 15.0d0*ramp_weight + 10.0d0)
      saha_log10_density = saha_log10_density - log10_density
      saha_dlnrho_dlnt = saha_dlnrho_dlnt - dlnrho_dlnt
      saha_dlnrho_dlnp = saha_dlnrho_dlnp - dlnrho_dlnp
      saha_specific_heat_cp = saha_specific_heat_cp - specific_heat_cp
      saha_adiabatic_gradient = saha_adiabatic_gradient - adiabatic_gradient
      saha_specific_heat_cp_dt = saha_specific_heat_cp_dt - specific_heat_cp_dt
      saha_specific_heat_cp_dp = saha_specific_heat_cp_dp - specific_heat_cp_dp
      saha_adiabatic_gradient_dt = saha_adiabatic_gradient_dt - adiabatic_gradient_dt
      saha_adiabatic_gradient_dp = saha_adiabatic_gradient_dp - adiabatic_gradient_dp
      log10_density = log10_density + ramp_factor*saha_log10_density
      density = exp(ln10*log10_density)
      specific_gas_constant = beta*pressure/(density*temperature)
      electron_mean_weight_inverse = specific_gas_constant/gas_constant - &
           ion_mean_weight_inverse
      ion_fraction(1) = 1.0d0 + ramp_factor*(ion_fraction(1) - 1.0d0)
      ion_fraction(2) = ramp_factor*ion_fraction(2)
      ion_fraction(3) = 1.0d0 + ramp_factor*(ion_fraction(3) - 1.0d0)
      dlnrho_dlnt = dlnrho_dlnt + ramp_factor*saha_dlnrho_dlnt
      dlnrho_dlnp = dlnrho_dlnp + ramp_factor*saha_dlnrho_dlnp
      specific_heat_cp = specific_heat_cp + ramp_factor*saha_specific_heat_cp
      adiabatic_gradient = adiabatic_gradient + &
           ramp_factor*saha_adiabatic_gradient
      if (want_derivatives) then
         dlnrho_dlnt_dt = dlnrho_dlnt_dt + ramp_factor*saha_dlnrho_dlnt_dt
         dlnrho_dlnp_dt = dlnrho_dlnp_dt + ramp_factor*saha_dlnrho_dlnp_dt
         specific_heat_cp_dt = specific_heat_cp_dt + &
              ramp_factor*saha_specific_heat_cp_dt
         specific_heat_cp_dp = specific_heat_cp_dp + &
              ramp_factor*saha_specific_heat_cp_dp
         adiabatic_gradient_dt = adiabatic_gradient_dt + &
              ramp_factor*saha_adiabatic_gradient_dt
         adiabatic_gradient_dp = adiabatic_gradient_dp + &
              ramp_factor*saha_adiabatic_gradient_dp
      end if
      end if
      end if

!     1995 OPAL eqos
      if (eos_settings%use_opal95_eos) then
      if (temperature >= opal95_t_min .and. log10_temperature <= opal95_log10_t_max .and. &
           log10_density <= opal95_log10_rho_max) then

         call compute_opal95_eos(log10_temperature, temperature, log10_pressure, &
              pressure, opal_log10_density, opal_density, &
              hydrogen_fraction, metal_fraction, opal_beta, &
              opal_beta_inverse, opal_beta14, opal_specific_gas_constant, &
              opal_ion_mean_weight_inverse, &
              opal_electron_mean_weight_inverse, opal_dlnrho_dlnt, &
              opal_dlnrho_dlnp, opal_specific_heat_cp, &
              opal_adiabatic_gradient, ierr, opal_covers_point)
         if (ierr /= 0 .or. .not. opal_covers_point) return

         call find_opal95_ramp(temperature, opal_log10_density, ramp_factor, &
              in_opal_table, needs_ramp, ierr)
         if (ierr /= 0) return

         if (.not.in_opal_table) return

         call adopt_or_blend_opal_result()
      end if
      end if

!     2001 OPAL eqos  LLP 6/17/03
      if (eos_settings%use_opal2001_eos) then
      if (temperature >= opal01_t_min .and. temperature <= opal01_t_max .and. &
           log10_density <= opal01_log10_rho_max) then
         call compute_opal2001_eos(log10_temperature, temperature, log10_pressure, &
              pressure, opal_log10_density, opal_density, &
              hydrogen_fraction, metal_fraction, opal_beta, &
              opal_beta_inverse, opal_beta14, opal_specific_gas_constant, &
              opal_ion_mean_weight_inverse, &
              opal_electron_mean_weight_inverse, opal_dlnrho_dlnt, &
              opal_dlnrho_dlnp, opal_specific_heat_cp, &
              opal_adiabatic_gradient, ierr, opal_covers_point)
         if (ierr /= 0 .or. .not. opal_covers_point) return

         call find_opal2001_ramp(temperature, opal_log10_density, ramp_factor, &
              in_opal_table, needs_ramp, ierr)
         if (ierr /= 0) return
!        find_opal2001_ramp determines whether or not the point is in the OPAL
!        2001 EOS table

         if (.not.in_opal_table) return
!        Use OPAL results if not in (rho,T) regime where ramp needed
         call adopt_or_blend_opal_result()
      end if
      end if

!     2006 OPAL eqos  LLP 10/13/2996
      if (eos_settings%use_opal2006_eos) then
      if (temperature >= opal06_t_min .and. temperature <= opal06_t_max .and. &
           log10_density <= opal06_log10_rho_max) then
         call compute_opal2006_eos(log10_temperature, temperature, log10_pressure, &
              pressure, opal_log10_density, opal_density, &
              hydrogen_fraction, metal_fraction, opal_beta, &
              opal_beta_inverse, opal_beta14, opal_specific_gas_constant, &
              opal_ion_mean_weight_inverse, &
              opal_electron_mean_weight_inverse, opal_dlnrho_dlnt, &
              opal_dlnrho_dlnp, opal_specific_heat_cp, &
              opal_adiabatic_gradient, ierr, opal_covers_point)
         if (ierr /= 0 .or. .not. opal_covers_point) return

         call find_opal2006_ramp(temperature, opal_log10_density, ramp_factor, &
              in_opal_table, needs_ramp)
!        find_opal2006_ramp determines whether or not the point is in the OPAL
!        2006 EOS table.
!        If in_opal_table is true, point is in Opal 2006 table.
!        If needs_ramp is true, point is in ramp area, and ramping is
!        required. ramp_factor is the weight of the point for ramping
!        purposes.
!        Also, to eliminate a point, one can set needs_ramp to true and
!        ramp_factor to zero.

         if (.not.in_opal_table) return

!        Use OPAL 2006 results only if not in (rho,T) regime where
!        ramping is needed
         call adopt_or_blend_opal_result()
      end if
      end if

      return

contains

! Take over the OPAL result when the point is fully inside the OPAL
! table (needs_ramp false), or blend it into the Yale/SCV result with
! weight ramp_factor when the point lies in the ramp zone at the table
! edge. Shared by the 1995, 2001 and 2006 branches above; every
! variable is host-associated from blend_eos.
subroutine adopt_or_blend_opal_result()
         if (.not.needs_ramp) then
            log10_density = opal_log10_density
            density = exp10(log10_density)
            beta = opal_beta
            beta_inverse = 1.0d0/beta
            beta14 = 1.0d0 - beta
            specific_gas_constant = opal_specific_gas_constant
            ion_mean_weight_inverse = opal_ion_mean_weight_inverse
            electron_mean_weight_inverse = opal_electron_mean_weight_inverse
            dlnrho_dlnt = opal_dlnrho_dlnt
            dlnrho_dlnp = opal_dlnrho_dlnp
            specific_heat_cp = opal_specific_heat_cp
            adiabatic_gradient = opal_adiabatic_gradient
         else
            log10_density = log10_density + &
                 ramp_factor*(opal_log10_density - log10_density)
            density = exp10(log10_density)
            beta = beta + ramp_factor*(opal_beta - beta)
            beta_inverse = 1.0d0/beta
            beta14 = 1.0d0 - beta
            specific_gas_constant = specific_gas_constant + &
                 ramp_factor*(opal_specific_gas_constant - specific_gas_constant)
            ion_mean_weight_inverse = ion_mean_weight_inverse + &
                 ramp_factor*(opal_ion_mean_weight_inverse - ion_mean_weight_inverse)
            electron_mean_weight_inverse = electron_mean_weight_inverse + &
                 ramp_factor*(opal_electron_mean_weight_inverse - electron_mean_weight_inverse)
            dlnrho_dlnt = dlnrho_dlnt + &
                 ramp_factor*(opal_dlnrho_dlnt - dlnrho_dlnt)
            dlnrho_dlnp = dlnrho_dlnp + &
                 ramp_factor*(opal_dlnrho_dlnp - dlnrho_dlnp)
            specific_heat_cp = specific_heat_cp + &
                 ramp_factor*(opal_specific_heat_cp - specific_heat_cp)
            adiabatic_gradient = adiabatic_gradient + &
                 ramp_factor*(opal_adiabatic_gradient - adiabatic_gradient)
         end if
end subroutine adopt_or_blend_opal_result

end subroutine blend_eos

end module blended_eos

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
      use eos_results, only: eos_result
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
! The blended EOS at (log T, log P, X, Z): blend_eos once. With
! want_derivatives and use_numerical_derivatives, the six second
! derivatives are instead taken by centred differences of blend_eos's
! first derivatives, four extra evaluations at log T +- dtl and
! log P +- dpl (LLP 10-22-06, 8-5-07):
!   dlnrho_dlnt_dt = d(dln rho/dln T)/dln T, dlnrho_dlnt_dp the same
!   in P; adiabatic_gradient_dt/_dp and specific_heat_cp_dt/_dp the log
!   derivatives of grad_ad and cp.
! The shifted evaluations share beta, the ionization fractions and the
! mean weights with r (they are overwritten by the final call), carry
! their own density guess from one to the next, and use a copy of
! saha_state.
subroutine compute_blended_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
     metal_fraction, r, want_derivatives, in_atmosphere, saha_state, ierr)
      real(dp), intent(in) :: log10_temperature, log10_pressure
      real(dp), intent(in) :: hydrogen_fraction, metal_fraction
      type(eos_result), intent(inout) :: r
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state
      integer, intent(out) :: ierr
! the steps in log T and log P: a hundredth of the tables' approximate
! intervals (0.03, 0.15), which worked best empirically
      real(dp), parameter :: dtl = 0.030d0*0.01d0, dpl = 0.150d0*0.01d0
      real(dp), parameter :: dtl2 = 2.0d0*dtl, dpl2 = 2.0d0*dpl
! the shifted evaluations: their record, their density guess, the first
! derivatives (dln rho/dln T, dln rho/dln P, cp, grad_ad) above and
! below, and their unused second derivatives
      type(eos_result) :: shifted
      real(dp) :: shifted_log10_density, shifted_density
      real(dp) :: above(4), below(4), unused(6)
      real(dp) :: kept(6)
      integer :: saha_state_local

      ierr = 0
      if (want_derivatives .and. eos_settings%use_numerical_derivatives) then
         saha_state_local = saha_state
         shifted_log10_density = r%log10_density
         shifted_density = r%density
         above = 0.0d0
         below = 0.0d0
         unused = 0.0d0
         call evaluate_shifted(log10_temperature + dtl, log10_pressure, above)
         if (ierr /= 0) return
         call evaluate_shifted(log10_temperature - dtl, log10_pressure, below)
         if (ierr /= 0) return
         r%dlnrho_dlnt_dt = (above(1) - below(1))/dtl2/ln10
         r%specific_heat_cp_dt = (log10(above(3)) - log10(below(3)))/dtl2
         r%adiabatic_gradient_dt = (log10(above(4)) - log10(below(4)))/dtl2
         call evaluate_shifted(log10_temperature, log10_pressure + dpl, above)
         if (ierr /= 0) return
         call evaluate_shifted(log10_temperature, log10_pressure - dpl, below)
         if (ierr /= 0) return
         r%dlnrho_dlnt_dp = (above(1) - below(1))/dpl2/ln10
         r%specific_heat_cp_dp = (log10(above(3)) - log10(below(3)))/dpl2
         r%adiabatic_gradient_dp = (log10(above(4)) - log10(below(4)))/dpl2
      end if

      if (want_derivatives .and. .not. eos_settings%use_numerical_derivatives) then
         call blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, r, .true., in_atmosphere, saha_state, ierr)
      else
! no analytic second derivatives: r keeps the ones it has (the
! numerical ones, or the caller's)
         kept = [r%dlnrho_dlnt_dt, r%dlnrho_dlnt_dp, r%adiabatic_gradient_dt, &
              r%adiabatic_gradient_dp, r%specific_heat_cp_dt, r%specific_heat_cp_dp]
         call blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
              metal_fraction, r, .false., in_atmosphere, saha_state, ierr)
         r%dlnrho_dlnt_dt = kept(1)
         r%dlnrho_dlnt_dp = kept(2)
         r%adiabatic_gradient_dt = kept(3)
         r%adiabatic_gradient_dp = kept(4)
         r%specific_heat_cp_dt = kept(5)
         r%specific_heat_cp_dp = kept(6)
      end if

contains

! One shifted evaluation (no second derivatives). first carries this
! side's first derivatives in and out.
subroutine evaluate_shifted(lt, lp, first)
      real(dp), intent(in) :: lt, lp
      real(dp), intent(inout) :: first(4)
      shifted = r
      shifted%log10_density = shifted_log10_density
      shifted%density = shifted_density
      shifted%dlnrho_dlnt = first(1)
      shifted%dlnrho_dlnp = first(2)
      shifted%specific_heat_cp = first(3)
      shifted%adiabatic_gradient = first(4)
      shifted%dlnrho_dlnt_dt = unused(1)
      shifted%dlnrho_dlnt_dp = unused(2)
      shifted%adiabatic_gradient_dt = unused(3)
      shifted%adiabatic_gradient_dp = unused(4)
      shifted%specific_heat_cp_dt = unused(5)
      shifted%specific_heat_cp_dp = unused(6)
      call blend_eos(lt, lp, hydrogen_fraction, metal_fraction, shifted, &
           .false., in_atmosphere, saha_state_local, ierr)
      r%temperature = shifted%temperature
      r%pressure = shifted%pressure
      r%beta = shifted%beta
      r%beta_inverse = shifted%beta_inverse
      r%beta14 = shifted%beta14
      r%ion_fraction = shifted%ion_fraction
      r%specific_gas_constant = shifted%specific_gas_constant
      r%ion_mean_weight_inverse = shifted%ion_mean_weight_inverse
      r%electron_mean_weight_inverse = shifted%electron_mean_weight_inverse
      r%electron_degeneracy_parameter = shifted%electron_degeneracy_parameter
      shifted_log10_density = shifted%log10_density
      shifted_density = shifted%density
      first = [shifted%dlnrho_dlnt, shifted%dlnrho_dlnp, shifted%specific_heat_cp, &
           shifted%adiabatic_gradient]
      unused = [shifted%dlnrho_dlnt_dt, shifted%dlnrho_dlnt_dp, shifted%adiabatic_gradient_dt, &
           shifted%adiabatic_gradient_dp, shifted%specific_heat_cp_dt, shifted%specific_heat_cp_dp]
end subroutine evaluate_shifted

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
subroutine blend_eos(log10_temperature, log10_pressure, hydrogen_fraction, &
     metal_fraction, r, want_derivatives, in_atmosphere, saha_state, ierr)

      real(dp), intent(in) :: log10_temperature
      type(eos_result), intent(inout) :: r
      real(dp), intent(in) :: log10_pressure
      real(dp), intent(in) :: hydrogen_fraction
      real(dp), intent(in) :: metal_fraction
      logical, intent(in) :: want_derivatives, in_atmosphere
      integer, intent(inout) :: saha_state
! --- locals ---
      integer :: species_idx
      integer, parameter :: num_species = 12
      real(dp), parameter :: inverse_atomic_weights(4) = [0.9921d0, 0.24975d0, &
           0.08322d0, 0.4995d0]
      real(dp), parameter :: saha_ramp_width = 0.500d0, saha_ramp_scale = 2.000d0
! atomic weights in the ix_* species order (Na, Al, Mg, Fe, Si, C, H,
! O, N, Ar, Ne, He)
      real(dp), parameter :: atomic_weights_full(n_mix_species) = [23.0d0, &
           26.99d0, 24.32d0, 55.86d0, 28.1d0, 12.015d0, 1.008d0, 16.0d0, &
           14.01d0, 39.96d0, 20.19d0, 4.004d0]

      logical :: need_saha_solution, skip_relativistic_eos
      real(dp) :: saha_mass_fractions(n_mix_species)
      real(dp) :: metal_ratio, amu_correction, h_excess, y_excess
      real(dp) :: envelope_amu_over_amu
      real(dp) :: dfx1, dfx12, dfx4, amu_inverse, envelope_amu_frac
      logical :: use_envelope_mixture
      logical :: valid_table_point
! the OPAL result, adopted or blended in (its T and P are copied to r
! right after each OPAL call, as the call used to write them there)
      type(eos_result) :: opal_r
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
      real(dp) :: saha_dlnrho_dlnt_dt, saha_dlnrho_dlnt_dp, &
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
         r%ion_mean_weight_inverse = eos_mix%amuenv
         do species_idx = ix_na, ix_c
            r%ion_mean_weight_inverse = r%ion_mean_weight_inverse + &
                 amu_correction*eos_mix%fxenv(species_idx)
         end do
         h_excess = (hydrogen_fraction - eos_mix%envelope_hydrogen_fraction)/ &
              atomic_weights_full(ix_h)
         r%ion_mean_weight_inverse = r%ion_mean_weight_inverse + h_excess
         do species_idx = ix_o, ix_ne
            r%ion_mean_weight_inverse = r%ion_mean_weight_inverse + &
                 amu_correction*eos_mix%fxenv(species_idx)
         end do
         y_excess = (eos_mix%envelope_hydrogen_fraction + eos_mix%envelope_metal_fraction - &
              hydrogen_fraction - metal_fraction)/atomic_weights_full(ix_he)
         r%ion_mean_weight_inverse = r%ion_mean_weight_inverse + y_excess
         envelope_amu_over_amu = metal_ratio*eos_mix%amuenv/ &
              r%ion_mean_weight_inverse
         if (need_saha_solution) then
            do species_idx = ix_na, ix_c
               saha_mass_fractions(species_idx) = envelope_amu_over_amu* &
                    eos_mix%fxenv(species_idx)
            end do
            saha_mass_fractions(ix_h) = (eos_mix%fxenv(ix_h)* &
                 eos_mix%amuenv + h_excess)/r%ion_mean_weight_inverse
            do species_idx = ix_o, ix_ne
               saha_mass_fractions(species_idx) = envelope_amu_over_amu* &
                    eos_mix%fxenv(species_idx)
            end do
            saha_mass_fractions(ix_he) = (eos_mix%fxenv(ix_he)* &
                 eos_mix%amuenv + y_excess)/r%ion_mean_weight_inverse
         end if
      else
!        Ion mean weight from the envelope mixture plus the H/he/metal
!        EXCESS (SHARED WITH compute_mean_weights), THEN THE SAHA NUMBER FRACTIONS
         call ion_mean_weight_excess(hydrogen_fraction, metal_fraction, &
              r%ion_mean_weight_inverse, dfx1, dfx12, dfx4, use_envelope_mixture)
         if (use_envelope_mixture) then
!           Use envelope abundances
            if (need_saha_solution) then
               do species_idx = 1, num_species
                  saha_mass_fractions(species_idx) = &
                       eos_mix%fxenv(species_idx)
               end do
            end if
         else
            amu_inverse = 1.0d0/r%ion_mean_weight_inverse
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
      r%temperature = exp(ln10*log10_temperature)
      r%pressure = exp(ln10*log10_pressure)
      r%beta14 = radiation_constant_over_3*(r%temperature**2)**2/r%pressure
      r%beta = 1 - r%beta14

      if (r%beta < 0d0) then
         write(run_log_unit,*) &
              'EQSTAT2: BETA is negative: TL,PL,BETA,1-BETA', &
              log10_temperature, log10_pressure, r%beta, 1d0 - r%beta
      end if

      r%beta_inverse = 1.0d0/r%beta
      r%beta14 = 4.0d0*r%beta_inverse*min(1.0d0, r%beta14)
      r%electron_degeneracy_parameter = 0.0d0
      if (.not.need_saha_solution) then
!        Compute values for fully ionized gas and return
         r%electron_mean_weight_inverse = hydrogen_fraction*inverse_atomic_weights(1) &
              + (1.0d0 - hydrogen_fraction)*inverse_atomic_weights(4)
         r%specific_gas_constant = gas_constant*(r%ion_mean_weight_inverse + &
              r%electron_mean_weight_inverse)
         r%ion_fraction(1) = 1.0d0
         r%ion_fraction(2) = 0.0d0
         r%ion_fraction(3) = 1.0d0
         call compute_fully_ionized_eos(log10_temperature, log10_pressure, r)
      else
!     Check if Saumon, Chabrier, and Van Horn equation of state needed.
!     THIS EOS REPLACES THE CALL TO compute_saha_eos, EXCEPT FOR DERIVATIVE PURPOSES.
      if (eos_settings%use_scv_eos) then
         call compute_saha_eos(saha_mass_fractions, log10_temperature, &
              log10_pressure, r, want_derivatives, in_atmosphere, &
              saha_state)
         call compute_scv_envelope_eos(log10_temperature, hydrogen_fraction, &
              metal_fraction, r, valid_table_point, ierr)
         if (ierr /= 0) return
! The SCV path takes its derivatives from compute_saha_eos (analytic).
      else
!        CALL TO PRATHER EOS - Either because SCV was not requested or
!        it has failed.
         call compute_saha_eos(saha_mass_fractions, log10_temperature, &
              log10_pressure, r, want_derivatives, in_atmosphere, &
              saha_state)
      end if
!     STORE compute_saha_eos VALUES FOR INTERPOLATION WITH compute_fully_ionized_eos VALUES
      saha_log10_density = r%log10_density
      saha_dlnrho_dlnt = r%dlnrho_dlnt
      saha_dlnrho_dlnp = r%dlnrho_dlnp
      saha_specific_heat_cp = r%specific_heat_cp
      saha_adiabatic_gradient = r%adiabatic_gradient
      if (want_derivatives) then
         saha_dlnrho_dlnt_dt = r%dlnrho_dlnt_dt
         saha_dlnrho_dlnt_dp = r%dlnrho_dlnt_dp
         saha_specific_heat_cp_dt = r%specific_heat_cp_dt
         saha_specific_heat_cp_dp = r%specific_heat_cp_dp
         saha_adiabatic_gradient_dt = r%adiabatic_gradient_dt
         saha_adiabatic_gradient_dp = r%adiabatic_gradient_dp
      end if
!     Bypass relativistic EOS (compute_fully_ionized_eos) if in low temperature region
      if (.not. skip_relativistic_eos) then
!     Compute values for fully ionized gas
      r%electron_mean_weight_inverse = hydrogen_fraction*inverse_atomic_weights(1) + &
           (1.0d0 - hydrogen_fraction)*inverse_atomic_weights(4)
      r%specific_gas_constant = gas_constant*(r%ion_mean_weight_inverse + &
           r%electron_mean_weight_inverse)
      r%ion_fraction(1) = 1.0d0
      r%ion_fraction(2) = 0.0d0
      r%ion_fraction(3) = 1.0d0
      call compute_fully_ionized_eos(log10_temperature, log10_pressure, r)
!     Interpolate values
      ramp_weight = saha_ramp_scale*(eos_settings%saha_log10t_cutoff - log10_temperature)
      ramp_weight_sq = ramp_weight*ramp_weight
      ramp_factor = ramp_weight_sq*ramp_weight* &
           (6.0d0*ramp_weight_sq - 15.0d0*ramp_weight + 10.0d0)
      saha_log10_density = saha_log10_density - r%log10_density
      saha_dlnrho_dlnt = saha_dlnrho_dlnt - r%dlnrho_dlnt
      saha_dlnrho_dlnp = saha_dlnrho_dlnp - r%dlnrho_dlnp
      saha_specific_heat_cp = saha_specific_heat_cp - r%specific_heat_cp
      saha_adiabatic_gradient = saha_adiabatic_gradient - r%adiabatic_gradient
      saha_specific_heat_cp_dt = saha_specific_heat_cp_dt - r%specific_heat_cp_dt
      saha_specific_heat_cp_dp = saha_specific_heat_cp_dp - r%specific_heat_cp_dp
      saha_adiabatic_gradient_dt = saha_adiabatic_gradient_dt - r%adiabatic_gradient_dt
      saha_adiabatic_gradient_dp = saha_adiabatic_gradient_dp - r%adiabatic_gradient_dp
      r%log10_density = r%log10_density + ramp_factor*saha_log10_density
      r%density = exp(ln10*r%log10_density)
      r%specific_gas_constant = r%beta*r%pressure/(r%density*r%temperature)
      r%electron_mean_weight_inverse = r%specific_gas_constant/gas_constant - &
           r%ion_mean_weight_inverse
      r%ion_fraction(1) = 1.0d0 + ramp_factor*(r%ion_fraction(1) - 1.0d0)
      r%ion_fraction(2) = ramp_factor*r%ion_fraction(2)
      r%ion_fraction(3) = 1.0d0 + ramp_factor*(r%ion_fraction(3) - 1.0d0)
      r%dlnrho_dlnt = r%dlnrho_dlnt + ramp_factor*saha_dlnrho_dlnt
      r%dlnrho_dlnp = r%dlnrho_dlnp + ramp_factor*saha_dlnrho_dlnp
      r%specific_heat_cp = r%specific_heat_cp + ramp_factor*saha_specific_heat_cp
      r%adiabatic_gradient = r%adiabatic_gradient + &
           ramp_factor*saha_adiabatic_gradient
      if (want_derivatives) then
         r%dlnrho_dlnt_dt = r%dlnrho_dlnt_dt + ramp_factor*saha_dlnrho_dlnt_dt
         r%dlnrho_dlnt_dp = r%dlnrho_dlnt_dp + ramp_factor*saha_dlnrho_dlnt_dp
         r%specific_heat_cp_dt = r%specific_heat_cp_dt + &
              ramp_factor*saha_specific_heat_cp_dt
         r%specific_heat_cp_dp = r%specific_heat_cp_dp + &
              ramp_factor*saha_specific_heat_cp_dp
         r%adiabatic_gradient_dt = r%adiabatic_gradient_dt + &
              ramp_factor*saha_adiabatic_gradient_dt
         r%adiabatic_gradient_dp = r%adiabatic_gradient_dp + &
              ramp_factor*saha_adiabatic_gradient_dp
      end if
      end if
      end if

!     1995 OPAL eqos
      if (eos_settings%use_opal95_eos) then
      if (r%temperature >= opal95_t_min .and. log10_temperature <= opal95_log10_t_max .and. &
           r%log10_density <= opal95_log10_rho_max) then

         call compute_opal95_eos(log10_temperature, log10_pressure, &
              hydrogen_fraction, metal_fraction, opal_r, ierr, &
              opal_covers_point)
         r%temperature = opal_r%temperature
         r%pressure = opal_r%pressure
         if (ierr /= 0 .or. .not. opal_covers_point) return

         call find_opal95_ramp(r%temperature, opal_r%log10_density, ramp_factor, &
              in_opal_table, needs_ramp, ierr)
         if (ierr /= 0) return

         if (.not.in_opal_table) return

         call adopt_or_blend_opal_result()
      end if
      end if

!     2001 OPAL eqos  LLP 6/17/03
      if (eos_settings%use_opal2001_eos) then
      if (r%temperature >= opal01_t_min .and. r%temperature <= opal01_t_max .and. &
           r%log10_density <= opal01_log10_rho_max) then
         call compute_opal2001_eos(log10_temperature, log10_pressure, &
              hydrogen_fraction, metal_fraction, opal_r, ierr, &
              opal_covers_point)
         r%temperature = opal_r%temperature
         r%pressure = opal_r%pressure
         if (ierr /= 0 .or. .not. opal_covers_point) return

         call find_opal2001_ramp(r%temperature, opal_r%log10_density, ramp_factor, &
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
      if (r%temperature >= opal06_t_min .and. r%temperature <= opal06_t_max .and. &
           r%log10_density <= opal06_log10_rho_max) then
         call compute_opal2006_eos(log10_temperature, log10_pressure, &
              hydrogen_fraction, metal_fraction, opal_r, ierr, &
              opal_covers_point)
         r%temperature = opal_r%temperature
         r%pressure = opal_r%pressure
         if (ierr /= 0 .or. .not. opal_covers_point) return

         call find_opal2006_ramp(r%temperature, opal_r%log10_density, ramp_factor, &
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
            r%log10_density = opal_r%log10_density
            r%density = exp10(r%log10_density)
            r%beta = opal_r%beta
            r%beta_inverse = 1.0d0/r%beta
            r%beta14 = 1.0d0 - r%beta
            r%specific_gas_constant = opal_r%specific_gas_constant
            r%ion_mean_weight_inverse = opal_r%ion_mean_weight_inverse
            r%electron_mean_weight_inverse = opal_r%electron_mean_weight_inverse
            r%dlnrho_dlnt = opal_r%dlnrho_dlnt
            r%dlnrho_dlnp = opal_r%dlnrho_dlnp
            r%specific_heat_cp = opal_r%specific_heat_cp
            r%adiabatic_gradient = opal_r%adiabatic_gradient
         else
            r%log10_density = r%log10_density + &
                 ramp_factor*(opal_r%log10_density - r%log10_density)
            r%density = exp10(r%log10_density)
            r%beta = r%beta + ramp_factor*(opal_r%beta - r%beta)
            r%beta_inverse = 1.0d0/r%beta
            r%beta14 = 1.0d0 - r%beta
            r%specific_gas_constant = r%specific_gas_constant + &
                 ramp_factor*(opal_r%specific_gas_constant - r%specific_gas_constant)
            r%ion_mean_weight_inverse = r%ion_mean_weight_inverse + &
                 ramp_factor*(opal_r%ion_mean_weight_inverse - r%ion_mean_weight_inverse)
            r%electron_mean_weight_inverse = r%electron_mean_weight_inverse + &
                 ramp_factor*(opal_r%electron_mean_weight_inverse - r%electron_mean_weight_inverse)
            r%dlnrho_dlnt = r%dlnrho_dlnt + &
                 ramp_factor*(opal_r%dlnrho_dlnt - r%dlnrho_dlnt)
            r%dlnrho_dlnp = r%dlnrho_dlnp + &
                 ramp_factor*(opal_r%dlnrho_dlnp - r%dlnrho_dlnp)
            r%specific_heat_cp = r%specific_heat_cp + &
                 ramp_factor*(opal_r%specific_heat_cp - r%specific_heat_cp)
            r%adiabatic_gradient = r%adiabatic_gradient + &
                 ramp_factor*(opal_r%adiabatic_gradient - r%adiabatic_gradient)
         end if
end subroutine adopt_or_blend_opal_result

end subroutine blend_eos

end module blended_eos

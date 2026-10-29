!----------------------------------------------------------------------
! eos_mixture
!----------------------------------------------------------------------
! The envelope mixture the EOS evaluates with (eos_mix): the envelope X
! and Z, the inverse mean atomic weight sum_i X_i/A_i (1/mu_ion of the
! envelope) and the 12 species' number fractions, read by the Saha and
! SCV paths of blend_eos and by compute_mean_weights.
! The star layer sets it through eos_set_mixture whenever it
! recomputes the envelope composition. Also here:
!     ion_mean_weight_excess  1/mu_ion of a mixture with other X and Z
!     compute_mean_weights    1/mu_ion, 1/mu_e and the gas constant of
!                             a zone
!     mean_weights_from_density  R/mu and 1/mu_e that a result's density
!                             implies
module eos_mixture
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use eos_results, only: eos_result
      use composition_layout, only: ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, &
           ix_h, ix_o, ix_n, ix_ar, ix_ne, ix_he, n_mix_species
      use physical_constants, only: gas_constant
      implicit none
      private

! The 12-species surface mixture layout (composition_layout's ix_*,
! re-exported): element per slot of number_fractions, of
! mixture_atomic_weights and of the saha_number_fractions vector
! blend_eos hands to eos/yale/compute_saha_eos.
      public :: ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, ix_h, ix_o, ix_n, &
           ix_ar, ix_ne, ix_he, n_mix_species

      type, public :: eos_mixture_state
            real(dp) :: envelope_hydrogen_fraction, &
                 envelope_metal_fraction
! sum_i X_i/A_i over the envelope's species [1/amu]
            real(dp) :: inverse_atomic_weight
! number_fractions(i): number fraction of species i in the envelope mixture,
! in the ix_* order above (the metals Na..Si, C, then H, O..Ne, He;
! slots ix_h and ix_he are set from X and Y directly, the metals from
! the mixture_weights_seed control scaled to Z).
            real(dp) :: number_fractions(n_mix_species)
      end type eos_mixture_state

      type(eos_mixture_state), public, save :: eos_mix

! the atomic weights of the mixture's species, in the ix_* order (Na,
! Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He): those core's starting model
! builds number_fractions with (X_i/A_i over sum_j X_j/A_j)
      real(dp), parameter, public :: mixture_atomic_weights(n_mix_species) = [23.0d0, &
           26.99d0, 24.32d0, 55.86d0, 28.1d0, 12.015d0, 1.008d0, 16.0d0, &
           14.01d0, 39.96d0, 20.19d0, 4.004d0]

      public :: ion_mean_weight_excess, compute_mean_weights, mean_weights_from_density

contains

! ion_mean_weight_excess: the inverse ion mean molecular weight of a
! mixture with hydrogen fraction x and metal fraction z, as the envelope
! value eos_mix%inverse_atomic_weight plus the excess of each of H, He and metals
! (excess metals taken to be carbon-12) over the envelope mixture.
! use_envelope is .true. when x and z are within 1e-5 of the envelope
! values: then mu_ion_inv = inverse_atomic_weight and dfx1/dfx12 hold the raw
! (unweighted) differences; otherwise dfx1/dfx12/dfx4 are the
! weighted H/metal/He excesses that blend_eos's zone_mixture also needs
! for the Saha number fractions (dfx4 is defined only in that branch).
! Shared by compute_mean_weights and zone_mixture.
pure subroutine ion_mean_weight_excess(x, z, mu_ion_inv, dfx1, dfx12, dfx4, &
     use_envelope)
      real(dp), intent(in) :: x, z
      real(dp), intent(out) :: mu_ion_inv, dfx1, dfx12, dfx4
      logical, intent(out) :: use_envelope
! inverse_atomic_weights: 1/A for H, He and C(12), and Z/A = 1/2 of He
! and the metals (blend_eos's 1/mu_e; unused here)
      real(dp), parameter :: inverse_atomic_weights(4) = [0.9921d0, 0.24975d0, 0.08322d0, 0.4995d0]

! Set up fractional abundances
      dfx1 = (x - eos_mix%envelope_hydrogen_fraction)
      dfx12 = (z - eos_mix%envelope_metal_fraction)
      use_envelope = abs(dfx1) + abs(dfx12) < 1.0d-5
      if (use_envelope) then
! Use envelope abundances
         mu_ion_inv = eos_mix%inverse_atomic_weight
      else
         dfx1 = dfx1*inverse_atomic_weights(1)
         dfx12 = dfx12*inverse_atomic_weights(3)
         dfx4 = (eos_mix%envelope_hydrogen_fraction + eos_mix%envelope_metal_fraction - &
              x - z)*inverse_atomic_weights(2)
! Assume excess Z(metals) is in the form of carbon(12)
         mu_ion_inv = eos_mix%inverse_atomic_weight + dfx1 + dfx4 + dfx12
      end if
end subroutine ion_mean_weight_excess

! compute_mean_weights: the ion and electron inverse mean molecular
! weights and the specific gas constant R*(1/mu_ion + 1/mu_e) of a
! zone for the OPAL and MHD paths: 1/mu_ion from ion_mean_weight_excess,
! and 1/mu_e backed out of the ideal-gas law with the given beta (the
! gas-pressure fraction). So the R it returns equals beta P/(rho T) for
! any input: MHD's comparison of the two catches only a table result
! inconsistent with its own P (a stale or empty row), not a wrong mu.
pure subroutine compute_mean_weights(temperature, pressure, density, hydrogen_fraction, &
     metal_fraction, specific_gas_constant, ion_mean_weight_inverse, &
     electron_mean_weight_inverse, beta)


      real(dp), intent(in) :: temperature, pressure, density, &
           hydrogen_fraction, metal_fraction, beta
      real(dp), intent(out) :: specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse

! dfx1/dfx12/dfx4: per-species excesses returned by
! ion_mean_weight_excess; only ion_mean_weight_inverse is needed here.
      real(dp) :: dfx1, dfx12, dfx4, ee
      logical :: use_envelope

! Ion mean weight from the envelope mixture plus the H/he/metal excess
      call ion_mean_weight_excess(hydrogen_fraction, metal_fraction, &
           ion_mean_weight_inverse, dfx1, dfx12, dfx4, use_envelope)
      ee = ((beta*pressure)/(density*temperature*gas_constant* &
           ion_mean_weight_inverse)) - 1.0d0
      electron_mean_weight_inverse = ee*ion_mean_weight_inverse
      specific_gas_constant = gas_constant*(ion_mean_weight_inverse + &
           electron_mean_weight_inverse)

end subroutine compute_mean_weights

! mean_weights_from_density: the mean weights of a result whose density
! comes from a table or a blend (SCV, the Saha ramp) rather than from
! the composition: R/mu = beta P/(rho T), the gas constant the
! ideal-gas law needs for r's density, and 1/mu_e = (R/mu)/R - 1/mu_ion,
! r's 1/mu_ion kept. (The Saha solution's own R/mu is this by
! construction; OPAL's and MHD's come from compute_mean_weights, which
! backs 1/mu_e out the same way.)
pure subroutine mean_weights_from_density(r)
      type(eos_result), intent(inout) :: r
      r%specific_gas_constant = r%beta*r%pressure/(r%density*r%temperature)
      r%electron_mean_weight_inverse = r%specific_gas_constant/gas_constant - &
           r%ion_mean_weight_inverse
end subroutine mean_weights_from_density

end module eos_mixture

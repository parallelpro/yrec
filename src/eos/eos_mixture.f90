!----------------------------------------------------------------------
! eos_mixture
!----------------------------------------------------------------------
! The envelope mixture the EOS evaluates with (eos_mix): the envelope X
! and Z, the mean atomic mass and the 12 species number fractions read
! by the Saha and SCV paths of blend_eos and by compute_mean_weights.
! The star layer sets it through eos_set_mixture whenever it
! recomputes the envelope composition. Also here:
!     ion_mean_weight_excess  1/mu_ion of a mixture with other X and Z
!     compute_mean_weights    1/mu_ion, 1/mu_e and the gas constant of
!                             a zone
module eos_mixture
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use composition_layout, only: ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, &
           ix_h, ix_o, ix_n, ix_ar, ix_ne, ix_he, n_mix_species
      use physical_constants, only: gas_constant
      implicit none
      private

! The 12-species surface mixture layout (composition_layout's ix_*,
! re-exported): element per slot of fxenv, blend_eos's
! atomic_weights_full table and the saha_mass_fractions vector handed
! to eos/yale/compute_saha_eos.
      public :: ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, ix_h, ix_o, ix_n, &
           ix_ar, ix_ne, ix_he, n_mix_species

      type, public :: eos_mixture_state
            real(dp) :: envelope_hydrogen_fraction, &
                 envelope_metal_fraction
            real(dp) :: amuenv
! fxenv(i): number fraction of species i in the envelope mixture,
! in the ix_* order above (the metals Na..Si, C, then H, O..Ne, He;
! slots ix_h and ix_he are set from X and Y directly, the metals from
! the mixture_weights_seed control scaled to Z).
            real(dp) :: fxenv(n_mix_species)
      end type eos_mixture_state

      type(eos_mixture_state), public, save :: eos_mix

      public :: ion_mean_weight_excess, compute_mean_weights

contains

! ion_mean_weight_excess: the inverse ion mean molecular weight of a
! mixture with hydrogen fraction x and metal fraction z, as the envelope
! value eos_mix%amuenv plus the excess of each of H, He and metals
! (excess metals taken to be carbon-12) over the envelope mixture.
! use_envelope is .true. when x and z are within 1e-5 of the envelope
! values: then mu_ion_inv = amuenv and dfx1/dfx12 hold the raw
! (unweighted) differences; otherwise dfx1/dfx12/dfx4 are the
! weighted H/metal/He excesses that compute_blended_eos also needs for its Saha
! number fractions (dfx4 is defined only in that branch).  Shared by
! compute_mean_weights and compute_blended_eos.
pure subroutine ion_mean_weight_excess(x, z, mu_ion_inv, dfx1, dfx12, dfx4, &
     use_envelope)
      real(dp), intent(in) :: x, z
      real(dp), intent(out) :: mu_ion_inv, dfx1, dfx12, dfx4
      logical, intent(out) :: use_envelope
! inverse_atomic_weights: 1/A for H, He, C(12), and a spare slot
      real(dp), parameter :: inverse_atomic_weights(4) = [0.9921d0, 0.24975d0, 0.08322d0, 0.4995d0]

! Set up fractional abundances
      dfx1 = (x - eos_mix%envelope_hydrogen_fraction)
      dfx12 = (z - eos_mix%envelope_metal_fraction)
      use_envelope = abs(dfx1) + abs(dfx12) < 1.0d-5
      if (use_envelope) then
! Use envelope abundances
         mu_ion_inv = eos_mix%amuenv
      else
         dfx1 = dfx1*inverse_atomic_weights(1)
         dfx12 = dfx12*inverse_atomic_weights(3)
         dfx4 = (eos_mix%envelope_hydrogen_fraction + eos_mix%envelope_metal_fraction - &
              x - z)*inverse_atomic_weights(2)
! Assume excess Z(metals) is in the form of carbon(12)
         mu_ion_inv = eos_mix%amuenv + dfx1 + dfx4 + dfx12
      end if
end subroutine ion_mean_weight_excess

! compute_mean_weights: the ion and electron inverse mean molecular
! weights and the specific gas constant R*(1/mu_ion + 1/mu_e) of a
! zone, as a cross-check on an EOS table's own bookkeeping (the OPAL
! and MHD paths call it). 1/mu_ion comes from ion_mean_weight_excess;
! 1/mu_e is backed out of the ideal-gas law with the given beta (the
! gas-pressure fraction).
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

end module eos_mixture

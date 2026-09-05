!----------------------------------------------------------------------
! eos_mixture_lib
!----------------------------------------------------------------------
! Added 2026 (physics-purity pass -- see ROADMAP.md "Decoupling the
! physics domains from the model"). The envelope mixture the eos
! domain evaluates with: surface hydrogen/metal fractions, the mean
! atomic mass, and the 12-species number fractions consumed by the
! Yale (eqrelv/eqsaha) and SCV branches of eqstat2 and by eos/mu.f90.
!
! Previously the eos internals read this straight from the model
! (star%...), the one place the eos domain touched
! star_info. Now the STAR LAYER pushes the mixture through the
! eos_set_mixture facade entry (eos_lib) whenever it recomputes the
! envelope composition (core/read_starting_model.f90's mixture blocks,
! setup/rescale_model.f90's rescale), and the eos domain reads only this,
! its own state -- same pattern as every other physics-domain state
! module. MESA-equivalently: composition is an INPUT of the eos,
! owned on the physics side of the boundary.
module eos_mixture_lib
      implicit none
      private

! The 12-species surface mixture layout: element per slot of fxenv,
! eqstat2's atomic_weights_full table and the saha_mass_fractions
! vector handed to eos/yale/saha_eos. Same values as star_info_lib's
! ix_* / n_mix_species (the star layer fills fxenv in that order in
! core/read_starting_model.f90's update_surface_mixture); a copy
! because this module must not depend on star_info_lib (physics-purity
! pass), like net_lib's n_reactions.
      integer, parameter, public :: &
           ix_na = 1, ix_al = 2, ix_mg = 3, ix_fe = 4, ix_si = 5, &
           ix_c = 6, ix_h = 7, ix_o = 8, ix_n = 9, ix_ar = 10, &
           ix_ne = 11, ix_he = 12, n_mix_species = 12

      type, public :: eos_mixture_state
            double precision :: envelope_hydrogen_fraction, &
                 envelope_metal_fraction
            double precision :: amuenv
! fxenv(i): number fraction of species i in the envelope mixture,
! in the ix_* order above (the metals Na..Si, C, then H, O..Ne, He;
! slots ix_h and ix_he are set from X and Y directly, the metals from
! the mixture_weights_seed control scaled to Z).
            double precision :: fxenv(n_mix_species)
      end type eos_mixture_state

      type(eos_mixture_state), public, save :: eos_mix

      public :: ion_mean_weight_excess

contains

! ion_mean_weight_excess: the inverse ion mean molecular weight of a
! mixture with hydrogen fraction x and metal fraction z, as the envelope
! value eos_mix%amuenv plus the excess of each of H, He and metals
! (excess metals taken to be carbon-12) over the envelope mixture.
! use_envelope is .true. when x and z are within 1e-5 of the envelope
! values: then mu_ion_inv = amuenv and dfx1/dfx12 hold the raw
! (unweighted) differences; otherwise dfx1/dfx12/dfx4 are the
! weighted H/metal/He excesses that eqstat also needs for its Saha
! number fractions (dfx4 is defined only in that branch).  Shared by
! eos/mu.f90 and eos/eqstat.f90 since 2026 W2; statement order is
! that of mu.f90.
subroutine ion_mean_weight_excess(x, z, mu_ion_inv, dfx1, dfx12, dfx4, &
     use_envelope)
      implicit none
      double precision, intent(in) :: x, z
      double precision, intent(out) :: mu_ion_inv, dfx1, dfx12, dfx4
      logical, intent(out) :: use_envelope
! inverse_atomic_weights: 1/A for H, He, C(12), and a spare slot
      double precision :: inverse_atomic_weights(4)
      data inverse_atomic_weights/0.9921d0, 0.24975d0, 0.08322d0, 0.4995d0/

! SET UP FRACTIONAL ABUNDANCES
      dfx1 = (x - eos_mix%envelope_hydrogen_fraction)
      dfx12 = (z - eos_mix%envelope_metal_fraction)
      use_envelope = dabs(dfx1) + dabs(dfx12) .lt. 1.0d-5
      if (use_envelope) then
! USE ENVELOPE ABUNDANCES
         mu_ion_inv = eos_mix%amuenv
      else
         dfx1 = dfx1*inverse_atomic_weights(1)
         dfx12 = dfx12*inverse_atomic_weights(3)
         dfx4 = (eos_mix%envelope_hydrogen_fraction + eos_mix%envelope_metal_fraction - &
              x - z)*inverse_atomic_weights(2)
! ASSUME EXCESS Z(METALS) IS IN THE FORM OF CARBON(12)
         mu_ion_inv = eos_mix%amuenv + dfx1 + dfx4 + dfx12
      end if
end subroutine ion_mean_weight_excess

end module eos_mixture_lib

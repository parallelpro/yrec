!----------------------------------------------------------------------
! mu
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original mu.f; only variable names, source form, and comment style
! were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! Mean molecular weight routine, called from meqos.f90/oeqos.f90/
! oeqos01.f90/oeqos06.f90 as a cross-check on the EOS's own mean-weight
! bookkeeping. Given T, P, density, X, Z, and beta (gas pressure
! fraction), returns the specific gas constant R*(mu_ion^-1+mu_e^-1),
! the ion mean-molecular-weight inverse, and the electron
! mean-molecular-weight inverse (backed out from the ideal-gas law
! using the supplied beta).
subroutine mu(temperature, pressure, density, hydrogen_fraction, &
     metal_fraction, specific_gas_constant, ion_mean_weight_inverse, &
     electron_mean_weight_inverse, beta)

      use phys_const_lib
      use eos_mixture_lib, only: eos_mix
      implicit none

      double precision, intent(in) :: temperature, pressure, density, &
           hydrogen_fraction, metal_fraction, beta
      double precision, intent(out) :: specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse

      double precision :: inverse_atomic_weights(4)
      data inverse_atomic_weights/0.9921d0, 0.24975d0, 0.08322d0, 0.4995d0/
      double precision :: dfx1, dfx12, dfx4, ee

! SET UP FRACTIONAL ABUNDANCES
      dfx1 = (hydrogen_fraction - eos_mix%envelope_hydrogen_fraction)
      dfx12 = (metal_fraction - eos_mix%envelope_metal_fraction)
      if (dabs(dfx1) + dabs(dfx12) .lt. 1.0d-5) then
! USE ENVELOPE ABUNDANCES
         ion_mean_weight_inverse = eos_mix%amuenv
      else
         dfx1 = dfx1*inverse_atomic_weights(1)
         dfx12 = dfx12*inverse_atomic_weights(3)
         dfx4 = (eos_mix%envelope_hydrogen_fraction + eos_mix%envelope_metal_fraction - &
              hydrogen_fraction - metal_fraction)*inverse_atomic_weights(2)
! ASSUME EXCESS Z(METALS) IS IN THE FORM OF CARBON(12)
         ion_mean_weight_inverse = eos_mix%amuenv + dfx1 + dfx4 + dfx12
      end if
      ee = ((beta*pressure)/(density*temperature*gas_constant* &
           ion_mean_weight_inverse)) - 1.0d0
      electron_mean_weight_inverse = ee*ion_mean_weight_inverse
      specific_gas_constant = gas_constant*(ion_mean_weight_inverse + &
           electron_mean_weight_inverse)

      return
end subroutine mu

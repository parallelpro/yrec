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
      use eos_mixture_lib, only: ion_mean_weight_excess
      implicit none

      double precision, intent(in) :: temperature, pressure, density, &
           hydrogen_fraction, metal_fraction, beta
      double precision, intent(out) :: specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse

! dfx1/dfx12/dfx4: per-species excesses returned by
! ion_mean_weight_excess; only ion_mean_weight_inverse is needed here.
      double precision :: dfx1, dfx12, dfx4, ee
      logical :: use_envelope

! ION MEAN WEIGHT FROM THE ENVELOPE MIXTURE PLUS THE H/He/METAL EXCESS
      call ion_mean_weight_excess(hydrogen_fraction, metal_fraction, &
           ion_mean_weight_inverse, dfx1, dfx12, dfx4, use_envelope)
      ee = ((beta*pressure)/(density*temperature*gas_constant* &
           ion_mean_weight_inverse)) - 1.0d0
      electron_mean_weight_inverse = ee*ion_mean_weight_inverse
      specific_gas_constant = gas_constant*(ion_mean_weight_inverse + &
           electron_mean_weight_inverse)

      return
end subroutine mu

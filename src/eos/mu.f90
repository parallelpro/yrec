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

      use const_lib
      implicit none

      double precision, intent(in) :: temperature, pressure, density, &
           hydrogen_fraction, metal_fraction, beta
      double precision, intent(out) :: specific_gas_constant, &
           ion_mean_weight_inverse, electron_mean_weight_inverse

! common/comp/: only envelope_hydrogen_fraction/envelope_metal_fraction/
! envelope_amu are used here. The remaining members are declared only
! to preserve the storage layout shared with every other file that
! references common/comp/. Naming matches getopac.f90/meqos.f90.
      double precision :: envelope_hydrogen_fraction, &
           envelope_metal_fraction, zenvm, envelope_amu, &
           envelope_species_fractions(12), xnew, znew, stotal, senv
      common/comp/ envelope_hydrogen_fraction, envelope_metal_fraction, &
           zenvm, envelope_amu, envelope_species_fractions, xnew, znew, &
           stotal, senv




! common/ccout2/: not used in this file; declared only to preserve
! layout. Naming matches meqos.f90.
      logical :: ldebug, lcorr, lmilne, ltrack, lstpch
      common/ccout2/ ldebug, lcorr, lmilne, ltrack, lstpch

! DATA NZP1/12/
      double precision :: atomic_weights(4)
      data atomic_weights/0.9921d0, 0.24975d0, 0.08322d0, 0.4995d0/
      save

      double precision :: dfx1, dfx12, dfx4, ee

! SET UP FRACTIONAL ABUNDANCES
      dfx1 = (hydrogen_fraction - envelope_hydrogen_fraction)
      dfx12 = (metal_fraction - envelope_metal_fraction)
      if (dabs(dfx1) + dabs(dfx12) .lt. 1.0d-5) then
! USE ENVELOPE ABUNDANCES
         ion_mean_weight_inverse = envelope_amu
      else
         dfx1 = dfx1*atomic_weights(1)
         dfx12 = dfx12*atomic_weights(3)
         dfx4 = (envelope_hydrogen_fraction + envelope_metal_fraction - &
              hydrogen_fraction - metal_fraction)*atomic_weights(2)
! ASSUME EXCESS Z(METALS) IS IN THE FORM OF CARBON(12)
         ion_mean_weight_inverse = envelope_amu + dfx1 + dfx4 + dfx12
      end if
      ee = ((beta*pressure)/(density*temperature*gas_constant* &
           ion_mean_weight_inverse)) - 1.0d0
      electron_mean_weight_inverse = ee*ion_mean_weight_inverse
      specific_gas_constant = gas_constant*(ion_mean_weight_inverse + &
           electron_mean_weight_inverse)

      return
end subroutine mu

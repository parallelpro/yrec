!----------------------------------------------------------------------
! eos_lib
!----------------------------------------------------------------------
! Added 2026 as part of the YREC readability refactor's phase two
! (disentangling the solver from the physics domains -- see
! GUIDELINES.md's "Physics domains still entangled with the solver").
!
! eos_get is the first facade of this phase: a single explicit-
! interface entry point that replaces the `if (use_mhd_eos) call meqos
! else call eqstat` dispatch that used to be duplicated at every call
! site, and additionally centralizes the Debye-Huckel composition
! setup (debye_huckel_x/y/z_total/z(3), read only by eos/yale/eqrelv.f90
! on the non-MHD path) that 6 of those call sites also used to
! duplicate verbatim. eos/eqstat.f90 (which also hosts eqstat2, its
! co-located pair), eos/yale/eqrelv.f90, and eos/mhd/meqos.f90 are all
! unchanged; this is a pure dispatch/setup wrapper.
!
! composition_at_zone is OPTIONAL: callers that carry full per-species
! composition (misc/coefft.f90, io/wrtout.f90, core/starin.f90,
! misc/physic.f90, mixing/hsubp.f90, mixing/sconvec.f90) pass
! composition(:,idx) and get the same Debye-Huckel setup they used to
! compute themselves. Callers that only ever carried bulk
! hydrogen_fraction/metal_fraction (atm/envint.f90, atm/qatm.f90,
! atm/qenv.f90, wind/massloss.f90) omit it -- exactly matching their
! prior behavior, since none of them set these fields before either.
!
! core/starin.f90 previously had a bug here (confirmed against the
! original F77 source via git history): a missing ELSE meant it called
! *both* meqos and eqstat when MHD was on, and *neither* when MHD was
! off. eos_get's if/else has the structurally-correct form, so
! migrating that call site fixes the bug by construction.
!
! mixing/hsubp.f90, mixing/sconvec.f90, and wind/massloss.f90
! previously called eqstat unconditionally, with no LMHD check at all
! -- confirmed authentic original YREC behavior (unchanged since the
! very first commit, not a modernization artifact). Migrating them to
! eos_get extends real MHD support to these three secondary/diagnostic
! calculations for the first time, per explicit user sign-off; this is
! an acknowledged numerics change for use_mhd_eos=.true. runs, which
! the Stage-0 regression suite cannot verify (no test case sets LMHD).
!
! wind/calcad.f90 is NOT migrated: it calls esac06 directly (bypassing
! eqstat2's boundary-ramping) under its own use_opal2006_eos check for
! a self-contained acoustic-depth diagnostic, and never checks
! use_mhd_eos at all -- confirmed to match the original F77 exactly,
! a deliberate design choice, not part of this dispatch pattern.
module eos_lib
      implicit none
contains

subroutine eos_get(log10_temperature, temperature, log10_pressure, &
     pressure, log10_density, density, hydrogen_fraction, metal_fraction, &
     beta, beta_inverse, beta14, ion_fraction, specific_gas_constant, &
     ion_mean_weight_inverse, electron_mean_weight_inverse, &
     electron_degeneracy_parameter, dlnrho_dlnt, dlnrho_dlnp, &
     specific_heat_cp, adiabatic_gradient, dlnrho_dlnt_dt, &
     dlnrho_dlnp_dt, adiabatic_gradient_dt, adiabatic_gradient_dp, &
     specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
     in_atmosphere, saha_state, composition_at_zone)

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
      double precision, intent(in), optional :: composition_at_zone(15)

      if (use_mhd_eos) then
         call meqos(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp)
      else
         if (present(composition_at_zone) .and. use_debye_huckel_correction) then
            debye_huckel_x = composition_at_zone(1)
            debye_huckel_y = composition_at_zone(2)+composition_at_zone(4)
            debye_huckel_z_total = composition_at_zone(3)
            debye_huckel_z(1) = composition_at_zone(5)+composition_at_zone(6)
            debye_huckel_z(2) = composition_at_zone(7)+composition_at_zone(8)
            debye_huckel_z(3) = composition_at_zone(9)+composition_at_zone(10)+ &
                 composition_at_zone(11)
         end if
         call eqstat(log10_temperature, temperature, log10_pressure, &
              pressure, log10_density, density, hydrogen_fraction, &
              metal_fraction, beta, beta_inverse, beta14, ion_fraction, &
              specific_gas_constant, ion_mean_weight_inverse, &
              electron_mean_weight_inverse, electron_degeneracy_parameter, &
              dlnrho_dlnt, dlnrho_dlnp, specific_heat_cp, &
              adiabatic_gradient, dlnrho_dlnt_dt, dlnrho_dlnp_dt, &
              adiabatic_gradient_dt, adiabatic_gradient_dp, &
              specific_heat_cp_dt, specific_heat_cp_dp, want_derivatives, &
              in_atmosphere, saha_state)
      end if

      return
end subroutine eos_get

end module eos_lib

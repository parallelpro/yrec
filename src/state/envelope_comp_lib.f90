!----------------------------------------------------------------------
! envelope_comp_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/comp/: envelope composition/mass
! state, recomputed every timestep (core/starin.f90, setup/hpoint.f90,
! setup/rscale.f90, wind/mdot.f90) and read broadly across eos/,
! mixing/, atm/, and elsewhere -- not set-once startup config despite
! initial appearances (envelope_hydrogen_fraction/zenvm/amuenv are
! recomputed from the current envelope composition in
! core/starin.f90; senv/stotal track the current envelope/total mass
! in log space).
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! oldmod_lib/scrtch_lib/etc: genuinely evolving per-model state, read
! from distant, unrelated points in the call graph. Field names are
! unchanged from the original COMMON member names, matching that same
! precedent.
module envelope_comp_lib
      implicit none
      private

      type, public :: envelope_composition_state
            double precision :: envelope_hydrogen_fraction, &
                 envelope_metal_fraction
            double precision :: zenvm, amuenv, fxenv(12)
            double precision :: xnew, znew, stotal, senv
      end type envelope_composition_state
! 2026 (phase six, step 1 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module envelope_comp_lib

!----------------------------------------------------------------------
! fluxes_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/fluxes/: solar neutrino flux
! diagnostics, zeroed and rebuilt from scratch every Newton iteration
! by misc/coefft.f90 (when lsnu is active) -- the same cadence as
! scrtch_lib -- and read by the output writers (io/wrtout.f90,
! io/wrtmonte.f90) and core/main.f90 (which also has a second,
! permanently-disabled `compute_neutrino_fluxes = .false.` block that
! recomputes neutrino_flux_total/neutrino_flux_zone; left converted
! rather than dropped since it still compiles and the surrounding live
! code in the same file reads these fields).
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! scrtch_lib/turnover_lib: genuinely evolving per-model state, read
! from distant, unrelated points in the call graph. Field names are
! unchanged from the original COMMON member names, matching that same
! precedent.
module fluxes_lib
      implicit none
      private

      type, public :: neutrino_flux_state
            double precision :: neutrino_flux(10), neutrino_flux_total(10)
            double precision :: cl37_snu_rate, ga71_snu_rate
      end type neutrino_flux_state
! 2026 (phase six, step 1 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module fluxes_lib

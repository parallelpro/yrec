!----------------------------------------------------------------------
! pulse_diag_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/pulse1/ and common/pulse2/:
! pulsation/diagnostic quantities used by the GYRE-format pulsation
! output feature (io/write_gyre_pulse.f90) and the envelope
! integrators. common/pulse1/'s 9 members are per-shell (json-indexed)
! arrays computed by misc/coefft.f90 every Newton iteration (same
! cadence as scrtch_lib) plus lpumod; common/pulse2/'s 20 members are
! scalars holding the current envelope point's diagnostics, set by
! atm/qatm.f90/atm/qenv.f90 and read by atm/atm_lib.f90/io/wrtmod.f90.
! Bundled into one type since every file using more than one of these
! two former blocks already handles them together.
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! scrtch_lib/mdphy_lib/temp_lib/envstruct_lib: genuinely evolving
! per-model state, read from distant, unrelated points in the call
! graph. Field names are unchanged from the original COMMON member
! names (pulse2's especially cryptic ones kept as-is, matching
! mdphy_lib's own precedent for names with no established readable
! equivalent), matching that same precedent.
module pulse_diag_lib
      implicit none
      private
      integer, parameter :: json = 5000

      type, public :: pulsation_diagnostics_state
! former common/pulse1/
            double precision :: pulse_dlnrho_dlnp(json), &
                 pulse_dlneps_dlnrho(json)
            double precision :: pulse_dlneps_dlnt(json), &
                 pulse_dlnkap_dlnrho(json)
            double precision :: pulse_dlnkap_dlnt(json), &
                 pulse_specific_heat(json)
            double precision :: pulse_mean_molecular_weight(json), &
                 pulse_dlnrho_dlnt(json)
            double precision :: pulse_electron_mean_molecular_weight(json)
            logical :: lpumod
! former common/pulse2/
            double precision :: qqdp, qqed, qqet, qqod, qqot, qdel, qdela, &
                 qqcp
            double precision :: qrmu, qtl, qpl, qdl, qo, qol, qt, qp
            double precision :: qqdt, qemu, qd, qfs
      end type pulsation_diagnostics_state
! 2026 (phase six, step 1 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module pulse_diag_lib

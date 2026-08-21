!----------------------------------------------------------------------
! engeb_diag_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/neweps/, common/be7/, and
! common/grab/: per-model nuclear energy-generation diagnostics, all
! set once per shell by nuclear/engeb.f90 and read from a handful of
! distant callers -- misc/coefft.f90 (alpha-capture energy, neutrino
! loss rate, and the He3+He3/He3+He4 diagnostic rate arrays) and
! util/ytime.f90 (alpha-capture energy, neutrino loss rate), plus
! core/main.f90 (be7_mass_fraction, saved into a per-zone history
! array).
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! oldmod_lib/scrtch_lib/turnover_lib/light_burn_lib: genuinely
! evolving per-model state, read from distant, unrelated points in the
! call graph rather than handed from one caller to one callee. Bundled
! into one type despite the three different original COMMON block
! names since they're all engeb.f90-computed diagnostics of the same
! character, some read together by the same callers (coefft.f90 reads
! both neweps and grab). Field names are unchanged from the original
! COMMON member names, matching that same precedent.
module engeb_diag_lib
      implicit none
      private
      integer, parameter :: json = 5000

      type, public :: engeb_diagnostics_state
! former common/neweps/
            double precision :: alpha_capture_energy, neutrino_loss_rate
! former common/be7/
            double precision :: be7_mass_fraction
! former common/grab/: He3+He3/He3+He4 luminosity and per-shell rate
! diagnostics (JVS 10/11).
            double precision :: he3_luminosity_placeholder, &
                 he3_total_placeholder
            double precision :: he3_he3_rate_placeholder(json), &
                 he3_he4_rate_placeholder(json)
      end type engeb_diagnostics_state

      type(engeb_diagnostics_state), public :: engeb_diag

end module engeb_diag_lib

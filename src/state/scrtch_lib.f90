!----------------------------------------------------------------------
! scrtch_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/scrtch/: per-shell physics
! diagnostics computed by the Henyey coefficient-builder
! (misc/coefft.f90) every Newton iteration -- opacity, the nuclear
! energy-generation breakdown, gradients, electron degeneracy,
! convective velocity, specific heat, ionization fractions -- and read
! by essentially every output writer (io/wrtout.f90, io/putstore.f90,
! io/write_gyre_pulse.f90, ...). The single most central, most-
! frequently-touched piece of per-shell state in the whole solver.
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same
! as common/oldmod/: genuinely evolving state (recomputed every model,
! not set once), but still module territory since it's read from many
! distant, unrelated points in the call graph, not handed cleanly from
! one caller to one callee.
!
! A real derived type, not bare module variables, for the same reason
! as oldmod_lib: unambiguous provenance at every use site
! (`shell_diag%so` vs. a bare `so`), and a bundle that's straightforward
! to iterate/serialize later for profile/pulse output. Field names are
! unchanged from the original COMMON member names (so/sesum/seg/etc)
! rather than renamed to something more readable -- this pass is about
! adding structure/type-safety, not renaming; that was a separate,
! earlier phase of this project. Single module-level instance
! (shell_diag), not multi-instance -- same reasoning as oldmod_lib.
!
! Distinct from oldmod_lib's prev_model_state: that type holds the
! primary structure variables snapshotted from the *previous*
! timestep; this one holds *derived diagnostics* of the *current*
! (in-progress or just-converged) model. Different fields, different
! purpose, hence its own type.
module scrtch_lib
      implicit none
      private
      integer, parameter :: json = 5000

      type, public :: shell_diagnostics_state
            double precision :: sesum(json), seg(7,json), sbeta(json), &
                 seta(json)
            logical :: locons(json)
            double precision :: so(json), del_grad(3,json), &
                 sfxion(3,json), svel(json), scp(json)
      end type shell_diagnostics_state
! 2026 (phase four, step 4 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module scrtch_lib

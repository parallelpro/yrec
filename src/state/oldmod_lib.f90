!----------------------------------------------------------------------
! oldmod_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/oldmod/: a snapshot of the previous
! converged model's structure, kept for timestep-to-timestep
! comparisons (structure-change-based timestep control, deuterium/
! light-element burning rate changes, remeshing, etc). Written from
! several distinct places (core/starin.f90, core/main.f90,
! setup/hpoint.f90) and read from many more -- genuinely shared,
! evolving global state, not fixed configuration -- but per
! GUIDELINES.md's module-vs-argument test, this is still module
! territory: nothing here flows cleanly from one caller to one callee
! within a single call sequence the way common/tridi/ did, it's global
! state touched from many distant, unrelated points in the call graph.
!
! Deliberately a real derived type (MESA's prev_star_info pattern)
! rather than bare module variables like const_lib/luout_lib: the
! `prev_model%` prefix makes it unambiguous at every use site that a
! given field is specifically previous-model state (not, say, a
! current-model quantity that merely happens to share a name), and
! bundling the fields into one type makes it straightforward to
! iterate/serialize the whole snapshot later (profile/history/pulse
! output). This does NOT introduce multi-instance support (YREC runs
! one star at a time, unlike MESA) -- prev_model below is a single
! module-level instance, accessed the same way common/oldmod/ was,
! just through a named, typed field instead of a bare positional slot.
!
! `json` (the max shell-count array bound) is declared private to this
! module purely to size the type's array components -- every caller
! already declares its own local `integer, parameter :: json = 5000`,
! which would collide with json if it were also exported here.
module oldmod_lib
      implicit none
      private
      integer, parameter :: json = 5000

      type, public :: prev_model_state
            double precision :: old_pressure(json), old_temperature(json), &
                 old_radius(json), old_luminosity(json), old_density(json)
            double precision :: old_composition(15,json)
            double precision :: old_shell_mass(json)
            logical :: old_convective_flag(json), old_cz_flag(json)
            double precision :: old_teff
            integer :: old_num_zones
      end type prev_model_state
! 2026 (phase four, step 4 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module oldmod_lib

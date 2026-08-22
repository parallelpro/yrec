!----------------------------------------------------------------------
! turnover_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces the state half of common/ovrtrn/: the
! convective-turnover-timescale and photospheric-pressure values used
! by the wind-torque and Rossby-number machinery (wind/amcalc.f90,
! wind/mwind.f90, rotation/getw.f90, ...), each held as a current and
! a previous-substep value so callers can linearly interpolate across
! a timestep via fracstep. The other two former common/ovrtrn/
! members, use_new_turnover_timescale and calc_envelope_flag, are
! genuine set-once NAMELIST /physics/ configuration and went to
! const_lib instead (Category 1, not this Category 1b type).
!
! Per GUIDELINES.md's module-vs-argument test this is case 1b, same as
! oldmod_lib/scrtch_lib: genuinely evolving per-model state, read from
! many distant, unrelated points in the call graph rather than handed
! cleanly from one caller to one callee. Field names are unchanged
! from the original COMMON member names, matching that same
! precedent.
module turnover_lib
      implicit none
      private

      type, public :: turnover_state
            double precision :: convective_turnover_timescale, &
                 convective_turnover_timescale_old
            double precision :: pphot, pphot0
            double precision :: fracstep
      end type turnover_state
! 2026 (phase six, step 1 -- ROADMAP.md): the instance moved into
! star_info (state/star_info_lib.f90); this module now only defines
! the type.
end module turnover_lib

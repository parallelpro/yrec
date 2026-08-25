!----------------------------------------------------------------------
! const_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/const1/, common/const2/, and
! common/const3/: physical/mixing-length constants that are set once
! (by setup/setups.f90, at run startup) and read broadly across the
! codebase, never varying per call. This is global configuration, not
! per-call data, so it becomes a module of plain (non-parameter)
! module-level variables rather than subroutine arguments -- matching
! MESA's own const_def/chem_def convention. Every file that used to
! declare any of these three COMMON blocks now does `use const_lib`
! instead; setup/setups.f90's existing assignment statements
! (unchanged) now set these module variables directly instead of the
! old COMMON slots.
!
! const1 members: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
! const2 members: gas_constant, radiation_constant_over_3, ca3l, csig,
!                 csigl, cgl, cmkh, cmkhn
! const3 members: cdelrl, cmixl, cmixl2, cmixl3, clndp,
!                 seconds_per_year -- cmixl (the mixing length) is the
!                 one member here that isn't a pure physical constant:
!                 it's copied from the per-kind-card namelist array
!                 cmixla(nk) each time a new kind card starts, rather
!                 than computed once at the very start of the run like
!                 the rest of const1-3 -- still "occasional
!                 configuration read broadly," not per-call data, so
!                 the same module treatment applies.
! ctlim members: atime, tcut, saha_log10t_cutoff, tenv0, tenv1,
!                tgcut -- NAMELIST /physics/ values (core/parmin.f90),
!                same "set once, read broadly" character as const1-3;
!                the computed member tenv (0.5*(tenv0+tenv1)) moved
!                to star%tenv in the 2026 phase-A eviction.
! 2026 (phase six, step 3 -- ROADMAP.md): const_lib is now a
! compatibility umbrella. Its members were split by ownership --
! physical constants into phys_const_lib, the namelist targets into
! controls_lib, and the table/working-state strays into their domain
! state modules (atm_table_lib, opacity_table_lib, scv_eos_lib,
! yale_eos_lib). The use-and-reexport below keeps every existing
! `use const_lib` working unchanged (bare names resolve through the
! umbrella); new code should use the specific modules, and the
! per-file migration off the umbrella can proceed incrementally.
! The evicted domain-state members are deliberately NOT re-exported:
! their users reference the domain modules directly.
module const_lib
      use phys_const_lib
      use controls_lib
      implicit none
      public
end module const_lib

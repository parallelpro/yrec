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
! ctlim members: atime, tcut, saha_log10t_cutoff, tenv0, tenv1, tenv,
!                tgcut -- NAMELIST /physics/ values (core/parmin.f90),
!                same "set once, read broadly" character as const1-3;
!                tenv is the one computed member (0.5*(tenv0+tenv1)),
!                still only once, not per call.
module const_lib
      implicit none

! former common/const1/
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! former common/const2/
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn

! former common/const3/. cmixl's default (1.4d0) was previously set by
! a DATA statement in core/parmin.f90 (data lkuthe,cmixl/.false.,
! 1.4d0/) -- moved here since parmin.f90 can no longer target a
! use-associated variable with DATA; main.f90 overwrites it with the
! per-kind-card mixing length (cmixla(nk)) before it's ever read for a
! real model, same as before.
      double precision :: cdelrl, cmixl2, cmixl3, clndp, seconds_per_year
      double precision :: cmixl = 1.4d0

! former common/ctlim/. Defaults (previously two DATA statements in
! core/parmin.f90, now illegal there since these are use-associated
! rather than locally declared) moved here as declaration-time
! initializers. tenv is the one member computed at runtime (see
! core/parmin.f90: tenv = 0.5d0*(tenv0+tenv1)), so it has no default.
      double precision :: atime(14) = (/1.0d-3,2.0d-2,5.0d-1,2.0d-2, &
           3.0d-1,1.5d-3,1.0d-1,2.0d-2,4.0d-2,2.0d-2,2.0d-2,0.25d0, &
           1.5d0,0.25d0/)
      double precision :: tcut(5) = (/6.5d0,6.5d0,6.82d0,7.7d0,7.5d0/)
      double precision :: saha_log10t_cutoff = 6.0d0
      double precision :: tenv0 = 3.0d0, tenv1 = 9.0d0, tgcut = 6.9d0
      double precision :: tenv

! former common/cross/, common/weak/. All 6 cross members and
! weak_screening_threshold get their real, final values from
! setup/remap.f90 (called from core/parmin.f90 right after the
! NAMELIST /physics/ read) -- cross_section_scale/qs0e_scale/
! qqs0ee_scale/o16_gamma_scale/c12_alpha_scale are fully recomputed
! there from other inputs regardless of any namelist value; only
! use_new_nuclear_rates and weak_screening_threshold are read
! (unchanged) from the namelist itself (lnewnuc/weakscreening in
! core/parmin.f90, which can't be renamed -- see the naming note at
! that file's top -- so parmin.f90 copies them into these canonical
! names right after the namelist read, before remap runs).
      double precision :: cross_section_scale(17), qs0e_scale(8), &
           qqs0ee_scale(8), o16_gamma_scale, c12_alpha_scale
      logical :: use_new_nuclear_rates
      double precision :: weak_screening_threshold

end module const_lib

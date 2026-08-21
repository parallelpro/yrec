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

end module const_lib

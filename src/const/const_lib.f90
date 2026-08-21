!----------------------------------------------------------------------
! const_lib
!----------------------------------------------------------------------
! New (2026) as part of the YREC module-modernization project (see
! GUIDELINES.md). Replaces common/const1/ and common/const2/: physical
! constants that are set once (by setup/setups.f90, at run startup)
! and read broadly across the codebase, never varying per call. This
! is global configuration, not per-call data, so it becomes a module
! of plain (non-parameter) module-level variables rather than
! subroutine arguments -- matching MESA's own const_def/chem_def
! convention. Every file that used to declare common/const1/ and/or
! common/const2/ now does `use const_lib` instead; setup/setups.f90's
! existing assignment statements (unchanged) now set these module
! variables directly instead of the old COMMON slots.
!
! const1 members: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi
! const2 members: gas_constant, radiation_constant_over_3, ca3l, csig,
!                 csigl, cgl, cmkh, cmkhn
module const_lib
      implicit none

! former common/const1/
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! former common/const2/
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn

end module const_lib

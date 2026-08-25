!----------------------------------------------------------------------
! phys_const_lib
!----------------------------------------------------------------------
! Added 2026 (phase six, step 3 -- ROADMAP.md). The genuinely physical
! and derived constants split out of const_lib: former common/const1/,
! /const2/, /const3/ (set once by setups, read broadly; cmixl's
! per-kind-card caveat documented below, carried from const_lib), and
! the version string. Everything here is immutable after startup
! except cmixl (see its block comment).
module phys_const_lib
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

! former common/version/: yrec_version_string/git_hash_string
! (originally yrecver/githash) are not namelist values -- genuinely
! used in core/parmin.f90, renamed in place there.
      character(len=10) :: yrec_version_string
      character(len=20) :: git_hash_string


end module phys_const_lib

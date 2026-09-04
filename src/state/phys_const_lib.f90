!----------------------------------------------------------------------
! phys_const_lib
!----------------------------------------------------------------------
! Added 2026 (phase six, step 3 -- ROADMAP.md). The genuinely physical
! and derived constants split out of const_lib: former common/const1/,
! /const2/, /const3/ (set once by setups, read broadly), and the
! version string. Everything here is immutable after startup: cmixl,
! the historical mutable exception, moved to
! star%mixing_length_alpha in the 2026 phase-A eviction. (cmixl2/
! cmixl3, despite the names, are MLT formula constants -- cc13 and
! 16*sqrt(2)*sigma -- not powers of the mixing length.)
module phys_const_lib
      implicit none

! former common/const1/
      double precision :: ln10, clni, c4pi, c4pil, c4pi3l, cc13, cc23, cpi

! former common/const2/
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn

! former common/const3/.
      double precision :: cdelrl, cmixl2, cmixl3, clndp, seconds_per_year

! former common/version/: yrec_version_string/git_hash_string
! (originally yrecver/githash) are not namelist values -- genuinely
! used in io/read_controls.f90, renamed in place there.
      character(len=10) :: yrec_version_string
      character(len=20) :: git_hash_string


end module phys_const_lib

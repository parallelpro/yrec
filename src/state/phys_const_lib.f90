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

! ---- named physical constants (2026 readability sweep) ----
! Each keeps exactly the literal value of the site it replaced; the
! *_legacy names mark values that differ from the ones setups.f90
! computes into the variables above (csig) or from the namelist
! defaults (solar_radius_cgs), so the two deliberately coexist.
!
! (10^9 sidereal years / 1 s) * (1 amu / 1 g), the C21 of engeb/rates/
! deutrate: converts a rate per second per gram into the per-Gyr
! per-amu units stored in star%reaction_rate and
! star%deuterium_burning_rate.
      double precision, parameter :: gyr_amu_per_sec_gram = 5.240358d-8
! proton mass [g] (microdiff electron-density estimates)
      double precision, parameter :: m_proton_cgs = 1.6726d-24
! atomic mass unit [g], the value viscos uses for the ion number
! densities (CODATA 1973; differs from the 1.6605402d-24 of the
! rates/engeb mass tables' documentation)
      double precision, parameter :: amu_cgs_legacy = 1.6605655d-24
! electron mass [amu] (the fourth species of the microdiff
! atomic-weight tables)
      double precision, parameter :: m_electron_amu = 5.486d-4
! solar radius [cm] hard-coded by the microdiff/gravitational-
! settling scales (Bahcall-Loeb units); the model itself uses the
! namelist value star%solar_radius_cgs (default also 6.9598d10)
      double precision, parameter :: rsun_cgs_legacy = 6.9598d10
! Stefan-Boltzmann constant [erg/cm^2/s/K^4] as used in viscos's
! radiative thermal diffusivity; differs from csig (5.67051d-5)
      double precision, parameter :: sigma_sb_cgs_legacy = 5.669d-5
! 4a/(15c) (a = radiation constant, c = speed of light, to the
! literal's own precision): coefficient of the Ledoux (1958)
! radiative viscosity nu_rad = 4aT^4/(15 c kappa rho^2) in viscos.f90
      double precision, parameter :: rad_viscosity_coeff_cgs = 6.7282653d-26

! former common/version/: yrec_version_string/git_hash_string
! (originally yrecver/githash) are not namelist values -- genuinely
! used in io/read_controls.f90, renamed in place there.
      character(len=10) :: yrec_version_string
      character(len=20) :: git_hash_string


end module phys_const_lib

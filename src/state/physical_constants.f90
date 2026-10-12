!----------------------------------------------------------------------
! physical_constants
!----------------------------------------------------------------------
! Added 2026 (phase six, step 3 -- docs/ROADMAP.md). The genuinely physical
! and derived constants split out of const_lib: former common/const1/,
! /const2/, /const3/ (set once by setups, read broadly), and the
! version string. The mathematical constants are parameters (2026-10:
! setups assigned them at run time, so any caller of a physics domain
! had to run setups first; each compile-time value is bit-identical to
! the run-time one). The rest is immutable after startup: cmixl,
! the historical mutable exception, moved to
! star%mixing_length_alpha in the 2026 phase-A eviction. (cmixl2/
! cmixl3, despite the names, are MLT formula constants -- cc13 and
! 16*sqrt(2)*sigma -- not powers of the mixing length.)
module physical_constants
      implicit none

! mathematical constants (former common/const1/, and clndp and
! seconds_per_year of /const3/); clndp is ln10's legacy twin
      double precision, parameter :: ln10 = log(10.0d0), clndp = ln10, &
           clni = 1.0d0/ln10, cc13 = 1.0d0/3.0d0, cc23 = cc13 + cc13, &
           cpi = 3.1415926535898d0, c4pi = 4.0d0*cpi, c4pil = log10(c4pi), &
           c4pi3l = log10(cc13*c4pi), seconds_per_year = 3.1558d7

! former common/const2/
      double precision :: gas_constant, radiation_constant_over_3, ca3l, &
           csig, csigl, cgl, cmkh, cmkhn

! former common/const3/.
      double precision :: cdelrl, cmixl2, cmixl3

! ---- named physical constants (2026 readability sweep) ----
! Each keeps exactly the literal value of the site it replaced; the
! *_legacy names mark values that differ from the ones run_setup's set_up_constants_and_tables
! computes into the variables above (csig) or from the namelist
! defaults (solar_radius_cgs), so the two deliberately coexist.
!
! atomic mass unit [g], the value of the network's isotope-mass tables
! (Bahcall 1997): turns a rate per second per gram into the per-second
! per-amu units of the composition rate vector and the deuterium burning
! rate. (2026-10: it replaced gyr_amu_per_sec_gram = 5.240358d-8, the C21
! of engeb/rates/deutrate -- the same factor times 10^9 sidereal years,
! which burned 4.7e-6 more than the timestep's 3.1558d7 s per year.)
      double precision, parameter :: amu_cgs = 1.6605402d-24
! proton mass [g] (microdiff electron-density estimates)
      double precision, parameter :: m_proton_cgs = 1.6726d-24
! atomic mass unit [g], the value viscos uses for the ion number
! densities (CODATA 1973; differs from the 1.6605402d-24 of the
! rates/compute_zone_energy mass tables' documentation)
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

contains

! log10 of the stellar radius (cm) from log10 L/Lsun and log10 Teff,
! L = 4 pi R**2 sigma Teff**4 (the star%log_R_surface expression, before
! the solar-radius offset; formerly wind_lib's and observables's).
! envelope_integration's surface_log10_radius_cm orders the terms
! differently and can differ in the last bit.
pure double precision function log10_radius_from_l_teff(log_luminosity_lsun, log_teff, &
     log10_solar_luminosity)
      double precision, intent(in) :: log_luminosity_lsun, log_teff, log10_solar_luminosity
      log10_radius_from_l_teff = 0.5d0*(log_luminosity_lsun + log10_solar_luminosity &
           - c4pil - csigl - 4.0d0*log_teff)
end function log10_radius_from_l_teff

end module physical_constants

!----------------------------------------------------------------------
! physical_constants
!----------------------------------------------------------------------
! The mathematical and physical constants: compile-time parameters
! (2026-10: setups assigned most of them at run time, so a physics
! domain gave the right answer only after core's setup had run; each
! compile-time value was checked to be bit-identical to the run-time
! one). The one run-time value left is cgl, log10 G, which the control
! G_cgs can set: run_setup assigns it, from log10_newton_g. Also the
! version strings. (cmixl2/cmixl3, despite the names, are MLT formula
! constants -- cc13 and 16*sqrt(2)*sigma -- not powers of the mixing
! length.)
module physical_constants
      implicit none

! mathematical constants (former common/const1/, and clndp and
! seconds_per_year of /const3/); clndp is ln10's legacy twin
      double precision, parameter :: ln10 = log(10.0d0), clndp = ln10, &
           clni = 1.0d0/ln10, cc13 = 1.0d0/3.0d0, cc23 = cc13 + cc13, &
           cpi = 3.1415926535898d0, c4pi = 4.0d0*cpi, c4pil = log10(c4pi), &
           c4pi3l = log10(cc13*c4pi), seconds_per_year = 3.1558d7

! The CODATA 1986 values (Cohen & Taylor 1987), cgs, that the setup has
! always used (former common/const2/ and /const3/): the Stefan-Boltzmann
! constant sigma [erg/(cm^2 s K^4)], the speed of light [cm/s] (exact),
! the molar gas constant R [erg/(mol K)], and the electron mass [g], the
! Boltzmann constant [erg/K] and the Planck constant [erg s] of the Saha
! constant below.
      double precision, parameter :: csig = 5.67051d-5, &
           speed_of_light_cgs = 2.99792458d10, gas_constant = 8.314510d7, &
           electron_mass_cgs = 9.1093897d-28, boltzmann_constant_cgs = 1.380658d-16, &
           planck_constant_cgs = 6.6260755d-27
! log10 sigma, and a/3 = 4 sigma/(3c) [erg/(cm^3 K^4)], with a the
! radiation constant (P_rad = (a/3) T^4)
      double precision, parameter :: csigl = log10(csig), &
           radiation_constant_over_3 = csig*4.0d0/(3.0d0*speed_of_light_cgs)
! the Saha-equation constant cmkh = log10[(2 pi m_e)^(3/2) k^(5/2)/h^3]
      double precision, parameter :: cmkh = 1.50d0*log10(2.0d0*cpi*electron_mass_cgs) + &
           2.5d0*log10(boltzmann_constant_cgs) - 3.0d0*log10(planck_constant_cgs)
! the radiative-gradient constant cdelrl = log10[3/(64 pi sigma)], so
! that grad_rad = 10^cdelrl kappa L P/(G m T^4); and the mixing-length
! constants 1/3 and 16 sqrt(2) sigma
      double precision, parameter :: cdelrl = -c4pi3l - csigl - log10(16.0d0), &
           cmixl2 = cc13, cmixl3 = 16.0d0*sqrt(2.0d0)*csig
! log10 G [cm^3/(g s^2)]: the CODATA 1986 G = 6.67259e-8, its log10
! rounded to five decimals (G = 6.67252e-8); the value unless the
! control G_cgs sets another (log10_newton_g)
      double precision, parameter :: log10_newton_g_default = -7.17571d0
! log10 G for the run: the default until run_setup sets it from the
! control G_cgs (former common/const2/). Only the star layer reads it;
! the domains and the star-blind kernels take log10 G from their
! settings records (log10_gravitational_constant: convection's,
! rotation's, the envelope kernel's and the shell equations'), and
! tools/check_boundaries.py keeps them off this variable.
      double precision :: cgl = log10_newton_g_default

! ---- named physical constants (2026 readability sweep) ----
! Each keeps exactly the literal value of the site it replaced; the
! *_legacy names mark values that differ from the parameters above
! (csig) or from the namelist defaults (solar_radius_cgs), so the two
! deliberately coexist.
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
! (originally yrecver/githash) are not namelist values: set by
! core/run/yrec_version.f90, printed in the controls read's banner
! (io/controls_reader.f90, io/derived_controls.f90) and the history.
      character(len=10) :: yrec_version_string
      character(len=20) :: git_hash_string

contains

! log10 G (cgs) for a G_cgs control: log10 G_cgs when it is positive,
! else the CODATA 1986 value log10_newton_g_default (a non-positive
! G_cgs, the default, keeps the historical log10 G = -7.17571 bit for
! bit). run_setup sets cgl from it, and the settings builders the
! domains' copies.
pure double precision function log10_newton_g(G_cgs)
      use math_lib
      double precision, intent(in) :: G_cgs
      if (G_cgs > 0.0d0) then
         log10_newton_g = log10(G_cgs)
      else
         log10_newton_g = log10_newton_g_default
      end if
end function log10_newton_g

! log10 of the stellar radius (cm) from log10 L/Lsun and log10 Teff,
! L = 4 pi R**2 sigma Teff**4 (the star%log_R_surface_rsun expression, before
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

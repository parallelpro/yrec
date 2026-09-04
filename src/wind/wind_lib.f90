!----------------------------------------------------------------------
! wind_lib
!----------------------------------------------------------------------
! Added 2026 (readability sweep, batch R3): the two expressions that
! every wind routine used to spell out by hand. Each function body is
! token-for-token the expression it replaces, so results are
! bit-identical to the former inline code.
module wind_lib
      implicit none
contains

!----------------------------------------------------------------------
! log10_radius_from_l_teff
!----------------------------------------------------------------------
! log10 of the stellar radius (cm) from log10 L/Lsun and log10 Teff via
! L = 4 pi R**2 sigma Teff**4.
double precision function log10_radius_from_l_teff(log_luminosity_lsun, log_teff)
      use star_info_lib, only: star
      use phys_const_lib
      implicit none
      double precision, intent(in) :: log_luminosity_lsun, log_teff
      log10_radius_from_l_teff = 0.5d0*(log_luminosity_lsun+star%log10_solar_luminosity-c4pil- &
           csigl-4.d0*log_teff)
end function log10_radius_from_l_teff

!----------------------------------------------------------------------
! matt_centrifugal_factor
!----------------------------------------------------------------------
! Centrifugal reduction term of Matt+2012 ApJ 754, L26, implemented
! relative to the Sun: ((c_2**2 + fsun)/(c_2**2 + f))**excen with
! f = omega**2 R**3/(2 G M) evaluated at angular velocity omega. fsun is
! the same quantity for the Sun, computed by the caller.
double precision function matt_centrifugal_factor(omega, fsun, log10_radius, &
     total_mass_msun)
      use star_info_lib, only: star
      use phys_const_lib
      use math_lib
      implicit none
      double precision, intent(in) :: omega, fsun, log10_radius, total_mass_msun
      double precision :: fcorr_local
      fcorr_local = 0.5*omega**2*exp(ln10*(3.0*log10_radius-cgl))/ &
           total_mass_msun/star%solar_mass_cgs
      matt_centrifugal_factor = pow(((star%ctrl%c_2**2+fsun)/(star%ctrl%c_2**2+fcorr_local)), star%ctrl%excen)
end function matt_centrifugal_factor

end module wind_lib

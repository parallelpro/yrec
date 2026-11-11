!----------------------------------------------------------------------
! bahcall_loeb_units
!----------------------------------------------------------------------
! The units of the Bahcall & Loeb (1990) settling equations, in which
! both pipelines solve: a quantity in these units is the CGS value times
! its scale below.
!   radius       the solar radius (bl_radius_scale = 1/R_sun, in 1/cm);
!   mass         100 g/cm^3 R_sun^3 (bl_mass_scale = 1e-2/R_sun^3, in
!                1/g, so that the density unit is 100 g/cm^3);
!   temperature  1e7 K (bl_temp_scale = 1e-7, in 1/K);
!   time         2.7e13 yr (bl_time_scale, in s per unit): Bahcall &
!                Loeb's settling time with the Coulomb logarithm 2.2
!                folded in, so the coefficients divide by ln Lambda/2.2.
! The year is physical_constants' seconds_per_year (3.1558e7 s). Each
! scale is a compile-time parameter; the values are bit for bit those
! the setup once computed at run time (checked 2026-11-10).
module bahcall_loeb_units
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use physical_constants, only: rsun_cgs_legacy, seconds_per_year
      implicit none
      private
      real(dp), parameter, public :: &
           bl_radius_scale = 1.0d0/rsun_cgs_legacy, &
           bl_mass_scale = 1.0d-2*bl_radius_scale**3, &
           bl_temp_scale = 1.0d-7, &
           bl_time_scale = 2.7d13*seconds_per_year
end module bahcall_loeb_units

! bahcall_loeb_units_lib: the Bahcall & Loeb (1990) unit scales used by
! the diffusion/settling routines.  set_bahcall_loeb_scales stores them
! in star%bl_radius_scale / bl_mass_scale / bl_temp_scale / bl_time_scale.
! Written once per settling call from microdiff_setup and
! gravitational_settling_setup (which carried two copies of the same
! statements before 2026 W2); read back by microdiff_etm and
! equal_to_model to undo the in-place conversion of the callers' arrays.
module bahcall_loeb_units_lib
      implicit none
contains

subroutine set_bahcall_loeb_scales()
      use star_info_lib, only: star
      use phys_const_lib, only: rsun_cgs_legacy
      implicit none
!     solar_radius_bl = SOLAR RADIUS (CM)
!     seconds_per_year_bl = NUMBER OF SECONDS IN A YEAR.
      double precision, parameter :: solar_radius_bl = rsun_cgs_legacy
      double precision, parameter :: seconds_per_year_bl = 3.1558d7
!     star%bl_mass_scale=CONVERSION FACTOR FOR MASS.
!     star%bl_radius_scale=CONVERSION FACTOR FOR RADIUS.
!     star%bl_temp_scale=CONVERSION FACOTR FOR TEMPERATURE.
!     star%bl_time_scale=CONVERSION FACTOR FOR TIME.
      star%bl_radius_scale=1.0d0/solar_radius_bl
      star%bl_mass_scale=1.0d-2*star%bl_radius_scale**3
      star%bl_temp_scale=1.0d-7
!     INCLUDES FACTOR OF 2.2 FROM LN LAMBDA
      star%bl_time_scale=2.7d13*seconds_per_year_bl
end subroutine set_bahcall_loeb_scales

end module bahcall_loeb_units_lib

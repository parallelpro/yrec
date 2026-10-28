! bahcall_loeb_units: the Bahcall & Loeb (1990) unit scales of the
! settling routines. set_bahcall_loeb_scales stores them in diff_scr
! (bl_radius_scale, bl_mass_scale, bl_temp_scale, bl_time_scale) once per
! settling call, from microdiff_setup and gravitational_settling_setup;
! microdiff_etm and equal_to_model read them back to undo the in-place
! conversion of the callers' arrays.
module bahcall_loeb_units
      use element_diffusion_scratch, only: diffusion_scratch
      implicit none
contains

subroutine set_bahcall_loeb_scales(diff_scr)
      use physical_constants, only: rsun_cgs_legacy
      implicit none
      type(diffusion_scratch), intent(inout) :: diff_scr
!     solar_radius_bl = SOLAR RADIUS (CM)
!     seconds_per_year_bl = NUMBER OF SECONDS IN A YEAR.
      double precision, parameter :: solar_radius_bl = rsun_cgs_legacy
      double precision, parameter :: seconds_per_year_bl = 3.1558d7
!     diff_scr%bl_mass_scale=CONVERSION FACTOR FOR MASS.
!     diff_scr%bl_radius_scale=CONVERSION FACTOR FOR RADIUS.
!     diff_scr%bl_temp_scale=CONVERSION FACOTR FOR TEMPERATURE.
!     diff_scr%bl_time_scale=CONVERSION FACTOR FOR TIME.
      diff_scr%bl_radius_scale=1.0d0/solar_radius_bl
      diff_scr%bl_mass_scale=1.0d-2*diff_scr%bl_radius_scale**3
      diff_scr%bl_temp_scale=1.0d-7
!     INCLUDES FACTOR OF 2.2 FROM LN LAMBDA
      diff_scr%bl_time_scale=2.7d13*seconds_per_year_bl
end subroutine set_bahcall_loeb_scales

end module bahcall_loeb_units

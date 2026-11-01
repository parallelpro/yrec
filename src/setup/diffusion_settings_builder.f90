!----------------------------------------------------------------------
! diffusion_settings_builder
!----------------------------------------------------------------------
! element diffusion's settings record, built from the run's controls and
! the Monte-Carlo run's metal settling factor (setup/README.md).
! settings_install's install_card_settings installs it at the start of
! every run card (set_diffusion_settings):
!     diffusion_settings_from_star    diffusion_settings_t
module diffusion_settings_builder
      use element_diffusion_lib, only: diffusion_settings_t
      implicit none
      private
      public :: diffusion_settings_from_star

contains

!----------------------------------------------------------------------
! diffusion_settings_from_star: element diffusion's settings from the
! controls (the Monte-Carlo run's fgrz included).
!   star reads: physics%constant_mixing_coeff,
!     physics%constant_settling_reduction, physics%coulomb_log_choice,
!     physics%diffuse_lithium, physics%helium_diffusion_min,
!     physics%hydrogen_diffusion_floor, physics%settling_num_iterations,
!     physics%settling_tolerance, physics%settling_y_factor,
!     physics%use_diffusion_z, physics%use_new_diffusion_routines,
!     physics%use_thoul_diffusion, physics%use_thoul_fit,
!     settling_z_factor_in_use
!   star writes: none
pure function diffusion_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(diffusion_settings_t) :: settings

      settings%use_new_diffusion_routines = s%physics%use_new_diffusion_routines
      settings%use_diffusion_z = s%physics%use_diffusion_z
      settings%diffuse_lithium = s%physics%diffuse_lithium
      settings%fgry = s%physics%settling_y_factor
      settings%fgrz = s%settling_z_factor_in_use
      settings%hydrogen_diffusion_floor = s%physics%hydrogen_diffusion_floor
      settings%helium_diffusion_min = s%physics%helium_diffusion_min
      settings%settling_num_iterations = s%physics%settling_num_iterations
      settings%settling_tolerance = s%physics%settling_tolerance
      settings%constant_mixing_coeff = s%physics%constant_mixing_coeff
      settings%coulomb_log_choice = s%physics%coulomb_log_choice
      settings%use_thoul_fit = s%physics%use_thoul_fit
      settings%use_thoul_diffusion = s%physics%use_thoul_diffusion
      settings%constant_settling_reduction = s%physics%constant_settling_reduction
end function diffusion_settings_from_star

end module diffusion_settings_builder

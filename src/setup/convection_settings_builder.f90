!----------------------------------------------------------------------
! convection_settings_builder
!----------------------------------------------------------------------
! convection's settings record, built from the run's controls, the
! card's mixing length and the solar units (setup/README.md).
! settings_install's install_card_settings installs it at the start of
! every run card (set_convection_settings), and
! install_surface_mixture_settings passes its envelope_metals_no_cno
! again once the first card's model gives it:
!     convection_settings_from_star   convection_settings_t
module convection_settings_builder
      use convection_lib, only: convection_settings_t
      implicit none
      private
      public :: convection_settings_from_star

contains

!----------------------------------------------------------------------
! convection_settings_from_star: convection's settings from the
! controls and the run card's mixing length; log10 G is derived from
! G_cgs (log10_newton_g).
!   star reads: envelope_metals_no_cno, log10_solar_luminosity,
!     mixing_length_alpha, physics%adiabatic_overshoot_active,
!     physics%core_overshoot_active, physics%envelope_overshoot_active,
!     physics%G_cgs, physics%gradient_logt_cutoff,
!     physics%instability_transport_active,
!     physics%intermediate_overshoot_active,
!     physics%overshoot_alpha_core, physics%overshoot_alpha_envelope,
!     physics%overshoot_alpha_min, physics%overshoot_beta_core,
!     physics%overshoot_max_extent_flag, physics%rotation_active,
!     physics%spot_depth_varies, physics%spot_filling_factor,
!     physics%spot_temp_contrast, physics%use_extended_composition,
!     solar_luminosity_cgs
!   star writes: none
pure function convection_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      use physical_constants, only: log10_newton_g
      type(star_info), intent(in) :: s
      type(convection_settings_t) :: settings

      settings%mixing_length_alpha = s%mixing_length_alpha
      settings%log10_t_cut = s%physics%gradient_logt_cutoff
      settings%adiabatic_overshoot = s%physics%adiabatic_overshoot_active
      settings%solar_luminosity_cgs = s%solar_luminosity_cgs
      settings%log10_solar_luminosity = s%log10_solar_luminosity
      settings%log10_gravitational_constant = log10_newton_g(s%physics%G_cgs)
      settings%spot_filling_factor = s%physics%spot_filling_factor
      settings%spot_temp_contrast = s%physics%spot_temp_contrast
      settings%spot_depth_varies = s%physics%spot_depth_varies
      settings%core_overshoot = s%physics%core_overshoot_active
      settings%envelope_overshoot = s%physics%envelope_overshoot_active
      settings%intermediate_overshoot = s%physics%intermediate_overshoot_active
      settings%overshoot_alpha_core = s%physics%overshoot_alpha_core
      settings%overshoot_alpha_envelope = s%physics%overshoot_alpha_envelope
      settings%overshoot_alpha_intermediate = s%physics%overshoot_alpha_min
      settings%limit_core_overshoot = s%physics%overshoot_max_extent_flag
      settings%core_overshoot_limit = s%physics%overshoot_beta_core
      settings%rotation_active = s%physics%rotation_active
      settings%instability_transport_active = s%physics%instability_transport_active
      settings%use_extended_composition = s%physics%use_extended_composition
      settings%envelope_metals_no_cno = s%envelope_metals_no_cno
end function convection_settings_from_star

end module convection_settings_builder

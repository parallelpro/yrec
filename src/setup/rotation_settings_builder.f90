!----------------------------------------------------------------------
! rotation_settings_builder
!----------------------------------------------------------------------
! The rotation domain's settings record, built from the run's controls
! by the star layer; settings_install's install_card_settings installs
! it (set_rotation_settings) at the start of every run card, before the
! starting model is read. For a rotating run it also derives the wind
! law's coefficients (magnetic_braking's derive_wind_law), clamps the
! convection zone's rotation law to [-2, 0] and makes
! force_solid_body_rotation imply solid_body_mode_flag = 1; a run without
! rotation keeps the controls as read and no wind law. log10 G is
! derived from G_cgs (log10_newton_g), as convection's.
module rotation_settings_builder
      use rotation_configuration, only: rotation_settings_t
      use magnetic_braking, only: wind_law_inputs, wind_law_coefficients, derive_wind_law
      implicit none
      private
      public :: rotation_settings_from_star

contains

!----------------------------------------------------------------------
! rotation_settings_from_star: rotation's settings from the controls.
!   star reads: log10_solar_luminosity, physics%chi_grid_scale,
!     physics%constant_background_diffusion_coeff,
!     physics%convergence_tolerance, physics%critical_reynolds,
!     physics%cz_rotation_law_exponent, physics%difad_velocity_scale,
!     physics%diffusion_timestep_factor, physics%disk_locking_active,
!     physics%disk_locking_age_gyr, physics%disk_omega_rad_s,
!     physics%es_mixing_scale, physics%es_velocity_scale,
!     physics%force_solid_body_rotation, physics%fpft_average_weight,
!     physics%fpft_iter1, physics%fpft_iter2, physics%G_cgs,
!     physics%gsf_inhibition_mode, physics%gsf_mixing_scale,
!     physics%gsf_velocity_scale, physics%helium_diffusion_min,
!     physics%hydrogen_diffusion_floor,
!     physics%instability_transport_active,
!     physics%intermediate_overshoot_active,
!     physics%max_diffusion_iters, physics%max_domega_global,
!     physics%mixing_velocity_scale, physics%no_am_transport_in_core,
!     physics%overshoot_alpha_envelope, physics%pmm_exponent_a,
!     physics%pmm_exponent_b, physics%pmm_exponent_c,
!     physics%pmm_exponent_d, physics%pmm_exponent_m,
!     physics%pmm_norm_jdot, physics%pmm_solar_omega,
!     physics%pmm_solar_pressure, physics%pmm_solar_turnover_timescale,
!     physics%rotation_active, physics%scale_by_b_field,
!     physics%scale_by_rossby_number,
!     physics%secular_shear_mixing_scale,
!     physics%secular_shear_velocity_scale,
!     physics%settling_timestep_fraction, physics%solar_omega_sat,
!     physics%solberg_hoiland_mixing_scale,
!     physics%solberg_hoiland_velocity_scale,
!     physics%solid_body_mode_flag,
!     physics%use_constant_background_diffusion,
!     physics%use_diffusion_advection_transport,
!     physics%use_extended_composition, physics%use_pmm_wind_law,
!     physics%use_wind_torque, physics%velocity_form_correction,
!     physics%wind_law_constant, physics%wind_law_index,
!     physics%wind_law_name, physics%wind_saturation_omega,
!     solar_luminosity_cgs, solar_mass_cgs, solar_radius_cgs
!   star writes: none
pure function rotation_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      use physical_constants, only: log10_newton_g
      use control_slots, only: ichi_dm_max, ichi_dl_max, ichi_dp_core_max
      use transport_mechanism_layout, only: i_mech_eddington_sweet, i_mech_gsf, &
           i_mech_shear, i_mech_solberg_hoiland
      type(star_info), intent(in) :: s
      type(rotation_settings_t) :: settings
      type(wind_law_inputs) :: wind_controls
      type(wind_law_coefficients) :: wind

      wind_controls = wind_law_inputs(use_pmm_wind_law = s%physics%use_pmm_wind_law, &
           law_name = s%physics%wind_law_name, law_index = s%physics%wind_law_index, &
           law_constant = s%physics%wind_law_constant, exponent_a = s%physics%pmm_exponent_a, &
           exponent_b = s%physics%pmm_exponent_b, exponent_c = s%physics%pmm_exponent_c, &
           exponent_d = s%physics%pmm_exponent_d, exponent_m = s%physics%pmm_exponent_m, &
           solar_torque = s%physics%pmm_norm_jdot, solar_omega = s%physics%pmm_solar_omega, &
           saturation_omega = s%physics%wind_saturation_omega, &
           saturation_omega_solar = s%physics%solar_omega_sat, &
           scale_by_b_field = s%physics%scale_by_b_field)
      if (s%physics%rotation_active) then
         wind = derive_wind_law(wind_controls)
      else
         wind = wind_law_coefficients(saturation_omega = s%physics%wind_saturation_omega, &
              scale_by_b_field = s%physics%scale_by_b_field)
      end if
      settings%instability_transport_active = s%physics%instability_transport_active
      settings%diffusion_timestep_factor = s%physics%diffusion_timestep_factor
      settings%max_domega_global = s%physics%max_domega_global
      settings%max_diffusion_iters = s%physics%max_diffusion_iters
      settings%convergence_tolerance = s%physics%convergence_tolerance
      settings%shape_outer_iterations = s%physics%fpft_iter1
      settings%shape_inner_iterations = s%physics%fpft_iter2
      settings%shape_tolerance = s%physics%fpft_average_weight
      settings%cz_rotation_law_exponent = s%physics%cz_rotation_law_exponent
      settings%force_solid_body_rotation = s%physics%force_solid_body_rotation
      settings%solid_body_mode_flag = s%physics%solid_body_mode_flag
      if (s%physics%rotation_active) then
         settings%cz_rotation_law_exponent = min(0.0d0, max(-2.0d0, s%physics%cz_rotation_law_exponent))
         if (s%physics%force_solid_body_rotation) settings%solid_body_mode_flag = 1
      end if
      settings%velocity_scale(i_mech_eddington_sweet) = s%physics%es_velocity_scale
      settings%velocity_scale(i_mech_gsf) = s%physics%gsf_velocity_scale
      settings%velocity_scale(i_mech_shear) = s%physics%secular_shear_velocity_scale
      settings%velocity_scale(i_mech_solberg_hoiland) = s%physics%solberg_hoiland_velocity_scale
      settings%critical_reynolds = s%physics%critical_reynolds
      settings%gsf_inhibition_mode = s%physics%gsf_inhibition_mode
      settings%use_zahn_coupling_factor = s%physics%velocity_form_correction
      settings%no_am_transport_in_core = s%physics%no_am_transport_in_core
      settings%mixing_velocity_scale = s%physics%mixing_velocity_scale
      settings%mixing_scale(i_mech_eddington_sweet) = s%physics%es_mixing_scale
      settings%mixing_scale(i_mech_gsf) = s%physics%gsf_mixing_scale
      settings%mixing_scale(i_mech_shear) = s%physics%secular_shear_mixing_scale
      settings%mixing_scale(i_mech_solberg_hoiland) = s%physics%solberg_hoiland_mixing_scale
      settings%use_constant_background_diffusion = s%physics%use_constant_background_diffusion
      settings%constant_background_diffusion_coeff = s%physics%constant_background_diffusion_coeff
      settings%use_diffusion_advection_transport = s%physics%use_diffusion_advection_transport
      settings%difad_velocity_scale = s%physics%difad_velocity_scale
      settings%chi_mass_scale = s%physics%chi_grid_scale(ichi_dm_max)
      settings%chi_luminosity_scale = s%physics%chi_grid_scale(ichi_dl_max)
      settings%chi_pressure_scale = s%physics%chi_grid_scale(ichi_dp_core_max)
      settings%use_wind_torque = s%physics%use_wind_torque
      settings%wind_constant = wind%constant
      settings%wind_law_omega_exponent = wind%omega_exponent
      settings%wind_mass_loss_exponent = wind%mass_loss_exponent
      settings%wind_radius_exponent = wind%radius_exponent
      settings%wind_mass_exponent = wind%mass_exponent
      settings%wind_saturation_omega = wind%saturation_omega
      settings%use_pmm_wind_law = s%physics%use_pmm_wind_law
      settings%wind_luminosity_exponent = wind%luminosity_exponent
      settings%wind_pressure_exponent = wind%pressure_exponent
      settings%wind_turnover_exponent = wind%turnover_exponent
      settings%wind_centrifugal_exponent = wind%centrifugal_exponent
      settings%wind_centrifugal_constant = wind%centrifugal_constant
      settings%pmm_solar_omega = s%physics%pmm_solar_omega
      settings%pmm_solar_pressure = s%physics%pmm_solar_pressure
      settings%pmm_solar_turnover_timescale = s%physics%pmm_solar_turnover_timescale
      settings%scale_by_rossby_number = s%physics%scale_by_rossby_number
      settings%scale_by_b_field = wind%scale_by_b_field
      settings%disk_locking_active = s%physics%disk_locking_active
      settings%disk_locking_age_gyr = s%physics%disk_locking_age_gyr
      settings%disk_omega_rad_s = s%physics%disk_omega_rad_s
      settings%hydrogen_diffusion_floor = s%physics%hydrogen_diffusion_floor
      settings%helium_diffusion_min = s%physics%helium_diffusion_min
      settings%settling_timestep_fraction = s%physics%settling_timestep_fraction
      settings%use_extended_composition = s%physics%use_extended_composition
      settings%intermediate_overshoot = s%physics%intermediate_overshoot_active
      settings%overshoot_alpha_envelope = s%physics%overshoot_alpha_envelope
      settings%solar_luminosity_cgs = s%solar_luminosity_cgs
      settings%log10_solar_luminosity = s%log10_solar_luminosity
      settings%solar_radius_cgs = s%solar_radius_cgs
      settings%solar_mass_cgs = s%solar_mass_cgs
      settings%log10_gravitational_constant = log10_newton_g(s%physics%G_cgs)
end function rotation_settings_from_star

end module rotation_settings_builder

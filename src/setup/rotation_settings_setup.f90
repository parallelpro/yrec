!----------------------------------------------------------------------
! rotation_settings_setup
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
module rotation_settings_setup
      use rotation_configuration, only: rotation_settings_t
      use magnetic_braking, only: wind_law_inputs, wind_law_coefficients, derive_wind_law
      implicit none
      private
      public :: rotation_settings_from_star

contains

!----------------------------------------------------------------------
! rotation_settings_from_star: rotation's settings from the controls.
!   star reads: ctrl%acfpft, ctrl%chi_grid_scale,
!     ctrl%constant_background_diffusion_coeff,
!     ctrl%convergence_tolerance, ctrl%critical_reynolds,
!     ctrl%difad_velocity_scale, ctrl%diffusion_timestep_factor,
!     ctrl%disk_locking_active, ctrl%disk_locking_age_gyr,
!     ctrl%disk_omega_rad_s, ctrl%es_mixing_scale,
!     ctrl%es_velocity_scale, ctrl%force_solid_body_rotation,
!     ctrl%G_cgs, ctrl%gsf_inhibition_mode, ctrl%gsf_mixing_scale,
!     ctrl%gsf_velocity_scale, ctrl%helium_diffusion_min,
!     ctrl%hydrogen_diffusion_floor, ctrl%instability_transport_active,
!     ctrl%itfp1, ctrl%itfp2, ctrl%lovstm, ctrl%lvfc,
!     ctrl%max_diffusion_iters, ctrl%max_domega_global,
!     ctrl%mixing_velocity_scale, ctrl%no_am_transport_in_core,
!     ctrl%overshoot_alpha_envelope, ctrl%pmm_exponent_a,
!     ctrl%pmm_exponent_b, ctrl%pmm_exponent_c, ctrl%pmm_exponent_d,
!     ctrl%pmm_exponent_m, ctrl%pmm_norm_jdot, ctrl%pmm_solar_omega,
!     ctrl%pmm_solar_pressure, ctrl%pmm_solar_turnover_timescale,
!     ctrl%rotation_active, ctrl%scale_by_b_field,
!     ctrl%scale_by_rossby_number, ctrl%secular_shear_mixing_scale,
!     ctrl%secular_shear_velocity_scale,
!     ctrl%settling_timestep_fraction, ctrl%solar_omega_sat,
!     ctrl%solberg_hoiland_mixing_scale,
!     ctrl%solberg_hoiland_velocity_scale, ctrl%solid_body_mode_flag,
!     ctrl%use_constant_background_diffusion,
!     ctrl%use_diffusion_advection_transport,
!     ctrl%use_extended_composition, ctrl%use_pmm_wind_law,
!     ctrl%use_wind_torque, ctrl%walpcz, ctrl%wind_law_constant,
!     ctrl%wind_law_index, ctrl%wind_law_name,
!     ctrl%wind_saturation_omega, log10_solar_luminosity,
!     solar_luminosity_cgs, solar_mass_cgs, solar_radius_cgs
!   star writes: none
function rotation_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      use physical_constants, only: log10_newton_g
      use control_slots, only: ichi_dm_max, ichi_dl_max, ichi_dp_core_max
      use transport_mechanism_layout, only: i_mech_eddington_sweet, i_mech_gsf, &
           i_mech_shear, i_mech_solberg_hoiland
      type(star_info), intent(in) :: s
      type(rotation_settings_t) :: settings
      type(wind_law_inputs) :: wind_controls
      type(wind_law_coefficients) :: wind

      wind_controls = wind_law_inputs(use_pmm_wind_law = s%ctrl%use_pmm_wind_law, &
           law_name = s%ctrl%wind_law_name, law_index = s%ctrl%wind_law_index, &
           law_constant = s%ctrl%wind_law_constant, exponent_a = s%ctrl%pmm_exponent_a, &
           exponent_b = s%ctrl%pmm_exponent_b, exponent_c = s%ctrl%pmm_exponent_c, &
           exponent_d = s%ctrl%pmm_exponent_d, exponent_m = s%ctrl%pmm_exponent_m, &
           solar_torque = s%ctrl%pmm_norm_jdot, solar_omega = s%ctrl%pmm_solar_omega, &
           saturation_omega = s%ctrl%wind_saturation_omega, &
           saturation_omega_solar = s%ctrl%solar_omega_sat, &
           scale_by_b_field = s%ctrl%scale_by_b_field)
      if (s%ctrl%rotation_active) then
         wind = derive_wind_law(wind_controls)
      else
         wind = wind_law_coefficients(saturation_omega = s%ctrl%wind_saturation_omega, &
              scale_by_b_field = s%ctrl%scale_by_b_field)
      end if
      settings%instability_transport_active = s%ctrl%instability_transport_active
      settings%diffusion_timestep_factor = s%ctrl%diffusion_timestep_factor
      settings%max_domega_global = s%ctrl%max_domega_global
      settings%max_diffusion_iters = s%ctrl%max_diffusion_iters
      settings%convergence_tolerance = s%ctrl%convergence_tolerance
      settings%shape_outer_iterations = s%ctrl%itfp1
      settings%shape_inner_iterations = s%ctrl%itfp2
      settings%shape_tolerance = s%ctrl%acfpft
      settings%cz_rotation_law_exponent = s%ctrl%walpcz
      settings%force_solid_body_rotation = s%ctrl%force_solid_body_rotation
      settings%solid_body_mode_flag = s%ctrl%solid_body_mode_flag
      if (s%ctrl%rotation_active) then
         settings%cz_rotation_law_exponent = min(0.0d0, max(-2.0d0, s%ctrl%walpcz))
         if (s%ctrl%force_solid_body_rotation) settings%solid_body_mode_flag = 1
      end if
      settings%velocity_scale(i_mech_eddington_sweet) = s%ctrl%es_velocity_scale
      settings%velocity_scale(i_mech_gsf) = s%ctrl%gsf_velocity_scale
      settings%velocity_scale(i_mech_shear) = s%ctrl%secular_shear_velocity_scale
      settings%velocity_scale(i_mech_solberg_hoiland) = s%ctrl%solberg_hoiland_velocity_scale
      settings%critical_reynolds = s%ctrl%critical_reynolds
      settings%gsf_inhibition_mode = s%ctrl%gsf_inhibition_mode
      settings%use_zahn_coupling_factor = s%ctrl%lvfc
      settings%no_am_transport_in_core = s%ctrl%no_am_transport_in_core
      settings%mixing_velocity_scale = s%ctrl%mixing_velocity_scale
      settings%mixing_scale(i_mech_eddington_sweet) = s%ctrl%es_mixing_scale
      settings%mixing_scale(i_mech_gsf) = s%ctrl%gsf_mixing_scale
      settings%mixing_scale(i_mech_shear) = s%ctrl%secular_shear_mixing_scale
      settings%mixing_scale(i_mech_solberg_hoiland) = s%ctrl%solberg_hoiland_mixing_scale
      settings%use_constant_background_diffusion = s%ctrl%use_constant_background_diffusion
      settings%constant_background_diffusion_coeff = s%ctrl%constant_background_diffusion_coeff
      settings%use_diffusion_advection_transport = s%ctrl%use_diffusion_advection_transport
      settings%difad_velocity_scale = s%ctrl%difad_velocity_scale
      settings%chi_mass_scale = s%ctrl%chi_grid_scale(ichi_dm_max)
      settings%chi_luminosity_scale = s%ctrl%chi_grid_scale(ichi_dl_max)
      settings%chi_pressure_scale = s%ctrl%chi_grid_scale(ichi_dp_core_max)
      settings%use_wind_torque = s%ctrl%use_wind_torque
      settings%wind_constant = wind%constant
      settings%wind_law_omega_exponent = wind%omega_exponent
      settings%wind_mass_loss_exponent = wind%mass_loss_exponent
      settings%wind_radius_exponent = wind%radius_exponent
      settings%wind_mass_exponent = wind%mass_exponent
      settings%wind_saturation_omega = wind%saturation_omega
      settings%use_pmm_wind_law = s%ctrl%use_pmm_wind_law
      settings%wind_luminosity_exponent = wind%luminosity_exponent
      settings%wind_pressure_exponent = wind%pressure_exponent
      settings%wind_turnover_exponent = wind%turnover_exponent
      settings%wind_centrifugal_exponent = wind%centrifugal_exponent
      settings%wind_centrifugal_constant = wind%centrifugal_constant
      settings%pmm_solar_omega = s%ctrl%pmm_solar_omega
      settings%pmm_solar_pressure = s%ctrl%pmm_solar_pressure
      settings%pmm_solar_turnover_timescale = s%ctrl%pmm_solar_turnover_timescale
      settings%scale_by_rossby_number = s%ctrl%scale_by_rossby_number
      settings%scale_by_b_field = wind%scale_by_b_field
      settings%disk_locking_active = s%ctrl%disk_locking_active
      settings%disk_locking_age_gyr = s%ctrl%disk_locking_age_gyr
      settings%disk_omega_rad_s = s%ctrl%disk_omega_rad_s
      settings%hydrogen_diffusion_floor = s%ctrl%hydrogen_diffusion_floor
      settings%helium_diffusion_min = s%ctrl%helium_diffusion_min
      settings%settling_timestep_fraction = s%ctrl%settling_timestep_fraction
      settings%use_extended_composition = s%ctrl%use_extended_composition
      settings%intermediate_overshoot = s%ctrl%lovstm
      settings%overshoot_alpha_envelope = s%ctrl%overshoot_alpha_envelope
      settings%solar_luminosity_cgs = s%solar_luminosity_cgs
      settings%log10_solar_luminosity = s%log10_solar_luminosity
      settings%solar_radius_cgs = s%solar_radius_cgs
      settings%solar_mass_cgs = s%solar_mass_cgs
      settings%log10_gravitational_constant = log10_newton_g(s%ctrl%G_cgs)
end function rotation_settings_from_star

end module rotation_settings_setup

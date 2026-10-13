!----------------------------------------------------------------------
! rotation_settings_setup
!----------------------------------------------------------------------
! The rotation domain's settings record, built from the run's controls
! by the star layer; run_yrec installs it (set_rotation_settings) at
! the start of every run card, before the starting model is read.
module rotation_settings_setup
      use rotation_configuration, only: rotation_settings_t
      implicit none
      private
      public :: rotation_settings_from_star

contains

!----------------------------------------------------------------------
! rotation_settings_from_star: rotation's settings from the controls.
!   star reads: ctrl%c_2, ctrl%chi_grid_scale,
!     ctrl%constant_background_diffusion_coeff, ctrl%constfactor,
!     ctrl%convergence_tolerance, ctrl%critical_reynolds,
!     ctrl%difad_velocity_scale, ctrl%es_mixing_scale,
!     ctrl%es_velocity_scale, ctrl%excen, ctrl%exl, ctrl%exm, ctrl%exmd,
!     ctrl%expr, ctrl%exr, ctrl%extau, ctrl%force_solid_body_rotation,
!     ctrl%gsf_inhibition_mode, ctrl%gsf_mixing_scale,
!     ctrl%gsf_velocity_scale, ctrl%helium_diffusion_min,
!     ctrl%hydrogen_diffusion_floor, ctrl%itfp1, ctrl%itfp2, ctrl%lvfc,
!     ctrl%max_diffusion_iters, ctrl%mixing_velocity_scale,
!     ctrl%no_am_transport_in_core, ctrl%overshoot_alpha_envelope,
!     ctrl%pmm_solar_omega, ctrl%pmm_solar_pressure,
!     ctrl%pmm_solar_turnover_timescale, ctrl%scale_by_b_field,
!     ctrl%scale_by_rossby_number, ctrl%secular_shear_mixing_scale,
!     ctrl%secular_shear_velocity_scale,
!     ctrl%settling_timestep_fraction, ctrl%solid_body_mode_flag,
!     ctrl%use_constant_background_diffusion,
!     ctrl%use_diffusion_advection_transport, ctrl%use_pmm_wind_law,
!     ctrl%walpcz, ctrl%wind_law_omega_exponent, job%acfpft,
!     job%diffusion_timestep_factor, job%disk_locking_active,
!     job%disk_locking_age_gyr, job%disk_omega_rad_s,
!     job%instability_transport_active, job%lovstm,
!     job%max_domega_global, job%use_extended_composition,
!     job%use_wind_torque, job%wind_saturation_omega,
!     log10_solar_luminosity, solar_luminosity_cgs, solar_mass_cgs,
!     solar_radius_cgs
!   star writes: none
function rotation_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(rotation_settings_t) :: settings

      settings%instability_transport_active = s%job%instability_transport_active
      settings%diffusion_timestep_factor = s%job%diffusion_timestep_factor
      settings%max_domega_global = s%job%max_domega_global
      settings%max_diffusion_iters = s%ctrl%max_diffusion_iters
      settings%convergence_tolerance = s%ctrl%convergence_tolerance
      settings%shape_outer_iterations = s%ctrl%itfp1
      settings%shape_inner_iterations = s%ctrl%itfp2
      settings%shape_tolerance = s%job%acfpft
      settings%cz_rotation_law_exponent = s%ctrl%walpcz
      settings%force_solid_body_rotation = s%ctrl%force_solid_body_rotation
      settings%solid_body_mode_flag = s%ctrl%solid_body_mode_flag
      settings%es_velocity_scale = s%ctrl%es_velocity_scale
      settings%gsf_velocity_scale = s%ctrl%gsf_velocity_scale
      settings%secular_shear_velocity_scale = s%ctrl%secular_shear_velocity_scale
      settings%critical_reynolds = s%ctrl%critical_reynolds
      settings%gsf_inhibition_mode = s%ctrl%gsf_inhibition_mode
      settings%use_zahn_coupling_factor = s%ctrl%lvfc
      settings%no_am_transport_in_core = s%ctrl%no_am_transport_in_core
      settings%mixing_velocity_scale = s%ctrl%mixing_velocity_scale
      settings%es_mixing_scale = s%ctrl%es_mixing_scale
      settings%gsf_mixing_scale = s%ctrl%gsf_mixing_scale
      settings%secular_shear_mixing_scale = s%ctrl%secular_shear_mixing_scale
      settings%use_constant_background_diffusion = s%ctrl%use_constant_background_diffusion
      settings%constant_background_diffusion_coeff = s%ctrl%constant_background_diffusion_coeff
      settings%use_diffusion_advection_transport = s%ctrl%use_diffusion_advection_transport
      settings%difad_velocity_scale = s%ctrl%difad_velocity_scale
      settings%chi_grid_scale = s%ctrl%chi_grid_scale
      settings%use_wind_torque = s%job%use_wind_torque
      settings%wind_constant = s%ctrl%constfactor
      settings%wind_law_omega_exponent = s%ctrl%wind_law_omega_exponent
      settings%wind_mass_loss_exponent = s%ctrl%exmd
      settings%wind_radius_exponent = s%ctrl%exr
      settings%wind_mass_exponent = s%ctrl%exm
      settings%wind_saturation_omega = s%job%wind_saturation_omega
      settings%use_pmm_wind_law = s%ctrl%use_pmm_wind_law
      settings%wind_luminosity_exponent = s%ctrl%exl
      settings%wind_pressure_exponent = s%ctrl%expr
      settings%wind_turnover_exponent = s%ctrl%extau
      settings%wind_centrifugal_exponent = s%ctrl%excen
      settings%wind_centrifugal_constant = s%ctrl%c_2
      settings%pmm_solar_omega = s%ctrl%pmm_solar_omega
      settings%pmm_solar_pressure = s%ctrl%pmm_solar_pressure
      settings%pmm_solar_turnover_timescale = s%ctrl%pmm_solar_turnover_timescale
      settings%scale_by_rossby_number = s%ctrl%scale_by_rossby_number
      settings%scale_by_b_field = s%ctrl%scale_by_b_field
      settings%disk_locking_active = s%job%disk_locking_active
      settings%disk_locking_age_gyr = s%job%disk_locking_age_gyr
      settings%disk_omega_rad_s = s%job%disk_omega_rad_s
      settings%hydrogen_diffusion_floor = s%ctrl%hydrogen_diffusion_floor
      settings%helium_diffusion_min = s%ctrl%helium_diffusion_min
      settings%settling_timestep_fraction = s%ctrl%settling_timestep_fraction
      settings%use_extended_composition = s%job%use_extended_composition
      settings%intermediate_overshoot = s%job%lovstm
      settings%overshoot_alpha_envelope = s%ctrl%overshoot_alpha_envelope
      settings%solar_luminosity_cgs = s%solar_luminosity_cgs
      settings%log10_solar_luminosity = s%log10_solar_luminosity
      settings%solar_radius_cgs = s%solar_radius_cgs
      settings%solar_mass_cgs = s%solar_mass_cgs
end function rotation_settings_from_star

end module rotation_settings_setup

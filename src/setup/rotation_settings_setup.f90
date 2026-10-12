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
!     job%use_wind_torque, job%wind_saturation_omega
!   star writes: none
function rotation_settings_from_star() result(s)
      use star_info_lib, only: star
      type(rotation_settings_t) :: s

      s%instability_transport_active = star%job%instability_transport_active
      s%diffusion_timestep_factor = star%job%diffusion_timestep_factor
      s%max_domega_global = star%job%max_domega_global
      s%max_diffusion_iters = star%ctrl%max_diffusion_iters
      s%convergence_tolerance = star%ctrl%convergence_tolerance
      s%shape_outer_iterations = star%ctrl%itfp1
      s%shape_inner_iterations = star%ctrl%itfp2
      s%shape_tolerance = star%job%acfpft
      s%cz_rotation_law_exponent = star%ctrl%walpcz
      s%force_solid_body_rotation = star%ctrl%force_solid_body_rotation
      s%solid_body_mode_flag = star%ctrl%solid_body_mode_flag
      s%es_velocity_scale = star%ctrl%es_velocity_scale
      s%gsf_velocity_scale = star%ctrl%gsf_velocity_scale
      s%secular_shear_velocity_scale = star%ctrl%secular_shear_velocity_scale
      s%critical_reynolds = star%ctrl%critical_reynolds
      s%gsf_inhibition_mode = star%ctrl%gsf_inhibition_mode
      s%use_zahn_coupling_factor = star%ctrl%lvfc
      s%no_am_transport_in_core = star%ctrl%no_am_transport_in_core
      s%mixing_velocity_scale = star%ctrl%mixing_velocity_scale
      s%es_mixing_scale = star%ctrl%es_mixing_scale
      s%gsf_mixing_scale = star%ctrl%gsf_mixing_scale
      s%secular_shear_mixing_scale = star%ctrl%secular_shear_mixing_scale
      s%use_constant_background_diffusion = star%ctrl%use_constant_background_diffusion
      s%constant_background_diffusion_coeff = star%ctrl%constant_background_diffusion_coeff
      s%use_diffusion_advection_transport = star%ctrl%use_diffusion_advection_transport
      s%difad_velocity_scale = star%ctrl%difad_velocity_scale
      s%chi_grid_scale = star%ctrl%chi_grid_scale
      s%use_wind_torque = star%job%use_wind_torque
      s%wind_constant = star%ctrl%constfactor
      s%wind_law_omega_exponent = star%ctrl%wind_law_omega_exponent
      s%wind_mass_loss_exponent = star%ctrl%exmd
      s%wind_radius_exponent = star%ctrl%exr
      s%wind_mass_exponent = star%ctrl%exm
      s%wind_saturation_omega = star%job%wind_saturation_omega
      s%use_pmm_wind_law = star%ctrl%use_pmm_wind_law
      s%wind_luminosity_exponent = star%ctrl%exl
      s%wind_pressure_exponent = star%ctrl%expr
      s%wind_turnover_exponent = star%ctrl%extau
      s%wind_centrifugal_exponent = star%ctrl%excen
      s%wind_centrifugal_constant = star%ctrl%c_2
      s%pmm_solar_omega = star%ctrl%pmm_solar_omega
      s%pmm_solar_pressure = star%ctrl%pmm_solar_pressure
      s%pmm_solar_turnover_timescale = star%ctrl%pmm_solar_turnover_timescale
      s%scale_by_rossby_number = star%ctrl%scale_by_rossby_number
      s%scale_by_b_field = star%ctrl%scale_by_b_field
      s%disk_locking_active = star%job%disk_locking_active
      s%disk_locking_age_gyr = star%job%disk_locking_age_gyr
      s%disk_omega_rad_s = star%job%disk_omega_rad_s
      s%hydrogen_diffusion_floor = star%ctrl%hydrogen_diffusion_floor
      s%helium_diffusion_min = star%ctrl%helium_diffusion_min
      s%settling_timestep_fraction = star%ctrl%settling_timestep_fraction
      s%use_extended_composition = star%job%use_extended_composition
      s%intermediate_overshoot = star%job%lovstm
      s%overshoot_alpha_envelope = star%ctrl%overshoot_alpha_envelope
end function rotation_settings_from_star

end module rotation_settings_setup

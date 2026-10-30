!----------------------------------------------------------------------
! microphysics_settings
!----------------------------------------------------------------------
! The settings records of the microphysics domains, built from the
! run's controls by the star layer (2026-10, Phase 3); each domain reads
! no other configuration. setup/settings_install.f90 installs them:
!     kap_settings_from_star          kap (kap_init, per run)
!     kap_table_paths_from_star       kap's table paths (kap_init)
!     eos_settings_from_star          eos (eos_init, per run)
!     eos_table_paths_from_star       eos's table paths (eos_init)
!     atm_settings_from_star          atm (atm_init, per run)
!     atm_table_paths_from_star       atm's table paths (atm_init)
!     convection_settings_from_star   convection (set_convection_settings,
!                                     per run card)
!     diffusion_settings_from_star    element diffusion
!                                     (set_diffusion_settings, per run card)
module microphysics_settings
      use kap_lib, only: kap_settings_t, kap_table_paths_t
      use eos_lib, only: eos_settings_t, eos_table_paths_t
      use atm_lib, only: atm_settings_t, atm_table_paths_t
      use convection_lib, only: convection_settings_t
      use element_diffusion_lib, only: diffusion_settings_t
      implicit none
      private
      public :: kap_settings_from_star, kap_table_paths_from_star, &
           eos_settings_from_star, eos_table_paths_from_star, &
           atm_settings_from_star, atm_table_paths_from_star, &
           convection_settings_from_star, diffusion_settings_from_star

contains

!----------------------------------------------------------------------
! kap_settings_from_star: kap's settings from the controls.
!   star reads: ctrl%alex_table_z1, ctrl%kurucz_table_z1,
!     ctrl%kurucz_table_z2, ctrl%laol_extrapolation_tolerance,
!     ctrl%laol_table_z1, ctrl%laol_table_z2,
!     ctrl%metal_fraction_match_tolerance,
!     ctrl%molecular_opacity_logt_max, ctrl%molecular_opacity_logt_min,
!     ctrl%opal95_single_table_z, ctrl%opal_table_z1,
!     ctrl%opal_table_z2, ctrl%use_alex06_tables,
!     ctrl%use_alex95_tables, ctrl%use_conductive_opacity,
!     ctrl%use_kurucz90_tables, ctrl%use_laol89_tables,
!     ctrl%use_opal92_tables, ctrl%use_opal95_tables,
!     ctrl%use_pure_z_table, use_two_z_tables
!   star writes: none
function kap_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(kap_settings_t) :: settings

      settings%use_opal95_tables = s%ctrl%use_opal95_tables
      settings%use_opal92_tables = s%ctrl%use_opal92_tables
      settings%use_laol89_tables = s%ctrl%use_laol89_tables
      settings%use_pure_z_table = s%ctrl%use_pure_z_table
      settings%use_alex06_tables = s%ctrl%use_alex06_tables
      settings%use_alex95_tables = s%ctrl%use_alex95_tables
      settings%use_kurucz90_tables = s%ctrl%use_kurucz90_tables
      settings%use_conductive_opacity = s%ctrl%use_conductive_opacity
      settings%use_two_z_tables = s%use_two_z_tables
      settings%opal95_single_table_z = s%ctrl%opal95_single_table_z
      settings%opal_table_z1 = s%ctrl%opal_table_z1
      settings%opal_table_z2 = s%ctrl%opal_table_z2
      settings%laol_table_z1 = s%ctrl%laol_table_z1
      settings%laol_table_z2 = s%ctrl%laol_table_z2
      settings%kurucz_table_z1 = s%ctrl%kurucz_table_z1
      settings%kurucz_table_z2 = s%ctrl%kurucz_table_z2
      settings%alex_table_z1 = s%ctrl%alex_table_z1
      settings%molecular_opacity_logt_min = s%ctrl%molecular_opacity_logt_min
      settings%molecular_opacity_logt_max = s%ctrl%molecular_opacity_logt_max
      settings%metal_fraction_match_tolerance = s%ctrl%metal_fraction_match_tolerance
      settings%laol_extrapolation_tolerance = s%ctrl%laol_extrapolation_tolerance
end function kap_settings_from_star

!----------------------------------------------------------------------
! kap_table_paths_from_star: the tables kap_init reads, from the
! controls.
!   star reads: ctrl%alex06_table_path, ctrl%alex95_table_paths,
!     ctrl%conductive_table_path, ctrl%kurucz_table2_path,
!     ctrl%kurucz_table_path, ctrl%laol_table2_path,
!     ctrl%laol_table_path, ctrl%opal92_table2_path,
!     ctrl%opal92_table_path, ctrl%opal95_table_path,
!     ctrl%pure_z_table_path
!   star writes: none
function kap_table_paths_from_star(s) result(paths)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(kap_table_paths_t) :: paths

      paths%opal95 = s%ctrl%opal95_table_path
      paths%opal92 = s%ctrl%opal92_table_path
      paths%opal92_z2 = s%ctrl%opal92_table2_path
      paths%laol = s%ctrl%laol_table_path
      paths%laol_z2 = s%ctrl%laol_table2_path
      paths%pure_z = s%ctrl%pure_z_table_path
      paths%alex06 = s%ctrl%alex06_table_path
      paths%alex94 = s%ctrl%alex95_table_paths
      paths%kurucz = s%ctrl%kurucz_table_path
      paths%kurucz_z2 = s%ctrl%kurucz_table2_path
      paths%conductive = s%ctrl%conductive_table_path
end function kap_table_paths_from_star

!----------------------------------------------------------------------
! eos_settings_from_star: eos's settings from the controls (the OPAL
! vintage is derive_controls' opal_eos_vintage).
!   star reads: ctrl%debye_huckel_eta_max, ctrl%debye_huckel_eta_min,
!     ctrl%opal_eos_vintage, ctrl%saha_log10t_cutoff,
!     ctrl%use_debye_huckel_correction, ctrl%use_diffusion_z,
!     ctrl%use_mhd_eos, ctrl%use_numerical_derivatives, ctrl%use_scv_eos
!   star writes: none
function eos_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(eos_settings_t) :: settings

      settings%use_mhd_eos = s%ctrl%use_mhd_eos
! the OPAL vintage derive_controls chose (the newest selected)
      settings%opal_vintage = s%ctrl%opal_eos_vintage
      settings%use_scv_eos = s%ctrl%use_scv_eos
      settings%use_numerical_derivatives = s%ctrl%use_numerical_derivatives
      settings%use_diffusion_z = s%ctrl%use_diffusion_z
      settings%use_debye_huckel_correction = s%ctrl%use_debye_huckel_correction
      settings%debye_huckel_eta_min = s%ctrl%debye_huckel_eta_min
      settings%debye_huckel_eta_max = s%ctrl%debye_huckel_eta_max
      settings%saha_log10t_cutoff = s%ctrl%saha_log10t_cutoff
end function eos_settings_from_star

!----------------------------------------------------------------------
! eos_table_paths_from_star: the tables eos_init reads, from the
! controls.
!   star reads: ctrl%centre1_table_path, ctrl%centre2_table_path,
!     ctrl%centre3_table_path, ctrl%centre4_table_path,
!     ctrl%centre5_table_path, ctrl%fermi_table_path,
!     ctrl%opal_eos_table_path, ctrl%scv_h_table_path,
!     ctrl%scv_he_table_path, ctrl%scv_z_table_path,
!     ctrl%zams_a_table_path, ctrl%zams_b_table_path,
!     ctrl%zams_c_table_path
!   star writes: none
function eos_table_paths_from_star(s) result(paths)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(eos_table_paths_t) :: paths

      paths%fermi = s%ctrl%fermi_table_path
      paths%scv_h = s%ctrl%scv_h_table_path
      paths%scv_he = s%ctrl%scv_he_table_path
      paths%scv_z = s%ctrl%scv_z_table_path
      paths%mhd_zams = [s%ctrl%zams_a_table_path, s%ctrl%zams_b_table_path, &
           s%ctrl%zams_c_table_path]
      paths%mhd_centre = [s%ctrl%centre1_table_path, s%ctrl%centre2_table_path, &
           s%ctrl%centre3_table_path, s%ctrl%centre4_table_path, s%ctrl%centre5_table_path]
      paths%opal = s%ctrl%opal_eos_table_path
end function eos_table_paths_from_star

!----------------------------------------------------------------------
! atm_settings_from_star: atm's settings from the controls.
!   star reads: ctrl%allard_target_alpha, ctrl%allard_target_feh,
!     ctrl%allard_use_tau100, job%atm_choice
!   star writes: none
function atm_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(atm_settings_t) :: settings

      settings%atm_choice = s%job%atm_choice
      settings%allard_use_tau100 = s%ctrl%allard_use_tau100
      settings%allard_target_feh = s%ctrl%allard_target_feh
      settings%allard_target_alpha = s%ctrl%allard_target_alpha
end function atm_settings_from_star

!----------------------------------------------------------------------
! atm_table_paths_from_star: the tables atm_init reads, from the
! controls.
!   star reads: ctrl%allard_table_path, ctrl%atm_table_path
!   star writes: none
function atm_table_paths_from_star(s) result(paths)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(atm_table_paths_t) :: paths

      paths%surface_pressure = s%ctrl%atm_table_path
      paths%allard = s%ctrl%allard_table_path
end function atm_table_paths_from_star

!----------------------------------------------------------------------
! convection_settings_from_star: convection's settings from the
! controls and the run card's mixing length; log10 G is derived from
! G_cgs (log10_newton_g).
!   star reads: ctrl%adiabatic_overshoot_active,
!     ctrl%core_overshoot_active, ctrl%envelope_overshoot_active,
!     ctrl%G_cgs, ctrl%gradient_logt_cutoff,
!     ctrl%instability_transport_active,
!     ctrl%intermediate_overshoot_active, ctrl%overshoot_alpha_core,
!     ctrl%overshoot_alpha_envelope, ctrl%overshoot_alpha_min,
!     ctrl%overshoot_beta_core, ctrl%overshoot_max_extent_flag,
!     ctrl%rotation_active, ctrl%spot_depth_varies,
!     ctrl%spot_filling_factor, ctrl%spot_temp_contrast,
!     ctrl%use_extended_composition, envelope_metals_no_cno,
!     log10_solar_luminosity, mixing_length_alpha, solar_luminosity_cgs
!   star writes: none
function convection_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      use physical_constants, only: log10_newton_g
      type(star_info), intent(in) :: s
      type(convection_settings_t) :: settings

      settings%mixing_length_alpha = s%mixing_length_alpha
      settings%log10_t_cut = s%ctrl%gradient_logt_cutoff
      settings%adiabatic_overshoot = s%ctrl%adiabatic_overshoot_active
      settings%solar_luminosity_cgs = s%solar_luminosity_cgs
      settings%log10_solar_luminosity = s%log10_solar_luminosity
      settings%log10_gravitational_constant = log10_newton_g(s%ctrl%G_cgs)
      settings%spot_filling_factor = s%ctrl%spot_filling_factor
      settings%spot_temp_contrast = s%ctrl%spot_temp_contrast
      settings%spot_depth_varies = s%ctrl%spot_depth_varies
      settings%core_overshoot = s%ctrl%core_overshoot_active
      settings%envelope_overshoot = s%ctrl%envelope_overshoot_active
      settings%intermediate_overshoot = s%ctrl%intermediate_overshoot_active
      settings%overshoot_alpha_core = s%ctrl%overshoot_alpha_core
      settings%overshoot_alpha_envelope = s%ctrl%overshoot_alpha_envelope
      settings%overshoot_alpha_intermediate = s%ctrl%overshoot_alpha_min
      settings%limit_core_overshoot = s%ctrl%overshoot_max_extent_flag
      settings%core_overshoot_limit = s%ctrl%overshoot_beta_core
      settings%rotation_active = s%ctrl%rotation_active
      settings%instability_transport_active = s%ctrl%instability_transport_active
      settings%use_extended_composition = s%ctrl%use_extended_composition
      settings%envelope_metals_no_cno = s%envelope_metals_no_cno
end function convection_settings_from_star

!----------------------------------------------------------------------
! diffusion_settings_from_star: element diffusion's settings from the
! controls (the Monte-Carlo run's fgrz included).
!   star reads: ctrl%constant_mixing_coeff,
!     ctrl%constant_settling_reduction, ctrl%coulomb_log_choice,
!     ctrl%diffuse_lithium, ctrl%helium_diffusion_min,
!     ctrl%hydrogen_diffusion_floor, ctrl%settling_num_iterations,
!     ctrl%settling_tolerance, ctrl%settling_y_factor,
!     ctrl%use_diffusion_z, ctrl%use_new_diffusion_routines,
!     ctrl%use_thoul_diffusion, ctrl%use_thoul_fit,
!     job%settling_z_factor
!   star writes: none
function diffusion_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(diffusion_settings_t) :: settings

      settings%use_new_diffusion_routines = s%ctrl%use_new_diffusion_routines
      settings%use_diffusion_z = s%ctrl%use_diffusion_z
      settings%diffuse_lithium = s%ctrl%diffuse_lithium
      settings%fgry = s%ctrl%settling_y_factor
      settings%fgrz = s%job%settling_z_factor
      settings%hydrogen_diffusion_floor = s%ctrl%hydrogen_diffusion_floor
      settings%helium_diffusion_min = s%ctrl%helium_diffusion_min
      settings%settling_num_iterations = s%ctrl%settling_num_iterations
      settings%settling_tolerance = s%ctrl%settling_tolerance
      settings%constant_mixing_coeff = s%ctrl%constant_mixing_coeff
      settings%coulomb_log_choice = s%ctrl%coulomb_log_choice
      settings%use_thoul_fit = s%ctrl%use_thoul_fit
      settings%use_thoul_diffusion = s%ctrl%use_thoul_diffusion
      settings%constant_settling_reduction = s%ctrl%constant_settling_reduction
end function diffusion_settings_from_star

end module microphysics_settings

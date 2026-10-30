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
!   star reads: control%alex_table_z1, control%kurucz_table_z1,
!     control%kurucz_table_z2, control%laol_table_z1,
!     control%laol_table_z2, control%opal95_single_table_z,
!     control%opal_table_z1, control%opal_table_z2,
!     control%use_alex06_tables, control%use_alex95_tables,
!     control%use_kurucz90_tables, control%use_laol89_tables,
!     control%use_opal92_tables, control%use_opal95_tables,
!     control%use_pure_z_table, physics%laol_extrapolation_tolerance,
!     physics%metal_fraction_match_tolerance,
!     physics%molecular_opacity_logt_max,
!     physics%molecular_opacity_logt_min,
!     physics%use_conductive_opacity, use_two_z_tables
!   star writes: none
function kap_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(kap_settings_t) :: settings

      settings%use_opal95_tables = s%control%use_opal95_tables
      settings%use_opal92_tables = s%control%use_opal92_tables
      settings%use_laol89_tables = s%control%use_laol89_tables
      settings%use_pure_z_table = s%control%use_pure_z_table
      settings%use_alex06_tables = s%control%use_alex06_tables
      settings%use_alex95_tables = s%control%use_alex95_tables
      settings%use_kurucz90_tables = s%control%use_kurucz90_tables
      settings%use_conductive_opacity = s%physics%use_conductive_opacity
      settings%use_two_z_tables = s%use_two_z_tables
      settings%opal95_single_table_z = s%control%opal95_single_table_z
      settings%opal_table_z1 = s%control%opal_table_z1
      settings%opal_table_z2 = s%control%opal_table_z2
      settings%laol_table_z1 = s%control%laol_table_z1
      settings%laol_table_z2 = s%control%laol_table_z2
      settings%kurucz_table_z1 = s%control%kurucz_table_z1
      settings%kurucz_table_z2 = s%control%kurucz_table_z2
      settings%alex_table_z1 = s%control%alex_table_z1
      settings%molecular_opacity_logt_min = s%physics%molecular_opacity_logt_min
      settings%molecular_opacity_logt_max = s%physics%molecular_opacity_logt_max
      settings%metal_fraction_match_tolerance = s%physics%metal_fraction_match_tolerance
      settings%laol_extrapolation_tolerance = s%physics%laol_extrapolation_tolerance
end function kap_settings_from_star

!----------------------------------------------------------------------
! kap_table_paths_from_star: the tables kap_init reads, from the
! controls.
!   star reads: control%alex06_table_path, control%alex95_table_paths,
!     control%conductive_table_path, control%kurucz_table2_path,
!     control%kurucz_table_path, control%laol_table2_path,
!     control%laol_table_path, control%opal92_table2_path,
!     control%opal92_table_path, control%opal95_table_path,
!     control%pure_z_table_path
!   star writes: none
function kap_table_paths_from_star(s) result(paths)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(kap_table_paths_t) :: paths

      paths%opal95 = s%control%opal95_table_path
      paths%opal92 = s%control%opal92_table_path
      paths%opal92_z2 = s%control%opal92_table2_path
      paths%laol = s%control%laol_table_path
      paths%laol_z2 = s%control%laol_table2_path
      paths%pure_z = s%control%pure_z_table_path
      paths%alex06 = s%control%alex06_table_path
      paths%alex94 = s%control%alex95_table_paths
      paths%kurucz = s%control%kurucz_table_path
      paths%kurucz_z2 = s%control%kurucz_table2_path
      paths%conductive = s%control%conductive_table_path
end function kap_table_paths_from_star

!----------------------------------------------------------------------
! eos_settings_from_star: eos's settings from the controls (the OPAL
! vintage is derive_controls' opal_eos_vintage).
!   star reads: physics%debye_huckel_eta_max,
!     physics%debye_huckel_eta_min, physics%opal_eos_vintage,
!     physics%saha_log10t_cutoff, physics%use_debye_huckel_correction,
!     physics%use_diffusion_z, physics%use_mhd_eos,
!     physics%use_numerical_derivatives, physics%use_scv_eos
!   star writes: none
function eos_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(eos_settings_t) :: settings

      settings%use_mhd_eos = s%physics%use_mhd_eos
! the OPAL vintage derive_controls chose (the newest selected)
      settings%opal_vintage = s%physics%opal_eos_vintage
      settings%use_scv_eos = s%physics%use_scv_eos
      settings%use_numerical_derivatives = s%physics%use_numerical_derivatives
      settings%use_diffusion_z = s%physics%use_diffusion_z
      settings%use_debye_huckel_correction = s%physics%use_debye_huckel_correction
      settings%debye_huckel_eta_min = s%physics%debye_huckel_eta_min
      settings%debye_huckel_eta_max = s%physics%debye_huckel_eta_max
      settings%saha_log10t_cutoff = s%physics%saha_log10t_cutoff
end function eos_settings_from_star

!----------------------------------------------------------------------
! eos_table_paths_from_star: the tables eos_init reads, from the
! controls.
!   star reads: control%centre1_table_path, control%centre2_table_path,
!     control%centre3_table_path, control%centre4_table_path,
!     control%centre5_table_path, control%fermi_table_path,
!     control%scv_h_table_path, control%scv_he_table_path,
!     control%scv_z_table_path, control%zams_a_table_path,
!     control%zams_b_table_path, control%zams_c_table_path,
!     physics%opal_eos_table_path
!   star writes: none
function eos_table_paths_from_star(s) result(paths)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(eos_table_paths_t) :: paths

      paths%fermi = s%control%fermi_table_path
      paths%scv_h = s%control%scv_h_table_path
      paths%scv_he = s%control%scv_he_table_path
      paths%scv_z = s%control%scv_z_table_path
      paths%mhd_zams = [s%control%zams_a_table_path, s%control%zams_b_table_path, &
           s%control%zams_c_table_path]
      paths%mhd_centre = [s%control%centre1_table_path, s%control%centre2_table_path, &
           s%control%centre3_table_path, s%control%centre4_table_path, s%control%centre5_table_path]
      paths%opal = s%physics%opal_eos_table_path
end function eos_table_paths_from_star

!----------------------------------------------------------------------
! atm_settings_from_star: atm's settings from the controls.
!   star reads: atm_choice_in_use, physics%allard_target_alpha,
!     physics%allard_target_feh, physics%allard_use_tau100
!   star writes: none
function atm_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(atm_settings_t) :: settings

      settings%atm_choice = s%atm_choice_in_use
      settings%allard_use_tau100 = s%physics%allard_use_tau100
      settings%allard_target_feh = s%physics%allard_target_feh
      settings%allard_target_alpha = s%physics%allard_target_alpha
end function atm_settings_from_star

!----------------------------------------------------------------------
! atm_table_paths_from_star: the tables atm_init reads, from the
! controls.
!   star reads: control%allard_table_path, control%atm_table_path
!   star writes: none
function atm_table_paths_from_star(s) result(paths)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(atm_table_paths_t) :: paths

      paths%surface_pressure = s%control%atm_table_path
      paths%allard = s%control%allard_table_path
end function atm_table_paths_from_star

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
function convection_settings_from_star(s) result(settings)
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
function diffusion_settings_from_star(s) result(settings)
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

end module microphysics_settings

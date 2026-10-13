!----------------------------------------------------------------------
! microphysics_settings
!----------------------------------------------------------------------
! The settings records of the microphysics domains, built from the
! run's controls by the star layer (2026-10, Phase 3). Each domain
! installs its record at init and reads no other configuration:
!     opacity_settings   kap (kap_init)
!     equation_of_state_settings       eos (eos_init)
!     atmosphere_settings   atm (atm_init)
!     convection_settings_from_star   convection (set_convection_settings,
!                        at the start of every run card)
!     diffusion_settings_from_star   element diffusion
!                        (set_diffusion_settings, at the start of every
!                        run card)
module microphysics_settings
      use kap_lib, only: kap_settings_t
      use eos_lib, only: eos_settings_t
      use atm_lib, only: atm_settings_t
      use convection_lib, only: convection_settings_t
      use element_diffusion_lib, only: diffusion_settings_t
      implicit none
      private
      public :: opacity_settings, equation_of_state_settings, atmosphere_settings, &
           convection_settings_from_star, diffusion_settings_from_star

contains

!----------------------------------------------------------------------
! opacity_settings: kap's settings from the controls.
!   star reads: ctrl%alex_table_z1, ctrl%kurucz_table_z1,
!     ctrl%kurucz_table_z2, ctrl%laol_table_z1, ctrl%laol_table_z2,
!     ctrl%metal_fraction_match_tolerance,
!     ctrl%molecular_opacity_logt_max, ctrl%molecular_opacity_logt_min,
!     ctrl%opal95_single_table_z, ctrl%opal_table_z1,
!     ctrl%opal_table_z2, ctrl%tollaol, ctrl%use_alex06_tables,
!     ctrl%use_alex95_tables, ctrl%use_conductive_opacity,
!     ctrl%use_kurucz90_tables, ctrl%use_laol89_tables,
!     ctrl%use_opal92_tables, ctrl%use_opal95_tables,
!     ctrl%use_pure_z_table, use_two_z_tables
!   star writes: none
function opacity_settings(s) result(settings)
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
      settings%laol_extrapolation_tolerance = s%ctrl%tollaol
end function opacity_settings

!----------------------------------------------------------------------
! equation_of_state_settings: eos's settings from the controls.
!   star reads: ctrl%debye_huckel_eta_max, ctrl%debye_huckel_eta_min,
!     ctrl%saha_log10t_cutoff, ctrl%use_debye_huckel_correction,
!     ctrl%use_mhd_eos, ctrl%use_numerical_derivatives,
!     ctrl%use_opal2001_eos, ctrl%use_opal2006_eos, ctrl%use_opal95_eos,
!     ctrl%use_scv_eos, job%use_diffusion_z
!   star writes: none
function equation_of_state_settings(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(eos_settings_t) :: settings

      settings%use_mhd_eos = s%ctrl%use_mhd_eos
      settings%use_opal95_eos = s%ctrl%use_opal95_eos
      settings%use_opal2001_eos = s%ctrl%use_opal2001_eos
      settings%use_opal2006_eos = s%ctrl%use_opal2006_eos
      settings%use_scv_eos = s%ctrl%use_scv_eos
      settings%use_numerical_derivatives = s%ctrl%use_numerical_derivatives
      settings%use_diffusion_z = s%job%use_diffusion_z
      settings%use_debye_huckel_correction = s%ctrl%use_debye_huckel_correction
      settings%debye_huckel_eta_min = s%ctrl%debye_huckel_eta_min
      settings%debye_huckel_eta_max = s%ctrl%debye_huckel_eta_max
      settings%saha_log10t_cutoff = s%ctrl%saha_log10t_cutoff
end function equation_of_state_settings

!----------------------------------------------------------------------
! atmosphere_settings: atm's settings from the controls.
!   star reads: ctrl%allard_target_alpha, ctrl%allard_target_feh,
!     ctrl%allard_use_tau100, job%atm_choice
!   star writes: none
function atmosphere_settings(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(atm_settings_t) :: settings

      settings%atm_choice = s%job%atm_choice
      settings%allard_use_tau100 = s%ctrl%allard_use_tau100
      settings%allard_target_feh = s%ctrl%allard_target_feh
      settings%allard_target_alpha = s%ctrl%allard_target_alpha
end function atmosphere_settings

!----------------------------------------------------------------------
! convection_settings_from_star: convection's settings from the
! controls and the run card's mixing length.
!   star reads: ctrl%alpham, ctrl%betac, ctrl%ladov, ctrl%lovmax,
!     ctrl%overshoot_alpha_core, ctrl%overshoot_alpha_envelope,
!     ctrl%spot_depth_varies, ctrl%spot_filling_factor,
!     ctrl%spot_temp_contrast, ctrl%tgcut, job%core_overshoot_active,
!     job%envelope_overshoot_active, job%instability_transport_active,
!     job%lovstm, job%rotation_active, job%use_extended_composition,
!     log10_solar_luminosity, mixing_length_alpha, solar_luminosity_cgs,
!     zenvm
!   star writes: none
function convection_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(convection_settings_t) :: settings

      settings%mixing_length_alpha = s%mixing_length_alpha
      settings%log10_t_cut = s%ctrl%tgcut
      settings%adiabatic_overshoot = s%ctrl%ladov
      settings%solar_luminosity_cgs = s%solar_luminosity_cgs
      settings%log10_solar_luminosity = s%log10_solar_luminosity
      settings%spot_filling_factor = s%ctrl%spot_filling_factor
      settings%spot_temp_contrast = s%ctrl%spot_temp_contrast
      settings%spot_depth_varies = s%ctrl%spot_depth_varies
      settings%core_overshoot = s%job%core_overshoot_active
      settings%envelope_overshoot = s%job%envelope_overshoot_active
      settings%intermediate_overshoot = s%job%lovstm
      settings%overshoot_alpha_core = s%ctrl%overshoot_alpha_core
      settings%overshoot_alpha_envelope = s%ctrl%overshoot_alpha_envelope
      settings%overshoot_alpha_intermediate = s%ctrl%alpham
      settings%limit_core_overshoot = s%ctrl%lovmax
      settings%core_overshoot_limit = s%ctrl%betac
      settings%rotation_active = s%job%rotation_active
      settings%instability_transport_active = s%job%instability_transport_active
      settings%use_extended_composition = s%job%use_extended_composition
      settings%envelope_metals_no_cno = s%zenvm
end function convection_settings_from_star

!----------------------------------------------------------------------
! diffusion_settings_from_star: element diffusion's settings from the
! controls (the Monte-Carlo run's fgrz included).
!   star reads: ctrl%constant_mixing_coeff,
!     ctrl%constant_settling_reduction, ctrl%coulomb_log_choice,
!     ctrl%diffuse_lithium, ctrl%fgry, ctrl%helium_diffusion_min,
!     ctrl%hydrogen_diffusion_floor, ctrl%settling_num_iterations,
!     ctrl%settling_tolerance, ctrl%use_new_diffusion_routines,
!     ctrl%use_thoul_diffusion, ctrl%use_thoul_fit, job%fgrz,
!     job%use_diffusion_z
!   star writes: none
function diffusion_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(diffusion_settings_t) :: settings

      settings%use_new_diffusion_routines = s%ctrl%use_new_diffusion_routines
      settings%use_diffusion_z = s%job%use_diffusion_z
      settings%diffuse_lithium = s%ctrl%diffuse_lithium
      settings%fgry = s%ctrl%fgry
      settings%fgrz = s%job%fgrz
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

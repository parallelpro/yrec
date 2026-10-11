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
module microphysics_settings
      use kap_lib, only: kap_settings_t
      use eos_lib, only: eos_settings_t
      use atm_lib, only: atm_settings_t
      use convection_lib, only: convection_settings_t
      implicit none
      private
      public :: opacity_settings, equation_of_state_settings, atmosphere_settings, &
           convection_settings_from_star

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
function opacity_settings() result(s)
      use star_info_lib, only: star
      type(kap_settings_t) :: s

      s%use_opal95_tables = star%ctrl%use_opal95_tables
      s%use_opal92_tables = star%ctrl%use_opal92_tables
      s%use_laol89_tables = star%ctrl%use_laol89_tables
      s%use_pure_z_table = star%ctrl%use_pure_z_table
      s%use_alex06_tables = star%ctrl%use_alex06_tables
      s%use_alex95_tables = star%ctrl%use_alex95_tables
      s%use_kurucz90_tables = star%ctrl%use_kurucz90_tables
      s%use_conductive_opacity = star%ctrl%use_conductive_opacity
      s%use_two_z_tables = star%use_two_z_tables
      s%opal95_single_table_z = star%ctrl%opal95_single_table_z
      s%opal_table_z1 = star%ctrl%opal_table_z1
      s%opal_table_z2 = star%ctrl%opal_table_z2
      s%laol_table_z1 = star%ctrl%laol_table_z1
      s%laol_table_z2 = star%ctrl%laol_table_z2
      s%kurucz_table_z1 = star%ctrl%kurucz_table_z1
      s%kurucz_table_z2 = star%ctrl%kurucz_table_z2
      s%alex_table_z1 = star%ctrl%alex_table_z1
      s%molecular_opacity_logt_min = star%ctrl%molecular_opacity_logt_min
      s%molecular_opacity_logt_max = star%ctrl%molecular_opacity_logt_max
      s%metal_fraction_match_tolerance = star%ctrl%metal_fraction_match_tolerance
      s%laol_extrapolation_tolerance = star%ctrl%tollaol
end function opacity_settings

!----------------------------------------------------------------------
! equation_of_state_settings: eos's settings from the controls.
!   star reads: ctrl%debye_huckel_eta_max, ctrl%debye_huckel_eta_min,
!     ctrl%saha_log10t_cutoff, ctrl%use_debye_huckel_correction,
!     ctrl%use_mhd_eos, ctrl%use_numerical_derivatives,
!     ctrl%use_opal2001_eos, ctrl%use_opal2006_eos, ctrl%use_opal95_eos,
!     ctrl%use_scv_eos, job%use_diffusion_z
!   star writes: none
function equation_of_state_settings() result(s)
      use star_info_lib, only: star
      type(eos_settings_t) :: s

      s%use_mhd_eos = star%ctrl%use_mhd_eos
      s%use_opal95_eos = star%ctrl%use_opal95_eos
      s%use_opal2001_eos = star%ctrl%use_opal2001_eos
      s%use_opal2006_eos = star%ctrl%use_opal2006_eos
      s%use_scv_eos = star%ctrl%use_scv_eos
      s%use_numerical_derivatives = star%ctrl%use_numerical_derivatives
      s%use_diffusion_z = star%job%use_diffusion_z
      s%use_debye_huckel_correction = star%ctrl%use_debye_huckel_correction
      s%debye_huckel_eta_min = star%ctrl%debye_huckel_eta_min
      s%debye_huckel_eta_max = star%ctrl%debye_huckel_eta_max
      s%saha_log10t_cutoff = star%ctrl%saha_log10t_cutoff
end function equation_of_state_settings

!----------------------------------------------------------------------
! atmosphere_settings: atm's settings from the controls.
!   star reads: ctrl%allard_target_alpha, ctrl%allard_target_feh,
!     ctrl%allard_use_tau100, job%atm_choice
!   star writes: none
function atmosphere_settings() result(s)
      use star_info_lib, only: star
      type(atm_settings_t) :: s

      s%atm_choice = star%job%atm_choice
      s%allard_use_tau100 = star%ctrl%allard_use_tau100
      s%allard_target_feh = star%ctrl%allard_target_feh
      s%allard_target_alpha = star%ctrl%allard_target_alpha
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
function convection_settings_from_star() result(s)
      use star_info_lib, only: star
      type(convection_settings_t) :: s

      s%mixing_length_alpha = star%mixing_length_alpha
      s%log10_t_cut = star%ctrl%tgcut
      s%adiabatic_overshoot = star%ctrl%ladov
      s%solar_luminosity_cgs = star%solar_luminosity_cgs
      s%log10_solar_luminosity = star%log10_solar_luminosity
      s%spot_filling_factor = star%ctrl%spot_filling_factor
      s%spot_temp_contrast = star%ctrl%spot_temp_contrast
      s%spot_depth_varies = star%ctrl%spot_depth_varies
      s%core_overshoot = star%job%core_overshoot_active
      s%envelope_overshoot = star%job%envelope_overshoot_active
      s%intermediate_overshoot = star%job%lovstm
      s%overshoot_alpha_core = star%ctrl%overshoot_alpha_core
      s%overshoot_alpha_envelope = star%ctrl%overshoot_alpha_envelope
      s%overshoot_alpha_intermediate = star%ctrl%alpham
      s%limit_core_overshoot = star%ctrl%lovmax
      s%core_overshoot_limit = star%ctrl%betac
      s%rotation_active = star%job%rotation_active
      s%instability_transport_active = star%job%instability_transport_active
      s%use_extended_composition = star%job%use_extended_composition
      s%envelope_metals_no_cno = star%zenvm
end function convection_settings_from_star

end module microphysics_settings

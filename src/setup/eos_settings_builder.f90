!----------------------------------------------------------------------
! eos_settings_builder
!----------------------------------------------------------------------
! eos's settings record and its tables' paths, built from the run's
! controls (setup/README.md). settings_install's install_run_settings
! installs them once per run, with eos_init, which also loads the tables:
!     eos_settings_from_star          eos_settings_t
!     eos_table_paths_from_star       eos_table_paths_t
module eos_settings_builder
      use eos_lib, only: eos_settings_t, eos_table_paths_t
      implicit none
      private
      public :: eos_settings_from_star, eos_table_paths_from_star

contains

!----------------------------------------------------------------------
! eos_settings_from_star: eos's settings from the controls (the OPAL
! vintage is derive_controls' opal_eos_vintage).
!   star reads: physics%debye_huckel_eta_max,
!     physics%debye_huckel_eta_min, physics%opal_eos_vintage,
!     physics%saha_log10t_cutoff, physics%use_debye_huckel_correction,
!     physics%use_diffusion_z, physics%use_mhd_eos,
!     physics%use_numerical_derivatives, physics%use_scv_eos
!   star writes: none
pure function eos_settings_from_star(s) result(settings)
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
pure function eos_table_paths_from_star(s) result(paths)
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

end module eos_settings_builder

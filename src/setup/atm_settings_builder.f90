!----------------------------------------------------------------------
! atm_settings_builder
!----------------------------------------------------------------------
! atm's settings record and its tables' paths, built from the run's
! controls (setup/README.md). settings_install's install_run_settings
! installs them once per run, with atm_init, which also loads the tables:
!     atm_settings_from_star          atm_settings_t
!     atm_table_paths_from_star       atm_table_paths_t
module atm_settings_builder
      use atm_lib, only: atm_settings_t, atm_table_paths_t
      implicit none
      private
      public :: atm_settings_from_star, atm_table_paths_from_star

contains

!----------------------------------------------------------------------
! atm_settings_from_star: atm's settings from the controls; atm_choice
! is the deck's (physics%atm_choice), which selects the tables atm_init
! loads, not the run's working choice atm_choice_in_use, which the
! surface boundary switches.
!   star reads: physics%allard_target_alpha, physics%allard_target_feh,
!     physics%allard_use_tau100, physics%atm_choice
!   star writes: none
pure function atm_settings_from_star(s) result(settings)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(atm_settings_t) :: settings

      settings%atm_choice = s%physics%atm_choice
      settings%allard_use_tau100 = s%physics%allard_use_tau100
      settings%allard_target_feh = s%physics%allard_target_feh
      settings%allard_target_alpha = s%physics%allard_target_alpha
end function atm_settings_from_star

!----------------------------------------------------------------------
! atm_table_paths_from_star: the tables atm_init reads, from the
! controls.
!   star reads: control%allard_table_path, control%atm_table_path
!   star writes: none
pure function atm_table_paths_from_star(s) result(paths)
      use star_info_lib, only: star_info
      type(star_info), intent(in) :: s
      type(atm_table_paths_t) :: paths

      paths%surface_pressure = s%control%atm_table_path
      paths%allard = s%control%allard_table_path
end function atm_table_paths_from_star

end module atm_settings_builder

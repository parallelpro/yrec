!----------------------------------------------------------------------
! settings_install
!----------------------------------------------------------------------
! Where the star layer installs the physics domains' settings records:
! every star-layer call of a domain's installer (set_<domain>_settings,
! or <domain>_init for the domains that also load their tables) is here,
! and tools/check_boundaries.py (INSTALL_FILES) keeps it so. Each
! procedure builds the records with setup/'s builders
! (<domain>_settings_from_star) and installs them, in a fixed order:
!   install_run_settings              per run, after the controls read
!                                     (run_setup's set_up_constants_and_tables):
!                                     kap, eos and atm, with their tables
!   install_monte_carlo_settings      per Monte-Carlo run, after
!                                     apply_monte_carlo_parameters
!                                     (run_yrec): net
!   install_card_settings             per run card, after the card's
!                                     mixing length, before its model read
!                                     (begin_kind_card): convection,
!                                     rotation, element diffusion
!   install_surface_mixture_settings  after the first card's model read
!                                     (read_starting_model's
!                                     update_surface_mixture): convection's
!                                     envelope_metals_no_cno, known only then
! A new domain's record is installed by adding its call to the procedure
! of its lifetime (setup/README.md, "Add a physics domain").
module settings_install
      use star_info_lib, only: star_info
      implicit none
      private
      public :: install_run_settings, install_monte_carlo_settings, &
           install_card_settings, install_surface_mixture_settings

contains

!----------------------------------------------------------------------
! install_run_settings: the opacities (kap_init), the equation of state
! (eos_init: the Fermi-Dirac degenerate-electron table always, and the
! OPAL, MHD and SCVH tables when the settings select them) and the
! tabulated surface boundary conditions (atm_init: the Kurucz (atm_choice
! 3 and 4), Allard (4) and Kurucz/Castelli (5) tables), each with its
! table paths. The table paths are controls, so the controls read must
! have run. Errors return through ierr.
!   star reads: envelope_hydrogen_fraction, envelope_metal_fraction
!   star writes: none
subroutine install_run_settings(s, ierr)
      use kap_lib, only: kap_init
      use eos_lib, only: eos_init
      use atm_lib, only: atm_init
      use microphysics_settings, only: kap_settings_from_star, kap_table_paths_from_star, &
           eos_settings_from_star, eos_table_paths_from_star, atm_settings_from_star, &
           atm_table_paths_from_star
      type(star_info), intent(in) :: s
      integer, intent(out) :: ierr

      ierr = 0
      call kap_init(kap_settings_from_star(s), kap_table_paths_from_star(s), &
           s%envelope_hydrogen_fraction, s%envelope_metal_fraction, ierr)
      if (ierr /= 0) return
      call eos_init(eos_settings_from_star(s), eos_table_paths_from_star(s), ierr)
      if (ierr /= 0) return
      call atm_init(atm_settings_from_star(s), atm_table_paths_from_star(s), ierr)
end subroutine install_run_settings

!----------------------------------------------------------------------
! install_monte_carlo_settings: the network's record, which takes the
! Monte-Carlo run's rate factors; set_net_settings checks the reaction
! table first and returns ierr.
!   star reads: none
!   star writes: none
subroutine install_monte_carlo_settings(s, ierr)
      use net_lib, only: set_net_settings
      use net_settings_setup, only: net_settings_from_star
      type(star_info), intent(in) :: s
      integer, intent(out) :: ierr

      call set_net_settings(net_settings_from_star(s), ierr)
end subroutine install_monte_carlo_settings

!----------------------------------------------------------------------
! install_card_settings: convection's record (the card's mixing
! length), rotation's and element diffusion's.
!   star reads: none
!   star writes: none
subroutine install_card_settings(s)
      use convection_lib, only: set_convection_settings
      use rotation_configuration, only: set_rotation_settings
      use element_diffusion_lib, only: set_diffusion_settings
      use microphysics_settings, only: convection_settings_from_star, &
           diffusion_settings_from_star
      use rotation_settings_setup, only: rotation_settings_from_star
      type(star_info), intent(in) :: s

      call set_convection_settings(convection_settings_from_star(s))
      call set_rotation_settings(rotation_settings_from_star(s))
      call set_diffusion_settings(diffusion_settings_from_star(s))
end subroutine install_card_settings

!----------------------------------------------------------------------
! install_surface_mixture_settings: convection's envelope_metals_no_cno,
! the first card's surface Z without C, N and O, which the starting
! model's surface mixture gives (install_card_settings installed the
! value the star had before: 0 on the first card).
!   star reads: envelope_metals_no_cno
!   star writes: none
subroutine install_surface_mixture_settings(s)
      use convection_lib, only: set_convection_envelope_metals
      type(star_info), intent(in) :: s

      call set_convection_envelope_metals(s%envelope_metals_no_cno)
end subroutine install_surface_mixture_settings

end module settings_install

!----------------------------------------------------------------------
! settings_install
!----------------------------------------------------------------------
! Where the star layer installs the physics domains' configuration:
! every star-layer call of a domain's installer (set_<domain>_settings,
! <domain>_init for the domains that also load their tables, and the
! setters of the state that follows the surface mixture:
! set_convection_envelope_metals, eos_set_mixture,
! kap_update_surface_tables) is here, and tools/check_boundaries.py
! (INSTALL_FILES) keeps it so. Each procedure builds the records with
! setup/'s builders (<domain>_settings_from_star) or takes the star's
! surface mixture, and installs them, in a fixed order, one procedure
! (or, per surface mixture, three) per lifetime:
!   install_run_settings              per run, after the controls read
!                                     (run_setup's set_up_constants_and_tables):
!                                     kap, eos and atm, with their tables
!   install_monte_carlo_settings      per Monte-Carlo run, after
!                                     begin_monte_carlo_run (run_yrec): net
!   install_card_settings             per run card, after the card's
!                                     mixing length, before its model read
!                                     (begin_kind_card): convection,
!                                     rotation, element diffusion
!   per surface mixture, whenever the star's surface mixture changes:
!   install_surface_mixture_settings  the whole mixture, after the first
!                                     card's model read
!                                     (update_surface_mixture): convection's
!                                     envelope_metals_no_cno, known only
!                                     then, and the EOS mixture
!   install_surface_z_settings        its Z (run_setup's run-plan Z, before
!                                     kap_init; a Z rescale,
!                                     rescale_interior_z): the EOS mixture
!   install_surface_x_settings        its X (each card's model read; a
!                                     rezone that moved the surface X,
!                                     refresh_surface_opacity): kap's
!                                     surface opacity tables
! A new domain's record is installed by adding its call to the procedure
! of its lifetime (setup/README.md, "Add a physics domain").
module settings_install
      use star_info_lib, only: star_info
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: install_run_settings, install_monte_carlo_settings, &
           install_card_settings, install_surface_mixture_settings, &
           install_surface_z_settings, install_surface_x_settings

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
! install_surface_mixture_settings: the surface mixture the first card's
! model gives (update_surface_mixture): convection's
! envelope_metals_no_cno, its surface Z without C, N and O
! (install_card_settings installed the value the star had before: 0 on
! the first card), and the EOS mixture (install_surface_z_settings).
!   star reads: envelope_metals_no_cno
!   star writes: none
subroutine install_surface_mixture_settings(s)
      use convection_lib, only: set_convection_envelope_metals
      type(star_info), intent(in) :: s

      call set_convection_envelope_metals(s%envelope_metals_no_cno)
      call install_surface_z_settings(s)
end subroutine install_surface_mixture_settings

!----------------------------------------------------------------------
! install_surface_z_settings: the EOS mixture (eos_set_mixture) for the
! star's envelope X, Z, inverse mean atomic weight and number fractions:
! with SCV on, the EOS builds its envelope table for that X and Z. Called
! whenever the envelope Z changes: the run plan's Z (run_setup), a Z
! rescale (rescale_interior_z), and with the whole mixture.
!   star reads: envelope_hydrogen_fraction,
!     envelope_inverse_atomic_weight, envelope_metal_fraction,
!     envelope_number_fractions
!   star writes: none
subroutine install_surface_z_settings(s)
      use eos_lib, only: eos_set_mixture
      type(star_info), intent(in) :: s

      call eos_set_mixture(s%envelope_hydrogen_fraction, &
           s%envelope_metal_fraction, s%envelope_inverse_atomic_weight, &
           s%envelope_number_fractions)
end subroutine install_surface_z_settings

!----------------------------------------------------------------------
! install_surface_x_settings: kap's surface opacity tables (the fixed-X
! slices of the OPAL95, OPAL92 and ALEX94 tables) rebuilt for the
! surface X hydrogen_fraction: each card's model read passes the first
! card's envelope X (update_surface_mixture), a rezone the current
! surface X once it has moved (refresh_surface_opacity). Table errors
! return through ierr.
subroutine install_surface_x_settings(hydrogen_fraction, ierr)
      use kap_lib, only: kap_update_surface_tables
      real(dp), intent(in) :: hydrogen_fraction
      integer, intent(out) :: ierr

      call kap_update_surface_tables(hydrogen_fraction, ierr=ierr)
end subroutine install_surface_x_settings

end module settings_install

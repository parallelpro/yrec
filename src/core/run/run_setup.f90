!----------------------------------------------------------------------
! run_setup
!----------------------------------------------------------------------
! The setup before the run cards: the run plan's envelope Z, the
! physical constants and the microphysics tables
! (set_up_constants_and_tables), and the Monte-Carlo run range
! (prepare_run).
module run_setup
      use star_info_lib, only: star_info
      use control_slots, only: i_target_z
      implicit none
      private
      public :: prepare_run, set_up_constants_and_tables

contains

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! prepare_run (formerly star_setup)
!----------------------------------------------------------------------
! The setup step before the run cards: the run plan's envelope Z
! (set_run_plan_envelope_z), the constants and the table loads
! (set_up_constants_and_tables, which passes ierr through the kap, eos
! and atm initializations), then the Monte-Carlo run range and each
! run's sampled parameters. Errors return through ierr.
!   star reads: none
!   star writes: none
subroutine prepare_run(s, ierr)

      use monte_carlo, only: setup_monte_carlo_runs
      implicit none
      type(star_info), intent(inout) :: s

      integer, intent(out) :: ierr

      ierr = 0
! The envelope Z the run plan asks for, which the opacity tables are
! read for, and the EOS mixture with it.
      call set_run_plan_envelope_z(s)
! Set the constants and read the tables.
      call set_up_constants_and_tables(s, ierr)
      if (ierr /= 0) return
! The Monte-Carlo run range and each run's sampled parameters
! (core/run/monte_carlo.f90).
      call setup_monte_carlo_runs(s)

      return
end subroutine prepare_run

! ---------------------------------------------------------------
! The envelope Z of the run plan: each card with a positive rescale Z
! sets it (0, the default, keeps it), so the last such card's wins, as
! PARMIN's did; kap_init reads it, and the EOS mixture takes it with the
! star's envelope mixture as it stands (settings_install's
! install_surface_z_settings; the starting model sets both again).
!   star reads: card_rescale_targets, num_cards
!   star writes: envelope_metal_fraction
subroutine set_run_plan_envelope_z(s)
      use settings_install, only: install_surface_z_settings
      type(star_info), intent(inout) :: s
      integer :: n

      do n = 1, s%num_cards
         if (s%card_rescale_targets(i_target_z,n) > 0.0d0) &
              s%envelope_metal_fraction = s%card_rescale_targets(i_target_z,n)
      end do
      call install_surface_z_settings(s)
end subroutine set_run_plan_envelope_z

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! set_up_constants_and_tables (legacy setups.f)
!----------------------------------------------------------------------
! One-time, model-independent setup at the start of a run: log10 G
! (physical_constants' cgl) and the solar constants on star, and the
! per-run settings records with their tables
! (settings_install's install_run_settings: the opacities, the equation
! of state and the tabulated surface boundary conditions).
!
! The table paths are controls (star%control), so the caller must run
! read_controls first. Errors return through ierr.
!   star reads: physics%G_cgs, physics%solar_luminosity_cgs,
!     physics%solar_mass_cgs, physics%solar_radius_cgs
!   star writes: ln10_over_lsun, log10_solar_luminosity,
!     log10_solar_radius, solar_luminosity_cgs, solar_mass_cgs,
!     solar_radius_cgs
subroutine set_up_constants_and_tables(s, ierr)
      use settings_install, only: install_run_settings
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

! --- Constants -------------------------------------------------------
! (the mathematical and physical constants are parameters of
! physical_constants; only log10 G depends on a control.)
! The solar units, from the controls (a Monte-Carlo run rescales Lsun).
! Solar luminosity: log10 Lsun, and ln10/Lsun (not a log; s/erg).
      s%solar_luminosity_cgs = s%physics%solar_luminosity_cgs
      s%log10_solar_luminosity = log10(s%solar_luminosity_cgs)
      s%ln10_over_lsun = ln10/s%solar_luminosity_cgs
! Solar mass (a namelist control; default 1.9891d33 g).
      s%solar_mass_cgs = s%physics%solar_mass_cgs
! Solar radius.
      s%solar_radius_cgs = s%physics%solar_radius_cgs
      s%log10_solar_radius = log10(s%solar_radius_cgs)
! log10 G: the control G_cgs when positive, else the CODATA 1986 value.
      cgl = log10_newton_g(s%physics%G_cgs)
! --- The per-run settings records and their tables ------------------
      call install_run_settings(s, ierr)
end subroutine set_up_constants_and_tables

end module run_setup

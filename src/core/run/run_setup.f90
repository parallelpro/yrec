!----------------------------------------------------------------------
! run_setup
!----------------------------------------------------------------------
! The setup before the run cards: the physical constants and the
! microphysics tables (set_up_constants_and_tables), and the Monte-Carlo
! run range (prepare_run).
module run_setup
      use star_info_lib, only: star_info
      implicit none
      private
      public :: prepare_run, set_up_constants_and_tables

contains

! ---------------------------------------------------------------
!----------------------------------------------------------------------
! prepare_run (formerly star_setup)
!----------------------------------------------------------------------
! The setup step before the run cards: the constants and the table
! loads (set_up_constants_and_tables, which passes ierr through the
! kap, eos and atm initializations), then the Monte-Carlo run range and
! each run's sampled parameters. Errors return through ierr.
!   star reads: none
!   star writes: none
subroutine prepare_run(s, ierr)

      use monte_carlo, only: setup_monte_carlo_runs
      implicit none
      type(star_info), intent(inout) :: s

      integer, intent(out) :: ierr

      ierr = 0
! Set the constants and read the tables.
      call set_up_constants_and_tables(s, ierr)
      if (ierr /= 0) return
! The Monte-Carlo run range and each run's sampled parameters
! (core/run/monte_carlo.f90).
      call setup_monte_carlo_runs(s)

      return
end subroutine prepare_run

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
! The table paths are controls (star%ctrl), so the caller must run
! read_controls first. Errors return through ierr.
!   star reads: ctrl%G_cgs, ctrl%solar_luminosity_cgs,
!     ctrl%solar_mass_cgs, ctrl%solar_radius_cgs
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
      s%solar_luminosity_cgs = s%ctrl%solar_luminosity_cgs
      s%log10_solar_luminosity = log10(s%solar_luminosity_cgs)
      s%ln10_over_lsun = ln10/s%solar_luminosity_cgs
! Solar mass (a namelist control; default 1.9891d33 g).
      s%solar_mass_cgs = s%ctrl%solar_mass_cgs
! Solar radius.
      s%solar_radius_cgs = s%ctrl%solar_radius_cgs
      s%log10_solar_radius = log10(s%solar_radius_cgs)
! log10 G: the control G_cgs when positive, else the CODATA 1986 value.
      cgl = log10_newton_g(s%ctrl%G_cgs)
! --- The per-run settings records and their tables ------------------
      call install_run_settings(s, ierr)
end subroutine set_up_constants_and_tables

end module run_setup

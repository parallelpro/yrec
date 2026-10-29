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
! (physical_constants' cgl) and the solar constants on star, the HSRA
! T(tau) offset, and the tables, loaded through the domain
! initializations: the opacities (kap_init); the equation of state
! (eos_init: the Fermi-Dirac degenerate-electron table, and the OPAL,
! MHD and SCVH tables when selected); and the tabulated surface
! boundary conditions (atm_init: the Kurucz, Allard or Kurucz/Castelli
! tables for atm_choice 3-5).
!
! The table paths and the mixture weights are read from star%job, so
! the caller must fill star%job (read_controls) first. The paths are
! passed as character arguments rather than through common blocks,
! where character strings can cause problems (MHP 8/25). Errors return
! through ierr.
subroutine set_up_constants_and_tables(s, ierr)
      use eos_lib
      use kap_lib
      use microphysics_settings, only: opacity_settings, equation_of_state_settings, &
           equation_of_state_table_paths, atmosphere_settings
      use atm_lib, only: atm_init, hsra_t_tau_offset
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

! --- Constants -------------------------------------------------------
! (the mathematical and physical constants are parameters of
! physical_constants; only log10 G depends on a control.)
! Solar luminosity: log10 Lsun, and ln10/Lsun (not a log; s/erg).
      s%log10_solar_luminosity = log10(s%solar_luminosity_cgs)
      s%ln10_over_lsun = ln10/s%solar_luminosity_cgs
! Solar mass (a namelist control; default 1.9891d33 g).
      s%solar_mass_cgs = s%ctrl%solar_mass_cgs
! Solar radius.
      s%log10_solar_radius = log10(s%solar_radius_cgs)
! Bolometric magnitude of the Sun.
      s%solar_bolometric_magnitude = 4.79d0
! log10 G: the control G_cgs when positive, else the CODATA 1986 value.
      cgl = log10_newton_g(s%ctrl%G_cgs)
! --- HSRA T(tau) offset -----------------------------------------------
! The HSRA T(tau) fit evaluated at tau = 2/3, used by atm_choice = 2.
      s%atm_hras = hsra_t_tau_offset()
! --- Opacity tables ---------------------------------------------------
      ierr = 0
      call kap_init(opacity_settings(s), s%envelope_hydrogen_fraction, &
           s%envelope_metal_fraction, s%job%mixture_weights, &
           s%job%alex06_table_path,s%job%kurucz_table_path,s%job%kurucz_table2_path, &
           s%job%laol_table_path,s%job%laol_table2_path, &
           s%job%opal95_table_path,s%job%opal92_table_path,s%job%opal92_table2_path, &
           s%job%pure_z_table_path,s%job%alex95_table_paths, &
           s%job%conductive_table_path, ierr)
      if (ierr /= 0) return
! --- Equation-of-state tables -----------------------------------------
! eos_init reads the Fermi-Dirac degenerate-electron table always, and
! the OPAL, MHD and SCVH tables when the settings select them.
      call eos_init(equation_of_state_settings(s), equation_of_state_table_paths(s), ierr)
      if (ierr /= 0) return
! --- Surface-pressure tables for the surface boundary condition -------
! atm_init reads the Kurucz (atm_choice 3 and 4), Allard (atm_choice 4)
! and Kurucz/Castelli (atm_choice 5) tables.
      call atm_init(atmosphere_settings(s), s%job%atm_table_path, &
           s%job%allard_table_path, ierr)
      if (ierr /= 0) return

      return
end subroutine set_up_constants_and_tables

end module run_setup

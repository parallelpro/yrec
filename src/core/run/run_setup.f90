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
! One-time, model-independent setup at the start of a run: the physical
! and derived constants (in physical_constants and star%), the HSRA
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
      use microphysics_settings, only: opacity_settings, equation_of_state_settings, atmosphere_settings
      use atm_lib, only: atm_init, hsra_t_tau_offset
      use physical_constants
      use math_lib
      implicit none
      type(star_info), intent(inout) :: s
      integer, intent(out) :: ierr

! --- locals ---
      double precision :: speed_of_light, electron_mass, boltzmann_constant, &
           planck_constant

! --- Constants -------------------------------------------------------
! (ln10, cc13, cpi and the other mathematical constants are parameters
! of physical_constants.)
! Solar luminosity: log10 Lsun, and ln10/Lsun (not a log; s/erg).
      s%log10_solar_luminosity = log10(s%solar_luminosity_cgs)
      s%ln10_over_lsun = ln10/s%solar_luminosity_cgs
! Solar mass (a namelist control; default 1.9891d33 g).
      s%solar_mass_cgs = s%ctrl%solar_mass_cgs
! Solar radius.
      s%log10_solar_radius = log10(s%solar_radius_cgs)
! Bolometric magnitude of the Sun.
      s%solar_bolometric_magnitude = 4.79d0
! Speed of light.
      speed_of_light = 2.99792458d10
! Stefan-Boltzmann constant sigma, and its log.
      csig = 5.67051d-5
      csigl = log10(csig)
! a/3 = 4 sigma/(3 c), with a the radiation constant, and its log.
      radiation_constant_over_3 = csig*4.0d0/(3.0d0*speed_of_light)
      ca3l = log10(radiation_constant_over_3)
! Molar gas constant.
      gas_constant = 8.314510d7
! log10 G. A non-positive G_cgs (the default) keeps the historical
! hard-coded log10 G = -7.17571 (G = 6.6726d-8) bit for bit; a
! positive G_cgs (a namelist control) overrides it.
      if (s%ctrl%G_cgs > 0.0d0) then
         cgl = log10(s%ctrl%G_cgs)
      else
         cgl = -7.17571d0
      end if
! Electron mass.
      electron_mass = 9.1093897d-28
! Boltzmann constant.
      boltzmann_constant = 1.380658d-16
! Planck constant.
      planck_constant = 6.6260755d-27
! The Saha-equation constant:
! cmkh = log10[ (2 pi m_e)^(3/2) k^(5/2) / h^3 ].
      cmkh = 1.50d0*log10(2.0d0*cpi*electron_mass)+2.5d0* &
           log10(boltzmann_constant)- &
           3.0d0*log10(planck_constant)
! The radiative-gradient constant, cdelrl = log10[ 3/(64 pi sigma) ],
! so that grad_rad = 10^cdelrl kappa L P / (G m T^4).
      cdelrl = -c4pi3l - csigl - log10(16.0d0)
! Mixing-length-theory constants: 1/3 and 16 sqrt(2) sigma.
      cmixl2 = cc13
      cmixl3 = 16.0d0*dsqrt(2.0d0)*csig
! (The Debye-Hueckel constants are set in eos_init.)
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
      call eos_init(equation_of_state_settings(s), s%job%fermi_table_path,s%job%scv_h_table_path,s%job%scv_he_table_path, &
           s%job%scv_z_table_path,s%job%zams_a_table_path,s%job%zams_b_table_path, &
           s%job%zams_c_table_path,s%job%centre1_table_path,s%job%centre2_table_path, &
           s%job%centre3_table_path,s%job%centre4_table_path,s%job%centre5_table_path, &
           s%job%opal_eos_table_path, ierr)
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

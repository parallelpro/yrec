!----------------------------------------------------------------------
! setups
!----------------------------------------------------------------------
! Modernized (free-form, readable names) 2026 as part of the YREC
! readability refactor. Logic and numerics are unchanged from the
! original setups.f; only variable names, source form, and comment
! style were updated. Validated against the Stage 0 regression suite
! (examples/run_standard_solar_model).
!
! MHP 8/25 Passed file names directly instead of through common, as
! character strings in common blocks can cause problems.
!
! One-time model-independent setup performed at program start:
! defines the physical/mathematical constants used throughout the
! code (physical_constants and star%), loads the opacity tables (setupopac/mhdtbl), the degenerate-electron
! (Fermi-Dirac) equation-of-state table used by fully_ionized_eos.f90, the
! Kurucz/Castelli surface-pressure table used for the T-tau surface
! boundary condition, and (if enabled) the SCVH envelope
! equation-of-state tables.
subroutine setups(ierr)
      use eos_lib
      use kap_lib
      use microphysics_settings, only: opacity_settings, equation_of_state_settings
      use atm_lib
      use atm_tables
      use star_info_lib, only: star
      use physical_constants
      use math_lib
      use t_tau_relations, only: hsra_t_tau_offset
      implicit none
      integer, intent(out) :: ierr

! 2026 de-tramp: the 24 table paths and the mixture work array come
! from star%job directly (they were only ever star%job members passed
! positionally through star_setup) -- callers set star%job first.
! --- locals ---
      double precision :: speed_of_light, electron_mass, boltzmann_constant, &
           planck_constant, hydrogen_atom_mass, electron_charge_esu

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SETUP CONSTANTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! (ln10, cc13, cpi and the other mathematical constants are parameters
! of physical_constants)
! Luminosity of Sun
      star%log10_solar_luminosity = log10(star%solar_luminosity_cgs)
      star%ln10_over_lsun = ln10/star%solar_luminosity_cgs
! Mass of Sun (namelist control since 2026; default 1.9891d33)
      star%solar_mass_cgs = star%ctrl%solar_mass_cgs
      star%log10_solar_mass = log10(star%solar_mass_cgs)
! Radius of Sun
      star%log10_solar_radius = log10(star%solar_radius_cgs)
! Bolometric magnitude of the sun
      star%solar_bolometric_magnitude = 4.79d0
! Speed of Light
      speed_of_light = 2.99792458d10
! Stefan-Boltzmann constant
      csig = 5.67051d-5
      csigl = log10(csig)
! radiation constant/3
      radiation_constant_over_3 = csig*4.0d0/(3.0d0*speed_of_light)
      ca3l = log10(radiation_constant_over_3)
! molar gas constant
      gas_constant = 8.314510d7
! log of the gravitational constant. G_cgs < 0 (default) keeps the
! historical hard-coded log10 G = -7.17571 (G=6.6726D-8) bit-for-bit;
! a positive G_cgs (namelist control since 2026) overrides it.
      if (star%ctrl%G_cgs > 0.0d0) then
         cgl = log10(star%ctrl%G_cgs)
      else
         cgl = -7.17571d0
      end if
! mass of the electron
      electron_mass = 9.1093897d-28
! Boltzmann constant
      boltzmann_constant = 1.380658d-16
! Planck's constant
      planck_constant = 6.6260755d-27
      cmkh = 1.50d0*log10(2.0d0*cpi*electron_mass)+2.5d0* &
           log10(boltzmann_constant)- &
           3.0d0*log10(planck_constant)
      cdelrl = -c4pi3l - csigl - log10(16.0d0)
      cmixl2 = cc13
      cmixl3 = 16.0d0*dsqrt(2.0d0)*csig
! (the Debye-Hueckel constants are eos_init's, 2026-10)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     EVALUATE TAU = 2/3 TEMPERATURE FOR HRA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      star%atm_hras = hsra_t_tau_offset()
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SET UP OPACITY TABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ierr = 0
      call kap_init(opacity_settings(), star%envelope_hydrogen_fraction, &
           star%envelope_metal_fraction, star%job%mixture_weights, &
           star%job%alex06_table_path,star%job%kurucz_table_path,star%job%kurucz_table2_path, &
           star%job%laol_table_path,star%job%laol_table2_path, &
           star%job%opal95_table_path,star%job%opal92_table_path,star%job%opal92_table2_path, &
           star%job%pure_z_table_path,star%job%alex95_table_paths, &
           star%job%conductive_table_path, ierr)
      if (ierr /= 0) return
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     READ IN MHD EOS TABLES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MHP 8/25 Passed file names directly instead of through common, as
! character strings in common blocks can cause problems
! 2026 (docs/ROADMAP.md stage 1): the use_mhd_eos gate moved inside
! eos_lib's eos_init; a no-op when MHD is off, exactly as before.
! eos_init also reads the Fermi-Dirac degenerate-electron table
! (the F-tables block that used to sit inline just below).
      call eos_init(equation_of_state_settings(), star%job%fermi_table_path,star%job%scv_h_table_path,star%job%scv_he_table_path, &
           star%job%scv_z_table_path,star%job%zams_a_table_path,star%job%zams_b_table_path, &
           star%job%zams_c_table_path,star%job%centre1_table_path,star%job%centre2_table_path, &
           star%job%centre3_table_path,star%job%centre4_table_path,star%job%centre5_table_path, &
           star%job%opal_eos_table_path, ierr)
      if (ierr /= 0) return
!
!
! (The F-tables read for the degenerate equation of state that lived
! here moved into eos_lib's eos_init above -- 2026, docs/ROADMAP.md
! stage 1 -- along with its state's relocation from atm_tables to
! yale_eos_tables.)
! JMH 8/18/91
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! INPUT PRESSURE TABLE FOR SURFACE BOUNDARY CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 2026 (docs/ROADMAP.md stage 1): the Kurucz (atm_choice 3/4), Allard
! (atm_choice 4), and Kurucz/Castelli (atm_choice 5) surface table
! loads that lived inline here moved into atm_lib's atm_init.
      call atm_init(star%job%atm_table_path,star%job%allard_table_path, ierr)
      if (ierr /= 0) return

! (The SCV equation-of-state table reads that lived here moved into
! eos_lib's eos_init -- 2026, docs/ROADMAP.md stage 1 -- alongside the
! Fermi-Dirac table read.)

      return
end subroutine setups

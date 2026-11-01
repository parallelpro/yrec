!----------------------------------------------------------------------
! eos_configuration
!----------------------------------------------------------------------
! The settings record of eos/: everything the equation of state needs
! from the run's configuration, installed by eos_init (from the star
! layer's eos_settings_from_star, setup/eos_settings_builder.f90)
! or set_eos_settings; nothing else writes it. Its defaults are the
! controls' (defaults/controls_registry.tsv), so a star-free caller that
! sets only what it changes gets what a deck that leaves the controls
! alone gets.
module eos_configuration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use mhd_eos_tables, only: mhd_n_zams, mhd_n_centre
      implicit none
      private
! The OPAL table vintages eos_settings_t%opal_vintage selects:
! opal_vintage_conflict is the builder's answer to controls that select
! more than one, an error at eos_init.
      integer, parameter, public :: opal_vintage_none = 0, opal_vintage_1995 = 1, &
           opal_vintage_2001 = 2, opal_vintage_2006 = 3, opal_vintage_conflict = -1

      type, public :: eos_settings_t
! the MHD tables instead of everything below; the OPAL table vintage
! (opal_vintage_*); the SCV EOS at low temperature
         logical :: use_mhd_eos = .false.
         integer :: opal_vintage = opal_vintage_none
         logical :: use_scv_eos = .false.
! derivatives by finite differences rather than analytically; the metal
! fraction varies (diffusion of Z), so the mixture is rescaled per zone
         logical :: use_numerical_derivatives = .false., use_diffusion_z = .false.
! the Debye-Hueckel correction and the degeneracy (eta) range over which
! it ramps from its nondegenerate to its degenerate form (eta_max >
! eta_min, which eos_init checks)
         logical :: use_debye_huckel_correction = .false.
         real(dp) :: debye_huckel_eta_min = -1.0d0, debye_huckel_eta_max = 1.0d0
! log10 T (K) at and above which the gas is fully ionized; the Saha
! solution below it, blended into the fully ionized one over the 0.5 dex
! under it
         real(dp) :: saha_log10t_cutoff = 6.0d0
      end type eos_settings_t
      type(eos_settings_t), public, protected :: eos_settings

! The tables eos_init reads, by path: the Fermi-Dirac table (always), the
! SCV tables (use_scv_eos), the MHD tables (use_mhd_eos: the three
! ZAMS-type compositions and the five centre-type ones) and the OPAL
! table of the selected vintage. The star layer fills it from the
! controls (eos_table_paths_from_star, setup/eos_settings_builder.f90).
      type, public :: eos_table_paths_t
         character(len=256) :: fermi = ''
         character(len=256) :: scv_h = '', scv_he = '', scv_z = ''
         character(len=256) :: mhd_zams(mhd_n_zams) = '', mhd_centre(mhd_n_centre) = ''
         character(len=256) :: opal = ''
      end type eos_table_paths_t
      public :: set_eos_settings

contains

! set_eos_settings: install the run's EOS settings.
subroutine set_eos_settings(settings)
      type(eos_settings_t), intent(in) :: settings
      eos_settings = settings
end subroutine set_eos_settings

end module eos_configuration

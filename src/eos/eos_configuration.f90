!----------------------------------------------------------------------
! eos_configuration
!----------------------------------------------------------------------
! The settings record of eos/: everything the equation of state needs
! from the run's configuration, installed by eos_init (from the star
! layer's equation_of_state_settings, setup/microphysics_settings.f90)
! or set_eos_settings; nothing else writes it.
module eos_configuration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      type, public :: eos_settings_t
! the MHD tables instead of everything below; the OPAL table vintage
! (1995, 2001 or 2006; at most one); the SCV EOS at low temperature
         logical :: use_mhd_eos = .false.
         logical :: use_opal95_eos = .false., use_opal2001_eos = .false., &
              use_opal2006_eos = .false.
         logical :: use_scv_eos = .false.
! derivatives by finite differences rather than analytically; the metal
! fraction varies (diffusion of Z), so the mixture is rescaled per zone
         logical :: use_numerical_derivatives = .false., use_diffusion_z = .false.
! the Debye-Hueckel correction and the degeneracy (eta) range it ramps
! over
         logical :: use_debye_huckel_correction = .false.
         real(dp) :: debye_huckel_eta_min = 0.0d0, debye_huckel_eta_max = 0.0d0
! log10 T below which the Saha solution is used (ramped)
         real(dp) :: saha_log10t_cutoff = 0.0d0
      end type eos_settings_t
      type(eos_settings_t), public, protected :: eos_settings
      public :: set_eos_settings

contains

! set_eos_settings: install the run's EOS settings.
subroutine set_eos_settings(settings)
      type(eos_settings_t), intent(in) :: settings
      eos_settings = settings
end subroutine set_eos_settings

end module eos_configuration

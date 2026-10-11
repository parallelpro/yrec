!----------------------------------------------------------------------
! convection_configuration
!----------------------------------------------------------------------
! The convection domain's settings: everything it reads besides its
! arguments, installed by the star layer through convection_lib's
! set_convection_settings at the start of every run card (the mixing
! length is per card) and updated by set_convection_envelope_metals
! when the starting model is read.
module convection_configuration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: convection_settings_t, convection_settings

      type :: convection_settings_t
! mixing-length theory
         real(dp) :: mixing_length_alpha = 0.0d0
! skip mixing-length theory above this log T (cores: the gradient is
! adiabatic)
         real(dp) :: log10_t_cut = 0.0d0
! the adiabatic gradient in the overshoot region below the envelope
! (ctrl%ladov)
         logical :: adiabatic_overshoot = .false.
! the Sun's luminosity, cgs and log10
         real(dp) :: solar_luminosity_cgs = 0.0d0, log10_solar_luminosity = 0.0d0
! starspots (G Somers 9/14): filling factor, temperature contrast, and
! whether the contrast varies with depth
         real(dp) :: spot_filling_factor = 0.0d0, spot_temp_contrast = 0.0d0
         logical :: spot_depth_varies = .false.
! overshoot of the core, the envelope and intermediate zones (the job
! flags lovstm etc.) and its extent in pressure scale heights; lovmax
! limits core overshoot to betac times the radius (Woo & Demarque 2001)
         logical :: core_overshoot = .false., envelope_overshoot = .false., &
              intermediate_overshoot = .false.
         real(dp) :: overshoot_alpha_core = 0.0d0, overshoot_alpha_envelope = 0.0d0, &
              overshoot_alpha_intermediate = 0.0d0
         logical :: limit_core_overshoot = .false.
         real(dp) :: core_overshoot_limit = 0.0d0
! rotating models with instability transport keep at least one
! radiative point in an overshoot region
         logical :: rotation_active = .false., instability_transport_active = .false.
! homogenize the extended composition (15 species, not 11)
         logical :: use_extended_composition = .false.
! the envelope's metal fraction outside C, N and O (star%zenvm), the
! floor of species 9 after homogenizing
         real(dp) :: envelope_metals_no_cno = 0.0d0
      end type convection_settings_t

      type(convection_settings_t), save :: convection_settings
end module convection_configuration

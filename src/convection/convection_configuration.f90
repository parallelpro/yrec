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
      use physical_constants, only: log10_newton_g_default
      implicit none
      private
      public :: convection_settings_t, convection_settings

! The codes convection's routines return in ierr, besides the EOS's and
! the opacity's, which they pass up as they are: the mixing-length
! square root of a negative number (compute_temperature_gradients), and
! more convection zones than the caller's bounds array has rows
! (find_convection_zones).
      integer, parameter, public :: convection_unstable_radicand = 1, &
           convection_too_many_zones = 2

      type :: convection_settings_t
! mixing-length theory
         real(dp) :: mixing_length_alpha = 0.0d0
! skip mixing-length theory above this log T (cores: the gradient is
! adiabatic)
         real(dp) :: log10_t_cut = 0.0d0
! the adiabatic gradient in the overshoot region below the envelope
! (physics%adiabatic_overshoot_active)
         logical :: adiabatic_overshoot = .false.
! the Sun's luminosity, cgs and log10
         real(dp) :: solar_luminosity_cgs = 0.0d0, log10_solar_luminosity = 0.0d0
! log10 G (cgs): the control G_cgs, or the CODATA 1986 value
         real(dp) :: log10_gravitational_constant = log10_newton_g_default
! starspots (Somers 2014): filling factor, temperature contrast, and
! whether the contrast varies with depth
         real(dp) :: spot_filling_factor = 0.0d0, spot_temp_contrast = 0.0d0
         logical :: spot_depth_varies = .false.
! overshoot of the convective core, the surface convection zone and the
! intermediate zones (the controls core_overshoot_active,
! envelope_overshoot_active and intermediate_overshoot_active), each by
! alpha pressure scale heights (overshoot_alpha_core,
! overshoot_alpha_envelope, overshoot_alpha_min); limit_core_overshoot
! (overshoot_max_extent_flag) caps the core's at core_overshoot_limit
! (overshoot_beta_core) times the core's radius (Woo & Demarque 2001)
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
! the envelope's metal fraction outside C, N and O (star%envelope_metals_no_cno), the
! floor of species 9 after homogenizing
         real(dp) :: envelope_metals_no_cno = 0.0d0
      end type convection_settings_t

      type(convection_settings_t), save :: convection_settings
end module convection_configuration

!----------------------------------------------------------------------
! rotation_configuration
!----------------------------------------------------------------------
! The rotation domain's settings: every control it reads, installed by
! the star layer through set_rotation_settings at the start of every
! run card (setup/rotation_settings_setup builds the record from
! star%ctrl and star%job). Rotation reads no control directly
! (tools/check_boundaries.py); its state lives on star% and in
! rotation_scratch_lib.
module rotation_configuration
      use, intrinsic :: iso_fortran_env, only: dp => real64
      implicit none
      private
      public :: rotation_settings_t, rotation_settings, set_rotation_settings

      type :: rotation_settings_t
! instability transport (secular_transport) on. The step is cut into
! sub-steps that each change omega by at most diffusion_timestep_factor
! (fractionally), counting the step's largest change of omega, or
! max_domega_global when the model has none; the wind's sub-steps use
! max_domega_global.
         logical :: instability_transport_active = .false.
         real(dp) :: diffusion_timestep_factor = 0.0d0, max_domega_global = 0.0d0
! the transport iterations: their limit and the convergence tolerance
! on the angular momentum
         integer :: max_diffusion_iters = 0
         real(dp) :: convergence_tolerance = 0.0d0
! the shape of the rotating star (rotating_shape): the iteration limits
! on the equipotential radius (outer and inner loops) and the relative
! tolerance on r0^3 (the former job%acfpft; check_angular_momentum
! tightens it for one call)
         integer :: shape_outer_iterations = 0, shape_inner_iterations = 0
         real(dp) :: shape_tolerance = 0.0d0
! the rotation law in convection zones: omega ~ r^cz_rotation_law_exponent
! for exponents between -2 and 0, solid body at 0 and above, constant
! specific angular momentum at -2 and below (the former walpcz); the
! whole star a solid body (force_solid_body_rotation, or
! solid_body_mode_flag = 1). With a negative exponent, the surface
! convection zone follows the law and solid_body_mode_flag = 2 rotates
! the interior as a separate solid body, 3 as a solid body at the
! zone's base rate (see enforce_rotation_profile)
         real(dp) :: cz_rotation_law_exponent = 0.0d0
         logical :: force_solid_body_rotation = .false.
         integer :: solid_body_mode_flag = 0
! the circulation and instability velocities: their scales, the GSF
! inhibition mode, the critical Reynolds number, the Zahn (1992)
! coupling factor (the former lvfc), no transport in a convective core
         real(dp) :: es_velocity_scale = 0.0d0, gsf_velocity_scale = 0.0d0, &
              secular_shear_velocity_scale = 0.0d0, critical_reynolds = 0.0d0
         integer :: gsf_inhibition_mode = 0
         logical :: use_zahn_coupling_factor = .false.
         logical :: no_am_transport_in_core = .false.
! the mixing of the composition: the overall scale and the scales of
! the three mechanisms; a constant background diffusion coefficient
         real(dp) :: mixing_velocity_scale = 0.0d0, es_mixing_scale = 0.0d0, &
              gsf_mixing_scale = 0.0d0, secular_shear_mixing_scale = 0.0d0
         logical :: use_constant_background_diffusion = .false.
         real(dp) :: constant_background_diffusion_coeff = 0.0d0
! the advection-diffusion transport (Maeder & Zahn 1998), commented out:
! a rotating run that sets it stops; and the overall scale of the
! diffusion coefficients, D = 4 pi difad_velocity_scale (U_ES + U_GSF +
! U_SS) r (legacy FW)
         logical :: use_diffusion_advection_transport = .false.
         real(dp) :: difad_velocity_scale = 0.0d0
! the grid scales (the former hpttol, slots ichi_* in controls_buffer):
! the chi coordinate's mass, luminosity and pressure scales are slots
! ichi_dm_max, ichi_dl_max and ichi_dp_core_max
         real(dp) :: chi_grid_scale(12) = 0.0d0
! the wind (magnetic_braking): applied at all; the Kawaler law's
! constant and exponents (omega, mass-loss rate, radius, mass), and the
! saturation omega
         logical :: use_wind_torque = .false.
         real(dp) :: wind_constant = 0.0d0, wind_law_omega_exponent = 0.0d0, &
              wind_mass_loss_exponent = 0.0d0, wind_radius_exponent = 0.0d0, &
              wind_mass_exponent = 0.0d0, wind_saturation_omega = 0.0d0
! the Matt/PMM law: on; its exponents of luminosity, photospheric
! pressure, turnover time and the centrifugal term, and the centrifugal
! constant; the solar reference values; the Rossby-number and field
! scalings
         logical :: use_pmm_wind_law = .false.
         real(dp) :: wind_luminosity_exponent = 0.0d0, wind_pressure_exponent = 0.0d0, &
              wind_turnover_exponent = 0.0d0, wind_centrifugal_exponent = 0.0d0, &
              wind_centrifugal_constant = 0.0d0
         real(dp) :: pmm_solar_omega = 0.0d0, pmm_solar_pressure = 0.0d0, &
              pmm_solar_turnover_timescale = 0.0d0
         logical :: scale_by_rossby_number = .false., scale_by_b_field = .false.
! disk locking (MHP 9/94): on, until this age (Gyr), at this omega (rad/s)
         logical :: disk_locking_active = .false.
         real(dp) :: disk_locking_age_gyr = 0.0d0, disk_omega_rad_s = 0.0d0
! helium settling within the rotation sub-steps: the hydrogen and
! helium floors and the sub-step as a fraction of the settling time
         real(dp) :: hydrogen_diffusion_floor = 0.0d0, helium_diffusion_min = 0.0d0, &
              settling_timestep_fraction = 0.0d0
! the extended (15-species) composition; overshoot of intermediate
! convection zones (job%lovstm) and the envelope's overshoot extent
         logical :: use_extended_composition = .false., intermediate_overshoot = .false.
         real(dp) :: overshoot_alpha_envelope = 0.0d0
      end type rotation_settings_t

      type(rotation_settings_t), save :: rotation_settings

contains

! Installs the run card's rotation settings.
subroutine set_rotation_settings(settings)
      type(rotation_settings_t), intent(in) :: settings
      rotation_settings = settings
end subroutine set_rotation_settings

end module rotation_configuration

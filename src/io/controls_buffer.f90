!----------------------------------------------------------------------
! controls_buffer
!----------------------------------------------------------------------
! The run controls under their canonical names: the buffer the input
! decks are read into. read_controls seeds it from the pristine star%ctrl
! defaults (seed_controls_buffer), the namelists read over it, and
! adopt_canonical_names and map_user_inputs fill the members the decks
! set under other names; store_controls_to_star then copies it into
! star%ctrl and star%job, where the rest of the code reads its controls.
! The declarations below are the defaults: tools/gen_controls_state.py
! generates star%ctrl's type and the two copies (io/controls_sync.f90)
! from them, so regenerate after changing a member, and add the control
! to defaults/controls_registry.tsv, which maps it to its legacy and
! new-inlist names. Each control's default is stated once in the code:
! here; set_legacy_defaults (io/legacy_controls.f90) gives each legacy
! name its canonical twin's default, and holds only the defaults of the
! legacy names that have no canonical home. The legacy namelists are
! frozen: a new control is read only by the new format. Only the
! controls reader, map_user_inputs,
! controls_sync and test_net use this module; any file may import its
! named constants (max_runs, and the slot names itime_* and ichi_*,
! which are state/control_slots.f90's, re-exported here).
! tools/check_boundaries.py enforces both.

module controls_buffer
      use, intrinsic :: iso_fortran_env, only: dp => real64
      use control_slots
      implicit none
! Everything here is public: the buffer is the controls read's shared
! state (tools/check_boundaries.py's CONTROLS_BUFFER_IMPORTERS lists the
! files that may use it), and the named constants max_runs, itime_* and
! ichi_* are for any file.
      public
      private :: dp

! Maximum number of runs (kind cards) in one job: the dimension of
! every per-run array below (the calibration protocols in
! core/run/calibration.f90 / setup_solar_calibration.f90 fill
! them up to this limit). tools/gen_controls_state.py substitutes the
! literal value into the generated star%ctrl / star%job type bodies,
! which star_info_lib includes without seeing this module.
      integer, parameter :: max_runs = 50

! MESA-style output (new-format inlists only). The run-log, profile and
! pulsation-output switches below apply to both output modes.
      character(len=256) :: star_history_name = 'history.data'
      character(len=256) :: history_columns_file = ' '
      character(len=256) :: profile_columns_file = ' '
      integer :: profile_interval = 50
! Structure-limit stop conditions, per run like the other stopping
! criteria (target_end_age, central_*_stop); a limit at its +-1d99
! sentinel is disabled.
      real(dp) :: log_L_upper_limit(max_runs) = 1.0d99
      real(dp) :: log_L_lower_limit(max_runs) = -1.0d99
      real(dp) :: Teff_upper_limit(max_runs) = 1.0d99
      real(dp) :: Teff_lower_limit(max_runs) = -1.0d99
      real(dp) :: log_g_upper_limit(max_runs) = 1.0d99
      real(dp) :: log_g_lower_limit(max_runs) = -1.0d99
      real(dp) :: nu_max_upper_limit(max_runs) = 1.0d99
      real(dp) :: nu_max_lower_limit(max_runs) = -1.0d99
! The solar mass and Newton's G. G_cgs < 0 (the default) keeps the
! historical log10 G = -7.17571 (G = 6.6726e-8) bit-for-bit; set a
! positive value (e.g. 6.6743d-8, CODATA 2018) to override.
      real(dp) :: solar_mass_cgs = 1.9891d33
      real(dp) :: G_cgs = -1.0d0
! solar reference values for the asteroseismic scaling relations:
      real(dp) :: nu_max_sun = 3090.0d0
      real(dp) :: delta_nu_sun = 135.1d0
      real(dp) :: Teff_sun = 5772.0d0
      character(len=256) :: profile_data_prefix = 'profile'
! the copy of the deck; blank: inlist_used in the output directory
! (the run log's)
      character(len=256) :: inlist_used_file = ' '
      character(len=8) :: pulse_format = 'GYRE'
! Independent toggles: profiles only, pulse only, both, or neither
! (cadence stays with profile_interval / pulse_gyre_interval):
      logical :: write_profile_flag = .false.
      logical :: write_pulse_flag = .false.
! Run-log verbosity: the compact per-model progress line prints every
! terminal_interval models (0 = off; the final model always prints),
! and report_solver_diagnostics restores the verbose solver forensics
! (Henyey iteration trace, envelope/atmosphere integration statistics,
! triangle/rezoning bookkeeping).
      integer :: terminal_interval = 50
! model numbering for fresh start-model loads: the first converged
! model after a first_call_flag load is numbered this (default 1);
! <= 0 keeps the counter stored in the model file (restart
! continuation, the pre-2026 behavior).
      integer :: set_initial_model_number = 1
      logical :: report_solver_diagnostics = .false.
! diagnostic: per-zone solar neutrino production table (compute_zone_energy per
! shell) written to the short/log stream at the start of each run
      logical :: compute_neutrino_fluxes = .false.

! Timestep factors (atime) and log T cutoffs: nuclear burning, Saha
! ionization, the envelope's temperature window (tenv0, tenv1; star%envelope_store_log_temperature
! is their mean) and the simplified temperature gradient (tgcut).
      real(dp) :: atime(14) = (/1.0d-3,2.0d-2,5.0d-1,2.0d-2, &
           3.0d-1,1.5d-3,1.0d-1,2.0d-2,4.0d-2,2.0d-2,2.0d-2,0.25d0, &
           1.5d0,0.25d0/)
      real(dp) :: nuclear_logT_cutoffs(5) = (/6.5d0,6.5d0,6.82d0,7.7d0,7.5d0/)
      real(dp) :: saha_log10t_cutoff = 6.0d0
      real(dp) :: tenv0 = 3.0d0, tenv1 = 9.0d0, tgcut = 6.9d0

! Nuclear rate selection and weak screening.
      logical :: use_new_nuclear_rates = .false.
      real(dp) :: weak_screening_threshold = 0.03d0

! Overshoot and semiconvection. The overshoot range is model state
! (star%envelope_overshoot_bottom_index/envelope_overshoot_top_index).
      real(dp) :: dpenv = 1.0d0
      logical :: core_overshoot_active = .false.
      real(dp) :: overshoot_alpha_core = 0.0d0
      logical :: envelope_overshoot_active = .false.
      real(dp) :: overshoot_alpha_envelope = 0.0d0
      logical :: lovstm = .false.
      real(dp) :: alpham = 0.0d0
      logical :: ladov = .false.
      logical :: lovmax = .false.
      real(dp) :: betac = 0.15d0
      logical :: use_semiconvection = .false.

! Rotation. The reader clamps walpcz to [-2, 0].
      real(dp) :: acfpft = 1.0d-36
      integer :: itfp1 = 5, itfp2 = 20
      logical :: rotation_active = .false.
      real(dp) :: walpcz = 0.0d0
      logical :: instability_transport_active = .false., set_initial_omega = .false.
      real(dp) :: initial_omega = 0.0d0

! Mass accretion and the Reimers wind. core/mass_change switches
! star%job%use_mass_accretion off when the disk is exhausted.
      real(dp) :: mass_accretion_rate = -1.0d-14
      real(dp) :: fczdmdt = 0.1d0, ftotdmdt = 1.0d-2
      real(dp) :: accreted_composition(15) = (/0.71668d0, &
           0.265721d0,0.01757d0,2.9d-5,3.013d-3,3.385d-5,9.346d-4, &
           0.0d0,8.462d-3,0.0d0,1.696d-5,0.0d0,2.0d-9,2.0d-9,3.0d-11/)
      real(dp) :: reimers_scaling_factor = -4.0d-13
      logical :: use_reimers_wind = .false., use_mass_accretion = .false.

! The level-4 Henyey iteration limit, the first-guess extrapolation and
! the solar neutrino fluxes.
      integer :: max_iter_level4 = 0
      logical :: improved_first_guess_flag = .false., calc_neutrinos = .false.

! Light-element rate scales: the reader divides the decks' xsli6 etc by
! their reference values (sli6 = xsli6/5.5d3).
      real(dp) :: li6_rate_scale = 1.0d0
      real(dp) :: li7_rate_scale = 1.0d0
      real(dp) :: be9_pg_rate_scale = 1.0d0
      real(dp) :: be9_pd_rate_scale = 1.0d0
      real(dp) :: be9_palpha_rate_scale = 1.0d0

! Itoh et al. 1996 neutrino losses (net_lib's compute_zone_energy).
      logical :: use_itoh_neutrino_loss = .false.

! Whether to run the full envelope integration.
      logical :: calc_envelope_flag = .true.

! Extended (15-species) or default composition array.
      logical :: use_extended_composition = .false.

! Starspots.
      real(dp) :: spot_filling_factor = 0.00d0
      real(dp) :: spot_temp_contrast = 1.00d0
      logical :: spot_depth_varies = .false.

! Disk locking.
      real(dp) :: disk_locking_age_gyr = 0.0d0
      real(dp) :: disk_omega_rad_s = 7.2722d-6
      logical :: disk_locking_active = .false.

! Per-run stopping criteria and timestep override; the two flags are
! derived from them by the reader.
      real(dp) :: target_end_age(max_runs) = 0.0d0
      real(dp) :: timestep_override(max_runs) = 0.0d0
      real(dp) :: central_deuterium_stop(max_runs) = 0.0d0
      real(dp) :: central_hydrogen_stop(max_runs) = 0.0d0
      real(dp) :: central_helium_stop(max_runs) = 0.0d0
      logical :: end_age_stop_active(max_runs), timestep_override_active(max_runs)

! The MHD EOS, and the units of the ZAMS and central tables (set by the
! reader).
      logical :: use_mhd_eos = .false.
      integer :: unit_zams_a, unit_zams_b, unit_zams_c
      integer :: unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

! Opacity-table metallicity matching.
      real(dp) :: metal_fraction_match_tolerance = 1.0d-8

! Envelope profile print count.
      integer :: npenv = 0


! Unit numbers, set by the reader. opal92_table_unit is the unit
! kap/opal92/read_opal92_tables.f90 opens the OPAL92 table on.
      integer :: first_unit, run_unit, standard_unit, fermi_unit
      integer :: dynamics_unit, opal92_table_unit
      integer :: composition_unit, kurucz_table_unit

! Monte Carlo runs.
      logical :: monte_carlo_active = .false.
      integer :: mc_run_start = 1, mc_run_end = 1

! Newton iteration: tolerances, corrections, zoning scales and iteration
! limits. fcorr has no input; the run's copy, star%job%fcorr, is set by
! the Newton iteration.
      real(dp) :: htoler(5,2) = reshape((/6.0d-5,4.5d-5,3.0d-5, &
           9.0d-5,3.0d-5,9.0d-1,5.0d-1,5.0d-1,2.0d0,2.5d-6/), (/5,2/))
      real(dp) :: fcorr0 = 0.8d0, fcorri = 0.1d0, fcorr
      real(dp) :: chi_grid_scale(12) = (/ 1.0d-8, 8.0d-2, 5.0d-2, &
           5.0d-2, 1.0d0, 1.0d0, 0.0d0, 5.0d-2, 2.0d-2, 5.0d-2, 5.0d-2, &
           1.0d-1 /)
      integer :: max_iter_level1 = 2, max_iter_level2 = 20, max_iter_level3 = 2

! Angular momentum transport: the sub-timestep factor, iteration limits
! and convergence tolerance.
      real(dp) :: diffusion_timestep_factor = 1.0d-2
      real(dp) :: convergence_tolerance = 1.0d-4
      integer :: num_rotation_structure_iters = 1, max_diffusion_iters = 10

! Gravitational settling.
      real(dp) :: settling_tolerance = 1.0d-8
      integer :: coulomb_log_choice = 1
      integer :: settling_num_iterations = 10
      logical :: diffuse_helium_active = .false.
      real(dp) :: settling_timestep_fraction = 0.1d0
      real(dp) :: hydrogen_diffusion_floor = 1.0d-3
      real(dp) :: helium_diffusion_min = 1.0d-3
      logical :: use_thoul_fit = .false.
      real(dp) :: fgry = 1.0d0, fgrz = 1.0d0
      logical :: use_thoul_diffusion = .false.
      logical :: use_diffusion_z = .false.
      logical :: use_new_diffusion_routines = .false.
      logical :: diffuse_lithium = .false.

! Atmosphere and envelope integration steps and tolerances.
      real(dp) :: atm_error_tol = 3.0d-4
      real(dp) :: atm_step_initial = 1.0d-10
      real(dp) :: atm_step_begin = 1.0d-1
      real(dp) :: atm_step_min = 1.0d-1
      real(dp) :: atm_step_max = 5.0d-1
      real(dp) :: env_error_tol = 3.0d-4
      real(dp) :: env_step_begin = 1.0d-1
      real(dp) :: env_step_min = 1.0d-1
      real(dp) :: env_step_max = 5.0d-1
! The Bulirsch-Stoer integrator of both: the tolerance on the fitting
! mass, as a fraction of |log_fit_point_mass_fraction|, which is also
! the atmosphere's step tolerance; the number of stages; and the
! extrapolation order (numerics/bulirsch_stoer.f90 takes 1 to 11 stages
! and orders 2 to 7; the controls read checks them).
      real(dp) :: tolerance_fraction = 1.0d-3
      integer :: max_stage_index = 11
      integer :: extrap_order = 7

! The atmosphere choice as read. The run's working copy is
! star%job%atm_choice (the surface boundary condition and the
! calibration cycles change it); star%atm_choice_initial keeps this value.
      integer :: atm_choice = 0

! Initial abundances of a 12-species set (Na/Al/Mg/Fe/Si/C/H/O/N/Ar/Ne/He,
! a G&N93 solar mixture).
      real(dp) :: mixture_weights_seed(12) = (/0.001999d0, 0.003238d0, &
           0.037573d0, 0.071794d0, 0.040520d0, 0.173285d0, 0.000000d0, &
           0.482273d0, 0.053152d0, 0.005379d0, 0.098668d0, 0.000000d0/)

! Monte Carlo input files and their units (set by the reader). The paths
! have no default.
      character(len=256) :: monte_carlo_file1_path = repeat(achar(0), 256)
      character(len=256) :: monte_carlo_file2_path = repeat(achar(0), 256)
      integer :: monte_carlo_unit1, monte_carlo_unit2

! GYRE pulsation-output cadence.
      integer :: pulse_gyre_interval = 0

! The envelope-fit triangle (its limits; lnew0 forces a fresh one at run
! start) and the requested envelope mass (set from the starting model).
      real(dp) :: tri_delta_teffl = 1.0d-2
      real(dp) :: tri_delta_logl = 8.0d-2
      logical :: lnew0 = .false.
      real(dp) :: requested_envelope_mass
      logical :: change_envelope_mass_flag

! The run list: per-run model counts, first-call flags and rescaling
! (rescale_params and rescale_kind are derived from rsclm/rsclx/rsclz/
! rsclcm and kindrn by the reader).
      real(dp) :: rescale_params(4,max_runs)
      integer :: rescale_kind(max_runs)
      integer :: num_models(max_runs) = 0
      integer :: num_runs = 1
      logical :: first_call_flag(max_runs) = .true.

! Timestep limits from rotation and structure changes.
      real(dp) :: max_domega_global = 1.0d1
      logical :: use_structure_dt_limits = .true.

! Envelope and atmosphere generation step sizes.
      real(dp) :: atm_step_size = 0.5d0
      real(dp) :: envelope_step_size = 0.5d0
      logical :: envelope_generation_flag = .false.

! Helium flash.
      logical :: helium_flash_active = .false.

! The starting envelope composition, used by the reader.
      real(dp) :: initial_envelope_x = 0.7d0
      real(dp) :: initial_envelope_z = 0.02d0

! Rescaling one species of the starting model.
      real(dp) :: new_species_value = 1.3d1
      integer :: new_species_index
      logical :: rescale_species_active = .false.
      logical :: value_relative_to_h

! Changing the starting model's CNO mixture and isotope ratios. The two
! flags are set by the reader's validation of these inputs.
      integer :: mixture_change_mode = 0
      integer :: isotope_change_mode = 0
      logical :: change_cno_mixture_active, change_isotope_ratios_active
      real(dp) :: target_carbon_cno_fraction = 0.172148d0
      real(dp) :: target_nitrogen_cno_fraction = 0.050426d0
      real(dp) :: target_oxygen_cno_fraction = 0.468195d0
      real(dp) :: c12_to_c13_ratio = 88.27d0
      real(dp) :: n14_to_n15_ratio = 411.9d0
      real(dp) :: o16_to_o17_ratio = 471.4d0
      real(dp) :: o16_to_o18_ratio = 2693.0d0
      real(dp) :: target_metal_fraction = 0.02292d0
      real(dp) :: initial_h2_fraction = 2.781d-5
      real(dp) :: initial_he3_fraction = 3.461d-5
      real(dp) :: initial_li6_fraction = 7.187d-10
      real(dp) :: initial_li7_fraction = 1.025d-8
      real(dp) :: initial_be9_fraction = 1.595d-10
      real(dp) :: initial_b10_fraction = 1.002d-9
      real(dp) :: initial_b11_fraction = 4.405d-9

! Rotational mixing velocity and coefficient scales (the legacy ies, imu
! and fmu are read and then unused).
      real(dp) :: fo = 1.0d0
      real(dp) :: difad_velocity_scale = 1.0d0
      real(dp) :: mixing_velocity_scale = 1.0d0
      real(dp) :: es_velocity_scale = 1.0d0
      real(dp) :: gsf_velocity_scale = 1.0d0
      real(dp) :: secular_shear_velocity_scale = 1.0d0
      real(dp) :: critical_reynolds = 1.0d3
      integer :: ies = 1
      integer :: imu = 1
      real(dp) :: es_mixing_scale = 1.0d0
      real(dp) :: secular_shear_mixing_scale = 1.0d0
      real(dp) :: gsf_mixing_scale = 1.0d0
! The Solberg-Hoiland instability's velocity scale (0: off) and mixing
! efficiency; new in 2026, with no legacy names.
      real(dp) :: solberg_hoiland_velocity_scale = 0.0d0
      real(dp) :: solberg_hoiland_mixing_scale = 1.0d0
      integer :: gsf_inhibition_mode = 1

! The burning network's integration tolerances.
      real(dp) :: min_abundance = 1.0d-20
      real(dp) :: absolute_tolerance = 1.0d-5
      real(dp) :: relative_tolerance = 1.0d-4
      integer :: max_burn_iterations = 50

! OPAL95 table unit (set by the reader).
      integer :: opal95_table_unit


! Extending the core inward.
      logical :: extend_core_inward = .false.
      integer :: num_core_shells_added = 0
      real(dp) :: core_mass_reduction_factor = 1.0d0

! Isochrone output and the short file.
      logical :: rewind_short_file = .false.
      logical :: isochrone_output_active = .false.
      integer :: isochrone_file_unit

! Per-run starting composition, mixing length and envelope mass.
      real(dp) :: initial_x_array(max_runs) = 0
      real(dp) :: initial_z_array(max_runs) = 0
      real(dp) :: mixing_length_array(max_runs) = 0
      real(dp) :: senv0_array(max_runs) = 1.26d-4
      logical :: has_senv0_array(max_runs) = .false.

! Solar calibration.
      real(dp) :: luminosity_tolerance = 1.0d-5
      real(dp) :: radius_tolerance = 1.0d-4
      real(dp) :: zx_tolerance = 1.0d-3
      real(dp) :: target_solar_zx = 0.02292d0
      real(dp) :: target_solar_age = 4.57d9
      logical :: calibrate_solar_model = .false.
      logical :: calibrate_solar_zx = .false.

! The metallicity ramp, per run.
      real(dp) :: rsclzc(max_runs) = -1.0d0, rsclzm1(max_runs) = -1.0d0, &
           rsclzm2(max_runs) = -1.0d0
      logical :: use_z_ramp = .false.

! Calibrating a star to a target luminosity, Teff or radius.
      real(dp) :: target_luminosity_lsun = 0
      real(dp) :: target_star_luminosity_tolerance = 0
      real(dp) :: target_teff = 0
      real(dp) :: target_radius_rsun = 0
      logical :: specify_teff_flag = .false.
      logical :: calibrate_star_flag = .false.

! The OPAL EOS tables and numerical EOS derivatives (iopale is a unit).
      logical :: use_opal95_eos = .false.
      logical :: use_opal2001_eos = .false.
      logical :: use_opal2006_eos = .false.
      logical :: use_numerical_derivatives = .false.
      integer :: iopale
! the SCV low-temperature EOS (namelist lscv), and the Debye-Hueckel
! correction with its eta limits (namelist ldh, etadh0, etadh1)
      logical :: use_scv_eos = .false.
      logical :: use_debye_huckel_correction = .false.
      real(dp) :: debye_huckel_eta_min = -1.0d0
      real(dp) :: debye_huckel_eta_max = 1.0d0

! Opacity tables: their metallicities, the molecular-opacity range and
! which tables to use.
      real(dp) :: laol_table_z1 = 0
      real(dp) :: laol_table_z2 = 0
      real(dp) :: opal_table_z1 = 0
      real(dp) :: opal_table_z2 = 0
      real(dp) :: opal95_single_table_z = 0
      real(dp) :: alex_table_z1 = 0
      real(dp) :: kurucz_table_z1 = 0
      real(dp) :: kurucz_table_z2 = 0
      real(dp) :: molecular_opacity_logt_min = 4.0d0
      real(dp) :: molecular_opacity_logt_max = 4.1d0
      logical :: use_alex06_tables = .false.
      logical :: use_laol89_tables = .false.
      logical :: use_opal92_tables = .false.
      logical :: use_opal95_tables = .false.
      logical :: use_kurucz90_tables = .false.
      logical :: use_alex95_tables = .false.

! The LAOL89 interpolators extrapolate linearly within log(tollaol) of a
! table edge.
      real(dp) :: tollaol = 10.0d0
      logical :: use_pure_z_table = .false.

! Conductive opacities, and the second Kurucz table's unit (icondopacp is
! a unit too).
      integer :: kurucz_table2_unit, icondopacp
      logical :: use_conductive_opacity = .false.

! Alexander table units (set by the reader). alex94_table_unit reads the
! Alexander 1994 low-T tables (kap/alex94/alex94_opacity.f90); the
! namelist names lalex95, use_alex95_tables and alex95_table_z_list keep
! their historical "95".
      integer :: alex94_table_unit
      integer :: alex06_table_unit

! Angular momentum transport options: the alternative velocity form in
! the diffusion coefficients (lvfc) and diffusion-advection transport.
      logical :: lvfc = .false.
      logical :: use_diffusion_advection_transport = .false.

! No angular momentum transport in the core.
      logical :: no_am_transport_in_core = .false.

! SCV EOS table units (set by the reader).
      integer :: scv_h_unit, scv_he_unit, scv_z_unit

! Allard atmosphere tables: the target [Fe/H] and [alpha/Fe], tau=100
! tables, and the unit (set by the reader).
      real(dp) :: allard_target_feh = 0d0
      real(dp) :: allard_target_alpha = 0d0
      logical :: allard_use_tau100 = .false.
      integer :: allard_table_unit

! Solid-body rotation.
      logical :: force_solid_body_rotation = .false.
      integer :: solid_body_mode_flag = 0

! Constant mixing and settling coefficients.
      real(dp) :: constant_mixing_coeff = 1.0d0, constant_settling_reduction = 1.0d0

! The timestep limit from the envelope-fit triangle.
      logical :: use_envelope_triangle_dt = .false.

! Magnetic braking: the Matt/PMM wind law (K97 and V13 set their own
! exponents; rotation/magnetic_braking.f90's derive_wind_law).
      real(dp) :: pmm_exponent_a = 2.0d0
      real(dp) :: pmm_exponent_b = 1.0d0
      real(dp) :: pmm_exponent_c = 2.0d0
      real(dp) :: pmm_exponent_d = 0.0d0
      real(dp) :: pmm_exponent_m = 0.5d0
      real(dp) :: pmm_norm_jdot = 1.32d30
      real(dp) :: pmm_norm_mdot = 1.27d12
      real(dp) :: pmm_solar_pressure = 4.9304d0
      real(dp) :: pmm_solar_omega = 2.836d-6
      real(dp) :: pmm_solar_turnover_timescale = 9.4805d5
      logical :: use_pmm_wind_law = .false.
      logical :: scale_by_rossby_number = .true.
      logical :: scale_by_b_field = .false.
      character(len=3) :: wind_law_name = 'K97'

! Magnetic braking: the Kawaler (1988) law's index alfa and the constant
! factor fk of both laws; the saturation omega (rad/s), or, below 1000,
! the same in units of the solar omega.
      real(dp) :: wind_law_index = 1.5d0
      real(dp) :: wind_law_constant = 1.0d0
      real(dp) :: wind_saturation_omega = 3.0d-4
      real(dp) :: solar_omega_sat = 1000.0d0
      logical :: use_wind_torque = .true.

! Constant background diffusion.
      real(dp) :: constant_background_diffusion_coeff = 2.5d4
      logical :: use_constant_background_diffusion = .false.

! The timestep and zoning tolerances under readable names; when lnewvars
! is set, map_user_inputs fills atime and chi_grid_scale from them.
      real(dp) :: flag_dx = 0.05d0, flag_dw = 0.10d0, flag_dz = 0.05d0
      real(dp) :: time_core_min = 1.0d-3, time_dl = 2.0d-2, &
           time_dp = 4.0d-2, time_dr = 2.0d-2, time_dt = 2.0d-2, &
           time_dw_global = 8.0d-2, time_dw_mix = 8.0d-2, &
           time_dx_core_frac = 0.5d0, time_dx_core_tot = 2.0d-2, &
           time_dx_shell = 0.1d0, time_dx_total = 1.5d-3, &
           time_dy_core_frac = 0.5d0, time_dy_core_tot = 2.0d-2, &
           time_dy_shell = 0.1d0, time_dy_total = 1.5d-3, &
           time_max_dt_frac = 1.5d0
      real(dp) :: tol_czbase_fine_width = 0.0d0, tol_dl_max = 0.02d0, &
           tol_dm_max = 0.08d0, tol_dm_min = 1.0d-8, &
           tol_dp_core_max = 0.05d0, tol_dp_czbase_max = 0.05d0, &
           tol_dp_env_max = 0.05d0, tol_dx_max = 1.0d0, tol_dz_max = 1.0d0
      logical :: lstruct_time = .false., lnewvars = .false.

! The nuclear S-factors (keV b) and their first and second derivatives
! (b, b/keV), used with use_new_nuclear_rates: Solar Fusion II
! (Adelberger et al. 2011). pep is relative to pp, and Be7 + e- is the
! constant of the linear term; the code uses T9, so a coefficient for
! T/1e6 K is divided by sqrt(1000). map_user_inputs scales them against
! the Bahcall & Ulrich (1988) S0 and the Solar Fusion I derivatives.
      real(dp) :: s0_pp = 4.01d-22
      real(dp) :: s0_he3he3 = 5.21d3
      real(dp) :: s0_he3he4 = 5.6d-1
      real(dp) :: s0_p_c12 = 1.34d0
      real(dp) :: s0_p_c13 = 7.6d0
      real(dp) :: s0_p_n14 = 1.66d0
      real(dp) :: s0_p_o16 = 1.06d1
      real(dp) :: s0_be7_electron = 1.7709d-10
      real(dp) :: s0_be7_p = 0.0208d0
      real(dp) :: s0_n15_p_c12_branch = 7.3d4
      real(dp) :: s0_n15_p_o16_branch = 3.6d1
      real(dp) :: s0p_pp = 4.49d-24
      real(dp) :: s0p_he3he3 = -4.9d0
      real(dp) :: s0p_he3he4 = -3.6d-4
      real(dp) :: s0p_p_c12 = 2.6d-3
      real(dp) :: s0p_p_c13 = -7.83d-3
      real(dp) :: s0p_p_n14 = -3.3d-3
      real(dp) :: s0p_p_o16 = -5.4d-2
      real(dp) :: s0pp_p_c12 = 8.3d-5
      real(dp) :: s0pp_p_c13 = 7.29d-4
      real(dp) :: s0pp_p_o16 = 0.0d0
      real(dp) :: s0p_be7_p = -3.12d-5
      real(dp) :: s0pp_be7_p = -2.288d-7
      real(dp) :: s0_pep = 3.5734d-6, s0_hep = 8.6d-20

end module controls_buffer

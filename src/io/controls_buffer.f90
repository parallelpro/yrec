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
! generates star%ctrl's type and the two copies (state/controls_sync.f90)
! from them, so regenerate after changing a member, and add the control
! to defaults/controls_registry.tsv, which maps it to its legacy and
! new-inlist names. Each control's default is stated once in the code:
! here, or, for a control read under a legacy name with no default here,
! in set_legacy_defaults (io/legacy_controls.f90). The legacy namelists
! are frozen: a new control is read only by the new format. Only the
! controls reader, map_user_inputs,
! controls_sync and test_net use this module; any file may import its
! named constants (max_runs, itime_*, ichi_*). tools/check_boundaries.py
! enforces both.

module controls_buffer
      implicit none

! Maximum number of runs (kind cards) in one job: the dimension of
! every per-run array below (the calibration protocols in
! setup/setup_star_calibration.f90 / setup_solar_calibration.f90 fill
! them up to this limit). tools/gen_controls_state.py substitutes the
! literal value into the generated star%ctrl / star%job type bodies,
! which star_info_lib includes without seeing this module.
      integer, parameter, public :: max_runs = 50

! MESA-style output (new-format inlists only). The run-log, profile and
! pulsation-output switches below apply to both output modes.
      character(len=256), public :: star_history_name = 'history.data'
      character(len=256), public :: history_columns_file = ' '
      character(len=256), public :: profile_columns_file = ' '
      integer, public :: profile_interval = 50
! Structure-limit stop conditions, per run like the other stopping
! criteria (target_end_age, central_*_stop); a limit at its +-1d99
! sentinel is disabled.
      double precision, public :: log_L_upper_limit(max_runs) = 1.0d99
      double precision, public :: log_L_lower_limit(max_runs) = -1.0d99
      double precision, public :: Teff_upper_limit(max_runs) = 1.0d99
      double precision, public :: Teff_lower_limit(max_runs) = -1.0d99
      double precision, public :: log_g_upper_limit(max_runs) = 1.0d99
      double precision, public :: log_g_lower_limit(max_runs) = -1.0d99
      double precision, public :: nu_max_upper_limit(max_runs) = 1.0d99
      double precision, public :: nu_max_lower_limit(max_runs) = -1.0d99
! The solar mass and Newton's G. G_cgs < 0 (the default) keeps the
! historical log10 G = -7.17571 (G = 6.6726e-8) bit-for-bit; set a
! positive value (e.g. 6.6743d-8, CODATA 2018) to override.
      double precision, public :: solar_mass_cgs = 1.9891d33
      double precision, public :: G_cgs = -1.0d0
! solar reference values for the asteroseismic scaling relations:
      double precision, public :: nu_max_sun = 3090.0d0
      double precision, public :: delta_nu_sun = 135.1d0
      double precision, public :: Teff_sun = 5772.0d0
      character(len=256), public :: profile_data_prefix = 'profile'
      character(len=256), public :: inlist_used_file = '{YREC_OUTPUT}/inlist_used'
      character(len=8), public :: pulse_format = 'GYRE'
! Independent toggles: profiles only, pulse only, both, or neither
! (cadence stays with profile_interval / pulse_gyre_interval):
      logical, public :: write_profile_flag = .false.
      logical, public :: write_pulse_flag = .false.
! Run-log verbosity: the compact per-model progress line prints every
! terminal_interval models (0 = off; the final model always prints),
! and report_solver_diagnostics restores the verbose solver forensics
! (Henyey iteration trace, envelope/atmosphere integration statistics,
! triangle/rezoning bookkeeping).
      integer, public :: terminal_interval = 50
! model numbering for fresh start-model loads: the first converged
! model after a first_call_flag load is numbered this (default 1);
! <= 0 keeps the counter stored in the model file (restart
! continuation, the pre-2026 behavior).
      integer, public :: set_initial_model_number = 1
      logical, public :: report_solver_diagnostics = .false.
! diagnostic: per-zone solar neutrino production table (compute_zone_energy per
! shell) written to the short/log stream at the start of each run
      logical, public :: compute_neutrino_fluxes = .false.

! Timestep factors (atime) and log T cutoffs: nuclear burning, Saha
! ionization, the envelope's temperature window (tenv0, tenv1; star%tenv
! is their mean) and the simplified temperature gradient (tgcut).
      double precision :: atime(14) = (/1.0d-3,2.0d-2,5.0d-1,2.0d-2, &
           3.0d-1,1.5d-3,1.0d-1,2.0d-2,4.0d-2,2.0d-2,2.0d-2,0.25d0, &
           1.5d0,0.25d0/)
! Slot names for atime (the per-phase timestep-control factors). The
! meaning of each slot is fixed by setup/map_user_inputs.f90, which
! fills it from the like-named time_* namelist inputs when lnewvars
! is set; readers subscript with these names.
      integer, parameter, public :: &
           itime_core_min = 1, &      ! time_core_min: min central abundance for the core-burning criterion
           itime_dx_core_tot = 2, &   ! time_dx_core_tot: max absolute delta-X burned in core per step
           itime_dx_core_frac = 3, &  ! time_dx_core_frac: max fraction of core X burned per step
           itime_dy_core_tot = 4, &   ! time_dy_core_tot: max absolute delta-Y burned in core per step
           itime_dy_core_frac = 5, &  ! time_dy_core_frac: max fraction of core Y burned per step
           itime_dx_total = 6, &      ! time_dx_total: max Msun of X burned per step
           itime_dx_shell = 7, &      ! time_dx_shell: max fraction of X burned in the shell per step
           itime_dt = 8, &            ! time_dt: max change in T at any shell per step
           itime_dp = 9, &            ! time_dp: max change in P at any shell per step
           itime_dr = 10, &           ! time_dr: max change in R at any shell per step
           itime_dl = 11, &           ! time_dl: max change in L at any shell per step
           itime_dy_shell = 12, &     ! time_dy_shell: max fraction of Y burned in the shell per step
           itime_max_dt_frac = 13, &  ! time_max_dt_frac: max timestep change factor model to model
           itime_dy_total = 14        ! time_dy_total: max Msun of Y burned in the shell per step
      double precision :: nuclear_logT_cutoffs(5) = (/6.5d0,6.5d0,6.82d0,7.7d0,7.5d0/)
      double precision :: saha_log10t_cutoff = 6.0d0
      double precision :: tenv0 = 3.0d0, tenv1 = 9.0d0, tgcut = 6.9d0

! Nuclear rate selection and weak screening.
      logical :: use_new_nuclear_rates
      double precision :: weak_screening_threshold

! Overshoot and semiconvection. iov1/iov2, the overshoot range, are
! model state (star%iov1/iov2).
      double precision :: dpenv = 1.0d0
      logical :: core_overshoot_active = .false.
      double precision :: overshoot_alpha_core = 0.0d0
      logical :: envelope_overshoot_active = .false.
      double precision :: overshoot_alpha_envelope = 0.0d0
      logical :: lovstm = .false.
      double precision :: alpham = 0.0d0
      logical :: ladov = .false.
      logical :: lovmax = .false.
      double precision :: betac = 0.15d0
      logical :: use_semiconvection = .false.

! Rotation. The reader clamps walpcz to [-2, 0].
      double precision :: acfpft = 1.0d-36
      integer :: itfp1 = 5, itfp2 = 20
      logical :: rotation_active = .false.
      double precision :: walpcz = 0.0d0
      logical :: instability_transport_active = .false., set_initial_omega = .false.
      double precision :: initial_omega = 0.0d0

! Mass accretion and the Reimers wind. core/mass_change switches
! star%job%use_mass_accretion off when the disk is exhausted.
      double precision :: mass_accretion_rate = -1.0d-14
      double precision :: fczdmdt = 0.1d0, ftotdmdt = 1.0d-2
      double precision :: accreted_composition(15) = (/0.71668d0, &
           0.265721d0,0.01757d0,2.9d-5,3.013d-3,3.385d-5,9.346d-4, &
           0.0d0,8.462d-3,0.0d0,1.696d-5,0.0d0,2.0d-9,2.0d-9,3.0d-11/)
      double precision :: reimers_scaling_factor = -4.0d-13
      logical :: use_reimers_wind = .false., use_mass_accretion = .false.

! The level-4 Henyey iteration limit, the first-guess extrapolation and
! the solar neutrino fluxes.
      integer :: max_iter_level4 = 0
      logical :: improved_first_guess_flag = .false., calc_neutrinos = .false.

! Light-element rate scales: the reader divides the decks' xsli6 etc by
! their reference values (sli6 = xsli6/5.5d3).
      double precision :: li6_rate_scale, li7_rate_scale, &
           be9_pg_rate_scale, be9_pd_rate_scale, be9_palpha_rate_scale

! Itoh et al. 1996 neutrino losses (net_lib's compute_zone_energy).
      logical :: use_itoh_neutrino_loss = .false.

! Whether to run the full envelope integration.
      logical :: calc_envelope_flag = .true.

! Extended (15-species) or default composition array.
      logical :: use_extended_composition = .false.

! Starspots.
      double precision :: spot_filling_factor = 0.00d0
      double precision :: spot_temp_contrast = 1.00d0
      logical :: spot_depth_varies = .false.

! Disk locking.
      double precision :: disk_locking_age_gyr = 0.0d0
      double precision :: disk_omega_rad_s = 7.2722d-6
      logical :: disk_locking_active = .false.

! Per-run stopping criteria and timestep override; the two flags are
! derived from them by the reader.
      double precision :: target_end_age(max_runs) = 0.0d0
      double precision :: timestep_override(max_runs) = 0.0d0
      double precision :: central_deuterium_stop(max_runs) = 0.0d0
      double precision :: central_hydrogen_stop(max_runs) = 0.0d0
      double precision :: central_helium_stop(max_runs) = 0.0d0
      logical :: end_age_stop_active(max_runs), timestep_override_active(max_runs)

! The MHD EOS, and the units of the ZAMS and central tables (set by the
! reader).
      logical :: use_mhd_eos = .false.
      integer :: unit_zams_a, unit_zams_b, unit_zams_c
      integer :: unit_centre1, unit_centre2, unit_centre3, unit_centre4, &
           unit_centre5

! Opacity-table metallicity matching.
      double precision :: metal_fraction_match_tolerance
      double precision :: zsi = 0.0d0

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
      double precision :: htoler(5,2) = reshape((/6.0d-5,4.5d-5,3.0d-5, &
           9.0d-5,3.0d-5,9.0d-1,5.0d-1,5.0d-1,2.0d0,2.5d-6/), (/5,2/))
      double precision :: fcorr0 = 0.8d0, fcorri = 0.1d0, fcorr
      double precision :: chi_grid_scale(12)
! Slot names for chi_grid_scale (the former hpttol vector). The
! meaning of each slot is fixed by setup/map_user_inputs.f90, which
! fills it from the like-named tol_*/flag_* namelist inputs when
! lnewvars is set; readers subscript with these names.
      integer, parameter, public :: &
           ichi_dm_min = 1, &            ! tol_dm_min: min log-m spacing
           ichi_dm_max = 2, &            ! tol_dm_max: max log-m spacing
           ichi_flag_dx = 3, &           ! flag_dx: X jump that flags a shell
           ichi_flag_dz = 4, &           ! flag_dz: Z jump that flags a shell
           ichi_dx_max = 5, &            ! tol_dx_max: max X spacing
           ichi_dz_max = 6, &            ! tol_dz_max: max Z spacing
           ichi_czbase_fine_width = 7, & ! tol_czbase_fine_width: log-P width of fine zoning at CZ base
           ichi_dp_env_max = 8, &        ! tol_dp_env_max: max log-P spacing, convective envelope
           ichi_dl_max = 9, &            ! tol_dl_max: max L/Lsurf spacing
           ichi_dp_czbase_max = 10, &    ! tol_dp_czbase_max: max log-P spacing near CZ base
           ichi_dp_core_max = 11, &      ! tol_dp_core_max: max log-P spacing, radiative interior
           ichi_flag_dw = 12             ! flag_dw: max log-omega spacing (rotating models)
      integer :: max_iter_level1 = 2, max_iter_level2 = 20, max_iter_level3 = 2

! Angular momentum transport: the sub-timestep factor, iteration limits
! and convergence tolerance.
      double precision :: diffusion_timestep_factor = 1.0d-2
      double precision :: convergence_tolerance
      integer :: num_rotation_structure_iters = 1, max_diffusion_iters = 1

! Gravitational settling.
      double precision :: settling_tolerance
      integer :: coulomb_log_choice, settling_num_iterations
      logical :: diffuse_helium_active
      double precision :: settling_timestep_fraction, &
           hydrogen_diffusion_floor, helium_diffusion_min
      logical :: use_thoul_fit
      double precision :: fgry = 1.0d0, fgrz = 1.0d0
      logical :: use_thoul_diffusion = .false.
      logical :: use_diffusion_z
      logical :: use_new_diffusion_routines
      logical :: diffuse_lithium = .false.

! Atmosphere and envelope integration steps and tolerances.
      double precision :: atm_error_tol, atm_step_initial, atm_step_begin
      double precision :: atm_step_min, atm_step_max
      double precision :: env_error_tol, env_step_begin, env_step_min
      double precision :: env_step_max

! The atmosphere choice as read. The run's working copy is
! star%job%atm_choice (the surface boundary condition and the
! calibration cycles change it); star%atm_choice_initial keeps this value.
      integer :: atm_choice

! Initial abundances of a 12-species set (Na/Al/Mg/Fe/Si/C/H/O/N/Ar/Ne/He,
! a G&N93 solar mixture).
      double precision :: mixture_weights_seed(12) = (/0.001999d0, 0.003238d0, &
           0.037573d0, 0.071794d0, 0.040520d0, 0.173285d0, 0.000000d0, &
           0.482273d0, 0.053152d0, 0.005379d0, 0.098668d0, 0.000000d0/)

! Monte Carlo input files and their units (set by the reader). The paths
! have no default.
      character(len=256) :: monte_carlo_file1_path, monte_carlo_file2_path
      integer :: monte_carlo_unit1, monte_carlo_unit2

! GYRE pulsation-output cadence.
      integer :: pulse_gyre_interval = 0

! The envelope-fit triangle (its limits; lnew0 forces a fresh one at run
! start) and the requested envelope mass (set from the starting model).
      double precision :: tri_delta_teffl, tri_delta_logl
      logical :: lnew0 = .false.
      double precision :: requested_envelope_mass
      logical :: change_envelope_mass_flag

! The run list: per-run model counts, first-call flags and rescaling
! (rescale_params and rescale_kind are derived from rsclm/rsclx/rsclz/
! rsclcm and kindrn by the reader).
      double precision :: rescale_params(4,max_runs)
      integer :: rescale_kind(max_runs)
      integer :: num_models(max_runs), num_runs
      logical :: first_call_flag(max_runs)

! Timestep limits from rotation and structure changes.
      double precision :: max_domega_global
      logical :: use_structure_dt_limits

! Envelope and atmosphere generation step sizes.
      double precision :: atm_step_size, envelope_step_size
      logical :: envelope_generation_flag

! Helium flash.
      logical :: helium_flash_active

! The starting envelope composition, used by the reader.
      double precision :: initial_envelope_x = 0.7d0
      double precision :: initial_envelope_z = 0.02d0

! Rescaling one species of the starting model.
      double precision :: new_species_value
      integer :: new_species_index
      logical :: rescale_species_active
      logical :: value_relative_to_h

! Changing the starting model's CNO mixture and isotope ratios. The two
! flags are set by the reader's validation of these inputs.
      integer :: mixture_change_mode, isotope_change_mode
      logical :: change_cno_mixture_active, change_isotope_ratios_active
      double precision :: target_carbon_cno_fraction, &
           target_nitrogen_cno_fraction, target_oxygen_cno_fraction
      double precision :: c12_to_c13_ratio, n14_to_n15_ratio, &
           o16_to_o17_ratio, o16_to_o18_ratio
      double precision :: target_metal_fraction
      double precision :: initial_h2_fraction, initial_he3_fraction, &
           initial_li6_fraction, initial_li7_fraction, &
           initial_be9_fraction, initial_b10_fraction, &
           initial_b11_fraction

! Rotational mixing velocity and coefficient scales (ies and imu are
! read and then unused).
      double precision :: fo = 1.0d0
      double precision :: difad_velocity_scale, mixing_velocity_scale, &
           es_velocity_scale, gsf_velocity_scale, mu_gradient_scale, &
           secular_shear_velocity_scale, critical_reynolds
      integer :: ies = 1
      integer :: imu = 1
      double precision :: es_mixing_scale, secular_shear_mixing_scale, &
           gsf_mixing_scale
      integer :: gsf_inhibition_mode

! The burning network's integration tolerances.
      double precision :: min_abundance, absolute_tolerance, &
           relative_tolerance
      integer :: max_burn_iterations

! OPAL95 table unit (set by the reader).
      integer :: opal95_table_unit


! Extending the core inward.
      logical :: extend_core_inward
      integer :: num_core_shells_added
      double precision :: core_mass_reduction_factor

! Isochrone output and the short file.
      logical :: rewind_short_file, isochrone_output_active
      integer :: isochrone_file_unit

! Per-run starting composition, mixing length and envelope mass.
      double precision :: initial_x_array(max_runs), initial_z_array(max_runs), &
           mixing_length_array(max_runs), senv0_array(max_runs)
      logical :: has_senv0_array(max_runs)

! Solar calibration.
      double precision :: luminosity_tolerance, radius_tolerance, &
           zx_tolerance, target_solar_zx, target_solar_age
      logical :: calibrate_solar_model, calibrate_solar_zx

! The metallicity ramp, per run.
      double precision :: rsclzc(max_runs) = -1.0d0, rsclzm1(max_runs) = -1.0d0, &
           rsclzm2(max_runs) = -1.0d0
      logical :: use_z_ramp

! Calibrating a star to a target luminosity, Teff or radius.
      double precision :: target_luminosity_lsun, target_star_luminosity_tolerance, target_teff, target_radius_rsun
      logical :: specify_teff_flag, calibrate_star_flag

! The OPAL EOS tables and numerical EOS derivatives (iopale is a unit).
      logical :: use_opal95_eos, use_opal2001_eos, use_opal2006_eos, &
           use_numerical_derivatives
      integer :: iopale
! the SCV low-temperature EOS (namelist lscv), and the Debye-Hueckel
! correction with its eta limits (namelist ldh, etadh0, etadh1)
      logical :: use_scv_eos, use_debye_huckel_correction
      double precision :: debye_huckel_eta_min, debye_huckel_eta_max

! Opacity tables: their metallicities, the molecular-opacity range and
! which tables to use.
      double precision :: laol_table_z1, laol_table_z2, opal_table_z1, &
           opal_table_z2, opal95_single_table_z, alex_table_z1, &
           kurucz_table_z1, kurucz_table_z2, &
           molecular_opacity_logt_min, molecular_opacity_logt_max
      logical :: use_alex06_tables, use_laol89_tables, use_opal92_tables, use_opal95_tables, use_kurucz90_tables, use_alex95_tables

! The LAOL89 interpolators extrapolate linearly within log(tollaol) of a
! table edge.
      double precision :: tollaol = 10.0d0
      logical :: use_pure_z_table = .false.

! Conductive opacities, and the second Kurucz table's unit (icondopacp is
! a unit too).
      integer :: kurucz_table2_unit, icondopacp
      logical :: use_conductive_opacity

! Alexander table units (set by the reader). alex94_table_unit reads the
! Alexander 1994 low-T tables (kap/alex94/alex94_opacity.f90); the
! namelist names lalex95, use_alex95_tables and alex95_table_z_list keep
! their historical "95".
      integer :: alex94_table_unit
      integer :: alex06_table_unit

! The composition the Alexander tables are read at.
      double precision :: alex_mixture_x = 0.7d0, alex_mixture_z = 0.02d0

! Angular momentum transport options: the alternative velocity form in
! the diffusion coefficients (lvfc) and diffusion-advection transport.
      logical :: lvfc = .false.
      logical :: use_diffusion_advection_transport

! No angular momentum transport in the core.
      logical :: no_am_transport_in_core

! SCV EOS table units (set by the reader).
      integer :: scv_h_unit, scv_he_unit, scv_z_unit

! Allard atmosphere tables: the target [Fe/H] and [alpha/Fe], tau=100
! tables, and the unit (set by the reader).
      double precision :: allard_target_feh, allard_target_alpha
      logical :: allard_use_tau100
      integer :: allard_table_unit

! Solid-body rotation.
      logical :: force_solid_body_rotation
      integer :: solid_body_mode_flag

! Constant mixing and settling coefficients.
      double precision :: constant_mixing_coeff = 1.0d0, constant_settling_reduction = 1.0d0

! The timestep limit from the envelope-fit triangle.
      logical :: use_envelope_triangle_dt

! Magnetic braking: the Matt/PMM wind law (the reader overrides its
! parameters for the K97 and V13 laws).
      double precision :: pmm_exponent_a, pmm_exponent_b, pmm_exponent_c, &
           pmm_exponent_d, pmm_exponent_m, pmm_norm_jdot, pmm_norm_mdot, &
           pmm_solar_pressure, pmm_solar_omega, pmm_solar_turnover_timescale
      logical :: use_pmm_wind_law, scale_by_rossby_number, scale_by_b_field
      character(len=3) :: wind_law_name

! Magnetic braking: the Kawaler (1988) wind law and its saturation.
      double precision :: exmd, extau, exr, exm, exl, expr, constfactor, &
           excen, c_2, wind_law_omega_exponent
      double precision :: wind_saturation_omega
      logical :: use_wind_torque = .true.

! Constant background diffusion.
      double precision :: constant_background_diffusion_coeff
      logical :: use_constant_background_diffusion

! The timestep and zoning tolerances under readable names; when lnewvars
! is set, map_user_inputs fills atime and chi_grid_scale from them.
      double precision :: flag_dx = 0.05d0, flag_dw = 0.10d0, flag_dz = 0.05d0
      double precision :: time_core_min = 1.0d-3, time_dl = 2.0d-2, &
           time_dp = 4.0d-2, time_dr = 2.0d-2, time_dt = 2.0d-2, &
           time_dw_global = 8.0d-2, time_dw_mix = 8.0d-2, &
           time_dx_core_frac = 0.5d0, time_dx_core_tot = 2.0d-2, &
           time_dx_shell = 0.1d0, time_dx_total = 1.5d-3, &
           time_dy_core_frac = 0.5d0, time_dy_core_tot = 2.0d-2, &
           time_dy_shell = 0.1d0, time_dy_total = 1.5d-3, &
           time_max_dt_frac = 1.5d0
      double precision :: tol_czbase_fine_width = 0.0d0, tol_dl_max = 0.02d0, &
           tol_dm_max = 0.08d0, tol_dm_min = 1.0d-8, &
           tol_dp_core_max = 0.05d0, tol_dp_czbase_max = 0.05d0, &
           tol_dp_env_max = 0.05d0, tol_dx_max = 1.0d0, tol_dz_max = 1.0d0
      logical :: lstruct_time = .false., lnewvars = .false.

! User-supplied nuclear S-factors and their first and second derivative
! ratios (relative to the Adelberger et al. 1998 Solar Fusion I values).
      double precision :: s0_pp, s0_he3he3, s0_he3he4, s0_p_c12, s0_p_c13, &
           s0_p_n14, s0_p_o16, s0_be7_electron, s0_be7_p, &
           s0_n15_p_c12_branch, s0_n15_p_o16_branch
      double precision :: s0p_pp, s0p_he3he3, s0p_he3he4, s0p_p_c12, &
           s0p_p_c13, s0p_p_n14, s0p_p_o16, s0pp_p_c12, s0pp_p_c13, &
           s0pp_p_o16, s0p_be7_p, s0pp_be7_p
      double precision :: s0_pep = 3.5734d-6, s0_hep = 8.6d-20

end module controls_buffer

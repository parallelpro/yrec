!----------------------------------------------------------------------
! controls_reset_lib
!----------------------------------------------------------------------
! GENERATED (2026, phase five step C -- ROADMAP.md): pristine-snapshot
! support for controls_lib AND phys_const_lib, one mirror variable per
! member.
!
! Fresh-process semantics for repeated run_yrec calls in one process:
! parmin's namelist reads only OVERWRITE what the input file provides,
! so a second call would inherit run-1 values for any control the
! second namelist omits. phys_const is included because of cmixl, the
! documented const3 caveat member (copied from cmixla(nk) per kind
! card): the acceptance test caught it leaking run-1's mixing length
! into run 2's parmin echo. The other phys_const members are
! recomputed by setups each run; snapshotting them too is harmless
! and keeps the rule simple.
!
! controls_capture runs once, at the first run_yrec entry, recording
! the declaration-time default state; controls_restore rewinds to it
! at every later entry, so parmin always reads over pristine defaults
! -- exactly what a new process sees. Regenerate by re-running the
! generation script in the phase-five step C commit when either
! module's members change (the build fails loudly on drift).
module controls_reset_lib
      use controls_lib
      use phys_const_lib
      implicit none

      double precision :: snap_atime(14)
      double precision :: snap_tcut(5)
      double precision :: snap_saha_log10t_cutoff
      double precision :: snap_tenv0
      double precision :: snap_tenv1
      double precision :: snap_tgcut
      double precision :: snap_tenv
      double precision :: snap_cross_section_scale(17)
      double precision :: snap_qs0e_scale(8)
      double precision :: snap_qqs0ee_scale(8)
      double precision :: snap_o16_gamma_scale
      double precision :: snap_c12_alpha_scale
      logical :: snap_use_new_nuclear_rates
      double precision :: snap_weak_screening_threshold
      double precision :: snap_dpenv
      logical :: snap_lovstc
      double precision :: snap_alphac
      logical :: snap_envelope_overshoot_active
      double precision :: snap_alphae
      logical :: snap_lovstm
      double precision :: snap_alpham
      logical :: snap_ladov
      logical :: snap_lovmax
      double precision :: snap_betac
      logical :: snap_lsemic
      integer :: snap_iov1
      integer :: snap_iov2
      integer :: snap_iovim
      double precision :: snap_acfpft
      integer :: snap_itfp1
      integer :: snap_itfp2
      logical :: snap_rotation_active
      double precision :: snap_walpcz
      logical :: snap_instability_transport_active
      logical :: snap_lwnew
      double precision :: snap_wnew
      double precision :: snap_mass_accretion_rate
      double precision :: snap_fczdmdt
      double precision :: snap_ftotdmdt
      double precision :: snap_accreted_composition(15)
      double precision :: snap_creim
      logical :: snap_lreimer
      logical :: snap_use_mass_accretion
      integer :: snap_niter4
      logical :: snap_lnews
      logical :: snap_lsnu
      double precision :: snap_li6_rate_scale
      double precision :: snap_li7_rate_scale
      double precision :: snap_be9_pg_rate_scale
      double precision :: snap_be9_pd_rate_scale
      double precision :: snap_be9_palpha_rate_scale
      logical :: snap_use_itoh_neutrino_loss
      logical :: snap_use_new_turnover_timescale
      logical :: snap_calc_envelope_flag
      double precision :: snap_solar_luminosity_cgs
      double precision :: snap_log10_solar_luminosity
      double precision :: snap_ln_solar_luminosity
      double precision :: snap_solar_mass_cgs
      double precision :: snap_log10_solar_mass
      double precision :: snap_solar_radius_cgs
      double precision :: snap_log10_solar_radius
      double precision :: snap_solar_bolometric_magnitude
      logical :: snap_use_extended_composition
      double precision :: snap_spot_filling_factor
      double precision :: snap_spot_temp_contrast
      logical :: snap_spot_depth_varies
      double precision :: snap_disk_temperature
      double precision :: snap_disk_pressure
      logical :: snap_disk_locking_active
      double precision :: snap_disk_lifetime
      double precision :: snap_target_end_age(50)
      double precision :: snap_timestep_override(50)
      double precision :: snap_central_deuterium_stop(50)
      double precision :: snap_central_hydrogen_stop(50)
      double precision :: snap_central_helium_stop(50)
      logical :: snap_end_age_stop_active(50)
      logical :: snap_timestep_override_active(50)
      logical :: snap_use_mhd_eos
      integer :: snap_unit_zams_a
      integer :: snap_unit_zams_b
      integer :: snap_unit_zams_c
      integer :: snap_unit_centre1
      integer :: snap_unit_centre2
      integer :: snap_unit_centre3
      integer :: snap_unit_centre4
      integer :: snap_unit_centre5
      double precision :: snap_metal_fraction_match_tolerance
      double precision :: snap_zsi
      integer :: snap_idt
      integer :: snap_idd(4)
      logical :: snap_lstore
      logical :: snap_lstatm
      logical :: snap_lstenv
      logical :: snap_lstmod
      logical :: snap_lstphys
      logical :: snap_lstrot
      logical :: snap_lscrib
      logical :: snap_lstch
      logical :: snap_lphhd
      integer :: snap_npenv
      integer :: snap_nprtmod
      integer :: snap_npoint
      integer :: snap_print_point_interval
      logical :: snap_ldebug
      logical :: snap_lcorr
      logical :: snap_lmilne
      logical :: snap_ltrack
      logical :: snap_lstpch
      integer :: snap_first_unit
      integer :: snap_run_unit
      integer :: snap_standard_unit
      integer :: snap_fermi_unit
      integer :: snap_opal_model_unit
      integer :: snap_opal_envelope_unit
      integer :: snap_opal_atm_unit
      integer :: snap_dynamics_unit
      integer :: snap_laol_table_unit
      integer :: snap_neutrino_unit
      integer :: snap_composition_unit
      integer :: snap_kurucz_table_unit
      logical :: snap_lmonte
      integer :: snap_imbeg
      integer :: snap_imend
      double precision :: snap_htoler(5,2)
      double precision :: snap_fcorr0
      double precision :: snap_fcorri
      double precision :: snap_fcorr
      double precision :: snap_chi_grid_scale(12)
      integer :: snap_niter1
      integer :: snap_niter2
      integer :: snap_niter3
      double precision :: snap_dtdif
      double precision :: snap_convergence_tolerance
      integer :: snap_itdif1
      integer :: snap_itdif2
      double precision :: snap_settling_tolerance
      integer :: snap_coulomb_log_choice
      integer :: snap_settling_num_iterations
      logical :: snap_diffuse_helium_active
      double precision :: snap_settling_timestep_fraction
      double precision :: snap_hydrogen_diffusion_floor
      double precision :: snap_helium_diffusion_min
      logical :: snap_use_thoul_fit
      double precision :: snap_fgry
      double precision :: snap_fgrz
      logical :: snap_lthoul
      logical :: snap_use_diffusion_z
      logical :: snap_use_new_diffusion_routines
      logical :: snap_ldifli
      double precision :: snap_atm_error_tol
      double precision :: snap_atm_step_initial
      double precision :: snap_atm_step_begin
      double precision :: snap_atm_step_min
      double precision :: snap_atm_step_max
      double precision :: snap_env_error_tol
      double precision :: snap_env_step_begin
      double precision :: snap_env_step_min
      double precision :: snap_env_step_max
      logical :: snap_pulsation_output_active
      integer :: snap_pulsation_file_version
      double precision :: snap_pulsation_mass_msun
      integer :: snap_atm_choice
      integer :: snap_atm_choice_initial
      logical :: snap_use_ttau_relation
      double precision :: snap_atm_hras
      double precision :: snap_vnew(12)
      character(len=256) :: snap_monte_carlo_file1_path
      character(len=256) :: snap_monte_carlo_file2_path
      integer :: snap_monte_carlo_unit1
      integer :: snap_monte_carlo_unit2
      integer :: snap_pulse_gyre_interval
      double precision :: snap_tri_delta_teffl
      double precision :: snap_tri_delta_logl
      logical :: snap_lnew0
      double precision :: snap_requested_envelope_mass
      logical :: snap_change_envelope_mass_flag
      double precision :: snap_rescale_params(4,50)
      integer :: snap_rescale_kind(50)
      integer :: snap_num_models(50)
      integer :: snap_num_runs
      logical :: snap_first_call_flag(50)
      double precision :: snap_max_domega_global
      logical :: snap_use_structure_dt_limits
      double precision :: snap_atm_step_size
      double precision :: snap_envelope_step_size
      logical :: snap_envelope_generation_flag
      logical :: snap_helium_flash_active
      double precision :: snap_initial_envelope_x
      double precision :: snap_initial_envelope_z
      double precision :: snap_new_species_value
      integer :: snap_new_species_index
      logical :: snap_rescale_species_active
      logical :: snap_value_relative_to_h
      integer :: snap_mixture_change_mode
      integer :: snap_isotope_change_mode
      logical :: snap_change_cno_mixture_active
      logical :: snap_change_isotope_ratios_active
      double precision :: snap_target_carbon_cno_fraction
      double precision :: snap_target_nitrogen_cno_fraction
      double precision :: snap_target_oxygen_cno_fraction
      double precision :: snap_c12_to_c13_ratio
      double precision :: snap_n14_to_n15_ratio
      double precision :: snap_o16_to_o17_ratio
      double precision :: snap_o16_to_o18_ratio
      double precision :: snap_target_metal_fraction
      double precision :: snap_initial_h2_fraction
      double precision :: snap_initial_he3_fraction
      double precision :: snap_initial_li6_fraction
      double precision :: snap_initial_li7_fraction
      double precision :: snap_initial_be9_fraction
      double precision :: snap_initial_b10_fraction
      double precision :: snap_initial_b11_fraction
      double precision :: snap_fo
      double precision :: snap_difad_velocity_scale
      double precision :: snap_mixing_velocity_scale
      double precision :: snap_es_velocity_scale
      double precision :: snap_gsf_velocity_scale
      double precision :: snap_mu_gradient_scale
      double precision :: snap_secular_shear_velocity_scale
      double precision :: snap_critical_reynolds
      integer :: snap_ies
      integer :: snap_imu
      double precision :: snap_es_mixing_scale
      double precision :: snap_secular_shear_mixing_scale
      double precision :: snap_gsf_mixing_scale
      integer :: snap_gsf_inhibition_mode
      double precision :: snap_min_abundance
      double precision :: snap_absolute_tolerance
      double precision :: snap_relative_tolerance
      integer :: snap_max_burn_iterations
      integer :: snap_opal95_table_unit
      double precision :: snap_po_weight_l
      double precision :: snap_po_weight_teff
      double precision :: snap_po_weight_age
      double precision :: snap_po_max_len_sq
      logical :: snap_po_output_enabled
      integer :: snap_track_file_version
      logical :: snap_extend_core_inward
      integer :: snap_num_core_shells_added
      double precision :: snap_core_mass_reduction_factor
      logical :: snap_lrwsh_placeholder
      logical :: snap_isochrone_output_active
      integer :: snap_isochrone_file_unit
      double precision :: snap_initial_x_array(50)
      double precision :: snap_initial_z_array(50)
      double precision :: snap_mixing_length_array(50)
      double precision :: snap_senv0_array(50)
      logical :: snap_has_senv0_array(50)
      double precision :: snap_luminosity_tolerance
      double precision :: snap_radius_tolerance
      double precision :: snap_zx_tolerance
      double precision :: snap_target_solar_zx
      double precision :: snap_target_solar_age
      logical :: snap_calibrate_solar_model
      logical :: snap_calibrate_solar_zx
      double precision :: snap_rsclzc(50)
      double precision :: snap_rsclzm1(50)
      double precision :: snap_rsclzm2(50)
      integer :: snap_iolaol2
      integer :: snap_ioopal2
      integer :: snap_nk
      logical :: snap_use_z_ramp
      double precision :: snap_target_luminosity_lsun
      double precision :: snap_target_star_luminosity_tolerance
      double precision :: snap_target_teff
      double precision :: snap_target_radius_rsun
      double precision :: snap_log_l_prev_model
      double precision :: snap_log_r_prev_model
      double precision :: snap_age_at_target_radius
      double precision :: snap_log_l_at_target_radius
      double precision :: snap_log_l_at_target_radius_prev_run
      double precision :: snap_age_prev_model
      logical :: snap_star_found_flag
      logical :: snap_specify_teff_flag
      logical :: snap_just_passed_target_radius_flag
      logical :: snap_calibrate_star_flag
      logical :: snap_use_opal95_eos
      logical :: snap_use_opal2001_eos
      logical :: snap_use_opal2006_eos
      logical :: snap_use_numerical_derivatives
      integer :: snap_iopale
      double precision :: snap_laol_table_z1
      double precision :: snap_laol_table_z2
      double precision :: snap_opal_table_z1
      double precision :: snap_opal_table_z2
      double precision :: snap_opal95_single_table_z
      double precision :: snap_alex_table_z1
      double precision :: snap_kurucz_table_z1
      double precision :: snap_kurucz_table_z2
      double precision :: snap_molecular_opacity_logt_min
      double precision :: snap_molecular_opacity_logt_max
      logical :: snap_use_alex06_tables
      logical :: snap_use_laol89_tables
      logical :: snap_use_opal92_tables
      logical :: snap_use_opal95_tables
      logical :: snap_use_kurucz90_tables
      logical :: snap_use_alex95_tables
      logical :: snap_use_two_z_tables
      integer :: snap_ikur2
      integer :: snap_icondopacp
      logical :: snap_use_conductive_opacity
      integer :: snap_alex95_table_unit
      integer :: snap_alex06_table_unit
      double precision :: snap_alex_mixture_x
      double precision :: snap_alex_mixture_z
      double precision :: snap_vfc(5000)
      logical :: snap_lvfc
      logical :: snap_use_diffusion_advection_transport
      logical :: snap_no_am_transport_in_core
      integer :: snap_scv_h_unit
      integer :: snap_scv_he_unit
      integer :: snap_scv_z_unit
      double precision :: snap_allard_target_feh
      double precision :: snap_allard_target_alpha
      logical :: snap_allard_use_tau100
      integer :: snap_allard_table_unit
      logical :: snap_force_solid_body_rotation
      integer :: snap_solid_body_mode_flag
      double precision :: snap_cstmixing
      double precision :: snap_cstdiffmix
      double precision :: snap_taucz_placeholder
      double precision :: snap_deladj_placeholder(5000)
      double precision :: snap_tauhe_placeholder
      double precision :: snap_tnorm_placeholder
      double precision :: snap_tcz_placeholder
      double precision :: snap_whe_placeholder
      double precision :: snap_acatmr_placeholder(5000)
      double precision :: snap_acatmd_placeholder(5000)
      double precision :: snap_acatmp_placeholder(5000)
      double precision :: snap_acatmt_placeholder(5000)
      double precision :: snap_tatmos_placeholder
      integer :: snap_iacat_placeholder
      integer :: snap_ijlast_placeholder
      logical :: snap_laoly_placeholder
      integer :: snap_ijvs_placeholder
      integer :: snap_ijent_placeholder
      integer :: snap_ijdel_placeholder
      logical :: snap_compute_acoustic_depth
      double precision :: snap_ageout_placeholder(5)
      logical :: snap_lclcd_placeholder
      logical :: snap_ljlast_placeholder
      logical :: snap_ljwrt_placeholder
      logical :: snap_acoustic_depth_output
      integer :: snap_iclcd_placeholder
      logical :: snap_use_envelope_triangle_dt
      double precision :: snap_pmm_exponent_a
      double precision :: snap_pmm_exponent_b
      double precision :: snap_pmm_exponent_c
      double precision :: snap_pmm_exponent_d
      double precision :: snap_pmm_exponent_m
      double precision :: snap_pmm_norm_jdot
      double precision :: snap_pmm_norm_mdot
      double precision :: snap_pmm_solar_pressure
      double precision :: snap_pmm_solar_omega
      double precision :: snap_pmm_solar_turnover_timescale
      logical :: snap_use_pmm_wind_law
      logical :: snap_scale_by_rossby_number
      logical :: snap_scale_by_b_field
      character(len=3) :: snap_wind_law_name
      double precision :: snap_exmd
      double precision :: snap_extau
      double precision :: snap_exr
      double precision :: snap_exm
      double precision :: snap_exl
      double precision :: snap_expr
      double precision :: snap_constfactor
      double precision :: snap_structfactor
      double precision :: snap_excen
      double precision :: snap_c_2
      double precision :: snap_wind_law_omega_exponent
      double precision :: snap_wind_saturation_omega
      logical :: snap_ljdot0
      double precision :: snap_constant_background_diffusion_coeff
      logical :: snap_use_constant_background_diffusion
      double precision :: snap_flag_dx
      double precision :: snap_flag_dw
      double precision :: snap_flag_dz
      double precision :: snap_time_core_min
      double precision :: snap_time_dl
      double precision :: snap_time_dp
      double precision :: snap_time_dr
      double precision :: snap_time_dt
      double precision :: snap_time_dw_global
      double precision :: snap_time_dw_mix
      double precision :: snap_time_dx_core_frac
      double precision :: snap_time_dx_core_tot
      double precision :: snap_time_dx_shell
      double precision :: snap_time_dx_total
      double precision :: snap_time_dy_core_frac
      double precision :: snap_time_dy_core_tot
      double precision :: snap_time_dy_shell
      double precision :: snap_time_dy_total
      double precision :: snap_time_max_dt_frac
      double precision :: snap_tol_czbase_fine_width
      double precision :: snap_tol_dl_max
      double precision :: snap_tol_dm_max
      double precision :: snap_tol_dm_min
      double precision :: snap_tol_dp_core_max
      double precision :: snap_tol_dp_czbase_max
      double precision :: snap_tol_dp_env_max
      double precision :: snap_tol_dx_max
      double precision :: snap_tol_dz_max
      logical :: snap_lstruct_time
      logical :: snap_lnewvars
      double precision :: snap_s0_pp
      double precision :: snap_s0_he3he3
      double precision :: snap_s0_he3he4
      double precision :: snap_s0_p_c12
      double precision :: snap_s0_p_c13
      double precision :: snap_s0_p_n14
      double precision :: snap_s0_p_o16
      double precision :: snap_s0_be7_electron
      double precision :: snap_s0_be7_p
      double precision :: snap_s0_n15_p_c12_branch
      double precision :: snap_s0_n15_p_o16_branch
      double precision :: snap_s0p_pp
      double precision :: snap_s0p_he3he3
      double precision :: snap_s0p_he3he4
      double precision :: snap_s0p_p_c12
      double precision :: snap_s0p_p_c13
      double precision :: snap_s0p_p_n14
      double precision :: snap_s0p_p_o16
      double precision :: snap_s0pp_p_c12
      double precision :: snap_s0pp_p_c13
      double precision :: snap_s0pp_p_o16
      double precision :: snap_s0p_be7_p
      double precision :: snap_s0pp_be7_p
      double precision :: snap_s0_pep
      double precision :: snap_s0_hep
      double precision :: snap_ln10
      double precision :: snap_clni
      double precision :: snap_c4pi
      double precision :: snap_c4pil
      double precision :: snap_c4pi3l
      double precision :: snap_cc13
      double precision :: snap_cc23
      double precision :: snap_cpi
      double precision :: snap_gas_constant
      double precision :: snap_radiation_constant_over_3
      double precision :: snap_ca3l
      double precision :: snap_csig
      double precision :: snap_csigl
      double precision :: snap_cgl
      double precision :: snap_cmkh
      double precision :: snap_cmkhn
      double precision :: snap_cdelrl
      double precision :: snap_cmixl2
      double precision :: snap_cmixl3
      double precision :: snap_clndp
      double precision :: snap_seconds_per_year
      double precision :: snap_cmixl
      character(len=10) :: snap_yrec_version_string
      character(len=20) :: snap_git_hash_string

contains

subroutine controls_capture
      snap_atime = atime
      snap_tcut = tcut
      snap_saha_log10t_cutoff = saha_log10t_cutoff
      snap_tenv0 = tenv0
      snap_tenv1 = tenv1
      snap_tgcut = tgcut
      snap_tenv = tenv
      snap_cross_section_scale = cross_section_scale
      snap_qs0e_scale = qs0e_scale
      snap_qqs0ee_scale = qqs0ee_scale
      snap_o16_gamma_scale = o16_gamma_scale
      snap_c12_alpha_scale = c12_alpha_scale
      snap_use_new_nuclear_rates = use_new_nuclear_rates
      snap_weak_screening_threshold = weak_screening_threshold
      snap_dpenv = dpenv
      snap_lovstc = lovstc
      snap_alphac = alphac
      snap_envelope_overshoot_active = envelope_overshoot_active
      snap_alphae = alphae
      snap_lovstm = lovstm
      snap_alpham = alpham
      snap_ladov = ladov
      snap_lovmax = lovmax
      snap_betac = betac
      snap_lsemic = lsemic
      snap_iov1 = iov1
      snap_iov2 = iov2
      snap_iovim = iovim
      snap_acfpft = acfpft
      snap_itfp1 = itfp1
      snap_itfp2 = itfp2
      snap_rotation_active = rotation_active
      snap_walpcz = walpcz
      snap_instability_transport_active = instability_transport_active
      snap_lwnew = lwnew
      snap_wnew = wnew
      snap_mass_accretion_rate = mass_accretion_rate
      snap_fczdmdt = fczdmdt
      snap_ftotdmdt = ftotdmdt
      snap_accreted_composition = accreted_composition
      snap_creim = creim
      snap_lreimer = lreimer
      snap_use_mass_accretion = use_mass_accretion
      snap_niter4 = niter4
      snap_lnews = lnews
      snap_lsnu = lsnu
      snap_li6_rate_scale = li6_rate_scale
      snap_li7_rate_scale = li7_rate_scale
      snap_be9_pg_rate_scale = be9_pg_rate_scale
      snap_be9_pd_rate_scale = be9_pd_rate_scale
      snap_be9_palpha_rate_scale = be9_palpha_rate_scale
      snap_use_itoh_neutrino_loss = use_itoh_neutrino_loss
      snap_use_new_turnover_timescale = use_new_turnover_timescale
      snap_calc_envelope_flag = calc_envelope_flag
      snap_solar_luminosity_cgs = solar_luminosity_cgs
      snap_log10_solar_luminosity = log10_solar_luminosity
      snap_ln_solar_luminosity = ln_solar_luminosity
      snap_solar_mass_cgs = solar_mass_cgs
      snap_log10_solar_mass = log10_solar_mass
      snap_solar_radius_cgs = solar_radius_cgs
      snap_log10_solar_radius = log10_solar_radius
      snap_solar_bolometric_magnitude = solar_bolometric_magnitude
      snap_use_extended_composition = use_extended_composition
      snap_spot_filling_factor = spot_filling_factor
      snap_spot_temp_contrast = spot_temp_contrast
      snap_spot_depth_varies = spot_depth_varies
      snap_disk_temperature = disk_temperature
      snap_disk_pressure = disk_pressure
      snap_disk_locking_active = disk_locking_active
      snap_disk_lifetime = disk_lifetime
      snap_target_end_age = target_end_age
      snap_timestep_override = timestep_override
      snap_central_deuterium_stop = central_deuterium_stop
      snap_central_hydrogen_stop = central_hydrogen_stop
      snap_central_helium_stop = central_helium_stop
      snap_end_age_stop_active = end_age_stop_active
      snap_timestep_override_active = timestep_override_active
      snap_use_mhd_eos = use_mhd_eos
      snap_unit_zams_a = unit_zams_a
      snap_unit_zams_b = unit_zams_b
      snap_unit_zams_c = unit_zams_c
      snap_unit_centre1 = unit_centre1
      snap_unit_centre2 = unit_centre2
      snap_unit_centre3 = unit_centre3
      snap_unit_centre4 = unit_centre4
      snap_unit_centre5 = unit_centre5
      snap_metal_fraction_match_tolerance = metal_fraction_match_tolerance
      snap_zsi = zsi
      snap_idt = idt
      snap_idd = idd
      snap_lstore = lstore
      snap_lstatm = lstatm
      snap_lstenv = lstenv
      snap_lstmod = lstmod
      snap_lstphys = lstphys
      snap_lstrot = lstrot
      snap_lscrib = lscrib
      snap_lstch = lstch
      snap_lphhd = lphhd
      snap_npenv = npenv
      snap_nprtmod = nprtmod
      snap_npoint = npoint
      snap_print_point_interval = print_point_interval
      snap_ldebug = ldebug
      snap_lcorr = lcorr
      snap_lmilne = lmilne
      snap_ltrack = ltrack
      snap_lstpch = lstpch
      snap_first_unit = first_unit
      snap_run_unit = run_unit
      snap_standard_unit = standard_unit
      snap_fermi_unit = fermi_unit
      snap_opal_model_unit = opal_model_unit
      snap_opal_envelope_unit = opal_envelope_unit
      snap_opal_atm_unit = opal_atm_unit
      snap_dynamics_unit = dynamics_unit
      snap_laol_table_unit = laol_table_unit
      snap_neutrino_unit = neutrino_unit
      snap_composition_unit = composition_unit
      snap_kurucz_table_unit = kurucz_table_unit
      snap_lmonte = lmonte
      snap_imbeg = imbeg
      snap_imend = imend
      snap_htoler = htoler
      snap_fcorr0 = fcorr0
      snap_fcorri = fcorri
      snap_fcorr = fcorr
      snap_chi_grid_scale = chi_grid_scale
      snap_niter1 = niter1
      snap_niter2 = niter2
      snap_niter3 = niter3
      snap_dtdif = dtdif
      snap_convergence_tolerance = convergence_tolerance
      snap_itdif1 = itdif1
      snap_itdif2 = itdif2
      snap_settling_tolerance = settling_tolerance
      snap_coulomb_log_choice = coulomb_log_choice
      snap_settling_num_iterations = settling_num_iterations
      snap_diffuse_helium_active = diffuse_helium_active
      snap_settling_timestep_fraction = settling_timestep_fraction
      snap_hydrogen_diffusion_floor = hydrogen_diffusion_floor
      snap_helium_diffusion_min = helium_diffusion_min
      snap_use_thoul_fit = use_thoul_fit
      snap_fgry = fgry
      snap_fgrz = fgrz
      snap_lthoul = lthoul
      snap_use_diffusion_z = use_diffusion_z
      snap_use_new_diffusion_routines = use_new_diffusion_routines
      snap_ldifli = ldifli
      snap_atm_error_tol = atm_error_tol
      snap_atm_step_initial = atm_step_initial
      snap_atm_step_begin = atm_step_begin
      snap_atm_step_min = atm_step_min
      snap_atm_step_max = atm_step_max
      snap_env_error_tol = env_error_tol
      snap_env_step_begin = env_step_begin
      snap_env_step_min = env_step_min
      snap_env_step_max = env_step_max
      snap_pulsation_output_active = pulsation_output_active
      snap_pulsation_file_version = pulsation_file_version
      snap_pulsation_mass_msun = pulsation_mass_msun
      snap_atm_choice = atm_choice
      snap_atm_choice_initial = atm_choice_initial
      snap_use_ttau_relation = use_ttau_relation
      snap_atm_hras = atm_hras
      snap_vnew = vnew
      snap_monte_carlo_file1_path = monte_carlo_file1_path
      snap_monte_carlo_file2_path = monte_carlo_file2_path
      snap_monte_carlo_unit1 = monte_carlo_unit1
      snap_monte_carlo_unit2 = monte_carlo_unit2
      snap_pulse_gyre_interval = pulse_gyre_interval
      snap_tri_delta_teffl = tri_delta_teffl
      snap_tri_delta_logl = tri_delta_logl
      snap_lnew0 = lnew0
      snap_requested_envelope_mass = requested_envelope_mass
      snap_change_envelope_mass_flag = change_envelope_mass_flag
      snap_rescale_params = rescale_params
      snap_rescale_kind = rescale_kind
      snap_num_models = num_models
      snap_num_runs = num_runs
      snap_first_call_flag = first_call_flag
      snap_max_domega_global = max_domega_global
      snap_use_structure_dt_limits = use_structure_dt_limits
      snap_atm_step_size = atm_step_size
      snap_envelope_step_size = envelope_step_size
      snap_envelope_generation_flag = envelope_generation_flag
      snap_helium_flash_active = helium_flash_active
      snap_initial_envelope_x = initial_envelope_x
      snap_initial_envelope_z = initial_envelope_z
      snap_new_species_value = new_species_value
      snap_new_species_index = new_species_index
      snap_rescale_species_active = rescale_species_active
      snap_value_relative_to_h = value_relative_to_h
      snap_mixture_change_mode = mixture_change_mode
      snap_isotope_change_mode = isotope_change_mode
      snap_change_cno_mixture_active = change_cno_mixture_active
      snap_change_isotope_ratios_active = change_isotope_ratios_active
      snap_target_carbon_cno_fraction = target_carbon_cno_fraction
      snap_target_nitrogen_cno_fraction = target_nitrogen_cno_fraction
      snap_target_oxygen_cno_fraction = target_oxygen_cno_fraction
      snap_c12_to_c13_ratio = c12_to_c13_ratio
      snap_n14_to_n15_ratio = n14_to_n15_ratio
      snap_o16_to_o17_ratio = o16_to_o17_ratio
      snap_o16_to_o18_ratio = o16_to_o18_ratio
      snap_target_metal_fraction = target_metal_fraction
      snap_initial_h2_fraction = initial_h2_fraction
      snap_initial_he3_fraction = initial_he3_fraction
      snap_initial_li6_fraction = initial_li6_fraction
      snap_initial_li7_fraction = initial_li7_fraction
      snap_initial_be9_fraction = initial_be9_fraction
      snap_initial_b10_fraction = initial_b10_fraction
      snap_initial_b11_fraction = initial_b11_fraction
      snap_fo = fo
      snap_difad_velocity_scale = difad_velocity_scale
      snap_mixing_velocity_scale = mixing_velocity_scale
      snap_es_velocity_scale = es_velocity_scale
      snap_gsf_velocity_scale = gsf_velocity_scale
      snap_mu_gradient_scale = mu_gradient_scale
      snap_secular_shear_velocity_scale = secular_shear_velocity_scale
      snap_critical_reynolds = critical_reynolds
      snap_ies = ies
      snap_imu = imu
      snap_es_mixing_scale = es_mixing_scale
      snap_secular_shear_mixing_scale = secular_shear_mixing_scale
      snap_gsf_mixing_scale = gsf_mixing_scale
      snap_gsf_inhibition_mode = gsf_inhibition_mode
      snap_min_abundance = min_abundance
      snap_absolute_tolerance = absolute_tolerance
      snap_relative_tolerance = relative_tolerance
      snap_max_burn_iterations = max_burn_iterations
      snap_opal95_table_unit = opal95_table_unit
      snap_po_weight_l = po_weight_l
      snap_po_weight_teff = po_weight_teff
      snap_po_weight_age = po_weight_age
      snap_po_max_len_sq = po_max_len_sq
      snap_po_output_enabled = po_output_enabled
      snap_track_file_version = track_file_version
      snap_extend_core_inward = extend_core_inward
      snap_num_core_shells_added = num_core_shells_added
      snap_core_mass_reduction_factor = core_mass_reduction_factor
      snap_lrwsh_placeholder = lrwsh_placeholder
      snap_isochrone_output_active = isochrone_output_active
      snap_isochrone_file_unit = isochrone_file_unit
      snap_initial_x_array = initial_x_array
      snap_initial_z_array = initial_z_array
      snap_mixing_length_array = mixing_length_array
      snap_senv0_array = senv0_array
      snap_has_senv0_array = has_senv0_array
      snap_luminosity_tolerance = luminosity_tolerance
      snap_radius_tolerance = radius_tolerance
      snap_zx_tolerance = zx_tolerance
      snap_target_solar_zx = target_solar_zx
      snap_target_solar_age = target_solar_age
      snap_calibrate_solar_model = calibrate_solar_model
      snap_calibrate_solar_zx = calibrate_solar_zx
      snap_rsclzc = rsclzc
      snap_rsclzm1 = rsclzm1
      snap_rsclzm2 = rsclzm2
      snap_iolaol2 = iolaol2
      snap_ioopal2 = ioopal2
      snap_nk = nk
      snap_use_z_ramp = use_z_ramp
      snap_target_luminosity_lsun = target_luminosity_lsun
      snap_target_star_luminosity_tolerance = target_star_luminosity_tolerance
      snap_target_teff = target_teff
      snap_target_radius_rsun = target_radius_rsun
      snap_log_l_prev_model = log_l_prev_model
      snap_log_r_prev_model = log_r_prev_model
      snap_age_at_target_radius = age_at_target_radius
      snap_log_l_at_target_radius = log_l_at_target_radius
      snap_log_l_at_target_radius_prev_run = log_l_at_target_radius_prev_run
      snap_age_prev_model = age_prev_model
      snap_star_found_flag = star_found_flag
      snap_specify_teff_flag = specify_teff_flag
      snap_just_passed_target_radius_flag = just_passed_target_radius_flag
      snap_calibrate_star_flag = calibrate_star_flag
      snap_use_opal95_eos = use_opal95_eos
      snap_use_opal2001_eos = use_opal2001_eos
      snap_use_opal2006_eos = use_opal2006_eos
      snap_use_numerical_derivatives = use_numerical_derivatives
      snap_iopale = iopale
      snap_laol_table_z1 = laol_table_z1
      snap_laol_table_z2 = laol_table_z2
      snap_opal_table_z1 = opal_table_z1
      snap_opal_table_z2 = opal_table_z2
      snap_opal95_single_table_z = opal95_single_table_z
      snap_alex_table_z1 = alex_table_z1
      snap_kurucz_table_z1 = kurucz_table_z1
      snap_kurucz_table_z2 = kurucz_table_z2
      snap_molecular_opacity_logt_min = molecular_opacity_logt_min
      snap_molecular_opacity_logt_max = molecular_opacity_logt_max
      snap_use_alex06_tables = use_alex06_tables
      snap_use_laol89_tables = use_laol89_tables
      snap_use_opal92_tables = use_opal92_tables
      snap_use_opal95_tables = use_opal95_tables
      snap_use_kurucz90_tables = use_kurucz90_tables
      snap_use_alex95_tables = use_alex95_tables
      snap_use_two_z_tables = use_two_z_tables
      snap_ikur2 = ikur2
      snap_icondopacp = icondopacp
      snap_use_conductive_opacity = use_conductive_opacity
      snap_alex95_table_unit = alex95_table_unit
      snap_alex06_table_unit = alex06_table_unit
      snap_alex_mixture_x = alex_mixture_x
      snap_alex_mixture_z = alex_mixture_z
      snap_vfc = vfc
      snap_lvfc = lvfc
      snap_use_diffusion_advection_transport = use_diffusion_advection_transport
      snap_no_am_transport_in_core = no_am_transport_in_core
      snap_scv_h_unit = scv_h_unit
      snap_scv_he_unit = scv_he_unit
      snap_scv_z_unit = scv_z_unit
      snap_allard_target_feh = allard_target_feh
      snap_allard_target_alpha = allard_target_alpha
      snap_allard_use_tau100 = allard_use_tau100
      snap_allard_table_unit = allard_table_unit
      snap_force_solid_body_rotation = force_solid_body_rotation
      snap_solid_body_mode_flag = solid_body_mode_flag
      snap_cstmixing = cstmixing
      snap_cstdiffmix = cstdiffmix
      snap_taucz_placeholder = taucz_placeholder
      snap_deladj_placeholder = deladj_placeholder
      snap_tauhe_placeholder = tauhe_placeholder
      snap_tnorm_placeholder = tnorm_placeholder
      snap_tcz_placeholder = tcz_placeholder
      snap_whe_placeholder = whe_placeholder
      snap_acatmr_placeholder = acatmr_placeholder
      snap_acatmd_placeholder = acatmd_placeholder
      snap_acatmp_placeholder = acatmp_placeholder
      snap_acatmt_placeholder = acatmt_placeholder
      snap_tatmos_placeholder = tatmos_placeholder
      snap_iacat_placeholder = iacat_placeholder
      snap_ijlast_placeholder = ijlast_placeholder
      snap_laoly_placeholder = laoly_placeholder
      snap_ijvs_placeholder = ijvs_placeholder
      snap_ijent_placeholder = ijent_placeholder
      snap_ijdel_placeholder = ijdel_placeholder
      snap_compute_acoustic_depth = compute_acoustic_depth
      snap_ageout_placeholder = ageout_placeholder
      snap_lclcd_placeholder = lclcd_placeholder
      snap_ljlast_placeholder = ljlast_placeholder
      snap_ljwrt_placeholder = ljwrt_placeholder
      snap_acoustic_depth_output = acoustic_depth_output
      snap_iclcd_placeholder = iclcd_placeholder
      snap_use_envelope_triangle_dt = use_envelope_triangle_dt
      snap_pmm_exponent_a = pmm_exponent_a
      snap_pmm_exponent_b = pmm_exponent_b
      snap_pmm_exponent_c = pmm_exponent_c
      snap_pmm_exponent_d = pmm_exponent_d
      snap_pmm_exponent_m = pmm_exponent_m
      snap_pmm_norm_jdot = pmm_norm_jdot
      snap_pmm_norm_mdot = pmm_norm_mdot
      snap_pmm_solar_pressure = pmm_solar_pressure
      snap_pmm_solar_omega = pmm_solar_omega
      snap_pmm_solar_turnover_timescale = pmm_solar_turnover_timescale
      snap_use_pmm_wind_law = use_pmm_wind_law
      snap_scale_by_rossby_number = scale_by_rossby_number
      snap_scale_by_b_field = scale_by_b_field
      snap_wind_law_name = wind_law_name
      snap_exmd = exmd
      snap_extau = extau
      snap_exr = exr
      snap_exm = exm
      snap_exl = exl
      snap_expr = expr
      snap_constfactor = constfactor
      snap_structfactor = structfactor
      snap_excen = excen
      snap_c_2 = c_2
      snap_wind_law_omega_exponent = wind_law_omega_exponent
      snap_wind_saturation_omega = wind_saturation_omega
      snap_ljdot0 = ljdot0
      snap_constant_background_diffusion_coeff = constant_background_diffusion_coeff
      snap_use_constant_background_diffusion = use_constant_background_diffusion
      snap_flag_dx = flag_dx
      snap_flag_dw = flag_dw
      snap_flag_dz = flag_dz
      snap_time_core_min = time_core_min
      snap_time_dl = time_dl
      snap_time_dp = time_dp
      snap_time_dr = time_dr
      snap_time_dt = time_dt
      snap_time_dw_global = time_dw_global
      snap_time_dw_mix = time_dw_mix
      snap_time_dx_core_frac = time_dx_core_frac
      snap_time_dx_core_tot = time_dx_core_tot
      snap_time_dx_shell = time_dx_shell
      snap_time_dx_total = time_dx_total
      snap_time_dy_core_frac = time_dy_core_frac
      snap_time_dy_core_tot = time_dy_core_tot
      snap_time_dy_shell = time_dy_shell
      snap_time_dy_total = time_dy_total
      snap_time_max_dt_frac = time_max_dt_frac
      snap_tol_czbase_fine_width = tol_czbase_fine_width
      snap_tol_dl_max = tol_dl_max
      snap_tol_dm_max = tol_dm_max
      snap_tol_dm_min = tol_dm_min
      snap_tol_dp_core_max = tol_dp_core_max
      snap_tol_dp_czbase_max = tol_dp_czbase_max
      snap_tol_dp_env_max = tol_dp_env_max
      snap_tol_dx_max = tol_dx_max
      snap_tol_dz_max = tol_dz_max
      snap_lstruct_time = lstruct_time
      snap_lnewvars = lnewvars
      snap_s0_pp = s0_pp
      snap_s0_he3he3 = s0_he3he3
      snap_s0_he3he4 = s0_he3he4
      snap_s0_p_c12 = s0_p_c12
      snap_s0_p_c13 = s0_p_c13
      snap_s0_p_n14 = s0_p_n14
      snap_s0_p_o16 = s0_p_o16
      snap_s0_be7_electron = s0_be7_electron
      snap_s0_be7_p = s0_be7_p
      snap_s0_n15_p_c12_branch = s0_n15_p_c12_branch
      snap_s0_n15_p_o16_branch = s0_n15_p_o16_branch
      snap_s0p_pp = s0p_pp
      snap_s0p_he3he3 = s0p_he3he3
      snap_s0p_he3he4 = s0p_he3he4
      snap_s0p_p_c12 = s0p_p_c12
      snap_s0p_p_c13 = s0p_p_c13
      snap_s0p_p_n14 = s0p_p_n14
      snap_s0p_p_o16 = s0p_p_o16
      snap_s0pp_p_c12 = s0pp_p_c12
      snap_s0pp_p_c13 = s0pp_p_c13
      snap_s0pp_p_o16 = s0pp_p_o16
      snap_s0p_be7_p = s0p_be7_p
      snap_s0pp_be7_p = s0pp_be7_p
      snap_s0_pep = s0_pep
      snap_s0_hep = s0_hep
      snap_ln10 = ln10
      snap_clni = clni
      snap_c4pi = c4pi
      snap_c4pil = c4pil
      snap_c4pi3l = c4pi3l
      snap_cc13 = cc13
      snap_cc23 = cc23
      snap_cpi = cpi
      snap_gas_constant = gas_constant
      snap_radiation_constant_over_3 = radiation_constant_over_3
      snap_ca3l = ca3l
      snap_csig = csig
      snap_csigl = csigl
      snap_cgl = cgl
      snap_cmkh = cmkh
      snap_cmkhn = cmkhn
      snap_cdelrl = cdelrl
      snap_cmixl2 = cmixl2
      snap_cmixl3 = cmixl3
      snap_clndp = clndp
      snap_seconds_per_year = seconds_per_year
      snap_cmixl = cmixl
      snap_yrec_version_string = yrec_version_string
      snap_git_hash_string = git_hash_string
end subroutine controls_capture

subroutine controls_restore
      atime = snap_atime
      tcut = snap_tcut
      saha_log10t_cutoff = snap_saha_log10t_cutoff
      tenv0 = snap_tenv0
      tenv1 = snap_tenv1
      tgcut = snap_tgcut
      tenv = snap_tenv
      cross_section_scale = snap_cross_section_scale
      qs0e_scale = snap_qs0e_scale
      qqs0ee_scale = snap_qqs0ee_scale
      o16_gamma_scale = snap_o16_gamma_scale
      c12_alpha_scale = snap_c12_alpha_scale
      use_new_nuclear_rates = snap_use_new_nuclear_rates
      weak_screening_threshold = snap_weak_screening_threshold
      dpenv = snap_dpenv
      lovstc = snap_lovstc
      alphac = snap_alphac
      envelope_overshoot_active = snap_envelope_overshoot_active
      alphae = snap_alphae
      lovstm = snap_lovstm
      alpham = snap_alpham
      ladov = snap_ladov
      lovmax = snap_lovmax
      betac = snap_betac
      lsemic = snap_lsemic
      iov1 = snap_iov1
      iov2 = snap_iov2
      iovim = snap_iovim
      acfpft = snap_acfpft
      itfp1 = snap_itfp1
      itfp2 = snap_itfp2
      rotation_active = snap_rotation_active
      walpcz = snap_walpcz
      instability_transport_active = snap_instability_transport_active
      lwnew = snap_lwnew
      wnew = snap_wnew
      mass_accretion_rate = snap_mass_accretion_rate
      fczdmdt = snap_fczdmdt
      ftotdmdt = snap_ftotdmdt
      accreted_composition = snap_accreted_composition
      creim = snap_creim
      lreimer = snap_lreimer
      use_mass_accretion = snap_use_mass_accretion
      niter4 = snap_niter4
      lnews = snap_lnews
      lsnu = snap_lsnu
      li6_rate_scale = snap_li6_rate_scale
      li7_rate_scale = snap_li7_rate_scale
      be9_pg_rate_scale = snap_be9_pg_rate_scale
      be9_pd_rate_scale = snap_be9_pd_rate_scale
      be9_palpha_rate_scale = snap_be9_palpha_rate_scale
      use_itoh_neutrino_loss = snap_use_itoh_neutrino_loss
      use_new_turnover_timescale = snap_use_new_turnover_timescale
      calc_envelope_flag = snap_calc_envelope_flag
      solar_luminosity_cgs = snap_solar_luminosity_cgs
      log10_solar_luminosity = snap_log10_solar_luminosity
      ln_solar_luminosity = snap_ln_solar_luminosity
      solar_mass_cgs = snap_solar_mass_cgs
      log10_solar_mass = snap_log10_solar_mass
      solar_radius_cgs = snap_solar_radius_cgs
      log10_solar_radius = snap_log10_solar_radius
      solar_bolometric_magnitude = snap_solar_bolometric_magnitude
      use_extended_composition = snap_use_extended_composition
      spot_filling_factor = snap_spot_filling_factor
      spot_temp_contrast = snap_spot_temp_contrast
      spot_depth_varies = snap_spot_depth_varies
      disk_temperature = snap_disk_temperature
      disk_pressure = snap_disk_pressure
      disk_locking_active = snap_disk_locking_active
      disk_lifetime = snap_disk_lifetime
      target_end_age = snap_target_end_age
      timestep_override = snap_timestep_override
      central_deuterium_stop = snap_central_deuterium_stop
      central_hydrogen_stop = snap_central_hydrogen_stop
      central_helium_stop = snap_central_helium_stop
      end_age_stop_active = snap_end_age_stop_active
      timestep_override_active = snap_timestep_override_active
      use_mhd_eos = snap_use_mhd_eos
      unit_zams_a = snap_unit_zams_a
      unit_zams_b = snap_unit_zams_b
      unit_zams_c = snap_unit_zams_c
      unit_centre1 = snap_unit_centre1
      unit_centre2 = snap_unit_centre2
      unit_centre3 = snap_unit_centre3
      unit_centre4 = snap_unit_centre4
      unit_centre5 = snap_unit_centre5
      metal_fraction_match_tolerance = snap_metal_fraction_match_tolerance
      zsi = snap_zsi
      idt = snap_idt
      idd = snap_idd
      lstore = snap_lstore
      lstatm = snap_lstatm
      lstenv = snap_lstenv
      lstmod = snap_lstmod
      lstphys = snap_lstphys
      lstrot = snap_lstrot
      lscrib = snap_lscrib
      lstch = snap_lstch
      lphhd = snap_lphhd
      npenv = snap_npenv
      nprtmod = snap_nprtmod
      npoint = snap_npoint
      print_point_interval = snap_print_point_interval
      ldebug = snap_ldebug
      lcorr = snap_lcorr
      lmilne = snap_lmilne
      ltrack = snap_ltrack
      lstpch = snap_lstpch
      first_unit = snap_first_unit
      run_unit = snap_run_unit
      standard_unit = snap_standard_unit
      fermi_unit = snap_fermi_unit
      opal_model_unit = snap_opal_model_unit
      opal_envelope_unit = snap_opal_envelope_unit
      opal_atm_unit = snap_opal_atm_unit
      dynamics_unit = snap_dynamics_unit
      laol_table_unit = snap_laol_table_unit
      neutrino_unit = snap_neutrino_unit
      composition_unit = snap_composition_unit
      kurucz_table_unit = snap_kurucz_table_unit
      lmonte = snap_lmonte
      imbeg = snap_imbeg
      imend = snap_imend
      htoler = snap_htoler
      fcorr0 = snap_fcorr0
      fcorri = snap_fcorri
      fcorr = snap_fcorr
      chi_grid_scale = snap_chi_grid_scale
      niter1 = snap_niter1
      niter2 = snap_niter2
      niter3 = snap_niter3
      dtdif = snap_dtdif
      convergence_tolerance = snap_convergence_tolerance
      itdif1 = snap_itdif1
      itdif2 = snap_itdif2
      settling_tolerance = snap_settling_tolerance
      coulomb_log_choice = snap_coulomb_log_choice
      settling_num_iterations = snap_settling_num_iterations
      diffuse_helium_active = snap_diffuse_helium_active
      settling_timestep_fraction = snap_settling_timestep_fraction
      hydrogen_diffusion_floor = snap_hydrogen_diffusion_floor
      helium_diffusion_min = snap_helium_diffusion_min
      use_thoul_fit = snap_use_thoul_fit
      fgry = snap_fgry
      fgrz = snap_fgrz
      lthoul = snap_lthoul
      use_diffusion_z = snap_use_diffusion_z
      use_new_diffusion_routines = snap_use_new_diffusion_routines
      ldifli = snap_ldifli
      atm_error_tol = snap_atm_error_tol
      atm_step_initial = snap_atm_step_initial
      atm_step_begin = snap_atm_step_begin
      atm_step_min = snap_atm_step_min
      atm_step_max = snap_atm_step_max
      env_error_tol = snap_env_error_tol
      env_step_begin = snap_env_step_begin
      env_step_min = snap_env_step_min
      env_step_max = snap_env_step_max
      pulsation_output_active = snap_pulsation_output_active
      pulsation_file_version = snap_pulsation_file_version
      pulsation_mass_msun = snap_pulsation_mass_msun
      atm_choice = snap_atm_choice
      atm_choice_initial = snap_atm_choice_initial
      use_ttau_relation = snap_use_ttau_relation
      atm_hras = snap_atm_hras
      vnew = snap_vnew
      monte_carlo_file1_path = snap_monte_carlo_file1_path
      monte_carlo_file2_path = snap_monte_carlo_file2_path
      monte_carlo_unit1 = snap_monte_carlo_unit1
      monte_carlo_unit2 = snap_monte_carlo_unit2
      pulse_gyre_interval = snap_pulse_gyre_interval
      tri_delta_teffl = snap_tri_delta_teffl
      tri_delta_logl = snap_tri_delta_logl
      lnew0 = snap_lnew0
      requested_envelope_mass = snap_requested_envelope_mass
      change_envelope_mass_flag = snap_change_envelope_mass_flag
      rescale_params = snap_rescale_params
      rescale_kind = snap_rescale_kind
      num_models = snap_num_models
      num_runs = snap_num_runs
      first_call_flag = snap_first_call_flag
      max_domega_global = snap_max_domega_global
      use_structure_dt_limits = snap_use_structure_dt_limits
      atm_step_size = snap_atm_step_size
      envelope_step_size = snap_envelope_step_size
      envelope_generation_flag = snap_envelope_generation_flag
      helium_flash_active = snap_helium_flash_active
      initial_envelope_x = snap_initial_envelope_x
      initial_envelope_z = snap_initial_envelope_z
      new_species_value = snap_new_species_value
      new_species_index = snap_new_species_index
      rescale_species_active = snap_rescale_species_active
      value_relative_to_h = snap_value_relative_to_h
      mixture_change_mode = snap_mixture_change_mode
      isotope_change_mode = snap_isotope_change_mode
      change_cno_mixture_active = snap_change_cno_mixture_active
      change_isotope_ratios_active = snap_change_isotope_ratios_active
      target_carbon_cno_fraction = snap_target_carbon_cno_fraction
      target_nitrogen_cno_fraction = snap_target_nitrogen_cno_fraction
      target_oxygen_cno_fraction = snap_target_oxygen_cno_fraction
      c12_to_c13_ratio = snap_c12_to_c13_ratio
      n14_to_n15_ratio = snap_n14_to_n15_ratio
      o16_to_o17_ratio = snap_o16_to_o17_ratio
      o16_to_o18_ratio = snap_o16_to_o18_ratio
      target_metal_fraction = snap_target_metal_fraction
      initial_h2_fraction = snap_initial_h2_fraction
      initial_he3_fraction = snap_initial_he3_fraction
      initial_li6_fraction = snap_initial_li6_fraction
      initial_li7_fraction = snap_initial_li7_fraction
      initial_be9_fraction = snap_initial_be9_fraction
      initial_b10_fraction = snap_initial_b10_fraction
      initial_b11_fraction = snap_initial_b11_fraction
      fo = snap_fo
      difad_velocity_scale = snap_difad_velocity_scale
      mixing_velocity_scale = snap_mixing_velocity_scale
      es_velocity_scale = snap_es_velocity_scale
      gsf_velocity_scale = snap_gsf_velocity_scale
      mu_gradient_scale = snap_mu_gradient_scale
      secular_shear_velocity_scale = snap_secular_shear_velocity_scale
      critical_reynolds = snap_critical_reynolds
      ies = snap_ies
      imu = snap_imu
      es_mixing_scale = snap_es_mixing_scale
      secular_shear_mixing_scale = snap_secular_shear_mixing_scale
      gsf_mixing_scale = snap_gsf_mixing_scale
      gsf_inhibition_mode = snap_gsf_inhibition_mode
      min_abundance = snap_min_abundance
      absolute_tolerance = snap_absolute_tolerance
      relative_tolerance = snap_relative_tolerance
      max_burn_iterations = snap_max_burn_iterations
      opal95_table_unit = snap_opal95_table_unit
      po_weight_l = snap_po_weight_l
      po_weight_teff = snap_po_weight_teff
      po_weight_age = snap_po_weight_age
      po_max_len_sq = snap_po_max_len_sq
      po_output_enabled = snap_po_output_enabled
      track_file_version = snap_track_file_version
      extend_core_inward = snap_extend_core_inward
      num_core_shells_added = snap_num_core_shells_added
      core_mass_reduction_factor = snap_core_mass_reduction_factor
      lrwsh_placeholder = snap_lrwsh_placeholder
      isochrone_output_active = snap_isochrone_output_active
      isochrone_file_unit = snap_isochrone_file_unit
      initial_x_array = snap_initial_x_array
      initial_z_array = snap_initial_z_array
      mixing_length_array = snap_mixing_length_array
      senv0_array = snap_senv0_array
      has_senv0_array = snap_has_senv0_array
      luminosity_tolerance = snap_luminosity_tolerance
      radius_tolerance = snap_radius_tolerance
      zx_tolerance = snap_zx_tolerance
      target_solar_zx = snap_target_solar_zx
      target_solar_age = snap_target_solar_age
      calibrate_solar_model = snap_calibrate_solar_model
      calibrate_solar_zx = snap_calibrate_solar_zx
      rsclzc = snap_rsclzc
      rsclzm1 = snap_rsclzm1
      rsclzm2 = snap_rsclzm2
      iolaol2 = snap_iolaol2
      ioopal2 = snap_ioopal2
      nk = snap_nk
      use_z_ramp = snap_use_z_ramp
      target_luminosity_lsun = snap_target_luminosity_lsun
      target_star_luminosity_tolerance = snap_target_star_luminosity_tolerance
      target_teff = snap_target_teff
      target_radius_rsun = snap_target_radius_rsun
      log_l_prev_model = snap_log_l_prev_model
      log_r_prev_model = snap_log_r_prev_model
      age_at_target_radius = snap_age_at_target_radius
      log_l_at_target_radius = snap_log_l_at_target_radius
      log_l_at_target_radius_prev_run = snap_log_l_at_target_radius_prev_run
      age_prev_model = snap_age_prev_model
      star_found_flag = snap_star_found_flag
      specify_teff_flag = snap_specify_teff_flag
      just_passed_target_radius_flag = snap_just_passed_target_radius_flag
      calibrate_star_flag = snap_calibrate_star_flag
      use_opal95_eos = snap_use_opal95_eos
      use_opal2001_eos = snap_use_opal2001_eos
      use_opal2006_eos = snap_use_opal2006_eos
      use_numerical_derivatives = snap_use_numerical_derivatives
      iopale = snap_iopale
      laol_table_z1 = snap_laol_table_z1
      laol_table_z2 = snap_laol_table_z2
      opal_table_z1 = snap_opal_table_z1
      opal_table_z2 = snap_opal_table_z2
      opal95_single_table_z = snap_opal95_single_table_z
      alex_table_z1 = snap_alex_table_z1
      kurucz_table_z1 = snap_kurucz_table_z1
      kurucz_table_z2 = snap_kurucz_table_z2
      molecular_opacity_logt_min = snap_molecular_opacity_logt_min
      molecular_opacity_logt_max = snap_molecular_opacity_logt_max
      use_alex06_tables = snap_use_alex06_tables
      use_laol89_tables = snap_use_laol89_tables
      use_opal92_tables = snap_use_opal92_tables
      use_opal95_tables = snap_use_opal95_tables
      use_kurucz90_tables = snap_use_kurucz90_tables
      use_alex95_tables = snap_use_alex95_tables
      use_two_z_tables = snap_use_two_z_tables
      ikur2 = snap_ikur2
      icondopacp = snap_icondopacp
      use_conductive_opacity = snap_use_conductive_opacity
      alex95_table_unit = snap_alex95_table_unit
      alex06_table_unit = snap_alex06_table_unit
      alex_mixture_x = snap_alex_mixture_x
      alex_mixture_z = snap_alex_mixture_z
      vfc = snap_vfc
      lvfc = snap_lvfc
      use_diffusion_advection_transport = snap_use_diffusion_advection_transport
      no_am_transport_in_core = snap_no_am_transport_in_core
      scv_h_unit = snap_scv_h_unit
      scv_he_unit = snap_scv_he_unit
      scv_z_unit = snap_scv_z_unit
      allard_target_feh = snap_allard_target_feh
      allard_target_alpha = snap_allard_target_alpha
      allard_use_tau100 = snap_allard_use_tau100
      allard_table_unit = snap_allard_table_unit
      force_solid_body_rotation = snap_force_solid_body_rotation
      solid_body_mode_flag = snap_solid_body_mode_flag
      cstmixing = snap_cstmixing
      cstdiffmix = snap_cstdiffmix
      taucz_placeholder = snap_taucz_placeholder
      deladj_placeholder = snap_deladj_placeholder
      tauhe_placeholder = snap_tauhe_placeholder
      tnorm_placeholder = snap_tnorm_placeholder
      tcz_placeholder = snap_tcz_placeholder
      whe_placeholder = snap_whe_placeholder
      acatmr_placeholder = snap_acatmr_placeholder
      acatmd_placeholder = snap_acatmd_placeholder
      acatmp_placeholder = snap_acatmp_placeholder
      acatmt_placeholder = snap_acatmt_placeholder
      tatmos_placeholder = snap_tatmos_placeholder
      iacat_placeholder = snap_iacat_placeholder
      ijlast_placeholder = snap_ijlast_placeholder
      laoly_placeholder = snap_laoly_placeholder
      ijvs_placeholder = snap_ijvs_placeholder
      ijent_placeholder = snap_ijent_placeholder
      ijdel_placeholder = snap_ijdel_placeholder
      compute_acoustic_depth = snap_compute_acoustic_depth
      ageout_placeholder = snap_ageout_placeholder
      lclcd_placeholder = snap_lclcd_placeholder
      ljlast_placeholder = snap_ljlast_placeholder
      ljwrt_placeholder = snap_ljwrt_placeholder
      acoustic_depth_output = snap_acoustic_depth_output
      iclcd_placeholder = snap_iclcd_placeholder
      use_envelope_triangle_dt = snap_use_envelope_triangle_dt
      pmm_exponent_a = snap_pmm_exponent_a
      pmm_exponent_b = snap_pmm_exponent_b
      pmm_exponent_c = snap_pmm_exponent_c
      pmm_exponent_d = snap_pmm_exponent_d
      pmm_exponent_m = snap_pmm_exponent_m
      pmm_norm_jdot = snap_pmm_norm_jdot
      pmm_norm_mdot = snap_pmm_norm_mdot
      pmm_solar_pressure = snap_pmm_solar_pressure
      pmm_solar_omega = snap_pmm_solar_omega
      pmm_solar_turnover_timescale = snap_pmm_solar_turnover_timescale
      use_pmm_wind_law = snap_use_pmm_wind_law
      scale_by_rossby_number = snap_scale_by_rossby_number
      scale_by_b_field = snap_scale_by_b_field
      wind_law_name = snap_wind_law_name
      exmd = snap_exmd
      extau = snap_extau
      exr = snap_exr
      exm = snap_exm
      exl = snap_exl
      expr = snap_expr
      constfactor = snap_constfactor
      structfactor = snap_structfactor
      excen = snap_excen
      c_2 = snap_c_2
      wind_law_omega_exponent = snap_wind_law_omega_exponent
      wind_saturation_omega = snap_wind_saturation_omega
      ljdot0 = snap_ljdot0
      constant_background_diffusion_coeff = snap_constant_background_diffusion_coeff
      use_constant_background_diffusion = snap_use_constant_background_diffusion
      flag_dx = snap_flag_dx
      flag_dw = snap_flag_dw
      flag_dz = snap_flag_dz
      time_core_min = snap_time_core_min
      time_dl = snap_time_dl
      time_dp = snap_time_dp
      time_dr = snap_time_dr
      time_dt = snap_time_dt
      time_dw_global = snap_time_dw_global
      time_dw_mix = snap_time_dw_mix
      time_dx_core_frac = snap_time_dx_core_frac
      time_dx_core_tot = snap_time_dx_core_tot
      time_dx_shell = snap_time_dx_shell
      time_dx_total = snap_time_dx_total
      time_dy_core_frac = snap_time_dy_core_frac
      time_dy_core_tot = snap_time_dy_core_tot
      time_dy_shell = snap_time_dy_shell
      time_dy_total = snap_time_dy_total
      time_max_dt_frac = snap_time_max_dt_frac
      tol_czbase_fine_width = snap_tol_czbase_fine_width
      tol_dl_max = snap_tol_dl_max
      tol_dm_max = snap_tol_dm_max
      tol_dm_min = snap_tol_dm_min
      tol_dp_core_max = snap_tol_dp_core_max
      tol_dp_czbase_max = snap_tol_dp_czbase_max
      tol_dp_env_max = snap_tol_dp_env_max
      tol_dx_max = snap_tol_dx_max
      tol_dz_max = snap_tol_dz_max
      lstruct_time = snap_lstruct_time
      lnewvars = snap_lnewvars
      s0_pp = snap_s0_pp
      s0_he3he3 = snap_s0_he3he3
      s0_he3he4 = snap_s0_he3he4
      s0_p_c12 = snap_s0_p_c12
      s0_p_c13 = snap_s0_p_c13
      s0_p_n14 = snap_s0_p_n14
      s0_p_o16 = snap_s0_p_o16
      s0_be7_electron = snap_s0_be7_electron
      s0_be7_p = snap_s0_be7_p
      s0_n15_p_c12_branch = snap_s0_n15_p_c12_branch
      s0_n15_p_o16_branch = snap_s0_n15_p_o16_branch
      s0p_pp = snap_s0p_pp
      s0p_he3he3 = snap_s0p_he3he3
      s0p_he3he4 = snap_s0p_he3he4
      s0p_p_c12 = snap_s0p_p_c12
      s0p_p_c13 = snap_s0p_p_c13
      s0p_p_n14 = snap_s0p_p_n14
      s0p_p_o16 = snap_s0p_p_o16
      s0pp_p_c12 = snap_s0pp_p_c12
      s0pp_p_c13 = snap_s0pp_p_c13
      s0pp_p_o16 = snap_s0pp_p_o16
      s0p_be7_p = snap_s0p_be7_p
      s0pp_be7_p = snap_s0pp_be7_p
      s0_pep = snap_s0_pep
      s0_hep = snap_s0_hep
      ln10 = snap_ln10
      clni = snap_clni
      c4pi = snap_c4pi
      c4pil = snap_c4pil
      c4pi3l = snap_c4pi3l
      cc13 = snap_cc13
      cc23 = snap_cc23
      cpi = snap_cpi
      gas_constant = snap_gas_constant
      radiation_constant_over_3 = snap_radiation_constant_over_3
      ca3l = snap_ca3l
      csig = snap_csig
      csigl = snap_csigl
      cgl = snap_cgl
      cmkh = snap_cmkh
      cmkhn = snap_cmkhn
      cdelrl = snap_cdelrl
      cmixl2 = snap_cmixl2
      cmixl3 = snap_cmixl3
      clndp = snap_clndp
      seconds_per_year = snap_seconds_per_year
      cmixl = snap_cmixl
      yrec_version_string = snap_yrec_version_string
      git_hash_string = snap_git_hash_string
end subroutine controls_restore

end module controls_reset_lib

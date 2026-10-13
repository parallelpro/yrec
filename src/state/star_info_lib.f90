!----------------------------------------------------------------------
! star_info_lib
!----------------------------------------------------------------------
! Added 2026 (phase four, step 3 -- docs/ROADMAP.md "Phase four: the star
! layer"). The single owned representation of the stellar model,
! adapted from MESA's star_info to YREC's single-star-per-process
! design (phase-one decision: one module-level instance, no handles --
! the same pattern as every state/ module).
!
! Fields are exactly the model arrays that used to be declared in
! `program main` and threaded through 20-60-argument positional
! signatures (crrect took 60, starin 50, midmod/wrtout 43, hpoint 35,
! getw 29, mix 21) -- the measured root cause of the star-layer
! entanglement, and the source of all three intent bugs phases 2-3
! surfaced. Names and shapes are main's own, unchanged.
!
! What is deliberately NOT here:
!  * run control (timesteps, iteration counters, convergence
!    tolerances): driver policy, stays in arguments by design.
!  * physics-domain state: eos/kap/atm/nuclear keep their table state
!    in their own modules and take per-point inputs as plain
!    arguments. They still READ their configuration from star (about
!    60 star%ctrl/star%job members plus the rate scales, envelope X/Z
!    and use_two_z_tables) -- no per-zone model state. Passing each
!    domain its own controls instead is planned
!    (docs/TARGET_ARCHITECTURE.md section 4.6).
! (The model-level scalars and the former-COMMON state/ modules that
! earlier revisions of this note listed as pending have since been
! absorbed -- see the "model-level scalars" block below and the
! "former common/..." markers throughout the type.)
!
! Mid-timestep working copies (the *_mid arrays getw builds for
! seculr/mid_timestep_model) are working state of the rotation pipeline, not the
! model; routines that receive those keep explicit arguments.
module star_info_lib
! the slots of the composition and neutrino-flux vectors, re-exported
! (they live in constants-only modules so that the physics domains can
! share them without using this module)
      use composition_layout
      use envelope_points, only: envelope_structure_state, &
           atmosphere_structure_state, atmosphere_point_state
      use stitched_model_state, only: stitched_model_record
      use neutrino_flux_layout
      implicit none
      private
      public :: json

! ---- named indices for the positional array slots (2026) ----
! The composition array xa(15,json) and its species counts:
! composition_layout.
      public :: i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, i_n14, i_n15, &
           i_o16, i_o17, i_o18, i_h2, i_li6, i_li7, i_be9, n_species_basic, &
           n_species_extended

! The 12-species surface mixture slots: composition_layout.
      public :: ix_na, ix_al, ix_mg, ix_fe, ix_si, ix_c, ix_h, ix_o, ix_n, &
           ix_ar, ix_ne, ix_he, n_mix_species

! Monte-Carlo sample capacity of the job%s11_rate .. job%age_target
! arrays.
      integer, parameter, public :: max_mc_runs = 1000

! Capacity of the convective (mixed) and radiative zone-bound lists
! mixed_zone_bounds(max_convective_zones,2) and
! radiative_zone_bounds(max_radiative_zones,2) (a model of n
! convective zones has at most n+1 radiative ones).
      integer, parameter, public :: &
           max_convective_zones = 12, max_radiative_zones = 13

! historical gradient-row indices: DEL_GRAD(3,json) is now the three
! flat arrays gradr/gradT/grada; these constants remain for the
! other 3-row gradient carriers (current_gradients, env_gradients).
      integer, parameter, public :: &
           i_grad_rad = 1, i_grad_actual = 2, i_grad_ad = 3

! eps_channels(7,json): specific energy generation channels [erg/g/s].
      integer, parameter, public :: &
           i_eps_pp1 = 1, i_eps_pp2 = 2, i_eps_pp3 = 3, &
           i_eps_cno = 4, i_eps_he3 = 5, i_eps_neu = 6, &
           i_eps_grav = 7

! luminosity_breakdown(n_lum_channels): integrated luminosity per
! channel [Lsun].
      integer, parameter, public :: &
           i_lum_pp1 = 1, i_lum_pp2 = 2, i_lum_pp3 = 3, &
           i_lum_cno = 4, i_lum_3alpha = 5, i_lum_neu = 6, &
           i_lum_grav = 7, i_lum_he_c = 8, n_lum_channels = 8

! The neutrino-flux slots: neutrino_flux_layout.
      public :: i_nu_pp, i_nu_pep, i_nu_hep, i_nu_be7, i_nu_b8, i_nu_n13, &
           i_nu_o15, i_nu_f17, i_nu_he3he3, i_nu_he3he4, n_nu_fluxes

! 2026 state consolidation: the sub-struct TYPE definitions moved
! here from their former one-type-per-file modules in state/ --
! the instances have lived inside star_info since the phase-6 and
! MESA-convention folds; this completes the merge (one root, one
! file). Physics-domain table state stays with its domain.



! (2026: the pulse q* values moved out of star; they are the envelope
! integration's point_physics, core/envelope/envelope_integrands.f90.)

! ---- from state/star_job_lib.f90 ----
! 2026 controls->star% campaign, phase B: the authoritative home of
! every namelist control. Component list GENERATED from the read
! buffer's declarations (io/controls_buffer.f90) by
! tools/gen_controls_state.py -- regenerate on any member change.
! Every component is default-initialized, so
!   star%ctrl = controls_state()
! is the structural "reset controls to pristine defaults" (used by
! read_controls before every read; replaces the retired
! controls_reset_lib snapshot machinery). Immutable after
! read_controls stores into it; consumers migrate from the buffer's
! bare names to star%ctrl%... in phase C.
      type, public :: controls_state
            include 'controls_state_def.inc'
      end type controls_state

      type, public :: star_job
            character(len=256) :: alex06_table_path, allard_table_path, &
                 atm_table_path, fermi_table_path, kurucz_table_path, &
                 kurucz_table2_path, laol_table_path, laol_table2_path, &
                 opal95_table_path, opal92_table_path
            character(len=256) :: zams_a_table_path, zams_b_table_path, &
                 zams_c_table_path, centre1_table_path, centre2_table_path, &
                 centre3_table_path, centre4_table_path, centre5_table_path
            character(len=256) :: opal92_table2_path, pure_z_table_path, &
                 scv_h_table_path, scv_he_table_path, scv_z_table_path, &
                 conductive_table_path, opal_eos_table_path
            character(len=256) :: alex95_table_paths(7)
            double precision :: mixture_weights(n_mix_species)
! phase C flattening: the Monte-Carlo sample arrays (former
! common/monte2/), read from the MC input file by star_setup -- job
! configuration, moved here from the old run_diagnostics grab-bag.
            double precision :: s11_rate(max_mc_runs), s33_rate(max_mc_runs), &
                 s34_rate(max_mc_runs), s17_rate(max_mc_runs), &
                 metal_to_h_ratio(max_mc_runs), &
                 helium_fraction_param(max_mc_runs), &
                 diffusion_factor(max_mc_runs), &
                 luminosity_target(max_mc_runs), age_target(max_mc_runs)
! 2026 phase A batch 6: the current kind-card index (former
! common/zramp/ NK), the run list's cursor -- set by run_yrec's card
! loop, read broadly (evolve_step, the io writers, the calibration
! protocol arithmetic). A structure component cannot be a
! DO-variable, so run_yrec drives it from a local index and
! preserves the historical post-loop value (num_runs+1 on
! exhaustion) explicitly.
            integer :: nk
! phase C batch 2: the run-list / calibration-protocol card arrays
! and latches -- namelist/card-read (parmin fills the buffer, the
! generated sync stores them HERE, not into star%ctrl) but mutable:
! setcal/chkcal/setscal/chkscal rewrite next-cycle cards, the
! stop-disarm pass negates thresholds, the MC loop scales ages.
! GENERATED component list; regenerate via tools/gen_controls_state.py.
            include 'job_controls_def.inc'
      end type star_job

      type, public :: star_info
! structure / thermodynamics, per zone
            double precision :: log_mass(json), luminosity_lsun(json), &
                 logR(json), logP(json), &
                 logT(json), logRho(json)
            logical :: convective_flag(json)
            double precision :: xa(n_species_extended,json), m(json), &
                 dm(json), eps_grav(json)
! rotation corrections to the structure equations
            double precision :: fp_rot(json), &
                 ft_rot(json)
! Henyey solver work arrays and corrections
            double precision :: elim_coeff(4,2,json), elim_rhs(4,json)
            double precision :: log_pressure_delta(json), &
                 log_temperature_delta(json)
            double precision :: max_residual(4)
            integer :: max_correction_index(4)
! surface / envelope fit
            double precision :: surface_bc(6), stored_envelope_state(4), &
                 envelope_fit_coeffs(9), luminosity_breakdown(n_lum_channels)
            double precision :: triangle_log_l(3), &
                 triangle_log_teff(3), fit_point_pressure(3), &
                 fit_point_temperature(3), fit_point_radius(3)
! rotation state, per zone
            double precision :: omega(json), i_rot(json), &
                 j_rot(json), kinetic_energy_rot(json), &
                 kinetic_energy_rot_old(json), mean_radius(json), &
                 eta_squared(json), qiw(json), mean_gravity(json)
! omega_mid: omega (rad/s) of the latest rotation sub-step's model, the
! structure interpolated to the sub-step's end (rot_scr's *_mid model).
! mid_timestep_model, secular_transport and the wind write it;
! evolve_rotation copies it into omega at the end of the step. It
! outlives the step: the next step's first solid-body sub-step
! (rotate_midpoint_model) starts its Newton iteration from its mean.
            double precision :: omega_mid(json)
! am_transport_convective_flag (legacy LCZ): T where a zone is
! convective for angular-momentum purposes (overshoot included), from
! am_convective_regions on the current grid. Written by
! read_starting_model and rezone, i.e. before every structure solve;
! read by henyey_iterate's omega_from_j, which enforces solid-body
! rotation over those zones during the Newton iterations. Never
! assigned before 2026-10 (inherited from F77: main.f's LCZ was a
! local passed only to CRRECT/GETROT), so the iterations treated every
! zone as radiative.
            logical :: am_transport_convective_flag(json)
! nuclear / neutrino diagnostics, per zone
            double precision :: be7_mass_fraction_zone(json), &
                 neutrino_flux_zone(n_nu_fluxes,json)
! model-level scalars (2026, phase four follow-on): absorbed after the
! per-call-site audit the star_info header called for -- every
! converted routine receives exactly main's storage for these at every
! call site. NOT absorbed, still arguments, because call sites pass
! separate storage: composition_step's core_cz_edge/envelope_cz_edge
! (its re-runs at levels >= 3 pass their own locals) and its
! mixed_zone_bounds_no_overshoot (same reason, see step 3).
            integer :: nz, model_number
! Newton iterations the last converged model took (set by
! henyey_iterate at convergence; shown in the run-log model line).
            integer :: newton_iterations = 0
            integer :: core_cz_top_index, envelope_cz_bottom_index
            double precision :: log_total_mass, star_mass
            double precision :: log_Teff, log_L
! 2026 controls->star% campaign, phase A: computed/working state
! evicted from controls_buffer (they were never namelist values). Per
! the agreed shape these land FLAT on star%. atm_hras is set once by
! setups (Krishna-Swamy T(tau) at tau=2/3); tenv once by parmin
! (0.5*(tenv0+tenv1)); atm_choice_initial once by parmin (from
! kttau); use_ttau_relation is toggled at runtime by surfbc/envint
! (genuinely mutable working state -- the reason it cannot be a
! control).
            double precision :: atm_hras, tenv
            integer :: atm_choice_initial
            logical :: use_ttau_relation
! phase A batch 2: the nuclear cross-section scale family -- fully
! recomputed by setup/user_inputs.f90 from the s0_* namelist inputs after
! every controls read (never namelist values themselves), and
! cross_section_scale(1:3,16) is overwritten per Monte-Carlo run by
! apply_monte_carlo_parameters. Consumed by net_lib's compute_network_rates and
! burn_lib's compute_zone_energy.
            double precision :: cross_section_scale(17), qs0e_scale(8), &
                 qqs0ee_scale(8), o16_gamma_scale, c12_alpha_scale
! phase A batch 4: the working mixing length (former CMIXL, the one
! mutable member physical_constants carried) -- copied from the per-card
! control mixing_length_array(nk) at every kind-card start
! (run_yrec's begin_kind_card) and rewritten by solar calibration.
! MESA name; the 1.4d0 declaration default reproduces the historical
! parmin-echo value before the first kind card (and after a
! yrec_reset star0 restore).
            double precision :: mixing_length_alpha = 1.4d0
! phase A batch 5: the solar-unit octet (former common/const/).
! Computed at startup by core/run/run_setup.f90 -- solar_luminosity_cgs/
! solar_radius_cgs seeded from the NAMELIST /physics/ clsun/crsun
! values (copied by parmin), the rest derived from them -- and the
! luminosity trio is overwritten per Monte-Carlo run by
! apply_monte_carlo_parameters. Job-configured unit definitions,
! i.e. computed state, not controls. ln10_over_lsun is
! ln(10)/solar_luminosity_cgs (it is NOT a logarithm of Lsun): times
! the enclosed mass it scales the luminosity-equation terms in
! henyey_coefficients.
            double precision :: solar_luminosity_cgs, &
                 log10_solar_luminosity, ln10_over_lsun, &
                 solar_mass_cgs, solar_radius_cgs, &
                 log10_solar_radius, solar_bolometric_magnitude
! phase C batch 2: the rotation-mixing velocity-factor work array
! (former common/varfc/ VFC, written per zone by getfc/seculr every
! step) -- never a control, evicted from the buffer.
            double precision :: vfc(json)
! 2026 solver-scratch cleanup: the rot/circ members with PROPERTY
! evidence (read by the output writers or across domains) flatten
! here; the true solver workspace moved to
! rotation/rotation_scratch_lib.f90.
            double precision :: valfmlt(json), vphmlt(json), vcmxmlt(json)
            double precision :: es_circulation_velocity(json), &
                 secular_shear_velocity(json), gsf_circulation_velocity(json)
! phase C batch 3: working/diagnostic state evicted from the buffer
! (zero read-path references -- never namelist targets). The
! acoustic-depth family holds calcad's OUTPUTS (read by wrtout's
! track columns); the age/log_l/log_r group is chkscal's
! star-calibration bookkeeping (star_found_flag arms the final
! age-stopped run); compute_acoustic_depth is the calcad toggle the
! calibration verdict disables; iov1/iov2 are the overshoot zone
! indices find_convection_zones returns; use_two_z_tables is derived by
! setupopac from the table configuration; disk_lifetime is the
! disk-locking countdown evolve_step advances; pulsation_mass_msun
! is stamped per kind card by begin_kind_card.
            double precision :: age_at_target_radius, age_prev_model, &
                 log_l_at_target_radius, log_l_at_target_radius_prev_run, &
                 log_l_prev_model, r_rsun_prev_model
! run-log bookkeeping: why the (last) kind card ended -- the default
! set per card by init_stop_conditions, overridden by whichever stop
! fires; read by run_log's end-of-run summary.
            character(len=96) :: termination_reason = 'model budget exhausted'
! run-log bookkeeping: the settling routines have reported that
! settling is suspended (a fully convective model, or one without
! hydrogen or helium to settle); the message prints once per suspension,
! and "SETTLING RESUMED" clears it.
            logical :: settling_suspended_reported = .false.
! the mass change's state (core/mass_change; henyey_coefficients reads
! the structure shifts): the envelope's mean specific entropy, and the
! shifts of log P and log T accretion makes at the surface
            double precision :: envelope_specific_entropy, &
                 delta_log_pressure, delta_log_temperature
! (commented out with mass_change's rotation-scaled solar wind,
! 2026-10-12: its rate, reference and maximum omega, and its switch)
!            double precision :: solar_wind_mass_loss_rate_msun_yr, &
!                 wind_reference_omega, wind_max_omega
!            logical :: use_rotation_scaled_solar_wind
            logical :: star_found_flag = .false., &
                 just_passed_target_radius_flag = .false., &
                 use_two_z_tables = .false.
            double precision :: disk_gate_age_gyr, pulsation_mass_msun
            integer :: iov1, iov2
! mixed/radiative zone bookkeeping
            integer :: mixed_zone_bounds(max_convective_zones,2), &
                 mixed_zone_bounds_no_overshoot(max_convective_zones,2), &
                 radiative_zone_bounds(max_radiative_zones,2)
! 2026 (phase four, step 4): the former-COMMON model-state modules,
! folded in as components -- their types stay defined in their own
! state/ files, the single instances now live here. prev is the
! previous-model store (former oldmod_lib prev_model), diag the
! per-shell diagnostics (former scrtch_lib shell_diag), run the
! run-level diagnostics (former run_diag_lib run_diag), rot the
! rotation/diffusion working state (former rotdiff_lib rot_diff).
! Physics-domain files (atm, nuclear, wind internals) that read these
! now visibly reference star%... -- the remaining star-coupling inside
! physics domains is grep-able as star% under those directories.
! ---- 2026 phase-C flattening: the former prev/diag/run/turnover/
! flux/compute_zone_energy/env_comp/thermo/light_burn/evo sub-structs dissolved to
! flat members (MESA shape: state is flat, only job/ctrl input
! bundles and the rotation-solver scratch structs stay nested).
! Member declarations moved verbatim, per-group comments retained.
! -- former prev_model_state (prev (model snapshot)) --
! 2026 MESA-convention pass: start-of-step values carry MESA's _start
! suffix, mirroring the live star_info member they snapshot.
            double precision :: logP_start(json), logT_start(json), &
                 logR_start(json), luminosity_lsun_start(json), logRho_start(json)
            double precision :: xa_start(n_species_extended,json)
            logical :: convective_flag_start(json)
            double precision :: log_Teff_start
! The start-of-step surface luminosity (Lsun), L at the outer point
! before the rezone, for the HR-motion timestep limit.
            double precision :: surface_luminosity_lsun_start
! -- former shell_diagnostics_state (per-zone diagnostics), member
! micro-renames 2026: eps_total (SESUM, total specific energy
! generation), eps_channels (SEG, per-channel via i_eps_*), beta
! (SBETA, gas-pressure fraction), eta (SETA, electron degeneracy),
! conv_vel (SVEL, MESA name), opacity_zone (SO -- per-zone OPACITY,
! not O16: the 2026 o16_zone guess was wrong),
! fxion_zone (SFXION, ionization fractions), gradr/gradT/grada (the
! DEL_GRAD(3,:) rows, split into MESA-named arrays). scp keeps its
! name deliberately: it is the Henyey-solve-time specific heat that
! the profile writers read, distinct in fill time from cp (the
! rotation-pipeline copy) -- merging or renaming awaits a semantics
! audit.
            double precision :: eps_total(json), eps_channels(7,json), &
                 beta(json), eta(json)
! 2026 retire-legacy (.FULL): the model-grid rotational transport
! coefficients [cm^2/s], stored by secular_transport for the
! D_omega / D_mix_rot profile columns (zero when rotation is off;
! last substep of the step wins).
            double precision :: am_diffusion_coeff(json), &
                 mixing_diffusion_coeff(json)
            double precision :: opacity_zone(json), gradr(json), gradT(json), &
                 grada(json), fxion_zone(3,json), conv_vel(json), scp(json)
! -- former turnover_state (turnover) --
            double precision :: convective_turnover_timescale, &
                 convective_turnover_timescale_old
            double precision :: pphot, pphot0
            double precision :: fracstep
! -- former light_element_burn_state (light_burn) --
! former common/newrat/: Li6/Li7/Be9 burning rates at the end of the
! timestep, at the current (possibly overshoot-adjusted) depth.
            double precision :: rate_li6(json), rate_li7(json), &
                 rate_be9(json)
! former common/oldrat/: same, at the start of the timestep.
            double precision :: rate_li6_start(json), rate_li7_start(json), &
                 rate_be9_start(json)
! former common/liov/: pressure scale heights used to search downward
! from the CZ base for the true (overshoot-corrected) base location.
            double precision :: envelope_overshoot_depth_start, &
                 envelope_overshoot_depth_end
! former common/prevcz/: previous end-of-timestep values, used as the
! new beginning-of-timestep values.
            double precision :: cz_base_radius_prev
            integer :: envelope_cz_base_zone_prev
! former common/deuter/: deuterium burning rate per s (current/start
! of timestep; per Gyr until 2026-10), and the convection-zone base
! zone index used by the deuterium-burning routines.
            double precision :: deuterium_burning_rate(json), &
                 deuterium_burning_rate_start(json)
! accreted_mass: the mass accreted during the step, g (wind/mdot;
! named accreted_mass_fraction until 2026-10 despite being a mass)
            double precision :: accreted_mass
            integer :: jcz
! He3+He3 and He3+He3 + He3+He4 luminosity per shell (JVS 10/11;
! former common/grab/)
            double precision :: he3_he3_luminosity_zone(json), &
                 he3_burning_luminosity_zone(json)
! -- former neutrino_flux_state (flux) --
            double precision :: neutrino_flux_total(n_nu_fluxes)
            double precision :: cl37_snu_rate, ga71_snu_rate
! -- former envelope_composition_state (env_comp) --
            double precision :: envelope_hydrogen_fraction, &
                 envelope_metal_fraction
            double precision :: zenvm, amuenv, fxenv(n_mix_species)
            double precision :: xnew, znew, stotal, senv
! -- former shell_temp_state (thermo) --
! mu = MESA name (was mean_molecular_weight)
            double precision :: cp(json), mu(json)
            double precision :: qdt(json), thdif(json), visc(json)
! -- former run_diagnostics_state (run (observables + driver bookkeeping)) --
! 2026 MESA-style output: per-model history sources. Computed by
! compute_observables after all physics for the step is done;
! writers (write_history) only read them. Zero at run start
! (star0 snapshot / static zero), refreshed every converged model.
           double precision :: log_R_surface, log_g_surface
! asteroseismic observables (2026): nu_max from (log g, Teff);
! delta_nu_rho [uHz] from the mean density (a scaling estimate,
! hence the _rho suffix); delta_nu [uHz], the asymptotic p-mode
! large separation, from the sound-travel-time integral; delta_Pg
! [s], the asymptotic l=1 g-mode period spacing, from the
! Brunt-Vaisala integral -- computed by observables's seismic
! theme, consumed by the stop conditions and the optional history
! columns
           double precision :: nu_max, delta_nu_rho
           double precision :: delta_nu, delta_Pg
           double precision :: total_moment_of_inertia, cz_moment_of_inertia
           double precision :: rotation_period_days, surf_velocity_kms
           double precision :: h_shell_bot_mass, h_shell_mid_mass, &
                h_shell_top_mass, h_shell_bot_radius, h_shell_mid_radius, &
                h_shell_top_radius
! former common/entrop/
           double precision :: dlog_temperature_dt(json), &
                dlog_pressure_dt(json), &
                dln_luminosity_dt(json), &
                dlog_radius_dt(json)
! former common/theage/
           double precision :: age_gyr
! former common/calsun/
           double precision :: dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, &
                log_l_prev, log_r_prev, delta_x, delta_alpha
           logical :: solar_calibration_converged
! former common/sound/
           double precision :: adiabatic_index_gamma1(json)
! former common/cent/
           double precision :: central_log10_temperature, central_log10_pressure, &
                central_log10_density, envelope_mass, envelope_radius
! 2026 (phase four, step 5): observables formerly computed as
! locals inside the retired write_legacy_output, now filled by the star layer
! (core/diagnostics/observables.f90) and only READ by the writers.
           double precision :: central_beta, central_degeneracy_eta
! The Saha-table state of the central EOS call (compute_central_conditions),
! carried from model to model.
           integer :: central_saha_state = 0
           double precision :: core_cz_mass
           double precision :: envelope_cz_temperature, envelope_cz_density, &
                envelope_cz_pressure, envelope_cz_opacity
! former common/origstart/
           double precision :: orig_specific_angular_momentum(json), &
                orig_composition(n_species_extended,json)
! former common/comp2/
           double precision :: envelope_helium_fraction, envelope_he3_fraction
! (former common/envprt/'s current_* point values: now the envelope
! integration's point_physics, core/envelope/envelope_integrands.f90)
! former common/oldrot/
           double precision :: old_omega(json), &
                old_specific_angular_momentum(json), &
                old_hg(json), &
                old_mean_radius(json), old_eta_squared(json)
! The start-of-step physics on the step's mesh, from which the rotation
! sub-steps interpolate in time (mid_timestep_model): the rezone stores
! it (rerun_physics_on_new_grid; eps_total and eps_total + eps_neu +
! eps_grav splined by regrid_start_of_step_eps). The gradients, mu, the
! opacity, c_P, d ln rho/d ln T, the convective velocity, the viscosity,
! the thermal diffusivity and the two energy-generation sums, named
! after the members they copy (opacity: opacity_zone).
           double precision :: old_gradT(json), &
                old_grada(json), old_mu(json), &
                old_opacity(json), old_cp(json), &
                old_qdt(json), old_conv_vel(json), &
                old_visc(json), old_thdif(json), &
                old_eps_total(json), old_gradr(json), &
                old_eps_plus_neu_grav(json)
! former common/i2o/
           character(len=4) :: initial_composition_code
! -- former evolve_state (evo (driver-step state)) --
            logical :: has_h_shell, model_diverged_flag, punch_pending_flag, &
                 recompute_envelope_triangle, reset_triangle
            integer :: h_shell_end_index, h_shell_midpoint_zone, &
                 h_shell_zone_begin, ikut_flag, stored_vertex_index
! chosen_dt (former DELTSH): the timestep compute_timestep adopted
! (the minimum over all its limiters, seconds); read_starting_model
! reuses the slot transiently as |delta_time|. Nothing reads it as
! an H-burning governor.
            double precision :: convective_velocity, dt, &
                 dt_saved, dlnrho_dlnp, dlnrho_dlnt, chosen_dt, &
                 max_domega_frac, timestep_yr, total_angular_momentum, &
                 total_rotational_ke, triangle_orientation
! -- the extended-model pulse physics arrays (former pulse%) --
! pulse_dlnrho_*/pulse_dlnkap_* are logarithmic derivatives;
! pulse_deps_dlnrho/dlnt are ABSOLUTE d eps/dln rho, d eps/dlnT in
! erg/g/s (compute_zone_energy's burn%deps_*; named pulse_dlneps_* until 2026-10).
            double precision :: pulse_dlnrho_dlnp(json), &
                 pulse_deps_dlnrho(json)
            double precision :: pulse_deps_dlnt(json), &
                 pulse_dlnkap_dlnrho(json)
            double precision :: pulse_dlnkap_dlnt(json), &
                 pulse_specific_heat(json)
            double precision :: pulse_mean_molecular_weight(json), &
                 pulse_dlnrho_dlnt(json)
            double precision :: pulse_electron_mean_weight_inverse(json)
! 2026 (phase six, step 1): the nine remaining model-state modules,
! folded in the same way as prev/diag/run/rot -- types stay in their
! own files, the single instances live here. See docs/ROADMAP.md phase six
! for the per-module classification (envstruct/atmstruct deliberately
! stay atm-domain state and are NOT here).
! job configuration and driver-step state (2026 MESA-convention pass:
! folded in as nested sub-structs -- MESA's s% job precedent; one
! root means yrec_reset's star snapshot covers them automatically)
            type(star_job) :: job
! phase B: the controls bundle (see controls_state above). Nested
! like star%job -- the two input bundles are the only nested
! sub-structs in the target shape.
            type(controls_state) :: ctrl
! The outer layers of the latest envelope integration (envelope_points):
! the envelope and the atmosphere point by point, and the atmosphere
! integration's current point. Written by the integration kernel through
! atm_get; read by the stitched model, the envelope refit and rebuild,
! read_starting_model and the mass change.
            type(envelope_structure_state) :: env_points
            type(atmosphere_structure_state) :: atm_points
            type(atmosphere_point_state) :: atm_point
! The stitched model of the converged star (stitched_model_state).
            type(stitched_model_record) :: stitched
      end type star_info

! the one star this process evolves (no handles -- see header)
      type(star_info), public, save :: star


! Module-level (deliberately OUTSIDE the star snapshot/reset):
! evolve_step/output_diag reset flags (set by yrec_reset's
! prologue) and the libyrec namelist-path overrides.
      character(len=256), public, save :: control_nml_override = ' '
      character(len=256), public, save :: physics_nml_override = ' '
      logical, public, save :: evolve_step_reset_pending = .false.

end module star_info_lib

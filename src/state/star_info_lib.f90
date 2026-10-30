!----------------------------------------------------------------------
! star_info_lib
!----------------------------------------------------------------------
! The type of the stellar model, star_info (MESA's name), with the two
! input bundles nested in it (controls_state, star_job) and the
! capacities and slot constants its members are sized and indexed by.
! The module declares no star: run_yrec (core/run/run_driver.f90) owns
! the run's star, and every star-layer procedure takes it as its
! argument s (docs/GUIDELINES.md, "The star object: six rules"). The
! physics domains never use this module; they take the slot constants
! from the constants-only layout modules it re-exports, and their
! configuration from the settings records setup/ builds.
!
! The members of star_info are grouped by lifetime, longest-lived first:
! per run (the input bundles, the constants set at setup), per run card
! (bookkeeping), per model (the model, its rotation state and regions),
! per step (control, the start-of-step snapshot, working state), per
! Newton iteration (the solve, the diagnostics written during it), then
! the diagnostics written once per converged model, and last the members
! one module alone uses. Each
! class opens with a banner,
!   ! ==== <class> ====
! and a paragraph on what it holds, ending with "Valid:" and when its
! values can be trusted. Each member has a comment directly above it,
!   ! <meaning> [<unit>]; <writers>. <notes>
! the unit in brackets ([-] for a dimensionless number or an index;
! [logical], [text]), then the modules that write it (the main ones),
! then optional notes. tools/gen_star_docs.py turns these comments into
! the tables of state/README.md, with every writer and reader taken from
! the code, and fails (--check, run by test_generated_files.py) when a
! member has no comment or names a module that does not write it, or is
! in the wrong class: an earlier value (<member>_start, _prev,
! _iteration_start) belongs in the Step snapshot class, and a member one
! module alone uses in "Private to one module". Keep one member per
! declaration.
!
! A run starts from a static template star (core/run/yrec_reset.f90), so
! every member starts each run at zero, .false. or blank unless its
! declaration gives an initializer.
module star_info_lib
! the slots of the composition, neutrino-flux and transport-mechanism
! vectors, re-exported (they live in constants-only modules so that the
! physics domains can share them without using this module)
      use composition_layout
      use envelope_points, only: envelope_structure_state, &
           atmosphere_structure_state, atmosphere_point_state
      use stitched_model_state, only: stitched_model_record
      use neutrino_flux_layout
      use transport_mechanism_layout
! the sizes of the reaction-rate scales
      use reaction_layout, only: num_reactions, num_extra_scales
      implicit none
      private
      public :: json

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

! The rows of eps_channels(7,json): the chains' fractions of eps_total
! (rows 1-5, dimensionless), eps_nu and eps_grav (rows 6-7, erg/(g s)).
      integer, parameter, public :: &
           i_eps_pp1 = 1, i_eps_pp2 = 2, i_eps_pp3 = 3, &
           i_eps_cno = 4, i_eps_he3 = 5, i_eps_neu = 6, &
           i_eps_grav = 7

! luminosity_breakdown_lsun(n_lum_channels): integrated luminosity per
! channel [Lsun]; i_lum_extra is the extra energy sources'
! (core/structure/energy_sources.f90), zero without any.
      integer, parameter, public :: &
           i_lum_pp1 = 1, i_lum_pp2 = 2, i_lum_pp3 = 3, &
           i_lum_cno = 4, i_lum_3alpha = 5, i_lum_neu = 6, &
           i_lum_grav = 7, i_lum_he_c = 8, i_lum_extra = 9, &
           n_lum_channels = 9

! The transport-mechanism slots: transport_mechanism_layout.
      public :: i_mech_eddington_sweet, i_mech_gsf, i_mech_shear, &
           i_mech_solberg_hoiland, n_transport_mechanisms, transport_mechanism_name
! The neutrino-flux slots: neutrino_flux_layout.
      public :: i_nu_pp, i_nu_pep, i_nu_hep, i_nu_be7, i_nu_b8, i_nu_n13, &
           i_nu_o15, i_nu_f17, i_nu_he3he3, i_nu_he3he4, n_nu_fluxes

! The namelist controls the run reads, by their canonical names
! (star%ctrl; the input-only and unused ones stay in the read's buffer,
! registry home buffer). The component list is GENERATED from the read
! buffer's declarations (io/controls_buffer.f90) and the registry's home
! column by tools/gen_controls_state.py; regenerate on any change. read_controls resets star%ctrl (star%ctrl =
! controls_state()) before every read and stores the buffer into it;
! nothing changes it after the read. The controls are documented in
! src/defaults/controls_registry.md.
      type, public :: controls_state
            include 'controls_state_def.inc'
      end type controls_state

! The run's job (star%job): the controls the run changes
! (job_controls_def.inc, generated from the registry's home column: the
! read sets them, then the calibrations rewrite the next cycle's cards,
! the stop conditions disarm thresholds, the atmosphere switch changes
! atm_choice, the solver ramps fcorr), the Monte-Carlo samples and the
! run-card cursor. The table paths are controls (star%ctrl).
      type, public :: star_job
! the Monte-Carlo samples, one entry per Monte-Carlo run, from the sample
! file: the S-factor scales of p+p, 3He+3He, 3He+4He and 7Be+p, Z/X,
! the helium (diffusion) parameter and its copy diffusion_factor, the
! solar-luminosity and age factors [-]; monte_carlo
! (setup_monte_carlo_runs). Read by apply_monte_carlo_parameters and the
! Monte-Carlo summary.
            double precision :: s11_rate(max_mc_runs) = 0.0d0, &
                 s33_rate(max_mc_runs) = 0.0d0, s34_rate(max_mc_runs) = 0.0d0, &
                 s17_rate(max_mc_runs) = 0.0d0, metal_to_h_ratio(max_mc_runs) = 0.0d0, &
                 helium_fraction_param(max_mc_runs) = 0.0d0, &
                 diffusion_factor(max_mc_runs) = 0.0d0, &
                 luminosity_target(max_mc_runs) = 0.0d0, age_target(max_mc_runs) = 0.0d0
! the current run card's index, the run list's cursor [-]; run_driver
! (run_yrec's card loop runs a local index, as a structure component
! cannot be a DO variable, and leaves num_runs+1 after the last card).
! Read by the step, the writers and the calibrations' cycle arithmetic.
            integer :: nk = 0
! The run-card arrays and latches. GENERATED component list; regenerate
! with tools/gen_controls_state.py.
            include 'job_controls_def.inc'
      end type star_job

      type, public :: star_info
! ==== Input bundles ====================================================
! The controls (ctrl) and the run's job (job), read by read_controls at
! the start of every run_yrec call. Their members are documented in
! src/defaults/controls_registry.md and the *.defaults files, not here.
! Valid: from the controls read on; ctrl is not changed after it, job is
! changed by the run (see the type star_job).

! the namelist controls the run reads, by their canonical names [-];
! controls_reader, controls_sync. Reset to the declared defaults before each read; members:
! src/defaults/controls_registry.md, generated into controls_state_def.inc.
            type(controls_state) :: ctrl
! the run's job: the run-card arrays and latches (the controls the run
! changes), the Monte-Carlo samples and the run-card cursor nk [-];
! controls_sync, run_driver, calibration, monte_carlo. Changed by the
! run: the calibrations rewrite the next cards, the stop conditions
! disarm thresholds, the atmosphere switch changes atm_choice, the
! solver ramps fcorr.
            type(star_job) :: job

! ==== Constants set at setup ===========================================
! Values fixed for the whole run by the controls read, the controls'
! derivation or the setup before the first card: the solar units, the
! nuclear rate scales, a few atmosphere and opacity settings.
! Valid: after read_controls and prepare_run; a Monte-Carlo run rescales
! the solar luminosity and the first rate scales.

! Lsun, a copy of ctrl%solar_luminosity_cgs (legacy clsun) [erg/s];
! run_setup, monte_carlo. A Monte-Carlo run scales it.
            double precision :: solar_luminosity_cgs
! log10 Lsun [log10 erg/s]; run_setup, monte_carlo
            double precision :: log10_solar_luminosity
! ln(10)/Lsun, not a logarithm: times the enclosed mass it scales the
! luminosity-equation terms [s/erg]; run_setup, monte_carlo
            double precision :: ln10_over_lsun
! Msun, a copy of ctrl%solar_mass_cgs [g]; run_setup
            double precision :: solar_mass_cgs
! Rsun, a copy of ctrl%solar_radius_cgs (legacy crsun) [cm]; run_setup
            double precision :: solar_radius_cgs
! log10 Rsun [log10 cm]; run_setup
            double precision :: log10_solar_radius
! per-reaction scale of the cross section S0 against Bahcall & Ulrich
! (1988), by reaction_layout's slot r_* [-]; user_inputs,
! controls_reader, monte_carlo. Those of pp, He3 + He3 and He3 + He4 are
! rescaled per Monte-Carlo run.
            double precision :: cross_section_scale(num_reactions)
! scale of each reaction's S'/S against Solar Fusion I, by reaction slot
! r_* [-]; user_inputs
            double precision :: qs0e_scale(num_reactions)
! scale of each reaction's S''/S against Solar Fusion I, by reaction slot
! r_* [-]; user_inputs
            double precision :: qqs0ee_scale(num_reactions)
! scale of the rates outside the reaction table, by reaction_layout's
! extra slot xs_* (pep, Be7 + e, Be7 + p, hep) [-]; user_inputs,
! controls_reader, monte_carlo. Be7 + p's is rescaled per Monte-Carlo run.
            double precision :: extra_cross_section_scale(num_extra_scales)
! scale of Be7 + p's S'/S against Solar Fusion I [-]; user_inputs
            double precision :: be7p_qs0e_scale
! scale of Be7 + p's S''/S against Solar Fusion I [-]; user_inputs
            double precision :: be7p_qqs0ee_scale
! scale of the N15(p,gamma)O16 branch [-]; user_inputs
            double precision :: o16_gamma_scale
! scale of the N15(p,alpha)C12 branch [-]; user_inputs
            double precision :: c12_alpha_scale
! log T at which an envelope integration stores its state for the surface
! boundary (stored_envelope_state): the mean of ctrl%tenv0 and
! ctrl%tenv1 [log10 K]; derived_controls
            double precision :: envelope_store_log_temperature
! atm_choice as read from the deck [-]; controls_reader. The run may
! change job%atm_choice (see use_ttau_relation); this keeps the original.
            integer :: atm_choice_initial
! load a second opacity-table set at another Z (the Z ramp or Z
! diffusion) [logical]; run_driver
            logical :: use_two_z_tables = .false.

! ==== Run and card bookkeeping =========================================
! Per-run and per-card state: the card's mixing length and mass, the
! atmosphere switch, the surface mixture of the first card, the solar
! and star calibration loops, the reason the run ended.
! Valid: from the card (or run) that sets them.

! mixing-length parameter of the current run card [-]; run_driver. From
! job%mixing_length_array(nk) at each card start (the solar calibration
! changes the later cards' entries); in .mod; 1.4 before the first card.
            double precision :: mixing_length_alpha = 1.4d0
! star_mass at the start of the card, for the profile and pulsation
! headers [Msun]; run_driver
            double precision :: pulsation_mass_msun
! the run has switched atm_choice to the gray atmosphere because Teff left
! the Kurucz or Allard table [logical]; surface_boundary,
! envelope_integration, controls_reader. Switches back when Teff returns.
            logical :: use_ttau_relation
! why the last run card ended, for the run summary [text];
! stop_conditions, run_driver. Set to the default per card by
! init_stop_conditions, overridden by whichever stop fires.
            character(len=96) :: termination_reason = 'model budget exhausted'
! the settling routines have reported that settling is suspended (a fully
! convective model, or no hydrogen or helium to settle) [logical];
! composition_step, mid_timestep_model. The message prints once per
! suspension; "SETTLING RESUMED" clears it.
            logical :: settling_suspended_reported = .false.
! helium settles in this run card [logical]; run_driver, composition_step,
! mid_timestep_model. Armed from the control ctrl%diffuse_helium_active
! (never changed after the read) at the start of every card, and so of
! every Monte-Carlo run, by begin_kind_card; the settling routines switch
! it off once the central X falls below ctrl%hydrogen_diffusion_floor.
            logical :: settling_active = .false.
! surface X of the first run card: the reference X of the EOS mixture and
! of 1/mu_ion in compute_shell_physics [-]; starting_model. Set once
! (update_surface_mixture) and kept for the run on purpose; the controls
! read and kap_init use it earlier, while it is still 0.
            double precision :: envelope_hydrogen_fraction
! surface Z of the EOS mixture: the first card's, or the last Z
! rescale's [-]; starting_model, model_rescaling, run_cards. The controls
! read (interpret_kind_cards) also sets it from a card's rescale Z.
            double precision :: envelope_metal_fraction
! surface Y = 1 - X - Z - X(He3) of the first card [-]; starting_model
            double precision :: envelope_helium_fraction
! surface He3 of the first card [-]; starting_model
            double precision :: envelope_he3_fraction
! the part of the first card's surface Z not in C, N and O, passed to
! convection's settings (envelope_metals_no_cno) [-]; starting_model
            double precision :: envelope_metals_no_cno
! inverse mean atomic weight of the surface mixture, sum X_i/A_i
! [1/amu]; starting_model
            double precision :: envelope_inverse_atomic_weight
! number fractions of the 12 surface-mixture elements, slots ix_*, passed
! to the EOS (eos_set_mixture) [-]; starting_model
            double precision :: envelope_number_fractions(n_mix_species)
! d log L/dX of the last solar-calibration Newton step: calibration's
! parameter dlogl_dx_fixed, -3.78 [-]; calibration. 0 before the first
! step; printed in the Monte-Carlo summary, as the three below.
            double precision :: dlum_dx
! d log R/dX of the last Newton step, dlogr_dx_fixed, -0.89 [-];
! calibration
            double precision :: drad_dx
! d log L/d alpha of the last Newton step, dlogl_dalpha_fixed, 0.0139
! [-]; calibration
            double precision :: dlum_dalpha
! d log R/d alpha of the last Newton step, dlogr_dalpha_fixed, -0.050
! [-]; calibration
            double precision :: drad_dalpha
! the solar calibration has converged [logical]; calibration
            logical :: solar_calibration_converged
! the star calibration found a track through the target L and R, so the
! final age-stopped run is armed [logical]; calibration
            logical :: star_found_flag = .false.
! this model has just crossed the star calibration's target radius
! [logical]; calibration
            logical :: just_passed_target_radius_flag = .false.

! ==== Model ============================================================
! The converged model: the mesh, the structure, the composition, the
! global scalars, and the envelope fit that is its outer boundary
! condition. Persists from step to step; the core of it is written to
! the .mod file (io/model_files.f90; "in .mod" below). Point 1 is the
! innermost point, point nz the fitting point.
! Valid: between steps, and inside a step for the model being solved; the
! Newton corrections overwrite the structure in place, so during a step
! it holds the current iterate.

! number of points of the interior model, out to the fitting point [-];
! rezoning, starting_model, model_files. In .mod.
            integer :: nz
! counter of converged models in the run; the step ages the model only
! when it is >= 0 [-]; step_driver, run_driver, model_files. In .mod;
! reset per card from ctrl%set_initial_model_number.
            integer :: model_number
! age of the model [Gyr]; step_driver, model_files. In .mod.
            double precision :: age_gyr
! total mass M [Msun]; starting_model, mass_change, model_files. In .mod;
! MESA's history name, without the _msun suffix.
            double precision :: star_mass
! log10 of the total mass M [log10 g]; starting_model, mass_change,
! model_files. In .mod.
            double precision :: log_total_mass
! log10 of the luminosity at the fitting point, L(nz) [log10 Lsun];
! henyey_iteration, model_files. In .mod; reset after every Newton
! correction.
            double precision :: log_L_lsun
! log10 effective temperature, from the envelope fit (envelope_fit_coeffs
! 7-9) [log10 K]; henyey_iteration, model_files. In .mod; reset after
! every Newton correction.
            double precision :: log_Teff
! log10 of the mass m inside point k, the mesh coordinate [log10 g];
! rezoning, starting_model, mass_change, model_files. In .mod.
            double precision :: log_mass(json)
! mass inside point k [g]; rezoning, starting_model, mass_change. Rebuilt
! from log_mass after every rezone or read.
            double precision :: m(json)
! mass of the shell around point k, half the span to the neighbouring
! points; dm(nz) reaches the total mass [g]; rezoning, starting_model,
! mass_change
            double precision :: dm(json)
! luminosity L(m) through point k [Lsun]; henyey_iteration, rezoning,
! starting_model, model_files. In .mod; a Newton variable (linear, not
! log).
            double precision :: luminosity_lsun(json)
! log10 radius [log10 cm]; henyey_iteration, rezoning, starting_model,
! model_files. In .mod; a Newton variable.
            double precision :: logR(json)
! log10 pressure [log10 dyn/cm^2]; henyey_iteration, rezoning,
! starting_model, model_files. In .mod; a Newton variable.
            double precision :: logP(json)
! log10 temperature [log10 K]; henyey_iteration, rezoning,
! starting_model, model_files. In .mod; a Newton variable.
            double precision :: logT(json)
! log10 density, from the EOS at (P, T, X, Z) [log10 g/cm^3];
! structure_equations, rezoning, starting_model, model_files. In .mod;
! rewritten from the EOS on every Newton iteration, at the structure that
! iteration linearized, so it can lag the last correction.
            double precision :: logRho(json)
! true where the zone is convective by the temperature-gradient test
! [logical]; structure_equations, rezoning, starting_model, model_files.
! In .mod; rewritten on every Newton iteration and by
! compute_shell_physics.
            logical :: convective_flag(json)
! mass fractions; row = species slot i_* of composition_layout, row
! i_metals the total Z [-]; burn_step, composition_step, rezoning,
! starting_model, model_files. In .mod; rows 12-15 only with
! ctrl%use_extended_composition.
            double precision :: xa(n_species_extended,json)
! surface X for the envelope integrations and the surface opacity tables;
! follows xa(i_h1,nz) [-]; rezoning, starting_model, model_rescaling.
! Refreshed when it moves by more than 1e-8 (refresh_surface_opacity), or
! set by a rescale.
            double precision :: surface_hydrogen_fraction
! surface Z for the envelope integrations, updated with
! surface_hydrogen_fraction [-]; rezoning, starting_model,
! model_rescaling
            double precision :: surface_metal_fraction
! log10(m(nz)/M), the mass fraction inside the fitting point [log10];
! henyey_iteration, rezoning, starting_model, mass_change,
! envelope_rebuild. Recomputed at every solve level (refresh_surface_bc).
            double precision :: log_fit_point_mass_fraction
! (log P, log T, log r, log m/M) where a vertex's envelope integration
! first reached envelope_store_log_temperature [log10 cgs];
! henyey_iteration. The rezone uses it to add a point at the surface.
            double precision :: stored_envelope_state(4)
! coefficients c(1:9) of the envelope fit at the fitting point: log R,
! log L and log Teff each = c log P + c log T + c [-]; henyey_iteration,
! model_files. Rebuilt by compute_surface_boundary when the triangle
! changes; read from YREC7/MODEL2 files, not written to .mod.
            double precision :: envelope_fit_coeffs(9)
! log10 L of the three vertices of the envelope-fit triangle
! [log10 Lsun]; henyey_iteration, model_files
            double precision :: triangle_log_l_lsun(3)
! log10 Teff of the three vertices of the envelope-fit triangle
! [log10 K]; henyey_iteration, model_files
            double precision :: triangle_log_teff(3)
! log10 P at the fitting point of each vertex's envelope integration
! [log10 dyn/cm^2]; henyey_iteration, model_files
            double precision :: fit_point_pressure(3)
! log10 T at the fitting point of each vertex [log10 K];
! henyey_iteration, model_files
            double precision :: fit_point_temperature(3)
! log10 r at the fitting point of each vertex [log10 cm];
! henyey_iteration, model_files
            double precision :: fit_point_radius(3)
! which triangle vertex stored stored_envelope_state (0: none) [-];
! step_driver, run_driver
            integer :: stored_vertex_index
! orientation (+1 or -1) of the envelope-fit triangle [-]; step_driver,
! run_driver
            double precision :: triangle_orientation

! ==== Rotation state ===================================================
! The rotation curve and the shape of the rotating star. The specific
! angular momentum j_rot is conserved; omega, the moments of inertia and
! the shape factors follow from it and the structure.
! Valid: between steps. Meaningful only with ctrl%rotation_active;
! without rotation the shape factors are 1.

! angular velocity of each shell [rad/s]; henyey_iteration, rezoning,
! starting_model, evolve_rotation, model_files. In .mod.
            double precision :: omega(json)
! specific angular momentum of each shell, the conserved quantity; omega
! follows from it (omega_from_j) [cm^2/s]; rezoning, starting_model,
! evolve_rotation, secular_transport
            double precision :: j_rot(json)
! moment of inertia of each shell [g cm^2]; henyey_iteration, rezoning,
! starting_model, evolve_rotation
            double precision :: i_rot(json)
! rotational kinetic energy of each shell, omega j dm/2 [erg];
! henyey_iteration, rezoning, starting_model
            double precision :: kinetic_energy_rot(json)
! rotational correction f_P to the hydrostatic equation, 1 without
! rotation [-]; henyey_iteration, rezoning, starting_model. Set from omega
! by rotation_shape_factors.
            double precision :: fp_rot(json)
! rotational correction f_T to the radiative-transfer equation, 1 without
! rotation [-]; henyey_iteration, rezoning, starting_model
            double precision :: ft_rot(json)
! r0, the mean radius of the shell's level surface [cm];
! henyey_iteration, rezoning, starting_model, evolve_rotation
            double precision :: mean_radius(json)
! Radau's eta_2 of each shell, the shape of its level surface [-];
! henyey_iteration, rezoning, starting_model, evolve_rotation
            double precision :: eta_squared(json)
! effective gravity averaged over the level surface, <g> [cm/s^2];
! henyey_iteration, rezoning, starting_model
            double precision :: mean_gravity(json)
! true where a zone is convective for angular-momentum transport,
! overshoot included [logical]; rezoning, starting_model. From
! am_convective_regions at the model read and every rezone, i.e. before
! every structure solve; henyey_iterate's omega_from_j holds those zones
! in solid-body rotation during the Newton iterations.
            logical :: am_transport_convective_flag(json)
! total angular momentum J of the star [g cm^2/s]; step_driver,
! run_driver, evolve_rotation. Zeroed at each card start, then totalled by
! the model read and every rezone.
            double precision :: total_angular_momentum
! total rotational kinetic energy [erg]; step_driver, run_driver,
! evolve_rotation
            double precision :: total_rotational_ke

! ==== Regions of the current model =====================================
! Point indices and bounds of the convective core, the surface convection
! zone, the hydrogen-burning shell and the mixed regions.
! Valid: for the model they were located on; located again at each stage
! that needs them (the composition step, the rezone, after convergence).

! outermost point of the convective core (1: none) [-]; step_driver,
! rezoning, run_driver. Set by the composition step's predictor pass, the
! rezone, burn_light_elements and, after convergence,
! locate_shell_boundaries.
            integer :: core_cz_top_index
! innermost point of the surface convection zone (nz: none) [-];
! step_driver, rezoning, run_driver. Set as core_cz_top_index is.
            integer :: envelope_cz_bottom_index
! first point of the overshoot region below the surface convection zone
! (-1: none) [-]; composition_step, rezoning, starting_model,
! evolve_rotation. From find_convection_zones; with ctrl%ladov the
! radiative points envelope_overshoot_bottom_index to
! envelope_overshoot_top_index take the adiabatic gradient.
            integer :: envelope_overshoot_bottom_index
! last point of that overshoot region (-1: none) [-]; composition_step,
! rezoning, starting_model, evolve_rotation
            integer :: envelope_overshoot_top_index
! first and last point of each mixed (convective plus overshoot) region,
! up to max_convective_zones of them [-]; composition_step, step_driver.
! Rewritten by every composition-step pass and by burn_light_elements.
            integer :: mixed_zone_bounds(max_convective_zones,2)
! the model burns hydrogen in a shell rather than in the core [logical];
! step_driver, run_driver
            logical :: has_h_shell
! first point of the H-burning shell, where X > ctrl%atime(1) [-];
! step_driver, run_driver. Valid when has_h_shell.
            integer :: h_shell_bot_index
! first point where X exceeds half its surface value [-]; step_driver,
! run_driver. Valid when has_h_shell.
            integer :: h_shell_mid_index
! last point of the H-burning shell [-]; step_driver, run_driver. Valid
! when has_h_shell.
            integer :: h_shell_top_index

! ==== Step control =====================================================
! The timestep and the flags that steer one step's solve.
! Valid: during and between steps.

! timestep of the step being taken (the next step's after
! compute_timestep) [s]; step_driver, run_driver, mass_change. Negative
! means solve without aging (a rescale or relaxation); made positive
! after convergence.
            double precision :: dt
! copy of dt: when the end-age stop clips the next timestep to zero, the
! timestep just taken is restored from it [s]; step_driver, run_driver
            double precision :: dt_saved
! dt in years [yr]; step_driver, run_driver. In .mod.
            double precision :: timestep_yr
! largest fractional change of omega over the step [-]; step_driver,
! run_driver. Written by compute_timestep and evolve_angular_momentum.
            double precision :: max_domega_frac
! the structure solve of this step diverged; a set flag stops the run on
! the next pass [logical]; step_driver, run_driver
            logical :: model_diverged_flag
! a model has been computed, so the next helium-flash card saves its
! starting model [logical]; run_driver, step_driver
            logical :: punch_pending_flag
! build a new envelope-fit triangle at the next solve [logical];
! step_driver, run_driver, mass_change. Set by a model read, the mass
! change, the envelope rebuild or ctrl%lnew0; cleared by
! compute_surface_boundary.
            logical :: recompute_envelope_triangle
! the rezone moved the outer point: redo the triangle at level 2
! [logical]; step_driver, run_driver
            logical :: reset_triangle

! ==== Step snapshot ====================================================
! Copies of model quantities at the start of the step (<member>_start),
! of the previous model (<member>_prev) or of the structure-rotation
! iteration (<member>_iteration_start), on the step's mesh. The rezone
! interpolates the old model onto the new mesh into the *_start arrays
! and then copies them into the live arrays, so after the rezone they
! hold the start-of-step model on the new mesh; without a rezone
! (ctrl%helium_flash_active) rezone_or_snapshot copies the structure
! into the *_start arrays itself. Valid: from the start-of-step copy to
! the end of the step.

! log P at the start of the step, on the current mesh [log10 dyn/cm^2];
! rezoning, starting_model, step_driver
            double precision :: logP_start(json)
! log T at the start of the step [log10 K]; rezoning, starting_model,
! step_driver
            double precision :: logT_start(json)
! log r at the start of the step [log10 cm]; rezoning, starting_model,
! step_driver
            double precision :: logR_start(json)
! L at the start of the step [Lsun]; rezoning, starting_model,
! step_driver
            double precision :: luminosity_lsun_start(json)
! log rho at the start of the step [log10 g/cm^3]; rezoning,
! step_driver. Not set by snapshot_step_start at a card start.
            double precision :: logRho_start(json)
! composition at the start of the step [-]; rezoning, step_driver,
! evolve_rotation. Copied before the composition step; the corrector
! passes and the structure-rotation iteration restart from it.
            double precision :: xa_start(n_species_extended,json)
! convective_flag at the start of the step, rotating models only
! [logical]; rezoning
            logical :: convective_flag_start(json)
! log Teff at the start of the step, for the HR-motion timestep limit
! [log10 K]; rezoning, starting_model, step_driver
            double precision :: log_Teff_start
! L(nz) at the start of the step, before the rezone, for the HR-motion
! timestep limit [Lsun]; rezoning, starting_model, step_driver
            double precision :: surface_luminosity_lsun_start
! kinetic_energy_rot at the start of the step, for d(E_rot)/dt in the
! energy equation of rotating models [erg]; step_driver
            double precision :: kinetic_energy_rot_start(json)
! omega at the start of the step [rad/s]; rezoning, starting_model. Read
! by the omega timestep limit and the transport.
            double precision :: omega_start(json)
! j_rot at the start of the step, the rezone's interpolation buffer
! [cm^2/s]; rezoning
            double precision :: j_rot_start(json)
! mean_gravity at the start of the step [cm/s^2]; rezoning
            double precision :: mean_gravity_start(json)
! mean_radius at the start of the step, the rezone's interpolation buffer
! [cm]; rezoning
            double precision :: mean_radius_start(json)
! eta_squared at the start of the step, the rezone's interpolation buffer
! [-]; rezoning
            double precision :: eta_squared_start(json)
! gradT at the start of the step, from which the rotation sub-steps
! interpolate in time (mid_timestep_model) [-]; rezoning. The rezone
! stores the twelve physics arrays from here on: the two energy
! sums in regrid_start_of_step_eps, the rest in rerun_physics_on_new_grid.
            double precision :: gradT_start(json)
! grada at the start of the step [-]; rezoning
            double precision :: grada_start(json)
! mu at the start of the step [-]; rezoning
            double precision :: mu_start(json)
! opacity_zone at the start of the step [cm^2/g]; rezoning
            double precision :: opacity_zone_start(json)
! cp at the start of the step [erg/(g K)]; rezoning
            double precision :: cp_start(json)
! dlnrho_dlnt (d ln rho/d ln T at constant P) at the start of the step [-];
! rezoning
            double precision :: dlnrho_dlnt_start(json)
! conv_vel at the start of the step [cm/s]; rezoning
            double precision :: conv_vel_start(json)
! kinematic_viscosity at the start of the step [cm^2/s]; rezoning
            double precision :: kinematic_viscosity_start(json)
! thermal_diffusivity at the start of the step [cm^2/s]; rezoning
            double precision :: thermal_diffusivity_start(json)
! eps_total at the start of the step, splined onto the step's mesh
! (regrid_start_of_step_eps) [erg/(g s)]; rezoning
            double precision :: eps_total_start(json)
! gradr at the start of the step [-]; rezoning
            double precision :: gradr_start(json)
! eps_total + eps_nu + eps_grav at the start of the step [erg/(g s)];
! rezoning
            double precision :: eps_plus_neu_grav_start(json)
! Li6 proton-capture rate rho X lambda at the start of the step; zero
! outward of log T = 6 [1/s]; burn_step, evolve_rotation. With
! ctrl%use_extended_composition.
            double precision :: rate_li6_start(json)
! Li7 proton-capture rate at the start of the step [1/s]; burn_step,
! evolve_rotation
            double precision :: rate_li7_start(json)
! Be9 proton-capture rate at the start of the step [1/s]; burn_step,
! evolve_rotation
            double precision :: rate_be9_start(json)
! deuterium burning rate at the start of the step, stored by the
! composition step's predictor pass [1/s]; burn_step, rezoning
            double precision :: deuterium_burning_rate_start(json)
! envelope-overshoot depth alpha_ov H_P below the surface convection zone
! at the start of the step [cm]; step_driver, evolve_rotation
            double precision :: envelope_overshoot_depth_start
! convective_turnover_timescale of the previous model [s]; observables,
! starting_model. Copied after each converged model (and at a model
! read); the wind interpolates between it and the current value at
! substep_time_fraction.
            double precision :: convective_turnover_timescale_prev
! photosphere_logP of the previous model, copied with
! convective_turnover_timescale_prev [log10 dyn/cm^2]; observables,
! starting_model
            double precision :: photosphere_logP_prev
! j_rot at the start of the structure-rotation iteration [cm^2/s];
! step_driver. Restored on each repeat pass when
! job%num_rotation_structure_iters > 1.
            double precision :: j_rot_iteration_start(json)
! xa_start at the start of the structure-rotation iteration [-];
! step_driver
            double precision :: xa_iteration_start(n_species_extended,json)
! base point of the surface convection zone as the deuterium-burning cap
! uses it (nz: none) [-]; step_driver, run_driver, mass_change. A copy of
! envelope_cz_bottom_index taken at the card start and each rezone.
            integer :: envelope_cz_bottom_index_start

! ==== Step working state ===============================================
! Values one stage computes for a later stage of the same or the next
! step: rates of change, burning rates, the depth of the surface mixed
! region, the turnover time, the shell physics that rotation uses, the
! outer layers of the latest envelope integration, the last rotation
! sub-step's omega. A few are printed too (the history prints the
! turnover time and the photospheric pressure).
! Valid: from the stage that writes them; each member says which.

! d log T/dt over the last step, for the next step's first guess [1/s];
! step_driver, rezoning, run_driver
            double precision :: dlog_temperature_dt(json)
! d log P/dt over the last step [1/s]; step_driver, rezoning, run_driver
            double precision :: dlog_pressure_dt(json)
! d ln L/dt over the last step, 2 (L - L_start)/(L + L_start)/dt [1/s];
! step_driver, rezoning, run_driver
            double precision :: dln_luminosity_dt(json)
! d log r/dt over the last step [1/s]; step_driver, rezoning, run_driver
            double precision :: dlog_radius_dt(json)
! age of the last converged model, for the disk-locking tests [Gyr];
! step_driver
            double precision :: disk_gate_age_gyr
! convective turnover time one pressure scale height above the base of
! the surface convection zone [s]; turnover_timescale, run_driver,
! mass_change. From the stitched model, every converged model; zero at a
! card start and without a surface convection zone.
            double precision :: convective_turnover_timescale
! log10 photospheric pressure of the last atmosphere integration that
! reached the photosphere [log10 dyn/cm^2]; envelope_integration. Set
! inside atm_get, i.e. by any envelope integration (the triangle's too).
            double precision :: photosphere_logP
! fraction of the step elapsed at the current rotation sub-step [-];
! evolve_rotation, starting_model. 0.5 after a model read.
            double precision :: substep_time_fraction
! omega of the latest rotation sub-step's model, the structure
! interpolated to the sub-step's end [rad/s]; mid_timestep_model,
! secular_transport, evolve_rotation. evolve_rotation copies it into
! omega at the end of the step; it outlives the step: the next step's
! first solid-body sub-step starts its Newton iteration from its mean.
            double precision :: omega_substep(json)
! Li6 proton-capture rate rho X lambda at the end of the step, on the
! converged model [1/s]; burn_step. With ctrl%use_extended_composition.
            double precision :: rate_li6(json)
! Li7 proton-capture rate at the end of the step [1/s]; burn_step
            double precision :: rate_li7(json)
! Be9 proton-capture rate at the end of the step [1/s]; burn_step
            double precision :: rate_be9(json)
! deuterium burning rate of each shell above the nuclear cutoff [1/s];
! burn_step. Zero while the surface D is below 1e-11.
            double precision :: deuterium_burning_rate(json)
! envelope-overshoot depth alpha_ov H_P at the end of the step [cm];
! step_driver, evolve_rotation
            double precision :: envelope_overshoot_depth_end
! radius of the base of the surface mixed region found by the last call
! within the step [cm]; convection_geometry, step_driver. With
! ctrl%use_extended_composition; zeroed at each step start to mean none
! yet.
            double precision :: cz_base_radius_prev
! mass accreted during the step [g]; mass_change. Set only when the star
! accretes; its readers check that it does.
            double precision :: accreted_mass
! shift of log P of the convection zone for the entropy of accreted
! matter [log10]; mass_change. Zero at present: change_mass zeroes
! envelope_specific_entropy first.
            double precision :: delta_log_pressure
! shift of log T for the same [log10]; mass_change
            double precision :: delta_log_temperature
! specific heat at constant pressure, for rotation (the viscosity, the
! transport) [erg/(g K)]; shell_physics
            double precision :: cp(json)
! mean molecular weight, for rotation [-]; shell_physics
            double precision :: mu(json)
! d ln rho/d ln T at constant P, for rotation [-]; shell_physics
            double precision :: dlnrho_dlnt(json)
! radiative thermal diffusivity K, from viscos, for rotational mixing
! [cm^2/s]; shell_physics
            double precision :: thermal_diffusivity(json)
! kinematic viscosity (radiative plus ionic), from viscos, for rotational
! mixing [cm^2/s]; shell_physics
            double precision :: kinematic_viscosity(json)
! the envelope of the latest envelope integration, point by point: log P,
! log T, log m, log rho, log r, X, Z, the gradients, the opacity and its
! derivatives [cgs, logs base 10]; envelope_integration, envelope_rebuild,
! starting_model. Written by the integration kernel through atm_get.
            type(envelope_structure_state) :: env_points
! the atmosphere of the latest envelope integration, point by point
! [cgs, logs base 10]; envelope_integration
            type(atmosphere_structure_state) :: atm_points
! the atmosphere integration's current point: log tau, log P, log T,
! log rho, the opacity, the ionization fractions [logs base 10];
! envelope_integration, mass_change. The mass change falls back on its
! last values when the Allard lookup fails.
            type(atmosphere_point_state) :: atm_point

! ==== Newton ===========================================================
! Arrays of the Henyey solve, rewritten on every iteration.
! Valid: during the solve; log_pressure_delta and log_temperature_delta,
! the changes over the step, until the next step's solve starts.

! block-elimination coefficients c(k,j,i) of henyey_elimination (its
! header gives the relations) [-]; henyey_elimination,
! structure_equations
            double precision :: elim_coeff(4,2,json)
! right-hand sides e(k,i): the residuals, then (after henyey_solve) the
! corrections dlog P, dlog T, dlog r, dL [log10, Lsun];
! henyey_elimination, henyey_iteration, structure_equations
            double precision :: elim_rhs(4,json)
! change of log P since the start of the step [log10]; step_driver,
! henyey_iteration. Zeroed (or set to the first-guess extrapolation)
! before the solve, then accumulated by every correction; feeds eps_grav.
            double precision :: log_pressure_delta(json)
! change of log T since the start of the step [log10]; step_driver,
! henyey_iteration. As log_pressure_delta.
            double precision :: log_temperature_delta(json)
! per equation (P, T, r, L): the largest |residual|, then the largest
! |correction|, then the signed correction there [log10, Lsun];
! henyey_elimination, henyey_iteration. Misnamed: holds residuals or
! corrections depending on the stage of the iteration.
            double precision :: max_residual(4)
! Newton iterations of the step's solve, counted over its levels up to
! the last convergence, for the run-log model line [-]; henyey_iteration
            integer :: newton_iterations = 0

! ==== Diagnostics written during the solve =============================
! Per-zone physics stored for output by compute_henyey_coefficients at
! every shell on every Newton iteration (energy channels, EOS and opacity
! values, gradients, pulsation derivatives). Rule 5 of the six
! (diagnostics written once, after convergence) is the planned
! replacement.
! Valid: after convergence they hold the last iteration's values, the
! structure that iteration linearized, before its correction when one was
! applied. Some are also rewritten by compute_shell_physics; the last
! writer wins.

! eps_grav = -T ds/dt of each zone [erg/(g s)]; structure_equations,
! step_driver
            double precision :: eps_grav(json)
! the extra energy sources' sum at each zone
! (core/structure/energy_sources.f90) [erg/(g s)]; structure_equations.
! Zero without a registered source.
            double precision :: eps_extra(json)
! nuclear energy generation, the sum of the chains' eps: pp1 + pp2 +
! pp3 + CNO + triple alpha + the alpha captures [erg/(g s)];
! structure_equations. Zero below the burning cutoff; the neutrino losses
! are not in it.
            double precision :: eps_total(json)
! per channel i_eps_*: rows 1-5 the fractions of eps_total of pp1, pp2,
! pp3, CNO and triple alpha, row 6 eps_nu, row 7 eps_grav [rows 1-5: -;
! rows 6-7: erg/(g s)]; structure_equations. The alpha captures have no
! row, so rows 1-5 sum to less than 1 where they burn.
            double precision :: eps_channels(7,json)
! luminosity per channel i_lum_* (pp1, pp2, pp3, CNO, 3 alpha, neutrinos,
! gravothermal, alpha captures, the extra energy sources), scaled to sum
! to L(nz) [Lsun]; structure_equations, observables, henyey_iteration,
! model_files. Summed every iteration, rescaled again after convergence
! (renormalize_luminosity_breakdown); compute_timestep reads it; read from
! YREC7/MODEL2 files, zeroed by the .mod reader.
            double precision :: luminosity_breakdown_lsun(n_lum_channels)
! gas-pressure fraction P_gas/P [-]; structure_equations
            double precision :: beta(json)
! electron degeneracy parameter [-]; structure_equations. 0 unless the
! EOS treats the gas as fully ionized.
            double precision :: eta(json)
! Rosseland mean opacity kappa [cm^2/g]; structure_equations,
! shell_physics. Also written by compute_shell_physics (the rezone,
! rotation, the card start); the last writer wins.
            double precision :: opacity_zone(json)
! radiative temperature gradient grad_rad [-]; structure_equations,
! shell_physics. As opacity_zone.
            double precision :: gradr(json)
! actual temperature gradient d ln T/d ln P [-]; structure_equations,
! shell_physics. As opacity_zone.
            double precision :: gradT(json)
! adiabatic gradient grad_ad [-]; structure_equations, shell_physics. As
! opacity_zone.
            double precision :: grada(json)
! convective velocity from mixing-length theory [cm/s];
! structure_equations, shell_physics. As opacity_zone.
            double precision :: conv_vel(json)
! first adiabatic exponent Gamma_1 [-]; structure_equations
            double precision :: adiabatic_index_gamma1(json)
! solar-neutrino fluxes at the Earth by slot i_nu_* (slots 9-10: the
! He3+He3 and He3+He4 rates in the same units) [10^10/(cm^2 s)];
! structure_equations, observables, neutrino_flux_output. Only with
! ctrl%calc_neutrinos; also overwritten by the card-start neutrino table.
            double precision :: neutrino_flux_total(n_nu_fluxes)
! d eps/d ln rho, absolute (not logarithmic) [erg/(g s)];
! structure_equations. Zero below the burning cutoff.
            double precision :: pulse_deps_dlnrho(json)
! d eps/d ln T, absolute [erg/(g s)]; structure_equations. Zero below
! the burning cutoff.
            double precision :: pulse_deps_dlnt(json)
! d ln kappa/d ln rho [-]; structure_equations
            double precision :: pulse_dlnkap_dlnrho(json)
! d ln kappa/d ln T [-]; structure_equations
            double precision :: pulse_dlnkap_dlnt(json)
! specific heat at constant pressure c_P [erg/(g K)]; structure_equations
            double precision :: pulse_specific_heat(json)
! mean molecular weight mu = R/(specific gas constant) [-];
! structure_equations
            double precision :: pulse_mean_molecular_weight(json)
! d ln rho/d ln T at constant P [-]; structure_equations
            double precision :: pulse_dlnrho_dlnt(json)
! 1/mu_e [-]; structure_equations
            double precision :: pulse_electron_mean_weight_inverse(json)

! ==== Diagnostics written after convergence ============================
! Observables computed once per model by compute_observables (after the
! stitched model is built), and the rotational transport's coefficients
! and velocities kept for the profiles.
! Valid: from compute_observables (or the rotation step) until the next
! model; zero at a run start.

! the stitched model: the converged star (interior, envelope and
! atmosphere) on one extended grid, with its profile and pulse columns
! [cgs]; stitched_model. Allocated, and zeroed, by the first build.
            type(stitched_model_record) :: stitched
! each rotational transport mechanism's velocity at interface k, slots
! i_mech_* of transport_mechanism_layout [cm/s]; secular_transport. From
! the last iteration of the step's last rotation sub-step.
            double precision :: mechanism_velocity(n_transport_mechanisms, json)
! angular-momentum diffusion coefficient D_omega on the model grid
! [cm^2/s]; secular_transport. The step's last sub-step wins; 0 without
! rotation.
            double precision :: am_diffusion_coeff(json)
! rotational mixing coefficient D_mix on the model grid [cm^2/s];
! secular_transport. As am_diffusion_coeff.
            double precision :: mixing_diffusion_coeff(json)
! log10 of the photospheric radius from L = 4 pi R^2 sigma Teff^4
! [log10 Rsun]; observables
            double precision :: log_R_surface_rsun
! log10 surface gravity G M/R^2 [log10 cm/s^2]; observables
            double precision :: log_g_surface
! frequency of maximum oscillation power, scaled from the Sun's by g and
! Teff [uHz]; observables
            double precision :: nu_max_uhz
! large frequency separation scaled from the mean density [uHz];
! observables
            double precision :: delta_nu_rho_uhz
! asymptotic p-mode large separation, from the sound-travel time [uHz];
! observables
            double precision :: delta_nu_uhz
! asymptotic l = 1 g-mode period spacing of the inner cavity, from the
! Brunt-Vaisala integral [s]; observables. 0 without a cavity.
            double precision :: delta_Pg
! moment of inertia of the star [g cm^2]; observables
            double precision :: total_moment_of_inertia
! moment of inertia of the surface convection zone [g cm^2]; observables
            double precision :: cz_moment_of_inertia
! surface rotation period 2 pi/omega(nz), capped at 9999 [d];
! observables. 0 without rotation.
            double precision :: rotation_period_days
! surface equatorial velocity omega(nz) R [km/s]; observables
            double precision :: surf_velocity_kms
! mass coordinate of the bottom of the H-burning shell [Msun];
! observables. 0 without a shell, as the five below.
            double precision :: h_shell_bot_mass_msun
! mass coordinate of the middle of the H-burning shell [Msun];
! observables
            double precision :: h_shell_mid_mass_msun
! mass coordinate of the top of the H-burning shell [Msun]; observables
            double precision :: h_shell_top_mass_msun
! radius of the bottom of the H-burning shell [r/R]; observables
            double precision :: h_shell_bot_radius
! radius of the middle of the H-burning shell [r/R]; observables
            double precision :: h_shell_mid_radius
! radius of the top of the H-burning shell [r/R]; observables
            double precision :: h_shell_top_radius
! central temperature, extrapolated from the first point [log10 K];
! observables
            double precision :: central_log10_temperature
! central pressure [log10 dyn/cm^2]; observables
            double precision :: central_log10_pressure
! central density, from the EOS [log10 g/cm^3]; observables
            double precision :: central_log10_density
! central P_gas/P [-]; observables
            double precision :: central_beta
! central electron degeneracy parameter [-]; observables
            double precision :: central_degeneracy_eta
! mass of the surface convection zone [Msun]; observables
            double precision :: envelope_mass_msun
! radius of the base of the surface convection zone [Rsun]; observables,
! starting_model. Written twice per model: by locate_surface_cz_base,
! then by the turnover-time walk.
            double precision :: envelope_radius_rsun
! mass of the convective core [Msun]; observables. 0 without one.
            double precision :: core_cz_mass_msun
! temperature at the base of the surface convection zone, linear [K];
! observables
            double precision :: envelope_cz_temperature
! density at the base of the surface convection zone [g/cm^3];
! observables
            double precision :: envelope_cz_density
! pressure at the base of the surface convection zone [dyn/cm^2];
! observables
            double precision :: envelope_cz_pressure
! opacity at the base of the surface convection zone [cm^2/g];
! observables
            double precision :: envelope_cz_opacity
! chlorine (37Cl) detector capture rate [SNU]; observables. 0 without
! ctrl%calc_neutrinos.
            double precision :: cl37_snu_rate
! gallium (71Ga) detector capture rate [SNU]; observables. 0 without
! ctrl%calc_neutrinos.
            double precision :: ga71_snu_rate

! ==== Private to one module ===========================================
! Values kept on the star for one module alone, which writes and reads
! them: the state it carries between its calls (the calibrations'
! previous checks, models and corrections, the central EOS call's Saha
! state, the Newton solve's boundary rows and convergence log, the
! surface-mixing search's starting zone). tools/gen_star_docs.py places
! every member that one module alone uses here, and nothing else.
! Valid: only inside that module; each member says when.

! log(L/Lsun) at the previous solar-calibration check, only printed
! [log10 Lsun]; calibration
            double precision :: log_l_lsun_prev_check
! log(R/Rsun) at the previous solar-calibration check, only printed
! [log10 Rsun]; calibration
            double precision :: log_r_rsun_prev_check
! the solar calibration's Newton correction to X [-]; calibration
            double precision :: delta_x
! the solar calibration's Newton correction to the mixing length [-];
! calibration
            double precision :: delta_alpha
! age interpolated at the star calibration's target radius [Gyr];
! calibration
            double precision :: age_gyr_at_target_radius
! age of the previous model, for the star calibration [Gyr]; calibration
            double precision :: age_gyr_prev_model
! log(L/Lsun) interpolated at the star calibration's target radius
! [log10 Lsun]; calibration
            double precision :: log_l_lsun_at_target_radius
! log(L/Lsun) at the target radius in the previous star-calibration
! run [log10 Lsun]; calibration
            double precision :: log_l_lsun_at_target_radius_prev_run
! log(L/Lsun) of the previous model, for the star calibration
! [log10 Lsun]; calibration
            double precision :: log_l_lsun_prev_model
! R/Rsun of the previous model, for the star calibration [Rsun];
! calibration
            double precision :: r_rsun_prev_model
! zone of the base of the surface mixed region found by the last call
! within the step (cz_base_radius_prev), the starting zone of the
! next call [-]; convection_geometry
            integer :: envelope_cz_base_zone_prev
! the Saha state the central EOS call carries from one model to the next
! [-]; observables
            integer :: central_saha_state = 0
! point of the largest correction of each variable, for the convergence
! log [-]; henyey_iteration
            integer :: max_correction_index(4)
! outer boundary-condition rows from the envelope fit, in henyey_solve's
! form [-]; henyey_iteration
            double precision :: surface_bc(6)
      end type star_info

! Module-level (deliberately outside the star, which run_yrec builds
! afresh for every run): the namelist-path overrides of the C API
! (yrec_capi) and the flag yrec_reset sets so that evolve_step
! re-initializes its saved locals.
      character(len=256), public, save :: control_nml_override = ' '
      character(len=256), public, save :: physics_nml_override = ' '
      logical, public, save :: evolve_step_reset_pending = .false.

end module star_info_lib

!======================================================================
! Sources and history
!
! star_info was added in 2026 (docs/ROADMAP.md, "Phase four: the star
! layer"), adapted from MESA's star_info. Its first members were the
! model arrays once declared in the F77 `program main` and threaded
! through 20-60-argument positional signatures (CRRECT took 60, STARIN
! 50, MIDMOD and WRTOUT 43, HPOINT 35, GETW 29, MIX 21). The COMMON
! blocks and the one-instance state modules (prev_model, shell_diag,
! run_diag, rot_diff, turnover, light_burn, flux, env_comp, thermo, evo,
! pulse) were then folded in as flat members, and the computed values
! once kept in the controls buffer moved here. The module-level star was
! removed in 2026-10. The members' F77 homes:
!   /calsun/    dlum_dx, drad_dx, dlum_dalpha, drad_dalpha, log_l_lsun_prev_check,
!               log_r_rsun_prev_check, delta_x, delta_alpha,
!               solar_calibration_converged
!   /cent/      central_log10_temperature, central_log10_pressure,
!               central_log10_density, envelope_mass_msun, envelope_radius_rsun
!   /comp2/     envelope_helium_fraction, envelope_he3_fraction
!   /const/     the solar units (solar_luminosity_cgs ...; CMBOL, the
!               Sun's bolometric magnitude, is physical_constants'
!               solar_bolometric_magnitude)
!   CMIXL       mixing_length_alpha (the one mutable value the constants
!               module carried)
!   /deuter/    deuterium_burning_rate, deuterium_burning_rate_start
!               (per Gyr until 2026-10, now per s),
!               envelope_cz_bottom_index_start (JCZ)
!   /entrop/    dlog_temperature_dt, dlog_pressure_dt, dln_luminosity_dt,
!               dlog_radius_dt
!   /liov/      envelope_overshoot_depth_start, envelope_overshoot_depth_end
!   /monte2/    job%s11_rate ... job%age_target
!   /newrat/, /oldrat/  rate_li6, rate_li7, rate_be9 and their _start
!               copies
!   /oldrot/    omega_start, j_rot_start, mean_gravity_start (HG),
!               mean_radius_start, eta_squared_start
!   /origstart/ j_rot_iteration_start, xa_iteration_start
!   /prevcz/    cz_base_radius_prev, envelope_cz_base_zone_prev
!   /sound/     adiabatic_index_gamma1
!   /theage/    age_gyr (DAGE)
!   /varfc/     zahn_coupling (VFC)
!   /zramp/     job%nk (NK)
! and the per-zone diagnostics of the former shell_diagnostics_state:
! eps_total (SESUM), eps_channels (SEG), beta (SBETA), eta (SETA),
! conv_vel (SVEL), opacity_zone (SO, the opacity, not O16), and gradr,
! gradT, grada (the rows of DEL_GRAD(3,:)). The legacy flags and
! indices: LCZ (am_transport_convective_flag), LPUNCH
! (punch_pending_flag), IOV1 and IOV2 (envelope_overshoot_bottom_index,
! envelope_overshoot_top_index). am_transport_convective_flag was never
! assigned before 2026-10: main.f's LCZ was a local passed only to
! CRRECT and GETROT, so the Newton iterations treated every zone as
! radiative.
!
! Renamed in 2026-10 (old name first): gravitational_luminosity ->
! eps_grav; temperature_, pressure_, luminosity_, radius_entropy_term ->
! dlog_temperature_dt, dlog_pressure_dt, dln_luminosity_dt,
! dlog_radius_dt; pressure_scale_height_start, _end ->
! envelope_overshoot_depth_start, _end; trial_log_temperature,
! trial_log_luminosity, trial_sign_flag -> triangle_log_teff,
! triangle_log_l, triangle_orientation; istore_flag ->
! stored_vertex_index; log_r_prev_model -> r_rsun_prev_model (it is
! linear); solar_calibration_active -> solar_calibration_converged; dage
! -> age_gyr; accreted_mass_fraction -> accreted_mass (it is a mass);
! omega_mid (rot_scr%omega_mid before) -> omega_substep (2026-10-28);
! pulse_dlneps_dlnrho, _dlnt -> pulse_deps_dlnrho, _dlnt (absolute
! derivatives); mean_molecular_weight -> mu; the module variables
! env_struct, atmo_struct and observables' ksaha_center -> env_points,
! atm_points and central_saha_state; es_circulation_velocity,
! gsf_circulation_velocity, secular_shear_velocity -> the slots of
! mechanism_velocity. state/README.md section 5 lists the members
! removed in 2026-10.
!
! The rotation-scaled solar wind of core/step/mass_change.f90, commented
! out on 2026-10-12, read four members that were removed with it: its
! rate solar_wind_mass_loss_rate_msun_yr [Msun/yr], its reference and
! maximum omega wind_reference_omega and wind_max_omega [rad/s] (double
! precision), and its switch use_rotation_scaled_solar_wind (logical).
! Reviving it means declaring them again, in the constants set at setup,
! and filling them from the controls.
!
! Renamed after the 2026-10-29 review (old name first; the earlier values
! as <member>_start / _prev, the legacy names, the unit suffixes):
! old_omega -> omega_start,
! old_specific_angular_momentum -> j_rot_start,
! old_hg -> mean_gravity_start, old_mean_radius -> mean_radius_start,
! old_eta_squared -> eta_squared_start, old_gradT -> gradT_start,
! old_grada -> grada_start, old_gradr -> gradr_start,
! old_mu -> mu_start, old_opacity -> opacity_zone_start,
! old_cp -> cp_start, old_qdt -> dlnrho_dlnt_start,
! old_conv_vel -> conv_vel_start, old_visc -> kinematic_viscosity_start,
! old_thdif -> thermal_diffusivity_start,
! old_eps_total -> eps_total_start,
! old_eps_plus_neu_grav -> eps_plus_neu_grav_start,
! kinetic_energy_rot_old -> kinetic_energy_rot_start,
! convective_turnover_timescale_old -> convective_turnover_timescale_prev,
! pphot0 -> photosphere_logP_prev,
! orig_specific_angular_momentum -> j_rot_iteration_start,
! orig_composition -> xa_iteration_start,
! jcz -> envelope_cz_bottom_index_start,
! iov1 -> envelope_overshoot_bottom_index,
! iov2 -> envelope_overshoot_top_index, vfc -> zahn_coupling,
! qdt -> dlnrho_dlnt, thdif -> thermal_diffusivity,
! visc -> kinematic_viscosity, xnew -> surface_hydrogen_fraction,
! znew -> surface_metal_fraction, senv -> log_fit_point_mass_fraction,
! zenvm -> envelope_metals_no_cno,
! amuenv -> envelope_inverse_atomic_weight,
! fxenv -> envelope_number_fractions, pphot -> photosphere_logP,
! tenv -> envelope_store_log_temperature, atm_hras -> hsra_log_t_offset,
! fracstep -> substep_time_fraction,
! h_shell_zone_begin -> h_shell_bot_index,
! h_shell_midpoint_zone -> h_shell_mid_index,
! h_shell_end_index -> h_shell_top_index,
! envelope_mass -> envelope_mass_msun,
! core_cz_mass -> core_cz_mass_msun,
! h_shell_bot_mass -> h_shell_bot_mass_msun,
! h_shell_mid_mass -> h_shell_mid_mass_msun,
! h_shell_top_mass -> h_shell_top_mass_msun,
! envelope_radius -> envelope_radius_rsun,
! log_R_surface -> log_R_surface_rsun, nu_max -> nu_max_uhz,
! delta_nu -> delta_nu_uhz, delta_nu_rho -> delta_nu_rho_uhz,
! age_prev_model -> age_gyr_prev_model,
! age_at_target_radius -> age_gyr_at_target_radius,
! log_l_at_target_radius -> log_l_lsun_at_target_radius,
! log_l_at_target_radius_prev_run -> log_l_lsun_at_target_radius_prev_run,
! log_l_prev_model -> log_l_lsun_prev_model,
! log_l_prev -> log_l_lsun_prev_check,
! log_r_prev -> log_r_rsun_prev_check,
! luminosity_breakdown -> luminosity_breakdown_lsun,
! log_L -> log_L_lsun, triangle_log_l -> triangle_log_l_lsun
! and stotal, a copy of log_total_mass, merged into it.
!======================================================================

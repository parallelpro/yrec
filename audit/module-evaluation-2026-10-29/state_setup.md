# state/ + setup/ review (2026-10-29, read-only evaluation; paths under src/)

Scores vs rotation/: state Readability 2.5, Organization 2, Maintainability 2.5, Extensibility 2.5, Docs 2.5;
setup 4 / 3 / 3.5 / 3 / 3.

The layout modules, envelope_points and stitched_model types and the settings builders meet the rotation bar;
state/README.md's tables (meaning, unit, writers, readers, lifetime of every member) are the best newcomer
resource in the star layer. The weak point is star_info: order and comments from the F77 COMMON blocks (not
lifetime or owner); ctrl/job split mislabelled; declared control defaults not the real ones for 107 controls;
the README that rescues it is not generated and has lost 12 members. physical_constants and
integration_parameters hide configuration in module variables set at run time, read by star-free eos/,
convection/ and a star-blind kernel (Code style rule 3). The settings-record convention covers 6 domains
(kap, eos, atm, convection, element diffusion, rotation), not net (builder in core), the envelope kernel, or
the constants. setup/ also holds map_user_inputs, which is really part of the controls read.

## Possible bugs
1. The settling switch never turns back on for later cards -- s%job%diffuse_helium_active is a user control
   that two domains overwrite as a latch (core/step/composition_step.f90:263, rotation/mid_timestep_model.f90:466;
   both say settling stops "for good"); only derive_controls sets it true, once per run (from read_controls).
   VERIFIED by the integrator. A later card (or Monte Carlo run) that starts a hydrogen-rich model after an
   earlier card's centre fell below hydrogen_diffusion_floor runs without settling. (atm_choice is saved and
   restored per cycle, core/run/run_driver.f90:230/:263; job%use_mass_accretion has the same latch,
   core/step/mass_change.f90:125, documented as "for the rest of the run".) Fix: keep the control immutable
   and put the latch on the star, re-armed per card. Measure (no shipped deck is expected to move).
2. Declared control defaults are not the real ones: 107 of 274 controls declare a default that never takes
   effect; the real one comes from set_legacy_defaults (io/legacy_controls.f90:340-666) through
   adopt_canonical_names. Examples (declared vs effective): es/gsf/secular_shear velocity and mixing scales 0
   vs 1; s0_pp 0 vs 4.01e-22; atm_error_tol 0 vs 3e-4; settling_tolerance 0 vs 1e-8; scale_by_rossby_number
   .false. vs .true.; wind_law_name ' ' vs 'K97'. controls_state() is not the pristine defaults;
   test_rotation.f90:351 sets mu_gradient_scale by hand because of this; tools/gen_inlist.py's "defaults stay
   defined in one place" is false.
3. A rotation kernel test runs with an unset constant (test only): test_rotation_kernels never sets csigl;
   spindown_check gets log10_radius_from_l_teff ~2.12 dex low (radius ~1/130 of true); self-consistent, passes.
4. Controls stored but never read: zsi, alex_mixture_x, alex_mixture_z (io/controls_buffer.f90:195, :415; not in
   the registry; read nowhere; generated into s%ctrl). mu_gradient_scale has no effect (registry says so; O15),
   status 'canonical', a test still sets it.
5. Dead constants: cmkhn (state/physical_constants.f90:27) never assigned or read; ca3l assigned
   (core/run/run_setup.f90:94) never read.

## Findings (ranked)
1. state/controls_state_def.inc, job_controls_def.inc (generated from io/controls_buffer.f90 initializers):
   two default sources disagree for 107 controls. One default source: move set_legacy_defaults' literals into
   the controls_buffer initializers and seed the legacy names from the canonical ones; a test that
   controls_state() equals defaults/*.defaults. A for decks; M.
2. tools/gen_controls_state.py:37-67 JOB_MEMBERS: the ctrl/job split is mislabelled -- 33 of the 70 job members
   are never written after the read (rotation_active, use_extended_composition, core_overshoot_active,
   envelope_overshoot_active, lovstm, use_semiconvection, instability_transport_active, disk_*, atm/env_step_*,
   acfpft, pulse_format, the 8 *_limit arrays, ...); the "model-restore set (getyrec7/getmodel2 write them
   back)" does not exist. Make job exactly the 37 actually mutated, enforced by a check_boundaries rule (a
   member is in job iff code outside the read assigns it). A, M.
3. state/star_info_lib.f90:166-529: member order follows COMMON provenance, not README sec.2's 12 lifetime
   classes (Model class spread over 168-173, 187-191, 224-230, 406, 438; calibration over 296, 318, 440).
   Reorder by lifetime class with one line per member (meaning, unit, writer); history to a Sources block. A, M.
   (Later wave.)
4. state/README.md tables drifted; the script that made them (README:53 "derived at e9c00d7f") is not in the
   repo: the 12 old_* members (56aa2e6f) have no rows; :23 describes the retired star0 snapshot (now
   yrec_reset.f90's star_template); :20 says atmosphere and envelope "are not in star" (contradicts the rows for
   env_points, atm_points, stitched); module table (:5-12) omits envelope_points, stitched_model_state,
   transport_mechanism_layout; :159 lists disk_gate_age_gyr's reader as rotation_scratch_lib (now types only)
   and misses henyey_iteration, rezoning, envelope_rebuild, mid_timestep_model, transport_checks; :43 says 228
   members (237 besides ctrl, job); :21 cites age_gyr as an exception to the unit rule (it has the suffix).
   Generate the written-by/read-by columns from check_star_contracts' analysis via generated_files.py, with a
   --check test that every member has a row. A, M.
5. state/physical_constants.f90:25-30 with core/run/run_setup.f90:88-121: run-time "constants" gas_constant,
   csig, csigl, radiation_constant_over_3, cmkh, cdelrl, cmixl2, cmixl3 are module variables assigned from
   literals in core (CODATA 1986: sigma 5.67051e-5, R 8.314510e7, m_e 9.1093897e-28, k 1.380658e-16,
   h 6.6260755e-27, G 6.6726e-8 -- no file says so). Star-free domains read them: eos/ (gas_constant, cmkh,
   radiation_constant_over_3), convection/ (cdelrl, cmixl2, cmixl3, cgl), so eos/test/test_eos.f90:203 and
   convection/test/test_convection.f90:86 must call core's set_up_constants_and_tables(star); cause of bug 3.
   Parameters with sources, each compile-time value checked against the run-time one (rule 3); cgl (depends on
   ctrl%G_cgs) into the settings records; delete cmkhn and ca3l. A once bit-identity verified. M.
6. state/integration_parameters.f90: hidden configuration -- tolerance_fraction, max_stage_index, extrap_order
   in module variables, set in io/read_controls.f90:252-254, read only by
   core/envelope/envelope_integration_kernel.f90 (460, 525, 676, 830), which is in STAR_BLIND_FILES ("takes
   its configuration explicitly"). Header names bsstep.f90 and atm_lib.f90 as readers and "a DATA statement"
   as the default (both wrong; default legacy_controls.f90:408-410). Fields of the kernel's cfg record (filled
   in envelope_integration.f90); delete the module. A, S.
7. state/controls_sync.f90 is in the wrong directory: uses controls_buffer (io), used only by
   io/read_controls.f90, so state/ depends on io/; not private, re-exports the whole buffer. Generator emits
   it into io/, private with a public list. A, S.
8. Star fields that are constants or derived controls: dlum_dx, drad_dx, dlum_dalpha, drad_dalpha set to the
   same literals every call (core/run/calibration.f90:76-79); solar_bolometric_magnitude = 4.79
   (run_setup.f90:86); atm_hras = atm's hsra_t_tau_offset(); use_two_z_tables derived at run_driver.f90:116,
   read only by opacity_settings; tenv = 0.5*(tenv0+tenv1). Named parameters or computed in the builder.
   (Later wave.)
9. Duplicated quantities: stotal == log_total_mass (5 writers keep it equal: starting_model.f90:346,
   mass_change.f90:540, model_rescaling.f90:212/285/316); timestep_yr == dt/seconds_per_year; solar_mass_cgs ==
   ctrl%solar_mass_cgs; three pairs written by two procedures (cp/pulse_specific_heat,
   mu/pulse_mean_molecular_weight, qdt/pulse_dlnrho_dlnt; compute_shell_physics vs
   compute_henyey_coefficients); jcz a lagged copy of envelope_cz_bottom_index (run_driver.f90:361,
   step_driver.f90:333) -- rename. (Later wave.)
10. Names: six conventions for an earlier value (*_start, old_*, *_old, pphot0, *_prev, orig_*); ~24 cryptic
   legacy names (jcz, iov1/iov2, vfc, qdt, thdif, visc, xnew, znew, stotal, senv, zenvm, amuenv, fxenv, pphot,
   tenv, atm_hras, old_hg, qs0e_scale, qqs0ee_scale, ...). _start for start-of-step snapshots and _prev for the
   previous model. (Later wave.)
11. Unit-suffix rule broken by ~15 members: envelope_mass, core_cz_mass, h_shell_{bot,mid,top}_mass (Msun);
   envelope_radius (Rsun); log_R_surface (log Rsun); nu_max, delta_nu, delta_nu_rho (uHz); age_prev_model,
   age_at_target_radius (Gyr); luminosity_breakdown, log_L, triangle_log_l (Lsun). Rename; output columns keep
   MESA names. (Later wave.)
12. ~104 ctrl members only feed map_user_inputs or derive_controls (tol_*, time_*, flag_*, s0_*,
   pmm_exponent_*, unit numbers); their s%ctrl copies are never read. Two pairs alias: ctrl%tol_dm_max vs
   ctrl%chi_grid_scale(ichi_dm_max) agree only with lnewvars; ctrl%s0_pp vs cross_section_scale only with
   use_new_nuclear_rates. An INPUT_ONLY set in the generator, or document them; delete zsi and
   alex_mixture_x/z; mark mu_gradient_scale 'unused'.
13. setup/ organization: map_user_inputs (setup/user_inputs.f90) writes the buffer (why setup/ is in
   CONTROLS_BUFFER_IMPORTERS) -> io/ beside derive_controls (later wave); net's builder configure_network lives
   in core/step/burn_step.f90:52, takes intent(inout) although its contract writes nothing, builds and installs
   in one call; two naming schemes (opacity_settings / equation_of_state_settings / atmosphere_settings vs
   *_settings_from_star); user_inputs.f90 reference S-factors in DATA statements (45-92, should be
   parameters), .gt., bare use controls_buffer (34).
14. setup/README.md stale: :23 "element_diffusion/ still reads star%ctrl and star%job directly" (false; in
   STAR_FREE_DIRS since 926f34ef); :11 "or star in the owners" (no global star since c35778aa); :21 omits
   element_diffusion/test; neither README has an "add a field / control / domain" section.
15. state/star_info_lib.f90 stale comments: :13-18 names "main's own, unchanged"; :20-22 run control "stays in
   arguments" yet dt, dt_saved, newton_iterations, max_residual are members; :33-35 *_mid copies (gone since
   f0e8b4b6); :100-104 and :506-510 "one root, one file" and "envstruct/atmstruct NOT here", :515-517 ctrl and
   job "the only nested sub-structs" (env_points, atm_points, stitched are nested); :233-234 atm_hras
   "Krishna-Swamy" (it is the HSRA fit, run_setup.f90:124); :246-247 burn_lib's compute_zone_energy (it is
   configure_network); :290 compute_acoustic_depth / acoustic-depth family (no such members); :292-293
   use_two_z_tables "derived by setupopac" (run_driver sets it); :313-317 commented-out declarations; :325-334
   prev/diag/run/rot sub-structs and star coupling "grep-able" (gone; domains star-free); :111
   "state/star_job_lib.f90" (gone); :484 dangling "former common/i2o/"; 11 legacy routine names (parmin,
   setupopac, surfbc, envint, getfc, calcad, wrtout).

## Domains still reading s%ctrl / s%job directly
net, kap, eos, atm, numerics, convection, element_diffusion: 0. rotation: mid_timestep_model.f90 3 refs (445,
466, 556; job%diffuse_helium_active; 466 writes it). Star-layer core (allowed): mass_change 13,
structure_equations 7, henyey_iteration 22, surface_boundary 10, envelope_integration 15 (cfg built from
ctrl/job every call), envelope_rebuild 9, composition_step 13, convection_geometry 1, burn_step 27
(configure_network), timestep_control 37, rezoning 33, model_rescaling 55, starting_model 69,
shell_boundaries 2, observables 7, stitched_model 8, neutrino_flux_output 1. s%ctrl is never assigned after
the read in production.

## Structural recommendations
1. One source of truth for controls: a registry row carries name, type, default, unit, target domain;
   generate controls_buffer declarations, defaults, controls_state, the job/ctrl split (job derived from code
   that writes the member); retire set_legacy_defaults' literals.
2. Single-source, grouped star documentation: reorder star_info by lifetime class with a comment per member;
   generate README sec.3 from the code with a --check test.
3. Finish the no-hidden-configuration rule and the directory split: constants as parameters;
   integration_parameters and cgl into settings records; every builder in setup/ under one name
   <domain>_settings_from_star, installed from one place (install_run_settings(s), install_card_settings(s));
   controls_sync and map_user_inputs to io/, leaving state/ with types, layouts, constants and setup/ with
   builders only.

## Newcomer traces
(a) Add a star field: declare in state/star_info_lib.f90 (no rule for which block; non-zero start needs an
initializer since the run starts from a static template); a row in state/README.md sec.3 by hand (nothing
checks; how the 12 old_* rows went missing); assign through s and update the star reads/writes lines
(check_star_contracts.py --fix; only 26 of 137 procedures carry them); optionally io/model_files.f90,
stitched_model.f90 and column lists. Rebuilds everything.
(b) Add a control into a domain's settings record: io/controls_buffer.f90 (legacy twin's default in
io/legacy_controls.f90); registry row; namelist group in io/legacy_controls.f90 for two-file decks;
gen_inlist, gen_controls_state (JOB_MEMBERS first if run-written), gen_controls_dump; the domain's settings
type; the copy in setup/*_settings_setup.f90 + check_star_contracts --fix; read it; YREC_REGEN_CONTROLS=1
pytest test_controls.py. 5 hand edits, 3 generators, 1 regeneration.
(c) Add a physics domain: new directory with <dom>_lib (facade) and <dom>_configuration (settings type,
module-level record, setter); a builder in setup/ (two naming schemes); install by hand in
core/run/run_setup.f90 (per run) or begin_kind_card in core/run/run_driver.f90 ~:312 (per card) -- no install
list; tools/check_boundaries.py (DOMAINS, PUBLIC, STAR_FREE_DIRS, FACADE_ONLY_DOMAINS); Makefile and
gen_makefile_deps.py; controls by (b); README and tests (a unit test must call core's
set_up_constants_and_tables on a star if the domain uses gas_constant, cgl or cdelrl).

Fine as is: composition_layout, neutrino_flux_layout, transport_mechanism_layout; envelope_points and
stitched_model_state; the generation approach (generated_files.py --check, GENERATED headers); the setup
builders (pure copies, intent(in), accurate contracts, every settings field read); no star member
unreferenced or write-only (all 237 checked); s%ctrl immutable after the read; README sec.2 and the
per-member tables (keep, but generate).

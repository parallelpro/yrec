# state/ — the star's type and the layouts it is built from

This directory defines the star and the records and slot constants it is made of. `star_info_lib` defines the type `star_info` and the two input bundles nested in it, `star%ctrl` (type `controls_state`) and `star%job` (type `star_job`); it declares no instance: `run_yrec` (`core/run/run_driver.f90`) owns the run's star and hands it to every star-layer procedure as the argument `s` (`docs/GUIDELINES.md`, "The star object: six rules"). The physics domains never use `star_info_lib`: they take the slot constants from the constants-only layout modules, and their configuration from the settings records `setup/` builds. The components of `controls_state` (`controls_state_def.inc`) and the controls of `star_job` (`job_controls_def.inc`) are generated from `io/controls_buffer.f90` and the registry's `home` column by `tools/gen_controls_state.py`, which also writes the copy module `io/controls_sync.f90`; this README does not list them one by one. `star%ctrl` never changes after the controls read; `star%job`'s generated part (`job_controls_def.inc`) holds exactly the controls that code outside the read assigns (the run-card arrays the calibration rewrites, the stop thresholds the stop conditions disarm, `atm_choice`, `fcorr`, the latches `use_structure_dt_limits` and `use_mass_accretion`, the Monte-Carlo range), as the registry's `home` column says and `tools/check_boundaries.py` checks against the code; its hand-declared part holds the table paths, the surface-mixture weights, the Monte-Carlo samples and the run-card cursor `nk`, which no check covers yet. A latch that the run switches and every card re-arms is a star member armed from its control (`settling_active`, from `ctrl%diffuse_helium_active`). The controls are documented in `src/defaults/controls_registry.md` and, with their defaults, in `src/defaults/physics.defaults` and `src/defaults/control.defaults`.

| Module (file) | What it holds |
|---|---|
| `star_info_lib` | The types `star_info`, `star_job` and `controls_state`; the capacities `max_convective_zones`, `max_radiative_zones` and `max_mc_runs`; the slots of the energy channels (`i_eps_*`) and of the luminosity channels (`i_lum_*`); it re-exports the layouts below. Outside the type: the namelist-path overrides of the C API (`yrec_capi`) and the flag `yrec_reset` sets to re-initialize `evolve_step`'s saved locals. |
| `composition_layout` | The size `json` = 5000 of every per-zone array; the rows of `xa`, `i_h1` ... `i_be9` (row `i_metals` is the total Z), and the 12 surface-mixture slots `ix_*`. |
| `neutrino_flux_layout` | The slots `i_nu_*` of the neutrino-flux vectors. |
| `reaction_layout` | The slots of net/'s per-reaction arrays: the reactions `r_*` (`num_reactions`), the composition rate vector's branching rows `rv_*` (`num_rate_rows`), the extra cross-section scales `xs_*` (`num_extra_scales`) and `num_sstandard`. net/'s `net_configuration` re-exports them; `star_info_lib` sizes the star's scale arrays with them. |
| `transport_mechanism_layout` | The slots `i_mech_*` of rotation's transport mechanisms, with each one's name, the slot it corrects and its two controls (`star%mechanism_velocity`, `rotation/README.md`). |
| `envelope_points` | The records of the outer layers the envelope integration fills: the envelope and the atmosphere point by point, and the atmosphere integration's current point (`star%env_points`, `star%atm_points`, `star%atm_point`). |
| `stitched_model_state` | The record of the stitched model, the converged star on one grid from the centre to the atmosphere (`star%stitched`). |
| `physical_constants` | The mathematical and physical constants, compile-time parameters (CODATA 1986, with their sources), except log10 G (`cgl`): it starts at the CODATA value, and `core/run/run_setup.f90` sets it once per run from the control `G_cgs` (`log10_newton_g`). Only the star layer reads `cgl`: the physics domains, rotation and the star-blind kernels take log10 G from their settings records (`log10_gravitational_constant`), which `tools/check_boundaries.py` checks. The version strings are `core/run/yrec_version.f90`'s. |
| `control_slots` | The slot names of the controls that are vectors: `atime(itime_*)`, the timestep-control factors, and `chi_grid_scale(ichi_*)`, the zoning tolerances (`io/controls_buffer.f90` re-exports them for the controls read). |
| `controls_state_def.inc`, `job_controls_def.inc` | Generated (see above): the components of `star%ctrl` and the controls of `star%job`. |

The rest of this README describes the members of `star_info`: what each one is, its unit, who writes and reads it, and when its value can be trusted. Sections 2 and 3 are generated (below).

---

## 1. Conventions

- **Points.** A per-zone array has `json` = 5000 entries, of which `1:s%nz` are the model. Point 1 is the innermost point (not the centre itself), point `nz` the fitting point where the envelope integration takes over. The envelope and atmosphere above it are in the star only as the records of the latest envelope integration (`env_points`, `atm_points`) and in the stitched model (`stitched`), which `core/diagnostics/stitched_model.f90` builds for output.
- **Units.** CGS unless the name carries its unit (`_lsun`, `_msun`, `_rsun`, `_gyr`, `_yr`, `_days`, `_kms`, `_uhz`, `_cgs`); `log` means log₁₀. The structure variables are log₁₀ P, log₁₀ T, log₁₀ r and L in L☉ against log₁₀ m. The unit follows the quantity and a lifetime qualifier follows the unit (`luminosity_lsun_start`, `r_rsun_prev_model`); `star_mass` (M☉) keeps MESA's name; the tables give every unit, with `-` for a dimensionless quantity or an index and `log10 g` for log₁₀ of a quantity in grams.
- **Lifetime.** Most members have no explicit lifetime: the same storage is rewritten at different stages of a step. `star_info_lib.f90` declares the members grouped by lifetime class (section 2), longest-lived first, and each class's banner says when its values are valid.
- **A second run in one process.** `yrec_reset` gives every `run_yrec` call a copy of `star_template`, a static star that nothing writes, so every member starts a run at its declaration default (zero, `.false.` or blank where none is given).
- **One source.** Each member's meaning, unit and notes are the comment above its declaration in `star_info_lib.f90`, in the form `! <meaning> [<unit>]; <writers>. <notes>`; `tools/gen_star_docs.py` builds sections 2 and 3 from those comments and from the code. Edit the comments, not the tables.

## 2. Lifetime classes

<!-- BEGIN GENERATED classes: tools/gen_star_docs.py from state/star_info_lib.f90 -->
| Class | Members | What it is | When its value is valid |
|---|---|---|---|
| [Input bundles](#31-input-bundles) | 2 | The controls (`ctrl`) and the run's job (job), read by `read_controls` at the start of every `run_yrec` call. Their members are documented in `src/defaults/controls_registry.md` and the `*.defaults` files, not here. | From the controls read on; `ctrl` is not changed after it, job is changed by the run (see the type `star_job`). |
| [Constants set at setup](#32-constants-set-at-setup) | 17 | Values fixed for the whole run by the controls read, the controls' derivation or the setup before the first card: the solar units, the nuclear rate scales, a few atmosphere and opacity settings. | After `read_controls` and `prepare_run`; a Monte-Carlo run rescales the solar luminosity and the first rate scales. |
| [Run and card bookkeeping](#33-run-and-card-bookkeeping) | 20 | Per-run and per-card state: the card's mixing length and mass, the atmosphere switch, the surface mixture of the first card, the solar and star calibration loops, the reason the run ended. | From the card (or run) that sets them. |
| [Model](#34-model) | 29 | The converged model: the mesh, the structure, the composition, the global scalars, and the envelope fit that is its outer boundary condition. Persists from step to step; the core of it is written to the .mod file (`io/model_files.f90`; "in .mod" below). Point 1 is the innermost point, point `nz` the fitting point. | Between steps, and inside a step for the model being solved; the Newton corrections overwrite the structure in place, so during a step it holds the current iterate. |
| [Rotation state](#35-rotation-state) | 12 | The rotation curve and the shape of the rotating star. The specific angular momentum `j_rot` is conserved; `omega`, the moments of inertia and the shape factors follow from it and the structure. | Between steps. Meaningful only with `ctrl%rotation_active`; without rotation the shape factors are 1. |
| [Regions of the current model](#36-regions-of-the-current-model) | 9 | Point indices and bounds of the convective core, the surface convection zone, the hydrogen-burning shell and the mixed regions. | For the model they were located on; located again at each stage that needs them (the composition step, the rezone, after convergence). |
| [Step control](#37-step-control) | 8 | The timestep and the flags that steer one step's solve. | During and between steps. |
| [Step snapshot](#38-step-snapshot) | 37 | Copies of model quantities at the start of the step (&lt;member&gt;_start), of the previous model (&lt;member&gt;_prev) or of the structure-rotation iteration (&lt;member&gt;_iteration_start), on the step's mesh. The rezone interpolates the old model onto the new mesh into the `*_start` arrays and then copies them into the live arrays, so after the rezone they hold the start-of-step model on the new mesh; without a rezone (`ctrl%helium_flash_active`) `rezone_or_snapshot` copies the structure into the `*_start` arrays itself. | From the start-of-step copy to the end of the step. |
| [Step working state](#39-step-working-state) | 26 | Values one stage computes for a later stage of the same or the next step: rates of change, burning rates, the depth of the surface mixed region, the turnover time, the shell physics that rotation uses, the outer layers of the latest envelope integration, the last rotation sub-step's `omega`. A few are printed too (the history prints the turnover time and the photospheric pressure). | From the stage that writes them; each member says which. |
| [Newton](#310-newton) | 6 | Arrays of the Henyey solve, rewritten on every iteration. | During the solve; `log_pressure_delta` and `log_temperature_delta`, the changes over the step, until the next step's solve starts. |
| [Diagnostics written during the solve](#311-diagnostics-written-during-the-solve) | 22 | Per-zone physics stored for output by `compute_henyey_coefficients` at every shell on every Newton iteration (energy channels, EOS and opacity values, gradients, pulsation derivatives). Rule 5 of the six (diagnostics written once, after convergence) is the planned replacement. | After convergence they hold the last iteration's values, the structure that iteration linearized, before its correction when one was applied. Some are also rewritten by `compute_shell_physics`; the last writer wins. |
| [Diagnostics written after convergence](#312-diagnostics-written-after-convergence) | 34 | Observables computed once per model by `compute_observables` (after the stitched model is built), and the rotational transport's coefficients and velocities kept for the profiles. | From `compute_observables` (or the rotation step) until the next model; zero at a run start. |
| [Private to one module](#313-private-to-one-module) | 14 | Values kept on the star for one module alone, which writes and reads them: the state it carries between its calls (the calibrations' previous checks, models and corrections, the central EOS call's Saha state, the Newton solve's boundary rows and convergence log, the surface-mixing search's starting zone). `tools/gen_star_docs.py` places every member that one module alone uses here, and nothing else. | Only inside that module; each member says when. |

236 members in all, 6 of them records (the two input bundles among them).
<!-- END GENERATED classes -->

## 3. The members

How to read the tables:

- **member** gives the declaration's shape; `json` is the zone dimension (use `1:nz`).
- **written by** names the procedures in which `s%member` is the target of an assignment, of a `read` statement or of an `allocate`, or an actual argument to a dummy argument declared `intent(out)`, or `intent(inout)` that the callee (directly, in its internal procedures or through its own callees) assigns. So a member written inside a callee appears under its caller when the caller passes `s%member` itself: `stored_vertex_index` is written by `compute_surface_boundary` through `refresh_surface_bc`'s argument, and the table names `step_driver` (`solve_level`). With more than six procedures only the modules are listed.
- **read by** names the modules where the value of `s%member` is used: in an expression, a subscript or a loop bound, or passed to a dummy argument that the callee reads. Passing it to an `intent(out)` dummy, or to an `intent(inout)` one that the callee never reads, is not a read. With more than seven modules the count is given by directory.
- Module names are the names in the files' `module` statements, which are the file names (`core/step/step_driver.f90` holds `step_driver`). Tests (`*/test/`) and `retired/` code are not counted.
- **meaning**, **unit** and **notes** are the members' comments in `star_info_lib.f90`. "In .mod" marks what the `.mod` model file (`io/model_files.f90`) holds.
- The tables are generated by `python3 tools/gen_star_docs.py` (in `src/`); `test_generated_files.py` runs it with `--check`, which fails when the tables are stale, when a member has no comment in the form above, names a writer that does not write it, or is never written or never read outside the tests. To check one row by hand, `grep -rn '%member\b' src --include='*.f90'`.

<!-- BEGIN GENERATED members: tools/gen_star_docs.py from state/star_info_lib.f90 -->
### 3.1 Input bundles

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `ctrl` (type `controls_state`) | the namelist controls the run reads, by their canonical names | - | `controls_reader` (read_controls); `controls_sync` (store_controls_to_star) | 33 modules (core/ 23, io/ 7, setup/ 3) | Reset to the declared defaults before each read; members: `src/defaults/controls_registry.md`, generated into `controls_state_def.inc` |
| `job` (type `star_job`) | the run's job: the run-card arrays and latches (the controls the run changes), the Monte-Carlo samples and the run-card cursor nk | - | `envelope_integration`, `model_rescaling`, `starting_model`, `calibration`, `monte_carlo`, `run_driver`, `mass_change`, `step_driver`, `henyey_iteration`, `surface_boundary`, `controls_reader`, `controls_sync` (22 procedures) | 20 modules (core/ 14, io/ 5, setup/ 1) | Changed by the run: the calibrations rewrite the next cards, the stop conditions disarm thresholds, the atmosphere switch changes `atm_choice`, the solver ramps fcorr |

### 3.2 Constants set at setup

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `solar_luminosity_cgs` | Lsun, a copy of `ctrl%solar_luminosity_cgs` (legacy clsun) | erg/s | `monte_carlo` (apply_monte_carlo_parameters); `run_setup` (set_up_constants_and_tables) | 12 modules (core/ 6, io/ 3, setup/ 2, rotation/ 1) | A Monte-Carlo run scales it |
| `log10_solar_luminosity` | log10 Lsun | log10 erg/s | `monte_carlo` (apply_monte_carlo_parameters); `run_setup` (set_up_constants_and_tables) | `observables`, `envelope_integration`, `starting_model`, `run_driver`, `mass_change`, `microphysics_settings`, `rotation_settings_setup` |  |
| `ln10_over_lsun` | ln(10)/Lsun, not a logarithm: times the enclosed mass it scales the luminosity-equation terms | s/erg | `monte_carlo` (apply_monte_carlo_parameters); `run_setup` (set_up_constants_and_tables) | `structure_equations` |  |
| `solar_mass_cgs` | Msun, a copy of `ctrl%solar_mass_cgs` | g | `run_setup` (set_up_constants_and_tables) | 10 modules (core/ 7, io/ 1, rotation/ 1, setup/ 1) |  |
| `solar_radius_cgs` | Rsun, a copy of `ctrl%solar_radius_cgs` (legacy crsun) | cm | `run_setup` (set_up_constants_and_tables) | 8 modules (core/ 6, io/ 1, setup/ 1) |  |
| `log10_solar_radius` | log10 Rsun | log10 cm | `run_setup` (set_up_constants_and_tables) | `observables`, `run_driver`, `yrec_output` |  |
| `cross_section_scale(num_reactions)` | per-reaction scale of the cross section S0 against Bahcall & Ulrich (1988), by `reaction_layout`'s slot `r_*` | - | `monte_carlo` (apply_monte_carlo_parameters); `controls_reader` (adopt_canonical_names); `user_inputs` (set_cross_section_scales) | `net_settings_setup` | Those of pp, He3 + He3 and He3 + He4 are rescaled per Monte-Carlo run |
| `qs0e_scale(num_reactions)` | scale of each reaction's S'/S against Solar Fusion I, by reaction slot `r_*` | - | `user_inputs` (reset_rate_scale_defaults, set_first_derivative_scales) | `net_settings_setup` |  |
| `qqs0ee_scale(num_reactions)` | scale of each reaction's S''/S against Solar Fusion I, by reaction slot `r_*` | - | `user_inputs` (reset_rate_scale_defaults, set_second_derivative_scales) | `net_settings_setup` |  |
| `extra_cross_section_scale(num_extra_scales)` | scale of the rates outside the reaction table, by `reaction_layout`'s extra slot `xs_*` (pep, Be7 + e, Be7 + p, hep) | - | `monte_carlo` (apply_monte_carlo_parameters); `controls_reader` (adopt_canonical_names); `user_inputs` (set_cross_section_scales) | `net_settings_setup` | Be7 + p's is rescaled per Monte-Carlo run |
| `be7p_qs0e_scale` | scale of Be7 + p's S'/S against Solar Fusion I | - | `user_inputs` (reset_rate_scale_defaults, set_first_derivative_scales) | `net_settings_setup` |  |
| `be7p_qqs0ee_scale` | scale of Be7 + p's S''/S against Solar Fusion I | - | `user_inputs` (reset_rate_scale_defaults, set_second_derivative_scales) | `net_settings_setup` |  |
| `o16_gamma_scale` | scale of the N15(p,gamma)O16 branch | - | `user_inputs` (reset_rate_scale_defaults, set_cross_section_scales) | `net_settings_setup` |  |
| `c12_alpha_scale` | scale of the N15(p,alpha)C12 branch | - | `user_inputs` (reset_rate_scale_defaults, set_cross_section_scales) | `net_settings_setup` |  |
| `envelope_store_log_temperature` | log T at which an envelope integration stores its state for the surface boundary (`stored_envelope_state`): the mean of `ctrl%tenv0` and `ctrl%tenv1` | log10 K | `derived_controls` (derive_controls) | `envelope_integration` |  |
| `atm_choice_initial` | `atm_choice` as read from the deck | - | `controls_reader` (adopt_canonical_names) | `surface_boundary` | The run may change `job%atm_choice` (see `use_ttau_relation`); this keeps the original |
| `use_two_z_tables` | load a second opacity-table set at another Z (the Z ramp or Z diffusion) | logical | `run_driver` (run_yrec) | `microphysics_settings` |  |

### 3.3 Run and card bookkeeping

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `mixing_length_alpha` | mixing-length parameter of the current run card | - | `run_driver` (begin_kind_card) | `starting_model`, `history_output`, `model_files`, `yrec_output`, `microphysics_settings` | From `job%mixing_length_array(nk)` at each card start (the solar calibration changes the later cards' entries); in .mod; 1.4 before the first card |
| `pulsation_mass_msun` | `star_mass` at the start of the card, for the profile and pulsation headers | Msun | `run_driver` (begin_kind_card) | `profile_output` |  |
| `use_ttau_relation` | the run has switched `atm_choice` to the gray atmosphere because Teff left the Kurucz or Allard table | logical | `envelope_integration` (atm_get); `surface_boundary` (switch_atmosphere_choice); `controls_reader` (adopt_canonical_names) | `surface_boundary` | Switches back when Teff returns |
| `termination_reason` | why the last run card ended, for the run summary | text | `run_driver` (end_of_card_calibration); `stop_conditions` (abundance_stop_triggered, converged_model_is_nan, init_stop_conditions, structure_limit_stop_triggered) | `stop_conditions`, `run_log` | Set to the default per card by `init_stop_conditions`, overridden by whichever stop fires |
| `settling_suspended_reported` | the settling routines have reported that settling is suspended (a fully convective model, or no hydrogen or helium to settle) | logical | `composition_step` (settle); `mid_timestep_model` (settle_helium) | `composition_step`, `mid_timestep_model` | The message prints once per suspension; "SETTLING RESUMED" clears it |
| `settling_active` | helium settles in this run card | logical | `run_driver` (begin_kind_card); `composition_step` (settle); `mid_timestep_model` (settle_helium) | `composition_step`, `mid_timestep_model` | Armed from the control `ctrl%diffuse_helium_active` (never changed after the read) at the start of every card, and so of every Monte-Carlo run, by `begin_kind_card`; the settling routines switch it off once the central X falls below `ctrl%hydrogen_diffusion_floor` |
| `envelope_hydrogen_fraction` | surface X of the first run card: the reference X of the EOS mixture and of 1/`mu_ion` in `compute_shell_physics` | - | `starting_model` (update_surface_mixture) | `model_rescaling`, `starting_model`, `shell_physics`, `run_cards`, `transport_checks`, `settings_install` | Set once (`update_surface_mixture`) and kept for the run on purpose; the controls read and `kap_init` use it earlier, while it is still 0 |
| `envelope_metal_fraction` | surface Z of the EOS mixture: the first card's, or the last Z rescale's | - | `model_rescaling` (rescale_interior_z); `starting_model` (update_surface_mixture); `run_cards` (interpret_kind_cards) | `model_rescaling`, `starting_model`, `shell_physics`, `run_cards`, `transport_checks`, `settings_install` | The controls read (`interpret_kind_cards`) also sets it from a card's rescale Z |
| `envelope_helium_fraction` | surface Y = 1 - X - Z - X(He3) of the first card | - | `starting_model` (update_surface_mixture) | `shell_physics`, `transport_checks` |  |
| `envelope_he3_fraction` | surface He3 of the first card | - | `starting_model` (update_surface_mixture) | `shell_physics`, `transport_checks` |  |
| `envelope_metals_no_cno` | the part of the first card's surface Z not in C, N and O, passed to convection's settings (`envelope_metals_no_cno`) | - | `starting_model` (update_surface_mixture) | `microphysics_settings`, `settings_install` |  |
| `envelope_inverse_atomic_weight` | inverse mean atomic weight of the surface mixture, sum `X_i`/`A_i` | 1/amu | `starting_model` (update_surface_mixture) | `model_rescaling`, `starting_model`, `shell_physics`, `run_cards`, `transport_checks` |  |
| `envelope_number_fractions(n_mix_species)` | number fractions of the 12 surface-mixture elements, slots `ix_*`, passed to the EOS (`eos_set_mixture`) | - | `starting_model` (update_surface_mixture) | `model_rescaling`, `starting_model`, `run_cards` |  |
| `dlum_dx` | d log L/dX of the last solar-calibration Newton step: calibration's parameter `dlogl_dx_fixed`, -3.78 | - | `calibration` (check_solar_calibration) | `monte_carlo` | 0 before the first step; printed in the Monte-Carlo summary, as the three below |
| `drad_dx` | d log R/dX of the last Newton step, `dlogr_dx_fixed`, -0.89 | - | `calibration` (check_solar_calibration) | `monte_carlo` |  |
| `dlum_dalpha` | d log L/d alpha of the last Newton step, `dlogl_dalpha_fixed`, 0.0139 | - | `calibration` (check_solar_calibration) | `monte_carlo` |  |
| `drad_dalpha` | d log R/d alpha of the last Newton step, `dlogr_dalpha_fixed`, -0.050 | - | `calibration` (check_solar_calibration) | `monte_carlo` |  |
| `solar_calibration_converged` | the solar calibration has converged | logical | `calibration` (check_solar_calibration, setup_solar_calibration) | `monte_carlo`, `run_driver` |  |
| `star_found_flag` | the star calibration found a track through the target L and R, so the final age-stopped run is armed | logical | `calibration` (check_star_calibration, setup_star_calibration) | `run_driver`, `stop_conditions` |  |
| `just_passed_target_radius_flag` | this model has just crossed the star calibration's target radius | logical | `calibration` (check_star_calibration, setup_star_calibration) | `stop_conditions` |  |

### 3.4 Model

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `nz` | number of points of the interior model, out to the fitting point | - | `neutrino_flux_output`, `rezoning`, `starting_model`, `step_driver`, `model_files` (10 procedures) | 23 modules (core/ 16, io/ 4, rotation/ 3) | In .mod |
| `model_number` | counter of converged models in the run; the step ages the model only when it is &gt;= 0 | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step); `model_files` (read_mod_model, read_model2, read_yrec7) | 9 modules (io/ 6, core/ 3) | In .mod; reset per card from `ctrl%set_initial_model_number` |
| `age_gyr` | age of the model | Gyr | `step_driver` (advance_composition_and_age); `model_files` (read_mod_model, read_model2, read_yrec7) | 9 modules (io/ 6, core/ 3) | In .mod |
| `star_mass` | total mass M | Msun | `starting_model` (rescale_and_refit_envelope); `mass_change` (update_masses); `model_files` (read_mod_model, read_model2, read_yrec7) | 10 modules (core/ 5, io/ 3, rotation/ 2) | In .mod; MESA's history name, without the _msun suffix |
| `log_total_mass` | log10 of the total mass M | log10 g | `starting_model` (rescale_and_refit_envelope); `mass_change` (update_masses); `model_files` (read_mod_model, read_model2, read_yrec7) | 17 modules (core/ 12, rotation/ 3, io/ 2) | In .mod |
| `log_L_lsun` | log10 of the luminosity at the fitting point, L(nz) | log10 Lsun | `henyey_iteration` (apply_henyey_corrections, refresh_surface_bc); `model_files` (read_mod_model, read_model2, read_yrec7) | 15 modules (core/ 9, io/ 4, rotation/ 2) | In .mod; reset after every Newton correction |
| `log_Teff` | log10 effective temperature, from the envelope fit (`envelope_fit_coeffs` 7-9) | log10 K | `henyey_iteration` (apply_henyey_corrections); `model_files` (read_mod_model, read_model2, read_yrec7) | 19 modules (core/ 11, io/ 6, rotation/ 2) | In .mod; reset after every Newton correction |
| `log_mass(json)` | log10 of the mass `m` inside point k, the mesh coordinate | log10 g | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `model_files` (11 procedures) | 14 modules (core/ 10, rotation/ 3, io/ 1) | In .mod |
| `m(json)` | mass inside point k | g | `rezoning` (set_new_shell_masses); `starting_model` (build_shell_masses); `composition_step` (settle); `mass_change` (update_masses); `step_driver` (evolve_step); `mid_timestep_model` (advance_rotation_substep) | 11 modules (core/ 7, rotation/ 3, io/ 1) | Rebuilt from `log_mass` after every rezone or read |
| `dm(json)` | mass of the shell around point k, half the span to the neighbouring points; `dm(nz)` reaches the total mass | g | `rezoning` (set_new_shell_masses); `starting_model` (build_shell_masses); `mass_change` (update_masses); `step_driver` (evolve_step) | 13 modules (core/ 10, rotation/ 3) |  |
| `luminosity_lsun(json)` | luminosity L(m) through point k | Lsun | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `model_files` (11 procedures) | 15 modules (core/ 11, io/ 2, rotation/ 2) | In .mod; a Newton variable (linear, not log) |
| `logR(json)` | log10 radius | log10 cm | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `henyey_iteration`, `model_files` (12 procedures) | 17 modules (core/ 13, io/ 2, rotation/ 2) | In .mod; a Newton variable |
| `logP(json)` | log10 pressure | log10 dyn/cm^2 | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `henyey_iteration`, `model_files` (12 procedures) | 17 modules (core/ 13, io/ 2, rotation/ 2) | In .mod; a Newton variable |
| `logT(json)` | log10 temperature | log10 K | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `henyey_iteration`, `model_files` (12 procedures) | 18 modules (core/ 14, io/ 2, rotation/ 2) | In .mod; a Newton variable |
| `logRho(json)` | log10 density, from the EOS at (P, T, X, Z) | log10 g/cm^3 | `rezoning`, `starting_model`, `composition_step`, `step_driver`, `structure_equations`, `model_files` (12 procedures) | 17 modules (core/ 13, io/ 2, rotation/ 2) | In .mod; rewritten from the EOS on every Newton iteration, at the structure that iteration linearized, so it can lag the last correction |
| `convective_flag(json)` | true where the zone is convective by the temperature-gradient test | logical | `rezoning`, `starting_model`, `composition_step`, `step_driver`, `structure_equations`, `model_files`, `evolve_rotation` (16 procedures) | 11 modules (core/ 8, rotation/ 2, io/ 1) | In .mod; rewritten on every Newton iteration and by `compute_shell_physics` |
| `xa(n_species_extended,json)` | mass fractions; row = species slot `i_*` of `composition_layout`, row `i_metals` the total Z | - | `rezoning`, `starting_model`, `burn_step`, `composition_step`, `mass_change`, `step_driver`, `model_files`, `evolve_rotation`, `mid_timestep_model`, `secular_transport` (30 procedures) | 21 modules (core/ 14, io/ 4, rotation/ 3) | In .mod; rows 12-15 only with `ctrl%use_extended_composition` |
| `surface_hydrogen_fraction` | surface X for the envelope integrations and the surface opacity tables; follows `xa(i_h1,nz)` | - | `model_rescaling` (rescale_surface_x); `rezoning` (refresh_surface_opacity); `starting_model` (interpolate_new_fitting_point, reset_state_after_read) | `model_rescaling`, `rezoning`, `starting_model`, `henyey_iteration` | Refreshed when it moves by more than 1e-8 (`refresh_surface_opacity`), or set by a rescale |
| `surface_metal_fraction` | surface Z for the envelope integrations, updated with `surface_hydrogen_fraction` | - | `model_rescaling` (rescale_interior_z); `rezoning` (refresh_surface_opacity); `starting_model` (interpolate_new_fitting_point, reset_state_after_read) | `model_rescaling`, `starting_model`, `henyey_iteration` |  |
| `log_fit_point_mass_fraction` | log10(`m(nz)`/M), the mass fraction inside the fitting point | log10 | `envelope_rebuild`, `rezoning`, `starting_model`, `mass_change`, `henyey_iteration` (7 procedures) | `envelope_integration`, `envelope_rebuild`, `envelope_refit`, `rezoning`, `starting_model` | Recomputed at every solve level (`refresh_surface_bc`) |
| `stored_envelope_state(4)` | (log P, log T, log r, log `m`/M) where a vertex's envelope integration first reached `envelope_store_log_temperature` | log10 cgs | `henyey_iteration` (refresh_surface_bc) | `rezoning` | The rezone uses it to add a point at the surface |
| `envelope_fit_coeffs(9)` | coefficients c(1:9) of the envelope fit at the fitting point: log R, log L and log Teff each = c log P + c log T + c | - | `henyey_iteration` (refresh_surface_bc); `model_files` (read_model2, read_yrec7) | `henyey_iteration` | Rebuilt by `compute_surface_boundary` when the triangle changes; read from YREC7/MODEL2 files, not written to .mod |
| `triangle_log_l_lsun(3)` | log10 L of the three vertices of the envelope-fit triangle | log10 Lsun | `henyey_iteration` (refresh_surface_bc); `model_files` (read_model2, read_yrec7) | `henyey_iteration` |  |
| `triangle_log_teff(3)` | log10 Teff of the three vertices of the envelope-fit triangle | log10 K | `henyey_iteration` (refresh_surface_bc); `model_files` (read_model2, read_yrec7) | `henyey_iteration` |  |
| `fit_point_pressure(3)` | log10 P at the fitting point of each vertex's envelope integration | log10 dyn/cm^2 | `henyey_iteration` (refresh_surface_bc); `model_files` (read_model2, read_yrec7) | `rezoning`, `henyey_iteration` |  |
| `fit_point_temperature(3)` | log10 T at the fitting point of each vertex | log10 K | `henyey_iteration` (refresh_surface_bc); `model_files` (read_model2, read_yrec7) | `henyey_iteration` |  |
| `fit_point_radius(3)` | log10 r at the fitting point of each vertex | log10 cm | `henyey_iteration` (refresh_surface_bc); `model_files` (read_model2, read_yrec7) | `rezoning`, `henyey_iteration` |  |
| `stored_vertex_index` | which triangle vertex stored `stored_envelope_state` (0: none) | - | `run_driver` (begin_kind_card); `step_driver` (reload_model_if_diverged, solve_level) | `step_driver` |  |
| `triangle_orientation` | orientation (+1 or -1) of the envelope-fit triangle | - | `run_driver` (begin_kind_card); `step_driver` (reload_model_if_diverged, solve_level) | `step_driver` |  |

### 3.5 Rotation state

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `omega(json)` | angular velocity of each shell | rad/s | `envelope_rebuild`, `rezoning`, `starting_model`, `henyey_iteration`, `model_files`, `evolve_rotation` (14 procedures) | 14 modules (core/ 10, io/ 3, rotation/ 1) | In .mod |
| `j_rot(json)` | specific angular momentum of each shell, the conserved quantity; `omega` follows from it (`omega_from_j`) | cm^2/s | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `henyey_iteration`, `evolve_rotation`, `mid_timestep_model`, `secular_transport` (17 procedures) | 8 modules (core/ 5, rotation/ 3) |  |
| `i_rot(json)` | moment of inertia of each shell | g cm^2 | `rezoning` (respin_on_new_grid); `starting_model` (initialize_rotation_state); `step_driver` (evolve_step); `henyey_iteration` (update_rotation_after_corrections); `evolve_rotation` (evolve_angular_momentum, spin_down_without_transport) | `observables`, `stitched_model`, `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `evolve_rotation` |  |
| `kinetic_energy_rot(json)` | rotational kinetic energy of each shell, `omega` j `dm`/2 | erg | `rezoning` (retally_angular_momentum); `starting_model` (initialize_rotation_state); `step_driver` (evolve_step); `henyey_iteration` (update_rotation_after_corrections) | `starting_model`, `step_driver`, `structure_equations` |  |
| `fp_rot(json)` | rotational correction `f_P` to the hydrostatic equation, 1 without rotation | - | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration` (7 procedures) | `stitched_model`, `rezoning`, `starting_model`, `henyey_iteration`, `structure_equations`, `evolve_rotation` | Set from `omega` by `rotation_shape_factors` |
| `ft_rot(json)` | rotational correction `f_T` to the radiative-transfer equation, 1 without rotation | - | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration` (7 procedures) | `stitched_model`, `rezoning`, `starting_model`, `henyey_iteration`, `structure_equations`, `evolve_rotation` |  |
| `mean_radius(json)` | r0, the mean radius of the shell's level surface | cm | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `evolve_rotation` (8 procedures) | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration` |  |
| `eta_squared(json)` | Radau's `eta_2` of each shell, the shape of its level surface | - | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `evolve_rotation` (8 procedures) | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `evolve_rotation` |  |
| `mean_gravity(json)` | effective gravity averaged over the level surface, &lt;g&gt; | cm/s^2 | `rezoning` (respin_on_new_grid); `starting_model` (initialize_rotation_state); `step_driver` (converge_with_rotation, evolve_step); `henyey_iteration` (update_rotation_after_corrections) | `rezoning`, `mid_timestep_model` |  |
| `am_transport_convective_flag(json)` | true where a zone is convective for angular-momentum transport, overshoot included | logical | `rezoning` (respin_on_new_grid); `starting_model` (read_starting_model) | `rezoning`, `henyey_iteration` | From `am_convective_regions` at the model read and every rezone, i.e. before every structure solve; `henyey_iterate`'s `omega_from_j` holds those zones in solid-body rotation during the Newton iterations |
| `total_angular_momentum` | total angular momentum J of the star | g cm^2/s | `run_driver` (begin_kind_card); `step_driver` (evolve_step, reload_model_if_diverged, rezone_or_snapshot); `evolve_rotation` (total_rotation_state) | `run_driver`, `step_driver`, `history_output`, `evolve_rotation` | Zeroed at each card start, then totalled by the model read and every rezone |
| `total_rotational_ke` | total rotational kinetic energy | erg | `run_driver` (begin_kind_card); `step_driver` (evolve_step, reload_model_if_diverged, rezone_or_snapshot); `evolve_rotation` (total_rotation_state) | `run_driver`, `step_driver`, `history_output`, `evolve_rotation` |  |

### 3.6 Regions of the current model

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `core_cz_top_index` | outermost point of the convective core (1: none) | - | `rezoning` (locate_new_cz_edges); `run_driver` (begin_kind_card); `step_driver` (advance_composition_and_age, burn_light_elements, evolve_step) | `observables`, `rezoning`, `run_driver`, `step_driver` | Set by the composition step's predictor pass, the rezone, `burn_light_elements` and, after convergence, `locate_shell_boundaries` |
| `envelope_cz_bottom_index` | innermost point of the surface convection zone (`nz`: none) | - | `rezoning` (locate_new_cz_edges); `run_driver` (begin_kind_card); `step_driver` (advance_composition_and_age, burn_light_elements, evolve_step) | `observables`, `rezoning`, `starting_model`, `run_driver`, `mass_change`, `step_driver` | Set as `core_cz_top_index` is |
| `envelope_overshoot_bottom_index` | first point of the overshoot region below the surface convection zone (-1: none) | - | `rezoning`, `starting_model`, `composition_step`, `step_driver`, `evolve_rotation`, `mid_timestep_model` (7 procedures) | `shell_physics`, `structure_equations` | From `find_convection_zones`; with `ctrl%ladov` the radiative points `envelope_overshoot_bottom_index` to `envelope_overshoot_top_index` take the adiabatic gradient |
| `envelope_overshoot_top_index` | last point of that overshoot region (-1: none) | - | `rezoning`, `starting_model`, `composition_step`, `step_driver`, `evolve_rotation`, `mid_timestep_model` (7 procedures) | `shell_physics`, `structure_equations` |  |
| `mixed_zone_bounds(max_convective_zones,2)` | first and last point of each mixed (convective plus overshoot) region, up to `max_convective_zones` of them | - | `composition_step` (advance_composition); `step_driver` (burn_light_elements) | `burn_step`, `composition_step`, `step_driver` | Rewritten by every composition-step pass and by `burn_light_elements` |
| `has_h_shell` | the model burns hydrogen in a shell rather than in the core | logical | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `observables`, `step_driver`, `yrec_output` |  |
| `h_shell_bot_index` | first point of the H-burning shell, where X &gt; `ctrl%atime(1)` | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `observables`, `run_driver`, `step_driver`, `yrec_output` | Valid when `has_h_shell` |
| `h_shell_mid_index` | first point where X exceeds half its surface value | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `observables`, `run_driver`, `step_driver` | Valid when `has_h_shell` |
| `h_shell_top_index` | last point of the H-burning shell | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `observables` | Valid when `has_h_shell` |

### 3.7 Step control

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `dt` | timestep of the step being taken (the next step's after `compute_timestep`) | s | `run_driver` (begin_kind_card); `mass_change` (limit_mass_change_timestep); `step_driver` (evolve_step, reload_model_if_diverged) | `run_driver`, `stop_conditions`, `mass_change`, `step_driver` | Negative means solve without aging (a rescale or relaxation); made positive after convergence |
| `dt_saved` | copy of `dt`: when the end-age stop clips the next timestep to zero, the timestep just taken is restored from it | s | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `step_driver` |  |
| `timestep_yr` | `dt` in years | yr | `run_driver` (begin_kind_card); `step_driver` (advance_composition_and_age, evolve_step, reload_model_if_diverged) | `run_driver`, `step_driver`, `model_files`, `profile_output`, `run_log` | In .mod |
| `max_domega_frac` | largest fractional change of `omega` over the step | - | `run_driver` (begin_kind_card); `step_driver` (converge_with_rotation, evolve_step) | `step_driver` | Written by `compute_timestep` and `evolve_angular_momentum` |
| `model_diverged_flag` | the structure solve of this step diverged; a set flag stops the run on the next pass | logical | `run_driver` (begin_kind_card); `step_driver` (converge_with_rotation, solve_level, solve_structure) | `run_driver`, `step_driver` |  |
| `punch_pending_flag` | a model has been computed, so the next helium-flash card saves its starting model | logical | `run_driver` (run_yrec); `step_driver` (reload_model_if_diverged) | `run_driver` |  |
| `recompute_envelope_triangle` | build a new envelope-fit triangle at the next solve | logical | `run_driver` (begin_kind_card); `mass_change` (update_masses); `step_driver` (evolve_step, reload_model_if_diverged, solve_level, solve_structure) | `step_driver` | Set by a model read, the mass change, the envelope rebuild or `ctrl%lnew0`; cleared by `compute_surface_boundary` |
| `reset_triangle` | the rezone moved the outer point: redo the triangle at level 2 | logical | `run_driver` (begin_kind_card); `step_driver` (rezone_or_snapshot, solve_level) | `step_driver` |  |

### 3.8 Step snapshot

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `logP_start(json)` | log P at the start of the step, on the current mesh | log10 dyn/cm^2 | `rezoning` (regrid_structure_variables, seed_new_point_run); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `rezoning`, `timestep_control`, `mid_timestep_model` |  |
| `logT_start(json)` | log T at the start of the step | log10 K | `rezoning` (regrid_structure_variables); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `rezoning`, `burn_step`, `timestep_control`, `mid_timestep_model` |  |
| `logR_start(json)` | log r at the start of the step | log10 cm | `rezoning` (regrid_structure_variables); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `rezoning`, `convection_geometry`, `step_driver`, `timestep_control`, `mid_timestep_model` |  |
| `luminosity_lsun_start(json)` | L at the start of the step | Lsun | `rezoning` (regrid_structure_variables, seed_new_point_run); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `rezoning`, `step_driver`, `timestep_control`, `mid_timestep_model` |  |
| `logRho_start(json)` | log rho at the start of the step | log10 g/cm^3 | `rezoning` (regrid_structure_variables); `step_driver` (rezone_or_snapshot) | `rezoning`, `burn_step`, `mid_timestep_model` | Not set by `snapshot_step_start` at a card start |
| `xa_start(n_species_extended,json)` | composition at the start of the step | - | `rezoning` (regrid_composition_terms); `step_driver` (advance_composition_and_age, converge_with_rotation); `evolve_rotation` (burn_light_elements_over_substep, evolve_angular_momentum) | `rezoning`, `burn_step`, `step_driver`, `timestep_control`, `evolve_rotation`, `mid_timestep_model` | Copied before the composition step; the corrector passes and the structure-rotation iteration restart from it |
| `convective_flag_start(json)` | `convective_flag` at the start of the step, rotating models only | logical | `rezoning` (retally_angular_momentum) | `mid_timestep_model` |  |
| `log_Teff_start` | log Teff at the start of the step, for the HR-motion timestep limit | log10 K | `rezoning` (seed_new_point_run); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `timestep_control` |  |
| `surface_luminosity_lsun_start` | L(nz) at the start of the step, before the rezone, for the HR-motion timestep limit | Lsun | `rezoning` (seed_new_point_run); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `timestep_control` |  |
| `kinetic_energy_rot_start(json)` | `kinetic_energy_rot` at the start of the step, for d(E_rot)/`dt` in the energy equation of rotating models | erg | `step_driver` (rezone_or_snapshot) | `structure_equations` |  |
| `omega_start(json)` | `omega` at the start of the step | rad/s | `rezoning` (regrid_rotation_variables, retally_angular_momentum); `starting_model` (snapshot_step_start) | `rezoning`, `timestep_control`, `evolve_rotation` | Read by the `omega` timestep limit and the transport |
| `j_rot_start(json)` | `j_rot` at the start of the step, the rezone's interpolation buffer | cm^2/s | `rezoning` (regrid_rotation_variables, retally_angular_momentum) | `rezoning` |  |
| `mean_gravity_start(json)` | `mean_gravity` at the start of the step | cm/s^2 | `rezoning` (retally_angular_momentum) | `mid_timestep_model` |  |
| `mean_radius_start(json)` | `mean_radius` at the start of the step, the rezone's interpolation buffer | cm | `rezoning` (regrid_rotation_variables, retally_angular_momentum) | `rezoning` |  |
| `eta_squared_start(json)` | `eta_squared` at the start of the step, the rezone's interpolation buffer | - | `rezoning` (regrid_rotation_variables, retally_angular_momentum) | `rezoning` |  |
| `gradT_start(json)` | `gradT` at the start of the step, from which the rotation sub-steps interpolate in time (`mid_timestep_model`) | - | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` | The rezone stores the twelve physics arrays from here on: the two energy sums in `regrid_start_of_step_eps`, the rest in `rerun_physics_on_new_grid` |
| `grada_start(json)` | `grada` at the start of the step | - | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` |  |
| `mu_start(json)` | `mu` at the start of the step | - | `rezoning` (rerun_physics_on_new_grid) | `evolve_rotation`, `mid_timestep_model` |  |
| `opacity_zone_start(json)` | `opacity_zone` at the start of the step | cm^2/g | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` |  |
| `cp_start(json)` | `cp` at the start of the step | erg/(g K) | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` |  |
| `dlnrho_dlnt_start(json)` | `dlnrho_dlnt` (d ln rho/d ln T at constant P) at the start of the step | - | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` |  |
| `conv_vel_start(json)` | `conv_vel` at the start of the step | cm/s | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` |  |
| `kinematic_viscosity_start(json)` | `kinematic_viscosity` at the start of the step | cm^2/s | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` |  |
| `thermal_diffusivity_start(json)` | `thermal_diffusivity` at the start of the step | cm^2/s | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` |  |
| `eps_total_start(json)` | `eps_total` at the start of the step, splined onto the step's mesh (`regrid_start_of_step_eps`) | erg/(g s) | `rezoning` (regrid_start_of_step_eps) | `mid_timestep_model` |  |
| `gradr_start(json)` | `gradr` at the start of the step | - | `rezoning` (rerun_physics_on_new_grid) | `mid_timestep_model` |  |
| `eps_plus_neu_grav_start(json)` | `eps_total` + `eps_nu` + `eps_grav` at the start of the step | erg/(g s) | `rezoning` (regrid_start_of_step_eps) | `mid_timestep_model` |  |
| `rate_li6_start(json)` | Li6 proton-capture rate rho X lambda at the start of the step; zero outward of log T = 6 | 1/s | `burn_step` (store_li_be_rates); `evolve_rotation` (burn_light_elements_over_substep) | `burn_step` | With `ctrl%use_extended_composition` |
| `rate_li7_start(json)` | Li7 proton-capture rate at the start of the step | 1/s | `burn_step` (store_li_be_rates); `evolve_rotation` (burn_light_elements_over_substep) | `burn_step` |  |
| `rate_be9_start(json)` | Be9 proton-capture rate at the start of the step | 1/s | `burn_step` (store_li_be_rates); `evolve_rotation` (burn_light_elements_over_substep) | `burn_step` |  |
| `deuterium_burning_rate_start(json)` | deuterium burning rate at the start of the step, stored by the composition step's predictor pass | 1/s | `rezoning` (regrid_composition_terms); `burn_step` (compute_burn_rates) | `rezoning`, `burn_step`, `mid_timestep_model` |  |
| `envelope_overshoot_depth_start` | envelope-overshoot depth `alpha_ov` `H_P` below the surface convection zone at the start of the step | cm | `step_driver` (rezone_or_snapshot); `evolve_rotation` (burn_light_elements_over_substep) | `convection_geometry` |  |
| `convective_turnover_timescale_prev` | `convective_turnover_timescale` of the previous model | s | `observables` (refresh_turnover_timescale); `starting_model` (read_starting_model) | `evolve_rotation`, `secular_transport` | Copied after each converged model (and at a model read); the wind interpolates between it and the current value at `substep_time_fraction` |
| `photosphere_logP_prev` | `photosphere_logP` of the previous model, copied with `convective_turnover_timescale_prev` | log10 dyn/cm^2 | `observables` (refresh_turnover_timescale); `starting_model` (read_starting_model) | `evolve_rotation`, `secular_transport` |  |
| `j_rot_iteration_start(json)` | `j_rot` at the start of the structure-rotation iteration | cm^2/s | `step_driver` (converge_with_rotation) | `step_driver` | Restored on each repeat pass when `job%num_rotation_structure_iters` &gt; 1 |
| `xa_iteration_start(n_species_extended,json)` | `xa_start` at the start of the structure-rotation iteration | - | `step_driver` (converge_with_rotation) | `step_driver` |  |
| `envelope_cz_bottom_index_start` | base point of the surface convection zone as the deuterium-burning cap uses it (`nz`: none) | - | `run_driver` (begin_kind_card); `mass_change` (measure_convective_envelope); `step_driver` (rezone_or_snapshot) | `burn_step`, `structure_equations` | A copy of `envelope_cz_bottom_index` taken at the card start and each rezone |

### 3.9 Step working state

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `dlog_temperature_dt(json)` | d log T/`dt` over the last step, for the next step's first guess | 1/s | `rezoning` (regrid_composition_terms); `run_driver` (begin_kind_card); `step_driver` (solve_structure) | `rezoning`, `step_driver` |  |
| `dlog_pressure_dt(json)` | d log P/`dt` over the last step | 1/s | `rezoning` (regrid_composition_terms); `run_driver` (begin_kind_card); `step_driver` (solve_structure) | `rezoning`, `step_driver` |  |
| `dln_luminosity_dt(json)` | d ln L/`dt` over the last step, 2 (L - `L_start`)/(L + `L_start`)/`dt` | 1/s | `rezoning` (regrid_composition_terms); `run_driver` (begin_kind_card); `step_driver` (solve_structure) | `rezoning`, `step_driver` |  |
| `dlog_radius_dt(json)` | d log r/`dt` over the last step | 1/s | `rezoning` (regrid_composition_terms); `run_driver` (begin_kind_card); `step_driver` (solve_structure) | `rezoning`, `step_driver` |  |
| `disk_gate_age_gyr` | age of the last converged model, for the disk-locking tests | Gyr | `step_driver` (converge_with_rotation) | 8 modules (core/ 4, rotation/ 4) |  |
| `convective_turnover_timescale` | convective turnover time one pressure scale height above the base of the surface convection zone | s | `turnover_timescale` (compute_turnover_timescale, time_bottom_cell_crossing, time_scale_height_above_base, turnover_from_interior_new); `run_driver` (begin_kind_card); `mass_change` (measure_convective_envelope) | 9 modules (core/ 6, rotation/ 2, io/ 1) | From the stitched model, every converged model; zero at a card start and without a surface convection zone |
| `photosphere_logP` | log10 photospheric pressure of the last atmosphere integration that reached the photosphere | log10 dyn/cm^2 | `envelope_integration` (atm_get) | `observables`, `starting_model`, `history_output`, `evolve_rotation`, `secular_transport` | Set inside `atm_get`, i.e. by any envelope integration (the triangle's too) |
| `substep_time_fraction` | fraction of the step elapsed at the current rotation sub-step | - | `starting_model` (read_starting_model); `evolve_rotation` (evolve_angular_momentum, spin_down_without_transport) | `evolve_rotation`, `secular_transport` | 0.5 after a model read |
| `omega_substep(json)` | `omega` of the latest rotation sub-step's model, the structure interpolated to the sub-step's end | rad/s | `evolve_rotation` (brake_surface_convection_zone); `mid_timestep_model` (rotate_substep_model); `secular_transport` (check_iterate, restart_iteration, spin_down_uncoupled_surface) | `evolve_rotation`, `mid_timestep_model`, `secular_transport` | `evolve_rotation` copies it into `omega` at the end of the step; it outlives the step: the next step's first solid-body sub-step starts its Newton iteration from its mean |
| `rate_li6(json)` | Li6 proton-capture rate rho X lambda at the end of the step, on the converged model | 1/s | `burn_step` (store_li_be_rates) | `burn_step`, `evolve_rotation` | With `ctrl%use_extended_composition` |
| `rate_li7(json)` | Li7 proton-capture rate at the end of the step | 1/s | `burn_step` (store_li_be_rates) | `burn_step`, `evolve_rotation` |  |
| `rate_be9(json)` | Be9 proton-capture rate at the end of the step | 1/s | `burn_step` (store_li_be_rates) | `burn_step`, `evolve_rotation` |  |
| `deuterium_burning_rate(json)` | deuterium burning rate of each shell above the nuclear cutoff | 1/s | `burn_step` (compute_burn_rates) | `burn_step`, `mid_timestep_model` | Zero while the surface D is below 1e-11 |
| `envelope_overshoot_depth_end` | envelope-overshoot depth `alpha_ov` `H_P` at the end of the step | cm | `step_driver` (burn_light_elements); `evolve_rotation` (burn_light_elements_over_substep) | `convection_geometry`, `evolve_rotation` |  |
| `cz_base_radius_prev` | radius of the base of the surface mixed region found by the last call within the step | cm | `convection_geometry` (locate_surface_mixing_base); `step_driver` (rezone_or_snapshot) | `convection_geometry` | With `ctrl%use_extended_composition`; zeroed at each step start to mean none yet |
| `accreted_mass` | mass accreted during the step | g | `mass_change` (mix_in_accreted_composition) | `burn_step`, `mid_timestep_model` | Set only when the star accretes; its readers check that it does |
| `delta_log_pressure` | shift of log P of the convection zone for the entropy of accreted matter | log10 | `mass_change` (adjust_cz_for_accreted_entropy) | `mass_change`, `structure_equations` | Zero at present: `change_mass` zeroes `envelope_specific_entropy` first |
| `delta_log_temperature` | shift of log T for the same | log10 | `mass_change` (adjust_cz_for_accreted_entropy) | `mass_change`, `structure_equations` |  |
| `cp(json)` | specific heat at constant pressure, for rotation (the viscosity, the transport) | erg/(g K) | `shell_physics` (compute_shell_physics) | `rezoning`, `shell_physics`, `mid_timestep_model` |  |
| `mu(json)` | mean molecular weight, for rotation | - | `shell_physics` (compute_shell_physics) | `rezoning`, `evolve_rotation`, `mid_timestep_model` |  |
| `dlnrho_dlnt(json)` | d ln rho/d ln T at constant P, for rotation | - | `shell_physics` (compute_shell_physics) | `rezoning`, `mid_timestep_model` |  |
| `thermal_diffusivity(json)` | radiative thermal diffusivity K, from viscos, for rotational mixing | cm^2/s | `shell_physics` (compute_shell_physics) | `rezoning`, `mid_timestep_model` |  |
| `kinematic_viscosity(json)` | kinematic viscosity (radiative plus ionic), from viscos, for rotational mixing | cm^2/s | `shell_physics` (compute_shell_physics) | `rezoning`, `mid_timestep_model` |  |
| `env_points` (type `envelope_structure_state`) | the envelope of the latest envelope integration, point by point: log P, log T, log `m`, log rho, log r, X, Z, the gradients, the opacity and its derivatives | cgs, logs base 10 | `envelope_integration` (run_kernel); `envelope_rebuild` (shift_envelope_points); `starting_model` (append_envelope_as_shells) | `stitched_model`, `envelope_integration`, `envelope_rebuild`, `envelope_refit`, `starting_model` | Written by the integration kernel through `atm_get` |
| `atm_points` (type `atmosphere_structure_state`) | the atmosphere of the latest envelope integration, point by point | cgs, logs base 10 | `envelope_integration` (run_kernel) | `stitched_model` |  |
| `atm_point` (type `atmosphere_point_state`) | the atmosphere integration's current point: log tau, log P, log T, log rho, the opacity, the ionization fractions | logs base 10 | `envelope_integration` (run_kernel); `mass_change` (find_accreted_surface_state) | `envelope_integration`, `mass_change` | The mass change falls back on its last values when the Allard lookup fails |

### 3.10 Newton

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `elim_coeff(4,2,json)` | block-elimination coefficients c(k,j,i) of `henyey_elimination` (its header gives the relations) | - | `henyey_elimination` (henyey_eliminate); `structure_equations` (eliminate_shell_pair) | `henyey_elimination`, `henyey_iteration` |  |
| `elim_rhs(4,json)` | right-hand sides e(k,i): the residuals, then (after `henyey_solve`) the corrections dlog P, dlog T, dlog r, dL | log10, Lsun | `henyey_elimination` (henyey_eliminate); `henyey_iteration` (henyey_iterate); `structure_equations` (eliminate_shell_pair) | `henyey_elimination`, `henyey_iteration` |  |
| `log_pressure_delta(json)` | change of log P since the start of the step | log10 | `step_driver` (solve_structure); `henyey_iteration` (apply_henyey_corrections) | `step_driver`, `henyey_iteration`, `structure_equations` | Zeroed (or set to the first-guess extrapolation) before the solve, then accumulated by every correction; feeds `eps_grav` |
| `log_temperature_delta(json)` | change of log T since the start of the step | log10 | `step_driver` (solve_structure); `henyey_iteration` (apply_henyey_corrections) | `step_driver`, `henyey_iteration`, `structure_equations` | As `log_pressure_delta` |
| `max_residual(4)` | per equation (P, T, r, L): the largest \|residual\|, then the largest \|correction\|, then the signed correction there | log10, Lsun | `henyey_elimination` (henyey_eliminate); `henyey_iteration` (find_largest_corrections, linearize_structure_equations, test_henyey_convergence) | `henyey_elimination`, `henyey_iteration` | Misnamed: holds residuals or corrections depending on the stage of the iteration |
| `newton_iterations` | Newton iterations of the step's solve, counted over its levels up to the last convergence, for the run-log model line | - | `henyey_iteration` (accept_small_residuals, set_henyey_correction_factor) | `run_log` |  |

### 3.11 Diagnostics written during the solve

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `eps_grav(json)` | `eps_grav` = -T ds/`dt` of each zone | erg/(g s) | `step_driver` (solve_structure); `structure_equations` (add_shell_energy_budget) | `structure_equations` |  |
| `eps_extra(json)` | the extra energy sources' sum at each zone (`core/structure/energy_sources.f90`) | erg/(g s) | `structure_equations` (add_shell_energy_budget) | `stitched_model` | Zero without a registered source |
| `eps_total(json)` | nuclear energy generation, the sum of the chains' eps: pp1 + pp2 + pp3 + CNO + triple alpha + the alpha captures | erg/(g s) | `structure_equations` (store_zone_energy_channels) | `stitched_model`, `rezoning`, `structure_equations`, `pulsation_output`, `mid_timestep_model` | Zero below the burning cutoff; the neutrino losses are not in it |
| `eps_channels(7,json)` | per channel `i_eps_*`: rows 1-5 the fractions of `eps_total` of pp1, pp2, pp3, CNO and triple alpha, row 6 `eps_nu`, row 7 `eps_grav` | rows 1-5: -; rows 6-7: erg/(g s) | `structure_equations` (store_zone_energy_channels) | `stitched_model`, `rezoning`, `mid_timestep_model` | The alpha captures have no row, so rows 1-5 sum to less than 1 where they burn |
| `luminosity_breakdown_lsun(n_lum_channels)` | luminosity per channel `i_lum_*` (pp1, pp2, pp3, CNO, 3 alpha, neutrinos, gravothermal, alpha captures, the extra energy sources), scaled to sum to L(nz) | Lsun | `observables`, `henyey_iteration`, `structure_equations`, `model_files` (8 procedures) | 8 modules (core/ 6, io/ 2) | Summed every iteration, rescaled again after convergence (`renormalize_luminosity_breakdown`); `compute_timestep` reads it; read from YREC7/MODEL2 files, zeroed by the .mod reader |
| `beta(json)` | gas-pressure fraction `P_gas`/P | - | `structure_equations` (store_zone_thermodynamics) | `stitched_model` |  |
| `eta(json)` | electron degeneracy parameter | - | `structure_equations` (store_zone_thermodynamics) | `stitched_model` | 0 unless the EOS treats the gas as fully ionized |
| `opacity_zone(json)` | Rosseland mean opacity kappa | cm^2/g | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `observables`, `stitched_model`, `rezoning`, `shell_physics`, `pulsation_output`, `mid_timestep_model` | Also written by `compute_shell_physics` (the rezone, rotation, the card start); the last writer wins |
| `gradr(json)` | radiative temperature gradient `grad_rad` | - | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `observables`, `stitched_model`, `rezoning`, `step_driver`, `evolve_rotation`, `mid_timestep_model` | As `opacity_zone` |
| `gradT(json)` | actual temperature gradient d ln T/d ln P | - | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `observables`, `stitched_model`, `rezoning`, `composition_step`, `pulsation_output`, `mid_timestep_model` | As `opacity_zone` |
| `grada(json)` | adiabatic gradient `grad_ad` | - | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `observables`, `stitched_model`, `rezoning`, `step_driver`, `pulsation_output`, `evolve_rotation`, `mid_timestep_model` | As `opacity_zone` |
| `conv_vel(json)` | convective velocity from mixing-length theory | cm/s | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `stitched_model`, `rezoning`, `mid_timestep_model` | As `opacity_zone` |
| `adiabatic_index_gamma1(json)` | first adiabatic exponent `Gamma_1` | - | `structure_equations` (store_zone_thermodynamics) | `stitched_model`, `pulsation_output` |  |
| `neutrino_flux_total(n_nu_fluxes)` | solar-neutrino fluxes at the Earth by slot `i_nu_*` (slots 9-10: the He3+He3 and He3+He4 rates in the same units) | 10^10/(cm^2 s) | `neutrino_flux_output` (write_neutrino_flux_table); `observables` (compute_snu_rates); `structure_equations` (add_shell_energy_budget, start_energy_bookkeeping) | `neutrino_flux_output`, `observables`, `monte_carlo`, `structure_equations`, `history_output` | Only with `ctrl%calc_neutrinos`; also overwritten by the card-start neutrino table |
| `pulse_deps_dlnrho(json)` | d eps/d ln rho, absolute (not logarithmic) | erg/(g s) | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` | Zero below the burning cutoff |
| `pulse_deps_dlnt(json)` | d eps/d ln T, absolute | erg/(g s) | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` | Zero below the burning cutoff |
| `pulse_dlnkap_dlnrho(json)` | d ln kappa/d ln rho | - | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` |  |
| `pulse_dlnkap_dlnt(json)` | d ln kappa/d ln T | - | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` |  |
| `pulse_specific_heat(json)` | specific heat at constant pressure `c_P` | erg/(g K) | `structure_equations` (store_pulse_derivatives) | `stitched_model` |  |
| `pulse_mean_molecular_weight(json)` | mean molecular weight `mu` = R/(specific gas constant) | - | `structure_equations` (store_pulse_derivatives) | `stitched_model` |  |
| `pulse_dlnrho_dlnt(json)` | d ln rho/d ln T at constant P | - | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` |  |
| `pulse_electron_mean_weight_inverse(json)` | 1/`mu_e` | - | `structure_equations` (store_pulse_derivatives) | `stitched_model` |  |

### 3.12 Diagnostics written after convergence

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `stitched` (type `stitched_model_record`) | the stitched model: the converged star (interior, envelope and atmosphere) on one extended grid, with its profile and pulse columns | cgs | `stitched_model` (build_stitched_model, compute_seismic_columns, fill_stitched_arrays) | `observables`, `stitched_model`, `turnover_timescale`, `yrec_output` | Allocated, and zeroed, by the first build |
| `mechanism_velocity(n_transport_mechanisms, json)` | each rotational transport mechanism's velocity at interface k, slots `i_mech_*` of `transport_mechanism_layout` | cm/s | `secular_transport` (compute_velocities) | `stitched_model`, `secular_transport` | From the last iteration of the step's last rotation sub-step |
| `am_diffusion_coeff(json)` | angular-momentum diffusion coefficient `D_omega` on the model grid | cm^2/s | `secular_transport` (set_diffusion_coefficients) | `stitched_model` | The step's last sub-step wins; 0 without rotation |
| `mixing_diffusion_coeff(json)` | rotational mixing coefficient `D_mix` on the model grid | cm^2/s | `secular_transport` (set_diffusion_coefficients) | `stitched_model` | As `am_diffusion_coeff` |
| `log_R_surface_rsun` | log10 of the photospheric radius from L = 4 pi R^2 sigma Teff^4 | log10 Rsun | `observables` (compute_surface_globals) | `observables`, `history_output`, `run_log`, `yrec_output` |  |
| `log_g_surface` | log10 surface gravity G M/R^2 | log10 cm/s^2 | `observables` (compute_surface_globals) | `observables`, `stop_conditions`, `history_output`, `yrec_output` |  |
| `nu_max_uhz` | frequency of maximum oscillation power, scaled from the Sun's by g and Teff | uHz | `observables` (compute_seismic_observables) | `stop_conditions`, `history_output` |  |
| `delta_nu_rho_uhz` | large frequency separation scaled from the mean density | uHz | `observables` (compute_seismic_observables) | `history_output` |  |
| `delta_nu_uhz` | asymptotic p-mode large separation, from the sound-travel time | uHz | `observables` (compute_seismic_observables) | `history_output` |  |
| `delta_Pg` | asymptotic l = 1 g-mode period spacing of the inner cavity, from the Brunt-Vaisala integral | s | `observables` (compute_seismic_observables) | `history_output` | 0 without a cavity |
| `total_moment_of_inertia` | moment of inertia of the star | g cm^2 | `observables` (compute_moment_of_inertia) | `observables`, `history_output` |  |
| `cz_moment_of_inertia` | moment of inertia of the surface convection zone | g cm^2 | `observables` (compute_rotation_observables) | `observables`, `history_output` |  |
| `rotation_period_days` | surface rotation period 2 pi/`omega(nz)`, capped at 9999 | d | `observables` (compute_rotation_observables) | `history_output` | 0 without rotation |
| `surf_velocity_kms` | surface equatorial velocity `omega(nz)` R | km/s | `observables` (compute_rotation_observables) | `history_output` |  |
| `h_shell_bot_mass_msun` | mass coordinate of the bottom of the H-burning shell | Msun | `observables` (compute_h_shell_boundaries) | `history_output` | 0 without a shell, as the five below |
| `h_shell_mid_mass_msun` | mass coordinate of the middle of the H-burning shell | Msun | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `h_shell_top_mass_msun` | mass coordinate of the top of the H-burning shell | Msun | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `h_shell_bot_radius` | radius of the bottom of the H-burning shell | r/R | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `h_shell_mid_radius` | radius of the middle of the H-burning shell | r/R | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `h_shell_top_radius` | radius of the top of the H-burning shell | r/R | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `central_log10_temperature` | central temperature, extrapolated from the first point | log10 K | `observables` (compute_central_conditions) | `observables`, `monte_carlo`, `history_output` |  |
| `central_log10_pressure` | central pressure | log10 dyn/cm^2 | `observables` (compute_central_conditions) | `observables`, `monte_carlo`, `history_output` |  |
| `central_log10_density` | central density, from the EOS | log10 g/cm^3 | `observables` (compute_central_conditions) | `observables`, `monte_carlo`, `history_output` |  |
| `central_beta` | central `P_gas`/P | - | `observables` (compute_central_conditions) | `history_output` |  |
| `central_degeneracy_eta` | central electron degeneracy parameter | - | `observables` (compute_central_conditions) | `history_output` |  |
| `envelope_mass_msun` | mass of the surface convection zone | Msun | `observables` (locate_surface_cz_base) | `monte_carlo`, `history_output` |  |
| `envelope_radius_rsun` | radius of the base of the surface convection zone | Rsun | `observables` (locate_surface_cz_base, refresh_turnover_timescale); `starting_model` (read_starting_model) | `monte_carlo`, `history_output` | Written twice per model: by `locate_surface_cz_base`, then by the turnover-time walk |
| `core_cz_mass_msun` | mass of the convective core | Msun | `observables` (locate_core_cz) | `history_output` | 0 without one |
| `envelope_cz_temperature` | temperature at the base of the surface convection zone, linear | K | `observables` (locate_surface_cz_base) | `history_output` |  |
| `envelope_cz_density` | density at the base of the surface convection zone | g/cm^3 | `observables` (locate_surface_cz_base) | `history_output` |  |
| `envelope_cz_pressure` | pressure at the base of the surface convection zone | dyn/cm^2 | `observables` (locate_surface_cz_base) | `history_output` |  |
| `envelope_cz_opacity` | opacity at the base of the surface convection zone | cm^2/g | `observables` (locate_surface_cz_base) | `history_output` |  |
| `cl37_snu_rate` | chlorine (37Cl) detector capture rate | SNU | `observables` (compute_snu_rates) | `observables`, `monte_carlo`, `history_output` | 0 without `ctrl%calc_neutrinos` |
| `ga71_snu_rate` | gallium (71Ga) detector capture rate | SNU | `observables` (compute_snu_rates) | `observables`, `monte_carlo`, `history_output` | 0 without `ctrl%calc_neutrinos` |

### 3.13 Private to one module

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `log_l_lsun_prev_check` | log(L/Lsun) at the previous solar-calibration check, only printed | log10 Lsun | `calibration` (check_solar_calibration) | `calibration` |  |
| `log_r_rsun_prev_check` | log(R/Rsun) at the previous solar-calibration check, only printed | log10 Rsun | `calibration` (check_solar_calibration) | `calibration` |  |
| `delta_x` | the solar calibration's Newton correction to X | - | `calibration` (check_solar_calibration) | `calibration` |  |
| `delta_alpha` | the solar calibration's Newton correction to the mixing length | - | `calibration` (check_solar_calibration) | `calibration` |  |
| `age_gyr_at_target_radius` | age interpolated at the star calibration's target radius | Gyr | `calibration` (check_star_calibration) | `calibration` |  |
| `age_gyr_prev_model` | age of the previous model, for the star calibration | Gyr | `calibration` (check_star_calibration) | `calibration` |  |
| `log_l_lsun_at_target_radius` | log(L/Lsun) interpolated at the star calibration's target radius | log10 Lsun | `calibration` (check_star_calibration) | `calibration` |  |
| `log_l_lsun_at_target_radius_prev_run` | log(L/Lsun) at the target radius in the previous star-calibration run | log10 Lsun | `calibration` (check_star_calibration) | `calibration` |  |
| `log_l_lsun_prev_model` | log(L/Lsun) of the previous model, for the star calibration | log10 Lsun | `calibration` (check_star_calibration) | `calibration` |  |
| `r_rsun_prev_model` | R/Rsun of the previous model, for the star calibration | Rsun | `calibration` (check_star_calibration, setup_star_calibration) | `calibration` |  |
| `envelope_cz_base_zone_prev` | zone of the base of the surface mixed region found by the last call within the step (`cz_base_radius_prev`), the starting zone of the next call | - | `convection_geometry` (locate_surface_mixing_base) | `convection_geometry` |  |
| `central_saha_state` | the Saha state the central EOS call carries from one model to the next | - | `observables` (compute_central_conditions) | `observables` |  |
| `max_correction_index(4)` | point of the largest correction of each variable, for the convergence log | - | `henyey_iteration` (find_largest_corrections) | `henyey_iteration` |  |
| `surface_bc(6)` | outer boundary-condition rows from the envelope fit, in `henyey_solve`'s form | - | `henyey_iteration` (form_surface_bc_rows) | `henyey_iteration` |  |
<!-- END GENERATED members -->

## 4. Members removed

Removed in 2026-10, with what became of each. `tools/gen_star_docs.py --check` fails when the first column names a live member, so a name a later rename took again is described, not named (the d ln ρ/d ln T row).

| removed | why | what became of it |
|---|---|---|
| `old_shell_mass`, `converged_zone` | never referenced | deleted |
| `old_moment_of_inertia`, `log10_solar_mass` | written, never read | deleted |
| `fxion_zone`, `scp`, `he3_he3_luminosity_zone`, `he3_burning_luminosity_zone`, `pulse_dlnrho_dlnp` | diagnostics no writer emitted (`scp` was a copy of `pulse_specific_heat`) | deleted |
| `valfmlt`, `vphmlt`, `vcmxmlt` | always zero, in the F77 code too: their source scalars were only ever zeroed, for a pulsation format since retired | deleted |
| `qiw` | `omega_from_j`'s dI/dω, recomputed by every user | a local of each user |
| `chosen_dt`, `ikut_flag`, `convective_velocity` | outputs no caller used | locals of `compute_timestep` and `read_starting_model` |
| `mixed_zone_bounds_no_overshoot`, `radiative_zone_bounds` | outputs no caller used | locals of the convection-zone search |
| `initial_composition_code` | read with a Model2 file, never used | a local of the reader |
| `envelope_cz_log_radius` | used by one procedure | a local of `locate_surface_cz_base` |
| `dlnrho_dlnp` and the first d ln ρ/d ln T member | reached the next coefficient pass only as an EOS starting value every EOS family overwrites (`docs/bug_sweep.md` O23) | locals of `compute_henyey_coefficients`; the member now named `dlnrho_dlnt` is the former `qdt`, `compute_shell_physics`'s value for rotation |
| `solar_wind_mass_loss_rate_msun_yr`, `wind_reference_omega`, `wind_max_omega`, `use_rotation_scaled_solar_wind` | read only by the rotation-scaled solar wind, commented out on 2026-10-12 | deleted; `star_info_lib.f90`'s Sources section says how to revive them |
| `stotal` | a copy of `log_total_mass`, kept equal by the same writers | merged into `log_total_mass` |
| `solar_bolometric_magnitude` | a constant, 4.79 | a parameter of `physical_constants` |
| `hsra_log_t_offset` | the HSRA T(τ) fit at τ = 2/3, a constant | the envelope kernel's configuration takes atm's pure `hsra_t_tau_offset()` |
| `be7_mass_fraction_zone`, `neutrino_flux_zone` | used by one procedure, the neutrino table's | locals of `write_neutrino_flux_table` |
| `zahn_coupling` | rotation's working value within a sub-step | a field of rotation's `circulation_velocity_state` (`circ_scr`) |
| `envelope_specific_entropy` | passed from `change_mass` to `apply_mass_change` | a local of `change_mass`, passed as an argument |

The generator keeps the count of never-referenced members at zero.

## 5. Names

Renamed in 2026-10 (the old name in parentheses), because each invited a wrong reading: `eps_grav` (`gravitational_luminosity`; eps_grav in erg/(g s), not a luminosity); `dlog_temperature_dt`, `dlog_pressure_dt`, `dln_luminosity_dt`, `dlog_radius_dt` (the four `*_entropy_term`; rates of change for the next step's first guess); `envelope_overshoot_depth_start`, `_end` (`pressure_scale_height_start`, `_end`; alpha_ov H_P in cm); `triangle_log_teff`, `triangle_log_l_lsun`, `triangle_orientation` (`trial_log_temperature`, `trial_log_luminosity`, `trial_sign_flag`); `stored_vertex_index` (`istore_flag`); `r_rsun_prev_model` (`log_r_prev_model`; linear R/R☉); `solar_calibration_converged` (`solar_calibration_active`); `age_gyr` (`dage`). `star_info_lib.f90`'s Sources section lists the other renames and the members' F77 homes.

Still to watch:

- Earlier values: `<member>_start` is the start of the step, `<member>_prev` the previous model (`cz_base_radius_prev` and `envelope_cz_base_zone_prev`: the previous call within the step; `log_l_lsun_prev_check`, `log_r_rsun_prev_check`: the previous calibration check), `<member>_iteration_start` the start of the structure-rotation iteration.
- `max_residual` holds residuals or corrections, depending on the stage of the Newton iteration.
- `ln10_over_lsun` is ln(10)/L☉, not a logarithm.
- `star_mass` is in M☉ without a unit suffix: MESA's name.

## 6. Checks

state/ holds types and constants, so its checks are about how the rest of the code uses them:
- `test_generated_files.py` runs `tools/gen_star_docs.py --check`: it fails when sections 2 and 3 are stale, when a member has no comment of the form `<meaning> [<unit>]; <writers>. <notes>`, is declared with another member on one line, names a writer that does not write it, or is never written or never read outside the tests; when the first column of section 4's table names a live member; and when a member breaks a class rule: an earlier value (`<member>_start`, `_prev`, `_iteration_start`) is in the Step snapshot class and that class holds nothing else, and any other member that one module alone writes and reads is in the class "Private to one module" and that class holds nothing else. Each of these was confirmed to fail on a mutated copy of `star_info_lib.f90` or of this README.
- `test_boundaries.py` runs `tools/check_boundaries.py` (the physics domains in `STAR_FREE_DIRS` never use `star_info_lib`; `STAR_ACCESS_BUDGET` keeps every directory at zero references to a module-level star) and `tools/check_star_contracts.py` (the declared `star reads:` / `star writes:` lists match the code; every procedure of `core/run/` and `core/step/` that takes the star declares them, and procedures are told apart by module and host, so two modules' `wind_now` are each checked against their own code).
- `test_reentry.py`: a second `run_yrec` call in one process writes the same files as a fresh process, i.e. the template star of `yrec_reset` is a clean start.
- The decks: moving a member changes the layout of `star_info`, never a result; the reordering by lifetime class was byte-identical over the full deck set.

## 7. How to extend

**Add a star field** (traced on 2026-10-29 with a throwaway diagnostic written by `observables` and read by the history writer):
1. Choose its lifetime class (section 2) and declare it in that class's block of `type star_info` in `state/star_info_lib.f90`, one member per declaration, with its comment directly above: `! <meaning> [<unit>]; <the modules that write it>. <notes>`. Two classes follow from the code, and `gen_star_docs.py` says so when a member is elsewhere: an earlier value, named `<member>_start` (the start of the step), `<member>_prev` (the previous model) or `<member>_iteration_start`, goes in the Step snapshot class; a member that one module alone writes and reads goes in "Private to one module" (or is better a local, or a field of that domain's scratch record). If it must start a run at something other than zero, `.false.` or blank, give the declaration an initializer: every run starts from the static template star of `core/run/yrec_reset.f90`. If every run card must start it afresh (a latch the run switches off, as `settling_active`), arm it from its control in `begin_kind_card` (`core/run/run_driver.f90`), which runs at the start of every card and so of every Monte-Carlo run. A value that never changes during a run is a named parameter (`docs/GUIDELINES.md`, Code style rule 3), not a member; a value one procedure computes and uses is a local, and a domain's working arrays are a scratch record its driver passes down (star rule 5).
2. Write and read it as `s%<name>` in star-layer code (`core/`, `io/`, `setup/`, rotation's drivers), with the code rules `tools/check_boundaries.py` checks (a real power is `exp10()` or `pow()`, and `exp`, `log` and the like need `use math_lib` in scope). The physics domains (`STAR_FREE_DIRS` in `tools/check_boundaries.py`) never see the star: pass the value to them as an argument, or put a configuration value in the domain's settings record (`setup/README.md`).
3. If a procedure that writes or reads it declares its star contract (the `! star reads:` / `! star writes:` lines above its `subroutine` line), run `python3 tools/check_star_contracts.py --fix` in `src/`; `test_boundaries.py` checks the contracts.
4. Run `python3 tools/gen_star_docs.py` in `src/` to add its row here; it refuses a member without a comment in the form above, one whose comment names a module that does not write it, or one that nothing outside the tests reads.
5. If the field must survive a restart or a reload between run cards, write and read it in `io/model_files.f90` (`io/README.md` §5, "Add a field to the `.mod` model": one named line per global value; readers ignore names they do not know). If it is output, add a history or profile column (`io/README.md` §5, "Add a history column" or "Add a profile column"; a profile column's table row and values are in `core/diagnostics/stitched_model.f90`).
6. `make -j8` in `src/` rebuilds everything that uses `star_info_lib` (the dependencies are already there unless a file gains a `use`; then `python3 tools/gen_makefile_deps.py`). Verify with the unit tests and the decks; moving a member does not change any result.

**Add a control:** a control is not a star field: it is a member of `controls_buffer` (`io/controls_buffer.f90`) with its default as the initializer, and `tools/gen_controls_state.py` generates its component of `star%ctrl` (or of `star%job`, when its registry row's `home` is `job`) into `controls_state_def.inc` or `job_controls_def.inc`. Nothing in `star_info_lib.f90` is edited by hand. The recipe is `io/README.md` §7, "Add a control": the declaration, its row in `defaults/controls_registry.tsv` (`canonical`, `home`, `used`), the three generators, its copy into a domain's settings record (`setup/README.md`), and `YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py`.

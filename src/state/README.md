# state/ — the model state: `star` and the constants it is built from

This directory holds YREC's global state. `star_info_lib` defines the type `star_info`; it declares no instance (since 2026-10): `run_yrec` owns the run's star and hands it to every star-layer routine as the argument `s` (`docs/GUIDELINES.md`, "The star object: six rules"); it also defines the two input bundles nested in it, `star%ctrl` (type `controls_state`) and `star%job` (type `star_job`). `composition_layout` and `neutrino_flux_layout` hold the slot constants of the composition and neutrino-flux vectors, in modules of their own so that the physics domains can use them without `star`. `physical_constants` holds the mathematical and physical constants. The components of `controls_state`, the card arrays of `star_job` and the copy module `controls_sync` are generated from `io/controls_buffer.f90` by `tools/gen_controls_state.py`; this README does not list them one by one. The controls are documented in `src/defaults/controls_registry.md` and, with their defaults, in `src/defaults/physics.defaults` and `src/defaults/control.defaults`.

| Module (file) | What it holds |
|---|---|
| `star_info_lib` | The types `star_info`, `star_job` and `controls_state`; the size `json` = 5000 of every per-zone array and the other capacities (`max_convective_zones`, `max_radiative_zones`, `max_mc_runs`); the index constants of the luminosity channels (`i_lum_*`), the energy channels (`i_eps_*`) and the gradient rows; and, outside `star`, the namelist-path overrides of `pyyrec` and the reset flags of `yrec_reset`. |
| `composition_layout` | The rows of `xa`: `i_h1` ... `i_be9` (row `i_metals` is the total Z), and the 12 surface-mixture slots `ix_*`. |
| `neutrino_flux_layout` | The slots `i_nu_*` of the neutrino-flux vectors. |
| `physical_constants` | `ln10`, `cc13`, `cpi`, `seconds_per_year` and the other compile-time constants; `csig`, `cgl` (log G), `gas_constant` and the other constants that `core/run/run_setup.f90` sets once per run. |
| `integration_parameters` | The Bulirsch-Stoer integrator's controls (`tolerance_fraction`, `max_stage_index`, `extrap_order`), set by the controls read. |
| `controls_state_def.inc`, `job_controls_def.inc`, `controls_sync` | Generated (see above): the components of `star%ctrl`, the run-card arrays of `star%job`, and the copies between `controls_buffer` and `star`. |

The rest of this README describes the members of `star`: what each one is, its unit, who writes and reads it, and when its value can be trusted.

---

## 1. Conventions

- **Points.** A per-zone array has `json` = 5000 entries, of which `1:star%nz` are the model. Point 1 is the innermost point (not the centre itself), point `nz` the fitting point where the envelope integration takes over; the atmosphere and envelope above it are not in `star` (`core/diagnostics/stitched_model.f90` builds them for output).
- **Units.** CGS unless the name carries its unit (`_lsun`, `_msun`, `_gyr`, `_yr`, `_days`, `_kms`); `log` means log₁₀. The structure variables are log₁₀ P, log₁₀ T, log₁₀ r and L in L☉ against log₁₀ m. A few old names break the rule (`age_gyr` is in Gyr, `star_mass` in M☉); the tables give every unit, with `-` for a dimensionless quantity and `log10 g` for log₁₀ of a quantity in grams.
- **Lifetime.** Most members have no explicit lifetime: the same storage is rewritten at different stages of a step. Section 2 sorts them into classes; the tables give the exceptions.
- **A second run in one process.** `yrec_reset` copies the whole of `star` at the first `run_yrec` call and restores that copy at every later one, so every member starts a run at its declaration default (zero, `.false.` or blank where none is given).

## 2. Lifetime classes

| Class | Members | What it is | When its value is valid |
|---|---|---|---|
| Model | 30 | The converged model: the mesh, the structure, the composition, the global scalars, and the envelope fit that is its outer boundary condition. Persists from step to step; the core of it is written to the `.mod` file (`io/model_files.f90`, marked "in .mod"). | Between steps, and inside a step for the model being solved. The structure is overwritten in place by the Newton corrections, so during a step it holds the current iterate. |
| Rotation state | 12 | The rotation curve and the shape of the rotating star. The specific angular momentum `j_rot` is conserved; omega, the moments of inertia and the shape factors are derived from it and the structure. | Between steps. Meaningful only with `job%rotation_active`; without rotation the shape factors are 1. |
| Regions of the current model | 10 | Point indices and bounds of the convective core, the surface convection zone, the hydrogen-burning shell and the mixed regions. | Located again at each stage that needs them (composition step, rezone, after convergence); valid for the model they were located on. |
| Step snapshot | 24 | Copies of model quantities at the start of the step (`*_start`, `old_*`, `*_old`, `pphot0`) or of the structure-rotation iteration (`orig_*`), mapped onto the step's mesh. | From the start-of-step copy (the rezone, or `rezone_or_snapshot` without one) to the end of the step. |
| Step working state | 28 | Values one stage computes for a later stage of the same or the next step, and that are not output: rates of change, burning rates, the depth of the surface mixed region, the turnover time, the shell physics that rotation uses, the outer layers of the latest envelope integration, the last rotation sub-step's omega. | From the stage that writes them; each table row says which. |
| Newton | 10 | Arrays of the Henyey solve, rewritten on every iteration. | During the solve. `log_pressure_delta` and `log_temperature_delta`, the changes over the step, stay valid until the next step's solve starts. |
| Diagnostics written during the solve | 21 | Per-zone physics stored for output by `compute_henyey_coefficients` on every Newton iteration (energy channels, EOS and opacity values, gradients, pulsation derivatives). | After convergence they hold the last iteration's values: the structure that iteration linearized, before its correction when one was applied. Some are also rewritten by `compute_shell_physics`; the last writer wins. |
| Diagnostics written after convergence | 36 | Observables computed once per model by `compute_observables` (after the stitched model is built), and the rotational transport's coefficients and velocities kept for the profiles. | From `compute_observables` (or the rotation step) until the next model; zero at a run start. |
| Step control | 8 | The timestep and the flags that steer one step's solve. | During and between steps. |
| Run and card bookkeeping | 29 | Per-run and per-card state: the calibration loops, the card's mixing length and mass, the surface mixture of the first card, the atmosphere switch, the reason the run ended. | From the card (or run) that sets them. |
| Constants set at setup | 16 | Values fixed for the whole run during the controls read or the setup before the first card: the solar units, the nuclear rate scales, a few atmosphere and opacity settings. | After `read_controls` and `prepare_run`; a Monte-Carlo run rescales some of them. |
| Scratch | 4 | Values stored in `star` but used only inside the module (often the procedure) that writes them. | Only inside that module. |
| Never referenced | 0 | Declared and not used (the two there were removed in 2026-10). | — |

228 members in all, besides `ctrl` and `job` (section 4).

## 3. The members

How to read the tables:

- **member** gives the shape; `json` is the zone dimension (use `1:nz`).
- **written by** names the procedures in which `star%member` is the target of an assignment or of a `read` statement, or an actual argument to a dummy argument declared `intent(out)`, or `intent(inout)` that the callee (directly or through its own callees) assigns. So a member written inside a callee appears under its caller: `stored_vertex_index` is written by `compute_surface_boundary` through `refresh_surface_bc`'s argument, and the table names `step_driver (solve_level)`. With more than six procedures only the modules are listed.
- **read by** names the modules where the value of `star%member` is used: in an expression, or passed to a dummy argument that the callee reads. Passing it to an `intent(out)` dummy, or to an `intent(inout)` one that the callee never reads, is not a read. With more than seven modules the count is given by directory.
- Module names are file names (`core/step/step_driver.f90` holds `step_driver`), except `controls_reader`, in `io/read_controls.f90`. Tests (`*/test/`) are not counted.
- The two columns were derived from the code at commit e9c00d7f by a script, not by hand. To check one, `grep -rn 'star%member' src --include='*.f90'`.

### 3.1 Model

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `log_mass(json)` | log10 of the mass m enclosed by point i (the mesh coordinate) | log10 g | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `model_files` (11 procedures) | 14 modules (core/ 10, rotation/ 3, io/ 1) | in .mod |
| `luminosity_lsun(json)` | luminosity L(m) through point i | Lsun | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `model_files` (11 procedures) | 15 modules (core/ 11, io/ 2, rotation/ 2) | in .mod; a Newton variable (linear, not log) |
| `logR(json)` | log10 radius | log10 cm | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `henyey_iteration`, `model_files` (12 procedures) | 17 modules (core/ 13, io/ 2, rotation/ 2) | in .mod; a Newton variable |
| `logP(json)` | log10 pressure | log10 dyn/cm^2 | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `henyey_iteration`, `model_files` (12 procedures) | 17 modules (core/ 13, io/ 2, rotation/ 2) | in .mod; a Newton variable |
| `logT(json)` | log10 temperature | log10 K | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `henyey_iteration`, `model_files` (12 procedures) | 18 modules (core/ 14, io/ 2, rotation/ 2) | in .mod; a Newton variable |
| `logRho(json)` | log10 density, from the EOS at (P, T, X, Z) | log10 g/cm^3 | `rezoning`, `starting_model`, `composition_step`, `step_driver`, `structure_equations`, `model_files` (12 procedures) | 17 modules (core/ 13, io/ 2, rotation/ 2) | in .mod; rewritten from the EOS on every Newton iteration (`store_shell_structure`), at the structure that iteration linearized, so it can lag the last correction |
| `convective_flag(json)` | true where the zone is convective by the temperature-gradient test | logical | `rezoning`, `starting_model`, `composition_step`, `step_driver`, `structure_equations`, `model_files`, `evolve_rotation` (16 procedures) | 11 modules (core/ 8, rotation/ 2, io/ 1) | in .mod; rewritten on every Newton iteration and by `compute_shell_physics` |
| `xa(15,json)` | mass fractions, row = species slot `i_*` of `composition_layout`; row `i_metals` is the total Z | - | `rezoning`, `starting_model`, `burn_step`, `composition_step`, `mass_change`, `step_driver`, `model_files`, `evolve_rotation`, `mid_timestep_model`, `secular_transport` (29 procedures) | 21 modules (core/ 14, io/ 4, rotation/ 3) | in .mod; rows 12-15 only with `job%use_extended_composition` |
| `m(json)` | mass enclosed by point i | g | `rezoning` (set_new_shell_masses); `starting_model` (build_shell_masses); `composition_step` (settle); `mass_change` (rescale_cz_masses); `step_driver` (evolve_step); `mid_timestep_model` (advance_rotation_substep) | 11 modules (core/ 7, rotation/ 3, io/ 1) | rebuilt from log_mass after every rezone or read |
| `dm(json)` | mass of the shell around point i: half the span to the neighbouring points; dm(nz) reaches the total mass | g | `rezoning` (set_new_shell_masses); `starting_model` (build_shell_masses); `mass_change` (update_total_mass); `step_driver` (evolve_step) | 13 modules (core/ 10, rotation/ 3) |  |
| `stored_envelope_state(4)` | state (log P, log T, log r, log m/M) where a vertex's envelope integration first reached tenv | log10 cgs | `henyey_iteration` (refresh_surface_bc) | `rezoning` | `rezone` uses it to add a point at the surface |
| `envelope_fit_coeffs(9)` | coefficients c(1:9) of the envelope fit at the fitting point: log R, log L and log Teff each = c log P + c log T + c | - | `henyey_iteration` (refresh_surface_bc); `model_files` (read_yrec7, read_model2) | `henyey_iteration` | rebuilt by `compute_surface_boundary` when the triangle changes; read from YREC7/MODEL2 files, not written to .mod |
| `triangle_log_l(3)` | log10 L of the three triangle vertices | log10 Lsun | `henyey_iteration` (refresh_surface_bc); `model_files` (read_yrec7, read_model2) | `henyey_iteration` | was `trial_log_luminosity` until 2026-10 |
| `triangle_log_teff(3)` | log10 Teff of the three vertices of the envelope-fit triangle | log10 K | `henyey_iteration` (refresh_surface_bc); `model_files` (read_yrec7, read_model2) | `henyey_iteration` | was `trial_log_temperature` until 2026-10 |
| `fit_point_pressure(3)` | log10 P at the fitting point of each triangle vertex's envelope integration | log10 dyn/cm^2 | `henyey_iteration` (refresh_surface_bc); `model_files` (read_yrec7, read_model2) | `rezoning`, `henyey_iteration` |  |
| `fit_point_temperature(3)` | log10 T at the fitting point of each vertex | log10 K | `henyey_iteration` (refresh_surface_bc); `model_files` (read_yrec7, read_model2) | `henyey_iteration` |  |
| `fit_point_radius(3)` | log10 r at the fitting point of each vertex | log10 cm | `henyey_iteration` (refresh_surface_bc); `model_files` (read_yrec7, read_model2) | `rezoning`, `henyey_iteration` |  |
| `nz` | number of points of the interior model (out to the fitting point) | - | `rezoning`, `starting_model`, `step_driver`, `model_files` (9 procedures) | 22 modules (core/ 15, io/ 4, rotation/ 3) | in .mod |
| `model_number` | counter of converged models in the run; the step is aged only when it is >= 0 | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step); `model_files` (read_mod_model, read_yrec7, read_model2) | 9 modules (io/ 6, core/ 3) | in .mod; reset per card from `ctrl%set_initial_model_number` |
| `log_total_mass` | log10 of the total mass M | log10 g | `starting_model` (rescale_and_refit_envelope); `mass_change` (update_total_mass); `model_files` (read_mod_model, read_yrec7, read_model2) | 14 modules (core/ 9, rotation/ 3, io/ 2) | in .mod |
| `star_mass` | total mass M | Msun | `starting_model` (rescale_and_refit_envelope); `mass_change` (update_total_mass); `model_files` (read_mod_model, read_yrec7, read_model2) | 10 modules (core/ 5, io/ 3, rotation/ 2) | in .mod; the name lacks the `_msun` suffix (MESA's history name) |
| `log_Teff` | log10 effective temperature, from the envelope fit (`envelope_fit_coeffs` 7-9) | log10 K | `henyey_iteration` (apply_henyey_corrections); `model_files` (read_mod_model, read_yrec7, read_model2) | 19 modules (core/ 11, io/ 6, rotation/ 2) | in .mod; reset after every Newton correction |
| `log_L` | log10 of the luminosity at the fitting point, L(nz) | log10 Lsun | `henyey_iteration` (refresh_surface_bc, apply_henyey_corrections); `model_files` (read_mod_model, read_yrec7, read_model2) | 15 modules (core/ 9, io/ 4, rotation/ 2) | in .mod; reset after every Newton correction |
| `xnew` | surface X used for the envelope integrations and the surface opacity tables | - | `model_rescaling` (rescale_surface_x); `rezoning` (refresh_surface_opacity); `starting_model` (reset_state_after_read, interpolate_new_fitting_point) | `model_rescaling`, `rezoning`, `starting_model`, `henyey_iteration` | follows xa(i_h1,nz); refreshed when it moves by more than 1e-8 (`refresh_surface_opacity`), or set by a rescale |
| `znew` | surface Z used for the envelope integrations | - | `model_rescaling` (rescale_interior_z); `rezoning` (refresh_surface_opacity); `starting_model` (reset_state_after_read, interpolate_new_fitting_point) | `model_rescaling`, `starting_model`, `henyey_iteration` | updated with xnew |
| `stotal` | copy of log_total_mass | log10 g | `model_rescaling` (rescale_main_sequence_mass, rescale_hb_total_mass); `starting_model` (read_starting_model_file); `mass_change` (update_total_mass) | `observables`, `envelope_integrands`, `envelope_rebuild`, `envelope_refit`, `starting_model` | kept equal to log_total_mass by its writers; read by the envelope integrations and observables |
| `senv` | log10(m(nz)/M): the mass fraction inside the fitting point | log10 | `envelope_rebuild`, `rezoning`, `starting_model`, `mass_change`, `henyey_iteration` (7 procedures) | `envelope_integration`, `envelope_rebuild`, `envelope_refit`, `rezoning`, `starting_model` | recomputed at every solve level (`refresh_surface_bc`) |
| `age_gyr` | age of the model | Gyr | `step_driver` (advance_composition_and_age); `model_files` (read_mod_model, read_yrec7, read_model2) | 9 modules (io/ 6, core/ 3) | was `dage` until 2026-10 |
| `stored_vertex_index` | which triangle vertex stored `stored_envelope_state` (0: none) | - | `run_driver` (begin_kind_card); `step_driver` (reload_model_if_diverged, solve_level) | `step_driver` | was `istore_flag` until 2026-10 |
| `triangle_orientation` | orientation (+1 or -1) of the envelope-fit triangle | - | `run_driver` (begin_kind_card); `step_driver` (reload_model_if_diverged, solve_level) | `step_driver` | was `trial_sign_flag` until 2026-10 |

### 3.2 Rotation state

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `fp_rot(json)` | rotational correction f_P to the hydrostatic equation | - | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration` (7 procedures) | `stitched_model`, `rezoning`, `starting_model`, `henyey_iteration`, `structure_equations`, `evolve_rotation` | 1 without rotation; set from omega by `rotation_shape_factors` |
| `ft_rot(json)` | rotational correction f_T to the radiative-transfer equation | - | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration` (7 procedures) | `stitched_model`, `rezoning`, `starting_model`, `henyey_iteration`, `structure_equations`, `evolve_rotation` | as fp_rot |
| `omega(json)` | angular velocity of each shell | rad/s | `envelope_rebuild`, `rezoning`, `starting_model`, `henyey_iteration`, `model_files`, `evolve_rotation` (14 procedures) | 14 modules (core/ 10, io/ 3, rotation/ 1) | in .mod |
| `i_rot(json)` | moment of inertia of each shell | g cm^2 | `rezoning` (respin_on_new_grid); `starting_model` (initialize_rotation_state); `step_driver` (evolve_step); `henyey_iteration` (update_rotation_after_corrections); `evolve_rotation` (evolve_angular_momentum, spin_down_without_transport) | `observables`, `stitched_model`, `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `evolve_rotation` |  |
| `j_rot(json)` | specific angular momentum of each shell | cm^2/s | `rezoning`, `starting_model`, `mass_change`, `step_driver`, `henyey_iteration`, `evolve_rotation`, `mid_timestep_model`, `secular_transport` (17 procedures) | 8 modules (core/ 5, rotation/ 3) | the conserved quantity; omega follows from it (`omega_from_j`) |
| `kinetic_energy_rot(json)` | rotational kinetic energy of each shell, omega j dm / 2 | erg | `rezoning` (retally_angular_momentum); `starting_model` (initialize_rotation_state); `step_driver` (evolve_step); `henyey_iteration` (update_rotation_after_corrections) | `rezoning`, `starting_model`, `step_driver`, `structure_equations` |  |
| `mean_radius(json)` | r0, the mean radius of the shell's level surface | cm | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `evolve_rotation` (8 procedures) | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration` |  |
| `eta_squared(json)` | Radau's eta_2 of each shell (shape of the level surface) | - | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `evolve_rotation` (8 procedures) | `rezoning`, `starting_model`, `step_driver`, `henyey_iteration`, `evolve_rotation` |  |
| `mean_gravity(json)` | effective gravity averaged over the level surface, <g> | cm/s^2 | `rezoning` (respin_on_new_grid); `starting_model` (initialize_rotation_state); `step_driver` (evolve_step, converge_with_rotation); `henyey_iteration` (update_rotation_after_corrections) | `rezoning`, `mid_timestep_model` |  |
| `am_transport_convective_flag(json)` | true where a zone is convective for angular-momentum transport (overshoot included) | logical | `rezoning` (respin_on_new_grid); `starting_model` (read_starting_model) | `rezoning`, `henyey_iteration` | from `am_convective_regions` at the model read and at every rezone; `omega_from_j` holds those zones in solid-body rotation during the Newton iterations |
| `total_angular_momentum` | total angular momentum J of the star | g cm^2/s | `run_driver` (begin_kind_card); `step_driver` (evolve_step, reload_model_if_diverged, rezone_or_snapshot) | `run_driver`, `step_driver`, `history_output` | zeroed at each card start, then totalled by the model read and every rezone |
| `total_rotational_ke` | total rotational kinetic energy | erg | `run_driver` (begin_kind_card); `step_driver` (evolve_step, reload_model_if_diverged, rezone_or_snapshot) | `run_driver`, `step_driver`, `history_output` | as total_angular_momentum |

### 3.3 Regions of the current model

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `core_cz_top_index` | outermost point of the convective core (1: no convective core) | - | `rezoning` (locate_new_cz_edges); `run_driver` (begin_kind_card); `step_driver` (evolve_step, advance_composition_and_age, burn_light_elements) | `observables`, `rezoning`, `run_driver`, `step_driver` | set by the composition step's predictor pass, the rezone, `burn_light_elements` and, after convergence, `locate_shell_boundaries` |
| `envelope_cz_bottom_index` | innermost point of the surface convection zone (nz: none) | - | `rezoning` (locate_new_cz_edges); `run_driver` (begin_kind_card); `step_driver` (evolve_step, advance_composition_and_age, burn_light_elements) | `observables`, `rezoning`, `starting_model`, `run_driver`, `mass_change`, `step_driver` | as core_cz_top_index |
| `iov1` | first point of the overshoot region below the surface convection zone (-1: none) | - | `rezoning`, `starting_model`, `composition_step`, `step_driver`, `evolve_rotation`, `mid_timestep_model` (7 procedures) | `shell_physics`, `structure_equations` | from `find_convection_zones`; with `ctrl%ladov` the radiative points from iov1 to iov2 take the adiabatic gradient |
| `iov2` | last point of that overshoot region (-1: none) | - | `rezoning`, `starting_model`, `composition_step`, `step_driver`, `evolve_rotation`, `mid_timestep_model` (7 procedures) | `shell_physics`, `structure_equations` | as iov1 |
| `mixed_zone_bounds(12,2)` | first and last point of each mixed (convective plus overshoot) region | - | `composition_step` (advance_composition); `step_driver` (burn_light_elements) | `burn_step`, `composition_step`, `step_driver` | rows up to max_convective_zones; rewritten by every composition-step pass and by `burn_light_elements` |
| `jcz` | base point of the surface convection zone, as the deuterium-burning cap uses it | - | `run_driver` (begin_kind_card); `mass_change` (measure_convective_envelope); `step_driver` (rezone_or_snapshot) | `burn_step`, `structure_equations` | a copy of envelope_cz_bottom_index taken at the card start and each rezone; nz without a surface convection zone |
| `has_h_shell` | true if the model burns hydrogen in a shell rather than in the core | logical | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `observables`, `step_driver`, `yrec_output` |  |
| `h_shell_end_index` | last point of the H-burning shell | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `observables` | valid when has_h_shell |
| `h_shell_midpoint_zone` | first point where X exceeds half its surface value | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `observables`, `run_driver`, `step_driver` | valid when has_h_shell |
| `h_shell_zone_begin` | first point of the H-burning shell (X > ctrl%atime(1)) | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `observables`, `run_driver`, `step_driver`, `yrec_output` | valid when has_h_shell |

### 3.4 Step snapshot

`rezone` interpolates the old model onto the new mesh into the `*_start` and `old_*` arrays and then copies them into the live arrays, so after the rezone these arrays hold the start-of-step model on the new mesh. Without a rezone (`ctrl%helium_flash_active`), `rezone_or_snapshot` copies the structure into the `*_start` arrays itself.

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `kinetic_energy_rot_old(json)` | kinetic_energy_rot at the start of the step | erg | `step_driver` (rezone_or_snapshot) | `structure_equations` | rotating models; enters the energy equation as d(E_rot)/dt |
| `logP_start(json)` | log P at the start of the step, on the current mesh | log10 dyn/cm^2 | `rezoning` (seed_new_point_run, regrid_structure_variables); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `rezoning`, `timestep_control`, `mid_timestep_model` | `rezone` splines the old model onto the new mesh into the *_start arrays; without a rezone `rezone_or_snapshot` copies them |
| `logT_start(json)` | log T at the start of the step | log10 K | `rezoning` (regrid_structure_variables); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `rezoning`, `burn_step`, `timestep_control`, `mid_timestep_model` | as logP_start |
| `logR_start(json)` | log r at the start of the step | log10 cm | `rezoning` (regrid_structure_variables); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `rezoning`, `convection_geometry`, `step_driver`, `timestep_control`, `mid_timestep_model` | as logP_start |
| `luminosity_lsun_start(json)` | L at the start of the step | Lsun | `rezoning` (seed_new_point_run, regrid_structure_variables); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `rezoning`, `step_driver`, `timestep_control`, `mid_timestep_model` | as logP_start |
| `logRho_start(json)` | log rho at the start of the step | log10 g/cm^3 | `rezoning` (regrid_structure_variables); `step_driver` (rezone_or_snapshot) | `rezoning`, `burn_step`, `mid_timestep_model` | as logP_start; not set by `snapshot_step_start` at a card start |
| `xa_start(15,json)` | composition at the start of the step | - | `rezoning` (regrid_composition_terms); `step_driver` (advance_composition_and_age, converge_with_rotation); `evolve_rotation` (evolve_angular_momentum, burn_light_elements_over_substep) | `rezoning`, `burn_step`, `step_driver`, `timestep_control`, `evolve_rotation`, `mid_timestep_model` | copied before the composition step; the corrector passes and the structure-rotation iteration restart from it |
| `convective_flag_start(json)` | convective_flag at the start of the step | logical | `rezoning` (retally_angular_momentum) | `mid_timestep_model` | rotating models only (`retally_angular_momentum`) |
| `log_Teff_start` | log Teff at the start of the step | log10 K | `rezoning` (seed_new_point_run); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `timestep_control` | for the HR-motion timestep limit |
| `surface_luminosity_lsun_start` | L(nz) at the start of the step, before the rezone | Lsun | `rezoning` (seed_new_point_run); `starting_model` (snapshot_step_start); `step_driver` (rezone_or_snapshot) | `timestep_control` | for the HR-motion timestep limit |
| `convective_turnover_timescale_old` | convective_turnover_timescale of the previous model | s | `observables` (refresh_turnover_timescale); `starting_model` (read_starting_model) | `magnetic_braking` | copied after each converged model (and at a model read); the wind interpolates between it and the current value at fracstep |
| `pphot0` | pphot of the previous model | log10 dyn/cm^2 | `observables` (refresh_turnover_timescale); `starting_model` (read_starting_model) | `magnetic_braking` | copied with convective_turnover_timescale_old; pphot itself changes during the step |
| `rate_li6_start(json)` | Li6 proton-capture rate rho X lambda at the start of the step | 1/s | `burn_step` (store_li_be_rates); `evolve_rotation` (burn_light_elements_over_substep) | `burn_step` | with `job%use_extended_composition`; zero outward of log T = 6 |
| `rate_li7_start(json)` | Li7 proton-capture rate at the start of the step | 1/s | `burn_step` (store_li_be_rates); `evolve_rotation` (burn_light_elements_over_substep) | `burn_step` | as rate_li6_start |
| `rate_be9_start(json)` | Be9 proton-capture rate at the start of the step | 1/s | `burn_step` (store_li_be_rates); `evolve_rotation` (burn_light_elements_over_substep) | `burn_step` | as rate_li6_start |
| `envelope_overshoot_depth_start` | envelope-overshoot depth alpha_ov H_P below the surface convection zone at the start of the step | cm | `step_driver` (rezone_or_snapshot); `evolve_rotation` (burn_light_elements_over_substep) | `convection_geometry` | was `pressure_scale_height_start` until 2026-10 |
| `deuterium_burning_rate_start(json)` | deuterium burning rate at the start of the step | 1/s | `rezoning` (regrid_composition_terms); `burn_step` (compute_burn_rates) | `rezoning`, `burn_step`, `mid_timestep_model` | stored by the composition step's predictor pass |
| `orig_specific_angular_momentum(json)` | j_rot at the start of the structure-rotation iteration | cm^2/s | `step_driver` (converge_with_rotation) | `step_driver` | restored on each repeat pass when `job%num_rotation_structure_iters` > 1 |
| `orig_composition(15,json)` | xa_start at the start of the structure-rotation iteration | - | `step_driver` (converge_with_rotation) | `step_driver` | as orig_specific_angular_momentum |
| `old_omega(json)` | omega at the start of the step | rad/s | `rezoning` (regrid_rotation_variables, retally_angular_momentum); `starting_model` (snapshot_step_start) | `rezoning`, `timestep_control`, `evolve_rotation`, `transport_checks` | rezone splines omega onto the new mesh into it, then copies it back; read by the omega timestep limit and the transport |
| `old_specific_angular_momentum(json)` | j_rot at the start of the step | cm^2/s | `rezoning` (regrid_rotation_variables, retally_angular_momentum) | `rezoning` | read only by `rezone` itself (interpolation buffer) |
| `old_hg(json)` | mean_gravity at the start of the step | cm/s^2 | `rezoning` (retally_angular_momentum) | `mid_timestep_model` | legacy name (HG) |
| `old_mean_radius(json)` | mean_radius at the start of the step | cm | `rezoning` (regrid_rotation_variables, retally_angular_momentum) | `rezoning` | read only by `rezone` itself (interpolation buffer) |
| `old_eta_squared(json)` | eta_squared at the start of the step | - | `rezoning` (regrid_rotation_variables, retally_angular_momentum) | `rezoning` | read only by `rezone` itself (interpolation buffer) |

### 3.5 Step working state

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `delta_log_pressure` | shift of log P of the convection zone for the entropy of accreted matter | log10 | `mass_change` (adjust_cz_for_accreted_entropy) | `mass_change`, `structure_equations` | zero at present: `change_mass` zeroes envelope_specific_entropy first |
| `delta_log_temperature` | shift of log T for the same | log10 | `mass_change` (adjust_cz_for_accreted_entropy) | `mass_change`, `structure_equations` | as delta_log_pressure |
| `disk_gate_age_gyr` | age of the last converged model, for the disk-locking tests | Gyr | `step_driver` (converge_with_rotation) | `mass_change`, `rotation_scratch_lib` | a copy of dage taken after convergence |
| `convective_turnover_timescale` | convective turnover time one pressure scale height above the base of the surface convection zone | s | `turnover_timescale` (compute_turnover_timescale, turnover_from_interior_new, time_scale_height_above_base, time_bottom_cell_crossing); `run_driver` (begin_kind_card); `mass_change` (measure_convective_envelope) | 8 modules (core/ 6, io/ 1, rotation/ 1) | from the stitched model, every converged model; zero at a card start and without a surface convection zone |
| `pphot` | log10 photospheric pressure of the last atmosphere integration that reached the photosphere | log10 dyn/cm^2 | `envelope_integration` (atm_get) | `observables`, `starting_model`, `history_output`, `magnetic_braking` | set inside `atm_get`, i.e. by any envelope integration (the triangle's too) |
| `fracstep` | fraction of the step elapsed at the current rotation sub-step | - | `starting_model` (read_starting_model); `evolve_rotation` (evolve_angular_momentum, spin_down_without_transport) | `magnetic_braking` | 0.5 after a model read |
| `rate_li6(json)` | Li6 proton-capture rate rho X lambda at the end of the step | 1/s | `burn_step` (store_li_be_rates) | `burn_step`, `evolve_rotation` | with `job%use_extended_composition`; computed on the converged model for the Li/Be burn |
| `rate_li7(json)` | Li7 proton-capture rate at the end of the step | 1/s | `burn_step` (store_li_be_rates) | `burn_step`, `evolve_rotation` | as rate_li6 |
| `rate_be9(json)` | Be9 proton-capture rate at the end of the step | 1/s | `burn_step` (store_li_be_rates) | `burn_step`, `evolve_rotation` | as rate_li6 |
| `envelope_overshoot_depth_end` | envelope-overshoot depth alpha_ov H_P at the end of the step | cm | `step_driver` (burn_light_elements); `evolve_rotation` (burn_light_elements_over_substep) | `convection_geometry`, `evolve_rotation` | was `pressure_scale_height_end` until 2026-10 |
| `cz_base_radius_prev` | radius of the base of the surface mixed region found by the last call within the step | cm | `convection_geometry` (locate_surface_mixing_base); `step_driver` (rezone_or_snapshot) | `convection_geometry` | with `job%use_extended_composition`; zeroed at each step start to mean 'none yet' |
| `envelope_cz_base_zone_prev` | zone of that base, the starting zone of the next call within the step | - | `convection_geometry` (locate_surface_mixing_base) | `convection_geometry` |  |
| `deuterium_burning_rate(json)` | deuterium burning rate of each shell above the nuclear cutoff | 1/s | `burn_step` (compute_burn_rates) | `burn_step`, `mid_timestep_model` | zero while the surface D is below 1e-11 |
| `accreted_mass` | mass accreted during the step | g | `mass_change` (mix_in_accreted_composition) | `burn_step`, `mid_timestep_model` | set only when the star accretes; its readers check that it does |
| `cp(json)` | specific heat at constant pressure, from `compute_shell_physics` | erg/(g K) | `shell_physics` (compute_shell_physics) | `rezoning`, `mid_timestep_model`, `viscosity` | for rotation (viscosity, transport); compare scp |
| `mu(json)` | mean molecular weight | - | `shell_physics` (compute_shell_physics) | `rezoning`, `evolve_rotation`, `mid_timestep_model` | from `compute_shell_physics`, for rotation |
| `qdt(json)` | dln rho/dln T at constant P | - | `shell_physics` (compute_shell_physics) | `rezoning`, `mid_timestep_model` | from `compute_shell_physics`, for rotation |
| `thdif(json)` | radiative thermal diffusivity K | cm^2/s | `viscosity` (viscos) | `rezoning`, `mid_timestep_model` | from `viscos`, for rotational mixing |
| `visc(json)` | kinematic viscosity (radiative plus ionic) | cm^2/s | `viscosity` (viscos) | `rezoning`, `mid_timestep_model` | from `viscos`, for rotational mixing |
| `dlog_temperature_dt(json)` | d log T/dt over the last step, for the next step's first guess | 1/s | `rezoning` (regrid_composition_terms); `run_driver` (begin_kind_card); `step_driver` (solve_structure) | `rezoning`, `step_driver` | was `temperature_entropy_term` until 2026-10 |
| `dlog_pressure_dt(json)` | d log P/dt over the last step | 1/s | `rezoning` (regrid_composition_terms); `run_driver` (begin_kind_card); `step_driver` (solve_structure) | `rezoning`, `step_driver` | was `pressure_entropy_term` until 2026-10 |
| `dln_luminosity_dt(json)` | d ln L/dt over the last step, 2 (L - L_start)/(L + L_start)/dt | 1/s | `rezoning` (regrid_composition_terms); `run_driver` (begin_kind_card); `step_driver` (solve_structure) | `rezoning`, `step_driver` | was `luminosity_entropy_term` until 2026-10 |
| `dlog_radius_dt(json)` | d log r/dt over the last step | 1/s | `rezoning` (regrid_composition_terms); `run_driver` (begin_kind_card); `step_driver` (solve_structure) | `rezoning`, `step_driver` | was `radius_entropy_term` until 2026-10 |
| `env_points` (type `envelope_structure_state`, `envelope_points`) | the envelope of the latest envelope integration, point by point: log P, log T, log m, log rho, log r, X, Z, the gradients, the opacity and its derivatives, ... | cgs; logs base 10 | `envelope_integration` (atm_get, through the integration kernel) | `stitched_model`, `envelope_refit`, `envelope_rebuild`, `starting_model` | the module variable `env_struct` until 2026-10 |
| `atm_points` (type `atmosphere_structure_state`) | the atmosphere of the latest integration, point by point | cgs; logs base 10 | `envelope_integration` (atm_get) | `stitched_model` | the module variable `atmo_struct` until 2026-10 |
| `atm_point` (type `atmosphere_point_state`) | the atmosphere integration's current point: log tau, log P, log T, log rho, the opacity, the ionization fractions | logs base 10 | `envelope_integration` (atm_get) | `mass_change` | the mass change falls back on its last values when the Allard lookup fails; a module variable until 2026-10 |
| `central_saha_state` | the Saha state the central EOS call carries from one model to the next | - | `observables` (compute_observables) | `observables` | observables' module variable `ksaha_center` until 2026-10 |
| `omega_mid(json)` | omega of the latest rotation sub-step's model, the structure interpolated to the sub-step's end | rad/s | `mid_timestep_model` (rotate_midpoint_model); `secular_transport`; `evolve_rotation` (through `enforce_rotation_profile`) | `evolve_rotation`, `secular_transport`, `mid_timestep_model` | `rot_scr%omega_mid` until 2026-10. It outlives the step: the next step's first solid-body sub-step starts its Newton iteration from its mean |

### 3.6 Newton

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `elim_coeff(4,2,json)` | block-elimination coefficients c(k,j,i) of `henyey_elimination` | - | `henyey_elimination` (henyey_eliminate); `structure_equations` (eliminate_shell_pair) | `henyey_elimination`, `henyey_iteration` | see `henyey_elimination`'s header for the relations |
| `elim_rhs(4,json)` | right-hand sides e(k,i): the residuals, then (after `henyey_solve`) the corrections dlog P, dlog T, dlog r, dL | log10, Lsun | `henyey_elimination` (henyey_eliminate); `henyey_iteration` (henyey_iterate); `structure_equations` (eliminate_shell_pair) | `henyey_elimination`, `henyey_iteration` |  |
| `log_pressure_delta(json)` | change of log P since the start of the step | log10 | `step_driver` (solve_structure); `henyey_iteration` (apply_henyey_corrections) | `step_driver`, `henyey_iteration`, `structure_equations` | zeroed (or set to the first-guess extrapolation) before the solve, then accumulated by every correction; feeds eps_grav |
| `log_temperature_delta(json)` | change of log T since the start of the step | log10 | `step_driver` (solve_structure); `henyey_iteration` (apply_henyey_corrections) | `step_driver`, `henyey_iteration`, `structure_equations` | as log_pressure_delta |
| `max_residual(4)` | per equation (P, T, r, L): the largest \|residual\|, then the largest \|correction\|, then the signed correction there | log10, Lsun | `henyey_elimination` (henyey_eliminate); `henyey_iteration` (linearize_structure_equations, find_largest_corrections, test_henyey_convergence) | `henyey_elimination`, `henyey_iteration` | misnamed: holds residuals or corrections depending on the stage of the iteration |
| `max_correction_index(4)` | point of the largest correction of each variable | - | `henyey_iteration` (find_largest_corrections) | `henyey_iteration` | only for the convergence log |
| `surface_bc(6)` | outer boundary-condition rows from the envelope fit, in `henyey_solve`'s form | - | `henyey_iteration` (form_surface_bc_rows) | `henyey_iteration` |  |
| `newton_iterations` | Newton iterations of the step's solve, counted over its levels up to the last convergence | - | `henyey_iteration` (accept_small_residuals, set_henyey_correction_factor) | `run_log` | for the run-log model line |

### 3.7 Diagnostics written during the solve

`compute_henyey_coefficients` stores these at every shell on every Newton iteration. Rule 5 of the six (diagnostics written once, after convergence) is the planned replacement.

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `eps_grav(json)` | eps_grav = -T ds/dt of each zone | erg/(g s) | `step_driver` (solve_structure); `structure_equations` (add_shell_energy_budget) | `structure_equations` | was `gravitational_luminosity` until 2026-10 |
| `luminosity_breakdown(9)` | luminosity per channel `i_lum_*` (pp1, pp2, pp3, CNO, 3 alpha, neutrinos, gravothermal, alpha captures, the extra energy sources), scaled to sum to L(nz) | Lsun | `observables`, `henyey_iteration`, `structure_equations`, `model_files` (7 procedures) | 8 modules (core/ 6, io/ 2) | summed every iteration, rescaled again after convergence (`renormalize_luminosity_breakdown`); `compute_timestep` reads it; read from YREC7/MODEL2 files, zeroed by the .mod reader |
| `eps_extra(json)` | the extra energy sources' sum at each zone (`core/structure/energy_sources.f90`) | erg/(g s) | `structure_equations` (add_shell_energy_budget) | `stitched_model` (profile column `eps_extra`) | zero without a registered source (2026-10) |
| `eps_total(json)` | nuclear energy generation pp1 + pp2 + pp3 + CNO + 3 alpha | erg/(g s) | `structure_equations` (store_zone_energy_channels) | `stitched_model`, `rezoning`, `structure_equations`, `pulsation_output`, `mid_timestep_model` | zero below the burning cutoff |
| `eps_channels(7,json)` | per channel `i_eps_*`: rows 1-5 the fractions of eps_total, row 6 eps_nu, row 7 eps_grav | rows 1-5: -; 6-7: erg/(g s) | `structure_equations` (store_zone_energy_channels) | `stitched_model`, `rezoning`, `mid_timestep_model` | mixed units by row |
| `beta(json)` | gas-pressure fraction P_gas/P | - | `structure_equations` (store_zone_thermodynamics) | `stitched_model` |  |
| `eta(json)` | electron degeneracy parameter | - | `structure_equations` (store_zone_thermodynamics) | `stitched_model` | 0 unless the EOS treats the gas as fully ionized |
| `opacity_zone(json)` | Rosseland mean opacity kappa | cm^2/g | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `observables`, `stitched_model`, `rezoning`, `pulsation_output`, `mid_timestep_model`, `viscosity` | also written by `compute_shell_physics` (rezone, rotation, card start); the last writer wins |
| `gradr(json)` | radiative temperature gradient grad_rad | - | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `observables`, `stitched_model`, `rezoning`, `convection_geometry`, `mid_timestep_model` | as opacity_zone |
| `gradT(json)` | actual temperature gradient d ln T/d ln P | - | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `observables`, `stitched_model`, `rezoning`, `composition_step`, `pulsation_output`, `mid_timestep_model` | as opacity_zone |
| `grada(json)` | adiabatic gradient grad_ad | - | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `observables`, `stitched_model`, `rezoning`, `convection_geometry`, `pulsation_output`, `mid_timestep_model` | as opacity_zone |
| `conv_vel(json)` | convective velocity from mixing-length theory | cm/s | `shell_physics` (compute_shell_physics); `structure_equations` (store_zone_thermodynamics) | `stitched_model`, `rezoning`, `mid_timestep_model` | as opacity_zone |
| `neutrino_flux_total(10)` | solar-neutrino fluxes at the Earth by slot `i_nu_*` (slots 9-10: He3+He3 and He3+He4 rates in the same units) | 10^10 /(cm^2 s) | `neutrino_flux_output` (write_neutrino_flux_table); `observables` (compute_snu_rates); `structure_equations` (start_energy_bookkeeping, add_shell_energy_budget) | `neutrino_flux_output`, `observables`, `monte_carlo`, `structure_equations`, `history_output` | only with `ctrl%calc_neutrinos`; also overwritten by the card-start neutrino table; unit from `neutrino_flux_output`'s header |
| `adiabatic_index_gamma1(json)` | first adiabatic exponent Gamma_1 | - | `structure_equations` (store_zone_thermodynamics) | `stitched_model`, `pulsation_output` |  |
| `pulse_deps_dlnrho(json)` | d eps/dln rho (absolute, not logarithmic) | erg/(g s) | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` | zero below the burning cutoff |
| `pulse_deps_dlnt(json)` | d eps/dln T (absolute) | erg/(g s) | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` | as pulse_deps_dlnrho |
| `pulse_dlnkap_dlnrho(json)` | dln kappa/dln rho | - | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` |  |
| `pulse_dlnkap_dlnt(json)` | dln kappa/dln T | - | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` |  |
| `pulse_specific_heat(json)` | c_P | erg/(g K) | `structure_equations` (store_pulse_derivatives) | `stitched_model` |  |
| `pulse_mean_molecular_weight(json)` | mean molecular weight mu = R/(specific gas constant) | - | `structure_equations` (store_pulse_derivatives) | `stitched_model` |  |
| `pulse_dlnrho_dlnt(json)` | dln rho/dln T at constant P | - | `structure_equations` (store_pulse_derivatives) | `stitched_model`, `pulsation_output` | for the pulsation output |
| `pulse_electron_mean_weight_inverse(json)` | 1/mu_e | - | `structure_equations` (store_pulse_derivatives) | `stitched_model` |  |

### 3.8 Diagnostics written after convergence

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `stitched` (type `stitched_model_record`, `stitched_model_state`) | the stitched model: the converged star (interior, envelope and atmosphere) on one extended grid, with its profile and pulse columns | cgs | `stitched_model` (build_stitched_model) | `observables`, `turnover_timescale`, `yrec_output`; the profile and pulse writers through their argument | stitched_model's module arrays until 2026-10; allocated, and zeroed, by the first build |
| `es_circulation_velocity(json)` | Eddington-Sweet circulation velocity at each interface | cm/s | `circulation` (circulation_velocities, clear_velocities, apply_velocity_scales, average_with_previous, zero_transport_in_core) | `stitched_model`, `circulation` | written by the rotation transport, last iteration of the last sub-step |
| `secular_shear_velocity(json)` | secular (or dynamical) shear-instability velocity | cm/s | `circulation` (circulation_velocities, clear_velocities, apply_velocity_scales, average_with_previous, zero_transport_in_core) | `stitched_model`, `circulation` | as es_circulation_velocity |
| `gsf_circulation_velocity(json)` | GSF-instability velocity | cm/s | `circulation` (circulation_velocities, clear_velocities, apply_velocity_scales, average_with_previous, zero_transport_in_core) | `stitched_model`, `circulation` | as es_circulation_velocity |
| `am_diffusion_coeff(json)` | angular-momentum diffusion coefficient D_omega | cm^2/s | `secular_transport` (set_diffusion_coefficients) | `stitched_model` | written by the rotation transport (last sub-step of the step); 0 without rotation |
| `mixing_diffusion_coeff(json)` | rotational mixing coefficient D_mix | cm^2/s | `secular_transport` (set_diffusion_coefficients) | `stitched_model` | as am_diffusion_coeff |
| `cl37_snu_rate` | chlorine (37Cl) detector capture rate | SNU | `observables` (compute_snu_rates) | `observables`, `monte_carlo`, `history_output` | 0 without `ctrl%calc_neutrinos` |
| `ga71_snu_rate` | gallium (71Ga) detector capture rate | SNU | `observables` (compute_snu_rates) | `observables`, `monte_carlo`, `history_output` |  |
| `log_R_surface` | log10 of the photospheric radius from L = 4 pi R^2 sigma Teff^4 | log10 Rsun | `observables` (compute_surface_globals) | `observables`, `history_output`, `run_log`, `yrec_output` |  |
| `log_g_surface` | log10 surface gravity G M/R^2 | log10 cm/s^2 | `observables` (compute_surface_globals) | `observables`, `stop_conditions`, `history_output`, `yrec_output` |  |
| `nu_max` | frequency of maximum oscillation power, scaled from the Sun's by g and Teff | uHz | `observables` (compute_seismic_observables) | `stop_conditions`, `history_output` |  |
| `delta_nu_rho` | large frequency separation scaled from the mean density | uHz | `observables` (compute_seismic_observables) | `history_output` |  |
| `delta_nu` | asymptotic large separation from the sound-travel time | uHz | `observables` (compute_seismic_observables) | `history_output` |  |
| `delta_Pg` | asymptotic l = 1 g-mode period spacing of the inner cavity | s | `observables` (compute_seismic_observables) | `history_output` | 0 without a cavity |
| `total_moment_of_inertia` | moment of inertia of the star | g cm^2 | `observables` (compute_moment_of_inertia) | `observables`, `history_output` |  |
| `cz_moment_of_inertia` | moment of inertia of the surface convection zone | g cm^2 | `observables` (compute_rotation_observables) | `observables`, `history_output` |  |
| `rotation_period_days` | surface rotation period 2 pi/omega(nz), capped at 9999 | d | `observables` (compute_rotation_observables) | `history_output` | 0 without rotation |
| `surf_velocity_kms` | surface equatorial velocity omega(nz) R | km/s | `observables` (compute_rotation_observables) | `history_output` |  |
| `h_shell_bot_mass` | mass coordinate of the bottom of the H-burning shell | Msun | `observables` (compute_h_shell_boundaries) | `history_output` | 0 without a shell |
| `h_shell_mid_mass` | mass coordinate of its middle | Msun | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `h_shell_top_mass` | mass coordinate of its top | Msun | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `h_shell_bot_radius` | radius of the bottom of the H-burning shell | r/R | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `h_shell_mid_radius` | radius of its middle | r/R | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `h_shell_top_radius` | radius of its top | r/R | `observables` (compute_h_shell_boundaries) | `history_output` |  |
| `central_log10_temperature` | central temperature, extrapolated from the first point | log10 K | `observables` (compute_central_conditions) | `observables`, `monte_carlo`, `history_output` |  |
| `central_log10_pressure` | central pressure | log10 dyn/cm^2 | `observables` (compute_central_conditions) | `observables`, `monte_carlo`, `history_output` |  |
| `central_log10_density` | central density, from the EOS | log10 g/cm^3 | `observables` (compute_central_conditions) | `observables`, `monte_carlo`, `history_output` |  |
| `envelope_mass` | mass of the surface convection zone | Msun | `observables` (locate_surface_cz_base) | `monte_carlo`, `history_output` |  |
| `envelope_radius` | radius of the base of the surface convection zone | Rsun | `observables` (locate_surface_cz_base, refresh_turnover_timescale); `starting_model` (read_starting_model) | `monte_carlo`, `history_output` | written twice per model: `locate_surface_cz_base`, then overwritten by the turnover-time walk (`compute_turnover_timescale`) |
| `central_beta` | central P_gas/P | - | `observables` (compute_central_conditions) | `history_output` |  |
| `central_degeneracy_eta` | central electron degeneracy parameter | - | `observables` (compute_central_conditions) | `history_output` |  |
| `core_cz_mass` | mass of the convective core | Msun | `observables` (locate_core_cz) | `history_output` | 0 without one |
| `envelope_cz_temperature` | temperature at the base of the surface convection zone | K | `observables` (locate_surface_cz_base) | `history_output` | linear, not log |
| `envelope_cz_density` | density there | g/cm^3 | `observables` (locate_surface_cz_base) | `history_output` |  |
| `envelope_cz_pressure` | pressure there | dyn/cm^2 | `observables` (locate_surface_cz_base) | `history_output` |  |
| `envelope_cz_opacity` | opacity there | cm^2/g | `observables` (locate_surface_cz_base) | `history_output` |  |

### 3.9 Step control

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `model_diverged_flag` | the structure solve of this step diverged | logical | `run_driver` (begin_kind_card); `step_driver` (solve_structure, converge_with_rotation, solve_level) | `run_driver`, `step_driver` | a set flag stops the run on the next pass |
| `punch_pending_flag` | a model has been computed, so the next helium-flash card saves its starting model | logical | `run_driver` (run_yrec); `step_driver` (reload_model_if_diverged) | `run_driver` | legacy LPUNCH |
| `recompute_envelope_triangle` | build a new envelope-fit triangle at the next solve | logical | `run_driver` (begin_kind_card); `mass_change` (update_total_mass); `step_driver` (evolve_step, reload_model_if_diverged, solve_structure, solve_level) | `step_driver` | set by a model read, mass change, the envelope rebuild or `ctrl%lnew0`; cleared by `compute_surface_boundary` |
| `reset_triangle` | rezone moved the outer point; redo the triangle at level 2 | logical | `run_driver` (begin_kind_card); `step_driver` (rezone_or_snapshot, solve_level) | `step_driver` |  |
| `dt` | timestep of the step being taken (next step's after `compute_timestep`) | s | `run_driver` (begin_kind_card); `mass_change` (limit_mass_change_timestep); `step_driver` (evolve_step, reload_model_if_diverged, solve_structure) | `run_driver`, `stop_conditions`, `mass_change`, `step_driver` | negative means 'solve, do not age' (a rescale or relaxation); made positive after convergence |
| `dt_saved` | copy of dt; when the end-age stop clips the next timestep to zero, the timestep just taken is restored from it | s | `run_driver` (begin_kind_card); `step_driver` (evolve_step) | `step_driver` |  |
| `max_domega_frac` | largest fractional change of omega over the step | - | `run_driver` (begin_kind_card); `step_driver` (evolve_step, converge_with_rotation) | `step_driver` | written by `compute_timestep` and `evolve_angular_momentum` |
| `timestep_yr` | dt in years | yr | `run_driver` (begin_kind_card); `step_driver` (evolve_step, reload_model_if_diverged, advance_composition_and_age) | `run_driver`, `step_driver`, `model_files`, `profile_output`, `run_log` | in .mod |

### 3.10 Run and card bookkeeping

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `use_ttau_relation` | the run has switched atm_choice to the gray atmosphere because Teff left the Kurucz or Allard table | logical | `envelope_integration` (atm_get); `surface_boundary` (switch_atmosphere_choice); `controls_reader` (adopt_canonical_names) | `surface_boundary` | switches back when Teff returns to the table |
| `mixing_length_alpha` | mixing-length parameter of the current run card | - | `run_driver` (begin_kind_card) | `starting_model`, `history_output`, `model_files`, `yrec_output`, `microphysics_settings` | in .mod; from `job%mixing_length_array(nk)` at each card start (the solar calibration changes the later cards' entries) |
| `age_at_target_radius` | age interpolated at the target radius | Gyr | `calibration` (check_star_calibration) | `calibration` |  |
| `age_prev_model` | age of the previous model | Gyr | `calibration` (check_star_calibration) | `calibration` |  |
| `log_l_at_target_radius` | log(L/Lsun) interpolated at the target radius | log10 Lsun | `calibration` (check_star_calibration) | `calibration` |  |
| `log_l_at_target_radius_prev_run` | the same, from the previous calibration run | log10 Lsun | `calibration` (check_star_calibration) | `calibration` |  |
| `log_l_prev_model` | log(L/Lsun) of the previous model | log10 Lsun | `calibration` (check_star_calibration) | `calibration` |  |
| `r_rsun_prev_model` | R/Rsun of the previous model | Rsun | `calibration` (check_star_calibration, setup_star_calibration) | `calibration` | was `log_r_prev_model` until 2026-10 |
| `termination_reason` | why the last run card ended | text | `run_driver` (end_of_card_calibration); `stop_conditions` (converged_model_is_nan, structure_limit_stop_triggered, abundance_stop_triggered, init_stop_conditions) | `stop_conditions`, `run_log` | default 'model budget exhausted'; read by the run summary |
| `settling_suspended_reported` | the settling routines have reported that settling is suspended (fully convective, or no hydrogen or helium to settle) | logical | `composition_step` (settle), `mid_timestep_model` (burn_and_settle), through `diffuse_elements` | the same | the message prints once per suspension; "SETTLING RESUMED" clears it. `diff_scr`'s until 2026-10 |
| `star_found_flag` | the star calibration found a track through the target L and R; the final run is armed | logical | `calibration` (check_star_calibration, setup_star_calibration) | `run_driver`, `stop_conditions` |  |
| `just_passed_target_radius_flag` | this model has just crossed the target radius | logical | `calibration` (check_star_calibration, setup_star_calibration) | `stop_conditions` |  |
| `pulsation_mass_msun` | star_mass at the start of the card | Msun | `run_driver` (begin_kind_card) | `profile_output` | for the profile and pulsation headers |
| `envelope_hydrogen_fraction` | surface X of the first run card: the reference X of the EOS mixture and of 1/mu_ion in `compute_shell_physics` | - | `starting_model` (update_surface_mixture) | `model_rescaling`, `starting_model`, `run_setup`, `shell_physics`, `run_cards`, `transport_checks` | set once (`update_surface_mixture`, first card) and kept for the run on purpose; the controls read (`interpret_kind_cards`) and `kap_init` use it earlier, while it is still 0 |
| `envelope_metal_fraction` | surface Z of the EOS mixture: the first card's, or the last Z rescale's | - | `model_rescaling` (rescale_interior_z); `starting_model` (update_surface_mixture); `run_cards` (interpret_kind_cards) | `model_rescaling`, `starting_model`, `run_setup`, `shell_physics`, `run_cards`, `transport_checks` | also set from a card's rescale Z by the controls read (`interpret_kind_cards`) and by `rescale_interior_z` |
| `zenvm` | part of the surface Z not in C, N and O | - | `starting_model` (update_surface_mixture) | `starting_model`, `microphysics_settings` | passed to convection's settings (`envelope_metals_no_cno`) |
| `amuenv` | inverse mean atomic weight of the surface mixture, sum X_i/A_i | 1/amu | `starting_model` (update_surface_mixture) | `model_rescaling`, `starting_model`, `shell_physics`, `run_cards`, `transport_checks` |  |
| `fxenv(12)` | number fractions of the 12 surface-mixture elements, slot `ix_*` | - | `starting_model` (update_surface_mixture) | `model_rescaling`, `starting_model`, `run_cards` | passed to the EOS (`eos_set_mixture`) |
| `dlum_dx` | fixed d log L/dX of the solar calibration's Newton step | - | `calibration` (check_solar_calibration) | `calibration`, `monte_carlo` | a constant (-3.78) set on every check; printed in the Monte-Carlo summary |
| `drad_dx` | fixed d log R/dX | - | `calibration` (check_solar_calibration) | `calibration`, `monte_carlo` | constant |
| `dlum_dalpha` | fixed d log L/d alpha | - | `calibration` (check_solar_calibration) | `calibration`, `monte_carlo` | constant |
| `drad_dalpha` | fixed d log R/d alpha | - | `calibration` (check_solar_calibration) | `calibration`, `monte_carlo` | constant |
| `log_l_prev` | log(L/Lsun) at the previous calibration check | log10 Lsun | `calibration` (check_solar_calibration) | `calibration` | only printed |
| `log_r_prev` | log(R/Rsun) at the previous calibration check | log10 Rsun | `calibration` (check_solar_calibration) | `calibration` | only printed |
| `delta_x` | the Newton correction to X | - | `calibration` (check_solar_calibration) | `calibration` |  |
| `delta_alpha` | the Newton correction to the mixing length | - | `calibration` (check_solar_calibration) | `calibration` |  |
| `solar_calibration_converged` | the solar calibration has converged | logical | `calibration` (check_solar_calibration, setup_solar_calibration) | `monte_carlo`, `run_driver` | was `solar_calibration_active` until 2026-10 |
| `envelope_helium_fraction` | surface Y = 1 - X - Z - X(He3) of the first card | - | `starting_model` (update_surface_mixture) | `shell_physics`, `transport_checks` |  |
| `envelope_he3_fraction` | surface He3 of the first card | - | `starting_model` (update_surface_mixture) | `shell_physics`, `transport_checks` |  |

### 3.11 Constants set at setup

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `atm_hras` | the HSRA T(tau) fit at tau = 2/3 | log10 K | `run_setup` (set_up_constants_and_tables) | `envelope_integrands`, `envelope_integration` | a constant set at run time; the declaration comment calls it Krishna-Swamy, the code uses the HSRA fit |
| `tenv` | log T at which an envelope integration stores its state for the surface boundary (stored_envelope_state): the mean of ctrl tenv0 and tenv1 | log10 K | `derived_controls` (derive_controls) | `envelope_integration` | set by `derive_controls` |
| `atm_choice_initial` | atm_choice as read from the deck | - | `controls_reader` (adopt_canonical_names) | `surface_boundary` | the run may change `job%atm_choice`; this keeps the original |
| `cross_section_scale(17)` | per-reaction scale of the cross section S0 against Bahcall & Ulrich (1988) | - | `monte_carlo` (apply_monte_carlo_parameters); `controls_reader` (adopt_canonical_names); `user_inputs` (set_cross_section_scales) | `burn_step` | slots `r_*`/`xs_*` of net_lib; slots 1-3 and 16 rescaled per Monte-Carlo run |
| `qs0e_scale(8)` | scale of each reaction's S'/S against Solar Fusion I | - | `user_inputs` (reset_rate_scale_defaults, set_first_derivative_scales) | `burn_step` | slots `r_*`, `iq_be7p` |
| `qqs0ee_scale(8)` | scale of each reaction's S''/S against Solar Fusion I | - | `user_inputs` (reset_rate_scale_defaults, set_second_derivative_scales) | `burn_step` |  |
| `o16_gamma_scale` | scale of the N15(p,gamma)O16 branch | - | `user_inputs` (reset_rate_scale_defaults, set_cross_section_scales) | `burn_step` |  |
| `c12_alpha_scale` | scale of the N15(p,alpha)C12 branch | - | `user_inputs` (reset_rate_scale_defaults, set_cross_section_scales) | `burn_step` |  |
| `solar_luminosity_cgs` | Lsun | erg/s | `monte_carlo` (apply_monte_carlo_parameters); `controls_reader` (adopt_canonical_names) | 13 modules (core/ 6, io/ 3, rotation/ 3, setup/ 1) | from the namelist clsun; a Monte-Carlo run scales it |
| `log10_solar_luminosity` | log10 Lsun | log10 erg/s | `monte_carlo` (apply_monte_carlo_parameters); `run_setup` (set_up_constants_and_tables) | `observables`, `envelope_integration`, `starting_model`, `run_driver`, `mass_change`, `magnetic_braking`, `microphysics_settings` | as solar_luminosity_cgs |
| `ln10_over_lsun` | ln(10)/Lsun | s/erg | `monte_carlo` (apply_monte_carlo_parameters); `run_setup` (set_up_constants_and_tables) | `structure_equations` | not a logarithm |
| `solar_mass_cgs` | Msun | g | `run_setup` (set_up_constants_and_tables) | 11 modules (core/ 8, rotation/ 2, io/ 1) | copy of `ctrl%solar_mass_cgs` |
| `solar_radius_cgs` | Rsun | cm | `controls_reader` (adopt_canonical_names) | 8 modules (core/ 6, io/ 1, rotation/ 1) | from the namelist crsun |
| `log10_solar_radius` | log10 Rsun | log10 cm | `run_setup` (set_up_constants_and_tables) | `observables`, `run_driver`, `yrec_output` |  |
| `solar_bolometric_magnitude` | M_bol of the Sun, 4.79 | mag | `run_setup` (set_up_constants_and_tables) | `yrec_output` | a constant set at run time |
| `use_two_z_tables` | load a second opacity-table set at another Z (Z ramp or Z diffusion) | logical | `run_driver` (run_yrec) | `microphysics_settings` |  |

### 3.12 Scratch

| member | meaning | unit | written by | read by | notes |
|---|---|---|---|---|---|
| `be7_mass_fraction_zone(json)` | Be7 mass fraction of each shell in the neutrino table | - | `neutrino_flux_output` (write_neutrino_flux_table) | `neutrino_flux_output` | used only inside `write_neutrino_flux_table` |
| `neutrino_flux_zone(10,json)` | each shell's contribution to the neutrino fluxes, then its normalized profile | 10^10 /(cm^2 s) | `neutrino_flux_output` (write_neutrino_flux_table) | `neutrino_flux_output` | used only inside `write_neutrino_flux_table` |
| `vfc(json)` | Zahn (1992) coupling factor between angular-momentum transport and mixing, <= 1 | - | `circulation` (zahn_coupling_factor); `secular_transport` (compute_unperturbed_velocities) | `circulation` | rotation transport only (`circulation`, `secular_transport`) |
| `envelope_specific_entropy` | mass-weighted mean specific entropy of the convective envelope (ideal gas, up to a constant) | erg/(g K) | `mass_change` (change_mass, average_envelope_energies) | `mass_change` | used only inside `mass_change`; `change_mass` zeroes it right after computing it, so the accretion shifts it drives are zero |

### 3.13 Never referenced

None: `old_shell_mass` and `converged_zone` were removed in 2026-10.

## 4. The input bundles

| member | meaning | written by | read by | notes |
|---|---|---|---|---|
| `ctrl` (type `controls_state`) | The namelist controls, by their canonical names (the run-card arrays and latches are in `job`). | `controls_reader` (`read_controls`, which resets it to the defaults before each read) and `controls_sync` (`store_controls_to_star`) | 36 modules | Not changed after the controls read. Members: `src/defaults/controls_registry.md`, `src/defaults/physics.defaults`; the component list is generated into `controls_state_def.inc`. |
| `job` (type `star_job`) | The run's job: the table paths, the surface-mixture weights, the Monte-Carlo samples (`s11_rate` ... `age_target`), the run-card cursor `nk`, and the run-card arrays and latches of `job_controls_def.inc`. | the controls read, then 17 modules in all (calibration, Monte Carlo, run and step drivers, the solver, the atmosphere switch, ...) | 37 modules | Mutable: the calibration rewrites the next cards, the stop conditions disarm thresholds, the atmosphere switch changes `atm_choice`, the solver ramps `fcorr`. Members: `src/defaults/control.defaults` and the registry. |

## 5. Candidates for removal

None of these has been removed; each needs a look before it goes.

**Removed in 2026-10**: `old_shell_mass` and `converged_zone` (never referenced), `old_moment_of_inertia` and `log10_solar_mass` (written, never read); the diagnostics no writer emitted, `fxion_zone`, `scp` (a copy of `pulse_specific_heat`), `he3_he3_luminosity_zone`, `he3_burning_luminosity_zone` and `pulse_dlnrho_dlnp`, and the always-zero `valfmlt`, `vphmlt` and `vcmxmlt` (zero in the F77 code too: their source scalars were only ever zeroed, for a pulsation format since retired); the outputs no caller used, now locals of the callers: `qiw` (`omega_from_j`'s di_domega, recomputed by every user), `chosen_dt`, `ikut_flag` and `convective_velocity` (of `compute_timestep` and `read_starting_model`), `mixed_zone_bounds_no_overshoot` and `radiative_zone_bounds` (of the convection-zone search), and `initial_composition_code` (read with a Model2 file); `envelope_cz_log_radius` became a local of `locate_surface_cz_base`; `dlnrho_dlnp` and `dlnrho_dlnt`, which reached the next coefficient pass only as an EOS starting value every EOS family overwrites, became locals of `compute_henyey_coefficients` (`docs/bug_sweep.md` O23).

## 6. Names that mislead

Renamed in 2026-10 (the old name in parentheses), because each invited a wrong reading: `eps_grav` (`gravitational_luminosity`; eps_grav in erg/(g s), not a luminosity); `dlog_temperature_dt`, `dlog_pressure_dt`, `dln_luminosity_dt`, `dlog_radius_dt` (the four `*_entropy_term`; rates of change for the next step's first guess); `envelope_overshoot_depth_start`, `_end` (`pressure_scale_height_start`, `_end`; alpha_ov H_P in cm); `triangle_log_teff`, `triangle_log_l`, `triangle_orientation` (`trial_log_temperature`, `trial_log_luminosity`, `trial_sign_flag`); `stored_vertex_index` (`istore_flag`); `r_rsun_prev_model` (`log_r_prev_model`; linear R/R☉); `solar_calibration_converged` (`solar_calibration_active`); `age_gyr` (`dage`).

Still to watch:

- `max_residual` holds residuals or corrections, depending on the stage of the Newton iteration.
- `ln10_over_lsun` is ln(10)/L☉, not a logarithm.
- `star_mass` is in M☉ without a unit suffix: MESA's name, kept.

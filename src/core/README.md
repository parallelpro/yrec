# core/ — the star layer: the run, the step and the structure solve

This directory drives a YREC run and evolves the model. It owns the model state (`star%`, declared in `state/star_info_lib.f90`) and calls the physics domains (eos, kap, atm, net, convection, rotation, element_diffusion) and the numerics through their public modules.

Its subdirectories follow the layers of `docs/TARGET_ARCHITECTURE.md` §4.1:

```
  run/          the run: controls, tables, run cards, calibration, Monte Carlo, stops
    step/       one model: the stages in order, and the next timestep
      structure/    the Henyey solve of the structure equations
      envelope/     the outer layers: atmosphere and envelope integrations
      model/        the model itself: the starting model, rescaling, rezoning
    diagnostics/  the converged star: the stitched model and derived quantities
```

`main.f90` is the command-line entry point; `run/yrec_capi.f90` is the same engine for C and Python (`pyyrec`). Paths below are under `src/` unless they start at the repository root.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `main` | The command-line program: calls `run_yrec` and exits 1 on a positive `ierr` (§3). |
| `core_messages` | The routes of core's messages (§3): `core_report`, `core_terminal`, `core_log`, `core_stdout`, `core_stdout_and_log`; the lines that keep their trailing blanks (`core_terminal_values`, `core_log_values`, `core_stdout_values`, `log_envelope_start_retry`, `log_mixed_zone_bounds`); `rewind_core_log`. |
| **run/** | |
| `run_driver` | `run_yrec`: the whole run (§2.1). Reads the controls, opens the run's files, sets up the tables, and goes through the run cards inside the Monte-Carlo loop, evolving each card's model with `evolve_step`; its steps `begin_calibration`, `begin_kind_card`, `end_kind_card` and `end_of_card_calibration`. |
| `run_setup` | `prepare_run` (the run plan's envelope Z, then `set_up_constants_and_tables`, then the Monte-Carlo samples) and `set_up_constants_and_tables` (the solar units on the star and log10 G, the one run-time constant; then the kap/eos/atm settings and table loads, `setup/settings_install`'s `install_run_settings`). |
| `calibration` | The solar and star calibrations: their setup before the run (`setup_solar_calibration`, `setup_star_calibration`), the check after each triple or at each crossing of R* (`check_solar_calibration`, `check_star_calibration`), and the star calibration's end-of-card verdict (`end_star_calibration_card`). |
| `monte_carlo` | The Monte-Carlo runs: the run range (`first_mc_run`, `last_mc_run`), the samples (`setup_monte_carlo_runs`), the per-run hook (`begin_monte_carlo_run`), and the summaries (`write_run_summaries`). |
| `stop_conditions` | The stops after each converged model (`check_stop_conditions`): the end age (`reached_end_age`), the central-abundance table (disarmed at a card's start where the starting model has passed a threshold, `init_stop_conditions`), the structure limits (log L, Teff, log g, nu_max) and the star calibration's target radius; the `step_status` protocol; the NaN guard (`converged_model_is_nan`) and the rotation check (`check_rotation_initialised`). |
| `yrec_capi`, `yrec_reset` | The C entry point (`yrec_run`: a path over 256 characters is refused; the run's own units are closed after it), and the re-entry prologue (`yrec_run_prologue`, `close_run_units`) that makes a second run in one process start from a fresh star. |
| `yrec_version` | The version stamp of the run log. |
| **step/** | |
| `step_driver` | `evolve_step`: one model (§2.2). Its locals live for the step only; what outlives it is in the star. `refresh_turnover_timescale`: the turnover time of the converged model, for the next step. |
| `mass_change` | Mass loss or accretion at the surface (`change_mass`): the rate, the timestep limits, the entropy shift; gathers the convection zone from the star, calls the kernels, and stores the results. |
| `mass_change_kernel` | The kernels of the mass change, with strict in/out arguments and no star: the convection zone's mass and new masses (`convection_zone_mass`, `change_shell_masses`), the angular momentum a thermal wind carries off (`wind_specific_angular_momentum`, `remove_wind_angular_momentum`), and the accreted composition mixed into the zone (`mix_accreted_composition`). |
| `composition_step` | `advance_composition`: the operator-split composition change over the step at fixed structure: burning, convective mixing, semiconvection, settling. |
| `burn_step` | The star-layer side of nuclear burning: reads the model, calls net/, stores the results; the light elements (`burn_li_be`, `store_start_li_be_rates`, `store_end_li_be_rates`). |
| `shell_physics` | `compute_shell_physics(s, ierr)`: every zone's EOS, opacity, gradients and viscosity, for rotation and the convection-zone edges. |
| `convection_geometry` | Where the convection zones end: the fraction of a zone at which del_rad - del_ad changes sign (`cz_base_fraction`), the lowest zone above a radius (`first_zone_above`), and the base of the surface mixed region for the Li/Be burn (`locate_surface_mixing_base`). |
| `timestep_control` | `compute_timestep(s)`: the next timestep. Gathers the limits' inputs from the star, calls the kernels of `timestep_limits` and takes the tightest (`take_tightest_limit`); the fixed-timestep override and the end age. |
| `timestep_limits` | The timestep limits as star-blind kernels behind a settings record (`timestep_settings`): `structure_change_limit`, `omega_change_limit`, `hydrogen_burning_limit`, `helium_burning_limit`, `hr_motion_limit`, and `tightest_timestep`. |
| **structure/** | |
| `henyey_iteration` | `refresh_surface_bc` and `henyey_iterate`: the surface boundary of one solve level, then its Newton iterations, as named steps with explicit arguments (linearize, renormalize the luminosity channels, accept small residuals, the surface rows, solve, the largest corrections, the convergence test, the correction factor, the log lines, apply, rotation). |
| `structure_equations` | `compute_henyey_coefficients`: across the model, gathers each shell from the star, calls the kernel, eliminates, and stores the shells' output values; `luminosity_channel_sum`. |
| `shell_equations` | `evaluate_shell_equations`: the kernel -- the four equations' terms and partials at one shell (eos, opacity, gradients, nuclear, gravitational and extra energy), with strict in/out arguments and no star. Its header holds the equations and the partials' conventions. |
| `henyey_elimination` | The block elimination (`henyey_eliminate`) and back-substitution (`henyey_solve`): kernels on arrays, which return `numerics_singular` for a zero pivot; `eliminate_star_shell_pair` hands the elimination the star's arrays. |
| `surface_boundary` | `compute_surface_boundary`: the outer boundary condition from the triangle of envelope integrations: switches the atmosphere, keeps the star's triangle around the model, integrates the envelopes at its vertices, and stores the fit. |
| `triangle_fit` | The kernels of the triangle -- `start_triangle`, the point-in-triangle test (`edge_side`, `vertex_to_reflect`), `reflect_vertex`, `reflect_until_inside`, and `fit_triangle_coefficients` (log R, log L and log Teff as planes in the fitting point's log P and log T) -- with strict in/out arguments and no star. |
| `energy_sources` | The extra energy sources' registry (`add_energy_source`, `clear_energy_sources`, `evaluate_energy_sources`): a source is a subroutine of an `energy_point`; their sum and partials enter the luminosity equation (§6, "Add an energy source"). |
| **envelope/** | |
| `envelope_integration` | `atm_get`: the driver of one envelope integration: fills the kernel's settings from the star, calls it, applies its two effects on the star (the photospheric pressure, the switch to a gray atmosphere), and reports its events (`report_envint_events`). |
| `envelope_integration_kernel` | `integrate_envelope_atmosphere`: the atmosphere and envelope integrations from the surface to the fitting point, star-blind; it prints nothing and returns its notes, solver diagnostics and failures as events (`envint_events`) beside `ierr`. |
| `envelope_integrands` | Their right-hand sides, in optical depth and in pressure (`compute_atmosphere_derivs`, `compute_envelope_derivs`; star-blind). |
| `envelope_refit`, `envelope_rebuild` | Appending the envelope's points inside a new fitting mass as shells (`append_envelope_points(s, ...)`), and moving the fitting point outward when the star accretes (`rebuild_envelope(s, target, ierr)`). |
| (state/envelope_points) | The outer layers of the latest envelope integration, point by point: records the kernel fills and the star holds (`s%env_points`, `s%atm_points`, `s%atm_point`). |
| **model/** | |
| `starting_model` | `read_starting_model(s, ierr)`: each card's starting model, read or reused, rescaled if asked, its envelope refitted; its steps take the star and a `starting_model_state` record. |
| `model_rescaling` | `rescale_model`: a new X, Z, species, mass or core mass. |
| `rezoning` | `rezone(s, ierr)`: the new mesh of each step and the model on it, as steps that take the star and a `rezone_state` record (the flagged points, the new mesh); one step per mesh criterion (§6). |
| `regridding` | The kernels of the rezone (no star): the shell masses of a mesh (`compute_shell_masses`), the profiles splined onto new points (`regrid_profile`, `regrid_composition`, `regrid_structure`, `regrid_rotation`), the angular momentum totals (`compute_angular_momentum_totals`). |
| `shell_boundaries` | `find_shell_boundaries` (the kernel, no star) and `locate_shell_boundaries(s)` (on the star): the convective core, the surface convection zone, the hydrogen-burning shell. |
| **diagnostics/** | |
| `stitched_model` | `build_stitched_model`: the whole converged star (interior, envelope, atmosphere) and every quantity the output needs; the profile columns. |
| `observables` | `compute_observables`: the per-model observables, one procedure per theme; nothing it computes feeds the physics (§2.2). |
| `turnover_timescale` | `compute_turnover_timescale`: the convective turnover time near the base of the surface convection zone, from the stitched model (called by the step and by `read_starting_model`). |
| `neutrino_flux_output` | `write_neutrino_flux_table`: the opt-in neutrino-flux table of each card's starting model. |

Kept commented-out code, not compiled, each file with a header saying where it came from and how to revive it, all three commented out with the advection-diffusion transport on 2026-10-12 (with `rotation/retired/advection_diffusion.f90.txt`): `model/retired/rezoning_theta_term.f90.txt` (the transport's theta term on the new mesh), `step/retired/rotation_scaled_solar_wind.f90.txt` (the rotation-scaled solar wind) and `structure/retired/rotation_stability_copies.f90.txt` (the structure's copies for the rotation stability coefficients).

## 2. The flow

### 2.1 The run

`run_yrec` (`run/run_driver.f90`), called by `main.f90` and by the C entry `yrec_run`, runs one job:

1. `yrec_run_prologue` (`run/yrec_reset.f90`): the run's star, a copy of the process-start star; at a re-entry, the previous run's log and final-model units closed first (`close_run_units`).
2. `read_controls` (io): the controls into `s%control` and `s%physics`, and the run plan.
3. `open_run_files` (io): the run log, the copy of the deck, the starting and final models, the Monte-Carlo and isochrone files; `output_open` (io): the history and profile output.
4. `prepare_run` (`run/run_setup.f90`): the run plan's envelope Z; the solar units, log10 G, the kap/eos/atm settings and tables; the Monte-Carlo samples. Then `warn_inconsistent_controls` (io).
5. For each Monte-Carlo run (`first_mc_run` to `last_mc_run`; one run without Monte Carlo):
   1. `begin_monte_carlo_run` (`run/monte_carlo.f90`): every run copy of a control re-armed, the state a run leaves cleared, the run's samples applied; `install_monte_carlo_settings` (setup).
   2. `begin_calibration`: `setup_solar_calibration` or `setup_star_calibration` (`run/calibration.f90`) lay out the run cards.
   3. For each run card: `begin_kind_card` (the card's composition and mixing length, `install_card_settings`, `read_starting_model`, `check_rotation_initialised`, `init_stop_conditions`, the neutrino table, the output headers, `locate_shell_boundaries`, `compute_timestep`); then `evolve_step` (§2.2) up to the card's model count, until a stop (`step_kind_card_done`) or the star calibration's crossing of R* (`step_leave_run_loop`, which goes on to the next card); `end_kind_card`; `end_of_card_calibration` (`check_solar_calibration` after each triple, `end_star_calibration_card`), which may end the card list.
   4. `log_run_summary` (io): why the run ended (`s%termination_reason`) and the final model; `write_run_summaries`: the Monte-Carlo summary.
6. `output_close` (io).

The stops (`run/stop_conditions.f90`, `check_stop_conditions`, after each converged model is written): the end age, the central D/X/Y table, the structure limits, then the star calibration's target radius.

### 2.2 A step

`evolve_step` (`step/step_driver.f90`) advances the model by one timestep:

1. **Mass and composition** (`advance_composition_and_age`): the mass change (`change_mass`), then the composition step over the timestep (`advance_composition`'s predictor pass), and the new age.
2. **The mesh** (`rezone_or_snapshot`): the rezone, and the start-of-step snapshot.
3. **The structure** (`solve_structure`): the Henyey solve in levels. Level 1 uses the previous model's envelope triangle; level 2 checks the surface boundary first; level 3 follows, and from level 3 on the composition step is redone before each solve (`solve_level`: `refresh_surface_bc`, the corrector pass, then `henyey_iterate` on the equations of `compute_henyey_coefficients`).
4. **Level 4 and rotation** (`converge_with_rotation`): level 4 runs for every model, inside the structure-rotation iteration, which has one pass without rotation; then the convergence verdict, the convection zones mixed fully and, with rotation, the angular momentum (`evolve_angular_momentum`) and the shape factors.
5. **After convergence:** the shell boundaries (`locate_shell_boundaries`), the light elements (`burn_light_elements`), the envelope rebuild after accretion (`rebuild_envelope`), the next timestep (`compute_timestep`), the stitched model (`build_stitched_model`), the observables (`compute_observables`), the turnover time (`refresh_turnover_timescale`), the output (`output_write_model`), and the stop conditions (`check_stop_conditions`).

The observables only fill the output's quantities. The turnover time, which the next step's deuterium cap and wind read, is a stage of the step: `refresh_turnover_timescale` walks the stitched model (`compute_turnover_timescale`) and copies the turnover time and `photosphere_logP` to their `_prev` members, the start values the wind interpolates from during the next step. The walk also writes the radius of the convection zone's base it finds into `envelope_radius_rsun`, over `compute_observables`' value, so the output has the walk's.

A step that diverges stops the run: nothing cuts the timestep.

## 3. Conventions

- **The star comes in as an argument:** a core procedure that touches it takes `type(star_info), intent(inout) :: s` first and reads `s%member`. `run_yrec` owns the run's star (there is no global one since 2026-10) and passes it down; `yrec_reset`'s prologue resets it to its process-start state at the start of every run. A driver reads the model from the star; it takes the model's arrays as arguments only when its callers pass different models (`burn_li_be` and `store_end_li_be_rates`: the converged model or a rotation sub-step's; `advance_composition`'s zone outputs at the corrector passes), and its header says so.
- **Kernels** take no star (star rule 6) and print nothing. `shell_equations`, `triangle_fit`, `mass_change_kernel`, `timestep_limits`, `regridding`, `envelope_integration_kernel` and `envelope_integrands` are star-blind files (`tools/check_boundaries.py`'s `STAR_BLIND_FILES`); `find_shell_boundaries`, `henyey_eliminate` and `henyey_solve` are star-free procedures in files that also hold a star driver. Their drivers gather from the star and store back. A new kernel's file goes in `STAR_BLIND_FILES`, with ok/FAIL checks in a test program and their mutations (§4).
- **Steps with arguments:** `rezone`, `read_starting_model` and `henyey_iterate` are drivers whose steps are module procedures with explicit arguments (the star, the driver's state record where it has one, `rezone_state` or `starting_model_state`, documented at its type, and `ierr`); everything else a step needs is its own local. `compute_timestep`'s steps are too, but for `take_tightest_limit`, the one line that lists its host's limits. The other long drivers' steps (`evolve_step`'s stages among them) are internal procedures that share their host's locals.
- **Units:** CGS unless a name carries its unit (`_msun`, `_lsun`, `_yr`, `_gyr`); `log` is log₁₀.
- **The structure variables** are log₁₀ P, log₁₀ T, log₁₀ R and L in L☉, against log₁₀ m. The difference equations and their partials follow `docs/TARGET_ARCHITECTURE.md` §4.4 and `shell_equations`' header.
- **The aging flag** (`evolve_model_flag` in `evolve_step`) marks a model solved but not aged (on a rescale card); `s%dt` is negative only when a model file stores a negative timestep, and `compute_timestep` keeps its magnitude.
- **Errors:** a procedure that can fail returns `ierr` up to `run_yrec`, which returns it: `main` exits 1 on a positive `ierr`, and `yrec_run` returns it to C. The one negative code is `envelope_termination` (-2, `envelope_integration_kernel`): a Bulirsch-Stoer step of the envelope integration failed, the historical "solution diverged" end, after which the outputs up to that model are valid and `main` exits 0. The kernels' own codes (`numerics_bad_input`, `numerics_singular`, the eos and kap codes) pass through unchanged. `yrec_run` refuses a namelist path longer than 256 characters with 1, before anything runs.
- **Messages** go through `core_messages`, never through a `write` of core's own. The caller formats into an array of lines (`character(len=message_length) :: lines(n)`, set to `unwritten_line`, then an internal write) and calls one route: `core_report` (the terminal and the run log: an error the run stops on, or a warning it goes on with), `core_terminal`, `core_log` (the run log: the per-model notes and the solver diagnostics), `core_stdout` (standard output, what the F77 `write(*)` meant; it differs from `core_terminal` only before io assigns the units) or `core_stdout_and_log`. A line whose form a route would change has its own procedure: the list-directed lines go through `core_terminal_values`, `core_log_values` and `core_stdout_values`, which keep their trailing blanks, and two lines of the envelope integration and the composition step through `log_envelope_start_retry` and `log_mixed_zone_bounds`. The kernels print nothing: `integrate_envelope_atmosphere` returns events and `atm_get` reports them; `henyey_eliminate` and `henyey_solve` return `numerics_singular` and their drivers report it. `check_boundaries`' `check_no_io` keeps the files of `STAR_BLIND_FILES` and core's entries in `KERNEL_FILES` (`shell_boundaries`, `model_rescaling`, `calibration`, `observables`, `turnover_timescale`, `stitched_model`) free of I/O statements. The I/O core keeps is file handling: the starting model's `rewind` (`starting_model`), the Monte-Carlo sample read, summary writes and rewinds (`monte_carlo`), and `close_run_units`' inquire and close (`yrec_reset`).
- **Controls:** core is the star layer, so its drivers read `s%control` and `s%physics` directly (io/README.md §5, "Add a control the star layer reads directly"); its kernels take what they need as arguments or in a settings record (`timestep_settings`, which `timestep_settings_of` builds from the star). The physics domains read only the settings records that `setup/` builds.
- **Kept by decision:** the extra energy sources' registry is the process's configuration, not model state, and outlives `run_yrec` until `clear_energy_sources` (`energy_sources`' header); `run_yrec`'s locals, its star among them, are static (a blanket `save`), and the prologue resets the star; a diverging step stops the run (§2.2); `envelope_termination` exits 0, as the F77's "solution diverged" did.

## 4. Tests

Build the programs with `make tests` in `src/`, then run each root `test_*.py` with pytest from the repository root (CI runs them exactly and with `YREC_TEST_TOLERANT=1`). The nine programs with a manifest print ok/FAIL checks and exit nonzero on a FAIL; the pinned outputs are `core/test/expected_test_*.out`, regenerated from `src/<program>` after a deliberate change, with the reason in the commit.

| Test (driver) | What it checks | Mutations |
|---|---|---|
| `test_structure` (`test_structure.py`) | Every partial of `evaluate_shell_equations` against central differences at six shells of a TAMS model (OPAL tables, 3e-3) and four of a birthline model (Yale EOS, 1e-4, with the eps_grav partials); the `extra` case: a registered test source's partials and the size of its term. The partials marked `known` are documented approximations: the convective gradient's dependence on L is not linearized (O1), and its P and T derivatives near the adiabat rest on the EOS's derivatives of grad_ad. | 10, on the Yale case with the extra source; among them 18869257's 1/chi of the eps_grav partials reverted |
| `test_triangle_fit` (`test_triangle_fit.py`) | The triangle's kernels on hand-made triangles: the planes at the vertices and inside, the slopes, the degenerate fits refused with `ierr`, the point-in-triangle test, the reflection, the walks and their cut-off after `max_reflections`. | 15 |
| `test_mass_change` (`test_mass_change.py`) | The mass-change kernels' conservation on a synthetic 40-shell model: the shell masses, the zone's scaling, dm(1) of a zone that reaches the centre, the wind's angular momentum, and accretion of every basic species into a zone of uniform composition (the one case the mixing conserves, O7). | 12; among them bc874267 reverted and J_cz measured after the change (91a6b567) |
| `test_henyey_elimination` (`test_henyey_elimination.py`) | `henyey_eliminate`/`henyey_solve` on random systems of 2, 3 and 40 shells against a dense LU solve of the Newton system of the module header (7e-15), the residual maxima exactly, and the five zero-pivot contracts. | 15 |
| `test_envelope_integration` (`test_envelope_integration.py`) | The integrands against their equations; the integration's status, events and end points; every Bulirsch-Stoer step against RK4 between the kernel's points; the vertex interpolation; the thin envelope; the retried start; the error contracts (maxstp in the atmosphere and the envelope, the Bulirsch-Stoer limits, a step that underflows). | 29 |
| `test_timestep_limits` (`test_timestep_limits.py`) | Each timestep limit against its formula, with its switches, clipping and floors, and the tightest within the growth limit. | 33 |
| `test_regrid` (`test_regrid.py`) | The rezone's regridding kernels on a synthetic 1 Msun model of 60 points regridded onto 65: the shell masses sum to M and follow the centre, interior and outer rules; profiles linear in log m are reproduced; the points the meshes share keep their values; each species' mass, X + Y + Z, the angular momentum and the rotational kinetic energy change by their pinned amounts (the spline conserves none of them); a failed spline returns osplin's code. 19 checks, byte-compared with `expected_test_regrid.out` (residuals below 1e-12 compared as zero). | 19 |
| `test_model_kernels` (`test_model_kernels.py`) | `find_shell_boundaries` on synthetic 40-point profiles: the convective core's edge, the surface convection zone's base, a fully convective model, and the hydrogen shell's beginning, middle and end (by the luminosity's step and by the surface X), on either side of `x_core_min`. 9 checks. | 9 |
| `test_run_contracts` (`test_run_contracts.py`) | The star calibration on a star set up in the test: a crossing of R* before the last card sets up the next pair; a crossing on the last card (card 50) writes no card past it and ends with a verdict; `end_star_calibration_card`'s verdicts. The C entry called as a host would, with the host's own unit open: a missing deck returns a positive status, the host's unit stays open across calls and re-entries, the run's log and final-model units are closed after a run, a path of exactly 256 characters runs (the cm_env_overshoot deck), and a control or physics path one longer is refused with 1. 17 checks. | 17 |
| `test_reentry` (`test_reentry.py`) | Two `run_yrec` calls in one process: the second run's output files are byte-identical to a fresh process's (the solar deck, and a carry-over case where the first run sets a control the second omits). | none (a byte comparison of whole runs) |
| `test_stop_conditions.py` | The four structure-limit stops (log L, Teff, log g, nu_max) end the solar deck's evolving card at its first model with their STOP line; the seismic history columns are opt-in. | none (deck runs) |
| `test_monte_carlo.py` | Two Monte-Carlo runs of identical samples are the same run (M169: no state carries over between runs). | none |
| `test_pyyrec.py` | The C entry through pyyrec (ctypes): a run, and two runs in one process, byte-identical to a fresh process's (`make lib` first). | none |
| `test_other_energy.py` | The energy-source hook end to end (`make other_energy` first): a registered source that adds nothing changes nothing, and a positive one brightens the star and fills the history column `lum_extra`. | none |
| `test_error_exit.py` | Error paths end the process with status 1 and an ERROR line; among them core's inputs: a malformed or missing starting model, the envelope triangle's size and the envelope's Bulirsch-Stoer controls (refused at the controls read). | none |
| `test_physics.py` | The standard solar model is still a solar model at the solar age, within loose bounds. | none |

Every deck runs core, so the deck pins are its integration test: the 35 pinned run logs (`*/standard/`) that `src/tools/compare_decks.py --set all` and `test_all.py` compare, and `test_mesa_baselines.py`'s two MESA-format baselines. No pinned deck calibrates to a target star, rescales a horizontal-branch model or writes the neutrino table; check a change there by running such a deck against the previous commit's binary.

Each check of the nine programs with a manifest was confirmed to fail under a mutation of the code it covers: the manifests `core/test/<program>.mutations.toml` hold 159 mutations, and `python3 src/tools/mutate.py core` (after `make tests` in `src/`) reported all 159 caught on 2026-11-13 (test_envelope_integration 29, test_henyey_elimination 15, test_mass_change 12, test_model_kernels 9, test_regrid 19, test_run_contracts 17, test_structure 10, test_timestep_limits 33, test_triangle_fit 15).

Under the development flags (`make devel_tests`, `python3 src/tools/devel_tests.py`) the programs run clean but for one known failure, O4: the starting model's stitched atmosphere reads `env_log10_radius(0)` (`diagnostics/stitched_model.f90:475`), a bounds error that stops `test_structure`'s three cases and `test_run_contracts` (`devel_tests.py`'s `KNOWN_FAILURES`).

## 5. Coverage

75.5% of core/'s 5865 executable lines run (4426; 2026-11-13, at r2-docs e7551d23), under the seven argument-free test programs (`test_triangle_fit`, `test_mass_change`, `test_henyey_elimination`, `test_envelope_integration`, `test_timestep_limits`, `test_regrid`, `test_model_kernels`) and 30 decks: the 16 configuration-matrix decks, the 12 pinned testsuite solar decks, cm_accretion and o24_rescale_rotating. Measured with
```
python3 src/tools/coverage.py --work DIR --dir core --test test_triangle_fit --test test_mass_change \
    --test test_henyey_elimination --test test_envelope_integration --test test_timestep_limits \
    --test test_regrid --test test_model_kernels \
    --deck examples/run_config_matrix/cm_<each>.inlist --deck testsuite/Test_solar<each>.inlist \
    --deck src/tools/compare_decks/cm_accretion.inlist --deck src/tools/compare_decks/o24_rescale_rotating.inlist
```
(a -O0 --coverage build in DIR; one `--deck` per deck).

| Directory | Lines | Run | Coverage |
|---|---|---|---|
| step/ | 1147 | 1014 | 88.4% |
| envelope/ | 909 | 788 | 86.7% |
| structure/ | 829 | 699 | 84.3% |
| diagnostics/ | 916 | 752 | 82.1% |
| model/ | 1411 | 922 | 65.3% |
| run/ | 593 | 224 | 37.8% |
| `main`, `core_messages` | 60 | 27 | 45.0% |

The kernels: `henyey_elimination`, `mass_change_kernel`, `shell_physics` and `triangle_fit` 100%, `envelope_integrands` 98.0%, `timestep_limits` 97.7%, `shell_boundaries` 97.4%, `regridding` 97.0%, `envelope_integration_kernel` 94.6%, `shell_equations` 86.2%. What does not run:
- outside the measure: `calibration` and `yrec_capi` (0%) run under `test_run_contracts`, and the solar calibration's decks are .nml1/.nml2 pairs; `energy_sources` (5.7%) runs under `test_other_energy` and `test_structure`'s extra case; `coverage.py` takes neither programs with arguments nor deck pairs;
- `monte_carlo` (34.1%): cm_monte_carlo stops at once in the measure, because its parameter file is not copied beside it (bug_sweep M313);
- `neutrino_flux_output` (0%): no deck asks for the table;
- `starting_model` (56.4%) and `model_rescaling` (38.9%): the horizontal-branch and core rescales, the core extension, the envelope refits of the ZAMS decks, and the error paths; agent M's run with all 41 testsuite inlists (2026-11-13, at r2-core-m 80d07e6c) reached 76.7% and 64.7%, and model/ 77.2%;
- `turnover_timescale` (68.7%): the bottom cell's crossing time when several cells lie in the envelope, and the spline failures (90.0% in agent M's run);
- `convection_geometry` (42.2%): the Li/Be burn's envelope-overshoot branches, which no deck reaches (bug_sweep M303);
- `stop_conditions` (58.7%), `core_messages` (37.7%) and the drivers' remaining lines: error paths and the stops and messages no deck in the measure meets.

## 6. How to extend

**Add an energy source** (tidal heating, magnetic heating, an injected energy), with no edit of YREC's source:
1. Write the source as a subroutine with `energy_sources`' interface, `subroutine my_source(pt, eps, deps_dlnrho, deps_dlnt, ierr)` with `type(energy_point), intent(in) :: pt` (m, dm, log T, log P, log rho, log r, L, xa, dt): eps in erg/g/s, its partials at constant T and at constant rho, `ierr` nonzero to stop the solve.
2. In a program, register it before the run, `call add_energy_source('name', my_source, ierr)`, then `call run_yrec(ierr)`: `examples/other_energy/yrec_extra_heat.f90` is the template, and `make other_energy` in `src/` the build rule to copy.
3. The term enters the luminosity equation at every interior shell (`evaluate_shell_equations`' step `form_extra_energy_terms`), its partials the Jacobian, and its luminosity the history column `lum_extra`. A source above the fitting point has no effect (the envelope is integrated at constant L).
4. Test it: register it in a copy of `test_structure`'s `extra` case (module `test_structure_source`) to check its partials against finite differences, and run a deck as `test_other_energy.py` does.

Traced (2026-11-12, in a scratch copy of the tree): a program that registers a constant 1 erg/g/s source before `run_yrec`, linked as `make other_energy` links `yrec_extra_heat`; cm_base's history column `lum_extra` is 0.5191 Lsun, M/Lsun for 1 erg/g/s.

**Add a timestep limit:**
1. In `core/step/timestep_limits.f90`, a kernel: a function of the settings record and the arrays it needs, returning the limit in seconds, `dt_unlimited` when it does not apply. A setting it needs goes in `timestep_settings` and in `timestep_settings_of` (`timestep_control.f90`), from a control (a new control is io/README.md §5's recipe).
2. In `compute_timestep` (`timestep_control.f90`), a local for its result, the call after the other limits, and the result in `take_tightest_limit`'s list (one line).
3. In `core/test/test_timestep_limits.f90`, a check of its formula and of its switch, and a mutation of each in `test_timestep_limits.mutations.toml`; `python3 src/tools/mutate.py core --program test_timestep_limits`.
4. Re-pin what moves (every deck where it binds), in a measured commit.

Traced (2026-11-12, in a scratch copy of the tree) with a cap kernel (`cap_limit`: a settings field, the kernel, its call and one entry in `take_tightest_limit`, a check in `test_timestep_limits`): with a 1e10 s cap, cm_base's largest timestep went from 7882 yr to 187 yr, all below the cap.

**Add a stop condition.** A stop ends the current run card when the converged model crosses a limit. The structure limits are a table in `run/stop_conditions.f90`'s `structure_limit_stop_triggered`: one row per quantity, with an upper and a lower limit per card (a limit at its ±1d99 default is off).
1. The control (io/README.md §5, "Add a control"): declare it in `io/controls_buffer.f90` beside `Teff_upper_limit`, per card and off by default, with its description comment:
   ```fortran
   ! Per kind card: stop the card when log10 of the central temperature [K]
   ! exceeds this (1d99 = off)
         real(dp) :: log_Tc_upper_limit(max_runs) = 1.0d99
   ```
   add its row to `defaults/controls_registry.tsv` (a copy of `Teff_upper_limit`'s: a `control`, used by the run), and its name to the `&control` namelist list in `io/legacy_controls.f90`.
2. The row: in `structure_limit_stop_triggered`, raise `nlim` by one and add the quantity's name to `qname` (7 characters), its value to `qval` (here `s%logT(1)`), and its limits to `qup` and `qlo` (`-1.0d99` when it has no lower limit).
3. `make` in `src/` reruns the controls generators (io/README.md, "Generated files"); `python3 -I tools/check_star_contracts.py --fix` adds the member to the function's contract; `make check`.
4. A case in `test_stop_conditions.py`: the solar deck with the limit set so that its evolving card stops early, and the `STOP:` line in the run log. Re-pin `io/test/expected_controls/` (`YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py`: the pristine dump gains the control).

The check runs after the observables and the output (`check_stop_conditions` is the last stage of `evolve_step`), so `s%` holds the converged model's values. A stop that is not a limit (a new kind of event) is a further test inside `check_stop_conditions` that writes `s%termination_reason` and returns `step_kind_card_done`, before the star calibration's test.

Traced on a throwaway branch (2026-11-12): the four hand edits above, `make`, `check_star_contracts --fix`, `make check`; the solar deck with `log_Tc_upper_limit(2) = 7.0d0` stops its evolving card at its first model with log T_c > 7, printing `STOP: log_Tc =  7.0027E+00 above log_Tc_upper_limit =  7.0000E+00` and ending with `run finished: log_Tc above log_Tc_upper_limit`.

**Add a mesh criterion.** The rezone keeps the flagged points and spaces the new mesh evenly in chi between them (`model/rezoning.f90`). Each criterion is a step of `flag_fixed_points` that adds the points it wants kept with `append_flag_point(rz, zone)`; `sort_flag_points` then sorts the list and removes duplicates.
1. Write the step beside the others, taking the star and the rezone's state:
   ```fortran
   ! ---------------------------------------------------------------
   ! The first point where L reaches half the surface luminosity.
   !   star reads: none
   !   star writes: none
   subroutine flag_half_luminosity(s, rz)
         type(star_info), intent(inout) :: s
         type(rezone_state), intent(inout) :: rz
         integer :: i
         do i = 2, s%nz
            if (s%luminosity_lsun(i) >= 0.5_dp*s%luminosity_lsun(s%nz)) then
               call append_flag_point(rz, i)
               exit
            end if
         end do
   end subroutine flag_half_luminosity
   ```
2. Call it in `flag_fixed_points` before `sort_flag_points`, and name it in the list in that procedure's header. The list holds `max_flag_points` (100); `append_flag_point` drops and reports (once per rezone) what does not fit, so call a criterion that must win before the scan of jumps (`flag_jumps`).
3. A threshold is a setting: read it from `s%physics%chi_grid_scale(ichi_...)` (a new slot is a control, io/README.md §5) or name it as a parameter of the module.
4. `python3 -I tools/check_star_contracts.py --fix` fills the contract; `make check`.
5. Check it on a deck with `report_solver_diagnostics = .true.` in `&physics`: each rezone writes its `FLAG-POINTS` line to the run log. A new criterion changes the mesh, so compare_decks shows which decks move: it is a measured commit, with its pins re-seeded.

Traced on a throwaway branch (2026-11-12): the step and its call, `check_star_contracts --fix`, `make check`; cm_base with the diagnostics on keeps point 397 at its first rezone (`FLAG-POINTS 397 969`, against `FLAG-POINTS 969` without the step).

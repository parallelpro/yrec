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

`main.f90` is the command-line entry point; `run/yrec_capi.f90` is the same engine for C and Python (`pyyrec`).

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| **run/** | |
| `run_driver` | `run_yrec`: the whole run. Reads the controls, sets up the tables, and goes through the run cards, evolving each card's model with `evolve_step`; the calibration verdict at the end of each card. |
| `run_setup` | `set_up_constants_and_tables` (the physical constants, and the kap/eos/atm table loads through their `*_init`) and `prepare_run` (that, then the Monte-Carlo run range). |
| `calibration` | The solar and star calibrations: their setup before the run (`setup_*_calibration`) and the check after each card (`check_*_calibration`). |
| `monte_carlo` | The Monte-Carlo runs: the sampled parameters and the summaries. |
| `stop_conditions` | The stops: end age, the central-abundance table, the star-calibration target radius. |
| `yrec_capi`, `yrec_reset` | The C entry point, and the reset that makes a second run in one process start from a fresh state. |
| `yrec_version` | The version stamp of the run log. |
| **step/** | |
| `step_driver` | `evolve_step`: one model (§2). |
| `mass_change` | Mass loss or accretion at the surface: the rate, the timestep limits, the entropy shift; gathers the convection zone from the star, calls the kernels, and stores the results. |
| `mass_change_kernel` | The kernels of the mass change, with strict in/out arguments and no star: the convection zone's mass and new masses (`convection_zone_mass`, `change_shell_masses`), the angular momentum a thermal wind carries off (`wind_specific_angular_momentum`, `remove_wind_angular_momentum`), and the accreted composition mixed into the zone (`mix_accreted_composition`). |
| `composition_step` | The operator-split composition change over the step at fixed structure: burning, convective mixing, semiconvection, settling. |
| `burn_step` | The star-layer side of nuclear burning: reads the model, calls net/, stores the results; the light elements. |
| `shell_physics` | `compute_shell_physics`: every zone's EOS, opacity, gradients and viscosity, for rotation and the convection-zone edges. |
| `convection_geometry` | Where the convection zones end (the base of the surface mixed region, pressure scale heights). |
| `timestep_control` | `compute_timestep`: the next timestep, the most stringent of the structure, ω, burning and HR-diagram limits. |
| **structure/** | |
| `henyey_iteration` | `henyey_iterate`: the Newton iterations of one solve level. |
| `structure_equations` | `compute_henyey_coefficients`: across the model, gathers each shell from the star, calls the kernel, eliminates, and stores the shells' output values. |
| `shell_equations` | `evaluate_shell_equations`: the kernel -- the four equations' terms and partials at one shell (eos, opacity, gradients, nuclear and gravitational energy), with strict in/out arguments and no star. |
| `henyey_elimination` | The block elimination (`henyey_eliminate`) and back-substitution (`henyey_solve`). |
| `surface_boundary` | `compute_surface_boundary`: the outer boundary condition from the triangle of envelope integrations: switches the atmosphere, keeps the triangle around the model, integrates the envelopes at its vertices, and stores the fit. |
| `triangle_fit` | The kernels of the triangle -- `start_triangle`, the point-in-triangle test (`edge_side`, `vertex_to_reflect`), `reflect_vertex`, `reflect_until_inside`, and `fit_triangle_coefficients` (log R, log L and log Teff as planes in the fitting point's log P and log T) -- with strict in/out arguments and no star. |
| `point_scratch_lib` | The per-point physics the integrands and `compute_henyey_coefficients` share. |
| **envelope/** | |
| `envelope_integration`, `envelope_integration_kernel` | The atmosphere and envelope integrations from the surface to the fitting point (the kernel is star-blind). |
| `envelope_integrands` | Their right-hand sides, in optical depth and in pressure. |
| `envelope_refit`, `envelope_rebuild` | Refitting the envelope, and moving the fitting point outward when the star accretes (the envelope above it thickens). |
| `envelope_structure`, `atmosphere_structure` | The outer layers of the current model, point by point. |
| **model/** | |
| `starting_model` | `read_starting_model`: each card's starting model, read or reused, rescaled if asked. |
| `model_rescaling` | `rescale_model`: a new mass, composition or core. |
| `rezoning` | `rezone`: the new mesh of each step and the model on it. |
| `shell_boundaries` | `locate_shell_boundaries`: the convective core, the surface convection zone, the hydrogen-burning shell. |
| **diagnostics/** | |
| `stitched_model` | The whole converged star (interior, envelope, atmosphere) and every quantity the output needs. |
| `observables` | The per-model observables. |
| `turnover_timescale` | The convective turnover time near the base of the surface convection zone. |
| `neutrino_flux_output` | The opt-in neutrino-flux table of the starting model. |

## 2. A step

`evolve_step` (`step/step_driver.f90`) advances the model by one timestep:

1. **Mass and composition** (`advance_composition_and_age`): the mass change, then the composition step over the timestep.
2. **The mesh** (`rezone_or_snapshot`): the rezone, and the start-of-step snapshot.
3. **The structure** (`solve_structure`): the Henyey solve in levels. Level 1 uses the previous model's envelope triangle; level 2 checks the surface boundary first; level 3 follows, and from level 3 on the composition step is redone before each solve (`solve_level`). Each level is `henyey_iterate`'s Newton iterations on the equations of `compute_henyey_coefficients`.
4. **Rotation** (`converge_with_rotation`): with rotation on, level 4 inside the structure⇄rotation iteration, which evolves the angular momentum (`evolve_angular_momentum`) and the shape factors.
5. **After convergence:** the shell boundaries, the light elements, the envelope rebuild after accretion, the next timestep (`compute_timestep`), the stitched model, the observables, the output, and the stop conditions.

A step that diverges stops the run: nothing cuts the timestep.

## 3. Conventions

- CGS unless a name carries its unit (`_msun`, `_lsun`, `_yr`, `_gyr`); `log` is log₁₀.
- The structure variables are log₁₀ P, log₁₀ T, log₁₀ R and L in L☉, against log₁₀ m. The difference equations and their partials follow `docs/TARGET_ARCHITECTURE.md` §4.4 (and `structure_equations`' header).
- A negative `star%dt` means "solve, do not age" (a rescale, a relaxation); it is made positive after convergence.
- Errors return `ierr` up to `run_yrec`, which stops the run.

## 4. Tests

Every deck runs core, so the deck pins are its test: the 37 pins of `test_all.py`/`test_mesa_baselines.py`, the legacy byte comparison of every example's log, and `test_reentry` (two runs in one process). No pinned deck rescales a model, calibrates to a target star, or writes the neutrino table; check changes there by running such a deck against the previous commit's binary.

The kernels have unit tests of their own (`docs/GUIDELINES.md`, star rule 6):

- `test_structure` (`core/test/test_structure.f90`, run by `test_structure.py`) compares every partial of `evaluate_shell_equations` with central differences of its own terms at six shells of a TAMS model (OPAL tables) and four of a birthline model (Yale EOS). The partials it marks `known` are documented approximations: the convective gradient's dependence on L is not linearized, and its P and T derivatives near the adiabat rest on the EOS's derivatives of grad_ad.
- `test_triangle_fit` (`core/test/test_triangle_fit.f90`, run by `test_triangle_fit.py`) checks the triangle's kernels on hand-made triangles: the fitted planes reproduce linear functions of (log Teff, log L) at the vertices and inside, for log R, log L and the reverse relation log Teff; degenerate fitting points (collinear, or vertices 1 and 2 at the same log T) give non-finite coefficients, unchecked; the point-in-triangle test inside, beyond, on the edges and at the vertices, both orientations; the reflection is the half turn about the opposite edge's midpoint; and walks of `reflect_until_inside` end with the point inside (they would never end for a triangle of zero area, or one whose orientation is the opposite of its vertices', as from a zero `tri_delta_teffl` or `tri_delta_logl`, or one of them negative).
- `test_mass_change` (`core/test/test_mass_change.f90`, run by `test_mass_change.py`) checks the mass-change kernels' conservation on a synthetic 40-shell model, with a partial and a fully convective envelope: after a loss or gain the shell masses sum to the new total and the zone's masses scale about its base (dm(1) from the center when the zone reaches it); a thermal wind leaves the zone J_cz + dJ for walpcz = 0, -1 and -2; accretion conserves every basic species.

# rotation/ — rotating stars

This directory evolves a star's rotation over each timestep:
- the shape of the rotating star, which enters the structure equations (fp, ft);
- the rotation curve implied by the angular momentum;
- the transport of angular momentum by secular instabilities (diffusion, and advection-diffusion);
- the mixing of the composition that transport drives;
- the wind's spin-down of the surface convection zone (magnetic braking).

Rotation has no facade: by the author's decision, core and setup use its modules by name. `tools/check_boundaries.py` lists the six procedures they may call: `evolve_angular_momentum`, `omega_from_j`, `rotation_shape_factors`, `zone_moments_of_inertia`, `am_convective_regions` and `viscos`.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `evolve_rotation` | `evolve_angular_momentum`, the step's driver (§2). |
| `mid_timestep_model` | One rotation sub-step: `advance_rotation_substep` interpolates the structure in time, and `burn_and_settle` burns and settles helium over the sub-step. |
| `omega_profile` | The rotation curve from the angular momentum: `omega_from_j`, with convection zones as solid bodies (`enforce_rotation_profile`) and radiative shells one at a time (`solid_body_omega`). |
| `am_regions` | `am_convective_regions`: the convection zones, with their overshoot, for angular momentum transport. |
| `viscosity` | `viscos`: the microscopic viscosity and the thermal diffusivity. |
| `magnetic_braking` | The wind torque on the surface convection zone, `brake_convection_zone`, and its solid-body form. Covers the Kawaler (1988) and Matt/PMM laws (K97, V13, or custom), with saturation, Rossby and field scaling. |
| `rotation_scratch_lib` | The working state: `rot_scr` (coefficients, grids, the structure at the middle of the sub-step), `mid_scr` (the physics interpolated to the sub-step) and `circ_scr`. |
| **shape/** | |
| `rotating_shape` | `rotation_shape_factors`: each shell's fp, ft, distortion, mean radius and mean gravity (Kippenhahn & Thomas / Endal & Sofia). |
| `moments_of_inertia` | `zone_moments_of_inertia`: I and dI/dω per zone (Law 1980). |
| `shell_quadrature` | The two quadratures of the distortion integral. |
| **seculr/** | |
| `secular_transport` | `apply_secular_transport`: one sub-step of transport, and the tridiagonal angular momentum solve. |
| `am_transport_coefficients` | The diffusion and advection-diffusion coefficients, and the banded solver. |
| `circulation` | The circulation and instability velocities, their scales, the Zahn coupling factor, the quadrupole moment. |
| `rotational_mixing` | The mixing of the composition: its grid, coefficients and solve. |
| `am_grid`, `equal_grid` | The equally spaced grids and the maps back to the model. |
| `transport_checks` | The angular momentum and composition checks, which can ask for a redo. |

`burn_settle_mix.f90`, `burn_mix_extrapolated.f90` and `burn_li_be_midpoint.f90` hold the Bulirsch–Stoer burn-and-mix path. It is commented out by the author's decision; each file's header says how to revive it.

## 2. A step

`evolve_angular_momentum` (core's `evolve_step` calls it when `rotation_active`):
1. Finds the regions for angular momentum transport.
2. **Without instability transport** (`instability_transport_active` false): the rotation curve, plus the wind on the surface convection zone over the whole step.
3. **With it:**
   - computes the rotation curve, the number of sub-steps (the ω-change limit) and the shell physics (`shell_physics` in core);
   - then, for each sub-step:
     - carry θ forward;
     - advance the mid-sub-step model;
     - apply the wind (fully convective models);
     - apply the secular transport; a redo repeats the sub-step from its start;
     - store the composition at the sub-step's end, and burn the light elements;
   - finally sets ω from the new angular momentum.

Rotation calls back into core (`shell_physics`, burn_step's light-element burn) and setup (`rotation_stability_setup`). Element diffusion is its own domain, reached through `element_diffusion_lib`.

## 3. Conventions

- CGS, and log is log₁₀. ω is in rad/s, and `j_rot` is the specific angular momentum. Zones run from the center out.
- Convection zones rotate as solid bodies unless an imposed law (`walpcz`) says otherwise. `force_solid_body_rotation` (or `solid_body_mode_flag = 1`) makes the whole star one.
- Failures (a singular solve, a step cut too often) return `ierr`.

## 4. Tests and extending

`test_rotation` (`make tests` in `src/`, then `pytest test_rotation.py`) covers:
- magnetic braking, every law and branch;
- the shape quadratures;
- the angular momentum solve;
- the rotating structure on a synthetic model: the moments of inertia, ω from J, the shape factors and the viscosity.

The decks pin the full step: the rotating solar deck, the five rotating configuration-matrix decks, and the three testsuite rotating-mixing decks (differential and solid body). A deliberate change means regenerating `rotation/test/expected_test_rotation.out`.

**Add a braking law:** give it a control, set its coefficients in `magnetic_braking`'s law setup, and add its cases to `test_rotation`'s braking checks.

# rotation/ — rotating stars

This directory evolves a star's rotation over each timestep:
- the shape of the rotating star, which enters the structure equations (fp, ft);
- the rotation curve implied by the angular momentum;
- the transport of angular momentum by secular instabilities, as diffusion;
- the mixing of the composition that transport drives;
- the wind's spin-down of the surface convection zone (magnetic braking).

Rotation has no facade: by the author's decision, core and setup use its modules by name. `tools/check_boundaries.py` lists the seven procedures they may call: `evolve_angular_momentum`, `omega_from_j`, `rotation_shape_factors`, `zone_moments_of_inertia`, `am_convective_regions`, `viscos` and `set_rotation_settings`.

Rotation reads its controls, the solar units and log10 G from a settings record, `rotation_settings` (module `rotation_configuration`), never from `s%ctrl`, `s%job` or physical_constants' run-time `cgl` (the checker enforces both). The star layer builds the record from the controls (`setup/rotation_settings_setup`), and installs it at the start of every run card (`setup/settings_install`'s `install_card_settings`, from `run_yrec`'s `begin_kind_card`).

The star comes in as an argument: the drivers (`evolve_angular_momentum`; the sub-step's `advance_rotation_substep` and `burn_and_settle`; `apply_secular_transport`; the checks `check_angular_momentum` and `check_composition`) take `type(star_info) :: s`, as core's procedures do, and nothing in rotation names a module-level star. Everything else is star-free: the drivers pass what it needs as arguments (the wind laws, for one, get the step's turnover times and photospheric pressures in a `wind_conditions` record). The kernels that touch no scratch either are listed in `check_boundaries.py`'s `STAR_BLIND_FILES` and tested by `test_rotation_kernels` (§4). `rotation_scratch_lib` is rotation's alone: no file outside `rotation/` may use it.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `evolve_rotation` | `evolve_angular_momentum`, the step's driver (§2). |
| `mid_timestep_model` | One rotation sub-step: `advance_rotation_substep` interpolates the structure in time, and `burn_and_settle` burns and settles helium over the sub-step. |
| `omega_profile` | The rotation curve from the angular momentum: `omega_from_j`, with convection zones as solid bodies (`enforce_rotation_profile`) and radiative shells one at a time (`solid_body_omega`). |
| `am_regions` | `am_convective_regions`: the convection zones, with their overshoot, for angular momentum transport. |
| `viscosity` | `viscos`: the microscopic viscosity and the thermal diffusivity. |
| `magnetic_braking` | The wind torque on the surface convection zone, `brake_convection_zone`, and its solid-body form; `outer_point_wind_terms`, the wind's boundary terms for a transport region that reaches the surface; `derive_wind_law`, the law's constant, exponents and saturation from its controls (a pure function `rotation_settings_from_star` calls). Covers the Kawaler (1988) and Matt/PMM laws (K97, V13, or custom), with saturation, Rossby and field scaling. |
| `rotation_messages` | `rotation_report` and `rotation_note`, rotation's two kinds of message (§3), and `rotation_numerics_report`, a numerics routine's failure as a one-line report (where, and `numerics_status_text` of its code). |
| `rotation_configuration` | The settings record `rotation_settings` (`rotation_settings_t`): every control rotation reads, with readable names (`cz_rotation_law_exponent` for `walpcz`, ...), and the wind law's coefficients from `derive_wind_law`. |
| `rotation_scratch_lib` | The types of the working state, one record per stage: `sub_model` (the structure interpolated to the sub-step, its convection zones, composition snapshot and starting ω), `sub_physics` (the physics interpolated to it), `stab` (the structure at the interfaces and the stability coefficients), `chi_grid` (the transport region's equally spaced grid) and `circ_scr` (the mechanisms' previous velocities and the Zahn coupling factor). `evolve_angular_momentum` allocates them for the step and passes each routine the records it uses; the module declares no instance. Nothing in them outlives a step; what does (the last sub-step's omega, `star%omega_substep`) is on the star. |
| **shape/** | |
| `rotating_shape` | `rotation_shape_factors`: each shell's fp, ft, distortion, mean radius and mean gravity (Kippenhahn & Thomas / Endal & Sofia). |
| `moments_of_inertia` | `zone_moments_of_inertia`: I and dI/dω per zone (Law 1980). |
| `shell_quadrature` | The two quadratures of the distortion integral. |
| **seculr/** | |
| `instability_velocities` | The velocity of each transport mechanism, as star-blind kernels: the Eddington–Sweet circulation and its μ-gradient inhibition, the GSF instability, the secular and dynamical shear, and the Solberg–Høiland instability. |
| `transport_mechanisms` | `compute_mechanism_velocity`: one slot's mechanism, its kernel called with what it reads from the interface record (star-blind). The slots, with each mechanism's name, the slot it corrects and its two controls, are `state/transport_mechanism_layout`'s. |
| `transport_regions` | `next_transport_region`: the runs of interfaces where the total velocity U is positive. |
| `region_transport` | `diffuse_region_angular_momentum`: one region's angular momentum diffusion over a sub-step, on its grid: the system, its solve and the map back to j/m (star-blind). |
| `secular_transport` | `apply_secular_transport`: one sub-step of transport, the iteration; its steps (the interfaces' setup, the velocities, the restart, the uncoupled spin-down, the region loop, the checks, the mixing) are module procedures with explicit arguments. |
| `am_transport_coefficients` | The angular momentum diffusion equation's tridiagonal system (`am_diffusion_coeffs`) and its solve (`solve_angular_momentum_transport`). |
| `circulation` | `circulation_velocities`: the structure at the interfaces, the mechanisms' velocities (`instability_velocities`) with their scales and averaging, their sum U; the diffusion coefficients from them (`diffusion_velocity_scales`); the Zahn coupling factor; the quadrupole moment. |
| `rotational_mixing` | The mixing of the composition: its grid, coefficients and solve. |
| `am_grid`, `equal_grid` | The equally spaced grids and the maps back to the model. |
| `stability_coefficients` | `rotation_stability_setup`: the structure at the zone interfaces, the zones' radii, the dynamical-shear limit and the coefficients of the circulation and instability velocities, once per sub-step (into `stab`). |
| `transport_checks` | The angular momentum and composition checks, which can ask for a redo. |

The commented-out code, kept by the author's decision, lives in `retired/` (not compiled): the Bulirsch–Stoer burn-and-mix path (`burn_settle_mix.f90.txt`, `burn_mix_extrapolated.f90.txt`, `burn_li_be_midpoint.f90.txt`, and its call sites in `bulirsch_stoer_calls.f90.txt`) and the advection-diffusion transport (`advection_diffusion.f90.txt`, every fragment with the place it came from). Each file's header says how to revive it.

The advection-diffusion transport (`use_diffusion_advection_transport`, legacy `LDIFAD`; Maeder & Zahn 1998) is commented out too (2026-10-12, the author's decision): it never worked, because most of its coefficient arrays were never filled. its banner in `retired/advection_diffusion.f90.txt` lists the pieces to restore; a rotating run that sets the flag stops.

## 2. A step

`evolve_angular_momentum` (core's `evolve_step` calls it when `rotation_active`):
1. Finds the regions for angular momentum transport.
2. **Without instability transport** (`instability_transport_active` false): the rotation curve, plus the wind on the surface convection zone over the whole step.
3. **With it:**
   - computes the rotation curve, the number of sub-steps (the ω-change limit) and the shell physics (core's `shell_physics`);
   - then, for each sub-step:
     - record the sub-step's starting ω (for the wind);
     - advance the sub-step's model (interpolated in time to the sub-step's end);
     - apply the wind (fully convective models);
     - apply the secular transport; a redo repeats the sub-step from its start;
     - store the composition at the sub-step's end, and burn the light elements;
   - finally sets ω from the new angular momentum.

Rotation calls back into core for two services: `shell_physics` (the per-shell EOS, opacity, gradients and viscosity, which every run also needs at its start and at rezoning) and burn_step's light-element burn. Element diffusion is its own domain, reached through `element_diffusion_lib`.

## 3. Conventions

- CGS, and log is log₁₀. ω is in rad/s, and `j_rot` is the specific angular momentum. Zones run from the center out.
- Convection zones rotate as solid bodies unless an imposed law (`cz_rotation_law_exponent`, the control `walpcz`) says otherwise. `force_solid_body_rotation` (or `solid_body_mode_flag = 1`) makes the whole star one; test it with `solid_body_everywhere()` (`rotation_configuration`).
- Composition slots come from `composition_layout`: the transport mixes slots 1 to `last_mu_species` (X, Y, Z, He3, which set the μ gradient) at every iteration, the rest up to `num_tracked_species()` once.
- Failures (a singular solve, a spline that cannot be evaluated, a step cut too often) return `ierr`; a numerics routine's code is passed up as it is, after a `rotation_numerics_report`. Messages go through `rotation_messages`: an error or a cut of the sub-step to the terminal and the run log (`rotation_report`), a per-sub-step note (the dynamical-shear interfaces, the transport's convergence) to the terminal, and to the run log too when solver diagnostics are on (`rotation_note`). Kernels print nothing.
- Kept by decision (2026-10-12): the μ gradient enters the Eddington–Sweet inhibition twice (`inhibit_by_mu_gradient`), so it cuts U_ES by at most a half, as in the F77 code. `docs/HOW_YREC_WORKS.md` §4.6 explains it.

## 4. Tests and extending

`test_rotation` (`make tests` in `src/`, then `pytest test_rotation.py`) covers:
- magnetic braking, every law and branch;
- the shape quadratures;
- the angular momentum solve;
- the rotating structure on a synthetic model: the moments of inertia, ω from J, the shape factors and the viscosity.

`test_rotation_kernels` (`pytest test_rotation_kernels.py`) checks the kernels against the laws they keep, on a synthetic model:
- the moments of inertia and the shape factors in the spherical limit, and their growth with ω (1 − fp as ω², 1 − ft as ω⁴);
- `omega_from_j` gives each radiative zone and the convective envelope back their angular momentum, the envelope as a solid body;
- the angular momentum diffusion conserves J without a wind, loses exactly the wind's surface term with one, and keeps a solid body;
- the composition mixing conserves each species and keeps a uniform composition;
- the Kawaler law's ω exponents below and above saturation;
- the wind law's coefficients (`derive_wind_law`): K97 has the Kawaler law's R and M exponents at any index and ω exponent 3, which the Kawaler law has at index 1.5; a Matt law's torque at the solar ω is fk times the solar torque; `solar_omega_sat` below 1000 is in solar units; K97 and V13 set the field scaling;
- the wind's spin-down against the analytic laws (exponential when saturated, Ω₀(1 + 2KΩ₀²t)^(-1/2) below), with the error falling as 1/N² in N calls, the implicit midpoint rule's order;
- the angular momentum diffusion against its eigenmodes: the cosine modes are the implicit step's eigenvectors exactly, and the decay converges to the continuum's exp(−κ(π/L)²t) as h² in space and dt in time;
- the instability velocities: U_ES and U_GSF ~ ω², no shear term and no GSF for a solid body, the μ-gradient factor between 1/2 and 1, no shear velocity for a solid body or below the secular threshold, a quadratic shear onset above it, and a μ gradient's stabilization; the Solberg–Høiland criterion N² + κ² ≥ 0 on single interfaces (a solid body and a uniform j stable, Ω ~ r⁻⁴ unstable with its velocity, a μ gradient's stabilization, Ω reversed, convective pairs skipped) and, on a grid of N_T, N_μ, Ω and dΩ/d ln r, that it is unstable only where the shear is dynamical, with U_SH ≤ U_DSI;
- the region scan on a synthetic U profile;
- the transport grid and the map back: equally spaced chi points, grid masses that add up to the region's, J conserved in a redistribution and a surface loss taken within 1e-5, the convection zone below the region taking its first change, a spline failure returning osplin's code with j unchanged;
- one region's transport (`diffuse_region_angular_momentum`): J conserved without a wind while j moves, the wind's surface term taken within 1e-5, a solid body unchanged; and the wind's terms for it (`outer_point_wind_terms`) in each case: below the surface, disk-locked, braked, no torque;
- the stability coefficients against their formulas, on a structure for which the interpolation to the interfaces is exact;
- the checks' contracts: convergence, ω from j and the stop on a negative j (`check_angular_momentum`); a small negative abundance repaired, a larger one a halved sub-step, the fourth cut a stop, and μ of a fully ionized gas (`check_composition`).

Each of these checks was confirmed to fail under a mutation of the code it covers.

`docs/ROTATION_VALIDATION.md` collects the validation: the analytic checks' measured orders, a sub-step convergence study (`src/tools/substep_convergence.py`), and a comparison with MESA for a solid body under J conservation (`sample_cases/rotation_vs_MESA`).

`python3 src/tools/coverage.py --work DIR --deck ...` measures the line coverage of rotation/ under the two test programs and the decks given (a -O0 --coverage build in DIR). With the configuration-matrix rotating decks, the rotating-variant decks (radiative cores, 50 Myr) and the rescale-only deck it is 96.9% of 2168 lines (2026-10-29, with the Solberg–Høiland kernel; 92.7% before the checks' unreachable branches and the gradient-reversal safeguard were removed, O15).

The decks pin the full step: the rotating solar deck, the five rotating configuration-matrix decks, and the three testsuite rotating-mixing decks (differential and solid body). A deliberate change means regenerating `rotation/test/expected_test_rotation.out`.

**Add a control:** add the field to `rotation_settings_t`, copy it from `star%ctrl` in `rotation_settings_from_star`, and read `rotation_settings%<field>`.

**Add a transport mechanism** (a new instability or circulation that transports angular momentum and mixes as a diffusion):
1. Declare it in `state/transport_mechanism_layout.f90`: the next slot (`i_mech_<name>`, `n_transport_mechanisms` + 1), its name, the slot it corrects (0 for none; an earlier one) and the names of its two controls. The coefficients sum the mechanisms in slot order, so the next slot leaves the present ones' sums unchanged bit for bit.
2. Write its velocity as a kernel in `seculr/instability_velocities.f90`: arguments in, its velocity at interfaces `zone_min..zone_max` out, star-free. Give it a check in `test_rotation_kernels` that a mistake would break (a limit, a scaling, a stability criterion).
3. Add its case to `compute_mechanism_velocity` (`seculr/transport_mechanisms.f90`): the call of its kernel with what it reads from `stab`, `sub_physics` and ω, into `mechanism_velocity(i_mech_<name>,:)`. Clearing, scaling, averaging between iterations, summing into U and the diffusion coefficients loop over the slots already; `test_rotation_kernels` checks that it writes no slot it has not declared.
4. Add its two controls to the registry (`defaults/controls_registry.tsv`, then the generators) and copy them into `rotation_settings%velocity_scale(i_mech_<name>)` and `%mixing_scale(i_mech_<name>)` in `rotation_settings_from_star`; `test_rotation_kernels.py` checks both against the names the layout declares.
5. For output, add a profile column that reads `s%mechanism_velocity(i_mech_<name>,k)` (`core/diagnostics/stitched_model.f90`, `io/profile_output.f90`).

**Add any other kernel:** keep it star-free (arguments in, results out), add its file to `STAR_BLIND_FILES`, and give it a check in `test_rotation_kernels` that a mistake would break (a conservation law, a limit, a scaling).

**Add a braking law:** give it a control (as above), set its coefficients in `magnetic_braking`'s `derive_wind_law` (with a check in `test_rotation_kernels`' `wind_coefficients_check`), and add its cases to `test_rotation`'s braking checks (which set `rotation_settings` directly).

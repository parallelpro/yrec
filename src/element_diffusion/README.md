# element_diffusion/ — gravitational settling and microscopic diffusion

This directory moves the elements through the radiative interior over a timestep: the gravitational settling and microscopic diffusion of hydrogen and helium, the metals and the light elements Li6, Li7 and Be9, in the formalism of Bahcall & Loeb (1990) with the diffusion coefficients of Thoul, Bahcall & Loeb (1994). Each species obeys

    dX/dt = (1/(ρ r²)) d/dr [ D1 + D2 dX/dr ],

a drift term D1 (the pressure, temperature and hydrogen-concentration gradients) and a diffusion term D2, solved on equally spaced grids in radius between any convective core and the surface convection zone.

Code outside element_diffusion/ uses it only through `element_diffusion_lib` (enforced by `tools/check_boundaries.py`, `FACADE_ONLY_DOMAINS`). Like the microphysics domains it never touches the star (`STAR_FREE_DIRS`): its controls come from a settings record, `diffusion_settings` (module `element_diffusion_configuration`), which the star layer builds (`setup/diffusion_settings_builder`'s `diffusion_settings_from_star`) and installs at the start of every run card (`set_diffusion_settings`, in `setup/settings_install`'s `install_card_settings`). The one value that changes within a card comes in through `diffuse_elements`: the run log's latch (`suspension_reported`, the star's `settling_suspended_reported`).

There are two pipelines. **microdiff** (`use_new_diffusion_routines`) is the one every shipped deck runs (354 of the 515 shipped decks settle, all with it). The **original routines** (`gravitational_settling`) are selected by no shipped deck; whether to fix and test them or retire them to `retired/` is the authors' decision (`audit/module-evaluation-7-2026-11-10/B_LISTS.md`, element_diffusion 1). Both share the region, the grids, the Coulomb logarithms, the scheme, the Thoul solve and the messages.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `element_diffusion_lib` | **The public interface.** `diffuse_elements`, one settling step: the region (`find_settling_region`), its note, then microdiff or the original routines; `plan_settling_substeps` and its `settling_substep_plan`, the sub-steps the callers take; `set_diffusion_settings` and `diffusion_settings_t`. |
| `element_diffusion_configuration` | The settings record `diffusion_settings` (`diffusion_settings_t`): every control element diffusion reads (§3 says which pipeline each acts in). |
| `element_diffusion_messages` | Every text the domain writes, with its destination (§3): `note_settling_region`, `note_solve_pass`, `report_unconverged`, `report_numerics_failure`. |
| `settling_region` | `find_settling_region`: the zones that settle, and the verdict when there are none (`region_found`, `region_fully_convective`, `region_no_hydrogen`, `region_no_helium`). |
| `settling_substeps` | `plan_settling_substeps`: the sub-steps of a step, none longer than `settling_timestep_fraction` of helium's settling timescale at the outer boundary (Michaud et al. 1984), and their count. |
| `microdiff_driver` | microdiff: `microdiff` (the stages of §2), `convert_to_bahcall_loeb`, and `report_species_solve` (a solve's messages). |
| `microdiff_mapping` | microdiff's grids: the `microdiff_grid` record (the structure and the mass fractions at one grid's points), `map_microdiff_to_equal_grid`, and `microdiff_etm`, the changes back on the model (through `apply_grid_changes`). |
| `microdiff_run` | `run_microdiff`, one species' step (the drift by Lax-Wendroff, the diffusion term by `solve_diffusion_term`), and the `species_solve_report` it returns. |
| `microdiff_coefficients` | `compute_microdiff_coefficients`: D1 and D2 at one grid's points from the Thoul solve, with the hydrogen gradient. |
| `diffusion_grid` | The equal grids and their maps, for both pipelines: `size_equal_grid`, `place_equal_grids`, `build_equal_grid_map` (the `equal_grid_map` record) and `map_field` onto the grids; `apply_grid_changes` and `change_zone`, the changes back on the model's zones; `locate_stencil`, the 4-point search both directions use. |
| `diffusion_scheme` | The scheme both pipelines use: `lax_wendroff_step1`, `lax_wendroff_step2` and `implicit_diffusion_coeffs`. |
| `coulomb_logarithm` | `compute_concentrations` and `compute_coulomb_logarithms`: the pairwise screened Coulomb logarithms, for microdiff and the original routines' `coulomb_log_choice` 4. |
| `thoul_solver` | `thoul_diffusion`: Burgers' equations solved for the diffusion coefficients (A. A. Thoul's routine). |
| `bahcall_loeb_units` | The four unit scales of the Bahcall & Loeb equations (radius, mass, temperature, time), compile-time parameters. |
| `species_table` | The species' weights and charges: the four Thoul columns (`thoul_col_*`; H, He, Fe, electrons) and the light elements as a table (`light_diffused`, their composition rows, `num_light_diffused`). |
| `original_settling` | The original routines: `gravitational_settling` (the stages), `iterate_hydrogen_diffusion` (the implicit step's iteration), `gravitational_settling_setup` (the units in place, and the coefficients: Thoul et al.'s, their fit, or the older analytic ones) and `coulomb_log_of_choice`. |
| `settling_grid` | The original routines' maps: `model_to_equal` and `equal_to_model` (the changes back, and the units back), both on `diffusion_grid`. |
| `element_diffusion_scratch` | The type of the original routines' working state (`diff_scr`: the metals' coefficients and changes). `diffuse_elements` allocates the record for their step and passes it down; the module declares no instance. |

## 2. The flow

The callers, each with its own outer boundary and convective flags:
- core's `composition_step` (`settle`) settles every model while the star's `settling_active` holds, with the deep-mixing flags (the convective flags, without overshoot);
- rotation's `burn_and_settle` (`settle_helium`, in `mid_timestep_model`) settles within each rotation sub-step, with the mixed regions and their overshoot.

Each caller (1) stops settling for the rest of the card once the central X falls below `hydrogen_diffusion_floor` (`settling_active` false), (2) needs a surface convection zone (core skips with a run-log warning, rotation stops the run), (3) takes its outer boundary (core: the first zone inward from the envelope's base with Y above `helium_diffusion_min`; rotation: the top of the outermost radiative region), (4) cuts the step with `plan_settling_substeps`, and (5) calls `diffuse_elements` once per sub-step.

`diffuse_elements`:
1. `find_settling_region`: from the top of any convective core (zone 1 without one), moved out past a hydrogen-exhausted core to the first zone with X above `hydrogen_diffusion_floor`, to the base of the surface convection zone, moved in to the first zone with Y above `helium_diffusion_min`.
2. `note_settling_region` (messages): with no region, the run-log note, once per suspension; with one after a suspension, "SETTLING RESUMED". No region, no settling.
3. microdiff (or the original routines, below).

`microdiff`:
1. `convert_to_bahcall_loeb`: the timestep, d ln P/dr, the enclosed mass and the total mass in Bahcall & Loeb units, on copies (the model's own stay exactly as they were); the radius in solar radii; the density and temperature unlogged.
2. `map_microdiff_to_equal_grid`: the grids over the region (`size_equal_grid`: the span over the narrowest of the region's outer 20 zones, rounded up, is the point count, so the intervals are a little wider than that zone), their points (`place_equal_grids`), each point's stencil (`build_equal_grid_map`), and the structure and the mass fractions on them (`map_field`).
3. `run_microdiff`, for hydrogen (iron as the heavy species; helium follows from it), then the metals (`use_diffusion_z`, as iron), then each light element (`diffuse_lithium`):
   - `compute_microdiff_coefficients` at the centres: each point's concentrations and Coulomb logarithms (`coulomb_logarithm`), Burgers' equations (`thoul_diffusion`), D1 and D2; the hydrogen call computes d ln C_H/dr, which the others read;
   - `lax_wendroff_step1`: the drift over half the step at the midpoints, and the coefficients there from the abundance it gives;
   - `lax_wendroff_step2`: the drift over the step at the centres;
   - `solve_diffusion_term`: the diffusion term, one tridiagonal solve (`implicit_diffusion_coeffs`, numerics' `solve_tridiagonal`);
   - `report_species_solve`: the solve's messages.
4. `microdiff_etm`: the changes back on the model's zones (`apply_grid_changes`, `change_zone`): the core's zones take the first centre's change, the envelope's the last centre's, the zones between the change interpolated at their radius; then Y = 1 − X − Z − He3.

`gravitational_settling` (the original routines): `gravitational_settling_setup` converts the timestep, d ln P/dr, the enclosed mass and the total mass to Bahcall & Loeb units in place and sets each zone's coefficients and their derivatives in X (the analytic ones, or Thoul et al.'s with `coulomb_log_choice`'s logarithm, or their fit); `model_to_equal` maps them to the grids; the two Lax-Wendroff steps for X and Z; an iteration of hydrogen's implicit step (`iterate_hydrogen_diffusion`; the metals have no diffusion term here) that corrects D2 by its derivative until the largest change is below `settling_tolerance`; `equal_to_model` maps the changes back and converts the units back (the round trip can change the last bit of the enclosed mass).

Calls into other packages: numerics (`interp`, `lagrange4`, `solve_tridiagonal`, `ludcmp`, `lubksb`, `numerics_status_text`), math (`pow`, `exp`, `log` under CRMATH), state (`composition_layout`, `physical_constants`, `run_messages`).

## 3. Conventions

- **Units.** CGS unless a name says otherwise, log is log₁₀. Inside the pipelines the radius, mass, timestep and d ln P/dr are in Bahcall & Loeb units (`_bl`; `bahcall_loeb_units`): the solar radius, 100 g/cm³ R☉³, 10⁷ K and 2.7 × 10¹³ yr, in which the coefficients carry ln Λ/2.2. Zones run from the centre out.
- **Errors.** A singular matrix (Burgers' system, `ludcmp`; the implicit step, `solve_tridiagonal`) returns numerics' status in `ierr`, after `report_numerics_failure` names the routine; the callers stop the run on it. An implicit solve that misses `settling_tolerance` writes `report_unconverged` and the run goes on.
- **Messages** go through `element_diffusion_messages`, the only file with an I/O statement; the kernels print nothing (`check_boundaries.py`: the six state-free kernels are in `STAR_BLIND_FILES`, the four that read the settings in `KERNEL_FILES`). The texts and their destinations are the F77's, kept until the B phase rewords them (`check_style.py` counts no I/O elsewhere: element_diffusion's ten style counts are all 0):
  - the run log: a model that cannot settle ("FULLY CONVECTIVE MODEL - NO SETTLING", "X BELOW ... IN WHOLE MODEL-NO SETTLING", "Y BELOW ..."), once per suspension and every model under solver diagnostics, and "SETTLING RESUMED";
  - the run log under solver diagnostics: each pass of an implicit solve ("ITERATION ... DXMAX ... IMAX");
  - the terminal and the run log: "MICRODIFF (GRSETT) FAILED TO CONVERGE ...";
  - standard output only: "<routine>: <status text>" for a singular matrix.
- **Controls** come from `diffusion_settings`, never from `s%control` or `s%physics`; the callers read `settling_timestep_fraction` from their own records and pass it to `plan_settling_substeps`. Which setting acts where:

  | Setting (control) | microdiff | original routines |
  |---|---|---|
  | `use_new_diffusion_routines` | selects microdiff | selects the original routines when false |
  | `use_diffusion_z` | the metals settle (as iron) | the metals settle (as iron) |
  | `diffuse_lithium` | Li6, Li7 and Be9 settle (as column 3) | not read |
  | `fgry` (`settling_y_factor`) | factor on hydrogen's (so helium's) coefficients | factor on the coefficients |
  | `fgrz` (`settling_z_factor`, Monte Carlo) | factor on the metals' **and the light elements'** coefficients (bug_sweep M250) | factor on the metals' coefficients |
  | `hydrogen_diffusion_floor`, `helium_diffusion_min` | the region (and the callers' stop and outer boundary) | the same |
  | `settling_num_iterations`, `settling_tolerance` | only the messages: the implicit step is solved once (M253) | the iteration's passes and tolerance |
  | `coulomb_log_choice` | not read: each pair's screened logarithm | 2 Noerdlinger's, 3 Loeb's, 4 each pair's screened logarithm, else 2.2 |
  | `use_thoul_diffusion`, `use_thoul_fit` | not read: always Thoul et al.'s solution | Thoul et al.'s solution, their fit, or (off) the older analytic coefficients |
  | `constant_mixing_coeff` | not read | with Thoul's coefficients, a factor on D1 and D2 (the metals' too) |
  | `constant_settling_reduction` | not read | with Thoul's coefficients, a factor on the metals' D1 |
  | `settling_timestep_fraction` | the callers' sub-steps | the same |

  `diffuse_lithium` brings five other settings with it when the controls are read (io's `derive_diffusion_switches`, bug_sweep M181): `use_diffusion_z`, `diffuse_helium_active`, `use_thoul_diffusion` and `use_new_diffusion_routines` on, `use_thoul_fit` off and `coulomb_log_choice` 4.
- **Kept by decision** (until the B phase): microdiff divides by ln Λ = 2.2 as the F77's single-precision literal (2.2000000476837158) where the original routines use 2.2d0 (M257); the first midpoint of both grids is searched over model zones 2..num_points as in the F77, which extrapolates from inside a convective core when the region starts at or beyond that zone (M247: the fix moves a pinned deck); each map back keeps its own guard for grids of fewer than four points (O21).

## 4. Tests

`test_element_diffusion` (`make tests` in `src/`, then `pytest test_element_diffusion.py` from the repository root; `YREC_TEST_TOLERANT=1` compares the numbers to a tolerance) runs one 10⁸-year step on a synthetic 200-zone solar-like model and prints the changes, then its checks.

The pinned cases (the expected output's first part):
- the original routines with helium only, with metals, each Coulomb logarithm, the mixing and settling factors, the scale factors, the Thoul fit, the older coefficients, and a solve stopped short of convergence;
- microdiff with helium only, with metals, and with the light elements;
- for both, the abundance floors (the outer edge where Y first exceeds `helium_diffusion_min`, and the notes when Y or X is below its floor everywhere).

The case "old: without the Thoul et al. coefficients" reads the original routines' stale metal coefficients (bug_sweep M33), so its pin holds one draw and differs in about 4 runs of 100 (lines 237-239 only; measured 2026-11-10); the checks below do not use it.

The checks (`# checks`), each printed `ok` or `FAIL` (a FAIL stops the program with status 1):
- microdiff on the model (helium, metals, light elements): X + Y + Z + He3 = 1 in every zone to 10⁻¹⁵; each species' mass conserved to 3 × 10⁻⁴ of the mass it moves, in cells centred on the zones (the scheme conserves exactly on its grid; the map back by interpolation keeps this much); helium, the metals, Li6, Li7 and Be9 sink and hydrogen rises at the surface, helium gathers at the centre; a zero timestep changes nothing, bit for bit; the change is linear in a small timestep (twice the step, twice the change, to 10⁻⁴ at 10⁴ yr); with a convective core the region starts at its top and the core's zones all take one change.
- `find_settling_region`: a core's top to the envelope's base, moved in to Y above its floor; the three verdicts without a region.
- `plan_settling_substeps`: Michaud's timescale and its scalings in q and T; the count is max(1, ceiling(dt/longest)), an exact multiple taking no extra sub-step (O16).
- `diffusion_grid`: the point count, the span over the narrowest zone rounded up, and the even spacing; `map_field` exact on a linear field at every point and on a cubic at the interior points; `locate_stencil`'s two points each side, within the table.
- `diffusion_scheme`: both Lax-Wendroff steps and the implicit step conserve the species' mass on the grid (flux form).
- `compute_coulomb_logarithms`: symmetric, and the H–H logarithm from k_B and e (the Debye length) to 10⁻⁴; the fit's limits, 0.97494 ln x − 1.3592 for weak coupling and 0.15249 x^1.2 for strong.
- `thoul_diffusion`: no current and no net mass flow for any gradient; the coefficients scale as 1/ln Λ; hydrogen rises, helium and iron sink.

Each check has a mutation in `test/test_element_diffusion.mutations.toml` that makes it print FAIL; `python3 src/tools/mutate.py element_diffusion` (after `make tests`) reports 21 of 21 mutations caught (2026-11-10). The program runs clean under the development flags (`make devel_tests`, `python3 src/tools/devel_tests.py`).

A published fit of Thoul et al.'s coefficients was not checked: the original routines' `use_thoul_fit` polynomials agree with the solution only at an effective ln Λ of 2.6–2.8 (2–4 %), not at the 2.2 they are used with (20 %), and the paper was not at hand; the checks above test what the solution must satisfy instead.

The decks pin microdiff: the solar GS decks, the settling and rotating cases of the configuration matrix (cm_settle_he, cm_settle_he_z, cm_rot_*), and the CF10 solar pair. A deliberate change of the pinned cases means regenerating `test/expected_test_element_diffusion.out` and saying why in the commit.

## 5. Coverage

97.1% of 829 lines (2026-11-10, at 5628d302), from `test_element_diffusion` and the settling decks cm_settle_he, cm_settle_he_z and cm_rot_settle_he (the configuration matrix) and rv_ref and rv_nocore (`src/tools/compare_decks/`, rotating with radiative cores, light elements on):

    python3 src/tools/coverage.py --work DIR --dir element_diffusion --test test_element_diffusion \
        --deck examples/run_config_matrix/cm_settle_he.inlist --deck examples/run_config_matrix/cm_settle_he_z.inlist \
        --deck examples/run_config_matrix/cm_rot_settle_he.inlist --deck src/tools/compare_decks/rv_ref.inlist \
        --deck src/tools/compare_decks/rv_nocore.inlist

Every file is above 90% but the messages module (80.5%). The 24 lines not run: the singular-matrix paths (`ludcmp`'s and `solve_tridiagonal`'s failures in microdiff_run and original_settling, and `report_numerics_failure`), the solver-diagnostics note (`note_solve_pass`: no deck or test turns diagnostics on), microdiff's unconverged report (its solve converges in its passes in every case), and the second midpoint of a two-point grid in `place_equal_grids`.

## 6. How to extend

The settings record is the only way a control reaches element diffusion (`STAR_FREE_DIRS`); the callers' sub-step policy is `plan_settling_substeps`, which takes what it needs as arguments.

**Add a control** (traced 2026-11-10 on a throwaway branch, with the treatment below):
1. The control itself, as `io/README.md` §5 ("Add a control") says: its declaration and description in `io/controls_buffer.f90`, and its row in `defaults/controls_registry.tsv` with the settings cell `diffusion_settings%<field>`. `make -j8` in `src/` rewrites the generated files (`defaults/controls_registry.md`, `defaults/physics.defaults`, `io/controls_buffer_public.inc`, `io/controls_sync.f90`, `io/inlist_decl.inc`, `io/test/controls_dump.inc`, `state/physics_group_def.inc` for a `physics` control); commit them.
2. The field in `diffusion_settings_t` (`element_diffusion_configuration.f90`), with what it does, and its copy in `diffusion_settings_from_star` (`setup/diffusion_settings_builder.f90`); then `python3 tools/check_star_contracts.py --fix` (the builder's `star reads:` line) and `python3 tools/check_boundaries.py` (it checks the settings cell against the builder).
3. Read `diffusion_settings%<field>` where it acts, and add its row to §3's table.
4. `YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py` from the repository root: `src/io/test/expected_controls/pristine.txt` gains the control's line; then `python3 -m pytest -q test_generated_files.py test_boundaries.py`.

**Add a settling treatment** (a pipeline beside microdiff, or a new set of coefficients):
1. Its module, `element_diffusion/<name>.f90`: star-free, arguments in and results out, over the region `diffuse_elements` hands it (`zone_begin`, `zone_end`); reuse `diffusion_grid`, `diffusion_scheme`, `coulomb_logarithm` and `thoul_solver`, and write any message through `element_diffusion_messages` (a kernel prints nothing).
2. Its control (above), and its branch in `diffuse_elements` (`element_diffusion_lib.f90`), after the region is found.
3. `python3 tools/gen_makefile_deps.py` in `src/` (deps.mk; `make check` fails until it is current), and the file in `tools/check_boundaries.py`'s `STAR_BLIND_FILES` (no module state) or `KERNEL_FILES` (it reads the settings record).
4. Its checks in `test/test_element_diffusion.f90` (a conservation law, a limit, a scaling) through `report_check`, each with a mutation in `test/test_element_diffusion.mutations.toml`; `make tests`, then append the new lines to `test/expected_test_element_diffusion.out` (the pinned cases above them must not move) and run `python3 src/tools/mutate.py element_diffusion`.
5. Its rows here: §1's table, §2's flow, §3's settings.

**Add a light element**: a row of `light_diffused` and of `light_diffused_row` (its weight, charge and composition slot) in `species_table.f90`, and `num_light_diffused` one more; microdiff's grids, its loop over the light elements and the map back take the count from there (bug_sweep M252). The slot must exist in `state/composition_layout.f90` (state's README) and in the starting models; then a check in `test_element_diffusion` that the new element sinks, and its mutation.

**Add a kernel**: as the treatment's steps 1, 3 and 4: star-free, in `STAR_BLIND_FILES` or `KERNEL_FILES`, with a check a mistake would break and its mutation.

# convection/ — convection physics

This directory decides where the star convects and what its temperature gradient is there:
- the gradient from mixing-length theory, with starspots;
- the convection zones from the shells' convective flags, extended by step overshoot and by semiconvection;
- the homogenization of the zones' composition;
- the pressure scale height and the envelope overshoot depth that core's light-element burn and rotation's sub-steps use.

Code outside convection/ uses it only through `convection_lib` (`tools/check_boundaries.py`, `FACADE_ONLY_DOMAINS`). convection/ reads nothing from the star (`STAR_FREE_DIRS`): its state comes in as arguments and its configuration as one settings record, `convection_settings`, which the star layer builds in setup/ and installs at the start of every run card. Its kernels print nothing: `temperature_gradients`, `convection_zones`, `overshoot` and `semiconvection` are in `STAR_BLIND_FILES` and `KERNEL_FILES`, and every message goes through `convection_messages` (§3).

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `convection_lib` | **The public interface.** It re-exports `compute_temperature_gradients`, `find_convection_zones`, `homogenize_convection_zones`, `extend_by_semiconvection`, `pressure_scale_height`, `envelope_overshoot_depth`, the settings type and the two error codes; `set_convection_settings` installs the settings, and `set_convection_envelope_metals` passes the starting model's envelope metal fraction outside C, N and O. |
| `convection_configuration` | The settings record `convection_settings_t` and the installed copy `convection_settings`; the error codes `convection_unstable_radicand` and `convection_too_many_zones`. |
| `temperature_gradients` | `compute_temperature_gradients` (formerly TPGRAD), as named steps: the radiative gradient (`radiative_gradient_of`), the adiabatic shortcut above `log10_t_cut` (`take_adiabatic`), the starspot correction (`spot_factor`, Somers 2014), the mixing-length coefficients and cubic (`mixing_length_coefficients`, `solve_mlt_cubic`), the convective velocity, and the gradient's log derivatives for the Henyey Jacobian (`mlt_derivatives`). The module header writes the equations, the derivatives included. |
| `convection_zones` | `find_convection_zones` (formerly CONVEC): the convection zones from the flags (`next_convective_run`, `list_convection_zones`), extended by overshoot and merged (`merge_overlapping_zones`), and the radiative regions between them (`list_radiative_regions`). `homogenize_convection_zones` (formerly MIXCZ) mixes each zone (`homogenize_region`) and floors the composition (`apply_composition_floors`). |
| `overshoot` | Step overshoot: `extend_by_overshoot` (formerly OVERSH) moves the edges of the mixed regions by α H_P (`extend_upward`, `extend_downward`, H_P at the edge with the EOS's density in `edge_scale_height`, formerly HSUBP); `pressure_scale_height`, H_P = P r²/(G m ρ); `envelope_overshoot_depth`, α_env H_P below the surface convection zone (H_P from `edge_scale_height`): the extent `extend_by_overshoot` moves that zone's edge by and the one its callers outside take. |
| `semiconvection` | `extend_by_semiconvection` (formerly SCONVEC): semiconvection at the upper edge of the convective core (Castellani et al. 1971), the edge in `extend_core_edge`, each shell's ∇_rad in `radiative_gradient_at`. |
| `convection_messages` | Every message of the directory, with its text and destination (§3). |

The settings builder is `setup/convection_settings_builder.f90` (`convection_settings_from_star`). core keeps the geometry of the zones' edges (`core/step/convection_geometry.f90`) and, in `composition_step`, its own homogenization of the mixed regions with overshoot (`homogenize_mixed_regions`, a byte-identical copy of `homogenize_region`, which `convection_lib` can export for core's round).

## 2. The flow

convection/ is called, never calls back into the star layer; it calls eos_lib (`eos_get`) and kap_lib (`kap_get`) for the EOS and the opacity it needs.

1. **Every run card:** `settings_install`'s `install_card_settings` (from run_driver's `begin_kind_card`) installs `convection_settings_from_star`; `install_surface_mixture_settings` sets the envelope's metals once the first card's model is read.
2. **Each shell's gradient** (`compute_temperature_gradients`), wherever a shell's physics is evaluated: core's `shell_physics` (every shell's gradient, convective flag and ∇_rad, at the start of a run, after a rezone and in rotation's sub-steps), `shell_equations` (each shell's structure equations and their derivatives in the Henyey iteration), `envelope_integrands` (the envelope and atmosphere integration, star-blind) and `starting_model` (the outermost shell's convection); within convection/, semiconvection's `radiative_gradient_at`. The steps are in `temperature_gradients`' header: ∇_rad; a radiative shell (∇ = ∇_rad, or ∇_ad in the overshoot region below the envelope with `adiabatic_overshoot`); a convective shell above `log10_t_cut` (∇ = ∇_ad); the spot correction; the mixing-length cubic; the velocity; the derivatives.
3. **Each composition step** (core's `advance_composition`, in the predictor pass and in each corrector pass of the iteration):
   1. `find_convection_zones` with the deep-mixing flags: the zones (runs of more than one convective shell), then, with an overshoot switch on, `extend_by_overshoot` and `merge_overlapping_zones`, then the radiative regions;
   2. core burns the composition (its rates and network);
   3. with `use_semiconvection`, on the passes that run it, `extend_by_semiconvection`: the convective core's upper edge (`extend_core_edge`; nothing without a mixed region, without a convective core or for a fully convective model), then `merge_overlapping_zones`;
   4. core homogenizes the mixed regions and settles the elements (element_diffusion);
   5. after deep mixing, `find_convection_zones` again with the ordinary flags.
4. **After the step converges:** core's `evolve_step` mixes the convection zones (`homogenize_convection_zones` with the convective flags); for the Li/Be burn it finds the end-of-step zones (`find_convection_zones`) and the overshoot depth below the surface convection zone at the start and the end of the step (`envelope_overshoot_depth`).
5. **Rotation** (rotating models): the regions for angular momentum transport (`am_convective_regions` → `find_convection_zones`, at the start, after a rezone and in each sub-step), the sub-step's mixing (`homogenize_convection_zones` with the transport's flags) and its Li/Be burn (`find_convection_zones`, then `envelope_overshoot_depth` when the surface is convective).

## 3. Conventions

- **Units:** CGS unless the name has a unit suffix; log is log₁₀. `luminosity_lsun` is L/L☉, linear. The derivatives are d ln ∇/d ln T, d ln P, d ln r.
- **Zone lists:** shells run from the centre (1) out. A zone list holds, in row k, the first and the last shell of the k-th region from the centre (`bounds(k,1)`, `bounds(k,2)`); `mixed_zone_bounds` has `max_convective_zones` (12) rows. `find_convection_zones` and `homogenize_convection_zones` set `convective_flag(num_zones+1)` to `.false.`, a sentinel callers' own scans read.
- **Errors:** a routine that can fail returns `ierr` (positive): `convection_unstable_radicand` (1) when the mixing-length square root would take a negative number, `convection_too_many_zones` (2) when a model has more convection zones than the bounds array has rows, or the EOS's or the opacity's code, passed up as it is. The callers in core and rotation pass it up.
- **Messages:** through `convection_messages`, each with its text and destination since the F77: to the run log, the zone limit, the mixed regions with and without overshoot (`extend_by_overshoot`), the two merges, the semiconvective edge and the cubic's non-convergence; to the terminal, the mixing-length square-root warning (before the `ierr` return). The run log is `run_messages`' `run_log_unit`, written directly. Kernels print nothing (`check_boundaries.py`'s `KERNEL_FILES`). The rewording, and the formats that overflow (bug_sweep M254), wait for the B phase's re-pin.
- **Controls:** read from `convection_settings` only, never from `s%control`/`s%physics` (the directory is star-free). A new control goes through the settings record (§6).
- **Kept by decision:**
  - semiconvection acts at the convective core only (the YREC team's core-only switch of 2025, a8e60674): nothing below an envelope or an intermediate zone (bug_sweep M241, the authors' choice of 2026-11-11);
  - the cubic's start A3^-0.333333333 (the F77's literal for -1/3; results depend on it);
  - the semiconvection search evaluates every shell at the first outside shell's radius (the F77's; the radius enters only the cubic);

## 4. Tests

`test_convection` (`make tests` in `src/`, then `pytest test_convection.py` from the repository root; plain, and with `YREC_TEST_TOLERANT=1` as CI runs it). It first prints, for byte-comparison, the gradient on hand-set EOS and opacity results (radiative points and the adiabatic-overshoot gate, a convective core above `log10_t_cut`, mixing-length convection with and without derivatives, both starspot modes, the error path), the zones of a synthetic 40-zone model with overshoot of every kind and a merge (Yale EOS scale heights), and homogenization. Then the ok/FAIL checks, each confirmed to fail under a mutation:
- **gradient**, on a synthetic EOS and opacity (power laws with exact derivatives):
  - ∇_rad = 10^cdelrl κ L P/(G m T⁴) f_T/f_P;
  - a radiative shell: ∇ = ∇_rad and v_c = 0, and ∇_ad in the overshoot region with `adiabatic_overshoot`;
  - above `log10_t_cut`: ∇ = ∇_ad and ∇_ad's derivatives;
  - the mixing-length ∇ and v_c against a solution computed in the test from the header's equations, its cubic by bisection, over 18 decades of density;
  - ∇_ad ≤ ∇ ≤ ∇_rad from ρ = 10⁻²¹ to 10³;
  - the limits' scalings: (∇ − ∇_ad)/Δ ∝ ρ^-1/3 (efficient) and (∇_rad − ∇)/Δ ∝ ρ² (inefficient);
  - the three derivatives against centred differences (relative error below 10⁻⁶);
  - with spots (constant and depth-varying contrast): ∇ and v_c against the computed solution, and the three derivatives against centred differences, the T and P ones with ∇_rad/S and d ln S/d ln T (bug_sweep M242);
  - the error contract: a negative radicand returns `convection_unstable_radicand`;
- **zones:** runs of one shell are not zones; the mixed and radiative regions tile the model; twelve zones fill the twelve rows; a thirteenth returns `convection_too_many_zones`;
- **overshoot:** `pressure_scale_height` = P r²/(G m ρ); `envelope_overshoot_depth` = α_env H_P with the EOS's density, and 0 when off; each kind of region's moved edge is the last shell within α H_P (H_P from the Yale EOS, computed in the test), the core's limited to β r; the overshoot region below the envelope, and none without a convection zone, for a fully convective model or from an intermediate zone's overshoot (bug_sweep M249); the rotating rule's one shell, below the envelope and above the core; the merge of regions that meet, their union when the lower one's overshoot reaches past the upper one's top (bug_sweep M285);
- **homogenization:** each species' mass in each zone is kept; each zone is uniform; the radiative shells are untouched; a uniform composition stays; the floors (X ≥ 0, Z ≤ 1 − X, O16 ≥ 0.99×10⁻³ of the CNO);
- **semiconvection**, on a 40-shell model with a convective core under a hydrogen gradient (X = 0.30 in the core, 0.45 to 0.70 above), with the OPAL95 opacities (`kap_init`) and the Yale EOS: the core's edge moves out two shells; f ∇_rad ≥ ∇_ad in the shells it adds and not in the next (one shell has ∇_rad > ∇_ad but f ∇_rad < ∇_ad); a 1 s timestep adds only the first shell (the propagation limit); a region the extension overlaps is merged, one that only shares the edge shell is not; a stable edge stays; nothing moves without a convective core (on the same structure, a surface convection zone, X = 0.70 above 0.45, whose lower edge would move with its composition, and a region off the centre) or without a mixed region (a row an earlier model left) (bug_sweep M241).

The mutations are `convection/test/test_convection.mutations.toml`; `python3 src/tools/mutate.py convection` (after `make tests`) reported 43 of 43 caught (2026-11-10). A check that fails makes the program exit 1. The program runs clean under the development flags (`make devel_tests`, `python3 src/tools/devel_tests.py`). A deliberate change to the output means regenerating `convection/test/expected_test_convection.out` in a commit that says why; `test_convection.py` compares the checks' verdicts, not their error digits (round-off differs between platforms).

The decks pin the rest: the configuration matrix's overshoot decks and the testsuite's solar decks log the mixed regions. No pinned deck reaches semiconvection (`cm_semiconvection`'s pre-main-sequence star is fully convective for all its models, bug_sweep M255); the UpperMS sample decks that use it are not pinned.

## 5. Coverage

`python3 src/tools/coverage.py --work DIR --dir convection --test test_convection --deck examples/run_config_matrix/cm_core_overshoot.inlist --deck examples/run_config_matrix/cm_env_overshoot.inlist --deck examples/run_config_matrix/cm_semiconvection.inlist --deck examples/run_config_matrix/cm_base.inlist` (a -O0 --coverage build in DIR) measured **97.8% of 465 lines** on 2026-11-10: convection_lib, convection_zones and overshoot 100%, semiconvection 98.7%, temperature_gradients 97.4%, convection_messages 88.5%. The lines not run are the mixing-length cubic's non-convergence (its run-log message and the adiabatic fallback: no test point makes 25 Chebyshev steps fail), the intermediate region's case of semiconvection (both edges; each edge's own path runs, the core's and the envelope's), and temperature_gradients' `end module` line, where gfortran puts the derived type's generated procedures. Since then semiconvection's downward branch and its intermediate case are gone (bug_sweep M241).

## 6. How to extend

Traced on 2026-11-10 on a throwaway branch (a least distance for core overshoot, a new `physics` control `overshoot_min_distance_core` read as `convection_settings%core_overshoot_min_distance`): the files below are the files the branch touched, and `make check`, `test_controls`, `test_generated_files`, `test_boundaries` and `test_convection` passed; the new check's mutation was caught.

**Add a setting** (a control convection reads; the settings record is the only way in):
1. The control, as `io/README.md` §5 ("Add a control") says: its declaration and description in `io/controls_buffer.f90` and its row in `defaults/controls_registry.tsv`, with the settings cell `convection_settings%<field>`. `make -j8` in `src/` regenerates `defaults/controls_registry.md`, `defaults/physics.defaults`, `io/controls_buffer_public.inc`, `io/controls_sync.f90`, `io/inlist_decl.inc`, `io/test/controls_dump.inc` and `state/physics_group_def.inc` (`control.defaults` and `control_group_def.inc` for a `&control` one; commit them).
2. The field in `convection_settings_t` (`convection_configuration.f90`) and its copy in `convection_settings_from_star` (`setup/convection_settings_builder.f90`); then `python3 tools/check_star_contracts.py --fix` in `src/`, which adds the control to the builder's `star reads:` line.
3. Read `convection_settings%<field>` where the physics needs it.
4. From the repository root: `make -C src tests`, `YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py` (`src/io/test/expected_controls/pristine.txt` gains the control's line), then `python3 -m pytest -q test_generated_files.py test_boundaries.py test_convection.py` and `make -C src check`.

**Add an overshoot prescription** (a different extent, or a cap):
1. Its setting(s), as above.
2. The extent in `overshoot.f90`: the region's branch in `extend_by_overshoot` (the distance, or the moved shell in `extend_upward`/`extend_downward`); for the surface convection zone, `envelope_overshoot_depth`, which that branch, core's Li/Be burn and rotation's sub-steps take.
3. For an envelope rule, the caller that applies its own envelope overshoot: the rezone's overshoot base (`core/model/rezoning.f90`, α taken as a depth in log P).
4. A check in `test_convection`'s `overshoot_checks` (the moved edge against the new extent, as the α H_P checks do), its mutation in `convection/test/test_convection.mutations.toml` (`python3 src/tools/mutate.py convection --only <id>`), and `convection/test/expected_test_convection.out` regenerated for the lines it adds.
5. With the new setting at its default nothing moves; a deck that sets it moves (a measured commit for a changed default, `docs/GUIDELINES.md` "Changes that move results").

**Add a kernel** (a new step of the gradient, a zone rule): keep it star-free (arguments in, results out, settings from `convection_settings`) and `pure` where it can be; put its messages in `convection_messages`; a new file goes in `check_boundaries.py`'s `STAR_BLIND_FILES` and `KERNEL_FILES`, and `python3 tools/gen_makefile_deps.py` in `src/` rewrites `deps.mk`; export from `convection_lib` what core or rotation call; give it an ok/FAIL check with a mutation.

## Sources

The routines are the Yale code's TPGRAD, CONVEC, OVERSH, HSUBP, MIXCZ and SCONVEC (2026 conversion: mixing/ became convection/, the settings record replaced the COMMON blocks, the 33-argument gradient took the EOS and opacity result arrays). Mixing-length theory after Böhm-Vitense (1958); semiconvection after Castellani et al. (1971, Ap&SS 10, 340); starspots after G. Somers (2014); the limit on small cores' overshoot after Woo & Demarque (2001). Each file's Sources section has its own history.

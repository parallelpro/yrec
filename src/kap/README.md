# kap/ — the opacities

This directory computes the Rosseland mean opacity of a zone, κ, and its logarithmic derivatives, from tables. Code outside kap/ uses it only through `kap_lib` (`tools/check_boundaries.py` enforces it): `kap_init` once per run, `kap_get` per zone, `kap_update_surface_tables` when the envelope X changes. kap reads nothing from the star: its configuration is one settings record and one paths record, built in `setup/kap_settings_builder.f90` and installed by `kap_init` at the start of every run.

The opacity tables come in **families** (OPAL 1995, OPAL 1992, LAOL 1989, Alexander 2006, Alexander 1994, Kurucz 1990, the pure-Z table, the Potekhin conductive table). Each family is a directory with two modules, its lookup (`<family>_opacity`) and its table reader (`<family>_reader`), and describes itself to kap_eval by a `kap_source` descriptor. `opacity_families` orders the families of each role in one precedence list, and `kap_eval` uses the first selected family of each role.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `kap_lib` (`kap_lib.f90`) | **The public interface**, and kap's state. `kap_init` checks the settings (`check_kap_settings`), installs them, chooses the families (`select_families`) and loads their tables; `kap_get` returns a zone's opacity in a result array (`i_kap`, `i_log10_kap`, `i_dlnkap_dlnrho`, `i_dlnkap_dlnt`); `kap_update_surface_tables`; `set_kap_settings` (installs settings without loading tables: for test programs). The private `kap_eval` evaluates the opacity in named steps (section 2). It holds the one instance of the tables (`tables`) and a never-written copy of their initial state (`initial_tables`). |
| `opacity_tables` (`opacity_tables.f90`) | The types: one table set per family (`opal95_table_set`, `opal92_table_set`, `laol_table_set`, `pure_z_table_set`, `alex06_table_set`, `alex94_table_set`, `kurucz_table_set`, `conductive_table_set`) collected in `opacity_table_state`; the dimensions, empty-cell markers and composition-cache tolerances; the settings record `kap_settings_t` (and `kap_settings`, protected, installed by `install_kap_settings`) and the paths record `kap_table_paths_t`; the four-point stencil locators. |
| `opacity_source` (`opacity_source.f90`) | What a family is to kap_eval: the `kap_source` descriptor (its name, whether the settings select it, its window in log T, whether it covers any Z itself or blends two tables in Z, and its `load`, `refresh`, `lookup` and `core_lookup` procedures), the point it is evaluated at (`kap_point`) and its result (`kap_result`). |
| `opacity_families` (`opacity_families.f90`) | The precedence lists: `interior_families` (OPAL 1995, OPAL 1992, LAOL 1989), `molecular_families` (Alexander 2006, Alexander 1994, Kurucz 1990), and `select_families`, the first selected of each with the pure-Z table, the conductive opacity and the ramp's windows. The one place that orders the families. |
| `opacity_blending` (`opacity_blending.f90`) | The formulas kap_eval combines opacities with, each once: `blend_in_z` and `blend_results_in_z` (linear in Z between two tables), `ramp_to_molecular`, `hubbard_lampe_applies` and `hubbard_lampe_opacity` (Hubbard & Lampe's conductive fit, its coefficients named), `add_conduction` (the harmonic sum). |
| `kap_messages` (`kap_messages.f90`) | Every message kap writes (section 3): `kap_report` (terminal and run log), `kap_terminal`, `kap_log`, and `kap_terminal_values`/`kap_log_values` for a list-directed record that ends in a number. |
| `opacity_table_files` (`opacity_table_files.f90`) | The readers' two helpers: `open_table` (status `old`, a failure reported) and `read_failed` (a read's iostat checked, the unit closed). |
| **opal95/** | |
| `opal95_opacity` | The OPAL and OP opacities in the OPAL 1995 layout (`GS98.OP17` and the others in `input/opacity/atomic/legacy/opal95/`): `build_opal95_fixed_z_table` (the tables at `opal95_single_table_z`), `build_opal95_surface_table` (at the envelope X too), `compute_opal95_opacity`, with 4-point Lagrange interpolation in 2D (T, log R: the surface table), 3D (X, T, log R: the fixed-Z table) or 4D (Z, X, T, log R: the full set); its interpolators share `opal95_rows_at_logr` and `extrapolate_in_logr`. The descriptor `opal95_source`. |
| `opal95_reader` | `read_opal95_tables`: the file's 126 tables, checked against the X and Z grids. |
| **opal92/** | |
| `opal92_opacity` | The 1992 OPAL tables, one file per Z (two when blending): `prepare_opal92_tables` (splines in log R), `build_opal92_surface_table`, `compute_opal92_opacity` (splines in log R and log T, linear in X); `opal92_source`. |
| `opal92_reader` | `read_opal92_tables` (3 X tables of 50 T rows per file; a shorter file is an error). |
| **laol89/** | |
| `laol_opacity` | The Los Alamos LAOL89 tables (one or two Z) and the pure-Z table (Z = 1, a pure mixture of equal numbers of C and O, `PURECO.DBGLAOL`) that helium-burning cores blend toward: `build_laol_splines`, `compute_laol_opacity`, `build_pure_z_splines`, `compute_pure_z_opacity`, both lookups through one column kernel, `laol_column_opacity`; `laol89_source` and `pure_z_source`. |
| `laol_reader` | `read_laol_tables`, `read_pure_z_table`. |
| **alex06/** | |
| `alex06_opacity` | The low-temperature molecular opacities of Ferguson et al. (2005), "Alexander 2006" (`alexmol06gs98.tab`): `build_alex06_composition_table` (the table at a composition, rebuilt when X or Z moves), `compute_alex06_opacity`; `alex06_source`. |
| `alex06_reader` | `read_alex06_tables`, every header checked. |
| **alex94/** | |
| `alex94_opacity` | The molecular opacities of Alexander & Ferguson (1994), seven files, one per X (the controls call them alex95: `use_alex95_tables`, `alex95_table_paths`): `build_alex94_fixed_z_table`, `build_alex94_surface_table`, `compute_alex94_opacity`; `alex94_source`. |
| `alex94_reader` | `read_alex94_tables`. |
| **kurucz90/** | |
| `kurucz_opacity` | The molecular opacities of Kurucz (1990), one table per [M/H] (two when blending): `build_kurucz_splines`, `compute_kurucz_opacity` (no value outside its table: the interior opacity is used alone there); `kurucz90_source`. |
| `kurucz_reader` | `read_kurucz_tables`. |
| **conductive/** | |
| `conductive_opacity` | Electron-conduction opacities from Potekhin et al. (1999): `compute_conductive_opacity` (the mixture of H, He and an O proxy for the metals, weighted by their ionized fractions), with Potekhin's interpolation routines under their own names (`condopacp`, `cinterp3`; the search is numerics_lib's `hunt`); `conductive_source`. |
| `conductive_reader` | `read_conductive_table` (`condall06.d`). |

`retired/conductive_opacity.f90.txt` holds the formulas the F77 conductive mixture kept as comments (an earlier mixing rule by mass), with how to revive them. Each file ends with a Sources section: the tables' references and the code's history.

## 2. The flow

**At the start of each run** (setup's `install_run_settings`, from `run_yrec`), `kap_init(settings, paths, X_env, Z_env, ierr)`:
1. `check_kap_settings` refuses settings that cannot give an opacity (section 3).
2. `set_kap_settings`: installs the settings (`install_kap_settings`) and chooses the families (`opacity_families`' `select_families`): the first selected interior and molecular family of each precedence list, the pure-Z table and the conductive opacity, with the windows of the ramp.
3. `tables = initial_tables`: every table set, cursor and cache starts as in a fresh process, so a second run in one process (pyyrec, `test_reentry`) starts where the first did.
4. Each chosen family's `load`: the conductive table, the interior family's tables (OPAL 1995: read, then the fixed-Z and surface tables; OPAL 1992 and LAOL 1989: read and spline, with the second Z table under `use_two_z_tables`), the pure-Z table, the molecular family's.

**When the envelope X changes** (setup's `install_surface_x_settings`, from core's starting model and rezoning), `kap_update_surface_tables(X, ierr)` calls the interior and molecular families' `refresh` (the surface tables of OPAL 1995, OPAL 1992 and Alexander 1994).

**For each zone** (core's shell physics, envelope integration and starting model, convection's semiconvection), `kap_get(log ρ, log T, X, Z, res, ion_fraction, ierr)` builds a `kap_point` and calls `kap_eval`, whose steps are:
1. `molecular_opacity`: the molecular family's `lookup` at log T ≤ `molecular_opacity_logt_max` (Kurucz 1990 may have no value there).
2. `interior_opacity`, unless log T < `molecular_opacity_logt_min` and the molecular family has a value: the interior family's `lookup`. In a helium-burning core (Z > 0.1 and log T > 7, `helium_core_opacity`) its `core_lookup` (OPAL 1995 at its top Z, 0.1; OPAL 1992 and LAOL 1989 at their first table's Z) is blended linearly in Z toward the pure-Z table, placed at Z = 1. At Z > 0.12 below 10⁷ K, or at a Z off the envelope's for a family that has no table there (OPAL 1992 or LAOL 1989 without `use_two_z_tables`), there is no interior opacity: ierr = 1.
3. The ramp (`ramp_to_molecular`): between the two windows' ends, κ = w κ_int + (1 − w) κ_mol with w linear in log T, and κ's own log derivatives (the κ-weighted blend of the two families' plus (dw/dln T)(κ_int − κ_mol)/κ in T); below `molecular_opacity_logt_min`, κ_mol.
4. Conduction: the Potekhin opacity when `use_conductive_opacity` and its table covers the point, else Hubbard & Lampe's fit where it applies (log T ≥ 4.2 and log ρ ≥ 2 log T − 13), added harmonically (`add_conduction`: 1/κ = 1/κ_rad + 1/κ_cond).

kap calls numerics_lib (`interp`, `lagrange4`, the splines, `locate`, `findex`, `hunt`) and math_lib; it uses state's `physical_constants` (ln 10), `composition_layout` (the LAOL header's metal slots) and `run_messages`' run-log unit (for kap_messages); nothing else.

## 3. Conventions

- **Units:** CGS; κ in cm²/g; log is log₁₀ (`log10_opacity`), except where a name or a comment says ln. The derivatives are logarithmic: dln κ/dln ρ at constant T and dln κ/dln T at constant ρ. The tables' density coordinate is log R = log ρ − 3 log T₆.
- **ion_fraction** (the ionized fractions of H, He⁺ and He⁺⁺, from the EOS) enters the conductive opacity only; kap does not change it.
- **Errors:** every routine that can fail returns `ierr`: 1 for an error kap finds, or a numerics routine's code (`numerics_status`) passed up as it is. Nothing in kap stops the run. `kap_init` returns 1 for refused settings (the pure-Z table without an interior table; two Z tables at one Z for a family that blends them; `molecular_opacity_logt_min` not below `_max` with a molecular family on), a table that cannot be opened (a missing path is never created), and a short or malformed table. `kap_get` returns ierr for a point off its tables (and leaves `res` as it was), its callers return on it.
- **Messages** go through `kap_messages`, and every one is written where the failure is found: errors `kap_init` finds go to the terminal and the run log (`kap_report`); the older diagnostics keep the destination and the words they have always had (the terminal, `kap_terminal`, or the run log, `kap_log`; their rewording, sentence case and every error a report, is on the measured-change list). The table readers and the lookups' failure paths are kap's only writers; the lookups are listed in `check_boundaries.py`'s `STAR_BLIND_FILES`, kap_lib and the type modules in `KERNEL_FILES`, so `check_no_io` keeps them free of I/O statements.
- **Controls:** from `kap_settings_t` and `kap_table_paths_t` only, never `s%control`/`s%physics`; a new control goes through io/README.md §5.
- **State:** kap_lib's `tables` is the only state, and only kap_lib writes it: the family routines get their table set as an argument, OPAL95's per-call stencil is a local, and the cursors and composition caches start from `initial_tables` at every `kap_init`.
- **Kept by decision, or open:**
  - LAOL89's derivatives are the secants of the knots around the point (within 10% of the derivatives of the interpolated κ; `test_kap` checks that);
  - the conductive opacity's branch below log ρ = −6 takes the mixture term with the opposite sign and extrapolates the wrong way (M222, κ moves by about 10⁻⁶ in the photosphere: a measured fix);
  - the LAOL tables' extrapolation distance is ln(`laol_z_tolerance`) on a log₁₀ ρ axis (M227: 2.3 dex for the default 10);
  - OPAL 1995's low-X, low-T corner (X < 0.2, the T stencil's first row at or below the fifth) is not checked against the empty part of the tables (the F77 code's postponed check);
  - a family selected but not first in its role is neither loaded nor used;
  - the conductive mixture weights each species by its ionized number fraction without renormalizing: a partly ionized gas conducts more than a fully ionized one in this formula (the F77 code's).

## 4. Tests

`test_kap` (`make tests` in `src/`, then `pytest test_kap.py` from the repository root, exact and with `YREC_TEST_TOLERANT=1`) prints the opacity on a grid of points for every family (the pinned part: the production combination at and away from the surface table and on a helium core, OPAL92 and LAOL89 with two Z tables, Kurucz and Alexander 1994, and two error paths), then 35 ok/FAIL checks; a FAIL makes the program exit 1. The checks:
- **exactness at the table nodes:** OPAL 1995's 2D, 3D and 4D paths each return the file's value at a node (the tables loaded at the node Z = 0.02 and X = 0.70, the test reading the file's entry itself); Alexander 2006 returns a 3-decimal table value at a node;
- **derivatives against centred differences** of log κ (step 10⁻⁴ in log ρ and log T), family by family: OPAL 1995 on each path, Alexander 2006, Potekhin conduction with ionized fractions, the helium core's pure-Z blend, the ramp between Alexander 2006 and OPAL 1995 (at log T 4.07), OPAL 1992 with two Z, LAOL89 (its secants, to 10%), Kurucz 1990 with two Z, Alexander 1994; and Hubbard & Lampe's fit;
- **continuity:** κ across the ramp's two ends, across OPAL 1995's lower log R edge (extrapolated below), and from the 2D to the 4D path in Z; the surface table rebuilt by `kap_update_surface_tables` at a new X gives the 3D path's value there;
- **the formulas:** `blend_in_z` returns each table's value at its own Z, and OPAL 1992's two-Z blend at `opal_table_z1` is its first table; `add_conduction` halves κ for equal opacities and averages their slopes; the Potekhin mixture's κ_c is Σ wᵢ κᵢ for a partly ionized mixture;
- **the error contracts:** `kap_init` refuses each of the three bad settings, a missing table (not creating it) and a short OPAL 1992 table; `kap_get` returns ierr off the tables (with `res` unchanged), at Z > 0.12 below 10⁷ K and for a helium core without the pure-Z table; a refused Alexander 2006 composition fails again on a retry; a second `kap_init` gives the first run's first lookup bit for bit (at a node, where the warm-start cursors matter).

Each check was confirmed to fail under a mutation of the code it covers: `kap/test/test_kap.mutations.toml` holds 38 mutations, and `python3 src/tools/mutate.py kap` (after `make tests` in `src/`) reported all of them caught on 2026-11-10. `test_kap` runs clean under the development flags (`make devel_tests`, `python3 src/tools/devel_tests.py`).

The pinned output is `kap/test/expected_test_kap.out`; a deliberate change means regenerating it from `src/test_kap` (run from the repository root with `YREC_INPUT=input`) and saying why in the commit. The decks pin the production combination end to end (`compare_decks.py`; every shipped deck uses OPAL 1995, Alexander 2006, the pure-Z table and Potekhin conduction): `test_kap` is the only coverage of OPAL92, LAOL89, Kurucz 1990 and Alexander 1994.

## 5. Coverage

79.4% of kap/'s 2175 executable lines run (2026-11-10), under `test_kap` and four configuration-matrix decks (`cm_base`, `cm_settle_he_z` for the 4D path under Z diffusion, `cm_atm_krishnaswamy`, `cm_core_overshoot`), measured with
```
python3 src/tools/coverage.py --work DIR --dir kap --test test_kap \
    --deck examples/run_config_matrix/cm_base.inlist --deck examples/run_config_matrix/cm_settle_he_z.inlist \
    --deck examples/run_config_matrix/cm_atm_krishnaswamy.inlist --deck examples/run_config_matrix/cm_core_overshoot.inlist
```
(a -O0 --coverage build in DIR). By file: `opacity_blending`, `opacity_families` and the conductive reader 100%; `kap_lib` 96.6%; `opacity_tables` 97.2%; the lookups 72-89% (OPAL 1995 85.8%, Alexander 2006 88.9%, conductive 84.2%, Kurucz 76.9%, LAOL 74.8%, OPAL 1992 74.0%, Alexander 1994 71.7%); the readers 54-87%. The lines not run are the readers' and lookups' failure branches for tables YREC does not ship (a mismatched header, a malformed row, a row of fewer than 4 points), the messages of those branches, Alexander 1994's lookup off the surface X, and OPAL 1992's and LAOL's extrapolation past their rows.

## 6. How to extend

**Add a control** that kap reads: io/README.md §5 ("Add a control"; "Add a table path" for a table's file), with its field in `kap_settings_t` (or `kap_table_paths_t`) in `opacity_tables.f90` and its copy in `kap_settings_from_star` (or `kap_table_paths_from_star`) in `setup/kap_settings_builder.f90`.

**Add a table family** (traced on 2026-11-10 on a throwaway branch with a toy molecular family, `flat`, a constant opacity read from a one-number file: the files below are the ones it touched):
1. Its directory, `kap/<family>/`, with two modules:
   - `<family>_reader.f90`: `read_<family>_table(tbl, path, ierr)`, opening with `opacity_table_files`' `open_table` and checking every read with `read_failed`;
   - `<family>_opacity.f90`: the lookup (`compute_<family>_opacity`, the table set and the point as arguments, no I/O statement), and the descriptor `<family>_source(settings)` with its procedures behind `opacity_source`'s interfaces: `load` (the reader, then whatever it builds), `lookup`, and where needed `refresh` (a surface table) and `core_lookup` (an interior family, for helium-burning cores). Set `covers_any_z`, or `blends_two_z` with `table_z` and `table_z_controls`, as the family is.
2. `kap/opacity_tables.f90`: its `<family>_table_set` type, its component in `opacity_table_state` (in order of size: on arm64 Mach-O a component must start within 8 MB of the state's start), its `use_<family>_tables` field in `kap_settings_t` and its path field in `kap_table_paths_t`.
3. `kap/opacity_families.f90`: its `use` line and its entry in its role's list, at its place in the precedence (the list's length grows by one).
4. The controls, by io/README.md §5 ("Add a control", "Add a table path"): their declarations in `io/controls_buffer.f90`, their rows in `defaults/controls_registry.tsv` (settings cells `kap_settings%use_<family>_tables` and `kap_table_paths%<family>`), and their copies in `kap_settings_from_star` and `kap_table_paths_from_star` (`setup/kap_settings_builder.f90`).
5. `src/tools/check_boundaries.py`: the lookup in `STAR_BLIND_FILES`. `src/tools/style_counts.tsv`: the new directory's `io_statement` count, the reader's opens and reads (a new directory starts at 0, and the reader may keep them; a raise is a hand edit the commit explains).
6. A check in `test_kap` (its nodes, its derivatives against differences), its mutation in `kap/test/test_kap.mutations.toml`, and this README's row in section 1.
7. From `src/`: `python3 tools/gen_makefile_deps.py` (two new modules); `make -j8` (it reruns the controls generators: `defaults/control.defaults`, `defaults/controls_registry.md`, `io/controls_buffer_public.inc`, `io/controls_sync.f90`, `io/expand_paths.inc`, `io/inlist_decl.inc`, `io/test/controls_dump.inc` and `state/control_group_def.inc` change); `python3 tools/check_star_contracts.py --fix` (the builder's `star reads:` lines); `make check`; `make tests`; from the repository root `YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py` (`io/test/expected_controls/pristine.txt` gains the two controls) and `pytest test_kap.py` (`kap/test/expected_test_kap.out` gains the check's line); `python3 src/tools/mutate.py kap`.

So a family is its directory (a lookup and a reader), its table-set type, one list entry, two settings fields and two controls: eight hand-edited files besides the test and the README, and the generated ones make writes.

**Add a check:** a call in `test_kap`'s `run_checks` (the helpers `node_check`, `derivative_check`, `continuity_check` and `report` print the line), its mutation in `kap/test/test_kap.mutations.toml` (`python3 src/tools/mutate.py kap --check`, then the run), and its line in `kap/test/expected_test_kap.out`.

# kap/ — the opacities

This directory computes the Rosseland mean opacity of a zone, κ, and its logarithmic derivatives, from tables. Code outside kap/ uses it only through `kap_lib` (`tools/check_boundaries.py` enforces it): `kap_init` once per run, `kap_get` per zone, `kap_update_surface_tables` when the envelope X changes, `kap_set_envelope_metal_fraction` when the envelope Z changes. kap reads nothing from the star: its configuration is one settings record and one paths record, built in `setup/kap_settings_builder.f90` and installed by `kap_init` at the start of every run.

The opacity tables come in **families** (OPAL 1995, OPAL 1992, LAOL 1989, Alexander 2006, Alexander 1994, Kurucz 1990, the pure-Z table, the Potekhin conductive table). Each family is a directory with two modules, its lookup (`<family>_opacity`) and its table reader (`<family>_reader`), and describes itself to kap_eval by a `kap_source` descriptor. `opacity_families` orders the families of each role in one precedence list, and `kap_eval` uses the first selected family of each role.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `kap_lib` (`kap_lib.f90`) | **The public interface**, and kap's state. `kap_init` checks the settings (`check_kap_settings`), installs them, chooses the families (`select_families`) and loads their tables; `kap_get` returns a zone's opacity in a result array (`i_kap`, `i_log10_kap`, `i_dlnkap_dlnrho`, `i_dlnkap_dlnt`); `kap_update_surface_tables`; `kap_set_envelope_metal_fraction`; `set_kap_settings` (installs settings without loading tables: for test programs). The private `kap_eval` evaluates the opacity in named steps (section 2). It holds the installed settings record (`kap_settings`), the one instance of the tables (`tables`) and a never-written copy of their initial state (`initial_tables`). |
| `opacity_tables` (`opacity_tables.f90`) | The types: one table set per family (`opal95_table_set`, `opal92_table_set`, `laol_table_set`, `pure_z_table_set`, `alex06_table_set`, `alex94_table_set`, `kurucz_table_set`, `conductive_table_set`) collected in `opacity_table_state`; the dimensions, empty-cell markers and composition-cache tolerances; the settings record's type `kap_settings_t` and the paths record `kap_table_paths_t`; the four-point stencil locators; `alex_table_at_point`, the two Alexander lookups' interpolation in one table (in log R at four T rows, with the linear extrapolation above the top log R, then in log T). |
| `opacity_source` (`opacity_source.f90`) | What a family is to kap_eval: the `kap_source` descriptor (its name, whether the settings select it, its window in log T, whether it covers any Z itself or blends two tables in Z, and its `load`, `refresh`, `lookup` and `core_lookup` procedures), the point it is evaluated at (`kap_point`) and its result (`kap_result`). |
| `opacity_families` (`opacity_families.f90`) | The precedence lists: `interior_families` (OPAL 1995, OPAL 1992, LAOL 1989), `molecular_families` (Alexander 2006, Alexander 1994, Kurucz 1990), and `select_families`, the first selected of each with the pure-Z table, the conductive opacity and the ramp's windows. The one place in kap that orders the families (io switches off the molecular controls that lose before kap sees them: section 6). |
| `opacity_blending` (`opacity_blending.f90`) | The formulas kap_eval combines opacities with, each once: `blend_in_z` and `blend_results_in_z` (linear in Z between two tables), `ramp_to_molecular`, `hubbard_lampe_applies` and `hubbard_lampe_opacity` (Hubbard & Lampe's conductive fit, its coefficients named), `add_conduction` (the harmonic sum). |
| `kap_messages` (`kap_messages.f90`) | Every message kap writes (section 3), through `kap_report`: to the terminal and, through `run_messages`' `log_line`, to the run log. |
| `opacity_table_files` (`opacity_table_files.f90`) | The readers' helpers: `open_table` (status `old`, a failure reported), `read_failed` (a read's iostat checked) and `table_z_differs` (the Z a file states against its control's); a refused table's unit is closed and its reason reported in one place (`refuse_table`). |
| **opal95/** | |
| `opal95_opacity` | The OPAL and OP opacities in the OPAL 1995 layout (`GS98.OP17` and the others in `input/opacity/atomic/legacy/opal95/`): `build_opal95_fixed_z_table` (the tables at `opal95_single_table_z`), `build_opal95_surface_table` (at the envelope X too), `compute_opal95_opacity`, with 4-point Lagrange interpolation in 2D (T, log R: the surface table), 3D (X, T, log R: the fixed-Z table) or 4D (Z, X, T, log R: the full set). The lookup is five named steps (`locate_opal95_t_and_logr`, `locate_opal95_composition`, `shift_rows_off_empty_corner`, `opal95_logr_weights`, then `interpolate_opal95_2d`, `_3d` or `_4d`, which share `opal95_rows_at_logr` and `extrapolate_in_logr`); `opal95_x_tables_at_z` is the one rule for the X tables that exist at a Z. The descriptor `opal95_source`. |
| `opal95_reader` | `read_opal95_tables`: the file's 126 tables, each checked against the X and Z grids; `opal95_row_columns`, the number of log R columns a row has. |
| **opal92/** | |
| `opal92_opacity` | The 1992 OPAL tables, one file per Z (two when blending): `prepare_opal92_tables` (`spline_opal92_rows`: splines in log R), `build_opal92_surface_table` (`build_opal92_surface_rows`), `compute_opal92_opacity` (splines in log R and log T, `interpolate_opal92_2d`; linear in X); `opal92_source`. |
| `opal92_reader` | `read_opal92_tables` (`read_opal92_table`: 3 X tables of 50 T rows per file; a shorter file, or one that states another Z than its control's, is an error). |
| **laol89/** | |
| `laol_opacity` | The Los Alamos LAOL89 tables (one or two Z) and the pure-Z table (Z = 1, a pure mixture of equal numbers of C and O, `PURECO.DBGLAOL`) that helium-burning cores blend toward: `build_laol_splines` and `build_pure_z_splines` (both through `spline_row_in_log_rho`), `compute_laol_opacity`, `compute_pure_z_opacity`, both lookups through one column kernel, `laol_column_opacity`; `laol89_source` and `pure_z_source`. |
| `laol_reader` | `read_laol_tables` (`read_laol_table`; a file that states another Z than its control's is an error), `read_pure_z_table`. |
| **alex06/** | |
| `alex06_opacity` | The low-temperature molecular opacities of Ferguson et al. (2005), "Alexander 2006" (`alexmol06gs98.tab`): `build_alex06_composition_table` (the table at a composition, rebuilt when X or Z moves), `compute_alex06_opacity` (the stencils' weights, then `alex_table_at_point`); `alex06_source`. |
| `alex06_reader` | `read_alex06_tables` (`read_alex06_table`, every header and row checked). |
| **alex94/** | |
| `alex94_opacity` | The molecular opacities of Alexander & Ferguson (1994), seven files, one per X (the controls call them alex95: `use_alex95_tables`, `alex95_table_paths`): `build_alex94_fixed_z_table`, `build_alex94_surface_table`, `compute_alex94_opacity` (the stencils' weights, then `alex_table_at_point` in the surface table, or in four fixed-Z tables and in X); `alex94_source`. |
| `alex94_reader` | `read_alex94_tables`. |
| **kurucz90/** | |
| `kurucz_opacity` | The molecular opacities of Kurucz (1990), one table per [M/H] (two when blending): `build_kurucz_splines`, `compute_kurucz_opacity` (`kurucz_rows_at_density`, `locate_kurucz_density`; no value outside its table: the interior opacity is used alone there); `kurucz90_source`. |
| `kurucz_reader` | `read_kurucz_tables` (`read_kurucz_table`). |
| **conductive/** | |
| `conductive_opacity` | Electron-conduction opacities from Potekhin et al. (1999): `compute_conductive_opacity` (the mixture of H, He and an O proxy for the metals by Potekhin et al.'s rule: each species weighted by its share of the electrons), `species_conductivities` (each species' table read at the density where it alone has the mixture's electron density), with Potekhin's interpolation routines under their own names (`condopacp`, `cinterp3`; the search is numerics_lib's `hunt`); `conductive_source`. |
| `conductive_reader` | `read_conductive_table` (`condall06.d`). |

`retired/conductive_opacity.f90.txt` holds the formulas the F77 conductive mixture kept as comments (an earlier mixing rule by mass), with how to revive them. Each file ends with a Sources section: the tables' references and the code's history (the F77 routines' names are there).

## 2. The flow

**At the start of each run** (setup's `install_run_settings`, from `run_yrec`), `kap_init(settings, paths, X_env, Z_env, ierr)`:
1. `check_kap_settings` refuses settings that cannot give an opacity (section 3).
2. `set_kap_settings`: installs the settings and chooses the families (`opacity_families`' `select_families`): the first selected interior and molecular family of each precedence list, the pure-Z table and the conductive opacity, with the windows of the ramp.
3. `tables = initial_tables`: every table set, cursor and cache starts as in a fresh process, so a second run in one process (pyyrec, `test_reentry`) starts where the first did. The envelope Z starts at `Z_env`, the run plan's.
4. Each chosen family's `load`: the conductive table, the interior family's tables (OPAL 1995: read, then the fixed-Z and surface tables; OPAL 1992 and LAOL 1989: read and spline, with the second Z table under `use_two_z_tables`), the pure-Z table, the molecular family's. A family checks at its load what only it knows (section 3).

**When the envelope X changes** (setup's `install_surface_x_settings`, from core's starting model and rezoning), `kap_update_surface_tables(X, ierr)` calls the interior and molecular families' `refresh` (the surface tables of OPAL 1995, OPAL 1992 and Alexander 1994).

**When the envelope Z changes** (setup's `install_surface_z_settings`: the run plan's Z, the first card's model, a Z rescale), `kap_set_envelope_metal_fraction(Z)`: an interior family with a single table that does not interpolate in Z serves the zones at that Z.

**For each zone** (core's shell physics, shell equations, envelope integration and starting model; convection's temperature gradients and semiconvection), `kap_get(log ρ, log T, X, Z, res, ion_fraction, ierr)` builds a `kap_point` and calls `kap_eval`, whose steps are:
1. `molecular_opacity`: the molecular family's `lookup` at log T ≤ `molecular_opacity_logt_max` (Kurucz 1990 may have no value there).
2. `interior_opacity`, unless log T < `molecular_opacity_logt_min` and the molecular family has a value: the interior family's `lookup`. In a helium-burning core (Z > 0.1 and log T > 7, `helium_core_opacity`) its `core_lookup` (OPAL 1995 at its top Z, 0.1; OPAL 1992 and LAOL 1989 at their first table's Z) is blended linearly in Z toward the pure-Z table, placed at Z = 1. At Z > 0.12 below 10⁷ K, without an interior family, or at a Z off the envelope's for a family that has no table there (OPAL 1992 or LAOL 1989 without `use_two_z_tables`), there is no interior opacity: ierr = 1.
3. The ramp (`ramp_to_molecular`), where the molecular family has a value: between the two windows' ends, κ = w κ_int + (1 − w) κ_mol with w linear in log T, and κ's own log derivatives (the κ-weighted blend of the two families' plus (dw/dln T)(κ_int − κ_mol)/κ in T); below `molecular_opacity_logt_min`, κ_mol.
4. Conduction: the Potekhin opacity when `use_conductive_opacity` and its table covers the point (and the gas has free electrons), else Hubbard & Lampe's fit where it applies (log T ≥ 4.2 and log ρ ≥ 2 log T − 13), added harmonically (`add_conduction`: 1/κ = 1/κ_rad + 1/κ_cond).

kap calls numerics_lib (`interp`, `lagrange4`, the splines, `locate`, `findex`, `hunt`) and math_lib; it uses state's `physical_constants` (ln 10, and σ for the conductive opacity) and `run_messages`' `log_line` (for kap_messages); nothing else.

## 3. Conventions

- **Units:** CGS; κ in cm²/g; log is log₁₀ (`log10_opacity`), except where a name or a comment says ln. The derivatives are logarithmic: dln κ/dln ρ at constant T and dln κ/dln T at constant ρ. The tables' density coordinate is log R = log ρ − 3 log T₆.
- **ion_fraction** (the ionized fractions of H, He⁺ and He⁺⁺, from the EOS) enters the conductive opacity only; it is `intent(in)`.
- **Errors:** every routine that can fail returns `ierr`: 1 for an error kap finds, or a numerics routine's code (`numerics_status`) passed up as it is. Nothing in kap stops the run.
  - `kap_init` returns 1 for refused settings (the pure-Z table without an interior table; two Z tables at one Z for a family that blends them; `molecular_opacity_logt_min` not below `_max` with a molecular family on), for what a family checks at its load (an OPAL 1992 or LAOL 1989 file that states another Z than its control's, beyond 10⁻⁴, since the files round it; `laol_z_tolerance` not positive; a Kurucz 1990 file without rows), a table that cannot be opened (a missing path is never created), and a short or malformed table (a header or a row that is not the grid's).
  - `kap_get` returns ierr (and leaves `res` as it was) for log T off the interior tables, for Z > 0.12 below 10⁷ K, for a helium core without the pure-Z table, for a composition off the Alexander 2006 tables, and for a point that OPAL 1992, LAOL 1989 or the pure-Z table do not cover. **It does not for a point beyond the log R range of OPAL 1995 or of the Alexander tables, or below the Alexander tables' lowest log T: those are extrapolated** (the list below says how); Kurucz 1990 hands such a point to the interior family. Its callers return on ierr.
- **Messages** are errors only, one line each in sentence case after the entry point that met them (`kap_init: ...`, `kap_get: ...`, `kap_update_surface_tables: ...`), and every one goes to the terminal and the run log (`kap_messages`' `kap_report`; the run log's line through `log_line`, so a message from before the log opens waits for it). The table readers report their own errors at `kap_init`. The lookups and the tables' builds write nothing: each is `pure` and returns ierr with a one-line reason in `message`, through the `kap_source` procedures, and kap_lib reports it (`kap_eval`, `kap_init`'s loads, `kap_update_surface_tables`). The lookups are listed in `check_boundaries.py`'s `STAR_BLIND_FILES`, kap_lib and the type modules in `KERNEL_FILES`, so `check_no_io` keeps them free of I/O statements.
- **Controls:** from `kap_settings_t` and `kap_table_paths_t` only, never `s%control`/`s%physics`; a new control goes through io/README.md §5.
- **State:** kap_lib's `kap_settings`, `tables` and chosen families are the only state, and only kap_lib writes them: the family routines get the settings and their table set as arguments, OPAL95's per-call stencil is a local, and the cursors and composition caches start from `initial_tables` at every `kap_init`.
- **10^x** is written three ways, each family's own since the F77: math_lib's `exp10` (OPAL 1992, LAOL, Kurucz, conduction, Hubbard & Lampe), `pow(10, x)` (OPAL 1995) and `exp(ln10 x)` (the two Alexander lookups). They differ at round-off, so one spelling for all would move every deck; kept.
- **Kept by decision, or open:**
  - **the lookups that extrapolate with ierr 0** (M342; the F77 code's; a refusal would stop a shipped deck, so they stay):
    - OPAL 1995 in log R: linearly from the table's first column below it, without limit; and at the empty high-R, high-T corner, past a row's last filled column: by the cubic of the row's last four columns up to one column beyond it, linearly from that column further out. log κ jumps where a row passes from the one to the other (measured at the fixed Z 0.016232, X = 0.70: −0.007 at log T 7.45 and log R 0, −0.017 at log T 7.95 and log R −0.5);
    - Alexander 2006 and Alexander 1994 in log R and log T: from their last four nodes (the cubic) below the first log R and below the first log T, and linearly from the last two columns above the top log R (Alexander 1994: at the surface X; off it, the cubic). Above the top log R the returned derivatives are the cubic's at the edge, not the line's (at log T 3.9954: dln κ/dln ρ 0.415 returned for 0.388, dln κ/dln T 7.07 for 7.15);
    - OPAL 1992 in log R, linearly from a row's ends, without limit; LAOL 1989 and the pure-Z table up to their tolerance (the next line);
    - where the shipped decks go off the tables (counted 2026-11-18 by a build that reports the Alexander 2006 and OPAL 1995 lookups' points off their tables, on 16 of the 56 decks of `compare_decks --set all`, at least one per kind of run): only the 0.03 M☉ chain (`m0030gs98z018804_prems`, `m0030gs98z018804_Dbl`, `Test_m0030_feh0_noGS_ZAMS`) does, with Alexander 2006 above its top log R (1.0) at log T 3.37–4.00, up to log R 3.79 (2.8 dex past the table) in the ZAMS deck; no counted deck goes below Alexander 2006's first log R or log T, below OPAL 1995's first log R, or past a row's last filled column of OPAL 1995. A refusal there would stop those three decks, so the extrapolation stays and is stated here.
  - LAOL89's derivatives are the secants of the knots around the point (within 10% of the derivatives of the interpolated κ; `test_kap` checks that; bug_sweep M280);
  - the LAOL tables' extrapolation distance is ln(`laol_z_tolerance`) on a log₁₀ ρ axis (M227: 2.3 dex for the default 10);
  - OPAL 1995's low-X, low-T corner (X < 0.2, the T stencil's first row at or below the fifth) is not checked against the empty part of the tables (the F77 code's postponed check);
  - OPAL 1995's X = 0.95 table at Z = 0.04 is in the files and not used: the stencils treat Z ≥ 0.04 as without it;
  - a family selected but not first in its role is neither loaded nor used (io warns for two molecular families and switches the others off; nothing warns for two interior families: M344);
  - below log ρ = −6 the conductive opacity is extrapolated in log ρ with the mixture's slope at −6, its dln κ/dln T the mixture's there (round 3B; the H table's alone before, the F77 code's: M299), and it stops where log R crosses 0 (M299's second half, the authors' decision);
  - the conductive mixture's atomic weights are the F77 code's (1.008, 4.004, 16);
  - the conductive reader opens its table itself, without the compiler's reason in its message, where the other readers open through `open_table`.

## 4. Tests

`test_kap` (`make tests` in `src/`, then `pytest test_kap.py` from the repository root, exact and with `YREC_TEST_TOLERANT=1`) prints the opacity on a grid of points for every family (the pinned part: the production combination at and away from the surface table and on a helium core, OPAL92 and LAOL89 with two Z tables, Kurucz and Alexander 1994, and two error paths), then 86 ok/FAIL checks; a FAIL makes the program exit 1. The checks, by kind:
- **exactness at the table nodes:** OPAL 1995's 2D, 3D and 4D paths, its X = 1 − Z table, and OPAL 1995 alone below the ramp's top each return the file's value at a node (the test reads the file's entry itself); Alexander 2006 returns a 3-decimal table value at a node;
- **derivatives against centred differences** of log κ (step 10⁻⁴ in log ρ and log T), family by family: OPAL 1995 on each path, where its rows shift at the empty high-R corner and where they are extrapolated past it, Alexander 2006, Potekhin conduction with ionized fractions, the helium core's pure-Z blend, the ramp between Alexander 2006 and OPAL 1995 (at log T 4.07), OPAL 1992 with two Z, LAOL89 (its secants, to 10%), Kurucz 1990 with two Z, Alexander 1994; and Hubbard & Lampe's fit;
- **continuity and consistency:** κ across the ramp's two ends, across OPAL 1995's lower log R edge (extrapolated below), and from the 2D to the 4D path in Z (at the run's Z and at Z = 0.04); the surface table rebuilt by `kap_update_surface_tables` gives the 3D path's value (at X = 0.65, 0.92, and 0.85 at Z = 0.06); the Potekhin κ_c across log ρ = −6 and the extrapolation's slope below it; OPAL 1995's high-R shift on a table whose empty cells differ by X (a doctored copy): independent of the previous lookup's X, no empty cell read, and at the real corner log κ within 0.05 of its neighbours';
- **the formulas:** `blend_in_z` returns each table's value at its own Z, and OPAL 1992's two-Z blend at `opal_table_z1` is its first table; `add_conduction` halves κ for equal opacities and averages their slopes; the Potekhin mixture's 1/K is Σ xᵢ/Kᵢ(ρᵢ) for a partly ionized mixture, and a gas without free electrons has no conductive opacity; the precedence lists (two families of a role selected: the first alone, bit for bit); Kurucz 1990 without a value: the interior opacity alone;
- **the error contracts:** `kap_init` refuses each of the three bad settings, a table file not at its control's Z (OPAL 1992, LAOL 1989), a `laol_z_tolerance` that is not positive, a missing table through each of the two opens (not creating it), and the readers' refusals on doctored copies the test writes (a short OPAL 1992 table; an OPAL 1995 file with a wrong header, cut inside a table, cut after one, with a row that is not numbers; an Alexander 2006 table with a wrong X, log R column or log T row; an Alexander 1994 table with a wrong Z or row length; a LAOL 1989 or pure-Z file with too many X; a Kurucz 1990 file without rows or with an unreadable row); `kap_get` returns ierr off the tables (with `res` unchanged), at Z > 0.12 below 10⁷ K, for a helium core without the pure-Z table, and without an interior family (by its message); a refused Alexander 2006 composition fails again on a retry; a single-Z table serves the envelope's Z after `kap_set_envelope_metal_fraction` moves it; a second `kap_init` gives the first run's first lookup bit for bit; a message before the run log opens waits for it.

Each check was confirmed to fail under a mutation of the code it covers: `kap/test/test_kap.mutations.toml` holds 93 mutations, and `python3 src/tools/mutate.py kap` (after `make tests` in `src/`) reported all of them caught on 2026-11-18. `test_kap` runs clean under the development flags (`make devel_tests`, `python3 src/tools/devel_tests.py`). The checks, in the order `test_kap` prints them, each with the mutations that make it FAIL (the manifest's ids):

| Check (the line test_kap prints) | Mutations that make it FAIL |
|---|---|
| OPAL95 2D path (surface X, Z): exact at a table node | `opal95-2d-node` |
| OPAL95 3D path (X off the surface): exact at a table node | `opal95-3d-node` |
| OPAL95 4D path (Z off the fixed-Z table): exact at a node | `opal95-4d-node` |
| Alexander 2006: exact at a table node (log T 3.40, log R -3) | `alex06-node` |
| OPAL95 3D path: exact at a node of the X = 1 - Z table | `opal95-x-1mz-node` |
| OPAL95 2D: dln kappa/dln rho, dln T against differences | `opal95-2d-dt` |
| OPAL95 3D: dln kappa/dln rho, dln T against differences | `opal95-3d-drho` |
| OPAL95 4D: dln kappa/dln rho, dln T against differences | `opal95-4d-drho` |
| Alexander 2006: dln kappa/dln rho, dln T against differences | `alex06-dt` |
| Potekhin conduction (ionized): derivatives against differences | `conductive-dt`, `conduction-sum-drho` |
| helium core (pure Z blend): derivatives against differences | `helium-core-blend` |
| ramp: kappa continuous at molecular_opacity_logt_min | `ramp-min` |
| ramp: kappa continuous at molecular_opacity_logt_max | `ramp-max` |
| ramp: dln kappa/dln rho, dln T against differences | `ramp-dweight`, `ramp-drho-weights` |
| OPAL95: kappa continuous at the table's lowest log R | `opal95-extrapolation` |
| OPAL95: the 2D path continuous with the 4D path in Z | `opal95-fixed-z` |
| kap_update_surface_tables: the new 2D path is the old 3D path | `opal95-surface-x` |
| OPAL95 surface table at X = 0.92: the 2D path is the 3D path | `opal95-surface-x-above-0.9` |
| OPAL95 surface table at Z = 0.06, X = 0.85: 2D is 3D | `opal95-surface-x-weights` |
| OPAL95 at Z = 0.04, X = 0.85: 2D continuous with 4D in Z | `opal95-surface-no-x095` |
| OPAL95 high-R shift: independent of the previous lookup's X | `opal95-shift-stale-x` |
| OPAL95 high-R shift: no empty cell of the stencil's X tables read | `opal95-shift-first-x-table`, `opal95-shift-surface-tables` |
| OPAL95 high-R shift (rows end at different log R): derivatives | `opal95-shift-derivative` |
| OPAL95 high-R shift: log kappa within 0.05 of its neighbours' | `opal95-shift-row-end`, `opal95-no-shift` |
| OPAL95 past a row's last log R (extrapolated): derivatives | `opal95-extrapolation-high` |
| precedence: OPAL95 with OPAL92 also selected is OPAL95 alone | `precedence-interior` |
| precedence: Alexander 2006 with 1994 also selected is 2006's | `precedence-molecular` |
| no molecular family: OPAL95 exact at a node below log T 4.1 | `ramp-without-molecular` |
| Potekhin mixture: 1/K = sum x_i/K_i(rho_i), electron fractions | `conductive-he-weight`, `conductive-weights-normalized`, `conductive-electron-density` |
| Potekhin conduction: kappa_c continuous at log rho = -6 | `conductive-low-density-base` |
| Potekhin conduction below log rho = -6: dln kappa/dln rho | `conductive-low-density-direction` |
| Potekhin conduction: none for a gas without free electrons | `conductive-no-electrons` |
| OPAL92 (two Z): derivatives against differences | `opal92-dt` |
| OPAL92 two-Z blend: at opal_table_z1 it is the first table | `blend-at-z1` |
| LAOL89 (secant derivatives): against differences, to 10% | `laol-dt` |
| Kurucz 1990 (two Z): derivatives against differences | `kurucz-drho` |
| Alexander 1994: derivatives against differences | `alex94-dt` |
| Alexander 1994 off the surface X: derivatives against differences | `alex94-x-drho` |
| Alexander 1994 at another Z: its tables rebuilt, exact at a node | `alex94-rebuild-at-z` |
| Alexander 1994 surface table at X = 0.62: 2D is 3D | `alex94-surface-x` |
| Alexander 1994 above its top log R: linear from the last two | `alex94-top-logr` |
| Alexander 2006 above its top log R: linear from the last two | `alex06-top-logr` |
| Kurucz 1990: a density no row has returns ierr | `kurucz-no-value` |
| OPAL92 between two X tables: derivatives against differences | `opal92-x-interpolation` |
| OPAL92 at an X of its grid: the X interpolation beside it | `opal92-grid-x-table` |
| OPAL92: log T off the table returns ierr | `opal92-t-range` |
| helium core (OPAL92): the Z blend of its table and the pure-Z | `opal92-core-table-z` |
| helium core (LAOL89): the Z blend of its table and the pure-Z | `laol-core-table-z` |
| LAOL column: exact on a linear table, 1 dex past a row's ends | `laol-column-low-end` |
| LAOL column: beyond the tolerance returns ierr | `laol-column-tolerance` |
| LAOL column: log kappa above 35 is capped | `laol-cap` |
| LAOL column: a row of fewer than four points returns ierr | `laol-column-short-row` |
| blend_in_z: each table's value at its own Z | `blend-slope` |
| add_conduction: equal opacities halve kappa, average the slopes | `conduction-sum-slope` |
| hubbard_lampe_opacity: derivatives against differences | `hubbard-lampe-drho` |
| kap_init refuses the pure-Z table without an interior table | `settings-pure-z` |
| kap_init refuses two OPAL92 tables at one Z | `settings-two-z` |
| kap_init refuses molecular_opacity_logt_min >= _max | `settings-ramp` |
| kap_init: a missing table returns ierr (and is not created) | `missing-table` |
| kap_init: a missing OPAL92 table returns ierr (open_table) | `open-table-status` |
| kap_init: a short OPAL92 table returns ierr | `short-table` |
| kap_get: log T off the tables returns ierr, res unchanged | `opal95-t-range`, `res-on-error` |
| kap_get: Z above 0.12 below 10^7 K returns ierr | `z-limit` |
| Alexander 2006: a refused composition fails again on a retry | `alex06-cache` |
| kap_init again: the first lookup is the first run's, bit for bit | `kap-init-state` |
| kap_get: a helium core without the pure-Z table returns ierr | `helium-core-without-pure-z` |
| kap_init refuses an OPAL92 table that is not at opal_table_z1 | `opal92-table-z` |
| kap_init refuses a LAOL89 table that is not at laol_table_z1 | `laol-table-z` |
| kap_init refuses a laol_z_tolerance that is not positive | `laol-tolerance` |
| kap_init refuses a Kurucz 1990 table without rows | `kurucz-no-rows` |
| kap_init refuses an OPAL95 table at another X than its slot's | `opal95-header` |
| kap_init refuses an Alexander 2006 table at another X | `alex06-header` |
| kap_init refuses an Alexander 1994 table at another Z | `alex94-header` |
| kap_init refuses a LAOL89 table with more X than fit | `laol-size` |
| kap_init refuses a Kurucz 1990 table with an unreadable row | `kurucz-unreadable-row` |
| kap_init: an OPAL95 file that ends inside a table says so | `opal95-ends-inside` |
| kap_init: an OPAL95 file that ends after a whole table says so | `opal95-ends-after` |
| kap_init: an OPAL95 row that is not numbers is refused | `opal95-unreadable-row` |
| kap_init refuses an Alexander 2006 table with other log R columns | `alex06-logr-columns` |
| kap_init refuses an Alexander 2006 row at another log T | `alex06-row-logt` |
| kap_init refuses an Alexander 1994 row of another length | `alex94-row-length` |
| kap_init refuses a pure-Z table with more than one X | `pure-z-size` |
| Kurucz 1990 without a value: the interior opacity alone | `kurucz-out-of-table` |
| kap_get without an interior family says so, at any Z | `no-interior-family-first` |
| kap_report: a message before the run log opens waits for it | `run-log-through-log-line` |
| kap_set_envelope_metal_fraction: one Z table serves the new Z only | `envelope-z-follows` |

The pinned output is `kap/test/expected_test_kap.out`; a deliberate change means regenerating it from `src/test_kap` (run from the repository root with `YREC_INPUT=input`) and saying why in the commit. Several checks print an error at round-off; `test_kap.py` lets two numbers below 10⁻¹² count as equal in the tolerant mode. The decks pin the production combination end to end (`compare_decks.py`; every shipped deck uses OPAL 1995, Alexander 2006, the pure-Z table and Potekhin conduction): `test_kap` is the only coverage of OPAL92, LAOL89, Kurucz 1990 and Alexander 1994.

## 5. Coverage

92.4% of kap/'s 2061 executable lines run (2026-11-18), under `test_kap` and four configuration-matrix decks (`cm_base`, `cm_settle_he_z` for the 4D path under Z diffusion, `cm_atm_krishnaswamy`, `cm_core_overshoot`), measured with
```
YREC_INPUT=$PWD/input YREC_START=$PWD/startmodels python3 src/tools/coverage.py --work DIR --dir kap --test test_kap \
    --deck examples/run_config_matrix/cm_base.inlist --deck examples/run_config_matrix/cm_settle_he_z.inlist \
    --deck examples/run_config_matrix/cm_atm_krishnaswamy.inlist --deck examples/run_config_matrix/cm_core_overshoot.inlist
```
from the repository root (a -O0 --coverage build in DIR; `test_kap` needs `YREC_INPUT`: without it the test exits 1 and only the decks count). By file:

| File | Lines | Run | Coverage |
|---|---|---|---|
| `kap_lib` | 159 | 154 | 96.9% |
| `kap_messages` | 9 | 9 | 100% |
| `opacity_blending` | 41 | 41 | 100% |
| `opacity_families` | 25 | 25 | 100% |
| `opacity_table_files` | 33 | 32 | 97.0% |
| `opacity_tables` | 50 | 49 | 98.0% |
| `opal95_opacity` / `opal95_reader` | 295 / 108 | 285 / 89 | 96.6% / 82.4% |
| `opal92_opacity` / `opal92_reader` | 250 / 39 | 218 / 34 | 87.2% / 87.2% |
| `laol_opacity` / `laol_reader` | 225 / 90 | 209 / 90 | 92.9% / 100% |
| `alex06_opacity` / `alex06_reader` | 108 / 50 | 99 / 50 | 91.7% / 100% |
| `alex94_opacity` / `alex94_reader` | 108 / 41 | 104 / 39 | 96.3% / 95.1% |
| `kurucz_opacity` / `kurucz_reader` | 144 / 60 | 125 / 51 | 86.8% / 85.0% |
| `conductive_opacity` / `conductive_reader` | 201 / 25 | 177 / 25 | 88.1% / 100% |

(The figure of 2026-11-10 was 79.4% of 2175 lines, with four lookups between 72% and 77%.) The lines not run are failure branches and corners that no shipped or doctored table reaches: a spline that fails; a row that does not start at its first density or log R; a table that cannot be closed; the conductive kernel's own range refusals and its last-interval branch (Potekhin's `condopacp` and `cinterp3`); OPAL 1995's fixed-Z build at a Z at or above its top nodes, its reader's unreadable headers and a file cut before a table's log R grid; OPAL 1992's refusals of an X off its grid (unreachable: `findex` returns no index 0); Alexander 2006 at a Z above its top node or at X > 0.9 with Z ≥ 0.1; a Kurucz density exactly at a row's last; `kap_update_surface_tables`' error return (no table that `kap_init` accepts makes a refresh fail); and two messages for settings that `kap_init` refuses earlier.

## 6. How to extend

**Add a control** that kap reads: io/README.md §5 ("Add a control"; "Add a table path" for a table's file), with its field in `kap_settings_t` (or `kap_table_paths_t`) in `opacity_tables.f90` and its copy in `kap_settings_from_star` (or `kap_table_paths_from_star`) in `setup/kap_settings_builder.f90`.

**Add a table family** (traced on 2026-11-10 on a throwaway branch with a toy molecular family, `flat`, a constant opacity read from a one-number file: the files below are the ones it touched; steps 4's last line and 6's messages added on 2026-11-18 from the eighth evaluation's reading, M344):
1. Its directory, `kap/<family>/`, with two modules:
   - `<family>_reader.f90`: `read_<family>_table(tbl, path, ierr)`, opening with `opacity_table_files`' `open_table` and checking every read with `read_failed` (and the file's Z with `table_z_differs`, if a control names it);
   - `<family>_opacity.f90`: the lookup (`compute_<family>_opacity`, `pure`, the table set and the point as arguments, no I/O statement, array dummies assumed-shape), and the descriptor `<family>_source(settings)` with its procedures behind `opacity_source`'s interfaces: `load` (the reader, then whatever it builds), `lookup`, and where needed `refresh` (a surface table) and `core_lookup` (an interior family, for helium-burning cores). Set `covers_any_z`, or `blends_two_z` with `table_z` and `table_z_controls`, as the family is.
2. `kap/opacity_tables.f90`: its `<family>_table_set` type, its component in `opacity_table_state` (in order of size: on arm64 Mach-O a component must start within 8 MB of the state's start), its `use_<family>_tables` field in `kap_settings_t` and its path field in `kap_table_paths_t`.
3. `kap/opacity_families.f90`: its `use` line and its entry in its role's list, at its place in the precedence (the list's length grows by one).
4. The controls, by io/README.md §5 ("Add a control", "Add a table path"): their declarations in `io/controls_buffer.f90`, their rows in `defaults/controls_registry.tsv` (settings cells `kap_settings%use_<family>_tables` and `kap_table_paths%<family>`), and their copies in `kap_settings_from_star` and `kap_table_paths_from_star` (`setup/kap_settings_builder.f90`). **A molecular family also goes in `io/derived_controls.f90`'s `select_eos_and_opacity_tables`**, which counts the molecular controls that are set, warns for more than one and switches off the ones that lose, in the precedence order: its count, its warning's text and its switch-offs.
5. `src/tools/check_boundaries.py`: the lookup in `STAR_BLIND_FILES`. `src/tools/style_counts.tsv`: the new directory's `io_statement` count, the reader's opens and reads (a new directory starts at 0, and the reader may keep them; a raise is a hand edit the commit explains).
6. **An interior family's control also goes in kap_lib's two messages that list the interior controls by name** (`check_kap_settings`' "use_pure_z_table needs an interior table (...)" and `no_interior_family`'s "no interior opacity table chosen (...)"; both are pinned in `expected_test_kap.out`).
7. A check in `test_kap` (its nodes, its derivatives against differences, a doctored copy of its table for its reader's refusals), its mutation in `kap/test/test_kap.mutations.toml`, and this README's rows in sections 1 and 4.
8. From `src/`: `python3 tools/gen_makefile_deps.py` (two new modules); `make -j8` (it reruns the controls generators: `defaults/control.defaults`, `defaults/controls_registry.md`, `io/controls_buffer_public.inc`, `io/controls_sync.f90`, `io/expand_paths.inc`, `io/inlist_decl.inc`, `io/test/controls_dump.inc` and `state/control_group_def.inc` change); `python3 tools/check_star_contracts.py --fix` (the builder's `star reads:` lines); `make check`; `make tests`; from the repository root `YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py` (`io/test/expected_controls/pristine.txt` gains the two controls) and `pytest test_kap.py` (`kap/test/expected_test_kap.out` gains the check's line); `python3 src/tools/mutate.py kap`.

So a family is its directory (a lookup and a reader), its table-set type, one list entry, two settings fields and two controls: nine hand-edited files besides the test and the README (the two modules, `opacity_tables`, `opacity_families`, `controls_buffer`, the registry, the builder, `check_boundaries`, `style_counts.tsv`), a tenth for a molecular family (`io/derived_controls.f90`) or an interior one (`kap_lib`'s messages), and the generated ones make writes.

**Add a check:** a call in `test_kap`'s `run_checks` (the helpers `node_check`, `opal95_node_check`, `derivative_check`, `continuity_check`, `surface_table_check` and `report` print the line; `copy_table` writes a copy of a table with one line edited, `kap_init_with` loads one configuration with shipped tables but the one named, and `last_log_line` reads the last message), its mutation in `kap/test/test_kap.mutations.toml` (`python3 src/tools/mutate.py kap --check`, then the run), its line in `kap/test/expected_test_kap.out` and its row in section 4's table. A check of a refusal should read the refusal's message (`last_log_line`): a short or doctored table is often refused for a second reason once the first is mutated away.

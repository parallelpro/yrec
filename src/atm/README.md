# atm/ — the atmosphere

This directory gives the stellar model its surface boundary condition: the pressure at the photosphere, either from a table of model atmospheres or from an analytic T-tau relation that the star layer's envelope integrator integrates. Code outside atm/ uses it only through `atm_lib` (enforced by `tools/check_boundaries.py`, `FACADE_ONLY_DOMAINS`). atm/ reads nothing from `star`: its configuration is one settings record and one paths record, installed at the start of every run.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `atm_lib` (`atm_lib.f90`) | **The public interface**, and the run's settings and tables (private module variables that only `atm_init` writes). `atm_init` checks the settings, installs them and reads the table they select. `atm_get_surface_pt` looks up log P (and log T) at (log Teff, log g) in the table of the caller's current choice and reports what the lookup returns. Its other entries: the choice queries (`atm_choice_is_known`, `atm_choice_is_tabulated`, `atm_choice_name`, the named choices and `atm_last_choice`), the switch to gray (`atm_leaves_table`, `atm_reenters_table`, `atm_note_choice_switch`), `atm_allard_teff_max`, `atm_table_is_loaded` (whether `atm_init` read a choice's table), and the T-tau relations below. |
| `atm_choices` (`atm_choices.f90`) | **What each `atm_choice` is**: the named choices (`atm_eddington_gray` = 0 ... `atm_castelli_kurucz` = 5), and for each its name in the run log, whether it is a table, and whether the star layer replaces it by the gray relation above a log Teff (`kurucz_gray_switch_log_teff` = 3.95 for Kurucz; the Allard tables' last log Teff for Allard; never for the others). |
| `atm_tables` (`atm_tables.f90`) | The types, no state: the settings record `atm_settings_t`, the paths record `atm_table_paths_t`, a surface-pressure table `surface_pressure_table` and the Allard tables `allard_table`, with their sizes. |
| `atm_messages` (`atm_messages.f90`) | Every message atm writes, one procedure each, with its destinations (§3). |
| `t_tau_relations` (`t_tau_relations.f90`) | The analytic T-tau relations, Eddington gray, Krishna-Swamy (1966) and the HSRA fit: `ttau_log10_temperature`, `ttau_start_log10_temperature` (the start of the integration), `ttau_photosphere_x_limit` (its end) and `hsra_t_tau_offset`, with their constants named (a star-blind kernel). |
| **tables/** | |
| `kurucz_atmosphere` | `interpolate_surface_pressure`: log P at (log Teff, log g) from a Kurucz (1993) or Castelli & Kurucz (2003) table: natural cubic splines in log g along four log Teff rows (3-point Lagrange or linear in the top columns above log Teff 4.5), then in log Teff. The table in, log P and a status out (a star-blind kernel). |
| `kurucz_table_reader` | `read_surface_pressure_table`: a surface-pressure file into a `surface_pressure_table`, and each row's first and last valid log g column. |
| `allard_atmosphere` | `interpolate_allard_surface`: the Allard tables (NextGen, AMES-Cond, AMES-Dusty) at (log Teff, log g), in named steps: the tables' range, the 4x4 stencil, 4-point Lagrange in log g along each row (`polint`), then across the rows in log Teff. The tables in, an `allard_lookup_result` out: the result, or which limit the point is beyond (a star-blind kernel). |
| `allard_table_reader` | `read_allard_tables`: an Allard file in either layout (the old NextGen 1 file; the 2007 layout at the deck's [Fe/H] and [α/Fe], in two passes), then `find_table_ranges`, the lookups' limits and the table's check. |

The atmosphere integration itself (`core/envelope/envelope_integration_kernel.f90`, `core/envelope/envelope_integrands.f90`) belongs to the star layer. It calls these services, and its results (`atmosphere_structure`, `envelope_structure`, the current point `atm_point`) live in core/.

## 2. The flow

**At the start of a run** (`setup/settings_install.f90`'s `install_run_settings`, called by `core/run/run_setup.f90`):
1. `atm_init(atm_settings_from_star(s), atm_table_paths_from_star(s), ierr)`; `setup/atm_settings_builder.f90` builds both records from the controls.
2. `atm_init` refuses a choice that is not one of the six (`atm_choice_is_known`), installs the settings, returns every table to its unread state (so that a second run in one process starts clean, M283) and reads the choice's table: 3 `read_surface_pressure_table` (the Kurucz table, `atm_table_file`), 4 `read_allard_tables` (`allard_atm_table_file`, at `allard_target_feh` and `allard_target_alpha`), 5 `read_surface_pressure_table` (the Castelli/Kurucz table, `atm_table_file`). Choices 0-2 read nothing.

**At each envelope integration** (`core/envelope/envelope_integration_kernel.f90`):
1. `prepare_surface_boundary`: when `atm_choice_is_tabulated(choice)`, `atm_get_surface_pt` gives log P and log T at the surface. It calls the choice's kernel (`interpolate_surface_pressure` or `interpolate_allard_surface`) and reports what it returns through `atm_messages`. An Allard point above the tables in log Teff or outside them in log g is `lookup_failed`, and the kernel integrates the Eddington relation for that integration instead.
2. `integrate_atmosphere`, when no table gave the surface: the T-tau relation from `ttau_start_log10_temperature` to `ttau_photosphere_x_limit`, with `ttau_log10_temperature` at each step (`core/envelope/envelope_integrands.f90`). `hsra_t_tau_offset` goes into every integration's configuration (`core/envelope/envelope_integration.f90`).

**At each surface boundary** (`core/structure/surface_boundary.f90`'s `switch_atmosphere_choice`): `atm_leaves_table` says when the run's working choice (`star%atm_choice_in_use`) leaves its table for the gray relation, `atm_reenters_table` when the deck's choice takes the surface back, and `atm_note_choice_switch` writes the terminal's note of each switch.

**At accretion** (`core/step/mass_change.f90`'s `find_accreted_surface_state`): when `atm_table_is_loaded(atm_allard)` (the deck's choice is 4), `atm_get_surface_pt(atm_allard, ...)` at the star's Teff and log g; otherwise, or when the point is outside the tables, core takes the photosphere of the last envelope integration (O19).

**In the run log's header** (`io/derived_controls.f90`): `atm_choice_name` gives the "surface boundary:" line.

## 3. Conventions

- **Units**: log is log₁₀; P in dyn/cm², T in K, g in cm/s². A table lookup returns log P at T = Teff with log T = log Teff, except Allard with `allard_use_tau100`, which returns log P and log T at tau = 100. The outputs of `atm_get_surface_pt` are set only when the lookup succeeds.
- **Errors** return `ierr` (positive) and never stop the run:
  - `atm_init`: 1 for a choice that is not one of the six, or a table that cannot be opened, read or used (a table without a 4-point stencil, a -999 inside a row, a row of fewer than four valid columns, an old NextGen 1 file with `allard_use_tau100`, whose tau = 100 columns it lacks: O19); 2 (`allard_composition_missing`) for an Allard file without the deck's [Fe/H] and [α/Fe], with a message that lists the file's;
  - `atm_get_surface_pt`: 1 for a point outside a Kurucz table (below log Teff 3.5 or log g -0.5, more than one step above the last row or column, or a stencil reaching a row without the columns it needs) or below the Allard tables, for a choice that is not a table, or for a table `atm_init` did not read (O19); a spline's or `polint`'s `numerics_status` when it fails. A point above the Allard tables in log Teff, or outside them or a stencil row in log g, is not an error: `lookup_failed`, and the caller falls back to gray;
  - the kernels return their status (`surface_pressure_outside`, the numerics codes, the `allard_lookup_result` flags), and `atm_lib` reports it.
- **Messages** go through `atm_messages`, one procedure each; kernels print nothing. The texts and destinations are still the F77's and the port's, byte for byte: capitals and F77 names ("ALSURFP: TEFFL greater than TeffLmax"), list-directed lines with their trailing blanks, the "Program Terminated" banners. The terminal is standard output, except the Kurucz range message and the 2026 checks' errors (`report_atm_error`, run_messages' `report_error`), which go to run_messages' `terminal_unit`; most messages also go to the run log. A successful lookup writes nothing (M235). Rewording them is the B phase's message commit.
- **Controls**: atm reads its settings record and its paths record, built in `setup/atm_settings_builder.f90`, never `s%control` or `s%physics` (`STAR_FREE_DIRS`).
- **Kept by decision** (the F77's, named; changing any of them moves decks):
  - the photospheric log tau to nine decimals (the ¼ log₁₀(3/4) of the Eddington and Krishna-Swamy relations is computed, not the F77's -0.031235, which put T 3.2×10⁻⁷ in log T off Teff at the photosphere);
  - the Allard reader's tolerances: [Fe/H] within 10⁻⁵ of the deck's, [α/Fe] exactly equal;
  - the Allard tables' low-gravity extension: four column widths below a row's first valid log g, one above its last (M233);
  - the Kurucz table's switch to gray at log Teff 3.95, below the table's top (MHP 9/01).

## 4. Tests

`test_atm` (`make tests` in `src/`, then `pytest test_atm.py` from the repository root; it reads the repository's tables through `YREC_INPUT`) has two parts.

**The pins**, through `atm_lib`, compared byte for byte with `atm/test/expected_test_atm.out`: the lookups of every table (Kurucz, Castelli/Kurucz, Allard NextGen 2 with and without the tau = 100 point, the old NextGen 1 layout, a lookup above the Allard tables), the three T-tau relations on a grid of optical depths, and a Kurucz lookup below the table.

**The checks**, 39 of them, each printed with its error (or 0/1 for a contract) and `ok` or `FAIL`; the program exits 1 on a FAIL. `test_atm.py` compares a check's line by its label and status, since its error is round-off that differs between platforms.
- **The T-tau relations**: T = Teff at each photosphere (Eddington's tau = 2/3 and HSRA's to round-off, Krishna-Swamy's 0.312156 to 10⁻⁹ in log T, its log tau's nine decimals); the photospheric log tau of Eddington and HSRA is log₁₀(2/3); Eddington's and Krishna-Swamy's start is their T as tau -> 0; Eddington's T⁴ is linear in tau + 2/3 (the gray law); T rises with tau from 10⁻³ to 10 in each relation (the HSRA fit has the solar temperature minimum near tau = 3×10⁻⁴); the HSRA fit's value at tau = 1; NaN for a table's choice.
- **The Kurucz and Castelli/Kurucz lookups** (`interpolate_surface_pressure` on tables read by `read_surface_pressure_table`): exact at every node below log Teff 4.5; continuous across the log Teff rows, where the stencil moves; exact on a plane in log Teff and log g on the Castelli grid in every branch (splines, 3-point Lagrange, linear); a point more than a step above the table refused (M232); a row ending below the top column (atmk1990m50) exact at its last valid node, and a stencil on an all -999 row (atmk1990p01's row 29) refused (M234); a Castelli/Kurucz stencil reaching the last row refused (M231), and a 3-point row without its top three columns; the reader's errors, a missing and a truncated file.
- **The Allard lookups and reader** (`interpolate_allard_surface`, `read_allard_tables`): exact at the nodes of Nextgen2; continuous across the rows; T(tau = 100) above Teff; the outcome's flags (above, outside, below). On synthetic files in the Allard layout that the test writes (rows out of order, a record of another composition to skip): exact on a polynomial cubic in log g and quadratic in log Teff, to the file's nine digits; the low-gravity extension of four column widths and one above (M233); a -999 inside a row and a row of three valid columns refused. The reader's errors: an [Fe/H] or [α/Fe] the file lacks is `allard_composition_missing` (M221); three rows; a missing, empty or unrecognised file; a record whose log P cannot be read; the old NextGen 1 file asked for the tau = 100 point (O19).
- **atm_lib's contracts**: a lookup in a table `atm_init` did not read returns ierr 1 (checked first, before any table is read; the Allard tables too, O19); `atm_table_is_loaded` is true for the table `atm_init` read and for no other, and for none before (O19); a choice that is none of the choices refused (M236); each `atm_init` starts clean, without the tables an earlier one read (M283); a choice that is not a table returns ierr 1; below the Allard tables ierr 1, and off a stencil row's log g `lookup_failed` with the outputs untouched; the descriptor's known, tabulated and named choices; Kurucz and Allard switch to gray at their limit and back below it, the other choices never; the four switch notes (pinned).

Each check was confirmed to fail under a mutation of the code it covers: `atm/test/test_atm.mutations.toml` holds 45, and `python3 src/tools/mutate.py atm` (after `make tests`) reports 45 caught, 0 survived (2026-11-11). `pytest test_atm.py` also checks that every mutation still finds its text. The program runs clean under the development flags (`make devel_tests`, `python3 -I src/tools/devel_tests.py`: test_atm ok, 2026-11-10).

A deliberate change to a pin or to a check's label means regenerating `atm/test/expected_test_atm.out` (run `src/test_atm` with `YREC_INPUT` set, in a scratch directory, and copy its output) in a commit that says why.

The decks pin the whole surface boundary: the gray decks (most of them), `cm_atm_krishnaswamy` (Krishna-Swamy), `Test_solarGS98_kurucz` and `Test_solarM22M_difkur` (Kurucz), `Test_solarGS98_allard` (Allard), and `cm_accretion` (an accreting gray deck, whose accreted surface is the envelope integration's photosphere since O19's fix: it asks atm for no table). No shipped deck uses HSRA or Castelli/Kurucz, and none crosses a table's switch to gray; `test_atm` is their only coverage.

## 5. Coverage

95.8% of atm/'s 708 executable lines (2026-11-10), from the repository root:

    python3 src/tools/coverage.py --work DIR --dir atm --test test_atm \
        --deck testsuite/Test_solarGS98_allard.inlist --deck testsuite/Test_solarGS98_kurucz.inlist \
        --deck src/tools/compare_decks/cm_accretion.inlist \
        --deck examples/run_config_matrix/cm_atm_krishnaswamy.inlist

By file: atm_choices, t_tau_relations and allard_atmosphere 100%, atm_lib 98.9%, kurucz_table_reader 98.0%, kurucz_atmosphere 97.1%, atm_messages 93.6%, allard_table_reader 92.9%. The 30 lines not run:
- a spline's or `polint`'s failure, and its message (`kurucz_atmosphere`, `atm_lib`, `atm_messages`): their only failure is coincident nodes, which a table's sorted, distinct grid cannot give;
- the Allard second pass's "Impossible failure" records (a Teff or log g the first pass did not collect) and their message: the two passes read the same records with the same filter, so the branches only guard against a file that changes between them;
- the read errors of the old NextGen 1 layout, of the Allard layout's header and of its first pass (`allard_table_reader`); the test's unreadable record fails in the second pass;
- one line of `read_surface_pressure_table`'s allocation, which gcov counts on a continuation line though the statement runs.

## 6. How to extend

The two recipes for a new choice were traced on throwaway branches on 2026-11-10: the change made, built, `make check` and the tests run, a deck run with the new choice; the files listed are the files the traces touched.

**Add a T-tau relation** (choice 6, say). Nothing in core/ or io/'s code changes: the envelope integrates any choice that is not tabulated, the banner takes the choice's name, and `atm_init` accepts it once it is known.
1. `atm/atm_choices.f90`: its named choice (`atm_<name> = 6`) as `atm_last_choice`, and its row in `choices`: its name for the run log, `.false.` (not a table), `.false.` (no switch to gray).
2. `atm/t_tau_relations.f90`: a case in each of the three functions (`ttau_log10_temperature`, `ttau_start_log10_temperature`, `ttau_photosphere_x_limit`), its constants named at the top, its formula in the header.
3. `atm/atm_lib.f90`: re-export `atm_<name>` (its `use atm_choices` list and a `public` line).
4. The control's description: the comment above `atm_choice` in `io/controls_buffer.f90`; `make` rewrites `defaults/physics.defaults` and `defaults/controls_registry.md` from it.
5. `atm/test/test_atm.f90`: a check that T = Teff at its photosphere (and any limit it keeps), its mutation in `atm/test/test_atm.mutations.toml`, then regenerate `atm/test/expected_test_atm.out` (it gains the check's line).
6. This README: §1's row and §4's list.

**Add a table** (choice 6, a surface-pressure table in the Kurucz layout with a path of its own and a switch to gray, as traced). Nothing in core/ changes: the envelope kernel asks `atm_choice_is_tabulated`, the surface boundary `atm_leaves_table` and `atm_reenters_table`, the banner `atm_choice_name`.
1. `atm/atm_choices.f90`: its named choice as `atm_last_choice`, its row (`.true.`, a table; `.true.` if the star layer switches it to gray), and its switch's log Teff as a named parameter if it has one.
2. `atm/atm_tables.f90`: a field in `atm_table_paths_t` for its file. (A table in another layout also needs its type here, with its sizes.)
3. `atm/atm_lib.f90`: re-export the choice (the `use atm_choices` list and a `public` line); a module variable for its table; its case in `atm_init` (the reader, with its row count) and in `atm_get_surface_pt` (`look_up_surface_pressure`, or a lookup procedure of its own that reports its kernel's outcome); its case in `gray_switch_log_teff` if it switches to gray.
4. `atm/atm_messages.f90`: its note in `note_choice_switch` if it switches to gray (and the messages of a table in another layout).
5. `setup/atm_settings_builder.f90`: the path's copy in `atm_table_paths_from_star`; then `python3 tools/check_star_contracts.py --fix` in `src/` (its `star reads:` line).
6. The path control, as `io/README.md` §5 "Add a table path" says: its declaration and description in `io/controls_buffer.f90` (and the new choice in `atm_choice`'s description), its row in `defaults/controls_registry.tsv` with the settings cell `atm_table_paths%<field>`. `make` rewrites the generated files (the trace's: `defaults/control.defaults`, `defaults/physics.defaults`, `defaults/controls_registry.md`, `io/controls_buffer_public.inc`, `io/controls_sync.f90`, `io/expand_paths.inc`, `io/inlist_decl.inc`, `io/inlist_read.inc`, `io/test/controls_dump.inc`, `state/control_group_def.inc`); then the re-pins that section lists (`test_controls`' `pristine.txt`).
7. `atm/test/test_atm.f90`: the choice in the descriptor check's list of tables, and checks of its lookup (exact at the nodes, its refusals) with their mutations; regenerate `atm/test/expected_test_atm.out`.
8. This README: §1, §2 and §4.

A table in another layout adds, before step 3, a reader in `atm/tables/<name>_table_reader.f90` (reporting through `atm_messages`) and a pure kernel in `atm/tables/<name>_atmosphere.f90` (the table in, the result and a status out, printing nothing) listed in `tools/check_boundaries.py`'s `STAR_BLIND_FILES` and `KERNEL_FILES`; then `python3 tools/gen_makefile_deps.py` in `src/`.

**Add a setting** (a control atm reads, as `allard_use_tau100`): the control as `io/README.md` §5 "Add a control" says, with the settings cell `atm_settings%<field>`; the field in `atm_settings_t` (`atm/atm_tables.f90`) and its copy in `atm_settings_from_star` (`setup/atm_settings_builder.f90`), then `check_star_contracts.py --fix`; `atm_lib` reads `installed_settings%<field>` and passes it to the reader or the kernel as an argument (they read no settings); a check in `test_atm`. The table recipe's path field traced the same steps for the paths record.

**Add a message**: a procedure in `atm/atm_messages.f90` that says it, with its destinations, called by `atm_lib` or a reader; never a write in a kernel (`check_boundaries.py`'s `check_no_io`).

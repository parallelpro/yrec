# eos/ — the equation of state

This directory computes the density of a zone, its thermodynamic quantities (beta, the ionization fractions, the mean molecular weights, cp, the adiabatic gradient) and their derivatives, from (log T, log P, X, Z). Code outside eos/ uses it only through `eos_lib` (`tools/check_boundaries.py` lists its three entries).

eos/ reads nothing from `star`: it is in `STAR_FREE_DIRS`, and its solvers in `STAR_BLIND_FILES`. Its configuration is one settings record and one record of table paths, both built by the star layer and handed to `eos_init`, and the envelope mixture the star layer hands `eos_set_mixture`. Its constants are compile-time parameters of `state/physical_constants.f90`, so it gives the right answer without any other setup having run.

The equation of state is the Yale gas (Saha below a cutoff temperature, fully ionized above it, blended across a ramp) with **table sources** applied over it: SCV at low temperature, an OPAL table, or MHD instead of everything. Each source is described once (`eos_table_source`) and listed once (`eos_source_list`), so a new source is a module, one list entry and its settings (§6).

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `eos_lib` | **The public interface.** `eos_init` installs the settings and reads the tables they select; `eos_get` evaluates a zone into an `eos_result`; `eos_set_mixture` takes the envelope mixture (and has SCV's envelope table built for it). It re-exports the settings and paths types and the `opal_vintage_*` constants. |
| `eos_configuration` | `eos_settings_t` (everything eos reads from the run's controls, with the controls' defaults) and the installed `eos_settings` (`set_eos_settings`, called by `eos_init`); `eos_table_paths_t`; the `opal_vintage_*` constants; `saha_ramp_width`. |
| `eos_results` | `eos_result`, the record every source fills: one named field per quantity, with its unit, and the fields that are also inputs (§3). |
| `eos_mixture` | The envelope mixture `eos_mix` (protected: `set_eos_mixture_state` writes it) and its 12 species' atomic weights; the composition terms' constants; `ion_mean_weight_excess`, `compute_mean_weights`, `mean_weights_from_density`; `debye_huckel_composition_of`, the Debye-Hückel correction's composition at a zone. |
| `eos_messages` | `eos_log`, `eos_report` and `eos_log_once`, eos's messages (§3), and the once-a-run warnings' flags (`reset_eos_warnings`, `eos_warned`). |
| `eos_table_source` | The table-source descriptor `eos_table_source_t` (name, stage, and the procedures `selected`, `read_tables`, `in_window`, `lookup`, `edge_weight`, `adopt`, `blend`), the stages, and the standard procedures a source names: `always_in_window`, `whole_weight_inside`, `adopt_whole_result`, `adopt_table_result`, `blend_table_result`. |
| `eos_source_list` | `table_sources()`, the ordered list of sources: SCV, OPAL, MHD. |
| `blended_eos` | The EOS of a run without a whole-EOS source (§2.1): `blend_eos` in named steps (`zone_mixture`, `radiation_share`, `gas_eos`, then the over-the-gas sources through `apply_table_source`), with the Saha ramp's `blend_saha_into` and `fully_ionized_gas`; `compute_blended_eos` adds the numerical second derivatives. |
| `yale_eos` (`yale/`) | The Yale EOS of Prather (1976): `compute_saha_eos` (the Saha equations of 11 species and He's two stages, by Newton iteration, in steps: `set_saha_constants`, `solve_saha_system` with `charge_balance` and `eliminate_helium`, the two derivative steps) and `compute_fully_ionized_eos` (degenerate, relativistic electrons from the Fermi-Dirac table, the optional Debye-Hückel correction; steps `start_fully_ionized_solution`, `iterate_fully_ionized_density`, `interpolate_fermi_table`, `differentiate_pressure`, `differentiate_density`, `compute_internal_energy_and_cp`). |
| `fermi_table_reader` (`yale/`) | `read_fermi_table`, the Fermi-Dirac table's reader. |
| `scv_eos` (`scv/`) | The Saumon-Chabrier-Van Horn EOS: `scv_table_source` (its descriptor entry), `build_scv_envelope_table`, `compute_scv_envelope_eos` (the envelope table, for a zone at the X and Z it was built for, to 1e-5) and `compute_scv_eos` (any X, Z, from the three tables), which share `find_scv_cell` and `scv_ionization_fractions`. |
| `scv_table_reader` (`scv/`) | `read_scv_tables`, the H, He and metal tables' reader. |
| `opal_eos` (`opal/`) | The OPAL EOS of the 1995, 2001 and 2006 vintages, written once and driven by the vintage's row of `opal_layout`: `opal_table_source`, `compute_opal_eos` (the density that gives P, and the derivatives there), `find_opal_edge_ramp`; under them `opal_density_search`, `interpolate_opal_eos` (OPAL's `esac`), its coverage sums (`clamped_coverage`; 1995: `unclamped_coverage`, `stencil_entry`, §3), `t6rinteos`, the radiation corrections `add_radiation` (2001, 2006) and `add_radiation_to_gammas` (1995), `opal_mass_per_mole`, the edge ramps `find_table_edge_ramp` (1995, 2001) and `find_ragged_edge_ramp` (2006), and the messages in each vintage's words (`report_*`). |
| `opal_table_reader` (`opal/`) | `open_opal_eos_table` (at `eos_init`) and `readcoeos`, the reader `interpolate_opal_eos` calls on a vintage's first use, with the reader's messages. |
| `mhd_eos` (`mhd/`) | The Mihalas-Hummer-Däppen EOS: `mhd_table_source`, `compute_mhd_eos`, Däppen's interpolation `mhdpx`, `mhdpx1`, `mhdpx2`, and `interpolate_mhd_table`, the 4x4 lookup in one table. |
| `mhd_table_reader` (`mhd/`) | `read_mhd_tables` and Däppen's readers `mhdst`, `mhdst1`, `rabu`, `rtab`. |
| `yale_eos_tables`, `scv_eos_tables`, `opal_eos_tables`, `mhd_eos_tables` | Each family's tables and state (§5). `yale_eos_tables`: the Fermi-Dirac table and its shape (`fermi_n_*`, the rows' log T, §3), the Debye-Hückel constants and composition record. `scv_eos_tables`: the three tables, the envelope table with the composition it was built for, the search start `last_scv_cell`. `opal_eos_tables`: the instances `opal_tables(iv_opal95:iv_opal06)` of one state type, `opal_shared` (the open table's unit, the loaded flags, each table's Z), `opal_layout` (one row per vintage: everything the code tells the vintages apart by), the empty-entry fill and its thresholds, and the helpers `interpolate_opal_t6_rho`, `find_opal_ramp`, `quad`. `mhd_eos_tables`: the MHD tables and the last interpolation. |
| `eos_test_entries` | What `test_eos` reaches beyond `eos_lib`, in one place: `gamma1_at_density` (Γ₁ and ∇_ad at T6, rho), the OPAL entry points `esac95`, `esac01`, `esac06`, and the internals the checks inspect, re-exported. No production module carries a procedure for the test. |

`eos/retired/opal95_radsub_notes.f90.txt` keeps the commented-out lines and the open question of the 1995 radiation correction (§3).

## 2. The flow

The star layer calls eos at three times:

1. **Once a run**, `install_run_settings` (`setup/settings_install.f90`, from `core/run/run_setup.f90`) calls `eos_init` with `eos_settings_from_star` and `eos_table_paths_from_star` (`setup/eos_settings_builder.f90`). `eos_init` checks the settings (one OPAL vintage, eta_max > eta_min, a Saha cutoff of at least 5.7), installs them, resets the once-a-run warnings, reads the Fermi-Dirac table, then each selected source's tables in the list's order: SCV's three tables, the OPAL table (opened here; its reader reads it on first use), the MHD tables.
2. **Whenever the envelope composition changes** (the run plan's envelope Z, the starting model, the rescaling of Z), `install_surface_z_settings` calls `eos_set_mixture`, which stores the mixture and, with SCV on, builds SCV's envelope table for its X and Z (once the tables are read; a table built for them already stays).
3. **For every zone evaluation**, `eos_get`: a selected whole-EOS source's lookup alone (MHD, §2.2); otherwise the Debye-Hückel composition of the zone (when the correction is on) and `compute_blended_eos`, which calls `blend_eos` once or, with `want_derivatives` and `use_numerical_derivatives`, four more times at log T ± 3e-4 and log P ± 1.5e-3 for centred-difference second derivatives.

### 2.1 blend_eos

For one (log T, log P, X, Z), with the cutoff T_c = `saha_log10t_cutoff` and the ramp width `saha_ramp_width` = 0.5 dex:

1. **The zone's mixture** (`zone_mixture`): 1/mu_ion, and the Saha equations' abundances, from the envelope mixture. With `use_diffusion_z` all the metals scale with Z/Z_env; otherwise (and when the envelope has no metals) the zone differs from the envelope by its excess H, He and metals (the metals counted as carbon).
2. **Radiation** (`radiation_share`): beta = P_gas/P = 1 − (a/3) T⁴/P; beta ≤ 0 is an error.
3. **The gas** (`gas_eos`), by log T:

   | log T | the gas |
   |---|---|
   | ≥ T_c | fully ionized (`fully_ionized_gas`: 1/mu_e from X, then `compute_fully_ionized_eos`) |
   | T_c − 0.5 < log T < T_c | the Saha result (with the in-place-of-Saha sources applied) blended into the fully ionized one, weight f = w³(6w² − 15w + 10), w = 2(T_c − log T) (`blend_saha_into`) |
   | ≤ T_c − 0.5 | the Saha equations (`compute_saha_eos`), with the in-place-of-Saha sources applied |

   The in-place-of-Saha source is SCV: where its tables cover the point its lookup, from a copy of the Saha result, replaces it whole (`adopt_whole_result`); it has no edge weight (M88, a measured change for a later round).
4. **The over-the-gas sources** (OPAL), each in turn (`apply_table_source`): inside the source's window, its lookup from a copy of the result; T and P become the lookup's; where the table holds the point its result replaces the result (`adopt_table_result`), and within one table cell of its edge it is blended in with the edge weight (`find_opal_edge_ramp`, `blend_table_result`); elsewhere the result stands.

So a deck with SCV and OPAL 2006 (304 of the 310 decks in testsuite/ and examples/, the reference solar decks among them) gets OPAL wherever its table covers the point (log T from 3.27 up), SCV at low T outside it, the Saha solution in between, and the fully ionized gas at high density and in the centre above the OPAL table.

**Which source fills which field** (the per-path table; "kept" is the value the step before gave):

| Field | Saha | fully ionized | Saha ramp | SCV | OPAL (adopted; blended linearly near the edge) | MHD |
|---|---|---|---|---|---|---|
| density | composition, P_gas/(R/mu T) | Fermi-Dirac iteration | blend in log rho | SCV table | OPAL's density search | table |
| beta, beta14 | radiation_share | radiation_share | kept | kept | OPAL (1 − P_rad/P), beta14 from it | table (P_gas/P) |
| R/mu, 1/mu_e | composition | composition | beta P/(rho T) | beta P/(rho T) | `compute_mean_weights`: beta P/(rho T) | `compute_mean_weights` |
| 1/mu_ion | `zone_mixture` | `zone_mixture` | kept | kept | excess metals as carbon (M129) | excess metals as carbon |
| ionization fractions | Saha | 1, 0, 1 | blended | SCV's H and He tables | kept (OPAL gives none) | table |
| eta | 0 | Fermi-Dirac | the fully ionized one's, not blended | kept | kept | table |
| first derivatives | analytic | analytic | blended | SCV table, with radiation | OPAL table | table |
| second derivatives | analytic (with `want_derivatives`, not `in_atmosphere`) | analytic | blended | Saha's kept | kept (Yale/SCV's) | from the table's slots, always |
| Z | the zone's (mixture) | the zone's | — | the zone's (envelope table: the envelope's) | **the table's own Z** (the zone's enters only 1/mu_ion) | the tables' |

With `use_numerical_derivatives` the second derivatives are centred differences of the blend's own first derivatives, on every blended path.

### 2.2 MHD

`compute_mhd_eos` interpolates the MHD tables (three ZAMS-type compositions, each in a lower and an upper part in T, and five centre-type ones) bicubically in (log P, log T) in each table of the family the temperature selects, then in X (quadratically across the ZAMS compositions, cubically across four of the centre ones). A point outside the tables returns ierr = 1. The tables are not shipped. Its input is the total log P, while Däppen's header names the gas pressure (an open question, M157: with real tables, 10^(slot 2) − 10^(slot 20) = aT⁴/3 at a node would settle it).

## 3. Conventions

- **Units** are cgs (T in K, P in dyn/cm², rho in g/cm³, cp in erg/(g K)); log is log₁₀. OPAL works internally in T6 = T/10⁶ and P in 10¹² dyn/cm².
- **Derivatives are logarithmic**: ∂lnρ/∂lnT at constant P, ∂lnρ/∂lnP at constant T, and the six second derivatives `dlnrho_dlnt_dt` (∂(∂lnρ/∂lnT)/∂lnT), `dlnrho_dlnt_dp`, `adiabatic_gradient_dt`/`_dp` (∂ln∇ad/∂lnT, ∂lnP), `specific_heat_cp_dt`/`_dp`.
- **`want_derivatives`** asks for the second derivatives; on the blended path, without it they keep the caller's values (MHD sets them always). **`in_atmosphere`** skips only the Saha solver's second derivatives.
- **Inputs carried in the result**: the density is the starting guess of the fully ionized iteration, and the `saha_*` fields the Saha solver's start (it returns its solution in them), so a result depends, within those iterations' tolerances, on the record it was given. A fresh record starts the fully ionized iteration at log rho 0, from which it did not converge at log T 8.4 to 9.35 and log P 21 to 23.5 where it was tried (an error now, M93; from a start near the answer it converges). `saha_iterations` is a counter nothing in eos reads.
- **The mean weights**: where the density comes from the composition, R/mu = R(1/mu_ion + 1/mu_e); where it comes from a table or a blend, R/mu = beta P/(rho T) and 1/mu_e is backed out of it; so R/mu steps by 7e-4 at the Saha cutoff (an open question to the authors, docs/bug_sweep.md).
- **The Fermi-Dirac table's range** is log T 5.2 to 9 (rows 5 + 0.2j, `fermi_log10t_min`/`_max`) and X = log rho − log mu_e − 1.5 log T from −16 to −3. `eos_init` refuses a Saha cutoff below 5.7, so the ramp never reaches below the table; a run's first point beyond it elsewhere (above log T 9, or past X = −3, very degenerate gas) gets one warning in the run log, and is extrapolated.
- **Errors** return ierr = 1 (SCV's spline failure returns numerics' status): a missing or corrupt table (every reader checks its reads, the places of the Fermi-Dirac cards, the SCV tables' common grid and the OPAL T6 grid, and closes a table it rejects), a bad configuration (`opal_vintage_conflict`, which `derive_controls` never produces since it lets the newest vintage win with a warning, an eta range with eta_max ≤ eta_min, a cutoff below 5.7), beta ≤ 0, the fully ionized iteration not converging, an MHD point outside its tables, an OPAL interpolation error. A point outside an OPAL or SCV table is not an error: the result so far covers it. `eos_get`'s caller decides what an error means.
- **Messages** go through `eos_messages`: a warning or diagnostic the run continues after to the run log (`eos_log`; `eos_log_once` for the Fermi-Dirac range), an error's message to the terminal and the run log (`eos_report`). Each is formatted where the condition is found, into lines (internal writes, so the formats are the original ones), and the table readers, which do the I/O, are apart from the kernels: `yale_eos`, `scv_eos`, `opal_eos` and `mhd_eos` hold no I/O statement (`check_no_io`). The Saha solver's non-convergence and OPAL's density search falling back are warnings: the run continues with the last iterate, or the result without OPAL. Some messages keep the F77's capitals; their rewording waits for a re-pin commit.
- **Controls** come only through `eos_settings` and `eos_table_paths_t`, built in setup/; eos never reads `s%control`/`s%physics`.
- **Kept by decision, with the reason** (each a measured change for a later round, docs/bug_sweep.md): the OPAL 1995 stencils are not clamped to its table (M94): `stencil_entry` gives a read past the table what the 1995 code's own storage held there, so the 1995 numbers stay those of its code until the fix, which sets `clamps_stencil` for 1995 and moves test_eos's pin at (log T 3.70, log P 4.5) by 2.4e-6 dex; SCV has no edge weight (M88); the second derivatives inside a table are the Yale solution's (O2); OPAL's and MHD's 1/mu_ion count a zone's extra metals as carbon (M129); SCV's H⁺ fraction (M128) and its general path's mixing entropy before the envelope table exists (M92); the SCV downward smoothing's reach into a shorter row (M206).

## 4. Configuration and state

**Configuration.** `eos_settings_t` holds the MHD switch, the OPAL vintage (`opal_vintage_none`, `_1995`, `_2001`, `_2006`), the SCV switch, numerical derivatives, diffusion of Z, the Debye-Hückel correction and its degeneracy range, and the Saha cutoff. Its defaults are the controls' (test_eos.py checks them against `defaults/physics.defaults`), so a star-free caller sets only what it changes.

`eos_table_paths_t` holds the tables' paths. The star layer fills both from the controls in `setup/eos_settings_builder.f90`. The OPAL vintage and its path are derived controls (`io/derived_controls.f90`'s `select_eos_and_opacity_tables`, from the three buffer switches `use_opal95_eos`, `use_opal2001_eos`, `use_opal2006_eos` and their paths; the newest selected wins, with a warning).

**State** that outlives a call, and who sets it:

| State | Module | Set by | Lifetime |
|---|---|---|---|
| `eos_settings` | `eos_configuration` | `eos_init` | the run |
| `eos_mix` | `eos_mixture` | `eos_set_mixture` | until the next mixture |
| the warnings' flags | `eos_messages` | `eos_log_once`; reset by `eos_init` | the run |
| the Fermi-Dirac table | `yale_eos_tables` | `eos_init` | the run |
| the SCV tables; the envelope table and its composition | `scv_eos_tables` | `eos_init` (which drops the envelope table); `eos_set_mixture` | the run; until the next mixture |
| SCV's search start `last_scv_cell` | `scv_eos_tables` | every SCV lookup; a run starts at the tables' first cell | the process: a lookup's cell depends on it only for a point exactly on a node |
| the OPAL tables, their loaded flags and Z, the open unit | `opal_eos_tables` | `eos_init` opens the unit and marks the vintage unread; its reader reads it on first use | until the next `eos_init` that names the vintage |
| OPAL's stencil, `quad` caches, last results, `t_row_index` | the vintage's instance | every interpolation; `find_ragged_edge_ramp` starts from the last stencil (its result does not depend on it), `find_opal_ramp` from its cached row | the process |
| the MHD tables and the last interpolation | `mhd_eos_tables` | `eos_init`; every MHD call | the run |

## 5. Tests and coverage

**test_eos** (`make tests` in `src/`, then `pytest test_eos.py` from the repository root) calls eos without a star, through `eos_lib` and `eos_test_entries`: it installs the settings itself and sets the envelope mixture (core's algorithm on the default mixture seed, X = 0.70, Z = 0.016492). Its standard output is byte-compared with `eos/test/expected_test_eos.out`; it writes its run log to `test_eos.short` in its working directory. The pinned values: the production EOS (SCV and OPAL 2006) on a (log T, log P) grid; `gamma1_at_density` on both branches; the error paths; every configuration a deck can choose (OPAL 1995 and 2001, SCV only, Yale only, numerical derivatives, diffusion of Z, Debye-Hückel); each configuration's derivatives and their finite differences; the MHD lookup on a synthetic table. With `YREC_MHD_TABLES` naming the real MHD tables it also evaluates the MHD path (the output then differs from the pin).

Its 65 ok/FAIL checks (a FAIL makes it exit 1):
- derivatives against centred finite differences where they must agree: the Saha solution's (also with the Debye-Hückel correction on), the fully ionized gas's, OPAL 2006's, 1995's and 2001's first, the numerical second derivatives;
- MHD: the table lookup's range and lir's status; the reader on eight synthetic tables (an ideal gas, so the interpolation is exact): the density in each family, the tabulated gradients, ierr outside the tables in T and in P, a missing and a truncated file;
- limits: a neutral ideal gas, a weakly degenerate fully ionized one (with the first Fermi-Dirac correction), grad_ad and cp of the gas plus radiation, grad_ad → 1/4; beta ≤ 0 is ierr; the fully ionized iteration's failure is ierr; a point beyond the Fermi-Dirac table is reported once;
- continuity of every blended field at both ends of the Saha ramp (Yale alone and with SCV), and at both ends of the OPAL 2006 edge ramp;
- `in_atmosphere`; the SCV envelope table for a new composition, and not for another one; SCV's general path before the envelope table; R/mu = beta P/(rho T) on both SCV paths and its continuity; SCV's cell independent of its search start;
- the Debye-Hückel correction: it raises rho, by the ratio of its charge-weighted factor between two mixtures, its heavier metals following the zone's Z, a zone's metals adding up to Z;
- `gamma1_at_density` below the OPAL table; diffusion of Z at the envelope's Z and with a metal-free envelope;
- `eos_init`'s setting checks (two vintages, the eta range, the cutoff) and the record's defaults; the table-source list;
- OPAL: the 1995 and 2001 trial calls; a truncated table on first use (after another table was read: M105), a missing one at `eos_init`; an entry only the second X quadratic reads; a huge entry that either X quadratic weights negatively, and one read alone at X = 0, each an error (M292); a 1995 table with a zero T6;
- the readers: SCV (a re-read, a missing, a truncated and an off-grid table), Fermi-Dirac (a missing, a truncated, a misplaced card, abscissae too wide).

Every check has a mutation in `eos/test/test_eos.mutations.toml` that makes it print FAIL, and `python3 src/tools/mutate.py eos` (after `make tests` in src/) reports all 70 caught (2026-11-11). `test_eos.py` also checks `eos_settings_t`'s defaults against the controls' (`test_eos_settings_defaults_are_the_controls`). Regenerate the pin only for a deliberate change: run `src/test_eos` with `YREC_INPUT` set in an empty directory, review the diff, and commit it saying why.

Under the development flags (`make devel_tests`, `python3 src/tools/devel_tests.py`) test_eos stops at M92 (SCV's general path takes log10 of a zero mixing entropy before the envelope table is built: a divide-by-zero trap), a measured fix for a later round; its OPAL 1995 reads past the arrays (M94) no longer stop it.

**The decks** pin the EOS in use: 304 of the 310 decks of testsuite/ and examples/ use OPAL 2006, SCV and the Debye-Hückel correction, two SCV and the correction without OPAL, four the Yale EOS with the correction alone; none uses MHD or OPAL 1995 or 2001, which test_eos alone covers.

**Coverage**: `python3 src/tools/coverage.py --work DIR --dir eos --test test_eos --deck examples/run_config_matrix/cm_base.inlist --deck testsuite/Test_solarGS98_base.inlist --deck testsuite/Test_solarGS98_yaleeos.inlist --deck examples/run_from_dbl_to_zams/Test_m0030_feh0_noGS_ZAMS.inlist` (a -O0 --coverage build in DIR; the four decks reach the OPAL 2006, SCV and Yale paths) gives 88.5% of 3045 lines (2026-11-10): `blended_eos`, `eos_lib`, `eos_mixture`, the descriptor, the list and the messages 100%, the Yale solver 99.7%, SCV 97.1%, MHD 92.9% (on the synthetic tables), the OPAL code 74.8% and its reader 67.2%. What is not run is mostly error branches and their messages: OPAL's for a corrupt table's other faults (a block out of order, a short density column, a mass fraction above 1, a Z that differs between the X tables, an empty region of the stencil), the interpolators' invalid-argument messages, MHD's reader's checks of a table's composition.

## 6. Extending

**Add a control**: the field in `eos_settings_t` (with the control's default), copied in `eos_settings_from_star` (`setup/eos_settings_builder.f90`; then `python3 tools/check_star_contracts.py --fix`), read as `eos_settings%<field>`; the control itself through the registry (io/README.md §5, "Add a control"); test_eos.py's defaults test then checks the default.

**Add a table source** (a table or a solver that covers part of the (T, rho) plane). Traced on 2026-11-10 on a throwaway branch: a toy over-the-gas source (steps 1 to 3 and 6) built, and test_eos kept its pin; `make check` then failed until step 4 sets the two new fields (check_boundaries' settings-field contract keeps that step from being forgotten). Steps 4 and 5 are the builder's and io's own recipes. The files:
1. `eos/<name>/<name>_eos.f90`: a module with
   - `<name>_table_source()`, returning its `eos_table_source_t`: a name, a stage (`stage_in_place_of_saha`, `stage_over_the_gas`, or `stage_whole_eos`), and the procedures below;
   - `selected` (a logical function of `eos_settings`), `read_tables(paths, ierr)`, `in_window(log10_temperature, r)` (or `always_in_window`);
   - `lookup(log10_temperature, log10_pressure, X, Z, table_r, covers, ierr)`. `table_r` arrives as a copy of the result so far: T, P (cgs), beta, beta14, 1/mu_ion of the zone, the Saha solver's fields, and the gas's density as a guess. It fills the density and its log, the four first derivatives, and the mean weights (`mean_weights_from_density(table_r)` gives R/mu and 1/mu_e from the density, T, P, beta and the 1/mu_ion it holds); it sets `covers` false where its table does not reach;
   - `edge_weight(temperature, log10_density, weight, in_table, near_edge, ierr)` (or `whole_weight_inside`); a weight across a range of log T should be a smoothstep, as `blend_saha_into`'s, so that `continuity_in_t` can check it;
   - `adopt` and `blend`: `adopt_table_result` and `blend_table_result` carry the density, beta, the mean weights and the first derivatives, and leave the ionization fractions, eta and the second derivatives the result's; `adopt_whole_result` takes the whole record (SCV's), the fields the lookup did not set included.

   Its tables in `eos/<name>_eos_tables.f90`, its reader in `eos/<name>/<name>_table_reader.f90` (open on a `newunit` with iostat, every read with iostat, a corrupt file ierr = 1 with a message through `eos_report`, the unit closed). A module that calls `exp`, `log10`, ... has `use math_lib`. The Makefile finds the files; run `python3 tools/gen_makefile_deps.py` in `src/`. The source module goes into `STAR_BLIND_FILES` in `tools/check_boundaries.py` (not its reader).
2. `eos/eos_source_list.f90`: its entry in `table_sources()`, in the order of application within its stage (a later over-the-gas source wins where two cover a point), and `n_table_sources`.
3. `eos/eos_configuration.f90`: its switch in `eos_settings_t` and its path in `eos_table_paths_t`.
4. `setup/eos_settings_builder.f90`: the control and the path, copied in `eos_settings_from_star` and `eos_table_paths_from_star`; then `python3 tools/check_star_contracts.py --fix`.
5. The control and the table path: io/README.md §5, "Add a control" and "Add a table path" (`io/controls_buffer.f90`, `defaults/controls_registry.tsv` with the `settings` cell `eos_settings%<switch>`, the generators, the controls dumps: `YREC_REGEN_CONTROLS=1 pytest test_controls.py`).
6. `eos/test/test_eos.f90`: an option of `configure_eos`, a configuration in `check_configurations` and `check_derivatives`, the source in `source_list_checks`, and ok/FAIL checks: a limit it must reproduce, its reader on a missing and a truncated file, continuity at both ends of its edge ramp (`ramp_end_is_continuous`); a mutation for each in `eos/test/test_eos.mutations.toml`; regenerate the pin.
7. This README: §1, §2.1 and its field table, §4.

Nothing else changes: `blend_eos`, `eos_get` and `eos_init` loop over the list.

**Add an OPAL vintage** (a table of the 2001/2006 kind). Traced on 2026-11-10 on a throwaway branch: a fourth vintage whose row copies 2006's (steps 1 and 2), read from the 2006 table, built and gave 2006's density to the last bit, its messages in 2006's words; steps 4 and 5 were checked against the 2006 vintage's own rows and branch, not traced:
1. `eos/eos_configuration.f90`: `opal_vintage_<year>`, the next number.
2. `eos/opal_eos_tables.f90`: the table's size (`n_eos<yy>_nr`, `_nt`; `n_eos_nr_max`/`_nt_max` if it is larger), its grids (X, the lowest T6 row per density column, the highest density column per T6 row), `n_opal_vintages` and `iv_opal<yy>`, its row of `opal_layout` (every field), and the upper bound of `opal_tables`.
3. `eos/opal/opal_eos.f90` and `eos/opal/opal_table_reader.f90`: nothing, unless its messages differ from 2006's words (the `report_*` procedures branch on 1995 and 2001 and take 2006's otherwise).
4. `io/derived_controls.f90`: the vintage's branch in `select_eos_and_opacity_tables` (the newest selected wins; it sets `opal_eos_vintage` and `opal_eos_table_path`, which the builder already copies), and its switch in the check for more than one vintage.
5. Its switch and its path are buffer controls that only that read uses: their rows in `defaults/controls_registry.tsv` have home `buffer`, used `read` and the `settings` cell blank, as `use_opal2006_eos` and `opal2006_eos_table_path` do (io/README.md §5, "Add a control", step 2), declared in `io/controls_buffer.f90`, then the generators. `setup/eos_settings_builder.f90` needs no edit.
6. `eos/test/test_eos.f90`: a configuration in `check_configurations` and `check_derivatives` with its table, the edge-ramp continuity for it, the reader on a truncated copy; regenerate the pin.
7. This README: §1, §2.1, §4.

**Add a blended region** (a ramp between two solutions, like the Saha ramp): in `eos/blended_eos.f90` keep a copy of the record with the first result (as `gas_eos`'s `saha = r`) before computing the second into it, and blend with a weight from 1 at one end to 0 at the other, as `blend_saha_into` does every field (with `want_derivatives` the second derivatives too); its position and width are named parameters with their units, or settings (`eos_configuration`, the builder, a control as above); `continuity_in_t` in test_eos checks both ends (`results_agree` names the fields that jump); this README §2.1.

## Sources

The Yale EOS: Prather (1976), with the Debye-Hückel correction of Guenther, Demarque, Kim & Pinsonneault (1992). SCV: Saumon, Chabrier & Van Horn (1995); the metal table was computed with the Prather EOS. OPAL: Rogers, Swenson & Iglesias (1996), Rogers & Nayfonov (2002), the 2005 update for the 2006 tables. MHD: Hummer & Mihalas (1988), Mihalas, Däppen & Hummer (1988), Däppen et al. (1988). Each module's Sources section says what its code follows and what it was in the F77.

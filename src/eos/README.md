# eos/ — the equation of state

This directory computes the density of a zone and its thermodynamic quantities (beta, the ionization fractions, the mean molecular weights, cp, the adiabatic gradient) and their derivatives, from (log T, log P, X, Z). Code outside eos/ uses it only through `eos_lib` (enforced by `tools/check_boundaries.py`, which lists eos's three entries). eos/ reads nothing from `star` (it is in `STAR_FREE_DIRS`, and its solvers in `STAR_BLIND_FILES`): its configuration is one settings record and one record of table paths, both built by the star layer and handed to `eos_init`, and the envelope mixture the star layer hands `eos_set_mixture`. Its constants are compile-time parameters of `state/physical_constants.f90`, so it gives the right answer without any other setup having run.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `eos_lib` (`eos_lib.f90`) | **The public interface.** `eos_init` installs the settings and reads the tables they select; `eos_get` evaluates a zone into an `eos_result`; `eos_set_mixture` takes the envelope mixture the star layer owns (and builds the SCV envelope table for it). It re-exports the settings and paths types and the `opal_vintage_*` constants. |
| `eos_configuration` | `eos_settings_t` (everything eos reads from the run's controls, with the controls' defaults) and the installed `eos_settings`, set by `eos_init` (or, in test_eos, `set_eos_settings`); `eos_table_paths_t` (the tables, by path); the `opal_vintage_*` constants. |
| `eos_results` | `eos_result`, the record every EOS family fills: one named field per quantity, each with its definition and unit, and the fields that are also inputs (§3). |
| `eos_mixture` | The envelope mixture (`eos_mix`) and the atomic weights of its 12 species; the composition terms' constants (`inverse_a_hydrogen`, ...); `ion_mean_weight_excess`, `compute_mean_weights` and `mean_weights_from_density` (the mean weights a density implies); `debye_huckel_composition_of`, the Debye-Hückel correction's composition at a zone (X, Y, Z, C, N and O the zone's; Ne, Na, Mg, Al, Si, Ar and Fe in the envelope mixture's proportions, their sum the zone's Z less its C, N and O, so that the metals add up to Z and the C and O a core makes stay C and O). |
| `eos_messages` | `eos_log` and `eos_report`, eos's two kinds of message (§3). |
| `blended_eos` | The EOS of every run without MHD (§2): `blend_eos` evaluates one point in named steps (`zone_mixture`, `radiation_share`, `gas_eos`, `opal_step`), with the blends `blend_saha_into`, `adopt_table_result` and `blend_table_result`; `compute_blended_eos` adds the numerical second derivatives; `gamma1_at_density` (test_eos only) gives Γ₁ and ∇_ad at (T, ρ). |
| `yale_eos` (`yale/`) | The Yale EOS of Prather (1976): `compute_saha_eos` (partial ionization: the Saha equations of 11 species, H and 10 metals, and He's two stages, by Newton iteration; the steps `charge_balance`, `eliminate_helium`, the derivative steps), `compute_fully_ionized_eos` (degenerate, relativistic electrons from the Fermi-Dirac table, optional Debye-Hückel correction) and `read_fermi_table`. |
| `scv_eos` (`scv/`) | The Saumon-Chabrier-Van Horn EOS: `read_scv_tables`, `build_scv_envelope_table` (the H, He and metal tables combined at the envelope composition), `compute_scv_envelope_eos` (that table, when the zone has the X and Z it was built for, to 1e-5) and `compute_scv_eos` (any X, Z, from the three tables), which share `find_scv_cell` and `scv_ionization_fractions`. |
| `opal_eos` (`opal/`) | The OPAL EOS: `compute_opal_eos` (density and derivatives at P, for any vintage) and `find_opal_edge_ramp` (whether the vintage's table covers the point, and its weight near the edge); the 2001 and 2006 tables' code, written once and driven by the vintage's row of `opal_layout`: `opal_density_search` (OPAL's `rhoofp01`/`rhoofp06`), `interpolate_opal_eos` (entered as OPAL's `esac01`/`esac06`), `add_radiation`, `readcoeos`, the edge ramps `find_table_edge_ramp` (2001) and `find_ragged_edge_ramp` (2006); and their messages in each vintage's words (`report_*`). |
| `opal95_eos` (`opal/`) | The OPAL 1995 tables, with their own state layout (see `opal_eos_tables`): `find_opal95_ramp`, `rhoofp`, `esac`, `t6rinterp`, `radsub`, `gmass`, `readco`; they read the 1995 row of `opal_layout`. |
| `mhd_eos` (`mhd/`) | The Mihalas-Hummer-Däppen EOS: `read_mhd_tables` and `compute_mhd_eos`. Däppen's routines keep their names (`mhdst`, `mhdst1`, `rabu`, `rtab`, `mhdpx`, `mhdpx1`, `mhdpx2`); `interpolate_mhd_table` is the 4x4 table lookup. |
| `opal_eos_tables`, `mhd_eos_tables`, `scv_eos_tables`, `yale_eos_tables` | Each family's tables and state. `opal_eos_tables` holds the instances (`opal95`; `opal_tables(iv_opal01:iv_opal06)`, one array, since 2001 and 2006 share one state type), `opal_shared` (the open table's unit, the loaded flags, each table's Z), and `opal_layout`, one row per vintage with everything the code tells the vintages apart by: the table's size and grids, the slots of `eos_output`, the radiation constants, the reader's row format, the density search's edge table, a/c and radiation handling, the kind of edge, `compute_opal_eos`'s T6 range and `opal_step`'s window; and the helpers the vintages share (`interpolate_opal_t6_rho`, `find_opal_ramp`, `quad`, `gmass01`, `open_opal_eos_table`). `yale_eos_tables` holds the Fermi-Dirac table and the Debye-Hückel constants and composition type, `scv_eos_tables` the SCV tables and the envelope table with the composition it was built for. |

## 2. What eos_get does

With `use_mhd_eos`, `eos_get` calls `compute_mhd_eos` (§2.2). Otherwise it builds the Debye-Hückel composition of the zone (when the correction is on) and calls `compute_blended_eos`, which calls `blend_eos` once, or, with `want_derivatives` and `use_numerical_derivatives`, four more times at log T ± 3e-4 and log P ± 1.5e-3 for centred-difference second derivatives (each shifted evaluation starts from the caller's record and the previous one's density, and leaves the record as it was).

### 2.1 blend_eos

For one (log T, log P, X, Z), with the cutoff T_c = `saha_log10t_cutoff` and the ramp width 0.5, four steps, each taking the result so far:

1. **The zone's mixture** (`zone_mixture`): 1/mu_ion, and the Saha equations' abundances, from the envelope mixture `eos_mix`. With `use_diffusion_z` all the metals scale with Z/Z_env; otherwise (and when the envelope has no metals) the zone differs from the envelope by its excess H, He and metals (counted as carbon). Both use 1/A = 0.9921 for H and 0.24975 for He.
2. **Radiation** (`radiation_share`): beta = P_gas/P = 1 − (a/3) T⁴/P. beta ≤ 0 is an error (ierr = 1).
3. **The gas** (`gas_eos`), by log T:

   | log T | source |
   |---|---|
   | ≥ T_c | fully ionized (`fully_ionized_gas`: 1/mu_e from X, then `compute_fully_ionized_eos`) |
   | T_c − 0.5 < log T < T_c | the Saha (or SCV) result blended into the fully ionized one (`blend_saha_into`), weight f = w³(6w² − 15w + 10), w = 2(T_c − log T) |
   | ≤ T_c − 0.5 | the Saha equations (`compute_saha_eos`) |

   With `use_scv_eos`, below T_c SCV's result (density, the mean weights it implies, the ionization fractions and the first derivatives; Saha's second derivatives stay) replaces Saha's where SCV's table covers the point; elsewhere Saha's result stands.
4. **OPAL** (`opal_step`), when a vintage is selected (`opal_vintage`) and the point is inside its (T, rho) window (its row of `opal_layout`): `compute_opal_eos` finds the density that gives P and the derivatives there; where the table does not cover the point the result so far stands; inside the table OPAL's density, beta, mean weights and first derivatives replace it (`adopt_table_result`), and near the table's edge (`find_opal_edge_ramp`: one table cell wide) they are blended in linearly (`blend_table_result`). The second derivatives stay the Yale/SCV ones.

So a deck with SCV and OPAL 2006 (298 of the 304 shipped decks, the reference solar decks among them) gets OPAL wherever its table covers the point (log T from 3.27 up), SCV at low T outside it, the Saha solution in between, and the fully ionized gas at high density (log rho > 7) and in the stellar centre above the OPAL table.

### 2.2 compute_mhd_eos

The MHD tables (three ZAMS-type compositions, each in a lower and an upper part in T, and five centre-type ones) are interpolated bicubically in (log P, log T) in each table of the family the temperature selects, then in X (quadratically across the ZAMS compositions, cubically across four of the centre ones). A point outside the tables returns ierr = 1, and a lookup the interpolation cannot do lir's status (numerics_status). The tables are not shipped.

## 3. Conventions

- **Units** are cgs: T in K, P in dyn/cm², rho in g/cm³, cp in erg/(g K), and log is log₁₀. The OPAL routines work internally in T6 = T/10⁶ and P in units of 10¹² dyn/cm².
- **Derivatives are logarithmic.** The first derivatives are ∂lnρ/∂lnT at constant P and ∂lnρ/∂lnP at constant T. The six second derivatives are (test_eos's derivative check compares each with finite differences, and checks them where they must agree, §6):

  | Output | Definition |
  |---|---|
  | `dlnrho_dlnt_dt` | ∂(∂lnρ/∂lnT)/∂lnT |
  | `dlnrho_dlnt_dp` | ∂(∂lnρ/∂lnT)/∂lnP |
  | `adiabatic_gradient_dt`, `_dp` | ∂ln∇ad/∂lnT, ∂ln∇ad/∂lnP |
  | `specific_heat_cp_dt`, `_dp` | ∂ln cp/∂lnT, ∂ln cp/∂lnP |

- **`want_derivatives`** asks for the second derivatives; on the blended path, without it they keep the caller's values (the MHD path sets them always). Inside an OPAL or SCV table they are the Yale (Saha) solution's, not the table's, so they can be far from those of the EOS actually used; `use_numerical_derivatives` makes them consistent.
- **`in_atmosphere`** skips the Saha solver's second derivatives (they keep the caller's values); everything else, the first derivatives included, is computed as without it.
- **beta** is P_gas/P, and `beta14` is 4(1 − beta)/beta = 4 P_rad/P_gas, the radiation term of the Yale solvers, on every path.
- **The mean weights**: where the density comes from the composition (the Saha solution, the fully ionized gas) R/mu = R(1/mu_ion + 1/mu_e); where it comes from a table or a blend (SCV, OPAL, MHD, the Saha ramp) R/mu = beta P/(rho T), the gas constant the ideal-gas law needs for that density, and 1/mu_e is backed out of it. The two differ by the degeneracy correction, so R/mu steps by 7e-4 at the Saha cutoff (an open question to the authors, docs/bug_sweep.md).
- **Inputs carried in the result**: the density is the starting guess of the fully ionized density iteration, and the `saha_*` fields (E and the He+ and He++ fractions) the Saha solver's start, which it returns its solution in; so a result depends, within those iterations' tolerances, on the record it was given (a record carried from zone to zone warm-starts both; a new one starts from E = 0.1). The second derivatives are inputs as above. Every other field is overwritten; log T and log P come back unchanged. `saha_iterations` is a counter the Saha solver adds its iterations to; nothing in eos reads it.
- **Errors** return ierr (positive): a missing or corrupt table (every reader checks its reads, and closes a table it rejects), a bad configuration (two OPAL vintages, a Debye-Hückel degeneracy range with eta_max ≤ eta_min), beta ≤ 0, an MHD point outside its tables, an OPAL interpolation error, an OPAL vintage with no table open to read, a failed spline in the SCV envelope table. A point outside an OPAL or SCV table is not an error: the result so far covers it.
- **Messages** go through `eos_messages`: a warning or a diagnostic the run continues after (a Saha or density iteration that did not converge, an OPAL lookup outside its table, a density search that falls back) to the run log (`eos_log`); every message of an error eos returns ierr for to the terminal and the run log (`eos_report`). Each is formatted by its caller into lines (an internal write, so the formats are the original ones) and printed without trailing blanks; no other eos module writes to a unit.

## 4. Configuration

`eos_settings_t` (in `eos_configuration`) is everything eos reads: MHD; the OPAL vintage (`opal_vintage`: `opal_vintage_none`, `_1995`, `_2001`, `_2006`); SCV; numerical derivatives; diffusion of Z; the Debye-Hückel correction and its degeneracy range; the Saha cutoff. Its defaults are the controls' (the cutoff log T = 6, the range −1 to 1, every switch off), so a star-free caller sets only what it changes. `eos_table_paths_t` holds the tables' paths (Fermi-Dirac, SCV, MHD, the OPAL table of the selected vintage). The star layer fills both from the controls (`eos_settings_from_star` and `eos_table_paths_from_star` in `setup/microphysics_settings.f90`; the controls let the newest vintage win, and settings that still select two are `opal_vintage_conflict`, which `eos_init` refuses) and passes them to `eos_init` once per run (`install_run_settings` in `setup/settings_install.f90`, called from `core/run/run_setup.f90`). eos opens every table itself with `newunit`; the OPAL table is read by its vintage's reader on first use, the others at `eos_init`. The envelope mixture comes from `eos_set_mixture`: the run plan's envelope Z in `prepare_run` (before the opacity tables are read for it), the starting model's mixture, and the rescaling of Z.

## 5. State

What eos keeps between calls, module by module, and who sets it:

| State | Module | Set by | Lifetime |
|---|---|---|---|
| `eos_settings` | `eos_configuration` | `eos_init` (`set_eos_settings` in test_eos) | the run |
| the envelope mixture `eos_mix` | `eos_mixture` | `eos_set_mixture` (`prepare_run`, the starting model, the rescaling of Z) | until the next call |
| the Fermi-Dirac table | `yale_eos_tables` | `eos_init` | the run |
| the SCV tables, and the envelope table | `scv_eos_tables` | `eos_init` reads the tables (and drops the envelope table); `eos_set_mixture` builds the envelope table for the envelope's X and Z. The table records the X and Z it was built for, and the lookup tests the zone against those, so a zone at another composition takes SCV's general path | the run; until the next mixture |
| SCV's warm-start cell (`idtt`, `idp`) | `scv_eos_tables` | every SCV lookup (`find_scv_cell`): the next search starts there; a lookup's cell depends on it only for a point exactly on a table node | never reset |
| the OPAL tables, their loaded flags and Z, the open table's unit | `opal_eos_tables` (`opal95`, `opal_tables`, `opal_shared`) | the vintage's reader on first use (a table counts as loaded once read); `eos_init` opens the unit | the process (a reset keeps the loaded tables) |
| OPAL's interpolation stencil, quadratic coefficients (`quad` caches) and last results (`eos_output`) | the vintage's instance | every interpolator call; `find_ragged_edge_ramp` reuses the last stencil, `find_opal_ramp` its cached row (`t_row_index`) | never reset |
| the MHD tables and the last interpolation (`mhd_output`) | `mhd_eos_tables` | `eos_init`; every MHD call | the run |

Nothing else outlives a call: the Saha solver's start and the Debye-Hückel composition travel with the call (the result record, and a record `eos_get` builds).

## 6. Tests and coverage

`test_eos` (`make tests` in `src/`, then `pytest test_eos.py` from the repository root) runs without a star: it calls `eos_init` itself and sets the envelope mixture (core's algorithm on the default mixture seed, X = 0.70, Z = 0.016492); set again after `eos_init`, as a run's starting model does, the mixture builds the SCV envelope table. Its output is byte-compared with `eos/test/expected_test_eos.out`:

- pinned values: the production EOS (SCV and OPAL 2006) on a (log T, log P) grid; `gamma1_at_density` on both branches; the error paths; every configuration a deck can choose (OPAL 1995 and 2001, SCV only, Yale only, numerical derivatives, diffusion of Z, Debye-Hückel); the derivatives against finite differences; the MHD table lookup on a synthetic table. With `YREC_MHD_TABLES` set to a directory with the real MHD tables it also evaluates the MHD path.
- ok/FAIL checks (49), each confirmed to fail under a mutation of the code it covers (a failed check makes the program exit 1):
  - the derivatives against centred finite differences where they must agree: the Saha solution's (to 1e-6 and 1e-4), the fully ionized gas's (1e-3 and 0.1), OPAL 2006's first (3e-3), the numerical second derivatives (1e-4);
  - the MHD reader on eight synthetic unformatted tables (an ideal gas in every row, so the interpolation is exact): the density in each family and at the ZAMS/centre boundary, the tabulated gradients, ierr outside the tables in T and in P, a missing and a truncated file; the table lookup's range and lir's status;
  - the ideal-gas limits: a neutral gas (rho = mu P_gas/(R T), grad_ad → 0.4, cp → 5R/(2 mu)), a weakly degenerate fully ionized one (rho = mu P_gas/(R T)·(1 − 2^(−5/2) e^η P_e/P_gas), to 2e-4), grad_ad and cp against the ideal gas plus radiation, grad_ad → 1/4 when radiation dominates; beta ≤ 0 returns ierr;
  - continuity of every blended field (log rho, rho, beta, R/mu, 1/mu_e, the ionization fractions, the four first and six second derivatives; at the cutoff without the mean weights, §3) at the foot of the Saha ramp (Yale alone and with SCV under it) and at the cutoff, and of every field the table blend sets at both ends of the OPAL 2006 edge ramp (each end bisected to where the blend first departs from its reference); `ramp_end_is_continuous` takes the log T and the source's references, so any table source can be checked so;
  - `in_atmosphere`: the density, ionization and first derivatives of a full evaluation;
  - SCV: the general path before the envelope table is built; the table `eos_set_mixture` builds for a new composition, and a zone at an earlier table's composition on the general path; R/mu = beta P/(rho T) on both paths and continuous at the ramp's foot; the reader (a re-read drops the envelope table; a missing or truncated table is ierr, no table left open);
  - the Debye-Hückel correction raises rho, the more with heavier metals in the envelope mixture, and its heavier metals scale with the zone's Z when all its metals do; a zone's record has its metals add up to its Z, and carbon added to a zone's Z leaves its Ne..Fe alone;
  - `gamma1_at_density` below the OPAL table takes the Yale/SCV value;
  - diffusion of Z: at the envelope's Z the default mixture's 1/mu_ion and rho, and an envelope without metals takes the default mixture;
  - `eos_init` refuses two OPAL vintages and a Debye-Hückel range with eta_max ≤ eta_min; the settings' defaults are the controls';
  - the OPAL 1995 and 2001 density searches' trial calls leave the other slots and give the full call's pressure;
  - the Fermi-Dirac and OPAL 2006 readers: a missing and a truncated table are ierr = 1, the truncated one closed.

A deliberate change to the output means regenerating `eos/test/expected_test_eos.out`. The decks pin the full EOS in use: 298 of the 304 shipped decks (testsuite and examples) use OPAL 2006, SCV and the Debye-Hückel correction, two SCV and the correction without OPAL (`Test_m0p3feh+0p0GS98_scv_ZAMS`, `Test_m0p08feh+0p0GS98_allardscv_15Gyr`), and four the Yale EOS with the correction alone (`Test_solarGS98_yaleeos` and three more `*yaleeos*` decks); none uses MHD, OPAL 1995 or 2001, which test_eos alone covers.

`python3 src/tools/coverage.py --work DIR --dir eos --test test_eos --deck ...` measures the line coverage of eos/ (a -O0 --coverage build in DIR). With test_eos and four decks that reach the OPAL, SCV and Yale paths (`examples/run_config_matrix/cm_base.inlist`, `testsuite/Test_solarGS98_base.inlist` and `Test_solarGS98_yaleeos.inlist`, `examples/run_from_dbl_to_zams/Test_m0030_feh0_noGS_ZAMS.inlist`) it is 85.5% of 3262 lines (2026-10-29): `eos_lib`, `eos_mixture` and the small modules 100%, `blended_eos` and the Yale solvers 97.9%, SCV 97.2%, MHD 87.5% (the synthetic tables), `opal_eos_tables` 86.1%, the OPAL 1995 and 2001/2006 code 71%; what is not run is mostly the readers' and interpolators' error branches and their messages.

## 7. Extending

**Add a control:** add the field to `eos_settings_t` (with the control's default), copy it from `star%ctrl` in `eos_settings_from_star` (`setup/microphysics_settings.f90`; then `tools/check_star_contracts.py --fix`), and read `eos_settings%<field>`. The control itself comes through the registry (`defaults/controls_registry.tsv` and the generators; io/README.md's controls recipe).

**Add an EOS source** (a table or a solver that covers part of the (T, rho) plane and is blended in like OPAL; this recipe was traced on 2026-10-29 by adding a throwaway source, steps 1 to 4 and 7, which built and passed its checks). A source is four things: a reader, a window (where it is consulted at all), a lookup that fills an `eos_result` and says whether it covers the point, and an edge weight (1 inside, falling continuously to 0 at its table's edge). The files:
1. `eos/<name>/<name>_eos.f90`, a module with
   - `read_<name>_tables(path, ierr)`: open on a `newunit` with iostat, every read with iostat; a missing or corrupt file is ierr = 1 with a message through `eos_messages` (`eos_report`), the unit closed;
   - `compute_<name>_eos(log10_temperature, log10_pressure, X, Z, r, covers, ierr)`: fills the density and log density, beta, the mean weights (`mean_weights_from_density` gives R/mu and 1/mu_e from the density), the ionization fractions and the four first derivatives, and leaves `r` as it was with `covers` false where its table does not reach;
   - its edge weight, and its window (a parameter, or a field of a layout row as `opal_layout`'s); a weight across a range of log T should be a smoothstep (as `blend_saha_into`'s), whose slope vanishes at the ends, so that `continuity_in_t` can check it; a weight across the table's edge in rho may be linear, as OPAL's (`ramp_end_is_continuous` bisects for its ends);

   its tables in `eos/<name>_eos_tables.f90`. A module that calls `exp`, `log10`, ... has `use math_lib` (check_boundaries' math-lib contract). The Makefile finds new files; run `python3 tools/gen_makefile_deps.py` in `src/`. A source that takes arguments in and returns results out, as the others do, goes into `STAR_BLIND_FILES` in `tools/check_boundaries.py`.
2. `eos/eos_configuration.f90`: its switch in `eos_settings_t` and its path in `eos_table_paths_t`.
3. `eos/eos_lib.f90`: `eos_init` reads its tables when the switch is on.
4. `eos/blended_eos.f90`: a step subroutine called from `blend_eos`, before `opal_step` if OPAL should win where both cover the point, after it otherwise; inside the window it calls `compute_<name>_eos` on a fresh record, then `adopt_table_result` where it covers the point, or `blend_table_result` with its edge weight near the edge. A source that replaces everything instead, as MHD does, is a branch in `eos_get`.
5. `setup/microphysics_settings.f90`: copy the control in `eos_settings_from_star` and the path in `eos_table_paths_from_star`; then `python3 tools/check_star_contracts.py --fix`.
6. The control and the table path: the recipes in `io/README.md` §5, "Add a control" and "Add a table path" (`io/controls_buffer.f90`, `defaults/controls_registry.tsv`, the generators, the controls dumps: `YREC_REGEN_CONTROLS=1 pytest test_controls.py`); a table path is a string control like `fermi_table_path`, expanded by `io/run_files.f90`'s `expand_paths`.
7. `eos/test/test_eos.f90`: an option of `configure_eos` (its switch and table path), a configuration in `check_configurations` and `check_derivatives`, and ok/FAIL checks that a mistake would break: a limit it must reproduce, its reader on a missing and a truncated file, continuity at both ends of its edge ramp (`ramp_end_is_continuous` with its two references: the source alone, and the blend without it); regenerate `eos/test/expected_test_eos.out`.
8. This README: §1, §2.1, §5.

**Add an OPAL vintage** (a table of the 2001/2006 kind; the 1995 tables have their own code):
1. `eos/eos_configuration.f90`: `opal_vintage_<year>` (the next number), and `eos_init`'s range check takes `n_opal_vintages`.
2. `eos/opal_eos_tables.f90`: the table's size (`n_eos<yy>_nr`, `_nt`; `n_eos_nr_max`/`n_eos_nt_max` if it is larger), its grids (X, the lowest T6 row per density column, the highest density column per T6 row, the density search's edge table), `n_opal_vintages` and `iv_opal<yy>`, its row of `opal_layout` (every field: sizes, grids, slots, a/c and the gas constant, the reader's row format and skipped columns, the routine id, the T6 range, the search's edge, a/c, radiation handling and name, the kind of edge, the window), and the upper bound of `opal_tables(iv_opal01:...)`.
3. `eos/opal/opal_eos.f90`: nothing, unless its messages differ from 2006's words (`report_*` branch on 2001 and take 2006's otherwise).
4. `io/derived_controls.f90`: the vintage's branch in `derive_controls`' choice (the newest selected wins; it sets `opal_eos_vintage` and `opal_eos_table_path`, which `eos_settings_from_star` and `eos_table_paths_from_star` in `setup/microphysics_settings.f90` copy), and its switch in the check for more than one vintage.
5. The control and the path, as for a source (steps 5 and 6 above).
6. `eos/test/test_eos.f90`: a configuration in `check_configurations` and `check_derivatives` with its table, the edge-ramp continuity for it, and the reader on a truncated copy; regenerate the pin.
7. This README: §1, §2.1, §4.

**Add a blended region** (a new ramp between two sources, like the Saha ramp):
1. `eos/blended_eos.f90`: in `gas_eos` (or a new step), keep a copy of the record with the first result (as `saha = r` does) before computing the second into it, and blend with a weight that runs from 1 at one end of the region to 0 at the other: `blend_saha_into` shows the smoothstep (continuous with its first two derivatives at both ends) and every field to blend (the log density, then the density and the mean weights from it; the ionization fractions; the four first derivatives; with `want_derivatives`, the six second derivatives).
2. The region's position and width: named parameters with their units, or fields of `eos_settings_t` (`eos/eos_configuration.f90`) copied from controls in `setup/microphysics_settings.f90` (the control as in step 6 above).
3. `eos/test/test_eos.f90`: `continuity_in_t` at both ends of a region in log T (`ramp_end_is_continuous` bisects for the ends of one in log P): it compares every blended field across each end (`results_agree`) and names the ones that jump.
4. This README: §2.1.

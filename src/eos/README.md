# eos/ — the equation of state

This directory computes the density of a zone and its thermodynamic quantities (beta, the ionization fractions, the mean molecular weights, cp, the adiabatic gradient) and their derivatives, from (log T, log P, X, Z). Code outside eos/ uses it only through `eos_lib` (enforced by `tools/check_boundaries.py`). eos/ reads nothing from `star` (it is in `STAR_FREE_DIRS`, and its solvers in `STAR_BLIND_FILES`): its configuration is one settings record and one record of table paths, both built by the star layer and handed to `eos_init`. Its constants are compile-time parameters of `state/physical_constants.f90`, so it gives the right answer without any other setup having run.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `eos_lib` (`eos_lib.f90`) | **The public interface.** `eos_init` installs the settings and reads the tables they select; `eos_get` evaluates a zone into an `eos_result`; `eos_set_mixture` and `eos_update_envelope_tables` take the envelope composition the star layer owns; `set_eos_settings` installs settings without reading tables; `eos_get_gamma1` returns gamma1 and the adiabatic gradient at (T, rho) (for test_eos). It re-exports the settings and paths types and the `opal_vintage_*` constants. |
| `eos_configuration` | `eos_settings_t` (everything eos reads from the run's controls) and the installed `eos_settings`; `eos_table_paths_t` (the tables, by path); the `opal_vintage_*` constants. |
| `eos_results` | `eos_result`, the record every EOS family fills: one named field per quantity, each with its definition and unit. |
| `eos_mixture` | The envelope mixture (`eos_mix`), the mixture's atomic weights, `ion_mean_weight_excess` and `compute_mean_weights`. |
| `eos_messages` | `eos_log` and `eos_report`, eos's two kinds of message (§3). |
| `blended_eos` | The EOS of every run without MHD (§2): `blend_eos` evaluates one point in named steps (`zone_mixture`, `fully_ionized_gas`, `blend_saha_into`, `adopt_opal`, `blend_opal_into`); `compute_blended_eos` adds the numerical second derivatives. |
| `yale_eos` (`yale/`) | The Yale EOS of Prather (1976): `compute_saha_eos` (partial ionization: the Saha equations for H, He and 11 metals by Newton iteration), `compute_fully_ionized_eos` (degenerate, relativistic electrons from the Fermi-Dirac table, optional Debye-Hückel correction) and `read_fermi_table`. |
| `scv_eos` (`scv/`) | The Saumon-Chabrier-Van Horn EOS: `read_scv_tables`, `build_scv_envelope_table` (the H, He and metal tables combined at the envelope composition), `compute_scv_envelope_eos` (that table, when the zone has the envelope's X and Z to 1e-5) and `compute_scv_eos` (any X, Z, from the three tables). |
| `opal_eos` (`opal/`) | The OPAL EOS: `compute_opal_eos` (density and derivatives at P, for any vintage) and `find_opal_edge_ramp` (whether the vintage's table covers the point, and its weight near the edge); the 2001 and 2006 tables' code, written once (`interpolate_opal_eos`, entered as OPAL's `esac01`/`esac06`; `rhoofp01`/`rhoofp06`, the density searches; `add_radiation`; `readcoeos`), and their messages in each vintage's words (`report_*`). |
| `opal95_eos` (`opal/`) | The OPAL 1995 tables, with their own state layout: `find_opal95_ramp`, `rhoofp`, `esac`, `t6rinterp`, `radsub`, `gmass`, `readco`. |
| `mhd_eos` (`mhd/`) | The Mihalas-Hummer-Däppen EOS: `read_mhd_tables` and `compute_mhd_eos`. Däppen's routines keep their names (`mhdst`, `mhdst1`, `rabu`, `rtab`, `mhdpx`, `mhdpx1`, `mhdpx2`); `interpolate_mhd_table` is the 4x4 table lookup. |
| `opal_eos_tables`, `mhd_eos_tables`, `scv_eos_tables`, `yale_eos_tables` | Each family's tables and state. `opal_eos_tables` also holds `opal_layout`, one row per OPAL vintage with everything the code tells the vintages apart by (table size and grids, the slots of `eos_output`, the radiation constants, the reader's row format, `compute_opal_eos`'s T6 range, `blend_eos`'s window), and the helpers the vintages share (`interpolate_opal_t6_rho`, `find_opal_ramp`, `quad`, `gmass01`, `open_opal_eos_table`). |

## 2. What eos_get does

With `use_mhd_eos`, `eos_get` calls `compute_mhd_eos` (§2.2). Otherwise it calls `compute_blended_eos`, which calls `blend_eos` once, or, with `want_derivatives` and `use_numerical_derivatives`, four more times at log T ± 3e-4 and log P ± 1.5e-3 for centred-difference second derivatives (each shifted evaluation starts from the previous one's density, and leaves `r` as it was).

### 2.1 blend_eos

For one (log T, log P, X, Z), with the cutoff T_c = `saha_log10t_cutoff` and the ramp width 0.5:

1. **The zone's mixture** (`zone_mixture`): 1/mu_ion, and the Saha equations' abundances, from the envelope mixture `eos_mix`. With `use_diffusion_z` all the metals scale with Z/Z_env; otherwise the zone differs from the envelope by its excess H, He and metals (counted as carbon). Both use 1/A = 0.9921 for H and 0.24975 for He.
2. **Radiation**: beta = P_gas/P = 1 − (a/3) T⁴/P. beta ≤ 0 is an error (ierr = 1).
3. **The gas**, by log T:

   | log T | source |
   |---|---|
   | ≥ T_c | fully ionized (`fully_ionized_gas`: 1/mu_e from X, then `compute_fully_ionized_eos`) |
   | T_c − 0.5 < log T < T_c | the Saha (or SCV) result blended into the fully ionized one (`blend_saha_into`), weight f = w³(6w² − 15w + 10), w = 2(T_c − log T) |
   | ≤ T_c − 0.5 | the Saha equations (`compute_saha_eos`) |

   With `use_scv_eos`, below T_c SCV's density and first derivatives replace Saha's where its table covers the point (Saha's second derivatives stay; a point outside the SCV table keeps Saha's result).
4. **OPAL**, when a vintage is selected (`opal_vintage`) and the point is inside its (T, rho) window (`opal_layout`): `compute_opal_eos` finds the density that gives P and the derivatives there; where the table does not cover the point the result so far stands; inside the table OPAL's density, beta, mean weights and first derivatives replace it (`adopt_opal`), and near the table's edge (`find_opal_edge_ramp`: one table cell wide) they are blended in linearly (`blend_opal_into`). The second derivatives stay the Yale/SCV ones.

So a deck with SCV and OPAL 2006 (the reference solar decks) gets OPAL wherever its table covers the point (log T from 3.27 up), SCV at low T outside it, the Saha solution in between, and the fully ionized gas at high density (log rho > 7) and in the stellar centre above the OPAL table.

### 2.2 compute_mhd_eos

The MHD tables (three ZAMS-type compositions, each in a lower and an upper part in T, and five centre-type ones) are interpolated bicubically in (log P, log T) in each table of the family the temperature selects, then in X (quadratically across the ZAMS compositions, cubically across four of the centre ones). A point outside the tables returns ierr = 1. The tables are not shipped.

## 3. Conventions

- **Units** are cgs: T in K, P in dyn/cm², rho in g/cm³, cp in erg/(g K), and log is log₁₀. The OPAL routines work internally in T6 = T/10⁶ and P in units of 10¹² dyn/cm².
- **Derivatives are logarithmic.** The first derivatives are ∂lnρ/∂lnT at constant P and ∂lnρ/∂lnP at constant T. The six second derivatives are (test_eos's derivative check compares each with finite differences):

  | Output | Definition |
  |---|---|
  | `dlnrho_dlnt_dt` | ∂(∂lnρ/∂lnT)/∂lnT |
  | `dlnrho_dlnt_dp` | ∂(∂lnρ/∂lnT)/∂lnP |
  | `adiabatic_gradient_dt`, `_dp` | ∂ln∇ad/∂lnT, ∂ln∇ad/∂lnP |
  | `specific_heat_cp_dt`, `_dp` | ∂ln cp/∂lnT, ∂ln cp/∂lnP |

- **`want_derivatives`** asks for the second derivatives; without it they keep the caller's values. Inside an OPAL or SCV table they are the Yale (Saha) solution's, not the table's, so they can be far from those of the EOS actually used; `use_numerical_derivatives` makes them consistent.
- **`in_atmosphere`** skips the Saha solver's second derivatives (they keep the caller's values); everything else, the first derivatives included, is computed as without it.
- **beta** is P_gas/P, and `beta14` is 4(1 − beta)/beta = 4 P_rad/P_gas, the radiation term of the Yale solvers, on every path.
- **Inputs carried in the results**: the density is the starting guess of the fully ionized density iteration, so a result depends on the previous call's density within the iteration's tolerance. The second derivatives are inputs as above. Every other field is overwritten; log T and log P come back unchanged. `saha_iterations` is a counter the Saha solver adds its iterations to; nothing in eos reads it.
- **Errors** return ierr (positive): a missing or corrupt table, a bad configuration (two OPAL vintages), beta ≤ 0, an MHD point outside its tables, an OPAL interpolation error. A point outside an OPAL or SCV table is not an error: the next source covers it.
- **Messages** go through `eos_messages`: a warning or a diagnostic the run continues after (a Saha or density iteration that did not converge, an OPAL lookup outside its table, a density search that falls back) to the run log (`eos_log`); an error to the terminal and the run log (`eos_report`). Each is formatted by its caller into lines (an internal write, so the formats are the original ones) and printed without trailing blanks; no other eos module writes to a unit.

## 4. Configuration

`eos_settings_t` (in `eos_configuration`) is everything eos reads: MHD; the OPAL vintage (`opal_vintage`: `opal_vintage_none`, `_1995`, `_2001`, `_2006`); SCV; numerical derivatives; diffusion of Z; the Debye-Hückel correction and its degeneracy range; the Saha cutoff. `eos_table_paths_t` holds the tables' paths (Fermi-Dirac, SCV, MHD, the OPAL table of the selected vintage). The star layer fills both from the controls (`equation_of_state_settings` and `equation_of_state_table_paths` in `setup/microphysics_settings.f90`; the controls let the newest vintage win, and settings that still select two are `opal_vintage_conflict`, which `eos_init` refuses) and passes them to `eos_init` once per run (`core/run/run_setup.f90`). eos opens every table itself with `newunit`; the OPAL table is read by its vintage's reader on first use, the others at `eos_init`.

## 5. State

What eos keeps between calls, module by module, and who sets it:

| State | Module | Set by | Lifetime |
|---|---|---|---|
| `eos_settings` | `eos_configuration` | `eos_init`, `set_eos_settings` | the run |
| the envelope mixture `eos_mix` | `eos_mixture` | `eos_set_mixture` (the starting model, and the rescaling of Z) | until the next call |
| the Debye-Hückel composition (X, Y, Z, the 18 metals) | `yale_eos_tables` | the metals beyond O by `eos_set_mixture`; X, Y, Z, C, N, O per zone by `eos_get`'s `composition_at_zone` (omitted, the last one stays) | until the next call |
| the Fermi-Dirac table | `yale_eos_tables` | `eos_init` | the run |
| the Saha solver's starting guess (E, x1, x2) | `yale_eos` (`saha_e_guess`, ...) | every Saha solve: the previous call's solution | never reset |
| the SCV tables, and the envelope table | `scv_eos_tables` | `eos_init`; `eos_update_envelope_tables` (the starting model). `eos_set_mixture` alone (the rescaling of Z) does not rebuild it, and the lookup tests the zone against `eos_mix`, so the two must be called as a pair | the run; until the next update |
| SCV's warm-start cell (`idtt`, `idp`) | `scv_eos_tables` | every SCV lookup (and the envelope table's build, as loop scratch) | never reset |
| the OPAL tables, their loaded flags and Z, the open table's unit | `opal_eos_tables` (`opal95`, `opal01`, `opal06`, `opal_shared`) | the vintage's reader on first use; `eos_init` opens the unit | the process (a reset keeps the loaded tables) |
| OPAL's interpolation stencil, quadratic coefficients (`quad` caches) and last results (`eos_output`) | the vintage's instance | every esac call; `find_opal2006_ramp` reuses the last stencil, `find_opal_ramp` its cached row (`t_row_index`) | never reset |
| the MHD tables and the last interpolation (`mhd_output`) | `mhd_eos_tables` | `eos_init`; every MHD call | the run |

## 6. Tests and coverage

`test_eos` (`make tests` in `src/`, then `pytest test_eos.py` from the repository root) runs without a star: it calls `eos_init` itself, sets the envelope mixture (core's algorithm on the default mixture seed, X = 0.70, Z = 0.016492) and builds the SCV envelope table. Its output is byte-compared with `eos/test/expected_test_eos.out`:

- pinned values: the production EOS (SCV and OPAL 2006) on a (log T, log P) grid; `eos_get_gamma1` on both branches; the error paths; every configuration a deck can choose (OPAL 1995 and 2001, SCV only, Yale only, numerical derivatives, diffusion of Z, Debye-Hückel); the derivatives against finite differences; the MHD table lookup on a synthetic table. With `YREC_MHD_TABLES` set to a directory with the real MHD tables it also evaluates the MHD path.
- ok/FAIL checks, each confirmed to fail under a mutation of the code it covers (a failed check makes the program exit 1):
  - the MHD reader on eight synthetic unformatted tables (an ideal gas in every row, so the interpolation is exact): the density in each family and at the ZAMS/centre boundary, the tabulated gradients, ierr outside the tables in T and in P, a missing and a truncated file;
  - the ideal-gas limits: a neutral gas (rho = mu P_gas/(R T), grad_ad → 0.4, cp → 5R/(2 mu)), a weakly degenerate fully ionized one (rho = mu P_gas/(R T)·(1 − 2^(−5/2) e^η P_e/P_gas), the first-order Fermi-Dirac correction, to 2e-4), grad_ad and cp against the ideal gas plus radiation, grad_ad → 1/4 when radiation dominates;
  - beta ≤ 0 (radiation pressure above P) returns ierr;
  - continuity at the foot of the Saha ramp and at the cutoff of every blended quantity it compares (rho, grad_ad, the ionization fractions, two second derivatives), and of rho and grad_ad at both ends of the OPAL 2006 edge ramp (each end bisected to 1e-11 in log P);
  - `in_atmosphere`: the density, ionization and first derivatives of a full evaluation;
  - the SCV envelope table rebuilt after a change of composition agrees with SCV's general path;
  - the Debye-Hückel correction raises rho, the more with heavier metals in the envelope mixture;
  - `eos_get_gamma1` below the OPAL table takes the Yale/SCV value;
  - with `use_diffusion_z`, at the envelope's Z, the default mixture's 1/mu_ion and rho;
  - `eos_init` refuses settings that select two OPAL vintages.

A deliberate change to the output means regenerating `eos/test/expected_test_eos.out`. The decks pin the full EOS in use (every shipped deck uses OPAL 2006, SCV and the Debye-Hückel correction; none uses MHD, OPAL 1995 or 2001).

`python3 src/tools/coverage.py --work DIR --dir eos --test test_eos --deck ...` measures the line coverage of eos/ (a -O0 --coverage build in DIR). With test_eos and four decks that reach the OPAL, SCV and Yale paths (`examples/run_config_matrix/cm_base.inlist`, `testsuite/Test_solarGS98_base.inlist` and `Test_solarGS98_yaleeos.inlist`, `examples/run_from_dbl_to_zams/Test_m0030_feh0_noGS_ZAMS.inlist`) it is 85.1% of 3365 lines (2026-10-29): `blended_eos` and `eos_mixture` 100%, the Yale and SCV solvers 96.9%, MHD 87.4% (the synthetic tables), the OPAL code 71-83%; what is not run is mostly the readers' and interpolators' error branches and their messages.

## 7. Extending

**Add a control:** add the field to `eos_settings_t`, copy it from `star%ctrl` in `equation_of_state_settings` (`setup/microphysics_settings.f90`; then `tools/check_star_contracts.py --fix`), and read `eos_settings%<field>`. The control itself comes through the registry (`defaults/controls_registry.tsv` and the generators, see `io/README.md`).

**Add an EOS source** (a table or a solver that covers part of the (T, rho) plane and is blended in like OPAL):
1. `eos/<name>/<name>_eos.f90`: a module with `read_<name>_tables(path, ierr)` (open on a newunit with iostat; a missing or corrupt file is ierr = 1, with a message through `eos_messages`) and `compute_<name>_eos(log10_temperature, log10_pressure, X, Z, r, ierr, covers)`, which fills an `eos_result` (density, beta, the mean weights, the four first derivatives; `covers` false where its table does not reach), its tables in `eos/<name>_eos_tables.f90`. The Makefile finds new files; run `python3 tools/gen_makefile_deps.py` in `src/`.
2. `eos/eos_configuration.f90`: its switch in `eos_settings_t` and its path in `eos_table_paths_t`.
3. `eos/eos_lib.f90`: `eos_init` reads its tables when the switch is on.
4. `eos/blended_eos.f90`: a step in `blend_eos` (before the OPAL step if OPAL should win where both cover the point, after it otherwise): inside the source's window, call `compute_<name>_eos`, then `adopt_table_result` where it covers the point, or `blend_table_result` with its edge weight (1 inside, falling continuously to 0 at the table's edge, as `find_opal_edge_ramp`'s does). A source that replaces everything instead, as MHD does, is a branch in `eos_get`.
5. `setup/microphysics_settings.f90`: copy the control in `equation_of_state_settings` and the path in `equation_of_state_table_paths`; then `python3 tools/check_star_contracts.py --fix`.
6. The control and the path: the recipe in `io/README.md` §5 (`io/controls_buffer.f90`, `defaults/controls_registry.tsv`, the generators, the controls dumps); a table path is also a file-name control (as `ffermi`, `io/legacy_controls.f90`), a `star%job` member (`star_job` in `state/star_info_lib.f90`, beside `fermi_table_path`) and its copy in `io/read_controls.f90` (as `fermi_table_path`'s).
7. `eos/test/test_eos.f90`: an option of `configure_eos` (its switch and table path), a configuration in `check_configurations`, and ok/FAIL checks that a mistake would break: a limit it must reproduce, continuity at both ends of its edge ramp (`ramp_end_is_continuous`); regenerate `eos/test/expected_test_eos.out`.
8. This README: §1, §2.1, §5.

**Add a blended region** (a new ramp between two sources, like the Saha ramp):
1. `eos/blended_eos.f90`: in `blend_eos`'s gas step, keep a copy of the record with the first result (as `saha = r` does) before computing the second into it, and blend with a weight that runs from 1 at one end of the region to 0 at the other: `blend_saha_into` shows the smoothstep (continuous with its first two derivatives at both ends) and every field to blend (the log density, then the density, R/mu and 1/mu_e from it; the ionization fractions; the four first derivatives; with `want_derivatives`, the six second derivatives).
2. The region's position and width: named parameters with their units, or fields of `eos_settings_t` (`eos/eos_configuration.f90`) copied from controls in `setup/microphysics_settings.f90` (the control as in step 6 above).
3. `eos/test/test_eos.f90`: `continuity_in_t` at both ends of a region in log T (`ramp_end_is_continuous` bisects for the ends of one in log P): it compares the density, grad_ad, the ionization fractions and two second derivatives across each end, and fails when a field is left out of the blend.
4. This README: §2.1.

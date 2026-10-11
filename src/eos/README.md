# eos/ — the equation of state

This directory computes the density of a zone and its thermodynamic quantities (beta, the ionization fractions, the mean molecular weights, cp, the adiabatic gradient) and their derivatives, from (log T, log P, X, Z). Code outside eos/ uses it only through `eos_lib` (enforced by `tools/check_boundaries.py`). eos/ reads nothing from `star`: its configuration is one settings record, installed at the start of every run.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `eos_lib` (`eos_lib.f90`) | **The public interface.** `eos_init` installs the settings and loads the tables they select. `eos_get` evaluates a zone into an `eos_result` record. `eos_set_mixture`, `eos_set_debye_huckel_z`, `eos_set_debye_huckel_composition` and `eos_update_envelope_tables` take the compositions the star layer owns. `eos_get_gamma1` returns gamma1 and the adiabatic gradient at (T, rho), for test_eos. |
| `eos_results` | `eos_result`, the record every EOS family fills: one named field per quantity, each with its definition and unit. |
| `eos_configuration` | The settings record (`eos_settings_t`, `eos_settings`). |
| `eos_mixture` | The envelope mixture (`eos_mix`), `ion_mean_weight_excess` and `compute_mean_weights`. |
| `blended_eos` (`blended_eos.f90`) | The EOS of every run without MHD (section 2): `blend_eos` evaluates one point; `compute_blended_eos` adds the numerical derivatives. |
| `yale_eos` (`yale/`) | The Yale EOS of Prather (1976): `compute_saha_eos` (partial ionization), `compute_fully_ionized_eos` (degenerate, relativistic electrons, optional Debye-Hückel correction) and `read_fermi_table`. |
| `scv_eos` (`scv/`) | The Saumon-Chabrier-Van Horn EOS: `read_scv_tables`, `build_scv_envelope_table` (the H, He and metal tables combined at the envelope composition), `compute_scv_envelope_eos` and `compute_scv_eos`. |
| `opal_eos` (`opal/`) | The OPAL EOS: `compute_opal_eos` (density and derivatives at P, for any vintage) and `find_opal_edge_ramp` (whether the vintage's table covers the point, and its weight near the edge), and the 2001 and 2006 tables, which share one state layout. Where the two differ (table size, slot numbers, messages, the reader's columns, the density search and radiation correction), the code branches on the vintage. OPAL's own routines keep their names (`esac01`/`esac06` are entry points to `interpolate_opal_eos`; `rhoofp01/06`, `radsub01/06`). |
| `opal95_eos` (`opal/`) | The OPAL 1995 tables, with their own state layout: `find_opal95_ramp`, `rhoofp`, `esac`, `readco`, `radsub`, ... |
| `mhd_eos` (`mhd/`) | The Mihalas-Hummer-Däppen EOS: `read_mhd_tables` and `compute_mhd_eos`. Däppen's routines keep their names (`mhdst`, `mhdpx`, ...). |
| `opal_eos_tables`, `mhd_eos_tables`, `scv_eos_tables`, `yale_eos_tables` | Each family's tables and state. `opal_eos_tables` also holds the helpers the OPAL vintages share (`interpolate_opal_t6_rho`, `find_opal_ramp`, `quad`, `gmass01`, `open_opal_eos_table`). |

## 2. What eos_get does

With `use_mhd_eos` it calls `compute_mhd_eos`. Otherwise it calls `compute_blended_eos`, which calls `blend_eos` once, or four more times for finite-difference derivatives when `use_numerical_derivatives` is set. For one (log T, log P, X, Z), `blend_eos`:
1. **Below the ionization cutoff** (`saha_log10t_cutoff`): solves the Saha equations (`compute_saha_eos`), or, with `use_scv_eos`, looks the point up in the SCV tables and takes from Saha only what SCV lacks.
2. **Above it**: treats the gas as fully ionized (`compute_fully_ionized_eos`). Near the cutoff the two results are interpolated.
3. **OPAL**: if an OPAL vintage is selected and its table covers the point, adopts the OPAL result, or blends it in with the weight `find_opal_edge_ramp` returns near the table's edge.

## 3. Conventions

- **Units** are cgs: T in K, P in dyn/cm², rho in g/cm³, cp in erg/(g K). The OPAL routines work internally in T6 = T/10⁶ and P in units of 10¹² dyn/cm².
- **Derivatives are logarithmic.** The first derivatives are ∂lnρ/∂lnT at constant P and ∂lnρ/∂lnP at constant T. The six second derivatives are (test_eos's derivative check confirms each against finite differences on the Yale paths):

  | Output | Definition |
  |---|---|
  | `dlnrho_dlnt_dt` | ∂(∂lnρ/∂lnT)/∂lnT |
  | `dlnrho_dlnt_dp` | ∂(∂lnρ/∂lnT)/∂lnP (named `dlnrho_dlnp_dt` until 2026-10) |
  | `adiabatic_gradient_dt`, `_dp` | ∂ln∇ad/∂lnT, ∂ln∇ad/∂lnP |
  | `specific_heat_cp_dt`, `_dp` | ∂ln cp/∂lnT, ∂ln cp/∂lnP |

  They are set only with `want_derivatives`; otherwise they keep the caller's values. Inside an OPAL table they are the Yale solution's, not OPAL's: OPAL supplies density and the first derivatives only, so the second derivatives can be far from those of the EOS actually used (up to the wrong sign). `use_numerical_derivatives` makes them consistent.
- **beta** is P_gas/P. `beta14` is 4(1−β)/β on the Yale path but 1−β after OPAL or MHD (an open item in docs/TARGET_ARCHITECTURE.md).
- **Inputs carried in the results**: the density is the starting guess of the fully ionized density iteration, so a result depends on the previous call's density within the iteration's tolerance (relative changes up to ~1e-4 in test_eos). `saha_state` carries the Saha solver's state. Every other result is overwritten; log T and log P come back unchanged.

## 4. Configuration

`eos_settings_t` (in `eos_configuration`) is everything eos reads: which EOS (MHD; the OPAL vintage; SCV), numerical derivatives, diffusion of Z, the Debye-Hückel correction and its degeneracy range, and the Saha cutoff. The star layer fills it from the controls (`equation_of_state_settings` in `setup/microphysics_settings.f90`) and passes it to `eos_init`, together with the table paths from `star%job`. eos opens every table itself with `newunit`. The OPAL tables are read by their vintage's reader on first use; the others are read at `eos_init`.

## 5. Tests and extending

`test_eos` (`make tests` in `src/`, then `pytest test_eos.py` from the repository root) pins the production EOS (SCV and OPAL 2006) on a (log T, log P) grid, `eos_get_gamma1` on both branches, the error paths, and every configuration a deck can choose: OPAL 1995 and 2001, SCV only, Yale only, numerical derivatives, diffusion of Z and the Debye-Hückel correction. The MHD tables are not shipped; with `YREC_MHD_TABLES` set to their directory, test_eos also evaluates the MHD path. A deliberate change to its output means regenerating `eos/test/expected_test_eos.out`.

**Add an EOS family:** a module in its own directory with `read_*_tables` (called from `eos_init`) and `compute_*_eos`; its state in a `*_tables` module; a `use_*` flag in `eos_settings_t` (and in `equation_of_state_settings`); a branch in `eos_get` or `blend_eos`; a `compute_*_eos` that fills an `eos_result`; and a `check_configurations` section in `test_eos`.

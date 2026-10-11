# kap/ — the opacities

This directory computes the Rosseland mean opacity of a zone, and its logarithmic derivatives, from tables. Code outside kap/ uses it only through `kap_lib` (enforced by `tools/check_boundaries.py`), and kap/ reads nothing from `star`: its configuration is one settings record, installed at the start of every run.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `kap_lib` (`kap_lib.f90`) | **The public interface.** `kap_init` installs the settings and loads the tables they select; `kap_get` and `kap_eval` return the opacity of a zone; `kap_update_surface_tables` rebuilds the surface-composition tables when the envelope X changes; `set_kap_settings` installs settings without reloading. `kap_eval` holds the blending logic (section 2). |
| `opacity_tables` (`opacity_tables.f90`) | The state: every family's tables (`opacity_table`), the table dimensions, the composition-cache tolerances, the shared four-point locators, and the settings record (`kap_settings_t`, `kap_settings`). |
| `opal95_opacity` (`opal95/`) | Atomic opacities in the OPAL 1995 layout (the OPAL and OP tables, e.g. `GS98.OP17`): 126 tables in X and Z. Interpolation in 2D on the surface table, or 3D (X, T, R) and 4D (Z, X, T, R) away from it. |
| `opal92_opacity` (`opal92/`) | The 1992 OPAL tables, one file per Z (two when blending). |
| `laol_opacity` (`laol89/`) | The Los Alamos LAOL89 tables (one or two Z), and the pure-Z table (X = 0, Y + Z = 1) that helium-burning cores blend toward. |
| `alex06_opacity` (`alex06/`) | Low-temperature molecular opacities of Ferguson et al. (2005), the "Alexander 2006" tables. |
| `alex94_opacity` (`alex94/`) | Molecular opacities of Alexander & Ferguson (1994/95): seven files, one per X. |
| `kurucz_opacity` (`kurucz90/`) | Molecular opacities of Kurucz (1990), one table per Z. |
| `conductive_opacity` (`conductive/`) | Electron-conduction opacities from Potekhin et al. (1999), with Potekhin's interpolation routines under their own names (`condopacp`, `cinterp3`, `hunt`). |

Each family module has the same three kinds of routine: `read_*_tables` (at `kap_init`), `build_*_surface_table` where the family has one (at init and when the envelope X changes), and `compute_*_opacity` (per zone).

## 2. What kap_eval does

For one zone (log ρ, log T, X, Z):
1. **Molecular opacity** below `molecular_opacity_logt_max`: Alexander 2006, else Alexander 1994/95, else Kurucz 1990 (blended between two Z tables if `use_two_z_tables`).
2. **Interior opacity** unless log T < `molecular_opacity_logt_min` and a molecular value exists: OPAL95-format, else OPAL92, else LAOL89 (two Z tables blended when set). In a helium-burning core (Z > 0.1, log T > 7) the interior value is blended in Z toward the pure-Z table, which plays the second table at Z = 1. Outside the tables, or with no table chosen, `ierr = 1`.
3. **The ramp**: between `molecular_opacity_logt_min` and `_max` the opacity is the linear blend of the two.
4. **Conduction**: the Potekhin opacity when `use_conductive_opacity`, else Hubbard & Lampe where it applies, added harmonically (1/κ = 1/κ_rad + 1/κ_cond).

## 3. Conventions

- **Opacity** κ in cm²/g; `log10_opacity` is log₁₀ κ.
- **Derivatives are logarithmic:** dln κ/dln ρ and dln κ/dln T.
- **log R** = log ρ − 3 log T₆ is the tables' density coordinate.
- **ion_fraction** (3 numbers, the ionization fractions of H, He and the metal proxy, from the EOS) only affects the conductive opacity.

## 4. Configuration

`kap_settings_t` (in `opacity_tables`) is everything kap reads: which tables, the tables' Z values, the ramp limits, the Z tolerance of the surface table and the LAOL89 extrapolation limit. The star layer fills it from the controls (`opacity_settings` in `setup/microphysics_settings.f90`) and passes it to `kap_init`, together with the table paths from `star%job`. Every table file, the conductive one included, is opened by kap itself with `newunit`.

## 5. Tests and extending

`test_kap` (`make tests` in `src/`, then `pytest test_kap.py` from the repository root) pins the OPAL95 interior on a (log T, log ρ) grid, the error paths, and every family: the production combination (Alexander 2006, OPAL95-format, pure Z, Potekhin) at and away from the surface table and on a helium core; OPAL92 and LAOL89 with two Z tables; Kurucz and Alexander 1994/95. A deliberate change to its output means regenerating `kap/test/expected_test_kap.out`.

**Add a table family:** a module in its own directory with `read_*_tables`, `compute_*_opacity` and, if it caches a composition slice, `build_*_surface_table`; its tables in `opacity_table`; a `use_*` flag in `kap_settings_t` (and in `opacity_settings`); a branch in `load_opacity_tables`, `refresh_surface_tables` and `kap_eval`; and a `check_families` section in `test_kap`.

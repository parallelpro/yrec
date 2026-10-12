# atm/ — the atmosphere

This directory gives the stellar model its surface boundary condition: the pressure at the photosphere, either from a table of model atmospheres or from an analytic T-tau relation that the envelope integrator integrates. Code outside atm/ uses it only through `atm_lib` (enforced by `tools/check_boundaries.py`). atm/ reads nothing from `star`: its configuration is one settings record, installed at the start of every run.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `atm_lib` (`atm_lib.f90`) | **The public interface.** `atm_init` installs the settings and reads the tables they select. `atm_get_surface_pt` returns log P (and log T) at (log Teff, log g) from the table of the caller's current `atm_choice`. `atm_allard_teff_max` gives the Allard tables' Teff limit. The T-tau relations below are re-exported here. |
| `atm_tables` (`atm_tables.f90`) | The settings record (`atm_settings_t`, `atm_settings`) and the tables (`atm_table`): the Kurucz and Castelli/Kurucz surface-pressure tables (`surface_p_table`) and the Allard tables. |
| `kurucz_atmosphere` (`tables/`) | `read_surface_pressure_table` and `interpolate_surface_pressure`, for the Kurucz (1993) and Castelli & Kurucz (2003) tables of log P at (log Teff, log g). |
| `allard_atmosphere` (`tables/`) | The Allard model atmospheres (NextGen, AMES-Cond, AMES-Dusty): `read_allard_tables` (both file layouts), `find_allard_table_bounds` and `interpolate_allard_surface`, which can also return P and T at tau = 100. |
| `t_tau_relations` (`t_tau_relations.f90`) | The analytic T-tau relations: Eddington gray, Krishna-Swamy (1966) and the HSRA fit. Each gives log T(tau), the starting temperature of the atmosphere integration and its photospheric limit. |

The atmosphere integration itself (`core/envelope/envelope_integration_kernel.f90`, `core/envelope/envelope_integrands.f90`) belongs to the star layer. It calls these services, and its results (`atmosphere_structure`, `envelope_structure`, the current point `atm_point`) live in core/.

## 2. The atmosphere choices

| `atm_choice` | Surface boundary condition |
|---|---|
| 0 | Eddington gray T-tau relation, integrated to tau = 2/3 |
| 1 | Krishna-Swamy T-tau relation, integrated to tau = 0.312 |
| 2 | HSRA T-tau relation (legacy) |
| 3 | Kurucz (1993) table |
| 4 | Allard tables (with the Kurucz table loaded as well) |
| 5 | Castelli & Kurucz (2003) table |

The star layer switches a run between its tabulated choice and gray as Teff leaves and re-enters the table (`core/structure/surface_boundary.f90`), so `atm_get_surface_pt` takes the current choice as an argument. When an Allard lookup falls outside the tables, `lookup_failed` is set and the integrator falls back to gray for that model.

## 3. Conventions

- log is log₁₀; P in dyn/cm², T in K, g in cm/s².
- A table lookup returns log P at T = Teff and log T = log Teff, except Allard with `allard_use_tau100`, which returns log P and log T at tau = 100.
- The outputs of `atm_get_surface_pt` are set only when the lookup succeeds.

## 4. Configuration

`atm_settings_t` (in `atm_tables`) is everything atm reads: which tables to load (`atm_choice` at startup), whether the Allard lookups return the tau = 100 point, and the [Fe/H] and [α/Fe] of the Allard table to read. The star layer fills it from the controls (`atmosphere_settings` in `setup/microphysics_settings.f90`) and passes it to `atm_init`, together with the table paths from `star%job`. atm opens every table itself with `newunit`.

## 5. Tests and extending

`test_atm` (`make tests` in `src/`, then `pytest test_atm.py` from the repository root) pins the lookups of every table: Kurucz, Castelli/Kurucz, Allard NextGen 2 with and without the tau = 100 point, the old NextGen 1 layout, and a lookup outside the Allard tables. It also pins the three T-tau relations on a grid of optical depths, and the error paths. A deliberate change to its output means regenerating `atm/test/expected_test_atm.out`.

**Add a table of model atmospheres:** a module in `tables/` with a reader (called from `atm_init`) and an interpolator; its table in `atm_tables`; a new `atm_choice` value with a branch in `atm_init` and `atm_get_surface_pt`; and a section in `test_atm`. **Add a T-tau relation:** a branch in each of the three functions of `t_tau_relations`, and a line in `test_atm`'s T-tau loop.

# element_diffusion/ — gravitational settling and microscopic diffusion

This directory moves the elements through the radiative interior over a timestep: gravitational settling and the microscopic diffusion of helium, the metals and the light elements, in the Bahcall & Loeb (1990) formulation with the Thoul, Bahcall & Loeb (1994) coefficients.

Code outside element_diffusion/ uses it only through `element_diffusion_lib` (enforced by `tools/check_boundaries.py`). It has two callers, and they call it the same way:
- core's `composition_step` settles every model with settling on, rotating or not;
- rotation's `burn_and_settle` (in `mid_timestep_model`) settles within each rotation sub-step.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `element_diffusion_lib` | **The public interface.** `diffuse_elements` runs one settling step: the newer routines (`microdiff`) when `use_new_diffusion_routines` is set, else the original ones (`gravitational_settling`). `reset_element_diffusion` restores the domain's working state at the start of each run (called by yrec_reset). |
| `original_settling` | The original routines: `gravitational_settling` and its setup of the diffusion coefficients (the Thoul et al. coefficients, their fit, or the older analytic ones). |
| `settling_grid` | The original routines' mapping of the model onto an equally spaced grid in radius and back. |
| `microdiff_driver` | The newer routines: `microdiff` and its setup (the diffusion region, the conversion to Bahcall & Loeb units). |
| `microdiff_mapping` | microdiff's equally spaced grids (zone centers and midpoints, the `microdiff_grid` type) and the mapping of the abundance changes back. |
| `microdiff_run`, `microdiff_coefficients` | microdiff's solve for one species, and its diffusion coefficients. |
| `diffusion_scheme` | The numerical scheme both pipelines use: the two-step Lax–Wendroff advection of the first term and the coefficients of the implicit solve of the second. |
| `thoul_solver` | `thoul_diffusion`: the solution of the Burgers equations for the coefficients (Thoul et al. 1994). |
| `species_table`, `bahcall_loeb_units` | The species' weights and charges; the Bahcall & Loeb unit scales. |
| `element_diffusion_scratch` | The working state (`diff_scr`): coefficient arrays, the unit scales, the "settling suspended" latch. |

## 2. Conventions

- Arrays run from the center (zone 1) out. The structure arrays are the caller's: log radius, density and temperature, the enclosed mass in grams (unlogged), d ln P/dr, the temperature gradient and the convective flags.
- The routines convert the timestep, the total mass, the enclosed mass and d ln P/dr to Bahcall & Loeb units in place and convert them back before returning. The round trip can change the last bit of the enclosed mass (`star%m`); that is the existing behaviour.
- Settling acts between any convective core and the surface convection zone, inside the outermost zone where Y exceeds `helium_diffusion_min`. With Y or X below its floor everywhere, there is no settling, with one run-log note per suspension.
- `use_thoul_fit` and `use_thoul_diffusion` act only in the original routines; microdiff always uses the Thoul et al. coefficients. Light elements settle only in microdiff (`diffuse_lithium`).
- A singular matrix (the Thoul solve or the tridiagonal solve) returns `ierr`. An iteration that misses its tolerance writes a warning and continues.

## 3. Tests and extending

`test_element_diffusion` (`make tests` in `src/`, then `pytest test_element_diffusion.py` from the repository root) runs one 10⁸-year step of both routines on a synthetic 200-zone solar-like model:
- the original routines with helium only, with metals, each Coulomb logarithm, the mixing and settling factors, the scale factors, the Thoul fit, the older coefficients, and a solve stopped short of convergence;
- microdiff with helium only, with metals, and with the light elements;
- for both, the abundance floors.

No deck uses the original routines, so this test is their only pin. The decks pin microdiff: the solar GS decks and the settling and rotating cases of the configuration matrix. A deliberate change means regenerating `element_diffusion/test/expected_test_element_diffusion.out`.

**Add a settling treatment:** put it in its own module beside the two pipelines, give it a control, branch on it in `diffuse_elements`, and add its cases to `test_element_diffusion`.

# convection/ — convection physics

This directory decides where the star convects and what its temperature gradient is there. It holds:
- the gradient from mixing-length theory;
- the convection zones from the shells' convective flags, with overshoot and semiconvection;
- the homogenization of the zones' composition.

Code outside convection/ uses it only through `convection_lib` (enforced by `tools/check_boundaries.py`). convection/ reads nothing from `star`. Its configuration is one settings record, installed by the star layer at the start of every run card.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `convection_lib` (`convection_lib.f90`) | **The public interface.** It re-exports the four procedures below. `set_convection_settings` installs the settings, and `set_convection_envelope_metals` passes the starting model's envelope metal fraction outside C, N and O. |
| `convection_configuration` | The settings record (`convection_settings_t`) and the installed copy. |
| `temperature_gradients` | `compute_temperature_gradients`: the radiative gradient and, where it exceeds the adiabatic one, the convective gradient from Böhm-Vitense mixing-length theory (the MLT cubic), the convective velocity, and the gradient's log derivatives for the Henyey Jacobian. It also handles starspots (Somers 2014), the adiabatic gradient in the overshoot region below the envelope, and the adiabatic shortcut above `log10_t_cut` (cores). |
| `convection_zones` | `find_convection_zones`: the zones from the convective flags, extended by step overshoot (private `extend_by_overshoot`, with the scale height from the EOS in `compute_scale_height`), merged where overshoot joins them, and the radiative regions between them. Also `homogenize_convection_zones`, which mixes each zone's composition. |
| `semiconvection` | `extend_by_semiconvection`: semiconvection at the edge of the convective core (Castellani et al. 1971). |

The burn-and-settle sub-stepping of rotating models (`mid_timestep_model`'s `burn_and_settle`, and the commented-out Bulirsch–Stoer path) lives in `rotation/`. core's `composition_step` drives the zones, semiconvection and homogenization for each timestep, and `core/convection_geometry.f90` turns the zone flags into boundary positions.

## 2. Settings and state

`convection_settings_t` holds:
- the mixing length α;
- `log10_t_cut`;
- the starspot parameters;
- the overshoot switches (core, envelope, intermediate zones) and their extents in pressure scale heights;
- the Woo & Demarque core-overshoot limit;
- the rotation flags, for the rule that keeps a radiative point in the overshoot region of rotating models;
- the extended-composition flag;
- L☉;
- the envelope's non-CNO metal fraction.

The star layer builds it in `convection_settings_from_star` (`setup/microphysics_settings.f90`). run_yrec installs it in `begin_kind_card`, after the card's mixing length and the Monte Carlo run's L☉ are set.

The overshoot region below the envelope is model state. `find_convection_zones` returns it (`overshoot_bottom`, `overshoot_top`, which the star layer keeps as `star%iov1`/`iov2`). The callers of `compute_temperature_gradients` pass `in_overshoot_region` for the zone they evaluate.

## 3. Conventions

- log is log₁₀; cgs units. Luminosities in the gradient formula are in solar units, linear.
- Model arrays run from the centre (zone 1) out. Zone pairs `(k,1)`/`(k,2)` are the first and last shell of the k-th region from the centre.
- Failures (an unstable square root in the MLT cubic, an EOS failure) return `ierr`.

## 4. Tests and extending

`test_convection` (`make tests` in `src/`, then `pytest test_convection.py` from the repository root) covers the gradient on hand-set EOS and opacity results:
- radiative points and the adiabatic-overshoot gate;
- a convective core;
- mixing-length convection with and without derivatives;
- both starspot modes;
- the error path.

It also covers the zones on a synthetic 40-zone model, with overshoot of every kind, a merge and the rotating rule (scale heights from the Yale EOS, Fermi table only), and homogenization. Semiconvection needs the opacity tables, so the configuration matrix's semiconvection deck pins it. A deliberate change to the output means regenerating `convection/test/expected_test_convection.out`.

**Add a convection treatment** (for example a different overshoot prescription):
1. Add its parameters to `convection_settings_t` and fill them in `convection_settings_from_star`.
2. Branch on them in the module that owns the physics.
3. Add a case to `test_convection`.

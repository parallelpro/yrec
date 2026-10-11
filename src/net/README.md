# Nuclear burning in YREC: `net/` and `core/burn_step.f90`

This note explains how YREC computes nuclear reaction rates, energy generation and composition changes: which routine does what, who calls it during a timestep, and how data flows between them. It ends with recipes for the most common extensions. File references are relative to `src/`.

Two layers are involved:

- **`net/` is the physics.** Routines here compute rates, energies and abundance changes from the values they are given. They never touch `star`: the model comes in as arguments, and the run's configuration as one record, `net_settings` (in `net_lib`), which the star layer fills at the start of every run (`configure_network`). `tools/check_boundaries.py` enforces this.
- **`core/burn_step.f90` is the star-layer driver.** It loops over the model's shells and regions, reads `star`, calls `net/`, and stores the results back in `star`.

Callers elsewhere in the code (the composition step `core/composition_step.f90`, the Henyey solver, the rotation sub-steps, the timestep limiters) go through `burn_step`, or call a `net/` routine directly with explicit arguments.

---

## 1. What lives where

### `net/net_lib.f90` (module `net_lib`)

The domain's public entry point. It also re-exports `thermal_neutrino_lib`.

| Routine | What it does |
|---|---|
| `nuclear_reaction_rates` | **The one rate evaluation.** For one zone (log ρ, log T, composition), it returns a `nuclear_rates` record: the rate of each reaction slot `r_*`, per gram per second with the abundance factors omitted (screened and scaled); dln rate/dlnT and dln rate/dln ρ; the Be7 and N15 branchings; the two Be7 channel rates and the screened hep factor. Covers screening (Salpeter weak, Graboske et al. intermediate and strong), the S-factor fits of reactions 1–7, the Caughlan & Fowler α captures and triple-α. |
| `network_rate_vector` | Calls `nuclear_reaction_rates`, zeroes any rate at or below 10⁻⁵, and repacks the result as the **composition rate vector** `rate_vec(rr_*)`, per s per amu. This is what the composition solve uses. |
| `deuterium_proton_capture_rate` | The d(p,γ)³He fit, shared by the energy and deuterium code. |
| `net_setup` | Installs the run's settings (`net_settings_t`: temperature cutoffs, screening threshold, cross-section and S-factor scales, Li/Be rate scales, neutrino switches). Called once per run by `configure_network`. |

It also holds the index families and data the network uses:
- `r_*`: reaction slots 1–13;
- `rr_*`: rows of the composition rate vector;
- `iso_*`: isotopes in the rate and screening composition;
- `s_*`: extra cross-section scale slots;
- the fit coefficients, screening factors, isotope masses and charges.

### `net/nuclear_energy.f90` (module `nuclear_energy_lib`)

| Routine | What it does |
|---|---|
| `nuclear_energy` | Energy generation of one zone, returned as a `nuclear_energy_result`. It contains ε (erg/g/s) and its **absolute** derivatives dε/dln ρ and dε/dln T for the Henyey Jacobian, plus the shares of the pp chain, Be7 branches, CNO and triple-α. It also holds the α-capture and thermal-neutrino terms, the He3 energy rates, and (when `ctrl%calc_neutrinos`) the solar neutrino fluxes and Be7 abundance. Internal steps: `fill_mass_fractions` → `nuclear_reaction_rates` → the deuterium term → `sum_energy_generation` → `chain_energy_shares` → `solar_neutrino_fluxes` → `thermal_neutrino_loss`. |

### `net/network_lib.f90` (module `network_lib`)

| Routine | What it does |
|---|---|
| `implicit_burn_step` | **The implicit network.** One backward-Euler step for one zone or region: solves X_i − X0_i − dt·A_i·Σ_j ν_ij R_j(X) = 0 by Newton iteration, with equations and Jacobian built from the reaction table. |
| `equilibrium_burn_rates` | The explicit alternative used at iteration level 1 and by the hydrogen-burning timestep limiter: dX/dt, dY/dt, dC/dt and dO/dt with He3 and the CN cycle in equilibrium. |
| `network_table_errors`, `write_network_table` | The table's self-check (every row conserves baryon number) and a printout. Both are pinned by `test_net`. |

The **reaction table** (`reactions(:)`) has one row per reaction. Each row gives:
- the rate row it reads (`rr_*`);
- a branching factor (none, f, or 1 − f for the two N15 outcomes);
- the rate law (the mass fractions that multiply the rate);
- the change list of (species, ν) pairs;
- `emitted_mass` for particles the network does not track (the neutron of C13(α,n)O16).

The solved species are, in this order, H, He3, He4, C12, C13, N14 and O16. O18 is an explicit product.

### `net/light_elements.f90` (module `light_element_lib`)

Light elements are burned in the trace approximation, separately from the network.

| Routine | What it does |
|---|---|
| `deuterium_rate` | The deuterium burning rate per shell, per s, with the cap that keeps deuterium in the surface convection zone from burning faster than the convective overturn. |
| `burn_deuterium_region` | Deuterium over a step for one shell or one convection zone: exponential sub-steps with the rate interpolated between its start- and end-of-step values, the D → He3 bookkeeping, and mass-weighted blending of accreted gas. |
| `li_be_rates` | The Li6, Li7 and Be9 proton-capture rates at T9 (Caughlan & Fowler 1988), scaled by the controls `li6_rate_scale`, `li7_rate_scale` and `be9_*_rate_scale`. |

### `net/thermal_neutrinos.f90` (module `thermal_neutrino_lib`)

Thermal neutrino losses from the fits of Itoh et al. (1996), as coded by F. X. Timmes and adapted by G. Newsham (2006).
- **YREC's wrappers:** `thermal_neutrino_losses` and `neutrino_losses_for_mixture`.
- **Timmes' routines, under their public names:** `sneut` (the emissivities), `azbar` (mean A and Z), `ifermi12` and `zfermim12` (inverse Fermi integrals).

### `core/burn_step.f90` (module `burn_step`): the star-layer driver

| Routine | Called from | What it does |
|---|---|---|
| `configure_network` | `run_yrec`, at the start of every run | Fills `net_settings` from `star%ctrl` and the rate factors (after the Monte Carlo factors are applied) and calls `net_setup`. |
| `compute_burn_rates` | `composition_step` | Each shell's rate vector (`network_rate_vector`) and deuterium rate (`deuterium_rate`), stored in `rate_by_zone`, `star%deuterium_burning_rate` (and `_start` in the predictor pass) and `rot_scr%reaction_rate_by_zone`. |
| `burn_network` | `composition_step` | `burn_network_region` per radiative shell and per mixed region. |
| `burn_network_region` | `burn_network`, `mixing/rotmix` | Model bookkeeping around `implicit_burn_step`: convection-zone averaging, the temperature cutoff, the first guess, the metallicity rule, writing `star%xa` back. |
| `burn_equilibrium` | `composition_step`, predictor pass | `equilibrium_burn_rates` per shell or region, overwriting X (and Z, C12, O16 under helium burning). |
| `burn_deuterium` | `composition_step` | `burn_deuterium_region` per shell or region. |
| `surface_accretion_now` | `burn_deuterium`, `mid_timestep_model` | The accreted gas over a step (`surface_accretion`: on, inflow, its deuterium fraction, the mass accreted) for `burn_deuterium_region`. |
| `compute_zone_energy` | `henyey_coefficients`, `neutrino_flux_table`, `timestep_limit_heburn` | `nuclear_energy` with the two inputs that depend on the shell's place in the model filled in (inside the surface convection zone? its turnover time). Returns the `nuclear_energy_result`; stores nothing. |
| `store_start_li_be_rates` | `evolve_step` | Li/Be rates of the start-of-step model → `star%rate_*_start`. |
| `store_end_li_be_rates` | `evolve_step`, `evolve_angular_momentum` | Li/Be rates of the model passed in (converged, or a rotation sub-step's midpoint) → `star%rate_*`. |
| `burn_li_be` | `evolve_step`, `evolve_angular_momentum` | Li/Be depletion over a step: radiative shells, then the surface convection zone with its base (overshoot included) located at the start and end of the step. |
| `burn_li_be_midpoint` | (compiled out) | Single-midpoint variant, used only by the Bulirsch–Stoer path; built only with `-DYREC_BULIRSCH_STOER_BURN`. |

---

## 2. When each piece runs during a timestep

`core/evolve_step.f90` advances one model. The calls that involve burning, in order:

```
evolve_step
 advance_composition_and_age
   star%xa_start := star%xa                      start-of-step composition
   composition_step(predictor_pass)              core/composition_step.f90
     find_convection_zones                       radiative shells and mixed regions
     compute_burn_rates      [burn_step]         per shell: network_rate_vector, deuterium_rate
     burn_network            [burn_step]         per region: burn_network_region -> implicit_burn_step
     burn_equilibrium        [burn_step]         predictor pass only: equilibrium_burn_rates
     (semiconvection), homogenize mixed regions, settling, renormalize
     burn_deuterium          [burn_step]         per region: burn_deuterium_region
 rezone_or_snapshot
   store_start_li_be_rates   [burn_step]         if extended composition
 solve_structure, levels 1-3: solve_level, each level
   refresh_surface_bc                            surface boundary condition
   (levels >= 3) star%xa := star%xa_start; composition_step(corrector_pass)
                                                 re-burn with the updated structure
   henyey_iterate, the Newton loop
     henyey_coefficients, every shell, every Newton iteration:
       compute_zone_energy   [burn_step]         -> nuclear_energy -> nuclear_reaction_rates
                                                 eps and d eps/dln rho, d eps/dln T enter the luminosity equation
 converge_with_rotation
   level 4 (solve_level again)
   evolve_angular_momentum (rotating runs), per sub-step:
     mid_timestep_model -> rotmix -> burn_network_region; burn_deuterium_region (the sub-step's accreted mass)
     store_end_li_be_rates (midpoint model), burn_li_be
 burn_light_elements                             non-rotating, or no rotational mixing
   store_end_li_be_rates (converged model), burn_li_be
 compute_timestep
   timestep_limit_hburn: network_rate_vector + equilibrium_burn_rates at the H-shell midpoint
   timestep_limit_heburn: compute_zone_energy at the hottest zone (triple-alpha share)
 output
   neutrino_flux_table (when solar neutrino output is on): compute_zone_energy per shell
```

Two things often surprise newcomers:
- **Composition and structure are operator-split.** Burning runs at fixed structure in `composition_step`, then the structure is solved with the new composition. At iteration levels ≥ 3, `solve_level` runs the composition step again (the corrector pass) to re-burn with the updated structure, between the surface boundary condition and the Newton loop.
- **There are two rate paths from one evaluation.** The composition path uses `network_rate_vector`: per s per amu, no derivatives, rates at or below 10⁻⁵ dropped. The energy path uses `nuclear_energy`: erg/g/s with derivatives, and only reactions 1–7 below 10⁻³⁰ dropped before ε is summed. Both call `nuclear_reaction_rates`, so a rate is defined in exactly one place.

---

## 3. Units and conventions

- **Composition.** `star%xa(i, k)` holds mass fractions in the slots `i_h1`=1, `i_he4`=2, `i_metals`=3, `i_he3`=4, `i_c12`=5, `i_c13`=6, `i_n14`=7, `i_n15`=8, `i_o16`=9, `i_o17`=10, `i_o18`=11, `i_h2`=12, `i_li6`=13, `i_li7`=14, `i_be9`=15. The light elements (12–15) are tracked with the extended composition.
- **Rates.** `nuclear_rates%rate` is reactions per gram per second with the abundance factors omitted. The composition rate vector is that times the atomic mass unit in grams, per s per amu: with Y = X/A, dY/dt = Σ ν R. Every time in the burning code is in seconds (2026-10; it used to mix years and Gyr, with a year length that differed by 4.7×10⁻⁶ between the timestep and the rates). Integer mass numbers convert X to Y.
- **Energy derivatives are absolute:** dε/dln T, not dln ε/dln T. That's what the luminosity equation needs (see the header of `core/henyey_eliminate.f90`).
- **Temperature cutoffs.** `net_settings%logT_cutoffs` (from `ctrl%nuclear_logT_cutoffs`): `(1)` gates all burning, `(2)` and `(3)` gate parts of `equilibrium_burn_rates`, `(4)` gates helium burning, and `(5)` the solar neutrino fluxes and thermal neutrino losses.

---

## 4. Extending the network

Run `test_net` after any change (`make tests` in `src/`, then `pytest test_net.py` from the repository root). It pins:
- the composition rate vector at four temperatures (log T = 6.9 to 7.6);
- the cross-section scales;
- `sneut` and `azbar`;
- the deuterium rate;
- a finite-difference check of every rate's derivatives;
- the reaction table with its baryon balance.

A deliberate change to its output means regenerating `net/test/expected_test_net.out` and committing it.

**Change a rate fit or coefficient.** Edit it in `nuclear_reaction_rates` (`net_lib.f90`), the only copy. If the functional form changes, update the matching `dlnrate_dlnt` expression next to it; `test_net`'s derivative check will tell you if they disagree. For the S-factor reactions 1–7, the fit coefficients are the `fit_q*` constants at the top of `net_lib`.

**Add a reaction to the implicit network.** Today this touches three places:
1. **The rate.** In `net_lib`, add a slot `r_*` (and raise `num_reactions`) and its formula in `nuclear_reaction_rates`, with its derivative. Add a row `rr_*` (and raise `num_rate_rows`), filled in `network_rate_vector`.
2. **The energy.** In `sum_energy_generation` (`nuclear_energy.f90`), add `eg(r_new)` and its Q-value term.
3. **The composition.** In `network_lib`'s `reactions(:)`, add a row with its rate row, rate law and change list, plus `emitted_mass` for any untracked particle. `network_table_errors` must stay 0.

A planned step (`docs/TARGET_ARCHITECTURE.md` §4.5(b)) puts the rate, Q-value and stoichiometry in a single table entry.

**Add a tracked isotope.** Today this needs:
- a slot in `state/star_info_lib.f90` (the species count is still a literal 15 in many declarations, which Phase 5 removes);
- an entry in `network_lib`'s species list (`s_*`, `mass_number`) if it is solved implicitly, and in `burn_network_region`'s `xa_slot` map;
- any reactions, as above.

**Add a light element in the trace approximation.** Follow deuterium: a rate function and a region integrator in `light_elements.f90` (pure, everything as arguments), driven per region from `burn_step`.

**Add an energy source that is not nuclear** (tidal or magnetic heating, for example). That belongs to the structure equations, not here. Phase 2 of `docs/TARGET_ARCHITECTURE.md` adds an energy-source list with a user hook. Until then the energy term is assembled in `core/henyey_coefficients.f90`; see the conventions in `core/henyey_eliminate.f90` before touching it.

---

## 5. Known limitations and open items

These are tracked in `ROADMAP.md` §13.
- **Special-purpose operators.** Deuterium, Li/Be and the equilibrium burn are not entries in the reaction table. Making deuterium a network species would change the pp-chain model (p + p → d would no longer be instantaneous).
- **Li/Be depletion mixes three jobs.** `burn_li_be` still contains model geometry (locating the convection-zone base with overshoot) and averaging alongside the depletion integration. The integration will move to `light_elements.f90` and the geometry to a shared convection-geometry module.
- **Legacy cross sections.** `LNEWNUC = .false.` takes the cross sections from the legacy `SSTANDARD` ratios. Every shipped deck uses `.true.`.

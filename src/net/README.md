# Nuclear burning in YREC: `net/` and `core/burn_step.f90`

This note explains how YREC computes nuclear reaction rates, energy generation and composition changes: which routine does what, who calls it during a timestep, and how data flows between them. It ends with recipes for the most common extensions. File references are relative to `src/`.

Two layers are involved:

- **`net/` is the physics.** Routines here compute rates, energies and abundance changes from the values they are given. They never touch `star`: the model comes in as arguments, and the run's configuration as one record, `net_settings` (in `reaction_rates`), which the star layer fills at the start of every run (`configure_network`). `tools/check_boundaries.py` enforces this.
- **`core/burn_step.f90` is the star-layer driver.** It loops over the model's shells and regions, reads `star`, calls `net/`, and stores the results back in `star`.

Callers elsewhere in the code (the composition step `core/composition_step.f90`, the Henyey solver, the rotation sub-steps, the timestep limiters) go through `burn_step`, or call a `net/` routine directly with explicit arguments. Either way they `use net_lib`, net/'s **public interface**: a module with no code of its own that re-exports what the rest of YREC uses from the five modules below, and nothing more. The individual reaction rates, the reaction table and the neutrino routines stay internal. Read `net/net_lib.f90` to see net/'s whole interface in one page; `tools/check_boundaries.py` rejects a `use` of any other net/ module from outside net/.

Names follow the project rules (`docs/GUIDELINES.md`): modules are nouns named for what they hold, procedures are verbs (`compute_…`, `burn_…`), and functions returning a quantity are named for it (`deuterium_rate`).

---

## 1. What lives where

### `net/reaction_rates.f90` (module `reaction_rates`)

The rate of every reaction in a zone, the composition rate vector, and the run's settings.

| Routine | What it does |
|---|---|
| `compute_reaction_rates` | **The one rate evaluation.** For one zone (log ρ, log T, composition), it returns a `nuclear_rates` record: the rate of each reaction slot `r_*`, per gram per second with the abundance factors omitted (screened and scaled); dln rate/dlnT and dln rate/dln ρ; the Be7 and N15 branchings; the two Be7 channel rates and the screened hep factor. In steps: composition moments → screening (Salpeter weak, Graboske et al. intermediate and strong) → the S-factor rates of reactions 1–7 → the Be7 and N15 branchings → the Caughlan & Fowler α captures and triple-α → hep. |
| `compute_network_rates` | Calls `compute_reaction_rates`, zeroes any rate at or below 10⁻⁵, and repacks the result as the **composition rate vector** `rate_vec`, per s per amu: the reactions at their `r_*` slots, then the two branching fractions at `rv_frac_n15_pa` and `rv_frac_be7_ecap`. This is what the composition solve uses. |
| `deuterium_proton_capture_rate` | The d(p,γ)³He fit, shared by the energy and deuterium code. |
| `set_net_settings` | Installs the run's settings (`net_settings_t`: temperature cutoffs, screening threshold, cross-section and S-factor scales, Li/Be rate scales, neutrino switches). Called once per run by `configure_network`. |

It also holds the index families and data the network uses:
- `r_*`: reaction slots 1–13, used by the rate record, the composition rate vector and the cross-section scales alike;
- `rv_frac_n15_pa`, `rv_frac_be7_ecap`: the two branching fractions that follow the reactions in the composition rate vector;
- `xs_pep`, `xs_be7e`, `xs_be7p`, `xs_hep`: the four cross-section scale slots after the reactions;
- `iso_*`: isotopes in the rate and screening composition;
- the fit coefficients, screening factors, isotope masses and charges.

### `net/nuclear_energy.f90` (module `nuclear_energy`)

| Routine | What it does |
|---|---|
| `compute_nuclear_energy` | Energy generation of one zone, returned as a `nuclear_energy_result`. It contains ε (erg/g/s) and its **absolute** derivatives dε/dln ρ and dε/dln T for the Henyey Jacobian, plus the shares of the pp chain, Be7 branches, CNO and triple-α. It also holds the α-capture and thermal-neutrino terms, the He3 energy rates, and (when `ctrl%calc_neutrinos`) the solar neutrino fluxes and Be7 abundance. It takes the zone's composition as one array (the `star%xa` layout) and is a short sequence of steps with explicit arguments: `compute_reaction_rates` → `energy_path_rates` → `compute_deuterium_burning` → `compute_reaction_energies` (the table) → `sum_energy_generation` → `store_chain_energies` → `compute_solar_neutrino_fluxes` → `add_thermal_neutrino_loss`. |

### `net/network.f90` (module `network`)

| Routine | What it does |
|---|---|
| `burn_implicitly` | **The implicit network.** One backward-Euler step for one zone or region: solves X_i − X0_i − dt·A_i·Σ_j ν_ij R_j(X) = 0 by Newton iteration, with equations and Jacobian built from the reaction table. |
| `compute_equilibrium_burn_rates` | The explicit alternative used at iteration level 1 and by the hydrogen-burning timestep limiter: dX/dt, dY/dt, dC/dt and dO/dt with He3 and the CN cycle in equilibrium. |
| `compute_reaction_energies` | The energy generation of the table's reactions, Σ b_j·rate·ΠX·Q_j, by rate slot and by chain. `nuclear_energy` calls it, so a reaction's composition change and its energy come from the same row. |
| `network_table_errors`, `write_network_table` | The table's self-check and a printout, both pinned by `test_net`. The check: every row conserves baryon number, and no row deposits more than the energy its atomic mass excesses release. The printout shows each row's Q from the masses and the neutrino energy that implies. |

The **reaction table** (`reactions(:)`) has one row per reaction. Each row gives:
- the rate slot it reads (`r_*`);
- a branching factor (none; f or 1 − f for the two N15 outcomes; f_e or 1 − f_e for the two Be7 branches);
- the rate law (the mass fractions that multiply the rate);
- the change list of (species, ν) pairs;
- `emitted_mass` for particles the network does not track (the neutron of C13(α,n)O16);
- `q_mev`, the energy deposited per reaction net of neutrino losses, and the `chain` it counts under (pp I, pp II, pp III, CNO, α captures, triple-α).

He3(He4,γ)Be7 appears twice, once per Be7 branch: the same abundance change, with the two branch fractions and their own Q-values.

The solved species are, in this order, H, He3, He4, C12, C13, N14 and O16. O18 is an explicit product.

### `net/light_elements.f90` (module `light_elements`)

Light elements are burned in the trace approximation, separately from the network.

| Routine | What it does |
|---|---|
| `deuterium_rate` | The deuterium burning rate per shell, per s, with the cap that keeps deuterium in the surface convection zone from burning faster than the convective overturn. |
| `burn_deuterium_region` | Deuterium over a step for one shell or one convection zone: exponential decay with the rate linear in time between its start- and end-of-step values (`surviving_fraction`), the D → He3 bookkeeping, and mass-weighted blending of accreted gas (`accreted_survival`). |
| `burn_exposure`, `surviving_fraction`, `accreted_survival` | **The trace-element integrator** shared by deuterium and Li/Be: a species destroyed at rate λ per unit abundance, with λ interpolated over the step linearly in time (`linear_in_time`) or in log λ (`linear_in_log`). The exposure Λ = ∫λ dt has a closed form for both; the species keeps exp(−Λ) of itself, and gas accreted steadily over the step keeps `accreted_survival` of its own. |
| `compute_li_be_rates` | The Li6, Li7 and Be9 proton-capture rates at T9 (Caughlan & Fowler 1988), scaled by the controls `li6_rate_scale`, `li7_rate_scale` and `be9_*_rate_scale`. |

### `net/thermal_neutrinos.f90` (module `thermal_neutrinos`)

Thermal neutrino losses from the fits of Itoh et al. (1996), as coded by F. X. Timmes and adapted by G. Newsham (2006).
- **YREC's wrappers:** `compute_thermal_neutrino_losses` and `compute_mixture_neutrino_losses`.
- **Timmes' routines, under their public names:** `sneut` (the emissivities), `azbar` (mean A and Z), `ifermi12` and `zfermim12` (inverse Fermi integrals).

### `core/burn_step.f90` (module `burn_step`): the star-layer driver

| Routine | Called from | What it does |
|---|---|---|
| `configure_network` | `run_yrec`, at the start of every run | Fills `net_settings` from `star%ctrl` and the rate factors (after the Monte Carlo factors are applied) and calls `set_net_settings`. |
| `network_rates_of_model` | `compute_burn_rates`, `mid_timestep_model` | The composition rate vector of every shell of the model passed in that is hotter than the nuclear cutoff; zero below. |
| `compute_burn_rates` | `advance_composition` | Each shell's rate vector (`network_rates_of_model`) and deuterium rate (`deuterium_rate`), stored in `rate_by_zone`, `star%deuterium_burning_rate` (and `_start` in the predictor pass). |
| `burn_regions` | `burn_network`, `burn_equilibrium`, `burn_deuterium`, `mixing/rotmix`, `mid_timestep_model` | The regions a burn visits, in order: each radiative shell above the nuclear cutoff, then each mixed region as a unit. It takes the temperatures and region tables as arguments, so the rotation sub-steps pass their midpoint model's. |
| `burn_network` | `advance_composition` | `burn_network_region` for each of the `burn_regions`. |
| `burn_network_region` | `burn_network`, `mixing/rotmix` | Model bookkeeping around `burn_implicitly`: convection-zone averaging, the temperature cutoff, the first guess, the metallicity rule, writing `star%xa` back. |
| `burn_equilibrium` | `advance_composition`, predictor pass | `compute_equilibrium_burn_rates` for each of the `burn_regions`, overwriting X (and Z, C12, O16 under helium burning). |
| `burn_deuterium` | `advance_composition` | `burn_deuterium_region` for each of the `burn_regions`. |
| `surface_accretion_now` | `burn_deuterium`, `mid_timestep_model` | The accreted gas over a step (`surface_accretion`: on, inflow, its deuterium fraction, the mass accreted) for `burn_deuterium_region`. |
| `compute_zone_energy` | `henyey_coefficients`, `neutrino_flux_table`, `timestep_limit_heburn` | `nuclear_energy` with the two inputs that depend on the shell's place in the model filled in (inside the surface convection zone? its turnover time). Returns the `nuclear_energy_result`; stores nothing. |
| `store_start_li_be_rates` | `evolve_step` | Li/Be rates of the start-of-step model → `star%rate_*_start`. |
| `store_end_li_be_rates` | `evolve_step`, `evolve_angular_momentum` | Li/Be rates of the model passed in (converged, or a rotation sub-step's midpoint) → `star%rate_*`. |
| `burn_li_be` | `evolve_step`, `evolve_angular_momentum` | Li/Be depletion over a step, with the trace-element integrator (log interpolation): radiative shells, then the surface convection zone as one unit on its mass-averaged rates (`average_li_be_over_cz`), then the shells a retreating CZ left. The CZ base at the start and end of the step, overshoot included, comes from `locate_surface_mixing_base` in `core/convection_geometry.f90`. |

---

The compiled-out Bulirsch–Stoer path keeps its own single-midpoint Li/Be variant, `burn_li_be_midpoint`, in `mixing/burn_li_be_midpoint.f90` beside the other two files of that path; build with `-DYREC_BULIRSCH_STOER_BURN` to include it.

## 2. When each piece runs during a timestep

`core/evolve_step.f90` advances one model. The calls that involve burning, in order:

```
evolve_step
 advance_composition_and_age
   star%xa_start := star%xa                      start-of-step composition
   advance_composition(predictor_pass)           core/composition_step.f90
     find_convection_zones                       radiative shells and mixed regions
     compute_burn_rates      [burn_step]         per shell: compute_network_rates, deuterium_rate
     burn_network            [burn_step]         per region: burn_network_region -> burn_implicitly
     burn_equilibrium        [burn_step]         predictor pass only: compute_equilibrium_burn_rates
     (semiconvection), homogenize mixed regions, settling, renormalize
     burn_deuterium          [burn_step]         per region: burn_deuterium_region
 rezone_or_snapshot
   store_start_li_be_rates   [burn_step]         if extended composition
 solve_structure, levels 1-3: solve_level, each level
   refresh_surface_bc                            surface boundary condition
   (levels >= 3) star%xa := star%xa_start; advance_composition(corrector_pass)
                                                 re-burn with the updated structure
   henyey_iterate, the Newton loop
     henyey_coefficients, every shell, every Newton iteration:
       compute_zone_energy   [burn_step]         -> compute_nuclear_energy -> compute_reaction_rates
                                                 eps and d eps/dln rho, d eps/dln T enter the luminosity equation
 converge_with_rotation
   level 4 (solve_level again)
   evolve_angular_momentum (rotating runs), per sub-step:
     mid_timestep_model: network_rates_of_model (midpoint model) -> rotmix -> burn_network_region;
                         burn_deuterium_region (the sub-step's accreted mass)
     store_end_li_be_rates (midpoint model), burn_li_be
 burn_light_elements                             non-rotating, or no rotational mixing
   store_end_li_be_rates (converged model), burn_li_be
 compute_timestep
   timestep_limit_hburn: compute_network_rates + compute_equilibrium_burn_rates at the H-shell midpoint
   timestep_limit_heburn: compute_zone_energy at the hottest zone (triple-alpha share)
 output
   neutrino_flux_table (when solar neutrino output is on): compute_zone_energy per shell
```

Two things often surprise newcomers:
- **Composition and structure are operator-split.** Burning runs at fixed structure in the composition step (`advance_composition`), then the structure is solved with the new composition. At iteration levels ≥ 3, `solve_level` runs the composition step again (the corrector pass) to re-burn with the updated structure, between the surface boundary condition and the Newton loop.
- **There are two rate paths from one evaluation.** The composition path uses `compute_network_rates`: per s per amu, no derivatives, rates at or below 10⁻⁵ dropped. The energy path uses `nuclear_energy`: erg/g/s with derivatives, and only reactions 1–7 below 10⁻³⁰ dropped before ε is summed. Both call `compute_reaction_rates`, so a rate is defined in exactly one place.

---

## 3. Units and conventions

- **Composition.** `star%xa(i, k)` holds mass fractions in the slots `i_h1`=1, `i_he4`=2, `i_metals`=3, `i_he3`=4, `i_c12`=5, `i_c13`=6, `i_n14`=7, `i_n15`=8, `i_o16`=9, `i_o17`=10, `i_o18`=11, `i_h2`=12, `i_li6`=13, `i_li7`=14, `i_be9`=15. The light elements (12–15) are tracked with the extended composition.
- **Rates.** `nuclear_rates%rate` is reactions per gram per second with the abundance factors omitted. The composition rate vector is that times the atomic mass unit in grams, per s per amu: with Y = X/A, dY/dt = Σ ν R. Every time in the burning code is in seconds (2026-10; it used to mix years and Gyr, with a year length that differed by 4.7×10⁻⁶ between the timestep and the rates). Integer mass numbers convert X to Y.
- **Energy derivatives are absolute:** dε/dln T, not dln ε/dln T. That's what the luminosity equation needs (see the header of `core/henyey_elimination.f90`).
- **Temperature cutoffs.** `net_settings%logT_cutoffs` (from `ctrl%nuclear_logT_cutoffs`): `(1)` gates all burning, `(2)` and `(3)` gate parts of `compute_equilibrium_burn_rates`, `(4)` gates helium burning, and `(5)` the solar neutrino fluxes and thermal neutrino losses.

---

## 4. Extending the network

Run `test_net` after any change (`make tests` in `src/`, then `pytest test_net.py` from the repository root). It pins:
- the composition rate vector at four temperatures (log T = 6.9 to 7.6);
- the cross-section scales;
- `sneut` and `azbar`;
- the deuterium rate;
- a finite-difference check of every rate's derivatives;
- the reaction table, with its baryon balance and each row's Q against the atomic mass excesses (the difference is the row's mean neutrino energy);
- `burn_implicitly` against backward Euler's closed form with only pp on, then with only C12(α,γ) on, and its mass conservation with every rate on;
- `surviving_fraction` and `accreted_survival` against closed forms and an independent quadrature;
- `compute_nuclear_energy` at the solar centre, a CNO core and a helium core, with a finite-difference check of its derivatives.

The checks print `ok` or `FAIL`.

A deliberate change to its output means regenerating `net/test/expected_test_net.out` and committing it.

**Change a rate fit or coefficient.** Edit it in `compute_reaction_rates` (`reaction_rates.f90`), the only copy. If the functional form changes, update the matching `dlnrate_dlnt` expression next to it; `test_net`'s derivative check will tell you if they disagree. For the S-factor reactions 1–7, the fit coefficients are the `fit_q*` constants at the top of `reaction_rates`.

**Add a reaction to the implicit network.** This touches two places:
1. **The rate.** In `reaction_rates`, add a slot `r_*` (and raise `num_reactions`) and its formula in `compute_reaction_rates`, with its derivative, and copy it into the rate vector in `compute_network_rates` (raise `num_rate_rows` and move the two `rv_*` fraction slots up).
2. **The table row.** In `network`'s `reactions(:)`, add a row with its rate slot, rate law, change list, `emitted_mass` for any untracked particle, Q-value and chain. The composition step and the energy generation both read it. `network_table_errors` must stay 0.

The rate formula itself still lives in `compute_reaction_rates`; TARGET_ARCHITECTURE §4.5(b) would also move it into the table entry as a rate function.

**Add a tracked isotope.** Today this needs:
- a slot in `state/star_info_lib.f90` (the species count is still a literal 15 in many declarations, which Phase 5 removes);
- an entry in `network`'s species list (`s_*`, `mass_number`) if it is solved implicitly, and in `burn_network_region`'s `xa_slot` map;
- any reactions, as above.

**Add a light element in the trace approximation.** Follow deuterium: a rate function and a region integrator in `light_elements.f90` (pure, everything as arguments), driven per region from `burn_step`.

**Add an energy source that is not nuclear** (tidal or magnetic heating, for example). That belongs to the structure equations, not here. Phase 2 of `docs/TARGET_ARCHITECTURE.md` adds an energy-source list with a user hook. Until then the energy term is assembled in `core/henyey_coefficients.f90`; see the conventions in `core/henyey_elimination.f90` before touching it.

---

## 5. Known limitations and open items

These are tracked in `docs/ROADMAP.md` §13.
- **Special-purpose operators.** Deuterium, Li/Be and the equilibrium burn are not entries in the reaction table. Making deuterium a network species would change the pp-chain model (p + p → d would no longer be instantaneous).
- **Surface Li has a numerical noise floor of about 0.005 dex.** The start-of-step CZ rate average uses the integer zone index of the CZ base (only the end of the step is interpolated within the zone), so a roundoff-level change in the structure can move it by one zone and shift surface Li7 by ~1% (ROADMAP §13).
- **Legacy cross sections.** `LNEWNUC = .false.` takes the cross sections from the legacy `SSTANDARD` ratios. Every shipped deck uses `.true.`.

# net/ — nuclear burning

This directory computes, for one zone (log ρ, log T, composition):
- the hydrogen- and helium-burning reaction rates and their derivatives, screened;
- the composition change over a step, by an implicit network built from a reaction table;
- the nuclear energy generation ε with its absolute derivatives for the Henyey Jacobian, by chain, with the thermal neutrino losses and the solar neutrino fluxes;
- deuterium and Li6, Li7, Be9 in the trace approximation, with the closed-form integrator they share.

net/ never touches `star`: the model comes in as arguments, and the run's configuration as one record, `net_settings` (`net_configuration`), which the star layer builds from the controls (`net_settings_from_star`, `setup/net_settings_setup.f90`) and installs with `set_net_settings` at the start of every Monte-Carlo run. `tools/check_boundaries.py` enforces this (`STAR_FREE_DIRS`).

Code outside net/ uses only `net_lib`, its public interface: a module with no code of its own that re-exports what the rest of YREC uses. Read `net/net_lib.f90` to see net/'s whole interface in one page; `check_boundaries.py` rejects a `use` of any other net/ module from outside (`FACADE_ONLY_DOMAINS`) and lists the procedures others may call (`PUBLIC["net"]`). The star-layer driver is `core/step/burn_step.f90` (§1.2); the Henyey solve calls `compute_nuclear_energy` itself, from the shell kernel `core/structure/shell_equations.f90`.

---

## 1. What lives where

### 1.1 net/

| Module (file) | What it holds |
|---|---|
| `net_lib` | The facade (above). |
| `net_configuration` | The layout and the settings: the reaction slots `r_*` (1 to `num_reactions`, 11), the two branching rows of the composition rate vector after them (`rv_frac_n15_pa`, `rv_frac_be7_ecap`; `num_rate_rows`), the extra cross-section scales `xs_*` (pep, Be7 + e, Be7 + p, hep: the rates outside the table; `num_extra_scales`), the settings record `net_settings_t` and the run's `net_settings` (protected), and `cross_section_scales_from_sstandard`, the legacy `sstandard` vector onto the scales. |
| `network` | The two tables and the code that reads them: the reaction table `reactions(:)` and the species table `species(:)` (`s_*` order); `compute_reaction_rates` (the one rate evaluation), `compute_network_rates` (the composition rate vector), `burn_implicitly` (the implicit network, with its codes `burn_not_converged`, `burn_singular_jacobian`), `compute_abundance_rates`, `compute_reaction_energies` (ε by slot and by chain from the table), `reaction_counts`, `reference_s0` (a reaction's S0 for the input scales), `s_factor_slots`; the chains `chain_*` with `chain_name`. |
| `network_checks` | The table's self-check (`network_table_errors`, with `table_row_problem` and `slot_problem`) and printout (`write_network_table`); `set_net_settings`, which runs the check and installs the settings, and `inert_scale_count`. |
| `reaction_rates` | The rate physics the table rows name: `compute_zone_conditions` (density, electron density, T9 powers, the screening plasma), `compute_screened_rate` (one row's rate and derivatives from its rate law, screened), the rate laws (`compute_s_factor_rate` and the Caughlan & Fowler `compute_c13_alpha_n_rate`, `compute_c12_alpha_gamma_rate`, `compute_n14_alpha_gamma_rate`, `compute_triple_alpha_rate`), `compute_be7_branching`, `n15_alpha_fraction`, `hep_rate_factor`, the d(p,γ) fit `compute_deuterium_proton_capture` and its overturn cap `cap_deuterium_rate`; the records `zone_conditions`, `nuclear_rates`, `be7_branching`, `t9_powers`; the rate-law names (`s_factor_fit`, `cf88_*`, `shares_rate`; `num_rate_laws`). |
| `screening` | The electron screening: the isotopes `iso_*` of the screening composition with their bare masses and charges, `screening_mass_fractions`, the plasma of a zone (`compute_composition_moments`, `compute_screening_plasma`), and one reaction's exponent with its derivatives (`compute_reaction_screening`: weak, intermediate or strong; `weak_intermediate_screening` for Be7 + p and hep), from its charges `screening_t` and the Graboske constants they give, `graboske_z`. |
| `nuclear_energy` | `compute_nuclear_energy`, the energy generation of one zone, as a `nuclear_energy_result` (§2.2). |
| `light_elements` | `deuterium_rate`, `burn_deuterium_region` (with `surface_accretion`), `compute_li_be_rates`, and the trace-element integrator `burn_exposure`, `surviving_fraction`, `accreted_survival` (`linear_in_time`, `linear_in_log`). |
| `thermal_neutrinos` | The Itoh et al. (1996) thermal neutrino losses as coded by F. X. Timmes: YREC's wrapper `compute_thermal_neutrino_losses` (G. Newsham 2006) and Timmes' `sneut`, `azbar`, `ifermi12`, `zfermim12` under their own names. Only `nuclear_energy` uses it. |
| `net_messages` | `net_report`: net's one kind of message (an error or a configuration warning), to the terminal and to the run log once it is open (§3). |

`retired/f77_rate_comments.txt` (not compiled) keeps the F77 comments of `rates.f` and `engeb.f` verbatim; each of `reaction_rates.f90` and `nuclear_energy.f90` ends with a Sources section that replaces them.

### 1.2 Outside net/

| File | What it does for net/ |
|---|---|
| `setup/net_settings_setup.f90` | `net_settings_from_star`: the settings from the controls and the star's cross-section scales (pure, `intent(in)` star). |
| `setup/user_inputs.f90` | `map_user_inputs`, during the controls read: the S0, S′ and S″ inputs (`s0_*`, `s0p_*`, `s0pp_*`) into `star%cross_section_scale`, `qs0e_scale`, `qqs0ee_scale` (by reaction slot), `extra_cross_section_scale` (by `xs_*`), `be7p_qs0e_scale`, `be7p_qqs0ee_scale` and the N15 + p branch scales. |
| `io/controls_reader.f90` | With `use_new_nuclear_rates` off, the legacy `sstandard` through `cross_section_scales_from_sstandard`. |
| `core/run/monte_carlo.f90` | A Monte-Carlo run rescales pp, He3 + He3, He3 + He4 (`r_*`) and Be7 + p (`xs_be7p`). |
| `core/run/run_driver.f90` | `run_yrec` installs the settings: `set_net_settings(net_settings_from_star(s), ierr)`, once per Monte-Carlo run. |
| `core/step/burn_step.f90` | The star-layer driver (table below). |
| `core/structure/shell_equations.f90` | `form_nuclear_energy_terms`: `compute_nuclear_energy` at every shell and Newton iteration of the Henyey solve. |
| `core/structure/structure_equations.f90` | The chains into `star%luminosity_breakdown_lsun` (`luminosity_channel_of_chain`) and `star%eps_channels` (`eps_channel_of_chain`); `eps_total`, the sum of every chain's ε. |

`core/step/burn_step.f90`'s procedures:

| Routine | Called from | What it does |
|---|---|---|
| `network_rates_of_model` | `compute_burn_rates`, `mid_timestep_model` | The composition rate vector of every shell of the model passed in that is hotter than the nuclear cutoff; zero below. |
| `compute_burn_rates` | `advance_composition` | Each shell's rate vector and deuterium rate (`deuterium_rate`), stored in `rate_by_zone`, `star%deuterium_burning_rate` (and `_start` in the predictor pass). |
| `burn_regions` | `burn_network`, `burn_deuterium`, `mid_timestep_model` | The regions a burn visits, in order: each radiative shell above the nuclear cutoff, then each mixed region as a unit. |
| `burn_network`, `burn_network_region` | `advance_composition`; `mid_timestep_model`'s `burn_and_settle` | Model bookkeeping around `burn_implicitly`: convection-zone averaging, the temperature cutoff, the first guess, the metallicity rule, writing `star%xa` back. |
| `burn_deuterium`, `surface_accretion_now` | `advance_composition`; `mid_timestep_model` | `burn_deuterium_region` for each region, with the step's accreted gas. |
| `compute_zone_energy` | `timestep_limit_heburn`, `write_neutrino_flux_table` | `compute_nuclear_energy` with the shell's place in the model filled in (in the surface convection zone? its turnover time); stores nothing. |
| `store_start_li_be_rates`, `store_end_li_be_rates`, `burn_li_be`, `average_li_be_over_cz` | `evolve_step`, `evolve_angular_momentum` | Li/Be rates per shell (`compute_li_be_rates`) and their depletion over a step with the trace-element integrator: radiative shells, then the surface convection zone on its mass-averaged rates, then the shells a retreating zone left. |

## 2. The flow

### 2.1 A burn (the composition step)

1. `compute_burn_rates` → `network_rates_of_model` → per shell `compute_network_rates`, which calls `compute_reaction_rates`:
   1. below `logT_cutoffs(1)` everything is zero;
   2. `compute_zone_conditions`: the density, the electron density, T9 and its powers, the screening plasma (composition moments, Salpeter's degeneracy factor, Λ₀, z̃, the mean charge), from `screening_mass_fractions(xa)`;
   3. each table row that computes a rate, when it burns (hydrogen-burning rows need hydrogen, helium-burning rows log T ≥ `logT_cutoffs(4)`): `compute_screened_rate` → its screening and its rate law, scaled at the row's slot;
   4. with hydrogen, the Be7 branching and the N15(p,α) fraction.

   `compute_network_rates` then drops rates at or below `min_composition_rate` (10⁻⁵) and returns the rate vector per s per amu: the reactions at `r_*`, the two fractions at `rv_*`.
2. `burn_network` → per region `burn_network_region` → `burn_implicitly`: Newton iteration on X − X₀ − dt A Σ ν R(X) = 0, its equations and Jacobian from the table rows (`compute_abundance_rates`, `reaction_rate_derivative`); the explicit products (O18) from the converged abundances. `ierr` is `burn_not_converged` when it runs out of iterations and `burn_singular_jacobian` when the Newton system is singular; `burn_network_region` writes a message for each (`report_error`: terminal and run log) and the run stops.
3. `burn_deuterium` → per region `burn_deuterium_region`: exponential decay with the rate linear in time between its start- and end-of-step values (`surviving_fraction`), D → He3, and the accreted gas (`accreted_survival`).

### 2.2 `compute_nuclear_energy`

1. Below `logT_cutoffs(1)` it returns `burning = .false.` and zeros.
2. `compute_reaction_rates` (as above).
3. `energy_path_rates`: the rates of the S-factor law (the slots `s_factor_slots` names, from the rows' `rate_law`) below `min_energy_path_rate` (10⁻³⁰) dropped, with their derivatives.
4. `compute_deuterium_burning`: the d(p,γ) energy and its d ln ε/d ln ρ, d ln ε/d ln T; capped in the surface convection zone (`cap_deuterium_rate`), where both derivatives are zero.
5. `reaction_counts` (for the fluxes) and `compute_reaction_energies`: Σ b_j R_j ΠX Q_j by slot and by chain, from the table, so a reaction's composition change and its energy come from one row.
6. `sum_energy_generation`: ε and its absolute derivatives Σ ε_j d ln r_j/d ln (ρ, T), plus deuterium's; zero derivatives at or below `min_eps_with_derivatives` (10⁻¹² erg/g/s).
7. `store_chain_energies`: `eps_by_chain(chain_*)`, d(p,γ) in pp I.
8. With `calc_neutrinos`, `compute_solar_neutrino_fluxes` (with the hep rate factor, `hep_rate_factor`, and the equilibrium Be7 abundance).
9. Above `logT_cutoffs(5)`, `add_thermal_neutrino_loss`: Beaudet, Petrosian & Salpeter (1967), or Itoh et al. (1996) with `use_itoh_neutrino_loss`, subtracted from ε with its derivatives.

### 2.3 When each piece runs during a timestep

```
evolve_step                                       core/step/step_driver.f90
 advance_composition_and_age
   advance_composition(predictor_pass)            core/step/composition_step.f90
     compute_burn_rates      [burn_step]          per shell: compute_network_rates, deuterium_rate
     burn_network            [burn_step]          per region: burn_network_region -> burn_implicitly
     (semiconvection), homogenize mixed regions, settling, renormalize
     burn_deuterium          [burn_step]          per region: burn_deuterium_region
 rezone_or_snapshot
   store_start_li_be_rates   [burn_step]          with the extended composition
 solve_structure, levels 1-3: solve_level
   refresh_surface_bc
   (levels >= 3) star%xa := star%xa_start; advance_composition(corrector_pass)
   henyey_iterate, the Newton loop
     compute_henyey_coefficients, every shell, every iteration:
       evaluate_shell_equations -> form_nuclear_energy_terms -> compute_nuclear_energy
                                                  eps and its derivatives in the luminosity equation
 converge_with_rotation
   evolve_angular_momentum (rotating runs), per sub-step:
     mid_timestep_model: network_rates_of_model (midpoint model) -> burn_and_settle -> burn_network_region;
                         burn_deuterium_region (the sub-step's accreted mass)
     store_end_li_be_rates (midpoint model), burn_li_be
 burn_light_elements                              non-rotating, or no rotational mixing
   store_end_li_be_rates (converged model), burn_li_be
 compute_timestep
   timestep_limit_hburn: compute_network_rates + compute_abundance_rates at the H-shell midpoint
   timestep_limit_heburn: compute_zone_energy at the hottest zone (the triple-alpha chain)
 output
   write_neutrino_flux_table (compute_neutrino_fluxes): compute_zone_energy per shell
```

Two things often surprise newcomers:
- **Composition and structure are operator-split.** Burning runs at fixed structure in the composition step, then the structure is solved with the new composition; from iteration level 3 the composition step runs again (the corrector pass) with the updated structure.
- **There are two rate paths from one evaluation.** The composition path (`compute_network_rates`: per s per amu, no derivatives, rates at or below 10⁻⁵ dropped) and the energy path (`compute_nuclear_energy`: erg/g/s with derivatives, the S-factor rates below 10⁻³⁰ dropped) both call `compute_reaction_rates`, so a rate is defined in one place.

## 3. Conventions

- **Units** are CGS and log is log₁₀; T9 = T/10⁹ K. Every time is in seconds. Energy derivatives are **absolute**, d ε/d ln T and d ε/d ln ρ, as the luminosity equation and GYRE's eps_eps_T want.
- **Composition.** `xa(i, k)` holds mass fractions in the slots of `state/composition_layout.f90` (net/ takes them from there, never from `star_info_lib`); the network species `s_*` map onto them with `network_species_slot`. The light elements (`i_h2`…`i_be9`) come with the extended composition.
- **Rates.** `nuclear_rates%rate(r_*)` is reactions per gram per second with the abundance factors omitted. The rate laws turn N_A⟨σv⟩ into it with ρ N_A/(A₁ A₂): the S-factor law through its fit (q7 holds −ln(A₀A₁) with bare nuclear masses and −ln(1 + δ₀₁) for identical nuclei); the Caughlan & Fowler laws through `per_gram_*` with the reactants' neutral atomic masses, and ρ² N_A/(3! A³) for triple α. A row's `reactants` list the mass fractions that multiply the rate (repeat a species for a square or a cube). The composition rate vector is the rate times `amu_cgs`, per s per amu, so that with Y = X/A, dY/dt = Σ ν R (integer mass numbers convert X to Y).
- **Q values** (`q_mev`) are the energy deposited per reaction, net of the neutrinos it emits: never more than the atomic mass excesses release. The self-check enforces that, to the table's 1 keV rounding (`q_rounding_mev`), and `write_network_table` prints the neutrino energy each row implies.
- **Screening** (`screening`). Each row carries the charges of its two nuclei, `screening_t(Z₁, Z₂)`; the weak form takes Z₁Z₂, and the intermediate and strong forms Graboske et al.'s (1973, table 4) z_p = (Z₁+Z₂)^p − Z₁^p − Z₂^p for p = 1.86 (intermediate) and 5/3, 4/3, 2/3 (strong), which `graboske_z` holds for every pair of charges from 1 to `max_screening_charge` (14), computed at compile time. The exponent is weak (Salpeter) below `weak_screening_threshold`, intermediate up to `max_intermediate_weak_exponent` (2), and above that the smaller of strong and intermediate, strong always from `min_strong_weak_exponent` (5).
- **Scales.** Only the S-factor law reads the per-reaction scales (`cross_section_scale`, `s_prime_scale`, `s_double_prime_scale` at the row's slot); a scale set on another reaction has no effect, and `set_net_settings` warns about it.
- **Temperature cutoffs** (`logT_cutoffs`, the control `nuclear_logT_cutoffs`): (1) gates everything, (4) helium burning, (5) the thermal neutrino losses; (2) and (3) are read and unused. The solar neutrino fluxes depend on `calc_neutrinos` only.
- **Errors and messages.** The kernels print nothing and are `pure`; failures return `ierr` (positive): `burn_implicitly` (`burn_not_converged`, `burn_singular_jacobian`), `set_net_settings` (a wrong table). The configuration's messages go through `net_messages`' `net_report`: the table's errors, and the warnings for a scale without effect or a legacy `sstandard` entry that belongs to no reaction.

## 4. Tests

`test_net` (`make tests` in `src/`, then `pytest test_net.py`) prints values for byte-comparison with `net/test/expected_test_net.out`, and 20 checks that print `ok` or `FAIL`:
- the configuration: `set_net_settings` installs `net_settings_from_star`'s record; the table's self-check; the legacy `sstandard` map; the scale reported as without effect;
- the rates' derivatives against centred differences (2×10⁻⁴ in T; 2×10⁻³ in ρ, where the screening plasma is held fixed);
- `burn_implicitly` against backward Euler's closed form with only pp, then only C12(α,γ); its mass conservation with every rate on; `burn_not_converged` when it runs out of iterations and `burn_singular_jacobian` for a singular Newton system;
- `surviving_fraction` and `accreted_survival` against closed forms and an independent quadrature;
- `compute_nuclear_energy`'s derivatives at the solar centre, a CNO core and a helium core (5×10⁻³); the thermal neutrino losses alone, on a C/O mixture where nothing burns, with both fits: negative, and their derivatives (BPS 10⁻⁶, Itoh 10⁻⁵); deuterium burning outside and inside the surface convection zone with the cap binding, with its derivatives; nothing below the burning cutoff;
- strong screening at log T 7, log ρ 6: the regime, the derivatives, and the exponent against Graboske et al.'s formula;
- `deuterium_rate`'s cap, exactly N_A/(m_d τ).

It also pins the composition rate vector at four temperatures, the cross-section scales, `sneut` and `azbar`, the reaction table with each row's Q against the mass excesses, `compute_abundance_rates`, `burn_deuterium_region` on a shell, a region, an accreting surface region and a region below the threshold, `compute_li_be_rates`, and the solar neutrino fluxes.

Each check was confirmed to fail under a mutation of the code it covers (26 mutations, 2026-10-29: the table, a dead slot, a slot computed twice, the derivatives of each rate law and of the neutrino losses, the deuterium cap, the cutoff, the screening regime and constants, the burn's equation and error, the integrator, the legacy map, the no-effect report). A deliberate change to the output means regenerating `expected_test_net.out` in a commit that says why. The decks pin the rest: every solar and testsuite deck runs the network, the energy and the light elements.

`python3 src/tools/coverage.py --work DIR --dir net --test test_net --deck ...` measures the line coverage of net/ under `test_net` and the decks given (a -O0 --coverage build in DIR). With `test_net` and five decks (`cm_accretion`, `cm_monte_carlo`, the 0.3 M☉ pre-main-sequence deck, the 1 M☉ ZAHB deck and a 0.05 M☉ Allard-atmosphere deck) it is 90.6% of 1582 lines (2026-10-29): 96.7% outside Timmes' `thermal_neutrinos` (84.7%, the regimes YREC does not reach). The lines not run are the self-check's error reports (run only under the mutations), the branch where the intermediate form is the smaller one between weak exponents 2 and 5, an unknown rate law's zero, and `surviving_fraction`'s underflow.

## 5. Extending

**Change a rate fit or coefficient.** An S-factor reaction's coefficients are the `fit` of its row in `network`'s table. A Caughlan & Fowler form is its rate-law routine in `reaction_rates.f90`; change its `dlnrate_dlnt` next to it, and `test_net`'s derivative check tells you whether they agree. Either moves results: a measured commit.

**Add a reaction** (among the tracked species, without input scales):
1. `net/net_configuration.f90`: name its slot (`r_<name> = num_reactions + 1`) and raise `num_reactions`. The composition rate vector's `rv_*` rows, the star's scale arrays and every per-reaction array follow.
2. `net/network.f90`: one row in `reactions(:)` (and `n_rows` + 1), importing the slot: its name, slot, branch (`unbranched`), reactants, change list (species, ν), `emitted_mass` for an untracked particle, `q_mev` net of neutrinos and its `chain`, `rate_law`, `burns_in` and `screening` (§3). The rate evaluation, the burn, the energy and the reaction counts all read the row; `set_net_settings` refuses a row that breaks the self-check.
3. A rate of a new functional form: in `net/reaction_rates.f90`, a rate-law constant (`num_rate_laws` + 1), a pure routine like `compute_c12_alpha_gamma_rate` with its exact derivatives and the prefactor convention of §3, and its `case` in `compute_screened_rate`.
4. Regenerate the pins: `net/test/expected_test_net.out` (the table, the derivative errors per slot) and the controls dumps (`YREC_REGEN_CONTROLS=1 pytest test_controls.py`: the star's scale arrays gain a slot).

**Add an S-factor reaction with input S0, S′ and S″ scales:** steps 1, 2 and 4 above, with `rate_law = s_factor_fit` and its `fit` (q1…q8: the S-factor law, `compute_s_factor_rate`, and the Sources at the end of `reaction_rates.f90`); its scales are the settings' at its slot already. Then:
1. `net/net_lib.f90`: export its slot.
2. Its three controls (`s0_<name>`, `s0p_<name>`, `s0pp_<name>`), as `io/README.md` §7 ("Add a control") says: `io/controls_buffer.f90` (the default as the initializer), `defaults/controls_registry.tsv` (with `canonical`, `home` `ctrl`, `used` `run`), the generators.
3. `setup/user_inputs.f90`: the fit's reference S0 (`s0_<name>_ref`) and its Solar Fusion S′, S″ (`qs0e_<name>_a98`, `qqs0ee_<name>_a98`, with `s0_<name>_a98`), and the three lines of `set_cross_section_scales`, `set_first_derivative_scales` and `set_second_derivative_scales` for its slot.
`net_settings_from_star`, `star_info` and Monte Carlo need nothing.

**Add a chain:**
1. `net/network.f90`: `chain_<name>` (`n_chains` + 1) and its `chain_name`; the rows that count under it. `nuclear_energy_result%eps_by_chain` follows.
2. `net/net_lib.f90`: export it.
3. `state/star_info_lib.f90`: its luminosity channel `i_lum_<name>` (`n_lum_channels` + 1), and an `i_eps_*` slot if it gets a fraction column (`eps_channels` is sized 7).
4. `core/structure/structure_equations.f90`: its case in `luminosity_channel_of_chain` (and `eps_channel_of_chain`, else 0). `eps_total` sums every chain already.
5. `core/diagnostics/observables.f90`: add the channel to `renormalize_luminosity_breakdown`'s sum.
6. For output, a history column that reads `s%luminosity_breakdown_lsun(i_lum_<name>)` (`io/history_output.f90`, `defaults/history_columns.list`).

**Add an isotope** (say Ne20, solved implicitly). This is the expensive extension:
1. `state/composition_layout.f90`: a slot `i_<name>` (`n_species_extended` + 1). The species count is still a literal 15 in some forty declarations across the code (Phase 5 removes them), and `io/model_files.f90` writes the 15 slots (a model-format change).
2. `net/network.f90`: the species `s_<name>`, before `s_o18` if solved (which renumbers `s_o18` and changes the pivot order: a measured change), `n_solved`/`n_network`, and its row in `species(:)` (name, xa slot, mass number, atomic mass excess), from which `network_species_slot`, `mass_number` and the mass excesses follow; then the reactions that make or burn it (above).
3. `net/screening.f90`: `screening_mass_fractions` (Ne20 and Mg24 have screening isotopes already); `nuclear_energy`'s thermal-neutrino mixture if it matters there.
4. Outside net/: the starting models (`core/model/starting_model.f90`), element diffusion's species (`element_diffusion/species_table.f90`), rotation's mixed species (`rotation_configuration`'s `num_tracked_species`), the outputs, the tests' compositions.

**Add a light element in the trace approximation.** Follow deuterium: a rate function and a region integrator in `light_elements.f90` (pure, everything as arguments, the integrator from `burn_exposure`), driven per region from `core/step/burn_step.f90`.

**Add a control.** Add the field to `net_settings_t` (`net/net_configuration.f90`) with the control's default as its declared default, copy it in `net_settings_from_star` (`setup/net_settings_setup.f90`), and read `net_settings%<field>`; the control itself as `io/README.md` §7 ("Add a control") says.

**Add an energy source that is not nuclear** (tidal or magnetic heating). That belongs to the structure equations, not here: `core/structure/energy_sources.f90` holds a list of extra sources a program registers (`add_energy_source`; `examples/other_energy` is one).

## 6. Known limitations and open items

These are tracked in `docs/ROADMAP.md` §13.
- **Special-purpose operators.** Deuterium and Li/Be are not entries in the reaction table. Making deuterium a network species would change the pp-chain model (p + p → d would no longer be instantaneous).
- **Surface Li has a numerical noise floor of about 0.005 dex.** The start-of-step CZ rate average uses the integer zone index of the CZ base, so a roundoff-level change in the structure can move it by one zone and shift surface Li7 by ~1%.
- **Legacy cross sections.** `use_new_nuclear_rates = .false.` takes the cross sections from the legacy `sstandard` ratios; its entries for the helium-burning reactions and for O16 + α and C12 + C12 have no effect (a warning says so). Every shipped deck uses the new path.
- **The resonant terms of the Caughlan & Fowler rates** matter only above T9 ≈ 1, which YREC does not reach; no check exercises them.

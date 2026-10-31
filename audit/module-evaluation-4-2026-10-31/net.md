# net/ -- fourth evaluation (2026-10-31, yrec-next e8c847ba, read-only)

Scope: src/net/ (9 modules, README, test_net, retired/), setup/net_settings_setup.f90, state/reaction_layout.f90, and the files README §1.2 lists (the 14 outside net/ that use net_lib or reaction_layout, and the others that hold a piece of net's layout), read against src/rotation/README.md, GUIDELINES "Naming" / "The star object" / "Declarations and interfaces" / "Code style" and TARGET_ARCHITECTURE §4.5(b). Nothing under /Applications/yrec was changed. The scratch builds and scripts are in SP/review4_net/: copies of net/, numerics/, math/ and state/ with an io_units stub, the drivers drive_nan*.f90, and defre*.py.

## Scores (1-10)

Earlier scores (out of 5, doubled) -> now. The code has not changed since the third evaluation apart from the controls rename (`s%ctrl` -> `s%physics` in the builder and two comments). A score that drops here does so because this review found new items, not because the code got worse.

- **Readability: 8 -> 9 -> 9 -> 9.** The table-driven network, the named steps (`compute_nuclear_energy` reads as its own list of steps) and the one-line headers are still clear. Below 10: the Beaudet-Petrosian-Salpeter block (nuclear_energy.f90:333-417) is the one block outside Timmes' code that keeps opaque F77 names (`el, eli, ez, polx10..polx32, qedn, qetnx, ex1..ex3`) and has no equations; network.f90:42-50 keeps provenance in the module header; the argument header at nuclear_energy.f90:100-103 is garbled.
- **Organization: 8 -> 9 -> 9 -> 8.5.** Dependencies run one way, there is a single facade, and the builder is in setup/. Below 9: the Li/Be reference cross sections are in io/, away from the fits they normalize (finding 1); the two thermal-neutrino fits are split, BPS inline in nuclear_energy and Itoh in thermal_neutrinos; network.f90 still holds the tables, the rate evaluation and the Newton solver; net_lib has two dead exports.
- **Maintainability: 6 -> 8 -> 9 -> 8.5.** test_net's 25 checks are each confirmed by a mutation, set_net_settings runs a self-check, and the builder is checked component by component. Weighing against: four open defects that change results (M89-M91, plus finding 1, which went unnoticed since 2014 because the reference and its fit live in different directories); a hole in burn_implicitly's ierr contract (finding 2); a boundary check that sees 2 of net's 9 public entries (finding 4); and M102, M114 and M118.
- **Extensibility: 6 -> 7 -> 8 -> 7.5.** Adding a reaction is still one row plus a slot: 2 files, a 3rd for a new rate law, then the pins. Below 8: following the chain recipe leaves a new hydrogen-burning chain out of the timestep limiter's L_H (finding 3); following the isotope recipe can count a tracer twice in the Itoh mixture (finding 6); `n_rows`, `n_network` and `n_chains` are hand-kept literals; and relative to §4.5(b), scales are addressed by slot rather than name, each rate law takes the screening (it is not a bare fit), and there is no way to register a law.
- **Docs accuracy: 6 -> 7 -> 8 -> 7.5.** README §1.1 to §4 check out claim by claim (users, flow, conventions, the 25 checks), and so do the recipes' listed steps. Newly false or stale:
  - §3's "one edit" claim, for Li/Be;
  - §1 claims an enforcement (`PUBLIC["net"]`) that does not happen;
  - §6 points to the wrong ROADMAP section and never mentions bug_sweep §9;
  - the chain and isotope recipes each miss a file;
  - the slot location is stale in three module headers and in state/README;
  - GLOSSARY has a rename slip.

  Three known traps (M118, M121, M90's comment) are still in place.

## Summary

- The README's map holds. All 14 files outside net/ that use net_lib or reaction_layout are listed, and every call in §1.2 and §2.3 is where it says. §3's conventions match the code: units, the Avogadro convention, the Q-value check, the screening regimes and the cutoffs, except M118. §4's 25 ok/FAIL checks are all present.
- **Recipe trace, counting hand-edited files.**
  - Add a reaction: 2 files (reaction_layout, network), a 3rd for a new rate law (reaction_rates), then the pins (expected_test_net.out and the controls dumps).
  - S-factor reaction with inputs: 6 files (reaction_layout, network, net_lib, controls_buffer, controls_registry.tsv, user_inputs) plus 3 generators and the pins.
  - Chain: the 4 files listed (network, net_lib, star_info_lib, structure_equations) plus io for the history column. A 5th is needed and not listed: timestep_control (finding 3).
  - Isotope: 3 files on net's side (composition_layout, network, nuclear_energy). nuclear_energy needs two lists, and the recipe names one. Then Phase 5's literal-15 declarations (~29 lines in ~28 files).
- **New physics item (inherited from G. Somers 6/14).** The Li6 rate scale is taken against 5.5 MeV b (Fowler et al. 1967), but the CF88 fit that net computes has S(0) ≈ 2.69 MeV b. The 1160 repository inlists that set `scale_li6_rate = .true.` and `li6_rate_value = 3100` (keV b) therefore burn Li6 at 0.49× the rate they ask for.
- **New latent contract hole.** `burn_implicitly` returns ierr 0 with all-NaN abundances when only the Newton right-hand side is NaN. This was confirmed on a scratch build at -O0 and -O2; at -O3 it reports "singular" instead, which is M102's mislabel.
- **New tooling/docs gap.** check_boundaries' `DEF_RE` cannot see `pure` procedures, so `PUBLIC["net"]` enforces 2 of its 9 entries (README:11 claims more). Only the facade check protects net today.
- **Known open items** (bug_sweep §9: M89-M91 change results, plus M102, M114, M118, M121) and the new ones above keep net below 9 on four of the five criteria. Readability stays at 9.

## New findings (not in bug_sweep), ranked

1. **io/legacy_controls.f90:222-227, io/derived_controls.f90:74-85 against net/light_elements.f90:321-368.**
   - Issue: the Li/Be override references are not the S(0) of the fits that net computes.
     - `derive_rate_scales` sets `li6_rate_scale = xsli6/5.5e3` ("Li6 5500 (Fowler et al. 1967)"). `compute_li_be_rates` uses Caughlan & Fowler 1988, whose Li6(p,α) prefactor 3.73e10 gives S(0) = 2689 keV b through N_A⟨σv⟩ = 7.8324e9 (Z0Z1/Â)^{1/3} S/(Â T9^{2/3}) e^{-τ}.
     - The conversion was checked two ways. It gives CF88 pp's S(0) as 4.06e-22 keV b, against net's 4.07e-22. The τ it predicts is 8.395, against the fit's 8.413.
     - In 1160 repository inlists (grid/, sample_cases/, startmodels/, examples/, testsuite/), plus 18 legacy .nml decks, `scale_li6_rate = .true.` with `li6_rate_value = 3100` ("Cross-section for p+6Li (keV barns)"). The applied scale is 0.564 where the intended one is 1.15, so Li6 burns at 0.49× the requested rate.
     - Li7 (0.87×) and Be9(p,γ) (0.85×) are also off, but within those fits' energy dependence. Be9(p,d) and Be9(p,α) agree to 2%.
     - The F77 has the same code (parmin.f:800-808 with lirate88.f:46-52).
     - README §3 says of the references that "each is defined once, so changing a fit's normalization and its reference is one edit". For Li/Be that is false: the references are in io/, and README §1.2 does not list those files.
   - Fix: put each Li/Be fit's S(0) beside its fit in light_elements, exported through net_lib as `reference_s0` is, and have `derive_rate_scales` scale against those values. Re-check Li7 and Be9(p,γ) against their fits. Measure the change in surface Li6 on the pre-MS and solar decks.
   - **B** (Li6 depletion in every extended-composition deck that sets the override) -- S (code) / M (measurement).
2. **net/network.f90:584-604.**
   - Issue: `burn_implicitly` can report a success it did not have.
     - `simeqc` tests only for NaN pivots, so a NaN confined to the right-hand side gets through; that happens when the solved species of `x_start` hold a NaN.
     - The corrections then come out NaN. `max(max_abs_change, abs(correction))` with a NaN argument is processor-dependent, and the convergence test can pass on the first iteration.
     - Confirmed on scratch copies (SP/review4_net/src/drive_nan.f90 and drive_nan2.f90, gfortran 15.2): with a NaN in x_start(He3) or x_start(C13), the routine returns ierr 0 and all-NaN x at -O0 and -O2. The coverage and debug builds use -O0. At -O3 the same input returns ierr 2, which is M102's "singular" message. A NaN in a rate gives ierr 2 at every level.
     - burn_network_region then writes the NaNs into `star%xa` with no message, against GUIDELINES code style 7.
   - Fix: test `all(ieee_is_finite(correction))` (and the residual) before the convergence test, and return a named code of its own. M102 proposes one for the NaN pivot, so one new code covers both. Give it a case and a check in test_net.
   - **A** -- S.
3. **core/step/timestep_control.f90:70, :170-173, :186-188; core/structure/henyey_iteration.f90:363-366.** README:182-187 ("Add a chain") misses these.
   - Issue: the step's hydrogen-burning luminosity is `luminosity_components(1)+...+(4)`, by literal index, and the helium limiter keys on component 5. The dummy is declared `(8)` while the actual has `n_lum_channels` = 9. The Henyey log line hand-sums pp1..CNO and 3α + He-C.
     - The recipe says that `eps_total` and both renormalizations follow a new chain. But a new hydrogen-burning chain is left out of L_H: the trace's own hep, or CNO split into CN and NO. The shell-burning limiter then allows longer steps.
   - Fix: give each chain a burning group in network's chain data (hydrogen or helium), and add one star-layer sum over the channels of a group, used by both the timestep and the log line. Use `i_lum_*` names, size the dummy with `n_lum_channels`, and add the step to the recipe.
   - **A** (today's chains); B for a run that follows the recipe -- S.
4. **src/tools/check_boundaries.py:116-121 (DEF_RE), :86-90 (`PUBLIC["net"]`).** README:11 relies on these.
   - Issue: `DEF_RE` accepts only `[recursive] subroutine|function|double precision/integer/logical function`. It sees 5 of net's 63 procedures and 2 of the 9 `PUBLIC["net"]` entries. The other seven, `compute_network_rates`, `burn_implicitly`, `compute_nuclear_energy` and the rest, are `pure`, so a cross-domain `call` of any pure net procedure passes; atm has 3 of 8 visible. Checked with a copy of the regex (SP/review4_net/defre_all.py). Only the facade check protects net today.
   - Fix: accept the prefixes `pure|elemental|impure|recursive|module` and a type-spec, as M116 proposes for the star tools, then re-run the check.
   - **A** -- S.
5. **Stale slot location in four places.**
   - net/network.f90:21 says "slots r_* (net_configuration)".
   - net/network.f90:64-66 is the module's own recipe: "naming its rate slot r_* in net_configuration".
   - net/reaction_rates.f90:24-25 says "the reaction slots ... are net_configuration's".
   - state/reaction_layout.f90:5 says "net_configuration re-exports them" (it is private and exports only its own names), and :9 says "setup's map_user_inputs" (it has been io's since 9cef61be). state/README.md:10 repeats the re-export claim.
   - All of these predate ec2eaaf6's move to reaction_layout. A newcomer who reads the module header instead of the README edits the wrong file.
   - Fix: point each to state/reaction_layout.f90.
   - **A** -- S.
6. **README:191 (isotope recipe) against net/nuclear_energy.f90:92-94, :363.**
   - Issue: the recipe offers to add the species to `itoh_slots` (with its nucleon number and charge) "if it should count ... on its own rather than as part of Z". It does not say that a metal tracer must also go into `metal_tracer_slots`, the list subtracted from Z for the Ne20 remainder; otherwise the tracer is counted twice in the Itoh mean A and Z. The list also restates composition_layout's "slots 5..11 and 13..15 are tracers inside Z" by hand.
   - Fix: name the tracer set once in composition_layout and use it in nuclear_energy; add the line to the recipe.
   - **A** -- S.
7. **README:202 (§6).**
   - Issue: it says "The open items are listed in docs/ROADMAP.md's last section ("Still open")". The last section is now "The controls follow the inlist groups" (ROADMAP:3241); net's still-open list is item 3 of round 2's "Still open" (ROADMAP:3086-3110). The open defects M89-M91, M102, M114, M118 and M121 (bug_sweep §9) are not pointed to from net/ at all.
   - Fix: point to the named subsection and to bug_sweep §9 by M number, as rotation's docs name the O items.
   - **A** -- S.
8. **net/network_checks.f90:132-176.** The self-check has two gaps.
   - A row's `screening_t` is range-checked only, never compared with its reactants' charges, although the species table already gives them (`screening_isotope` -> `isotope_charge`). O16 + p typed (1,7) would pass. Triple α's (2,6) is the one convention.
   - A `shares_rate` row's reactants are never compared with its owner's. `reaction_counts` (network.f90:345-360) uses the owner's reactants, while the energies (:475-498) and the burn use each row's own, so a mismatch makes the neutrino fluxes and the energies disagree silently.
   - Fix: derive the charges from the first two reactants unless the row overrides them, and check that a sharing row's reactants equal its owner's.
   - **A** -- S.
9. **core/diagnostics/neutrino_flux_output.f90:20-37, :83-96, with net/nuclear_energy.f90:137-138.**
   - Issue: with `compute_neutrino_fluxes = .true.` and `calc_neutrinos = .false.`, the per-shell flux table and its totals are all zeros, because net computes fluxes only with `calc_neutrinos`. Nothing warns: neither check_controls nor warn_inconsistent_controls has a test, and the registry row (controls_registry.tsv:388) does not mention the dependence.
   - Fix: warn (or require `calc_neutrinos`) at the read, and say it in the registry and in the routine's header.
   - **A** -- S.
10. **Small drifts.**
    - net_lib.f90:41-43 still re-exports `num_reactions` and `num_extra_scales`. Nothing outside net/ takes them from net_lib (star_info_lib uses reaction_layout), against the header's "and nothing more" (:7-8).
    - README:157 says "a new species or reaction changes the pins, not the program". The rate-grid pin prints 11 named slots (test_net.f90:221-226), so a new reaction's composition rate is not pinned until the program changes.
    - docs/GLOSSARY.md:255 says "n15_pg_scale and n15_pa_scale until 2026-10". The old names were `o16_gamma_scale` and `c12_alpha_scale`; 6882c0eb's rename caught its own history note.
    - nuclear_energy.f90:100-103 runs the `in_surface_cz` and `turnover_timescale` descriptions together.
    - **A** -- S.

## What would raise each score

- **Readability (9 -> 10):**
  - Rewrite `add_thermal_neutrino_loss`'s BPS part with named terms (pair, photo, plasma), the fit's equations in the header and readable locals.
  - Move network.f90:42-50's provenance into a Sources section at the end of the file.
  - Fix compute_nuclear_energy's argument header.
  - Break README §1.1's longest rows (network, reaction_rates) into short lists.
- **Organization (8.5 -> 9.5):**
  - Put the Li/Be reference S(0) values beside their fits (finding 1); io keeps only the inputs.
  - Move the BPS fit into thermal_neutrinos beside Itoh, as one module with two fits, each returning loss and absolute derivatives.
  - Move `burn_implicitly` and `compute_abundance_rates` into a module of their own, so network holds the tables and the rate evaluation.
  - Drop net_lib's dead exports.
  - Route `net_report` through io_units' report_* (M114), so net has one message path.
- **Maintainability (8.5 -> 9.5):**
  - Fix M89, M90, M91 and finding 1, each as a measured commit.
  - Add the NaN guard and its named code (finding 2, with M102).
  - Make `DEF_RE` see pure procedures (finding 4).
  - Close the self-check gaps (finding 8).
  - Raise network_checks' 82% coverage with a test that feeds `table_row_problem` and `slot_problem` a variant row.
  - Add a test that each light-element reference equals its fit's S(0), so the drift in finding 1 cannot recur.
- **Extensibility (7.5 -> 9):**
  - Use implied-shape tables: `reactions(*)`, `species(*)`, `chain_name(*)`, with `n_rows`, `n_network` and `n_chains` = `size(...)`. gfortran 15.2 accepts this; checked in SP/review4_net/src/implied.f90. That removes a hand step from three recipes and makes "the species table is the single source of the species count" (§4.5(b)) literal.
  - Add a burning group per chain, so L_H and the log line follow a new chain (finding 3).
  - Name the metal-tracer set once (finding 6).
  - Derive screening charges from the reactants.
  - Address rate scales by reaction name for the inputs and Monte Carlo (§4.5(b)).
  - Do Phase 5's literal-15 declarations, which isotopes need.
  - Either give rate laws a registration path (a procedure pointer), or record in §4.5(b) that integer dispatch is the decision; Phase 1's status says so, §4.5(b) does not.
- **Docs accuracy (7.5 -> 9):**
  - Fix findings 3, 5, 6, 7 and 10.
  - Correct README §3's "one edit" claim (finding 1) and README:11's `PUBLIC` claim (finding 4).
  - Fix the known traps M118, M121 and M90's "the screening's stand-in" comment.
  - Re-trace the chain and isotope recipes through the timestep limiter and the Itoh mixture.
  - Add a Status line to TARGET_ARCHITECTURE §4.5(b) saying what of the type/procedure-pointer design was not taken: rate laws take the screening, scales are by slot, light elements stay special.

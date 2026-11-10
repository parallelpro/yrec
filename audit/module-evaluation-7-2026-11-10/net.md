# net -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

**Scope.**
- **Read whole:** every hand-written .f90 of src/net/ except Timmes' `sneut` body (thermal_neutrinos.f90: the header, the wrapper and `azbar` read, the rest skimmed for I/O and F77 constructs); net/README.md; net/test/test_net.f90; test_net.py; the header of net/retired/f77_rate_comments.txt; setup/net_settings_builder.f90; state/reaction_layout.f90; core/step/burn_step.f90.
- **Read for net's claims:** the net parts of setup/settings_install.f90, io/controls_buffer.f90, io/user_inputs.f90, io/derived_controls.f90, defaults/controls_registry.tsv, io/test/controls_dump.inc, the Makefile's controls rules, state/run_messages.f90, eos_messages and rotation_messages, star_info_lib, timestep_control, structure_equations, check_boundaries.py (PUBLIC, DEF_RE, STAR_FREE_DIRS, FACADE_ONLY_DOMAINS, ONLY_LIST_MODULES), ci.yml, docs/GUIDELINES.md (Naming, the star object, Declarations, Code style), TARGET_ARCHITECTURE §4.5(b) and Phase 1, bug_sweep §9-§13 (every net, burn and light-element item), and ROADMAP from round 2's "Still open" to the end.
- **What changed since the fourth evaluation (e8c847ba).** `git log e8c847ba..823326a7 -- src/net` lists six commits, none of them net's own. They brought:
  - the builder's rename (net_settings_setup to net_settings_builder);
  - the import of `run_log_unit` from run_messages;
  - test_net's switch to the controls' run copies;
  - M161's header fixes (a43d1efa).
  - net/ code is otherwise byte-for-byte as the fourth evaluation read it.
- **Counts.** Scratch Python over the tree (`python3 -I`), in SP/review7_net/:
  - `netlib_uses.py`: which net_lib exports anything outside net/ takes;
  - `defre.py`: a copy of check_boundaries' `DEF_RE` against every procedure header in net/.
- **Method.** The other counts are greps quoted with their results. One recipe ("Add a reaction") was traced statically through make's generators. Nothing was built, run or written in MAIN.

## Scores (1-10)

Earlier scores (out of 5, doubled) -> fourth evaluation -> now. Nothing in net/ changed, so a score that moves does so because the code around net/ moved under its documents, or because this review found something the earlier ones did not.

- **Readability 8 -> 9 -> 9 -> 9 -> 9.** It still reads as a list of named steps over two tables, every kernel `pure`, with no commented-out code. Left:
  - the Beaudet-Petrosian-Salpeter block (nuclear_energy.f90:333-417) keeps opaque F77 names (`el, eli, ez, polx10..polx32, qedn, qetnx, ex1..ex3`), no equations, and inline polynomial and derivative literals (:381, :403);
  - provenance in module headers: network.f90:42-50, reaction_rates.f90:2, light_elements.f90:5-7, reaction_rates.f90:447-448;
  - compute_nuclear_energy's garbled argument header (nuclear_energy.f90:100-103);
  - the intermediate screening formula written twice, with bare Graboske coefficients (finding 4);
  - "ABSOLUTE" in capitals (nuclear_energy.f90:358);
  - three README rows or paragraphs over 1,000 characters (README:25, :159, :161; rotation's README has none).
- **Organization 8 -> 9 -> 9 -> 8.5 -> 8.5.** Dependencies run one way, there is one facade, the builder is in setup/ and net/ is in `STAR_FREE_DIRS` and `FACADE_ONLY_DOMAINS`. Left:
  - the Li/Be reference cross sections live in io/ (controls_buffer.f90:46-48), away from the fits they normalize (M127);
  - the two thermal-neutrino fits are split (BPS inline in nuclear_energy, Itoh in thermal_neutrinos);
  - network.f90 holds the tables, the rate evaluation and the Newton solver;
  - net_lib still re-exports `num_reactions` and `num_extra_scales`, which nothing outside net/ takes from it (SP/review7_net/netlib_uses.py: 0 users);
  - net_report keeps a message path of its own beside run_messages' queue (M114);
  - net's driver, burn_step, keeps F77 shapes and literal thresholds (finding 3).
- **Maintainability 6 -> 8 -> 9 -> 8.5 -> 8.5.** test_net's 25 ok/FAIL checks are each confirmed by a mutation. Two gains:
  - CI now builds and runs test_net (ci.yml:131), though its first GitHub run is still pending;
  - the registry's `settings` column ties every net control to its `net_settings` field, and check_boundaries holds the builder to it (ecaf66ff, 79c54d0b).

  Against that, every open defect of net and its driver is still open:
  - results-changing: M89, M90, M91 and M127, plus M176 (burn_step's start-of-step deuterium rates, found since);
  - latent: M102, M118, M131 and M139;
  - docs and tools: M114, M121, M141, M146, M147 and M162.

  The test gaps are unchanged too (no NaN case, no check of the eps_nuc derivatives, light-element results pinned only).
- **Extensibility 6 -> 7 -> 8 -> 7.5 -> 7.5.** Adding a reaction is still one table row plus a slot. I traced it through make: reaction_layout.f90 is one of the controls generators' inputs, so the dump and the controls follow on their own. Left:
  - `n_rows`, `n_network`, `n_chains`, `num_itoh_slots` and `metal_tracer_slots(10)` are hand-kept literals;
  - a new hydrogen chain is left out of L_H (M146);
  - a tracer can be counted twice in the Itoh mixture (M147);
  - the light-element recipe names 2 of the dozen files deuterium touches (finding 2);
  - scales are by slot, rate laws take the screening, and there is no registration path (§4.5(b), unrecorded);
  - Phase 5's literal 15: 47 lines in 26 files outside net/.
- **Docs 6 -> 7 -> 8 -> 7.5 -> 7.** §1.1, §2, §3 and §4 still check out claim by claim, but for the cutoff wording (finding 5), and M161 is fixed. Lowered because round 3 moved the rate scales under net's documents and they did not follow:
  - README §1.2 has 3 stale rows and 1 missing file, and the reaction recipe's pin step names the old dump (finding 1);
  - the light-element recipe is untraced and thin (finding 2);
  - the cutoff control's description names no element (finding 5);
  - the fourth report's docs items are all still open except M161: §3's "one edit" (M127), README:11's PUBLIC claim (M141), §6's pointer (M162, now further off: the last ROADMAP section is the sixth-evaluation pass), the chain and isotope recipes (M146, M147), the traps M90 (comment), M118 and M121, and GLOSSARY:255.

**Average 8.1** (8.2 at the fourth evaluation).

## Summary

- **net/ has had no round since the fourth evaluation.** Of that report's 10 findings:
  - 1 is done: M161, the slot location;
  - 8 are open: M127, M131, M146, M141, M147, M162, M139, and the self-check gaps, which bug_sweep left out;
  - finding 10's small drifts are open, except README:157, which bug_sweep judged true.

  Of its ~30 "what would raise" items, only the slot-location fix and the builder's tooling moved.
- **The code is still the strongest physics domain after rotation.**
  - A table-driven network with a self-check.
  - One rate evaluation behind both paths.
  - Every kernel pure and star-free, failures as named ierr codes, messages through one module.
  - 25 mutation-confirmed checks, and CI now runs them.
- **Docs fell behind round 3.** 3a44ae33 (2026-10-31) turned six of the eight rate scales into `physics%` controls and the other two into run copies, and moved the legacy `sstandard` map from controls_reader to user_inputs. state/README was updated (its rows 41 and 129 are right); net/README §1.2, its reaction recipe, GLOSSARY:255 and HOW_YREC_WORKS:445 were not.
- **New this time:**
  - the light-element recipe names 2 files where following deuterium touches about 12;
  - the intermediate screening formula is written twice;
  - burn_step's light-element bookkeeping has 11 literal thresholds and 23 F77 comparisons;
  - the user-facing description of `nuclear_logT_cutoffs` says nothing of what each of its five elements gates.

  No new bug: I re-derived the S-factor law's derivatives and scales, the BPS ρ and T derivatives, the Be7 and flux units and the deuterium bookkeeping, and found them right.

## Previous report's items: status

**New findings (fourth evaluation):**
1. Li/Be references are not the fits' S(0), and README §3's "one edit" claim (M127). **OPEN.**
   - The references are now named parameters in io/controls_buffer.f90:46-48 (`li6_reference_xs = 5.5d3`, ...), used by derive_rate_scales (io/derived_controls.f90:109-115).
   - The fits are unchanged (light_elements.f90:321-368), and README:141 is unchanged.
2. `burn_implicitly` can report success with all-NaN abundances (M131). **OPEN.**
   - network.f90:583-604 is unchanged; test_net has no NaN case.
3. The timestep limiter's L_H sums channels 1-4 by number (M146). **OPEN.**
   - timestep_control.f90:70 (`luminosity_components(8)` against `n_lum_channels` = 9), :172-173, :188; README:182-187 is unchanged.
4. `DEF_RE` cannot see `pure` procedures (M141). **OPEN.**
   - DEF_RE is now at check_boundaries.py:110-115, `PUBLIC["net"]` at :80-84.
   - Recount (SP/review7_net/defre.py): 5 of net's 63 procedures and 2 of the 9 PUBLIC entries visible, as before.
5. Stale slot location in four places (M161). **DONE** in a43d1efa.
   - network.f90:21 and :65, reaction_rates.f90:24-25, reaction_layout.f90:4-13 and state/README.md:41 all name state/reaction_layout correctly.
6. The isotope recipe omits `metal_tracer_slots` (M147). **OPEN.**
   - nuclear_energy.f90:83-94 and README:191 are unchanged.
7. README §6 points to the wrong ROADMAP section, and GLOSSARY's rename slip (M162). **OPEN.**
   - README:202 still says "the last section ("Still open")"; the last section is now "Sixth-evaluation pass" (ROADMAP:3755), with no such subsection. Net's open list is ROADMAP:3108-3113.
   - GLOSSARY:255 still says "`n15_pg_scale` and `n15_pa_scale` until 2026-10".
8. The self-check does not tie screening charges to the reactants, nor a `shares_rate` row's reactants to its owner's. **OPEN** (bug_sweep §10 left it out). network_checks.f90:132-176 is unchanged.
9. `compute_neutrino_fluxes` without `calc_neutrinos` gives a table of zeros (M139). **OPEN.** The control's description (controls_buffer.f90:440-443) still does not mention the dependence.
10. Small drifts.
    - net_lib's dead exports: **OPEN** (net_lib.f90:41-43).
    - README:157 "pins, not the program": **dropped**; bug_sweep §10 judged it true.
    - GLOSSARY: **OPEN**, under M162.
    - nuclear_energy.f90:100-103's argument header: **OPEN**.

**What would raise readability:**
- the BPS rewrite with named terms: **OPEN**;
- move network.f90:42-50's provenance to Sources: **OPEN**;
- fix compute_nuclear_energy's header: **OPEN**;
- shorten README §1.1's rows: **OPEN** (row 25 is 1,203 characters).

**Organization:**
- Li/Be S(0) beside the fits: **OPEN**;
- BPS into thermal_neutrinos: **OPEN**;
- `burn_implicitly` in a module of its own: **OPEN**;
- drop the dead exports: **OPEN**;
- net_report through run_messages (M114): **OPEN**, and now a smaller change. run_messages already queues a line written before the run log opens (log_line, flush_pending_log, run_messages.f90:111-137). net_report (net_messages.f90:22-35) still writes the terminal and an open log itself, so the controls read's `sstandard` warning never reaches the run log.

**Maintainability:**
- M89, M90, M91 and M127 as measured commits: **OPEN**;
- the NaN guard (M131 with M102): **OPEN**;
- `DEF_RE`: **OPEN**;
- the self-check gaps: **OPEN**;
- network_checks coverage: **OPEN**;
- a test that each Li/Be reference is its fit's S(0): **OPEN**.
- Not on the list, done since: CI runs test_net (ci.yml:131, f66688ca); the registry's `settings` column and check_boundaries' builder checks cover net's builder (ecaf66ff, 79c54d0b).

**Extensibility:**
- implied-shape tables: **OPEN** (network.f90:94, :144, :192);
- a burning group per chain: **OPEN**;
- one named tracer set: **OPEN**;
- screening charges from the reactants: **OPEN**;
- scales by name: **OPEN**;
- Phase 5's literal 15: **OPEN**;
- record the integer-dispatch decision in §4.5(b): **OPEN**.

**Docs:**
- findings 3, 6, 7 and 10: **OPEN**; finding 5: **DONE**;
- README §3's "one edit" and README:11: **OPEN**;
- the traps M118, M121 and M90's comment (nuclear_energy.f90:79-80): **OPEN**;
- re-trace the chain and isotope recipes: **OPEN**;
- a Status line in TARGET_ARCHITECTURE §4.5(b): **OPEN**.
  - §4.5(b) still says "Screening is applied by the network driver, using the charges in the species table, so rate functions stay bare fits".
  - It also says "Light elements become ordinary species".
  - Neither is the code. Phase 1's status records only the integer dispatch and the special light elements.

**ROADMAP "Still open" item 3 (net):**
- test_net has no finite-difference check of `deps_nuc_dlnt`/`deps_nuc_dlnrho`: **OPEN**;
- Phase 5: **OPEN**.

## New findings (not in bug_sweep)

No new bug. Docs and structure, in order of what a newcomer would act on.

1. **net/README.md:41, :43, :44, :173 (and a missing row); docs/GLOSSARY.md:255; docs/HOW_YREC_WORKS.md:445.**
   - **Issue:** net's documents still describe the rate scales as the star's members, written by the controls read through controls_reader. 3a44ae33 (2026-10-31, round 3) changed that, and a43d1efa updated state/README but not net/README. In detail:
     - **README:41 (star_info_lib).** "The star's scale arrays ... `cross_section_scale`, `s_prime_scale`, `s_double_prime_scale` ..., `extra_cross_section_scale`, `be7p_s_prime_scale`, `be7p_s_double_prime_scale`, `n15_pg_scale` and `n15_pa_scale` ..., under `net_settings_t`'s names."
       - star_info_lib now holds only the two run copies, `cross_section_scale_in_use` (:228) and `extra_cross_section_scale_in_use` (:233).
       - The other six are derived `physics` controls (io/controls_buffer.f90:660-673; registry rows 397-402).
       - The builder reads them as `s%physics%...` (setup/net_settings_builder.f90:42-48).
     - **README:43 (io/user_inputs).** It says map_user_inputs maps the S0 inputs "into the star's scale arrays". map_user_inputs takes no star since 3a44ae33. It writes the buffer's derived controls (user_inputs.f90:106-119); the read stores them on `star%physics`, and `arm_run_copies` copies two of them into the run copies (controls_sync.f90:911-912).
     - **README:44 (io/controls_reader).** "With `use_new_nuclear_rates` off, the legacy `sstandard` through `cross_section_scales_from_sstandard`": controls_reader no longer calls it. The call is user_inputs.f90:114, on `cross_section_ratios` (the control the legacy `sstandard` reads into, controls_buffer.f90:647).
     - **Missing row.** io/controls_buffer.f90 uses reaction_layout (:31) and declares every net control, with default and description:
       - the scale controls (:645-699);
       - the Li/Be references of M127 (:46-48);
       - `nuclear_logT_cutoffs` (:608) and the burn's tolerances (:632-644).

       A newcomer who asks where a rate scale's default lives is not sent there.
     - **README:173 ("Add a reaction", step 4).** It says "every dump gains the star's `cross_section_scale(<slot>)`". In fact:
       - the dump writes `physics%cross_section_scale`, `physics%s_prime_scale`, `physics%s_double_prime_scale` and `star%cross_section_scale_in_use` (io/test/controls_dump.inc:258-269, :608-610);
       - with M180's per-deck difference hashes only `pristine.txt` moves;
       - the step does not say that make regenerates io/test/controls_dump.inc: its `do i1 = 1, 11` loops are `num_reactions` written out by gen_controls_dump, and state/reaction_layout.f90 is a `CONTROLS_INPUTS` (Makefile:87-89). The regenerated file has to be committed, or test_generated_files.py fails.

       I traced the recipe statically through these rules. It otherwise holds: rv_*, the star's run copies, the controls and the probe in test_net all follow num_reactions.
     - **README:178.** The S-factor recipe's step 2 still has the newcomer run the three generators by hand. make has run them since 9906f760, and io/README §5's recipe, re-traced on 2026-11-01, says so. Harmless, but untraced since round 3.
     - **GLOSSARY:255 and HOW_YREC_WORKS:445** still name `star%cross_section_scale` and its siblings as star members.
   - **In the F77?** No; these are 2026 documents, true until 3a44ae33.
   - **Fix:**
     - rewrite rows 41, 43 and 44 from state/README.md:41 and :129, which have it right;
     - add an io/controls_buffer.f90 row (the scale controls and their defaults, the Li/Be references, the cutoffs and the tolerances);
     - in step 4, name the dump's lines, `pristine.txt` alone, and the regenerated controls_dump.inc;
     - update GLOSSARY:255 and HOW_YREC_WORKS:445;
     - re-trace the S-factor recipe through make.
   - **Kind** A. **Size** S.
2. **net/README.md:194 ("Add a light element in the trace approximation").**
   - **Issue:** the recipe reads "Follow deuterium: a rate function and a region integrator in `light_elements.f90` ..., driven per region from `core/step/burn_step.f90`". It names 2 files. Following deuterium touches these:
     - net/light_elements.f90: `deuterium_rate`, `burn_deuterium_region`, the thresholds;
     - net/reaction_rates.f90: the d(p,γ) fit and its overturn cap, :162-207;
     - net/nuclear_energy.f90: its energy, `compute_deuterium_burning`, added to pp I at :198 and :262;
     - net/net_lib.f90: five exports;
     - state/composition_layout.f90: the slot;
     - state/star_info_lib.f90: `deuterium_burning_rate` and `_start`;
     - core/step/burn_step.f90: `compute_burn_rates`, `burn_deuterium`, `surface_accretion_now`;
     - core/step/composition_step.f90: the call, gated on `use_extended_composition`, :183;
     - core/model/rezoning.f90: the regrid of the start-of-step rates, :774-810;
     - core/run/monte_carlo.f90: the per-run reset of the start rates, :133-148 (M169);
     - rotation/mid_timestep_model.f90 and rotation/rotation_scratch_lib.f90: the sub-step's rates and burn, :266-277 and :69-73;
     - test_net;
     - and Phase 5's 15-slot declarations for the new slot.

     A newcomer who follows the recipe leaves the new element unburned in the rotation sub-steps, unregridded at a rezone and stale across Monte-Carlo runs (M169's own bug, which this recipe would re-create). README:167 lists the four traced recipes, and this is not one of them; the standard asks for every file.
   - **In the F77?** No; a 2026 recipe.
   - **Fix:** trace it on a throwaway branch, as the isotope recipe was, and list the files in order, stopping at Phase 5 where the trace must. Otherwise, say that it is untraced and give the deuterium call graph above.
   - **Kind** A. **Size** S.
3. **core/step/burn_step.f90 (net's star-layer driver): :243, :244, :245, :256, :263, :289, :321, :333, :544; :90-91.**
   - **Issue:**
     - **Literal thresholds.** The light-element and network bookkeeping keeps 11 unnamed thresholds:
       - the burning-rate floor `1.0d-32` five times (:243 twice, :263 twice, :321);
       - the abundance floor `1.0d-24` four times (:244, :256, :289, :333);
       - the log T above which Li/Be are cleared, `7.0d0` (:245);
       - the hydrogen exhaustion that switches the metallicity rule, `5.0d-7` (:544).

       The standard asks for named thresholds with unit and meaning, and net/ meets it (`min_deuterium_step`, `min_composition_rate`, ...).
     - **F77 shapes.** The file also has 23 `.lt./.le./.gt./.ge./.eq.` comparisons, 40 `double precision` declarations and explicit-shape `(json)` arguments.
     - **Loose intents.** Five read-only drivers declare the star `intent(inout)`: surface_accretion_now :49, compute_zone_energy :74, average_li_be_over_cz :114, burn_regions :573 and network_rates_of_model :609. Their contracts say `star writes: none`.
     - **Stale path.** :90-91 points to "rotation/burn_li_be_midpoint.f90", which d51f2ea7 (2026-10-28) moved to rotation/retired/burn_li_be_midpoint.f90.txt.
   - **In the F77?** The literals, yes. `git show 6cd56739:src/liburn.f` has them at :134 (`RBE9(I).LE.1.0D-32`), :135-136 and :220-224 (`.LT.1.0D-24`) and :137 (`HT(I).GT.7.0D0`); kemcom.f:347 has `HCOMPP(1,IEND).LT.5.0D-7`. The path came with the 2026 move.
   - **Fix:**
     - name the four thresholds as parameters (`min_li_be_rate`, `min_li_be_abundance`, `li_be_clear_log_temperature`, `hydrogen_exhausted_fraction`), each with its meaning;
     - convert the comparisons and declarations in the byte-identical way net/ was;
     - tighten the five intents to `intent(in)`;
     - fix the path.
   - **Kind** A. **Size** S.
4. **net/screening.f90:204-209, :222-224, :243-245.**
   - **Issue:** the intermediate-screening exponent, `0.38*lambda0_086*charge_moment_158*graboske_z(z86,...)/(inv_mu_i*z_tilde_058*mean_charge_028)`, is written twice: in compute_reaction_screening (:204-205) and in weak_intermediate_screening (:243-245). Its derivative pair, `0.43*` and `-1.29*`, is written twice inside compute_reaction_screening (:207-209 and :222-224).
     - This is against code style 2, "each formula exists once".
     - Graboske et al.'s coefficients are bare literals (0.38, 0.43, 1.29; and in the strong form 0.624, 0.316, 0.737, 0.208), against code style 4.
     - test_net's independent copy of the strong formula (test_net.f90:775-777, :803-805) has no named source constant to cite.
     - A coefficient corrected in one copy would leave Be7 + p and hep screened by the other.
   - **In the F77?** Yes. `git show 6cd56739:src/rates.f` has `UINT=0.38*XXL8*XTR*Z86(I)/(AMU*Z58*Z28)` at :407 and again for Be7 + p at :509 (`Z86BE7P`); engeb.f:624-628 repeats both, with `0.43*UINT` and `-1.29*UINT`.
   - **Fix:** one pure `intermediate_exponent(plasma, s)` that both routines call; the derivative pair set once; the seven coefficients as named parameters citing Graboske et al. (1973), eq. 19 and table 4. Byte-identical if the expressions keep their order.
   - **Kind** A. **Size** S.
5. **io/controls_buffer.f90:607 (shown to users as defaults/physics.defaults:433); net/README.md:76, :87; net/network.f90:282; net/reaction_rates.f90:143-144.**
   - **Issue:** what each element of `nuclear_logT_cutoffs` does is stated only in code.
     - The control's description, the one users read, is "log T cutoffs (5) below which burning groups are switched off". It names no element.
     - It does not say that (2) and (3) do nothing (net_configuration.f90:37-39; a deck that sets them gets no warning).
     - Element (5) gates the thermal neutrino losses, not a burning group.
     - Element (1) gates at or below (network.f90:292, nuclear_energy.f90:122), not below.
     - net's own texts disagree: README:76, :87, network.f90:282 and reaction_rates.f90:143-144 say "below", while README:142, net_configuration.f90:36 and nuclear_energy.f90:43 say "at or below".
     - M118 covers burn_step's own `<` and the helium comment, not these.
   - **In the F77?** The control's vagueness, yes: the F77 TCUT was documented only in the code (`git show 6cd56739:src/rates.f`:308 `IF(TL.LE.TCUT(1))`, :574 `IF(TL.LT.TCUT(4))`). The wording drift is 2026's.
   - **Fix:** a description that says, per element: (1) nothing burns at or below; (2), (3) unused; (4) helium burning off below; (5) thermal neutrino losses off at or below. "At or below" in the four comments. Optionally a `warn_inconsistent_controls` warning for (2) or (3) set to another value than the default.
   - **Kind** A. **Size** S.

## What would raise each score

- **Readability (9 -> 10):**
  - BPS with named terms (pair, photo, plasma), its equations in the header and the polynomial coefficients as named arrays (with finding 4's screening constants);
  - provenance out of the module headers into Sources (network.f90:42-50, reaction_rates.f90:2 and :447-448, light_elements.f90:5-7);
  - fix nuclear_energy.f90:100-103 and :358;
  - README rows and paragraphs under 1,000 characters.
- **Organization (8.5 -> 9.5):**
  - the Li/Be S(0) beside their fits (M127);
  - BPS into thermal_neutrinos, one module with two fits;
  - burn_implicitly and compute_abundance_rates in a module of their own;
  - drop net_lib's two dead exports;
  - net_report through run_messages' report_warning and report_error (M114: now a few lines);
  - burn_step's thresholds and shapes (finding 3);
  - net's derived constants (erg_per_mev, m_e c²/k, which screening has as 5.92986 and BPS as 5.9302, the Fermi-momentum and Λ₀ coefficients, 4π AU²) derived from physical_constants, or sourced beside it.
- **Maintainability (8.5 -> 9.5):**
  - fix M89, M90, M91, M127 and M176, each as a measured commit;
  - add M131's NaN guard with M102's code, and a test_net case for it;
  - make `DEF_RE` see `pure` (M141);
  - close the self-check gaps;
  - add test_net checks for:
    - the eps_nuc derivatives (ROADMAP:3108-3113);
    - each light-element reference against its fit's S(0);
    - `burn_deuterium_region`'s conservation (ΔX + ΔHe3 + ΔD = 0) and closed form (D₀ e^{-Λ} on one shell), so the light elements have ok/FAIL checks, not only pins;
    - `net_settings_t()`'s defaults against the controls' defaults, which net_configuration.f90:32-34 claims;
  - a first green CI run on GitHub.
- **Extensibility (7.5 -> 9):**
  - implied-shape tables (`reactions(*)`, `species(*)`, `chain_name(*)`, `itoh_slots(*)`), with the counts from `size`;
  - a burning group per chain (M146);
  - one named tracer set in composition_layout (M147);
  - screening charges derived from the reactants unless overridden;
  - scales addressed by reaction name;
  - Phase 5's literal 15;
  - a traced light-element recipe (finding 2);
  - either a registration path for rate laws or §4.5(b) amended to the integer-dispatch decision.
- **Docs (7 -> 9):**
  - findings 1, 2 and 5;
  - the open docs items: M127's "one edit", M141's README:11, M146, M147, M162 (point §6 at ROADMAP:3108-3113 and name the open M numbers, as rotation's README names its O items), M118, M121, the M90 comment, GLOSSARY:255;
  - a Status line in TARGET_ARCHITECTURE §4.5(b) for the screening, the light elements and the slot-addressed scales;
  - re-trace every recipe after round 3, ending each with its checks passing, as io and state did.

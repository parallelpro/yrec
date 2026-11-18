# convection_diffusion -- eighth evaluation (2026-11-17, yrec-next 1179afbf, read-only)

Scope: every production .f90 of src/convection/ (7 files) and src/element_diffusion/ (18 files) read whole, both READMEs, setup/convection_settings_builder.f90 and setup/diffusion_settings_builder.f90, both mutation manifests, test_convection.py, and the two test programs (test_convection.f90 lines 1-760 whole and the rest by its check list; test_element_diffusion.f90's head, its case driver and its pipeline checks whole, the rest by its check list). Outside the two directories I read the call sites the READMEs' flow sections name: core/step/composition_step.f90 (advance_composition, homogenize_mixed_regions, settle), setup/settings_install.f90, and by grep the callers in core/step/step_driver.f90, core/model/rezoning.f90, rotation/evolve_rotation.f90, rotation/am_regions.f90 and rotation/mid_timestep_model.f90. From the tools: check_boundaries.py's facade, STAR_BLIND_FILES and KERNEL_FILES entries for the two directories, style_counts.tsv's rows and check_style.py's exemptions. From the docs: GUIDELINES "Naming", "Declarations and interfaces", "Code style" and "A package README"; ROADMAP rounds 0, 1 and 1B; bug_sweep's table rows O1, O10, O21, M33, M96, M97, M112, M241-M258, M285-M289, M308 and the entries of O1, O21, M285 and M308; B_LISTS.md's convection and element_diffusion sections; FINDING_MAP.md's CD rows. Every F77 claim was checked with `git show 6cd56739:src/<file>.f`. Nothing in the tree was built or run. In the scratch directory (review8/review8_convection_diffusion/): style.py and its output style_counts.txt (lines, procedures and F77 shapes of the production code, with rotation beside), and counts.txt (checks, mutations, commits).

The seventh report scored both directories for the first time (convection 6.4, element_diffusion 5.9). Since then convection had round 1 (9 commits, no result changes) and round 1B (8 measured commits); element_diffusion had round 1 only. The owner froze element_diffusion's measured fixes on 2026-11-11 ("don't do anything to element_diffusion"), so its open defects are a decision, not neglect. They still weigh in its maintainability score, as the seventh report weighed them, because the score says what the code is.

Counts, production code only (style.py):

| | convection | element_diffusion | rotation (the bar) |
|---|---|---|---|
| code lines / procedures (pure) | 881 / 33 (13) | 1625 / 39 (20) | 4502 / 110 (10) |
| longest procedure, median (code lines) | 76, 18 | 151, 24 | 165, 22.5 |
| ALL-CAPS comment lines (seventh: 132, 353) | 0 | 0 | 0 |
| `double precision`, old relational operators | 0, 0 | 0, 0 | 344, 110 |
| modules without `private` | 0 | 0 | 1 |
| I/O statements outside the messages module | 0 | 0 | 0 |
| ok/FAIL checks = mutations in the manifest | 49 | 21 | (its own harness) |
| README coverage figure | 97.8% (2026-11-10) | 97.1% (2026-11-10) | 96.9% |

## Scores (1-10)

### convection/

- **Readability 6 -> 8.5.** The prose is rewritten in sentences, the gradient's equations and all three derivatives are written once in temperature_gradients' header, the routines are short named steps, and every threshold is a named parameter with its meaning.

  Left:
  - mlt_derivatives keeps the F77's `qa1t`, `qa3p`, `qvr`, `qddelt` behind a one-line legend (temperature_gradients.f90:307-308, :317-349);
  - long positional interfaces: compute_temperature_gradients takes 20 arguments (7 results and ierr), find_convection_zones 19 (10 results), where Code style rule 1 asks for a typed result; call sites read `..., .false., is_convective, 1.0d0, 1.0d0, log_teff, .false., ierr` (semiconvection.f90:194-197);
  - extend_core_edge's header says the edge moves "while ... the distance covered ... stays within v_p dt", but the first shell joins whatever the distance (finding 6);
  - `log_mixed_regions`' argument `num_zones` is a region count (convection_messages.f90:59-61);
  - the test program is in `double precision` with literal slots (`xa(1,i)`, `comp(13,i)`) and labels in the F77's names ("tgcut", "lovmax limit", "iov1 iov2").
- **Organization 7 -> 8.5.** One zone scan, one merge, one EOS-opacity-gradient procedure, overshoot.f90 owning H_P and the envelope overshoot extent for core and rotation, a messages module, the four kernels in STAR_BLIND_FILES and KERNEL_FILES.

  Left:
  - core's `homogenize_mixed_regions` (core/step/composition_step.f90:220-247) is still a copy of `homogenize_region`; the README says so and says core's round could take the export, and core's round has passed;
  - `pressure_scale_height` is exported and has no caller outside the directory, while core/diagnostics/turnover_timescale.f90:341-342, :357 and rotation/seculr/stability_coefficients.f90:187 compute their own H_P (finding 7);
  - the kernels write `.false.` into the caller's `convective_flag(num_zones+1)`, a side channel that no caller reads and that is out of bounds at the largest model (finding 1);
  - the messages are written to `run_log_unit` directly (M112's open part), and `log_mixed_regions` writes two lines on every call with an overshoot switch on, several times per model, with no gate.
- **Maintainability 6 -> 8.5.** 49 ok/FAIL checks, each with a mutation in the manifest: the formula of ∇_rad, the mixing-length solution against a bisection, both limits' scalings, the three derivatives against centred differences with and without spots, the zone limit's error, every kind of overshoot against α H_P, conservation in homogenization, and semiconvection on a model with a core. Round 1B fixed M241, M242, M249, M254, M285 and O10 as measured commits, and round 2B M96 and M97.

  Left:
  - O1 (open by the owner's decision): above `log10_t_cut` and after an unconverged cubic the derivatives are ∇_ad's, and no L derivative is returned;
  - M286 (noted);
  - findings 1, 2 and 6;
  - no pinned deck reaches semiconvection (README.md:74 says so); the unit test is its only pin;
  - the spot settings have no range check: a contrast above 1 makes S > 1, and a shell with ∇_rad > ∇_ad > ∇_rad/S then returns `convection_unstable_radicand` (temperature_gradients.f90:153-168; io/controls_check.f90 checks the overshoot settings only).
- **Extensibility 6.5 -> 8.5.** A setting follows the record, the builder and the registry's settings cell (17 registry rows carry a `convection_settings%` cell), and an overshoot prescription is one branch of one file that core and rotation take through `envelope_overshoot_depth`.

  Left:
  - a new result of the gradient (the L derivative O1 wants, say) changes a 20-argument positional interface at its five call sites;
  - a change to how a region is mixed must be made twice (the copy in composition_step);
  - there is no recipe for a second gradient law (a different mixing-length variant): compute_temperature_gradients is one law with no slot.
- **Docs 6.5 -> 9.** The README has all six parts of the standard. I checked §1's table against the modules, §2's flow against the call sites (settings_install, shell_physics, shell_equations, envelope_integrands, starting_model, composition_step's eight steps, step_driver, rezoning, evolve_rotation, am_regions, secular_transport), §3's error codes and message destinations against convection_messages, and traced both recipes' file lists (the registry cells, the builder's `star reads`, the three outside users of `envelope_overshoot_depth`).

  Left:
  - the coverage figure predates round 1B's eight code commits (finding 8);
  - README.md:7 says core and rotation use the pressure scale height (finding 7);
  - docs/TARGET_ARCHITECTURE.md:508 and docs/bug_sweep.md's M308 are stale (finding 9).

Average 8.6 (8.5 / 8.5 / 8.5 / 8.5 / 9).

### element_diffusion/

- **Readability 5.5 -> 8.** Every header now states its equations and units; the legacy names, typos, `data` statements, explicit-shape dummies and missing `private` are gone; the scratch fields are named for what they hold.

  Left:
  - the message texts are the F77's capitals ("MICRODIFF FAILED TO CONVERGE ... IN  10ITERATIONS", element_diffusion_messages.f90:56-101), kept for the frozen B phase (B list item 6);
  - gravitational_settling_setup is 131 code lines: the Thoul fit's polynomials and their derivatives are unnamed literals (original_settling.f90:317-323, :336-343), and `coulomb_log_choice` is tested as bare 2, 3, 4 (:298, :357, :377-380);
  - thoul_diffusion is 151 lines of the paper's index arithmetic (`i-num_species` 18 times, thoul_solver.f90:147-169);
  - solve_diffusion_term still computes a midpoint change and a zero derivative that it multiplies together, and replays the F77's pass count for the messages (microdiff_run.f90:165-173, :190-197; finding 10);
  - the settings record keeps `fgry` and `fgrz`.
- **Organization 6.5 -> 8.** One region finder, one Coulomb logarithm, one grid builder and map back, one unit set, the sub-step policy in the domain for both callers, a messages module, six star-blind kernels and four settings-reading ones in the checker's lists.

  Left:
  - two pipelines, by the frozen decision: the original one makes `diffuse_elements` take the timestep, d ln P/dr, the enclosed mass and the total mass as `intent(inout)`, converted in place and back (element_diffusion_lib.f90:56-64), and keeps the scratch type and settling_grid alive for code no deck runs;
  - `apply_grid_changes` carries both pipelines' short-grid guards behind a flag (diffusion_grid.f90:245-246; O21);
  - a numerics failure is reported to standard output only (element_diffusion_messages.f90:109-116; B list item 6).
- **Maintainability 5.5 -> 7.** The seventh report's main gap is closed: 21 ok/FAIL checks, each with a mutation, cover conservation on the grid and through the map back, closure of the mass fractions, the sign of each species' change, linearity in a small timestep, the region, the sub-step plan, the grids' exactness, the Coulomb fit's limits and the Burgers solution's constraints. What holds the score at 7 is the list of open defects, every one frozen by the owner's decision, and the thin mutation cover of the coefficients.

  Left (a follow-up round's work list, if the owner lifts the freeze):
  - B list item 1: retire the original routines, or fix and test them (M243, M244, M245, M246, M288, and finding 4);
  - B list item 2, M33: the pin "old: without the Thoul et al. coefficients" holds one draw of unset memory and differs in about 4 runs of 100 (README.md:104);
  - B list item 3, O21: a 2-point grid's stencil reads outside the grid;
  - B list item 7, M247: the first midpoint's extrapolation, which acts in a pinned deck;
  - B list items 8 and 9, M287 and M257's 2.2;
  - M289 and M250 (decisions);
  - findings 3 and 5 (two guards that change no shipped deck);
  - the microdiff coefficients, the code every settling deck runs, have one mutation (microdiff_coefficients.f90) and microdiff_run one; no check tests D1's r² T^(5/2) scaling, the hydrogen-gradient term, or the `fgry`/`fgrz` factors;
  - several checks bundle assertions under one mutation ("helium, the metals and Li6, Li7, Be9 sink; hydrogen rises" has seven);
  - of the expected output's 19 pinned cases (560 of its 582 lines), 13 are the original routines'.
- **Extensibility 6 -> 8.** A light element is a row of species_table (I traced the recipe: the grid's rows, microdiff's loop and the map back all take `num_light_diffused`), a control follows the settings path (13 registry rows carry a `diffusion_settings%` cell), and the README's four recipes name their files and commands.

  Left:
  - a heavy species other than "the metals as iron" has no recipe: column 3 and `thoul_fe` are fixed in microdiff (microdiff_driver.f90:112-134);
  - `fgrz` scales the metals and the light elements together (M250);
  - a new treatment must sit beside two pipelines until the decision on the original one.
- **Docs 6 -> 8.5.** The README has the six parts; its "which setting acts where" table answers the seventh report's finding 17, and its flow matches both callers (composition_step.f90:258-300 read against §2).

  Left:
  - the two control descriptions M258 leaves open (io/controls_buffer.f90:1163-1166; frozen);
  - the table's row for `settling_num_iterations` omits M289's case (finding 10);
  - docs/TARGET_ARCHITECTURE.md:553 (finding 9).

Average 7.9 (8 / 8 / 7 / 8 / 8.5).

## Summary

- **Both directories now meet the standard in shape.** A messages module each, kernels that hold no I/O statement and are in the checker's lists, style counts at zero, a README with flow, conventions, checks, coverage and traced recipes, and test programs whose checks can fail. I checked the ROADMAP's claims for rounds 1 and 1B against the code (one scan and one merge, the triplet, overshoot.f90's ownership, M241, M242, M248, M249, M285's `max`, the dropped `diffuse_helium`, the single implicit solve, the table-sized light rows, `plan_settling_substeps` in both callers) and found each as described.
- **convection is at 8.6.** What is left is small: three inherited behaviours that bug_sweep does not have (a write one element past the flag arrays at 5000 shells, a merge that is not the union when the upper region reaches below the lower one, semiconvection's first shell joining whatever the timestep), a copy of the homogenization in core, long positional interfaces, and a coverage figure one round old.
- **element_diffusion is at 7.9, with maintainability at 7.** Its readability, structure and tests rose as far as byte-identical work could take them. Its open defects are the owner's frozen B list (M33, M243-M247, M250, M257, M287-M289, O21); they weigh here as they did in the seventh report. Three more inherited defects are new: a 0/0 at Z = 0, the original routines' derivative factors, and a degenerate one-zone region.
- **Root docs lag in two places:** TARGET_ARCHITECTURE's routine table for both directories, and bug_sweep's M308, which is fixed in the code.

## Previous report's items: status

Findings (the seventh report's numbers; FINDING_MAP's bug_sweep numbers):

1. Semiconvection's `only_check_core` (M241): **DONE**, 5a370511. semiconvection.f90:61-62 returns without a mixed region, without a core, and for a fully convective model; the switch and the downward branch are gone; test_convection checks all three.
2. The starspot Jacobian (M242): **DONE**, e6f70a34. temperature_gradients.f90:151-157 and :325-326 use ∇_rad/S and d ln S/d ln T; finite-difference checks with both spot modes.
3. The original routines' four defects (M243-M246): **OPEN** (frozen). Registered, and named in original_settling.f90:27-32, :140-142, :355-356.
4. The first midpoint's search (M247): **OPEN** (frozen; its fix moves a pinned deck). Documented at diffusion_grid.f90:103-111 and README.md:93.
5. The twelfth zone (M248): **DONE**, 70916805. convection_zones.f90:196-200 returns `convection_too_many_zones`; homogenization has no limit; two checks.
6. The stale overshoot region (M249): **DONE**, a21950d9. convection_zones.f90:90-92; overshoot.f90:104-105, :133-134 (the surface zone only).
7. `fgrz` on Li and Be (M250): **OPEN** (decision); documented in microdiff_coefficients.f90:17-20, the settings record and README.md:83. The control's description is unchanged.
8. `diffuse_helium` always true (M251): **DONE**, 81f856bd.
9. `light(3,json)` (M252): **DONE**, 60248a57 (microdiff_mapping.f90:30).
10. microdiff's repeated solve (M253): **DONE**, 8aac54fa; the pass count is replayed for the messages (see finding 10).
11. Overflowing formats (M254): **DONE**, 76a19fc7 (`i5`, `1pe13.5`, a format that stops after the last region).
12. Nothing reaches semiconvection (M255): **DONE** for the test (eight semiconvection checks) and the README's sentence; no pinned deck reaches it still.
13. Thin recipes (M256): **DONE**, 875dc67b and aefa6f84; traced here.
14. Formulas more than once (CD-14): **DONE** in convection (scan, merge, triplet, H_P and the envelope extent in overshoot.f90, with the measured commits 4c1c730a, 86b54640, 4cf4e0ac) except core's homogenization copy. **DONE** in element_diffusion (region, Coulomb logarithm with named constants, grids, units, the year, the sub-step policy) except the single-precision 2.2 (M257, frozen).
15. Messages (CD-15): **DONE** for both modules and the kernel lists; convection's square-root warning reaches the run log. **PARTLY**: element_diffusion's numerics failure goes to standard output only (frozen), `log_mixed_regions` is ungated, and both modules write the run log's unit directly (M112).
16. Stale text (M257): **DONE** (the README's units paragraph, microdiff_etm's five unused arguments, the D1 formula, 1990, the legacy file names); the 2.2 is open.
17. The controls' reach (M258): **DONE** in the README (§3's table, §2's two callers); the two descriptions in io are **OPEN** (frozen).

The scores' "Left" items of the seventh report:

- convection: the ALL-CAPS prose, history comments, operators and intrinsics, the long routines, the derivative block's algebra, the two misleading names, the dead branch and unreachable test, the literal slots and counts, the unnamed numbers: all **DONE**. O10, M96 and M97: **DONE** (rounds 1B and 2B). O1: **OPEN** (decision). The flow section, the messages convention, the coverage figure: **DONE**.
- element_diffusion: the prose, typos, legacy names, `private`, `only:`, assumed shape, the scratch names, the `data` statements, the duplicated stages, the sub-step policy, the flow section, the settings table, the coverage figure, the recipe: all **DONE**. O21 and M33: **OPEN** (frozen). The second pipeline: **OPEN** (frozen).

## New findings (not in bug_sweep)

Bugs first (1-6), then docs and structure (7-11).

1. **convection/convection_zones.f90:93, :298 (callers' arrays: state/star_info_lib.f90:542, :633; core/step/composition_step.f90:113)**
   - **Issue:** find_convection_zones and homogenize_convection_zones set `convective_flag(num_zones+1) = .false.`. Every caller's flag array has `json` (5000) elements, and a model may have exactly `json` shells (io/model_files.f90:714 accepts `nz <= json`; core/model/rezoning.f90:631, :736 and :829 refuse only more). For a 5000-shell model the write lands one element past the array: in the star, on whatever follows `convective_flag`; in composition_step, past the local `deep_mix_flag`. A devel build stops on it. The dummy is assumed-shape, so nothing in the kernel says the array must be longer than the model.
     - No reader needs it: the directory's own scans stop at `num_zones` (next_convective_run), and I found no caller that reads element `nz + 1` (the `flag(k+1)` reads in turnover_timescale.f90:196-198, :276 and starting_model.f90:890-892 stay below `nz`). README.md:45 and convection_zones.f90:47-48 call it a sentinel "callers' own scans read".
   - **F77:** yes (convec.f:55, `LC(M+1) = .FALSE.`, with `LC(JSON)`).
   - **Fix:** drop the two writes, or guard them with `num_zones < size(convective_flag)`; correct the README and the header; a check with a flag array of exactly `num_zones` elements under the devel flags.
   - **Kind:** A. **Size:** S.
2. **convection/convection_zones.f90:228-235 (header :208-216), with overshoot.f90:128-134 and :144-146**
   - **Issue:** the merge gives the merged region region k's first shell and the higher last shell, so it is the union only when region k+1 starts at or above region k's first shell. Downward overshoot can put region k+1's first shell below it: a surface convection zone whose α_env H_P reaches past a thin intermediate zone that does not itself overshoot (`envelope_overshoot` on, `intermediate_overshoot` off, or a smaller α for it).
     - Example: zones [18-19] and [22-40], the envelope's edge moved to 14. The merge gives [18-40]; shells 14-17 leave the mixed region, though `overshoot_bottom` = 14 still marks them for the adiabatic gradient, and they are listed as radiative.
     - If the envelope's edge passes below region k−1's top as well, that pair is never compared again (the loop moves forward only).
     - M285 fixed the same defect at the top (`max` of the two last shells) and its entry's "a merge is the union of the two regions" is the rule this case breaks. No test case has an upper region reaching below a lower one.
   - **F77:** yes (convec.f:129-141, `MXZONE(J+1,1) = MXZONE(J,1)`).
   - **Fix:** `merged = [min(bounds(k,1), bounds(k+1,1)), max(bounds(k,2), bounds(k+1,2))]`, stepping k back one row after a merge that lowered the first shell; a check and its mutation.
   - **Kind:** B (models with a thin zone within the envelope's overshoot; no pinned deck known to reach it, not measured). **Size:** S.
3. **element_diffusion/diffusion_grid.f90:298-300**
   - **Issue:** change_zone scales the CNO tracers by `z_ratio = z_new/xa(i_metals)`. With Z = 0 in a zone and the metals settling (`use_diffusion_z`, which `diffuse_lithium` also turns on, M181), the metals' change is exactly zero (the coefficients carry a factor of the abundance), so `z_new = 0`, the ratio is 0/0, and the seven tracer slots `i_c12:i_o18` become NaN in every zone at the first settling step. A devel build stops on the division.
     - Latent: a zero-metallicity model with settling. No shipped starting model is named for Z = 0; I did not scan the decks.
   - **F77:** yes (microdiff_etm.f:28, :87, :118 and equal_to_model.f:29, :80, :99, `ZZ2 = ZZ/HCOMP(3,I)`).
   - **Fix:** leave the tracers alone when `xa(i_metals) <= 0`; a check on a Z = 0 model. Byte-identical for every model with Z > 0.
   - **Kind:** A. **Size:** S.
4. **element_diffusion/original_settling.f90:320-329, :344-351**
   - **Issue:** two more inherited inconsistencies of the original routines, beside M243-M246 and M288.
     - `constant_mixing_coeff` multiplies D1 and D2 (:324-326) but not their X derivatives (:327-329), and `constant_settling_reduction` multiplies the metals' D1 (:344-346) but not its Z derivative (:350-351). With the factors at 0.5 and 0.8 (the test's case "constant mixing 0.5, settling reduction 0.8") the Lax-Wendroff midpoint correction and the implicit step's D2 correction are 2 and 1.25 times too large.
     - With Burgers' equations solved (`use_thoul_fit` off), only A_P and A_T come from the solution; D2's A_H and all three X derivatives are still the fit's polynomials (:320-323), so D2 ignores `coulomb_log_choice` 4's pairwise logarithms.
     - No deck runs these routines; only test_element_diffusion's pins see it.
   - **F77:** yes (setup_grsett.f:236-244, `COD1(I) = CSTMIXING*FAC*...` beside `QCOD1X(I) = FAC*...`; :264-269, `COD1Z(I) = CSTDIFFMIX*CSTMIXING*...` beside `QCOD1Z(I) = CSTMIXING*FAC*...`).
   - **Fix:** enter it in bug_sweep with M243-M246 for the same decision (retire, or fix and test).
   - **Kind:** B (test pins only). **Size:** S.
5. **element_diffusion/settling_region.f90:72-81, with diffusion_grid.f90:68-77**
   - **Issue:** find_settling_region moves `zone_end` inward to the first zone with Y above `helium_diffusion_min` and returns `region_found` without comparing it with `zone_begin`. When Y exceeds the minimum at the core's top zone and nowhere above it (a helium-rich convective core under an envelope whose Y is below a large `helium_diffusion_min`; the test's own floor case uses 0.5), the region is one zone. size_equal_grid then has `span = 0` and `narrowest = 0`, takes `int(0/0)` for the point count and `mod(0, 0)`, and the grids are built from NaN. A search that passes below `zone_begin` gives a negative span instead.
     - Latent: the default minimum is 1e-3, which every envelope exceeds.
   - **F77:** yes (microdiff_setup.f:64-101: `IBEGIN` from the core and the X search, `IEND` from the envelope and `DO 45 I=IEND,1,-1 / IF(HCOMP(2,I).GT.YMIN) GOTO 47`, with no comparison of the two).
   - **Fix:** return a no-region verdict when `zone_end <= zone_begin`; a check.
   - **Kind:** A. **Size:** S.
6. **convection/semiconvection.f90:133-158 (header :78-85)**
   - **Issue:** the propagation limit never applies to the first shell above the core. The loop starts at `outside + 1` and tests the distance before it evaluates a shell, so its earliest exit leaves `new_edge = outside`: once the first shell is unstable with the core's composition, it joins the core whatever v_p dt is. test_convection states it ("a 1 s timestep adds only the first shell"). The header reads as if a small enough timestep moved nothing. Semiconvection runs on several passes of a step, each from the zones found afresh, so the core's mixed edge is at least one shell out on every pass that runs it.
   - **F77:** yes (sconvec.f, the same loop from `JMR + 1`).
   - **Fix:** the authors' choice: test the first shell's weighted distance against v_p dt as well, or say in the header and in README §3's "Kept by decision" that the first shell is exempt.
   - **Kind:** B to change (decks with semiconvection; none pinned); A to document. **Size:** S.
7. **convection/convection_lib.f90:16, :23; convection/README.md:7**
   - **Issue:** `pressure_scale_height` is exported, and README.md:7 lists "the pressure scale height ... that core's light-element burn and rezone and rotation's sub-steps use". No code outside the directory calls it (the three users take `envelope_overshoot_depth`; test_convection is its only outside caller). H_P is still computed by hand at core/diagnostics/turnover_timescale.f90:341-342 and :357 and at rotation/seculr/stability_coefficients.f90:187 (with their own gravity). check_boundaries' facade list does not name it, since function references are not checked (M304).
   - **Fix:** have the two outside formulas call it where the gravity agrees (a measured commit where it does not), or make it private to overshoot.f90 and correct the README's sentence.
   - **Kind:** A for the export and the text; B to merge a formula that differs. **Size:** S.
8. **convection/README.md:78**
   - **Issue:** the coverage figure (97.8% of 465 lines, 2026-11-10) was measured before round 1B's eight commits changed the code it covers: semiconvection lost its downward and intermediate branches, overshoot.f90 gained `envelope_overshoot_depth`'s gate, the messages were rewritten. The paragraph notes the semiconvection change but keeps the old number and line count, and its per-file percentages (semiconvection 98.7%, convection_messages 88.5%) describe code that is gone.
   - **Fix:** rerun the command the paragraph gives and replace the figures.
   - **Kind:** A. **Size:** S.
9. **docs/TARGET_ARCHITECTURE.md:508, :553; docs/bug_sweep.md:480 and M308's entry**
   - **Issue:**
     - TARGET_ARCHITECTURE.md:508 puts `oversh` and `hsubp` in convection_zones.f90 as `extend_by_overshoot` and `compute_scale_height`. They are in overshoot.f90 since round 1, and `compute_scale_height` no longer exists (`edge_scale_height`, `pressure_scale_height`).
     - :553 gives microdiff's four routines two files; `microdiff_setup` is `convert_to_bahcall_loeb`, the grids are diffusion_grid's, and the table has no row for `model_to_equal`/`equal_to_model` (settling_grid.f90), `microdiff_run` or `microdiff_cod`. ROADMAP's round 1 notes the table's eos and kap rows as stale, not these.
     - bug_sweep lists M308 as "Open (test)". It is fixed: ebf4d507 made `boot_yale_eos` call `eos_set_mixture` (convection/test/test_convection.f90:124), and the expected output has no Saha line.
   - **Fix:** the three rows and the M308 status.
   - **Kind:** A. **Size:** S.
10. **element_diffusion/README.md:85; microdiff_run.f90:164-173, :190-197**
    - **Issue:**
      - The settings table says `settling_num_iterations` and `settling_tolerance` decide "only the messages" in microdiff. With `settling_num_iterations < 1` there is no implicit solve at all, so the diffusion term D2 is dropped (microdiff_run.f90:164; M289). The table should say so while M289 is open.
      - solve_diffusion_term still fills `change_mid` with the drift's midpoint change (stored one slot off, as M245 describes for the original routines) and a `zero_derivative` array, and passes both to implicit_diffusion_coeffs, which multiplies them: dead work that reads as if the coefficients were corrected. The loop at :190-197 then replays the F77's pass count so that the messages match.
    - **Fix:** the table's cell; pass no derivative to implicit_diffusion_coeffs from microdiff (an optional argument, byte-identical for finite changes), and drop the replay when the messages are reworded.
    - **Kind:** A. **Size:** S.
11. **convection/temperature_gradients.f90:106-111; convection_zones.f90:70-75**
    - **Issue:** the two main entries return their results as long tails of `intent(out)` arguments (7 and 10), against Code style rule 1 ("a routine with several results returns a typed record"). The three derivatives are documented as "set only for a convective shell with want_derivatives", so a caller that reads them after a radiative shell reads undefined values, and nothing in the interface says which results belong together. The five call sites of the gradient (four in core, one in semiconvection) pass rotation factors, log Teff and the overshoot flag by position.
    - **Fix:** a `temperature_gradient_result` record (gradient, ∇_rad, velocity, convective flag, the derivatives with a flag that says they are set) and a zone-list record for find_convection_zones; byte-identical.
    - **Kind:** A. **Size:** M.

## What would raise each score

**convection/**
- **Readability, to 9+:** finding 11's records; mlt_derivatives' names written out (`dlna1_dlnt`, ...); the extend_core_edge header (finding 6); the test program in `real(dp)` with named slots.
- **Organization, to 9+:** core's homogenization through `homogenize_region`; finding 7; finding 1's sentinel removed; `log_mixed_regions` under the solver-diagnostics gate or once per model (a re-pin of the logs that carry it); the messages through run_messages' routines when M112 closes.
- **Maintainability, to 9+:** findings 1 and 2 with their checks; a pinned deck with a convective core and semiconvection; a range check of the spot settings in io/controls_check.f90; O1's convection part when the owner takes it up (the L derivative, with a finite-difference check).
- **Extensibility, to 9:** the result record (a new result touches no call site); a recipe for a second gradient law.
- **Docs, to 9.5:** findings 8 and 9; README.md:7.

**element_diffusion/**
- **Readability, to 9:** the messages reworded (B list item 6); the fit's polynomials as named coefficient arrays with their derivatives computed from them; named Coulomb-logarithm choices; finding 10's dead work removed.
- **Organization, to 9:** the decision on the original routines (retired, `diffuse_elements` loses its `intent(inout)` round trip, the scratch type and settling_grid go to retired/); one short-grid rule (O21).
- **Maintainability, to 8.5-9 (each a measured commit, or a decision, after the freeze is lifted):**
  - M33, so that every pin is deterministic;
  - O21 and M247, with a convective-core case that reaches M247's path and a 2- and 3-point grid case;
  - M243-M246, M288 and finding 4 fixed with a check each, or the routines retired with their pins;
  - M287 and M257's 2.2;
  - findings 3 and 5 (guards, byte-identical on the shipped decks);
  - checks of the microdiff coefficients with their mutations: D1 ∝ r² T^(5/2), `fgry` and `fgrz` linear, the hydrogen-gradient term's sign, a settling timescale against Michaud's estimate;
  - the bundled checks split so that each assertion has its mutation.
- **Extensibility, to 9:** a recipe for a separately settled heavy species; M250's decision.
- **Docs, to 9:** M258's two descriptions; finding 10's table cell; finding 9's row.

Nothing was run.

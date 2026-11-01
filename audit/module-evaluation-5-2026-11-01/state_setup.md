# state/ + setup/ -- fifth evaluation (2026-11-01, yrec-next 6a75b72f, read-only)

What I read whole: everything under src/state/ and src/setup/ (.f90, .inc, both READMEs), src/tools/gen_star_docs.py and src/tools/check_star_contracts.py, and check_boundaries.py's state and setup rules (lines 1-120 and 255-1000). For context I also read core/run/monte_carlo.f90, core/run/run_driver.f90, the relevant parts of run_setup, starting_model, surface_boundary, envelope_integration, io/controls_reader and io/user_inputs, and the fourth report, ROADMAP "Round 3", bug_sweep §10-11, TARGET_ARCHITECTURE 4.2/4.6 and GUIDELINES (Naming, six rules, Declarations, Code style).

Checks I ran (`python3 -B`, tree unchanged before and after by `git status`):
- `gen_star_docs.py --check`: "outputs up to date".
- `check_star_contracts.py`: "Star contracts OK (92 declared)".
- Both tools' `--self-test`: all ok.

Sweeps (grep, plus one scratch script over the registry):
- Nothing outside `io/controls_sync.f90` assigns `s%control` or `s%physics`.
- No STAR_FREE dir and not rotation reads either record.
- None of the controls read's files writes a star member.
- `arm_run_copies` has 33 lines. Of the 33 controls with a run copy, only `physics%atm_choice` is read outside it, by the atm builder (the tables) and by surface_boundary's restore. Both are by design (M140).
- Contract coverage by directory (star-taking procedures with a contract / all of them):

| Directory | With a contract / all |
|---|---|
| core/run | 24/24 |
| core/step | 25/25 |
| setup | 15/15 |
| io | 26/29 (the 3 without are controls_sync's, exempt) |
| core/model | 1/4 |
| core/structure | 1/5 |
| core/envelope | 0/5 |
| core/diagnostics | 0/26 |
| rotation | 0/16 |
| **Total** | **92/149** |

## Scores (1-10)

The bar is rotation/'s. The first three earlier scores were on a 1-5 scale and are doubled here. The order is first -> second -> third -> fourth -> now.

**state/**
- Readability: 5 -> 7 -> 8 -> 8 -> **8.5**.
  - Better: real(dp) throughout, honest Monte-Carlo sample names, the README opening in four labelled paragraphs.
  - Left: one of those paragraphs is still 1,867 characters, the "Where the checks are" sentence is 980, and §6/§7 bullets run to 1,000-1,560.
- Organization: 4 -> 6 -> 7 -> 7.5 -> **8**.
  - Better: a Monte-Carlo class with a per-run re-arm hook; atm_choice has one home; the pseudo-controls are gone; the read writes only its records and the run copies.
  - Left: the per-Monte-Carlo-run values are split across two classes (new 5); a latch that is not reset with its run copy (new 1); run_messages has landed in a "types and constants" directory without being described there; the flat 265-member type.
- Maintainability: 5 -> 7 -> 8 -> 8 -> **8.5**.
  - Better: contracts are required in core/run, core/step, setup and io; the scanner gaps M116/M168 are fixed, with self-tests; the read's own copies are checked; mc_samples' fields are documented and checked; max_runs has one source.
  - Left: check_max_runs is now vacuous (new 4); 57 of 149 star-taking procedures have no contract; M169 has no test.
- Extensibility: 5 -> 7 -> 8 -> 8 -> **8.5**.
  - Better: all three state recipes are traced, including per-Monte-Carlo run, and the run-copy recipe is complete with the dump step.
  - Left: the per-Monte-Carlo-run recipe points to one class while the existing values of that lifetime sit in another; there is no rule for a latch that is paired with a run copy.
- Docs accuracy: 5 -> 7 -> 7 -> 7.5 -> **8**.
  - Better: the generated tables are right by construction, and M163's four stale statements are fixed.
  - Left: the README has no row for the new run_messages module, and §6 calls state/ "types and constants" (new 2); control_slots' header makes three false statements (new 4); a few member comments are stale (new 7); M161 is open.

**setup/**
- Readability: 8 -> 8 -> 9 -> 9 -> **9.5**.
  - Ten pure one-line-copy builders, one naming scheme (`<x>_settings_builder`), and settings_install's header is a table of lifetimes.
  - Left: the README's §1 carries a 1,222-character row about io's map_user_inputs.
- Organization: 6 -> 7 -> 8 -> 8 -> **9**.
  - Every install is in settings_install, one procedure per lifetime, including the three per-surface-mixture installs.
  - The atm builder takes the deck's choice, and the kap builder derives use_two_z_tables.
  - Left: M85/M96's shared recompute.
- Maintainability: 7 -> 7 -> 8 -> 8 -> **8.5**.
  - Holding: purity is compiler-checked, INSTALLER_CALL_RE is complete, check_settings_fields holds, and the contracts are required.
  - Left: check_settings_fields skips the three paths records, and nothing requires a builder to be `pure` (new 6); no check ties a field to its control (known).
- Extensibility: 6 -> 7 -> 8 -> 8 -> **9**.
  - The field, domain and per-surface-mixture recipes are traced, each naming the check that fails until its step is done.
- Docs accuracy: 6 -> 8 -> 8 -> 8 -> **8**.
  - The tables of builders and installs are right.
  - Round 3's merge left the README describing READ3's predecessor in four places (new 3):
    - the star members the builders read;
    - map_user_inputs' targets;
    - a read step that no longer exists, twice;
    - `s%envelope_store_log_temperature`.

## Summary

- **Round 3 did most of what the fourth report asked**, all checked in the code:
  - a Monte-Carlo class and a per-run hook (`begin_monte_carlo_run` calls the generated `arm_run_copies`);
  - one atm_choice home;
  - honest sample names;
  - contracts required in four directories;
  - real(dp);
  - one max_runs;
  - pure, uniformly named builders;
  - every install in settings_install.
- **The rules hold everywhere:**
  - only controls_sync writes the two records, and no domain reads them;
  - each run-copied control is read only by arm_run_copies (atm_choice's two readers are by design);
  - the generated tables and the 92 contracts pass.
- **What is left is mostly merge leftovers.** The two packages' texts were each true on their own branch. Four of the five problems are documentation:
  - setup/README still describes the read before READ3 (adopt_canonical_names; the rate scales as star members);
  - control_slots' header and check_max_runs describe a MAX_RUNS that "repeats" max_runs, while it now reads it, so the check cannot fail;
  - state/README never mentions the new state/run_messages.f90;
  - star_info_lib keeps the per-Monte-Carlo-run solar luminosity and rate-scale run copies in "Constants set at setup".
- **One behaviour item is new:** the gray-atmosphere latch `use_ttau_relation` is not reset when M95's re-arm (or a calibration cycle) restores `atm_choice_in_use`. The next run or cycle then prints a switch-back it did not make and forces a new triangle (probably result-neutral, not measured).
- **The fourth report's item 3 stays open by decision:** bug_sweep §10 judged `num_rotation_structure_passes` a real run copy.
- **Known open defects still weigh:** M85, M96, M97, M115, M117, M126, M161, M169, and TARGET_ARCHITECTURE 4.6's leftovers (the diagnostics gate and the read's buffer as module state).

## Fourth report's items: status

New findings 1-10:

| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | atm_choice in three homes | DONE | `atm_choice_initial` gone (grep: 0 hits); `atm_settings_from_star` copies `s%physics%atm_choice` (setup/atm_settings_builder.f90:30); surface_boundary restores from the control (:150-158); the read no longer resets `use_ttau_relation`; check_control_copies covers the read's own copies (check_boundaries.py:701-733). |
| 2 | settings_install not the one install place | DONE | `install_surface_mixture_settings`, `install_surface_z_settings` and `install_surface_x_settings` are in settings_install.f90:120-160; core calls only those (starting_model:1012,1015; model_rescaling:158; rezoning:1067; run_setup:66); INSTALLER_CALL_RE lists eos_set_mixture, kap_update_surface_tables and set_convection_envelope_metals (check_boundaries.py:301-303); README has the lifetime "per surface mixture". M85/M96's shared recompute is not done (both open). |
| 3 | num_rotation_structure_passes a pseudo run copy | OPEN, by decision | Still a run copy (controls_sync.f90:914; step_driver.f90:513-514). bug_sweep §10 "Left out" calls it a real run copy. The same pattern holds for `first_mc_run`/`last_mc_run` (new 9). |
| 4 | mc_samples names | DONE | `s11_factor` … `settling_z_factor`, `luminosity_factor`, `age_factor`; `diffusion_factor` dropped (star_info_lib.f90:135-162); the .mc1 format kept (monte_carlo.f90:256-262). |
| 5 | builders not pure | DONE | All ten builder functions are `pure function`. Not enforced by a check (new 6). |
| 6 | contract-coverage claim | DONE | REQUIRED_DIRS = core/run, core/step, setup, io (check_star_contracts.py:457), failing check plus self-test; docstring rewritten, no global star. |
| 7 | run-copy recipe missing steps | DONE | state/README.md:488 runs gen_controls_state, gen_inlist and gen_controls_dump; :490 runs `--fix` and gen_star_docs. |
| 8 | max_runs in four places | DONE, with leftovers | Registry cells spell `real(max_runs)` (32 cells, 0 `(50)`); generated `.inc` spells `max_runs`; controls_registry.py:332 reads it from control_slots. But control_slots.f90:7-12 and :20-25 and check_boundaries.py:571-583 still describe the old design, and check_max_runs compares the value with itself (new 4). |
| 9 | stale cross-references and banners | DONE | star_info_lib.f90:1302 says section 4; README:15 states the re-exports exactly; GUIDELINES:208-209 points to the counts; the Constants banner's Valid names the Monte-Carlo rescale; mc_samples has left the Private class. |
| 10 | setup/README gaps | PARTLY | The coverage figure is restated (100% of 226 lines, user_inputs noted as io's) and use_two_z_tables is derived in the kap builder (§3 says so). But the rate-scale list added by M164 was overtaken at the merge: six of its eight names are now `physics%` controls and two are now `_in_use` run copies (new 3). |

What would raise each score:

| Criterion | Item | Status |
|---|---|---|
| state Readability | split README line 3 | PARTLY: four labelled paragraphs, but line 7 is still 1,867 characters |
| state Readability | rename mc_samples fields | DONE |
| state Readability | real(dp) | DONE (MEMBER_DECL_RE and star_member_types accept it) |
| state Organization | drop the pseudo copies | PARTLY (1 DONE, 3 declined) |
| state Organization | retire the five no-input controls | DONE (M152) |
| state Organization | rate scales and envelope_store_log_temperature as derived controls with run copies | DONE (physics_group_def.inc:40,48-60; controls_sync:911-912) |
| state Organization | one Monte-Carlo class with a per-run hook | DONE for the four members named. The other per-Monte-Carlo-run values stayed in Constants (new 5). |
| state Organization | Phase 5 typed parts | OPEN |
| state Maintainability | M116 | DONE |
| state Maintainability | enforce contract coverage | DONE for 4 directories (92/149 overall) |
| state Maintainability | check_control_copies on the read's copies | DONE |
| state Maintainability | gen_star_docs on mc_samples | DONE (parse_record_fields, field_access) |
| state Maintainability | one max_runs | DONE (stale comments and a vacuous check left) |
| state Extensibility | complete §7 | DONE |
| state Extensibility | per-Monte-Carlo-run recipe | DONE (traced, state/README.md:494-499) |
| state Docs | item 9 | DONE |
| state Docs | item 6's claim | DONE |
| state Docs | intent(inout) semantics gap between the two tools stated | DONE (README §3 and both docstrings) |
| setup Readability | pure builders | DONE |
| setup Readability | history parentheticals | DONE (none left in the builders' and install's headers) |
| setup Readability | one naming scheme | DONE (`<x>_settings_builder`; no `microphysics_settings`, `net_settings_setup` or `rotation_settings_setup` left outside audit/ and history) |
| setup Organization | item 2 | DONE |
| setup Organization | atm builder copies the control | DONE |
| setup Organization | use_two_z_tables | DONE |
| setup Maintainability | item 5 | DONE |
| setup Maintainability | INSTALLER_CALL_RE | DONE |
| setup Maintainability | control-to-field check | OPEN (design in setup/README §6) |
| setup Extensibility | per-surface-mixture lifetime with its procedure and recipe | DONE |
| setup Docs | item 10 | PARTLY (see above) |
| setup Docs | net_settings_setup header | DONE (file renamed, header new) |
| setup Docs | §3's true scope | DONE (the rule now covers the surface-mixture setters and is enforced) |

## New findings (not in bug_sweep), ranked

1. **The gray-atmosphere latch is not reset when its run copy is restored.**
   - Where: core/run/monte_carlo.f90:107; core/run/run_driver.f90:299; core/structure/surface_boundary.f90:149-162; core/envelope/envelope_integration.f90:187-193; state/star_info_lib.f90:422-426.
   - The issue:
     - `use_ttau_relation` means "the run has switched `atm_choice_in_use` to 0".
     - `begin_monte_carlo_run` re-arms `atm_choice_in_use` from the control. That is M95's fix, so this case is new in round 3; before it, both values carried over together.
     - The solar calibration's per-cycle restore (run_driver:299) does the same.
     - Neither resets the latch.
     - So a Monte-Carlo run or calibration cycle that ended on the gray fallback leaves the next one with the deck's atmosphere and the latch still set. At the first surface-boundary call with Teff inside the table, `switch_atmosphere_choice` prints "SWITCH BACK TO KURUCZ/ALLARD ATMOSPHERE", which never happened in that run, and sets `start_new_triangle`.
     - `settling_suspended_reported` carries over the same way (message only).
     - This is a concrete case under M169's heading, not among its candidates.
   - Fix: reset the latch wherever `atm_choice_in_use` is restored, or derive it, `s%atm_choice_in_use == 0 .and. s%physics%atm_choice /= 0`, and drop the member.
   - B? Probably result-neutral, because the card that follows reads its model, which already asks for a new triangle; a spurious terminal line. Not measured. Size S.

2. **state/README does not describe run_messages.f90, new in state/ in round 3.**
   - Where: state/README.md:3, :13-24, :466; state/run_messages.f90.
   - The issue:
     - The module table has no row for it. Line 3 ("the star and the records and slot constants it is made of") and §6 ("state/ holds types and constants") are no longer true: run_messages holds procedures and module state (the pending run-log lines, the two units, the solver-diagnostics gate).
     - Its two units have no initial value until io's `assign_unit_numbers` runs. test_atm and test_element_diffusion set them by hand, and any other program writes to unit 0 by accident.
   - Fix:
     - add a row (what it holds, that every layer may use it, its module state, the TARGET_ARCHITECTURE 4.6 item) and §6's checks (check_physics_no_io, ONLY_LIST_MODULES), and say which test covers it;
     - optionally initialize `terminal_unit = output_unit` and `run_log_unit = output_unit` (the log then counts as closed).
   - A (result-neutral), S.

3. **setup/README describes the controls read as it was before READ3 (a merge leftover).**
   - Where: setup/README.md:46, :50, :57, :64.
   - The issue:
     - :46 lists, as star members the builders read, `cross_section_scale`, `extra_cross_section_scale`, `s_prime_scale`, … `n15_pa_scale`. The net builder reads the run copies `cross_section_scale_in_use` and `extra_cross_section_scale_in_use`, and the other six are `physics%` controls (net_settings_builder.f90:41-48).
     - :50 says map_user_inputs writes "the star's scale arrays (`s%cross_section_scale`, …)". It writes the buffer (io/user_inputs.f90:47-50).
     - :50 and :57 name `adopt_canonical_names`, which exists nowhere in the code. :57's step list also leaves out expand_paths and interpret_kind_cards (controls_reader.f90:60-81).
     - :64 names `s%envelope_store_log_temperature`, which is `s%physics%envelope_store_log_temperature`.
   - Fix: the list becomes "the run copies `cross_section_scale_in_use` and `extra_cross_section_scale_in_use` (net)"; the read's steps as in read_controls; :50 points to io/README, which owns the read.
   - A, S.

4. **control_slots' header and check_max_runs describe the pre-merge max_runs (a merge leftover), and the check is vacuous.**
   - Where: state/control_slots.f90:7-12, :20-25; tools/check_boundaries.py:571-583.
   - The issue: the header says three things that are no longer true:
     - that controls_buffer re-exports the slots (it imports only `max_runs`; its users take the slots from control_slots, as state/README:23 says);
     - that gen_controls_state "writes its value" into the records (they spell `max_runs`);
     - that controls_registry's MAX_RUNS and the registry cells "repeat it" (MAX_RUNS is read from this file, controls_registry.py:332, and the cells spell `real(max_runs)`).
     - As a result, check_max_runs compares the parameter with itself and can fail only on a parse difference.
   - Fix: rewrite the two comments, and replace the check with one that fails on a literal run-card size in the registry or the generated files.
   - A, S.

5. **The per-Monte-Carlo-run values live in two classes.**
   - Where: state/star_info_lib.f90:21-23, :186-219, :221-234.
   - The issue:
     - `cross_section_scale_in_use` and `extra_cross_section_scale_in_use` are run copies that every Monte-Carlo run re-arms and rescales. `solar_luminosity_cgs`, `log10_solar_luminosity` and `ln10_over_lsun` are set from every run's sample.
     - All five sit in "Constants set at setup", whose banner says "fixed for the whole run" and which the header files under "per run".
     - The Monte-Carlo class then calls them "constants of the class above".
     - state/README §7's per-Monte-Carlo-run recipe puts the next such value in "Monte-Carlo runs", so one lifetime has two homes.
   - Fix: move the five into "Monte-Carlo runs" (or split Constants into per run, the solar mass and radius, and per Monte-Carlo run); regenerate.
   - A (layout only), S.

6. **check_settings_fields skips the three paths records, and no check requires a builder to be pure.**
   - Where: tools/check_boundaries.py:335-382; setup/README.md:71-72.
   - The issue:
     - The check matches only `\w+_settings_t` types and `result(settings)` functions.
     - So a field added to `kap_table_paths_t`, `eos_table_paths_t` or `atm_table_paths_t` by io/README's "Add a table path" step 3, and not copied in `<x>_table_paths_from_star`, keeps `''` silently until the domain's open fails.
     - §3's "a builder sets every field of its record" and "is a pure function … so the compiler rejects a side effect" are conventions: a new builder without `pure` passes every check.
   - Fix: match `\w+_(settings|table_paths)_t` and `result\((settings|paths)\)`, and fail on a builder without the `pure` prefix.
   - A, S.

7. **Four stale or overclaiming comments.**
   - Where:
     - star_info_lib.f90:442-446;
     - star_info_lib.f90:352-355;
     - state/README.md:11;
     - state/README.md:476.
   - The issue:
     - `envelope_hydrogen_fraction`'s comment says "the controls read and kap_init use it earlier, while it is still 0". The read no longer touches it (M120 took it out of the dump). Its early users are `install_surface_z_settings` from `set_run_plan_envelope_z` (the EOS mixture, with 1/μ_ion and the number fractions also 0) and kap_init.
     - The Run-and-card-bookkeeping banner says the members "armed from a control" are armed by the read and by begin_monte_carlo_run. `settling_active` is armed from a control by begin_kind_card.
     - README:11 says the type-and-shape check covers "every copy of a control into a star member". It covers run copies only; §6 says it right.
     - The star-field recipe's per-card latch advice does not say that the comment must name the control and run_driver, which check_control_copies requires.
   - Fix: comments and docs only.
   - A, S.

8. **The solar calibration's saved restores duplicate the run-copy re-arm.**
   - Where: core/run/run_driver.f90:265-266, :298-299; star_info_lib.f90:383-384, :388-389.
   - The issue:
     - begin_calibration saves `atm_choice_in_use` and `structure_dt_limits_active` right after begin_monte_carlo_run has re-armed them, so the saved values are always the controls'.
     - The comments' "the value its Monte-Carlo run started with" is the deck's value.
   - Fix: restore from `s%physics%atm_choice` and `s%physics%use_structure_dt_limits` (run_driver is a named writer, so check_control_copies allows it); drop the two `saved_*` locals. Combine with 1.
   - A, S.

9. **begin_monte_carlo_run re-arms and re-derives the run range every run.**
   - Where: core/run/monte_carlo.f90:90-92, :108; controls_sync.f90:919-920.
   - The issue:
     - `first_mc_run` and `last_mc_run` are run copies whose only change is `set_monte_carlo_range`'s derivation from three controls fixed for the call (`monte_carlo_active`, `mc_run_start`/`mc_run_end`, `max_mc_runs`).
     - Each run undoes the clamp and redoes it, and a comment has to explain that the loop's bounds are fixed.
     - This is the declined item 3's pattern.
   - Fix: a function `monte_carlo_range(s)` used by run_yrec and setup_monte_carlo_runs; the two run copies dropped.
   - A, S. Low priority.

## What would raise each score (concrete)

**state/**
- **Readability 8.5 -> 9:**
  - break README line 7 into a list: what a control is, what a run copy is, and the run copies by class (pointing to §2);
  - turn line 11 into one bullet per check, each pointing to §6;
  - keep the §6/§7 bullets under about 600 characters;
  - new 5 removes the "constants … each run sets" wording.
- **Organization 8 -> 9:**
  - new 5 (one class per Monte-Carlo-run lifetime);
  - new 1 (a latch derived from its run copy, or reset with it);
  - one sentence on why run_messages is in state/, or a directory of its own for the shared services;
  - Phase 5's typed parts are what takes it past 9.
- **Maintainability 8.5 -> 9:**
  - new 4 and new 6;
  - M169's proposed test: two identical samples, diff run 2 against run 1. It would also have caught new 1;
  - extend REQUIRED_DIRS to core/model, core/structure, core/envelope, core/diagnostics and rotation: 57 procedures, two `none` lines each, then `--fix`.
- **Extensibility 8.5 -> 9:**
  - the per-Monte-Carlo-run recipe pointing to the one class (new 5);
  - a line in "Add a value the run changes" for a latch tied to a run copy (reset where the copy is restored);
  - the star-field recipe's per-card latch step naming the control and run_driver (new 7).
- **Docs 8 -> 9:**
  - new 2, new 4, new 7;
  - M161 (state/README:18 and reaction_layout's header still name net_configuration and setup's map_user_inputs).

**setup/**
- **Readability 9.5 -> 10:** move the map_user_inputs row (io's code) into io/README, leaving a one-line pointer, and shorten §1's table cells.
- **Organization 9 -> 9.5:**
  - M85/M96's shared recompute as a settings_install procedure (the surface mixture for a new Z, convection's metals included);
  - state as a rule that `install_surface_z_settings` may run before `eos_init` (the run-plan Z), or pass that Z to kap_init alone.
- **Maintainability 8.5 -> 9:**
  - new 6;
  - the control-to-field check (setup/README §6's registry column).
- **Extensibility 9 -> 9.5:**
  - the same registry column, which removes one hand edit from "Add a field";
  - the paths records under check_settings_fields, so the table-path recipe's builder step is checked.
- **Docs 8 -> 9:** new 3. With it, every statement I traced in setup/README holds.

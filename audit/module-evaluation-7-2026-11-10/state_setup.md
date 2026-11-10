# state/ + setup/ -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

**Scope.**
- Read whole: every hand-written .f90 under state/ and setup/; both READMEs (the generated tables of state/README §2-3 checked row by row where a finding needed it); tools/check_star_contracts.py; and check_boundaries.py's star and settings rules (install sites, `check_settings_fields`, `unset_slots`, `check_settings_lifetimes`, the controls checks, `check_max_runs`/`literal_run_sizes`, `check_named_rows`, `check_control_copies`).
- Read in part: gen_star_docs.py (its header, the class and comment parser, and the write/read analysis).
- Read for the drivers' use of the star and the settings: monte_carlo.f90, run_driver.f90, run_setup.f90 and yrec_reset.f90 whole; the relevant parts of calibration.f90, stop_conditions.f90, starting_model's update_surface_mixture, mid_timestep_model's rotate_substep_model, omega_profile's solid_body_omega, mass_change and envelope_rebuild.
- Context: ROADMAP's "Sixth-evaluation pass", bug_sweep §13 (M181-M195), the sixth report, GUIDELINES (Naming, the six star rules, Declarations, Code style) and TARGET_ARCHITECTURE 4.2/4.6.

**Checks.** Each was run as `python3 -I -B`, with the tools directory put on sys.path by hand. `git status` and `tools/__pycache__` were the same before and after.
- `check_boundaries.py`: "Domain boundaries OK".
- `check_star_contracts.py`: "Star contracts OK (149 declared)".
- `gen_star_docs.py --check`: "outputs up to date".

**Scratch scripts** (SP/review7_state_setup/). Both import check_star_contracts' scan and only read:
- `unresolved.py`: contract coverage by directory, and the star-member arguments to callees the scan cannot resolve.
- `nointent.py`: the star-member arguments passed to a dummy that declares no intent.

**Sweeps:**
- Contract coverage (procedures with a contract / procedures that take the star): core/run 26/26, core/step 25/25, core/model 4/4, core/structure 6/6, core/envelope 5/5, core/diagnostics 26/26, setup 15/15, io 26/29 (the three are controls_sync's, which are exempt), rotation 16/16. **149/152** (sixth: 132/151).
- Every call that passes a star member resolves to its callee, and every such dummy declares an intent (0 unresolved, 0 without an intent).
- `card_rescale_targets` and `rescale_params` subscripted by a literal row: 0 (sixth: 58).
- Whole-member copies of a control outside the read: 8. They are the same 8 as at the sixth evaluation, each named in its member's comment.
- Literal subscripts on star vectors that have no slot names (new 4): 91. Literal species rows of `s%xa`: 2.
- Non-generated README lines over 600 characters: state/README :38 (776) and :46 (819), both table cells. setup/README has none. The longest code line in setup/ is 104 characters.

## Scores (1-10)

The bar is rotation/'s. Fourth -> fifth -> sixth -> now.

**state/**
- **Readability 8 -> 8.5 -> 9 -> 9.5.** The sixth report listed three items for 9.5, and all three are done:
  - the rescale targets' rows are named (`i_target_*`) at every site, and a check keeps them named;
  - the run_messages row is split, and the comments are rewrapped;
  - the Sources section cites the rename commits instead of repeating the lists.

  Left:
  - physical_constants still exports 15 cryptic legacy names, two of which are twins of other constants (new 5);
  - six star vectors still have no slot names (new 4);
  - two table cells are over 600 characters (:38, :46).
- **Organization 7.5 -> 8 -> 9 -> 9.** Changes since the sixth evaluation:
  - the Monte-Carlo range is now a pair of pure functions of the controls;
  - the solar luminosity's exception to the run-copy mechanism is stated;
  - so is the home of "state a run leaves that is not a run copy".

  Not raised. Left:
  - the solar luminosity still has its own mechanism (documented, not unified; M97);
  - the per-run reset is still a hand list, and a reading finds one more member it misses (new 2);
  - "Run and card bookkeeping" mixes four lifetimes, and a new check now relies on the class (new 3);
  - Phase 5's typed parts.
- **Maintainability 8 -> 8.5 -> 9 -> 9.** Changes since the sixth evaluation:
  - a literal run-card size now fails `check_max_runs`;
  - rotation's 16 procedures declare their contracts, and rotation/ is required;
  - `check_named_rows` refuses a literal rescale row.

  Not raised: the sixth report's other two items stayed open (M190's fix and a two-run disk/rotation test; M176). Left:
  - both star tools count an internal WRITE into a star member as a read. One contract and one README row are wrong, and both tools pass (new 1);
  - `check_settings_lifetimes` judges change frequency by lifetime class, which is too coarse (new 3);
  - one more member carries over from run to run (new 2).
- **Extensibility 8 -> 8.5 -> 9 -> 9.** A new rescale target is now a named row: `n_rescale_targets` and `i_target_*` in control_slots, spelled through to the buffer and the registry. Not raised: the snapshot that would end the hand list was measured but not built (M190's option 2). Left:
  - the per-Monte-Carlo-run recipe leaves out test_controls' regeneration for a run copy (new 7);
  - "Add a value the run changes", step 4, can lead to a builder that holds a stale value, and the lifetime check passes it (new 3).
- **Docs accuracy 7.5 -> 8 -> 8.5 -> 9.** The sixth report's docs items are all done:
  - new 2: contracts cover their own body (README :101, GUIDELINES rule 4);
  - new 4: the four stale statements;
  - new 1: the reset is worded as "the members found so far".

  Left (small):
  - the bookkeeping banner says "Most are" run copies; 7 of 31 are (new 6);
  - `termination_reason`'s generated row and `check_stop_conditions`' contract are wrong (new 1);
  - a non-member is listed among the members' F77 homes, and three small stale statements (new 7).

Average **9.1** (sixth: 8.9).

**setup/**
- **Readability 9 -> 9.5 -> 9.5 -> 10.** Both items the sixth report listed for 10 are done:
  - §3's "Every field set" and §5's long steps are now sub-bullets;
  - settings_install's header is rewrapped.

  The builders are one copy per line under a contract, and the README is a table, the order of the calls, the rules and the recipes. Nothing material is left.
- **Organization 8 -> 9 -> 9 -> 9.** Unchanged. Left:
  - M85/M96's shared recompute (open);
  - install_run_settings hands two star members straight to kap_init, outside any builder. This is documented at the member, but no check sees it (new 3).
- **Maintainability 8 -> 8.5 -> 9 -> 9.5.** Both sixth-report items are done:
  - an array field counts as set only when the builder sets it whole or every slot (M193);
  - a builder is installed at least as often as the members it copies change (check_settings_lifetimes).

  Left:
  - the lifetime check judges a member's change by its class, and two classes hold members that change within a card (new 3);
  - the rotation-law clamp `[-2, 0]` is written as literals (rotation_settings_builder.f90:99).
- **Extensibility 8 -> 9 -> 9.5 -> 9.5.** "Add a field" now names `check_settings_lifetimes`. Not raised to 10: the same gap (new 3) lets a per-card builder copy a latch the run switches within a card, and it passes every check.
- **Docs accuracy 8 -> 8 -> 9 -> 9.5.** The sixth report's new 8 is done (§3 lists `solid_body_mode_flag = 1`; "Add a field" starts with its control). Left:
  - §3's lifetime rule implies that the bookkeeping class changes at most per card (new 3);
  - rotation_configuration.f90:147-149 says that the controls read sets `solid_body_mode_flag`. Setup's builder does (new 7).

Average **9.5** (sixth: 9.2).

## Summary

- **The pass did the sixth report's lists, and I checked each item in the code.**
  - Named rescale rows, with a check that refuses a literal row.
  - check_max_runs refuses a literal size (M192).
  - Array fields count as set only when every slot is (M193).
  - A lifetime check for installs.
  - Rotation's contracts are required (149/152).
  - The contract text is fixed (M194), and so are the stale statements (M195).
  - Monte Carlo's range is a pair of functions.
  - The solar luminosity's exception and the per-run reset's limits are documented.

  The three behaviour questions (M190, M191, and io's M181) are the authors' decisions, with measured options.
- **The three tools pass, and the tree was unchanged after running them.**
- **New 1 is the most useful new item: both star tools count an internal WRITE as a read.**
  - `check_stop_conditions` writes `s%termination_reason` with `write(s%termination_reason, ...)`.
  - Its contract lists the member under reads and says `star writes: none`.
  - README §3.5's row leaves out the writer, and names stop_conditions as a reader, which it is not.
  - `--fix` would not change any of this.
- **New 2: a reading finds one more member the per-run reset misses.**
  - `omega_substep` is the first guess of the first solid-body sub-step.
  - A later Monte-Carlo run, or a card that re-reads its starting model, starts that sub-step from the previous run's last ω instead of the template's 0.
  - The effect is at the 1e-14 Newton tolerance. The F77 had it too, because GETW SAVEd OMEGAM.
  - It supports building the snapshot rather than extending the hand list.
- **New 3: `check_settings_lifetimes` takes a member's lifetime class as its change frequency.**
  - "Run and card bookkeeping" holds latches the run switches within a card (`atm_choice_in_use`, `structure_dt_limits_active`, `mass_change_active`).
  - It also holds the surface mixture, written after the per-card install.
  - A per-card builder that copies any of these passes the check.
- **The other new items are about structure and wording:**
  - six star vectors without slot names (91 literal subscripts), and two literal `xa` species rows;
  - physical_constants' legacy names;
  - "Most are run copies" (7 of 31);
  - stale lines.
- **Still weighing (bug_sweep or ROADMAP):** M85, M96, M97, M115, M126, M176, M190, M191; the solver-diagnostics gate and the read's buffer as module state (TARGET_ARCHITECTURE 4.6); Phase 5.

## Previous report's items: status

New findings 1-9 of the sixth report:

| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | `disk_gate_age_gyr` stale across Monte-Carlo runs and cards | PARTLY | M190 is open as the authors' decision, and its options are measured (ROADMAP "Measured for the authors' decisions"). The wording is done: `clear_previous_run_state`'s header (monte_carlo.f90:125-150), the member's comment (star_info_lib.f90:883-889), README :26 and §7 :571-575 now say "the members found so far". No two-run accreting, disk-locked test was added. |
| 2 | contracts called a promise of "everything the procedure may change" | DONE | M194. README :101 and GUIDELINES rule 4 now say that the lists cover the body, not the reach, with begin_monte_carlo_run as the example. The `star calls:` line was not added, by decision. |
| 3 | `card_rescale_targets` rows unnamed (58 literal subscripts) | DONE | control_slots.f90:72-82 (`i_target_*`, `n_rescale_targets`). check_boundaries.py:836-857 (`check_named_rows`). 2b2f7975 declares `rescale_params(n_rescale_targets,max_runs)` in the buffer and the registry. 0 literal rows. |
| 4 | four stale statements | DONE | M195: README §6 :507, run_messages.f90:26-29, control_slots.f90:33-36 and :54-57. |
| 5 | check_max_runs blind to a literal size | DONE | M192: `literal_run_sizes` (check_boundaries.py:792-807), with mutation cases. |
| 6 | Monte Carlo sets cards 2 and 3's stop ages without the calibration | PARTLY | M191 is open as the authors' decision. Documented at monte_carlo.f90:201-211 and README §7 :585. |
| 7 | readability leftovers | DONE | The run_messages row is short (README :47), with its state as a list (:50-55). The banners are rewrapped. settings_install's header reads well. Two other cells are now over 600 characters (:38, :46). |
| 8 | setup doc gaps | DONE | §3 :71 lists `solid_body_mode_flag = 1`. "Add a field" starts with its control (:111). |
| 9 | solar luminosity outside the run-copy mechanism | PARTLY | It is named as "One exception" (README :22) and in the per-Monte-Carlo-run recipe (:584). It is still a second mechanism (M97 open). |

What would raise each score:

| Criterion | Item | Status |
|---|---|---|
| state Readability | new 3 (named rows) | DONE |
| state Readability | new 7 (split the run_messages row; rewrap) | DONE |
| state Readability | Sources drops the rename lists | DONE (star_info_lib.f90:1298-1305 cites b3a6675b, a568e084, 0197827c) |
| state Organization | first_mc_run/last_mc_run as functions of the controls | DONE (0ec5ed1d; monte_carlo.f90:48-60; README §4 :480) |
| state Organization | new 9 (one mechanism for values the run changes) | PARTLY (the exception is documented, not removed) |
| state Organization | a stated home for per-run state that is not a run copy | DONE as text (README §7 :571-575). The snapshot was measured, not built. |
| state Organization | Phase 5 typed parts | OPEN |
| state Maintainability | new 1's fix and a two-run disk/rotation test | OPEN (M190, the authors' decision) |
| state Maintainability | new 5 | DONE |
| state Maintainability | M176 decided | OPEN |
| state Maintainability | rotation's 16 contracts | DONE (b621c57e; REQUIRED_DIRS includes rotation, check_star_contracts.py:461) |
| state Extensibility | new 1's snapshot | OPEN (measured as M190's option 2) |
| state Extensibility | new 3 | DONE |
| state Docs | new 2, new 4, new 1's wording | DONE |
| setup Readability | sub-bullets for §3 and §5 step 5; rewrap settings_install :12 | DONE |
| setup Organization | M85/M96's shared recompute | OPEN |
| setup Maintainability | array field set only when every slot is | DONE (M193; `unset_slots`, check_boundaries.py:403-427) |
| setup Maintainability | lifetime check for a builder that copies a per-Monte-Carlo-run member | DONE (`check_settings_lifetimes`, :484-528), with a granularity gap (new 3) |
| setup Extensibility | the lifetime check in "Add a field" | DONE (setup/README :113) |
| setup Docs | new 8 | DONE |

## New findings (not in bug_sweep)

1. **tools/check_star_contracts.py:305-334 (`scan`) and tools/gen_star_docs.py:340-353: an internal WRITE into a star member counts as a read. Where it happens: core/run/stop_conditions.f90:132-136 and :145; state/README.md:187.**
   - **Issue:**
     - Both scans treat a `read(` statement's list items as written. They have no case for `write(`, so `write(s%termination_reason,'(a,es10.3,a)') ...` falls through to "every `s%` here is a read". This is the end-age stop's message, in `check_stop_conditions`.
     - So `check_stop_conditions` declares `termination_reason` under `star reads:` and `star writes: none`, and check_star_contracts passes it (`--fix` would keep it).
     - gen_star_docs' row for `termination_reason` (README :187) leaves `check_stop_conditions` out of "written by" and lists `stop_conditions` under "read by". No procedure of stop_conditions reads the member: every other reference there is an assignment (stop_conditions.f90:118, 221, 229, 268, 290).
     - The same scan feeds `check_controls_writer` and `check_control_homes`. An internal write into a character control (`write(s%control%x, ...)`) would therefore pass the rule that only controls_sync writes the records.
     - One site in the tree.
   - **Fix:**
     - In `scan` and in gen_star_docs' statement analysis, count an internal file unit as written: the first item of a `write(` control list when it is `s%<member>`.
     - Add a self-test case (`write(s%msg,'(a)') 'x'` is a write).
     - Run `--fix` and `gen_star_docs.py`: the contract gains `star writes: termination_reason`, and the row gains the writer and loses the reader.
   - **Kind** A. **Size** S.

2. **state/star_info_lib.f90:902-907 (`omega_substep`); rotation/mid_timestep_model.f90:366-375; rotation/omega_profile.f90:175-189; core/run/monte_carlo.f90:146-150: a member the per-run reset misses.**
   - **Issue:**
     - In the solid-body branch of `rotate_substep_model` (`solid_body_everywhere()`: `force_solid_body_rotation`, or `solid_body_mode_flag = 1`), `solid_body_omega` takes the mean of `s%omega_substep` as its first guess.
     - That is the previous step's last sub-step ω, as the member's comment says ("it outlives the step").
     - Nothing writes `omega_substep` at a model read, at a card start or in `clear_previous_run_state`.
     - So the first step of every run's first card starts the Newton iteration from the template's 0. The first step of Monte-Carlo run 2 or later, and of a card with `card_from_starting_model` after an evolved card, starts from the previous run's or card's last ω.
     - `solid_body_omega` stops within `default_inertia_tolerance` = 1e-14 of J. It then sets j from the final ω (omega_profile.f90:229-232), so the two starts give j_rot that differ at about that level, and the runs need not be bit-identical.
     - Whether the difference reaches a printed history digit was not measured (read-only).
     - Latent: no shipped deck makes two Monte-Carlo runs, and test_monte_carlo's deck does not rotate.
     - README §7 :573-575 lists the known members (the deuterium rates, the report flag, `disk_gate_age_gyr`, M176) and not this one.
     - The solar calibration's printed "Old BL/Old RL" (`log_l_lsun_prev_check`, `log_r_rsun_prev_check`; calibration.f90:132-135) likewise show the previous Monte-Carlo run's last check at a later run's first check. This is terminal output only.
   - **In F77: yes.** GETW keeps OMEGAM as a SAVEd local (`git show 6cd56739:src/getw.f`, lines 27 and 83). MIDMOD passes it to SOLID in the solid-body branch (midmod.f:455-458). SOLID makes "A FIRST GUESS AT OMEGA BY AVERAGING THE PREVIOUS VALUES" (solid.f:24-36).
   - **Fix:**
     - Set `omega_substep = omega` at every model read (with M190's option 1, at the same place), or zero it in `clear_previous_run_state`.
     - Add it to README §7's list.
     - Better: build M190's snapshot, which ends the hand list for runs.
     - A rotating, solid-body two-run variant of test_monte_carlo would cover both.
   - **Kind** B (round-off level; measured commit). **Size** S.

3. **tools/check_boundaries.py:484-528 (`INSTALL_LIFETIMES`, `check_settings_lifetimes`); setup/README.md:77-83; setup/settings_install.f90:72-73: the lifetime check equates a member's class with how often it changes.**
   - **Issue:**
     - `install_card_settings` may copy any member of "Run plan" or "Run and card bookkeeping", and §3 :79 states the rule as "per card".
     - But the bookkeeping class holds latches the run switches within a card:
       - `atm_choice_in_use`: surface_boundary.f90:145-168 and envelope_integration.f90:205;
       - `structure_dt_limits_active`: step_driver.f90:230;
       - `mass_change_active`: mass_change.f90:130.
     - It also holds the surface mixture written at the model read, after the per-card install. This is why `install_surface_mixture_settings` re-passes `envelope_metals_no_cno` by hand.
     - The run plan's stops are disarmed by `init_stop_conditions` after the install.
     - So a builder installed per card that copies any of these (say, after "Add a value the run changes", step 4) passes the check and holds a stale value for the rest of the card.
     - The check also reads only the `*_from_star` bodies. `install_run_settings` passes `s%envelope_hydrogen_fraction` and `s%envelope_metal_fraction` (bookkeeping members) straight to `kap_init` once per run, which the rule would refuse. That case is documented at the member, but no check sees it.
     - Latent: no builder copies a latch today.
   - **Fix:**
     - Give the check member-level frequencies: a set of members that change within a card (the latches, the stops, `termination_reason`, the calibration loop's state) that no install may copy, and the `envelope_*` surface mixture as "per surface mixture".
     - Better still: split the bookkeeping class into card-start values and run latches, which also answers new 6.
     - Scan the install procedures' own `s%` actual arguments too.
     - Restate §3's per-card bullet.
   - **Kind** A. **Size** S.

4. **state/star_info_lib.f90:578, 583, 975, 987, 1229, 1232 (declarations); core/structure/henyey_iteration.f90, henyey_elimination.f90, core/model/rezoning.f90:225-233; core/step/mass_change.f90:264-265: star vectors without slot names.**
   - **Issue:** GUIDELINES Code style rule 4 ("slots and species are named") is followed by every other star vector (`i_eps_*`, `i_lum_*`, `i_nu_*`, `i_mech_*`, and now `i_target_*`). These six are not, and are subscripted by literals at 91 sites:
     - `elim_rhs` (33);
     - `max_residual` (30);
     - `envelope_fit_coeffs` (13: c(1:9), log R, log L and log Teff, each in log P, log T and 1);
     - `surface_bc` (6);
     - `max_correction_index` (5);
     - `stored_envelope_state` (4: log P, log T, log r, log m/M).
     - Separately, mass_change.f90:264-265 reads `s%xa(1,nz)` and `s%xa(3,nz)`, where `i_h1` and `i_metals` exist. The comment at star_info_lib.f90:687 also writes `physics%timestep_factors(1)` for `itime_core_min`.
     - `NAMED_ROWS` (check_boundaries.py:836-839) is built to take such arrays once their literal rows are gone.
   - **Fix:**
     - Equation slots (`i_eq_p`, `i_eq_t`, `i_eq_r`, `i_eq_l`) for `elim_rhs`, `max_residual` and `max_correction_index`.
     - Fit-coefficient and stored-state slots, beside `i_target_*` in a constants-only module.
     - `i_h1` and `i_metals` in mass_change.
     - Then add the arrays to NAMED_ROWS. Byte-identical. Most sites are core/structure's.
   - **Kind** A. **Size** M.

5. **state/physical_constants.f90:17-28, 32-58: legacy names for the shared constants.**
   - **Issue:**
     - The module exports `cc13` (1/3), `cc23`, `clndp` (a second name for `ln10`), `clni`, `cpi`, `c4pi`, `c4pil` (log10 4π), `c4pi3l`, `csig`, `csigl`, `cmkh`, `cdelrl`, `cmixl2` (equal to `cc13`), `cmixl3` and `cgl`. They are used in 2-19 files each.
     - The header has to warn that `cmixl2`/`cmixl3` "despite the names" are not powers of the mixing length.
     - Two of them are twins of other constants, against Code style rules 2 and 4 and TARGET_ARCHITECTURE 4.7 ("no name is reused", readable names).
     - `cpi` = 3.1415926535898d0 is a 14-digit literal, about 2e-15 off the double-precision π, and nothing says that it is kept on purpose.
     - The newer names in the same module read well (`seconds_per_year`, `radiation_constant_over_3`, `amu_cgs`).
   - **Fix:**
     - Readable names (`one_third`, `log10_four_pi`, `log10_sigma_sb`, `log10_saha_constant`, `log10_g_run`, ...), with the legacy names in the comments and GLOSSARY.
     - Drop `clndp` and `cmixl2` for `ln10` and `one_third`.
     - Say why `cpi` keeps its literal.
     - Byte-identical. The renames touch the domains' files.
   - **Kind** A. **Size** M.

6. **state/star_info_lib.f90:349-351 (and the generated README :78): "Most are the run's working copies of controls it changes (run copies)".**
   - **Issue:** 7 of the class's 31 members are run copies: `atm_choice_in_use`, `structure_dt_limits_active`, `mass_change_active`, `core_mass_reduction`, `num_core_points_added`, and the two calibration targets. Most of the class is card-start values, the first card's surface mixture and the calibration loops' state. The sentence came in a43d1efa (2026-11-01).
   - **Fix:** "Seven are run copies of controls (...); the rest are set at a card start, by the first card's model read, or by the calibrations". Or split the class (new 3).
   - **Kind** A. **Size** S.

7. **Small stale statements:**
   - **Locations:**
     - state/star_info_lib.f90:1285;
     - state/README.md:507;
     - core/run/monte_carlo.f90:265-267, 297, 301;
     - rotation/rotation_configuration.f90:147-149;
     - state/README.md:38 and :46;
     - state/README.md:583 and :588.
   - **Issue:**
     - **The Sources section** lists `zahn_coupling (VFC)` among "the members' F77 homes". It is no longer a member: it is a field of rotation's `circulation_velocity_state` (README §4).
     - **README §6** gives `log10_newton_g`'s caller as "(`run_setup`, every run)". Convection's and rotation's builders also call it, every card. It gives `log10_radius_from_l_teff`'s as "(`observables`, every model)". run_driver, mass_change and rotation's magnetic_braking also call it.
     - **write_monte_carlo_model's comments** name `monte_carlo_unit1`/`monte_carlo_unit2`. The code's units are `mc_summary_unit`/`mc_model_unit`.
     - **rotation_configuration's comment** says that `solid_body_mode_flag = 1` is set by "the controls read". Setup's `rotation_settings_from_star` sets it (rotation_settings_builder.f90:100; derived_controls.f90:160-163 says so).
     - **Two state README table cells** are 776 and 819 characters long.
     - **The per-Monte-Carlo-run recipe:** step 1 sends a run copy to "steps 1-3" of the run-copy recipe. Its step 5 leaves out `YREC_REGEN_CONTROLS=1 python3 -m pytest -q test_controls.py`, which the new run copy's dump line needs (io/test/controls_dump.inc writes every run copy; run-copy recipe, step 6). A newcomer finds test_controls failing.
   - **Fix:** comments and docs only. Split the two cells into a short cell and a list below the table. Send step 1 to "steps 1-3 and 6", or add test_controls to step 5.
   - **Kind** A. **Size** S.

## What would raise each score

**state/**
- **Readability 9.5 -> 10:**
  - new 5 (readable constant names, no twins);
  - new 4 (slot names for the six vectors);
  - the two long cells (new 7).
- **Organization 9 -> 9.5:**
  - M190's snapshot (one mechanism for runs, so no hand list), with new 2;
  - split "Run and card bookkeeping" by lifetime: card-start values, run latches, the first card's surface mixture, calibration loop state (new 3, new 6);
  - the solar luminosity as a registered run copy, with M97.
- **Organization 9.5 -> 10:** Phase 5's typed parts.
- **Maintainability 9 -> 9.5:**
  - new 1 (internal writes in both tools, with a self-test case);
  - new 2;
  - M190 decided, with a two-run rotating or disk-locked test;
  - M176 decided.
- **Extensibility 9 -> 9.5:**
  - the snapshot (a new per-run member needs no bisection);
  - new 3's member-level lifetimes, so that the run-copy recipe's step 4 cannot produce a stale per-card copy;
  - new 7's recipe step.
- **Docs 9 -> 9.5:**
  - new 6;
  - new 7;
  - new 1's contract and row, which follow from the tool fix.

**setup/**
- **Readability:** 10; nothing material left.
- **Organization 9 -> 9.5:** M85/M96's shared recompute as one settings_install procedure.
- **Maintainability 9.5 -> 10:**
  - new 3 (member-level change frequency, and the installs' own star arguments scanned);
  - named parameters for the rotation-law clamp (rotation_settings_builder.f90:99), in rotation_configuration beside the law's branches.
- **Extensibility 9.5 -> 10:** new 3.
- **Docs 9.5 -> 10:** §3's per-card bullet restated (new 3); rotation_configuration's comment (new 7).

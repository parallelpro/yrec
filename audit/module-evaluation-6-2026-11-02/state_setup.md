# state/ + setup/ -- sixth evaluation (2026-11-02, yrec-next 2632f06c, read-only)

What I read whole: every .f90 and .inc under src/state/ and src/setup/ (the two generated .inc files for their approach only), both READMEs (the generated tables of state/README §2-3 for their approach only), src/tools/gen_star_docs.py, src/tools/check_star_contracts.py, and check_boundaries.py's state and setup rules (lines 1-120, 255-1000, and the main's list of checks). For context I read the fifth report, ROADMAP's "Cleanup round" section, bug_sweep §1 and §11-12, GUIDELINES' six rules, Naming, Declarations and Code style, rotation/README.md, and the core code the rules reach: monte_carlo.f90 (whole), run_driver.f90 (whole), run_setup.f90 (whole), yrec_reset.f90, surface_boundary's gray switch, envelope_integration's fallback, starting_model's update_surface_mixture, calibration's star-calibration setup and check, the controls_sync arm_run_copies, controls_registry.py's builder_copies/named_sizes/type_of, and the domains' message modules.

Checks run (`python3 -B`, `git status` unchanged before and after):
- `gen_star_docs.py --check`: "outputs up to date".
- `check_star_contracts.py`: "Star contracts OK (132 declared)".
- A scratch script (SP/review6_state/coverage.py) importing check_star_contracts' scan, read-only.

Sweeps:
- Contract coverage (star-taking procedures with a contract / all): core/run 25/25, core/step 25/25, core/model 4/4, core/structure 6/6, core/envelope 5/5, core/diagnostics 26/26, setup 15/15, io 26/29 (the three are controls_sync's, exempt), rotation 0/16. **132/151** (fifth: 92/149).
- Whole-member copies of a control outside the read: 8 statements (run_setup's three solar units, run_driver's two calibration restores and settling_active, surface_boundary's two switch-backs); each member's comment names the module and the control.
- 25 calls pass the whole star from a procedure whose contract says `star writes: none` to a callee that writes the star (e.g. begin_monte_carlo_run -> arm_run_copies, 33 members; read_controls -> store_controls_to_star, 301).
- No literal run-card size left in io/controls_buffer.f90 (34 `max_runs`); the registry has no `(50)` cell.
- `card_rescale_targets(<1-4>, ...)` literal-row subscripts: 58 (model_rescaling 28, calibration 25, monte_carlo 2, run_setup 2, plus one comment each in state/ and io/).
- Non-generated README lines over 600 characters: state/README :37 (755), :45 (618), :46 (1,357), :480 (833); setup/README :62, :75, :78, :107, :111 (651-905). None in state §6/§7.

## Scores (1-10)

The bar is rotation/'s. fourth -> fifth -> now.

**state/**
- Readability: 8 -> 8.5 -> **9**. The README's opening is now labelled bullets (controls and run copies, who arms what, where the checks are), and §6/§7 bullets are all under 600 characters. Left: the run_messages row is one 1,357-character table cell; ragged comment wraps in star_info_lib's header and Monte-Carlo banner; card_rescale_targets' four rows have no slot names (new 3).
- Organization: 7.5 -> 8 -> **9**. One class per Monte-Carlo-run lifetime, a derived gray-fallback latch with its rule, run_messages' place explained, calibration restores from the controls. Left: first_mc_run/last_mc_run as run copies (open); the solar luminosity is a run copy outside the run-copy mechanism (new 9); the flat 264-member type (Phase 5).
- Maintainability: 8 -> 8.5 -> **9**. check_max_runs now checks every registry cell and generated size; check_settings_fields covers the paths records and requires `pure`; all of core/, setup/ and io/ must declare contracts; test_monte_carlo.py is in CI. Left: the per-run reset is a list that a bisection found, with a one-deck test, and a reading finds another member it misses, with a larger effect (new 1); check_max_runs cannot see a literal size in the hand-written buffer (new 5); M176 open.
- Extensibility: 8 -> 8.5 -> **9**. All six recipes and the latch rule are traced. The per-Monte-Carlo-run recipe points to the one class, and the per-card latch step names the control and run_driver. Left: the "state a run leaves" step depends on finding the member by bisection; adding a rescale target means adding another literal row.
- Docs accuracy: 7.5 -> 8 -> **8.5**. The fifth's new 2, 4 and 7 and M161 are fixed. Left: the round's own rewrite added two false sentences: README §3 says a contract "names everything the procedure may change" (new 2), and §6 calls run_messages "the one module here with procedures". There are also three smaller stale statements (new 4), and the claim that clear_previous_run_state clears "the state the previous run left that the next would read" (new 1).

**setup/**
- Readability: 9 -> 9.5 -> **9.5**. The map_user_inputs row is now a one-line pointer to io/README, and the §1 cells are short. Left: §3 rules and §5 recipe steps of 650-905 characters (e.g. :111); a ragged line in settings_install's header (:12).
- Organization: 8 -> 9 -> **9**. The rule that install_surface_z_settings may run before eos_init is stated (§3). Left: M85/M96's shared recompute (open).
- Maintainability: 8 -> 8.5 -> **9**. Purity is enforced, the paths records are checked, and the registry's settings column with check_settings_column checks both ends of every copy. Left: an array field counts as set when any one element is (rotation's velocity_scale/mixing_scale; test_rotation_kernels covers that one); nothing ties a builder's inputs' lifetime to its install's.
- Extensibility: 8 -> 9 -> **9.5**. The paths records and the settings column are in the field, domain and table-path recipes, each step naming the check that fails until it is done.
- Docs accuracy: 8 -> 8 -> **9**. The fifth's new 3 is fixed: the builders' star inputs, the read's steps and map_user_inputs' owner are all right. Left: §3's list of rotation's derivations leaves out `solid_body_mode_flag = 1` (§1 has it); "Add a field" sends you to "Add a control" at step 5, after the copy it depends on.

## Summary

- **The cleanup round did what the fifth report asked, item for item, and I checked each in the code.** All nine of the fifth's new findings but item 9 (which stays open by ROADMAP) are done, and every "what would raise" item but Phase 5 and M85/M96 is done:
  - the Monte-Carlo values sit in one class;
  - the gray latch is derived from its run copy (on_gray_fallback; equivalent within a run, checked against the envelope integration's old latch-setting);
  - the calibration restores read the controls;
  - check_max_runs and check_settings_fields check what they say;
  - every core/, setup/ and io/ procedure that takes the star has a contract (132/151);
  - M169 has a test;
  - setup/README describes the read as it is.
- **Both tools pass, and the tree was unchanged after running them.**
- **The most important new item is a member the per-run reset misses (new 1).**
  - `disk_gate_age_gyr` (legacy SAGE) is written only after a step converges. So the first step of every later Monte-Carlo run, and of a card that starts again from the starting model, takes the disk-locking verdicts from the previous run's or card's last model.
  - With accretion and disk locking on, a stale age past `disk_locking_age_gyr` gives that step a negative timestep, and the mass change is switched off for the whole run.
  - Latent: the shipped decks with a finite t_disk make one run, carry each card on from the last, and start younger than t_disk. It shows that clear_previous_run_state's list, and the README's statement of what it covers, rest on one bisection of five decks.
- **The round's doc rewrite introduced two false statements:**
  - contracts are per body, while state/README §3 calls them a promise of "everything the procedure may change". begin_monte_carlo_run declares `star writes: none` and re-arms 40-odd members;
  - §6 calls run_messages "the one module here with procedures".
- **Other new items are small:**
  - card_rescale_targets' rows are unnamed (58 literal subscripts);
  - check_max_runs is blind to a literal size in the hand-written buffer;
  - Monte Carlo overwrites cards 2 and 3's stop ages even without the solar calibration. The shipped cm_monte_carlo's `max_age(2) = 4.568D9` becomes 1e8 yr (F77 behaviour, latent there);
  - three stale comments.
- **Still weighing (in bug_sweep or ROADMAP):** M85, M96, M97, M115, M126, M176; first_mc_run/last_mc_run; rotation's 16 contracts; the solver-diagnostics gate and the read's buffer as module state (TARGET_ARCHITECTURE 4.6); Phase 5's typed parts.

## Fifth report's items: status

New findings 1-9:

| # | Item | Status | Evidence |
|---|---|---|---|
| 1 | gray latch not reset with its run copy | DONE | `use_ttau_relation` gone; `on_gray_fallback` (surface_boundary.f90:20-33) derives it; envelope_integration.f90:202-205 now sets only the run copy (it also set the latch before, 2765341c^:191-192, so within a run the two are equivalent); README §7 :553 states the rule; §4 :473 lists the member removed. M175 fixed. |
| 2 | run_messages not described | DONE | README :1, :3 (why it is in state/), :46 (row: module state, checks, tests, TARGET_ARCHITECTURE 4.6); units start at `error_unit` (run_messages.f90:42). Leftover: §6 :491 "the one module here with procedures" (physical_constants has two functions); run_messages.f90:26 "always 6" (new 4). |
| 3 | setup/README describes the pre-READ3 read | DONE | §1 :42-48 (the run copies `cross_section_scale_in_use`, `extra_cross_section_scale_in_use`; net's other scales `physics%`; map_user_inputs is io's); §2 :55 matches read_controls' calls (controls_reader.f90: read_input, map_user_inputs, expand_paths, check_controls, warn_inert_controls, derive_controls, interpret_kind_cards, store, arm); :62 `physics%envelope_store_log_temperature`. |
| 4 | control_slots header / vacuous check_max_runs | DONE | control_slots.f90:7-12, :17-26 true; check_boundaries.py:588-679 checks each per-run registry cell spells max_runs and each generated size (four files, each must yield sizes) equals its value. Gap: a literal in the hand-written buffer (new 5). |
| 5 | per-Monte-Carlo-run values in two classes | DONE | star_info_lib.f90:186-254: Constants holds the solar mass and radius only; the two scales and the solar luminosity with its logs are in "Monte-Carlo runs". |
| 6 | paths records unchecked; purity unenforced | DONE | SETTINGS_TYPE `\w+_(?:settings|table_paths)_t`, `result(settings|paths)`, and a check that the builder is `pure` (check_boundaries.py:338-399). |
| 7 | four stale comments | DONE | envelope_hydrogen_fraction (:447-453) names kap_init and set_run_plan_envelope_z's install; bookkeeping banner (:354-363) names begin_kind_card for settling_active; README :31 says "run copy" vs "any other copy" right; the star-field recipe (:527) names the control and run_driver. |
| 8 | calibration's saved restores | DONE | run_driver.f90:294-297 restores from `physics%use_structure_dt_limits` and `physics%atm_choice`; no `saved_*` locals. |
| 9 | first_mc_run/last_mc_run as run copies | OPEN | controls_sync arm_run_copies (last two lines) and set_monte_carlo_range re-derive them every run; ROADMAP "Still open" lists it. |

What would raise each score:

| Criterion | Item | Status |
|---|---|---|
| state Readability | README line 7 into a list | DONE (:12-21) |
| state Readability | line 11 into one bullet per check | DONE (:29-33) |
| state Readability | §6/§7 bullets under ~600 chars | DONE (longest in §6/§7 under 600) |
| state Readability | new 5 removes the "constants each run sets" wording | DONE |
| state Organization | new 5 | DONE |
| state Organization | new 1 | DONE |
| state Organization | why run_messages is in state/ | DONE (:3) |
| state Organization | Phase 5 typed parts | OPEN |
| state Maintainability | new 4, new 6 | DONE |
| state Maintainability | M169's identical-samples test | DONE (test_monte_carlo.py, in CI per a4d16b05; fails with the reset removed). One deck, no disk locking or rotation-specific state (new 1). |
| state Maintainability | REQUIRED_DIRS to the rest of core/ and rotation | PARTLY: all of core/ (check_star_contracts.py:459), 132/151; rotation's 16 stay optional by decision |
| state Extensibility | per-MC-run recipe pointing to one class | DONE (:558) |
| state Extensibility | latch-with-run-copy line | DONE (:553-555) |
| state Extensibility | per-card latch step naming control and run_driver | DONE (:527) |
| state Docs | new 2, 4, 7 | DONE (two leftovers, new 4) |
| state Docs | M161 | DONE (README :40, reaction_layout.f90:4-24; bug_sweep marks it fixed) |
| setup Readability | map_user_inputs row to io/README; shorter §1 cells | DONE (:48; §1 cells now ~100-450 chars) |
| setup Organization | M85/M96's shared recompute | OPEN (both open in bug_sweep) |
| setup Organization | rule for install_surface_z_settings before eos_init | DONE (§3 :77; true: eos_set_mixture builds the SCV table only once the tables are loaded, scv_eos.f90's `scv_tables_loaded` guard) |
| setup Maintainability | new 6 | DONE |
| setup Maintainability | control-to-field check | DONE (registry `settings` column, check_settings_column, controls_registry.builder_copies) |
| setup Extensibility | registry column in "Add a field" | DONE (it adds a cell to write, which the check prints) |
| setup Extensibility | paths records under check_settings_fields | DONE |
| setup Docs | new 3 | DONE |

## New findings (not in bug_sweep)

1. **core/step/step_driver.f90:557; core/step/mass_change.f90:465-471; rotation/omega_profile.f90:185-187; state/star_info_lib.f90:886-888; core/run/monte_carlo.f90:121-141** -- `disk_gate_age_gyr` goes stale across Monte-Carlo runs and across cards that read the starting model.
   - **The issue:**
     - It is the age the disk-locking tests read, legacy SAGE. Its only writer is the converged-step assignment, `s%disk_gate_age_gyr = s%age_gyr`. Nothing sets it at a model read, at a card start or in clear_previous_run_state.
     - So the first step of Monte-Carlo run 2 and later, and of a card with `card_from_starting_model` that follows an evolved card, sees the previous run's or card's last age. That age reaches:
       - the first rezone's and solve's `omega_from_j` (`disk_locking_engaged(stale age)`: the convective envelope is not held at the disk's ω when it should be);
       - `limit_mass_change_timestep`. With accretion and disk locking on and a stale age past `disk_locking_age_gyr`, it sets `dt = (t_disk - stale age)*1e9*yr < 0` and `disk_exhausted`. A negative dt is the code's sign for "do not age" (what the rest of that step does with it was not traced), and `change_mass` then switches `mass_change_active` off for the rest of the run, so that run never accretes.
     - The first step of every run sees 0 (the template's value), not the starting model's age.
     - In F77: yes. SAGE is in COMMON /DISK/ and set only after convergence (main.f:791).
     - Latent in the shipped decks:
       - those with a finite t_disk (92 sample_cases/Gyrochrones decks, three testsuite rotating decks, cm_rot_disk) make one run, carry each card on from the last, and start from deuterium-birthline models younger than t_disk;
       - the rest (run_from_start_to_dbl/*, cm_monte_carlo) set t_disk = 1e5 Gyr, which no age reaches;
       - cm_rot_base has disk locking off.
     - So test_monte_carlo (cm_monte_carlo, two runs) and M169's bisection could not see it. Not run here (read-only).
   - **Fix:**
     - set it from `age_gyr` at every model read (read_starting_model's reset after the read), which covers all three cases;
     - give test_monte_carlo a two-run accreting, disk-locked case;
     - make clear_previous_run_state's comment, the Monte-Carlo banner (star_info_lib.f90:210-213), state/README :68 and §7 :551 say "the members found so far", not "the state the previous run left that the next would read".
     - Alternatively, snapshot the star at the first begin_monte_carlo_run and restore it at each later one (the Monte-Carlo class excepted). That ends the hand list for Monte-Carlo runs, but not for cards.
   - B, S.
2. **src/state/README.md:92 (and GUIDELINES rule 4); core/run/monte_carlo.f90:106-107** -- the contracts are per body, but the README calls them a full promise.
   - **The issue:** §3 says "(the contract is the stricter promise: it names everything the procedure may change)". The checker's own rule is the procedure's body, with "calls … carry their own contracts" (check_star_contracts.py:15-17; GUIDELINES rule 4 says "own body"). A call that passes `s` whole adds nothing. Examples:
     - `begin_monte_carlo_run` declares `star writes: none` and re-arms 33 run copies, the range, the two cleared members and seven sampled values;
     - `prepare_run`, `compute_observables` and `read_controls` are the same (25 such calls in all).
   - **Fix:**
     - say "everything its own body may change; a call that passes `s` whole carries the callee's contract";
     - optionally, have `--fix` add a third line naming the star-taking callees (`!   star calls: arm_run_copies, …`) so that a driver's header shows its reach.
   - A, S.
3. **state/star_info_lib.f90:302-306 (rows), with 58 uses in core/** -- `card_rescale_targets(4,max_runs)` (and `control%rescale_params`) has no slot names.
   - **The issue:** the rows are the total mass, X, Z and the core mass. They are subscripted by literals 1-4 at 58 sites: model_rescaling 28, calibration 25, monte_carlo 2, run_setup 2 (e.g. `card_rescale_targets(3,n)` in set_run_plan_envelope_z). That breaks Code style rule 4 ("slots and species are named"), which every other vector of the star follows (`i_eps_*`, `i_lum_*`, `i_mech_*`, `itime_*`, `ichi_*`).
   - **Fix:** `i_target_mass, i_target_x, i_target_z, i_target_core_mass` in control_slots (which already holds the vector controls' slots), used everywhere; byte-identical.
   - A, M.
4. **state/README.md:491; state/run_messages.f90:26; state/control_slots.f90:31, :50** -- three stale or false statements.
   - **The issue:**
     - §6, new in the round, says "`run_messages`, the one module here with procedures"; physical_constants has two (`log10_newton_g`, `log10_radius_from_l_teff`).
     - run_messages' header says "terminal_unit -- stdout (always 6 …)", while :42 starts it at `error_unit` (the round's M112 initializer) and test_atm sets it to the log's unit.
     - control_slots names the legacy `lnewvars` ("when lnewvars is set") and "the deck's hpttol" for `use_new_rate_variables` and `zoning_tolerances`/`chi_grid_scale`.
   - **Fix:** comments and docs only.
   - A, S.
5. **tools/check_boundaries.py:603-628 (per_run_controls); tools/controls_registry.py:191-205 (type_of)** -- check_max_runs is blind to a literal run-card size in the hand-written buffer.
   - **The issue:**
     - A control counts as per-run only if its buffer dimensions spell `max_runs`, and `type_of` turns `max_runs` into its value before every type comparison.
     - So a new per-run control declared `real(dp) :: foo(50)` in io/controls_buffer.f90, with a registry cell `real(50)`, passes check_max_runs, check_rows_against_buffer and, if it has one, its run copy's type check.
     - When max_runs changes, it stays at 50, and `s%control%foo(current_card)` reads out of bounds past card 50.
     - This is the one hand-written place a size is typed, the case the fifth's item 4 asked to catch.
   - **Fix:** fail on any buffer or registry dimension equal to max_runs' value that is not spelled `max_runs` (no other named size is 50).
   - A, S.
6. **core/run/monte_carlo.f90:191-195** -- under Monte Carlo, cards 2 and 3's stop ages are overwritten whether or not the solar calibration is on.
   - **The issue:**
     - `apply_monte_carlo_parameters` sets `card_max_age_yr(2) = 1e8` and `(3) = 4.57e9` for every Monte-Carlo deck. That is not from the sample, not scaled by its age factor, and ignores `target_solar_age`.
     - Only setup_solar_calibration, when on, overwrites both again.
     - The shipped `cm_monte_carlo` (no calibration) sets `max_age(2) = 4.568D9` and runs with 1e8 yr (latent there: its card 2 stops at 10 models).
     - The per-Monte-Carlo-run recipe calls this procedure the place for "what a run takes from its sample", which these two lines are not.
     - In F77: yes (main.f:243-244).
   - **Fix:** drop the two lines (setup_solar_calibration sets both when it runs), or limit them to `control%calibrate_solar_model`. Measured commit.
   - B, S.
7. **state/README.md:46; state/star_info_lib.f90:21-29, :213-214, :425, :872; setup/settings_install.f90:12** -- readability leftovers.
   - **The issue:**
     - the run_messages row is a single 1,357-character table cell (the ROADMAP's "no README paragraph is over 1,000 characters" holds only if cells do not count);
     - star_info_lib's header and the Monte-Carlo banner (moved in the round) have ragged wraps ("Each / class opens", "they keep / the controls' values");
     - a few comment lines run past 80 columns.
   - **Fix:**
     - split the row: what it holds and who may use it in the table, its module state and tests as a short list under it;
     - rewrap.
   - A, S.
8. **setup/README.md:68-71; :97-103** -- two small setup doc gaps.
   - §3's list of derivations leaves out rotation's `solid_body_mode_flag = 1` under `force_solid_body_rotation`, which §1 :31 and the builder's header have.
   - "Add a field" puts "If the control is new, add it first" at step 5, after the copy (step 2) that needs it. "Add a control" step 4 points back into "Add a field".
   - **Fix:** add the flag to §3; move step 5 to the top as a precondition.
   - A, S.
9. **core/run/run_setup.f90:99-101; core/run/monte_carlo.f90:185-187; state/README.md:16, :560** -- the solar luminosity is a run copy outside the run-copy mechanism.
   - **The issue:**
     - It is copied from `physics%solar_luminosity_cgs` by run_setup and changed by every Monte-Carlo run. README :16's rule is that a value the run changes is the control's run copy (registry `run_copy`, armed by arm_run_copies).
     - The per-Monte-Carlo-run recipe (:560) carves out the exception, but the "Controls and run copies" list (:16-21) does not mention it.
     - So one kind of value has two mechanisms, and check_control_copies treats it as an ordinary copy.
   - **Fix:** make it a registered run copy (its logs derived where it is set), together with M97's fix (the sampled factor times the deck's value); or name the exception in :16-21.
   - A, S. Low priority.

## What would raise each score

**state/**
- **Readability 9 -> 9.5:**
  - new 3 (named rescale-target slots);
  - new 7 (split the run_messages row; rewrap);
  - the star_info_lib Sources section could drop the rename lists that README §5 and GLOSSARY repeat.
- **Organization 9 -> 9.5:**
  - first_mc_run/last_mc_run as a function of the controls (fifth's item 9);
  - new 9 (one mechanism for values the run changes);
  - a stated home for "per-run state that is not a run copy" (new 1's snapshot, or a lifetime class whose members the template resets per run).
  - Phase 5's typed parts take it to 10.
- **Maintainability 9 -> 9.5:**
  - new 1 (the fix, and a two-run test on an accreting, disk-locked, rotating deck so the reset list is tested beyond one non-rotating deck);
  - new 5;
  - M176 decided;
  - rotation's 16 contracts (two `none` lines each, then `--fix`).
- **Extensibility 9 -> 9.5:**
  - new 1's snapshot (so a new per-run member needs no bisection to be found);
  - new 3 (a new rescale target becomes a named slot, not a literal).
- **Docs 8.5 -> 9:**
  - new 2;
  - new 4;
  - new 1's wording of what clear_previous_run_state covers.

**setup/**
- **Readability 9.5 -> 10:** split §3 "Every field set" and §5 "Add a physics domain" step 5 (905 characters) into sub-bullets; rewrap settings_install's header line 12.
- **Organization 9 -> 9.5:** M85/M96's shared recompute as one settings_install procedure (the surface mixture for a new Z, convection's metals included).
- **Maintainability 9 -> 9.5:**
  - check_settings_fields counts an array field set only when every slot is (or names its slots);
  - a check that a builder copying a per-Monte-Carlo-run member is installed per Monte-Carlo run or per card, which README §7's per-MC recipe step 4 states in prose only.
- **Extensibility 9.5 -> 10:** the same lifetime check, so "Add a field" for a per-run value cannot install it at the wrong lifetime.
- **Docs 9 -> 9.5:** new 8.

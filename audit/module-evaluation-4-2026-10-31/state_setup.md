# state/ + setup/ -- fourth evaluation (2026-10-31, yrec-next e8c847ba, read-only)

Read whole: everything under src/state/ and src/setup/ (.f90, .inc, READMEs), src/tools/gen_star_docs.py, src/tools/check_star_contracts.py. Checked against: rotation/README.md, GUIDELINES (six star rules, Naming, Declarations, Code style, "The controls"), TARGET_ARCHITECTURE 4.2/4.6, ROADMAP's last sections, bug_sweep section 9. I ran `gen_star_docs.py --check` ("outputs up to date") and `check_star_contracts.py` ("Star contracts OK (64 declared)") with `python3 -B`. The tree was unchanged afterwards. I used grep for the writer and reader sweeps, and one scratch script that counts the star procedures without contracts (SP/review4_state/contracts_present.py).

What the sweeps show about the split:
- **The controls records:** nothing outside `io/controls_sync.f90` assigns `s%control`/`s%physics`: no assignment, no read target, and all 7 calls that pass a control are intent(in) uses. No STAR_FREE dir and not rotation reads either record (0 files). No `ctrl`/`job` names are left in code or docs.
- **The run copies:** `arm_run_copies` has exactly 36 lines. None of the 36 controls with a run copy is read anywhere outside `arm_run_copies`.
- **Contracts:** core/run (22) and core/step (25) declare them on every star-taking procedure, and setup on all 14.

## Scores (1-10)

The rotation/ sets are the bar. Earlier scores were on a 1-5 scale and are doubled here: first -> second -> third -> now.

**state/**
- Readability: 5 -> 7 -> 8 -> **8**.
  - The member comments are uniform and complete, and the records are now named for the inlist groups.
  - Held back by:
    - the README's 3,046-character intro paragraph;
    - star_info_lib is still all `double precision`;
    - the new `mc_samples` record has misleading names.
- Organization: 4 -> 6 -> 7 -> **7.5**.
  - Better: `job` is gone, and the Run plan class with checked run copies is coherent.
  - Left:
    - a flat 274-member type;
    - atm_choice in three star homes;
    - run copies that the run never changes;
    - derived values still split between the read and star members.
- Maintainability: 5 -> 7 -> 8 -> **8**.
  - Strong: check_control_copies, check_controls_writer/records and gen_star_docs' class rules.
  - Gaps:
    - M116 (known);
    - the "other copy" check sees only `s%m = s%<group>%x` statements outside the read;
    - contract coverage is not enforced;
    - max_runs lives in four places.
- Extensibility: 5 -> 7 -> 8 -> **8**.
  - The add-a-field, add-a-control and add-a-run-copy recipes exist and were traced.
  - The run-copy recipe for an existing control leaves out a generator, so a newcomer following it hits a failing test_generated_files.
- Docs accuracy: 5 -> 7 -> 7 -> **7.5**.
  - The generated tables are right by construction, and the prose was carried through the split.
  - About six small stale claims are left: a section cross-reference, the "re-exports" claim, a check that does not exist, a recipe step, a member count, and two banners.

**setup/**
- Readability: 8 -> 8 -> 9 -> **9**. Flat one-line copies, every builder with a contract, readable field names.
- Organization: 6 -> 7 -> 8 -> **8**.
  - There is one install module, but the domains' surface-mixture state (the EOS mixture, kap's surface tables) is installed from core at five call sites.
  - Module names follow two schemes, and `microphysics_settings` also holds convection and element diffusion.
- Maintainability: 7 -> 7 -> 8 -> **8**.
  - check_settings_fields and the INSTALL_FILES check hold.
  - 9 of the 10 builders are not `pure`, so "without side effects" is not compiler-checked.
  - No check ties a field to the right control (known).
- Extensibility: 6 -> 7 -> 8 -> **8**.
  - The domain, field and control recipes are traced and their names verified (DOMAINS, PUBLIC, DOMAIN_DIRS, YREC_LIBS_ORDER).
  - There is no lifetime "per mixture" for a domain whose state follows the surface mixture.
- Docs accuracy: 6 -> 8 -> 8 -> **8**.
  - Mostly true.
  - The list of star members the builders read leaves out net's eight rate scales.
  - The coverage figure still counts user_inputs, which has moved to io/.
  - §3's derivation rule does not name the use_two_z_tables exception.
  - net_settings_setup's header is stale (known, ROADMAP).

## Summary

- **The rules hold everywhere:**
  - only `controls_sync` writes `s%control`/`s%physics`, and no domain reads them;
  - every one of the 36 run copies is armed once, in `arm_run_copies`, and is the only reader of its control;
  - the generated tables and contracts pass.
- **The split's leftovers are small but real.** It changed state for the better, but:
  - `atm_choice_initial` now duplicates `physics%atm_choice`, and the atm builder copies the run copy where it wants the deck's choice;
  - `num_rotation_structure_passes` is a run copy whose only "change" is a per-step re-derivation from two controls;
  - the five no-input derived controls (known) exist only to give run copies a starting value.
- **The new `mc_samples` record** has a field named helium that sets the metal settling factor, a duplicate array, and "targets" that are factors.
- **setup/'s one-install-module rule covers settings records only.** `update_surface_mixture` installs convection's half of the surface mixture through settings_install and the EOS's and kap's halves directly.
- **Docs:** no new results-changing bug. The new items are documentation and checker gaps, plus four byte-identical simplifications.
- **Known open defects** still weigh on the scores: M95-M97, M115-M117, M126, and the TARGET_ARCHITECTURE 4.6 leftovers.

## New findings (not in bug_sweep), ranked

1. **atm_choice has three star homes after the split.**
   - Where: `state/star_info_lib.f90:217-219` (`atm_choice_initial`), `io/controls_reader.f90:190-192`, `core/structure/surface_boundary.f90:151-159`, `setup/microphysics_settings.f90:166`, `io/controls_buffer.f90:249-251`.
   - Issue:
     - `atm_choice_initial` is set in adopt_canonical_names from the buffer's `kttau`. It holds the same value `physics%atm_choice` now keeps for the whole run, and the control itself is read only by `arm_run_copies`.
     - `atm_settings_from_star` copies the run copy `atm_choice_in_use` for an init-time choice (which tables `atm_init` loads, `atm/atm_lib.f90:51-62`). The two are equal only because the install comes before any switch.
     - The read also resets the run latch `use_ttau_relation`, which the template star already starts at `.false.`.
     - check_control_copies did not see any of this, because the read's own files are exempt and the copy comes from a buffer name.
   - Fix:
     - surface_boundary and the atm builder read `s%physics%atm_choice`;
     - drop `atm_choice_initial`, the two lines in adopt_canonical_names and the buffer comment;
     - regenerate test_controls' dump labels and the star docs;
     - optionally extend check_control_copies to the read's buffer-to-star copies.
   - A (result-neutral) -- S.
2. **settings_install is not the one place a domain's configuration is installed.**
   - Where: `setup/settings_install.f90:4-7,99-110`; `core/model/starting_model.f90:982,1008-1013`; `core/run/run_setup.f90:68`; `core/model/model_rescaling.f90:158`; `core/model/rezoning.f90:1067`; `tools/check_boundaries.py:306-307`.
   - Issue:
     - `install_surface_mixture_settings` installs convection's `envelope_metals_no_cno`.
     - Thirty lines later, `update_surface_mixture` calls `eos_set_mixture` and `kap_update_surface_tables` itself.
     - The EOS mixture is set from three core sites in all, and kap's surface tables from two.
     - `INSTALLER_CALL_RE` does not list these calls, so the README's "every installer call ... is in settings_install.f90" holds only by definition.
   - Fix:
     - one `install_surface_mixture_settings(s, ierr)` that installs convection, eos and kap from the star's mixture, plus a rescale variant (also the natural home for M85/M96's shared recompute);
     - add both calls to `INSTALLER_CALL_RE`;
     - list "per surface mixture" as a lifetime in setup/README.
   - A -- M.
3. **`num_rotation_structure_passes` is a run copy the run never changes.**
   - Where: `state/star_info_lib.f90:674-677`, `core/step/step_driver.f90:512-514`, `io/controls_sync.f90:834`.
   - Issue:
     - Its only write outside the read is `= 1` when `.not. physics%rotation_active`, re-done at every step.
     - The value is a function of two run constants. It sits in the "Step control" class, and the registry's run_copy column then calls it "a value the run changes".
   - Fix: a local of converge_with_rotation, `merge(s%physics%num_rotation_structure_iters, 1, s%physics%rotation_active)`; clear the run_copy cell, regenerate.
   - A -- S.
4. **The Monte-Carlo sample record's names invite wrong readings.**
   - Where: `state/star_info_lib.f90:126-145`, `core/run/monte_carlo.f90:58,110-116`.
   - Issue:
     - `helium_fraction_param` ("the helium (diffusion) parameter") sets the metal settling factor fgrz, and the helium factor fgry keeps its control value.
     - `diffusion_factor` is an element-by-element copy of it.
     - `luminosity_target` and `age_target` are factors (their comments say so).
     - All of this is in a record created on 2026-10-30.
   - Fix: `settling_z_factor`, `luminosity_factor`, `age_factor`; drop `diffusion_factor`, or keep the legacy names in the comments.
   - A -- S.
5. **9 of the 10 builders are not `pure`.**
   - Where: `setup/microphysics_settings.f90:47,85,112,140,161,177,203,244`; `setup/rotation_settings_setup.f90:61`.
   - Issue:
     - setup/README §3 says the builders are "without side effects". GUIDELINES' Declarations rule 5 asks for `pure`.
     - Only `net_settings_from_star` has it.
     - Everything the others call is already pure: `log10_newton_g`, `derive_wind_law`, structure constructors.
   - Fix: mark them `pure`. The compiler then enforces the rule.
   - A -- S.
6. **state/README §6 lists a contract-coverage check that does not exist.**
   - Where: `state/README.md:453`; `tools/check_star_contracts.py:4-7,31-33`.
   - Issue:
     - The README lists "every procedure of `core/run/` and `core/step/` that takes the star declares them" among what test_boundaries checks. Nothing checks it: the checker is opt-in, and its docstring says so.
     - The claim is true today: 47 of 47 procedures, and setup 14 of 14.
     - The same docstring still says contracts cover reads "through the global `star`".
   - Fix: a `REQUIRED_DIRS` list in check_star_contracts (core/run, core/step, setup) that fails on a star-taking procedure without the two lines; update the docstring.
   - A -- S.
7. **The "Add a value the run changes" recipe leaves out steps.**
   - Where: `state/README.md:470`.
   - Issue:
     - After a control gets a `run_copy` cell, the recipe runs only gen_controls_state.py.
     - `defaults/controls_registry.md`, generated by gen_inlist.py, carries a run-copy column, so test_generated_files fails until `gen_inlist.py` runs.
     - Moving readers from the control to the member also changes their `star reads:` lines, so `check_star_contracts.py --fix` is needed too.
     - The trace was made with a new control, which io's recipe covers, not an existing one.
   - Fix: add both steps, as io/README's recipe has them.
   - A -- S.
8. **max_runs is written in four places, and the generator's reason for the literal is gone.**
   - Where: `state/control_slots.f90:16-20`; `tools/controls_registry.py:110` (`MAX_RUNS = 50`); `tools/gen_controls_state.py:65-68`; 31 registry type cells `(50)`; the literal `(50)` in `control_group_def.inc`.
   - Issue:
     - The generator still writes the value as a literal because "star_info_lib ... [has no] `use controls_buffer`".
     - Since 2026-10-30, star_info_lib has `use control_slots, only: max_runs` in scope.
     - Changing the capacity fails loudly, but takes 33 edits.
   - Fix: emit `max_runs` by name, spell the registry types `real(max_runs)`, or have controls_registry.py read the parameter from control_slots.
   - A -- S.
9. **Stale cross-references and incomplete banners.**
   - `state/star_info_lib.f90:1278-1279`: "state/README.md section 5 lists the members removed". It is section 4.
   - `state/README.md:7`: says star_info_lib "re-exports the layouts below". It re-exports three: composition_layout, neutrino_flux_layout, transport_mechanism_layout. It does not re-export reaction_layout, control_slots, envelope_points, stitched_model_state or physical_constants.
   - `docs/GUIDELINES.md:208`: "~240 members". There are 274.
   - The Constants banner's "Valid:" sentence (`star_info_lib.f90:173-174`) leaves out `settling_z_factor_in_use`, which changes per Monte-Carlo run.
   - The "Private to one module" banner (`:1144-1151`) does not mention `mc_samples`, which is run input, not state carried between calls.
   - Fix: a comments-only commit.
   - A -- S.
10. **setup/README is incomplete in three places.**
    - `setup/README.md:28` lists the star members the builders read but leaves out net's eight rate scales (`cross_section_scale` ... `n15_pa_scale`).
    - `:69`'s coverage figure (95.3% of 316 lines) includes user_inputs, now in io/.
    - `:51`'s "Any other derived quantity is computed by derive_controls ..." is contradicted by `use_two_z_tables`, which run_yrec derives (`core/run/run_driver.f90:150`) and the kap builder copies. That was a decided exception (ROADMAP round 2), but the rule does not name it.
    - Fix: complete the list, re-measure or restate the coverage, and name the exception (or compute the flag in `kap_settings_from_star` and drop the member).
    - A -- S.

## What would raise each score (concrete)

**state/**
- **Readability 8 -> 9:**
  - split README line 3 into four or five short paragraphs or bullets (what the star holds; controls vs run copies; who arms what; where the checks are), as rotation/README's opening does;
  - rename the mc_samples fields (item 4);
  - move star_info_lib to `real(dp)`, teaching gen_star_docs' `MEMBER_DECL_RE` and check_boundaries' type map the spelling first.
- **Organization 7.5 -> 8.5:**
  - drop the pseudo copies (items 1 and 3);
  - retire the five no-input derived controls (`solver_underrelax_factor`, `requested_envelope_mass`, `change_envelope_mass_flag`, `initial_envelope_x/z`) as plain run-state members with initializers;
  - make the rate scales and `envelope_store_log_temperature` derived controls with run copies, so the read writes nothing on the star but the records and `arm_run_copies` (TARGET_ARCHITECTURE 4.6 leftover);
  - gather the Monte-Carlo state (`first_mc_run`, `last_mc_run`, `settling_z_factor_in_use`, `mc_samples`) in one per-Monte-Carlo-run class, with a per-run re-arm hook (the M95 fix).
  - Beyond that, Phase 5's typed parts.
- **Maintainability 8 -> 9:**
  - fix M116;
  - enforce contract coverage (item 6);
  - extend check_control_copies to the read's buffer-to-star copies (item 1);
  - have gen_star_docs cover `mc_samples`' member comments (format, writers, never-read);
  - one source for max_runs (item 8).
- **Extensibility 8 -> 9:**
  - complete §7 (item 7);
  - a traced recipe for "a value the run changes per Monte-Carlo run", which today has no arming point.
- **Docs 7.5 -> 8.5:**
  - item 9;
  - item 6's claim;
  - state the semantics gap between the two tools for `intent(inout)`: gen_star_docs counts a write only when the callee assigns, while check_star_contracts counts every inout.

**setup/**
- **Readability 9 -> 9.5:**
  - `pure` builders (item 5);
  - drop the history parentheticals in headers ("(2026-10, Phase 3)");
  - one module-naming scheme (for example `<domain>_settings_builder`, or a single `settings_builders` module), avoiding the clash with the domains' `<domain>_settings` variables.
- **Organization 8 -> 9:**
  - item 2 (all surface-mixture state installed from settings_install);
  - the atm builder copies the control (item 1);
  - `use_two_z_tables` computed in the kap builder, or listed as the rule's exception.
- **Maintainability 8 -> 9:**
  - item 5;
  - add `eos_set_mixture` and `kap_update_surface_tables` to `INSTALLER_CALL_RE`;
  - a control-to-field mapping check (a registry column or a generated table), the known open item.
- **Extensibility 8 -> 9:** a "per surface mixture" lifetime with its install procedure, and the add-a-domain recipe saying where such state goes.
- **Docs 8 -> 9:** item 10; net_settings_setup's header (known, ROADMAP); §3 "One install module" stating its true scope until item 2 is done.

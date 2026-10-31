# io/ -- fourth evaluation (2026-10-31, yrec-next e8c847ba, read-only)

Scope: every .f90/.inc of src/io/ and src/io/test/, src/io/README.md, the controls generators (gen_controls_state, gen_inlist, gen_controls_dump, gen_default_columns), tools/controls_registry.py and defaults/controls_registry.tsv, all read whole; the relevant parts of check_boundaries.py, gen_star_docs.py, starting_model.f90 and run_driver.f90. The four generators' `--check` pass. The one experiment ran the generators on a scratch copy of src/ (SP/review4_io/sim). Paths are under src/.

## Scores (1-10)

Earlier scores (1-5 scale, doubled) -> now:

- **Readability 5 -> 7 -> 8 -> 7.5.** The split's `s%control`/`s%physics` and run-copy names (`card_max_age_yr`, `current_card_envelope_x`) read better at every use site, and the output half is at rotation's level. The read half still has F77 shapes: format labels in derived_controls, un-indented else-branches in controls_reader.f90:89-155, and multi-name declarations without per-member comments in controls_buffer. This reading also counts traps the third did not: misnamed header fields (finding 5), buffer names imported through legacy_controls (finding 8), and stale generated headers (finding 9). The code is the same as at the third evaluation; the half point comes from the new evidence, not from a regression.
- **Organization 5 -> 6 -> 7 -> 7.5.** Ownership is now clear. The two records follow the inlist groups and nothing changes them after the read: only io/controls_sync.f90 assigns them (grep and check_boundaries agree). The 36 values the run changes are run copies, armed once in the generated `arm_run_copies`. What the third evaluation said holds io back is unchanged: the triple copy (inlist name -> legacy name -> canonical), the read's state outliving the read (the run-file paths, `output_dir`, the diagnostics gate), and the physics domains importing io_units and run_log. New this time: the output's switches straddle both records (finding 3), and five no-input "controls" exist only to feed `arm_run_copies`.
- **Maintainability 4 -> 6 -> 7 -> 7.5.** Gains: five new enforced rules (17 planted violations caught), the read-twice re-arm check, all four generators fresh, test_io and test_controls mutation-confirmed, 94% coverage. Still held back by:
  - about 700 default-public module variables (controls_buffer 366, legacy_controls 334, which also re-exports the buffer and control_slots' 27 constants);
  - M110, plus the other unchecked writes (finding 6), and M112, M120;
  - CI (M125): io's two test programs, test_generated_files and three of its four generators are not in CI;
  - a generator gap that lets a malformed row through to a compile error (finding 1);
  - star contracts on only 1 of io's 35 star-taking procedures. Rotation has none either, so this is not scored against io, but GUIDELINES rule 4 asks for them.
- **Extensibility 4 -> 7 -> 8 -> 7.5.** Columns and .mod fields are at rotation's level: 1-2 files each, statically checked. A control now takes 4 hand edits and 5 tool runs plus a regen (the third counted 3 generator runs). A run copy adds 2 more edits, all checked. The table-path recipe, followed literally, gives code that does not compile, and no generator stops it (finding 1).
- **Docs 6 -> 7 -> 8 -> 7.5.** Most of the README is true:
  - every row of the §1 table, the §2 sequence step by step (against read_controls and run_yrec), the check_boundaries rules it lists, and the generated-files table all check out.
  - New inaccuracies: the records paragraph and the "same groups" sentence (finding 3), step 1's "in the group it belongs to" (controls_buffer is laid out by topic, the groups interleaved), the table-path recipe and the claim that the generators stop on disagreements (finding 1), and stale or false text in the registry and its tools (findings 2, 4, 9).
  - Known doc items still open: M109, M113, M117 and M120.

## Summary

- The split's rule holds wherever it is claimed. Nothing outside io/controls_sync.f90 assigns, reads into or passes `s%control`/`s%physics` to an intent(out/inout) dummy (tests aside). Every control is in its group's record, every run copy is armed exactly once, and rotation reads neither record.
- The README's description of how a run proceeds is accurate line by line. Its description of what each record holds is not (finding 3).
- What holds io back is still structural and unchanged since the third evaluation: three copies per value, about 700 public module variables, io state outliving the read, physics importing io. The open M98, M108-M113 and M120 remain, and CI covers none of io's own tests (M125).
- New, most important first:
  - the table-path recipe produces uncompilable code that the generators accept;
  - the registry describes seven controls by what their run copies do, five of them never set at all;
  - the records' contents are misdocumented;
  - the registry's legacy/live semantics have drifted;
  - model_file_header misnames MODEL2's disk fields, and most of the header is dead;
  - the history, profile and pulse writes are as unchecked as the .mod (M110);
  - a set of stale headers is left over from the split.
- Average 7.5, against 7.6 at the third evaluation: organization and maintainability gain half a point from the split, while readability, extensibility and docs lose half a point on fresh evidence.

## New findings (not in bug_sweep), ranked

1. **README.md:149, :157-162 (with :120); tools/gen_inlist.py:139-158, :167-169, :199-211.**
   - **Issue:** the table-path recipe has a trap. Step 2 of "Add a control" says `legacy` and `proposed_name` are both the new name, so a table path with inlist name `my_table_file` and code name `my_table_path` gets `legacy = proposed = my_table_file`. Run on a scratch copy, `gen_controls_state` and `gen_inlist` both exit 0. They emit `my_table_file` into the `&control` group (inlist_new_decl.inc), `adopt_canonical.inc` (`my_table_path = my_table_file`) and `legacy_mirror.inc`, but that name is declared nowhere. The `keep` class declares nothing and `legacy_decl` returning None is ignored, so the build fails under `implicit none`. This contradicts the README's "the generators stop when the registry and the declarations disagree". The working fill is `legacy = my_table_path`, `proposed_name = my_table_file` (class `newlocal`); I checked it on the same copy.
   - **Fix:** the recipe names the cells for each case. gen_inlist exits when a live row's legacy name has no declaration in the reader's scope.
   - A, S.
2. **defaults/controls_registry.tsv:2, :389, :391-394, :396-397; tools/controls_registry.py:23-25; tools/gen_inlist.py:339-341 (into controls_registry.md and the .defaults files).**
   - **Issue:** these docs describe what the run copies do, not the controls. "the solar calibration rewrites the later runs' entries", "rewritten by the calibration", "the run ramps it from fcorr0 and fcorri", "from senv0_array at the card start", "a rescale changes it"; derived means a member "that the run keeps (the Newton factor fcorr)"; "the read computes them (unit numbers ...) or the run keeps them". After the split the controls never change. `solver_underrelax_factor`, `requested_envelope_mass` and `change_envelope_mass_flag` are set by no code at all (controls_buffer.f90:203, :272-273): they stay 0 or .false. in the records and exist only as `arm_run_copies` sources.
   - **Fix:** each doc says what the control holds and names its run copy. Drop "unit numbers". This pairs with ROADMAP's open "the five derived controls could leave".
   - A, S.
3. **README.md:57, :91 (records paragraph), defaults/controls_registry.tsv:14, :366-388.**
   - **Issue:** "star%control (... the files, the output ...)" and "a legacy pair ... the same groups under the legacy names" are both wrong.
     - The output is split across the two records. `star_history_name`, `profile_interval`, `profile_data_prefix`, `pulse_format`, `write_profile_flag`, `write_pulse_flag`, `terminal_interval`, `compute_neutrino_fluxes`, `report_solver_diagnostics` and `log_output_file` are `&physics`. `history_columns_file`, `profile_columns_file`, `pulse_gyre_interval` and `isochrone_output_active` are `&control`. yrec_output.f90:85-86/:151/:190, history_output.f90:86-87 and profile_output.f90:40/:63 each read both records.
     - Seven controls sit in the legacy CONTROL namelist but in the inlist's `&physics`: fshort, profile_data_prefix, terminal_interval, nu_max_sun, delta_nu_sun, Teff_sun and report_solver_diagnostics (legacy_controls.f90:280-293).
   - **Fix:** say so, with a short table of the output's controls and their records. Moving groups would break shipped inlists, so this is docs only.
   - A, S.
4. **tools/controls_registry.py:11, :15-19, :28-29, :87; registry :366-367, :370, :384-388, :383; README.md:149; legacy_controls.f90:24-25 vs :326.**
   - **Issue:** the registry's semantics have drifted.
     - "live: read by both formats ... the legacy decks under legacy" is false for 8 inlist-only rows (history/profile columns files, star_history_name, profile_interval, pulse_format, write_profile/pulse_flag, compute_neutrino_fluxes). Their `legacy` cell holds the inlist name though no legacy deck reads it, which is what README step 2 tells a newcomer to do, against the column's own doc ("the name in the legacy decks, '-': none").
     - `canonical` "star%<x> for a star member" is documented and checked, but no row has it.
     - "The legacy namelists are frozen ... nothing is added here", yet solberg_hoiland_* were added on 2026-10-29 (73ca7732).
     - `report_solver_diagnostics` is `used = read`, but every `solver_diagnostics()` call during the run reads run_log's copy of it.
   - **Fix:** one meaning per column: a binding-variable column (or `legacy = -` with gen_inlist binding the code name), and a status for inlist-only rows. Drop the star% branch and correct the docs.
   - A, S/M.
5. **io/model_files.f90:52-63, :533; io/test/test_io.f90:743, :945-946.**
   - **Issue:** `model_file_header%disk_temperature` and `%disk_pressure` are MODEL2's TDISK and PDISK (putmodel2.f:64), that is the disk-locking age in Gyr and the disk omega in rad/s. The default 7.2722d-6 is `disk_omega_rad_s`'s. `lovstm` keeps its legacy spelling beside canonical neighbours. 21 of the header's 24 fields reach no caller; only test_io checks them, and read_yrec7 itself uses the two zone edges. `triangle_orientation` is among them: test_io checks it, but `reset_state_after_read` sets it to 1, against starting_model.f90:43's "come from the model file" (:351).
   - **Fix:** rename the fields to `disk_locking_age_gyr`, `disk_omega_rad_s` and `intermediate_overshoot_active`. Say in the header comment that the MODEL2-only fields are documentation of the format, and correct starting_model.f90:43.
   - A, S.
6. **io/history_output.f90:331-359, io/profile_output.f90:68-104, io/pulsation_output.f90:76-136, :178-196, :248-253, io/run_files.f90:270-279.**
   - **Issue:** M110's kind reaches beyond the .mod. The opens go through `check_io`, but the history rows and their flush, the profile rows, the GYRE and FGONG writes and the inlist copy have no `iostat`, so a full disk ends the run in a runtime error. `write_gsm_pulse` ignores `herr` after `h5fcreate_f` (attributes, datasets, close) and returns success. README §3 says `check_io` follows "an open, read or write".
   - **Fix:** `iostat`/`iomsg` on each write through one helper, and check `herr`.
   - A, S.
7. **io/run_cards.f90:30-33 (called at controls_reader.f90:65).**
   - **Issue:** `interpret_kind_cards(s)` takes the star `intent(inout)` and touches none of it; all its work is on legacy and buffer variables. The signature claims a star write that does not happen.
   - **Fix:** drop the argument.
   - A, S.
8. **io/derived_controls.f90:17-39, io/run_files.f90:12-24.**
   - **Issue:** about 30 controls_buffer members are imported "from" legacy_controls, which re-exports the whole buffer by its default `public`. Examples: `rotation_active`, `pulse_gyre_interval`, `opal_eos_table_path`, `use_thoul_fit`, and the 30 table paths in run_files. A reader looking for `rotation_active` in legacy_controls finds nothing.
   - **Fix:** take buffer names from controls_buffer and legacy names from legacy_controls. legacy_controls `private` with a public list (it need not re-export the buffer).
   - A, S.
9. **Stale text from the split and earlier moves:**
   - controls_sync.f90:13,16 (from gen_controls_state.py:216,219): "call read_input(s, ierr)" and "All take the star". read_input takes only `ierr`, and the list skips adopt, map, check and derive.
   - controls_check.f90:7: "the controls read calls it (read_input)". read_controls does.
   - controls_reader.f90:177-178 and legacy_controls.f90:11: a canonical home "controls_buffer or star". No row's home is the star.
   - legacy_controls.f90:239: cites ctl2_renames.py, which does not exist.
   - controls_buffer.f90:41-42: "new-format inlists only ... both output modes". There has been one mode since 2026-08-28.
   - controls_buffer.f90:93 and registry :388: "short/log stream".
   - registry :40 and :219: ALL-CAPS "DEFAULT IS ..." notes citing read_controls DATA statements that are gone; they are copied into controls_registry.md and the .defaults files.
   - registry :370, :384-387: "MESA mode:".
   - check_boundaries.py:31-39: names 4 of the 10 buffer importers.
   - check_boundaries.py:703-706: puts map_user_inputs in setup/ and controls_sync in state/.
   - **Fix:** comments only.
   - A, S.
10. **io/pulsation_output.f90:65-68.**
   - **Issue:** FGONG's glob(6) is written as 0 with the comment "mixing length alpha not tracked here", though `s%mixing_length_alpha` is tracked (the history header and the .mod write it). glob(4) and glob(5) are the outermost point's current Z and X, where the FGONG spec has the initial Z and X; I recall MESA's pulse_fgong.f90 does the same, but have not checked.
   - **Fix:** after checking pulse_fgong.f90, write alpha and the card's starting X and Z (`current_card_envelope_x`/`_z`).
   - B (FGONG files change), S.

## What would raise each score

- **Readability, to 8.5:**
  - indent controls_reader's else-branches;
  - turn derived_controls' and run_cards' format labels into named format parameters (as model_files does);
  - one declaration and one comment per controls_buffer member;
  - import buffer names from controls_buffer (finding 8);
  - rename the header's disk fields and `lovstm` (finding 5);
  - fix finding 9's text and the registry's ALL-CAPS notes.
- **Organization, to 8.5:**
  - bind the inlist groups to the code names, so an inlist value is copied once instead of three times (TARGET_ARCHITECTURE 4.6, still open);
  - arm the five no-input run copies from their real sources and drop them from the records;
  - put the run-file paths, `output_dir` and the diagnostics gate in a record that run_yrec owns;
  - move io_units and the gate out of io, so physics stops importing io;
  - lay controls_buffer out in `&control` and `&physics` sections. This is byte-identical; only the dump order and pins move.
- **Maintainability, to 8.5:**
  - controls_buffer and legacy_controls `private`, with generated public lists (the ~700 variables);
  - every write checked (M110, finding 6);
  - M112's single run-log path;
  - M120's dump, with the six scales and the 36 run copies;
  - CI running test_io, test_controls, test_generated_files (all six generators) and test_boundaries (M125);
  - gen_inlist's declaration check (finding 1);
  - star contracts on io's drivers (1 of 35 today).
- **Extensibility, to 8.5:**
  - one declaration plus one registry row plus one generator run per control (the binding above);
  - expand_paths generated from the registry's path rows (a `path` type or column replacing test_generated_files' hand-kept NOT_PATHS), so a table path needs no run_files edit;
  - a registry linter or gen_inlist message that says which cells to fill for each case (inlist name equals code name, differs, run copy, table path).
- **Docs, to 8.5:**
  - findings 2-4 and 9;
  - the open doc items M109 (the alias sentence), M113, M117 and M120 ("any change shows up");
  - the records paragraph as a short table (record, what it holds, who reads it);
  - recipes that name the registry cells per case;
  - a coverage paragraph remeasured after the split (it is dated 2026-10-29).

## Traces and counts (supporting)

- **Add a control** (read by a domain):
  - 4 hand edits: controls_buffer, the registry row, the domain's settings type, the setup/ builder;
  - tool runs: gen_controls_state, gen_inlist, gen_controls_dump, check_star_contracts --fix, gen_star_docs, then make, check_boundaries and the `YREC_REGEN_CONTROLS` test_controls run;
  - about 15 files change.
  - A control core reads takes 3 hand edits.
- **With a run copy:** add a star_info_lib member (control's shape, comment naming the control and controls_sync) and the code that changes it, so 5-6 hand edits. The `arm_run_copies` line is generated, and check_control_copies checks the declaration, the comment and the arming.
- **Add a table path:** 6 hand edits (buffer, registry, run_files' expand_paths, the domain's paths record, the builder, the domain's init) plus the same tool runs, and finding 1's trap.
- **Add a history column:** 2 files (history_output.f90 in 3 places, history_columns.list). A new value adds star_info_lib, observables.f90 and gen_star_docs. A default column means reseeding the baselines.
- **Add a profile column:** 2 files (stitched_model.f90 in 5 places, profile_columns.list), with the same additions.
- **Add a model-file field:** a global value is 1 file in 2 places. A per-shell column is model_files (3-4 places) plus test_io (3 places) and a tag bump.
- **Module variables:**
  - default-public: controls_buffer 366, legacy_controls 334 (+5 parameters). legacy_controls re-exports the buffer and control_slots' 27 constants, so about 700 in all, the same as at the third evaluation.
  - explicitly public and mutable: io_units' 8 unit numbers.
  - private: every other io module's state.
- **The two records:** 36 run copies. `s%control` is read in 13 modules and `s%physics` in 32 (state/README).

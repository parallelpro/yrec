# io -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

Scope: I read every hand-written .f90 of src/io/ and src/io/test/ whole, and src/io/README.md. For the generated files (controls_sync.f90, the *.inc and controls_dump.inc) I reviewed the approach and the sixth pass's changes to them. I read whole tools/controls_registry.py, gen_controls_state.py, gen_controls_dump.py, gen_inlist.py, gen_default_columns.py and generated_files.py, and the header of gen_makefile_deps.py. I also read src/defaults/ (the registry, the two .list files and both rename maps) and the five io test files at the repository root (test_controls.py, test_io.py, test_generated_files.py, test_mesa_output.py and test_mesa_baselines.py). Where io's claims touch other files, I read those parts: the Makefile's generator rules, ci.yml, check_boundaries.py, check_star_contracts.py, test_error_exit.py, setup/kap_settings_builder.f90, state/control_slots.f90, docs/ROADMAP.md §13 ("Sixth-evaluation pass"), docs/bug_sweep.md §13 (M181-M195), docs/TARGET_ARCHITECTURE.md §4.6 and docs/GUIDELINES.md (naming, declarations, code style). For the F77 comparisons I read 6cd56739's src/parmin.f.

I diffed every io file from 2632f06c, the sixth evaluation's commit, to 823326a7: 15 commits touch io's files. I ran these, all read-only:
- all six generators' `--check`, through a scratch runner that writes no bytecode;
- `check_star_contracts.py`: OK, 149 declared;
- `check_boundaries.py`: OK.

The tree was unchanged after these runs. Nothing was built and no test was run. The counts below come from scratch Python in SP/review7_io/. Paths are under src/ unless they start at the repository root; the test_*.py files and docs/ are at the root.

## Scores (1-10)

Fourth -> fifth -> sixth -> now:

- **Readability 7.5 -> 8.5 -> 9 -> 9.5.** The sixth report left three items, and all are fixed:
  - descriptions now use one naming space, which gen_inlist enforces;
  - legacy_controls' 13 history comments are now a Sources note;
  - no io Fortran line is over 110 characters.

  The code reads at rotation's level. Left:
  - controls_check's two capitalized messages (deferred to a re-pin round, M189);
  - one digit count is spelled in three constants (finding 4);
  - seven formats are still written inline, while the rest of io names its formats (yrec_output.f90:192 and :412, history_output.f90:374, profile_output.f90:86, pulsation_output.f90:93, :281 and :285);
  - four procedures have no header comment of their own: `parse_columns` (output_columns.f90:74-76), `read_controls` (controls_reader.f90:51-53), `write_model_line` (run_log.f90:58-61) and `map_user_inputs` (user_inputs.f90:107);
  - the reader uses `getarg` and list-directed `print *` / `write(*,*)` (controls_reader.f90:99-111, :147-151);
  - four tool and test docstring lines run to 125-141 characters (gen_inlist.py:52, gen_controls_dump.py:61, test_io.py:9, test_controls.py:386).
- **Organization 7.5 -> 8.5 -> 9 -> 9.** The sixth report left two items, and both are partly done:
  - The fixed units are now named parameters (run_files.f90:40-44).
  - derive_controls' overrides are documented where a reader looks: the descriptions, the module header and README §1.

  The score stays at 9 because io keeps module state that rotation's standard has none of:
  - six public, mutable unit variables (run_files.f90:21-22, :36-38);
  - the solver-diagnostics gate;
  - the buffer between reads.

  ROADMAP and TARGET_ARCHITECTURE §4.6 list all three as open. Also left:
  - the overrides are still silent when the run happens (M181, the authors' decision);
  - small new items: one read error bypasses run_messages (finding 2), an export nothing outside its module uses (finding 3), and defaults/ has no README and holds an unused, stale map (finding 6).
- **Maintainability 7.5 -> 8 -> 9 -> 9.5.** All the sixth report's tool and test items are fixed. I checked each in the code:
  - M185: the cells guide is printed on make's first refusal, a 3-D control dumps, and a deleted generated file is remade;
  - M186: stale hash rows and orphan expected files fail;
  - M187;
  - M183: GSM is written under the run's own name.

  M174's last window is closed and has a test that races two processes. test_io now pins 96 checks (91 before). Left:
  - CI's first GitHub run is still in progress (the branch was pushed today);
  - the coverage figure predates the sixth pass's code (finding 7);
  - the run plan's rescale-target line still depends on the rows' order (finding 1).
- **Extensibility 7.5 -> 8 -> 9 -> 9.5.**
  - Every refusal make shows now carries the cells guide (gen_controls_state.py:67-88).
  - The output-file recipe names the exact-file-set test (README.md:236).
  - I traced the table-path recipe against the tools that read it (Supporting, below).

  A control a domain reads still takes four hand edits by design. Left:
  - adding a rescale target needs a format edit that nothing points to (finding 1);
  - changing the interior file's digit count takes three edits (finding 4);
  - the optional `gen_inlist --new`.
- **Docs 7.5 -> 8.5 -> 8.5 -> 9.** The sixth report's user-facing inaccuracies are all fixed:
  - diffuse_lithium's overrides (M181, documented);
  - the descriptions' naming spaces (M188);
  - GSM (M183);
  - the stale count, generator claim, docstrings and recipe line (M189);
  - M130's calibration case (README.md:109).

  What is left is developer-facing and small:
  - the legacy Sources note undercounts the additions (finding 5);
  - defaults/ has no README (finding 6);
  - a stale gen_inlist docstring, the test_io header and the test_mesa_baselines docstring (finding 7);
  - a coverage figure dated before the code it describes (finding 7).

Average 9.3 (8.9 at the sixth evaluation, 8.3 at the fifth).

## Summary

- **The sixth-evaluation pass did what the sixth report asked.** The report had 24 items, 10 findings and 14 "what would raise" items. I count:
  - 17 done;
  - 6 partly done: the units and the gate; CI, which is pushed and running; and two items that each appear twice, M181's visibility and the two capitalized messages;
  - 1 open, by choice: `gen_inlist --new`.

  I checked each item in the code, not from the ROADMAP. The two capitalized controls_check messages stay as they were, deliberately: one is in 35 pinned logs, and test_error_exit asserts the other.
- **No open io bug.** io has no open defect of its own in bug_sweep apart from M181, which is documented and waits on the authors. Its open neighbours are M112 (direct run-log writes outside io), M114, M170 and M176.
- **The new code is careful:**
  - `run_log_path` is a pure function, tested at the edges.
  - `profile_pulse_path` keeps every name below 10,000 profiles. no_shared_name checks it across the digits' edges.
  - `finish_replacement_file` serves GSM and the files io writes itself alike.
  - `open_output_file` unlinks the path first, so another run's open cannot cut this run's history. The test_io.py test that races two processes covers it, as does a test_io check.
- **What keeps io below rotation's 10:**
  - the module state above;
  - the silent overrides (M181);
  - a handful of small structural and text items, none of which changes a result.
- **New findings: seven, all Kind A (no result changes) and size S.** The one closest to a bug is finding 1. The run plan's target line labels four values by position: a fifth rescale target would make the internal write overflow its record, and a reordering of the rows would mislabel the targets.

## Previous report's items: status

**What would raise readability:**
- Finding 2, one naming space in the descriptions, enforced: **DONE** (M188, ba1f6ef5). The rule is in `description_problems`, controls_registry.py:390-399 and :423-426. The seven descriptions now cite inlist names, or none: controls_buffer.f90:318, :488, :809-810, :1094, :1180, :1185 and :1193.
- Finding 10's legacy comments and messages: **PARTLY**.
  - The comments are done, in 5c9d9dc5: a Sources note, legacy_controls.f90:155-169. The note undercounts the later names (new finding 5).
  - The capitalized messages are deferred to a re-pin round (M189): controls_check.f90:75 and :173/:178.
- The 161-character banner line: **DONE**, controls_reader.f90:100-101.

**What would raise organization:**
- run_files' units and the gate in a record run_yrec owns, or named parameters for the fixed units: **PARTLY**.
  - The fixed units are named parameters: `terminal_output_unit = output_unit`, `last_model_file_unit = 11`, `run_log_file_unit = 20` (run_files.f90:40-44, 99d3908d).
  - The six public mutable unit variables and the gate remain (run_files.f90:21-22, :36-38). ROADMAP lists them as "Still open".
- derive_controls' overrides made visible (finding 1): **PARTLY**.
  - They are visible to a reader: the descriptions (controls_buffer.f90:890-892, :1115-1155), derived_controls' header (:9-10) and procedures (:160-166, :294-314), and README §1's row (README.md:26).
  - They are still silent at run time: M181 is open, the authors' decision, with options (b) and (c) measured.

**What would raise maintainability:**
- CI's first run on GitHub: **PARTLY**. The branch was pushed on 2026-11-10 and the run is in progress.
- Findings 3, 6 and 7: **DONE** (M185, M186, M187).
  - gen_controls_state.py:67-88;
  - gen_controls_dump.py:69-111, which declares `i1..i<rank>`;
  - Makefile:106-110;
  - test_controls.py:121-136;
  - test_io.f90:1484-1488, which fills `deck_lines` before the run-log checks at :1489.
- Finding 5's GSM rename: **DONE** (M183): pulsation_output.f90:199-206 and :231; output_files.f90:250-265.

**What would raise extensibility:**
- Finding 3, every refusal make shows carries the cells guide: **DONE** (gen_controls_state.py:85-87).
- Finding 9's output-file recipe line: **DONE** (README.md:236).
- `gen_inlist --new`: **OPEN** (optional; ROADMAP "Still open").

**What would raise docs:**
- Findings 1, 2, 5 and 9: **DONE** (819c55a9, ba1f6ef5, 06e9e578, b676fce0).
- The eight long README paragraphs: **DONE**. None is over 1,000 characters now. One table row, the gen_inlist row at README.md:49, is 1,031.
- README's M130 note: **DONE** (README.md:109).

**The sixth report's new findings:**
1. **PARTLY** (M181): documented; the behaviour waits on the authors.
2. **DONE** (M188).
3. **DONE** (M185), all four parts.
4. **DONE** (M182): `run_log_path`, run_files.f90:205-224; test_io.f90:549-560.
5. **DONE** (M183).
6. **DONE** (M186).
7. **DONE** (M187).
8. **DONE** (M184): `profile_pulse_path`, yrec_output.f90:344-362; test_io.f90:591-618.
9. **DONE** (M189). I checked each text:
   - README.md:189 gives 31 values, which I recount from the pins (Supporting);
   - README.md:44;
   - test_controls.py:10-13;
   - test_mesa_output.py:9-10;
   - README.md:236.
10. **PARTLY**: the comments and the README paragraphs are done; the two messages are deferred to a re-pin round.

## New findings (not in bug_sweep)

1. **io/run_cards.f90:44-45 and :135-136 (with state/control_slots.f90:78-82)**
   - **Issue:** the run plan's targets line labels its values by position, while the values come from `rescale_params` by row number.
     - The format is `'  M =',f9.6,'  X =',f9.6,'  Z =',f9.6,'  core mass =',f10.6`, written with `(rescale_params(i,n), i = 1, n_rescale_targets)`.
     - So the labels are right only while control_slots' `i_target_mass..i_target_core_mass` are 1..4 in that order. Reordered rows print mislabelled targets.
     - With a fifth row (`n_rescale_targets = 5`), the internal write of a rescaling card reverts the format into a second record of the scalar `line`, a runtime error ("End of record"), until someone edits the format. Nothing points to it: check_named_rows accepts the loop.
     - The ROADMAP's sixth pass calls 2b2f7975 "io's last reliance on the rows' order", but the echo still relies on it.
     - Latent: n_rescale_targets is 4. The F77 had the same positional echo (`WRITE(IOWR,452) NMODLS(NKIND),(RESCAL(I,NKIND),I = 1,4)`, 6cd56739 parmin.f:1346-1350).
   - **Fix:** write the four targets by name (`rescale_params(i_target_mass,n), ...`), as the store above it does.
   - **Kind:** A (the text is unchanged). **Size:** S.
2. **io/controls_reader.f90:127-132 (and :99-101, :111, :147, :150-151)**
   - **Issue:** an inlist with the old group names `&star_job`/`&controls` is refused through three `write(*,'(a)')` lines, not `report_error`.
     - Its error has no `ERROR:` head, unlike every other read error. README.md:95 and :106 say every error goes through run_messages as `ERROR: <text>`.
     - The reader's other terminal lines also go to `*` rather than `terminal_unit`: the banner, and the "inlist : ..." / "CONTROL namelist : ..." lines.
     - This was present at the sixth evaluation, but not reported then.
   - **Fix:**
     - `report_error(trim(control_nml_file) // ' names its groups &star_job and &controls; ... python3 src/tools/upgrade_inlist.py ' // trim(control_nml_file))`. test_error_exit.py:172 asserts only `upgrade_inlist.py deck.inlist`, which the line keeps.
     - The banner and file lines through `terminal_unit`.
   - **Kind:** A (terminal only; no pinned file and no test assertion moves). **Size:** S.
3. **io/output_files.f90:77-79 and :143-154**
   - **Issue:** `remove_file` is public, but nothing outside output_files uses it. A grep of src/ finds only its call in `open_output_file` (:137); test_io does not import it.
     - The module header's list of procedures (:7-25) and README §1's output_files row (README.md:37) leave it out.
     - Making it public is against GUIDELINES "Declarations" rule 1 (private by default, export what others use). It also offers an unlink that skips `open_output_file`'s pairing with the open.
   - **Fix:** drop it from the public list.
   - **Kind:** A. **Size:** S.
4. **io/yrec_output.f90:53-60 (used at :262-263 and :359)**
   - **Issue:** the interior GYRE file's digit count is written three times:
     - `interior_model_digits = 5`;
     - `first_shared_profile = 10000`, which is 10^(digits-1);
     - `interior_model_format = '(i5.5)'`.

     A change to one that misses the others reopens M184's collision, or pads to the wrong width (GUIDELINES "Code style" rule 2, each formula once). test_io's `no_shared_name` would catch a mismatch between the first two only if the edges it tests moved with them.
   - **Fix:** `first_shared_profile = 10**(interior_model_digits - 1)`, and the format built from the digit count, or `interior_model_format` derived with a comment tying it.
   - **Kind:** A. **Size:** S.
5. **io/legacy_controls.f90:155-169 (with :8-10 and README.md:138)**
   - **Issue:** the new Sources note says the two groups are "the F77 PARMIN's CONTROL and PHYSICS namelists (6cd56739, src/parmin.f) ... and two later additions", pulse_gyre_interval and impjmod.
     - Against 6cd56739's parmin.f, the groups hold 23 names its namelists lack. One of them, fstch, is in a later F77 PARMIN (the parent of df0eb85d). The other 22 are 2025-2026 additions.
     - In /control/: set_initial_model_number, terminal_interval, report_solver_diagnostics, the eight structure-limit stops (log_L_upper_limit .. nu_max_lower_limit), nu_max_sun, delta_nu_sun, Teff_sun, inlist_used_file, profile_data_prefix and pulse_gyre_interval.
     - In /physics/: solar_mass_cgs, G_cgs, the Solberg-Hoiland pair and impjmod.
     - The header's "the Solberg-Hoiland pair ... is the one exception" (to the 2026-10-11 freeze) is right. The note, read on its own, tells a reader that a legacy deck's names are the F77's, which is not so.
   - **Fix:** the note names the 2026 additions made before the freeze as a group, and the Solberg-Hoiland pair, beside the dated F77 additions.
   - **Kind:** A. **Size:** S.
6. **src/defaults/ (flatten_rename_map.tsv, file_rename_map.tsv; tools/controls_registry.py:345-349 and :400-404)**
   - **Issue:** defaults/ has no README, and io/README's tables describe only its registry, `.defaults` and `.list` files. Two files there are undocumented and stale:
     - **flatten_rename_map.tsv** (385 lines, the 2026-08 flattening review map) is read by nothing in the tree: a grep finds no reader. 99 of the 158 `star%` homes it gives name no member of state/star_info_lib.f90.
     - **file_rename_map.tsv** is read only for the F77 routine names in its third column (`description_problems`). 112 of its 117 target paths no longer exist, and its first line says "EXECUTED 2026-08-24".

     A newcomer opening defaults/ for the controls finds two maps that look like inputs and describe a tree that is gone.
   - **Fix:**
     - a defaults/README.md (or a paragraph in io/README §1): each file, who generates it and who reads it;
     - flatten_rename_map.tsv moved out of src/ (docs/history, or deleted, since its commits keep it);
     - file_rename_map.tsv's role stated in its header: it is the list of F77 routine names a description must write as "legacy X".
   - **Kind:** A. **Size:** S.
7. **Stale or unwrapped text in io's tools and tests**
   - **Issue:**
     - **tools/gen_inlist.py:50-60:** the docstring lists what gen_inlist refuses, but not M188's rule that a description names a control by its inlist name. README.md:49 and `description_problems`' own docstring give the rule. Line 52 is 137 characters: the M188 insertion was not rewrapped.
     - **io/test/test_io.f90:1-129:**
       - the header documents the `copy` mode (:125-129) but not the `history` mode (:368-402, `test_io history <path> <times> <rows>`), which test_io.py:104-126 runs;
       - its mutation-dates sentence (:115-119) ends at 2026-11-01. The run-log path, profile pulsation-name and M174-window checks of 2026-11-02 are missing, though README §4 says each check was confirmed.
     - **test_mesa_baselines.py:3-13:**
       - "Runs the standard solar cases in MESA mode": there are no modes since use_legacy_output was retired;
       - "the MESA-format counterpart of the legacy .short pin", "the retire-legacy campaign": the pinned run logs are `.log` (examples/run_standard_solar_model/standard/*.log);
       - test_mesa_output.py:100 still says "the legacy run log".
     - **io/README.md:204-208:** the coverage figure is "93.4% of 2357 lines (2026-11-01, ...)", and output_files is "94% (65 of 69 lines ...)". Since then the sixth pass added `remove_file`, `finish_replacement_file`, `run_log_path`, `profile_pulse_path` and the GSM rename. The per-file counts no longer describe the code, and the figure was not remeasured.
   - **Fix:** text only, and one coverage run (`tools/coverage.py` as README §4 gives it).
   - **Kind:** A. **Size:** S.

## What would raise each score

- **Readability, to 10:**
  - controls_check's two messages reworded, in the re-pin round that M189 defers them to;
  - finding 4;
  - the seven inline formats as named formats (listed under the score);
  - header comments for `parse_columns`, `read_controls`, `write_model_line` and `map_user_inputs`;
  - `get_command_argument` for `getarg`, and the reader's terminal lines through `terminal_unit` (with finding 2);
  - the long docstring lines rewrapped;
  - "GENERATED", "AUTHORITY" and "REBUILD" as plain words in the prose of output_columns.f90:12 and defaults/*_columns.list.
- **Organization, to 9.5-10:**
  - run_files' six unit variables and run_messages' gate in a record run_yrec owns (TARGET_ARCHITECTURE §4.6);
  - findings 2, 3 and 6;
  - output_files.f90 and output_columns.f90, which are star-free by construction, listed in check_boundaries' STAR_BLIND_FILES so they stay star-free;
  - `write_fgong_pulse` taking its six star values (age, `card_initial_x(1)`, `card_initial_z(1)`, log Teff, alpha, model number) as arguments, which would make pulsation_output star-free too;
  - M181 decided. With option (c), the overrides are visible at run time.
- **Maintainability, to 10:**
  - CI's first GitHub run green;
  - the coverage remeasured on today's code (finding 7);
  - finding 1;
  - M181 decided.
- **Extensibility, to 10:**
  - findings 1 and 4, so that a new rescale target or interior digit count is a one-place edit;
  - optionally `gen_inlist --new`, which would leave the declaration, the settings field and the builder's copy as the only hand edits for a control.
- **Docs, to 9.5-10:**
  - findings 5, 6 and 7;
  - README §1's output_files row listing exactly the module's public procedures (after finding 3);
  - the gen_inlist row of the generated-files table (README.md:49) as bullets, as the other long cells became.

## Traces, counts and README checks (supporting)

- **Tools, run read-only on 823326a7:**
  - all six generators' `--check` pass: gen_controls_state, gen_inlist, gen_controls_dump, gen_default_columns, gen_makefile_deps and gen_star_docs;
  - `check_star_contracts`: OK, 149 declared;
  - `check_boundaries`: OK.

  The tree was unchanged afterwards (SP/review7_io/run_check.py and run_tool.py, which add tools/ to the path and write no bytecode).
- **Registry:** 406 rows (380 live, 21 derived, 4 unused, 1 legacy).
  - Homes: physics 192, control 109, buffer 105.
  - 31 run copies: the sixth evaluation's 33, less `first_mc_run` and `last_mc_run` (0ec5ed1d).
  - 41 path rows, 142 settings cells and 6 aliases.
  - 342 of the 380 live controls are bound to their code names, so README.md:60 holds.
- **test_io:** 96 checks pinned in expected_test_io.out, all `ok`. That is the sixth evaluation's 91 plus five: two for the run log path, two for the profile pulsation names, and one for M174's window.
- **Read-twice (README.md:189):** the pinned dumps of cm_rot_base against cm_base, and of Test_solar_GS_rot... against Test_solar_noGS_norot, each give 31 values left at their defaults: 23 controls and 8 run copies.
- **Shared paths (README.md:103):** no two shipped decks in one directory of examples/, testsuite/ or sample_cases/ share a `.mod` or run-log path, counted by deck text with the defaults applied (SP/review7_io/shared_paths.py).
- **Recipe traced, "Add a table path" (README.md:272-276):** each file the recipe names is the one the tools read.
  - The builder `kap_table_paths_from_star` is `result(paths)` with `type(kap_table_paths_t) :: paths` (setup/kap_settings_builder.f90:74-84), which `builder_copies` turns into the record name `kap_table_paths` (controls_registry.py:494-514).
  - The conductive table's row carries `kap_table_paths%conductive`.
  - gen_inlist writes the path into `expand_paths.inc` from the type (gen_inlist.py:289-306), and adds a local for a differing inlist name (:180-193).
  - `kap_table_paths_t` is in kap/opacity_tables.f90, as step 3 says.
- **Code checks:**
  - No io Fortran line is over 110 characters.
  - No `.eq.`-style operators, no commented-out code and no `use` without an only-list, but for `math_lib` (GUIDELINES' stated exception) and HDF5's modules.
  - Every io procedure that takes the star declares its contract (the checker requires io).
  - README §1's rows, §2's steps and §3's bullets match the code as rechecked, including the new `run_log_path`, `profile_pulse_path`, `remove_file`/`open_output_file` and GSM paths, except for the items in findings 3 and 7.
- **F77:** the positional echo of finding 1 is PARMIN's (6cd56739 parmin.f:1346-1350). The 23 names of finding 5 are absent from 6cd56739's namelists (parmin.f:263-302 and :303-350).
- **Not rechecked:** the mutation claims, the 56-deck byte comparisons and the coverage, since nothing was built or run.

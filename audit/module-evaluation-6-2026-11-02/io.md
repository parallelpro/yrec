# io/ -- sixth evaluation (2026-11-02, yrec-next 2632f06c, read-only)

Scope: every hand-written .f90 of src/io/ and src/io/test/ read whole (the generated controls_sync.f90, *.inc and controls_dump.inc reviewed for their approach), src/io/README.md, state/run_messages.f90, tools/controls_registry.py, gen_controls_state.py, gen_inlist.py, gen_controls_dump.py, gen_default_columns.py and generated_files.py whole, gen_star_docs.py and gen_makefile_deps.py skimmed, defaults/controls_registry.tsv; and the parts of test_controls.py, test_io.py, test_generated_files.py, test_boundaries.py, test_mesa_output.py, the Makefile, ci.yml, run_driver, yrec_reset, calibration, shell_boundaries and bug_sweep/ROADMAP that io's claims touch. All six generators' `--check` pass, and the tree is unchanged after them. Nothing was built or run; the counts below come from scratch Python over the tree (SP/review6_io/). Paths are under src/ unless they start at the repository root.

## Scores (1-10)

Fourth -> fifth -> now:

- **Readability 7.5 -> 8.5 -> 9.** The code now reads at rotation's level: no F77 shapes left, named formats, one documented member per line, assumed-shape arguments, and descriptions in today's names, which gen_inlist enforces. Left:
  - the descriptions name other controls in two name spaces (finding 2);
  - legacy_controls' 13 history comment lines, 4 of them in capitals (finding 10);
  - a few lines over 110 characters.
- **Organization 7.5 -> 8.5 -> 9.** There is one recipe per task, one parser of the buffer, model_files no longer imports run_files, and the output's state is one record. Left:
  - run_files' six public mutable unit numbers and the solver gate (acknowledged in ROADMAP and TARGET_ARCHITECTURE 4.6);
  - derive_controls rewrites about a dozen deck controls in place and leaves no trace for the reader (finding 1).
- **Maintainability 7.5 -> 8 -> 9.** Every io bug the fifth report named is fixed, and I checked each in the code: M108, M111, M112's io part, M171-M173 and M178-M180. The tests now check what they claim (the per-deck difference hashes, a read-twice test across two decks). Left:
  - CI's first GitHub run;
  - two test holes (findings 6 and 7);
  - the tools' second registry check (finding 3).
- **Extensibility 7.5 -> 8 -> 9.** Adding a control a domain reads takes 4 hand edits (2 of them by design: the settings field and the builder's copy), then make (which reruns the three generators), `check_star_contracts --fix`, and a re-pin that moves pristine.txt only. A check now ties each settings field to its control. Left:
  - make's first refusal of a bad row comes from gen_controls_state, without the cells guide;
  - the output-file recipe omits test_mesa_output's exact file set.
- **Docs 7.5 -> 8.5 -> 8.5.** Every docs item of the fifth report is fixed, and the §1-§3 claims I rechecked hold. Not raised, because of new user-facing inaccuracies:
  - the documented meaning of diffuse_lithium leaves out that the read turns metal settling on, in 169 shipped decks that set `use_diffusion_z = .FALSE.` (finding 1);
  - the claim that "the pulsation files" are written under the run's own name is false for GSM (finding 5);
  - a stale count ("29 controls") and stale docstrings (finding 9);
  - 8 README paragraphs over 1,000 characters (rotation's README has none, and state and setup adopted a 1,000-character limit).

Average 8.9 (8.3 at the fifth evaluation).

## Summary

- **The cleanup round did what the fifth report asked.** Of its 30 "what would raise" items and 10 new findings: 37 done, 2 partly (M112, which is open only outside io; the units and the gate), 1 open (CI's first GitHub run). I checked each in the code, not from the ROADMAP.
  - Descriptions in today's names, enforced by `controls_registry.description_problems`.
  - One parser of controls_buffer.
  - make reruns the controls generators (stamp `.controls_generated`).
  - The registry's `settings` column, with `check_settings_column`.
  - Per-deck difference hashes, so a new control moves pristine.txt alone (M180).
  - A real read-twice test: a rotating deck, then a non-rotating one, in both deck forms. I count 31 leftover lines, which it would catch.
  - Writes checked on disk for every file io keeps open.
  - M177's write-then-rename for the files written in one go.
- **Owing to the cleanup round:** io has no open bug of its own in bug_sweep. Its open neighbours are M112's 186 direct run-log writes outside io, M114 (net), M170 (kap) and M176 (core).
- **What keeps io below rotation's 10**, apart from the module state already acknowledged:
  - **diffuse_lithium is documented by omission.** The read silently turns on metal settling, the Thoul routines without the fit, and Coulomb log 4. 169 shipped decks set `use_diffusion_z = .FALSE.` with `diffuse_lithium = .TRUE.`, among them all five rotating configuration-matrix decks, and run with metal settling; their pinned dumps show `physics%use_diffusion_z = T`. Nothing the user reads says so. It also means M129's "186 inlists with use_diffusion_z" undercounts the decks M129 reaches.
  - **Descriptions mix name spaces.** The user-facing .defaults files, keyed by inlist name, cite 5 controls by a code name that cannot be set.
  - **Small tool and test holes:**
    - gen_controls_state's guide-less refusal;
    - a removed deck's stale hash row passes;
    - a test_io check whose input is filled after use.
- **New possible bugs, all latent:**
  - a `.short` anywhere in the run log's path renames the log into another directory;
  - GSM files are still written in place (M177's class);
  - from 10,000 profiles on, an interior GYRE file and a profile's GYRE file can share a name.

## Fifth report's items: status

**What would raise readability:**
- Finding 1, the descriptions and the identifier check: **DONE**.
  - `controls_registry.description_problems` (controls_registry.py:371-413), run by gen_inlist, `--check` included.
  - A grep of controls_buffer.f90 finds none of the cited names (rscale/setcal, tpgrad, hpoint, mdot, starin, microdiff, seculr, checkj, itrot, wtime, xsbe91, sli6, cno_mixture_mode, isotope_mixture_mode, read_first_model).
- The capitals and history at controls_buffer.f90:1049-1053: **DONE** (now :1077-1087, plain text).
- user_inputs.f90:18: **DONE** ("With use_new_rate_variables (legacy LNEWVARS; ...").
- The F77 shapes (model_files' `1.0D-6`/`10D0` and `rewind iread`, run_files' `endif`): **DONE**. No `.eq.`, `double precision`, `endif`, `enddo` or COMMON left in io.
- Assumed-shape `sel(:)` and `sources(:)`: **DONE** (output_columns.f90:81, history_output.f90:406, run_files.f90:275).

**What would raise organization:**
- The gate and run_files' six units in a record run_yrec owns: **PARTLY**.
  - The units now start at `error_unit`, so nothing writes to an unconnected unit (run_messages.f90:42, 8ac4496c).
  - `last_model_unit` .. `isochrone_unit` are still public mutable module variables (run_files.f90:20-21, :35-37), and the gate is module state (run_messages.f90:54). Both are listed open in ROADMAP and TARGET_ARCHITECTURE 4.6.
- model_files not importing run_files: **DONE**. The DUMPED MODEL line goes through `log_line` with the unit (model_files.f90:224-232).
- One run-copy recipe (finding 9): **DONE**. README.md:215 hands over to state/README §7, which has the traced recipe.

**What would raise maintainability:**
- M108: **DONE** (controls_check.f90:119-127; test_controls' kind test).
- M111 and M171: **DONE** (`integer_field`, output_columns.f90:61-72; `interior_gyre_path`, yrec_output.f90:242-251).
- M112: **PARTLY**. io's run-log writes all go through `log_line` (a grep finds no direct write in io), and the open part is the 186 writes outside io.
- M172: **DONE** (profile_output.f90:81-88).
- M173's .mod and isochrone size check: **DONE** (model_files.f90:219; yrec_output.f90:382-397; `check_run_log` after each model).
- CI's first run on GitHub: **OPEN** (needs a push; ROADMAP "Still open").
- Findings 2 (whole-dump hashes), 3 (read-twice) and 8 (two parsers): **DONE**.
  - test_controls.py:84-86 and :111;
  - SECOND_READS at test_controls.py:340-375, which asserts at least 20 leftover controls; I count 31 lines in both forms from the pins;
  - controls_registry.buffer_members is the one parser, used by all three generators.

**What would raise extensibility:**
- make rules for the three controls generators: **DONE** (Makefile:85-104; generated_files.emit rewrites only changed files).
- The settings column with its check: **DONE** (142 rows filled; `builder_copies`, controls_registry.py:464-497; `check_settings_column`).
- Finding 2's hash scheme: **DONE**, with M180 (the pristine dump arms the run copies).

**What would raise docs:**
- Findings 1, 6, 7, 9 and 10: **DONE**.
  - 6: README §3's write bullets hold after M173, but for GSM (new finding 5).
  - 7: the .defaults and registry.md headers name controls_buffer.
  - 10: all seven items fixed:
    - test_generated_files' docstring;
    - test_io's mutation dates;
    - user_inputs:18;
    - an unset logical shown `.false.`;
    - "one line per element";
    - test_net among the importers;
    - the history header's initial_x/z wording at README.md:100.
- M113's `used` column: **DONE** (`none`). M117: **DONE**.
- README §4's read-twice sentence: **DONE** (README.md:163), with a stale count (finding 9).
- §3's "Every write" and §4's test_controls paragraphs as bullets: **DONE** (README.md:89-96, :155-167). Eight other paragraphs are still over 1,000 characters (finding 10).

**The fifth report's new findings:**
1. **DONE** (above).
2. **DONE** (above).
3. **DONE** (above).
4. **DONE** (M178): `run_description = ' '`, controls_buffer.f90:82.
5. **DONE** (M179): `expand_value` reports a path that does not fit, and the output paths are deferred-length.
6. **DONE**.
7. **DONE**.
8. **DONE**.
9. **DONE**.
10. **DONE**.

## New findings (not in bug_sweep)

1. **io/derived_controls.f90:294-307 (with io/controls_buffer.f90:1116-1118, :1121-1123, :1137-1139; derived_controls.f90:161-163 and controls_buffer.f90:889-890)**
   - **Issue:** `diffuse_lithium` silently sets `use_diffusion_z`, `diffuse_helium_active`, `use_thoul_diffusion` and `use_new_diffusion_routines` true, `use_thoul_fit` false and `coulomb_log_choice` = 4, whatever the deck says. This is PARMIN's LDIFLI block, ported. The documented meaning, "Include Li/Be in gravitational settling", in physics.defaults and controls_registry.md, says none of this, and no message does.
     - 169 shipped inlists set `diffuse_lithium = .TRUE.` with `use_diffusion_z = .FALSE.`: 84 in run_from_seed_to_start, 78 in run_from_start_to_dbl, and all five rotating configuration-matrix decks. They run with metal settling; cm_rot_base's pinned dump has `physics%use_diffusion_z = T` under a deck line that reads `use_diffusion_z = .FALSE. ! Compute gravitational settling of metals`.
     - M129 counts its reach as "the 186 inlists with use_diffusion_z" by deck text, so it undercounts.
     - The same pattern, smaller: every rotating run forces `new_triangle_init_flag` on, and its description ("Force a fresh envelope-fit triangle at run start") does not say so.
   - **Fix:**
     - the descriptions of diffuse_lithium, use_diffusion_z, use_thoul_fit and coulomb_log_choice (and new_triangle_init_flag) state the overrides;
     - the 169 decks' comments are corrected, or the decks set the values they get;
     - optionally a `report_note` when the override changes a value the deck set.
   - **Kind:** A for the descriptions; B for the note, since it moves the pinned run logs. **Size:** S.
2. **io/controls_buffer.f90:137-138, :308-321, :486-489, :811, :1091-1095, :1163-1179 (copied verbatim by tools/gen_inlist.py:322-323, :393)**
   - **Issue:** a description names another control by its code name in 5 places and by its inlist name in 10, where the two differ. One description, use_new_rate_variables', does both: `use_structure_time_limits` is a code name (inlist name `lstruct_time`) and `zoning_tolerances` an inlist name (code name `chi_grid_scale`).
     - physics.defaults lists controls by inlist name, so a reader is pointed at names a deck cannot set: "a constant rate (mass_accretion_rate)" for `mass_change`; "to new_species_value" for `change_species_value`.
     - controls_buffer, read as code, cites names that exist only in the inlist (`use_mass_change`, `laol_opacity_table_file`, ...).
     - README §5's "in today's names" does not say which names.
   - **Fix:** one rule (descriptions name the inlist name). `description_problems` flags a code name whose inlist name differs, or gen_inlist writes each control's inlist name in the .defaults files.
   - **Kind:** A. **Size:** S.
3. **tools/gen_controls_state.py:71-94 (with src/Makefile:99-104, tools/gen_controls_dump.py:85-87, :98, :55)**
   - **Issue:** make runs gen_controls_state first, and its own registry-versus-buffer check repeats `controls_registry.check_rows_against_buffer`, but more weakly.
     - It names a row by its legacy cell ("registry row -: canonical my_scale is not in controls_buffer" for every new control) and prints no cells guide. So the commonest newcomer slip, a row missing or misspelled, stops make without the guide that README §5 step 3 promises ("stops the build with the problem and the cells to fill").
     - **Also in the tools:**
       - gen_controls_dump's subroutine declares only `i1, i2`, though `writes_for` indexes up to `i3`, so a 3-D control breaks test_controls' build;
       - its header says it reads control_group_def.inc and physics_group_def.inc, which it does not;
       - the stamp rule's empty recipe (`$(CONTROLS_GENERATED): .controls_generated ;`) does not remake a generated file that was deleted.
   - **Fix:**
     - gen_controls_state calls `check_rows_against_buffer` and prints `CELLS_GUIDE` on failure;
     - declare `i3` (or the indices the dims need);
     - correct the header;
     - remove the stamp when a generated file is missing.
   - **Kind:** A. **Size:** S.
4. **io/run_files.f90:186-196**
   - **Issue:** `open_run_log` finds the legacy `.short` with `index(log_output_file, '.short')`, the first occurrence anywhere in the path, and cuts the path there. A run log at `runs.short_tests/cm.short` opens `runs.log` in the parent directory, while `output_directory` (and with it inlist_used, the history and the profiles) stays `runs.short_tests/`. The same goes for any name containing `.short` (`my.shorthand` gives `my.log`). Latent: the shipped decks end in `<deck>.short`.
   - **Fix:** take the suffix only (`n >= 6 .and. log_output_file(n-5:n) == '.short'`), and add a test_io check.
   - **Kind:** A. **Size:** S.
5. **io/pulsation_output.f90:193 (with README.md:94, :96; io/output_files.f90:14-16, :48-51)**
   - **Issue:** `write_gsm_pulse` creates the file at its path (`h5fcreate_f(..., H5F_ACC_TRUNC_F, ...)`), not under the run's own name. So:
     - README's "the files io writes in one go (the profiles, the pulsation files, the copy of the deck)... written under the run's own name... renamed over its path" is false for GSM;
     - decks that keep the default `profile_data_prefix` on an HDF5 build still race on `profile<N>.data.GSM`, M177's class.
     - GSM also has no size check (HDF5's status only), which README:92 does state.
   - **Fix:** create the file at `replacement_path(pulse_path)` and rename it (output_files exports the name and the rename for library-written files), or narrow README:94 and :96 to "GYRE and FGONG".
   - **Kind:** A. **Size:** S/M.
6. **test_controls.py:116-123**
   - **Issue:** when hashes.txt differs from the computed list, the test fails only for a current deck whose hash moved. A deck that was removed or renamed leaves a stale row, the test passes, and the pin keeps a row that nothing checks. The same goes for a representative deck's expected file that outlives its deck.
   - **Fix:** fail on any row of hashes.txt whose deck is no longer in `decks()`, and on any expected_controls file that no representative deck writes.
   - **Kind:** A (test only). **Size:** S.
7. **io/test/test_io.f90:1276, :1392-1394, :1405-1407 (filled only at :1419-1421)**
   - **Issue:** the run-log checks, "a run log cut short by a full disk" and "rewound and written shorter", log `deck_lines` before the program fills them. Under the Makefile's `-finit-local-zero` they are 100 lines of 200 NUL bytes, which is why the full-disk check exceeds its 2000-byte limit. Under ifort or ifx (supported by the Makefile, with no such flag) the lines are indeterminate, and blank lines would make the full-disk check fail spuriously.
   - **Fix:** fill `deck_lines` before the run-log checks (move :1419-1421 up).
   - **Kind:** A (test only). **Size:** S.
8. **io/yrec_output.f90:249 and :305-306**
   - **Issue:** the interior GYRE file is `<prefix><model>.data.GYRE`, padded to five digits, and a profile's GYRE file is `<prefix><N>.data.GYRE` with N unpadded. From N = 10000 on the two name spaces meet: profile 12345's pulsation file and model 12345's interior file are the same path, and the later write replaces the other. With `profile_interval = 1` both are written in the same `output_write_model`, the interior file last. Latent: it needs at least 10,000 profiles with GYRE pulse files and `pulse_gyre_interval > 0`.
   - **Fix:** give the interior files a distinct infix (`<prefix><model>.interior.GYRE`) or pad both, and add a test_io check.
   - **Kind:** A. **Size:** S.
9. **Stale or inexact text:**
   - README.md:163 says the read-twice test "leaves at their defaults 29 controls". The pins give 31 lines (23 of them controls, 8 run copies) in both the inlist and the legacy form.
   - README.md:44 says "a generator rewrites only the files whose text changes", but gen_default_columns.py writes default_columns.inc unconditionally (it does not use generated_files.emit).
   - test_controls.py:10-11 says a deck with a missing starting model "stops inside read_controls". It stops in `open_run_files`; the read opens no file.
   - test_mesa_output.py:9-12 still describes "MESA mode writes exactly {history.data, profile*.data, CASE.log}", while the check is .log, .mod, history.data, inlist_used and the profiles.
   - README.md:196, "Write another output file", omits that a file written by default changes test_mesa_output's exact file set and §2's list of a job's files.
   - **Fix:** text only. **Kind:** A. **Size:** S.
10. **io/legacy_controls.f90:49, :62-63, :67-69, :102, :105-106, :109, :111, :116-117 (and io/controls_check.f90:75, :173, :178; io/README.md:66, :119, :121, :141, :169, :173, :196, :207)**
    - **Issue:**
      - the frozen namelists carry 13 history comments ("MHP 10/24 ADDED END_DCEN,... USED IN MAIN", "!JNT 09/2025 FOR 05/15", "G Somers END"), against GUIDELINES "Code style" rule 5. Freezing binds the names, not the comments.
      - controls_check's messages shout ("SEMI-CONVECTION", "NOT applied").
      - eight README paragraphs run 1,162-1,724 characters (rotation's longest is under 1,000; the cleanup round held state and setup to 1,000).
    - **Fix:** a Sources note in legacy_controls; plain messages; split the eight paragraphs into bullets.
    - **Kind:** A for the comments and README. The messages are B if any pinned log carries them; none of the shipped decks trips them. **Size:** S.

## What would raise each score

- **Readability, to 9.5-10:**
  - finding 2 (one name space in the descriptions, enforced);
  - finding 10's legacy comments and messages;
  - controls_reader.f90:100's 161-character banner line through `trim`.
- **Organization, to 9.5-10:**
  - run_files' six units and the solver gate in a record run_yrec owns, or named parameters for the fixed ones (terminal 6, final model 11, run log 20);
  - derive_controls' overrides made visible (finding 1).
- **Maintainability, to 9.5-10:**
  - CI's first run on GitHub;
  - findings 3, 6 and 7;
  - finding 5's GSM rename.
- **Extensibility, to 9.5-10:**
  - finding 3, so that every refusal make shows carries the cells guide;
  - finding 9's output-file recipe line;
  - optionally, a `--new` mode of gen_inlist that writes the registry row from a declared member. That would leave the declaration, the field and the copy as the only hand edits.
- **Docs, to 9-9.5:**
  - findings 1, 2, 5 and 9;
  - finding 10's eight long README paragraphs;
  - README.md:100 saying, as M130 records, that in a solar calibration FGONG's glob(4)-(5) and the profile's initial_z are the first cycle's values.

## Traces, counts and README checks (supporting)

- **Add a control a domain reads** (README §5): 4 hand edits.
  - The hand edits:
    - the declaration;
    - the row with its settings cell;
    - the settings field;
    - the builder's copy.
  - Then:
    - make, which reruns gen_controls_state, gen_inlist and gen_controls_dump;
    - `check_star_contracts --fix`;
    - `gen_star_docs --check`;
    - check_boundaries, which prints the settings cell;
    - `make tests`;
    - `YREC_REGEN_CONTROLS=1` test_controls, where only pristine.txt moves;
    - pytest test_generated_files and test_boundaries.
  - About 8 commands and 9-10 changed files, against 15 at the fifth evaluation. The fifth's 3 hand-run generators and wholesale re-pin are gone.
- **A control the run changes:** steps 1-2 of the above, then state/README §7 (the member, its comment, its readers and writers, its settings cell moved), 6-7 hand edits. The arming and dump lines are generated, and pristine.txt moves by 2 lines.
- **A table path:** 5 hand edits (declaration, `path` row, paths-record field, builder copy, domain open), with nothing in the read. A path that does not fit stops the read, naming the control.
- **History column, profile column, .mod field, starting-model reader:** unchanged from the fifth (2, 2, 2-4 and 3 files). The output-file recipe is right but for finding 9.
- **Registry:** 406 rows (380 live, 21 derived, 4 unused, 1 legacy).
  - homes: physics 192, control 109, buffer 105;
  - 33 run copies (equal to arm_run_copies' 33 lines), 41 path rows, 142 settings cells;
  - 342 of the 380 live controls bound to their code names (README.md:60 holds).
- **Module state:**
  - io's mutable public variables are run_files' 6 units and run_messages' 2 units;
  - controls_buffer's 406 members are public by name, with importers held by `CONTROLS_BUFFER_IMPORTERS`. A grep matches it: the read's 9 files and includes, controls_dump.inc and test_net.
  - Every io procedure that takes the star has its contract lines (a scratch scan of io/*.f90, controls_sync exempt).
- **README claims checked and true:**
  - §1's rows, against the modules;
  - §2's steps, against run_driver.f90:123-246 and yrec_reset;
  - §3's write bullets, but for GSM (finding 5);
  - the paths paragraph (M179), the records tables and the alias rules (inlist_read.inc);
  - the generated-files table and the CI paragraph (ci.yml);
  - test_io's list (91 checks pinned in expected_test_io.out, all `ok`);
  - test_controls' listed cases (test_controls.py);
  - the "Add a control" and "Add a table path" cells (CELLS_GUIDE).
- **Not rechecked:** the coverage figures, the mutation claims and the 56-deck byte comparisons, since nothing was run.

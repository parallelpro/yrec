# io/ -- fifth evaluation (2026-11-01, yrec-next 6a75b72f, read-only)

Scope: every .f90/.inc of src/io/ and src/io/test/, src/io/README.md, state/run_messages.f90, the generators tools/gen_*.py (the controls ones read whole, gen_star_docs and gen_makefile_deps skimmed), tools/controls_registry.py, defaults/controls_registry.tsv; test_controls.py, test_io.py, test_generated_files.py, test_boundaries.py, the CI workflow, and the parts of run_driver, starting_model, model_rescaling and check_boundaries that io's claims touch. All six generators' `--check` pass. Nothing was built or run; the identifier and deck scans below are scratch Python over the tree. Paths are under src/ unless they start at the repository root.

## Scores (1-10)

Earlier scores -> now (the first three on a 1-5 scale, doubled):

- **Readability 5 -> 7 -> 8 -> 7.5 -> 8.5.** The read half now reads as well as the output half: indented branches, named formats, one member per line with its description above it, buffer names imported from controls_buffer. What is left is in those descriptions, now the only ones: 29 of the 406 cite F77 routine or variable names, or name controls that do not exist (finding 1).
- **Organization 5 -> 6 -> 7 -> 7.5 -> 8.5.** The triple copy is gone: 342 of the 380 inlist controls are bound to their code names, and the other 38 take one copy from a local. The read writes only its two records and the run copies. The run's paths are star controls, and no physics domain imports io. Left:
  - the buffer's 406 public module variables between reads, and legacy_controls' 330 private ones;
  - the diagnostics gate and eight unit numbers as public module state;
  - the run-copy recipe in two READMEs.
- **Maintainability 4 -> 6 -> 7 -> 7.5 -> 8.** Gains: no default-public module variables (there were about 700), star contracts on all 26 star-taking io procedures, every output file's writes checked, generators that name the cells to fill, and CI configured. Held back by:
  - the open M108, M111, M112 and M171-M173 (section 9's three were in the round's own goal);
  - CI's first real run, which has not happened;
  - a hash pin that every new control re-pins wholesale (finding 2);
  - a read-twice test that cannot see a missing reset (finding 3);
  - two parsers of controls_buffer in the tools (finding 8).
- **Extensibility 4 -> 7 -> 8 -> 7.5 -> 8.** Better:
  - A table path is 5 hand edits, none of them in the read: expand_paths is generated from the `path` rows.
  - `gen_inlist.py --cells` says what to fill for each kind of row.
  - Columns and .mod fields stay at rotation's level.

  Still the same: a control a domain reads takes 4 hand edits, 3 generators run by hand (make runs none of them, though it runs gen_default_columns), the contracts fix, gen_star_docs and a wholesale re-pin. No check ties a settings field to its control.
- **Docs 6 -> 7 -> 8 -> 7.5 -> 8.5.** Every §1 row, every §2 step (checked against run_driver), the records tables, the seven legacy/inlist group differences (recounted), the CI paragraph and the test paragraphs I could check hold. Wrong or stale:
  - "every write is checked" (finding 6);
  - the generated documents' "edit the registry" headers (finding 7);
  - the descriptions (finding 1);
  - the run-copy recipe (finding 9);
  - a stale test docstring, and the open M113 and M117.

Average 8.3, against 7.5 at the fourth evaluation.

## Summary

- Round 3 did what the fourth report asked, almost item for item. I checked each in the code: of 39 items, 31 are done, 5 partly and 3 open (M112, the one-step control, M117).
  - The controls read is one copy into the code names.
  - The read writes nothing on the star but the two records and the run copies.
  - All 18 modules in scope are `private`; default-public module variables 700 -> 0.
  - The read's state is gone from run_files, and physics imports no io module.
  - Every write to io's output files is checked, through one helper with a full-disk size check (M174's fix reviewed: sound).
- io is now close to rotation's bar in code. What keeps it below 9 is mostly docs and tests, and the open small bugs:
  - the control descriptions, which became the one source of the user-facing .defaults files and still carry F77 names and three non-existent controls;
  - a test_controls hash pin and a read-twice test that check less than they claim;
  - M108, M111, M112 and M171-M173.
- New possible bugs, both latent:
  - a deck without run_description writes NUL bytes into the run log;
  - expanded paths and output paths longer than 256 characters are cut without a message.
- Structural items still open, and acknowledged in TARGET_ARCHITECTURE 4.6:
  - the buffer as module state;
  - the gate and the unit numbers as public module variables;
  - a control still takes four hand edits and three hand-run generators.

## Fourth report's items: status

**What would raise readability:**
- Indent controls_reader's else-branches: **DONE** (controls_reader.f90 rewritten, 236 lines, clean).
- Named formats in derived_controls and run_cards: **DONE** (derived_controls.f90:52-68, run_cards.f90:28-42; no numeric format labels left in io).
- One declaration and one comment per controls_buffer member: **DONE**, and enforced (controls_registry.py:229-266 stops on several members per line, a member without a comment, and `/*`).
- Buffer names from controls_buffer: **DONE** (derived_controls.f90:16-40; run_files imports no buffer name at all).
- Header disk fields and `lovstm`: **DONE** (model_files.f90:61-76: disk_locking_age_gyr, disk_omega_rad_s, intermediate_overshoot_active, with a comment on who uses each field).
- Finding 9's text and the registry's ALL-CAPS notes: **PARTLY**.
  - Every item of finding 9 was rechecked and is fixed (below), and the "DEFAULT IS" notes are gone.
  - controls_buffer.f90:1049-1053 still has "UNITS: Gyr" and "ANGULAR VELOCITY" with history notes.
  - 24 descriptions are written in F77 routine names (finding 1).

**What would raise organization:**
- Bind the inlist groups to the code names: **DONE** (inlist_decl.inc: 342 bound, 38 locals copied once; legacy pairs through legacy_to_code.inc).
- The five no-input run copies: **DONE** (M152; gone from controls_buffer, and arm_run_copies arms 33 copies from real controls).
- The run-file paths, output_dir and the gate in a record run_yrec owns: **PARTLY**.
  - The paths are controls in the star, and output_dir is computed (run_files.f90:167-171) and kept in output_state.
  - The gate is still module state (run_messages.f90:45), and so are run_files' six unit numbers (run_files.f90:34-36).
- Move io_units and the gate out of io: **DONE** (state/run_messages.f90; a grep of eos, kap, atm, net, numerics, convection, element_diffusion, rotation and math finds no use of any io module).
- controls_buffer in &control and &physics sections: **DONE** (controls_buffer.f90:56, :403).

**What would raise maintainability:**
- controls_buffer and legacy_controls private, with generated public lists: **DONE**. All 18 modules are `private`. The buffer's 406 members are public by name (controls_buffer_public.inc), and its importers are held by check_boundaries' allow-list and only-list checks.
- Every write checked: **PARTLY**.
  - Done for io's output files: history, profile, GYRE, FGONG, GSM (herr checked), deck copy, .mod, isochrone (output_files.f90; yrec_output.f90:175-179, :348-351).
  - The .mod and isochrone files have no size check, and the run log's writes none (M173, M112, open).
- M112: **OPEN**.
- M120's dump: **DONE**. The dump holds every run copy (controls_dump.inc), and the six scales are now derived controls in star%physics.
- CI running io's tests and all six generators: **PARTLY**. It is configured (.github/workflows/ci.yml: the six `--check`s, test_io, test_controls, test_generated_files, test_boundaries) and was replayed locally, but has never run on GitHub (it needs a push).
- gen_inlist's declaration check: **DONE** (gen_inlist.py:131, :153-158; controls_registry.py:310-327).
- Star contracts on io's drivers: **DONE**: 26 of the 26 star-taking io procedures (controls_sync's 3 exempt). There were 35 star-taking procedures before; the read's no longer take the star.

**What would raise extensibility:**
- One declaration, one row, one generator run per control: **OPEN**. A control is still declaration + row + domain field + builder copy, plus 3 generators, `check_star_contracts --fix`, gen_star_docs and the re-pin.
- expand_paths generated, with a path type: **DONE** (gen_inlist.py:283-298; the registry type `path`, 41 rows).
- A message naming the cells to fill: **DONE** (`gen_inlist.py --cells`; CELLS_GUIDE is printed on every registry failure).

**What would raise docs:**
- Findings 2-4 and 9: **DONE** (below).
- M109: **DONE**. M113: **PARTLY** (the description is right; the `used` column stays `read`). M117: **OPEN** (controls_buffer.f90:904-905). M120: **DONE**.
- The records paragraph as a table: **DONE** (README.md:96-112, two tables).
- Recipes that name the cells per case: **DONE** (README.md:181, :193; CELLS_GUIDE).
- Coverage remeasured: **DONE** (README.md:147, 2026-11-01; I did not rerun it).

**The fourth report's findings:**
1. The table-path trap: **DONE** (M149). An undeclared code name or a legacy name outside the frozen groups stops gen_inlist, with the guide.
2. Descriptions written as what the run copies do: **DONE** (M152). Spot-checked: mixing_length_array, end_age_stop_active and rescale_params now say what they hold.
3. The records misdocumented: **DONE** (M153). README §2 step 1 names the seven controls; I recounted them against legacy_groups(): exactly 7.
4. The registry semantics drifted: **DONE** (M154). The inlist-only rows have legacy `-`, the star% branch is gone, and the freeze sentence names the Solberg-Hoiland pair (legacy_controls.f90:8-10).
5. model_file_header names: **DONE** (M155).
6. Unchecked history, profile and pulse writes: **DONE** (M138).
7. interpret_kind_cards' unused star argument: **DONE** (run_cards.f90:49, no argument).
8. Imports through legacy_controls: **DONE** (legacy_controls' only public name is read_legacy_deck).
9. Stale text: **DONE**, each item checked:
   - controls_sync's header and controls_check's caller;
   - no "or star" homes left;
   - no ctl2_renames, "both output modes", "short/log", "MESA mode" or "DEFAULT IS" anywhere in io/, tools/ or defaults/;
   - check_boundaries.py:31-39 lists all nine importers.
10. FGONG glob(4)-(6): **DONE** (M130; pulsation_output.f90:78-80).

## New findings (not in bug_sweep), ranked

1. **io/controls_buffer.f90:203, :212, :221, :243, :261, and 24 more descriptions** (:111-116, :465, :598-603, :661-666, :677, :793, :810, :905-913, :1063-1065, :1092-1094, :1121). **Issue:** since round 3 the comment above each member is a control's only description, and gen_inlist copies it into control.defaults, physics.defaults and controls_registry.md.
   - Five of them name controls that exist nowhere: `cno_mixture_mode = 1` (the control is mixture_change_mode), `isotope_mixture_mode = 1` (isotope_change_mode), and "when read_first_model is set" (first_call_flag). A scratch scan of every `word_with_underscore` in the descriptions against every code, inlist, legacy and Fortran name found exactly these.
   - 24 describe the control by an F77 routine or variable instead of what the value means: "rscale/setcal: Z-mix rescale parameter 1", "tpgrad:", "hpoint envelope rezoning", "mix: envelope pressure step control", "mdot: total accreted-mass fraction parameter", "(microdiff)", "starin's mixture algorithm", "(seculr/checkj)", "(itrot loop)", "(wtime/ytime governors)", "overridden by xsbe91 with lxbe91", "sets sli6 = value/5.5e3".

   **Fix:** rewrite them in code names and in terms of the deck's value (M152's rule). gen_inlist warns when a description names an identifier that is no control, inlist, legacy or code name. A, S/M.
2. **test_controls.py:107 (with README.md:185).** **Issue:** each deck's hash covers its whole dump, so adding any control, even at its default, changes all 314 hashes. The recipe's `YREC_REGEN_CONTROLS=1` re-pins them all, and README step 6 says so ("every line of hashes.txt"). For the 294 decks with no representative file, the claim "only the new member's lines change" is then checked by nothing: an accidental change in one of them is re-pinned with the deliberate one. **Fix:** hash each deck's lines that differ from the pristine dump (`differing(dump, pristine)`, already in the file). A new control at its default then moves only pristine.txt and the decks that set it. A (test only), S.
3. **test_controls.py:298-316, io/test/test_controls.f90:51-65 (with README.md:143).** **Issue:** the read-twice test reads the same deck twice. A second read of one deck gives the same records and buffer whether or not seed_controls_buffer resets them, since nothing writes the records between the reads. So it checks only that three scribbled run copies are re-armed, not "resets the controls" (its docstring) or "a second read into a used star starts from the defaults" (README). test_reentry's A-then-B form covers one control (s0_pp) and needs the full input tree. **Fix:** let YREC_TEST_READ_TWICE take a second deck: read cm_rot_base, then cm_base, and compare with cm_base read alone. A, S.
4. **io/derived_controls.f90:340-341, io/controls_buffer.f90:80.** **Issue:** run_description defaults to `repeat(achar(0), 256)`, and `trim` does not remove NUL.
   - A deck that leaves run_description(2) unset gets a run-log line of 256 NUL bytes (its len_trim is 256). One that leaves (1) unset gets them after "description:".
   - control.defaults documents the default as `''` (gen_inlist's `shown()`).
   - Latent: all 1,159 shipped inlists set both elements. A newcomer's deck, or test_controls' empty deck, writes a run log that grep calls binary.

   **Fix:** a blank default for run_description (only the pristine dump's two lines move), or skip a NUL description. A, S.
5. **io/run_files.f90:298; io/history_output.f90:66, :91; io/profile_output.f90:63, :70; io/yrec_output.f90:277, :287 (and :226, 320 characters).** **Issue:** `path_value = trim(value) // path_value(n+1:)` cuts an expansion longer than the 256-character path control without a word, and the writers build their output paths in 256-character buffers the same way. A long `YREC_INPUT` or output directory (an HPC scratch tree) then gives an open error on a cut name, or output written under a wrong name. **Fix:** expand_value reports an expansion that does not fit, and the writers check the length of the path they build (or the paths become 1024 characters). A, S/M.
6. **io/README.md:73 and :89.** **Issue:** "Every write of steps 2-7 is checked: a failed one stops the run with exit 1" and "Every write of a file is checked" are not true of all writes. Step 7's Monte Carlo summaries, the neutrino-flux table and every run-log write (progress lines, summary, `log_line`) have no iostat (M173 and M112 record the code; the README states the opposite). **Fix:** "every write of io's output files (history, profiles, pulsation files, .mod, isochrone, deck copy)", pointing to M173 for the rest. A, S.
7. **tools/gen_inlist.py:305-308 and :337-340 (into defaults/control.defaults:1-4, physics.defaults, controls_registry.md:3-6).** **Issue:** the generated documents say "GENERATED by tools/gen_inlist.py from defaults/controls_registry.tsv (edit that, then regenerate)". Since round 3 the descriptions, and the defaults, live in io/controls_buffer.f90, and the registry has no description column. A reader who fixes a description in a .defaults file is sent to the wrong file. The .inc headers (gen_inlist.py:67-68) name both sources. **Fix:** the headers say "descriptions and defaults: io/controls_buffer.f90; names and groups: the registry". A, S.
8. **tools/gen_controls_state.py:63-158 vs tools/controls_registry.py:190-300.** **Issue:** controls_buffer.f90 has two parsers in the tools:
   - gen_controls_state.parse() (its own TYPE_RE, joined_decl_lines and split_entities), which gen_controls_dump also uses for the buffer-only members;
   - controls_registry.buffer_members() (_DECL_RE, _joined, _entities), which gen_inlist and the row checks use.

   Their attribute regexes and their handling of a missing initializer differ (ZERO vs None), and `constants()` takes only real parameters into the public list. A declaration form one accepts and the other skips gives records and readers that disagree, caught only by the compiler. **Fix:** gen_controls_state and gen_controls_dump take buffer_members(). A, S.
9. **io/README.md:187 vs state/README.md:485-492.** **Issue:** "Add a control the run changes" says "steps 1-3 as above", which run gen_controls_state, gen_inlist and gen_controls_dump. It then declares the star member and runs all three again; none of them reads star_info_lib, so the second round changes nothing. state/README §7 has the traced recipe for the same task, in another order (member first). Two copies of one recipe is how M150 drifted. **Fix:** io's recipe points to state/README §7, as rotation's "Add a control" points to io's. A, S.
10. **Stale or inexact text left by round 3:**
    - test_generated_files.py:1-2 describes a "hand-kept list of path controls" that the round removed;
    - io/test/test_io.f90:93-94 dates the mutation confirmation 2026-10-29, before the writes, M174 and FGONG checks;
    - io/user_inputs.f90:18 reads "With use_new_rate_variables (use_new_rate_variables, ...", a rename-script artifact of 692bc883;
    - gen_inlist.py:98-99 shows an uninitialized logical's default as `0` in controls_registry.md's "Members no deck sets";
    - README.md:185's "one line in pristine.txt" is one per element for an array control;
    - README.md:21 says "only the read's files import" controls_buffer, but test_net.f90 does too;
    - README.md:93 and pulsation_output.f90:70-72 call FGONG's glob(4)-(5) "the history header's initial_z and initial_x". The header writes current_card_envelope_x and _z at the first row (history_output.f90:376-378), which a card-1 rescale sets to its target (model_rescaling.f90:114, :153). Two shipped decks have a card-1 target unlike their initial X or Z (sample_cases/UpperMS/Sample_m9p0feh+0p0GS98_os0p4_TAMS, examples/run_from_seed_to_start/m10000gs98z040000_prems).

    **Fix:** comments and docs only. A, S.

## What would raise each score

- **Readability, to 9:**
  - finding 1 (rewrite the 29 descriptions; the identifier check in gen_inlist);
  - controls_buffer.f90:1049-1053's ALL-CAPS and history into plain text, with a Sources note;
  - user_inputs.f90:18;
  - the F77 shapes left: model_files' `1.0D-6`/`10D0` and `rewind iread` statements, run_files' `endif`;
  - assumed-shape `sel(:)` and `sources(:)` in output_columns, history_output and run_files (GUIDELINES "Declarations" rule 4).
- **Organization, to 9:**
  - the solver-diagnostics gate and run_files' six units in a record run_yrec owns. The fixed units (terminal 6, run log 20, final model 11) could be named parameters, if test_io captures messages another way;
  - model_files not importing run_files for a log label (model_files.f90:178, :218);
  - one run-copy recipe (finding 9).

  The buffer as module state is close to inherent to namelist reads, so I would not score against it.
- **Maintainability, to 9:**
  - M108, M111, M112, M171, M172;
  - M173's .mod and isochrone size check (flush_output_file after each model's write);
  - CI's first run on GitHub;
  - findings 2, 3 and 8.
- **Extensibility, to 9:**
  - make rules that run gen_controls_state, gen_inlist and gen_controls_dump when controls_buffer.f90 or the registry changes, as for default_columns.inc;
  - the registry column naming each control's settings field, with its check (setup/README §6);
  - finding 2's hash scheme, so a new control re-pins only the decks that set it.

  With these, a domain control is declaration + row + field + copy, and nothing else by hand.
- **Docs, to 9:**
  - findings 1, 6, 7, 9 and 10;
  - M113's `used` column and M117;
  - README §4's read-twice sentence matched to the test (or finding 3's test);
  - README §3's "Every write" paragraph (about 1,300 characters) and §4's test_controls paragraph split into bullets, as rotation's README is.

## Traces, counts and README checks (supporting)

- **Add a control a domain reads** (README §5):
  - 4 hand edits: controls_buffer.f90, controls_registry.tsv, the domain's settings type, the setup/ builder;
  - tool runs: gen_controls_state, gen_inlist, gen_controls_dump (by hand: make runs none of them), `check_star_contracts --fix`, gen_star_docs, make, check_boundaries, `YREC_REGEN_CONTROLS=1` test_controls (re-pins all 314 hashes), then pytest test_generated_files and test_boundaries;
  - about 15 files change. Unchanged from the fourth evaluation.
- **With a run copy:** add a star_info_lib member with the checked comment, its writer and its readers: 6-7 hand edits. The arming line is generated. The recipe runs the generators twice (finding 9).
- **A control core reads:** 3 hand edits (declaration, row, use site), plus the same runs.
- **Add a table path:** 5 hand edits (declaration, row with type `path`, the domain's paths type, the builder's copy, the domain's open). No read edit, since expand_paths.inc is generated (was 6). The cells guide covers a differing inlist name.
- **Add a history column:** 2 files (history_output.f90 in 3 places, history_columns.list). A new value adds star_info_lib, observables.f90 and gen_star_docs; a default column means a baseline reseed.
- **Add a profile column:** 2 files (stitched_model.f90: index, table row and the three region subroutines; profile_columns.list).
- **Add a .mod field:** a global is model_files in 2 places plus test_io in 2. A per-shell column is model_files in 3-4 places, test_io in 3, and a tag bump.
- **Write another output file:** the writer (open_output_file, `iostat` writes, close_output_file), its call in output_write_model, and test_io checks copied from write_gyre_pulse's: 2-3 files, plus the control recipe if it has a switch.
- **Module variables:**
  - default-public: 0 (all 18 modules `private`; there were about 700);
  - public by name and mutable: controls_buffer's 406 members (plus 5 parameters; importers held by CONTROLS_BUFFER_IMPORTERS and the only-list check), run_files' 6 unit numbers and run_messages' 2, 414 in all;
  - private module state: legacy_controls' 330 legacy variables (reset at every read), run_messages' pending queue and gate.
- **Registry:** 406 rows (live 380, derived 21, unused 4, legacy 1). Homes: physics 192, control 109, buffer 105. 33 run copies, 41 path rows.
- **README claims checked and true:**
  - every §1 row;
  - §2 steps 1-7 against run_driver.f90:128-246;
  - the seven group differences;
  - 342 of 380;
  - §3's records tables against the registry homes;
  - the generated-files table;
  - the CI paragraph against ci.yml;
  - test_io's list against expected_test_io.out (79 checks);
  - test_controls' and test_error_exit's listed cases;
  - every §5 recipe's file and procedure names.

  Not rechecked: the coverage figures and the mutation claims (nothing was run).

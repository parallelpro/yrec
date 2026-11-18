# atm -- eighth evaluation (2026-11-17, yrec-next 1179afbf, read-only)

Scope: every file of src/atm/ read whole:
- README.md, atm_lib, atm_choices, atm_tables, atm_messages, t_tau_relations;
- tables/kurucz_atmosphere, tables/kurucz_table_reader, tables/allard_atmosphere, tables/allard_table_reader;
- test/test_atm.f90, test/test_atm.mutations.toml, test/expected_test_atm.out;
- also whole: setup/atm_settings_builder.f90 and the repository root's test_atm.py.

Read for the claims they touch:
- atm's callers: core/envelope/envelope_integration_kernel.f90 (`prepare_surface_boundary`, `integrate_atmosphere`), envelope_integration.f90 (`atm_get`, `gray_photosphere`, `report_envint_events`), envelope_integrands.f90, core/structure/surface_boundary.f90 (`switch_atmosphere_choice`, `on_gray_fallback`), core/structure/henyey_iteration.f90 (`refresh_surface_bc`), core/step/mass_change.f90 (`find_accreted_surface_state`), setup/settings_install.f90, io/derived_controls.f90's banner, and the controls' declarations and registry rows;
- numerics' `splint`, for what a spline does outside its nodes;
- the atm entries of tools/check_boundaries.py (PUBLIC, STAR_BLIND_FILES, KERNEL_FILES, STAR_FREE_DIRS, FACADE_ONLY_DOMAINS) and tools/style_counts.tsv;
- GUIDELINES "Naming", "Declarations and interfaces", "Code style" and "A package README (the standard)";
- docs/bug_sweep.md (every atm entry: O19, M112, M221, M231-M240, M283, M284) and docs/ROADMAP.md's rounds 0, 1 and 1B;
- the seventh report (audit/module-evaluation-7-2026-11-10/kap_atm.md) and its FINDING_MAP.md and B_LISTS.md.

How the evidence was gathered:
- every inherited defect was checked against `git show 6cd56739:src/<file>.f` (surfp.f, surfbc.f, envint.f, alfilein.f);
- nothing was built or run: the mutation count and the coverage figure are the README's, checked only as far as reading allows;
- the counts come from scratch Python over the tree (SP/review8_atm/): counts.py (style residue, I/O statements, the mutation manifest against the code and the test's labels), tables.py (the reader replayed on the shipped Kurucz and Castelli/Kurucz tables: each row's valid log g columns), allard.py (the Allard reader's two passes replayed on the four shipped Allard files), decks.py (which `atm_choice` and tables the 1,198 shipped deck files use).

Paths are under src/ unless they start at the repository root. Line numbers are at 1179afbf.

## Scores (1-10)

- **Readability 6.5 -> 9.** The 886-line F77-shaped Allard module is gone: two kernels of about 200 lines each in named steps (`find_stencil`, `interpolate_rows`, `interpolate_across_rows`; `outside_table`, `pressure_in_row`, `spline_in_log_g`), two readers, every threshold and format a named parameter with its meaning, headers that describe the code as it is, a Sources section per file, and no style residue (counts.py: 0 `endif`/`enddo`, 0 old relational operators, 0 lines over 110 characters, 0 commented-out statements, 0 capital-letter comment blocks). Left:
  - the arguments are `log_teff`, `log_g` in atm_lib.f90:119-123 and allard_atmosphere.f90:66, and `log10_teff`, `log10_gravity` in kurucz_atmosphere.f90:52-53 and atm_lib.f90:239, :246 (GUIDELINES "Naming": `log_x` is log10 of x);
  - six fields of `allard_table` that nothing reads, under a comment that says they are counted and printed (finding 4);
  - the choices listed by hand in two messages (finding 5);
  - small text items (finding 8).
- **Organization 6.5 -> 8.5.** One facade (enforced), no star, a settings and a paths record from setup/, the three kernels pure, argument-only and in STAR_BLIND_FILES and KERNEL_FILES, the run's tables private to atm_lib, every message in atm_messages (29 procedures, all called; no I/O statement elsewhere but the two readers' file I/O), and one descriptor (`atm_choices`) that the envelope kernel, the surface boundary and the banner ask. Left:
  - the switch back from gray asks atm about log Teff alone, while the lookup that sent the run to gray can have failed on log g: the working choice then flips at every surface boundary (finding 2);
  - a table's switch to gray has three homes inside atm: the flag in `choices` (atm_choices.f90:49-55), the log Teff in `gray_switch_log_teff` (atm_lib.f90:254-261) and the text in `note_choice_switch` (atm_messages.f90:356-379);
  - `atm_allard_teff_max` is a public entry that only test_atm calls (finding 4);
  - the Kurucz table's edges are constants of the kernel (log Teff 3.5, log g -0.5, the 4.5 and 4.55 rows; kurucz_atmosphere.f90:29-30), not properties of the table read: two shipped tables start below log Teff 3.5 (atmk1990m13 at 3.443 and m17 at 3.463; tables.py);
  - atm_messages writes to `run_log_unit` directly (atm_messages.f90:84, :90), not through `log_line`'s guard: M112's open part, as rotation's and eos's messages modules do.
- **Maintainability 6 -> 8.5.** test_atm has 39 ok/FAIL checks at the kernels, the readers and the facade (exactness at nodes, on a plane and on a polynomial, continuity, the photosphere of each relation, every refusal and error contract), 45 mutations, and pins of the messages. Every atm item of the seventh report is closed in the code. Left:
  - a Kurucz or Castelli/Kurucz lookup below a row's first valid log g is extrapolated without limit, with `ierr` 0, and no check covers it (finding 1);
  - finding 2's flip, and the Allard reader's unguarded paths (finding 3);
  - M284 (open, data): several shipped Kurucz tables hold rows the lookup refuses, and the README does not say which;
  - M112's direct run-log writes;
  - the coverage figure predates round 1B's rewrite of atm_messages (finding 7).
- **Extensibility 6 -> 8.5.** A new T-tau relation is one row in `choices`, a case in three functions and a re-export; a new table adds a path field, a module variable, two cases and a builder line; nothing in core/ changes, which the code confirms (core tests no choice number but `on_gray_fallback`'s 0, core/structure/surface_boundary.f90:33). README §6 has four recipes, two of them traced. Left:
  - the recipes were traced on 2026-11-10, and the next day's rewording added two hand-written lists of the choices that neither recipe names (finding 5);
  - the "table in the Kurucz layout" recipe works for a table of 57 rows, or of a multiple of 4 up to 76: the reader picks one of two fixed formats by row count (finding 6);
  - a table that switches to gray edits three files of atm for the switch alone.
- **Docs 6.5 -> 8.5.** The README has the standard's six items: the table matches the nine files, the flow names each caller with its procedure, the conventions give every `ierr` and where each message goes, the tests section says what each of the 39 checks verifies, and the kept-by-decision list is there. Left:
  - §2's sentence on accretion describes the code before round 2B (finding 7);
  - §5's coverage is dated before round 1B and quotes a message that no longer exists (finding 7);
  - §2 says a failed Allard lookup falls back "for that integration", while core latches the gray choice until the surface boundary switches back (finding 2);
  - M284's proposed note on the tables with holes is not in the README.

Average **6.3 -> 8.6**. No criterion is under 8.

## Summary

- **The round did what the ROADMAP says, and the code shows it.** Checked against the tree: the kernels take their table and return their result (kurucz_atmosphere.f90:52-95, allard_atmosphere.f90:66-94) and are in STAR_BLIND_FILES and KERNEL_FILES (tools/check_boundaries.py:195-199, :1236-1239); atm_lib keeps the tables privately (atm_lib.f90:60-62) and resets them at each `atm_init` (:88-90, M283); a successful Allard lookup writes nothing (:167-195, M235); the reader's tau = 100 check runs (allard_table_reader.f90:143-147) and accretion asks only loaded tables (core/step/mass_change.f90:270, O19); (1/4) log10(3/4) is computed (t_tau_relations.f90:57-60); the style ratchet is at 0 for atm and holds only the readers' 19 I/O statements for atm/tables (tools/style_counts.tsv:2-21).
- **The test claims hold as far as reading goes.** 39 `report`/`report_true` calls, 39 `ok` lines in the pinned output, 45 mutations whose `old` text each occurs exactly once in its file and whose `check` is one of the 39 labels, and no label without a mutation (counts.py). Not re-run.
- **Two inherited behaviours the round did not find.** Both are latent for the shipped decks (the 12 Kurucz decks are dwarfs and the Sun, no deck uses choice 5, the 116 Allard decks are low-mass main-sequence stars):
  - the Kurucz lookup has no per-row log g limit, so a point below a row's first valid column is a cubic (or linear) extrapolation of up to eight columns with `ierr` 0, where the Allard lookup refuses beyond four;
  - an Allard deck whose log g leaves the tables while log Teff stays inside flips between Allard and gray at every surface boundary, with three message lines and a rebuilt triangle each time.
- **The rest is small:** three unguarded paths in the Allard reader, a dead public entry and dead table fields, hand-written choice lists that the recipes miss, a reader format tied to two row counts, and three stale README passages.

## Previous report's items: status

The seventh report's atm items, by its "Left" bullets and findings (FINDING_MAP's KA numbers):

| Item | Status | Where |
|---|---|---|
| Readability: Allard headers in F77 names, typos, `fail_9999`, "Program Terminated" banners | DONE | the module was split and rewritten (allard_atmosphere.f90, allard_table_reader.f90); no hit for any of the seventh's strings but the Sources note at atm_messages.f90:384 |
| Readability: the misspelt `latmtptau100` | DONE | O19, 1f2b14fb; allard_table_reader.f90:143-147 |
| Readability: hard-coded 54 and 5 | DONE | `old_nextgen_num_teff`, `old_nextgen_num_g`, allard_table_reader.f90:43 |
| Readability: unnamed T-tau constants | DONE | t_tau_relations.f90:42-53, :57-60 |
| Readability: dead `print_flag`/`print_to_files` | DONE | gone |
| Readability: 53 `endif`/`enddo`, 21 long lines, history comments | DONE | 0 and 0 (counts.py); Sources sections |
| Organization: lookups return through module state (KA-12) | DONE | 92d62b61; results are arguments and `allard_lookup_result` |
| Organization: no messages module, per-lookup notes (KA-9, KA-5e) | DONE | atm_messages; M235 |
| Organization: an `atm_choice`'s meaning in four places outside atm | DONE | atm_choices.f90; envelope_integration_kernel.f90:368, surface_boundary.f90:138-147, io/derived_controls.f90:328-329. Residue: `on_gray_fallback`'s literal 0 (surface_boundary.f90:33), and finding 5 |
| Organization: nothing validates the choice (KA-6) | DONE | M236; atm_lib.f90:82-86; t_tau_relations' `case default` is NaN |
| Organization: the Kurucz table read for choice 4 (KA-12) | DONE | M239; atm_lib.f90:92-102 |
| Maintainability: byte-pin test (KA-10) | DONE | 39 checks, 45 mutations |
| Maintainability: O19, both halves, and its wider reach (KA-7) | DONE | 1f2b14fb; core's side in round 2B (mass_change.f90:259-279) |
| Maintainability: an Allard composition the file lacks (KA-3) | DONE | M221; allard_table_reader.f90:239-246 |
| Maintainability: KA-5a sentinels, 5b upper range, 5d gravity search | DONE | M231, M232, M234; kurucz_atmosphere.f90:99-110, :133-158 |
| Maintainability: KA-5c the extension width | DONE (doc; the code kept by decision) | M233; `low_gravity_widths`, allard_table_reader.f90:62-66 |
| Maintainability: M112's direct run-log writes | OPEN | atm_messages.f90:84, :90 (M112's open part) |
| Extensibility: recipes that miss files (KA-11) | DONE, with a new gap | M238; README.md:88-114; findings 5 and 6 |
| Docs: no error or message conventions, no coverage, no test claims | DONE | README §3-§5; §5 stale (finding 7) |
| Docs: kurucz_atmosphere's header promise; the Allard headers (KA-13) | DONE | M232, M233, M240 |
| KA-14: commented-out statements, history, unnamed constants in atm | DONE | 0 commented-out statements; no retired/ file is needed or present |

## New findings (not in bug_sweep)

1. **atm/tables/kurucz_atmosphere.f90:155-173 (and :133-152, :99-110): no per-row log g limit**
   - **Issue:** `outside_table` refuses log g < -0.5 for the whole table. Inside a row, a log g below the row's second valid column takes "the first four" valid columns (:164) and `spline_in_log_g` evaluates the natural spline there; numerics' `splint` continues its end cubic outside the nodes (numerics/cubic_spline.f90:120-155). So a point below the row's first valid column is extrapolated, down to log g -0.5, with `ierr` 0. Above log Teff 4.5 the linear and 3-point branches (:139-150) do the same from the top two or three columns.
     - In atmk1990p00.tab, the table of all 12 Kurucz decks, 40 of the 51 rows at or below log Teff 4.5 start above log g 0: rows 12-17 (log Teff 3.80-3.88) at 0.5, 18-21 at 1.0, 22-23 at 1.5, then 2.0 to 3.5 (tables.py). atmk2004p00.tab has 46 such rows.
     - Input that triggers it with choice 3: log Teff 3.93, log g 0.5 (a post-AGB crossing, or a massive star's yellow-supergiant phase). The stencil is rows 20-23, valid from log g 1.0, 1.0, 1.5 and 1.5, so every row is extrapolated by one or two columns. With choice 5, which never switches to gray, a row at log Teff 4.49 is valid from log g 3.5 and can be extrapolated by eight columns.
     - The Allard lookup has the limit (`row_log_g_min`, allard_atmosphere.f90:125-128); the Kurucz one does not, and M234's fix ("a log g inside a row's valid range is interpolated from its neighbours") did not add it. test_atm's checks stay inside the rows.
     - Latent for the shipped decks: the 12 Kurucz decks are 0.15-1.0 Msun dwarfs and two solar models; no deck uses choice 5 (decks.py).
   - **In the F77?** Yes: surfp.f:25 tests only `TEFFL.LT.3.5D0 .OR. GL.LT.-0.5D0`, under a comment that promises to quit "MORE THAN ONE TABLE POINT FROM THE DATA", and :104-116 is "DESIRED LOG G BELOW 2ND TABLE ENTRY -USE FIRST 4 POINTS".
   - **Fix:** in `pressure_in_row`, refuse (`surface_pressure_outside`) a log g more than one column width below the row's first valid column, in all three branches; a check in test_atm on atmk1990p00's row 22 with its mutation; the limit in the header and in README §3.
   - **Kind:** A for the shipped decks (to be confirmed on the 56-deck set); a deck that reaches it stops instead of extrapolating. **Size:** S.

2. **atm/atm_lib.f90:246-251 with core/structure/surface_boundary.f90:143-149 and core/envelope/envelope_integration.f90:176-180: the switch back to the Allard tables asks log Teff only**
   - **Issue:** `atm_get_surface_pt` returns `lookup_failed` for a point above the Allard tables in log Teff, outside them in log g, or outside a stencil row's log g range (atm_lib.f90:189). The kernel then integrates the gray relation, and `atm_get` latches `s%atm_choice_in_use = atm_eddington_gray`. At the next surface boundary `on_gray_fallback` holds and `atm_reenters_table(4, log Teff)` is true whenever log Teff is below the tables' last (4.0 for Nextgen2), whatever log g is. So for a failure on log g:
     - the surface boundary switches back to choice 4, writes "note: log Teff = ... is below 4.000, the Allard tables' last log Teff: switching back to the Allard tables", and sets `recompute_envelope_triangle`;
     - the first vertex's lookup fails again: atm writes "log g = ... is outside the Allard tables ...", core writes "note: Allard table lookup out of range; switching to a gray atmosphere" (envelope_integration.f90:354-359), and the choice is gray again;
     - this repeats at every `compute_surface_boundary`: three lines to the terminal and the run log, and three new envelope integrations where an unchanged triangle needs none.
     - Input that triggers it: an Allard deck at solar [Fe/H] whose log g falls to 1.5 or below with log Teff under 4.0 (the Nextgen2 rows start at log g 3.5 or 4.0, so the limit is 1.5 or 2.0; allard.py, and expected_test_atm.out's row range "1.500000 to 6.000000"): any star taken up the red giant branch.
     - Latent for the shipped decks: the 116 Allard decks are gyrochrones, lower-main-sequence and low-mass testsuite decks, and no pinned log holds either message.
     - README §2 says the kernel integrates the gray relation "for that integration"; core's comment says "for this and all future calls" (envelope_integration.f90:177-178). Neither describes the flip.
   - **In the F77?** The flip, yes: ENVINT set `KTTAU=0` and `LTTAU = .TRUE.` on `LAlFail` for "TeffL is above Allard max, or GL is out of range" (envint.f:146-153), and SURFBC switched back on `TEFFL.LT.AL_TEFFLmax` alone (surfbc.f:68-81). The rebuilt triangle at each flip is 2026: S12's fix (67d711a7) made the switch back start a new triangle.
   - **Fix:** give the switch back the point, not only log Teff: an `atm_reenters_table(choice, log_teff, log_g)` that also asks the kernel's range and stencil steps (no message), so that a point the lookup would refuse stays gray; or let a log g failure fall back for one integration without latching. A check in test_atm (a point inside in log Teff and outside in log g does not re-enter) and one sentence in README §2. The seam is core's as much as atm's.
   - **Kind:** A for the shipped decks; B for a deck that reaches it (its triangle history changes). **Size:** S/M.

3. **atm/tables/allard_table_reader.f90:193-199, :355-381, :226; atm_messages.f90:292, :313: three unguarded paths of the Allard reader**
   - **Issue:**
     - `skip_records` reads with no `iostat` (:197; the only such read of atm's 13, counts.py). A file in the Allard layout that is shorter than its header count (the first record says 8) ends the run in "Fortran runtime error: End of file" instead of `ierr`, in either pass (:182, :188). test_atm's "empty" and "neither layout" files do not reach it.
     - `find_table_ranges` uses each row's `first_column` and `last_column` (step 3, :373-381) before it checks them (step 4, :393-401). A row with no valid entry keeps both at their default 0, so `tbl%log_g(first)` is `log_g(0)` and `tbl%log_g(last-1)` is `log_g(-1)`: out of bounds, an abort under the development flags, and otherwise a refusal by step 4 with garbage ranges in the printout. It needs a file whose records at one Teff all hold -999 (no shipped file has an explicit -999; allard.py).
     - A blank record inside the data reads as Teff 0, log g 0, [Fe/H] 0, [α/Fe] 0 with `iostat` 0 (:226, :282; a formatted read pads a short record with blanks). At the default composition it becomes a row at log Teff = log10(0) with one valid column, and the table is refused with "the row at Teff = 0 K has 1 valid log g columns". No shipped file has one.
     - `log_allard_tables` formats a refused table with `20(...)` and `1P20E16.8` (atm_messages.f90:292, :313) while `allard_max_g` is 25: a refused table of 21 or more log g columns would end in a runtime error in the report itself (an internal write cannot start a second record). The shipped files have at most 13 columns.
     - `header_format = '(i2,x,a)'` (:51): `x` without a count is a compiler extension.
   - **In the F77?** Yes for the header skip (alfilein.f:140, :235, `read(ioatma,*)`) and the order of ALTABINIT's steps; the formats are the port's.
   - **Fix:** `skip_records` returns `ierr` through `report_allard_read_failure(path, 'its header')`; step 4's row check before step 3 (or a row without a valid entry refused in step 2); a record with Teff <= 0 refused; the formats sized by `allard_max_g`; `1x`. One synthetic file each in test_atm.
   - **Kind:** A (error paths). **Size:** S.

4. **atm/atm_lib.f90:204-207 and atm_tables.f90:70-71, :74-76, :88: a dead public entry and six dead fields**
   - **Issue:**
     - `atm_allard_teff_max` is exported (atm_lib.f90:47) and listed in `PUBLIC['atm']` (tools/check_boundaries.py:81), but since `atm_leaves_table` and `atm_reenters_table` took over the switch, only test_atm calls it. Its name and its result variable say "max", and it returns `log_teff_last` (4.000 for Nextgen2), while the table's `log_teff_max` is another number (4.008774, the lookups' limit).
     - `allard_table`'s `is_old_nextgen`, `num_feh`, `num_alpha`, `feh(:)`, `alpha(:)` and `log_teff_first` are never read; `feh` and `alpha` are never written. The reader sets `num_feh = 1` "for the table's printout" (allard_table_reader.f90:247-249), which does not print it, and atm_tables.f90:74-75 says the values are "counted".
   - **In the F77?** The fields are the COMMON blocks' (alatm01, alatm02); the entry is 2026.
   - **Fix:** remove the six fields and their assignments; make `atm_allard_teff_max` private to the test's needs (test_atm can read `tbl%log_teff_last` from a table it reads itself, as its other Allard checks do) or rename it `atm_allard_last_log_teff`.
   - **Kind:** A. **Size:** S.

5. **atm/atm_messages.f90:329 and :347, against README.md:92-98 and :104: the choices are listed by hand in two messages, and the recipes do not name them**
   - **Issue:** `report_unknown_choice` ends "(0 Eddington gray, 1 Krishna-Swamy, 2 HSRA, 3 Kurucz, 4 Allard, 5 Castelli/Kurucz)" and `report_not_a_table` "(3 Kurucz, 4 Allard, 5 Castelli/Kurucz)". Both came with the rewording of 2026-11-11 (12210cb9), one day after the recipes were traced. "Add a T-tau relation" lists no edit in atm_messages, and "Add a table" names only `note_choice_switch`. A newcomer who follows either gets a run whose seventh choice works while the message for an eighth says the choices are 0 to 5; test_atm pins that text for `atm_last_choice + 1` and would be regenerated with the stale list.
   - **In the F77?** No (2026).
   - **Fix:** build both lists from `atm_choices` (a loop from the first choice to `atm_last_choice` over the names, and over `atm_choice_is_tabulated`), so that the descriptor stays the one home; or name the two lines in both recipes.
   - **Kind:** A for the recipes; B for generated text (test_atm's pinned lines only). **Size:** S.

6. **atm/tables/kurucz_table_reader.f90:26-27, :68-72, against README.md:100-103: the reader's log Teff format is one of two, chosen by row count**
   - **Issue:** `read_surface_pressure_table` takes `kurucz_teff_format` when `num_teff == kurucz_num_teff` and `castelli_teff_format` otherwise. Read as formats: the Castelli one (`1p4e16.8,18(/1p4e16.8),/`) skips the "LOG G:" label only when the list ends on a full line, so it reads a table of 76 rows, or of any multiple of 4 below it; for another count (61, say) the next read meets the label and the table is refused as unreadable, and above 76 the format reverts. The recipe "Add a table ... in the Kurucz layout" says "the reader, with its row count" (step 3) and does not name the row-count parameter's home (atm_tables.f90:43) or the format.
   - **In the F77?** Yes: SETUPS had the two reads with their two formats.
   - **Fix:** one format for any count (`(1p4e16.8)` under format reversion, then one explicit skipped record), and the `<name>_num_teff` parameter named in recipe step 2. A check in test_atm on a synthetic table of 61 rows.
   - **Kind:** A. **Size:** S.

7. **atm/README.md:36, :75-86, :31: three stale passages**
   - **Issue:**
     - §2 "At accretion" (:36) says that without the Allard tables, "or when the point is outside the tables, core takes the photosphere of the last envelope integration (O19)". Since round 2B a failed lookup takes the gray atmosphere's photosphere at the star's Teff and log g (`gray_photosphere`, core/step/mass_change.f90:259-279; bug_sweep O19, "Round 2B, core's side").
     - §5 (:75-86) gives 95.8% of 708 lines on 2026-11-10. Round 1B's seven commits of 2026-11-11 came after it: 12210cb9 alone rewrote atm_messages (469 lines changed) and changed atm_lib and both readers, and the section was not touched (`git show 12210cb9 -- src/atm/README.md`). Its list of lines not run quotes the second pass's "Impossible failure" records, a text that commit removed, and its per-file figures (atm_messages 93.6%) are the old module's.
     - §2 (:31) says the kernel integrates the Eddington relation "for that integration" after a failed Allard lookup; core latches the choice (finding 2).
   - **In the F77?** Not applicable.
   - **Fix:** the accretion sentence as mass_change.f90:259-264 has it; the coverage re-measured with its date and its list; the latch and the switch back in §2.
   - **Kind:** A. **Size:** S.

8. **Small text and naming items**
   - **Issue:**
     - `log_teff`/`log_g` in atm_lib.f90:119-123 and allard_atmosphere.f90:66 against `log10_teff`/`log10_gravity` in kurucz_atmosphere.f90:52-53 and atm_lib.f90:239, :246, for the same two quantities;
     - `atm_init` returns on an unknown choice before it resets the tables (atm_lib.f90:82-90), so "each atm_init starts clean" does not hold after a refused choice (harmless: the run stops);
     - `atm_choice_name_length` is `public` (atm_choices.f90:43) and used only in its module;
     - `find_table_ranges` starts its extremes from the literals 999 and -999 (allard_table_reader.f90:371-372), beside the named `allard_invalid`;
     - test_atm's pinned section titles keep the F77's routine names, "surfp lookups" and "kcsurfp lookups" (test/test_atm.f90:104, :114-115);
     - atm_lib.f90:72-73's comment breaks mid-line;
     - atm_choices.f90 and atm_tables.f90 hold no I/O and are in neither STAR_BLIND_FILES nor KERNEL_FILES, so nothing keeps them so;
     - core's `on_gray_fallback` tests the literal 0 for `atm_eddington_gray` (core/structure/surface_boundary.f90:33).
   - **Fix:** text and names only.
   - **Kind:** A. **Size:** S.

## What would raise each score

- **Readability, to 9.5:** one spelling of the two arguments across the facade and both kernels; findings 4, 5 and 8.
- **Organization, to 9-9.5:** finding 2 (the switch back knows what the lookup knows); one home for a table's gray switch (its log Teff and its note's wording as data of the `choices` row, or of the table); finding 4's entry; the Kurucz table's edges from the table read; atm_messages' run-log lines through `log_line` (M112).
- **Maintainability, to 9-9.5:** findings 1 and 3 fixed, each with its check and mutation; finding 2's check; M284's tables listed in the README (or filled); the coverage re-measured after round 1B.
- **Extensibility, to 9-9.5:** findings 5 and 6; both recipes traced again at the current tree, since atm_messages changed after the trace; the gray switch as one row.
- **Docs, to 9-9.5:** finding 7; the per-row log g limit (once finding 1 is fixed) and the Kurucz tables' holes (M284) in §3; the date of the mutation run and of the coverage on the same tree.

## Traces, counts and README checks (supporting)

- **Decks** (decks.py over the 1,198 `.inlist`/`.nml1`/`.nml2` files of `git ls-files`): `atm_choice` 0 in 1,048 (and `kttau` 0 in 17 legacy files), 4 in 116 (92 gyrochrones, 18 lowerMS, 6 testsuite), 3 in 12 (10 lowerMS dwarfs, Test_solarGS98_kurucz, Test_solarM22M_difkur), 1 in 1; 2 and 5 in none. No deck sets `allard_target_feh`, `allard_target_alpha` or `allard_use_tau100`. So test_atm remains the only coverage of HSRA, Castelli/Kurucz and the tau = 100 point, as README §4 says.
- **Shipped tables** (tables.py, allard.py): 35 Kurucz and 16 Castelli/Kurucz tables read with the reader's formats; Nextgen2.all holds nine [Fe/H] at [α/Fe] 0 (62 rows and 6 log g columns at solar), AMES-Cond2.all five, AMES-Dusty2.all one; no Allard composition has a row of fewer than four columns or a hole.
- **Mutation manifest** (counts.py): 45 entries over 7 files (atm_lib 11, t_tau_relations 10, kurucz_atmosphere 8, allard_table_reader 8, allard_atmosphere 5, atm_choices 2, kurucz_table_reader 1); each `old` text found once; 39 labels, each with at least one mutation.
- **"Add a T-tau relation", traced by reading:** atm_choices (the row), t_tau_relations (three cases), atm_lib (re-export), io/controls_buffer.f90:853-856 (the description), test_atm and its manifest, the README; plus atm_messages.f90:329, which the recipe does not name (finding 5). Nothing in core/ or io/'s code: the kernel asks `atm_choice_is_tabulated` (envelope_integration_kernel.f90:368), the banner `atm_choice_name` (io/derived_controls.f90:328-329).
- **README claims checked and true:** §1's nine rows against the files and their public names; §2's start-of-run steps against atm_lib.f90:76-103 and setup/settings_install.f90:77; §2's envelope and surface-boundary steps against the kernel and `switch_atmosphere_choice`; §3's `ierr` values against atm_lib, both kernels and both readers; §3's destinations against atm_messages (`atm_report` for errors, fallbacks and notes, `atm_log` for the reader's findings); §4's 39 checks and 45 mutations; §4's deck list against decks.py.
- **Boundaries checked and true:** `PUBLIC['atm']` matches atm_lib's procedure exports; outside atm only core, setup, io's banner and test programs use atm_lib; `atm_settings_t`'s four fields and `atm_table_paths_t`'s two are each copied by the builder, which is pure and carries its `star reads:`/`star writes:` lines.
- **Open bug_sweep items of atm:** M284 (data) and M112's open part. O19, M221, M231-M240 and M283 are fixed in the code as their entries say.
- **Not checked:** anything that needs a build or a run: the mutation verdicts, the coverage figure, the development-flags run, and the size of finding 1's extrapolation error.

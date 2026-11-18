# eos -- eighth evaluation (2026-11-17, yrec-next 1179afbf, read-only)

Scope: every production file of src/eos/ read whole: README.md, eos_lib, eos_configuration, eos_results, eos_messages, eos_mixture, eos_table_source, eos_source_list, eos_test_entries, blended_eos, the four table modules, yale/yale_eos and fermi_table_reader, scv/scv_eos and scv_table_reader, opal/opal_eos and opal_table_reader, mhd/mhd_eos and mhd_table_reader, and retired/opal95_radsub_notes.f90.txt. Also read whole: setup/eos_settings_builder.f90, the repository root's test_eos.py, io/derived_controls.f90's `select_eos_and_opacity_tables`. Of eos/test/test_eos.f90 (1,991 lines) I read lines 1-925, 1290-1410, 1525-1700 and 1772-1832 and every `call check(` line; the rest was skimmed. Skimmed: expected_test_eos.out (its 76 ok lines), test_eos.mutations.toml (checked by script, below), tools/check_boundaries.py's eos entries, tools/check_style.py's header and tools/style_counts.tsv's eos rows, tools/devel_tests.py, .github/workflows/ci.yml, docs/GUIDELINES.md ("Naming", "Declarations and interfaces", "Code style", "A package README"), docs/ROADMAP.md rounds 0, 1 and 1B, docs/TARGET_ARCHITECTURE.md 4.5(c), and docs/bug_sweep.md's eos entries. The F77 was read with `git show 6cd56739:src/<file>.f` (eqbound06, eqstat2, rabu, mhdst1, mhdtbl). Since the seventh evaluation (823326a7) 35 commits touched eos, its builder or its test: round 1's byte-identical ones, round 1B's twelve measured ones, the integrator's, and M293's in round 2B. Nothing was built or run. The counts come from scratch Python over the tree (SP/review8_eos/: counts.py, manifest.py, unused.py, scv_rows.py, decks.py, with their .out files). Paths are under src/ unless they start at the repository root.

What I checked of the rounds' claims, against the code and not the ROADMAP:
- the manifest: all 82 mutations' `old` text is in the named file exactly once at 1179afbf, each of the pin's 76 ok lines has at least one mutation, and no mutation names a check that is not in the pin (manifest.out). I did not run them;
- the fixes of M93, M105, M106, M107, M136, M137, M144, M206-M211, M88, M92, M94, M128, M129, M290-M293, M305 and M306 are in the code where bug_sweep says (lines in the status section);
- `opal95_eos.f90` is gone and the three vintages run one code; `eos_mix` is protected; `blend_eos` is private; no eos entry is left in NO_IO_ALLOWANCE; devel_tests.py's KNOWN_FAILURES is empty and CI has a devel job;
- the README's deck figures hold: 310 decks, 304 with OPAL 2006 + SCV + Debye-Hueckel, 2 SCV only, 4 Yale only (decks.py);
- the shipped SCV H table's rows never get shorter with T (30 to 76 points, 63 rows, every step 0.2 in log P_gas and 0.08 in log T; scv_rows.py), which M88's edge weight and M206's guard rely on.

## Scores (1-10)

Fourth -> seventh -> now:

- **Readability 7.5 -> 7.5 -> 8.** The package no longer reads as F77: the style ratchet's ten counts are zero in every eos directory but the readers' I/O, the messages are sentences that start with where they come from, each solver file ends in a Sources section, the Saha solver's Jacobian has names (`jac_he2_x2`, `elim_he1_by_he2`), the Fermi-Dirac rows and the OPAL sentinels are named, and OPAL is written once. The kernels' insides are the point lost. Left:
  - 11 procedures over 100 lines, 1,908 lines in all (15 and 2,763 at the seventh): `interpolate_opal_eos` 318 (opal/opal_eos.f90:496-813, with loop bodies at three indentations), `compute_scv_envelope_eos` 307 (scv/scv_eos.f90:194-500, the four-way smoothing still written out at :298-389), `compute_scv_eos` 187, `readcoeos` 179, `opal_density_search` 153, `mhdst1` 146 with 36 arguments, `iterate_fully_ionized_density` 145;
  - the fully ionized solver's `pr6..pr9` and `cden1..3`, each reused for two meanings (yale/yale_eos.f90:858-944), `fermi_solution`'s 40 fields (`cl1`, `d8`, `emul8`, `ttdcu`; :141-185) and the Saha sensitivities `sqet`, `sx1qt` (:117-135), explained only by a comment;
  - about 15 thresholds still literals, most in the OPAL search and stencil: 1e-7 four times, 1e-6 twice, 1e-15, 1.25, 0.2 and 5, 1e-14 twice, 0.5e-7 (opal/opal_eos.f90:354-418, :560-676); `ion_mean_weight_excess`'s 1e-5 (eos_mixture.f90:108);
  - 17 of 140 procedures are `pure` (10 of 104 at the seventh); about 25 more have no side effects (the descriptor's standard procedures, `zone_mixture`, `blend_saha_into`, Yale's derivative steps, `find_scv_cell`, `scv_ionization_fractions`);
  - 76 production lines over 100 characters, 12 over 132 (9 of them in mhd/mhd_table_reader.f90:128-316);
  - the test keeps F77 declarations the ratchet does not count: 24 `double precision` lines, five `data` statements, `use eos_test_entries` without `only:`, and compositions filled by literal slot (test/test_eos.f90:35, :53-100, :404-409, :1867-1872);
  - findings 6 and 7 (stale comments, unused imports).
- **Organization 8 -> 8 -> 8.5.** The seventh's list is done but for one item: a table-source descriptor and one ordered list that `eos_init`, `eos_get` and `blend_eos` loop over; SCV a source with its own edge weight; OPAL 1995 on the shared code (about 1,000 duplicated lines gone); readers apart from the kernels; the test's entries in `eos_test_entries`; `eos_mix` protected; `blend_eos` private. Left:
  - the solvers still format and send their own messages: 23 `eos_report`/`eos_log` calls in yale_eos, scv_eos, opal_eos, mhd_eos, blended_eos and opal_eos_tables (13 of them in opal_eos). "Kernels print nothing" holds in check_no_io's sense (no I/O statement), not in the sense kap's lookups now meet (a reason returned, reported once by kap_lib). This is the seventh's "failures as codes, reported in `eos_get`";
  - the tables are public, writable module state (`tablex`, `tablenv`, `last_scv_cell`, `opal_tables`, `opal_shared`, `yale`, `mhd`): only convention keeps their writers to the readers, and the test writes `tablenv` and `last_scv_cell` (test/test_eos.f90:1368, :1392). `find_ragged_edge_ramp` still starts from the stencil the last interpolation left in the instance (opal/opal_eos.f90:953, :962);
  - `opal_eos_tables` holds data, three procedures and two messages (`find_opal_ramp`, :470-564);
  - the Yale base and the fully ionized gas are `blend_eos`'s own steps, not entries (TARGET_ARCHITECTURE 4.5(c) says so).
- **Maintainability 7 -> 7 -> 8.5.** Every latent defect the seventh listed is fixed and tested; non-convergence is an error on both Yale solvers; 76 ok/FAIL checks of behaviour (limits, continuity at both ends of four ramps, exactness on synthetic MHD tables, about twenty error contracts) with a manifest of 82 mutations that matches the code; test_eos runs under the development flags in CI. Left:
  - open in bug_sweep: O2 (second derivatives inside a table are Yale's), O18's OPAL 2001 search edge, M157's MHD pressure question, M307 (six shipped OPAL 2006 tables the reader refuses), M308; all but M307 wait for a decision;
  - the OPAL code's coverage is 74.8% and its reader's 67.2% (README.md:138); the figure is dated 2026-11-10, before round 1B's twelve commits;
  - findings 1 and 2: the MHD reader's one anonymous message and its read past the composition arrays, and two table-shape assumptions no reader checks;
  - finding 3: SCV's two paths give different numbers at the composition where one hands over to the other;
  - two OPAL diagnostics log and go on (`report_order_too_high`, `report_t6_search`; opal/opal_eos.f90:540, :650); no production call reaches them.
- **Extensibility 6.5 -> 6.5 -> 8.5.** A new table source is now a module, one list entry and its settings, and `blend_eos`, `eos_get` and `eos_init` need no edit; a vintage is a row of `opal_layout`. I traced both recipes by reading. Left:
  - finding 4: the source recipe gives `edge_weight` the signature it had before M88 and leaves out the style table's row for the new reader; the vintage recipe leaves out `eos_lib.f90`. Each stops the newcomer at the compiler or at `make check`, none silently;
  - a new base solver (a second Saha, another degenerate gas) is still an edit of `gas_eos`, with no recipe.
- **Docs 7.5 -> 7 -> 8.5.** The README meets the checklist: a table with a row for each of the 21 production files, the flow in three numbered times, conventions (units, errors, messages, controls, kept by decision), the 76 checks by kind, the mutation count with its date, a dated coverage command, three dated recipes, and the per-path table the fourth and seventh reports asked for. Every figure I re-counted holds (76 checks, 82 mutations, 304 of 310 decks, the Fermi range, the ramp width, the step sizes). Left:
  - finding 5: the per-path table's 1/mu_ion row still gives the rule M129 removed;
  - finding 4's three recipe gaps, and the recipes' trace date (2026-11-10) is before the interface change;
  - finding 6: about fifteen stale comments and two wrong section numbers;
  - "Kept by decision" names O2 only: finding 3's two SCV interpolations and O18's 2001 search edge are kept too;
  - three paragraphs over 900 characters (README.md:92, :118, :138; rotation's longest is 745).

Average **8.4** (7.2 at the seventh evaluation). No criterion is under 8; readability, at 8, is the lowest, and its "Left" list is a short round's work.

## Summary

- **eos is at the standard.** Thirty-five commits since the seventh evaluation did what that report's five "what would raise" lists asked, almost item by item: the source descriptor, OPAL 1995 on the shared layout, the readers out of the kernels, the test's entries in one module, the latent defects fixed (18 bug_sweep items in round 1, five of them texts; 11 measured ones in 1B; M293 after), messages reworded, history in Sources sections, a mutation manifest, the devel build in CI, and a README rewritten to the checklist.
- **The rounds' claims hold where I checked them.** The manifest's 82 mutations all apply to today's code and cover all 76 checks; each fix is where bug_sweep says; the README's counts are right.
- **What keeps it from 9:**
  - four kernels keep their F77 shape inside (`interpolate_opal_eos`, `compute_scv_envelope_eos`, the fully ionized density iteration, `mhdst1`);
  - the solvers send their own messages (23 calls);
  - the recipes and the per-path table were not re-traced after round 1B changed the descriptor's interface and the 1/mu_ion rule.
- **New, all latent or documentation:**
  - the MHD reader reads a composition record past its arrays before it checks the count, and its 17 failure sites share one message that names neither the file nor the check (F77);
  - no reader checks that the SCV rows' point counts never fall with T, or that the first OPAL block's density columns share one T6 grid, which the lookups assume (F77);
  - SCV's envelope-table path and its general path interpolate the same tables two ways and differ by up to 5e-4 dex in log rho, by the test's own tolerance, at the 1e-5 composition handover; it is in no list (F77);
  - three recipe gaps, one stale table row, about fifteen stale comments, 22 unused imports.

## Previous report's items: status

**The seventh report's new findings:**
1. SCV's downward smoothing reads unbuilt entries (M206): **DONE.** The guard, scv/scv_eos.f90:295-302; its check compares the lookup itself with the unbuilt entries filled (test/test_eos.f90:1360-1379).
2. OPAL 1995's pressure check prints two bare numbers (M207): **DONE** (opal/opal_eos.f90:192-201). The text hard-codes "OPAL 1995 EOS" and "5e-7" beside the row's `eos_name` and `pressure_check_tolerance` (finding 7).
3. The Fermi-Dirac table's top and degeneracy range unchecked (M208): **DONE.** The range is named (yale_eos_tables.f90:38-43), `compute_fully_ionized_eos` returns `extrapolated` (yale/yale_eos.f90:732-734), `fully_ionized_gas` logs it once a run (blended_eos.f90:410-416), README section 3 states it.
4. The second X quadratic checks the first's value (M209): **DONE** (opal/opal_eos.f90:752), and since M292 on the magnitude.
5. The vintage recipe's step 5 (M210): **DONE** (README.md:167). The recipe has a new gap in step 1 (finding 4).
6. test_eos.py's docstring (M211): **DONE** (test_eos.py:4-19); it now cites the wrong README section (finding 6).
7. One value written twice, the 1e35 family (plan item E-7): **DONE.** `opal_empty_entry` and its four thresholds are parameters with their meaning (opal_eos_tables.f90:92-105); the ramp test reads the layout (opal/opal_eos.f90:278); the search's sentinel is compared by name (:170).
8. Declarations and the test's slots (plan item E-8): **PARTLY.** Explicit-shape dummies are gone and `mhdst1`'s arrays have intents; 17 of 140 procedures are `pure`; the test's `double precision`, `data`, `use eos_test_entries` and literal slots are as they were (lines under Readability).

**The fourth report's findings, carried open by the seventh:** 1 (M148) **DONE** (README.md:148). 2 (M129) **DONE** (eos_mixture.f90:160-191; test/test_eos.f90:1167). 3 (M128) **DONE** (scv/scv_eos.f90:824). 4 (M136) **DONE** (scv/scv_table_reader.f90:85, :103; opal/opal_table_reader.f90:133). 5 (M137) **DONE** (eos_lib.f90:134). 6 and 7 (M157) **DONE** for the texts; the MHD pressure question is open, a decision. 8 (M156) **DONE** (README.md:63-75; finding 5 is one row of it). 9 (M144) **DONE** (test_eos.py:65-117). 10, style residue: **DONE** (retired/opal95_radsub_notes.f90.txt with its header; Sources sections in all six solver files and the four table modules); one "(was label 10)" is left (yale/yale_eos.f90:855).

**What would raise readability:**
- Name the Fermi-Dirac rows, split `fermi_solution`, replace `pr1..pr9`: **PARTLY** (rows `i_fd_*`, yale_eos_tables.f90:36-37; the rest open).
- Readable names for the Saha sensitivities: **OPEN** (yale/yale_eos.f90:117-135); the Jacobian and the residuals got names.
- `interpolate_opal_eos` in named steps; SCV's smoothing as one helper: **PARTLY** (`clamped_coverage`, the `report_*` procedures and `spline_cell` are out; the bodies are 318 and 307 lines).
- `mhdst1` with a record and assumed-shape arrays: **PARTLY** (assumed shape; still 36 arguments, mhd/mhd_table_reader.f90:177-182).
- History and commented-out code out of the bodies: **DONE.**

**What would raise organization:** SCV with an edge weight, MHD a source: **DONE** (scv/scv_eos.f90:96-141; mhd/mhd_eos.f90:49-60). The descriptor: **DONE** (eos_table_source.f90, eos_source_list.f90). Test-only entry points out of production modules: **DONE** (eos_test_entries.f90). `eos_mix` protected: **DONE** (eos_mixture.f90:53). Search state in a record (M107), OPAL 1995 on the shared layout (M94): **DONE** for both defects; the SCV cell is a typed module variable and OPAL's edge ramp still starts from the instance's last stencil.

**What would raise maintainability:** M92, M93, M94, M105, M106, M107 fixed: **DONE**, all six. Kernels' failures as codes, reported once: **PARTLY** (every failure returns ierr; the kernels still write the messages). One 1/mu_ion rule: **DONE.** test_eos under the devel flags in CI: **DONE.** OPAL coverage, M119, M87: **PARTLY** (M119 and M87 done; OPAL 71% -> 74.8%). A committed mutation manifest: **DONE.**

**What would raise extensibility:** the recipes corrected (M148, M123, finding 5): **DONE**, with finding 4's new gaps. The descriptor: **DONE.**

**What would raise docs:** the listed items: **DONE.** The per-path table: **DONE** (one stale row). The 1,256-character paragraph: **DONE.** Coverage re-measured with a pin change: **PARTLY** (re-measured in round 1, 2026-11-10; the pins moved again in 1B).

## New findings (not in bug_sweep)

1. **mhd/mhd_table_reader.f90:346-351** (`rabu`), **:45-51** (`read_mhd_tables`' one message), and the 17 `ierr = 1` sites of `mhdst`, `mhdst1`, `rabu` and `rtab` (:158-400)
   - **Issue:** two things about a bad MHD table.
     - `rabu` reads the element count and, in the same statement, that many triples: `read(...) num_chem, (atomic_weight(ic), number_abundance(ic), mass_fraction(ic), ic=1,num_chem), mean_molecular_weight`. The test `nchem0 < num_chem` comes after. A table whose header says 7 elements or more writes past the three 6-element columns it was handed (`mhd%zams_atomic_weight(:,ic)` and its neighbours, so into the next composition's column or the next component of `mhd`) before `ierr = 1` comes back; under `-fcheck=all` it is a bounds stop.
     - Whatever fails, the user reads "eos_init: an MHD table is corrupt or does not match its layout": not which of the eight files, and not which of the checks (variable count, table kind, composition flag, the X spacing of the companions, their T grids, the sizes, a short record). `read_mhd_tables` has the paths and loops over the files inside `mhdst`, so the name is one argument away.
   - **In the F77?** Yes to both: `READ(IR) NCHEM,(ATWT(IC),ABUN(IC),ABFRCS(IC),IC=1,NCHEM),GASMU` then `IF(NCHEM0.LT.NCHEM) STOP` (rabu.f), and bare `STOP`s in mhdst1.f whose FORMATs nothing writes. Latent: no MHD tables ship.
   - **Fix:** read the count alone, check it, `backspace` and read the record; `mhdst1`, `rabu` and `rtab` return a reason code (numerics' named-code pattern) and `read_mhd_tables` writes one line with the file and the reason. A test_eos check with a 7-element synthetic header.
   - **Kind:** A. **Size:** S.
2. **scv/scv_table_reader.f90:75-107** against scv/scv_eos.f90:752-753, :449-457, :572-574, :116; and **opal/opal_table_reader.f90:229-235, :140**
   - **Issue:** the lookups assume two things about a table's shape that no reader checks (M136 and M291 added the checks across tables and across X blocks; these are within one table).
     - SCV: `find_scv_cell` bounds a cell's columns by the count of its first row, "the one with the fewest points" (:752), and the lookups then read rows idtt+1 to idtt+3 at those columns, one more with the upward pressure shift. That holds only if the rows' point counts never fall as T rises. The reader checks 4 <= npts <= 76 and that the He and Z tables share the H table's grid, not that `nptsx` is non-decreasing, that `tlogx` increases, or that the step in log P_gas is the 0.2 `scv_edge_weight` takes from the first row. A table with a row shorter than the one below it gives a density from entries never read (zero, or an earlier table's after a second `eos_init`), with ierr = 0. The first column is also written twice: `scv_min_log10_gas_pressure = 4.0d0` (scv/scv_eos.f90:44) for coverage and `tablex(1,1,iscv_log10_p)` (:118) for the ramp.
     - OPAL: the interpolation uses one T6 grid, column 1's (`v%t6_grid(k) = v%t6_list(1,k)`, :140). A later block's rows are compared with the block before, column by column, but within the first block column j's T6 rows are never compared with column 1's.
   - **In the F77?** Yes (SETSCV and READCOEOS check neither). The shipped tables satisfy both (scv_rows.out for SCV); only a user's table can trigger it.
   - **Fix:** in `read_scv_rows`, ierr = 1 with the row named when `nptsx(t) < nptsx(t-1)`, `tlogx(t) <= tlogx(t-1)`, or a row's first log P_gas or its step differs from the first row's; coverage reads the table's first column. In `read_table_blocks`, compare each T6 read with `t6_list(1,t6_row)` once column 1 is in. One test_eos check each, as M136's and M291's.
   - **Kind:** A. **Size:** S.
3. **scv/scv_eos.f90:235-240** (`compute_scv_envelope_eos`'s handover to `compute_scv_eos`), with test/test_eos.f90:1793-1795 and :1816-1817
   - **Issue:** a zone within 1e-5 of the envelope's X and Z gets the envelope table through a natural spline with smoothing across cells; a zone farther off gets the three tables through 4-point Lagrange interpolation without smoothing. The two are not the same function of (T, P): the test that guards the envelope table accepts them 5e-4 dex apart in log rho and records 3e-4 at log T 4.13 ("agrees with SCV's general path to 5e-4 in log rho"). So the density steps by that much where a zone's composition crosses the tolerance, and cp and grad_ad are built differently on the two sides (the table's du/dT against the entropies). I did not measure it; the figure is the test's. It can matter where a zone under SCV differs from the surface: below a settling envelope in a deck without OPAL, or outside the OPAL table.
   - **In the F77?** Yes (EQSCVE hands over to EQSCVG at the same tolerance).
   - **Fix:** first say it: README section 3's "Kept by decision" and a bug_sweep entry, with the measured size over the SCV-only decks. Then a decision: one interpolation for both paths (the general path through the same spline and smoothing, the envelope table kept as its cache), a measured change.
   - **Kind:** A for the statement; B for the fix. **Size:** S; M.
4. **eos/README.md:149, :152, :163-166** (the "Add a table source" and "Add an OPAL vintage" recipes)
   - **Issue:** three places where a newcomer who follows the recipe is stopped.
     - Step 1 gives `edge_weight(temperature, log10_density, weight, in_table, near_edge, ierr)`. Since M88 (a0039581, two days after the recipe's trace) the descriptor's interface is `(log10_temperature, table_r, weight, in_table, near_edge, ierr)` (eos_table_source.f90:68-76), so the procedure as written does not match the pointer.
     - Step 1 puts the reader in a new directory, `eos/<name>/<name>_table_reader.f90`, with `open`, `read` and `close`. check_style allows a directory missing from `tools/style_counts.tsv` zero I/O statements ("a new directory starts clean"), so `make check` fails until a row `eos/<name>  io_statement  N` is added by hand. The recipe names the Makefile, the dependency generator and STAR_BLIND_FILES, not this file. The trace's toy source had no reader.
     - The vintage recipe's step 1 adds `opal_vintage_<year>` in `eos_configuration.f90` only. `io/derived_controls.f90:44-45` takes the constants from `eos_lib`, which re-exports them by name (eos_lib.f90:24-26, :41-42), so step 4's branch does not compile until `eos_lib.f90`'s two lists and derived_controls' `use` line have the new name. The recipe says steps 4 and 5 were "checked ..., not traced".
   - **Fix:** the three corrections, then trace both recipes again on a throwaway branch with a source that has a reader and a vintage that reaches `derive_controls`, and re-date them.
   - **Kind:** A. **Size:** S.
5. **eos/README.md:70** (the per-path table, row "1/mu_ion")
   - **Issue:** the OPAL cell reads "excess metals as carbon (M129)" and the MHD cell "excess metals as carbon". M129's fix (15b8b199) made both take `zone_mixture`'s rule through `compute_mean_weights` (the envelope's metals scaled under diffusion of Z, else the excess as carbon; eos_mixture.f90:178-184). That commit rewrote section 2.1's step 1 and the "Kept by decision" bullet and left this row, so the table now contradicts the paragraph above it and the code.
   - **Fix:** both cells "`compute_mean_weights`: `zone_mixture`'s rule".
   - **Kind:** A. **Size:** S.
6. **Stale comments and references** (GUIDELINES code style 5)
   - **Issue:**
     - eos_messages.f90:6-8: `eos_log` is "a Saha or density iteration that did not converge"; both are `eos_report` errors since M93 and M293.
     - opal/opal_eos.f90:261 and :889: "blend_eos's opal_step calls it"; there is no `opal_step` (it is `apply_table_source` through `opal_edge_weight`). :13-15 lists `readcoeos` among the module's own procedures.
     - opal_eos_tables.f90:64, :574, :580: `opal95%eos_output`, `opal95%quad` (the instance is `opal_tables(iv_opal95)`); :570: "esac*.f90/t6rint*.f90"; :304-308: one sentence written twice.
     - mhd_eos_tables.f90:9 and mhd/mhd_eos.f90:5-7: `read_mhd_tables` as mhd_eos's (it is mhd_table_reader's); mhd/mhd_table_reader.f90:358: "an unformatted opacity/EOS-style table".
     - yale/yale_eos.f90:1144: a comment cut in mid-parenthesis; :855: "(was label 10)".
     - test_eos.py:12 and test/test_eos.mutations.toml:2: "eos/README.md section 4" for the checks, which are section 5.
     - test/test_eos.f90:17-26 and :769-798: the headers list two and twelve check groups of the fourteen, and `settings_checks` as "eos_init's check of the OPAL vintage" alone.
     - tools/check_boundaries.py:67-69: the test-only entries "from their own modules" (they are eos_test_entries').
     - eos/README.md:91: only SCV's spline is said to return numerics' status; MHD's `lir` status comes back too (mhd/mhd_eos.f90:274, :325-326).
   - **Fix:** the texts.
   - **Kind:** A. **Size:** S.
7. **Unused imports, and message texts that repeat a named value**
   - **Issue:**
     - 22 names are imported and not used (unused.out): `eos_log` in blended_eos, yale_eos, scv_eos and mhd_eos; `debye_huckel_composition_of` in blended_eos; `n_mix_species` in yale_eos; `i_opal_p` and `opal_empty_entry` in opal_eos; 15 slots and sizes in mhd/mhd_eos.f90:17-27. An `only:` list is the module's stated dependence, so these overstate it.
     - Four texts carry a value that is a parameter beside them: "OPAL 1995 EOS" and "5e-7" (opal/opal_eos.f90:195-196, against the row's `eos_name` and `pressure_check_tolerance`); "log T 5.2 to 9 ... from -16 to -3" (blended_eos.f90:412-413, against `fermi_log10t_min`/`_max` and the table's abscissae); "beyond 1e20" (opal/opal_eos.f90:1325, against `opal_x_interp_overflow`). eos_lib.f90:120 says "more than one OPAL EOS vintage selected" for any vintage out of range.
   - **Fix:** drop the imports; write the values into the texts from the parameters (a re-pin of the test's run-log lines where they are pinned).
   - **Kind:** A. **Size:** S.

## What would raise each score

**Readability 8 -> 9**
- `interpolate_opal_eos` as its steps: the X stencil, the T6 and density stencils, the coverage and orders, the X interpolation of one variable (the two quadratics once), the scaling and radiation; one indentation.
- `compute_scv_envelope_eos`: the smoothing as one helper applied in T and then in P, carrying M206's guard; SCV's radiation correction (dln rho/dln P, dln P/dln T, du/dT, cp, grad_ad) once for both lookups.
- The fully ionized solver: names for `pr6..pr9` and `cden*` by what they hold (the inverse spacings, the Newton terms), and `fermi_solution` cut into its stencil, its pressures and its derivatives.
- `mhdst1` with a record for one table and one for the companions; its lines under 132 characters.
- Name the OPAL search's and stencil's tolerances in opal_eos_tables beside the sentinels.
- `pure` on the 25 side-effect-free procedures; the test in `real(dp)`, with `only:` and composition_layout's slots.
- Findings 6 and 7.

**Organization 8.5 -> 9**
- The solvers return a reason (a named code and the numbers) and `eos_get`, `apply_table_source` and `eos_init` report it once through eos_messages, as kap_lib does for its lookups; `opal_eos_tables` then needs no message.
- The table state private to each family's tables module with a reader entry and read-only access (`protected`), and the test's two pokes (`tablenv`, `last_scv_cell`) through eos_test_entries procedures.
- `find_ragged_edge_ramp` given the row and column to start from, or searching from the grid's ends.

**Maintainability 8.5 -> 9**
- Findings 1 and 2, each with its check and mutation.
- Finding 3's statement and measurement; then the decisions bug_sweep holds for eos (O2, O18's 2001 edge, M157, M307's six tables: remove them or mark them).
- OPAL's error branches under synthetic corrupt tables (a block out of order, a short column, a mass fraction above 1), to bring opal_eos and its reader past 85%; the coverage figure re-measured and re-dated at 1179afbf.
- `report_order_too_high` and `report_t6_search` return ierr.

**Extensibility 8.5 -> 9**
- Finding 4, with both recipes re-traced after round 1B.
- A fourth recipe, "Add a base solver", or the Yale base and the fully ionized gas as entries of the list.

**Docs 8.5 -> 9**
- Findings 4, 5 and 6.
- "Kept by decision": the two SCV interpolations (finding 3) and OPAL 2001's search edge (O18), each with its reason.
- The three paragraphs over 900 characters split.

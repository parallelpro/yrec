# kap -- eighth evaluation (2026-11-17, yrec-next 1179afbf, read-only)

Scope: every file of src/kap/ read whole:
- README.md, kap_lib, kap_messages, opacity_tables, opacity_source, opacity_families, opacity_blending, opacity_table_files;
- the seven family directories, each lookup and reader (opal95/, opal92/, laol89/, alex06/, alex94/, kurucz90/, conductive/), and retired/conductive_opacity.f90.txt;
- test/test_kap.f90, test/expected_test_kap.out, test/test_kap.mutations.toml, and the repository root's test_kap.py;
- setup/kap_settings_builder.f90 and setup/settings_install.f90.

Read for the claims they touch:
- kap's entries in tools/check_boundaries.py (PUBLIC, STAR_BLIND_FILES, KERNEL_FILES, STAR_FREE_DIRS, FACADE_ONLY_DOMAINS, the installer rule) and tools/style_counts.tsv;
- the callers (core/step/shell_physics, core/structure/shell_equations, core/envelope/envelope_integrands and envelope_integration_kernel, core/model/starting_model, model_rescaling and rezoning, core/run/run_setup, convection/semiconvection), io/derived_controls.f90's opacity parts, kap_messages' models (rotation_messages, eos_messages, state/run_messages);
- docs/GUIDELINES.md "Naming", "Declarations and interfaces", "Code style" and "A package README"; docs/ROADMAP.md rounds 0, 1, 1B and 2B; docs/TARGET_ARCHITECTURE.md 4.5(c) and its table of F77 homes; docs/bug_sweep.md (every kap entry, M219-M240, M277-M282, M298-M300, M112, O18);
- the seventh report (audit/module-evaluation-7-2026-11-10/kap_atm.md) and FINDING_MAP.md.

How the evidence was gathered:
- every defect below was read twice in the code and checked against `git show 6cd56739:src/<file>.f`;
- nothing was built or run, so coverage, the mutation verdicts and the devel run are the README's claims, not mine;
- counts come from scratch Python over the tree (SP/review8_kap/): counts.py (I/O sites, `kap_report` sites, long lines, procedure lengths, `pure`), checks.py (each ok line of the pinned output against the manifest's mutations) and unused.py (only-list names used once).

Paths are under src/ unless they start at the repository root. Another reviewer has src/atm.

## Scores (1-10)

- **Readability 6.5 -> 7.5.** kap_lib, opacity_source, opacity_families, opacity_blending and kap_messages read at the bar: `kap_eval` is four named steps (kap_lib.f90:110-146), each formula is once in opacity_blending with its coefficients named, and every file has a header and a Sources section. The style ratchet is at zero outside the readers' I/O. The seven lookups underneath are the F77 routines in modern syntax. Left:
  - `compute_opal95_opacity` is 282 lines with a nested procedure, and its nested `if` blocks are not indented (opal95_opacity.f90:270-551, :336-377, :419-431); `compute_kurucz_opacity` is 158 lines with two headers, one above and one inside (kurucz_opacity.f90:95-114); `compute_conductive_opacity` is 157 lines with the main branch written twice (conductive_opacity.f90:149-188 against :193-224);
  - opal95_reader is the one reader not on `read_failed`: 3 `goto` and 7 reads that jump to labels 9997-9999 by `end=`/`err=`, around an unindented `table_loop` (opal95_reader.f90:52-202); alex06_reader repeats its 40-line table body for the X = 1 - Z tables (alex06_reader.f90:53-100 against :107-145);
  - F77 names in the lookups (23 lines with `o0`, `qod0`, `chkd`, `kk2` or `ii1`; also `chko`, `mm1`, `ndss`, `ndf`, `im2`, `im3`, `grdnt`, `rhot3`, `jerr_gate`), and `grid_temp_k` for a T6 (opal92_reader.f90:37, :76-79);
  - internal procedures that work through the host's loop indices and work arrays (opal92_opacity.f90:72-112, :168-222; laol_opacity.f90:74-104; the readers);
  - three spellings of 10**x: `exp10` (opal92, laol, kurucz, conductive), `pow(1.0d1, x)` (opal95_opacity.f90:577, :601, :633) and `exp(ln10*x)` (alex06_opacity.f90:257; alex94_opacity.f90:208, :243);
  - unnamed values: OPAL92's surface slot 4 (opal92_opacity.f90:263, :359, :382; alex94 has `alex94_surface_slot`), its `-1.0d0` (:278), Kurucz's 25 (kurucz_opacity.f90:75), Alexander 2006's 0.1 (alex06_opacity.f90:92), `1.0d-6` (opal95_opacity.f90:419; conductive_opacity.f90:151-155, :195-199), Hubbard & Lampe's 0.6 (opacity_blending.f90:122), `np = 100` twice, `dummy(104)` (laol_reader.f90:162), and `len=256` 31 times;
  - 15 of 113 procedures are `pure`; most lookups could be (finding 10);
  - 9 of the 23 "(formerly ...)" banners give an October 2026 name, not the F77's (alex94_opacity.f90:35, :84, :116; opal92_opacity.f90:35, :117, :227, :322; opal95_opacity.f90:68, :183);
  - ragged indentation in opal92_reader.f90:63-83 and laol_opacity.f90:185-209; 14 lines over 110 characters.
- **Organization 7 -> 8.5.** One facade (enforced), star-free, one settings and one paths record from setup/, every family a `kap_source` descriptor in one precedence list, the state private to kap_lib and passed to the families as arguments, readers apart from lookups, one messages module, the seven lookups and opacity_blending in STAR_BLIND_FILES and kap_lib and the type modules in KERNEL_FILES. Left:
  - kap takes the envelope Z once, at `kap_init`, and has no entry to follow it (finding 1);
  - the molecular precedence is also coded in io (finding 8);
  - `opacity_table_state` is one central type that must list every family's table set, in order of size for a Mach-O limit (opacity_tables.f90:272-295): a family cannot own its tables;
  - conductive_reader opens its file itself, beside `open_table` (conductive_reader.f90:36-43);
  - laol_reader carries the header's metal mixture out through two argument lists into a local that `load_laol89` discards, and uses state's composition_layout for it (laol_reader.f90:30-33, :104-107; laol_opacity.f90:436-445);
  - `kap_settings` is a public, protected module variable of opacity_tables that only kap_lib needs (opacity_tables.f90:329).
- **Maintainability 5.5 -> 7.5.** test_kap has 37 ok/FAIL checks and a 42-mutation manifest in which every check has an entry (checks.py); the settings are checked at init; all 44 table reads and 8 opens are checked; failures return ierr with a reason that kap_lib reports. Left:
  - the envelope Z regression (finding 1), the OPAL95 surface table's X nodes (finding 2), the conductive opacity's zero-electron case and clamped derivative (findings 3 and 4);
  - open in bug_sweep: M298, M299, M280 (documented), M279 (data), and M112's part (kap_messages.f90:33, :42 write the run log's unit directly);
  - the Alexander lookups have no range and no ierr (finding 5); kap_init still accepts three things it could refuse (finding 6);
  - test gaps: `open_table`, OPAL95's high-R shift, the precedence list, five of the seven readers, and the single-Z contract have no check (finding 7);
  - coverage 79.4% (the README's figure, dated before round 1B and 2B's eight kap commits): the lookups 72-89%, the readers 54-87%; four families run only under test_kap, one derivative check each.
- **Extensibility 6 -> 8.** A family is a directory with a lookup and a reader, a table-set type, one list entry, two settings fields and two controls; the recipe was traced on a throwaway branch and names the generators, the checker lists, the manifest and the re-pins. I traced it again by reading: it holds but for two places (finding 8). Left:
  - io/derived_controls.f90's molecular selection and kap_lib's two messages that list the interior controls are not in the recipe;
  - nine hand-edited files (the README counts eight) against TARGET_ARCHITECTURE's "one module, one registration line, one control row": the table-set type and its component live in opacity_tables, the flag and path in two records and a builder.
- **Docs 6 -> 8.** The README has all six items of the standard, and its table, flow and conventions match the code where I checked (21 rows for 21 production files; the four steps against kap_lib.f90:122-144; the error and message rules against the 27 `kap_report` sites). Left:
  - README.md:76 still gives the conductive mixture as unnormalized, which M278 changed on 2026-11-14 and which README.md:39 and :84 describe correctly;
  - README.md:66's "`kap_get` returns ierr for a point off its tables" does not hold for the molecular families (finding 5);
  - the seven lookup modules' headers still list the readers' routines as their own; a formula's sign, a misplaced comment and other texts (finding 9);
  - the coverage figure predates the conductive rewrite.

Average **7.9** (seventh: 6.2).

## Summary

- **kap now meets the standard's structure.** What the seventh report called the missing middle layers are there and match the README: a messages module, behaviour tests with a mutation manifest, a descriptor per family in one precedence list, private state, checked settings, checked reads. I sampled the rounds' claims against the code rather than the ROADMAP: the fixes of M219, M222-M226, M228-M230, M236, M277, M278, O18's ramp and KA-14's constant are in the code as described; the manifest has 42 mutations and the pinned output 37 ok lines; I/O statements are only in the seven readers, opacity_table_files and kap_messages.
- **What keeps it under 8 on two criteria is below the facade.** The lookups are still long F77-shaped routines, and the tests check each family's derivative once but not the table edges, the readers' header checks or the shared open helper.
- **New defects, all latent in the shipped decks:**
  - kap's copy of the envelope Z is the run plan's at `kap_init` (0 when no card rescales Z) and never follows the model, so an OPAL92 or LAOL89 deck with one Z table is refused at its first interior lookup unless a card rescales Z to the model's own. The F77 compared with the live ZENV. No shipped deck uses those families.
  - OPAL95's surface table puts the X = 1 - Z table at X = 1, and at Z = 0.04 or 0.10 exactly reads an X = 0.95 slot that was never built, for an envelope X above 0.8 or 0.9. Round 1 fixed the same defect in Alexander 2006 (M225).
  - The conductive opacity returns NaN with ierr 0 for a gas with no electrons (Z = 0 and no ionized H or He), and since M278 its density derivative is not the derivative of the value where a species' shifted density is clamped at the table's edge.
- **One stale README line contradicts the code and the same README** (the conductive weights), and io still holds a second copy of the molecular precedence that the recipe does not name.

## Previous report's items: status

The seventh report's kap findings (KA numbers; bug_sweep numbers from FINDING_MAP.md). Its atm items are the other reviewer's.

| Item | Status | Where |
|---|---|---|
| KA-1 (M219) Alexander 2006 cache | DONE | the cache is set after a successful build (alex06_opacity.f90:193-198) and reset with the table (alex06_reader.f90:155-156); a check and a mutation |
| KA-2 (M220, M170, L10's gap) opens and reads | DONE | opacity_table_files.f90:23-55; every reader. The conductive reader keeps its own open (finding 7a) |
| KA-4a (M222) conductive low-density branch | DONE | conductive_opacity.f90:193-224; a continuity and a slope check |
| KA-4b (M223) OPAL95 low-X, low-T shortcut | DONE | flags reset per call (opal95_opacity.f90:333-336); the postponed check itself is listed as kept (README.md:73) |
| KA-4c (M224) OPAL95 fixed Z in (0.08, 0.10] | DONE | opal95_opacity.f90:93-106 |
| KA-4d (M225) Alexander 2006 fourth X node | DONE | alex06_opacity.f90:146-150. The same defect remains in OPAL95's surface table (finding 2) |
| KA-4e (M226) helium core without an interior table | DONE | kap_lib.f90:224-228, :355-357 |
| KA-4f (M227) LAOL tolerance ln against log10 | DONE (doc; code kept) | defaults/physics.defaults:461; laol_opacity.f90:15-18; README.md:72 |
| KA-4g (M228) second LAOL T grid | DONE | laol_opacity.f90:58-66 |
| KA-4h (M229) Kurucz density bound, skipped rows | DONE | kurucz_reader.f90:89-133 |
| KA-4i (M230) Alexander 2006 header checks | DONE | alex06_reader.f90:56-82, :109-131 (19 columns, every mismatch returns ierr) |
| KA-6 (M236) settings unchecked | DONE | `check_kap_settings` (kap_lib.f90:344-389): three checks, each with a test. Three more cases in finding 6 |
| KA-8 (M237) "pure C and N" | DONE | io/derived_controls.f90:330; defaults/control.defaults:176; README.md:27 |
| KA-9 messages | DONE | kap_messages; 27 `kap_report` sites, none in a lookup; sentence case. The run log's unit is written directly (M112, open) |
| KA-10 tests | DONE, with gaps | 37 checks, 42 mutations, a coverage figure. Gaps in finding 7. The first section still boots through the star layer (test_kap.f90:19-21, :114) and keeps two pinned headers of older wording (:123, :141) |
| KA-11 (M238) recipe | DONE, two omissions | README.md:105-116; finding 8 |
| KA-12 state | DONE | `tables` private in kap_lib (kap_lib.f90:60); OPAL95's stencil a local (opal95_opacity.f90:57-62, :288) |
| KA-13 (M240) stale text | DONE for every item listed | none of the quoted texts is left (grep); new ones in finding 9 |
| KA-14 style residue | PARTLY | commented-out code in retired/; history in Sources sections; operators; Hubbard & Lampe's, kap_eval's, OPAL95's and the conductive constants named, log10(16σ/3) computed. Left: the unnamed values and names under Readability |

The seventh report's other "Left" bullets:
- `kap_eval` a 310-line if-chain: DONE. `compute_opal95_opacity` 300 lines: PARTLY (282; the three interpolators share `opal95_rows_at_logr` and `extrapolate_in_logr`). LAOL's two copies: DONE (`laol_column_opacity`).
- Precedence coded three times: DONE in kap; a fourth copy in io was not counted then (finding 8). `kap_eval` exported: DONE (private).
- Names that misstate their contents: DONE (`kurucz_table_set%opacity`, `pure_z_table_set%rho_grid`); `use_alex95_tables` for the 1994 tables is documented (alex94_opacity.f90:13-14).
- O18's ramp: DONE (opacity_blending.f90:84-107, a derivative check inside the ramp). M104: DONE. Conventions, coverage and test claims in the README: DONE.

## New findings (not in bug_sweep)

1. **kap/kap_lib.f90:186-193, :293-294; kap/opacity_tables.f90:276-279; setup/settings_install.f90:72-73, :140-147; core/run/run_setup.f90:57-67**
   - **Issue:** `kap_init` stores its `envelope_metal_fraction` argument in `tables%envelope_metal_fraction`, and `interior_opacity` refuses a zone whose Z differs from it by more than `metal_fraction_match_tolerance` (default 1e-8, defaults/physics.defaults:419) when the interior family has one table and does not interpolate in Z (OPAL92 or LAOL89 without `use_two_z_tables`). Nothing updates that copy after `kap_init`.
     - At `kap_init` the star's `envelope_metal_fraction` is the run plan's: `set_run_plan_envelope_z` assigns it only for a card with a positive rescale Z (run_setup.f90:62-65), and the member has no initializer (state/star_info_lib.f90:450), so without such a card kap gets an unset member (0 in a fresh process; in a second run of one process, whatever it holds).
     - The model's Z arrives later: `update_surface_mixture` sets the star's envelope Z from the model (core/model/starting_model.f90:1201-1202) and calls only `install_surface_x_settings`; a Z rescale sets it again (core/model/model_rescaling.f90:200) and calls `install_surface_z_settings`, which tells the EOS alone.
     - So a deck with `use_opal92_tables` or `use_laol89_tables`, no Z ramp or Z diffusion, and no card that rescales Z to the model's Z fails at its first interior lookup: "kap_get: no interior opacity at Z = 0.01880: the OPAL 1992 table is at the envelope's Z = 0.00000 (use_z_ramp or use_diffusion_z reads two)". A deck whose second card rescales Z fails the same way after it.
     - No shipped deck uses OPAL92 or LAOL89, and OPAL95 skips the test (`covers_any_z`), so nothing shipped reaches it. test_kap passes the table's Z as the envelope Z in every configuration.
   - **In the F77?** No. GETOPAC compared with the live ZENV of COMMON /COMP/ (getopac.f:101-102, `ELSE IF((Z.GT.0.12D0) .OR. ((ABS(Z-ZENV) .GT. OPTOL) .AND..NOT.L2Z .AND. .NOT.LOPAL95))THEN`), which STARIN set from the model and RSCALE on a rescale. The copy came with 50f40772 (2026-08-21: "kap_init gains an envelope_metal_fraction argument; the surface-Z sanity check in kap_get uses the kap-side copy"); the same commit gave the EOS a setter for the three places the mixture changes, and kap none.
   - **Fix:** give kap the envelope Z where it changes: `kap_update_surface_tables(X, Z, ierr)` (or a `kap_set_envelope_z`) called from `install_surface_x_settings` and `install_surface_z_settings`. Then drop the argument from `kap_init`, which today reads two star members no one has set. A test_kap check: a single-Z OPAL92 load, the envelope Z moved, a lookup at the new Z succeeds and one at the old Z returns ierr.
   - **Kind:** A (no shipped deck). **Size:** S.

2. **kap/opal95/opal95_opacity.f90:203-233 (against :129, :485-549, and kap/opal95/opal95_reader.f90:65-69)**
   - **Issue:** `build_opal95_surface_table` interpolates the fixed-Z table in X with the nodes `tbl%grid_x(...)`. The tenth X table is X = 1 - Z (the reader checks it so, opal95_reader.f90:65-69), but `grid_x(10)` is 1.0 (opacity_tables.f90:153-155). Three cases, each for an envelope X above 0.8:
     - X > 0.9 at Z below 0.04: tables 7-10 with nodes 0.8, 0.9, 0.95, 1.0 (:207-218), where the fourth table is at 1 - Z (0.98 at Z = 0.02);
     - 0.04 < Z < 0.10: tables 6, 7, 8 and 10 with the fourth node `grid_x(10)` = 1.0 (:229-232), where the table is at 1 - Z (0.94 at Z = 0.06);
     - Z = 0.04 exactly (the first branch's `<=`, :203-204) or Z = 0.10 exactly (the "tables 6-9" branch, :219-224): the stencil takes `fixed_z_opacity(9,:,:)`, which `build_opal95_fixed_z_table` builds only below Z = 0.04 (:129), so it is the initial state's zeros.
     - The 3D and 4D paths have it right: `locate_opal95_x_stencil` sets the fourth node to 1 - Z and drops the tables that do not exist (:509-543). So at such a composition the 2D path (every zone at the surface composition) and the 3D path a hair off it disagree, with ierr 0. test_kap's "the new 2D path is the old 3D path" checks this at X = 0.65 only.
     - An envelope with X above 0.8 has Y below 0.2 - Z: a pure-hydrogen or test composition, no shipped deck.
   - **In the F77?** Yes: OP95XTAB took `A(I) = XXG(IDX+I-1)` and `A(4) = XXG(10)` with XXG(10) = 1.0 (op95xtab.f:35-72; ll95tbl.f:37-38), where GETOPAL95 took `A(4) = 1.0D0 - ZZ` (getopal95.f:385); its `ZTAB.LE.0.04D0` and "USE TABLES 6-9" branches are the same (op95xtab.f:46, :59-64). Round 1 fixed this defect in Alexander 2006 (M225's entry: "at X > 0.8 and Z <= 0.06 the node was 1.0, where the table is X = 1 - Z").
   - **Fix:** build the surface table with `locate_opal95_x_stencil`'s rule (one routine for the X stencil at a Z: the tables that exist and the node 1 - Z), and extend the 2D-against-3D check to X = 0.92 and to Z = 0.06 with X = 0.85.
   - **Kind:** A for the shipped decks (their envelope X is below 0.8). **Size:** S.

3. **kap/conductive/conductive_opacity.f90:125-134, :150-155, :173-174**
   - **Issue:** the weights are divided by `electron_sum`, and each species' density is shifted by `log10(electrons_per_amu)`, with no guard. For Z = 0 and `ion_fraction` = 0 (eos_results' default, eos/eos_results.f90:38; a neutral H and He gas) both are 0: the weights are 0/0, the shift is log10(0) (a trap under the devel flags), and the routine returns NaN for κ_c with `got_conductive_opacity` true and ierr 0, which `add_conduction` turns into a NaN opacity. With Z > 0 the metals' term keeps the sum positive, so the shipped decks do not reach it; a Z = 0 model with a cool, neutral surface would.
   - **In the F77?** In another form: its unnormalized weights gave log10(0), κ_c = 0 and so κ = 0 (condopacpint.f:32-34, :57, `CONDL = -LOG10(WtH1*CONDX + WtHe4*CONDY + WtOx*CONDZ)`). The 0/0 is M278's rewrite (95ae332b, 2026-11-14).
   - **Fix:** return no conduction (`got_conductive_opacity` false) when `electron_sum` is not positive; a test_kap check at X = 1, Z = 0, `ion_fraction` = 0.
   - **Kind:** A. **Size:** S.

4. **kap/conductive/conductive_opacity.f90:144-155, :175-185**
   - **Issue:** since M278 each species' table is read at its own density, log ρ + log10((A/Z) Y_e), clamped into the table. The comment says "d log rho_i/d log rho = 1, so the derivatives keep their form", which is false where a clamp binds: K_i is then constant in ρ, but the returned dln κ_c/dln ρ still subtracts the table's slope at the clamp. For hydrogen the shift is log10(Y_e) <= 0, so its density is clamped from log ρ = -6 up to -6 - log10(Y_e): about 0.1 dex for an ionized gas, several dex where hydrogen is partly ionized (Y_e of 1e-3 gives 3 dex), which is the outer envelope of every deck. Conduction is negligible there (κ_rad/κ_c of 1e-9 to 1e-7, M299's figures), so κ's derivative moves at round-off; the two derivative checks are at a fully ionized point (log ρ 2.93) and below -6, neither in the band.
   - **In the F77?** No (M278's change).
   - **Fix:** with M299's fix, which rewrites these lines: drop a clamped species' slope from the weighted mean (or state the band in "Kept by decision" beside M299's), and add a derivative check at log ρ = -5.5 with a partly ionized mixture.
   - **Kind:** B (round-off in κ's derivatives). **Size:** S.

5. **kap/alex06/alex06_opacity.f90:199-227; kap/alex94/alex94_opacity.f90:39-57, :124-170, :314-329; kap/README.md:66**
   - **Issue:** the README's error rule says "`kap_get` returns ierr for a point off its tables". The two Alexander families have no range test in log T or at low log R: the stencil is clamped to the last four nodes and the cubic Lagrange weights extrapolate (Alexander 2006 below log T 2.7 and below log R -8; only above the top log R is the extrapolation linear, :217-223, :240-247). `compute_alex94_opacity` has no ierr or message argument at all, `lookup_alex94` always returns 0, and `build_alex94_fixed_z_table` accepts any Z (:49-57). OPAL95 extrapolates in log R without limit on both sides, which the README states for the low side only (README.md:83).
   - **In the F77?** Yes: GETALEX06 and YALO3D clamp the index and interpolate (getalex06.f:25-47; yalo3d.f:33-44), with no range test and no STOP.
   - **Fix:** either refuse a point more than one node step off the table in log T and log R (measured first: a count over the 56 decks of how far the lookups go), or say in "Kept by decision" which families extrapolate, where and how, and correct README.md:66.
   - **Kind:** A for the README; the refusal is A only if no deck goes off the tables. **Size:** S.

6. **kap/kap_lib.f90:344-389, :186-198; kap/opal92/opal92_reader.f90:34-35, :65; kap/laol89/laol_reader.f90:51-52, :78-80, :104, :113, :167; kap/kurucz90/kurucz_reader.f90:83-139**
   - **Issue:** what `kap_init` still accepts:
     - (a) the Z a table file states is read and thrown away: OPAL92's header Z into `local_grid_z` ("never used afterwards"), LAOL89's `zlot` and `zhit`. `opal_table_z1/2` and `laol_table_z1/2`, the blend's abscissae, are never compared with them, so a deck that pairs `Z0200.OPALN` with `opal_table_z1 = 0.0188` blends in Z from the wrong node, with ierr 0. The OPAL95 and both Alexander readers check every header against their grids;
     - (b) `laol_extrapolation_tolerance <= 0` (the pure-Z table of every deck takes its ln, laol_reader.f90:167): a floating-point trap under the devel flags, and no extrapolation otherwise;
     - (c) a Kurucz 1990 file with its header lines and no row: `num_temps = 0`, `kap_init` returns 0, and the run fails at the first molecular lookup;
     - (d) no interior family at all (without the pure-Z table), where every lookup above the molecular window fails. With finding 1 its message is the wrong one: the Z test at :186 comes before `interior%selected` at :194, so the deck reads "the none table is at the envelope's Z" ('none' is the unselected descriptor's name, opacity_source.f90:95).
   - **In the F77?** Yes for (a)-(c) (setllo.f:54-55 reads GRIDZ and uses it for GRIDY alone; rdlaol.f; setkrz.f). (d)'s message is 2026.
   - **Fix:** compare the file's Z with the control's in the two readers (an error naming both); three lines in `check_kap_settings`; `num_temps == 0` an error in the Kurucz reader; test `interior%selected` first in `interior_opacity`.
   - **Kind:** A. **Size:** S.

7. **kap/test/test_kap.f90, kap/test/test_kap.mutations.toml: behaviour without a check**
   - **Issue:** every check has a mutation, but these have neither (checks.py: the manifest's files are the seven lookups, opacity_blending, kap_lib, conductive_reader and opal92_reader):
     - (a) `open_table` (opacity_table_files.f90:23-37), which seven readers open through. The "missing table ... (and is not created)" check passes a missing conductive path (test_kap.f90:480-483), and its mutation is in conductive_reader's own `open` (conductive_reader.f90:36-43; the manifest's `missing-table`), the one reader that does not use the helper. M220's case, a missing OPAL92 path created empty, is untested;
     - (b) OPAL95's high-R shift and extrapolation (opal95_opacity.f90:339-376, :404-410, :680-686): M277's fix and the open M298 are in this block. The manifest's ten OPAL95 mutations are the three nodes, the three derivatives, the low-R edge, the 2D-4D seam, the surface refresh and the T range. By my reading of the stencil the pinned helium-core line (test_kap.f90:194; log T 7.8, log R -1.4) takes the shift, as a number only, and nothing takes the extrapolation;
     - (c) the precedence list (opacity_families.f90 has no mutation): no check that two selected families of a role give the first;
     - (d) five of the seven readers: the header checks of OPAL95, Alexander 2006 and Alexander 1994, LAOL's size check and Kurucz's bounds (M229, M230's fixes), which the README's coverage paragraph lists as not run;
     - (e) the single-Z contract of finding 1, `kap_update_surface_tables`' error return, and Kurucz 1990 handing a point to the interior family (`in_table` false), which only a pinned line shows.
   - **Fix:** a check and a mutation for each: a missing OPAL92 path with the conductive reader moved to `open_table`; κ continuous across a row's last filled log R and a derivative check beyond it; two interior flags on, the lookup equal to OPAL95's; one doctored header per reader written by the test, as its short OPAL92 table is (test_kap.f90:484-496); (e)'s three.
   - **Kind:** A (tests). **Size:** M.

8. **io/derived_controls.f90:119-157; kap/opacity_families.f90:4-8, :17-18; kap/README.md:16, :105-116; kap/kap_lib.f90:239-240, :356-357**
   - **Issue:** opacity_families' header and the README call the precedence list "the one place that orders the families". io's `select_eos_and_opacity_tables` orders the molecular families a second time: it warns when more than one of `use_alex06_tables`, `use_alex95_tables`, `use_kurucz90_tables` is set and switches the losers off before kap sees them. For the interior families nothing warns: OPAL92 selected with OPAL95 is "neither loaded nor used", silently.
     - The "Add a table family" recipe, traced by reading, misses this: a new molecular family added by its seven steps is outside io's count and switch-offs, so selecting it with Alexander 2006 gives no warning. It also misses the two kap_lib messages that list the interior controls by name, for a new interior family.
     - The recipe's closing count says eight hand-edited files; its steps 1-5 name nine (two modules, opacity_tables, opacity_families, controls_buffer, the registry, the builder, check_boundaries, style_counts.tsv).
   - **In the F77?** The if-chains, yes; the two homes are 2026.
   - **Fix:** `select_families` returns the selected families it passed over, and `kap_init`'s caller (or io, from that list) warns for both roles; io's copy of the order goes. Until then, add io/derived_controls.f90 and the two messages to the recipe. Build the messages' control list from the families' names.
   - **Kind:** A (the interior warning adds a run-log line only for a deck that selects two). **Size:** S.

9. **Stale or inexact text**
   - **Issue:**
     - kap/README.md:76, "the conductive mixture weights each species by its ionized number fraction without renormalizing ... (the F77 code's)": fixed by M278 (95ae332b); README.md:39 and :84 and conductive_opacity.f90:55-65 say the opposite;
     - kap/README.md:93, the coverage figure and its per-file list are of 2026-11-10; eight kap commits followed, one of which rewrote the conductive mixture;
     - kap/README.md:54 lists kap_get's core callers as "shell physics, envelope integration and starting model"; core/structure/shell_equations.f90:187 calls it too;
     - kap/opacity_blending.f90:31-32, the fit's formula has "+ (c_t + c_t2 log T) log T"; the code (:124) and the F77 (getopac.f:186) have c_t - c_t2 log T;
     - the lookup modules' headers list the readers' routines as theirs, and some say what those do wrongly since the split: opal95_opacity.f90:7-8 ("read_opal95_tables: read a table file, then build the fixed-Z table"; `build_opal95_fixed_z_table` is not listed), opal92_opacity.f90:7 ("read ... and spline them"), alex94_opacity.f90:7-8, kurucz_opacity.f90:7, alex06_opacity.f90:7, laol_opacity.f90:8-14, conductive_opacity.f90:6;
     - kap/laol89/laol_opacity.f90:350-355, `cap_laol_opacity`'s comment sits above `laol_column_message`'s, and the routine itself (:375) has none;
     - kap/opacity_tables.f90:107-108, "table_number (2 for the second Z's, whose run-log lines say "#2" and "SECOND")": no such lines exist; it selects the name 'second LAOL 1989' in a message (laol_opacity.f90:150-151);
     - kap/opacity_tables.f90:62-63, the markers' tests quoted as ".ge./.gt."; :136-141, :145-146, :181-182, :209-210, :243-244 name the former COMMON blocks in the types' comments, which belong in a Sources section;
     - kap/conductive/conductive_opacity.f90:30-36, two comments run together without a break (the constant's, then the table's reach);
     - kap/test/test_kap.f90:58 cites core/read_input.f90, which does not exist; :62-64 name the F77 flags (LOPAL95, ZOPAL951, TMOLMIN); :135-140 is history.
   - **Fix:** text only; the coverage re-measured with the command the README gives.
   - **Kind:** A. **Size:** S.

10. **Code-style residue against GUIDELINES "Declarations and interfaces" 5 and "Code style" 1, 2 and 4**
    - **Issue:**
      - `pure`: 15 of 113 procedures. Lookups whose table set is `intent(in)` and that write only into `message` could be declared so and have the compiler hold "kernels print nothing": `compute_pure_z_opacity`, `laol_column_opacity`, `laol_column_message`, `select_four_point_window`, `compute_conductive_opacity`, `condopacp`, `cinterp3`, `alex_x_stencil_start` (already), and the builds (`build_opal95_surface_table`, `build_alex94_*`, `build_laol_splines`, `build_pure_z_splines`), as far as numerics' splines are pure;
      - `kap_get`'s `ion_fraction` is `intent(inout)` where the header says kap does not change it (kap_lib.f90:80), and is copied into a 3-element component with no size check (:85-86);
      - three names imported and not used (unused.py): `n_laol_x` (laol_opacity.f90:25), `n_opal95_tables` and `opal95_missing_opacity` (opal95_opacity.f90:24-26);
      - the conductive opacity's atomic weights are local literals, 1.008, 4.004 and 16 (conductive_opacity.f90:88-89; He4 is 4.0026), the F77's (condopacpint.f:14), with no source;
      - `log10(16 σ/3)` is computed at every call (conductive_opacity.f90:108) because a parameter cannot call math_lib's log10; a `protected` module value set once at `kap_init` would do;
      - the items listed under Readability: the three long routines, opal95_reader's labels, alex06_reader's repeated body, the three spellings of 10**x, the unnamed values, the F77 names, the host-associated work variables, the interim "(formerly ...)" names;
      - test_kap.f90 uses `double precision` throughout and `use kap_lib` and `use opacity_tables` without `only:` (:22-23); check_style does not cover test directories.
    - **Fix:** byte-identical commits: `pure` where the compiler accepts it; `intent(in)`; the imports; one `ten_to(x)` choice per family is a result change at round-off unless each family keeps its own, so name the three in one place and leave them; split `compute_opal95_opacity` into locate (T, R), check the high-R corner, weights in R, and the 2D/3D/4D choice; opal95_reader on `read_failed`; one `read_alex06_table(slot, x, z)`.
    - **Kind:** A (computing the atomic weights from a constants module would be B). **Size:** M.

## What would raise each score

- **Readability, to 8-8.5** (a short round's list, all byte-identical):
  - `compute_opal95_opacity` as four named steps, with the high-R corner test its own procedure (finding 10);
  - `compute_kurucz_opacity`: one header, the row search a procedure; `compute_conductive_opacity`: one `species_conductivities(logrho)` for the two branches;
  - opal95_reader on `read_failed` without labels; alex06_reader's table body once;
  - today's names for the 23 lines of F77 locals in opal92, alex06, alex94 and kurucz; `opal92_surface_slot`, and the other unnamed values with unit and meaning;
  - the internal procedures given their own locals; the interim "(formerly ...)" names dropped, the F77 ones moved to the Sources sections;
  - `pure` on the lookups and builds that qualify.
- **Organization, to 9:**
  - the envelope Z's entry (finding 1) and `kap_init` without the two unset star members;
  - the precedence in kap alone, with the warning for both roles (finding 8);
  - conductive_reader on `open_table`; LAOL's unused header mixture read into a local, without composition_layout;
  - `kap_settings` private to kap_lib.
- **Maintainability, to 8-8.5** (a short round's list):
  - findings 1, 2, 3 and 6 fixed, each with its check and mutation (all byte-identical for the shipped decks);
  - finding 7's checks: `open_table`, OPAL95's high-R edge, the precedence, one doctored header per reader, the single-Z contract;
  - M298 (byte-identical with the shipped tables) and, as one measured commit, M299 with finding 4;
  - finding 5 decided: a refusal after a count, or the README's statement;
  - M112's kap part (the run-log write through `log_line`);
  - the coverage re-measured, with the lookups above 85%.
- **Extensibility, to 8.5-9:** finding 8's two recipe lines now; then a family that owns its table-set type (an allocatable or class member of the state in place of the fixed aggregate), and settings the family reads from its own record, so that the hand edits fall from nine files toward TARGET_ARCHITECTURE's three.
- **Docs, to 9:** finding 9; README.md:66 and "Kept by decision" after finding 5; the coverage figure re-dated.

## Traces, counts and README checks (supporting)

- **The rounds' claims sampled against the code:**
  - one messages module: I/O statements are in the seven readers (13, 5, 5, 7, 22, 5, 16), opacity_table_files (2) and kap_messages (1), and nowhere else (counts.py); 27 `kap_report` call sites, none in a `*_opacity` file;
  - reads and opens: none of the readers' 44 reads and 8 opens is without `iostat`, `end=` or `err=`;
  - tests: 37 lines of expected_test_kap.out end in "ok"; the manifest has 42 `[[mutation]]` entries, at least one per check (checks.py). ROADMAP's 34/36 (round 1), 40 (1B) and the README's 42 (2B) are consistent with the commits that added checks;
  - the measured fixes as described: O18 (opacity_blending.f90:91-103), M222 (conductive_opacity.f90:214-221: the mixture's sign, the step toward the point), KA-14 (:108), M277 (opal95_opacity.f90:359-371 search `fixed_z_opacity`), M237 (io/derived_controls.f90:330), M278 (:125-134, :150-155);
  - check_boundaries: `PUBLIC['kap']` is kap_get, kap_init, kap_update_surface_tables; the seven lookups and opacity_blending are in STAR_BLIND_FILES; kap_lib, opacity_families, opacity_source and opacity_tables in KERNEL_FILES; kap's style_counts rows are 0 but `io_statement` in the reader directories.
- **README checked:** §1's 21 rows against the 21 production files and the procedures each names; §2's init order against kap_lib.f90:290-301 and the four steps against :122-144; §3's error and message rules; §4's count and manifest; §6's recipe (finding 8). Not checked: §5's coverage, the mutation run and the devel run (nothing was run).
- **Add a table family, traced by reading:** the family's two modules; opacity_tables (type, component, flag, path); opacity_families (use line, list entry, the list's length); io/controls_buffer and defaults/controls_registry.tsv; setup/kap_settings_builder (two copies); check_boundaries (STAR_BLIND_FILES); style_counts.tsv; test_kap, its manifest and pinned output; the README row; then the generators and re-pins of step 7. Not in the recipe: io/derived_controls.f90:147-157 (a molecular family) and kap_lib.f90:239-240, :356-357 (an interior family).
- **Open bug_sweep items of kap, weighed and not re-reported:** M298, M299, M280 (documented), M279 (data), M112 (kap_messages' direct write).
- **Longest procedures** (counts.py): compute_opal95_opacity 282, read_opal95_tables 174, condopacp 164 (Potekhin's, kept under its name), compute_kurucz_opacity 158, compute_conductive_opacity 157, read_alex06_tables 127, read_kurucz_tables 123, compute_alex94_opacity 122.

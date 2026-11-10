# eos -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

Scope: every file of src/eos/ read whole: README.md, eos_lib, eos_configuration, eos_results, eos_messages, eos_mixture, blended_eos, the four table modules, yale/yale_eos, scv/scv_eos, opal/opal_eos, opal/opal95_eos, mhd/mhd_eos and test/test_eos.f90. Also read whole: the repository root's test_eos.py, setup/eos_settings_builder.f90 and its install site (setup/settings_install.f90:63-75, :141-144), and io/derived_controls.f90's `select_eos_and_opacity_tables`. Skimmed: expected_test_eos.out (its 49 ok lines), the eos parts of tools/check_boundaries.py, tools/controls_registry.py's `builder_copies`, defaults/controls_registry.tsv's EOS rows, io/README §5's control recipes, .github/workflows/ci.yml, and docs/bug_sweep.md's eos entries (M85-M88, M92-M94, M105-M107, M119, M122, M123, M128, M129, M136, M137, M144, M148, M156, M157, O1, O2, O18). The F77 was read with `git show 6cd56739:src/<file>.f`. Since the fourth evaluation (e8c847ba), two commits touched src/eos/, 9 lines in all. 509cd96d makes eos_messages and test_eos take `run_log_unit` from state's run_messages instead of io_units, so eos imports no io module (`check_physics_no_io`). db9c44e6 renames the builder file in four comments. The code is otherwise the code the fourth evaluation scored. Nothing was built or run. The counts come from scratch Python over the tree (SP/review7_eos/counts.py and spline_weight.py, with their .out files). Paths are under src/ unless they start at the repository root.

## Scores (1-10)

First -> second -> third -> fourth -> now (the first three were scored out of 5 and are doubled here):

- **Readability 6 -> 7 -> 8 -> 7.5 -> 7.5.** The front half (eos_lib, blend_eos's four steps, eos_mixture, the settings and results records) still reads at rotation's level, and the kernels are unchanged. Left:
  - 15 procedures over 100 lines (2,763 lines in all): `esac` 381, `interpolate_opal_eos` 350, `compute_scv_envelope_eos` 294, `build_scv_envelope_table` 191, `compute_scv_eos` 171, `rhoofp` 165, the two OPAL readers 160 and 159, `opal_density_search` 153, `mhdst1` 144, `iterate_fully_ionized_density` 142, ...;
  - the Fermi-Dirac solver's unnamed rows (`ff(1..5)`) and reused `pr1..pr9`, and the Saha sensitivities `sqet`, `sx1qt`, ...;
  - about 35 unnamed thresholds in the kernels, among them one OPAL sentinel family written 18 times (finding 7);
  - 10 of 104 procedures are `pure`, and at least 20 more have no side effects; explicit-shape dummies are declared on 19 lines, mostly in MHD's readers (finding 8);
  - history in procedure bodies, commented-out code in `radsub`, and about 19 capitalized messages kept as the F77 had them.
- **Organization 8 -> 8.** One facade with three entries, one settings record, a module per family, star-free and star-blind; eos now imports no io module. Left:
  - no table-source descriptor (TARGET_ARCHITECTURE 4.5(c)). SCV is nested in `gas_eos` with no edge weight (M88), and MHD is a branch of `eos_get`;
  - OPAL 1995 is a second copy of the 2001/2006 code (`esac`, `rhoofp`, `readco`, `radsub`, `gmass`, `t6rinterp`: about 1,000 lines; M94);
  - test-only entry points public in production modules (`gamma1_at_density`, `blend_eos`, `set_eos_settings`, `esac01`/`esac06`, `esac`, `interpolate_mhd_table`). Writable public state: `eos_mix`, `opal_shared`, which test_eos.f90:824 resets by hand, and the SCV and OPAL search state (M107).
- **Maintainability 7 -> 7.** There are 49 ok/FAIL checks and a byte pin, but no eos defect has been fixed since 2026-10-30. Left:
  - open latent defects: M92-M94, M105-M107, M128, M136 and M137, and with them M129's two 1/mu_ion rules;
  - weak tests: M87, M119 and M144;
  - four new items (findings 1-4): a latent SCV read of unbuilt table entries, an error message that lost its text, an unchecked table range and a misdirected diagnostic;
  - the kernels format and print about 80 messages, and the Saha and OPAL non-convergences return success (GUIDELINES code style 7);
  - CI builds release only, so `-fcheck=all` never runs test_eos;
  - the "each confirmed by a mutation" claim cannot be re-checked from the tree: math/ and numerics/ ship a mutate_*.py, eos does not.
- **Extensibility 6.5 -> 6.5.** The recipes exist and the control recipe holds when traced: it hands over to io/README §5 for the `settings` cell, which `check_boundaries` verifies. Left:
  - no source descriptor, so a source is still edits in eight places;
  - the source recipe's mean-weight trap and dropped ion fractions (M148, M123);
  - a new trap in the OPAL-vintage recipe: its step 5 sends a newcomer to a builder copy of a control that is not in the star (finding 5).
- **Docs 7.5 -> 7.** The README is thorough and right on the main paths, and every claim in §1-§3 that I checked holds. It is lowered half a point because none of the fourth evaluation's documentation items was fixed (M119, M122, M123, M148, M156 and M157 are open), and because I found two more inaccuracies a reader would act on:
  - the vintage recipe (finding 5);
  - test_eos.py's docstring, which names a routine that no longer exists (finding 6).

  Also left: an undocumented validity range (finding 3), and §4's 1,256-character paragraph (README.md:73; rotation's README has none over 1,000).

Average **7.2** (7.3 at the fourth evaluation).

## Summary

- **Nothing has changed in eos since the fourth evaluation.** Two cross-cutting commits touched nine lines: the messages' run-log unit now comes from state/, and the builder file was renamed. eos is now the lowest-scored package. io, state and setup were raised past 8.9 in three rounds that did not include eos. The front half (facade, blend, mixture, settings, results) meets the rotation standard, but the family kernels (Fermi-Dirac, SCV, OPAL, MHD) keep their F77 shape. So do the 19 open eos defects in bug_sweep.
- **None of the fourth report's items is done** (status below). Its ten findings became M128, M129, M136, M137, M144, M148, M156 and M157, all open. The style residue (its finding 10) was left out of bug_sweep by that section's rule and is still in the code.
- **New, all latent or documentation:**
  - **SCV smoothing reads unbuilt table entries.** `compute_scv_envelope_eos` smooths toward the cell below at the columns of the cell's own first row, but that row below can be shorter (in 24 of the 62 row pairs of the shipped table). In practice zeros then enter the shifted cell's spline. This is inherited from the F77, the effect is at most about 0.5%, and it happens only at pressures above the lower row's highest, which no shipped deck is known to reach.
  - **The OPAL 1995 pressure check prints two bare numbers.** When the stop became an `ierr`, the F77's text "ERROR IN OEQOS PTOT" was dropped, so the message no longer says what failed.
  - **The Fermi-Dirac table's top and its degeneracy range are unchecked.** Above log T 9.0, and outside X = log rho - log mu_e - 1.5 log T in [-16, -3], the fully ionized gas extrapolates silently. M137 covers only the cutoff side.
  - **A misdirected OPAL diagnostic.** The second X quadratic's overflow check tests the first quadratic's value.
  - **Documentation and style:**
    - the vintage recipe's step 5;
    - test_eos.py's stale docstring;
    - the OPAL sentinel family as literals;
    - few `pure` procedures, explicit-shape dummies, and the test's literal composition slots.
- To reach 9, eos needs three things, all named in the fourth report: the source descriptor, the open latent defects fixed (most fixes are byte-identical), and the OPAL, SCV and Fermi-Dirac kernels split into named steps with OPAL 1995 on the shared layout.

## Previous report's items: status

**The fourth report's new findings:**
1. The source recipe's mean weights on a fresh record (now M148): **OPEN** (README.md:122, :128 unchanged).
2. Two 1/mu_ion rules under diffusion of Z (M129): **OPEN**. It is in bug_sweep's "Decisions needed": the Yale/SCV scaling for OPAL and MHD too, or the zone's composition everywhere. M181 (io) shows it reaches the 169 decks that `diffuse_lithium` puts on diffusion of Z.
3. SCV's H+ fraction (M128): **OPEN** (scv_eos.f90:649).
4. The SCV reader ignores the He/Z grids; OPAL 1995 accepts a zero T6 (M136): **OPEN** (scv_eos.f90:909-911; opal95_eos.f90:1007-1008).
5. The Saha cutoff is never range-checked (M137): **OPEN** (eos_lib.f90:102-118 checks the vintage and the Debye-Hueckel range only).
6. Deck counts "298 of the 304" (M157): **OPEN** (README.md:45, :111; there are still 310 decks in testsuite/ and examples/).
7. Stale comments (M157): **OPEN**:
   - yale_eos.f90:184, "Called from blend_eos.f90";
   - yale_eos_tables.f90:33, "eos_lib's debye_huckel_composition_of";
   - opal_eos_tables.f90:345, "48 MB apiece";
   - mhd_eos.f90:46, :127 and :160, "a message in the run log";
   - mhd_eos.f90:510, the PGL header.
8. Which source fills which field (M156): **OPEN** (README.md:43, :63; eos_results.f90:18-19, :39).
9. The defaults check against literals (M144): **OPEN** (test_eos.f90:1029-1033).
10. Style residue: **OPEN**. It is not in bug_sweep (section 10's rule).
    - opal95_eos.f90:826-836: the commented-out `cvt=`, `EN=EN*REVISE` and the open question "should this be ET=ET*REVISE???".
    - In-body history: scv_eos.f90:127-130 and :152, opal_eos.f90:245-247 and :758-761, opal95_eos.f90:630-633 and :739-743, test_eos.f90:274-284.
    - Only the four table modules have a Sources section.

**What would raise readability:**
- Name the Fermi-Dirac rows, split `fermi_solution`, replace `pr1..pr9`: **OPEN** (yale_eos.f90:126-170, :770-800, :937-946, :1014-1033).
- Readable names for the Saha sensitivities: **OPEN** (yale_eos.f90:100-120).
- `interpolate_opal_eos` in named steps; SCV's four-way smoothing as one helper: **OPEN** (opal_eos.f90:433-782; scv_eos.f90:110-243).
- `mhdst1` with a record and assumed-shape arrays: **OPEN** (mhd_eos.f90:275-418; I count 36 arguments, the fourth report 33).
- History and commented-out code out of the bodies: **OPEN** (above).

**What would raise organization:**
- SCV as a step with its own edge weight (M88), MHD as a source: **OPEN**.
- The table-source descriptor: **OPEN**.
- Test-only entry points out of production modules: **OPEN** (blended_eos.f90:36, eos_configuration.f90:54, opal_eos.f90:41-43, opal95_eos.f90:31, mhd_eos.f90:33).
- `eos_mix` protected: **OPEN** (eos_mixture.f90:49; only eos_lib.f90:153-156 writes it, but nothing enforces that).
- Search state in a record (M107); OPAL 1995 on the shared layout (M94): **OPEN**.

**What would raise maintainability:**
- Fix M92, M93, M94, M105, M106 and M107: **OPEN**, all six.
- Kernels' failures as codes, reported by `opal_step`/`eos_get`: **OPEN** (yale_eos.f90:342-348, :892-902; opal_eos.f90:372-379).
- One 1/mu_ion rule: **OPEN** (M129, decision pending).
- test_eos under `make devel` in CI: **OPEN** (ci.yml builds `make -C src -j4 tests` only).
- OPAL coverage (71%), M119, M87: **OPEN**.

**What would raise extensibility:**
- Correct the source recipe (M148, M123): **OPEN**.
- The descriptor, making a source one module plus one list entry: **OPEN**.

**What would raise docs:**
- Findings 1 and 6-9, M119, M122, M123: **OPEN**.
- The per-path field table in §3: **OPEN**.
- Coverage re-measured in a pin-change commit: **OPEN** (still "85.5% ... (2026-10-29)", README.md:113; no pin has moved since).

**Done elsewhere since:** eos no longer imports io. eos_messages and test_eos take `run_log_unit` from state/run_messages.f90 (509cd96d), and `check_physics_no_io` enforces it. The builder is `pure` in its own file, `setup/eos_settings_builder.f90` (db9c44e6), and the README's five references and eos_configuration's two follow it.

## New findings (not in bug_sweep)

1. **scv/scv_eos.f90:111-122, :159-163, :198-200, :297-310** (`compute_scv_envelope_eos` and its `spline_cell`), against `find_scv_cell` (:583-603) and `build_scv_envelope_table` (:784)
   - **Issue:** the downward temperature smoothing (`temp_smooth_direction = -1`, guarded only by `idtt > 1`) interpolates the cell shifted one row down, rows idtt-1..idtt+2, at columns idp..idp+3 (idp..idp+4 with the pressure shift as well). `find_scv_cell` bounds idp by the cell's own first row (idp <= nptsx(idtt) - 3). But the rows' point counts grow with T: the build fills row t only to nptsx(t) (:784), and the reader reads no more (:915).
     - **Where:** in 24 of the 62 row pairs of the shipped h_tab_i.dat, row idtt-1 is shorter, and then columns past nptsx(idtt-1) are never-built entries of `tablenv`. These are zero in practice: static storage, undefined by the standard. A row of log rho = 0, zero derivatives and zero du/dT enters the shifted cell's temperature spline.
     - **The input that triggers it:** a point in the band log T in [tlogx(idtt+1), tlogx(idtt+1) + 0.032) at a gas pressure above the lower row's highest columns. Examples: log T 3.700-3.732 with log P_gas >= 11.2, where row 3.54 has 38 points to log P 11.4 and row 3.62 has 49; and log T 4.18-4.212 with log P_gas >= 13.8.
     - **Size:** the natural spline's weight of the first node at x = 2 to 2.4 is 0 to 0.026, and the shifted cell's weight is 0.5(1 - d/tol). So the error is at most about 0.005 times each quantity (log rho in dex). It vanishes at both ends of the band, so it is a bump, not a step (SP/review7_eos/spline_weight.out). These are cool, dense points that OPAL 2006 does not cover; no shipped deck's star is known to reach them (not measured).
   - **In the F77?** Yes. EQSCVE sets ISMT = -1 under `IF(IDT.GT.1)` and reads `TABLENV(II,IDP+K-1,J+1)` with `II = IDT+I-1+ISMT`, unchecked against NPTSX(II) (`git show 6cd56739:src/eqscve.f`, :69-71, :176-182). SETSCV builds `DO J=1,NPTSX(I)` (setscv.f:101).
   - **Fix:** smooth downward only when the lower row holds the columns (`idp + 3 + max(press_smooth_direction, 0) <= nptsx(idtt-1)`, decided after the pressure smoothing), or else keep the cell's own result. `test_eos` gets a check in the band at log T 3.71, log P_gas 12. The same edge applies if the smoothing becomes one helper (the fourth report's readability item).
   - **Kind:** B (only where reached). **Size:** S.
2. **opal/opal_eos.f90:132-139** (`compute_opal_eos`'s pressure check, OPAL 1995's `checks_pressure`)
   - **Issue:** when the interpolated pressure differs from P by more than 5e-7, the error goes to the terminal and the run log as two bare numbers (`write(eos_lines,*) p_e12, output(i_p)`), and `ierr = 1` stops the run. Nothing says what failed or in which routine. eos_get's header says "the message is already out".
   - **In the F77?** No: the text is a 2026 loss. OEQOS wrote the same two numbers to the run log, then `STOP ' ERROR IN OEQOS PTOT'` (oeqos.f:68-71). 89b648ce ("Stage 3: eos domain converted to the ierr-not-stop discipline") dropped the stop text along with the stop.
   - **Fix:** one line that names the check, for example "OPAL 1995 EOS: the table's pressure at the density found differs from P by more than 5e-7 (P12, table):". Latent: no shipped deck uses OPAL 1995.
   - **Kind:** A. **Size:** S.
3. **yale/yale_eos.f90:733-737 and :805-808** (`start_fully_ionized_solution`'s T row, `iterate_fully_ionized_density`'s X bin)
   - **Issue:** the Fermi-Dirac table spans log T 5.2-9.0 (rows 5 + 0.2 j, j = 1 to 20) and X = log rho - log mu_e - 1.5 log T from -16 to -3 (the 43 abscissae, binned into 261 steps of 0.05). Both indices are clamped (`jt1` to 1..18, `kk` to 1..261), so a point beyond either edge is extrapolated quadratically, without a message or status. M137 covers only the low-T side, through the cutoff.
     - Above log T 9.0 the extrapolation is in T; past X = -3 (for example log rho >~ 8 at log T 7) it is in degeneracy.
     - README §2.1 gives "the fully ionized gas at high density (log rho > 7) and in the stellar centre above the OPAL table" and no range.
     - No shipped deck is likely to reach either edge: log T 9 is carbon ignition, and the shipped decks stop at the TAHB.
   - **In the F77?** Yes (eqrelv.f:56-57, :81-83).
   - **Fix:** name the table's range (with M137's first-row parameter). Then either an `eos_log` warning the first time a run leaves it, or `ierr`, which would make it a Kind B change. README §2.1 and §3 state the range.
   - **Kind:** A (documentation and warning). **Size:** S.
4. **opal/opal_eos.f90:729-731 and opal/opal95_eos.f90:592-601** (the second X quadratic's overflow message)
   - **Issue:** when X lies between x(2) and x(3) the interpolation blends two quadratics. The second one's check `if (v%x_interp_result(...) > 1.0d+20)` reads the first quadratic's value, already checked at opal_eos.f90:708 (opal95_eos.f90:563), not `x_interp_result_alt`. An empty entry (the 1e35 fill) in the fourth X table is never named by this diagnostic. The message prints the second stencil's three table values for a check made on the first. The bad value is still caught later by `t6rinteos`'s 1e15 test (an `ierr` with "Interpolation indices out of range"), so only the diagnostic is misdirected.
   - **In the F77?** Yes (esac06.f:447, `if(esk(it,ir) .gt. 1.D+20)` after computing `esk2`).
   - **Fix:** test `x_interp_result_alt`; better, once the shared code is one copy (M94), in one place.
   - **Kind:** A. **Size:** S.
5. **eos/README.md:139** (the "Add an OPAL vintage" recipe, step 5)
   - **Issue:** "The control and the path, as for a source (steps 5 and 6 above)". The source recipe's step 5 copies the control in `eos_settings_from_star`, and its step 6 hands over to io/README's "Add a control" and "Add a table path". Those recipes give the row home `control`/`physics`, `used` `run`, and `settings` the paths record's field (io/README.md:258, :274). None of this fits a vintage:
     - its switch and path are buffer-only controls that only `select_eos_and_opacity_tables` reads (io/derived_controls.f90:125-157, the recipe's own step 4). Their registry rows are home `buffer`, used `read`, `settings` blank, as `use_opal2006_eos` and `opal2006_eos_table_path` are;
     - the builder copies the derived `opal_eos_vintage` and `opal_eos_table_path` already and needs no edit.

     Followed literally, step 5's `s%physics%use_opal20xx_eos` does not compile (a buffer control is in neither record), and a row made to match the io recipe fails `check_boundaries`' `home`/`used` checks.
   - **Fix:** step 5 becomes: the switch and the path as buffer controls (io/README §5 step 2's "home `buffer` if only the read uses it", `used` `read`, `settings` blank); no builder edit. Re-trace the recipe.
   - **Kind:** A. **Size:** S.
6. **test_eos.py:7-10** (the repository root's test file docstring)
   - **Issue:** it says the test "exercises the eos facade entries (eos_get, eos_get_gamma1)" and names "eos_get_gamma1's Yale/SCV branch". `eos_get_gamma1` has not existed since 64500ba9 (2026-10-29): its successor is `blended_eos`'s test-only `gamma1_at_density`, and the facade's entries are `eos_init`, `eos_get` and `eos_set_mixture`. "Stage-0" is old jargon. The docstring also leaves out what makes the test fail besides the byte pin: the 49 ok/FAIL checks, through the program's exit status, which `assert result.returncode == 0` tests.
   - **Fix:** the facade's three entries, `gamma1_at_density` by name, and one sentence on the ok/FAIL checks.
   - **Kind:** A. **Size:** S.
7. **One value written twice, and an unnamed sentinel** (opal/opal95_eos.f90:70; opal/opal_eos.f90:102; the 1e35 family)
   - **Issue:**
     - `find_opal95_ramp` tests `log10_density < -14d0 .or. > 5.0d0` as literals. The 1995 row of `opal_layout` holds the same two values (`ramp_log10_rho_min`, `window_log10_rho_max`; opal_eos_tables.f90:227-229), which `find_table_edge_ramp` reads for 2001 (opal_eos.f90:218).
     - `compute_opal_eos` tests `density_cgs <= -998.0d0` against the named `opal_rho_not_found` (-999; opal_eos_tables.f90:95, whose comment repeats the -998).
     - The readers fill empty table entries with `1.0d+35` (opal_eos.f90:858, opal95_eos.f90:975). Four thresholds keyed to that fill, 1e30, 1e25, 1e20 and 1e15, plus a 2.0d+30 stand-in, are literals at 16 places across the two copies, 18 with the two fills (a grep of src/eos/opal/ and opal_eos_tables.f90). Changing the fill means finding them all.
   - **Fix:** `opal_empty_entry` and its thresholds as named parameters in opal_eos_tables, with their meaning; `find_opal95_ramp` reads `lay95`; `== opal_rho_not_found`.
   - **Kind:** A. **Size:** S.
8. **Declarations and the test's slots** (GUIDELINES "Declarations and interfaces" rules 2, 3 and 5, and "Code style" rule 4)
   - **Issue:** 10 of the 104 procedures in eos are `pure` (SP/review7_eos/counts.out). At least 20 others have no side effects and are not marked:
     - `blend_saha_into`, `adopt_table_result`, `blend_table_result`, `inside_opal_window` and `zone_mixture`;
     - `quad`, `gmass`, `gmass01` and `interpolate_opal_t6_rho`;
     - Yale's eleven derivative and interpolation steps, `set_saha_constants` through `compute_internal_energy_and_cp`;
     - `scv_ionization_fractions`.

     Explicit-shape dummies are declared on 19 lines (14 of them in mhd_eos.f90:129-753), not counting continuation lines. `mhdst1`'s eleven work-array dummies (:293-300) are declared without an intent. test_eos.f90 adds three more departures:
     - it has `use eos_lib` without an only-list (:32);
     - its main program declares `double precision` and `data` statements (:60-107);
     - it fills compositions by literal slot, `comp(1)`, `comp(5) = 0.17d0*z`, `comp(7)`, `comp(9)` (:412-417, :1435-1440), where composition_layout's `i_h1`, `i_c12`, `i_n14`, `i_o16` are meant.
   - **Fix:** mark the side-effect-free procedures `pure` (byte-identical); assumed shape in MHD with the `mhdst1` record; the test with `only:`, `real(dp)` and the named slots.
   - **Kind:** A. **Size:** M (mechanical).

## What would raise each score

**Readability 7.5 -> 9**
- Cut the F77-shaped kernels into named steps:
  - `interpolate_opal_eos`/`esac`: the X stencil, coverage and orders, the X interpolation, T6-rho, scaling and radiation;
  - `compute_scv_envelope_eos`'s smoothing as one helper applied in T and then P, with finding 1's guard;
  - `build_scv_envelope_table`: the mixing entropy, then the combination;
  - the Fermi-Dirac solver, with named table rows (`i_fd_*`) and named quantities for `pr1..pr9`/`cden*`.
- Name every threshold: finding 7's sentinel family; the density iteration's 1e-15, 4 dex and 1e-8 (yale_eos.f90:882-892); MHD's 5e-7 and 1e-4.
- Write SCV's radiation correction (dlnP/dlnT, du/dT, cp, grad_ad) once: it appears in both lookups (scv_eos.f90:251-273, :470-531) and in the build.
- Finding 8: `pure`, assumed shape, the `mhdst1` record.
- A Sources section in each of the six solver files, with the in-body history moved there. Drop the commented-out code or move it to eos/retired/ with a header. In a re-pin round, rewrite the capitalized messages and the misspelled "WARNNING... OPAL TBL FAIL".

**Organization 8 -> 9**
- The table-source descriptor (name, window, reader, lookup, edge weight), with `blend_eos` looping over an ordered list. SCV becomes a source with an edge weight, which also fixes M88.
- OPAL 1995 on the shared layout with clamped stencils (M94), which deletes about 1,000 duplicated lines.
- The test-only entry points behind a test module; `blend_eos` private; `eos_mix` protected and set only by `eos_set_mixture`; the OPAL loaded flags reset through an eos routine rather than poked by the test.
- The SCV and OPAL search state in a lookup record (M107; `find_ragged_edge_ramp` reads the last interpolation's stencil).

**Maintainability 7 -> 9**
- Fix the latent defects. M93, M105, M106, M107, M136, M137 and findings 1-4 are small, and most fixes are byte-identical. Fix M92, M94 and M128 with their tests.
- Return the kernels' failures as codes and report them once in `eos_get`: Saha non-convergence, M93, and the OPAL search's fallback flagged in the result.
- Run test_eos under `make devel` (`-fcheck=all -ffpe-trap`) in CI. That would catch M92, M94 and M107 mechanically.
- Commit a `mutate_eos.py` beside math's and numerics', so that "each confirmed by a mutation" can be re-run. Raise OPAL coverage (71%) with synthetic corrupt tables. Fix M87, M119 and M144.

**Extensibility 6.5 -> 9**
- The descriptor: a new source is then one module, one list entry and its settings.
- Correct the source recipe (M148, M123) and the vintage recipe (finding 5), and re-trace both from a fresh record on a throwaway branch.

**Docs 7 -> 9**
- Fix M148, M156, M157, M122, M123 and M119's test label, and findings 3, 5 and 6.
- Add the per-path table to §3: which source fills each `eos_result` field, and the Z and mixture rule each source uses.
- Split §4's 1,256-character paragraph.
- Date and re-measure the coverage figure with the next pin change.

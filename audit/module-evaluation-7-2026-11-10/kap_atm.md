# kap_atm -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

Scope: every file of src/kap/ and src/atm/ read whole:
- kap: README.md, kap_lib, opacity_tables, the seven family modules (opal95/, opal92/, laol89/, alex06/, alex94/, kurucz90/, conductive/), test/test_kap.f90 and expected_test_kap.out;
- atm: README.md, atm_lib, atm_tables, t_tau_relations, tables/kurucz_atmosphere, tables/allard_atmosphere, test/test_atm.f90 and expected_test_atm.out;
- also whole: the repository root's test_kap.py and test_atm.py, setup/kap_settings_builder.f90, setup/atm_settings_builder.f90 and setup/settings_install.f90.

Read for the claims they touch:
- the kap/atm entries of tools/check_boundaries.py (PUBLIC, STAR_BLIND_FILES, STAR_FREE_DIRS, FACADE_ONLY_DOMAINS);
- the callers in core/ (shell_physics, shell_equations, envelope_integration_kernel, envelope_integrands, surface_boundary, mass_change, starting_model, rezoning) and convection/semiconvection;
- io/derived_controls.f90's run-log banner, and the control descriptions in defaults/;
- GUIDELINES "Naming", "The star object", "Declarations and interfaces" and "Code style", and rotation/README.md as the bar;
- TARGET_ARCHITECTURE 1.1, 4.5(c) and 4.8;
- docs/bug_sweep.md (every kap/atm entry), and docs/ROADMAP.md sections 10-11 and MODERNIZATION_CHANGES section 5, where the pre-divergence sweeps' open items are kept;
- audit/bugsweep-2026-08-31/kap.md and atm_net_math.md, audit/bugsweep-2026-09-01-pass2/K_opacity.md and L_atm_controls.md, and audit/readability-sweep-2026-09-03/kapatm.md.

How the evidence was gathered:
- every inherited defect was checked against `git show 6cd56739:src/<file>.f`;
- nothing was built or run;
- the counts come from scratch Python over the tree (SP/review7_kap_atm/): counts.py (style residue, message sites, opens and reads), decks2.py (which families the 1,164 shipped decks use), ck_sentinels.py (the -999 entries the shipped atmosphere tables carry) and alex06_headers.py (read_alex06_tables' header checks replayed on the shipped table).

Paths are under src/ unless they start at the repository root. This is the first module evaluation of either directory. Before it there were only the bug sweeps of 2026-08-31 and 09-01 and the readability sweep of 09-03, with no scores.

## Scores (1-10)

No earlier scores; first score for each directory.

### kap/

- **Readability 6.5.** The interface is modern:
  - kap_lib's header, the settings and paths records, and the table-set types (`kurucz_table_set`, `opal92_table_set`, `laol_table_set`);
  - `blend_in_z` and the shared `stencil4_locate`;
  - verb-named family routines, `real(dp)`, `private` and `only:` lists.

  Below that the code keeps the F77's shape. Left:
  - `kap_eval` is a 310-line if-chain (kap_lib.f90:88-401). It does not read as named steps, the Hubbard & Lampe fit is inline, and its comments speak in legacy flags (`ldifz=T or lzramp=T`, `l2z=T`, :120-123, :256-258);
  - `compute_opal95_opacity` is 300 lines. The 2D/3D/4D interpolators repeat the same log R extrapolation block three times (opal95_opacity.f90:753-768, :835-850, :931-946);
  - `compute_laol_opacity` and `compute_pure_z_opacity` repeat the row and T-spline logic (laol_opacity.f90:277-377 against :527-608);
  - history in the code: 41 comment lines carry authors' initials and dates (MHP 7/98, DBG 5/94, JCZ 211125, YCK 3/91, llp), and 63 carry 2026 refactoring notes ("2026 wave 3 (R5): ...", "token-for-token");
  - 9 commented-out statements (3 more in atm; finding 14);
  - 8 all-capitals comment blocks; 93 `endif`/`enddo`; 64 lines over 110 characters;
  - unnamed thresholds and constants (kap_lib's 0.1, 7.0, 0.12, 4.2 and the Hubbard & Lampe coefficients; conductive's -3.5194, which is log10(16σ/3), and 3/9/9.75/-6);
  - names that misstate their contents: `kurucz_table_set%log10_opacity` holds κ, not its log (kurucz_opacity.f90:121); `zlaol_logrho_grid` holds ρ; `use_alex95_tables`/`alex95_table_paths` feed `alex94_opacity`/`paths%alex94`.
- **Organization 7.** One facade (enforced), star-free, one settings record and one paths record built in setup/, and one module per family. Left:
  - the domain's whole state is one public, writable module variable, `opacity_table` (opacity_tables.f90:267). It holds the tables, the warm-start cursors, the composition caches and OPAL95's per-call scratch (`opal95_weight_*`, `opal95_logr*`, `opal95_extrap_*`, :249-264), which the 2026 change made into arguments of the pure interpolators but still stores globally;
  - `kap_envelope_metal_fraction` is further module state in kap_lib (:55);
  - no messages module: 34 terminal writes (`write(*,*)`, not `terminal_unit`) and 24 run-log writes in 7 files (finding 9; M112 for the run-log part);
  - family precedence is coded three times as if-chains (`load_opacity_tables`, `refresh_surface_tables`, `kap_eval`); TARGET_ARCHITECTURE 4.5(c)'s `kap_source` is not started;
  - `kap_eval` is exported for the test alone (kap_lib.f90:41-49), and the README lists it as interface.
- **Maintainability 5.5.** Left:
  - the test is a byte pin with no ok/FAIL check. It never takes OPAL95's 3D path, the production interior path for a deck at its table's Z (finding 10);
  - open in bug_sweep: O18's ramp derivatives (kap_lib.f90:333-338), M104 (opal92_opacity.f90:226), M170 (alex06_opacity.f90:62) and L10's remaining gap (opal95_opacity.f90:169);
  - a missing table ends the run in a Fortran runtime error at 6 of kap's 9 opens (M170 records one), and OPAL92's two create the missing file and let kap_init succeed (finding 2);
  - the Alexander 2006 cache can hand back a table for another composition with ierr 0 (finding 1);
  - unchecked settings give NaN opacities with ierr 0 (finding 6);
  - nine pre-divergence defects are open and are recorded only in audit/ and ROADMAP §10-11, not in bug_sweep. One of them, the conductive opacity's extrapolation sign, is reached in the photosphere of all 1,164 shipped decks (finding 4).
- **Extensibility 6.** README §5 has a recipe for a table family, and the settings and paths records make a family's configuration a copy in one builder. Left:
  - the recipe, traced, names about half the hand edits (finding 11);
  - a new family is about 8 hand-edited files plus 6-7 regenerated or re-pinned ones; TARGET_ARCHITECTURE 1.1 counts "~18";
  - precedence must be threaded through three if-chains.
- **Docs 6.** The README's table and flow are right on the main paths. Left:
  - no conventions for errors and messages; no coverage figure; no test or mutation claims; an incomplete recipe;
  - user-facing text a reader would act on is wrong. Every run log says the pure-Z table is "pure C and N"; it is C and O (finding 8). The LAOL extrapolation tolerance is described as log10 while the code takes ln (finding 4);
  - about 20 stale comments (finding 13).

Average **6.2**.

### atm/

- **Readability 6.5.** t_tau_relations (pure, three functions, star-blind), atm_lib and atm_tables read at the bar, and kurucz_atmosphere names its thresholds. allard_atmosphere is 886 of atm's 1,520 lines, and it is the F77 in modern syntax. Left:
  - its headers are prose in the F77's names: KTTAU, LPRT, ISHORT, IMODPT, FALLARD, IOATMA, LAlFail, TEFFLs, GLs (allard_atmosphere.f90:30-111, :332-361, :704-724);
  - typos ("assoiated", "asjust", "teh", "Likewist", "weanticipate", "recoreds", "heve");
  - `fail_9999` and "Program Terminated" banners;
  - the kept misspelt `latmtptau100` (O19);
  - hard-coded 54 and 5 for the old layout (:429-430);
  - unnamed T-tau constants: -0.031235 = ¼ log10(3/4) and the two photospheric log τ (t_tau_relations.f90:40, :42, :59, :61, :76, :78);
  - `print_flag`/`print_to_files` are always `.false.` from atm_lib (dead parameters);
  - 53 `endif`/`enddo`, 21 lines over 110 characters, and history (G Somers 5/15, MHP 6/97, JNT 6/14).
- **Organization 6.5.** One facade (enforced), star-free, a settings and a paths record from setup/, and a star-blind T-tau kernel. Left:
  - the table lookups return their result through module state: `interpolate_surface_pressure` and `interpolate_allard_surface` write `atm_table%atm_log10_pressure`/`_temperature`, which `atm_get_surface_pt` copies out (atm_lib.f90:100-101, atm_tables.f90:57);
  - no messages module: 48 terminal writes and 66 run-log writes. The Allard kernel notes every successful lookup, 26% of the pinned Allard run log (finding 9);
  - what each `atm_choice` means is spread over four places outside atm: the envelope kernel's tabulated list `3 .or. 4 .or. 5` (core/envelope/envelope_integration_kernel.f90:306), surface_boundary's fallback switch with its literal 3.95 (core/structure/surface_boundary.f90:136-174), the run-log banner (io/derived_controls.f90:327-339) and the control's description;
  - nothing validates the choice (finding 6);
  - atm_init reads the Kurucz table for choice 4, and nothing uses it (finding 12).
- **Maintainability 6.** Left:
  - the test is a byte pin with no ok/FAIL check (finding 10);
  - O19 is open, both halves (allard_atmosphere.f90:375, :417; core/step/mass_change.f90:251-259), and it reaches further than recorded (finding 7);
  - a deck whose Allard [Fe/H] or [α/Fe] is not in the file ends in a runtime error (finding 3);
  - five pre-divergence defects are open and unregistered: the Castelli/Kurucz sentinels above about 48 kK, no upper range check, the Allard extension width, the gravity-search gap, and the per-lookup run-log notes (finding 5);
  - M112's direct run-log writes.
- **Extensibility 6.** README §5 has recipes for a table and a T-tau relation. Traced, the table recipe misses the core files that know the choice numbers, the path control and its record, and the fallback switch. A forgotten T-tau branch falls silently into HSRA through the functions' `else` (finding 11).
- **Docs 6.5.** README §1-§2 match the code (the rows, the choice table, the fallback to gray). Left:
  - §3 says nothing of ierr or where messages go; no coverage figure; no test claims beyond "pins"; incomplete recipes;
  - kurucz_atmosphere's header promises an error beyond one table step outside the table, which exists only at the low edges;
  - the Allard header gives one column's width where the code takes four, and says a lookup returns both P pairs (finding 13).

Average **6.3**.

## Summary

- **kap and atm meet the bar's boundaries, not its insides.** Both are behind a facade that check_boundaries enforces. Neither touches the star. Each takes one settings record and one paths record, built by pure copies in setup/ and installed from settings_install. The T-tau relations are a pure, star-blind kernel. Every family module was renamed and given `real(dp)`, `private` and `only:` lists in October, and the second-Z clones that the September sweep counted at about 1,100 lines became table-set types.
- **The standard's middle layers are missing.**
  - **No behaviour tests.** test_kap and test_atm are byte pins: 64 and 56 lines, 0 ok/FAIL checks, no mutation claims, no coverage figure. test_kap never runs OPAL95's 3D path, `kap_update_surface_tables`, or the conductive opacity with nonzero ion fractions.
  - **No messages module.** 58 kap and 115 atm direct writes are split between `*` and the run log, in capitals, some saying "RUN STOPPED" or "Program Terminated" on paths that return ierr.
  - **State.** OPAL95's per-call scratch sits in the public `opacity_table` record, and atm's lookups return through module state.
- **The defect register misses most of the domain's open defects.** bug_sweep has four kap items (O18's ramp, M104, M170, L10's gap) and one atm item (O19). Fourteen more found by the 2026-08-31 and 09-01 sweeps are open at 823326a7 and are not in bug_sweep. Some are in ROADMAP §10-11 or MODERNIZATION_CHANGES §5; some only in audit/. I re-verified each in the code and against the F77 (findings 4 and 5). The most consequential:
  - the conductive opacity's low-density branch flips the sign of log σ and extrapolates the wrong way, and the solar photosphere takes it in all 1,164 decks (negligible in κ, since conduction is negligible there);
  - choice 5 (Castelli/Kurucz) interpolates -999 sentinels above about 48 kK, silently;
  - OPAL95 reads an unset `density_shifted` for X < 0.2 at log T 4.0-4.05.
- **New possible bugs, all latent in the shipped decks:**
  - Alexander 2006's cache can serve a table for another composition with ierr 0, after a failed build or after a second kap_init in one process (pyyrec, test_reentry's two-deck form);
  - 6 of kap's 9 table opens end in a Fortran runtime error on a missing file, contrary to M170's "every other table reader reports" (M170 names one). OPAL92's two opens create the missing file, and kap_init succeeds;
  - an Allard [Fe/H] that the file lacks ends in a format runtime error;
  - unvalidated settings give silent NaN opacities, or the HSRA atmosphere for an unknown `atm_choice`;
  - accretion looks up Allard tables that atm_init never loaded unless `atm_choice = 4`.
- **User-facing text is wrong in two places.** Every run log of the 1,164 decks says "opacity: pure C and N tables enabled", but the table is pure C and O. The LAOL extrapolation tolerance is documented in log10 and applied in ln.
- **What it would take to reach rotation's bar:**
  - an ok/FAIL kernel test per domain (finite-difference derivatives, the ramp's continuity, exactness at table nodes, error contracts), each confirmed under a mutation, with a coverage figure;
  - kap_messages and atm_messages;
  - the open defects registered, and the byte-identical ones fixed (most are);
  - `kap_eval` split into named steps;
  - the README conventions and the recipes traced.

## New findings (not in bug_sweep)

1. **kap/alex06/alex06_opacity.f90:304-309 and :148-152 (cache: opacity_tables.f90:195)**
   - **Issue:** `compute_alex06_opacity` writes the requested (X, Z) into `alex06_cached_x/_z` before `build_alex06_composition_table` runs. The build can fail: X < 0, X > 1 − Z, Z < 0 or Z > 0.12 (:175-191) returns ierr = 1 before the table is touched. The cache then claims a composition the table does not hold. The next call at the same (X, Z) finds no difference (:304), skips the build, and returns the previous composition's opacity with ierr 0. That violates the ierr contract: a retry succeeds where the first call failed.
     - A second path: `read_alex06_tables` re-initializes the table to the (X = 0, Z = 0) slice (:148-152) without resetting the cache. After a second `kap_init` in one process (a second `run_yrec` call: pyyrec, `test_reentry`; `prepare_run` installs the settings once per call, core/run/run_driver.f90:146), the first lookup at the previous run's last surface (X, Z) uses the X = 0, Z = 0 molecular opacity.
     - The second path needs that composition to repeat within 1e-8, which a non-diffusing deck that starts past deuterium burning does. Alexander 2006 is in all 1,164 decks, and every one runs kap_init once per process, so no single run reaches it.
   - **In the F77?** No. GETALEX06 also set XE/ZE before ALEX06TAB (getalex06.f:18-22), but ALEX06TAB stopped on a bad composition (alex06tab.f:18-27: `WRITE(*,5) XE,ZE` / `STOP`). The failure became returnable with 073b6f92 (2026-08-21, "kap domain converted to the ierr-not-stop discipline"), and the re-read path came with 2026 re-entrancy.
   - **Fix:** set the cache after a successful build (or reset it to an impossible value on failure), and set `alex06_cached_x = 0`, `alex06_cached_z = 0` where `read_alex06_tables` loads the (0, 0) slice. Add test_kap checks: a failed composition fails twice, and a second kap_init followed by the same lookup gives the same κ.
   - **Kind:** A (byte-identical for every single-run deck). **Size:** S.

2. **kap/opal95/opal95_opacity.f90:61; opal92/opal92_opacity.f90:60, :66, :55, :87-107; laol89/laol_opacity.f90:92-93, :412-413; kurucz90/kurucz_opacity.f90:90; alex94/alex94_opacity.f90:67-68 (with M170's alex06_opacity.f90:62)**
   - **Issue:** M170 records the Alexander 2006 open alone and says "every other table reader reports and exits 1 (L5's rule)". At 823326a7 that holds only for the conductive table (conductive_opacity.f90:51) and atm's two readers. Five more opens have `status='OLD'` and no `iostat`, so a missing OPAL95, LAOL89, pure-Z, Kurucz or Alexander 1994/95 table ends the run in "Fortran runtime error: Cannot open file" (exit 2), as M170's Alexander 2006 does.
     - The two OPAL92 opens have no `status` either. A missing path is created empty where the directory allows. `read_opal92_table` meets end-of-file at once (`end=97`) and sets `num_x = 0` and `num_temps` from `num_temps_read`, an unset local (0 under `-finit-local-zero`, stack contents under ifort/ifx). `kap_init` then returns 0, and the run fails at its first lookup with "compute_opal92_opacity: T out of table".
     - A truncated OPAL92 file returns success with fewer temperature rows. LAOL89's 18 body reads, Alexander 2006's six and the conductive table's three have no `iostat` or `end=`.
     - No shipped deck uses OPAL92, LAOL89, Kurucz 1990 or Alexander 1994/95; all 1,164 use OPAL95, pure Z and Alexander 2006.
   - **In the F77?** Yes: setllo.f:50 `OPEN(UNIT=ILLDAT,FILE=FLLDAT)`, ll95tbl.f:48, rdlaol.f:32 and :80, rdzlaol.f:32, setkrz.f:45 and :76, and alxtbl.f:59, none with IOSTAT.
   - **Fix:** one `open_table(path, unit, label, ierr)` in kap that opens with `status='old', action='read', iostat, iomsg` and reports through kap's messages (finding 9). Every read gets `iostat`. `num_temps_read = 0` before the loop, and a short OPAL92 table is an error. M170's text should be corrected to name the eight unchecked opens.
   - **Kind:** A (error paths). **Size:** S/M.

3. **atm/tables/allard_atmosphere.f90:483-498, :738-742, :795-799**
   - **Issue:** `read_allard_tables` keeps only the records at the deck's `allard_target_feh` and `allard_target_alpha`. A value the file does not hold leaves `allard_num_teff = 0`. Nextgen2.all has [Fe/H] 0, -0.5, -1, -1.3, -1.5, -2, -3, -3.5, -4 and [α/Fe] 0 only, so -0.3, or an α of 0.4, is enough.
     - `find_allard_table_bounds` then reads `allard_teffl_grid(allard_num_teff)` and `(allard_num_teff-1)`, that is indices 0 and -1 (:739-742).
     - Its "Less than 4 rows" message is written with `910 format(a,i3)` against two character items and an integer (:797-799), so gfortran stops with a format runtime error instead of the diagnostic.
     - The same format error hides the message for any table with fewer than four rows. 116 shipped decks use Allard, all at the defaults [Fe/H] = 0 and [α/Fe] = 0 (no deck sets either).
   - **In the F77?** Yes: altabinit.f:104-106, `write(ISHORT,910) 'ALTABINIT: Bad input Allard Table: ',` with `910 format(A,I3)`.
   - **Fix:** after pass 1, `allard_num_teff == 0` is an error naming the [Fe/H] and [α/Fe] asked for and those the file has. Change the format to `(2a,i3)`.
   - **Kind:** A. **Size:** S.

4. **Nine pre-divergence kap defects, open at 823326a7 and not in bug_sweep** (found by audit/bugsweep-2026-08-31/kap.md and 2026-09-01-pass2/K_opacity.md; some carried into ROADMAP §10-11 or MODERNIZATION_CHANGES §5). Each was re-read in the code and against the F77.
   - **a. The conductive opacity's low-density branch** (conductive/conductive_opacity.f90:217-226 against :170 and :195).
     - The main branch takes `mix = -log10(Σ w_i/K_i)`, while the log ρ < -6 branch takes `+log10(...)`, so log κ_c changes sign in the conductivity term. It also extrapolates with `+QODC*(-6 - log ρ)`, the wrong way.
     - The branch runs at log ρ < -6 with log R ≥ 0, which includes the solar photosphere (log T 3.76, log ρ ≈ -6.6, log R ≈ 0.1) in all 1,164 decks (use_conductive_opacity in every one).
     - Conduction is negligible there either way (κ_c ≳ 1e6 against κ_rad ≲ 1), so κ moves by about 1e-6 relative.
     - F77: condopacpint.f:57 against :85 and :89-92.
     - Recorded in ROADMAP:663-664 and MODERNIZATION_CHANGES:420-421. **Kind B (tiny), S.**
   - **b. OPAL95's low-X, low-T shortcut** (opal95/opal95_opacity.f90:414, :461-505, :508; opacity_tables.f90:264).
     - For X < 0.2 and `opal95_index_t ≤ 5` (log T below about 4.05; the interior table is consulted there in the ramp, log T ≥ 4.0 by default), `density_shifted` is read unset at :508, and `opal95_extrap_hi`/`_row` keep the previous call's values.
     - `low_regime_flag` is set and never read. The check it postpones ("since the empty region depends on Z ... postpone") was never written.
     - F77: LSHIFT and LEXTH were SAVEd (getopal95.f:53, :134, :180), so the F77 read the previous call's values. The unset local is 2026.
     - Recorded in ROADMAP:441-442. **Kind A for `density_shifted = .false.` (GNU's `-finit-local-zero` gives that already); resetting `opal95_extrap_hi` changes that path only. S.**
   - **c. OPAL95's fixed-Z table rejects 0.08 < Z ≤ 0.10** (opal95_opacity.f90:211-225).
     - The search runs to `n_opal95_z - 1` (Z = 0.08), so an `opal95_single_table_z` in (0.08, 0.10] stops kap_init with "OUTSIDE OP95 TABLE RANGE ... RUN STOPPED". The tables reach 0.10 and the 4D path handles it.
     - The comment says "Desired Z > 0.1d0".
     - F77: op95ztab.f:30, :36. **Kind A for shipped decks, S.**
   - **d. Alexander 2006's fourth X node** (alex06/alex06_opacity.f90:263-267 against :237).
     - For Z > 0.06 the node is 1 − Z whatever `alex06_index_x` is, while the data at that position is the X = 0.8 table when X ≤ 0.8. The branch at :237 tests both indices; :263 tests only `iz`.
     - F77: alex06tab.f:76 against :100. Latent (Z > 0.06 in an atmosphere). **Kind A for shipped decks, S.**
   - **e. The helium-core branch has no `else`** (kap_lib.f90:205-233). With `use_pure_z_table` and none of OPAL95, OPAL92 or LAOL89 selected, `table_metal_fraction` and the interior κ are undefined, and the blend returns garbage with ierr 0. F77: getopac.f:81-92. **A, S** (with finding 6).
   - **f. The LAOL89 extrapolation tolerance** (laol89/laol_opacity.f90:262, :522).
     - `log(laol_extrapolation_tolerance)` is a natural log applied on a log10 ρ axis: 2.30 dex for the default 10.
     - The 2026 description says "log10 of this value" (defaults/physics.defaults:461, io/controls_buffer.f90:809-810).
     - F77: gtlaol.f:28 `TOLL = LOG(TOLLAOL)`. **A to correct the description, B to change the code; no shipped deck uses LAOL89. S.**
   - **g. The second LAOL89 table's T grid is logged over the first table's row count** (laol_opacity.f90:173-179). The code comment leaves it "an R6 decision". **A, S.**
   - **h. Kurucz 1990's density column is not bounded** against `kurucz_max_num_densities` (kurucz90/kurucz_opacity.f90:120-122). A malformed row is skipped without a word (:100). **A, S.**
   - **i. Alexander 2006's header check** (alex06_opacity.f90:73-78).
     - The X, Z check in the main loop prints "RUN STOPPED" and carries on (`!  STOP` commented out), while the same check at :117-123 returns ierr.
     - The log R check covers 16 of the 19 columns (:81, :125).
     - The shipped table passes every check (alex06_headers.py). F77: readalex06.f:69-71. **A, S.**
   - **Fix:** enter a-i in bug_sweep. The fixes of b-e and g-i are byte-identical for every shipped deck; a's moves κ by about 1e-6 in the photosphere, and f's code option moves LAOL89 decks only (none shipped).

5. **Five pre-divergence atm defects, open at 823326a7 and not in bug_sweep** (found by the 08-31 atm_net_math.md and the 09-01 L_atm_controls.md).
   - **a. Castelli/Kurucz sentinels** (atm/tables/kurucz_atmosphere.f90:89-109).
     - Rows above log Teff 4.5 interpolate the top two or three log g columns without `gmin_index`/`gmax_index`.
     - Every shipped Castelli/Kurucz table's last row (log Teff 4.699) holds -999 at log g 4.0 and 4.5 (ck_sentinels.py: atmk2004p00.tab and m10), so choice 5 returns garbage once the 4-row stencil reaches that row (Teff above about 48 kK). No fallback switch exists for choice 5.
     - F77: surfp.f:50-51. Recorded in MODERNIZATION_CHANGES:401-403. No shipped deck uses choice 5. **A for shipped decks, S.**
   - **b. No upper range check** (kurucz_atmosphere.f90:28-33, :49-53, :59-71).
     - The header promises ierr more than one table step outside the table, but only log Teff < 3.5 and log g < -0.5 are rejected. Above the top row or column the splines extrapolate without limit.
     - F77: surfp.f:25. **A, S.**
   - **c. The Allard extension width** (atm/tables/allard_atmosphere.f90:771 against :41-46 and :716-719). The code extends a row four column widths below its first valid log g; both headers say one. F77: altabinit.f:78. **A (doc), S.**
   - **d. The general-case gravity search** (kurucz_atmosphere.f90:135-173).
     - The loop starts at `gmax_index-3` and tests `log_g < logg(k+2)`, so a log g between `logg(gmax_index-1)` and the top falls through to the first four points and extrapolates from the low end.
     - A row of all -999 (`gmin_index = ng`, `gmax_index = 0`) skips the loop and reads columns ng..ng+3, out of bounds.
     - Latent with the shipped tables. **A, S.**
   - **e. The Allard kernel notes every successful lookup in the run log** (allard_atmosphere.f90:311-324). The `.not. print_to_files` branch, which is every call, writes two lines. That is 960 notes, 1,920 of the 7,383 lines, in testsuite/standard/Test_solarGS98_allard.log. **B (log text), S** (with finding 9).
   - **Fix:** enter a-e in bug_sweep. For a and d, use the row's valid range (or stop above the last full row); for b, add the upper checks.

6. **kap/kap_lib.f90:411-446 and atm/atm_lib.f90:43-67: the settings are installed unchecked**
   - **Issue:** `kap_init` and `atm_init` copy the record and load tables, but check nothing. GUIDELINES' model, `set_net_settings`, checks first and returns ierr. Consequences, all with ierr 0:
     - (a) `use_two_z_tables` with equal table Zs gives 0/0 in `blend_in_z` (kap_lib.f90:479-480), a NaN opacity. The controls' defaults are `opal_table_z1 = opal_table_z2 = 0` (likewise LAOL and Kurucz; defaults/control.defaults:311-347), so a deck that turns on the Z ramp or diffusion of Z with one of those families and leaves the Zs gets NaN;
     - (b) `molecular_opacity_logt_min >= _max` gives 0/0 in the ramp weight (:331-332);
     - (c) the helium-core branch without an interior table (finding 4e);
     - (d) `atm_choice` outside 0-5 loads no table, prints no banner line (io/derived_controls.f90:327-339), and runs the HSRA relation through t_tau_relations' `else` branches (t_tau_relations.f90:44-46, :62-64, :77-79).
   - **In the F77?** Partly. GETOPAC had the same if-chains, and an unknown KTTAU left TL stale in QATM and ENVINT (qatm.f:48-54, envint.f:166-171). The silent HSRA is the 2026 `else`.
   - **Fix:** a `check_kap_settings` and a `check_atm_settings`, called by the two inits, returning ierr with a message. One interior table when `use_pure_z_table` is set; distinct Zs under `use_two_z_tables`; min < max; atm_choice in 0..5. t_tau_relations' `else` becomes `atm_choice == 2`.
   - **Kind:** A (no shipped deck trips any of them). **Size:** S.

7. **core/step/mass_change.f90:251-259 with atm/atm_lib.f90:80-102 (O19's second half, wider than recorded)**
   - **Issue:** O19 says the accretion surface "uses a possibly failed Allard lookup". In fact the lookup always fails unless the deck's `atm_choice` is 4, because `find_accreted_surface_state` calls `atm_get_surface_pt(4, ...)` whatever the choice, and atm_init loads the Allard tables only for 4.
     - The unloaded record's limits are 0, so every accretion step writes "ALSURFP: TEFFL greater than TeffLmax ... 0.0 3.65 0.0" and "GL out of max range" (4 lines) to the run log and the terminal.
     - It then takes `s%atm_point`'s last values. The untracked examples/run_config_matrix/output/cm_accretion.log, written by 823326a7 on a gray-atmosphere accreting deck, has 24 such lines.
     - `atm_get_surface_pt` cannot tell "outside the table" from "table not loaded".
     - No shipped deck accretes.
   - **In the F77?** Yes: massloss.f:133-140 calls ALSURFP unconditionally ("valid for Allard atmospheres only"), and setups.f:250 loads the files for KTTAU = 4 only.
   - **Fix:** atm_init records which tables it read, and `atm_get_surface_pt` returns a named error for a choice whose tables are not loaded. core either requires `atm_choice = 4` for accretion or takes the envelope integration's photosphere. Note it on O19.
   - **Kind:** A (no shipped deck). **Size:** S.

8. **io/derived_controls.f90:341, defaults/control.defaults:176, defaults/controls_registry.md:81, io/controls_buffer.f90:323 (and kap/README.md:15, kap/laol89/laol_opacity.f90:5-6)**
   - **Issue:** the run-log banner says "opacity: pure C and N tables enabled", and the control's description says "Use the pure C and N opacity tables ... for metal-enriched mixtures". The shipped table, input/opacity/atomic/lanl/PURECO.DBGLAOL, is a "Pure mixture of equal numbers of C and O" (its first line; mass fractions C 0.4288, O 0.5712, N 1e-20). kap uses it only in helium-burning cores (Z > 0.1, log T > 7).
     - All 1,164 decks set `use_pure_z_table`, so every run log carries the wrong line, among them 35 pinned logs.
     - kap's README and laol_opacity's header describe the table as "(X = 0, Y + Z = 1)", while kap_eval places it at Z = 1.
   - **In the F77?** Yes: parmin.f:1263, `WRITE(ISHORT,*) ' USING PURE C AND N OPACITY TABLES'` (and :694, "PURE CN").
   - **Fix:** "pure C and O" in the banner and the description; "for helium-burning cores" in the description; "Z = 1 (C and O)" in the README and the header.
   - **Kind:** A for the descriptions; B for the banner (35 pinned logs move by one line). **Size:** S.

9. **kap and atm write their own messages (all family modules, allard_atmosphere and kurucz_atmosphere)**
   - **Issue:** GUIDELINES code style 7 asks for a domain's messages through one module and for kernels that print nothing; eos, net and rotation have one. kap has 34 `write(*,*)`/`print*` and 24 `write(run_log_unit,...)` in 7 files; atm has 48, 66 and 1 to `terminal_unit` (counts.py). M112 covers the 90 run-log writes; the rest is not recorded. Consequences:
     - `write(*,...)` bypasses `terminal_unit`, so test_atm's pinned stdout carries an Allard diagnostic (expected_test_atm.out:26-28), and an embedding host cannot redirect it;
     - the reason for most failures never reaches the run log. OPAL95's "LOG T ... OUT OF OPAL 95 TABLE RANGE", OPAL92's, Kurucz's spline checks, the conductive range and Alexander 2006's go to the terminal only;
     - messages in capitals say "RUN STOPPED" (kap_lib.f90:317-319, opal95_opacity.f90:220-221, alex94 :83 and :103, alex06 :76, :85, :100, :179, kurucz_atmosphere.f90:67-68) or "Program Terminated" (allard_atmosphere.f90:186-192, :648-657, :818-823) on paths that return ierr;
     - finding 5e's per-lookup notes.
   - **Fix:** kap_messages and atm_messages (report = terminal + run log, note), the kernels returning ierr and a short reason, and the domain's entry point reporting once. Plain-text messages.
   - **Kind:** B (run-log text). **Size:** M.

10. **test_kap.py / kap/test/test_kap.f90 and test_atm.py / atm/test/test_atm.f90: pins without behaviour checks**
    - **Issue:** both programs print numbers that the .py files byte-compare (64 and 56 lines). No line prints ok/FAIL, nothing claims a mutation, and neither README gives coverage, although tools/coverage.py takes `--dir kap` or `--dir atm`. Paths the production decks take, untested:
      - OPAL95's 3D interpolation (`interpolate_opal95_3d`, `locate_opal95_x_stencil(1, ...)`). Every point is either at the surface composition (2D) or at another Z (4D);
      - `kap_update_surface_tables`;
      - the conductive opacity with nonzero ion fractions (`fxion = 0` throughout);
      - the derivatives, against which O18's ramp error, finding 4a's extrapolation and the 2026-09 conductive derivative fix would all show.
    - Error contracts beyond two kap points and one atm point are untested: findings 1-3.
    - test_kap boots its first section through the star layer's `set_up_constants_and_tables` and a `star_info` (test_kap.f90:18-23, :113). Its "error paths via optional ierr" header is kept only because it is pinned (:132-137).
    - **Fix:** a kernel test per domain, rotation-style, each check confirmed under a mutation, plus a coverage figure in each README. For kap:
      - dlnκ/dlnρ and dlnκ/dlnT against centred differences of log κ for every family and across the ramp;
      - κ continuous at `molecular_opacity_logt_min`/`_max` and across the 2D/3D/4D switch (the surface composition against Z(1 + 2e-4));
      - exactness at table nodes (the production family already returns -4.322000000000000 at a node);
      - `blend_in_z` exact at z1 and z2;
      - the harmonic sum;
      - each error contract.

      For atm: T(τ_ph) = Teff for Eddington and Krishna-Swamy to 1e-6 (now 3.2e-7 and 2.5e-5 off, from the truncated constant and τ = 0.312); exactness at the Kurucz and Allard nodes; the error paths.
    - **Kind:** A (tests). **Size:** M.

11. **kap/README.md:46 and atm/README.md:46, the recipes, traced**
    - **Issue:** kap's "Add a table family" names the module, `opacity_table`, the `use_*` flag in `kap_settings_t` and `kap_settings_from_star`, the three branches and `check_families`. It omits:
      - the path field in `kap_table_paths_t` and `kap_table_paths_from_star`;
      - the controls themselves (io/README.md §5: the flag and the path, with their registry rows and settings cells);
      - `python3 tools/gen_makefile_deps.py` (eos/README.md:125 says it; kap's does not);
      - where the family goes in each precedence chain;
      - the README row and the re-pins (expected_test_kap.out, test_controls' pristine.txt).

      That is about 8 hand-edited files of 15 changed; the recipe names 4.

      atm's "Add a table of model atmospheres" omits:
      - the path control, `atm_table_paths_t` and `atm_table_paths_from_star`;
      - the envelope kernel's tabulated list (core/envelope/envelope_integration_kernel.f90:306), without which a new table is never consulted;
      - surface_boundary's fallback switch (core/structure/surface_boundary.f90:136-174);
      - the run-log banner (io/derived_controls.f90:327-339);
      - the control's description.

      "Add a T-tau relation" omits the banner and the description. A missed branch falls into HSRA silently (finding 6d).
    - **Fix:** list every file, with io/README §5 for the controls, and name the choice numbers' homes. Better: give atm one function that says whether a choice is tabulated, and use it in the kernel.
    - **Kind:** A. **Size:** S.

12. **kap/opacity_tables.f90:249-267 and atm/atm_tables.f90:55-57, atm_lib.f90:51-60**
    - **Issue:** `opacity_table` is `public` and writable. Besides the tables it holds OPAL95's per-call scratch (`opal95_weight_*`, `opal95_dweight_*`, `opal95_logr*`, `opal95_extrap_*`, "recomputed fresh each call"). `compute_opal95_opacity` computes these per call, stores them globally and passes them to the pure interpolators as arguments, so they could be locals.
      - In atm the two interpolators return their result through `atm_table%atm_log10_pressure`/`_temperature`.
      - `atm_init` reads the Kurucz table for choice 4 (:51-54), but nothing reads it: `atm_get_surface_pt(4)` uses Allard only, and the fallback from 4 is gray (surface_boundary.f90:149-158). An Allard deck must still name a Kurucz table that is never used, and fails if it is missing.
    - **In the F77?** Yes: COMMON /OP95FACT/, /OP95EXT/ and /ATMPRT/.
    - **Fix:** move the OPAL95 scratch into a local record of `compute_opal95_opacity`; make the interpolators return log P and log T in arguments; drop the Kurucz load for choice 4 (and the README's "with the Kurucz table loaded as well"); make `opacity_table` `protected` outside its writers, or split tables from cursors.
    - **Kind:** A (the Kurucz load: error path only). **Size:** M.

13. **Stale or inexact text**
    - **Issue:**
      - kap_lib.f90:191-193, "Z>0.10 T > 5 X 10^7 K", against the test `log10_temperature > 7.0d0` (:189). The threshold moved from 7.7 to 7.0 in the YREC team's 5b1a5ed6 (2025-12-17); the companion a8e60674 updated only the other message (:249);
      - opal95_opacity.f90:427-428, "the facade (kap_lib's kap_get) stops when its caller passes no ierr". `kap_get`'s ierr is required, and nothing stops;
      - opal95_opacity.f90:390-391, "cubic spline interploation", for a 4-point Lagrange scheme;
      - opal95_opacity.f90:218, "Desired Z > 0.1d0" (finding 4c);
      - opacity_tables.f90:239-243, "atm/turnover/acoustic_depths.f90 ... now maps onto ... this same type". No such file exists anywhere;
      - opal92_opacity.f90:39-41, "the first table on the fixed unit opal92_table_unit". Both are `newunit`;
      - kurucz_opacity.f90:37, "the controls' kurucz_table_unit";
      - routine-as-file references: opal92 :34, kurucz :34 and :141-143, laol :154-155 and :455, alex06 :160-162 and :283, conductive :80 and :83 (kap_get for kap_eval) and :238 (condall.d for condall06.d);
      - kap/README.md:19 lists `hunt` among the conductive module's routines; it is numerics_lib's;
      - kap/README.md:11 presents `kap_eval` as interface, while check_boundaries calls it domain-internal;
      - tools/check_boundaries.py:64-67, two contradictory comments ("kap_get_r is the ... variant of kap_get"; "kap_get is the ... (the former kap_get_r)");
      - kap README §2 omits the Z > 0.12 (or single-table Z) error branch, and that OPAL95's fixed-Z table is at `opal95_single_table_z`, not the envelope Z;
      - atm/tables/kurucz_atmosphere.f90:28-33, :49-53 (finding 5b);
      - allard_atmosphere.f90:30-33, "returns Log(P) ... as well as Log(P) and Log(T) at Tau = 100": it returns one pair;
      - allard_atmosphere.f90:41-46 and :716-719 (finding 5c);
      - allard_atmosphere.f90:111 (LPRT, ISHORT, IMODPT), :130-133 (KTTAU), :168 (LAlFail), :393 (FALLARD, IOATMA), the F77 unit and flag names.
    - **Fix:** text only.
    - **Kind:** A. **Size:** S.

14. **Code-style residue against GUIDELINES "Code style" 3-5 (no retired/ directory in either)**
    - **Issue:**
      - commented-out statements: alex06_opacity.f90:77 (`STOP`); conductive_opacity.f90:144-145, :172-173, :186-187, :192-193, :219-220, :311; opal95_opacity.f90:577-578; allard_atmosphere.f90:433, :438, :540;
      - history: 41 dated author notes and 63 "2026 ..." notes in kap, 7 in atm;
      - comparisons with `.le.`/`.ge.` (laol_opacity.f90:628, :631);
      - unnamed constants: the T-tau -0.031235, -0.505627854 and -0.176091259; conductive's -3.5194, which should come from σ in physical_constants; the Hubbard & Lampe fit; kap_eval's 0.1, 7.0, 0.12 and 4.2; the OPAL95 stencil's 0.04, 0.1, 0.8 and 0.2;
      - the literal surface slot 8 (alex94_opacity.f90:192, :264-277; `n_alex94_xt = 8` exists but is private);
      - the names in the Readability bullets.
    - **Fix:** a retired/ file per directory for the kept formulas, with their source; a Sources section per file for the history; named parameters with unit and source (naming is byte-identical; computing -3.5194 from σ is B at 1e-5); renames.
    - **Kind:** A (B for computing the constants). **Size:** M.

## What would raise each score

**kap/**
- **Readability, to 8-9:**
  - `kap_eval` as named steps: molecular, interior, helium-core blend, ramp, conduction, each a procedure returning a `kap_result` (TARGET_ARCHITECTURE 4.5(c)'s driver);
  - one log R extrapolation helper for OPAL95's three interpolators, and one row/T kernel for LAOL89's two;
  - findings 13 and 14;
  - Sources sections in place of the 104 history comments.
- **Organization, to 8.5-9:**
  - kap_messages (finding 9);
  - OPAL95's scratch made local, and `opacity_table` protected or split (finding 12);
  - one precedence list instead of three if-chains;
  - `kap_eval` private behind a test hook.
- **Maintainability, to 8-9:**
  - the ok/FAIL kernel test with mutations and a coverage figure (finding 10);
  - findings 1, 2, 4 and 6 fixed (all byte-identical but 4a and 4f's code option), with 4a-i entered in bug_sweep;
  - M104, M170 and L10's gap closed;
  - O18's ramp derivatives (a measured reseed).
- **Extensibility, to 8:** the recipe traced to every file (finding 11); a family descriptor (`kap_source`), so that a new family is one module, one list entry and two controls.
- **Docs, to 8-8.5:**
  - findings 8, 11 and 13;
  - README §3 conventions for ierr and messages;
  - a tests section that says what each check verifies, with coverage.

**atm/**
- **Readability, to 8-9:**
  - allard_atmosphere's headers rewritten in today's names and cut to what the code does;
  - `fail_9999` and the banners gone;
  - the old-layout dimensions named;
  - the T-tau constants named (¼ log10(3/4), log10 τ_ph);
  - the dead print parameters removed.
- **Organization, to 8.5-9:**
  - the lookups return through arguments (finding 12);
  - atm_messages, with no per-lookup notes;
  - one home for what a choice is (tabulated or not, its fallback) that the envelope kernel, surface_boundary and the banner ask;
  - settings checked at init (finding 6).
- **Maintainability, to 8-9:**
  - the ok/FAIL test (finding 10);
  - O19 fixed (both halves, and finding 7);
  - findings 3 and 5 fixed and registered.
- **Extensibility, to 8:** finding 11's recipes traced; the choice descriptor above.
- **Docs, to 8:** §3's error and message conventions; a tests section with coverage; findings 5b, 5c and 13.

## Traces, counts and README checks (supporting)

- **Decks** (decks2.py over git ls-files: 1,164 decks, each .inlist alone or .nml1 with its .nml2):
  - OPAL95, Alexander 2006, pure Z and the conductive opacity: 1,164 each;
  - OPAL92, LAOL89, Kurucz 1990, Alexander 1994/95: none;
  - `atm_choice` 0: 1,035, 4: 116, 3: 12, 1: 1; 2 and 5: none.

  So test_kap's family checks are the only coverage of four of the seven opacity families, and test_atm's the only coverage of choices 2 and 5.
- **Add a table family (kap), traced.** Hand edits:
  1. the module;
  2. opacity_tables (tables, flag, path);
  3. kap_lib (use line, three branches);
  4. kap_settings_builder (two copies);
  5. io/controls_buffer (two declarations);
  6. defaults/controls_registry.tsv (two rows with settings cells);
  7. test_kap;
  8. the README.

  Then gen_makefile_deps, make (which reruns the controls generators), `check_star_contracts --fix`, the re-pins of expected_test_kap.out and test_controls' pristine.txt. About 15 changed files; the recipe names 4.
- **Add a tabulated atmosphere (atm), traced.** The module, atm_tables, atm_lib (init, get), the paths record and builder, the path control and its registry row, envelope_integration_kernel.f90:306, surface_boundary (if the table has a limit), derived_controls' banner, the atm_choice description, test_atm: about 11 files, of which the recipe names 4.
- **Message sites** (counts.py, comments excluded):

  | | `write(*)`/`print*` | `write(run_log_unit)` | `write(terminal_unit)` | `report_*` |
  |---|---|---|---|---|
  | kap (7 files) | 34 | 24 | 0 | 0 |
  | atm (2 files) | 48 | 66 | 1 | 0 |

  Opens and reads:
  - kap: 9 opens, 1 with iostat; 34 reads without iostat or `end=`;
  - atm: 2 opens, both with iostat; 12 reads without.
- **Boundaries checked and true:**
  - `PUBLIC['kap']` = kap_get, kap_init, kap_update_surface_tables. Outside kap only test_kap calls `kap_eval` or `set_kap_settings`; core and convection call `kap_get`; setup calls `kap_init` and `kap_update_surface_tables`;
  - `PUBLIC['atm']` matches atm_lib's exports and its callers (core's envelope kernel, integrands, surface_boundary and mass_change; setup's settings_install);
  - STAR_FREE_DIRS has both directories, and STAR_BLIND_FILES has atm/t_tau_relations.f90 (the only kernel of the two that qualifies; kap's compute routines read and write the module's caches);
  - every one of `kap_settings_t`'s 21 fields and `atm_settings_t`'s 4 is copied by its builder, and both builders are pure with contract lines.
- **README claims checked:**
  - kap §1's rows, but for `hunt` and kap_eval (finding 13);
  - kap §2's four steps against kap_lib.f90:127-399, with the Z > 0.12 branch missing;
  - kap §4 against setup/;
  - kap §5's list against test_kap's sections;
  - atm §1's rows;
  - atm §2's choice table and the fallback against core's kernel and surface_boundary;
  - atm §3 ("set only when the lookup succeeds") against atm_lib.f90:84 and :99-101;
  - atm §5's list.
- **Open bug_sweep items of the two directories:** O18 (the ramp part), M104, M170, L10's gap and M112's 24 kap writes; O19 and M112's 66 atm writes. M128 reaches kap's conductive weight from the EOS side.
- **Not checked:** anything that needs a run. No coverage figure, no mutation, and the photospheric magnitude in 4a is estimated from the Potekhin table (log K ≈ 3.5 at log T 3.8, log ρ -6), not measured.

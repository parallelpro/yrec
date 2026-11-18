# core_b (core/model, core/run, core/diagnostics, main, core/README) -- eighth evaluation (2026-11-17, yrec-next 1179afbf, read-only)

Scope: the second half of core/: `run/` (8 files, 2,065 lines), `model/` (5 files, 3,544 lines), `diagnostics/` (4 files, 2,130 lines), `main.f90` (26), `core_messages.f90` (146, read because every message of the three directories goes through it), and core/README.md (236 lines), which this review owns as a whole. The other reviewer has structure/, step/ and envelope/.

- **Read whole:** every source file of run/, model/ and diagnostics/; main.f90; core_messages.f90; core/README.md; docs/ROADMAP.md's "Round 2 (core)" and "Round 2B" sections (:4451-4785); the seventh report (audit/module-evaluation-7-2026-11-10/core.md) and FINDING_MAP's core rows.
- **Sampled:** test_run_contracts.f90 (header and structure), test_model_kernels.f90, expected_test_regrid.out, the nine mutation manifests (counts per file), test_reentry.py and test_stop_conditions.py (headers and cases); step/step_driver.f90:180-300 and :640-660, step/timestep_control.f90:60-95, step/shell_physics.f90:60-150, step/burn_step.f90:70-85, envelope/envelope_integration.f90:60-90, :186-200 and :370-425 (what my directories call); io/README.md's generated-files table and its section 5 recipes; tools/check_star_contracts.py (header and scan), tools/check_style.py (metrics), tools/mutate.py (header).
- **Other inputs:** docs/bug_sweep.md (the table's open rows, and a grep before each finding); docs/GUIDELINES.md's star rules; the F77 at 6cd56739 (chkscal.f, setscal.f, main.f, starin.f, tauintnew.f, hpoint.f) for each bug's provenance.
- **Checks run (read-only; the tree is unchanged after them):** `check_boundaries.py` ("Domain boundaries OK"), `check_star_contracts.py` ("Star contracts OK (209 declared)"), `mutate.py core --check` ("179 of 179 mutations apply"). Nothing was built or run.
- **Counts:** scratch Python in /Users/yaguang/yrec-work/review8/review8_core_b/ (scan.py, scan2.py and their .txt outputs, lcore_decks.txt).
- Paths are under src/ unless they start at the repository root (docs/, audit/, test_*.py).

## Scores (1-10)

The seventh report scored core as a whole per criterion (readability 7, organization 7, maintainability 6.5, extensibility 6.5, docs 6.5; average 6.7) and each sub-directory by its average only (diagnostics 7.0, run 6.8, model 6.6). The per-criterion lines below therefore give "core-wide seventh -> this directory now", and each directory's last line gives its own average, seventh -> now.

### run/ (2,065 lines): 6.8 -> 8.0

- **Readability 7 -> 8.0.** No F77 shapes are left (0 `double precision`, 0 old operators, 0 `endif`/`enddo`; the seventh counted 24/23/20), the run's two calibration protocols are stated once at the call site (run_driver.f90:185-214), and the stops are a table. Left:
  - `run_yrec` is a 362-line host whose four steps are internal procedures sharing its static locals (run_driver.f90:251, :274, :329, :422); `begin_kind_card` and `check_star_calibration` keep ragged F77 indentation (run_driver.f90:330-416; calibration.f90:299-361);
  - unnamed values: the give-up count 15 and its 14 in the text (run_driver.f90:300-302), the 11 of the Monte-Carlo summary (monte_carlo.f90:271), `1.0d9` beside `years_per_gyr` (stop_conditions.f90:249), `mod(...,2)` for calibration's `star_calib_cards_per_cycle` (stop_conditions.f90:170), `do j = 1,4` for the rescale targets (calibration.f90:422, :526);
  - three lines over 132 columns (monte_carlo.f90:332, :349, :373).
- **Organization 7 -> 8.0.** Every message goes through core_messages, every procedure takes the star, every module is private, and the C entry no longer touches the host's units. Left:
  - file handling that belongs in io: monte_carlo's sample read, summary writes and rewinds (16 I/O statements, the style ratchet's only non-zero core counts with model/'s one), among them `rewind(run_log_unit)` (monte_carlo.f90:274, :280) beside core_messages' `rewind_core_log`;
  - events that reach the terminal and not the run log (finding 9): the solar calibration's Newton step, the central-abundance stop's line;
  - `run_yrec`'s blanket `save` (kept by decision, README section 3);
  - M312 (the run's `newunit` files stay open after a C call).
- **Maintainability 6.5 -> 7.5.** test_run_contracts (22 checks, 22 mutations over calibration, monte_carlo, yrec_capi and yrec_reset) covers the star calibration's last card, the Monte-Carlo run's start and the C entry; M268, M271, M34, M97, M115, M190 and M191 are fixed and I confirmed each in the code. Left:
  - findings 2 and 5 (both F77, one changes results);
  - no check with a mutation behind it for `stop_conditions` (deck tests only), `check_solar_calibration` or `run_driver`; run/ is 37.8% covered in the README's measure (calibration and yrec_capi lie outside it);
  - open in bug_sweep: O20 (decision), M312, M201 (tool: `check_stop_conditions`' contract says it writes nothing).
- **Extensibility 6.5 -> 8.0.** "Add a stop condition" is a traced recipe that I traced again against the tree (below): it names every hand-edited file. Left: an abundance stop is two `select case` edits with no recipe (stop_conditions.f90:313-335); a new run protocol means card-index arithmetic in three places (calibration, run_driver, monte_carlo's `run_index-2`).
- **Docs 6.5 -> 8.5.** README section 2.1 gives the run's flow in order and it matches `run_yrec`; headers give contracts and protocols. Left: finding 10's items in run/ (calibration's header omits `end_star_calibration_card`; M201's contract).

Average 8.0.

### model/ (3,544 lines): 6.6 -> 7.6

- **Readability 7 -> 7.5.** `rezone` (29 procedures, the longest 85 lines) and `read_starting_model` (21, the longest 107) are now drivers over steps with explicit arguments and a documented state record; no F77 shapes are left (the seventh counted 62/95/74). Left:
  - `rescale_model` is a 421-line host with 10 argument-less internal steps sharing some 25 locals (model_rescaling.f90:76-496), with five lines over 132 columns, up to 190 (:282, :335, :438, :473, :484);
  - the step bodies keep F77 layout and names: ragged indentation through starting_model.f90:725-763, :889-935, :968-1012 and rezoning.f90:169-224; `temp_scratch`, `iread`, `old_senv`; a once-through `do ... exit` loop (starting_model.f90:726-762); an `if` its own comment calls always true (:239-241);
  - `add_points_for_x_jumps` and `add_points_for_z_jumps` are one 85-line procedure written twice (finding 8);
  - literal slots and sizes (finding 11): `do species_idx = 5,11` three times, `(13,2)` and `(12,2)`.
- **Organization 7 -> 7.5.** regridding is star-blind, `find_shell_boundaries` and `compute_envelope_mixture` are pure kernels beside their drivers, and messages go through core_messages. Left:
  - model/ calls up into diagnostics/: `read_starting_model` builds the stitched model and the turnover time, and takes its shell masses from `observables`' copy of regridding's kernel (finding 8);
  - the mesh criteria read the star, so none is a kernel and none has a test;
  - `rescale_model`'s internal steps (above); the model file's `rewind` (starting_model.f90:233).
- **Maintainability 6.5 -> 7.0.** test_regrid (19 checks) and the new test_model_kernels (19) have a mutation per check (34 over regridding, shell_boundaries and model_rescaling's mixture kernel); O4's model side, M85, M86, M96, M126, M270 and O6's part A are fixed. Left:
  - finding 3 (the core extension: a NaN for one point, and a density profile it never shifts) and finding 6b;
  - the three drivers, 3,060 lines, have no ok/FAIL check: rezone's criteria and mesh, `read_starting_model`'s steps, `rescale_model`. Round 2's own regression there (M315, every deeper refit refused) was found by a count in round 2B, not by a test. Coverage: model/ 65.3%, starting_model 56.4%, model_rescaling 38.9%; no pinned deck rescales a horizontal-branch model;
  - open in bug_sweep: O6's and M176's parts B (decision), M202, M316, M317.
- **Extensibility 6.5 -> 8.0.** "Add a mesh criterion" is a traced recipe whose every name exists (`flag_fixed_points`, `append_flag_point`, `sort_flag_points`, `max_flag_points`, `flag_jumps`, the `FLAG-POINTS` line). Left: the species ranges by number (above); no recipe for a rescale target or a starting-model step.
- **Docs 6.5 -> 8.0.** Headers give each step's equations and the state records' fields. Left: finding 10's items in model/ (two stray contract blocks, three stale names, two comments that describe a defect instead of fixing it).

Average 7.6.

### diagnostics/ (2,130 lines): 7.0 -> 7.6

- **Readability 7 -> 7.5.** The observables are one procedure per theme with their formulas; the stitched model's columns have named indices and one-line kernels; no F77 shapes are left (the seventh counted 31/26/20). Left:
  - `turnover_from_interior_new` is a 354-line host: six argument-less internal steps sharing about 30 locals, driven by a `search: do ... cycle` loop (turnover_timescale.f90:118-471);
  - `build_stitched_model` keeps the F77's locals `b, gl, rl, ateffl, plim, lprt, lsbc0` (stitched_model.f90:116-118), and `build_pulse_points` five internal steps on shared scalars (:725-909);
  - unnamed values in observables (finding 11): `8.64d4`, `1.0d-5`, `9999.0d0`, `1.0d6`, the loop ends 8 and 10; the capture-rate tables are named `gasnuf_diag`/`clsnuf_diag` and carry no source (observables.f90:45-48).
- **Organization 7 -> 7.5.** The turnover refresh is now a stage of the step, the three diagnostic files are in `KERNEL_FILES` (no I/O), and the writers only read the stitched arrays. Left:
  - turnover_timescale.f90 is physics state (the wind and the deuterium cap read it), sits in diagnostics/, and returns its result through the star (star rule 2) from a walker that otherwise takes arrays;
  - `renormalize_luminosity_breakdown` changes model state inside the observables (finding 7);
  - `observables` exports a shell-mass kernel that model/ uses (finding 8); the stitch is the one writer of `photosphere_logP`, physics state, through `atm_get`.
- **Maintainability 6.5 -> 7.0.** O4 and O22 are fixed and a column without a value in a region is an error. Left:
  - findings 1 and 4 (a 2026 regression in the turnover time that a comment now describes; an envelope failure the stitch swallows) and 6a;
  - diagnostics/ has no ok/FAIL check and no mutation of its own (0 of core's 179): the turnover walk, the seismic columns, N^2, the scaling relations and the stitch's index map are covered by deck pins only (82.1%; neutrino_flux_output 0%);
  - open in bug_sweep: M139.
- **Extensibility 6.5 -> 8.0.** A history column and a profile column each have a recipe (io/README.md:212-226) that names stitched_model's five places, and a forgotten region is caught at run time (stitched_model.f90:948-957). Left: core/README does not point to them.
- **Docs 6.5 -> 8.0.** Headers with a Sources section for the stitch's fixes. Left: two statements that the code contradicts (finding 7: "nothing it computes feeds the physics"; turnover_timescale.f90:29-30: the deuterium limiter reads the current member, not the `_prev` copy).

Average 7.6.

### core's README and docs, whole package: 6.5 -> 8.5

README.md alone is at the standard, 9: the table's rows match the code (I checked every row of run/, model/ and diagnostics/, and the star-blind list against `STAR_BLIND_FILES`); section 2 gives the run and the step in order (the order of step 5, `compute_timestep` before the stitch and the observables, is step_driver.f90:180-221's); section 3 has the units, the aging flag, the ierr codes with the negative one, and the message routes; section 4 lists sixteen tests with what each checks and its mutations (179, the manifests' count; all 179 still apply at 1179afbf); section 5 has a coverage figure with its command and what does not run; section 6 has four recipes, each with its trace. Left for the README: finding 7's sentence; the message convention is stated more strictly than nine messages follow (finding 9); no pointer to io/README's column recipes or to a recipe for a kernel beyond section 3's one sentence; the coverage is of 2026-11-13, before round 2B. The in-code documentation in my half takes the whole to 8.5: finding 10.

**Recipe traced: "Add a stop condition"** (README.md:196-210), read against the tree for the existing `Teff_upper_limit`. The control occurs in 13 files. Six are hand-edited and the recipe names each: `io/controls_buffer.f90` (:152), `defaults/controls_registry.tsv` (:374), `io/legacy_controls.f90`'s namelist (:43), `core/run/stop_conditions.f90`'s four arrays (:194-210), `test_stop_conditions.py` (:69-72) and the re-pin of `io/test/expected_controls/`. The other seven are generated, each with a "GENERATED ... do not edit" head (`controls_sync.f90`, `inlist_decl.inc`, `control_group_def.inc`, `controls_buffer_public.inc`, `controls_dump.inc`, `control.defaults`, `controls_registry.md`), and the recipe's `make` step covers them (io/README.md:44-51). `make check` runs the contract checker the recipe names (Makefile:505-513). The recipe holds. I did not build it.

## Summary

- **What rounds 2 and 2B did here, confirmed in the code:** one message module and no stray write; `rezone` and `read_starting_model` as drivers over explicit-argument steps; the F77 shapes gone in all three directories; the star calibration's last card with a verdict and the C entry's path and unit rules, with test_run_contracts; test_model_kernels; the turnover refresh moved into the step; M34, M85/M96, M86, M97, M115, M126, M190, M191, M270, O4, O22, O6's part A and M31's narrow form. A sample of sixteen ROADMAP claims held; one is overstated ("the diagnostics feed no physics", finding 7).
- **What keeps model/ and diagnostics/ under 8:** the two long hosts that still share locals (`rescale_model`, `turnover_from_interior_new`); diagnostics/ with no test of its own and model/'s three drivers without one; the dependence of model/ on diagnostics/; and the defects below.
- **New defects, five:** a 2026 regression in the turnover time, recorded in the first bug sweep and never registered, now described by a comment (1); the run zeroes the turnover time the starting model has just computed, so each card's first step has no deuterium cap and a falling wind turnover time (2, F77); the core extension gives NaN for one added point, the default factor's case, and leaves the density profile unshifted (3, F77); the stitch swallows an envelope-integration failure after its ERROR line (4, 2026); the star calibration reports a crossing at the second model of every pass when the model starts above R* (5, F77).

## Previous report's items: status

All twelve items of the seventh report; items wholly in the other reviewer's half are marked "core_a" with what I saw of them.

1. **C-1 / M268, the star calibration's cards 51 and 52: DONE.** calibration.f90:287-298 sets up nothing past `num_cards` and writes a verdict; `end_star_calibration_card` (:556-568); test_run_contracts' star-calibration checks with 9 mutations.
2. **C-2 / M269, `henyey_iterate`'s undefined flag: DONE** (core_a; structure/henyey_iteration.f90:149-152 sets `model_diverged_flag` on every path).
3. **C-3 / M270, the starting model's warning: DONE.** starting_model.f90:321-344 prints the header's flag, worded as a warning.
4. **C-4 / M271, the C entry: DONE.** A path over 256 characters is refused (yrec_capi.f90:58-65); `close_run_units` closes the run log and the final model only (yrec_reset.f90:62-77). Left: M312 (open), and a NULL pointer from a C host reaches `strlen` (finding 6c).
5. **C-5, the messages module: DONE for routes, PARTLY for destinations.** core_messages exists and my three directories have no write of their own but file handling; the kernels of `KERNEL_FILES` are checked for I/O. The seventh's terminal-only events are still terminal-only in run/: the central-abundance stop (stop_conditions.f90:270) and the solar calibration's Newton step (calibration.f90:142-153) (finding 9). core_messages' `*_values` routes and `log_mixed_zone_bounds` write the run log's unit directly (M112's rest, recorded).
6. **C-6, the untested kernels: DONE** (core_a: test_henyey_elimination, test_envelope_integration). In my half the same gap is open for diagnostics/ and model/'s drivers (finding 12).
7. **C-7 / M272, `evolve_step`'s save: DONE** (core_a; no `save` in step_driver.f90).
8. **C-8, the star and its members as arguments; argument-less steps: PARTLY.** `read_starting_model(s, ierr)`, `rezone(s, ierr)`, `locate_shell_boundaries(s)` and `compute_timestep(s)` read the star themselves. Argument-less or host-sharing steps left in my half: model_rescaling 10, turnover_timescale 6, stitched_model 5, run_driver 4 (the seventh counted run 5, model 51, diagnostics 11).
9. **C-9, commented-out code: DONE.** model/retired/rezoning_theta_term.f90.txt with a header; `commented_code` is 0 in the style table. (Its header says the block sat in the "`if (s%physics%rotation_active)` branch"; the branch's condition is rezoning.f90:1312-1313's, with the overshoot case.)
10. **C-10 / M273, unnamed values: PARTLY.** Done: the calibration's loop ends from `num_calibration_runs` (calibration.f90:402-452), the species counts by name, `h1` passed to the shell kernel, the year and ln 10 (core_a). Left in my half: `do species_idx = 5,11` (model_rescaling.f90:183, :480, :488); the zone-bound dimensions (rezoning.f90:1205-1206); the atomic weights, now parameters but still two unsourced tables in core (model_rescaling.f90:33-34, step/shell_physics.f90:55) with copies in rotation and two tests (finding 11).
11. **C-11, the README: DONE** (the scores above).
12. **C-12 / M274, stale text: DONE** for every item in my half (README's six statements, calibration.f90's names, stop_conditions' header, yrec_capi's path, test_reentry.py's docstring; TARGET_ARCHITECTURE's Phase 2 status at :444-448). New stale text: finding 10.

## New findings (not in bug_sweep)

1. **core/diagnostics/turnover_timescale.f90:215-223, :345-353 (the comment at :346-348)**
   - **Issue:** when the convection zone's innermost point already lies more than one pressure scale height above the interpolated base (`val_test > 0` at `cz_base_index`), `time_scale_height_above_base` sets the turnover time to H_P/v_c there and leaves `spline_taucz_done` false; `time_bottom_cell_crossing` then runs and overwrites it with the bottom cell's width over v_c at its middle, another quantity. The code's comment says so ("This leaves spline_taucz_done false, so time_bottom_cell_crossing runs next and overwrites the result") instead of fixing it. The value feeds the wind's Rossby number and the deuterium cap of the next step, and the history column.
   - **Trigger:** a unique surface convection zone (`index_gap == 0`) whose base point is more than H_P above the interpolated base: a thin zone or a coarse mesh there.
   - **In the F77:** no. The bottom-cell code is inside the `ELSE` of that test, so the first branch kept its value:
     ```
     git show 6cd56739:src/tauintnew.f   # :142-147, :192, :230-232
              IF(VALTEST.GT.0.0)THEN
                 ...
                 TAUCZ = PSCA/CVEL
              ELSE
                 DO K = IMAX+1,M-1,1
                 ...
       110       CONTINUE
                 ...
                 TAUCZ = DELR/CVEL
       140       CONTINUE
              ENDIF
     ```
     The goto elimination (7b8c2e1d) made label 140 a flag that only the spline path sets. The first bug sweep found it (audit/bugsweep-2026-08-31/atm_net_math.md:65-82; docs/ROADMAP.md:457-460, "goto-elimination regression"); it has no bug_sweep number and was not fixed.
   - **Fix:** `spline_taucz_done = .true.` in the first branch (the flag's name then wants to be `timescale_done`). A check on a synthetic profile in a turnover test (finding 12), with this as its mutation.
   - **Kind:** B (every model that takes the branch; measure which decks do). **Size:** S.
2. **core/run/run_driver.f90:384-387 (with core/model/starting_model.f90:193-199; core/step/burn_step.f90:82-83, :699; core/structure/structure_equations.f90:138; rotation/evolve_rotation.f90:438-441; rotation/magnetic_braking.f90:582-583; net/reaction_rates.f90:201)**
   - **Issue:** `read_starting_model` builds the stitched model, computes the turnover time and copies it to `convective_turnover_timescale_prev`. `begin_kind_card` then sets `s%convective_turnover_timescale = 0` ("zero the convective turnover timescale (MHP 08/02)"), and nothing computes it again before the card's first step ends. During that step:
     - the deuterium cap reads the current member (burn_step and the shell kernel's config), finds 0, and is off (`turnover_timescale <= min_cap_turnover_time` returns);
     - the wind interpolates between the start value and the end value 0: `start + step_fraction*(0 - start)`, half the starting model's turnover time at the fraction 0.5.
     In every later step the two members are equal (`refresh_turnover_timescale` sets both from one model), so the first step of each card is the only one that differs, on every card, a card that continues the model in memory included.
   - **In the F77:** yes. The 2002 zero stayed in main after the 2017 initialization went into STARIN:
     ```
     git show 6cd56739:src/starin.f    # :904-908
     C INITIALIZE TAUCZ, PPHOT, AND FRACSTEP
           CALL GETTAU(...)
           TAUCZ0 = TAUCZ
     git show 6cd56739:src/main.f      # :290 CALL STARIN, then :412-414
     C MHP 08/02 STORE STARTING CZ PROPERTIES
              JCZ = JENV
              TAUCZ = 0.0D0
     ```
   - **Fix:** drop the assignment at run_driver.f90:387 (the starting model's value then stands in both members), and correct turnover_timescale.f90:29-30, which says the mid-step readers take the `_prev` copies.
   - **Kind:** B (the first step of each card: decks with a deuterium cap in a surface convection zone, and rotating decks whose wind law reads the turnover time). **Size:** S.
3. **core/model/starting_model.f90:605-606, :620-632, :639-648 (the comments at :620 and :635-638)**
   - **Issue, part a:** with one added point the extension divides 0 by 0. The shift copies old point 1 into slot 2 and leaves slot 1 as it was, so `actual_gradient = (logT(2)-logT(1))/(logP(2)-logP(1))` is 0/0; the new centre's log T is NaN. The comment says it ("for N = 1, slot 2 has already been overwritten and del is 0/0"). N is `int(log10(core_mass_reduction_factor)/chi_grid_scale(ichi_dm_max)) + 1`, so N = 1 for any factor below 10^0.08 = 1.20 at the default spacing, and the control's default is 1.0 (io/controls_buffer.f90:1203): `extend_core_inward = .true.` alone gives it. io checks neither control. M126's fix removed the case that reached it on a second card; this is the deck's own. All 213 shipped decks that extend the core set the factor to 10 (scratch lcore_decks.txt), so it is latent.
   - **Issue, part b:** the shift moves log m, log R, L, log P, log T, the flags, the composition and omega, and not log rho ("logRho is not moved"). After an extension by N points, `s%logRho(k)` for k > N is the density of the old point k, N points further out than the rest of point k, until the first structure solve stores a density (structure/structure_equations.f90:179). Meanwhile `density_estimate_offset` combines the old innermost point's P and T with old point N+1's density (:646-648); the starting model's stitched profile and turnover walk read the displaced density (diagnostics/stitched_model.f90:523); and when the card's first step ages the model, its predictor burn does too (step/burn_step.f90:682). N is 13 for the shipped decks' factor.
   - **In the F77:** yes, both:
     ```
     git show 6cd56739:src/starin.f    # :375, :390-400, :405-409
               MCORE = INT(DLOG10(FCORE)/HPTTOL(2))+1
               DO I=M,1, -1
                  HS(I+MCORE) = HS(I)  ...  HT(I+MCORE) = HT(I)   (no HD)
               RHOC = HD(1)
               DEL = (HT(2)-HT(1))/(HP(2)-HP(1))
               FACD = HP(MP1)-HD(MP1)-HT(MP1)
     ```
   - **Fix:** take the central density, the gradient and the luminosity from the two innermost old points before the shift (identical for N >= 2); refuse `core_mass_reduction_factor <= 1` with `extend_core_inward` in io's `check_controls`; shift log rho with the rest.
   - **Kind:** A for part a and the control check (only runs that now give NaN change); B for the density shift (measure the decks whose first card ages a model it reads with the extension). **Size:** S.
4. **core/diagnostics/stitched_model.f90:174-183 (with core/envelope/envelope_integration.f90:379-381, :416-419)**
   - **Issue:** `build_stitched_model` re-integrates the envelope and atmosphere at the converged model with the fixed output steps. When `atm_get` fails, the stitch keeps the interior alone and returns `ierr = 0`: `if (jerr /= 0) then; n_ie = n_ext; call fill_stitched_arrays(s, ierr); return`. `atm_get` has by then reported the failure ("ERROR: envelope integration ...: the envelope did not reach the fitting mass in N steps", or M132's microphysics error, or `envelope_termination`), so the run prints an ERROR line and goes on: the profile and pulse files of that model have no envelope or atmosphere, the seismic columns keep the last full build's values (:142-145), `photosphere_logP` keeps the previous model's, and the turnover walk runs on the interior. The procedure's header promises `ierr = 1` only for a missing column. README section 3: "a procedure that can fail returns ierr".
   - **Trigger:** any failure of that integration: the output step sizes (`atm_step_size`, `envelope_step_size`) needing more steps than the kernel allows, an EOS or opacity failure inside a step, an underflowing step.
   - **In the F77:** no (the stitched model is 2026; these lines are of 2026-08-25, b116224d).
   - **Fix:** return `jerr` (a positive code stops the run with exit 1, the negative one ends it as the solver's does), or, if an interior-only stitch is wanted for a model whose envelope cannot be output, report it as a warning that says so and document the fallback in the header.
   - **Kind:** A where no deck fails there (measure: the pins show whether any does); B for a run that now continues past the error. **Size:** S.
5. **core/run/calibration.f90:243-259 (with :507; core/run/stop_conditions.f90:170)**
   - **Issue:** `setup_star_calibration` sets `r_rsun_prev_model = 0`, and the check is skipped at each even card's first model, so nothing stores that model. At the second model the "previous radius" is 0: if the model's radius is above R*, the test `r_current > target .and. r_prev <= target` holds and the calibration reports "Just passed R*", interpolating log L and the age between this model and (R = 0, log L and age as the star's members were left, 0 on a first run), that is log L(R*) = log L * R*/R. The crossing path does not store the model either, so the next pass's second model does the same: every pass ends at its second model until card 50's "NOT converged". The header says the check works for a crossing "in either direction"; a track that starts above R* and contracts to it (a pre-main-sequence target) can never be calibrated. In a later pass the stored "previous model" is the previous pass's last model before its crossing, of another track.
   - **Trigger:** `calibrate_star_flag` with a starting model whose radius exceeds the target radius.
   - **In the F77:** yes:
     ```
     git show 6cd56739:src/setscal.f   # :58
           ALRI = 0
     git show 6cd56739:src/main.f      # :919-926
                 IF(MODELN.EQ.1) THEN
                    TEFFI = 10.0D0**TEFFL
                 ELSE
                    CALL CHKSCAL(BL, TEFFL, DAGE, NK)
     git show 6cd56739:src/chkscal.f   # the two tests after "Check if star has passed R*"
     ```
   - **Fix:** at the first model of an even card, store that model's L, age and radius as the previous model (a small `store_calibration_model(s)` that `check_stop_conditions` calls instead of skipping), so that each pass starts from its own track. A check in test_run_contracts: a pass that starts above R* and crosses it downward.
   - **Kind:** A (a calibration that starts below R* is unchanged: its second model is stored without a crossing either way). **Size:** S.
6. **Three reads or writes outside what was set**
   - **a. core/diagnostics/turnover_timescale.f90:196-203 (with :76-87).** The three-shell test reads `convective_flag(cz_base_index+1)` and `(+2)`. The caller fills its json-sized local arrays up to `num_points` only, so a base found at the walk's last or next-to-last point reads elements nothing set: `.false.` under gfortran's `-finit-local-zero` (the walk then ends as a radiative envelope), undefined under the Makefile's ifort flags, and outside the actual argument when `num_points` is json. **F77:** yes (`git show 6cd56739:src/tauintnew.f`, :90: `IF(.NOT.(LCZ(IMAX).AND.LCZ(IMAX+1).AND.LCZ(IMAX+2)))`). **Fix:** treat a point beyond `num_points` as radiative. **Kind:** A. **Size:** S.
   - **b. core/model/rezoning.f90:197-201.** A point is added outside the outer point (`s%nz = s%nz + 1`) with no test against json, where `place_points_along_chi` and the two jump steps test and stop. A model of exactly json points that is hotter than `rezone_env_logt1` at its surface writes element json + 1 of seven arrays. **F77:** yes (`git show 6cd56739:src/hpoint.f`, :104, `M = M + 1`). **Fix:** the same "needs more than json" error. **Kind:** A. **Size:** S.
   - **c. core/run/yrec_capi.f90:58-59.** `c_strlen` is called on both arguments before anything else; a C host that passes NULL for "use the default" crashes there. The header offers empty strings for that; pyyrec passes them. **Fix:** `c_associated` first, a NULL read as an empty path. **Kind:** A. **Size:** S.
7. **core/diagnostics/observables.f90:8-13, :111-135; README.md:66, :101**
   - **Issue:** the observables' header and the README say that nothing `compute_observables` computes feeds the physics ("The observables only fill the output's quantities"). `renormalize_luminosity_breakdown` rescales `star%luminosity_breakdown_lsun` in place (its own header: "This changes the model"), after `compute_timestep` has run for the step (step_driver.f90:181, :210). Its values are what the model file stores and what the next card's first `compute_timestep` reads (run_driver.f90:406; timestep_control.f90:79-85, :129: the hydrogen and helium burning limits and the sign of L_grav), where a step within a card reads the solver's values. The factor is 1 to within the solver's last correction, so the effect is small; the statement is still not true, and ROADMAP's "so that the diagnostics feed no physics" (docs/ROADMAP.md:4534-4535) rests on it.
   - **Fix:** move the renormalization into the step beside `refresh_turnover_timescale` and before `compute_timestep` (results change at the solver tolerance: measure it), or keep it and say in both places what it writes and who reads it.
   - **Kind:** A for the text; B for the move. **Size:** S.
8. **One formula in two places**
   - **core/diagnostics/observables.f90:88-108 and core/model/regridding.f90:79-100:** `shell_masses_from_log_mass` and `compute_shell_masses` are the same arithmetic in the same order (m from log m, dm from the midpoints, the last shell to M). The rezone calls regridding's, `read_starting_model` (starting_model.f90:1113-1116) and the neutrino table call observables', so model/ depends on diagnostics/ for a kernel model/ already has, and only one copy is tested (test_regrid's three shell-mass checks).
   - **core/model/rezoning.f90:672-756 and :767-849:** `add_points_for_x_jumps` and `add_points_for_z_jumps` differ in the species slot, the sign of the difference and the limit's index; 85 lines each.
   - **core/diagnostics/stitched_model.f90:162-168:** the spotted star's ambient Teff is written out again beside `ambient_log10_teff` (envelope/envelope_integration.f90:73-79, whose comment records the copy); the literals 0.25, 4.0 and 1.0 are exact in both precisions.
   - **Fix:** one shell-mass kernel (regridding's, pure); one `add_points_for_jumps(s, rz, slot, sign, limit, ierr)`; the stitch calls `ambient_log10_teff`. Each is byte-identical.
   - **Kind:** A. **Size:** S.
9. **Messages: where they go and how they read (README.md:114)**
   - **Terminal only, so the run log does not record them:** the solar calibration's Newton step, the lines "New X, Old X, DX" and their four companions (calibration.f90:142-153): a calibration's iterates are not in its log; the central-abundance stop's line (stop_conditions.f90:266-270; its reason does reach the log's summary); the Monte-Carlo samples' echo (monte_carlo.f90:105-108); the new interior points of an envelope refit (starting_model.f90:1092-1102).
   - **Against the README's convention** ("ERROR: " or "WARNING: ", the source, sentence case, "; run stopped" on an error that stops the run):
     - an error with no "ERROR: ": `' rezone: osplin: '//status` (rezoning.f90:963);
     - errors that stop the run without "; run stopped": the three spline failures of the rezone (rezoning.f90:650-651, :704-705, :799-800), the three of the turnover time (turnover_timescale.f90:387-389, :401-402, :462-463), the missing column (stitched_model.f90:950-953), and "; stopping" for the NaN guard (stop_conditions.f90:118);
     - a warning with no "WARNING: ", in capitals: "SPECIES ... MIX NOT MODIFIED IN EVOLVED MODEL" (starting_model.f90:438-440).
     - About 15 notes of model/ and run/ are still the F77's capitals ("OUTER POINTS DELETED", "NEW ENVELOPE MASS CALCULATED", "STARTING ... STOP DISABLED."); the convention is silent on notes.
   - **Fix:** `core_report` for the calibration's step and the abundance stop; the nine texts brought to the convention, in one message commit with its re-pin; a line in `check_style` that an "ERROR: " literal in core ends "; run stopped" would keep it.
   - **Kind:** B (run-log text only). **Size:** S.
10. **Stale or inexact text**
    - **Two contract blocks with a stray second pair** "star reads: none / star writes: none" under the real one: starting_model.f90:127-128 and rezoning.f90:115-116. The checker reads the first pair, so they pass.
    - **model_rescaling.f90:274:** "(locate_shell_boundaries, gn 10/05)"; the call is `find_shell_boundaries`. **:460-464:** the Z ramp's header uses the legacy names `rsclzc`, `rsclzm1`, `rsclzm2` for `card_z_ramp_*`.
    - **starting_model.f90:144:** "model_failed_flag (legacy LARGE)"; the member is `model_diverged_flag`. **:255, :307:** the file read "set[s] the timestep and star%log_total_mass"; it sets the timestep (the readers set the mass). **:239-241:** a test the comment calls always true is kept.
    - **turnover_timescale.f90:29-30:** "the mid-step consumers (the wind saturation, the deuterium limiter) read ... the _prev copies"; the deuterium limiter reads the current member (burn_step.f90:83, :699; finding 2).
    - **calibration.f90:4-7:** the module header lists four procedures; `end_star_calibration_card` is the fifth.
    - **README.md:66:** the `observables` row omits `shell_masses_from_log_mass`; the diagnostics rows do not point to io/README.md's "Add a history column" and "Add a profile column" (:212-226), the two recipes that end in diagnostics/.
    - **model/retired/rezoning_theta_term.f90.txt:** "the `if (s%physics%rotation_active)` branch" (it is the rotation-or-overshoot branch, rezoning.f90:1312-1313).
    - **stitched_model.f90:68:** `n_pulse_cols = 35` with 34 named columns (:74-82); say what column 35 is, or make it 34.
    - **Fix:** text only (the last item: check the writers first). **Kind:** A. **Size:** S.
11. **Unnamed values left (the seventh's finding 10, the part in this half)**
    - **Slots and sizes:** `do species_idx = 5,11` for `i_c12..i_o18` (model_rescaling.f90:183, :480, :488; starting_model uses the names); `radiative_zone_bounds(13,2)`, `convective_zone_bounds(12,2)` (rezoning.f90:1205-1206; starting_model.f90:135-136 uses `max_radiative_zones`, `max_convective_zones`); `do j = 1,4` for the rescale targets (calibration.f90:422, :526); the loop ends 8 and 10 for the neutrino slots (observables.f90:409, :416).
    - **Thresholds and conversions:** `1.0d-20` (model_rescaling.f90:140); `1.0D-8`, the surface-X change that rebuilds the opacity tables (rezoning.f90:1354); `0.250D0` and `8.0D0` of the added outer point (rezoning.f90:206, :214); `-1.0d-12` twice (starting_model.f90:782-783); `1.0d0` dex and `1.0d20` s in the turnover walk (turnover_timescale.f90:69, :188, :232); `8.64d4`, `1.0d-5`, `9999.0d0`, `1.0d6` (observables.f90:335, :439-442); the calibration's 15, 14 and 11 (run_driver.f90:300-302; monte_carlo.f90:271); the end-age test's 1 yr and 1e9 (stop_conditions.f90:249).
    - **Tables:** `mixture_atomic_weight` (model_rescaling.f90:33-34) and step/shell_physics.f90:55's four weights are two unsourced tables with different values for H and He (1.008 and 4.004 against 1.007825 and 4.002603); the capture-rate tables carry no source (observables.f90:45-48).
    - **Fix:** the names where they exist; module parameters with unit and meaning for the rest; a source line on each table. Byte-identical.
    - **Kind:** A. **Size:** S/M.
12. **Tests: diagnostics/ has none of its own, and model/'s drivers have none**
    - **Issue:** none of core's 179 mutations is in diagnostics/ (per file: regridding 19, calibration 9, shell_boundaries 9, yrec_capi 6, model_rescaling 6, monte_carlo 5, yrec_reset 2 in my half). The turnover walk is a procedure on arrays that takes the star for two solar constants and to store its result; findings 1 and 6a are in it, and the first has been known since August. `rezone`'s criteria and mesh, `read_starting_model`'s steps and `rescale_model` are checked by pins alone, and M315 passed them.
    - **Fix, each with its mutations:**
      - `test_turnover`: the walker star-free (the constants and the result as arguments; the file then qualifies for `STAR_BLIND_FILES`), on synthetic profiles: a unique zone timed at H_P/v_c, the innermost-point branch, two cells, a fully convective and a fully radiative model, a zone at the walk's last points;
      - `test_stitched`: the index map (no repeated radius at the two joins), the acoustic radius of a constant sound speed, a missing column's error, finding 4's error path;
      - in test_model_kernels: the mesh pieces that need no star once they take arrays (`accumulate_chi_spacing`, the jump insertion, `drop_crowded_points`, `sort_flag_points` with a triple duplicate), the core extension for N = 1 and N = 13, the horizontal-branch mass and core rescales (the masses the header's formulas give);
      - in test_run_contracts: `check_solar_calibration`'s Newton step against its 2x2 system, and the abundance stops' disarm and trigger.
    - **Kind:** A. **Size:** M.

## What would raise each score

- **run/ maintainability, 7.5 to 8.5:** findings 2 and 5 with their checks; finding 12's two run checks; M312; the Monte-Carlo files' handling moved to io (also organization).
- **model/ readability, 7.5 to 8.5:** `rescale_model`'s ten steps with explicit arguments and a state record, as `rezone` has (byte-identical); the step bodies re-indented and their F77 names replaced in one mechanical commit per file; finding 8's merged jump procedure; finding 11's names; model/'s eight lines over 132 columns wrapped (model_rescaling's five, rezoning.f90:1223, starting_model.f90:19-20).
- **model/ organization, 7.5 to 8.5:** one shell-mass kernel (finding 8), so that model/ calls diagnostics/ only for the starting model's stitch and turnover time, and those two calls made by the caller (`begin_kind_card`) after `read_starting_model`; the mesh criteria's arithmetic on arrays.
- **model/ maintainability, 7.0 to 8:** finding 3; finding 6b; finding 12's model checks; a pinned deck that rescales a horizontal-branch model; M202; O6's and M176's parts B when the owner schedules them.
- **diagnostics/ readability, 7.5 to 8.5:** `turnover_from_interior_new`'s steps with explicit arguments and a result record (and a name without "new"); `build_stitched_model`'s locals named; finding 11's values.
- **diagnostics/ organization, 7.5 to 8.5:** the turnover walk star-free and moved beside the step that owns it (step/ or a physics-state directory), returning its result; the luminosity renormalization moved or documented (finding 7); the shell-mass copy gone.
- **diagnostics/ maintainability, 7.0 to 8:** findings 1, 4 and 6a; `test_turnover` and `test_stitched` (finding 12); M139.
- **Docs, to 9:** finding 10; finding 7's sentences; finding 9's convention made true or stated as it is; a coverage figure taken after round 2B; the pointer to io/README's column recipes.

## Traces, counts and checks (supporting)

- **F77 shapes, seventh -> now** (lines with `double precision` / old operators / `endif`-`enddo`; scan.py): run 24/23/20 -> 0/0/0 (one `#endif` of the preprocessor); model 62/95/74 -> 0/0/0; diagnostics 31/26/20 -> 0/0/0. Whole-module `use` statements left are `use math_lib` only, which the style check exempts.
- **Message calls now** (`core_report` / `core_terminal` and `core_terminal_values` / `core_log` and `core_log_values`): run 6 / 14 / 3; model 24 / 2 / 13; diagnostics 4 / 0 / 3. Direct uses of the run log's unit: monte_carlo's two rewinds.
- **Long procedures** (scan2.py): `rescale_model` 421 lines, `run_yrec` 362, `turnover_from_interior_new` 354, `build_pulse_points` 185, `check_star_calibration` 172; nothing else in the three directories exceeds 110.
- **Controls read directly** (`s%control%`/`s%physics%` references): run 58, model 81, diagnostics 16; by decision for the star layer (README section 3).
- **Mutation manifests:** `mutate.py core --check` reports 179 of 179 applying at 1179afbf, the README's count and split.
- **ROADMAP claims sampled against the code, all found as stated:** M191 (monte_carlo.f90:233-236), M97 (:224), M34 (:343-344), M115 (:96-104), M126 (starting_model.f90:605-606), M190 (:372), M176's part A (:378), O9 (:381), M270 (:321-344), O6's part A (rezoning.f90:1023-1024), M31's narrow form (run_driver.f90:347-348), M268 (calibration.f90:289-298), M271 (yrec_capi.f90:58-65; yrec_reset.f90:62-77), O22 (stitched_model.f90:941-944), O4 (:470-474), the turnover refresh in the step (step/step_driver.f90:652-659).
- **Not checked:** anything that needs a build or a run (the tests' results, the pins, the coverage figure, the mutation results); the size of findings 1-3's effect on any deck.

# core_a -- eighth evaluation (2026-11-17, yrec-next 1179afbf, read-only)

Scope: the first half of core: `core/structure` (7 files, 2,402 lines), `core/step` (9 files, 3,583 lines) and `core/envelope` (5 files, 2,371 lines), plus `core/core_messages.f90` (146 lines), the parts of `core/README.md` (236 lines) that cover them, and their six test programs with manifests and pinned outputs (`test_structure`, `test_triangle_fit`, `test_mass_change`, `test_henyey_elimination`, `test_envelope_integration`, `test_timestep_limits`). Another reviewer has `core/model`, `run` and `diagnostics`.

- **Read whole:** all 21 production files of the three directories; core_messages.f90; core/README.md; the two `retired/*.f90.txt` headers; the six expected outputs; the headers of the six manifests; test_timestep_limits.py; mutate.py's header.
- **Sampled:** the test programs' sources (headers and check lists); atm_lib's lookup and switch functions and the Allard reader's range fields (for finding 3); convection's `homogenize_region`; numerics' bsstep statuses; io/controls_buffer for defaults.
- **Other inputs:** audit/module-evaluation-7-2026-11-10/core.md, FINDING_MAP.md and B_LIST_CORE.md (item 24); docs/ROADMAP.md "Round 2 (core)" and "Round 2B"; docs/GUIDELINES.md's four sections and "A package README"; docs/bug_sweep.md's table and the entries O1-O10, M89, M146, M303, M309-M318; `src/tools/check_boundaries.py` (STAR_BLIND_FILES, KERNEL_FILES, check_no_io), `style_counts.tsv`.
- **F77 evidence:** `git show 6cd56739:src/{mdot,ytime,htimer,envint,surfbc,getnewenv}.f`.
- **Counts:** scratch Python over the tree, in SP/review8_core_a/ (scan.py, unused.py, labels.py, mut_check.py and their .txt outputs). Nothing was built or run; no repository tool was run.
- **Checked against the rounds' claims (ROADMAP):** C-2, C-7, the messages route, `private`, the Henyey kernels and their test, C-8, `henyey_iterate`'s steps, the timestep kernels, the F77 shapes (style_counts.tsv: 0 in the three directories), retired/, the species counts, O5, O7, O8, O9, O19's core side, M132, M273's part B, item 24. All hold in the code except where the findings below say "in part". Every one of the 179 mutations' `old` strings occurs exactly once in its file at 1179afbf (mut_check.txt), so the manifests are in step with the code.
- Paths are under src/ unless they start at the repository root.

## Scores (1-10)

The seventh report scored core as a whole per criterion (readability 7, organization 7, maintainability 6.5, extensibility 6.5, docs 6.5) and the sub-directories by their average only (structure 7.7, step 7.0, envelope 6.2). The lines below give that average -> now, and the five criteria now.

### core/structure: 7.7 -> 8.1

- **Readability 8.0.** The headers derive the equations, the Henyey contract and the boundary rows with units; `henyey_iterate` reads as its eleven steps with explicit arguments; no F77 shapes are left. Left:
  - 21 argument-less internal steps that share their host's locals: `evaluate_shell_equations` 7 (shell_equations.f90:176-385), `compute_henyey_coefficients` 9 (structure_equations.f90:120-346), `compute_surface_boundary` 5 (surface_boundary.f90:137-284);
  - unnamed values (finding 9): the helium-flash factors, the tolerance slots `(5,1)`/`(5,2)`, channels by number, `1.0d-22`;
  - 13 unused names in `only:` lists (henyey_iteration.f90:61-62 has 11);
  - `s%max_residual` holds the residual maxima and then the correction maxima within one pass (henyey_iteration.f90:214-217, :338, :390).
- **Organization 8.0.** Two star-blind kernels (`shell_equations`, `triangle_fit`), two pure star-free kernels (`henyey_eliminate`, `henyey_solve`), every module private, messages through core_messages, contracts on every star procedure. Left:
  - the Newton scratch and about 25 output arrays are star members written on every iteration (`elim_coeff`, `elim_rhs`, `max_residual`, `surface_bc`; structure_equations.f90:48-55), against star rule 5 (TARGET_ARCHITECTURE's known debt);
  - `henyey_iterate` hands its divergence verdict back in `s%model_diverged_flag`, which `solve_structure` reads on the next line (henyey_iteration.f90:152, step_driver.f90:429), the pattern star rule 2 names;
  - finding 10: `energy_sources` reports from inside the shell kernel's call tree and is star-free but in no checker list; a nine-line star wrapper keeps henyey_elimination.f90 out of STAR_BLIND_FILES.
- **Maintainability 8.0.** 331 ok rows in three programs (test_structure 234, test_triangle_fit 83, test_henyey_elimination 14) and 40 mutations. Left:
  - open by decision: O1 (the Jacobian's gaps) and O3 (low-mass convergence);
  - finding 3 (the Allard fallback restarts the triangle at every model);
  - the temperature equation's radiative partials (`qt_dl`, `qt_dp`, `qt_dt`, shell_equations.f90:229-232) have ok rows but no mutation (finding 7);
  - `compute_surface_boundary`'s switch and triangle bookkeeping have no test of their own.
- **Extensibility 8.0.** The energy-source hook is complete, documented, traced and tested. Left:
  - TARGET_ARCHITECTURE §4.5a's two open items (the built-in terms as entries; a warning for a source above the fitting point);
  - a term that depends on R or L needs the terms record and the hand-written pivots extended (henyey_elimination.f90:38-41 says so; no recipe).
- **Docs 8.5.** The best headers of the half; the README rows match. Left: two stale header lines (henyey_elimination.f90:37, :87) and the message form (finding 6).

Average 8.1.

### core/step: 7.0 -> 7.7

- **Readability 7.5.** `timestep_limits`, `timestep_control` and `mass_change_kernel` are at the standard. Left:
  - 23 argument-less internal steps: `change_mass` 5 and `apply_mass_change` 7 (mass_change.f90:142-336, :435-614), `evolve_step` 7 sharing about 25 locals (step_driver.f90:233-636), `advance_composition` 4 (composition_step.f90:194-317);
  - `change_mass` spends about 200 lines on an entropy path whose results nothing reads (finding 5);
  - burn_step's thresholds are literals: `1.0d-32` five times, `1.0d-24` four times, `7.0d0`, `5.0d-7` (burn_step.f90:266-357, :571);
  - step_driver's card kind `2` and `6.6D0` (:199, :266) where `card_rescale` exists; five lines over 132 characters (:173, :443, :532, :553-554);
  - 17 unused names in `only:` lists (burn_step's `json` eight times);
  - three stale comments in step_driver (finding 9) and 76 lines of legacy and history notes, with no Sources section in the directory.
- **Organization 7.5.** Two star-blind kernel files and `convection_geometry`'s two pure functions. Left:
  - the Li/Be burn's algebra (`burn_li_be`, `average_li_be_over_cz`, 270 lines) and the composition step's mixing and closure take the star, so star rule 6's "something a reader could get wrong" is not a kernel there;
  - `homogenize_mixed_regions` (composition_step.f90:220-247) repeats convection's `homogenize_region`;
  - `burn_network_region` reports through run_messages, not core_messages (burn_step.f90:494, :565);
  - `shell_physics` keeps its own atomic-weight table, a copy of rotation's (shell_physics.f90:52-55).
- **Maintainability 7.5.** test_mass_change (78 rows, 14 mutations) and test_timestep_limits (33 rows, each with its mutation) are real conservation and formula tests; step/ is at 88.4% under the decks. Left:
  - findings 1 and 2 (latent, F77), finding 5;
  - open: M89, M146's rest;
  - no unit test of the Li/Be burn, the composition step's own arithmetic or `locate_surface_mixing_base` (convection_geometry was at 42.2% in the README's measure);
  - the shell branch's dX/dt limit (timestep_limits.f90:260-266) has no check.
- **Extensibility 8.0.** A timestep limit is a kernel, a call and one list entry, with a traced recipe. Left: the recipe leaves out four things a newcomer must touch (finding 8); no recipe for a stage of the step.
- **Docs 8.0.** Headers say what each stage does and why arguments are arguments. Left: finding 9's stale lines; mass_change.f90:8-9 still says the solar wind "is commented out in apply_mass_change" (it is in step/retired/).

Average 7.7.

### core/envelope: 6.2 -> 7.4

- **Readability 7.0.** The driver/kernel split and the events are clear. Left:
  - `integrate_envelope_atmosphere` is an 870-line host (envelope_integration_kernel.f90:221-1090) with 19 argument-less internal steps sharing about 70 locals, blocks not indented under their `if`/`do` (:403-435), and 44 lines that copy 22 EOS scalars into and out of `eos_res` (:461-511);
  - the same 25-line point store twice (:773-809, :951-999) and the Gamma1 formula three times (:628-631, :793-797, :972-976; a fourth in shell_equations.f90:380-384);
  - `point_physics`' names are the F77's (`qqdp`, `qqod`, `qdela`, `qrmu`; envelope_integrands.f90:23-34);
  - `atm_get` takes 24 arguments, 8 of them `intent(inout)` by habit (envelope_integration.f90:84-90, :109-118);
  - `rebuild_envelope` has 9 argument-less steps; `append_envelope_points` writes one interpolation twice (envelope_refit.f90:85-157).
- **Organization 7.5.** The kernel is star-blind, private, prints nothing and is in the no-I/O check. Left:
  - the kernel reads locals it never set and depends on the build's `-finit-local-zero` (finding 4);
  - its settings record is `intent(inout)` and the kernel changes `cfg%atm_choice` (:241, :380);
  - a dead sum, `taucz_env_accum`, with its flag and step (:773-774, :989-997);
  - the no-photosphere sentinel is `-1.0d99` in the kernel and `-1.0d98` in the driver (kernel :331, envelope_integration.f90:175);
  - the refit's and the rebuild's interpolation and shell-mass arithmetic take the star.
- **Maintainability 7.5.** test_envelope_integration: 35 checks (the integrands, every Bulirsch-Stoer step against RK4, the vertex, the thin envelope, the retry, six error contracts, the gray photosphere), 32 mutations; kernel 94.6%, integrands 98.0%. Left:
  - finding 4; open M309 and M310;
  - the tabulated surface (Kurucz, Allard, Castelli) and the fallback to gray have no check in the program (it asserts only "not switched to gray");
  - `append_envelope_points` and `rebuild_envelope` run under cm_accretion alone;
  - the `atm_retry` loop has no cap (:413-434): it ends only when the EOS fails on a vanishing density.
- **Extensibility 7.0.** `envint_step_config` lets a driver override the steps for one call. Left: no recipe for the envelope's likely changes: a quantity recorded at each point touches `point_physics`, both integrands, the kernel's two stores and its reversal, state's `envelope_points` and the stitched model; an event kind touches the kernel's list and `report_envint_events`.
- **Docs 8.0.** The kernel's header documents every event with its counts and values. Left: finding 6 (the messages' form, the rebuild's mislabelled table); envelope_integrands.f90:45 still speaks of `star%atm_choice_in_use` in a star-blind file.

Average 7.4.

### Package-level items, as they apply to this half

- **README (8.5):** all six parts of the standard are there, and every one of the 21 files and core_messages has a row that matches its public procedures. I traced "Add a timestep limit" against the code (finding 8). Left: finding 7 (the mutation sentence, the coverage figure's date, one stale coverage line); finding 6 (the message form it states); the Messages bullet is one paragraph of 1,653 characters (README.md:114); no recipe for a new test program or a new stage.
- **Messages (8.0):** one module, three routes, kernels silent, `check_no_io` over the star-blind files. Left: the stated form is not followed by seventeen messages of this half (finding 6); three procedures still write the run log's unit directly (M112's rest); core_messages has no test and ran at 37.7%.
- **Tests (8.0):** six programs, 477 ok rows, 119 mutations, all in step with the code. Left: finding 7; the untested drivers listed above.

## Summary

- **What the rounds did, and it holds:** the three directories have no F77 shapes left; the Henyey elimination and the timestep limits became kernels with tests whose every check has a mutation; the envelope kernel returns events instead of printing and has a 35-check test; the drivers read the star themselves; `henyey_iterate` is a list of steps; commented-out code is in retired/ with headers; the README has the standard's six parts with five traced recipes.
- **What keeps step and envelope under 8:** the long hosts with argument-less steps (72 in the half, against 90 at the seventh); the envelope kernel's shape and its dependence on a compiler flag; the dead accretion-entropy path; the Li/Be burn and the refit outside kernels and tests.
- **New possible bugs, each confirmed by reading the code twice and none by a run (read-only); whether a shipped deck meets the trigger of 2 or 3 is not checked:**
  - an accreting, disk-locked model older than the disk-locking age takes one step of negative length (F77);
  - a non-positive luminosity at the core edge gives a negative burning limit, which wins the minimum and becomes the next timestep by `abs` (F77);
  - an Allard deck whose gravity leaves the tables switches to gray and back at every model, with a new triangle each time (the switch is the F77's, the new triangle 2026's);
  - the envelope kernel reads unset locals, defined only by `-finit-local-zero` (2026: the F77 had `SAVE`).
- **Docs:** the README's "each check was confirmed to fail under a mutation" holds for two of this half's six programs; its coverage figure predates round 2B.

## Previous report's items: status

The seventh report's findings 1-12 (the plan's C-1 to C-12; bug_sweep M268-M274), as far as they touch this half:

1. **C-1, the star calibration's cards 51-52 (M268):** the other half (core/run).
2. **C-2, `henyey_iterate`'s undefined flag (M269): DONE.** The flag is set false at entry (henyey_iteration.f90:149-152). It is now the star member itself, not an argument (see structure's organization).
3. **C-3, the starting model's warning (M270):** the other half.
4. **C-4, the C entry (M271):** the other half.
5. **C-5, no messages module: DONE in the main, PARTLY in the detail.** `core_messages` exists and the three directories hold no I/O statement (style_counts.tsv). The kernel returns events (envelope_integration_kernel.f90:93-141) and `atm_get` reports them; the integrands' dead prints are gone; `check_no_io` covers STAR_BLIND_FILES. Left: `energy_sources` still reports for itself, now through core_report, when the shell kernel calls it (finding 10); the texts' form (finding 6).
6. **C-6, untested kernels: DONE.** test_envelope_integration and test_henyey_elimination, with manifests; `henyey_eliminate` and `henyey_solve` are pure and star-free (henyey_elimination.f90:104, :291).
7. **C-7, `evolve_step`'s vestigial save (M272): DONE.** No `save`, no reset block; the header says so (step_driver.f90:49-51).
8. **C-8, drivers take the star and its members: DONE for the members, PARTLY for the steps.** `compute_timestep(s)`, `rebuild_envelope(s, target, ierr)`, `compute_surface_boundary(s, ...)` and `compute_shell_physics(s, ierr)` read the star. Left: `atm_get`'s `log10_star_mass`, which every caller passes as `s%log_total_mass`. Of the argument-less steps, `henyey_iterate`'s and `compute_timestep`'s are done; 72 remain in this half (structure 21, step 23, envelope 28 with the rebuild's 9; scan.txt, which also counts `take_tightest_limit`, kept by design, and the module procedure `clear_energy_sources`), against 90 at the seventh.
9. **C-9, commented-out code: DONE.** structure/retired/rotation_stability_copies.f90.txt and step/retired/rotation_scaled_solar_wind.f90.txt, each with where from and how to revive; burn_step's path is corrected (burn_step.f90:93). mass_change's module header still says "commented out" (:8-9).
10. **C-10, unnamed values (M273 for part B): PARTLY.**
    - DONE: the year and its inverse from `seconds_per_year` (shell_equations.f90:148-149); ln 10 from physical_constants (timestep_limits.f90:28); the composition slots and species counts by name; `atm_eddington_gray` at the switches (surface_boundary.f90:139-147); the eps_grav floor's X named.
    - OPEN: the atomic-weight table (shell_physics.f90:52-55); `on_gray_fallback`'s literal 0 (surface_boundary.f90:33); the rest of finding 9.
11. **C-11, README: DONE.** The run's flow, the conventions with messages and the negative code, every test with its mutation count, the coverage section, five recipes with their traces, the energy_sources row. Left: findings 7 and 8.
12. **C-12, stale text (M274): DONE for this half's items.** The aging flag (README.md:112, timestep_control.f90:45-46), the triangle test's row, level 4 (README.md:98), the equations' home (README.md:111), `evolve_state`, the henyey-wave2 path (henyey_iteration.f90:237). New stale lines are in finding 9.

The seventh report's per-directory "Left":
- **structure:** `henyey_eliminate` takes the star and is untested: DONE. Finding 2: DONE. The commented-out block: DONE.
- **step:** timestep_control's F77 shape and 22-argument signature: DONE. The vestigial save: DONE. mass_change's dead accretion-entropy path: PARTLY: O5 and O7 are fixed, the dead path itself stays (finding 5).
- **envelope:** the untested kernel that prints: DONE. `if (.false.) print *`: DONE. Two modules public by default: DONE. M132: DONE (an integrand's failure inside a step returns its own positive code, envelope_integration_kernel.f90:906-909). O4's envelope side: fixed in round 2B (the other half's file).

## New findings (not in bug_sweep)

1. **core/step/mass_change.f90:477-483 (with core/step/step_driver.f90:283-297, :523; core/model/starting_model.f90:372)**
   - **Issue:** with accretion and disk locking, a step that would pass the disk-locking age is clipped to land on it: `s%dt = (disk_locking_age_gyr - disk_gate_age_gyr)*1.0d9*seconds_per_year`. Nothing tests that the model is younger than that age. `disk_gate_age_gyr` is the model's age (set at every model read and after every converged model).
   - **Trigger:** `use_mass_accretion`, `mass_accretion_rate > 0`, `disk_locking_active`, and a starting model older than `disk_locking_age_gyr` (default 0.0): for example a run restarted from a model saved after the disk-locking age, since the read re-arms `mass_change_active` (io/controls_sync.f90:915).
   - **What goes wrong:** the first step's `s%dt` is negative, by the model's age less the disk-locking age.
     - `delta_mass_cgs` (:419) is negative: the "accretion" removes rate x (age - disk age) from the convection zone, and with rotation takes the wind's angular momentum branch.
     - `evolve_step` then runs the composition step with the negative step (step_driver.f90:288), sets the age back to the disk-locking age (:297) and hands rotation the negative step (:523). The mass change switches off after it (:135).
   - **In the F77:** yes.
     ```
     git show 6cd56739:src/mdot.f   # :74-77
           IF(LDISK .AND. DMDT.GT.0.0D0)THEN
              TEST = SAGE + 1.0d-9*DELTS/CSECYR
              IF(TEST.GT.TDISK)THEN
                 DELTS = (TDISK-SAGE)*1.0D9*CSECYR
     ```
   - **Fix:** in `change_mass`, before the rate is used: when accretion is on with disk locking and `disk_gate_age_gyr >= disk_locking_age_gyr`, switch the mass change off and return without touching the timestep. A check in test_mass_change needs the clip as a small kernel (the limit and the clip: rate, masses, ages in, dt and `disk_exhausted` out).
   - **Kind:** B for the decks that trigger it (none pinned). **Size:** S.

2. **core/step/timestep_limits.f90:250-251, :320-321 (with :398-403; core/step/timestep_control.f90:92)**
   - **Issue:** the core branches of the hydrogen and helium limits divide by the luminosity at the core edge, `enclosed_mass(core_edge)/luminosity(core_edge)`, with no test of its sign. The shell branch of the helium limit has the test (`luminosity(shell_zone) > 0`, :325), and test_timestep_limits checks it there only ("helium: shell with L <= 0 -> unlimited").
   - **Trigger:** helium burning counted as on (triple-alpha luminosity above 1e-34 Lsun), central log rho below 5, and L <= 0 at the core edge, which is the centre when there is no convective core: a core whose neutrino losses exceed its heating, just below the flash's density threshold.
   - **What goes wrong:** the limit is negative. `tightest_timestep` takes `minval`, so the negative value beats every other limit and the growth limit, and `compute_timestep` makes it positive with `abs` (kept there for a stored model's negative step). The next timestep is |E_He dY m/L|, arbitrarily long for a small |L|. L = 0 gives +Inf, which loses the minimum (and traps in a devel build).
   - **In the F77:** yes.
     ```
     git show 6cd56739:src/ytime.f    # :68
     	    DELTSY = (5.85D17/CLSUN)*DELTSY*(HS1(JCORE)/HL(JCORE))
     git show 6cd56739:src/htimer.f   # :141, :178
              DELTS = MIN(DELTS,DELTSH,DELTSY,DELTSW,PDELT,EDELT)
           DELTS = ABS(DELTS)
     ```
   - **Fix:** `dt_unlimited` where the core edge's L is not positive, in both kernels, with a check and a mutation each; in `compute_timestep`, apply `abs` to the stored step only (the branch that skips the limits).
   - **Kind:** B only where it fires (byte-identical otherwise). **Size:** S.

3. **core/structure/surface_boundary.f90:31-34, :143-149 (with core/envelope/envelope_integration.f90:176-180; core/envelope/envelope_integration_kernel.f90:368-385; atm/atm_lib.f90:189, :246-251)**
   - **Issue:** two parts of the code decide when an Allard run is on the gray atmosphere, by different tests.
     - The kernel falls back to gray when a lookup fails: log Teff at or above the tables' maximum, or log g outside their range. `atm_get` then sets `atm_choice_in_use` to gray "for this and all future calls".
     - `switch_atmosphere_choice` switches back whenever the run is on the fallback and the model's log Teff is below the tables' last row. It does not look at log g, nor at the vertices' Teff. A switch back starts a new triangle (S12's fix).
   - **Trigger:** `atm_choice = 4` and a model whose surface gravity is outside the Allard tables' range while its Teff is inside (a low-gravity pre-main-sequence or giant model). Also a triangle vertex beyond the tables' maximum while the model is inside, which needs `tri_delta_teffl`/2 above the tables' last row spacing.
   - **What goes wrong:** at every model's surface check (level 2): the switch back to Allard with its note; a new triangle; the first vertex's lookup fails, with atm's message and "note: Allard table lookup out of range; switching to a gray atmosphere"; three envelope integrations instead of none. The run log gains three lines per model and the boundary condition is refitted from a fresh triangle every model. In the vertex case, vertex 1 is integrated with the Allard surface and vertices 2 and 3 with the gray one, so the planes are fitted through the jump between the two.
   - **In the F77:** the two tests, yes; the new triangle at each switch back, no (2026).
     ```
     git show 6cd56739:src/envint.f   # :144-150
              IF(LAlFail) then
                 KTTAU=0
                 LTTAU = .TRUE.
     git show 6cd56739:src/surfbc.f   # :75-77
              ELSE IF(KTTAU0.EQ.4.AND.TEFFL.LT.AL_TEFFLmax)THEN
                 KTTAU = KTTAU0
                 LTTAU = .FALSE.
     ```
   - **Fix:** one test for both directions, in atm: "the point is inside the tables" (Teff and log g), asked for the model and for the triangle's vertices. Leave the table when any vertex to be integrated is outside; return only when a fresh triangle about the model would be inside. The kernel's own fallback then never fires in a run and can become an error. A unit test of the switch on a synthetic table range.
   - **Kind:** B (Allard decks near the tables' edges). **Size:** M.

4. **core/envelope/envelope_integration_kernel.f90:461-482, :656-663, :685-686**
   - **Issue:** the kernel reads locals that nothing has set.
     - `find_atmosphere_start` copies 22 scalars into `eos_res` before `eos_get` (`pressure`, `log10_density`, `density`, `beta`, `beta_inverse`, `beta14`, `ion_fraction`, ...). On the first try of every integration only `temperature` has a value. `log10_density` is the density the EOS starts from.
     - `store_atmosphere_gradients` reads `prev_tau`, `prev_density` and `prev_opacity` at the first atmosphere point; its comment says the build's `-finit-local-zero` makes them zero, so that the first point's depth step is measured from tau = 0 with rho kappa = 0 there.
   - **What goes wrong:** with gfortran's flag the values are zero and every run is the same. The Makefile's ifort and ifx builds have no such flag (Makefile:201-210): there the EOS's starting density and the first `atmo_delta_depth`, which the stitched model's atmosphere radii are built from, are stack contents. The devel builds keep the flag, so no test sees it. bug_sweep has this class for kap and eos (M104, M107, M223), not for core.
   - **In the F77:** no in this form: ENVINT had a blanket `SAVE` (`git show 6cd56739:src/envint.f`, :108), so the values were the previous integration's; the SAVE went in 2026-08 (1c510cc5).
   - **Fix:** set `prev_tau`, `prev_density` and `prev_opacity` to zero at the start of `integrate_to_photosphere`; let `eos_res` keep its values between tries (it is the host's variable already) and delete the 44 lines of copies, or zero the scalars at entry. Both are byte-identical under the GNU builds; confirm with a build without `-finit-local-zero` and `-finit-real=snan`.
   - **Kind:** A. **Size:** S.

5. **core/step/mass_change.f90:118-128, :196-251, :258-302, :307-336, :490-541**
   - **Issue:** for accretion, `change_mass` computes the convection zone's mean thermal energy and entropy, looks up or integrates the surface the accreted matter joins (an Allard lookup, or since round 2B a gray atmosphere integration), calls the EOS there, and calls it again at Stahler's surface pressure and temperature. Then it sets the entropy to zero (:127), so every shift `adjust_cz_for_accreted_entropy` applies is zero, and `new_thermal_energy` is "not used further". The comments say so at each place; the second EOS call is kept "because the eos carries solver state".
   - **What goes wrong:** about 200 lines, an atmosphere integration and two EOS calls per accreting step whose results nothing reads, but whose failures end the run: an Allard lookup below the tables (:274), the gray integration (:279), either EOS call (:294, :334; a Saha non-convergence has been an error since M293). Stahler's temperature is far below any table for small rates: T_s^4 = (3/4) Mdot G M/(4 pi R^3 sigma) gives log T = 3.35 at 1e-9 Msun/yr on a solar-size star and 2.6 at 1e-12.
   - **In the F77:** yes (MASSLOSS's experimental entropy bookkeeping, "valid for Allard atmospheres only").
   - **Fix:** a measured commit that removes the path: `average_envelope_energies`' entropy sums, `find_accreted_surface_state`, `evaluate_stahler_surface`, `adjust_cz_for_accreted_entropy`, the members `delta_log_pressure` and `delta_log_temperature` with their use in structure_equations.f90:162-169, to step/retired/ with a header. Measure cm_accretion and the accretion decks (the EOS's carried state may move them at round-off).
   - **Kind:** B. **Size:** M.

6. **Messages' form: core/README.md:114 and core/core_messages.f90:17-19 against the texts**
   - **Issue:** both state the form: "an error begins 'ERROR: ', a warning 'WARNING: ', then where it comes from, in sentence case; an error the run stops on ends '; run stopped'". Item 24 reworded a list of messages; these in this half do not have the form:
     - no "ERROR: " on an error that stops the run: henyey_iteration.f90:195 and structure_equations.f90:260 (the two zero pivots), henyey_iteration.f90:463 (the corrections exceeded their limits), envelope_integration.f90:391-397 (a failed Bulirsch-Stoer step, an error when the integrands failed);
     - "ERROR: " without its origin or "; run stopped": surface_boundary.f90:193-196, :262-263; without "; run stopped": mass_change.f90:439-440, envelope_integration.f90:380-387, :417-418, envelope_rebuild.f90:100-101, :112-113, energy_sources.f90:161;
     - a warning in capitals without "WARNING: ": composition_step.f90:279 ("NO SURFACE CZ - DIFFUSION NOT MEANINGFUL", at every composition pass of a settling model without a surface convection zone);
     - notes in capitals to the terminal and the run log at every accreting model: envelope_rebuild.f90:150-151 and the table of :336-347, whose heading calls the L column "LOG L" and the log(m/M) column "LOG M" (the F77's heading, `git show 6cd56739:src/getnewenv.f`, :275-278);
     - a core message through run_messages: burn_step.f90:555-565.
   - **Fix:** one re-pin commit as item 24 was: the texts in the form; the rebuild's table to the run log with the solver diagnostics, its heading corrected; `burn_network_region` through core_report. Or say in both places which messages are exempt.
   - **Kind:** B (message text; cm_accretion's pin and the logs where the errors fire). **Size:** S.

7. **core/README.md:143, :149, :169, :175 (with the manifests of core/test/)**
   - **Issue:**
     - **:143** says "Each check of the nine programs with a manifest was confirmed to fail under a mutation of the code it covers". mutate.py counts a mutation caught when the one check it names prints FAIL. Of this half's six programs that covers every check in two (test_henyey_elimination 14 of 14, test_timestep_limits 33 of 33). In the others (labels.txt):
       - test_triangle_fit: 15 of 83 checks are named by a mutation;
       - test_structure: 8 of its 15 partials; `qt/lnP`, `qt/lnT` and `qt/L` have none, and in the manifest's Yale case they are `known` rows (convective shells), so the radiative branch's three partials (shell_equations.f90:229-232) are checked only without a mutation, in the TAMS case; ROADMAP's round-2 "Next" list says so, the README does not;
       - test_mass_change: 12 of 21 distinct checks (the others are the same check for other species);
       - test_envelope_integration: 32 of 35; the gray photosphere's "ierr = 0", "T = Teff" and "the star's atmosphere untouched" have none.
     - **:149** the coverage figure is of 2026-11-13 at e7551d23, before round 2B changed burn_step, mass_change, convection_geometry and the envelope driver.
     - **:175** "convection_geometry (42.2%): the Li/Be burn's envelope-overshoot branches, which no deck reaches (bug_sweep M303)": M303 is fixed, and cm_env_overshoot reaches them since round 2B.
     - **:169** lists `shell_physics` among "the kernels"; `compute_shell_physics` takes the star.
   - **Fix:** say what the manifests show (each mutation is caught by the check it names; give the checks covered per program); add mutations for the radiative `qt` partials, which needs mutate.py to take a second argument list per program, and for the three gray-photosphere checks; measure the coverage again at the current tree.
   - **Kind:** A. **Size:** S (the text), M (the mutations and the tool).

8. **core/README.md:188-194, "Add a timestep limit" (traced against core/step/timestep_limits.f90:31-33, core/step/timestep_control.f90:56-58, :143-147; test_timestep_limits.py; core/test/expected_test_timestep_limits.out)**
   - **Issue:** the recipe names the right three files and its trace is recorded, but four things a newcomer must also touch are not in it:
     - the kernel's name in `timestep_limits`' `public` list and in `compute_timestep`'s `use timestep_limits, only:` list (both modules are private with `only:` lists);
     - the pinned output core/test/expected_test_timestep_limits.out, which test_timestep_limits.py byte-compares: a new check line fails the pytest until it is regenerated;
     - `check_star_contracts --fix`, when the new setting or argument reads a star member (`timestep_settings_of`'s and `compute_timestep`'s contracts), and `make check`; the stop-condition recipe has this step;
     - the README's own counts (the test table's "33", the total of :143).
   - Two of the seventh report's five recipes are still missing: a new kernel with its test program (the Makefile's rule, `tests` and clean lists; devel_tests.py; the root `test_*.py`; the expected output; STAR_BLIND_FILES; README §4 and §5's command), and a new stage of the step. README.md:108 has one sentence on the first.
   - **Fix:** the four lines; the two recipes, each traced.
   - **Kind:** A. **Size:** S.

9. **Unnamed values, unused names and stale comments left in the three directories**
   - **Issue:**
     - **Thresholds and slots by number** (the bar: named parameters with unit and meaning):
       - burn_step.f90:266, :286, :344 (`1.0d-32`), :267, :279, :312, :356 (`1.0d-24`), :268 (`7.0d0`), :571 (`5.0d-7`);
       - henyey_iteration.f90:274-277 (`solver_tolerances(5,1)`, `(5,2)`), :412-413 (`5.0d-1`, `8.0d-1`, `5.0d-3`), :448 (`j=6,7` for the neutrino and gravothermal channels);
       - structure_equations.f90:282 (`do j = 1,6`), :292 (`1.0d-22`);
       - step_driver.f90:190 (`1.0D-3`), :199 and :266 (card kind `2`, where run/ and model/ use `card_rescale`), :266 (`6.6D0`);
       - envelope_refit.f90:93, :128 (`1.0D-14`); envelope_rebuild.f90:247 (`1.0D-10`); envelope_integration_kernel.f90:674 (`1.0d-6`), :274 (a parameter named `tiny`, the intrinsic's name); surface_boundary.f90:33 (`0` for `atm_eddington_gray`).
     - **Unused names in `only:` lists:** 30 (unused.txt): henyey_iteration.f90:61-62 (11 species slots), burn_step.f90 (`json` at :116, :218, :383, :447, :490, :634, :671, :721; `max_radiative_zones` at :721, :753), step_driver.f90:83-85 (`i_h1`, `i_h2`, `i_he4`), convection_geometry.f90:50, :87, mass_change.f90:31, composition_step.f90:85-86, shell_equations.f90:36, structure_equations.f90:62.
     - **Stale comments:**
       - step_driver.f90:136: the run stops with "MODEL FAILED TO CONVERGE"; the text is now "ERROR: starting model: the previous model failed to converge; run stopped" (model/starting_model.f90:148);
       - step_driver.f90:320: `timestep_limit_structure`, now `structure_change_limit`;
       - step_driver.f90:400: "the *_entropy_term arrays", now `dlog_temperature_dt` and its three siblings;
       - henyey_elimination.f90:37: "what compute_zone_energy returns"; the kernel calls `compute_nuclear_energy`; :87: "updates star%max_residual" in a star-free kernel;
       - mass_change.f90:8-9 (the solar wind "is commented out in apply_mass_change");
       - shell_physics.f90:16-18: a doubled banner that names the routine `shell_physics`;
       - core_messages.f90:36-39: the public list is grouped by round 2's two agents ("envelope/, structure/ and step/", "model/, run/ and diagnostics/"), not by what the routes do.
     - **Loose intents:** `log10_teff` is `intent(inout)` and only read in `compute_atmosphere_derivs` (envelope_integrands.f90:71); `atm_get`'s eight `intent(inout)` scalars (envelope_integration.f90:109-118), which force `compute_surface_boundary`'s two (surface_boundary.f90:94-95, its header says so).
   - **Fix:** names beside `store_li_be_rates`' `min_log_temperature`; the names dropped; the comments corrected; `intent(in)` down the chain. Byte-identical.
   - **Kind:** A. **Size:** S/M.

10. **core/structure/energy_sources.f90:43, :99-108, :160-162; core/structure/henyey_elimination.f90:75, :261-269; src/tools/check_boundaries.py:176-194**
    - **Issue:**
      - `evaluate_energy_sources` is called by the star-blind shell kernel (shell_equations.f90:354) and reports a failed source itself through core_report, so the kernel's call tree writes to the terminal and the run log; README.md:108 says the kernels print nothing. `energy_sources` uses no star, yet it is in neither STAR_BLIND_FILES nor KERNEL_FILES, so no check holds it.
      - `eliminate_star_shell_pair`, nine lines that pass the star's arrays, is the only reason henyey_elimination.f90 imports `star_info` and is not in STAR_BLIND_FILES; the README has to explain it as "star-free procedures in files that also hold a star driver".
    - **Fix:** `evaluate_energy_sources` returns the failing source's index beside `ierr`, and `compute_henyey_coefficients` reports it; the registry's two errors return `ierr` to the embedding program, which already tests it (or stay as the registry's own report, outside the kernel path); energy_sources.f90 and henyey_elimination.f90 into STAR_BLIND_FILES, with the wrapper moved into structure_equations.
    - **Kind:** A. **Size:** S.

11. **Formulas written more than once (GUIDELINES "Code style" 2)**
    - **Issue:**
      - Gamma1 from the EOS derivatives: shell_equations.f90:380-384 and envelope_integration_kernel.f90:628-631, :793-797, :972-976;
      - the envelope point's store, 25 lines twice: envelope_integration_kernel.f90:773-809 and :951-999;
      - the mass-weighted mean of a mixed region: composition_step.f90:220-247 against convection's `homogenize_region` (convection/convection_zones.f90:317);
      - the linear interpolation onto the fitting mass, twice with different operands: envelope_refit.f90:85-121 and :122-157;
      - the ambient Teff of a spotted surface: envelope_integration.f90:73-79 and a second copy in the stitched model, which the comment keeps apart as "not token-identical" (its single-precision literals 0.25, 4.0 and 1.0 are exact, so the two give the same bits);
      - the atomic weights: shell_physics.f90:52-55 and rotation's table (the seventh's item, open).
    - **Fix:** one `gamma1_from_eos` (eos_lib, or a small pure function in core); one `store_envelope_point(index)`; `homogenize_region` exported by convection_lib; one interpolation with its two end points as arguments; one ambient-Teff function; one sourced weights table in composition_layout. Each a byte-identical commit, checked on the decks.
    - **Kind:** A. **Size:** M.

## What would raise each score

- **structure, to 8.5-9:**
  - readability: the 21 steps of `evaluate_shell_equations`, `compute_henyey_coefficients` and `compute_surface_boundary` with explicit arguments (the shell kernel's can take a small working record); finding 9's names;
  - organization: finding 10; the divergence verdict as an argument of `henyey_iterate` again, stored by `solve_level`; the Newton scratch (`elim_coeff`, `elim_rhs`, `max_residual`, `surface_bc`) in a record `solve_structure` owns;
  - maintainability: finding 3 with a test of the switch; mutations for the radiative `qt` partials (finding 7); O1 when the owner schedules it;
  - extensibility: §4.5a's two items;
  - docs: the two stale header lines; finding 6.
- **step, to 8-8.5:**
  - readability: finding 5 (the dead path out); finding 9; `apply_mass_change`'s and `advance_composition`'s steps with arguments; `evolve_step`'s stages over a small step record (the level's flags, the zone outputs) instead of 25 shared locals;
  - organization: the Li/Be burn as a kernel (the rates and the base as arguments; `burn_li_be` gathers and stores), tested for the closed-form survival, the mixed average and the conservation of the accreted share; `homogenize_region` reused; one atomic-weight table; `burn_network_region` through core_messages;
  - maintainability: findings 1 and 2 with their checks and mutations; a check of the shell branch's dX/dt limit; M89 decided;
  - extensibility and docs: finding 8.
- **envelope, to 8-8.5:**
  - readability: the kernel's host split into the three parts its header names, as module procedures over one integration-state record (the 70 locals), with finding 11's single point store and Gamma1; `point_physics`' fields renamed; `atm_get`'s vertex group as one optional record and its inputs `intent(in)`;
  - organization: finding 4; `cfg` as `intent(in)` with the choice in use as a local; `taucz_env_accum` removed; one named no-photosphere sentinel; the refit's interpolation and the rebuild's shell masses as kernels;
  - maintainability: checks of a tabulated surface and of the fallback (a synthetic table through atm's test entry); a unit test of `append_envelope_points`' kernel (continuity at the fitting mass, the two branches, the degenerate interval); a cap on `atm_retry` with its event and error; M309;
  - extensibility: a recipe "record a quantity at each envelope point", traced;
  - docs: finding 6.
- **The package's items:**
  - README: finding 7's sentence and a coverage figure measured at the current tree; finding 8; the Messages bullet split;
  - messages: finding 6; a small test of core_messages (the unwritten lines, the trimmed lines, the three routes' destinations); M112's rest;
  - tests: finding 7's mutations.

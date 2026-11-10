# core -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

Scope: core/ is 38 source files (15,037 lines) in six sub-directories plus `main.f90`, five test programs (1,693 lines) with their pinned outputs, and one README (94 lines). This is core's first score.

- **Read whole:**
  - core/README.md; main.f90;
  - run/: run_driver, run_setup, yrec_reset, yrec_capi, calibration, monte_carlo and stop_conditions;
  - step/step_driver; structure/henyey_iteration and structure_equations;
  - step/timestep_control (lines 1-612) and step/mass_change (1-410);
  - the Python drivers test_structure.py, test_stop_conditions.py and test_monte_carlo.py.
- **Sampled by risk:**
  - henyey_elimination (its contract and the forward elimination), shell_equations (1-320), surface_boundary (100-305);
  - rezoning (the driver, the envelope window, the flagging, the json caps, the retired theta block);
  - starting_model (the driver, the model read, the controls check, rescale_and_refit_envelope);
  - model_rescaling (the driver and the X/Z rescale), composition_step (the driver), observables and energy_sources (their headers);
  - envelope_integration_kernel (every I/O statement and its error paths), envelope_integrands, the triangle_fit pivots and walk cap;
  - the headers and expected outputs of test_regrid, test_mass_change and test_triangle_fit; test_reentry.f90/.py, test_physics.py and test_error_exit.py (their headers and test list).
- **Skimmed through scans only:** stitched_model, turnover_timescale, burn_step, envelope_rebuild, envelope_refit, regridding, mass_change_kernel, shell_physics, convection_geometry, shell_boundaries and neutrino_flux_output.
- **Other inputs:**
  - docs/GUIDELINES.md's four sections; rotation/README.md (the standard); the sixth io report (the format);
  - docs/bug_sweep.md's table and the core items' entries;
  - docs/ROADMAP.md's core Stage A/B, design-debt (§8) and later core passages; docs/TARGET_ARCHITECTURE.md §4.3-4.5 and Phase 2.
- **Checks run:** tools/check_boundaries.py and tools/check_star_contracts.py, run read-only. Both pass ("Domain boundaries OK", "Star contracts OK (149 declared)"), and the tree is unchanged after them. Nothing was built or run.
- **Counts:** from scratch Python over the tree (SP/review7_core/: scan.py, subdir.py, literals.py, procs.py, paths.py and their .txt outputs).
- Paths are under src/ unless they start at the repository root (docs/, examples/, test_*.py).

## Scores (1-10)

First evaluation, so there is no history.

- **Readability 7.** The drivers read as their stages, and the headers are among the best in the tree: equations with units, the Henyey contract derived in full (henyey_elimination.f90:13-58), and star contracts on every driver. The code itself keeps much of its F77 shape. Left:
  - **F77 shapes:**
    - 459 `double precision` lines against 72 `real(dp)` (outside the tests);
    - 315 lines with `.lt.`/`.eq.`-style operators and 229 with `endif`/`enddo`;
    - 44 whole-module `use` statements (29 of physical_constants, 7 of eos_lib, 5 of kap_lib).
  - **Steps without arguments:** 157 internal subroutines share their host's locals and the star, against GUIDELINES "Code style" 1 (finding 8).
  - **Literals:** about 100 unnamed real literals in executable code, literal composition slots and `atm_choice` values (finding 10), and the bare card kinds the ROADMAP lists as open (docs/ROADMAP.md:3127-3128).
  - **Names explained instead of fixed:**
    - `dlogl_dlogr` and `dage_dlogr` are derivatives in linear R (calibration.f90:250-251);
    - `log_luminosity` holds L/Lsun (shell_physics.f90:28);
    - `corrections_too_large` is the caller's divergence flag.
  - **History in headers:** stop_conditions' header is a history narrative ("New (2026, core/ readability phase 2)", "2026 (bugsweep Batch 3)"), and 57 comment lines carry "2026" history. There are 2 Sources sections in 38 files.
- **Organization 7.** The layers follow TARGET_ARCHITECTURE §4.1, and core has six star-blind kernels in STAR_BLIND_FILES. There is no global star, and every star-taking procedure has its contract. Left:
  - core has no messages module (finding 5):
    - 69 direct terminal writes and M112's 70 run-log writes, against 5 calls to run_messages;
    - two star-blind kernels print.
  - drivers take the star and also its members as arguments, which every caller passes identically (finding 8; GUIDELINES star rule 3);
  - three modules are public by default: henyey_elimination, envelope_integration and envelope_integration_kernel;
  - process-level state:
    - evolve_step's vestigial blanket save (finding 7);
    - energy_sources' registry (acknowledged in its header).
  - diagnostics feed the physics: observables refreshes the turnover time that the deuterium cap and the mass change read (observables.f90:10-11).
- **Maintainability 6.5.** Four kernels have real ok/FAIL tests: the shell equations against finite differences (234 ok rows over three decks); the triangle (83 rows, including the ierr contracts); the mass change (56, conservation); the regridding (19). Two of them are recorded as failing under a mutation (docs/ROADMAP.md:1809, :1812). Left:
  - **Open bug_sweep items in core, about 22:**
    - the structure, solver and composition items: O1, O3-O10, O20, O22;
    - M31, M34, M85, M86, M96, M97, M112 (core's part), M115, M126, M132, M146, M176, M190 and M191.
  - **New latent bugs:** findings 1-4.
  - **Untested code:** the envelope integration kernel (1,052 lines) and the Henyey linear algebra (finding 6). Neither has a test of its own, and core has no coverage figure.
  - **Commented-out code** stays in the sources (finding 9).
- **Extensibility 6.5.** There are two documented extension points:
  - energy sources, a clean one (energy_sources.f90, `examples/other_energy`, test_structure's `extra` case, test_other_energy.py);
  - the stop table, whose header gives a one-line recipe (stop_conditions.f90:20-26).

  Left:
  - no README recipe at all (finding 11);
  - a new timestep limit edits a 22-argument signature at both of its callers (run_driver.f90:418-422, step_driver.f90:227-229);
  - a new mesh criterion goes into a 1,093-line rezone with 23 argument-less steps;
  - a new species meets the literal counts 11/15 (composition_step.f90:117-119, rezoning.f90:169-170) and `composition(15,json)` dummies (the ROADMAP's "literal 15").
- **Docs 6.5.** The in-code documentation is strong. README.md has §1's table (which matches every module but one) and a good §2 on the step. It has none of the standard's flow of the run, message conventions, coverage figure or recipes, and it has eight inaccurate or stale statements (findings 11 and 12). Left:
  - the README's missing sections;
  - the stale statements in README, TARGET_ARCHITECTURE and five headers (finding 12).

Average 6.7.

Per sub-directory (the five criteria averaged; they differ materially):

| Directory | Size | Score | Why |
|---|---|---|---|
| structure/ | 2,275 lines | **7.7** | The shell kernel and the triangle kernels are tested; energy sources are a clean extension point. Left: henyey_eliminate takes the star and is untested; finding 2; the commented-out block. |
| step/ | 3,607 | **7.0** | step_driver reads as its stages. Left: timestep_control's F77 shape (38 old-operator lines) and 22-argument signature; mass_change's dead accretion-entropy path (O5, O7); the vestigial save. |
| diagnostics/ | 2,124 | **7.0** | Observables come in themed procedures. Left: stitched_model at 994 lines; O4 and O22; the turnover time feeding back into the physics. |
| run/ | 1,929 | **6.8** | The protocols are stated at the call site. Left: findings 1 and 4; M190, M191, M97, M115, O20 and M34; terminal-only messages; history-heavy headers. |
| model/ | 2,978 | **6.6** | The regridding kernel is well tested. Left: rezone and read_starting_model at over 1,000 lines with 23 and 20 argument-less steps; 95 old-operator lines; literal slots; finding 3; O6, M85, M86, M96, M126 and M176. |
| envelope/ | 2,100 | **6.2** | Left: a star-blind 1,052-line kernel with no unit test that prints to the terminal and the run log; `if (.false.) print *`; two modules public by default; M132 and O4's envelope side. |

## Summary

- **What core does well.** These six points meet rotation's bar:
  - The run and the step are layered drivers with named stages, and the step's coupling iterations are written in the step driver (TARGET_ARCHITECTURE §4.3).
  - There is no global star. Every star-taking procedure declares what it reads and writes, and the checker passes.
  - The structure equations at one shell are a star-blind kernel with a finite-difference test of every partial. The approximations it marks `known` are documented.
  - The triangle, the mass change and the regridding are kernels with conservation or exactness tests, ok/FAIL and pinned.
  - The headers derive the Henyey system, the boundary rows and the timestep limits with their units.
  - The energy-source hook is a complete, documented and tested extension point.
- **What keeps core below 7 on three criteria:**
  - **Messages:** core writes its messages itself, 139 of them, while 5 go through run_messages. Some events reach the terminal but not the run log, and two star-blind kernels print.
  - **Tests:** the envelope integration, the largest kernel, has no test, nor has the Henyey elimination.
  - **README:** it lacks the flow of the run, the conventions for messages and for the negative ierr, a coverage figure, the tests outside the four kernels, and any recipe. It also carries stale claims about the negative timestep, the triangle test and the pinned decks.
  - **Style:** much of the F77 shape remains, and the steps are argument-less internal procedures.
- **Open defects already recorded:** about 22 bug_sweep items sit in core, all open. Among them:
  - the authors' decisions M176 (a card that reads the starting model takes the previous card's start-of-step deuterium rates), M190 (the disk-locking age across runs and cards; step_driver.f90:557) and M191 (monte_carlo.f90:210-211);
  - the later-card Z rescale's EOS metals and O16 floor (M85, M96);
  - the second core extension (M126);
  - the Jacobian, low-mass and stitched-atmosphere items (O1, O3, O4);
  - the rezone's deuterium gate (O6) and accretion mixing (O7).

  The ROADMAP's Stage B open list (its items 3-8: O6, O7, O1, the state items, O4, O3) is still open, item 6 in part (its dead stores were removed).
- **New possible bugs, all latent:**
  - a star calibration that reaches its 25th pair writes past the 50-slot run arrays and ends without a verdict (F77);
  - henyey_iterate leaves its `intent(out)` divergence flag undefined on its most common early return;
  - the starting model's mismatch warning prints a logical nothing sets (F77);
  - the C entry truncates paths over 256 characters silently.

## New findings (not in bug_sweep)

1. **core/run/calibration.f90:270-276, :288-292, :305-317 (with :468, :477, :489; core/run/stop_conditions.f90:171-178)**
   - **Issue:** `setup_star_calibration` lays out `max_runs` = 50 cards (25 rescale/evolve pairs), and `check_star_calibration` runs on every model of an even card. When the track crosses the target radius it writes the next pair's cards, `run_index+1` and `run_index+2`:
     - when L matches: `card_stops_active`, `card_max_age_yr`, `card_rescale_targets(i_target_x,:)` and `card_initial_x`;
     - otherwise: the next X guess in `card_rescale_targets(i_target_x,:)` and `card_initial_x`.
   - **Trigger:** a calibration that has not converged by card 50 and crosses R* there (24 passes without success, the first step and 23 secant steps, as an unreachable target or a too-tight tolerance gives). It writes cards 51 and 52 of arrays declared `(max_runs)` (state/star_info_lib.f90:282, :301, :315, :319), that is, into the star's next members (`card_initial_z(1:2)`, `card_z_ramp_core_z`, `card_stops_active(1:2)`, `card_central_deuterium_stop(1)`). A -fcheck=bounds build aborts there.
   - **What the user sees:**
     - On a match at card 50 the run prints "Have hit R* & L*, prepare final run to age" (:278-281), but no card is left. `step_leave_run_loop` ends the card loop, `end_of_card_calibration` never runs, and the run log's end-of-run reason is the card's "model budget exhausted" (stop_conditions.f90:290).
     - A calibration that runs out of cards ends the same way, with no verdict. The solar calibration instead says "solar calibration NOT converged after 14 cycles" (run_driver.f90:306-308).
   - **In the F77:** yes. CHKSCAL wrote `NK+1` and `NK+2` unconditionally, against SETSCAL's `NUMRUN = 50`:
     ```
     git show 6cd56739:src/chkscal.f   # :89-95, :109-111, :128-130
              LENDAG(NK+1) = .TRUE.
              ENDAGE(NK+1) = AGER*1.0D9
              LENDAG(NK+2) = .TRUE.
     git show 6cd56739:src/setscal.f   # :63
           NUMRUN = 50
     ```
   - **Fix:** in `check_star_calibration`, when `run_index + 2 > s%num_cards`, set a termination reason ("star calibration NOT converged after 24 passes", or "...: R* and L* reached on the last card, no final run") and leave the run list, as the solar calibration's verdict does. Use a test_error_exit-style deck with a target out of reach.
   - **Kind:** A (only runs that today write out of bounds change). **Size:** S.
2. **core/structure/henyey_iteration.f90:169, :191-192 (called at core/step/step_driver.f90:681-684)**
   - **Issue:** `corrections_too_large` is `intent(out)`, and its actual argument is `s%model_diverged_flag`. `henyey_iterate` returns at :191-192 when the surface boundary was rechecked without a new envelope and the previous level converged. It returns there without assigning the flag, so the flag is undefined on return.
     - That is level 2 of every step whose level 1 converged and whose triangle still holds the point, the commonest path through the solve. `solve_structure` then reads the flag (step_driver.f90:481).
     - gfortran today keeps the caller's `.false.` in practice. bug_sweep M135's entry records gfortran 15 at -O3 treating an unassigned `intent(out)` as clobbered (docs/bug_sweep.md:3458), and a compiler that poisons `intent(out)` dummies (NAG -C=undefined) would stop every run at its first step.
   - **In the F77:** no. CRRECT's `LARGE` had no intent and kept the caller's value on the same return:
     ```
     git show 6cd56739:src/crrect.f    # :113
           IF(LSBC .AND. .NOT.LFLP .AND. LFINI) RETURN
     ```
   - **Fix:** `corrections_too_large = .false.` at the top of `henyey_iterate`, beside `ierr = 0`.
   - **Kind:** A. **Size:** S.
3. **core/model/starting_model.f90:147-150, :365-376**
   - **Issue:** when an evolved model's mixing length or extended-composition flag differs from the deck's, the warning prints `lexcp0` for "MODEL USED", a local that nothing sets. Under the Makefile's -finit-local-zero it always prints F, whatever the model holds; under a compiler without that flag it prints garbage. The code comments the bug instead of fixing it ("a bug kept from starin"). The message also calls itself "ERROR IN SUBROUTINE STARIN" for a warning the run goes on after.
   - **In the F77:** yes. The 9/03 typo fix corrected the test, not the write:
     ```
     git show 6cd56739:src/starin.f    # :221-224
     C MHP 9/03 FIXED TYPO
              IF(.NOT.LMIXL .OR. LEXCOM0.NEQV.LEXCOM) THEN
                 WRITE(ISHORT,1040) CMIXL,CMIXL0,LEXCOM,LEXCP0
     ```
   - **Fix:** print `header%use_extended_composition` and drop `lexcp0`; word the message as a warning and send it through one report call.
   - **Kind:** B (the text of the run logs where the warning fires). **Size:** S.
4. **core/run/yrec_capi.f90:71-81 (with state/star_info_lib.f90:1239-1240; core/run/yrec_reset.f90:49-51; yrec_capi.f90:57-69)**
   - **Issue:** `c_string_value` copies at most 256 characters into the namelist override. A longer path from C or Python (pyyrec) is cut silently, and the run opens, or fails to open, the cut path with a message about a file the caller never named. M179 fixed the same class for control paths in io.
   - **Also:** both the C entry (after every run) and the re-entry prologue (before every run after the first) close every unit 7-99. A Fortran host that keeps a unit open there across `run_yrec` calls, such as a program like `examples/other_energy` calling it twice, has its file closed under it.
   - **In the F77:** no (2026).
   - **Fix:** take the C string's length (up to its NUL) into a deferred-length override and refuse a path the controls read cannot hold. Close only the units the run opened: run_files knows them, and `newunit` units are negative.
   - **Kind:** A. **Size:** S.
5. **Messages: core/structure/surface_boundary.f90:137-171; core/run/stop_conditions.f90:67-71, :265-267; core/run/calibration.f90:132-145; core/step/timestep_control.f90:122-125; core/run/monte_carlo.f90:81-84; core/envelope/envelope_integration_kernel.f90:320-323, :367-368, :512-513, :554, :810-811, :869, :896-899; core/structure/energy_sources.f90:95, :101, :153; core/envelope/envelope_integrands.f90:150-151, :297-299**
   - **Issue:** core has no messages module. It has 54 `write(*,...)` or `print` statements (two of them the dead `if (.false.)` lines below), 15 writes to `terminal_unit` and 70 to `run_log_unit` (the last are M112's), against 5 calls to run_messages (burn_step.f90:538, stitched_model.f90:944, neutrino_flux_output.f90:106-147). So:
     - **Events that reach the terminal but not the run log:** the atmosphere's switch to and from the gray boundary condition (a change of physics mid-run), the central-abundance stops, the "LROT set to TRUE, but OMEGA(1) = 0" configuration error, the calibration's Newton step, and the switch-off of the structure timestep limits. A run's log does not say they happened.
     - **Star-blind kernels print**, against GUIDELINES "Code style" 7:
       - the envelope integration kernel writes its notes and failures to the terminal and the run log itself;
       - energy_sources writes its registry errors with `write(*,...)`.
     - **Dead statements:** envelope_integrands keeps `if (.false.) print *, ...` to silence unused-dummy warnings.
     - check_boundaries' `check_no_io` covers numerics and math only, so nothing stops this.
   - **Fix:**
     - one core messages module (or run_messages' report_error/warning/note), through which every core message goes, to the terminal and the run log;
     - the envelope kernel and energy_sources return codes, and their drivers report;
     - extend `check_no_io` to STAR_BLIND_FILES;
     - drop the `if (.false.)` lines (unused dummies of a shared interface can be named in the interface's comment instead).
   - **Kind:** A where the text and its destinations stay; B for the run logs that gain the terminal-only lines. **Size:** M.
6. **Untested kernels: tools/check_boundaries.py:167, :171; core/envelope/envelope_integration_kernel.f90; core/envelope/envelope_integrands.f90; core/structure/henyey_elimination.f90:88, :240; README.md:89**
   - **Issue:** the envelope integration kernel (1,052 lines: atmosphere and envelope integrations, the Allard fallback, the stored envelope state) and its integrands are in STAR_BLIND_FILES, but no test program uses either; only the decks run them. The Henyey linear algebra (`henyey_eliminate`, `henyey_solve`) is tested the same way.
     - Star rule 6 asks a kernel for a unit test, and README:89 says "The kernels have unit tests of their own", which holds for four of six.
     - A wrong pivot in the elimination would show only as a slower or failed convergence. A wrong integrand would show as a shifted track, which a re-pin can absorb.
   - **In the F77:** no tests at all.
   - **Fix:**
     - `test_envelope`: a gray atmosphere against the Eddington T(τ) it integrates, a radiative envelope of constant opacity against its analytic P(T), an ideal-gas adiabatic one against P ∝ T^(5/2), and the ierr contracts (maxstp, bsstep bad input, envelope_termination);
     - `test_henyey_solve`: a random block-bidiagonal system of the documented form solved by `henyey_eliminate`/`henyey_solve` against a dense solve, which needs `henyey_eliminate` star-free first;
     - a mutation each.
   - **Kind:** A. **Size:** M.
7. **core/step/step_driver.f90:46-51, :138-168 (with state/star_info_lib.f90:1241; core/run/yrec_reset.f90:48)**
   - **Issue:** `evolve_step` keeps a blanket `save`, a block that resets 19 of its locals, and a process-wide flag that arms it. Its header says that "several [locals] carry state from one step to the next (the convergence bookkeeping, the surface convection zone's base zone of the previous step)... alongside evolve_state". None does, and nothing named `evolve_state` exists. Every local is assigned in the call before it is read:
     - `converged` and `iterations_done` in `solve_structure`;
     - `envelope_cz_zone_prev` in `rezone_or_snapshot` (:396), before `burn_light_elements` and rotation's sub-steps read it under the extended composition, the only case they do;
     - `target_envelope_mass` and `new_atmosphere_fit_needed` by `change_mass`;
     - the zone counts as `intent(out)` of `find_convection_zones`;
     - the rest in `solve_level`.
     - M169's bisection found likewise that "`evolve_step`'s saved locals did not matter" (docs/bug_sweep.md:3920).
     - ROADMAP §8 records the save and its reset list as debt, to be replaced by a derived-type state. The derived type is not needed: there is no state.
   - **Fix:** remove the `save`, the reset block, `evolve_step_reset_pending` and its arming in `yrec_reset`. This is byte-identical: every local is defined before use. Check it with the poison build of ROADMAP's method.
   - **Kind:** A. **Size:** S.
8. **Drivers take the star and its members (core/run/run_driver.f90:360-364, :413-415, :418-422; core/step/step_driver.f90:203-205, :212-214, :227-229, :292-296; core/structure/henyey_iteration.f90:109-118; calls of compute_shell_physics at core/model/starting_model.f90:232, core/model/rezoning.f90:1018 and rotation/evolve_rotation.f90:158)**
   - **Issue:** 22 calls pass `s` together with three or more of its own members:
     - `rebuild_envelope` gets `s` and 22 star members;
     - `compute_timestep` gets `s` and 19;
     - `compute_surface_boundary` gets `s` and 12;
     - `compute_shell_physics` gets `s` and 12, with its three callers passing the same members;
     - `locate_shell_boundaries` gets `s` and 10, and `read_starting_model` `s` and 9, each from two identical calls.

     GUIDELINES star rule 3 allows model arrays as arguments when callers pass different models, and asks the header to say why. Here every caller passes the same ones, and no header says why, so a reader cannot tell which arguments vary.
   - **The steps without arguments:** the long drivers' named steps are internal procedures with no arguments: 157 of them, sharing the host's locals and the star (rezone's 23, read_starting_model's 20, henyey_iterate's 10). GUIDELINES "Code style" 1 asks for steps "with explicit arguments, so the driver reads as the list of steps". Rotation has 30 such steps, net none.
   - **Fix:**
     - read from `s` what every caller passes from `s`;
     - make the timestep limits star-free kernels (arrays and a settings record in, a limit out), which also makes them testable;
     - give the steps of `henyey_iterate`, `compute_timestep` and `rezone` their inputs and outputs as arguments, one driver per byte-identical commit.
   - **Kind:** A. **Size:** M/L.
9. **Commented-out code: core/structure/structure_equations.f90:298-313; core/step/mass_change.f90:364, :391-404; core/model/rezoning.f90:154-156, :1044-1055**
   - **Issue:** three blocks and two declarations of superseded code stay in the sources: rotation's stability copies, the rotation-scaled solar wind and the theta term, all "commented out with the advection-diffusion transport (2026-10-12)". Rotation keeps such code in rotation/retired/*.f90.txt with a header saying where it came from and how to revive it. core has no retired/ directory, and mass_change's revival note asks the reader to "uncomment this block, its locals and its star members". burn_step.f90:90-91 points to "rotation/burn_li_be_midpoint.f90", which is now rotation/retired/burn_li_be_midpoint.f90.txt.
   - **Fix:** move the blocks to core/retired/*.f90.txt with rotation/retired's header form, leaving a one-line pointer; correct burn_step's path.
   - **Kind:** A. **Size:** S.
10. **Unnamed values: core/structure/surface_boundary.f90:32, :138-170; core/envelope/envelope_integration.f90:205; core/envelope/envelope_integration_kernel.f90:306, :315; core/step/mass_change.f90:258, :264-265; core/model/model_rescaling.f90:89, :121-122, :155-157, :195, :404, :426-436; core/model/shell_boundaries.f90:87-102; core/step/composition_step.f90:117-119; core/model/rezoning.f90:115, :169-170; core/run/calibration.f90:373, :392; core/structure/shell_equations.f90:146-147, :285; core/step/timestep_control.f90:331-332; core/step/shell_physics.f90:60-61; core/model/starting_model.f90:133-136**
    - **Issue:** GUIDELINES "Code style" 3-4 and the bar's "thresholds are named parameters, constants named and sourced":
      - **`atm_choice`:** compared and set as 0, 3, 4 and 5 at 13 places, with no named constants in atm_lib. The Kurucz switch temperature log Teff = 3.95 appears twice as a literal.
      - **Composition slots by number:** `composition(1,…)` and `composition(3,…)` in about 21 lines (model_rescaling, shell_boundaries, `s%xa(1,nz)`/`s%xa(3,nz)` in mass_change), where `i_h1`/`i_metals` exist. The species counts are written 11 and 15 where `n_species_basic`/`n_species_extended` exist, and the zone-bound dimensions 13 and 12 where `max_radiative_zones`/`max_convective_zones` exist.
      - **Solar calibration loops:** they end at a literal 46 and 47 while the run list is `max_runs - 2`, so the layout and `num_calibration_runs` disagree as soon as `max_runs` changes. check_max_runs (M192) does not see loop bounds.
      - **A second year in the shell kernel:**
        - `one_year_sec = 3.1558d7`, equal to physical_constants' `seconds_per_year`;
        - `one_year_sec_inv = 3.1688d-8`, which is not its inverse (3.16877e-8, a relative 7e-6);
        - the X > 0.01 threshold of the eps_grav floor is unnamed.

        All three are inherited (`git show 6cd56739:src/coefft.f`, :97-98 and :264: `TIME = 3.1558D7`, `TIMEI = 3.1688D-8`).
      - **ln 10 as 2.3026d0** in the structure timestep limit (F77 PTIME).
      - **Atomic weights** in two local `data` tables: shell_physics' is a copy of rotation/seculr/transport_checks.f90:303, and starting_model has its own mixture table.
    - **Fix:**
      - named `atm_*` constants exported by atm_lib;
      - the slot and size names at each site;
      - the loop bounds from `num_calibration_runs`;
      - `seconds_per_year` and its exact inverse, and `ln10`;
      - one atomic-weight parameter table in composition_layout or physical_constants, sourced.
    - **Kind:** A for the names, the loop bounds and the weights (byte-identical); B for the year's inverse and ln 10. **Size:** M.
11. **README.md (the standard's sections), :85-94**
    - **Issue:** the standard asks a README for six things, and core's falls short on four:
      - **The flow, what calls what in order:** §2 gives it for the step only. The run has none: read, open, setup, Monte Carlo runs, calibration set-up, cards, `begin_kind_card`'s eleven steps, the verdict, summaries. It is in run_driver's comments alone.
      - **Conventions:** §3 has no rule for messages (where they go), nothing on the negative ierr that main.f90:20-23 and yrec_capi.f90:19-23 define (`envelope_termination`: the run ends and the process exits 0).
      - **What each test checks:** §4 lists the four kernel tests and test_reentry. It does not list:
        - test_stop_conditions (the four structure-limit stops and the seismic columns);
        - test_error_exit (15 tests, 22 cases, several of them core's: the malformed starting model, the triangle size, the Bulirsch-Stoer limits);
        - test_monte_carlo (M169's two identical runs);
        - test_physics (a solar model within loose bounds);
        - test_other_energy;
        - test_structure's `extra` case.

        No test's mutation check is recorded (the ROADMAP records two). There is no coverage figure; tools/coverage.py takes `--dir core`.
      - **Recipes:** there is none. The extension points a newcomer meets are each documented only in a header, an example or another domain's README:
        - add an energy source;
        - add a stop condition;
        - add a timestep limit;
        - add a stage to the step;
        - add a kernel (STAR_BLIND_FILES, a test).

      §1's table has no row for structure/energy_sources.f90, the extension point TARGET_ARCHITECTURE §4.5a and net/README.md:198 send readers to.
    - **Fix:** a run section after §2; message and ierr conventions in §3; §4 with every test and a coverage figure; a §5 of recipes, each traced (at least "Add an energy source" and "Add a stop condition" end to end); the energy_sources row.
    - **Kind:** A. **Size:** M.
12. **Stale or inexact text:**
    - **README.md:82:** "A negative `s%dt` means 'solve, do not age' (a rescale, a relaxation); it is made positive after convergence."
      - The aging flag replaced the negative step on 2026-10-28 (88b32a90, docs/ROADMAP.md:2069; step_driver.f90:417-423). `s%dt` is now negative only if a model file stores a negative timestep.
      - timestep_control.f90:39-40 ("a model being relaxed") is likewise out of date.
    - **README.md:92**, test_triangle_fit:
      - It says degenerate fitting points "give non-finite coefficients, unchecked". Since d2ff6886 the fit refuses them with ierr (triangle_fit.f90:39-43, :196-221).
      - It says the walk "would never end" for a zero-area or mis-oriented triangle. The walk now gives up after `max_reflections` (:38, :171-172).
      - The test checks both (expected_test_triangle_fit.out:25-32).
    - **README.md:26:** `prepare_run` "(that, then the Monte-Carlo run range)". It reads the Monte-Carlo samples. The range has been `first_mc_run`/`last_mc_run`, functions of the controls, since 0ec5ed1d. The row also leaves out the run plan's envelope Z (run_setup.f90:39, :57-67).
    - **README.md:87:** "No pinned deck rescales a model". The pinned solar decks' card 1 is `KINDRN(1) = 3` with `RSCLX`, `RSCLZ` and `RSCLM` (examples/run_standard_solar_model/Test_solar_noGS_norot.nml1:22-27), and every card with `card_kind /= 1` rescales (starting_model.f90:655). What no pinned deck covers is the horizontal-branch rescale.
    - **README.md:93:** test_mass_change's "accretion conserves every basic species" holds on the test's zone of uniform composition (test_mass_change.f90:11-12), the one case O7's mixing conserves.
    - **README.md:72:** "with rotation on, level 4 inside the structure⇄rotation iteration". Level 4, the final convergence verdict and the full mixing of the convection zones run for every model (step_driver.f90:513-554); without rotation the iteration has one pass.
    - **README.md:81 and docs/TARGET_ARCHITECTURE.md:235:** both send the reader to structure_equations for the equations and `qt_dl`. They are in shell_equations' header, as structure_equations.f90:31-32 says.
    - **docs/TARGET_ARCHITECTURE.md:444:** Phase 2's status lists "the coupling iterations in the step driver, the `aging` flag" as open. Both are done (§4.3 at :213-215; 88b32a90).
    - **step_driver.f90:51:** "evolve_state" (nothing has that name; finding 7).
    - **henyey_iteration.f90:247:** cites "reports/henyey-wave2.md", which is audit/readability-sweep-2026-09-03/reports/henyey-wave2.md, as observables.f90:128 has it.
    - **calibration.f90:**
      - :202-204: "despite the log_r_* names of the model values". The model values are `r_rsun_*`.
      - :250-251: names the comment admits are wrong (finding 10's class; rename them `dlogl_dr`, `dage_dr`).
    - **stop_conditions.f90:129-131:** "iteration 1 only primes previous-model state". Nothing is called at iteration 1 (:172).
    - **yrec_capi.f90:6:** "the ctypes binding under /Applications/YREC/pyyrec", a path on one machine. It is pyyrec/ at the repository root.
    - **test_reentry.py:13-17:** the m0030 deck "ends via the BSSTEP numerics stop, which kills the process mid-call". bsstep returns a status (numerics has no stop, `check_no_io`), and the deck now stops at 50 models (docs/ROADMAP.md, core Stage A).
    - **Fix:** text only. **Kind:** A. **Size:** S.

## What would raise each score

- **Readability, to 8-8.5:**
  - the F77 shapes converted, a byte-identical commit per directory: `real(dp)` with the `iso_fortran_env` import, `<`/`==`, `end if`/`end do`, `only:` lists on physical_constants, eos_lib and kap_lib;
  - finding 10's names, and the card kinds from io/run_cards (the ROADMAP's open item);
  - finding 8's explicit-argument steps, at least for `henyey_iterate`, `compute_timestep` and `rezone`;
  - names that hold what they say (finding 12's calibration and shell_physics cases);
  - stop_conditions' and the other history headers as Sources sections.
- **Organization, to 8-8.5:**
  - finding 5: one route for messages, the kernels silent, and `check_no_io` extended to STAR_BLIND_FILES;
  - finding 8's star-member arguments removed;
  - `private` in henyey_elimination, envelope_integration and envelope_integration_kernel;
  - finding 7's vestigial save removed;
  - the turnover time's refresh moved out of observables into the step, so diagnostics only read.
- **Maintainability, to 7.5-8:**
  - findings 1-4;
  - finding 6's envelope and Henyey tests, with mutations;
  - a coverage figure for core under the unit tests and the pinned decks;
  - finding 9's retired/;
  - the open items the authors can settle without a decision, each a measured commit: O4's out-of-range read (Stage B item 7), O6 and O7 (items 3-4), M126;
  - M176, M190 and M191 decided.
- **Extensibility, to 7.5-8:**
  - finding 11's recipes, traced;
  - the timestep limits as kernels behind a settings record, so that a new limit is one kernel, one line in `take_tightest_limit` and a test;
  - the species counts and composition dimensions named (finding 10; the ROADMAP's literal 15);
  - TARGET_ARCHITECTURE §4.5a's open items for the energy sources: the built-in terms as entries, and the warning for a source above the fitting point.
- **Docs, to 8-8.5:**
  - finding 11's sections (the run's flow, the message and ierr conventions, every test with what it checks and its mutation, a coverage figure, the recipes, the energy_sources row);
  - finding 12's corrections;
  - the TARGET_ARCHITECTURE Phase 2 status brought up to date.

## Traces, counts and README checks (supporting)

- **Message routes:** per sub-directory, the terminal writes / the terminal_unit writes / the run_log_unit writes / the run_messages calls:

  | Directory | Terminal writes | `terminal_unit` | `run_log_unit` | run_messages |
  |---|---|---|---|---|
  | run | 12 | 6 | 10 | 0 |
  | step | 6 | 0 | 2 | 1 |
  | structure | 11 | 3 | 10 | 0 |
  | envelope | 10 | 3 | 9 | 0 |
  | model | 15 | 3 | 36 | 0 |
  | diagnostics | 0 | 0 | 3 | 4 |

- **F77 shapes, per sub-directory** (lines with `double precision` / old operators / `endif`-`enddo` / argument-less subroutines / json-sized dummies):

  | Directory | `double precision` | Old operators | `endif`/`enddo` | Argument-less subroutines | json-sized dummies |
  |---|---|---|---|---|---|
  | run | 24 | 23 | 20 | 5 | 0 |
  | step | 155 | 96 | 62 | 31 | 41 |
  | structure | 69 | 38 | 26 | 32 | 2 |
  | envelope | 118 | 37 | 27 | 27 | 8 |
  | model | 62 | 95 | 74 | 51 | 6 |
  | diagnostics | 31 | 26 | 20 | 11 | 5 |

  For comparison, rotation (the standard; 7,607 lines outside its tests): 344 `double precision` lines, 112 old-operator lines, 60 `endif`/`enddo`, 30 argument-less subroutines and 119 json-sized dummy lines. Core's density of `double precision` is lower than rotation's, of old operators and `endif`/`enddo` higher; io, state and setup have none.
- **Star-free vs star-taking:**
  - 91 core procedures take the star.
  - Six files are star-blind: triangle_fit, shell_equations, mass_change_kernel, regridding, envelope_integration_kernel and envelope_integrands. All six are in STAR_BLIND_FILES.
  - 25 procedures outside the tests are `pure`: the regridding, mass-change and triangle kernels (18), the channel sum and its two maps (3), the Monte-Carlo range (2), and one each in envelope_integration and surface_boundary. The shell and envelope kernels are not, as they call the EOS and opacity.
  - Core reads 119 distinct controls directly (298 references). Rotation reads none: it works from `rotation_settings`.
- **README checks:**
  - **Holds:**
    - every module but energy_sources has a §1 row, and every name the README cites exists in the code;
    - §2's five stages match `evolve_step` (but for finding 12's level-4 sentence);
    - the "37 pins" count is consistent with pytest.ini's directories (testsuite 12, config matrix 16, three examples with 1 each, the solar decks 4) plus test_mesa_baselines' 2.
  - **Inaccurate or stale:** the statements in finding 12. No paragraph exceeds 1,000 characters.
- **One recipe traced** (core has none; this is the stop_conditions header's "New stop options belong as further checks inside check_stop_conditions plus their registry rows", for a per-card upper limit on log T_c):
  - the hand edits:
    - the declaration in io/controls_buffer.f90 (as `Teff_upper_limit(max_runs) = 1.0d99`, :152);
    - its registry row;
    - `structure_limit_stop_triggered`'s `nlim`, `qname`, `qval`, `qup` and `qlo` (stop_conditions.f90:197-212);
    - a case in test_stop_conditions.py;
  - then:
    - make (the controls generators);
    - `check_star_contracts --fix`;
    - the pristine re-pin (io/README §5).

  The header's recipe is right but omits the declaration and the test. The value is current at the check, because observables runs first (step_driver.f90:259-268).
- **Not checked:** the coverage of core (nothing was run); the deck pins' byte claims; the mutation claims of the ROADMAP for test_structure and test_mass_change.

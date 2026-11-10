# convection_diffusion -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

Scope: every .f90 of src/convection/ and src/element_diffusion/ read whole (production and test programs), both READMEs, test_convection.py and test_element_diffusion.py, the two expected outputs (the convection one whole, the diffusion one in part), setup/convection_settings_builder.f90, setup/diffusion_settings_builder.f90 and setup/settings_install.f90's install procedures. I also read the call sites in core/step/composition_step.f90 (advance_composition, settle), rotation/mid_timestep_model.f90 (settle_helium), core/step/step_driver.f90 (envelope_overshoot_depth), core/step/convection_geometry.f90 and core/model/rezoning.f90 (the overshoot base), and the settling controls' descriptions in io/controls_buffer.f90. From the tools I read check_boundaries.py's STAR_FREE_DIRS, FACADE_ONLY_DOMAINS, STAR_BLIND_FILES and the scratch checks, and coverage.py's header. From the docs: GUIDELINES "Naming", "Declarations and interfaces" and "Code style"; the bug_sweep entries O1, O10, O16, O21, M33, M67, M96, M97, M112, M129 and M181; ROADMAP's mixing/ split and element_diffusion sections, and its agent-findings lists at ROADMAP.md:400-440 and :616-640. Every F77 claim was checked with `git show 6cd56739:src/<file>.f`, and two with the later F77 commits a8e60674 and 224c9f0d. Nothing in the tree was built or run. In the scratchpad (SP/review7_convection_diffusion/) I wrote a style counter (style.py), a deck scan (deck_settling.py, over the 515 shipped decks), the test model's mass split (cz_mass.py), and a three-line program that prints temperature_gradients' format 20. This is the first score for both directories, so there is no previous report.

## Scores (1-10)

First evaluation. The two directories are scored separately. Both had the architecture pass:
- a facade;
- a settings record built in setup/;
- star-free (STAR_FREE_DIRS);
- diffusion's scratch passed down as an argument.

Neither has had rotation's style, test or README pass. The counts below are production code only (test/ excluded), from SP/review7_convection_diffusion/style_counts.txt:

| | convection | element_diffusion | rotation (the bar) |
|---|---|---|---|
| code lines / procedures | 836 / 8 | 2024 / 17 | 4510 / 112 |
| ALL-CAPS comment lines | 132 | 353 | 0 |
| history comments (initials, dates) | 20 | 35 | 61 |
| direct `write`s in the code | 11 | 17 | 9 (in rotation_messages) |
| modules without `private` | 0 | 4 | 1 |
| ok/FAIL checks in the tests | 0 | 0 | every kernel check |
| README coverage figure | none | none | 96.9% |

### convection/

- **Readability 6.**
  - For: readable variable names, named loops, no GOTO.
  - Against:
    - the prose is the F77's (132 ALL-CAPS comment lines, "PSCALU", "ALPHAC", "LFLAGU", "JMC");
    - the history comments ("2026 de-tramp (ROADMAP item 3): 33 arguments -> 19", "G Somers 9/14 ... G Somers END");
    - the MLT derivative block (temperature_gradients.f90:229-257) is 30 lines of `qa1t`, `qa3p`, `qvt`, with no equation.

  Left:
  - 132 ALL-CAPS lines, 20 history comments, one line of commented-out code (convection_zones.f90:288), 52 `.eq.`-style operators, `dmax1`/`dmin1`/`dabs`/`dsqrt`, `double precision`;
  - every physics routine but compute_scale_height runs 83-292 lines (extend_by_semiconvection 292, compute_temperature_gradients 222), not a short sequence of named steps;
  - the derivative block's algebra, written down once;
  - misleading names: `boundary_mean_molecular_weight` is 1/μ (semiconvection.f90:188), and `max_overshoot_radius` starts as a gradient ratio (:191);
  - run-log formats that overflow (finding 11).
- **Organization 7.**
  - For: one facade (convection_lib), one settings record, installed once per card, and no star.
  - Against: the physics is not all here, and the routines repeat themselves.

  Left:
  - envelope overshoot is computed in five places, with two pressure-scale-height formulas (finding 14);
  - the zone scan and the merge loop are written twice each, the second merge with `>` where the first has `>=`;
  - semiconvection repeats its EOS-opacity-gradient call triplet three times;
  - the kernels write to the run log (M112) and to the terminal, and there is no messages module (finding 15);
  - semiconvection's `only_check_core` leaves a dead branch (finding 1);
  - extend_by_overshoot's fully convective test (convection_zones.f90:270-277) cannot be reached.
- **Maintainability 6.** The tests pin numbers and check nothing: test_convection prints 22 cases, and no line is an ok/FAIL check or was confirmed against a mutation. Semiconvection has no pin at all (finding 12).

  Left:
  - the open bug_sweep items: O1 (the convective gradient's Jacobian gaps), O10 (the overshoot gate and width), M96 (the O16 floor's metals after a Z rescale) and M97 (a Monte Carlo L☉);
  - new findings 1, 2, 5, 6 and 11;
  - literal slots and counts: `composition(1/3/9,…)`, 11 or 15 species, 24 slots in homogenize;
  - unnamed numbers: 1.0d-11 for the velocity, 1.0d+3 and -0.333333333 in the cubic's start, 0.99d-3 for the O16 floor, 10.0d0 in the Castellani speed, and `vtol` without a meaning.
- **Extensibility 6.5.**
  - For: a new parameter follows the settings path (record, builder, registry settings column).
  - Against: the README's recipe names 4 of the 6 to 10 files involved, and a new overshoot prescription, the README's own example, also needs four places outside convection/ (findings 13 and 14).

  Left:
  - the recipe;
  - the overshoot spread;
  - species counts and slots by literal.
- **Docs 6.5.** The README's table and settings list match the code, and its test list matches test_convection.

  Left:
  - no flow section (the order of the composition step's passes: zones, burn, semiconvection, homogenization, settling, zones again);
  - the conventions do not say where messages go;
  - no coverage figure;
  - the claim that cm_semiconvection pins semiconvection is false (finding 12);
  - the recipe is thin (finding 13);
  - the convective zone limit is not stated (finding 5).

Average 6.4.

### element_diffusion/

- **Readability 5.5.** The code is free-form with readable variable names, but it reads as the F77:
  - 353 ALL-CAPS comment lines (94 in original_settling, 90 in microdiff_driver);
  - typos ("COMPATABILITY", "PYHSICAL", "RECALCUALTED", "FACOTR", "CALUCLATE", "hydroden");
  - comments in names that no longer exist: GRTOL, LDOLI, ED_P, ELI_PRIME, EQCOD2X_H, SETUP_LISETT.F, CRSUN_BAH, `star%bl_*_scale`, `rot_scr%metal_*`;
  - "COMPUTE OVERSHOOT (TO BE ADDED)" four times.

  Left:
  - the above;
  - 4 modules without `private` (bahcall_loeb_units, microdiff_mapping, microdiff_run, microdiff_coefficients) and 10 `use` without `only:`;
  - 67 explicit-shape `(json)`/`(15,json)` dummy lines and `double precision` throughout;
  - names that say something else:
    - `metal_abundance_change(_mid)` holds Z until the end of the solve;
    - `eq_metal_diffusion_coeff1_mid` holds dD1/dZ;
    - the scratch type's comment calls both "coefficients";
  - `data` statements for constants (`num_species`);
  - microdiff_coefficients' "CRITICAL" header (finding 16).
- **Organization 6.5.**
  - For: one facade, a settings record, a scratch type with no instance (checked by check_scratch_types_only), the scheme and species table shared by both pipelines.
  - Against: the two pipelines duplicate whole stages, and the step's policy lives outside the domain.

  Left:
  - the region-finding block (about 78 lines) twice;
  - the Coulomb-logarithm block twice;
  - the grid mapping twice, each copy with its midpoint and centre loops written out twice;
  - the Bahcall & Loeb temperature unit and the year defined twice (finding 14);
  - the settling sub-step policy (Michaud's timescale, the count) copied in core's settle and rotation's settle_helium (ROADMAP.md:1687-1689);
  - kernels print, three of them to the terminal only (finding 15);
  - ten of its fourteen routines run 130-334 lines.
- **Maintainability 5.5.** The tests pin numbers only, and M33 makes one of the pins random. No line is an ok/FAIL check. The species' total mass the test prints is not checked (see "Traces").

  Left:
  - O21 (microdiff's 4-point stencil on a 2-point grid) and M33 (stale metal coefficients) open;
  - M181 (diffuse_lithium's overrides, a decision) and M129 weigh in through what the domain is asked to do;
  - findings 3, 4, 7, 8, 9 and 10;
  - the original routines are a second pipeline that no deck runs, carrying four inherited defects that bug_sweep does not track (finding 3).
- **Extensibility 6.**
  - For: species_table makes a light element a row.
  - Against: microdiff_grid hard-codes three of them (finding 9); settling_grid uses literal slots 1-11; the recipe names 3 of about 9 files (finding 13).

  Left:
  - the above.
- **Docs 6.** The README's table matches the code, its test list matches the program, and "no deck uses the original routines" holds (0 of the 354 settling decks; deck_settling.txt).

  Left:
  - §2's in-place unit conversion has been stale since O16 (finding 16);
  - the controls' reach is mis-stated or unstated (finding 17):
    - which controls act only in the original routines;
    - fgrz on Li and Be;
    - constant_mixing_coeff described as "added";
  - "they call it the same way" (finding 17);
  - no flow section (microdiff's pipeline is in microdiff_driver's header only);
  - no coverage figure;
  - the recipe (finding 13).

Average 5.9.

## Summary

- **The architecture is in place; the code and the tests are not.**
  - Both directories have their facade, a settings record built in setup/ and no star, all enforced (FACADE_ONLY_DOMAINS, STAR_FREE_DIRS). Diffusion's scratch has no module instance.
  - Inside, they read like the F77 with longer names: 485 ALL-CAPS comment lines between them (rotation has none), routines of 85-330 lines, and kernels that write to the run log and the terminal themselves.
- **The tests are pinned printouts with no ok/FAIL check and no mutation run.**
  - Semiconvection is not pinned at all: its deck stays fully convective for all 12 models, so the routine returns at its first test.
  - One diffusion pin is random (M33).
- **New bugs, most latent, most inherited:**
  - semiconvection takes the first mixed zone as the core whatever it is: below the envelope when there is no core, and from stale bounds when there is no zone at all (the stale read since the F77's 2025 core-only switch);
  - with starspots, the convective gradient's Jacobian uses the unspotted ∇_rad (three testsuite decks);
  - the original settling routines carry four inherited defects that ROADMAP's 2026 agent lists flagged and bug_sweep never took up:
    - the Noerdlinger derivative's sign;
    - a 2.2/lnΛ bias with Coulomb-log choices 2 and 3;
    - a shifted midpoint index;
    - a convergence test without `abs`;
  - both settling pipelines search the model's zones only up to the equal grid's point count;
  - a twelfth convection zone is dropped silently;
  - the overshoot region's indices go stale;
  - fgrz also scales Li and Be;
  - microdiff's documented "helium off" path reads an unset gradient;
  - `light(3,json)` is a literal;
  - microdiff's "iteration" repeats one identical solve;
  - run-log formats overflow.
- **Docs:** short, accurate READMEs that miss the bar's flow section, coverage figure, messages convention and complete recipes, with one false test claim (semiconvection) and one stale convention (the unit conversion).

## New findings (not in bug_sweep)

Bugs first (1-11), then docs and structure (12-17).

1. **convection/semiconvection.f90:69-71, :106-133, :141-143 (called from core/step/composition_step.f90:139-145)**
   - **Issue:** `only_check_core` sets the loop to zone 1 whatever zone 1 is, and whatever the zone count.
     - With `num_zones_mixed = 0` (a model with no convection zone), find_convection_zones has left `s%mixed_zone_bounds` from an earlier call (it returns at convection_zones.f90:118-127 without writing it). Semiconvection then evaluates, and may move, a stale edge. If the bounds are still zero, it indexes `luminosity_lsun(-1)`.
     - With no convective core, zone 1 is the envelope (or an intermediate zone). Semiconvection is then computed below the envelope's base, the downward branch, and the extended bounds are homogenized at :152-153. The header (:69-70, "only ever evaluated at the convective core") and README.md:20 say core only.
     - Example: cm_semiconvection run past its twelfth model, once a radiative core forms in its 1 Msun pre-MS star.
     - No pinned deck reaches it: cm_semiconvection is fully convective throughout (finding 12).
   - **F77:**
     - 6cd56739's SCONVEC looped over every zone (`DO 210 I = 1,NZONE`, sconvec.f:54), so it computed below envelopes by design and did nothing with NZONE = 0.
     - The core-only switch came with a8e60674 (2025-12-17, F77: `L_ONLY_CORE = .TRUE.` ... `ILOOP = 1`). That commit introduced the zero-zone read, and its "only the core" intent did not test that zone 1 is a core.
   - **Fix:** return when `num_zones_mixed < 1`, and with `only_check_core` when `mixed_zone_bounds(1,1) /= 1`. Or ask the authors whether envelope semiconvection is wanted and say so in the header and README.
   - **Kind:** B (decks with semiconvection and no core; no pinned one). **Size:** S.
2. **convection/temperature_gradients.f90:156-166, :234-235, :251-256**
   - **Issue:** with starspots, ΔV = ∇_rad/S − ∇_ad, with S = f x⁴ + 1 − f (x depending on T when `spot_depth_varies`). The Jacobian terms `qddelt` and `qddelp` still use the unspotted ∇_rad, `((-4 + tempot)*radiative_gradient - qdelat)/deldel`, and omit d ln S/d ln T = 4 f x³(1 − x)/S.
     - At the testsuite's f = 0.25, x = 0.8, the ∇_rad term is 1/S = 1.17 times too small.
     - This sits beside O1's gaps: the converged model is unchanged, but Newton's direction is off in every superadiabatic shell of the three spotted testsuite decks (Test_m0p1/m0p3/m1p0feh+0p0GS98_spot25_ZAMS: f = 0.25, depth-varying), low-mass stars where O3 says convergence is fragile.
   - **F77:** yes (tpgrad.f: `DELDEL = DELR/(SPOTF*DEEPX**4.0 + 1.0 - SPOTF) - DELA`, then `QDDELT = ((-4.0D0 + TEMPOT)*DELR - QDELAT)/DELDEL`).
   - **Fix:** use ∇_rad/S in both terms, add the d ln S/d ln T term for the depth-varying contrast, and add a finite-difference check of the three derivatives with spots to test_convection (test_structure already does this for the unspotted case). Or record it under O1.
   - **Kind:** B (iteration counts and round-off in the three decks). **Size:** S.
3. **element_diffusion/original_settling.f90:640-645, :503-516 with :586-597, :228-232, :258-266; microdiff_run.f90:154-158 and :186-194**
   - **Issue:** four inherited defects of the original routines. They were listed unverified in ROADMAP.md:432-435 and :633-636 (the 2026 agent findings) and never entered bug_sweep. I confirmed each.
     - (a) The Noerdlinger term of dD1/dX has the wrong sign. D1 ∝ 1/lnΛ and d lnΛ/dX = +1.5/(3X+1), so the term is −D1·1.5/((3X+1) lnΛ). The code adds +D1·1.5/(lnΛ(3X+1)). D2's X derivative lacks the term altogether.
     - (b) With `coulomb_log_choice` 2 or 3 and the full Thoul solution, lnΛ enters twice:
       - every pair's `coulomb_log(i,j) = ln_lambda` (:586-590), and the Thoul coefficients scale as 1/lnΛ;
       - the prefactor divides by `ln_lambda` again (:516);
       - the Bahcall & Loeb time unit carries only one 2.2 (bahcall_loeb_units.f90:27-28).

       So settling runs at 2.2/lnΛ of its rate, about 0.8 in the Sun's radiative interior (lnΛ ≈ 2.7). Choices 1 and 4, the fit and the analytic coefficients are consistent.
     - (c) In the iteration, the change at midpoint i−½ is stored at index i (`do eq_idx = 2, n: mid(eq_idx) = 0.5*(X(i)+X(i-1)-…)`), one slot off from the grid that `implicit_diffusion_coeffs` reads (mid(j) lies between points j and j+1). mid(1) keeps the Lax-Wendroff half-step's value. With the original routines' nonzero dD2/dX, each coefficient is corrected with its neighbour's change.
     - (d) Both convergence tests take the largest *signed* change (`if (dx > max)`), so an iteration whose corrections are all negative "converges" at once, and the warning cannot fire.

     No deck uses the original routines (README; 0 of 354 settling decks), so only test_element_diffusion's pins see (a)-(d). In microdiff, (c) and (d) are harmless (finding 10).
   - **F77:** yes, all four:
     - setup_grsett.f:274-279 (`QCOD1X(I)=QCOD1X(I) + COD1(I)*1.5D0/(LN_LAMBDA*(3.0D0*X+1.0D0))`) and :146-150;
     - grsett.f:159-161 (`EX_H(I) = 0.5D0*(EX(I)+EX(I-1)-EX_P(I)-EX_P(I-1))`) and :179-185;
     - microdiff_run.f:137-143 (`IF(DX.GT.DXMAX)`).
   - **Fix:** enter them in bug_sweep; fix (a)-(d) with a test check each, or retire the original routines with the authors' agreement.
   - **Kind:** B (test pins only). **Size:** S.
4. **element_diffusion/microdiff_mapping.f90:82-87; settling_grid.f90:91-96**
   - **Issue:** the first midpoint's linear interpolation searches the *model's* zones `iu = 2 .. num_eq_points`, the equal grid's point count. When the diffusion region starts above model zone `num_eq_points` (a convective core whose top index exceeds ⌈Δr/Δr_min⌉), the search runs out. `iu = num_eq_points` then extrapolates the mass, ρ, T, d ln P/dr, ∇, X, Y, Z (and Li/Be) of midpoint 1 linearly from two points deep inside the core.
     - Inside the range the search finds the same zone as a search from `zone_begin+1`, so the fix is byte-identical wherever the bug does not bite.
     - Latent: the shipped settling decks have more equal-grid points than core zones, as far as their structure suggests; I did not measure it.
   - **F77:** yes (microdiff_mte.f:59-62 and model_to_equal.f:55-58: `DO 15 IU=2,NPT ... IU=NPT`).
   - **Fix:** search `zone_begin+1 .. zone_end` in both.
   - **Kind:** A (byte-identical where it does not bite). **Size:** S.
5. **convection/convection_zones.f90:105-113 (max_convective_zones = 12, state/star_info_lib.f90:85); :492-502**
   - **Issue:** find_convection_zones records a zone only while `j_idx < size(mixed_zone_bounds,1)`, so the twelfth slot is never used.
     - The eleventh zone's end prints "TOO MANY MIX ZONES" even when no zone follows.
     - From the twelfth zone on, every zone, the surface one included, is dropped without an `ierr`. `envelope_cz_edge` then becomes `num_zones` (no envelope) for settling, the Li/Be burn and the rezone.
     - homogenize_convection_zones has its own limit, 12 zones (24 slots), and drops the rest without a word.
     - Latent: no tracked output carries the message.
   - **F77:** yes (convec.f: `IF(J.LT.12) GO TO 11 / WRITE(ISHORT,661)`).
   - **Fix:** record up to `size(...,1)` zones, return `ierr` (or merge the innermost pair) beyond it, and take homogenize's limit from the same parameter.
   - **Kind:** A (no deck reaches it). **Size:** S.
6. **convection/convection_zones.f90:116-127, :138-146, :161-165, :258-261, :355-357 (read at core/step/shell_physics.f90:122, core/structure/structure_equations.f90:169)**
   - **Issue:** the overshoot region that takes the adiabatic gradient (`adiabatic_overshoot_active`) is reset only inside extend_by_overshoot.
     - A model with no convection zone, or a fully convective one, keeps the previous model's `overshoot_bottom`/`overshoot_top`. Its radiative shells there then take ∇_ad.
     - The region is set by *any* zone's downward overshoot, the last one processed winning. With intermediate overshoot on and envelope overshoot off, it is an intermediate zone's region, described everywhere as "below the envelope" (README.md:39, star_info_lib.f90:670-675).
     - Latent: no shipped deck sets `adiabatic_overshoot_active` (520 decks set it false) or intermediate overshoot.
   - **F77:** yes (oversh.f:25-26 resets IOV1/IOV2; :108-109 sets them in every downward branch).
   - **Fix:** set both to −1 on find_convection_zones' early returns and when no overshoot flag is on; set them only for the zone that reaches the surface.
   - **Kind:** B in principle (no deck). **Size:** S.
7. **element_diffusion/microdiff_coefficients.f90:152-160 with microdiff_driver.f90:226-229**
   - **Issue:** microdiff scales column 3's coefficients by `fgrz`, the metals' settling factor. The light elements are diffused in column 3 too, so Li6, Li7 and Be9 settle at `settling_z_factor` times their rate.
     - The unscaled default at :154 is dead code: no call uses another column.
     - The control's description ("the factor on the metals' settling coefficients", io/controls_buffer.f90:1141-1143) and the README do not say this.
     - Every shipped deck sets the factor to 1.0; a Monte Carlo run that samples it (M151) moves Li and Be with it.
   - **F77:** not at 6cd56739 (microdiff_cod.f:112-113, no factor). It came with 224c9f0d (2026-01-26, F77, "Made FGRY and FGRZ work with new diffusion routines": `IF(J.EQ.3) FAC=FGRZ*...`), where J = 3 serves both.
   - **Fix:** the authors' choice. Either scale only the metal call (pass the factor in from microdiff), or document it in the description and README.
   - **Kind:** A to document; B to change (Monte Carlo runs only). **Size:** S.
8. **element_diffusion/element_diffusion_lib.f90:29-31, :48; microdiff_driver.f90:91, :106, :167; microdiff_coefficients.f90:215-217**
   - **Issue:** `diffuse_helium` is always true.
     - Both callers call diffuse_elements only while `s%settling_active` is true (composition_step.f90:161-163 and :262-265; mid_timestep_model.f90:462-465 and :483-486), and pass that same flag.
     - The original routines ignore it.
     - So the documented case "the metals and light elements settle regardless" (element_diffusion_lib.f90:30-31; README.md:5, "the two values that change within a card") never happens.
     - If it did, the metal and light-element calls would read `hydrogen_dlnc_dr`, which only the hydrogen call writes: an unset local, like M33's.
   - **F77:** yes (mix.f:377-379 and :420 call MICRODIFF only under LDIFY; microdiff.f:115 tests it again).
   - **Fix:** drop the argument (and the branch), or zero the gradient and test the case.
   - **Kind:** A. **Size:** S.
9. **element_diffusion/microdiff_mapping.f90:19-23 against species_table.f90:46-55**
   - **Issue:** `microdiff_grid%light` is dimensioned `light(3,json)`, a literal. species_table presents the light elements as a table (`num_light_diffused = 3`, one row each), so a newcomer adding a fourth (B11, say) bumps the count and the rows. map_microdiff_to_equal_grid and microdiff_etm then write `eq%light(4,i)` past the component's end, with no bounds check in a release build.
   - **F77:** the arrays were sized `ELI(NLIGHT,JSON)` from the argument (microdiff_mte.f); the literal is the 2026 type's.
   - **Fix:** `light(num_light_diffused, json)`.
   - **Kind:** A. **Size:** S.
10. **element_diffusion/microdiff_run.f90:151-210 (with diffusion_scheme.f90:160-166)**
    - **Issue:** microdiff's "iteration on the diffusion coefficients" (microdiff_driver.f90:36-37) recomputes nothing. It passes zero derivatives, so `implicit_diffusion_coeffs` adds 0 to D2, and every pass solves the same tridiagonal system for the same right-hand side.
      - Iteration 1's "change" is the implicit step itself, and iteration 2 always returns change 0.
      - So `settling_tolerance` and `settling_num_iterations`, documented as the solver's convergence controls, decide only whether a second identical solve runs: two per species per sub-step at the defaults.
      - With `settling_num_iterations = 1`, "MICRODIFF FAILED TO CONVERGE" prints whenever the step moves anything by more than the tolerance.
    - **F77:** yes (microdiff_run.f: `EQCOD2X_H(I) = 0.0`, `CALL GET_IMP_DIFFCO(...)`, `CALL TRIDIAG_GS(...)` in the loop).
    - **Fix:** solve once in microdiff (byte-identical) and say in the README that the controls act only in the original routines.
    - **Kind:** A. **Size:** S.
11. **convection/temperature_gradients.f90:210-211; convection_zones.f90:177-178, :387-392; semiconvection.f90:309-310, :323-325**
    - **Issue:** run-log formats that cannot hold their values.
      - `1p3e12.6` leaves no room for a sign or a separator. In "CUBIC NON-CONVERGENCE", dρ/dT (always negative) prints as `************`, and κ and c_P run together. In a scratch run of the format: `5.000000E+032.000000E+08************`.
      - The merge messages print zone indices with `i3`. A run of examples/run_from_zahb_to_tahb's ZAHB_m1p0feh+0p0_TAHB deck (its untracked output/ log here) prints `OLD[  1-***][***-***] NEW[  1-***]` 330 times, and semiconvection's edge message has the same `i3`.
      - The `4('[',i4,'-',i4,' ]')` formats end every line with a stray `[`.
    - **F77:** yes, all (tpgrad.f:135-136, convec.f's format 93, sconvec.f's 601 and 93).
    - **Fix:** `1p3e13.5` (or `es`), `i5`, and a format that stops at its last item. Do it when the messages move to one module (finding 15).
    - **Kind:** A (run-log text; test_convection's pin moves). **Size:** S.
12. **convection/README.md:56; convection/test/test_convection.f90:18-19**
    - **Issue:** "Semiconvection needs the opacity tables, so the configuration matrix's semiconvection deck pins it." cm_semiconvection rescales the 1 Msun deuterium-birthline model (`m1000gs98z018804_Dbl.first`) and evolves it for 12 models. The local reference log (examples/run_config_matrix/standard/cm_semiconvection.log; test_all.py's references are not tracked) shows log Teff 3.647, log L 0.80 and ages 1.6-1.9e5 yr, a Hayashi-track star that is fully convective for another million years or so.
      - So `extend_by_semiconvection` returns at semiconvection.f90:122 on every call (zone 1 spans the model). No test or deck computes a semiconvective extension, a merge, or the Castellani limit.
      - The deck's own comment says semiconvection is for "core He-burning models".
      - The two UpperMS sample decks that use it are not pinned.
    - **Fix:**
      - give cm_semiconvection a model with a convective core (a ZAHB or upper-MS start), or add a test_convection case on a synthetic model with the opacity tables test_kap already loads;
      - correct the sentence.
    - **Kind:** A (docs and test). **Size:** S/M.
13. **convection/README.md:58-61; element_diffusion/README.md:44 (recipes, traced)**
    - **Issue:** both recipes name too few of the files a newcomer must touch.
      - **convection, "Add a convection treatment (for example a different overshoot prescription)"**, traced:
        - the control: io/controls_buffer.f90 and defaults/controls_registry.tsv with its settings cell, per io/README §5;
        - the field (convection_configuration.f90) and the copy (setup/convection_settings_builder.f90), then `check_star_contracts --fix` for the builder's `star reads`;
        - the physics, then test_convection and its expected output;
        - for an envelope overshoot rule, also core/step/step_driver.f90:688-701 (`envelope_overshoot_depth`, for the Li/Be burn), core/step/convection_geometry.f90:126-…, core/model/rezoning.f90:309-318, :876, :1017 (the rezone's overshoot base, with α used as a depth in log P) and rotation/evolve_rotation.f90:390-391 (O10's gate).

        The recipe names 4 of these 6 to 10 files, and none of the steps after them.
      - **element_diffusion, "Add a settling treatment"**: the module, the control (two io files), the settings field and its copy, `diffuse_elements`, deps.mk (`python3 tools/gen_makefile_deps.py`; CI runs `--check`), the test and its expected output, and the README's table. The recipe names three.
      - Neither README says that the settings record is the only way in (rotation's "Add a control" line).
    - **Fix:** write both as rotation's numbered recipes, naming every file and command. For overshoot, first move the extent into convection/ (finding 14).
    - **Kind:** A. **Size:** S.
14. **Formulas and constants that exist more than once (GUIDELINES "Code style" rules 2-3)**
    - **Issue:**
      - **convection:**
        - the zone scan (convection_zones.f90:87-114 and :475-502);
        - the merge loop (:168-196 and semiconvection.f90:313-342, with `>=` against `>`);
        - the EOS-opacity-gradient triplet (semiconvection.f90:164-183, :197-215, :268-286);
        - H_P = P r²/(G m ρ) twice: compute_scale_height (convection_zones.f90:444-445, the EOS's ρ at the edge) and core/step/convection_geometry.f90:44-52 (the model's ρ, with physical_constants' `cgl` where convection takes G from its settings);
        - the envelope-overshoot extent in five places with three definitions: α H_P in radius in convection_zones, α H_P in step_driver and convection_geometry, α in log P in rezoning.
      - **element_diffusion:**
        - the diffusion region and the three suspension notes (microdiff_driver.f90:323-401 = original_settling.f90:404-480);
        - the Paquette Coulomb logarithm with its five unnamed, unsourced constants 0.23873 (3/4π), 6.9010, 2.3939e3, 0.81245, 0.18769 (original_settling.f90:542-584 = microdiff_coefficients.f90:84-142);
        - the equal-grid construction and the 4-point search (settling_grid.f90:58-200 against microdiff_mapping.f90:50-247, each with its midpoint and centre loops written twice);
        - the Bahcall & Loeb temperature unit (bahcall_loeb_units.f90:26 and microdiff_coefficients.f90:59, `bl_temp_scale_local`);
        - the year (bahcall_loeb_units.f90:19 `seconds_per_year_bl = 3.1558d7`, while physical_constants exports `seconds_per_year`);
        - 2.2 as `2.2d0` in one pipeline and as single-precision `2.2` in the other;
        - the settling sub-step policy in two callers (ROADMAP.md:1687-1689).
      - Unexplained: 2.7d13 in the time unit, 1.0d-24 for the mass-fraction floor, and the 20 layers of the grid spacing.
    - **Fix:** one procedure per formula (the region, the Coulomb logarithm, the grid, H_P from convection_lib), and named, sourced parameters. Move the overshoot extent and the sub-step policy into their domains. Merging copies that differ is a measured commit (GUIDELINES).
    - **Kind:** A for the identical copies; B where they differ (H_P, the overshoot extents, 2.2). **Size:** M.
15. **Messages: convection/temperature_gradients.f90:181-188, :208; convection_zones.f90:111, :174, :272-274, :385-392; semiconvection.f90:306, :320; element_diffusion/original_settling.f90:252, :269-281, :417-479; microdiff_run.f90:180, :197-206; microdiff_driver.f90:338-400; thoul_solver.f90:289**
    - **Issue:** both domains' kernels write their own messages, against GUIDELINES "Code style" rule 7 and rotation's `rotation_messages`.
      - The run-log writes are M112's.
      - Four failure messages go to the terminal only (`write(*,…)`, not the run log): the MLT square root, the two `solve_tridiagonal` failures and `ludcmp`'s. So a run that stops on them leaves no trace in its log.
      - extend_by_overshoot writes two "MIXED REGIONS" lines to the run log on every call, with no gate. find_convection_zones runs several times per model (composition step, step driver, rotation sub-steps), while composition_step's own zone line is gated by `solver_diagnostics()`.
    - **Fix:** a `convection_messages` and a `diffusion_messages` module (report to terminal and run log, note under the gate), the kernels returning `ierr` with the text built by their callers.
    - **Kind:** B for the run logs (the pins of the decks that print them move). **Size:** M.
16. **Stale or inexact text**
    - **Issue:**
      - element_diffusion/README.md:30 says the routines convert the timestep, masses and d ln P/dr "in place and convert them back... The round trip can change the last bit of star%m; that is the existing behaviour". Since O16 (c66901c9) microdiff, every deck's pipeline, converts copies; only the original routines still do the round trip.
        - In the same way, bahcall_loeb_units.f90:4-6 says microdiff_etm reads the scales back to undo it.
        - microdiff_etm's header (microdiff_mapping.f90:250-255) says it "converts DT/STOT/HS1/HQPR back", which :405-407 then denies.
        - microdiff_etm still takes five arguments it never reads (`diff_scr`, `timestep`, `dlnp_dr`, `enclosed_mass`, `total_mass`; :258-261).
      - microdiff_coefficients.f90:4-8: "0.0, 1.0, 0.5, 2.2 ... numerically different from a D0-suffixed literal". Only 2.2 differs; 0.0, 1.0 and 0.5 are exact in single precision.
      - original_settling.f90:494-499 states D1 with "X(1−X)" and "DEL*6*(X−0.32)", where the code (and its X derivative, 3.384 + 2.88X + 4.5X²) has X(1−X−Z−He3) and (X+0.32). The F77 comment has the same disagreement (setup_grsett.f:140).
      - "BAHCALL AND LOEB 1989" (microdiff_driver.f90:20, original_settling.f90:20) against the README's 1990.
      - temperature_gradients.f90:18-26 and :189-191: provenance history, and the legacy file names core/main, core/crrect, core/starin, setup/hpoint.
      - convection_configuration.f90:33-35: "the job flags lovstm etc.", "lovmax", "betac".
    - **Fix:** text, and drop the five arguments.
    - **Kind:** A. **Size:** S.
17. **The controls' reach (element_diffusion/README.md:5, :32; io/controls_buffer.f90:1139-1165)**
    - **Issue:** what each setting does is mis-stated or left out.
      - README.md:32 names `use_thoul_fit` and `use_thoul_diffusion` as acting only in the original routines. `coulomb_log_choice`, `constant_mixing_coeff` and `constant_settling_reduction` do too: microdiff always uses the pairwise Coulomb logarithm (microdiff_coefficients.f90:126-142) and neither factor.
      - `constant_mixing_coeff` is described as "Constant added diffusion coefficient (mixing)" (controls_buffer.f90:1162). It is a factor on D1 and D2, applied only in the Thoul branch (original_settling.f90:607-609, :630-635).
      - fgrz reaches Li and Be (finding 7).
      - diffuse_lithium's overrides (M181) are documented in io but not in this README.
      - README.md:5's "they call it the same way" does not hold for the flags the two callers pass:
        - core passes the deep-mix flags, the convective flags without overshoot (composition_step.f90:304-307);
        - rotation passes the mixed regions with overshoot (mid_timestep_model.f90:554-567).

        So the settling region's top differs with rotation.
    - **Fix:** a short "which settings act where" table in the README; correct the two descriptions.
    - **Kind:** A. **Size:** S.

## What would raise each score

**convection/**
- **Readability, to 8-9:**
  - the prose rewritten in sentences (no ALL-CAPS), the history moved to a Sources section, the commented-out line removed;
  - `<`/`==`, `real(dp)`, `max`/`min`/`abs`/`sqrt`;
  - compute_temperature_gradients split into named steps (radiative gradient, spot correction, cubic, velocity, derivatives) with the derivative equations written once;
  - the two misleading names;
  - finding 11's formats.
- **Organization, to 8.5-9:**
  - the overshoot extent and H_P owned by convection/ and called from core and rotation (finding 14);
  - one zone scan and one merge;
  - semiconvection's triplet as one procedure;
  - a `convection_messages` module (finding 15);
  - the dead switch and branch removed.
- **Maintainability, to 8.5-9:**
  - ok/FAIL kernel checks, each confirmed by a mutation:
    - ∇ between ∇_ad and ∇_rad;
    - the MLT efficient and inefficient limits;
    - the derivatives against finite differences, spots included;
    - homogenization conserving each species and leaving a uniform composition alone;
    - overshoot extents against α H_P;
    - a merge;
  - a pinned semiconvective extension (finding 12);
  - findings 1, 2, 5 and 6;
  - named slots and parameters for the literals;
  - O10 and M96.
- **Extensibility, to 8.5:** finding 13's recipe, traced; the overshoot spread removed; the species count taken from composition_layout.
- **Docs, to 8.5-9:**
  - a flow section;
  - the messages convention;
  - a coverage figure (`src/tools/coverage.py --dir convection` exists);
  - a list of what each test checks;
  - findings 12, 13 and 16.

**element_diffusion/**
- **Readability, to 8:**
  - the prose rewritten (353 ALL-CAPS lines), the legacy names in comments and the typos fixed;
  - `private` and `only:` everywhere;
  - assumed-shape arrays;
  - the two misleading scratch names (`metal_abundance_change` → `eq_metal_fraction`, …);
  - `data` statements as parameters;
  - finding 16.
- **Organization, to 8:**
  - one region finder, one Coulomb logarithm, one grid builder;
  - the sub-step policy in the domain;
  - a messages module;
  - a decision on the original routines (fix and test them, or retire them to `retired/` with the authors' agreement).
- **Maintainability, to 8-8.5:**
  - ok/FAIL checks confirmed by mutations:
    - each species' total mass conserved in a consistent cell convention;
    - X+Y+Z+He3 = 1;
    - the light elements sinking;
    - the change linear in a small timestep;
    - zero change at zero timestep;
    - Thoul's coefficients against the published fit at their sample points;
    - the Coulomb logarithm against Paquette's limits;
  - M33 fixed so that the pins are deterministic;
  - O21 and findings 3, 4, 8, 9 and 10;
  - a case with a convective core (zone_begin > 1), the path finding 4 is on.
- **Extensibility, to 8:** finding 9; named slots in settling_grid; finding 13's recipe.
- **Docs, to 8:**
  - a flow section (microdiff's stages);
  - "which settings act where" (finding 17);
  - the messages convention;
  - a coverage figure;
  - findings 16 and 7's documentation.

## Traces, counts and README checks (supporting)

- **Deck counts** (deck_settling.py over examples/, testsuite/ and sample_cases/: 515 decks):
  - 354 decks turn settling on, all with the new routines or diffuse_lithium, so 0 reach the original routines;
  - 3 decks use semiconvection (cm_semiconvection and two unpinned UpperMS 9 Msun samples);
  - 3 use starspots (the testsuite's spot25 decks, depth-varying);
  - 0 use adiabatic overshoot or intermediate overshoot;
  - all set `settling_z_factor = 1`.
- **README claims checked and true:**
  - convection §1's table, against the modules (convection_lib exports the four procedures, the two setters and the type);
  - §2's settings list, against convection_settings_t and the builder;
  - the install in `begin_kind_card` after the card's α (run_driver.f90:343-346, settings_install.f90:100-112);
  - §2's overshoot-region paragraph, against the arguments;
  - §4's gradient and zone case lists, against test_convection;
  - element_diffusion §1's table;
  - the install path;
  - the two callers;
  - "a singular matrix returns ierr";
  - the test list;
  - "no deck uses the original routines".
- **README claims not true:**
  - convection README.md:56 (finding 12);
  - element_diffusion README.md:30 (finding 16);
  - README.md:5 "the two values that change within a card" and "the same way" (findings 8 and 17).
- **Conservation in test_element_diffusion's own printout** (cz_mass.py, from build_model's formulas):
  - in "new: helium, metals and light elements" the surface convection zone (3.48% of the 1.25 Msun model, zones 143-200) gains ΔX = 1.177e-3, 4.10e-5 of the star's mass in hydrogen;
  - the printed mass-weighted change of X over the star is +2.67e-6, 6.5% of that; Z is −1.46e-7, 6.9% of its convection-zone loss;
  - one zone near the base weighs 1.6e-3 M, so 1.4 zones' worth of the convection zone's change, which the test's dm convection (dm(i) = m(i) − m(i−1)) could account for. So this is not a demonstrated defect, but it is the check the test should make, in a consistent cell convention.
- **Not checked:**
  - the coverage of either directory;
  - whether any shipped deck's diffusion region starts above model zone `num_eq_points` (finding 4);
  - the UpperMS semiconvection decks' zone 1 after the main sequence (finding 1).

Nothing was run.

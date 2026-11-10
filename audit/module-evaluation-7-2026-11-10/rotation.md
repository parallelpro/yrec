# rotation -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

Scope: every hand-written file of src/rotation/ read whole: the 24 live .f90 (7,607 lines: the top level, shape/, seculr/), rotation/README.md and retired/README.md, the headers and call-site notes of the five retired/*.f90.txt (1,695 lines, skimmed past their headers), both test programs (test_rotation.f90 whole; test_rotation_kernels.f90's header, its program block and the checks the findings touch) with their pins, and test_rotation.py, test_rotation_kernels.py and conftest's compare_output. Also read whole: setup/rotation_settings_builder.f90, state/transport_mechanism_layout.f90, the rotation entries of tools/check_boundaries.py (PUBLIC, STAR_BLIND_FILES, SCRATCH_*, CONTROLS_FREE_DIRS, the cgl check) and check_star_contracts.py's coverage rule; docs/GUIDELINES.md "Naming", "The star object", "Declarations and interfaces" and "Code style"; docs/bug_sweep.md section 6 (O11-O29), O10, O13, O16, M112, M117; docs/ROADMAP.md's rotation stages (2026-10-12 series, Stages B-E) and every later rotation mention; docs/TARGET_ARCHITECTURE.md's rotation parts; docs/HOW_YREC_WORKS.md §4.6; docs/ROTATION_VALIDATION.md §3-§5; io/README.md §5 (for the control recipes). F77 evidence from `git show 6cd56739:src/<file>.f`. Nothing was built or run; check_star_contracts.py and check_boundaries.py were run in check mode (`python3 -I -B`, no file written): both pass ("Star contracts OK (149 declared)", "Domain boundaries OK"). Counts come from scratch Python over the tree (SP/review7_rotation/scan.py). There is no previous module-evaluation report of rotation/, so section 4 is left out; the open bug_sweep and ROADMAP items are in the summary and the supporting section. Paths are under src/ unless they start at the repository root (docs/, test_*.py).

## Scores (1-10)

History: the 2026-10-12 rating (5.5/10 overall), the 2026-10-28 evaluation that started Stage C (out of 5, doubled: readability 7, organization 5, maintainability 5, extensibility 4; docs not scored), the re-evaluation before Stage D ("near 7 of 10"). rotation/ was then made the bar and never scored on it. The 10-28 readability was judged before the 2026-10 declaration and style rules were applied to rotation; now they are, as written.

- **Readability 7 -> 7.5.** The comments are the best in the tree (each routine's equations, units, sources and caller), and the drivers read as named steps; but the code has not had the declaration and style pass the other packages had. Left:
  - 344 `double precision` statements against 67 `real(dp)`; 117 `.lt.`/`.eq.`-style operators; 60 `endif`/`enddo`; `dfloat`, `dabs`, `dsqrt`, `dmin1`; `data` statements for constants; 42 lines over 110 characters (finding 15);
  - 86 lines of commented-out code in five live files, two of them inside a star-blind kernel, and 45 "(commented out ...)/(now in retired ...)" marker lines (finding 8);
  - all eight message formats are F77 capitals naming routines that no longer exist, one of them false (finding 9);
  - evolve_angular_momentum's header is F77 history (evolve_rotation.f90:21-37: "11/91 MHP revised subroutine (replaces old FINDW)", LINSTB, LJDOT0, LEXCOM, "11/91 JENV0 added to call"), and its local comment at :101-107 still describes a GOTO 9999 and MIDMOD;
  - names the code apologizes for instead of fixing: transport_checks.f90:308-310 ("Despite their names, ion_mean_weight_inverse and electron_mean_weight_inverse hold ... mu_ion and mu_e"); `shape` shadows the intrinsic.
- **Organization 5 -> 8.** The layout (shape/, seculr/, retired/), the settings record built in setup/, the types-only scratch, the star contracts and the mechanism slots are sound; the gaps are in how the rules are applied. Left:
  - five star-free files are kept out of STAR_BLIND_FILES only because they take `json` from star_info_lib, while three listed kernels reach star_info_lib through rotation_scratch_lib, and README.md:14 states a rule the list does not follow (finding 10);
  - substep_time_fraction is a star member used as a within-step channel, wind_now exists twice, and the settings record is public and not `protected`, as net's, kap's, eos's and atm's are (finding 11);
  - 37 steps are internal procedures sharing their host's locals (evolve_angular_momentum's 8 share about 30), where GUIDELINES "Code style" rule 1 asks for explicit arguments, as Stage D did for apply_secular_transport (finding 15);
  - acknowledged in ROADMAP (Stage D "Still open" 2-3): rotation and core use each other (shell_physics, burn_step, convection_geometry), and 22 procedures take 14 to 24 arguments (diffuse_composition_driver 24, am_transport_grid 21).
- **Maintainability 5 -> 8.** test_rotation_kernels' 92 ok/FAIL checks (limits, scalings, conservation, eigenmodes, convergence orders, contracts) are the strongest unit tests in the code, with a validation document and a MESA comparison. Left:
  - O29 and O10's gate (evolve_rotation.f90:390) open in bug_sweep;
  - new latent bugs: the composition write-back overwrites a neighbouring convection zone (finding 1), non-positive sub-step controls hang the step (finding 2), a redone sub-step advances the deuterium rates twice (finding 3), and two small physics defects (findings 4, 5);
  - each torque law written twice (finding 12; O17 was a divergence of these copies);
  - no ok/FAIL check covers viscos, the Zahn factor, diffusion_velocity_scales' bookkeeping, the power-law convection-zone modes 2 and 3 or disk locking in solid_body_omega (test_rotation pins them as numbers);
  - the coverage figure predates the error paths added on 2026-10-29 (finding 14).
- **Extensibility 4 -> 8.5.** A transport mechanism is a slot, a kernel and a case (Stage E added Solberg-Hoiland by the recipe in a day), and a control follows io/README §5. Left:
  - the mechanism recipe omits the control's declaration in io/controls_buffer.f90 and the re-pins (finding 13);
  - a braking law of a new form must be written into both copies of the torque (finding 12); the recipe (README.md:120) covers only new coefficients of the two existing forms;
  - explicit-shape `json` and `composition(15,json)` dummies everywhere (245 explicit-shape against 62 assumed-shape array dummies), so adding an isotope still edits rotation (Phase 5, acknowledged).
- **Docs (not scored) -> 8.** README.md has every part the bar asks for (what lives where, the step, conventions, what each test checks, coverage, recipes), and no paragraph over 1,000 characters (the longest prose paragraph, :14, is 745). Left:
  - claims the code does not keep: "Kernels print nothing" and "Messages go through rotation_messages" (:74; findings 6, 9), the STAR_BLIND rule (:14; finding 10), "two services" of core (:67), retired/README's "moved ... so that those read as the code that runs" and "each file's header says how to revive it" (README.md:48; finding 8);
  - stale comments in six files and stale lines in TARGET_ARCHITECTURE and HOW_YREC_WORKS (finding 14); GLOSSARY's 25 dead rotation paths (ROADMAP's still-open item 7).

Average 8.0 (7.5 / 8 / 8 / 8.5 / 8).

## Summary

- **What rotation still does best.** The README's shape (now every package's), the settings record and the controls-free check, the star contracts (required for rotation since b621c57e; all 16 star-taking procedures declare them and the check passes), the mechanism slots with their layout-to-registry test, and test_rotation_kernels: 92 behaviour checks against analytic laws, each documented as confirmed by a mutation (not re-run here). Of the open rotation items in bug_sweep only O29 (beyond break-up) and O10's gate (lovstm for lovste in the sub-steps, evolve_rotation.f90:390, rotation_settings_builder.f90:147) remain; M112's class covers rotation's direct run-log writes (rotation_messages.f90:35-42 through write_lines :62, mid_timestep_model.f90:538-548), and M32's inert Matt-law solar mass-loss rate is a rotation control. ROADMAP's Stage D "Still open" list is accurate (O29; the sub-step's composition through core; long argument lists; the constant-j law's physics).
- **Where the other packages now exceed it:**
  - io, state and setup carry no F77 shapes, and net, numerics and io follow the declaration rules; rotation keeps 344 `double precision`, 117 old operators and 245 explicit-shape dummies;
  - net, kap, eos and atm make their settings record `protected`; rotation's is public and writable;
  - io's recipes are traced and dated after each change, with a cells guide that make prints; rotation's mechanism recipe was last exercised by Stage E and now misses a file;
  - numerics and math hold no I/O, checked; rotation has a printing kernel.
- **New bugs, all latent in the shipped decks:**
  - O13's fix widened the composition write-back without changing it from a copy, so a convection zone beside a radiative end zone (and zone 1 or the surface zone, which the transport flags convective) takes the end zone's abundance, not its change; the test cannot see it (finding 1);
  - `diffusion_timestep_factor <= 0` or `max_domega_global < 0` in a rotating deck with instability transport sends the sub-step clock backward and hangs the step; nothing checks them (finding 2, F77);
  - a redone sub-step advances the sub-step deuterium rates twice (finding 3, F77);
  - the equipotential's theta-gradient term is dimensionally wrong (finding 4, F77, O(A^3) in g) and viscos' electron density is too high by 1/sum(X_i/A_i) (finding 5, F77, about 4% in nu deep inside).
- **Contract slips:** a star-blind kernel prints to the terminal (finding 6); advance_rotation_substep's contract says it writes `m` because burn_and_settle kept an `intent(inout)` that O16's fix made unnecessary, with a comment that is now false (finding 7).

## New findings (not in bug_sweep)

1. **rotation/seculr/rotational_mixing.f90:282-295 (with :265-276, secular_transport.f90:198, :229, :273-278)**
   - **Issue:** O13's fix (ed8f7af6) made diffuse_composition widen i0..i1 through the convection zones beside the region "whether or not the end zone is convective itself", but kept the write-back a copy: the convection zone below takes `composition(.,zone_begin)` (the end zone's new value), the one above `composition(.,zone_end)`. Before, the copy ran only when the end zone was itself in the convection zone, where both hold the same abundance, so a copy was the change. Now, with a radiative end zone whose X differs from the convection zone's, the zone's abundance is overwritten: the species' mass changes by m_cz (X_end - X_cz) more than the grid's change (the end cell's mass includes m_cz, composition_grid :433-470). The same holds at zone 1 and the surface zone, which flag_end_zones_convective marks convective during the transport: a region that starts at zone 2 or ends at zone nz-1, both radiative, copies their abundance into the centre or the surface zone. The angular momentum path is right: equal_grid_to_model multiplies by (1 + dJ/J) (am_grid.f90:396-420) and rescales. The kernel test (test_rotation_kernels.f90:1043-1057) gives the convection zone the end zone's composition before the call, so a copy and a change agree and it passes either way; ed8f7af6's message ("which now takes the end zone's change") describes the intent. O13 records that no shipped rotating deck reaches the case (the 26 were byte-identical).
   - **In the F77?** No: MIXCOM copied only under `IF(LCZ(IBEG).AND.IBEG.GT.1)` and `IF(LCZ(IEND).AND.IEND.LT.M)` (mixcom.f:97, :110, :137, :145), where the zone and the convection zone are one mixed region.
   - **Fix:** add the end point's change, `eval_y(1)` below and `eval_y(ntab)` above, to each widened zone (as equal_grid_to_model does for j), and give the test's convection zone an abundance different from the end zone's.
   - **Kind:** B (only runs that reach O13's case). **Size:** S.
2. **rotation/evolve_rotation.f90:278-286 (with :291-299, :167-208; io/controls_check.f90:91-98)**
   - **Issue:** count_substeps divides by `diffusion_timestep_factor` (inlist `time_dw_mix`) and counts `max_domega_global` (`time_dw_global`), and nothing checks either: check_controls validates only `max_diffusion_iters` among rotation's controls. In a deck with rotation and instability transport, a factor <= 0 or `max_domega_global < 0` gives a sub-step count <= 0 (or `int` of an infinity, undefined, followed by an overflowing `+ 1`), a non-positive or negative `sub_timestep`, and advance_substep_clock then moves `elapsed_substep_time` away from `full_timestep`: the loop's exit test (:207) never holds and the step hangs, interpolating backward in time. `max_domega_global = 0` gives a count of 0 and an infinite sub-step: one sub-step over the whole step, the omega-change limit silently off. Every value the tracked decks set is positive.
   - **In the F77?** Yes, unchecked in the same arithmetic (getw.f:170-177: `NSTEP = INT(DWMAX/DTDIF) ... DT = DELTS/DFLOAT(NSTEP)`).
   - **Fix:** in check_controls, for a rotating run with instability transport, require both > 0 (or document 0 as "no limit" and code it); one message naming the inlist controls.
   - **Kind:** A (no shipped deck changes). **Size:** S.
3. **rotation/mid_timestep_model.f90:248-263 (with evolve_rotation.f90:364-374)**
   - **Issue:** the sub-step deuterium rate is accumulated, `deuterium_rate_start = deuterium_rate; deuterium_rate += (sub_timestep/full_timestep)(r1 - r0)`. When the transport asks for a redo (a cut by check_composition or check_mixed_species), apply_instability_transport rewinds the clock, j, mu ("the caller takes the share back out on a redo", :175-179) and the composition, but not these rates, so the redone sub-step starts from the failed attempt's end rate and every later sub-step of the step burns deuterium ahead of the clock by the failed share. Reached only with surface deuterium above 1e-14 (pre-main sequence), instability transport and a cut after the first sub-step.
   - **In the F77?** Yes: MIDMOD's DRATEM accumulates the same way (midmod.f:283-291), and GETW's redo rewinds AMUM, HJM and HCOMP only (getw.f:283-293).
   - **Fix:** form both rates from the clock, r0 + f (r1 - r0) at the sub-step's start and end fractions (no accumulation, nothing to rewind), or rewind them with mu.
   - **Kind:** B (latent). **Size:** S.
4. **rotation/shape/shell_quadrature.f90:180-184**
   - **Issue:** `qphith = const*(p2*q*qr0t - 3.0d0*aint*cs) - w2*r*cs`. With Psi as the header writes it (:149) and qr0t = -(1/r) dr0/dtheta (:168), (1/r) dPsi/dtheta of the quadrupole term is -(4 pi/3)(r P2 q qr0t + 3 aint cos sin)/r^4: the code's first term lacks the factor r (q qr0t is cm^4 s^-2 where 3 aint is cm^5 s^-2) and has the opposite sign, as the comment at :182-183 itself notes. The term is in effect dropped (divided by r ~ 1e10 cm). The correct term is a fraction of order A of the aint term, so <g> and <1/g> move at O(A^3) relative: round-off for the Sun (A ~ 2e-5), of order 1e-3 for a star near break-up.
   - **In the F77?** Yes (func.f:32, `QPHITH = CONST*(P2*Q*QR0T-3.0D0*AINT*CS) - W2*R*CS`).
   - **Fix:** `const*(-r*p2*q*qr0t - 3.0d0*aint*cs)`, and a test of fp, ft at O(omega^6) or against a Roche model; a measured commit.
   - **Kind:** B (round-off for the shipped decks, larger only for fast rotators). **Size:** S.
5. **rotation/viscosity.f90:71-75 (comment :58-61)**
   - **Issue:** `electron_number_density = mean_charge*density_cgs/amu` is Z_bar rho/m_u; the electron density is Z_bar (rho/m_u) sum(X_i/A_i), as the comment says ("the two agree only when sum n_i = 1"), and the lumped metals count as neutral particles of weight 1 in Z_bar. For X = 0.70, Z = 0.02 the code's n_e is 1.26 times the fully ionized value, so the Coulomb logarithm is 0.12 lower; in the solar core, where it is about 2.6, the ion viscosity (which dominates nu there) is about 4% high. nu enters the secular-shear threshold (sqrt of the Prandtl number) and the GSF's viscous inhibition. Not in bug_sweep.
   - **In the F77?** Yes (viscos.f:43, `ELECN = ZF*RHO/AMU`).
   - **Fix:** n_e = (rho/m_u) sum(X_i Z_i/A_i), the metals' electrons included; a measured commit, or record it in bug_sweep as kept.
   - **Kind:** B. **Size:** S.
6. **rotation/seculr/instability_velocities.f90:85-91 (README.md:74; test/expected_test_rotation_kernels.out:83-86)**
   - **Issue:** compute_eddington_sweet_velocity, a star-blind kernel, writes to standard output (`write(*,303)`, "RAT ... POT SPH,CYL,QUA") when the quadrupole correction ratio is negative. The README says "Kernels print nothing", GUIDELINES "Code style" rule 7 says the same, and the line bypasses rotation_messages, so it never reaches the run log; the kernel cannot be `pure`, unlike its four siblings. Stage D moved the shear's report out of its kernel (the dispatcher reports it) but left this one. The kernel test's pin carries two such lines with numbers its masking does not cover.
   - **In the F77?** Yes (vcirc.f:178-183).
   - **Fix:** return the interfaces with a negative ratio (as compute_shear_velocities returns `dynamical`), report them from compute_mechanism_velocity through rotation_note, make the kernel pure.
   - **Kind:** A (terminal text; the kernel test's pin moves). **Size:** S.
7. **rotation/mid_timestep_model.f90:432-434 and :55-56 (element_diffusion/element_diffusion_lib.f90:39-48)**
   - **Issue:** burn_and_settle declares `enclosed_mass` `intent(inout)` because "the settling routines convert it to Bahcall & Loeb units in place and back, which can change its last bit". Since O16's fix (c66901c9) they convert copies: microdiff and gravitational_settling take it `intent(in)` (microdiff_driver.f90:84, :123-126; original_settling.f90:84), and only diffuse_elements still declares it, the timestep, dlnP/dr and the total mass `inout`. The comment is false, and the loose intent puts `m` in advance_rotation_substep's `star writes:` line, which GUIDELINES' star rule 4 says to fix by tightening the intent. The `inout` was introduced on 2026-10-12 for the in-place conversion (ROADMAP, the rotation series).
   - **Fix:** `intent(in)` for those four dummies of diffuse_elements and for burn_and_settle's `enclosed_mass`, the comment gone, `check_star_contracts --fix`.
   - **Kind:** A. **Size:** S.
8. **rotation/evolve_rotation.f90:308, :311, :317-321; seculr/am_transport_coefficients.f90:118-148; seculr/am_grid.f90:264-268; seculr/circulation.f90:189-217, :260-261; seculr/instability_velocities.f90:116-131, :173-178; seculr/stability_coefficients.f90:188 (retired/README.md:5-6, README.md:48, retired/burn_mix_extrapolated.f90.txt:1-8, retired/burn_settle_mix.f90.txt:1-8)**
   - **Issue:** 86 lines of commented-out code remain in five live files, most of them the advection-diffusion branches (they still name `rot_scr%` and `star%`, which no longer exist), two blocks inside a star-blind kernel. The standard puts kept code in retired/*.f90.txt; retired/README.md says it was moved "so that those read as the code that runs", and README.md:48 that advection_diffusion.f90.txt holds "every fragment with the place it came from". It holds one evolve_rotation fragment (:690-692), not these. 45 more lines are "(commented out ...)" / "(now in rotation/retired/...)" markers (17 in stability_coefficients.f90). Two retired files have no header saying where they came from or how to revive them, only the "Modernized ... logic unchanged" boilerplate GUIDELINES bans, against README.md:48's "each file's header says how to revive it"; burn_li_be_midpoint.f90.txt's header sends a reviver to "its tables in rotation_scratch_lib", which are now in advection_diffusion.f90.txt:835-.
   - **Fix:** move the fragments under headers in advection_diffusion.f90.txt, one marker line per site at most; a header for each retired file.
   - **Kind:** A. **Size:** S.
9. **rotation/seculr/transport_checks.f90:134-141, :237-238, :427-430, :442-444; mid_timestep_model.f90:490-493, :538-548; magnetic_braking.f90:494-496; seculr/instability_velocities.f90:90-91; seculr/transport_mechanisms.f90:95-96**
   - **Issue:**
     - every message is F77 capitals naming routines that no longer exist: "ERROR IN SR CHECKJ", "ERROR IN SR CHECKC", "STOPPED IN SUBROUTINE MIX", "ERROR IN WIND - TAUCZ NOT DEFINED ... STOPPED" (no space before STOPPED), "MAX D(J/M)/(J/M) ... AT PT.", "RAT ... POT SPH,CYL,QUA";
     - the negative-j stop says "AND 3 ATTEMPTS AT CUTTING TIMESTEP FAILED" though no attempt is made, as the comment at :124-127 records, and counts a cut it does not make (:134);
     - a recoverable cut is labelled "ERROR IN SR CHECKC";
     - "LIMIT" in the dynamical-shear note prints an omega (the inner omega at the limit);
     - the settling lines 'JMAX=', 'FM= ', 'TSCALE=', 'DTMAX=' go straight to the run log's unit, outside rotation_messages (M112's class).
     io was marked down for controls_check's shouting messages.
   - **Fix:** plain messages in today's names and units ("check_angular_momentum: j <= 0 in zone 20; the run stops"), the settling lines through rotation_note.
   - **Kind:** A for every model; the two rotation test pins and solver-diagnostics run logs move. **Size:** S.
10. **rotation/README.md:14; rotation/rotation_scratch_lib.f90:34; shape/rotating_shape.f90:41; seculr/region_transport.f90:13; seculr/transport_mechanisms.f90:21; seculr/circulation.f90:68, :321, :401, :528 (and am_grid.f90:60, :342; equal_grid.f90:42-180; rotational_mixing.f90:66-390; stability_coefficients.f90:70); tools/check_boundaries.py:164-200**
    - **Issue:** README.md:14 says "the kernels that touch no scratch either are listed in STAR_BLIND_FILES". Three listed files use rotation_scratch_lib (rotating_shape's `use` is unused), and the scratch module takes `json` from star_info_lib, so those "star-blind" kernels compile against star_info_lib through it; check_star_blind looks only at the file's text. Five files that take no star (circulation, am_grid, equal_grid, rotational_mixing, stability_coefficients: 2,084 lines) are unlisted only because they import `json` (and the slots) from star_info_lib, not composition_layout as the listed kernels do, so nothing guards them from growing a star reference.
    - **Fix:** take `json` and the slots from composition_layout in those five files and in rotation_scratch_lib; drop rotating_shape's unused `use`; list the five; README's rule "every file that takes no star".
    - **Kind:** A. **Size:** S.
11. **rotation/evolve_rotation.f90:173, :260, :437-444; seculr/secular_transport.f90:638-645; core/model/starting_model.f90:262; rotation/rotation_configuration.f90:137**
    - **Issue:**
      - `substep_time_fraction` is a star member that only rotation writes and reads within one step: written before each use (:173, :260) and read by the wind_now functions, so starting_model's 0.5 is never read. A one-step value on the star is what star rule 5 ("scratch stays with its driver") and rule 2 rule out.
      - wind_now is written twice (evolve_rotation's internal function and secular_transport's module function), against "each formula once".
      - `rotation_settings` is public and writable (test_rotation assigns its fields), where net, kap, eos and atm declare theirs `protected`.
    - **Fix:** pass the fraction (or the wind_conditions record) into apply_secular_transport, one wind_now, drop the member; make the record `protected`, with the tests building a record and calling set_rotation_settings.
    - **Kind:** A. **Size:** S.
12. **rotation/magnetic_braking.f90:319-372 against :447-472 and :519-540; :517-518**
    - **Issue:** each torque law is written twice: compute_solid_body_braking forms the Matt and Kawaler torques, the saturation and the Rossby/B-field scalings itself, while brake_convection_zone uses set_matt_law, set_kawaler_law and domega. The solar centrifugal term fsun is two identical lines (:353-354, :468-469). O17's missing 1/G was a divergence between exactly these copies, and a braking law of a new form must be written into both. The two forms are equal up to round-off (min(omega tau/tau_sun, w_sat) against min(omega tau/tau_sun, w_sat tau_sun/tau tau/tau_sun)), so merging them is a measured commit. domega's header still says tau_cz "is zero without Rossby scaling, as in the original", untrue since O17 (set_matt_law sets it with B-field scaling alone, :454-455).
    - **Fix:** compute_solid_body_braking calls domega at omega_old and at omega_surface with the law set_matt_law or set_kawaler_law builds; README's braking-law recipe then names one function per law.
    - **Kind:** B (round-off). **Size:** M.
13. **rotation/README.md:111-116 ("Add a transport mechanism")**
    - **Issue:** step 4, "Add its two controls to the registry (defaults/controls_registry.tsv, then the generators)", leaves out the controls' declarations in io/controls_buffer.f90, the only place a control's default and description live (io/README.md §5 step 1; gen_inlist refuses a row without one). Unlike the "Add a control" recipe above it, it does not point to io §5. It also omits the steps a newcomer needs to finish: `check_star_contracts --fix` for the builder's `star reads:` line, test_controls' pristine.txt re-pin, and the kernel test's pin (expected_test_rotation_kernels.out) gaining the new checks' lines. Traced against Stage E's addition (73ca7732): its hand edits were the layout, instability_velocities, transport_mechanisms, circulation's comments, controls_buffer, the registry, the builder, both test programs, stitched_model and profile_columns.list; the recipe names all but controls_buffer.f90.
    - **Fix:** step 4 "as io/README.md §5 says (the declaration, the row with `settings` = `rotation_settings%velocity_scale(i_mech_<name>)`, the builder's copy, the contracts and the re-pins)", and re-trace the recipe with a dummy mechanism.
    - **Kind:** A (text). **Size:** S.
14. **Stale comments and doc statements**
    - **Issue:**
      - seculr/transport_checks.f90:101-102: "max_diffusion_iters is not checked against it", where its own :17-19 and controls_check.f90:91-98 check it (O14).
      - rotation_configuration.f90:148-149: the controls read sets `solid_body_mode_flag = 1`; the builder does (rotation_settings_builder.f90:100). Its :66-70, :76-79 and :86-87 list the mechanisms' controls and D as three mechanisms, written before Solberg-Hoiland.
      - seculr/equal_grid.f90:18-22 describes the advection branch's interpolation as a block that "stays with" the callers; it is retired.
      - README.md:
        - :67 says rotation calls core for "two services". It uses three core modules and eight procedures: shell_physics; burn_step's store_end_li_be_rates, burn_li_be, network_rates_of_model, burn_regions, burn_network_region and surface_accretion_now; convection_geometry's pressure_scale_height.
        - :79-83's list for test_rotation omits its last case, the secular transport's stability factors, velocities and diffusion coefficients (test_rotation.f90:16-17).
        - :50 starts a sentence "its banner".
        - :105 dates the coverage (96.9% of 2168 lines) at 2026-10-29. It was measured at 2c88d896 (10-29 00:33), before c7d8beff and d4bcaf27 added the osplin error paths at nine calls and the polint report, which no deck reaches. docs/ROTATION_VALIDATION.md §4 carries the same figure.
      - docs/TARGET_ARCHITECTURE.md:
        - :507 names rotation/burn_settle_mix.f90 and burn_mix_extrapolated.f90 as files (now retired/*.f90.txt).
        - :542 puts solve_angular_momentum_transport in secular_transport.f90; it is in am_transport_coefficients.f90.
        - Phase 4 (:456-466) still lists the mechanism list and "one name per quantity" as work (done in Stages C-D), and says rotation has no truth-level test, beside ROTATION_VALIDATION's MESA comparison.
      - docs/HOW_YREC_WORKS.md:554-556: the sub-steps are a "dtdif fraction of the diffusion time". They are the step's largest omega change divided by dtdif (evolve_rotation.f90:276-286).
    - **Fix:** text only.
    - **Kind:** A. **Size:** S.
15. **The declaration and style rules (all of rotation/ but magnetic_braking.f90, rotation_configuration.f90, rotation_messages.f90 and instability_velocities.f90's dummies)**
    - **Issue:**
      - Against GUIDELINES "Declarations and interfaces": 344 `double precision` statements (real(dp) 67); 245 explicit-shape array dummies sized by `json`, `15`, `13`, `12` against 62 assumed-shape; one assumed-size `a(*)` (circulation.f90:275, acknowledged in ROADMAP); 22 whole-module `use`s besides math_lib (physical_constants 17, equal_grid 2, mid_timestep_model, secular_transport, rotation_scratch_lib); mid_timestep_model public by default (:22); 10 of about 118 procedures `pure`, though zone_moments_of_inertia, shape, solid_body_omega and most of seculr's helpers have no side effects.
      - Against "Code style": 117 `.lt.`-style operators; constants as run-time variables or `data` (rotating_shape.f90:81-86 and :211-214, where eps "cannot change" the result by its own comment; viscosity.f90:38-43; moments_of_inertia.f90:110-112; transport_checks.f90:302-303, whose He4 weight 4.002603 differs from viscosity's 4.0026); unnamed thresholds (evolve_rotation.f90:423 1e-24, mid_timestep_model.f90:141 1e-14, rotational_mixing.f90:206 1e-14, viscosity.f90:91 4.2e5 and :112 1e-38, magnetic_braking.f90:489-490's two 1.0 tests) and an inline coefficient (mid_timestep_model.f90:530, 4.348e21, sourced in the comment); slots by number (mid_timestep_model.f90:483, :508, :595-598; rotational_mixing.f90:204 `j_idx.ne.2`).
      - 37 internal procedures share their hosts' locals (evolve_angular_momentum 8, circulation_velocities 8, advance_rotation_substep 6, am_transport_grid 5, check_angular_momentum 4, rotation_stability_setup 3, burn_and_settle 2, shape 1).
    - **Fix:** the byte-identical conversions io and net had (`real(dp)`, operators, `end if`, standard intrinsics, private by default, only-lists, assumed shape, `pure`, parameters, named slots), file by file with objects compared. Then the host-associated steps as module procedures with explicit arguments, as Stage D Phase 3 did for apply_secular_transport.
    - **Kind:** A. **Size:** M (the conversions), L (the steps).

## What would raise each score

- **Readability, to 9:**
  - finding 15's conversions;
  - findings 8 and 9;
  - evolve_angular_momentum's header history into a Sources note, and its stale local comment (evolve_rotation.f90:101-107);
  - the names transport_checks.f90:308-310 apologizes for, renamed (mu_ion, mu_e), and `shape` renamed (compute_level_surfaces);
  - the 42 lines over 110 characters wrapped.
- **Organization, to 9:**
  - findings 10 and 11;
  - the host-associated steps of evolve_angular_momentum, circulation_velocities and advance_rotation_substep as module procedures with explicit arguments;
  - ROADMAP Stage D's items 2-3: a core service "advance the composition over an interval on a given model" (ending rotation's use of burn_step), and records instead of 14-24 arguments.
- **Maintainability, to 9:**
  - findings 1-7, and O29 and O10's gate;
  - finding 12;
  - ok/FAIL checks, each confirmed by a mutation, for:
    - viscos (the Spitzer and Ledoux limits);
    - the Zahn coupling factor;
    - diffusion_velocity_scales (slot order, the background term only in D_omega, the core's zeroing);
    - enforce_rotation_profile's modes 2 and 3 and the general law (J conserved, omega proportional to r^alpha);
    - disk locking in solid_body_omega;
  - the coverage re-measured.
- **Extensibility, to 9.5:**
  - finding 12, so that a braking law is one function and its coefficients;
  - finding 13, with the mechanism recipe re-traced;
  - a traced "add a braking law of a new form" recipe;
  - assumed-shape dummies (finding 15), so that rotation stops naming the zone and species capacities.
- **Docs, to 9:**
  - findings 13 and 14;
  - README.md:14, :48 and :74 made true (findings 6, 8, 10);
  - GLOSSARY's 25 dead rotation paths (ROADMAP still-open 7).

## Traces, counts and checks (supporting)

- **Add a control rotation reads** (README.md:109, via io/README.md §5): 4 hand edits, all named between the two READMEs.
  - The edits: the declaration in io/controls_buffer.f90, the registry row with `settings` = `rotation_settings%<field>`, the field in rotation_settings_t, the copy in rotation_settings_from_star.
  - Then make, `check_star_contracts --fix`, check_boundaries (which prints the `settings` cell) and test_controls' re-pin.
  - io traced the case on 2026-11-01 with "a physics scale rotation's settings read".
- **Add a transport mechanism** (README.md:111-116): traced against 73ca7732 (finding 13).
  - The layout's slot, name, correction and control names: one edit of state/transport_mechanism_layout.f90.
  - The kernel and its check, and the dispatcher case. Every bookkeeping loop already loops over the slots (circulation.f90:110-266; mechanism_sum, :289-296).
  - The controls' copy into velocity_scale/mixing_scale, which test_rotation_kernels.py's test_mechanism_controls checks against the layout and the registry.
  - The profile column in core/diagnostics/stitched_model.f90 (:572-585) and defaults/profile_columns.list.
  - Missing from the recipe: io/controls_buffer.f90 and the re-pins.
- **Counts, live code** (24 files, 7,607 lines):
  - 344 `double precision`, 67 `real(dp)`;
  - 117 old relational operators, 60 `endif`/`enddo`;
  - 245 explicit-shape against 62 assumed-shape array dummies;
  - about 118 procedures, 10 `pure`; 37 internal procedures in 8 hosts;
  - 22 procedures with 14 or more arguments (diffuse_composition_driver 24);
  - 86 lines of commented-out code and 45 marker lines; 42 lines over 110 characters;
  - 8 message formats, all in capitals.
- **Tests:**
  - test_rotation_kernels: 92 ok/FAIL checks in its pin (114 lines). Its pytest compares a check's label and status only, but the pin also carries unmasked numbers from the kernel's RAT print (finding 6) and the dynamical-shear notes.
  - test_rotation: 358 pinned lines, no ok/FAIL check.
  - test_rotation_kernels.py's test_mechanism_controls: layout to registry to builder.
  - CI runs both (.github/workflows/ci.yml:129).
- **Checked and true:**
  - README §1's 24 rows against the files and modules;
  - §2's step order against evolve_rotation.f90:127-216;
  - the seven public procedures (check_boundaries PUBLIC);
  - the settings install (setup/settings_install.f90:100-112, from run_driver.f90:346);
  - the controls-free and cgl checks; SCRATCH_TYPES_ONLY and SCRATCH_OWNERS;
  - the star contracts (the checker passes, 149 declared);
  - §3's slot ranges and the mu-twice note; §4's kernel-check list against test_rotation_kernels.f90's header;
  - the builder's copies (every rotation_settings_t field set, check_settings_fields);
  - ROADMAP Stage D's still-open list; bug_sweep's statuses for O11-O28 (fixed) and O29 (open).
- **Not rechecked:** the mutation claims, the coverage figure, the MESA comparison and the 56-deck byte comparisons (nothing was run).

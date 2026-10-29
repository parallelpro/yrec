# state/ + setup/ re-review (2026-10-29, yrec-next 7aa138b2, read-only)
## state/ + setup/: scores now (1-5, relative to rotation/)

**state/**
- Readability: 2.5 -> 3.5. Every star_info member now has its own comment in a fixed form (meaning, unit, writers). The 62 renames removed most of the cryptic names. What holds it back:
  - the hand-declared half of star_job has no comments;
  - misleading names are left: `i_eps_he3` is the triple-α channel; `fit_point_*` hold log10 values; `h_shell_*_radius` hold r/R; `qs0e_scale` and `qqs0ee_scale` remain;
  - `physical_constants` is neither private nor real(dp).
- Organization: 2 -> 3. The members are grouped into lifetime classes, controls_sync is out of state/, integration_parameters is gone, and the constants are parameters. What remains:
  - star_job still mixes read-only paths, a scratch array and run state with the controls the run changes;
  - constants and copies still live on the star;
  - the solar units are handled three different ways;
  - cgl is a global;
  - 97 of the 327 ctrl components are dead.
- Maintainability: 2.5 -> 3.5. The README tables are generated, and the --check catches a missing comment, a false writer or a dead member. There is one default source, and the job homes are checked. Unchecked:
  - comment prose and class placement (a stale eps_total passed the check);
  - the hand-declared half of star_job;
  - only 29 of 146 star procedures declare contracts.
- Extensibility: 2.5 -> 3.5. The "add a star field" recipe is accurate and the generator guides each step. Where a member goes, and when it is re-armed, is left to judgment. Contracts are optional.
- Docs accuracy: 2.5 -> 3.5. The writer and reader tables are right by construction (`gen_star_docs --check` passes at HEAD). About 8 prose errors came in with the merges and the renames (new issue 2).

**setup/**
- Readability: 4 -> 4. The builders are flat copies with contracts. user_inputs.f90 is unchanged: DATA statements, `.gt.`, a bare `use`.
- Organization: 3 -> 3.5. Net's builder moved in (pure, installed with ierr). Still open:
  - map_user_inputs stays here;
  - there are 2-3 naming schemes;
  - settings are installed from 4 places;
  - kap's and atm's table paths bypass setup;
  - the rotation builder imports io's controls_buffer.
- Maintainability: 3.5 -> 3.5. Contracts and 95% coverage, but no test that each field is copied from its control. Derivations are creeping into builders (EOS vintage, walpcz clamp, log10 G).
- Extensibility: 3 -> 3.5. There are now recipes to add a field, a domain and a control, and the checker and Makefile names they cite are correct. There is no install list and no check that a domain is installed, and the recipes were not traced.
- Docs accuracy: 3 -> 4. The rewritten README is mostly true. Small errors: the star-member list, "computes nothing" (three builders do), "the only place its default is written", "Two install points".

## Summary
The overnight work fixed most of the first report:
- one default source;
- job is exactly the 36 controls the run mutates, enforced by a check;
- star_info in 13 lifetime classes with one comment per member;
- the README tables generated and checked;
- constants as parameters;
- integration_parameters and controls_sync gone;
- the settling latch moved onto the star;
- net's builder moved into setup/.

What is left, and what the merges introduced:
- The hand-declared half of star_job (26 table paths, mixture_weights, the Monte-Carlo samples, nk) is outside the job rule, its check, the generator and the registry. mixture_weights is really dead scratch.
- About 8 doc errors came in with the merges and renames: eps_total after M1, `dlnrho_dlnt` listed as removed, damaged Sources.
- cgl is still a global read by star-blind kernels in core and rotation, and it is 0 until run_setup.
- 97 dead ctrl components are still generated.
- The setup naming and install plan is untouched (ROADMAP still-open item 4).

I found no new results-changing bug. Every possible bug below is latent.

## First report's findings: status
- **F1, two default sources (107 controls)**: FIXED (6085de4c). controls_state_def.inc carries the values that take effect: :111 settling_tolerance 1e-8, :124 atm_error_tol 3e-4, :176 es_velocity_scale 1, :267 scale_by_rossby_number .true., :269 'K97', :307 s0_pp 4.01e-22. set_legacy_defaults (io/legacy_controls.f90:345-438) now sets only the legacy-only names.
- **F2, the ctrl/job split mislabelled**: PARTLY.
  - Done: all 36 members of job_controls_def.inc:7-42 are assigned outside the read (I checked each, e.g. stop_conditions.f90:301-303, step_driver.f90:396/:483, mass_change.f90:126). ctrl is assigned only in io/controls_sync.f90. check_boundaries.py:315-356 enforces it.
  - Left: the hand-declared star_job members (star_info_lib.f90:115-139) are outside the rule and the check (new issue 1).
- **F3, star_info in COMMON order**: FIXED. Classes with banners at star_info_lib.f90:146-1018, history in Sources at :1031-1155. A few placements are arguable (new issue 5).
- **F4, README drifted, script missing**: FIXED. gen_star_docs.py, run with --check by test_generated_files.py:16-17 (passes at HEAD). Module table complete (README:5-14). star_template described (:25).
- **F5, run-time "constants"**: PARTLY.
  - Done: 8 parameters with sources (physical_constants.f90:18-48); cmkhn and ca3l gone; test_eos and test_convection boot no star.
  - Left: cgl is still a run-time global (physical_constants.f90:52, set at run_setup.f90:85) read by about 23 files, and only convection gets it through its record. The version strings (physical_constants.f90:92-93) are also run-time variables in this module.
- **F6, integration_parameters**: FIXED. Deleted (97c00d5a); the three values are ctrl members (controls_state_def.inc:133-135) copied into the kernel's cfg.
- **F7, controls_sync in state/**: FIXED (8398d475). Leftover: check_boundaries.py:310 still lists `state/controls_sync.f90`.
- **F8, constants and derived controls on the star**: PARTLY.
  - Done: the calibration's four derivatives are named parameters (calibration.f90:20).
  - Left: they are still copied into s%dlum_dx..drad_dalpha (calibration.f90:87-90; M34). Also unchanged: solar_bolometric_magnitude = 4.79 (run_setup.f90:83), hsra_log_t_offset (:88), use_two_z_tables (run_driver.f90:123), envelope_store_log_temperature (derived_controls.f90:43).
- **F9, duplicated quantities**: PARTLY.
  - Done: stotal merged into log_total_mass; jcz renamed.
  - Left: the solar_mass_cgs copy (run_setup.f90:78); star_mass/log_total_mass is still a pair kept in step (mass_change.f90:537/540); timestep_yr and the cp/mu/dlnrho pulse pairs were kept by decision.
- **F10, earlier-value conventions and cryptic names**: PARTLY.
  - Done: the 62 renames introduced `_start`, `_prev` and `_iteration_start`.
  - Left: `_prev` has four meanings (README:380); there is a lone `_end`; `envelope_cz_bottom_index_start` sits in the Regions class, not Step snapshot; still cryptic: qs0e_scale, qqs0ee_scale, be7p_*, o16_gamma_scale, c12_alpha_scale.
- **F11, unit-suffix rule**: PARTLY. 15 members renamed. Left: h_shell_*_radius (r/R), fit_point_* (log10), neutrino_flux_total/zone (1e10 cm⁻² s⁻¹).
- **F12, input-only ctrl members**: PARTLY.
  - Done: zsi and alex_mixture_x/z deleted, mu_gradient_scale retired (465bfc60); the registry's `used` column is checked.
  - Left: 97 components (30 marked `none`, 67 `read`) are still generated into s%ctrl.
- **F13, setup/ organization**: PARTLY.
  - Done: net's builder is setup/net_settings_setup.f90, installed at run_driver.f90:136.
  - Left: map_user_inputs is still in setup/; two naming schemes (setup/README.md:52); user_inputs.f90 still has DATA statements (:45-90), `.gt.` (:179-249), a bare `use controls_buffer` (:34) and a doubled banner (:14-15).
- **F14, setup README stale**: FIXED (rewritten, with recipes in §5). New small errors are listed under new issue 10.
- **F15, star_info_lib stale comments**: FIXED (all cited lines are gone). New staleness is listed under new issue 2.
- **Bug 1, settling latch**: FIXED (af946c0c). Armed at run_driver.f90:331, cleared at composition_step.f90:264 and mid_timestep_model.f90:467. The sibling use_mass_accretion latch is M31 in bug_sweep (open, mass_change.f90:126).
- **Bug 2, declared defaults not the real ones**: FIXED (with F1).
- **Bug 3, test with an unset csigl**: FIXED (d6736788); csigl is now a parameter.
- **Bug 4, controls stored and never read**: FIXED (465bfc60). M32's three remain, documented as `used = none`.
- **Bug 5, dead constants**: FIXED. cmkhn and ca3l appear nowhere.

## New issues (introduced overnight or newly found), ranked
Class A means the fix leaves results unchanged; class B means it changes results.

1. **star_info_lib.f90:106-139: star_job's hand-declared half breaks the job rule and escapes every check.**
   - The 26 table paths are assigned only by the read (controls_reader.f90:136-160, derived_controls.f90:75-85).
   - `mixture_weights` is not a control:
     - kap_init writes the LAOL header mixture into it (run_setup.f90:92 → kap_lib.f90:425, an intent(inout) dummy);
     - update_surface_mixture overwrites it (starting_model.f90:994-1027) before anything reads it, so that write is dead (the F77 `LLAOL` guard was never set; ce886f00).
   - None of these members has a comment, a README row or a registry home.
   - Two docs are wrong: star_info_lib.f90:107 says the mixture weights are "set by the read", and state/README:3 says job "holds exactly the controls that code outside the read assigns ... checks".
   - Fix:
     - paths into ctrl, or path records built in setup/ as eos's is;
     - mixture_weights becomes a local, and kap's laol array a kap local;
     - the Monte-Carlo samples and nk become documented star members.
   - A, M.
2. **Doc errors from the merges and renames** (gen_star_docs checks writers only, so all of these pass):
   - eps_total is described as "pp1+pp2+pp3+CNO+3 alpha" (star_info_lib.f90:839-840, README:300), but since M1 (ee268284) it includes the α captures.
   - README §4 (:372) says `dlnrho_dlnt` was removed and became a local, but the rename qdt→dlnrho_dlnt made it a live member (README:273). The local of the same name in compute_henyey_coefficients (structure_equations.f90:71) now sits beside `s%dlnrho_dlnt`.
   - The Step-working-state banner says these members are "not output" (star_info_lib.f90:688-689), but history_output prints convective_turnover_timescale and photosphere_logP (history_output.f90:287/:294).
   - stitched_model_state.f90:16 points to "io/yrec_output's registry"; the column table is now stitched_model's profile_columns (stitched_model.f90:300).
   - The rename script damaged the Sources section: :1068 now reads "zahn_coupling (zahn_coupling)" (the F77 name VFC is lost); at :1074-1075 the legacy flags IOV1/IOV2 were replaced by the new names; :1104 reads "controls.!" (a missing line break).
   - The README units list leaves out `_uhz`.
   - Fix: correct the text; add a gen_star_docs check that §4 names no live member.
   - A, S.
3. **physical_constants.f90:52: cgl has no initializer, so it is 0 (G = 1) until run_setup.f90:85 sets it.** It is read by star-blind kernels:
   - envelope_integration_kernel.f90:620, envelope_integrands.f90:244/:253, shell_equations.f90:215;
   - rotation's shape/, stability_coefficients.f90:138, magnetic_braking.f90:353/:467/:572.

   This is bug 3's trap again. Rotation's tests set the literal -7.17571 by hand (test_rotation.f90:59, test_rotation_kernels.f90:203). It also contradicts GUIDELINES Code style 3 ("a domain takes it from its settings record"). Fix: initialize `cgl = log10_newton_g_default` (results unchanged), then add log10 G fields to rotation_settings and the kernel cfgs. A, S/M.
4. **controls_state_def.inc: 97 of 327 ctrl components are dead.** The registry marks 30 `none` (22 unit numbers, fo, ies, imu, the canonical fmonte copies, M32's three) and 67 `read`; gen_controls_state emits them all. Fix: generate controls_state from `used = run` rows only, and leave the rest in the buffer. A, S-M.
5. **Leftovers of the naming and unit rules:**
   - h_shell_{bot,mid,top}_radius are r/R (observables.f90:459-467; star_info_lib.f90:956-961), while their `_msun` siblings were renamed.
   - fit_point_{pressure,temperature,radius} hold log10 values (:426-434).
   - neutrino_flux_* have no unit suffix.
   - i_eps_he3 is the triple-α slot (:76; structure_equations.f90:372), beside i_lum_3alpha.
   - The star's qs0e/qqs0ee/o16_gamma/c12_alpha names differ from net's record names for the same values (net_settings_setup.f90:36-42).
   - envelope_cz_bottom_index_start is in the Regions class, not Step snapshot.
   - Fix: a rename script, as a568e084 was. A, S.
6. **The renames stopped at the star's edge.** The same quantity keeps its old name one step away:
   - cfg%senv/tenv/atm_hras (envelope_integration.f90:155-159);
   - eos_mix%amuenv/fxenv (eos_lib.f90:221-222);
   - sub_physics%visc/thdif/qdt (mid_timestep_model.f90:186-189);
   - the dummy arguments h_shell_zone_begin and h_shell_midpoint_zone (timestep_control.f90:45-72, rezoning.f90:74).

   Comments still name the older `pressure_scale_height_start/_end` (step_driver.f90:360, convection_geometry.f90:133/:155). The test label 'star%tenv' remains (io/test/test_controls.f90:78), and controls_registry.md:263 still says "senv". A, S.
7. **check_star_contracts.py:140 keys procedures by name only (`procs.setdefault`).** A second procedure with the same name is never checked, and a callee's intents resolve to the first one found. There are two `wind_now` (evolve_rotation.f90:426, secular_transport.f90:611) and two `write_gsm_pulse`. gen_star_docs keys by (module, name), and check_control_homes reuses this footprint. It also counts `read(...) s%x` as a read where gen_star_docs counts a write. Fix: key by (module, name). A, S.
8. **microphysics_settings.f90:86-99: equation_of_state_settings derives the OPAL vintage, against the "pure copies" rule (setup/README.md:47).** Its conflict branch cannot be reached from a deck, because derive_controls already lets the newest vintage win, without a message (derived_controls.f90:76-86). There are two policies for one condition. Fix: check_controls rejects a deck that selects two vintages (B for such decks), or derive_controls sets one vintage control that the builder copies (A). S.
9. **The solar units are handled three ways.**
   - Msun is a control plus a star copy (run_setup.f90:78).
   - Lsun and Rsun are star members written by the read from legacy-only clsun/crsun (controls_reader.f90:248-249).
   - Teff_sun, nu_max_sun and delta_nu_sun are ctrl.

   Fix: make clsun and crsun canonical controls and drop the copy. A, M.
10. **setup's install and path plan.**
    - kap_init and atm_init take their table paths positionally from s%job (run_setup.f90:91-97, :107-108), while eos's paths go through setup's equation_of_state_table_paths.
    - Settings are installed from four places (run_setup, run_driver.f90:136 and :321-324, starting_model.f90:1014), and nothing checks that a domain is installed.
    - rotation_settings_setup.f90:56 imports io's controls_buffer for the ichi_* slots.
    - README errors: "solar_radius_cgs, solar_mass_cgs (convection, rotation)" (:21; convection reads neither); "computes nothing itself" (three builders derive); "the only place its default is written" (:87) against the domain-type initializers of step 1 (:70).

    A, M.
11. **Dead or duplicate star members kept:**
    - be7_mass_fraction_zone and neutrino_flux_zone are used inside one procedure;
    - zahn_coupling belongs in rotation's scratch record;
    - envelope_specific_entropy, delta_log_pressure and delta_log_temperature are always zero (mass_change.f90:118);
    - star_mass and log_total_mass are kept equal by the same writers.

    A, S.
12. **The checkers check less than the docs claim.** check_control_homes covers only registry rows that have a home (check_boundaries.py:318), and lists a deleted file and setup/user_inputs.f90 as part of the read (:309-310). The class banners (Scratch, the `_start` class) are not checked against the code. A, S.

## Domains still reading s%ctrl / s%job directly (code only, comments stripped)
- **Domains:** net, kap, eos, atm, numerics, math, convection, element_diffusion, rotation: 0. kap has 2 comment mentions (kurucz_opacity.f90:37, opal92_opacity.f90:41); rotation is enforced at 0 by CONTROLS_FREE_DIRS.
- **Configuration that bypasses the records:**
  - cgl in 8 rotation files and 4 core kernels;
  - net's burn tolerances passed per call from s%ctrl (burn_step.f90:515-517);
  - kap's and atm's paths passed positionally;
  - convection's metals through a second setter (starting_model.f90:1014).
- **Star layer (allowed):** calibration 173, starting_model 69, monte_carlo 65, model_rescaling 55, run_driver 47, step_driver 39, timestep_control 37, rezoning 33, stop_conditions 30, henyey_iteration 22, burn_step 18, envelope_integration 18, run_setup 16, mass_change 13, surface_boundary 10, composition_step 10, envelope_rebuild 9, stitched_model 8, structure_equations 7, observables 7, shell_boundaries 2, convection_geometry 1, neutrino_flux_output 1. Also io/ 869 and setup/ 157.
- **Writes:** s%ctrl only in io/controls_sync.f90 (320). s%job in 11 core files plus the read.

## Possible bugs (not in bug_sweep)
1. **load_table_paths** (io/controls_reader.f90:105-131). Its comment says a deck that names no table "keeps the previous read's path". That cannot happen: every run starts from star_template (yrec_reset.f90:47). Instead it copies the template's never-initialized job paths (NUL bytes) over set_legacy_defaults' blanks, so an omitted path holds 256 NULs. Latent: I found no code that tests a path for blank.
2. **job%mixture_weights:** kap_init's LAOL-header write is dead (new issue 1). Not a results bug: the F77 never used it either.
3. **cgl = 0 before run_setup** (new issue 3): any new kernel test of rotation's shape or the envelope/shell kernels silently runs with G = 1.
4. **A deck selecting two OPAL vintages** gets the newest without a word, and eos_init's conflict error never fires in a run (new issue 8).
5. **Check, not confirmed:** since M1, eps_channels rows 1-5 (the profile's eps_ppI..eps_he3) no longer sum to 1 where α captures burn, because no channel maps the α captures (structure_equations.f90:355-377). Decide whether that is intended, then document it.

## What would raise each score by one point
**state/**
- Readability: comment every star_job member (or move them out per new issue 1); run the new-issue-5 renames; make physical_constants private, real(dp) and without the version strings.
- Organization: split star_job per new issue 1; take the constants off the star (solar_bolometric_magnitude, hsra_log_t_offset, dlum_*); handle the solar units one way; generate ctrl from `used = run` only; move cgl into the records.
- Maintainability: extend gen_star_docs to star_job, to class rules (a single-module member must be Scratch; `_start`/`_prev` members in the snapshot class) and to §4 against live names; key check_star_contracts by (module, name); add contracts to every core/run and core/step driver.
- Extensibility: trace "add a star field" on a throwaway branch, as io's "add a control" was; add a rule for per-card re-arming of state.
- Docs: fix the new-issue-2 list.

**setup/**
- Readability: bring user_inputs.f90 to the code style (parameters with sources, `<`, only-lists, one banner), or move it to io/.
- Organization:
  - map_user_inputs to io/ beside derive_controls;
  - one builder name, `<dom>_settings_from_star`;
  - one install module (settings_install.f90: per run, per Monte-Carlo run, per card), with check_boundaries allowing `set_*` only there;
  - kap's and atm's paths as records built in setup/;
  - the ichi_*/itime_* slots into a state/ constants module.
- Maintainability: a static or mutation check that every settings field is assigned from its registry control; move the OPAL-vintage derivation to derive_controls so every builder is a pure copy.
- Extensibility: with the install module in place, trace "add a field" and "add a domain" on a throwaway branch.
- Docs: fix the new-issue-10 README errors and describe the install sites as they are.

All work was read-only. Besides grep and git, I ran `gen_star_docs.py --check` and `check_star_contracts.py` in check mode with `python3 -B`; both pass and the tree is unchanged.

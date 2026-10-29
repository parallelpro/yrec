# io/ review (2026-10-29, read-only evaluation; paths under src/)

Scores vs rotation/: Readability 2.5, Organization 2.5, Maintainability 2, Extensibility 2, Docs accuracy 3.

The controls reader is in good shape at its top level (read_input is a sequence of steps, the generators
are checked, test_controls pins all 313 decks). Underneath, an inlist value goes through three generated
copies (legacy -> new -> legacy -> canonical); controls_buffer mixes real controls with 29 unit numbers,
~10 derived wind-law coefficients and run-time latches; ~600 default-public module variables (~370 in
controls_buffer, ~230 in legacy_controls); the ctrl/job split is a Python list; column tables are parallel
integer-indexed arrays; output state in module variables. README "add a control" step 4 contradicts the
settings-record convention; "add a column" points to the wrong file for profiles.

## Possible bugs
1. read_controls.f90:249-482, controls_buffer.f90:217 -- with LNEWVARS / use_new_rate_variables = .false.
   (the default, controls_buffer.f90:481), chi_grid_scale stays zero: no initializer, and nothing copies
   hpttol (declared legacy_controls.f90:73, defaulted :412) into it; only map_user_inputs fills it, under
   lnewvars (setup/user_inputs.f90:116). VERIFIED by the integrator (every shipped deck sets it .TRUE.:
   19 LNEWVARS + 520 use_new_rate_variables lines). In F77, HPTTOL was the working vector (parmin.f:406-407)
   overwritten only under LNEWVARS (remap.f:196); the copy was dropped in a576ec06 ("lnewvars is true in
   every current example"), the N4 pattern. Effect: zero zoning tolerances, a division by zero at
   core/model/starting_model.f90:544-547 when the core is extended. The registry marks hpttol "unused".
   Fix: `if (.not. lnewvars) chi_grid_scale = hpttol` next to the N4 line (:478) and make hpttol live in the
   registry. Byte-identical for shipped decks.
2. run_cards.f90:64-66 -- format 370 labels end_dcen "central: log rho ="; it is the central deuterium
   mass-fraction stop (stop_conditions.f90:53); F77 printed "CENTRAL D" (parmin.f:1323). A 2026 regression
   in every run log of a card with a stop; the registry row (defaults/controls_registry.tsv:93,
   "mesa: log_center_density_limit") repeats it. Class B (pinned logs change).
3. run_cards.f90:96-118 (inherited) -- a kind-3 card (rescale and evolve) never sets end_age_stop_active or
   timestep_override_active, so max_age, timestep_override and central_*_stop on that card are ignored.
   ~200 decks use run_kind(1)=3; none sets a stop on that card. Lines 110-111 print "rescale the previous
   card's model" for a rescale-and-evolve card (the right format was dead in F77 too; commented-out at
   :112-113). Derive the flags for kinds 1 and 3, or warn.
4. model_files.f90:49 -- write_mod_model rewinds unconditionally; Monte Carlo writes each run's model through
   it (core/run/monte_carlo.f90:178-179), so the Monte Carlo model file keeps only the last run's model
   (F77 WRTLST rewound only ILAST, wrtlst.f:141). VERIFIED (rewind + per-run call) by the integrator.
   Came in 0066e7af. The pinned cm_monte_carlo deck does a single run. Class B.
5. run_files.f90:63-80 -- expand_paths expands 15 path controls; never expanded: fopale (OPAL95 EOS, used at
   derived_controls.f90:109), fopale01, fkur, fkur2, flaol, flaol2, fopal2, flldat, fmhd1-8, opecalex, fiso,
   fdyn, fmonte1/2 -- a {YREC_INPUT} there reaches the open literally. Expand every path control from one
   list (A for shipped decks).
6. yrec_output.f90:86-88 vs 213-227 -- the startup HDF5 check tests trim(pulse_format) == 'GSM' but the
   dispatcher also accepts 'gsm': pulse_format='gsm' on a non-HDF5 build passes startup, then warns at every
   profile and writes nothing. Any unknown format silently means GYRE.
7. derived_controls.f90:36-40 -- the semiconvection-plus-overshoot error goes to the run log only, as "ERROR
   IN SUBROUTINE PARMIN" in capitals; the two errors beside it (:47-64) go to the terminal too.
   test_error_exit.py:62 pins the text ("SEMI-CONVECTION"). B (failure-path text).
8. model_files.f90:318-379, 543-604 -- read_yrec7 and read_model2 never check s%nz against json
   (read_mod_model does, :188): an oversized header writes past the star's arrays. No iostat on their data
   reads; read_model2 has no ierr; read_mod_model's reads at :133, :154-174 lack iostat. A malformed model
   ends in a Fortran runtime error.
9. run_files.f90:95-106 -- no '/' in log_output_file runs system('mkdir -p ') with no operand; slice arithmetic
   at :103 drops the trailing '/' and overruns shell_cmd for a directory > ~247 chars; path unquoted;
   yrec_output.f90:72-77 computes the same directory a second way. (Prefer execute_command_line or no shell.)
10. run_cards.f90:119-124 (inherited, ZENV at parmin.f:1378) -- at read time the loop over all cards sets
   s%envelope_metal_fraction and calls eos_set_mixture, so the last card's positive rescale Z wins until
   starting_model sets the mixture again; io reaching into eos state at read time. Verify before changing.

## Findings (ranked)
1. derived_controls.f90:121-185: the wind law is computed in io -- Kawaler/PMM exponents, constfactor, c_2,
   the K97/V13 tables, the walpcz clamp and the solid-body flag are pseudo-controls (controls_buffer.f90:446-456)
   copied into rotation_settings (setup/rotation_settings_setup.f90:82-83). Move to a pure coefficient
   function in rotation/magnetic_braking called by the settings builder. A, M.
2. README.md:93-100, :8 "add a control" wrong/incomplete: step 4 ("read it from star%ctrl") contradicts
   CONTROLS_FREE_DIRS/STAR_FREE_DIRS (a domain needs its settings record and builder); omits that a
   new-format-only control needs legacy = proposed_name; which status value; that a run-mutated control must
   go into JOB_MEMBERS (tools/gen_controls_state.py:44-72); where validity checks go. A, S.
3. README.md:102 "add an output column" wrong for profiles: profile values live in
   core/diagnostics/stitched_model.f90 (n_prof_cols :38, profile_value :452, envelope and atmosphere cases
   :312-376 where a missing case silently gives 0); a history column needs n_hist_cols bumped, two parallel
   routines edited and is_int_hist_col checked; every column must be in defaults/*_columns.list (commented
   out if opt-in) or test_default_columns_lists_in_sync fails.
4. Inert controls presented as live: use_rotation_scaled_solar_wind = .TRUE. in 15 decks
   (examples/giant_differential_rotation/...) while its code is commented out (core/step/mass_change.f90:383-392;
   the registry says "copied to star%rot%..."); use_pmm_wind (lpmm) next to the real use_pmm_wind_law;
   pmm_omega_sat, opal95_z2, alex_z2, surface_composition_file used nowhere. Add controls_check warnings
   (terminal only, byte-safe -- check pinned logs) and correct the registry rows. A, S.
5. Registry quality (defaults/controls_registry.tsv): 163 of 383 live rows have no description; the type column
   is read by no generator or test (one row is "?"); status values 'canonical' and 'verified' undocumented and
   treated the same, 'keep' means something else inside gen_inlist; migration notes remain ("VERIFIED 2026:",
   "remap:", "oversh:"); wrong descriptions (control.defaults:408,411 "blank = all columns" -- blank selects the
   30-of-89 and 12-of-61 defaults; lnewvars "nuclear-rate variable set" -- user_inputs.f90:113 says the
   opposite; end_dcen). Add a test that types match controls_buffer; a home=ctrl|job column replacing
   JOB_MEMBERS. A, M.
6. run_files.f90:21-58: unit numbers stored as controls (26 fixed units assigned into controls_buffer
   variables :189-205, 276, 349, 359, 383, 404, 411-412, 426, 432; 29 star%ctrl members are unit numbers;
   dumped by test_controls; closed by number in yrec_reset's 7..99 sweep). newunit in the owning loaders.
   (Deferred to a later wave.)
7. history_output.f90:60-229: bare slot arithmetic vals(22+i), vals(32+i), vals(39+i), vals(47+i) where named
   slots exist (i_lum_pp1... in star_info_lib.f90:88-89, i_nu_* in neutrino_flux_layout); neut_diag_1/2
   columns opaque (He3+He3 and He3+He4 rates). Named slots A, S; renaming columns B.
8. Declarations: read_controls.f90 holds module controls_reader (Naming rule 1); legacy_controls,
   controls_buffer, io_units have no private and are imported whole; controls_buffer's scattered public
   attributes; duplicate use lines read_controls.f90:15/32-33 and 29/31; double precision, json-sized dummies
   in write_gyre_pulse. A, S.
9. derived_controls.f90 does ~eight jobs: the three fatal checks belong with controls_check; six
   terminal/log write pairs (:222-305) should be one report helper (mixes write(* and terminal_unit); an
   unknown aiso changes nothing silently (:313) where an unknown amix warns (:253); version_fmt is a reader
   global used only by the header (:370; built read_controls.f90:74-77). A, M.
10. run_cards.f90 F77 shape: ALL-CAPS comments and initials, commented-out code (:38, :112-113), .eq./.gt.,
   the rescale copy duplicated (:73-76, :100-103), card kinds 1/2/3 as bare integers (also yrec_output.f90:240-249).
   Named kinds, one rescale copy. A, S.
11. model_files.f90 header and duplicates: header (:7-11) says YREC7 and MODEL2 "share one layout" with the
   envelope records "skipped" (they are read :350-356, :583-588; MODEL2 is one row per shell; both readers take
   the same 26 arguments); read_mod_model's header says "three values" for four outputs (:98); "LOG J/M" at
   :457 is -log10 omega; writer and reader each spell out the slots list and the "columns:" literal
   (:45-47/118-120, :62-64/121-124) -> one module parameter each. A, S.
12. Stale comments: legacy_controls.f90:12 (reset_legacy_controls); derived_controls.f90:38 (PARMIN);
   profile_output.f90:12 ("stx_* arrays"); controls_check.f90:16 ("after star_setup"; it runs after
   prepare_run, run_driver.f90:117-121); controls_check.f90:52, :60-61, :68 (overshoot_boundaries.f90, a
   composition_step dispatch, setcal/setscal); controls_buffer.f90:479-480 calls the S-factors SF-I ratios
   (they are absolute keV b from SF-II, legacy_controls.f90:565-569); gen_inlist.py:176-180 writes "parmin"
   into the generated file; provenance notes at yrec_output.f90:84, :102-106, pulsation_output.f90:327. A, S.
13. Output state in module variables: yrec_output.f90:35-36, history_output.f90:24-26, profile_output.f90:23-25,
   run_log.f90:22-25; the history unit is a newunit found by inquire(file=) (history_output.f90:44, :244) so the
   7..99 sweep never closes it and a re-run into another directory leaks it. An output_state record owned by
   run_yrec. A, M.
14. Single-precision literals inherited from F77 DATA: sstandard at legacy_controls.f90:372-373; frac_o =
   0.468195 at :443 (vs 0.468195d0 at :448); pmmjd, pmmmd, pmmsolw at :536-539; m = 0.22 and c_2 = 0.0506 at
   derived_controls.f90:163, :180. B (moves PMM-wind decks ~1e-8), S.
15. intent(inout) on pure readers: every writer/logger takes s intent(inout) (history_output.f90:158,
   profile_output.f90:122, pulsation_output.f90:42, :152, run_log.f90:52, controls_check.f90:28,
   yrec_output.f90:182). intent(in). A, S.

## Structural recommendations
1. Rebuild the controls pipeline around the registry: namelist groups bind the canonical controls_buffer
   names directly; the legacy path becomes a generated legacy->canonical translation (replacing the
   hand-written adopt_canonical_names, 230 lines, and the mirror lines in set_legacy_defaults); the registry
   carries the home (ctrl|job), validated types and validation rules, so JOB_MEMBERS goes. Adding a control =
   one declaration, one registry row, one generator run. L. A if test_controls stays bitwise.
2. Named column tables, one per output: one array of (name, kind) records with named indices (ih_star_age,
   ip_logT); profile names next to profile_value in core, or both sides use shared ip_*;
   test_default_columns_lists_in_sync reads that table; shared helpers write the MESA header and a selected row
   for both writers. A column = one row plus one value case. M. A.
3. Split read, derive and open: read_controls only fills controls; each derivation goes to the settings
   builder that owns it (wind law to rotation, EOS/opacity selection to microphysics settings, diffusion
   implications to element_diffusion); one validation module holds errors and warnings; run_yrec opens run
   files and output after the read, keeping the order of log writes. Removes the mid-read
   store_controls_to_star (read_controls.f90:92-94); test_controls stops depending on starting-model files;
   io stops depending on eos. M-L. A if the write order is kept.

## Newcomer traces
Add a control: io/controls_buffer.f90 (declaration + initializer); registry row (legacy column must equal the
new name; status 'verified' or 'canonical'; type decorative); JOB_MEMBERS if run-mutated; gen_controls_state,
gen_inlist, gen_controls_dump (6-7 generated files); settings field + copy in setup/; optional check in
derived_controls/controls_check; YREC_REGEN_CONTROLS=1 (every hashes.txt line changes); gen_makefile_deps.
3-4 hand-edited files, 6-7 generated, 2 expected.
Add a history column: star member + core/diagnostics/observables.f90; history_output.f90 (n_hist_cols,
names(90), vals(90), is_int_hist_col); defaults/history_columns.list (as !name if opt-in); regenerate
default_columns.inc and reseed MESA baselines if default. Inserting anywhere but the end renumbers everything.
Add a profile column: stitched_model.f90 (n_prof_cols, interior/envelope/atmosphere cases -- a forgotten region
writes silent zeros); profile_output.f90 names(62) in the literal form the test's regex expects (integer
column: a col_* parameter and an edit at :168); the .list file.

Fine as is: io_units.f90 (bar default-public), output_columns.f90 (bar not using check_io at :53-56),
run_log.f90 (bar intents), test/test_controls.f90 and the generated-dumps approach, tools/generated_files.py,
gen_controls_dump.py, gen_default_columns.py, pulsation_output.f90 beyond known ROADMAP debt,
controls_check's structure.

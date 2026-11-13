# Round 2's measured fixes for core (Kind B), for the owner's approval (2026-11-13)

PLAN.md section 7: these run after round 2's byte-identical work (merged into yrec-next at 8d7feee2), one item per
commit, each with its measured change on the 56 decks and 35 pins and the pins reseeded in it. Numbers are
docs/bug_sweep.md's. The two core agents measured in scratch, against copies of the round's baseline and pins: agent
M (model/, run/, diagnostics/) built each fix in its own sparse clone at 350c9186 and ran `compare_decks --set all`
(56 decks, 35 pins), then the moved decks alone for their history files (its report:
`scratch_core_m/blist/BLIST_REPORT.md` in the work directory); agent E (envelope/, structure/, step/) built its two
cheap trials from scratch copies of its branch and ran only the decks they reach. Some decks of the set turn any
perturbation into a far larger change: the 12 testsuite solar decks, the `rv_` decks and o24_rescale_rotating take a
1e-9 change to 1e-6..1e-4 in their final log L, with one model more or fewer (rv_sbm2 more; docs/ROADMAP.md, round 1B).
For them an item gives the **first visible change**: the first progress line of the run log (printed every 50
models) that differs from the base, and by how much. That is the fix's own size; the end state after it is the
deck's amplification. "Decision" marks the choices that are the owner's; each carries a recommendation.

**The owner's choices** (items 10, 15, 5 and 16): M132, an envelope derivative failure as an error (exit 1) or the
historical clean stop; M176, part A only or also part B; M31, the narrow form or the broad; M303, a change of a
pinned deck. M190 and M191 were decided on 2026-11-10 (PLAN.md section 11: fix both).

## Byte-identical on the 56 decks

1. **O4** (the starting model's stitched atmosphere reads `env_log10_radius(0)`). Fix: in stitched_model's
   `atmosphere_radius_cm`, with no envelope point stored, the photosphere's radius from L and Teff (the expression
   `surface_log10_radius_cm` gives atm_get). **Recommendation:** fix. Moves: no output file. Measured: `--set all`
   ALL IDENTICAL, 35 pins; `-fcheck=all` runs of cm_base and both mass-loss decks, which stopped at
   stitched_model.f90:475, now run their 12 models; cm_base with every profile and pulse file written is
   byte-identical (the first profile written comes after the first step has rebuilt the stitched model). After it,
   devel_tests' four O4 known failures (test_structure's three cases, test_run_contracts) are to be run again and
   dropped if clean, and the battery can add a `-fcheck=all` run of the configuration matrix.
2. **O5** (the accretion term of the entropy equation reaches every zone at levels 1-2). Fix: in `solve_level`,
   `envelope_zone_index = s%envelope_cz_bottom_index` (the predictor's edge) at levels 1-2 instead of 0.
   **Recommendation:** fix. Moves: nothing today (the shifts it gates are zero). Measured: cm_accretion, the one
   accreting deck, .log and .mod identical; the other 55 decks do not accrete (not run).
3. **O6** (the rezone's deuterium gate reads the pre-rezone surface index). Fix, A: `regrid_composition_terms` gates
   on the new surface, `s%xa(i_h2,new_point_count)`; B: regrid the start rates unconditionally. **Recommendation:** A
   now; B with item 15's part B, after which the rates it regrids are defined. Moves: A, none; B, 18 decks and 5 pins
   (the three testsuite rotating solar decks, Test_solar_GS_rot_gray_OPALSCV_SF2_GS98_OP_CF10, Test_solarcal_noGS_norot)
   and 13 of the 15 rv_ decks. Measured: A, `--set all` ALL IDENTICAL, 35 pins; B, first visible change 1e-9 in 8
   decks, the end states at the decks' noise (Test_solarGS98_difrotmixper18p4d0p1 855 to 854 models, log L -3.4e-5;
   Test_solarcal_noGS_norot log L -7.2e-6; rv_sbm2 898 to 904 models).
4. **M126** (a second core extension in one run starts from the first one's log-mass shift). Fix:
   `extend_core_toward_center` computes the points to add from the deck's `core_mass_reduction_factor` at every
   extension, not from the stored shift. **Recommendation:** fix. Moves: none of the 56 (no shipped deck extends the
   core and reads its starting model on a second card). Measured: `--set all` ALL IDENTICAL, 35 pins;
   Test_solarcal_noGS_norot with `extend_core_inward` on: the base's second cycle adds one point (N = 1, "del is 0/0")
   and stops with a NaN, exit 1, no verdict; fixed, every cycle's model has 1020 zones and the calibration converges
   after 3 cycles.
5. **Decision, M31** (the structure-based timestep limits stay off for the rest of the run once switched off).
   Fix, narrow: `begin_kind_card` re-arms `structure_dt_limits_active` from its control on a card that reads the
   starting model; broad: on every card. **Recommendation:** the narrow form; the broad is a physics choice (should a
   main-sequence card restart with the pre-main-sequence limits?). Moves: narrow, none; broad, 3 decks, 2 pins.
   Measured: narrow, `--set all` ALL IDENTICAL, 35 pins, and a four-card copy of cm_settle_he that died at card 4
   ("blend_eos: beta = P_gas/P is not positive", M31's symptom) runs to its 1e8 yr stop; broad,
   Test_solar_GS_norot_gray_OPALSCV_SF2_GS98_OP_CF10 443 to 447 models (log L -5.1e-5), Test_solarcal_noGS_norot 442
   to 443 (-2.7e-6), o24_rescale_rotating's age at its model stop -14% (log L -2.3e-3).
6. **M191** (Monte Carlo sets cards 2 and 3's stop ages whether or not the solar calibration is on; decided: fix).
   Fix: `apply_monte_carlo_parameters` sets them only with `calibrate_solar_model` (dropping the two lines is the same
   run: `setup_solar_calibration` sets both afterwards). **Recommendation:** the gate. Moves: none of the 56.
   Measured: cm_monte_carlo .log, .mod and history identical; its copy with card 2 allowed 3000 models runs to the
   deck's max_age (350 models) instead of stopping at 1e8 yr (170).
7. **M34** (the Monte-Carlo summary prints zero calibration derivatives after a first-check convergence). Fix: record
   2 prints calibration's fixed derivatives; calibration's four star members go. **Recommendation:** fix. Moves: no
   deck (none calibrates with Monte Carlo). Not measured.
8. **M115** (the Monte-Carlo run range and parameter file are not checked). Fix: `setup_monte_carlo_runs(s, ierr)`
   reads with iostat and reports a short file (core_report, ierr = 1), `prepare_run` passes it on; io's
   controls_check rejects `mc_run_start` < 1, `mc_run_end` > `max_mc_runs` and a start after the end.
   **Recommendation:** fix. Moves: error paths only. Not measured.
9. **M112, core's part** (core_messages' run-log half bypasses log_line's guard). Fix: `core_report`, `core_log` and
   `core_stdout_and_log` send their run-log lines through run_messages' `log_line`. **Recommendation:** fix. Moves:
   no deck expected (the run log is open whenever core writes); a host or a test that writes before the log opens
   sees queued lines. Not measured (the routing it depends on landed in this round; one compare).

## Changes of results

10. **Decision, M132** (a derivative failure inside the envelope's bsstep ends the run as "solution diverged", exit
    0). Fix: the kernel returns bsstep's `derivatives_ierr` (the EOS's or opacity's code) as a positive error, and
    keeps `envelope_termination` for a step that underflows. The owner's choice: such a failure as an error (exit 1;
    pyyrec's status positive) or the historical clean stop. **Recommendation:** an error: it is a failure of the
    microphysics, not the solution diverging, and the same failure at the step's end already exits 1. Moves: no deck
    (none fails there). Not measured.
11. **M85 and M96** (after a later card's Z rescale, the EOS mixture and the convection zones' O16 floor keep the
    first card's metals). Fix: `rescale_interior_z` recomputes the envelope mixture with `update_surface_mixture`'s
    arithmetic (one routine they share) and installs it, on a card after the first. **Recommendation:** both in one
    commit (convection's B item 5 waits for it). Moves: Test_solarcal_noGS_norot (pinned), the decks with a later Z
    rescale to a new Z. Measured: `--set all`, one deck, Test_solarcal_noGS_norot: 442 models both, log L -5.4e-11,
    log Teff -2.0e-11, ages in the 9th digit; M85 alone gives the whole change, M96 alone is byte-identical (its floor
    binds in no deck run); testsuite/Test_m0p7feh-0p75_GN93_TAHB cut to 300 models: at most 1.3e-6 in a history
    column, log L +4.3e-9 at the end.
12. **M86** (the first card's envelope refit runs with the run plan's EOS mixture). Fix: `update_surface_mixture`
    called once, right after `rescale_model` in `rescale_and_refit_envelope`, not in `read_starting_model_file`.
    **Recommendation:** fix. Moves: the run logs of Test_m0030_feh0_noGS_ZAMS (pinned) and
    Test_m10000_feh0_noGSstdcno_ZAMS. Measured: `--set all`, those two logs only: their Saha non-convergence lines go,
    2,212 and 5 to 0, the 2,217 of M293; every other byte identical (the table EOS stages replace the Saha result at
    every refit point).
13. **M97** (a Monte-Carlo run's L☉ is 3.844e33 erg/s times its sample's factor, whatever the deck sets). Fix:
    `physics%solar_luminosity_cgs` times the factor; the constant goes. **Recommendation:** fix. Moves: every Monte
    Carlo deck; cm_monte_carlo (pinned). Measured: cm_monte_carlo's log L +1.809e-3 dex (log10(3.844/3.828)), log Teff
    -4.9e-7, same 12 models; its progress lines then equal cm_base's.
14. **M190** (a later Monte-Carlo run, or a card that re-reads its starting model, takes the previous one's last age
    as the disk-locking age; decided: fix). Fix: `s%disk_gate_age_gyr = s%age_gyr` at every model read, at the end of
    starting_model's `reset_state_after_read`. **Recommendation:** fix. Moves: testsuite/Test_solarGS98_difrotmixper18p4d0p1
    and its pin (it starts past its t_disk, and today its first step reads 0). Measured (that deck alone): 855 to 854
    models; at 4.568 Gyr log Teff 3.761337 to 3.761338, X_c 0.355849 to 0.355768, the run log's log L column -0.000042
    to -0.000025. A path-sensitive deck: its first visible change was not measured; the other 55 decks were not rerun
    (the entry's own 56-deck run found only this one).
15. **Decision, M176** (a card that reads the starting model takes the previous card's start-of-step deuterium
    rates). Fix, A: `reset_state_after_read` clears `deuterium_burning_rate_start`; B: `compute_burn_rates` clears
    the start rate wherever it clears the end rate. The owner's choice: A only, or also B. **Recommendation:** A now;
    B, the right invariant (it makes A, M169's reset and item 3's B redundant), in a round that reseeds the set
    anyway. Moves: A, Test_solarcal_noGS_norot (1 pin); B, 49 decks and 17 pins. Measured: A, `--set all` one deck:
    442 models both, log L -4.5e-10, the final timestep 1.585122e7 to 1.585140e7 yr, ages in the 8th digit; B, the
    configuration matrix at round-off (final log L within 1e-11), first visible change 1e-9 to 1.5e-8 in 22 decks,
    the testsuite -3.2e-5 to +2.7e-5 in log L, rv_iter1 652 to 639 models and rv_nocore 680 to 719 (log L -8.2e-3
    and -8.0e-3, first visible change already 1.6e-6 and 3.2e-6 at model 200), o24_rescale_rotating +7.2e-4.
16. **Decision, M303** (cm_env_overshoot never reaches envelope overshoot). Fix, a deck change: keep its starting
    model and give card 2 `max_model_number(2) = 600` and `max_age(2) = 1.0D7`. **Recommendation:** do it: no model in
    startmodels/ starts with a radiative core and runs (the 1 Msun TAMS model fails its first burn). Moves:
    cm_env_overshoot's log and .mod, and its pin. Measured: 78 models in 8.9 s (compare_decks' quick set about 7 s
    longer); envelope overshoot acts 127 times, the Li/Be overshoot depth 31 and 32 times, the rezone's overshoot base
    31 times; surface Li7 4.09e-10 against 4.27e-9 without overshoot.
17. **M273, the B part** (the shell kernel's inverse year is not its year's inverse; the structure timestep limit's ln
    10 is 2.3026). Fix: `1/seconds_per_year` for `one_year_sec_inv` (9.9e-6 relative; the entropy term's floor while
    dt < 1 yr and X > 0.01) and physical_constants' `ln10` in the L criterion of the structure limit (6.5e-6).
    **Recommendation:** both, in one measured commit. Moves: every deck where the L criterion sets the timestep (the
    pre-main-sequence and most evolved decks) and the first steps of any deck with dt < 1 yr: most pins, at the 1e-6
    level. Not measured.
18. **M270** (the starting model's mismatch warning prints an unset flag and calls itself an error). Fix: print
    `header%use_extended_composition`, drop `lexcp0`, one core_report warning ("WARNING: the starting model differs
    from the deck; the run goes on" / "  mixing length: deck d, model m" / "  extended composition: deck F, model T").
    **Recommendation:** fix, with the messages (item 24). Moves: cm_basic_composition.log (pinned) and rv_basic.log.
    Measured: those two logs only, the warning's 4 lines to 3; the model's flag is T where the old line printed F.
19. **O7** (accretion mixing assumes a uniform convection zone; the step uses two definitions of M_cz). Fix:
    `mix_accreted_composition` takes dm and mixes X_new = (Σ X dm + X_acc dM)/(Σ dm + dM); Σ dm is M_cz in the
    timestep limit, the entropy shift and the compression too. **Recommendation:** fix. Moves: accreting runs only:
    cm_accretion (no pin). Not measured.
20. **O19, core's side** (the accretion surface after a failed Allard lookup keeps the last envelope integration's
    point, which may be another Teff's and g's). Fix: integrate the gray atmosphere at the star's Teff and g for the
    fallback (atm_get, not saving the boundary), or return an error. **Recommendation:** the gray atmosphere: the
    state only seeds the accreted matter's EOS call, which an error would stop the run for. Moves: accreting runs
    with the Allard tables that leave them; none in the sets. Not measured.
21. **O8** (`store_li_be_rates` zeroes both rate sets past the first cool shell). Fix: clear only the set the call
    fills. **Recommendation:** fix, its own measured commit. Moves: surface Li and Be of the extended-composition
    decks (the testsuite's solar decks with light elements). Not measured this round (an earlier trial moved surface
    Li slightly).
22. **O9** (surface 7Li has a noise floor of about 0.005 dex). Fix: the start-of-step convection-zone average takes
    the fractional edge as the end-of-step one does. **Recommendation:** fix, after item 21, measured. Moves: the same
    decks, at the 0.005 dex level. Not measured.
23. **Decision, O1 (core's part) and O3.** O1: shell_equations' convective `qt_dl = 0` (marked `known` in
    test_structure) needs convection's d∇/dL; O3, the low-mass convergence, is an investigation.
    **Recommendation:** leave both to a later round: they change the Jacobian, so every deck's convergence path moves
    at the solver's tolerance (bug_sweep section 1's decision). Not measured.

## Messages

24. **The messages reworded**, one re-pin commit, text and destinations only; nothing else moves. Not measured. Logs
    move only where the lines fire: the mass-loss decks' standard output, a deck that reaches log Tc 7.1 with L_grav
    < 0 (the switch-off note), the solver-diagnostics runs, ZAHB_m1p0feh+0p0_TAHB's terminal (agent M's 10th below),
    and the logs of item 18. **Recommendation:** all of them, in one commit. The full list:
    - **The routes.** Fold `core_stdout` into `core_terminal`: they differ only before io assigns its units
      (standard output against standard error). Fold `core_stdout_and_log` into `core_report`. The `*_values` routes,
      which keep list-directed lines' trailing blanks, stay.
    - **envelope/, structure/, step/ (agent E):**
      1. mass_change's " TIMESTEP REDUCED FOR MASS LOSS", "MASS LOSS APPLIED", the wind's dJ line, " mass change:
         turnover time" and " mdot: envelope boundary at the surface (911)" reach standard output only, while their
         comments said the run log: through core_log or core_report;
      2. the structure limits' switch-off ("timestep control: structure-based limits disabled") and
         rebuild_envelope's notes reach standard output only: add the run log (they mark a change of the step's
         physics); energy_sources' registry errors likewise;
      3. the envelope kernel's legacy texts ("ATMOSPHERE INTEGRATION FAILED AFTER MAXSTP INTEGRATIONS.I QUIT.",
         "ENVELOPE INTEGRATION FAILED AFTER MAXSTP TRIES. I QUIT", the list-directed "ENVINT: X0>XLIM, X0,XLIM:",
         " bsstep: <status>" without the integration or the vertex) as sentences with their context;
      4. "-----CORRECTIONS EXCEEDED TOLERANCES" and "MODEL CONVERGED AFTER n ITERATIONS" in sentence case.
    - **model/, run/, diagnostics/ (agent M):**
      1. starting_model "STARIN: ***** RUN STOPPED ***** / MODEL FAILED TO CONVERGE" → "ERROR: starting model: the
         previous model failed to converge; run stopped" (core_report);
      2. "STARIN: ***** RUN TERMINATED, INVALID INPUT MODEL FILE. *****" → "ERROR: starting model: unknown model file
         layout '<tag>' in <path> (expected YMOD, NMOD or MOD2); run stopped";
      3. item 18's warning;
      4. "STARIN: Unable to extend core inward - JSON too small / ... / RUN TERMINATED" → "ERROR: core extension: the
         model needs N shells, more than json (M); run stopped";
      5. "ERROR IN SUBROUTINE STARIN / DESIRED ENVELOPE MASS ... TOO LARGE / ENVELOPE MASS NOT CHANGED" (a warning) →
         "WARNING: starting model: the requested envelope fit mass log(m/M) = ... lies inside the innermost shell;
         the envelope mass is not changed";
      6. "STARIN: degenerate envelope interpolation interval (was STOP 9998)" → "ERROR: starting model: a degenerate
         interval in mass while appending the envelope shells; run stopped";
      7. rezoning "ERROR IN HPOINT / MAX. STAR T LESS THAN MINIMUM ENVELOPE T.RUN STOPPED" → "ERROR: rezone: the
         hottest shell is cooler than the minimum envelope temperature tenv0; run stopped";
      8. rezoning's format 101 "DESIRED NUMBER OF POINTS N EXCEEEDS JSON" (misspelt, no blank, standard output only)
         and the two "EXCEEDS JSON" lines → "ERROR: rezone: the new mesh needs N points, more than json (M); run
         stopped" (core_report);
      9. model_rescaling "ERROR IN SUBROUTINE RSCALE ... X (Z) NOT RESCALED" (warnings) → "WARNING: rescale (card
         n): the requested X = x is above 1; X is not rescaled" (and Z); the species message, which ends "Z NOT
         RESCALED", → "... species k is not rescaled";
      10. model_rescaling's debugging prints "Entering old/new method", "SENV", "Envelope" (on ZAHB_m1p0feh+0p0_TAHB's
          terminal) → dropped, or one core_log note;
      11. "ERROR IN SUBROUTINE RSCALE / NEW ENVELOPE MASS LESS THAN ZERO ..." → "ERROR: rescale (card n): the new
          envelope mass outside the H shell would be negative (...); run stopped";
      12. stop_conditions "LROT set to TRUE, but OMEGA(1) = 0 ..." (F77 names, standard output only) → "ERROR:
          rotation_active is on but the model's central omega is 0: set set_initial_omega and initial_omega, or turn
          rotation_active off; run stopped" (core_report);
      13. calibration "NK=2, Y=Y+0.01, Setup next run, X=" → "star calibration, first pass: X lowered by 0.01 (Y
          raised), next X ="; " Setup next run, NK, X =" → "star calibration: secant step, card n, next X =", with
          the run log's "#";
      14. turnover_timescale "splint: <status>" (3 places) → "ERROR: turnover time: the spline of the convective
          velocity (or H_P) at the convection-zone base failed: <status>" (core_report).

## Not on either agent's list

25. **M146** (the timestep limiter's hydrogen-burning luminosity is channels 1-4 by number). Round 2 named them: the
    limiter takes `i_lum_pp1` .. `i_lum_cno` and `i_lum_3alpha` (136d1a74), and the explicit-shape dummy is gone. Left:
    a burning group per chain from net/, one sum the limiter and henyey_iterate's log line share (its channels 6-7
    are still by number), and the "Add a chain" recipe. **Recommendation:** do the rest with net's next round; it is
    byte-identical by its entry's account (Kind A). Moves: nothing. Not measured.
26. **O22** (profile column 54's N² disagrees with the pulse files'). Fix: one N² for both, the thermal form where it
    is negative. **Recommendation:** fix. Moves: column 54 of the profiles, in convection zones (no pinned profile).
    Not measured.
27. **Decision, O20** (the calibration's dL/dα is 0.0139, its comment's fit +0.139). Already the authors' question
    in bug_sweep section 1. **Recommendation:** keep it open until the fit is found; it moves only the solar
    calibration's path. Not measured.

## Round 1B's questions still open

28. **Decision, M278**, the conductive mixture's weights: B (each species weighted by its share of the electrons,
    each table read at its own density), A (the weights alone), or leave it. **Recommendation:** B, as round 1B
    researched it. Moves: every deck (the solar models by about -1.7e-5 in log L at 4.57 Gyr); the ZAHB trial is to
    be measured again on today's base.
29. **Decision, M285**, the overshoot merge: keep the merged fix (the union of the two regions), which lengthens
    ZAHB_m1p0feh+0p0_TAHB's core helium burning by 9.8% (1517 to 1917 models), or restore the F77's merge?
    **Recommendation:** keep it: the union is what a merge of two regions means.
30. **Decision, M293**, the Saha solver's non-convergence: a warning or an error. All 2,217 in the compare set come
    from M86's refit order, and item 12 removes every one of them. **Recommendation:** an error, after item 12 and a
    count over the shipped decks outside the set; measured then.

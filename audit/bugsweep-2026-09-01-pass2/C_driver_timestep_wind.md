# Bug sweep pass 2 -- Group C: evolution driver, stopping, timestep control, calibration, setup mapping, mass loss / wind

Date: 2026-09-01. Read-only review; no files edited, nothing built or run.
F77 reference: rev 6cd5673 (`git show 6cd5673:src/<name>.f`). Ancestor map used
for files that were renamed in the modernization: htimer.f/ytime.f/xtime.f/
ptime.f/wtime.f/entime.f -> util/compute_timestep + timestep_limit_*;
wind.f -> kawaler_wind; mwind.f -> matt_wind; cowind.f -> wind_spindown;
mcowind.f -> wind_spindown_matt; amcalc.f -> matt_structure_factor;
setcal.f/setscal.f -> setup_solar/star_calibration; rscale.f -> rescale_model;
remap.f -> map_user_inputs; findsh.f -> locate_shell_boundaries;
chkcal.f/chkscal.f -> check_solar/star_calibration; wrtmonte.f -> monte_carlo;
neutrino.f -> neutrino_flux_table.

Process disclosure: while grepping ROADMAP.md for an unrelated flag name early in
the pass I printed ~15 lines that turned out to be part of ROADMAP section 10
(first-sweep items 6-9, one of which concerns helium_dt). I stopped reading
immediately and did not consult ROADMAP or audit/bugsweep-2026-08-31/ further.
Finding 1 below had already been reached independently from the code before
that; everything else is from the code and the F77 alone.

## Files read in full

core/evolve_step.f90, core/run_yrec.f90, core/stop_conditions.f90, core/main.f90,
core/star_setup.f90, core/yrec_reset.f90, core/yrec_capi.f90,
core/check_solar_calibration.f90, core/check_star_calibration.f90,
core/monte_carlo.f90, core/neutrino_flux_table.f90,
util/compute_timestep.f90, util/timestep_limit_hburn.f90,
util/timestep_limit_heburn.f90, util/timestep_limit_hr.f90,
util/timestep_limit_omega.f90, util/timestep_limit_structure.f90,
util/version.f90, util/xrng4.f90, util/zero.f90,
setup/setups.f90, setup/setup_solar_calibration.f90,
setup/setup_star_calibration.f90, setup/map_user_inputs.f90,
setup/locate_shell_boundaries.f90, setup/rescale_model.f90,
wind/kawaler_wind.f90, wind/massloss.f90, wind/matt_structure_factor.f90,
wind/matt_wind.f90, wind/mdot.f90, wind/wcz.f90, wind/wind_spindown.f90,
wind/wind_spindown_matt.f90.

Also read for call-site / consumer checks (not in full unless noted):
rotation/seculr/secular_transport.f90 (wind call sites 262-300, 405-415),
rotation/evolve_angular_momentum.f90 (140-150, 276-290),
core/read_starting_model.f90 (450-465), io/read_controls.f90 (2000-2015,
Kawaler constant mapping), net/net_lib.f90 ~2059 and net/burn_lib.f90 ~1359
(consumers of qs0e_scale/qqs0ee_scale), state/star_info_lib.f90 (composition
index constants).

F77 read in full from 6cd5673: main.f, setups.f, massloss.f, mdot.f, wcz.f,
wind.f, cowind.f, mcowind.f, mwind.f, amcalc.f, wrtmonte.f, neutrino.f,
gmass.f, chkcal.f, chkscal.f, setcal.f, setscal.f, rscale.f, findsh.f,
getrot.f, locate.f, zero.f, xrng4.f, fpft.f, remap.f. (htimer.f/ytime.f etc.
were read via the modernized headers' quoted originals and main.f.)

Bare-external call-site checks (one `grep -rn "call <name>" /Applications/YREC/src`
each, no head): xrng4 (kap/laol89/gtpurz.f90:49, gtlaol.f90:64, gtlaol2.f90:62),
setups, map_user_inputs (net/test/test_net.f90:141, io/read_controls.f90:1757),
rescale_model (core/read_starting_model.f90:330, :690), locate_shell_boundaries
(evolve_step.f90:115, run_yrec.f90:341, rescale_model.f90:181), compute_timestep
(evolve_step.f90:139, run_yrec.f90:348), massloss (evolve_step.f90:255),
setup_solar_calibration (run_yrec.f90:190), setup_star_calibration
(run_yrec.f90:200), check_solar_calibration (run_yrec.f90:223),
check_star_calibration (stop_conditions.f90:98), timestep_limit_*
(compute_timestep.f90:92,107,120,129,140), wind_spindown
(wind_spindown_matt.f90:74), matt_structure_factor (wind_spindown_matt.f90:126,
matt_wind.f90:101), mdot (massloss.f90:314), kawaler_wind (matt_wind.f90:58),
matt_wind (evolve_angular_momentum.f90:176,287; secular_transport.f90:280,293),
wind_spindown_matt (secular_transport.f90:410), zero (eos/mhd/mhdpx2.f90:47),
wcz (no callers -- dead, documented as such). Actual/dummy counts and kinds
match at every site checked (massloss 22/22, mdot 23/23, compute_timestep 21/21,
matt_wind 12/12, kawaler_wind 12/12).

## Verified clean (compared line-by-line against the F77 and/or the cited physics)

- setups.f90 vs setups.f; setup_solar_calibration.f90 vs setcal.f;
  setup_star_calibration.f90 vs setscal.f; locate_shell_boundaries.f90 vs findsh.f;
  wcz.f90 vs wcz.f (plus the documented `iend.lt.num_zones` guard);
  wind_spindown.f90 vs cowind.f; matt_structure_factor.f90 vs amcalc.f;
  kawaler_wind.f90 vs wind.f; matt_wind.f90 iteration vs mwind.f;
  wind_spindown_matt.f90 vs mcowind.f; xrng4/zero vs xrng4.f/zero.f;
  monte_carlo vs wrtmonte.f; neutrino_flux_table vs neutrino.f;
  check_solar_calibration vs chkcal.f; check_star_calibration vs chkscal.f.
- run_yrec / evolve_step / stop_conditions control flow vs the main.f run loop
  (kind-card loop, LFIRST/rescale ordering, calibration re-entry, stop tests).
- Kawaler 1988 mapping (read_controls.f90:2008-2012): exmd=1-2a/3, exw=1+4a/3,
  exr=2-a, exm=-a/3 -- matches dJ/dt ~ Omega^{1+4an/3} R^{2-n} M^{-n/3} Mdot^{1-2n/3}
  with n=1.5 folded into the constant.
- Matt et al. 2012 centrifugal factor (matt_wind.f90:~118): fcorr = Omega^2 R^2/v_esc^2,
  fcen = ((c2^2+fsun)/(c2^2+fcorr))^{m}, c2=0.0506 -- correct form and exponent.
- Matt et al. 2015 Rossby scaling (extau/omega_saturation handling) consistent
  between matt_wind and wind_spindown_matt.
- Reimers (massloss.f90:157-161): Mdot[Msun/yr] = eta_fac * (L/Lsun)(R/Rsun)/(M/Msun);
  the 2026 fix is correct and the argument list to mdot is now 1:1.
- timestep_limit_structure / _omega / _hburn / _hr arithmetic; compute_timestep
  min-selection, atime(13) growth cap, end-age clamp.
- The F77 wind.f aliased its result onto COMMON/ROT/WNEW (the user's initial
  omega); the modern local `omega_iter_new` removes that clobber. This is a
  behavioral divergence from the F77 for multi-run jobs with set_initial_omega,
  but in the direction of correctness; noted, not reported as a bug.

---

## util/timestep_limit_heburn.f90:111 -- helium_dt read before it is written at core-He exhaustion -> zero timestep

- class: numerical
- severity: high
- confidence: high
- provenance: inherited (ytime.f: `DELTSY = (5.85D17/CLSUN)*DELTSY*(HS(JXBEG-1)/HL(JXBEG-1))`) aggravated by modernization (htimer.f's `SAVE` on DELTSY is gone; compute_timestep.f90:67 makes helium_dt a plain local)
- detail: In the non-flash branch, when the core helium fraction `1-Z(jcore)` drops below atime(1) (i.e. core He exhausted, Z_core -> 1), the routine computes `helium_dt = const*helium_dt*(M/L)` using helium_dt on the right-hand side. The dummy is intent(inout) and the header admits this "is preserved exactly, not fixed", but the actual argument is now an uninitialised local of compute_timestep (line 67), so with `-finit-local-zero` the result is exactly 0 (without it, stack garbage). compute_timestep then takes `min(..., helium_dt, ...)` at line 157/162 and returns previous_timestep = 0; the "fritz-out" guard at lines 98/113 protects only structure_dt and rotation_dt, not helium_dt. A zero dt feeds the Henyey gravitational-energy term (dS/dt) and the composition update as a division by zero, producing Inf/NaN in the next model. This is the most plausible root cause of the reported symptom "run_from_zahb_to_tahb Solar_m1p0 GN93 dies with NaN in the MLT cubic near model 950 at core He exhaustion": the failure is at exactly the transition where this branch first fires. In the F77 the SAVEd DELTSY made the expression a (dimensionally meaningless but nonzero) recurrence on the previous call's value, so it "worked". Fix direction: the shell branch should compute a fresh timestep from atime(4)/atime(5) and the shell's M/L exactly like the core branch (or from the He-shell burning rate), never from the incoming value.

## setup/map_user_inputs.f90:133 -- qs0e/qqs0ee S-factor derivative scalings normalised by the wrong reference S0 (pp instead of each reaction's own)

- class: physical
- severity: medium
- confidence: high
- provenance: inherited (remap.f:132 `QS0E(2) = (S0P_3_3/S0_3_3)/(QS0E_3_3_A98/S0_1_1_A98)` and the same pattern for reactions 3,4,5,7,8 and the QQS0EE entries)
- detail: The scaling of the first- and second-derivative S-factor terms is meant to be (S'/S)_user / (S'/S)_A98, i.e. each numerator ratio divided by the Adelberger-98 ratio for the same reaction. Only reaction 1 (pp) and reaction 8 use their own S0 in the denominator; lines 133, 138, 143, 148, 153, 158 (qs0e) and 174, 180, 187 (qqs0ee) divide `qs0e_<rxn>_a98` by `s0_pp_a98` (= 4e-22 MeV b) instead of `s0_<rxn>_a98`. With defaults (e.g. He3+He3: s0p=-4.9, s0=5.21e3, qs0e_a98=-4.1) qs0e_scale(2) comes out ~9e-26 instead of ~1.24. Consumers net_lib.f90:~2059 and burn_lib.f90:~1359 multiply the FCZ eq.(52) Q2/Q3 (and Q4/Q5) energy-derivative corrections by these factors, so the derivative corrections to S_eff for He3+He3, He3+He4, p+C12, p+C13, p+N14, p+O16 are silently switched off whenever use_new_nuclear_rates (lnewnuc) is true -- which it is in examples/run_standard_solar_model/*.nml2 and run_from_dbl_to_zams. Effect is a few-tenths-percent to percent-level bias in the pp-II/III and CNO rates at solar temperatures (the pp chain itself is unaffected, which is why solar calibration still converges). Note also line 36's `9.4E0` single-precision literal in the s0_p_o16 reference data (harmless at ~4e-8 relative).

## setup/rescale_model.f90:372 -- Z-ramp block scales CNO by (Zold-Zc)/Zold instead of Zc/Zold

- class: physical
- severity: medium
- confidence: high
- provenance: inherited (rscale.f:343 and :351 `ZFACT = (HCOMP(3,I)-RSCLZC(NK))/HCOMP(3,I)`)
- detail: Under use_z_ramp with rsclzc/rsclzm1/rsclzm2 > 0, the inner-core block (line 372) and ramp block (line 380) set Z to the target value and then multiply species 5..11 (C12..O18) by `z_rescale_factor = (Zold - Znew)/Zold`. The main Z-rescale block at line 97 correctly uses `(Zold+dZ)/Zold = Znew/Zold`. The ramp factor is the *removed* fraction of Z, not the retained fraction: for Zc = 0.5*Zold the CNO abundances are halved (coincidentally right), for Zc = 0.9*Zold they are cut to 10%, and for Zc > Zold they go negative. The resulting model has a total Z inconsistent with its C+N+O, so the CNO luminosity, opacity (via Z) and the summed metal abundance disagree from the first model of the run; negative CNO mass fractions would produce negative burning rates or an EOS/opacity failure. Also: the block comment at lines 92-93 says "CNO CYCLE ELEMENTS AND HE3 ARE MULTIPLIED" but every loop runs 5..11 and never touches He3 (index 4); comment/code mismatch (inherited).

## wind/mdot.f90:267 -- shell-mass recompute loop reads zone_mass_grams(0) when the star is fully convective

- class: logical
- severity: medium
- confidence: high
- provenance: inherited (mdot.f:183-186, identical loop bounds)
- detail: After changing the total mass, the loop `do zone_idx = envelope_boundary_zone, num_zones-1: shell_mass(i) = 0.5*(zone_mass(i+1) - zone_mass(i-1))` starts at the CZ base. The routine explicitly anticipates envelope_boundary_zone = 1 (lines 116-121 and ~182 special-case it), which is exactly the fully convective pre-main-sequence configuration in which accretion (use_mass_accretion with a positive rate) is normally used. With jenv = 1 the first iteration reads `zone_mass_grams(0)`: an array-bounds abort under `-fcheck=all`, otherwise a garbage shell mass for zone 1 that is then carried into the Henyey mass coordinates. The loop should start at max(envelope_boundary_zone, 2) with shell_mass(1) handled separately (as the centre point already is elsewhere in the code).

## wind/matt_wind.f90:57 -- PMM (Matt) wind path ignores wind_loss_active; two secular_transport call sites are not gated on use_wind_torque either

- class: logical
- severity: medium
- confidence: high
- provenance: inherited (mwind.f likewise never tests LJDOT; only wind.f does) -- interface half is modernization/unclear (the secular_transport call sites carry the flag but the callee drops it)
- detail: `wind_loss_active` (the per-run header flag LJDOT / use_wind_torque) is a dummy of matt_wind but is only forwarded to kawaler_wind (line 59), which tests it at kawaler_wind.f90:77; the `use_pmm_wind_law` branch never looks at it and always applies the torque. Callers in evolve_angular_momentum.f90:148 and :282 gate the whole call on the flag, so those paths are safe. However secular_transport.f90:270-272 (`surface_cz_active .and. diffusion_velocity(zone_max).eq.0 .and. .not.disk_lock_active`) and the solid-body branch at :286-289 only set `wind_loss_active = star%job%use_wind_torque` and then call matt_wind unconditionally. Consequence: a model with use_pmm_wind_law = .true. and the header wind flag off (LJDOT=F, e.g. a "no wind" control run or a disk-locking test with instability transport on) still loses angular momentum through the Matt torque on every secular sub-step, while the same input with the Kawaler law correctly applies no torque. Fix direction: test wind_loss_active at the top of matt_wind (or gate the two secular_transport sites the way evolve_angular_momentum does).

## wind/wind_spindown_matt.f90:135 -- fsun uses uninitialised local `gl` instead of the gravitational constant

- class: numerical
- severity: low
- confidence: high
- provenance: inherited (mcowind.f `FSUN = 0.5*PMMSOLW**2*CRSUN**3/EXP(CLN*GL)/CMSUN`, GL never set; the modern header documents it as a "PRESERVED BUG")
- detail: The solar centrifugal reference term should be 0.5*Omega_sun^2*Rsun^3/(G*Msun) ~ 1.05e-5; with `gl` = 0 (`-finit-local-zero`) the divisor is 1 instead of G = 6.67e-8, giving fsun ~ 7e-13. matt_wind.f90:118 computes the same quantity correctly with `cgl`. Because fsun is added to c2^2 = 2.56e-3 the numerical impact on fcen is ~0.4% relative, so the two wind routines used on the same star (matt_wind on the outer loop, wind_spindown_matt inside the secular sub-steps) apply slightly inconsistent torques. Without the init flag the value is stack garbage and could be anything. Trivial fix: use cgl.

## wind/mdot.f90:130 -- "add no more than 1/2 of the mass beyond the fitting point" limit is applied to mass *loss* as well

- class: logical
- severity: low
- confidence: medium
- provenance: inherited (mdot.f: `RATE = ABS(DMDT)` then `IF(RATE.GT.0.0D0)`)
- detail: mass_loss_rate_cgs is the absolute value of the rate (taken a few lines earlier), so the condition at line 130 is true for every nonzero rate and the accretion-only envelope limit `0.5*(M_total - M(fitting point))/rate` (DTENV) throttles the timestep of mass-loss runs too. For a Reimers RGB wind with a thin fitted envelope this can force much smaller steps than fczdmdt/ftotdmdt intend, and the comment/code disagree. Related: if use_mass_accretion is on with mass_accretion_rate = 0 and no Reimers wind, lines 124-126 divide by zero (FPE trap in the dev build; Inf limits in release, which are then harmless).

## Weak/uncertain observations

- util/timestep_limit_heburn.f90:112 -- `enclosed_mass(h_shell_zone_begin-1)` indexes 0 when locate_shell_boundaries reports no H shell (h_shell_zone_begin=1) while the core-He-exhausted branch fires; unreachable for ordinary stars but unguarded (inherited ytime.f).
- util/compute_timestep.f90:51/71 -- hydrogen_dt is intent(out) but never assigned when previous_timestep < 0 (relaxation case); caller receives an undefined value (modernization: F77 had it SAVEd).
- util/compute_timestep.f90:200-205 -- if the model's age already exceeds target_end_age, time_left_years < 0 is copied into the timestep, producing a negative dt for one step before stop_conditions ends the run (inherited).
- setup/rescale_model.f90:144 -- the error WRITE prints `icomp`, which is never assigned in this routine (F77 ICOMP likewise undefined); prints 0/garbage in the diagnostic.
- setup/rescale_model.f90:189 -- env_mass_old2 is recorded before the mass rescale and reused for the core rescale, so a run that rescales both mass and core uses a stale envelope-mass fraction (inherited rscale.f).
- setup/map_user_inputs.f90:36 -- `9.4E0` single-precision literal in a double-precision DATA-style initialisation (~4e-8 relative error; harmless).
- wind/massloss.f90 -- on the mass-loss path (rate < 0) mean_molecular_weight_local, accretion_specific_energy, mean_thermal_energy and new_atmosphere_fit_needed are left unset before the mdot call; harmless today because mdot only reads them under delta_mass > 0, but fragile.
- wind/massloss.f90:301/304 -- rot_scr%envelope_specific_entropy is unconditionally forced to 0 (historical toggle), so mdot's delta_log_radius correction is dead code; comment says the accretion entropy is used (inherited).
- wind/kawaler_wind.f90:58-75 -- uses the current star%convective_turnover_timescale while matt_wind interpolates between old/new turnover times with fracstep; the two laws see slightly different Rossby numbers on the same sub-step (inherited wind.f vs mwind.f).
- util/timestep_limit_hr.f90 header says log_teff "confirmed logged: main.f computes 10.0d0**TEFFL" -- header refers to a computation that lives in the caller; the argument is indeed log10(Teff) in both call paths, so the comment is merely misleading.
- util/timestep_limit_heburn.f90:75-78 -- h2_fraction passed to engeb is uninitialised when use_extended_composition is false (zero under -finit-local-zero; engeb ignores it in that mode).
- core/check_solar_calibration.f90 -- log_zx_mismatch is only assigned when the Z/X calibration is on but is written to the log unconditionally (prints 0).
- core/check_solar_calibration.f90 -- comment states dlum_dalpha uses a derivative "from Basu"; the coded constant is 0.0139 with no source; value is only a first-guess seed so the calibration converges regardless.
- setup/setup_solar_calibration.f90 -- the central_x/y/d stop criteria (end_xcen etc.) added in 10/24 are not propagated to runs > 3 the way endage/setdt are (inherited setcal.f); only matters if a user sets those on a calibration run.
- core/run_yrec.f90 -- `log_r_rsun_current` is actually log10(R/Rsun) of the *previous* model at the point where check_star_calibration uses it (name misleading; logic matches chkscal.f).
- wind/wcz.f90 -- dead code (no callers); retained per batch scope, matches wcz.f.

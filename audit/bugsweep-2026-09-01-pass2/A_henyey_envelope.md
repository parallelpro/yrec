# Sweep A (pass 2): structure solver, envelope integration, temperature gradients, observables

Reviewer scope: /Applications/YREC/src, branch yrec-modern. READ-ONLY sweep; F77
ancestors from `git show 6cd5673:src/<name>.f`. Findings written independently of
the 2026-08-31 pass.

## Files read in full

Assigned (all 15):
- core/henyey_coefficients.f90 (vs coefft.f)
- core/henyey_iterate.f90 (vs crrect.f)
- core/henyey_eliminate.f90 (vs reduce.f)
- core/henyey_solve.f90 (vs hsolve.f)
- core/envint_kernel.f90 (vs envint.f)
- core/envint_lib.f90 (vs envint.f driver part)
- core/envelope_derivs.f90 (vs qenv.f)
- core/atmosphere_derivs.f90 (vs qatm.f)
- core/surfbc.f90 (vs surfbc.f)
- core/rebuild_envelope.f90 (vs getnewenv.f, getrot.f)
- core/shell_physics.f90 (vs physic.f)
- core/turnover_timescale.f90 (vs gettau.f, tauintnew.f)
- core/observables_lib.f90 (vs wrtout.f blocks)
- mixing/temperature_gradients.f90 (vs tpgrad.f)
- mixing/compute_scale_height.f90 (vs hsubp.f)

F77 ancestors read in full: coefft.f, crrect.f, reduce.f, hsolve.f, envint.f,
qenv.f, qatm.f, surfbc.f, getnewenv.f, getrot.f, physic.f, tpgrad.f, hsubp.f,
gettau.f, tauintnew.f; wrtout.f (relevant blocks).

Supporting files read (in part, to trace flow): core/evolve_step.f90 (corrector
ladder, rebuild_envelope call site), core/stitched_model.f90 (atm_get call,
delta(1) handling), rotation/omega_from_j.f90 (full), atm/ttau_lib.f90
(dispatch), eos/eqstat.f90 (which result slots are always overwritten),
wind/mdot.f90 (where rot_scr%delta_log_temperature/pressure are set),
mixing/overshoot_boundaries.f90 and core/read_starting_model.f90 / core/rezone.f90
/ rotation/evolve_angular_momentum.f90 (call sites of bare externals),
Makefile (compiler flags).

Build-flag context that matters for every "uninitialized local" finding below:
Makefile line 85 puts `-finit-local-zero` in GNU_BASE_FFLAGS, so under gfortran
every plain local starts at zero (deterministic, but not necessarily the
intended value). The ifort/ifx flag sets (Makefile lines 162-171) do not have an
equivalent, so under those compilers the same reads are true garbage. The DEV
flags add `-ffpe-trap=invalid,zero,overflow`, which turns some of the masked
divisions below into hard traps.

## Verified clean (re-derived, found correct)

- henyey_eliminate.f90 and henyey_solve.f90: the forward elimination and back
  substitution are line-for-line equivalent to reduce.f / hsolve.f (index
  ranges, pivot ordering, the 4x4 block layout).
- henyey_coefficients.f90: the difference-equation coefficients and their
  partials with respect to (P,T,R,L) at both shell faces match coefft.f
  (including the dgrad/dP, dgrad/dT, dgrad/dR, dgrad/dL chain through
  temperature_gradients); the gravitational-energy (entropy) terms and the
  accretion branch sign conventions match.
- henyey_iterate.f90 driver flow vs crrect.f: convergence test, fcorr ramp,
  triangle rebuild gating, divergence flagging and mix() call placement match
  (except the JENV lifetime issue reported below).
- surfbc.f90: the triangle algebra (three envelope integrations, the linear
  interpolation of log P / log T / log R at the fitting point and the partials
  with respect to log L, log Teff) matches surfbc.f exactly.
- envelope_derivs.f90 / atmosphere_derivs.f90: the ODE right-hand sides
  (dlnP/dlnR, dlnT/dlnP, dtau) and the derivative-of-T(tau) terms match qenv.f /
  qatm.f.
- envint_lib.f90 atm_get driver: cfg snapshot, restart logic, post-return
  application of the gray-fallback flags match the envint.f driver section.
- shell_physics.f90 vs physic.f: interpolation of the shell quantities, the
  Endal-Sofia shear/Richardson criterion, the dropped IDT/IDD dead COMMON,
  and all three call sites (read_starting_model.f90:282, rezone.f90:947,
  evolve_angular_momentum.f90:213) match the 14 dummies.
- compute_scale_height.f90 vs hsubp.f: identical; four call sites in
  overshoot_boundaries.f90 match the 8 dummies.
- temperature_gradients.f90 vs tpgrad.f: MLT cubic (a1, coefficients, the
  Newton/bisection iteration, the `if (iter > 25)` post-loop test), the
  derivative terms, the tgcut switch, and the spot-modification of deldel
  match; the sqrt-of-negative now returns ierr instead of STOP 912.
- turnover_timescale.f90 vs gettau.f + tauintnew.f: the span selection
  (CHKPRS logic), the walker over the stitched model, the one-Hp-above-BCZ
  interpolation, and the JvS 02/2026 change are consistent with the F77; the
  F77 passed an undeclared scalar HS2 as an array, which the modern code
  correctly dropped.
- observables_lib.f90: central-condition extrapolation, nu_max, delta_nu from
  mean density, the delta_Pg = pi^2 sqrt(2) / integral(N/r dr) formula, and
  the moment-of-inertia sums all re-derive correctly.
- eqstat.f90: every path (Saha and fully-ionized) overwrites the
  dlnrho_dlnt / dlnrho_dlnp / cp / grada result slots, so the "historically
  inout" prepacking in henyey_coefficients is harmless.
- rebuild_envelope call site (evolve_step.f90:122-125): 25 actuals match the
  25 dummies, including the log_luminosity / log_luminosity_lsun pair.

## Findings

## core/atmosphere_derivs.f90:61 -- Allard-fallback integration mixes an HSRA T(tau) integrand with an Eddington start point and x-limit
- class: logical
- severity: medium
- confidence: high
- provenance: modernization (envint.f lines 146-153 set COMMON KTTAU=0 *before* the gray integration; modern envint_kernel.f90:260-262 sets only the local `cfg%atm_choice=0`, envint_lib.f90:93 writes `star%job%atm_choice = 0` only after the kernel returns)
- detail: The kernel snapshots star%job into `cfg` (envint_lib.f90:61) and, when the Allard table lookup fails for atm_choice=4, switches `cfg%atm_choice=0` and continues into the gray integration with the Eddington starting temperature (envint_kernel.f90:295-296) and Eddington photospheric x-limit (envint_kernel.f90:391). But the integrand `ttau_log10_temperature(optical_depth, log10_teff, star%job%atm_choice, star%atm_hras)` at atmosphere_derivs.f90:61-62 reads the LIVE star%job%atm_choice, which is still 4 during this call; ttau_lib dispatches 0 -> Eddington, 1 -> Krishna-Swamy, anything else -> HSRA, so the T(tau) integrand is HSRA while the start point and stop condition are Eddington. In the F77 the fallback wrote the COMMON KTTAU that QATM read, so the whole integration was consistently gray. Fires exactly once per Allard-table failure: the one atm_get call that fell out of range (one surfbc triangle vertex, or one photosphere evaluation) integrates an inconsistent atmosphere before the driver flips star%job%atm_choice for subsequent calls. Observable: a discontinuous jump in the fitted envelope (log P / log T at the fitting point) at the model where a KTTAU=4 run leaves the table range, i.e. an envelope-triangle vertex built with a different T(tau) than its neighbours.

## core/henyey_iterate.f90:133 -- envelope_zone_index (JENV) is an uninitialized plain local at corrector levels 1 and 2
- class: numerical
- severity: medium
- confidence: high
- provenance: modernization (crrect.f line 58 blanket SAVE keeps JENV across calls; MIX at crrect.f 103-104 sets it only at ITLVL>2, COEFFT at 121-123 receives it every level; coefft.f:253 `IF(IM.GE.JENV)`)
- detail: `envelope_zone_index` is declared as a plain local at henyey_iterate.f90:133, is written only by the `call mix(...)` at 214-216 (guarded by `iteration_level.gt.2 .and. delta_time.gt.0`), and is passed to henyey_coefficients at 240-242 on every iteration. At levels 1 and 2 (evolve_step.f90:397, 402) it is therefore 0 under gfortran (-finit-local-zero) and garbage under ifort/ifx. henyey_coefficients.f90:300 tests `if (im.ge.envelope_zone_index)` to decide whether to add rot_scr%delta_log_temperature / delta_log_pressure (the accretion entropy-term corrections, set in mdot.f90:179-180 only when mass is being accreted) -- with JENV=0 every zone is treated as belonging to the accreting envelope. In F77 the SAVE'd JENV carried the previous step's value into levels 1-2, which is approximately right; the modern code applies the accretion correction to the whole interior at levels 1-2. Level 3 then recomputes with the correct JENV, so the converged model is (mostly) recovered, but the level-2 solution is what is stored into the entropy-term predictor arrays at evolve_step.f90:407-414 (which feed the next model's first guess), and a bad level-1/2 iterate can trip the divergence flag and force a timestep cut. Only accretion runs (dM/dt > 0) are affected; non-accreting runs see zero deltas and are unaffected.

## core/rebuild_envelope.f90:96 -- comment claims a blanket SAVE that does not exist; omega(old_num_zones) read before set, so new envelope shells get j_rot = 0
- class: logical
- severity: medium
- confidence: high
- provenance: inherited + modernization (getnewenv.f line 46 `SAVE`, line 40 local OMEGA(JSON), lines 246/252 read OMEGA(MM) before this call's GETROT; modern file has no `save` statement at all -- grep finds only the env_*_saved scalars)
- detail: The header comment at rebuild_envelope.f90:92-103 says "Because of the blanket SAVE below, [omega] retains the previous call's values ... Preserved exactly as such", but no `save` exists in the modern file, so the local `omega(json)` is fresh on every call. Lines 352, 354 and 358 read `omega(old_num_zones)` before `call omega_from_j(...)` at 365-367 fills omega(1..num_zones) from specific_angular_momentum. Under gfortran that read is 0, so in the solid-body branch (352-354) every added envelope shell gets omega=0 and j_rot=0, and in the general branch (358) omega_ref=0 so j_rot=0 there too; under ifort/ifx it is garbage. The F77 was also wrong (stale value from the previous call, zero on the first call) but the comment says the modern code reproduces the stale-value behaviour, and it does not: it reproduces the first-call behaviour every time. In addition star%omega is never refreshed by rebuild_envelope (the F77 GETROT also only wrote the local OMEGA), yet the caller evolve_step.f90:122-128 immediately passes star%omega to rotation_shape_factors. Fires on every envelope rebuild in a mass-loss+rotation run where the number of zones changes; observable as a spurious j_rot=0 layer at the top of the model right after each rebuild, later smeared by the angular-momentum transport step.

## core/henyey_coefficients.f90:443 -- dlnepsilon partials and energy components below the nuclear cutoff are stored from stale/zero locals
- class: numerical
- severity: low
- confidence: high
- provenance: inherited by SAVE semantics (coefft.f: ENGEB called only inside `IF(TL.GT.TCUT(1))` at line 213; the derivative and component locals were SAVE'd statics); modern: plain locals
- detail: `zone_dlnepsilon_dlnrho`, `zone_dlnepsilon_dlnt` and `energy_gen_component(1:6)` are only assigned inside the `if (zone_log_temperature.gt.star%ctrl%nuclear_logT_cutoffs(1))` block (245-261), but lines 443-444 (`rot_scr%dlnepsilon_*(im)`) and 459-460 (`star%pulse_dlneps_*(im)`) store them for every zone, and the output branch at 379 tests `.lt.` (so the equality case takes the "else" path and reads them as well). In F77 the statics carried the last hot zone's derivatives into the cold envelope zones; in the modern code they carry the last hot zone's values from this call for zones above the cutoff within a call (same as F77) and zero at the first zone of a call before any hot zone has been seen. The Henyey system itself is unaffected (ql, ql_dt, ql_dp are zeroed at 242-244), so this only pollutes the rot_scr derivative arrays and the pulse-output dlneps columns for cool zones with the last hot zone's values -- the pulse output should carry 0 there.

## core/envint_kernel.f90:483 -- first-step delta_tau uses prev_tau/prev_density/prev_opacity that the dropped F77 "pulse initial point" block used to initialise
- class: numerical
- severity: low
- confidence: medium
- provenance: modernization (envint.f initialised TAUP, OP, DP at the "PULSE INITIAL POINT" block; the corresponding modern block at envint_kernel.f90:382-385 is empty and the comment at 474-479 acknowledges the first-step read)
- detail: `delta_tau_step = (tau_now - prev_tau)/(((density_now*opacity_now)+(prev_density*prev_opacity))/2)` at line 483 is evaluated on the first atmosphere step with prev_* never assigned in this call (zero under gfortran, garbage otherwise). The value goes into the stored atm structure's delta array; stitched_model.f90 skips delta(1) so the pulse output does not see it, but any other consumer of the first-point delta (or a future consumer) would. Under ifort/ifx the garbage could be NaN/Inf and, with -ffpe-trap, trap. Low severity because the only current consumer discards it.

## core/observables_lib.f90:430 -- rotation period divides by omega(nz), which is 0 for a non-rotating start or a fully braked surface
- class: numerical
- severity: low
- confidence: medium
- provenance: inherited (wrtout.f line 402/545 `ROTP = MIN(9999.0D0,0.5D0*C4PI/OMEGA(M)/8.64D4)`)
- detail: `star%rotation_period_days = min(9999.0d0, 0.5d0*c4pi/star%omega(star%nz)/8.64d4)` is guarded only by rotation_active. If omega(nz) is exactly 0 (rotation_active but the starting model has zero surface omega, or after rebuild_envelope zeroes the new surface shells per the finding above) the division yields +Inf, which min() masks to 9999 in a non-trapping build but raises a divide-by-zero trap under the DEV `-ffpe-trap=zero` flags. The rebuild_envelope finding makes this reachable in practice for mass-loss+rotation runs immediately after an envelope rebuild.

## mixing/temperature_gradients.f90:260 -- `convective_velocity = 1.0d0-11` evaluates to -10, not 1e-11
- class: physical
- severity: low
- confidence: high
- provenance: inherited (tpgrad.f line 153 `VEL = 1.0D0-11`)
- detail: The literal is a typo for `1.0d-11` (the same placeholder used elsewhere in the file for "no convection"). The expression parses as 1.0d0 minus 11 = -10.0. The branch is the `else` of the second stability test, which is algebraically unreachable: the second test `g*(-qdt)*presht*deldel .gt. 0` has the same sign as the first test `deldel*g*(-qdt)/presht .gt. 0` that must already have passed (presht > 0), so the code never assigns -10 in practice. Worth fixing because a future change to the surrounding tests would silently produce a negative convective velocity.

## Weak/uncertain observations

- observables_lib.f90:158-159: `core_boundary_fx2 = (grada(j+1)-gradr(j))/(grada(j+1)-gradr(j))` is identically 1 (documented FX/FX2 quirk in the comment at 155; value not consumed).
- observables_lib.f90:284: `call compute_turnover_timescale(star%envelope_radius, ierr)` overwrites envelope_radius with the walker's RBCZ; inherited (wrtout.f 504-505 passed ENVR the same way), but the name now lies.
- envint_kernel.f90:801-802: the inversion-sorting swap moves env_luminosity(idx1) <- env_luminosity(idx2) but never writes idx2 back (one line of the swap missing); harmless because env_luminosity is constant across the envelope.
- rebuild_envelope.f90:266-269: composition extrapolation `x(M) + f*(x(M) - envx)` has the opposite sign convention to the structure interpolation two lines up; inherited (getnewenv.f 179-180) and harmless because envint stores ENVX(NUMENV)=X constant.
- turnover_timescale.f90:239/248: for a fully convective model radius_at_bcz is set to 1. and then exponentiated to 10 cm; inherited quirk (gettau.f), only affects a diagnostic.
- turnover_timescale.f90:365-367: `if (tau .gt. 1.0d20) cycle search` with search_start_index = cz_base_index+1 could in principle loop forever if every candidate CZ returns tau>1e20 (needs conv_vel=1e-11 placeholder, i.e. logT>tgcut at one Hp above the base); inherited from tauintnew's label-50 retry; not reachable for normal models.
- henyey_iterate.f90:313: `dmin1(dabs(elim_rhs(4,i)),dabs(elim_rhs(4,i)/luminosity_lsun(i)))` divides by the zone luminosity, which is zero at the centre and can cross zero in a neutrino-cooled core; inherited (crrect.f); min() masks the Inf but it traps under -ffpe-trap.
- henyey_iterate.f90:160: `start_new_triangle.or.reset_triangle .and.iteration_level.eq.2` -- .and. binds tighter, so start_new_triangle forces a rebuild at any level; same precedence as crrect.f, so inherited intent.
- henyey_coefficients.f90:245 vs 379: the nuclear cutoff test is `.gt.` when computing and `.lt.` when zeroing the output, so a zone with logT exactly equal to the cutoff reports stale energy components; inherited (coefft.f 213/338); measure-zero.
- compute_scale_height.f90: latmo is passed .true. at an interior point (inherited hsubp.f); pscap is dead.
- henyey_coefficients.f90:150-189: `eos_res(i_dlnrho_dlnt) = dlnrho_dlnt` prepacking of the guess into the result array is harmless today only because every eqstat path overwrites it; the log10_density guess slot is the one that is genuinely read.
- NaN-in-MLT symptom (run_from_zahb_to_tahb Solar_m1p0 GN93, model ~950): temperature_gradients is only the messenger. A NaN deldel fails the `.le.1e-6` convergence test, then `test.gt.0` is false and ierr=1 is raised; the NaN must originate upstream in the EOS/opacity inputs (dlnrho_dlnt, cp, opacity) at core He exhaustion. No root cause was found in the assigned files.
- "UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL 1" (0.8 Msun): not in the assigned files (burn/net); nothing in the structure solver feeds that path except the converged T/rho.

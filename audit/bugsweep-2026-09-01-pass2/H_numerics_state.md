# Pass-2 sweep H: numerical library, math utilities, program state

Date: 2026-09-01. Read-only audit; nothing built or run. Originals
consulted at rev 6cd5673 via `git show`.

## Files read in full

- /Applications/YREC/src/numerics/numerics_lib.f90 (3243 lines, in chunks)
- /Applications/YREC/src/numerics/intpar_lib.f90
- /Applications/YREC/src/math/math_lib.f90
- /Applications/YREC/src/state/controls_sync_lib.f90
- /Applications/YREC/src/state/phys_const_lib.f90
- /Applications/YREC/src/state/star_info_lib.f90
- /Applications/YREC/src/state/controls_state_def.inc, job_controls_def.inc
- /Applications/YREC/src/io/controls_lib.f90 (the namelist buffer the sync copies)

Read in part (for provenance / consequence only): src/setup/setups.f90
lines 50-130 (where the phys_const values are assigned),
src/io/read_controls.f90 lines 1-120 and 1500-1510,
src/io/inlist_new_read.inc, src/eos/mhd/mhdpx1.f90, mhdpx2.f90,
src/core/envint_kernel.f90 (bsstep call sites), src/rotation/shape/
rotation_shape_factors.f90 (trapzd/polint/qgauss call sites),
src/core/burn_lib.f90 (ratext call sites), src/rotation/microdiff/
thoul_diffusion.f90 (ludcmp/lubksb), src/rotation/seculr/
diffuse_composition.f90 (osplin argument order).

Originals shown and compared line-by-line: lir.f, locate.f, intpt.f,
meval.f, search.f, cases.f, choose.f, slopes.f, simeqc.f, ysplin.f,
intpol.f, ratext.f, bsstep.f, trapzd.f, qgauss.f, setups.f.

## Verified clean (re-derived or diffed against original/NR)

- cspline / splinj / splinc / kspline / splnr: NR `spline` (natural and
  clamped ends), tridiagonal sweep and back-substitution correct.
- splint / splintd2 / ksplint: NR `splint`, bisection and cubic
  evaluation correct; splintd2 second-derivative form correct.
- inter3 / interp / intrp2: Lagrange weights and their derivatives
  re-derived, correct.
- findex, locate: bisection correct (locate's 0.99/1.01 edge fudge is
  verbatim from locate.f).
- ludcmp / lubksb: NR Crout with implicit pivoting, identical to NR
  including the `tiny` substitution.
- polint: NR Neville tableau, identical (stop -> ierr only).
- mmid: NR modified midpoint, identical to mmid.f.
- ratext: NR `rzextr` (matches ratext.f statement for statement).
- bsstep: matches bsstep.f statement for statement (see finding 3 for
  the one inherited defect); call sites in envint_kernel pass 25
  actuals matching the 25 dummies.
- qgauss: 10-point Gauss-Legendre nodes/weights checked against NR
  (w(2) is truncated at 11 digits, ~1e-13, harmless).
- tridiag_gs / ctridi / tridia: Thomas algorithm correct; tridia's
  `dj_n_seed` argument faithfully replaces the old COMMON smuggling.
- ysplin: slope-form cubic spline re-derived, correct for n >= 4.
- intpol: identical to intpol.f (bracketing, evaluation, derivative).
- simeqc: Gauss-Jordan on the flat 56-element array identical to
  simeqc.f index for index.
- lir / lir1 / lir_impl: the label 2-9 search web was re-traced against
  lir.f in both DIFF>0 and DIFF<0 directions, on-mesh-point branch,
  cubic pivot/closest selection, weight correction, and linear branch:
  equivalent. The SAVEd `N` (search_idx) reset arithmetic is identical.
- meval: the 90/100/110/120 pointer-advance web and the 150-180
  evaluation web were re-traced against meval.f branch by branch:
  equivalent, including the ERR=1/LCN=1 parameter-reuse shortcut and
  the label-210 tail shortcut. Spline parameters are plain locals now
  (were SAVEd) but every path that reads them sets them first.
- choose: the restructured decision tree reproduces every leaf of
  choose.f (labels 9-130), including the two eps-closeness routes to
  label 30/40.
- cases, search, spline (function), slopes: identical to the originals
  (slopes has the 2025 safedivide guards, semantically equal where the
  original did not divide by zero).
- controls_sync_lib: mechanically checked (script over the source) --
  364 assignments in seed_controls_buffer, 364 in store_controls_to_star,
  every seed target has a store source with the SAME name and vice
  versa, no duplicates, no crossed pairs; every declared member of
  controls_lib is synced, none extra; job-member routing (star%job vs
  star%ctrl) consistent in both directions.
- controls_state_def.inc / job_controls_def.inc defaults: every member's
  type and default matches the controls_lib declaration.
- phys_const values assigned in setups.f90: cpi, c, sigma, R, k, h,
  m_e, G (log), ln10 and derived logs all consistent with CODATA-1986
  and with each other; the only defect is finding 4.
- star_info defaults: newton_iterations=0, mixing_length_alpha=1.4,
  flags .false.; read_controls resets star%ctrl via controls_state()
  before seeding, star%job deliberately not reset (mutable run-list
  state). No default/assumption conflict found.
- intpar_lib: tolerance_fraction/max_stage_index/extrap_order have no
  declaration defaults but read_controls.f90:1009 DATA gives
  1.0d-3/11/7 and both inlist paths seed from those.

---

## src/numerics/numerics_lib.f90:1864-1868 -- intpt: loop-exit GOTO became RETURN; MHD table interpolation never runs for in-range T
- class: logical
- severity: high
- confidence: high
- provenance: modernization -- intpt.f lines `DO 100 N=1,NT / IF(TLOG(N).GE.TL) GOTO 101 / IT(1)=N / 100 CONTINUE / 101 IF(IT(1).GE.2) ...` (GOTO 101 exits the loop and continues). Pre-image at 9871628^ still had `goto 101` with `101 if(t_indices(1).ge.2)` after the loop; commit 9871628 ("goto passes 2-5") rewrote it as `if (...) then; continue; return; end if`.
- detail: The first loop scans the table's log T grid for the bracketing index. As soon as `table_log10t(n) >= log10_temperature` -- i.e. for every target temperature that lies inside the table -- the routine now RETURNs without touching `interp_vars` (intent(out)), whereas the original exited the loop and went on to do the 4x4 Lagrange interpolation. Only a target above the whole T grid gets past the loop, and then `t_indices(1)=num_t` is clamped to num_t-3 and the extrapolation proceeds. All eleven callers are in eos/mhd/mhdpx2.f90 (`call zero(interpolated_vars,...)` then `call intpt(...)`), so with the MHD EOS selected every table lookup returns all-zero thermodynamic variables (P, rho, derivatives, ...). Observable consequence: any run using the MHD equation of state is numerically dead (zeros -> divide-by-zero/NaN in eqstat/MHD path); the Stage-0 regression (OPAL/SCV EOS) never exercises intpt, which is why the byte-pinning did not catch it. Secondary: after the fix, `t_indices(1)` is still read uninitialised when `table_log10t(1) >= TL` (the F77 relied on SAVE and a stale value from the previous call) -- clamp it to 1 explicitly.

## src/numerics/numerics_lib.f90:1729-1732 -- trapzd: rho, omega^2 and eta2 are interpolated to a fixed point, not to the running midpoint y
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited -- trapzd.f lines 38-41: `RHOT = RHOP+DRHO*DEL / SMT = SMP+DM*(Y**2 - B1**2) / W2T = W2P + DW2*DEL / ETA22T = ETA22P + DETA2*DEL`.
- detail: In the refinement branch the loop runs over `it` new abscissae `y = b1 + (j-1/2)*del`, and `smt` is correctly evaluated at `y`. But `rhot`, `w2t` and `eta22t` use `drho*del` (constant, the value at b1+del) instead of `drho*(y-b1)`; only `r03t = y**3` and `smt` vary with y. With the caller's jmax=2 (rotation_shape_factors.f90:68) the only refinement is it=1, del=dr, so the single midpoint gets the *endpoint* density, omega^2 and eta2 of shell i while its r^7 and m are midpoint values -- and this inconsistent second trapezoid is then Romberg-extrapolated to h->0 with polint against the first. Observable consequence: a per-shell error of order (rho_i - rho_{i-1})/rho in the integral aint = int (rho/m) r^7 w^2 (5+eta)/(2+eta) dr that feeds <g>, hence fp/ft and the rotational structure corrections in every rotating run; percent-level in the integrand, smaller in the integral because the shell-to-shell density contrast is small, but it is a systematic (one-sided) bias, not noise.

## src/numerics/numerics_lib.f90:1431 -- bsstep: step-underflow test uses the hydrogen fraction instead of the independent variable
- class: numerical
- severity: low
- confidence: high
- provenance: inherited -- bsstep.f: dummy list `(Y,DYDX,NV,X0,HTRY,...,TEFFL,X,Z,KOUNT,KSAHA,ERR)` and line 51 `IF(X+H.EQ.X) THEN`; NR's test is on the integration variable X (here X0).
- detail: Numerical Recipes guards the step-halving loop with `if (x+h .eq. x)` on the integration variable. YREC renamed the independent variable X0 and used X for the hydrogen-fraction pass-through, so the test compares `hydrogen_fraction + h` with `hydrogen_fraction`. The integration variable is log P (or log tau) of order 1-10 while X is 0.7 (or 0 in a He-rich envelope), so the test fires at a different, and in the X=0 case only at the true denormal-underflow, point; the retry loop can therefore spin ~30 more quartering iterations than intended before returning ierr. No wrong result, but a delayed/mis-scaled failure detection in envint's atmosphere/envelope integrations.

## src/setup/setups.f90:118 (value of phys_const hydrogen_atom_mass) -- Avogadro digits transposed
- class: physical
- severity: low
- confidence: high
- provenance: inherited -- setups.f line 121 `HMASS = 1.0D0/6.0222137D23`.
- detail: The comment says "mass hydrogen atom (gm)" computed as 1/N_A, but the literal is 6.0222137e23; CODATA-1986 N_A is 6.0221367e23 (digits '13' and '22' swapped), a relative error of 1.3e-5. The value feeds only `debye_huckel_coefficient` (sqrt(pi/(k m_H^3)) e^3/3), so the Debye-Hueckel free-energy correction is off by ~2e-5 relative. Negligible for evolution, but it is a wrong constant and contradicts its comment; the rest of setups uses correctly-transcribed CODATA-1986 values.

## src/numerics/numerics_lib.f90:2051-2053 -- lir: degenerate table (z(1)==z(n)) returns with result_y unassigned
- class: interface
- severity: low
- confidence: high
- provenance: inherited -- lir.f `IF(DIFF) 4,102,3` ... `102 CONTINUE / RETURN` (the diagnostic WRITE at FORMAT 1002 is unreachable dead code there too).
- detail: `result_y` is intent(out) and the DIFF==0 path returns without writing it or `interp_flag`; the original had the same silent return (the FORMAT 1002 diagnostic was never issued). Callers (mhdpx1.f90:140 cubic-in-X branch, intpt) use the output unconditionally. Fires only when the 4 selected X nodes are identical (a badly built MHD table), so consequence is garbage rather than a crash in a rare path -- reporting because the modern intent(out) contract makes the silent path formally undefined.

---

## Weak/uncertain observations

- numerics_lib.f90:64-117 boole: `scalex = 1e-11`, `scaley = 1e7` single-precision literals (cancel on rescale) and the (n_grid-1)/4 truncation drop trailing points; no caller anywhere in src -- dead code.
- numerics_lib.f90:1099-1167 ysplin (inherited, ysplin.f identical): for n=3 the `i=2` and `i=n-1` boundary rows are the same row, second overwrites first with the wrong constants; callers (intpol from kurucz/opal92 with num_valid_temps/jt) are believed to always pass >= 4 points but nothing enforces it.
- numerics_lib.f90:1475/1482 intpol and ysplin: fixed scratch `spline_coeff(4,100)`, ysplin's F/H/D(100); no guard that n_grid <= 100 (opal92/kurucz temperature counts are < 100 today).
- numerics_lib.f90:822 polint declares `xa(20), ya(20)` explicit-shape; alsurfp.f90:241-257 passes array elements `allard_gl_grid(j1)` / `allard_teffl_grid(teffl_index)` by sequence association -- fine for n=4, but non-conforming if fewer than 20 elements remain and would trip -fcheck=bounds.
- numerics_lib.f90:2204 ratext reads `tableau(var_idx,k2)` at k2=num_use before it is ever written on the first pass (NR's rzextr guards with `IF(K.NE.NUSE)`); the value is discarded, so harmless. Also `fx(ncol=7)` / `substep_sequence(11)` are not checked against the namelist nuse/imax (defaults 7/11 are the maxima; any larger user value overruns).
- numerics_lib.f90:513-518 locate (inherited): the 0.99/1.01 out-of-range fudge assumes positive abscissae; for negative log tables the window is on the wrong side (simply never fires) -- all current callers pass positive grids.
- numerics_lib.f90:2254 simeqc (inherited): `write (5,1011)` writes the singular-matrix message to unit 5 (stdin) before returning ierr=1.
- numerics_lib.f90:2821-2843 search (inherited): with num_table_points==1 and no exact match `last-first` never equals 1 -> infinite loop; meval only calls it for >= 2 table points in practice.
- numerics_lib.f90 tridia (~1300): `fj = 1+(solution(n)-rhs_orig(n))/rhs_orig(n)` divides by the last RHS entry, which can be zero for a zero-flux boundary.
- numerics_lib.f90:1882-1884 intpt (inherited): P above the row maximum returns silently leaving the caller's zeroed output; there is no out-of-range flag to the MHD driver.
- math_lib.f90:90-93,113-116 crmath backend only: `pow_r(0d0, 0d0)` returns 0 (intrinsic gives 1) and `pow_r(0d0, -n)` returns 0 (intrinsic gives +Inf); the two backends therefore differ in semantics, not just last bits, at x==0. No current call site was found that passes a zero base with a non-positive exponent, so this is latent.
- src/defaults/controls.defaults:207/372/388 and controls_registry.tsv label imax/nuse/stolr0 as "0 (implicit)" defaults; the actual defaults are 11/7/1e-3 from read_controls.f90:1009 DATA -- documentation/registry mismatch only (both inlist paths seed correctly).
- phys_const_lib.f90:20 `cmkhn` is declared public, never assigned and never read anywhere in src -- dead.
- setups.f90:70 `star%ln_solar_luminosity = ln10/star%solar_luminosity_cgs` is ln(10)/Lsun, not ln(Lsun); the only consumer (henyey_coefficients.f90:341) uses it as the former, so misnamed, not wrong (inherited CLNSUN = CLN/CLSUN).
- numerics_lib.f90:1970 lir comment "MOST OF THE COMPUTATION IS PERFORMED IN SINGLE PRECISION" is stale (all double).

# Bug sweep: src/numerics/, src/setup/, src/state/

Files read in full: numerics/intpar_lib.f90, numerics/numerics_lib.f90 (3243 lines);
setup/setups.f90, locate_shell_boundaries.f90, map_user_inputs.f90, rescale_model.f90,
rezone.f90, rotation_stability_setup.f90, setup_solar_calibration.f90,
setup_star_calibration.f90; state/star_info_lib.f90, controls_sync_lib.f90,
phys_const_lib.f90 (+ controls_state_def.inc, job_controls_def.inc).
Provenance checked against git 6cd5673 originals (intpt.f, bsstep.f, mmid.f, trapzd.f,
lir.f, meval.f, choose.f, slopes.f, simeqc.f, findex.f, boole.f, hpoint.f, findsh.f,
rscale.f, remap.f, setupv.f, setcal.f). The controls_sync_lib seed/store pair was
verified programmatically to be exactly symmetric (364 mirrored assignments, no
mismatched targets).

## numerics_lib.f90:1864-1870 -- intpt: loop-exit GOTO converted to RETURN, killing the whole MHD-EOS interpolation path
- class: logical
- severity: high
- confidence: high
- provenance: modernization (original intpt.f: `DO 100 N=1,NT / IF(TLOG(N).GE.TL) GOTO 101 / IT(1)=N / 100 CONTINUE / 101 IF(IT(1).GE.2)...` -- GOTO 101 exits the loop and CONTINUES the interpolation)
- detail: The temperature-bracket search reads
  `do n=1,num_t; if (table_log10t(n).ge.log10_temperature) then; continue; return; end if; t_indices(1)=n; end do`.
  The original's `GOTO 101` (jump past the loop, then clamp t_indices(1) and do the
  4x4 Lagrangian interpolation) was turned into a `return`. Consequence: whenever the
  requested temperature is inside the table range -- i.e. in every normal call -- intpt
  returns immediately with `interp_vars` (intent(out)) never assigned; interpolation
  only ever runs in the off-table extrapolation case where TL exceeds the entire
  TLOG grid. intpt is the interpolation engine for all 11 call sites in
  eos/mhd/mhdpx2.f90 (ZAMS A/B/C and CENTRE1-5 tables), so any run with use_mhd_eos
  gets undefined EOS values. The file's own header admits this routine is "NOT
  exercised by that suite ... verified by build + code review only", which is how it
  survived Stage-0. The fix is `exit`, not `return`. (Secondary, inherited: the
  in-range `if(log10_pressure.gt.p_max) return` at line 1882 also leaves interp_vars
  unset, and t_indices(1) is undefined if table_log10t(1) >= TL -- the original
  papered over that with SAVE staleness.)

## numerics_lib.f90:1431 -- bsstep: step-underflow guard tests the hydrogen fraction, not the independent variable
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited (bsstep.f: signature `BSSTEP(Y,DYDX,NV,X0,HTRY,...,RL,TEFFL,X,Z,KOUNT,KSAHA,ERR)` -- X0 is the independent variable, X the pass-through hydrogen fraction -- yet the guard is `IF(X+H.EQ.X)`)
- detail: Numerical Recipes' BSSTEP guards against step underflow with
  `if (x+h .eq. x)` where x is the integration variable. When YREC added its opaque
  pass-through arguments, the independent variable became X0 but the guard kept
  testing X, which in YREC's calling convention is the envelope hydrogen fraction
  (~0.7). The modernized code faithfully reproduces this:
  `if(hydrogen_fraction+h.eq.hydrogen_fraction)`. The termination criterion therefore
  fires when h underflows relative to ~0.7 instead of relative to |indep_var| (log P
  during atm/envelope integration, typically O(1-10), but the scales differ and the
  sign of the effect varies). Failure mode when the scales diverge: either the
  step-halving loop runs long past the point where x0+h stops advancing (bsstep can
  then "accept" a step that leaves indep_var bit-identical, hanging the caller's
  integration loop), or it aborts earlier than NR intends. This is exactly the class
  of "convergence test vs standard form" defect the sweep targets.

## numerics_lib.f90:1726-1738 -- trapzd: shell-interpolation of rho/omega^2/eta2 is loop-invariant (uses del where y-b1 belongs)
- class: physical
- severity: medium
- confidence: high
- provenance: inherited (trapzd.f has the identical `RHOT = RHOP+DRHO*DEL`, `W2T = W2P + DW2*DEL`, `ETA22T = ETA22P + DETA2*DEL` inside the DO 10 loop)
- detail: In the refinement branch (n>1), the comment says "interpolate rho,m,omega,
  eta2+2 between shell i and shell i-1", and the mass indeed uses the midpoint
  ordinate (`smt = smp+dm*(y**2-b1**2)`), but rho, omega^2 and eta2+2 use
  `...p + slope*del` -- del = (b2-b1)/it is the current sub-step width, constant over
  the loop and shrinking with each refinement level. So every new midpoint at every
  refinement level gets the SAME rho/w2/eta2 (the b1-endpoint value plus a vanishing
  offset), instead of values interpolated to y. The correct factor is (y-b1). The
  refinement therefore converges, but to the integral of the wrong integrand (interior
  values pinned to the inner-shell endpoint); the trapezoid-refinement error estimate
  used by the caller (quad's aint for the rotating-shape <g> integrals) is computed
  against a systematically first-order-in-(b2-b1) biased integrand. Affects rotation
  shape factors/mean gravity accuracy in rotating models; exact for rigid uniform
  shells, biased where rho/omega/eta2 vary across a shell.

## setup/map_user_inputs.f90:133-190 -- S'/S and S''/S rate-derivative scales normalized by the pp S0 instead of each reaction's own A98 S0
- class: physical
- severity: medium
- confidence: medium
- provenance: inherited (remap.f: `QS0E(2) = (S0P_3_3/S0_3_3)/(QS0E_3_3_A98/S0_1_1_A98)` etc.; only QS0E(1) and QS0E(8)/QQS0EE(8) divide by their own S0_*_A98)
- detail: The intended scale is (S'/S)_user / (S'/S)_A98 per reaction; engeb/nrates
  multiply their hard-coded FCZ Seff Q2/Q3 (and Q4/Q5) terms by these scales
  (core/burn_lib.f90:1359, net/net_lib.f90:2059). For reactions 2-7 (he3he3, he3he4,
  p+C12, p+C13, p+N14, p+O16) and the second derivatives 4,5,7, the reference ratio
  in the denominator is built as qs0e_X_a98/s0_pp_a98 -- the A98 derivative of
  reaction X over the A98 S0 of the *pp* reaction (4.0e-22 keV b). The commented-out
  DATA block in the same file carries the correct per-reaction S0_*_A98 values
  (5.4e3, 5.3e-1, 1.34, ...) which the "KC 2025-05-31" cleanup deleted as unused.
  Numerically the computed scale is the intended one times S0_pp_A98/S0_X_A98
  (~1e-22 for he3he3, ~3e-22 for p+C12, ...), i.e. the user's S'-correction terms are
  effectively multiplied by zero for reactions 2-7 whenever use_new_nuclear_rates
  (LNEWNUC) is on -- the Seff derivative corrections silently vanish instead of
  scaling. The asymmetry with pp (index 1) and be7+p (index 8), which use their own
  S0, is strong evidence this is a defect rather than a convention; confidence is
  medium only because the namelist S0P inputs' intended normalization is documented
  ambiguously ("relative to the Solar Fusion I ... ratios") while example inputs use
  absolute S0 values.

## setup/rezone.f90:534-551, 590-604 -- gradient-driven point insertion has no json bound check
- class: logical
- severity: medium
- confidence: high
- provenance: inherited (hpoint.f: the `DO K = MM+NN,J+NN,-1 / HSS(K) = HSS(K-NN)` shift loops have no JSON guard either; only the chi-segment assignment checks `MNEW+NPT.GT.JSON`)
- detail: The chi-based point assignment stops with ierr when
  new_num_zones+segment_point_count > json (line 460), but the subsequent X-gradient
  and Z-gradient insertion passes shift `star%old_shell_mass` up by num_new_points per
  flagged jump (`do k = working_num_zones+num_new_points, j+num_new_points, -1`)
  without any check that working_num_zones+num_new_points <= json=5000. A model near
  the 5000-zone budget with steep X or Z gradients (post-dredge-up composition steps,
  pre-flash) writes past the end of old_shell_mass -- in the modernized layout that is
  a component of the module-level `star` structure, so the overrun corrupts the
  adjacent star_info member rather than trapping. The 2026 audit fixed the analogous
  flag_point overflow (known finding) but not this one.

## setup/rescale_model.f90:144 -- error message prints `icomp`, which is never assigned
- class: logical
- severity: low
- confidence: high
- provenance: inherited (rscale.f line 115: `WRITE(ISHORT,1004)ICOMP,NK,XNEWCP` with no assignment to ICOMP anywhere in the file)
- detail: In the single-species rescale error path (requested abundance >= 1), the
  log line "RESCALING OF SPECIES i3 IN KIND CARD #..." prints local `icomp`, which is
  declared but never set (the meaningful value would be
  star%ctrl%new_species_index). In the modernized version icomp is an ordinary
  uninitialized local (the F77 original at least had SAVE), so the report prints
  garbage for the species index on the diagnostic a user most needs to read. Error
  path only; no effect on evolution.

## numerics_lib.f90:2158-2166 -- ratext tableau/x_hist bounds tied to namelist extrap_order/max_stage_index with no guard
- class: logical
- severity: low
- confidence: medium
- provenance: inherited (ratext.f had the same IMAX=11/NMAX=15/NCOL=7 fixed dimensions fed by the /INTPAR/ namelist values)
- detail: `x_hist(11)` is indexed by est_index (= bsstep's stage index, up to
  max_stage_index) and `tableau(15,7)`'s second index runs to
  num_use = min(est_index, max_use), where max_use = extrap_order; both
  max_stage_index and extrap_order come straight from NAMELIST /physics/ (intpar_lib)
  with no validation anywhere on the path. A user setting extrap_order > 7 or
  max_stage_index > 11 gets silent out-of-bounds writes into SAVEd module state
  rather than an input error. With the shipped defaults (11/7) everything fits
  exactly; this only fires on nonstandard input.

## numerics_lib.f90:512-517 -- locate's out-of-range patch assumes a positive-valued ascending table
- class: numerical
- severity: low
- confidence: medium
- provenance: inherited (locate.f carries the same `IF ((J.EQ.0).AND.(X.GT.0.99D0*XX(1)))` / `X.LT.1.01D0*XX(N)` patches on the NR bisection)
- detail: The post-bisection tolerance patch re-admits x slightly outside the table by
  comparing against 0.99*xx(1) and 1.01*xx(n). That multiplicative tolerance has the
  intended sense only when xx(1), xx(n) > 0: for a table of negative values (e.g. log
  quantities below zero), 0.99*xx(1) > xx(1), so the "just below the low edge" case is
  never re-admitted while values genuinely inside the first 1% band above a negative
  xx(1)... are (the test degenerates to admitting in-range j==0 results that cannot
  occur). Harmless for the positive-valued tables it currently serves; a trap for any
  new caller with signed abscissas.

## numerics_lib.f90:2253-2254 -- simeqc singular-matrix diagnostic writes to unit 5 (stdin)
- class: logical
- severity: low
- confidence: high
- provenance: inherited (simeqc.f: `1010 WRITE (5,1011)` before `STOP29`)
- detail: The singular-pivot error message goes to Fortran unit 5, conventionally
  connected to standard input; on modern runtimes this typically lands in a stray
  fort.5 file or errors, so the "STOPPED AT 1010" breadcrumb is lost exactly when the
  composition solve (mixing/solve_composition) hits a singular system. The 2026 ierr
  conversion preserved the misdirected write. One-line fix: run_log/terminal unit.

## Weak/uncertain observations
- rezone.f90:258-317 (inherited, hpoint.f identical): pmax1-pmax5, including the ceiling-divided overshoot-region spacings pmax2/pmax3 and the appended flag points' companion spacing logic, are computed but never used -- the 7/02 chi-loop rewrite reads chi_grid_scale(8/10/11) directly, so the finely-zoned overshoot spacing the comment block describes is dead code; the flag-point appends for overshoot_base_zone/fine_zone_base still act.
- rezone.f90:546/600 (inherited): the insertion fill loop runs down to k=j-1, rewriting the pre-existing lower endpoint old_shell_mass(j-1) with an arithmetically-identical-but-rounded value; harmless roundoff churn.
- rezone.f90:609-625 (inherited): the too-close-point deletion pass does not protect flagged points (CZ edges, H-shell edge), so a flag point can be deleted if it lands within chi_grid_scale(1) of its neighbor; also the final point is always kept regardless of spacing.
- rezone.f90:989 (inherited): the surface-opacity-table refresh triggers only on a surface-X change > 1e-8; a Z-only surface change skips kap_update_surface_tables.
- numerics_lib.f90:590-618 (inherited from NR): ludcmp's scaling vector vv is fixed at 100 with no n<=100 check; the only current caller (thoul_diffusion, n<=42) fits.
- numerics_lib.f90:94-101 (inherited): boole silently ignores the trailing (n_grid-1) mod 4 grid points when n_grid is not 1+4k; caller contract only, no internal check.
- numerics_lib.f90:1707/1740 (inherited): trapzd's refinement counter `it` is SAVEd module state shared by all integrations -- two interleaved trapzd-based integrations would corrupt each other; single-threaded single-integral use today.
- numerics_lib.f90:2031-2033 (inherited): lir/lir1 with num_points<2 or a degenerate table (diff==0) returns with interp_flag (intent(out)) never set -- callers reading it get an undefined value; original left INTER equally unset.
- numerics_lib.f90:3067-3070 (inherited): meval with every eval point below the table calls search with eval_x(num_eval_points+1) -- an in-bounds (json-dimensioned) but semantically stale element; result discarded on the normal path.
- rescale_model.f90:366-386 (inherited): z-ramp divides by composition(3,zone) with no Z>0 guard, and z_ramp_slope divides by (rsclzm2-rsclzm1) with only >0 checks on each individually (equal values give a 0/0 slope).
- rescale_model.f90:135-136 (inherited): when value_relative_to_h is set, star%job%new_species_value is overwritten in place with the converted absolute abundance, so a second kind card re-converts an already-converted value.
- state/star_info_lib.f90 (documented in-file, known-adjacent): old_shell_mass despite its name holds pre-rezoning log-mass scratch, and scp/cp are same-quantity-different-fill-time aliases pending a semantics audit -- both flagged in the header comments; no new dimension inconsistency found (mixed 12 vs radiative 13 zone-bounds rows match the F77 MXZONE/MRZONE dims; the generated ctrl/job sync is exactly symmetric).
- numerics_lib.f90:66-67 (inherited, boole.f `PARAMETER(SCALEX=1e-11)`): single-precision e-notation scale constants -- harmless here because the identical constant multiplies in and divides out exactly; same literal-precision class as the already-reported burn/net/turnover findings.
- setup/rotation_stability_setup.f90: line-by-line comparison against setupv.f found a faithful conversion (interface Lagrangian weights, DDEL floors, FACT1-6, dynamical/diffusive shear coefficients, LDIFAD block all match); the FACT5 sign convention (CON-1 numerator against an outward-negative dlogP) is odd-looking but exactly the original's.
- locate_shell_boundaries.f90, setups.f90, setup_solar_calibration.f90, setup_star_calibration.f90, controls_sync_lib.f90, phys_const_lib.f90, intpar_lib.f90: no defects found beyond the above; constants in setups.f90 check out against CODATA-86-era values the code has always used (csig 5.67051e-5, G default log10 = -7.17571 for 6.6726e-8, etc.), consistent with the historical bit-pinning.

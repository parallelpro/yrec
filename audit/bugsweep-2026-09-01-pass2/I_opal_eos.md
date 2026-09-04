# Sweep I (pass 2): OPAL 1995 / 2001 / 2006 EOS interpolation packages

Scope: /Applications/YREC/src/eos/opal/*.f90 (27 files) and
/Applications/YREC/src/eos/opal_eos_lib.f90. F77 ancestors at 6cd5673
(src/esac.f, esac01.f, esac06.f, eqbound*.f, gmass*.f, oeqos*.f,
quad*.f, radsub*.f, readco*.f, rhoofp*.f, t6rint*.f).

## Files read in full

- src/eos/opal_eos_lib.f90
- src/eos/opal/esac.f90, esac01.f90, esac06.f90
- src/eos/opal/eqbound.f90, eqbound01.f90, eqbound06.f90
- src/eos/opal/gmass.f90, gmass01.f90, gmass06.f90
- src/eos/opal/oeqos.f90, oeqos01.f90, oeqos06.f90
- src/eos/opal/quad.f90, quadeos01.f90, quadeos06.f90
- src/eos/opal/radsub.f90, radsub01.f90, radsub06.f90
- src/eos/opal/readco.f90, readcoeos01.f90, readcoeos06.f90
- src/eos/opal/rhoofp.f90, rhoofp01.f90, rhoofp06.f90
- src/eos/opal/t6rinterp.f90, t6rinteos01.f90, t6rinteos06.f90
- Call sites / context (excerpts): src/eos/eqstat.f90:585-810,
  src/eos/eos_lib.f90:380-475, src/eos/mu.f90:17-60,
  src/eos/test/test_eos.f90:310-325, src/io/read_controls.f90:1865-1895,
  src/Makefile:85-171 (compiler flags).
- Originals compared: esac.f, esac06.f, eqbound.f, eqbound06.f, gmass.f,
  gmass01.f, gmass06.f, oeqos.f, oeqos06.f, radsub.f, radsub01.f,
  radsub06.f, readco.f, readcoeos01.f, readcoeos06.f, rhoofp.f,
  rhoofp01.f, rhoofp06.f, t6rinterp.f, t6rinteos06.f (plus grep-level
  checks of esac01.f, quad*.f, t6rinteos01.f).

## Verified clean (re-derived, found correct)

- quad/quadeos01/quadeos06: Lagrange-through-Newton quadratic
  coefficients (c2, c1, c0) and the recompute_flag/cache_slot scheme;
  slots are only reused with the same abscissae.
- X / T6 / rho bracketing searches in all three esac generations
  (descending T6 grid handled correctly; result k3 >= 3, l3 >= 3,
  x_index_lo >= 1), and the interpolation weights
  dix = (t6(k3)-slt)/(t6(k3)-t6(k2)), dix2 = (rho(l3)-slr)/(rho(l3)-rho(l2)),
  x_interp_weight = (x(x3)-X)/(x(x3)-x(x2)): each is 1 when the point
  sits at the centre of the first quadratic's stencil and 0 at the
  centre of the second, i.e. the overlapping-quadratic smoothing is
  correctly oriented in all three axes.
- The 2001/2006 edge guards (ipu/iqu pre-checks, `l3==nr .or. k3==nt`
  reset, x_index_lo reset for X<1e-6) and the goto-66 (Z mismatch)
  emulation via `exit` + `if (eos_var_idx > deriv_order)`.
- readcoeos01/readcoeos06 formats and variable counts (9 columns each;
  2001 layout 4=cv,5=chi_rho,6=chi_T,7=gamma1,8=gamma2/(gamma2-1);
  2006 layout 4=dE/drho,5=cv,6=chi_rho,7=chi_T,8=gamma1,9=gamma2/(gamma2-1)),
  and the nta/nra DATA lists moved into opal_eos_lib (element-by-element
  identical to esac06.f / rhoofp01.f / rhoofp.f).
- radsub family physics: P_rad = 4/3 (sigma/c) T^4 with sigma/c in
  Mb/T6^4 = 1.8914785e-3 (correct), E_rad = 3P_rad/rho, S_rad = 4/3 E_rad/T6,
  chi_rho = chi_rho,gas P_gas/P, chi_T = (P_gas chi_T,gas + 4P_rad)/P,
  cv = cv,gas + 4E_rad/T6, gamma3-1 = P chi_T/(rho T cv),
  gamma1 = chi_rho + chi_T (gamma3-1), gamma2/(gamma2-1) = gamma1/(gamma3-1)
  (Cox & Giuli 9.81-9.88). 2006 correctly skips the non-existent
  gamma3-1 column.
- oeqos/oeqos01/oeqos06 conversions: dlnrho/dlnP = 1/chi_rho,
  dlnrho/dlnT = -chi_T/chi_rho, cp = 1e6 cv gamma1/chi_rho
  (cp/cv = Gamma1/chi_rho; Mb-cc/(g T6) -> erg/g/K), grad_ad =
  1/(gamma2/(gamma2-1)); beta14 = 2.521971383e-3 T6^4/P12 equals
  (a/3)T^4/P and is consistent with radsub's 4/3*1.8914785e-3. The
  9-argument `call mu(...)` matches src/eos/mu.f90's dummy list.
- gmass/gmass01/gmass06: mole and mass bookkeeping, atomic weights,
  ionization energies all match the originals; only one OPAL year can be
  active per run (read_controls.f90:1875-1885), so the shared
  table_loaded_flag is safe.
- eqbound06 ramp factors (0 at the table edge, 1 one cell inside) and the
  linear row/column search seeded from esac06's k3/l2.

## Findings

## src/eos/opal/esac.f90:117 -- OPAL 1995 EOS is dead: goto-61 conversion made ierr=1/return unconditional on the table-loading call
- class: logical
- severity: high
- confidence: high
- provenance: modernization (esac.f: `IF(Z+XH-1.D-6 .GT. 1.0D0 ) GO TO 61`;
  converted in commit 9871628 "goto passes 2-5 ... all automated, all
  gated" to the block shown below, without THEN/END IF)
- detail: Lines 117-123 read `if (...) & write(...)` followed by
  un-nested `write`, `ierr = 1`, `return`. Only the first write is
  governed by the logical IF; the ierr=1/return executes unconditionally
  inside the `table_loaded_flag` first-call block. So the very first
  esac() call of any run with LOPALE=.TRUE. (use_opal95_eos) prints the
  Z/X values, returns ierr=1, and every rhoofp/oeqos caller propagates
  the error; the 1995 OPAL path can never produce a value. The gate that
  passed the automated goto pass only exercises LOPALE06 (all 19 example
  configs), which is why this was not caught. esac01.f90:122 and
  esac06.f90:122 kept a proper `if ... then` block. The same automated
  pattern was checked for in the other 27 files; this is the only
  instance.

## src/eos/opal/esac06.f90:406 (also esac.f90:382, esac01.f90 analogue, radsub*.f90) -- stale cv/gamma rescaling on every deriv_order=1 trial call; root cause of the "eos(5) growing without bound" crash that rhoofp06's tolerance change papers over
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited (esac06.f:475-477 `if (irad .eq. 1) call radsub06(...)
  else eos(iri(5))=eos(iri(5))*moles*aprop/tmass endif` after the
  `do 124 iv=1,iorder` loop, unconditional in iorder; radsub01.f/radsub.f
  likewise rewrite eos(2..9/10) from whatever is in the array)
- detail: rhoofp* call esac* with deriv_order=1 up to 14 times per
  inversion, so only eos(1) is freshly interpolated, yet the post-loop
  code always rescales cv (eos(5), or eos(4) for 2001) by
  moles*R/tmass (~1e2) and, for rad_flag=1 (1995/2001 trial calls),
  radsub/radsub01 recompute chi_rho, cv, gamma1..3 from the already
  rescaled values. The factor is only undone by the next full-order call
  in oeqos*. If rhoofp* fails (11-iteration cap, or an out-of-table trial
  point taking the *999 return) oeqos* takes `return 1`, eqstat silently
  keeps the Yale values (eqstat.f90:808) and no reset happens; ~11
  consecutive failing shells push cv to +Inf, and in radsub01/radsub
  gamma3-1 -> 0 and gamma2/(gamma2-1) = gamma1/0. With the DEV flags
  (`-ffpe-trap=invalid,zero,overflow`, src/Makefile:156) that aborts the
  run -- exactly the symptom recorded at rhoofp06.f90:174-176 ("eos(5)
  growing without bound and crashing"). The correct fix is to apply the
  rescaling/radsub only to the entries actually interpolated
  (deriv_order >= 5), not to loosen the inversion tolerance.

## src/eos/opal/rhoofp06.f90:176 -- convergence tolerance loosened to 1e-5 (0.5e-7 in the original and in rhoofp/rhoofp01): rho(P,T) from OPAL 2006 now carries ~1e-5 relative solver noise
- class: numerical
- severity: medium
- confidence: medium
- provenance: modernization-era divergence (rhoofp06.f:99 `IF (DABS((P3-PNR)/PNR) .LT. 0.5D-7)`;
  changed by commit 1bc94ae 2025-10-10; rhoofp.f90:183 and
  rhoofp01.f90:184 still use 0.5d-7)
- detail: The regula-falsi stop test is the only thing that fixes the
  returned density; at 1e-5 in P (chi_rho ~ 1) the density handed to
  eqstat has an iteration-dependent error of up to ~1e-5 relative, which
  is not a smooth function of (P,T). oeqos06's P-consistency check
  (`abs((p_e12-eos(1))/p_e12) > 0.5d-6`, still active in oeqos.f90:86) is
  commented out, so nothing detects it. Since every example config uses
  LOPALE06, this sets the noise floor of the density (and hence of the
  hydrostatic residual and its Henyey corrections) in all production
  runs at ~1e-5, 200x worse than the other two generations, and is a
  cross-generation divergence not explained by a table-format
  difference. The change was made to avoid the overflow described in the
  previous finding; once that is fixed the original 0.5d-7 can be
  restored.

## src/eos/opal/esac.f90:208,233-243 + t6rinterp.f90:71-90 -- OPAL 1995 has no low-T-edge guard: at k3 = nt the 3x4/4x4 table sums read row nt+1, ip=3 is selected, and t6rinterp uses t6_grid(nt+1) and x_interp_result(nt+1,:)
- class: logical
- severity: medium (OPAL95 only; unreachable until the first finding is fixed)
- confidence: high (out-of-bounds), medium (size of the numerical damage)
- provenance: inherited (esac.f:181 `K4=K3+1`, :205-210 `DO IT=K1,K1+3`,
  :236 `IF(L3 .EQ. NR) IQ=2` with no NT counterpart; esac01.f90's header
  documents the `k3==nt` reset as a 2001 addition)
- detail: The T6 bisection returns k3 = nt whenever t6 lies in
  [t6_grid(nt), t6_grid(nt-1)], which eqstat's OPAL95 gate
  (temperature >= 5000 K) admits. table_sum_3x4/4x4 then include
  eos_table(x,1,nt+1,ir) -- in memory the pressure of row 1 of the next
  density column, a finite number -- so t6_interp_order becomes 3
  (line 262), the X-interpolation loop writes x_interp_result(nt+1,ir)
  and t6rinterp evaluates the second quadratic with the abscissa
  t6_grid(nt+1), i.e. whatever follows t6_grid in the derived type. The
  blend weight (1-dix) of that garbage quadratic is 1 at the table's
  coolest row, so OPAL95 output in the lowest-T band is wrong (or trips
  the >1e15 ierr). Related inherited out-of-bounds reads in the same
  routine: line 219 `x_index_lo, x_index_lo+3` (esac.f:193 `DO M=MF,MF+3`)
  reads eos_table(6,...) for X > 0.6, and lines 238-243 read density
  column nr+1 when l3 = nr (harmless because line 268 then forces
  iq=2). 2001/2006 avoid all of these (`x_index_lo, x_index_hi`, ipu/iqu,
  `k3==nt`). A DEV build (`-fcheck=all`) would abort at the first one.

## src/eos/opal/t6rinterp.f90:78-90 -- comment says dix is carried over from a previous call via SAVE; it is now a plain local, so the density-only smoothing branch uses an uninitialized weight
- class: numerical
- severity: low
- confidence: high
- provenance: modernization (t6rinterp.f has `SAVE`; the 2026 "Save
  migration: eos and kap are fully save-free" commit 5cadb1e removed it
  and the comment was not updated)
- detail: When density_interp_order=3 but t6_interp_order/=3, line 90
  `esactq = esactq*dix + esactq2*(1-dix)` reads dix, which is never
  assigned in this call (zero under `-finit-local-zero`, arbitrary with
  another compiler). The result esactq is only consumed at line 98 inside
  `if (t6_interp_order.eq.3)`, so the value is dead today and the output
  is unaffected -- but the comment's stated intent (stale-but-defined
  SAVEd value) is false, and any future use of esactq in the iq=3/ip=2
  case would be garbage. The 2001/2006 versions nest the block inside
  the ip==3 test ("moved to loc a"), which is the correct structure.

## src/eos/opal/radsub.f90:76-82 -- 1995 unit revision applied to P and S but not to E, cv, dE/drho
- class: physical
- severity: low
- confidence: medium
- provenance: inherited (radsub.f: `REVISE=UNITF/UNITFOLD; PT=PT*REVISE;
  ST=ST*REVISE` only; the MHP 10/02 comments about EN/DEDRHOA already
  question it)
- detail: unit_ratio = 0.9648575/0.9652 = 0.99964 rescales the total
  pressure and entropy returned to oeqos, but the energy, cv and dE/drho
  stay in the old units and the chi/gamma set was computed from the
  unscaled pressure. For OPAL95 this makes P and cv mutually inconsistent
  at the 3.6e-4 level (cp via 1e6*cv*gamma1/chi_rho, and beta via
  P_rad/P, see the same P). The 2001/2006 radsubs dropped the revision
  entirely. Not reachable until the first finding is fixed.

## src/eos/opal/eqbound06.f90:94 (and eqbound.f90/eqbound01.f90 row search) -- index 0 read when the density equals the first grid point exactly
- class: logical
- severity: low
- confidence: medium
- provenance: inherited (eqbound06.f: `IF((D.LT.rho(1)) ...) GOTO 9999`
  then `40 IF(D .LE. rho(l-1)) l=l-1`)
- detail: The out-of-table test rejects density < density_grid_06(1)
  but admits equality; the downward search `do while
  (density.le.density_grid_06(density_row-1))` then decrements to
  density_row=1 and evaluates density_grid_06(0). Measure-zero for a
  computed 10**log10_density, but a DEV build would abort on it, and
  the 1995/2001 eqbound have the mirror-image edge (t_row_index=nt when
  t6 == t6_grid(nt), then t6_grid(nt+1)).

## src/eos/opal/readcoeos06.f90:119 -- stale loop index used as the log10_ne column; header claims it is SAVEd, it is not
- class: logical
- severity: low
- confidence: high
- provenance: inherited (readcoeos06.f `alogNe(jcs,j)` with j left over
  from the `do j=1,mv` init loop); the "var_idx is SAVE'd" comment is a
  modernization leftover
- detail: Every row's log10_ne_grid_06(density_row, var_idx) write goes
  to column var_idx = mv+1 = 11 (the init loops run in the same call),
  so the electron-density table is never populated and column 11 is
  overwritten nr*nt times. Harmless today because nothing reads
  log10_ne_grid_06, but if readcoeos06 were ever entered with
  readcoeos06_init_flag already set (e.g. a second table load), var_idx
  is an uninitialized local and the write lands at an arbitrary column.
  readcoeos01.f90 has the same shape.

## Weak/uncertain observations

- src/eos/eos_lib.f90:430-437 (call site, outside assignment): eos_get_gamma1 takes esac06's out-of-table alternate return to label 100 and returns the PREVIOUS point's eos_output_06(8),(9) with ierr=0; documented as "preserved verbatim from calcad", so pulse-output gamma1/grad_ad for shells outside the 2006 table silently repeat the neighbouring shell. Also no OPAL/SCV ramp here, unlike eqstat.
- src/eos/eqstat.f90:808: any rhoofp*/esac* failure inside the OPAL region falls back to the Yale/SCV values for that shell with no message and no ramp, producing EOS discontinuities between adjacent shells; with the 11-iteration cap and 0.5e-7 tolerance in rhoofp/rhoofp01 this is not a rare path (inherited design).
- src/eos/opal/rhoofp*.f90 bracket start (rhoofp06.f90:133-137): after `P1 > pnr` the lower guess is set to 0.2*rho1 without verifying P(0.2 rho1) < pnr, so regula falsi may start unbracketed (inherited from rhoofp.f).
- src/eos/opal/gmass01.f90:92 / gmass06.f90:89: electron-mass term (moles-1)*m_e dropped relative to gmass.f90:92 -> tmass, hence the cv scaling moles*R/tmass, differs by ~5e-4 between 1995 and 2001/2006 (inherited, upstream OPAL behaviour).
- src/eos/opal/esac.f90 (1995) discards the tabulated gamma1/gamma2 and returns the Landau-Lifshitz "DIRECT" values from cv, which the 2001 authors explicitly flag as not accurate enough (radsub01.f90:99-101); OPAL95 grad_ad/cp are therefore lower quality than 2001/2006 (inherited).
- src/eos/opal/esac06.f90:89-91 and esac01.f90:93-95: deriv_order > 9 only prints a warning and continues; column 10 is 1e35 in the 2001/2006 tables, so a caller passing 10 gets garbage rather than ierr (inherited; no current caller does).
- src/eos/opal/esac01.f90:231: `'(A, " ihi,ilo,imd",3I5)'` with only routine_id in the list; esac06 passes hi/lo/mid (cosmetic, unreachable since k3 >= 3).
- src/eos/opal/readcoeos01.f90:57 and readcoeos06.f90:74: unconditional `close(2)`; unit 2 is not used anywhere else in src, so harmless (inherited).
- src/eos/opal/radsub01.f90:39, esac01.f90/esac06.f90:55, rhoofp01.f90:39: single-precision literals for R (83.14510/83.14511) and sigma/c (1.8914785e-3), documented; the 1995 files use 83.1446304d0 -- the two R values differ by 8e-6 relative, a small cross-generation inconsistency in cv.
- src/eos/opal_eos_lib.f90: table_loaded_flag is shared by the three generations while eos_index_inverse*/x_grid_copy* are per-generation; safe only because read_controls forces a single OPAL year per run.
- src/eos/test/test_eos.f90:321 calls esac06 with rad_flag=2 deliberately to exercise the ierr path; fine, but note that path leaves eos_output_06 untouched.

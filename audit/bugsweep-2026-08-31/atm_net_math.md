# Bugsweep findings: atm/ (root + tables/), turnover/acoustic diagnostics, net/net_lib.f90, math/math_lib.f90

Files read in full: atm/atm_lib.f90, atm/atm_table_lib.f90, atm/atmstruct_lib.f90,
atm/envstruct_lib.f90, atm/ttau_lib.f90, atm/tables/{surfp,kcsurfp,alsurfp,altabinit,alfilein}.f90,
net/net_lib.f90, math/math_lib.f90. There is no atm/turnover/ directory in the tree; the
turnover/acoustic-depth diagnostics named in the assignment live in core/turnover_timescale.f90
(gettau/tauintnew successor) and core/observables_lib.f90 (acoustic radius / delta_nu / delta_Pg)
— both reviewed in full as part of this scope.

T-tau relations verified against the standard forms: Eddington offset 0.25*log10(3/4) =
-0.031235 (correct), Krishna-Swamy q(tau) = 1.39 - 0.815 e^{-2.54 tau} - 0.025 e^{-30 tau}
(correct), KS photosphere tau = 0.312 satisfies tau + q(tau) = 4/3 (x_limit values check out).

---

## net/net_lib.f90:343 -- sneut tfac2: `cap-cap` should be `cap*cap` (Itoh 1996 flavor coefficient)
- class: physical
- severity: high
- confidence: high
- provenance: inherited — `git show 6cd5673:src/sneut.f` line 125 has the identical
  `(cvp*cvp - cap-cap)`; the typo traces to a widely-circulated version of F. Timmes's
  sneut fit code.
- detail: The parameter block builds tfac2 = CV^2 - CA^2 + (n_nu-1)(CV'^2 - CA' - CA')
  instead of (CV'^2 - CA'^2) as required by Itoh et al. (1996) for the mu/tau-flavor
  contribution to the (CV^2 - CA^2)-type terms. With ca' = 0.5, `cap-cap` subtracts 1.0
  where 0.25 belongs: tfac2 = -1.318 instead of +0.182 — wrong sign and ~7x wrong
  magnitude. tfac2 propagates into tfac3 and tfac5, which enter the pair (line 574),
  photoneutrino (line 1052), and bremsstrahlung (lines 1233, 1389) loss rates, so all
  three processes carry a distorted qpair/qphot/gbrem correction whenever T > 1e7 K.
  MESA fixed the same typo in its copy of sneut5 years ago.

## net/net_lib.f90:2202 -- N14(alpha,gamma)F18 polynomial: 0.177 should be 0.117 (CF88)
- class: physical
- severity: medium
- confidence: high
- provenance: inherited — `git show 6cd5673:src/rates.f` line 594 has the same
  `+0.177D0*T9P13`.
- detail: The CF88 fit for N14(a,g)F18 is 7.78e9 T9^{-2/3} exp(-36.031 T9^{-1/3} -
  (T9/0.881)^2) * (1 + 0.012 T9^{1/3} + 1.45 T9^{2/3} + 0.117 T9 + 1.97 T9^{4/3} +
  0.406 T9^{5/3}) + ... (confirmed against the CF88 analytic-rate compilation at the
  China Nuclear Data Center mirror of Caughlan & Fowler 1988). In the code's
  T9^{-2/3}-factored form, the T9 coefficient slot (`0.177*t9p13` in r1) should be
  0.117 — a digit transposition. Overestimates rate(11) (He4+N14) by a few percent
  around T9 ~ 0.2-1, i.e. during He burning. All neighboring coefficients in this fit
  and in the C13(a,n) and C12(a,g) fits match CF88.

## net/net_lib.f90:2242-2243 (with 1917) -- frac outputs stored from uninitialized locals when logT <= cutoff
- class: logical
- severity: medium
- confidence: high
- provenance: inherited, aggravated by modernization — original rates.f line 308-312
  does `IF(TL.LE.TCUT(1)) ... GO TO 200` and label 200 (line 620) stores
  HF1(IU)=F3, HF2(IU)=F1 with F1/F3 never assigned on that path; the original's
  blanket `SAVE` made them stale-previous-call values, the modern automatic locals
  make them genuine stack garbage.
- detail: When log_temperature <= star%ctrl%nuclear_logT_cutoffs(1) (true for every
  cool envelope shell), `rates` zeroes rate(1..13) and jumps straight to the output
  block, where frac_c12_alpha(zone_idx) and frac_be7_electron(zone_idx) are assigned
  from c12_alpha_frac / be7_electron_frac, which no path has set. mixing/mix.f90:182
  calls rates for every zone of a mixed region including cool surface-CZ shells, and
  mixing/solve_composition.f90:210 mass-averages frac_c12_alpha over those zones, so
  garbage can enter the N15 branching fraction used in abundance updates (in the F77
  code it was a silently stale value from the previous hot shell).

## core/turnover_timescale.f90:258-263 + 315 -- innermost-point TAUCZ result overwritten by the averaging block
- class: logical
- severity: medium
- confidence: high
- provenance: modernization — original tauintnew.f puts the label-110 averaging code
  inside the ELSE of the `IF(VALTEST.GT.0.0)` test, so the "INNER MOST POINT PSCA IS
  BELOW BCZ" branch computes TAUCZ = PSCA/CVEL and control passes to the ENDIF,
  never reaching the averaging method.
- detail: In the goto-elimination the label-140 jump became the `spline_taucz_done`
  flag, but only the spline-exit path (line 308) sets it. The first-point branch
  (val_test > 0 at cz_base_index, lines 258-263) computes
  star%convective_turnover_timescale but leaves spline_taucz_done = .false., so the
  `if (.not. spline_taucz_done)` block at line 315 runs the multiple-cell averaging
  method and overwrites the timescale with cz_width/convective_velocity_bcz — a
  different quantity by construction. Fires whenever the first grid point of the
  surface CZ already lies more than one pressure scale height above the interpolated
  BCZ (coarse grids, thin CZs); silent numerical divergence from the F77 behavior in
  the wind-saturation/deuterium-limiter and history taucz.

## core/turnover_timescale.f90:192-198 -- core_cz_top_index read while uninitialized (guard tests garbage)
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited (tauintnew.f line 83 reads HP(JCORE) with JCORE unset on the
  same path; the F77 SAVE turned it into a stale value), aggravated by modernization
  (automatic local in a contained procedure = indeterminate every call).
- detail: The search loop assigns core_cz_top_index only when the scan crosses a
  convective-to-radiative transition. For the very common configuration "radiative
  center up to the base of the surface CZ" (no convective core; e.g. the Sun), the
  scan starts in_radiative_zone = .true. and exits at the first convective shell
  without ever setting core_cz_top_index. The KC 2025-05-30 guard
  `if (cz_base_index .gt. 0 .and. core_cz_top_index .gt. 0)` then evaluates an
  uninitialized integer: with nonzero garbage in [1,json] it can spuriously trigger
  (or wrongly suppress) the skip-shells-above-core `cycle search`, and with garbage
  > json line 193 reads log10_pressure out of the filled range. The guard papers over
  the runtime error seen with bounds checking instead of initializing the variable
  (e.g. to 0 before the search, plus skipping the trap when no core transition was
  found).

## atm/atm_lib.f90:220-227 -- atm_get_surface_pt swallows alsurfp's error (missing `ierr = jerr`)
- class: logical
- severity: medium
- confidence: high
- provenance: modernization (new-code facade, 2026; ierr conversion incomplete)
- detail: The routine sets ierr = 0, calls alsurfp(..., jerr), and then has
  `if (jerr == 0) then / return / end if` followed immediately by the end of the
  subroutine — the funnel body that every sibling facade has (`ierr = jerr; return`,
  cf. eos_lib.f90:191-194, envint_kernel.f90:236-239, and atm_init's own alfilein
  handling at atm_lib.f90:130-133) is missing. When alsurfp hits its fatal
  TEFFL < TEFFLmin path (which sets jerr = 1 after having already set
  lookup_failed = .false.), the caller (wind/massloss.f90:247) sees ierr = 0 and
  lookup_failed = .false. and proceeds with stale atm_table%atm_log10_pressure /
  atm_log10_temperature. The pre-conversion behavior was a hard STOP.

## atm/tables/alfilein.f90:81,132 -- latmtptau100 read while uninitialized; intended check never fires
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited typo (alfilein.f line 104 tests LATMTPTau100, a typo for the
  common-block flag LALTPTau100; IMPLICIT LOGICAL*4(L) made it a never-assigned
  static logical) — made worse by modernization: the modern explicit local
  `logical :: latmtptau100` has no SAVE and no initializer, so line 132 reads an
  indeterminate stack value each call.
- detail: The old-NextGen branch is supposed to abort when the user requests
  tau=100 boundary values (allard_use_tau100) with a 1999-format file that lacks
  them. Because of the inherited typo the test reads a distinct, never-assigned
  local instead of star%ctrl%allard_use_tau100: (a) the compatibility check can
  never fire as intended (tau100 request + old file proceeds against -999-filled
  tau100 tables), and (b) in the modern build the branch is formally UB and, if the
  stack garbage evaluates true, an old-NextGen run fails spuriously at table load.
  The 2026 header comment asserts the local is "always-default-valued", which is
  not a property Fortran gives an unsaved local — the comment's own justification
  for preserving it verbatim is unsound.

## atm/tables/surfp.f90:99-118 (same in kcsurfp.f90:102-121) -- high-gravity branch ignores the -999 edge (gmax index)
- class: physical
- severity: medium
- confidence: medium
- provenance: inherited — the G Somers 5/15 patch changed only the general-case
  search loop bound (NG -> IMINMAX, original surfp.f) and left the
  `IF(GL.GE.ATMGL(NG-1))` top-4-columns branch and the fall-through fallback
  untouched.
- detail: For rows where the high-gravity table entries are -999
  (kurucz_gmax_index(row) < ng — the very case the Somers patch handles), two
  requests go wrong: (1) log g >= logg_table(ng-1) takes the "top 4 log g" branch
  and splines over columns ng-3..ng, which include the -999 sentinels the patch was
  meant to exclude; (2) log g in [logg_table(gmax-1), logg_table(ng-1)) matches no k
  in the restricted search loop (largest tested interval tops out at
  logg_table(gmax-1)), so the loop completes and the `k < gmin` fallback fires,
  silently extrapolating from the FIRST four valid (low-gravity) points — the
  opposite end of the row. Both produce a badly wrong boundary pressure rather than
  a diagnostic. Fires only for tables carrying -999 high-g edges (the reason the
  gmin/gmax machinery exists).

## atm/tables/altabinit.f90:114-116 -- format 910 `(a,i3)` cannot print its two character args + integer
- class: logical
- severity: low
- confidence: high
- provenance: inherited — original altabinit.f line 104/106 has the same
  `write(ISHORT,910) 'ALTABINIT: ...', 'Less than 4 rows: nTeff = ', nTeff` against
  `910 format(A,I3)`.
- detail: The edit list is (character, character, integer) but the format supplies
  only one A before an I3; the second character actual meets the I3 descriptor,
  which is a runtime formatting error (gfortran: "Expected INTEGER ... found
  CHARACTER"). This is on the "fewer than 4 Teff rows" diagnostic path, so the one
  time a malformed Allard table needs reporting, the reporter itself crashes with
  an I/O error instead of the intended message. Formats 900 (`2a,f5.0,f7.2`) and
  920 match their edit lists; 910 is missing an `a` (should be `2a,i3`).

## atm/tables/altabinit.f90:88 -- low-gravity extended range is 4 column widths, comment and alsurfp doc say one
- class: physical (comment/intent mismatch)
- severity: low
- confidence: medium
- provenance: inherited — original altabinit.f line 78: `GLmin(i) = GLs(j1) -
  4D0*(GLs(j1+1) - GLs(j1))`, under the same "one column's width less" comment.
- detail: Step 3's header (lines 21-24, and alsurfp.f90's header lines 20-23)
  documents the permissible extrapolation range as one column width beyond each end
  of the row, and the max side implements exactly that (line 90); the min side
  subtracts FOUR column widths, allowing 4-point Lagrange extrapolation up to 4 grid
  spacings below the lowest tabulated gravity of a row before alsurfp's
  out-of-extended-range trap fires. Either the 4d0 is an unannotated deliberate
  widening (then three separate comments are wrong) or a genuine defect permitting
  far-field extrapolation of the surface pressure at the low-g edge. Reported as an
  intent/code mismatch per the sweep rules.

---

## Weak/uncertain observations (one line each)

- core/turnover_timescale.f90:202-203 — `convective_flag(cz_base_index+1/+2)` can index past num_points (reads unfilled elements) when the located CZ base is within 2 shells of the top; inherited (LCZ(IMAX+1/+2), tauintnew.f:90).
- core/turnover_timescale.f90:338-348 — averaging-branch edge guard only handles k==cz_base_index+1; k==cz_base_index (possible when the CZ-center radius is below the first point) gives kk-2 = cz_base_index-2, underflowing for a base near index 1; inherited.
- core/turnover_timescale.f90:244/326 — interp_fraction = dd2/(dd2-dd1) can divide by ~0 when gradr-grada is equal at adjacent points; inherited.
- core/turnover_timescale.f90:365-367 — the taucz>1e20 retry sets search_start_index = cz_base_index+1, which can exceed num_points; the next scan then reads convective_flag(search_start_index) beyond the filled range; inherited.
- atm/tables/surfp.f90:175-181 — modernization dropped the original's WRITE(ISTOR,...) copies of the log-P diagnostic (model-store stream loses two lines; output-only); same in kcsurfp.
- atm/tables/alfilein.f90:198/287 vs 322-323 — FeH skip tolerance 1d-5 but store-verify tolerance 1d-6: a record with |dFeH| in (1d-6,1d-5] is counted in pass 1 but silently not stored in pass 2, leaving -999 holes that altabinit then rejects; alpha uses `.gt.0d0` (exact equality) in the skip but 1d-6 in the verify; inherited.
- atm/tables/alfilein.f90:236-251 — FeH/alpha uniqueness scans bounded by allard_num_gl instead of their own counters (acknowledged in a code comment; preserved inherited bug; harmless while FeH/alpha grids are single-valued).
- atm/tables/altabinit.f90:63-77,86-91 — a row with no valid entries leaves gl_index_min/max unset (garbage j1/j2 used in step 3 before the step-4 validation runs); gl_grid(j1+1)/(j2-1) can go out of bounds when a row's single valid entry sits at an edge; inherited.
- atm/tables/alsurfp.f90:241-257 — passes trailing subarrays (from j1/teffl_index) to polint's fixed-size-20 dummies; fewer than 20 elements remain near the table top (sequence-association nonconformance, harmless since polint reads only n=4); inherited convention, documented at line 104.
- math/math_lib.f90:94,142 — `iy = floor(y)` executes before the |iy|<100 magnitude check, so a huge real exponent overflows the default integer (UB) before the guard; MESA's version tests `y == aint(y)` on the real first; new-code, theoretical for the exponents in use.
- net/net_lib.f90:2184 — C13(a,n)O16 linear coefficient 0.0129 vs CF88's published 0.013 (rounding-level; inherited, possibly from VandenBerg's notes).
- net/net_lib.f90:2035 — `hydrogen_fraction.eq.0.0` float-equality gate (KNOWN audit item class, labeled per instructions).
- net/net_lib.f90:1476-1479, 1950, 2094 etc. — single-precision literals (-2.25/-4.55 pow exponents, 1.017677E-4, 1.752E-10) (KNOWN finding class: single-precision e-notation literals in net; not re-reported).
- atm/atm_lib.f90:182 — Castelli failsafe uses `ng` not `ngc` (acknowledged in the file header as preserved and harmless, both 11; not re-reported).
- core/observables_lib.f90:158-161 — core_boundary_fx2 formula is numerator==denominator (identically 1) and the stale envelope_boundary_fx is used instead; inherited (wrtout.f:141-143), dead code, already documented in the module header as preserved.
- core/observables_lib.f90:124 — luminosity renormalization divides by the breakdown sum with no zero guard; inherited behavior.
- Verified clean: ttau_lib constants/forms; sneut's other Itoh-1996 fit coefficients and derivative algebra (spot-checked against the published sneut5 forms); rates' Graboske screening constants (lambda0 5.9426e-6, z86/z53 values, Fermi momentum 1.017677e-4, m_ec^2/k = 5.930), c21 = 5.2404e-8, CF88 C12(a,g)O16 structure; azbar/ifermi12/zfermim12 against Antia (1993) coefficient tables; the surfp/kcsurfp goto-to-do conversions are index-semantics-faithful (including zero-trip and fall-through loop-index values); compute_turnover_timescale's span selection matches gettau.f's CHKPRS logic; delta_Pg = 2 pi^2/sqrt(l(l+1)) and nu_max/delta-nu scalings in observables_lib are correct.

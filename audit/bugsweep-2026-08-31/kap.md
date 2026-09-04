# kap/ bug sweep — opacity table readers/interpolators (opal95, opal92, alex06, alex94, kurucz90, laol89, conductive)

All 39 assigned .f90 files read in full. Provenance checked against git revision
6cd5673 (`src/*.f`) for every finding.

## /Applications/YREC/src/kap/conductive/condopacpint.f90:119 -- conductive dlnkap_dlnT built from the dlnkap_dlnRHO species derivatives (copy-paste)
- class: logical
- severity: high
- confidence: high
- provenance: inherited — 6cd5673:src/condopacpint.f has `QOTC = CONDL * (WtH1*QODCX + WtHe4*QODCY + WtOx*QODCZ)`, i.e. the same QOD* (rho) inputs
- detail: Lines 114-120 compute both `conductive_dlnkap_dlnrho` AND
  `conductive_dlnkap_dlnt` from `dlnkap_dlnrho_h1/he4/ox`; the temperature
  derivatives `dlnkap_dlnt_h1/he4/ox` returned by the three condopacp calls
  (lines 95-102) are never used anywhere. The commented-out predecessor formula
  directly below (line 121-122, `QOTC = ... QOTCX ... QOTCY ... QOTCZ`) shows the
  intent was to use the T derivatives. Consequence: whenever Potekhin conductive
  opacity is active, dlnkap_dlnT of the conductive part equals
  3 - (rho-derivative mix) instead of 3 - (T-derivative mix), poisoning the
  blended dlnkap_dlnT handed to the Henyey solver (wrong Jacobian, not wrong
  equilibrium kappa). Inherited from the F77 original verbatim.

## /Applications/YREC/src/kap/conductive/condopacpint.f90:114 -- derivative mixing formula multiplies by the mixture log-conductivity instead of conductivity-weighted averaging
- class: physical
- severity: medium
- confidence: high
- provenance: inherited — same `QODC = CONDL * (WtH1*QODCX + ...)` in 6cd5673:src/condopacpint.f
- detail: With mix = -log10(S), S = sum_i w_i 10^(-K_i), the correct derivative is
  d(mix)/dx = [sum_i w_i 10^(-K_i) dK_i/dx] / S. The code computes
  `mix_log10_cond * sum_i(w_i dK_i/dx)` — it multiplies an unweighted derivative
  sum by the *value* of the mixture log conductivity (a dimensionally meaningless
  scale factor that can be any magnitude or sign). The commented-out older
  formula at least carried the 10^(-K_i) weights. Both conductive derivative
  outputs are wrong in magnitude whenever the Potekhin branch fires. (Again a
  Jacobian-quality problem, not a kappa-value problem.)

## /Applications/YREC/src/kap/conductive/condopacpint.f90:146-155 -- extrapolation branch (rho < 1e-6): CONDL sign flipped and extrapolation step applied with reversed sign
- class: physical
- severity: medium
- confidence: high (that the two branches are mutually inconsistent); medium (net impact)
- detail: In the main branch (rho >= 1e-6), `mix_log10_cond = -log10(S)` and
  `OCL = -3.5194 + 3*TL - DL - mix`. In the extrapolation branch the sign is
  dropped: `mix_log10_cond = +log10(S)` (line 146) yet OCL still subtracts it
  (line 151), so at the branch boundary DL = -6 the two formulas differ by
  2*log10(S) — a hard discontinuity proving one branch is wrong. Additionally
  the extrapolation step `OCL = OCL + conductive_dlnkap_dlnrho*(-6.0D0-DL)`
  (line 154-155) uses (-6 - DL) where linear extrapolation from the DL=-6 edge
  requires (DL - (-6)); the correction is applied with reversed sign, so kappa_cond
  *decreases* toward lower density instead of increasing (kappa_c ~ 1/(rho*lambda)).
  Fires for T < 10^4 K, rho < 1e-6 with log R >= 0 (cool-star photosphere
  conditions) when conductive opacity is enabled; in that regime the
  near-neutral gas usually makes kappa_cond enormous either way, so the blend
  1/k = 1/kr + 1/kc mostly hides it — hence medium not high.
- provenance: inherited — identical `CONDL = LOG10(...)` (no minus) and
  `OCL = OCL + QODC*(-6.0D0-DL)` in 6cd5673:src/condopacpint.f.

## /Applications/YREC/src/kap/opal95/getopal95.f90:142/183 -- density_shifted read uninitialized when the low-rho/low-T shortcut is taken
- class: numerical
- severity: medium
- confidence: high
- provenance: modernization — 6cd5673:src/getopal95.f has `SAVE` and `GOTO 60`, so LSHIFT was the (dubious but defined) stale value from the previous call; the modern local `logical :: density_shifted` has no save attribute
- detail: When `hydrogen_fraction < 0.2 .and. opal95_index_t <= 5` (line 136),
  the true branch sets only `low_regime_flag`; `density_shifted` is assigned only
  in the else branch (line 142). Line 183 then reads `.not.density_shifted` to
  choose between the uniform-rho-grid and per-temperature-row interpolation
  paths — an uninitialized-variable read (UB). The F77 original relied on SAVE to
  reuse the previous call's LSHIFT here; the conversion turned "stale state" into
  "undefined". On gfortran a zeroed stack usually gives .false. (the common
  case), but nothing guarantees it, and the paired state flag
  `opal95_extrap_hi(_row)` (kept in the type, so still effectively SAVEd) can now
  disagree with the local, applying stale high-edge extrapolation offsets.
  Fires exactly in the X<0.2, log T <~ 3.9 corner (cool He-rich material).

## /Applications/YREC/src/kap/opal95/ll95tbl.f90:80 (+ opacity_table_lib.f90:176) -- OPAL95 storage-slot collision: slots 110 and 119 double-booked, slots 109/118 never written
- class: logical
- severity: medium
- confidence: high
- provenance: inherited — identical `DATA NZ/0,10,...,100,109,118/` and `N = NZ(IZ)+IX` in 6cd5673:src/ll95tbl.f, and the same `NZ(JZ)+JX` addressing in op954d/op95ztab
- detail: Tables are stored at n = opal95_table_start_index(iz) + ix. With
  start_index = [0,10,...,90,100,109,118] and the sparse high-Z coverage
  (ix=1..8 for all 13 Z; ix=9 only Z<=0.04; ix=10 only Z<=0.08), the writes
  collide: (ix=1,iz=12: X=0, Z=0.08) and (ix=10,iz=11: X=0.94, Z=0.06) both map
  to slot 110; (ix=1,iz=13: X=0, Z=0.10) and (ix=10,iz=12: X=0.92, Z=0.08) both
  map to 119. Since ix=10 is read last, the X=0 tables at Z=0.08 and Z=0.10 are
  clobbered by X=1-Z tables. Any later read of slot start(12)+1 or start(13)+1 —
  the 4D interpolation (opal95_interp4d via getopal95) with jx=1 and a Z node at
  0.08/0.10 (i.e. X < 0.1 with Z >~ 0.04, metal-diffusion runs), or
  opal95_fixed_z_table's i=1 (X=0) plane for model Z > 0.04 — silently uses an
  X~0.93 table where an X=0 table belongs, a large kappa error. Slots 109/118
  are never written, but I verified no read can address them (the "absent
  tables at high X" remaps make jx=9 unreachable at those Z nodes). Dormant at
  solar Z; real for Z>0.04 / low-X regions.

## /Applications/YREC/src/kap/opal95/opal95_fixed_z_table.f90:34-48 -- Z in (0.08, 0.1] wrongly rejected as "outside OP95 table range"
- class: logical
- severity: medium
- confidence: high
- provenance: inherited — 6cd5673:src/op95ztab.f has the same `DO I = 3,NUMZ-1` search and STOP
- detail: The 4-nearest-tables search loops i = 3..num_z-1 = 12 looking for the
  first grid_z(i) >= Z, with grid_z(12) = 0.08 and grid_z(13) = 0.10. For any
  requested fixed-table Z in (0.08, 0.1] no i in range matches and the routine
  errors out ("DESIRED Z ... OUTSIDE OP95 TABLE RANGE" — the comment says
  "DESIRED Z > 0.1D0", which is wrong: it fires from Z > 0.08). The correct
  clamp is z_table_index = num_z-3 (nodes 0.04,0.06,0.08,0.10), exactly what
  getopal95's own 4D Z search does for the same value (getopal95.f90:300-324).
  Kills (or, post-ierr-conversion, error-returns) startup for a legal
  opal95_single_table_z in (0.08,0.1].

## /Applications/YREC/src/kap/alex06/alex06tab.f90:128 -- X interpolation node 4 set to 1-Z whenever iz==num_z-3, even when the 4th data column is a regular grid-X table
- class: physical
- severity: medium
- confidence: high
- provenance: inherited — 6cd5673:src/alex06tab.f has the same `IF(IZ.EQ.NUMZ-3)THEN QR(4) = 1.0D0-ZE`
- detail: The 4th X data column uses the special X=1-Z table only when
  `index_x == num_x-3 .AND. iz == num_z-3` (line 103); otherwise it is the
  ordinary table at grid_x(index_x+3) (line 110). But the matching abscissa
  override at line 128 tests only `iz == num_z-3`. So for Z > 0.06 (which forces
  iz = num_z-3) with moderate X (index_x < num_x-3, i.e. X <= ~0.5), the
  Lagrangian weights are computed for nodes (e.g.) {0.2, 0.35, 0.5, 1-Z~0.92}
  while the data columns are at {0.2, 0.35, 0.5, 0.7} — all four weights wrong,
  wrong low-T opacity. Condition should mirror line 103's compound test.

## /Applications/YREC/src/kap/opal95/getopal95.f90:169 -- high-rho shift search scans the full (Z,X) table array with an X-only index
- class: logical
- severity: low
- confidence: medium
- provenance: inherited — 6cd5673:src/getopal95.f: `IF(CAPPA(IXX,ITT,IDD+3)...)` then `IF(CAPPA2(IXX,ITT,J).LE.9.9D0)`
- detail: The per-temperature density-shift search checks the fixed-Z table
  (`opal95_fixed_z_opacity`, first index = X slot 1..10) at line 166 but scans
  `opal95_full_opacity` (first index = global table number 1..126) at line 169
  with the same X-slot index — addressing the Z=0 table block regardless of the
  actual Z. The 9.999 "no data" footprint is largely format-driven (rows 58-70)
  and thus identical across tables, which is why this mostly works; but any
  in-file 9.999 cells that differ between Z=0 and the working Z make the shifted
  rho base index wrong at the high-rho/high-T table edge. Should be
  `opal95_fixed_z_opacity(x_table_index,...)` (or the proper start_index+ix slot).

## /Applications/YREC/src/kap/alex94/alex94_interp3d.f90:53-64 -- first call after table load can read grid_logt(num_t+1) out of bounds (index clamp too weak for a 4-point stencil)
- class: logical
- severity: low
- confidence: medium
- provenance: inherited — 6cd5673:src/yalo3d.f lines 30-32 have the identical `IF((IT+2).gt.NUMT) IT=NUMT-2` insurance, and 6cd5673:src/alxtbl.f READs the file's leading row integer directly into the common /ALOT/ IT
- detail: read_alex94_tables.f90:84 reads each row's leading integer into the
  cached search index `alex95_index_t`; after loading, it holds the last row's
  label (up to 23 = num_t). The "insurance" clamp only enforces
  index_t <= num_t-2, but the 4-point stencil needs index_t <= num_t-3, and the
  downward search (`do i = index_t+1,2,-1`, sets index_t = i-1) can leave
  index_t = num_t-2 = 21 when log T is in (4.05, 4.10]. interp_nodes(4) then
  reads alex95_grid_logt(24) — past the 23-element array (in practice the
  adjacent type member, alex95_grid_x(1) = 0.0, giving garbage Lagrangian
  weights). Fires at most on the first alex95 evaluation after startup with
  log T in (4.05, 4.10]; subsequent calls have a sane cached index. Same clamp
  pattern exists in getalex06.f90:55-56 but there the index state is never
  clobbered above num_t-3, so it stays dead.

## /Applications/YREC/src/kap/conductive/condopacp.f90:58/93-95 -- hunt() warm-start indices lost their SAVE: read uninitialized on every call after the first
- class: numerical
- severity: low
- confidence: high (defect); low (numerical impact)
- provenance: modernization — 6cd5673:src/condopacp.f has a blanket `save`, so IZ/IT/IR persisted between calls; the modern locals `z_index, t_index, r_index` (line 58) have no save
- detail: The midpoint initialization at lines 93-95 happens only inside the
  one-time table-load branch; on every later call the three locals are
  uninitialized when passed as hunt()'s intent(inout) initial guess
  (lines 100/108/121). hunt resets an out-of-range guess and bisection finds the
  correct bracket for any in-range garbage, so results are still correct — but
  it is a formal uninitialized read (UB, and defeats the warm-start
  optimization the original SAVE provided). The table itself was migrated to
  conductive_table_lib during the save-migration campaign; these indices were
  missed.

## /Applications/YREC/src/kap/kap_lib.f90:354-389 -- Hubbard-Lampe conductive opacity applied even when use_conductive_opacity = .false.
- class: logical
- severity: low
- confidence: high (behavior); low (that it is unintended)
- provenance: inherited — 6cd5673:src/getopac.f lines 168-200: `IF (LcondOpacP) ... ELSE LCONDO=.FALSE.; IF (.NOT.LCONDO) ... Hubbard Lampe ...`
- detail: The use_conductive_opacity flag only selects Potekhin vs fallback;
  with the flag off, got_conductive_opacity=.false. sends control straight into
  the Hubbard-Lampe branch (line 369), which applies a conductive correction
  whenever log T >= 4.2 and log rho >= 2 log T - 13. There is no configuration
  that disables conductive opacity entirely. Faithful to the original (HL was
  the historical default), but the flag's name promises more than it does;
  worth knowing when comparing against codes with conduction truly off.

## /Applications/YREC/src/kap/alex06/readalex06.f90:67-72 -- X/Z header-mismatch check prints "RUN STOPPED" but neither stops nor sets ierr
- class: logical
- severity: low
- confidence: high
- provenance: inherited — 6cd5673:src/readalex06.f line 71 has the STOP commented out (`C STOP`) in the first block, while the second block's identical check (line 110) does STOP
- detail: In the main X-table loop the header consistency check writes the
  "ERROR ... RUN STOPPED" message and then continues reading (the STOP was
  commented out historically; the modern ierr-conversion faithfully preserved
  the no-op, unlike the R- and T-row checks in the same loop which do
  `ierr=1; return`). A mis-ordered or wrong-composition alex06 table file loads
  silently except for a log line claiming the run stopped. The parallel check
  for the X=1-Z block (line 111-117) does error out — an asymmetry between two
  blocks that should mirror each other.

## /Applications/YREC/src/kap/kap_lib.f90:180-235 -- pure-Z branch: table_metal_fraction used uninitialized if no OPAL95/OPAL92/LAOL interior table is selected
- class: logical
- severity: low
- confidence: medium
- provenance: inherited — 6cd5673:src/getopac.f pure-Z block has the same IF/ELSE IF chain with no final ELSE; ZTAB there was SAVE-stale rather than trapped
- detail: In the He-burning (Z>0.1, log T>7) branch, `table_metal_fraction` is
  set only inside the three interior-table arms. A (mis)configuration with only
  Kurucz/Alexander atmosphere tables plus a pure-Z table reaches line 225's
  `slope = (log10_opacity - purez...)/(table_metal_fraction - 1.0d0)` with both
  log10_opacity and table_metal_fraction unset. The "NO OPACITY TABLE CHOSEN"
  trap (line 322) guards only the non-pure-Z path. Requires a nonsensical but
  accepted configuration; the original at least reused a SAVEd stale ZTAB
  instead of a truly undefined local.

## /Applications/YREC/src/kap/laol89/gtlaol.f90:47 -- extrapolation tolerance takes the natural log of tollaol on a log10 axis
- class: physical
- severity: low
- confidence: medium
- provenance: inherited — 6cd5673:src/gtlaol.f line 28: `TOLL = LOG(TOLLAOL)` (same in gtlaol2.f, gtpurz.f)
- detail: `log_extrap_tolerance = log(tollaol)` with tollaol = 10.0 gives 2.303,
  which is then compared against log10(rho) offsets (lines 93/103): the code
  permits ~2.3 dex of linear extrapolation beyond the table edge, where a
  "factor of TOLLAOL" reading of the input (log10 -> 1.0 dex) was presumably
  intended. Comment-vs-code intent mismatch; consistently inherited across all
  three LAOL lookup routines.

## /Applications/YREC/src/kap/kurucz90/read_kurucz_tables.f90:67-69 -- no bound on density_index against the 50-column array
- class: logical
- severity: low
- confidence: medium
- provenance: inherited — 6cd5673:src/setkrz.f has the same unguarded ID=ID+1 per record (only the temperature count is bounds-checked, and that check is itself part of the modernization; the original had neither)
- detail: Within one temperature block, every data record increments
  density_index and stores into kurucz_log10_rho/opacity(num_read, density_index)
  with second dimension kurucz_max_num_densities = 50. A Kurucz-format file with
  more than 50 density rows per temperature silently overruns into the next
  row's storage (derived-type members, so no segfault, just corrupted tables).
  The temperature direction got an explicit guard during modernization; the
  density direction did not.

## Weak/uncertain observations

- getopal95.f90:136-139: `low_regime_flag` (original LCHK) is set and never read — the "postponed" low-rho/low-T empty-region check was never implemented (inherited dead code; the guard the comment promises does not exist).
- getopal95.f90:164-178: shift search can set opal95_index_rho(i) = j-3 <= 0 for j<=3, and leaves index_rho(i) stale if the inner scan finds nothing (inherited; original GOTO 50 with unset JD(I)).
- ll95tbl.f90:168: the between-tables header read (format 900) has no end= guard — a truncated OPAL95 file aborts with a runtime I/O error instead of a diagnostic (inherited).
- ll95tbl.f90:86, readalex06.f90:67/91, read_alex94_tables.f90:71/90, read_kurucz_tables.f90:49: float-equality comparisons between file-read values and DATA constants — safe in practice (identical decimal->binary conversion path), inherited idiom.
- read_opal92_tables.f90:58-65: num_temps_read is undefined if the very first header read hits EOF (jump to 97), and opal92_num_temps is taken from whichever table was read last, assuming all three X tables share a T grid (inherited).
- opal92_surface_table.f90:65/115: `rho_search_index` passed to opal92_interp2d(_z2) as the findex hint while uninitialized — benign (findex resets out-of-range hints and searches both directions); original relied on SAVE.
- build_kurucz_splines.f90:43 with read_kurucz_tables.f90:68: the reader stores exp10(logk), so the `chko <= 0` missing-data filter can only catch blank fields (read as 0.0 -> 1.0 after exp10, which passes!) via underflow; the filter as written can essentially never reject a value (inherited; actual Kurucz files appear to be dense so it never mattered).
- kurucz.f90:119-135: when local_logrho lands exactly on a row's last density knot, cubic coefficients c3/c4 of the final knot are used; ysplin never writes c(3,n)/c(4,n), so stale spline_work is multiplied by delta=0 — fine unless the stale garbage is NaN/Inf (inherited).
- kap_lib.f90:337-344: the atmosphere/interior ramp blends kappa linearly in linear space but blends dln-derivatives linearly too, so the returned derivatives are not the derivatives of the returned kappa (inherited approximation).
- kap_lib.f90:150-162/277-318: two-Z interpolations divide by (z1 - z2), and the pure-Z branch divides by (table_z - 1); equal configured Z values give division by zero (config-dependent, inherited).
- getalex06.f90:112-118: only high-side log R gets linear extrapolation; below grid_logr(1) the 4-point Lagrangian extrapolates cubically unguarded (inherited; OPAL95 handles the low side explicitly).
- opal95_fixed_z_table.f90:71-97: for metal_fraction exactly 0.04 the X=0.95 plane (index 9) is left stale; verified unused by the consistent index remaps in getopal95/opal95_surface_table (inherited, self-consistent).
- condopacpint.f90:77-79: the species weights are not renormalized by their sum, so partial ionization *raises* the inferred mixture conductivity (fewer ion scatterers model); crude but arguably intended — flagged only as model quirk (inherited).
- ll95tbl.f90:32-50: the out-of-sync 16- vs 17-member /newopac/ layout is pre-existing, harmless at the referenced offset, and already documented in-file; not re-reported as new.

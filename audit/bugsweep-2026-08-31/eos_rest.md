# EOS review: eos_lib.f90, eqstat.f90, mu.f90, eos/mhd/*, eos/scv/*, eos/yale/*

All assigned files were read in full and diffed against the F77 originals at
6cd5673 (eqstat.f, eqstat2.f, eqsaha.f, eqrelv.f, eqscve.f, eqscvg.f, setscv.f,
meqos.f, mhdpx.f, mhdpx1.f, mhdpx2.f, mhdst.f, mhdst1.f, rabu.f, rtab.f, mu.f,
setups.f, eqbound.f, physic.f). The modernizations are remarkably faithful;
most defects found are inherited from the F77 originals.

## /Applications/YREC/src/eos/eqstat.f90:587 -- ramp blend adds raw Saha QDTT/QDTP instead of the Saha-minus-fully-ionized difference
- class: physical
- severity: medium
- confidence: high
- provenance: inherited (6cd5673:src/eqstat2.f — lines 219-220 save `QDTT0 = QDTT`, `QDTP0 = QDTP`; the difference block at lines 240-248 subtracts DL0, QDT0, QDP0, QCP0, DELA0, QCPT0, QCPP0, QAT0, QAP0 but never `QDTT0 = QDTT0 - QDTT` / `QDTP0 = QDTP0 - QDTP`; lines 261-262 then apply `QDTT = QDTT + WT0*QDTT0`)
- detail: In the ionization-cutoff interpolation zone (saha_log10t_cutoff - 0.5 < logT < saha_log10t_cutoff, with want_derivatives), every blended quantity follows the pattern `save0 = saha_value; save0 = save0 - relv_value; result = relv_value + W*save0` — except the two mixed second derivatives dlnrho_dlnt_dt and dlnrho_dlnp_dt (modern lines 535-536 save, 587-588 apply, no subtraction at 565-572). The blend therefore yields `QDTT_relv + W*QDTT_saha` instead of `(1-W)*QDTT_relv + W*QDTT_saha`, over-counting by `W*QDTT_relv`. Fires on every derivative-mode EOS call in the ramp band; these feed the Henyey corrections (affecting convergence rate/step acceptance, not a converged model directly). The conspicuous asymmetry against the nine correctly-differenced siblings makes intent unambiguous.

## /Applications/YREC/src/eos/mhd/mhdpx1.f90:56 -- out-of-range temperature returns ierr=0 with stale table output
- class: logical
- severity: high (when it fires; MHD path has zero test coverage)
- confidence: high
- provenance: inherited (6cd5673:src/mhdpx1.f lines 34-37: `IF(TL.LT.TMINI .OR. TL.GT.TMAXI) ... GO TO 999`; line 108 `999 RETURN` — a plain return, no error signal)
- detail: When log10(T) falls outside [table_log10t_min, table_log10t_max], mhdpx1 returns immediately with ierr=0 and without touching mhd_eos%mhd_output. mhdpx then reads mhd_output(1) as the density and meqos builds the complete thermodynamic state from whatever the PREVIOUS call left in mhd_output (or zeros on the first call, giving rho = 10^0 = 1). The header comment in mhdpx.f90 (lines 13-14, inherited verbatim) promises "IERR = 1 SIGNALS PGL,TL OUTSIDE THE DOMAIN OF TABLES" — nothing ever sets it for this case, a direct comment/code contradiction. Silent wrong physics for any use_mhd_eos run that strays below/above the table in T (cool atmospheres are the likely place). The only diagnostic ever written for this ("OUT OF TABLE RANGE") sits after an unconditional `return` in mhdpx.f90 and is unreachable — dead code that was clearly meant to run.

## /Applications/YREC/src/eos/scv/scv_envelope_table.f90:290 -- clamped du/dt read is dead: immediately overwritten by unclamped read (and the whole cols 7-12 loop is dead anyway)
- class: logical
- severity: low (output columns never consumed)
- confidence: high
- provenance: inherited (6cd5673:src/setscv.f lines 215-216: `QR(K) = TABLENV(ii,jj,6)` followed immediately by `QR(K) = TABLENV(IDT+K-1,J,6)`)
- detail: In the "DERIVATIVES OF DU/DT" block, `interp_x(k_idx) = tablenv(ii,jj,6)` (jj clamped to min(nptsx(ii), p_idx), matching the adjacent d(ln cp)/d ln T block) is immediately overwritten by `interp_x(k_idx) = tablenv(idtt+k_idx-1, p_idx, 6)` with no clamp — the clamped read was clearly MEANT to be the fix (its twin in the cp block at lines 279-283 survived) but the older unclamped line was never deleted. Where a neighboring T-row is shorter than the current one (nptsx is ragged), the unclamped read pulls never-written (zero) entries. Mitigating fact discovered during review: tablenv columns 7-12 are written here and read NOWHERE in the tree (grep over src/*.f90: only scv_eos_lib declares tablenv, only this file writes cols 7-12; eqscve interpolates only cols 2-6), so the entire third loop of build_scv_envelope_table is dead code. That itself is worth knowing — if anyone ever wires those columns up, this bug is waiting.

## /Applications/YREC/src/eos/scv/scv_envelope_table.f90:202 -- dsmix_dlnt indexes neighbor T-rows at the current row's pressure index without clamping (feeds live cp column)
- class: logical
- severity: medium
- confidence: high (that out-of-row reads occur; medium that the numerical impact is material)
- provenance: inherited (6cd5673:src/setscv.f lines 136-137: `QSMIXT = FTD(1)*SMIX(IDT,J)+FTD(2)*SMIX(IDT+1,J)+FTD(3)*SMIX(IDT+2,J)`, no clamp; the clamped `jj = min(NPTSX(ii),j)` idiom appears only in the later derivative loop)
- detail: Unlike the dead cols-7-12 loop, this one matters: `dsmix_dlnt` (line 202-204) reads smix(idtt..idtt+2, p_idx) where p_idx runs to nptsx(t_idx) of the CURRENT row; a shorter neighbor row yields smix entries never assigned (0.0 in static storage). dsmix_dlnt enters dlns_dlnt -> cp_gas -> tablenv column 5, which eqscve interpolates and returns as specific_heat_cp for every SCV envelope evaluation near the ragged table edge. Same unclamped-neighbor-row pattern recurs in the read side: eqscve.f90 lines 219/264/322/380 (`ytab_work(k) = tablenv(ii, idp+k-1, j+1)`, nodes taken from row idtt but values from rows idtt..idtt+3) and the tablex/tabley reads at 528-535, 555-562; eqscvg.f90 lines 158-169 etc. — all inherited verbatim (eqscve.f lines 150/182/224/266). Fires only where nptsx varies across the 4-row stencil; produces silently wrong rho/cp/grada there rather than an error.

## /Applications/YREC/src/eos/yale/fully_ionized_eos.f90:209 -- non-convergence path returns with all intent(out) results unassigned
- class: numerical
- severity: medium (rare, but silently poisons the zone when it fires)
- confidence: high
- provenance: inherited (6cd5673:src/eqrelv.f lines 168-173: after NDEN>20, WRITE + `RETURN` before `DL = DL8` and all derivative computation), aggravated by modernization (intent(out) dummies are formally undefined on entry; the F77 version at least left the caller's previous-zone values in place by pass-by-reference)
- detail: If the Newton iteration on density fails to converge in 20 steps, the routine prints "EQRELV: Did not Converge" and returns without assigning log10_density, density, or any of the 11 derivative outputs. eqstat2 then blends/returns whatever those variables happened to hold (previous zone's values in practice). The run continues silently on stale thermodynamics. Under modern Fortran the outputs are formally undefined after the early return (intent(out)); in practice the storage still holds caller values, so behavior matches F77 — but it is now standard-UB and would break under a compiler that clobbers intent(out) args. An ierr here (the file already has the luout machinery) would be the stage-3-consistent fix.

## /Applications/YREC/src/eos/eqstat.f90:569 -- differencing of uninitialized saha_*_dt locals when want_derivatives is false (SAVE semantics lost in modernization)
- class: numerical
- severity: low
- confidence: high (uninitialized read is certain; consequences benign in practice)
- provenance: modernization (6cd5673:src/eqstat2.f has a blanket `SAVE`, so `QCPT0 = QCPT0-QCPT` at line 245 read stale-but-defined values; the modern eqstat2 locals saha_specific_heat_cp_dt/dp and saha_adiabatic_gradient_dt/dp — lines 366-368 — carry no SAVE and no initialization)
- detail: The subtractions at lines 569-572 run unconditionally in the ramp branch, but the saha_* derivative saves are only assigned under want_derivatives (lines 534-541). On a want_derivatives=.false. call the subtraction operates on never-set doubles. The results are dead (the apply block at 586-597 is guarded), so no numeric output is affected — but arithmetic on uninitialized doubles can trap under -ffpe-trap=invalid / -finit-real=snan debugging, and it is formally undefined. Moving lines 569-572 inside the existing want_derivatives guard is a zero-behavior-change fix.

## /Applications/YREC/src/eos/mhd/mhdpx1.f90:101 -- ZAMS-path X-interpolation runs quint over 25 variables when only 20 were filled
- class: numerical
- severity: low
- confidence: high
- provenance: inherited pattern, made formally undefined by modernization (original mhdpx1.f VAROUT was SAVE'd, so slots 21-25 held stale centre-table values; the modern `table_vars(ndimt,ivarx)` is an automatic local, uninitialized every call)
- detail: mhdpx2's six ZAMS branches copy only ivarc=20 variables into table_vars (lines 54-109) while the quint loop in mhdpx1 (line 101, `do iv=1,ivarx` with ivarx=25) interpolates all 25, reading table_vars(1:3, 21:25) that were never written this call. mhd_output(21:25) receive garbage; meqos reads only indices 1-20, so nothing downstream consumes it — but quint arithmetic on uninitialized doubles can FP-trap under debug builds, and the MHD path has zero test coverage to catch it. Restricting the ZAMS-branch loop to ivarc (or zero-filling table_vars) is behavior-preserving.

## /Applications/YREC/src/eos/scv/eqscve.f90:590 -- "fraction of hydrogen nuclei ionized" counts an H2 molecule as one nucleus
- class: physical
- severity: low
- confidence: medium
- provenance: inherited (6cd5673:src/eqscve.f lines 444: `XHP = XTF_HP / (XTF_H2 + XTF_H1 + XTF_HP)`)
- detail: xtf_h2/xtf_h1/xtf_hp are per-particle fractions (an H2 molecule is one particle carrying two H nuclei). The comment states xhp is "the fraction XHP of hydrogen nuclei that are singly ionized", which requires n(H+)/(2*n(H2)+n(H1)+n(H+)); the code weights H2 by 1, overestimating xhp wherever molecular hydrogen is significant (coolest SCV region). Same formula duplicated in eqscvg.f90:435. The helium expressions in the same block check out exactly (particle+charge conservation verified), so this is the one member of the family that deviates from its stated meaning. ion_fraction feeds opacity lookup (getopac in the original flow), so the impact is second-order. Possible it was an accepted approximation; the comment/code mismatch is why it is reported.

## /Applications/YREC/src/eos/eqstat.f90:648 -- beta14 changes meaning when an OPAL table takes over (1-beta vs 4*(1-beta)/beta)
- class: physical
- severity: low
- confidence: medium (real inconsistency; unclear any live consumer cares)
- provenance: inherited (6cd5673:src/eqstat2.f: main path `BETA14 = 4.0D0*BETAI*DMIN1(1.0D0,BETA14)`; every OPAL branch sets `BETA14=1.0D0-BETA`)
- detail: The Yale/Prather convention returned to callers is beta14 = 4(1-beta)/beta (capped), and saha_eos consumes it in that convention internally. But whenever an OPAL 95/01/06 branch replaces or ramps the result (lines 648, 699, 773, 789 etc.), beta14 is redefined as the plain radiation fraction 1-beta — a factor ~4/beta discrepancy handed to callers depending on which EOS regime covered the point. eqstat2 recomputes beta14 fresh at entry (line 463-473), so the EOS itself is self-consistent; the hazard is any external consumer of the i_beta14 result slot (envint_kernel currently just round-trips it). Inherited, so likely tolerated for decades — flagged because the discontinuity across the OPAL table edge is exactly the kind of thing the eqbound ramping tries to prevent for the other outputs (beta14 is NOT ramped, it is recomputed from the ramped beta as 1-beta, still in the wrong convention).

## /Applications/YREC/src/eos/eos_lib.f90:167 -- stale "OPTIONAL ierr" contract comments on mandatory ierr arguments; unreachable "error funnel" blocks
- class: logical (comment/code mismatch, dead code)
- severity: low
- confidence: high
- provenance: new-code
- detail: eos_eval, eos_init, and eos_get_gamma1 all document ierr as OPTIONAL with "when the caller omits ierr ... the same diagnostic messages, then a stop" (lines 167-174, 274, 398), but all three declare plain `integer, intent(out) :: ierr` — omission is impossible and no stop exists anywhere in these routines. The "error funnel" comments (lines 223-227, 362, 472) sit after unconditional `return` statements with no code. Similarly meqos.f90:110-114 keeps the original's unreachable "MHD TABLE FAIL" block (inherited dead there). Nothing misbehaves, but the documented error contract is not the implemented one, which will mislead the stage-3 ierr work.

## /Applications/YREC/src/eos/eos_lib.f90:317 -- eos_init error paths leak the open table file units
- class: logical
- severity: low
- confidence: high
- provenance: new-code (original setups.f STOPped, so the leak did not exist)
- detail: On the fermi-table glitch check (lines 311-320), eos_init now returns ierr=1 with fermi_unit still open (the close is at line 331); the original stopped the program so cleanup was moot. Callers that treat ierr as recoverable (the stated stage-3 direction) will leave the unit busy and a retry of eos_init will fail on the OPEN with status='OLD' on some systems. Same shape applies to mhdtbl.f90: mhdst failure returns before the eight close statements (lines 58-67).

## /Applications/YREC/src/eos/mhd/mhdst.f90:76 -- ZAMS A/B/C and centre 1-5 tables silently share one T-grid; last file read wins, no cross-file consistency check
- class: logical
- severity: low (input-validation gap, fires only with mismatched table files)
- confidence: high
- provenance: inherited (6cd5673:src/mhdst.f passes the same NT1,NR1,NT2,NR2,TLOG1,TLOG2 to all three ZAMS calls and the same NTX,NRX,TLOGX to all five centre calls)
- detail: Each of the three ZAMS mhdst1 calls overwrites zams_lower/upper_num_t/num_r and the shared log10t grids; each of the five centre calls overwrites centre_num_t/num_r/centre_log10t. If the files' grids differ, tables A/B are then interpolated on file C's grid (and centre1-4 on centre5's) with no diagnostic — mhdst1's up/down/centroid T-grid equality check operates only within one file. With zero MHD test coverage, a mismatched table set would produce silently wrong interpolation rather than an error. A cheap equality check when unit>1 tables are read would close it.

## /Applications/YREC/src/eos/mhd/meqos.f90:93 -- mu cross-check uses the carbon-excess mean-weight formula, inconsistent with eqstat2's metal-diffusion branch
- class: physical
- severity: low
- confidence: medium
- provenance: inherited (mu.f and meqos.f identical; eqstat2.f's LDIFZ branch computes AMU per-species while MU always uses the DFX1/DFX4/DFX12 carbon-excess update)
- detail: meqos validates the MHD table's specific gas constant against mu.f90 to 5e-7 relative and aborts (now ierr) on mismatch. mu.f90 assumes excess Z is carbon and no per-species metal scaling; a use_mhd_eos run with use_diffusion_z (where eqstat2 itself would use the per-species formula) can trip this check spuriously — MHD+Z-diffusion is effectively unsupported and fails with a misleading "CHECK MU" message rather than a compatibility diagnostic. mu.f90 itself was verified correct against its own convention (electron inverse-mu backed out of the ideal-gas law checks out algebraically; atomic_weights(4)=0.4995 = 1/2.002 is the fully-ionized-electron term, used only via eqstat2).

## /Applications/YREC/src/eos/mhd/mhdst1.f90:146 -- float equality check on file-read T grids (works only because reads are unformatted)
- class: numerical
- severity: low
- confidence: high (it is a float equality; low that it ever misfires)
- provenance: inherited (6cd5673:src/mhdst1.f lines 430: `IF(TLOG2(JT).NE.TLOGD(JT)) GOTO 600`)
- detail: The up/down/centroid grid-consistency checks compare doubles with .ne.; this is exact-bit comparison of values read from the same unformatted files, so it is safe as long as the tables remain binary and generated together, but becomes a false-failure trap if the tables are ever regenerated through any formatted/converted path. Related: `composition_tolerance = 0.05d0*abs(delta_x)` at line 121 becomes an exact-equality composition check if a table is written with delta_x = 0. Both are input-format landmines in a reader with no test coverage, per the extra-scrutiny request; neither misbehaves with the historical table files.

## Known-item sharpening (do not double-count)

- The known-but-undiagnosed ZAHB NaN in the MLT cubic (tpgrad sqrt-of-negative): two findings above are plausible contributors worth checking against the failing model's (T,P) track: (a) the eqstat2 QDTT/QDTP ramp blend (first finding) corrupts exactly the second derivatives (QDTT/QDTP -> QAT/QAP/QCPT/QCPP consistency) that tpgrad consumes, and only in the band saha_log10t_cutoff-0.5 < logT < cutoff — if shells near the H-burning shell/envelope boundary cross that band around model 950, the over-counted QDTT could push the cubic discriminant negative; (b) the eqrelv non-convergence silent return (fifth finding) would hand tpgrad a stale, inconsistent (QDT,QCP,DELA) set with only a log-file line as evidence — the run log of the failing case should be grepped for "EQRELV: Did not Converge" and "SAHA FAILURE" before hunting further.

## Weak/uncertain observations

- eqstat.f90:296 (eqstat numerical-derivative wrapper): faithful to original including interval choices; the original's KSSAHAx-assigned/KSAHAx-used typo (uninitialized KSAHA copy passed to probe calls) was silently FIXED in modernization (saha_state_local) — behavior improvement, not a bug.
- eqstat.f90:648 etc.: OPAL blends beta then sets beta_inverse=1/beta before checking beta>0 — a ramped beta could in principle be <=0 only if OPAL returned beta<0; negligible.
- saha_eos.f90:130 data-initialized locals (mean_electrons_per_ion, helium_ion_fraction_1/2) are implicitly SAVEd warm-starts across calls — intentional, matches original SAVE; already covered by the known "-Wsurprising static locals" audit item.
- saha_eos.f90:296 non-convergence prints "SAHA FAILURE" then continues with the unconverged E — inherited; only a log line marks the event.
- saha_eos.f90:321 derivative loops run nz1..nz0, skipping fully-ionized species with saha_ratio=1e16 — verified a deliberate negligible-term optimization, not an off-by-one.
- eqscve.f90:214 pressure interpolation nodes are always taken from row idtt while values come from rows idtt..idtt+3 — safe only because the SCV P-grid is column-aligned across T rows; inherited assumption, unchecked.
- eqscvg.f90:307 dlnsmix_dlnt takes log10 of interpolated smix, which is exactly 0 where the envelope helium fraction is 0 and 0 in the unwritten ragged region — -Inf would propagate into cp; unreachable for normal (Y>0) envelope compositions.
- eqscvg.f90 uses the smix table built for the ENVELOPE composition while weighting per-zone X/Y — inherited approximation for off-envelope compositions.
- scv_envelope_table.f90:158/166-169: log_p_work recomputed as log10(Pgas+Prad) and dlnp_dlnt_gas (line 195) are computed and never used — inherited dead statements.
- scv_envelope_table.f90 reads star%envelope_hydrogen_fraction directly rather than eos_mix — envelope-mixture consistency relies on both being synced by starin/eos_set_mixture; a physics-purity leftover, not a numeric bug today.
- setscv/tablex(:,:,12) and tabley(:,:,12) (du/dt on the native grid) are computed in the first loop and never read anywhere — inherited dead work.
- eos_lib.f90 eos_init: SCVHE/SCVZ header rows are read into dummies and ASSUMED to share tlogx/nptsx with the H table (inherited; comment admits it); the tablez read is deliberately reversed in P (inherited).
- eos_init fermi glitch check compares raw DVAL not INT(DVAL) — flagged by MHP 10/02 comment in the original, preserved verbatim; benign for the shipped table.
- fully_ionized_eos.f90:91 jt1 clamped to [1,18] silently extrapolates the Prather table beyond logT 9.0 — inherited.
- mhdpx2.f90: table_selector values map correctly (A/B/C to slots 1-3, centre 1-5 to slots 4-8; lower/upper tables paired) — checked for copy-paste transposition, none found.
- rtab.f90: bounds checks fire on the first record, before any out-of-bounds write — verified safe ordering.
- eqbound modern signature dropped the original's third argument DL0 — verified used only in commented-out debug writes in 6cd5673:src/eqbound.f; harmless.

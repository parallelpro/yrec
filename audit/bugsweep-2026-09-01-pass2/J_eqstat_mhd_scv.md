# Partition J -- EOS dispatcher/blending, MHD EOS tables, SCV EOS

Pass 2, 2026-09-01. Read-only review; no source edited, nothing built or run.

## Files read in full

- /Applications/YREC/src/eos/eqstat.f90 (810 lines: eqstat numerical-derivative wrapper + eqstat2 dispatcher)
- /Applications/YREC/src/eos/eos_lib.f90 (502 lines: eos_get / eos_eval / eos_init / eos_get_gamma1 / eos_set_mixture)
- /Applications/YREC/src/eos/mhd_eos_lib.f90 (98 lines)
- /Applications/YREC/src/eos/mhd/meqos.f90, mhdpx.f90, mhdpx1.f90, mhdpx2.f90, mhdst.f90, mhdst1.f90, mhdtbl.f90, rabu.f90, rtab.f90
- /Applications/YREC/src/eos/scv_eos_lib.f90 (27 lines)
- /Applications/YREC/src/eos/scv/eqscve.f90 (601), eqscvg.f90 (446), scv_envelope_table.f90 (319)
- Also read because they are the interpolation kernels the MHD files call and the assignment names "MHD table interpolation": numerics_lib.f90 `intpt` (lines 1839-1938), `quint` (885-897), `lir` entry (1971-1982); eos/mu.f90 (gas-constant check only)
- Call-site surveys done for every bare external: eqstat (eos_lib:205), eqstat2 (eos_lib:449 + 6 in eqstat.f90), meqos (eos_lib:182), mhdtbl (eos_lib:288), mhdst (mhdtbl:56), mhdst1 (8 calls in mhdst.f90), rabu (mhdst1:90-92), rtab (mhdst1:103-114), mhdpx (meqos:61), mhdpx1 (mhdpx:46), mhdpx2 (mhdpx1:65/71/78), intpt (8 calls in mhdpx2), eqscve (eqstat:503), eqscvg (eqscve:89), build_scv_envelope_table (core/read_starting_model:1185), eqbound/eqbound01/eqbound06 (eqstat.f90 vs eos/opal/eqbound*.f90 dummies). All actual lists match the dummy lists in count, order and type.
- F77 ancestors viewed at 6cd5673: eqstat.f, eqstat2.f, meqos.f, mhdpx.f, mhdpx1.f, mhdpx2.f, mhdst.f, mhdst1.f, rabu.f, rtab.f, intpt.f, quint.f, lir.f (labels only), eqscve.f, eqscvg.f, setscv.f, setups.f (SCV/Fermi read blocks), starin.f (call order).

## Verified clean (re-derived, found correct)

- eqstat numerical-derivative wrapper: step conventions and the returned quantities (cp_dt, grada_dt are dln x/dlnT; dlnrho_dlnt_dt is d(qdt)/dlnT) agree with fully_ionized_eos's analytic conventions and with F77 eqstat.f.
- eqstat2 Saha/relativistic ramp: weight `w = 2(cutoff - logT)` runs 0..1 over exactly the 0.5 dex where both branches are evaluated (need_saha_solution / skip_relativistic_eos brackets at eqstat.f90:379-381); quintic smoothstep has f(0)=0, f(1)=1, f'(0)=f'(1)=0; density, qdt, qdp, cp, grad_ad, cp_dt/dp, grada_dt/dp are blended as relv + f*(saha - relv), i.e. weights sum to 1 (but see finding 4 for the two that are not).
- Diffusion-Z mean-molecular-weight setup and saha_mass_fractions consistency between the diffusion and non-diffusion branches.
- eqbound / eqbound01 / eqbound06 call signatures (the F77 DL0 argument was deliberately dropped from the dummy lists too).
- SCV table read formats in eos_init (eos_lib.f90:334-358) vs setups.f 305-330 (tablez read in reversed column order, nptsx per row); Fermi-table read loop vs setups.f 161-191 (bin_end=261 at grid_idx=41).
- eqscve/eqscvg radiation corrections: qdp = qdp_gas/beta, qpt = qpt_gas*beta + 4(1-beta), qdt = -qpt*qdp, cv_rad = 12 Prad/(rho T), cp = cv + P qdt^2/(qdp rho T), grad_ad = -P qdt/(rho T cp); the entropy-sum form of cp_gas in eqscvg (dS/dlnT = X S_H dlnS_H/dlnT + Y S_He dlnS_He/dlnT + Smix dlnSmix/dlnT); ion-fraction bookkeeping xtf_h_e = 0.5(1 - xH2 - xH1), xtf_he_e = (2 - 2xHe - xHe+)/3.
- Column conventions used by eqscve (tablenv 1..6) and eqscvg (tablex/y cols 4,5,7,8,9,2,3; tablez cols 4,7,10,13) match the way build_scv_envelope_table / eos_init fill them and match the F77.
- Main SCV 4-point stencils never exceed nptsx: nptsx is monotone non-decreasing in logT (30 at 2.10 ... 76 at >= 5.14, verified against input/eos/scv/h_tab_i.dat), and idp <= nptsx(idtt)-3 is enforced on the lowest row of the stencil.
- eos_get_gamma1: cv = cp - P qdt^2/(qdp rho T), Gamma1 = chi_rho cp/cv.
- meqos derivative conversions: QDTT = QDT (QAT + 1 + QDT + QCPT) and QDPT = QDT (QAP - 1 + QDP + QCPP) follow from qdt = -grad_ad rho T cp / P; the dlog->dln factor cnvs on the grad_ad derivatives and its absence on the log cp derivatives are both correct.
- quint is the correct Newton forward-difference quadratic; mhdst/mhdst1/rtab/rabu are byte-faithful to the F77, and the sequence-association of the (79,21,20) scratch arrays onto (10,21,20) dummies for the centre tables is used consistently by reader (rtab) and consumer.
- SCV set-up order: rescale_and_refit_envelope can call eos_get (read_starting_model.f90:786) before update_surface_mixture builds tablenv (line 1185), same order as F77 starin.f (EQSTAT at 554, SETSCV at 827). This is harmless: with tablenv still all-zero, eqscve's pressure search runs off the end of the row and returns valid_table_point=.false., so the Saha result is kept.

---

## src/numerics/numerics_lib.f90:1864-1870 (intpt, called from eos/mhd/mhdpx2.f90) -- MHD table interpolator returns before interpolating whenever logT is inside the table
- class: logical
- severity: high
- confidence: high
- provenance: modernization -- F77 intpt.f lines `DO 100 N=1,NT / IF(TLOG(N).GE.TL) GOTO 101 / IT(1)=N / 100 CONTINUE / 101 IF(IT(1).GE.2) ...` (label 101 is the statement immediately after the loop, i.e. a loop EXIT). Commit 9871628 (2026-08-22, "goto passes 2-5 ... all automated, all gated") rewrote it as `if (table_log10t(n).ge.log10_temperature) then / continue / return / end if`.
- detail: The temperature bracket search in intpt now RETURNS from the subroutine at the first table temperature >= logT, i.e. for every point that lies inside the table (mhdpx1 has already rejected points outside [table_log10t_min, table_log10t_max]). `interp_vars` is never written; mhdpx2 had zeroed it (`call zero(interpolated_vars, ivarx)`), so every ZAMS/centre table returns all-zero variables, mhd_output(1:20) = 0, and meqos computes log10_density = 0, dlnrho_dlnp = 1/0, dlnrho_dlnt = 0/0, then fails its mean-weight consistency check (`specific_gas_constant` = 1/T vs mu's value) and returns ierr=1 -> the eos_get facade stops with "ERROR(MHD) IN MEAN WEIGHTS ... CHECK MU". Observable consequence: any run with lmhd=.true. (use_mhd_eos) dies at its first EOS evaluation; if the mu check were ever relaxed it would instead propagate Inf/NaN. The byte-identical gate that passed this commit only covers the two solar cases, which do not use MHD (no example/testsuite input sets lmhd=.true.), so the regression was invisible. The fix is `exit` instead of `return` (the sibling conversions in the same commit in `lir`, whose labels 101/102 really are `CONTINUE / RETURN`, are correct).

## src/eos/mhd/mhdpx1.f90:56-60 -- Out-of-range logT returns silently with ierr=0 and mhd_output untouched
- class: logical
- severity: medium
- confidence: high
- provenance: inherited -- F77 mhdpx1.f `IF(TL.LT.TMINI .OR. TL.GT.TMAXI) THEN GO TO 999 END IF` with `999 RETURN`; the diagnostic FORMAT 9001 was never referenced.
- detail: When logT lies outside the loaded tables, mhdpx1 returns without setting `mhd_eos%mhd_output` and without setting ierr. mhdpx then copies the STALE `mhd_output(1)` from the previous point into log10_density and meqos derives the whole thermodynamic set from the previous point's table values combined with the current T and P. The header of mhdpx.f90 (lines 13-14) promises "IERR = 1 SIGNALS PGL,TL OUTSIDE THE DOMAIN OF TABLES", and the unreachable block at mhdpx.f90:51-57 ("OUT OF TABLE RANGE") was clearly meant to fire. eos_eval routes EVERY point to meqos when use_mhd_eos is set (eos_lib.f90:181), including atmosphere/envelope points at logT ~3.5-3.9 that are below the ZAMS-table minimum, so the silent branch is reached routinely. In practice meqos's mu check will usually trip on the stale (rho,Pgas) pair and abort with the misleading "CHECK MU" message; the physically-right behaviour is either ierr=1 with the out-of-range message or a fallback to eqstat2. (Only relevant once finding 1 is fixed.)

## src/numerics/numerics_lib.f90:1882-1884 (intpt) -- Pressure above table maximum returns with outputs unset; no check below the minimum
- class: logical
- severity: medium
- confidence: high
- provenance: inherited -- F77 intpt.f `IF(PL.GT.PMAXI) THEN RETURN END IF`; `PMINI` is computed and never used in both versions.
- detail: For each of the 4 temperature rows, if logP exceeds the row's last tabulated pressure intpt returns silently; `interp_vars` keeps mhdpx2's zeros, so the affected composition table contributes all-zero variables to the X interpolation in mhdpx1 (a 1/3 or 1/4 weight of zeros in log rho, chi_rho, cp ...), with no error flag. Below the minimum pressure the code extrapolates from the first cell with no warning either. Consequence: silent corruption of the MHD EOS at the high-density edge (only after finding 1 is fixed).

## src/eos/eqstat.f90:564-572 and 587-588 -- Saha/relativistic ramp adds the raw Saha d(qdt)/dlnT and d(qdt)/dlnP instead of the difference
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited -- F77 eqstat2.f lines 240-248 difference DL0, QDT0, QDP0, QCP0, DELA0, QCPT0, QCPP0, QAT0, QAP0 but not QDTT0/QDTP0; lines 261-262 then do `QDTT = QDTT + WT0*QDTT0`.
- detail: In the 0.5 dex ramp below saha_log10t_cutoff every blended quantity is formed as relv + f*(saha - relv) EXCEPT `dlnrho_dlnt_dt` and `dlnrho_dlnp_dt`, whose stored Saha values are never differenced (lines 564-572 list nine subtractions; the two `saha_dlnrho_dln?_dt` are absent), so lines 587-588 produce relv + f*saha. The weights therefore sum to 1+f instead of 1: at the low end of the ramp (f=1) the value is relv+saha, whereas one step below (skip_relativistic_eos true) it is pure saha -- a discontinuity of the full relativistic value. These are the second derivatives that feed the Henyey Jacobian (want_derivatives=.true. path only), so the converged model is not changed but the Jacobian is wrong by O(100%) for zones with 5.5 < logT < 6.0, degrading convergence there. Fix: add `saha_dlnrho_dlnt_dt = saha_dlnrho_dlnt_dt - dlnrho_dlnt_dt` (and _dlnp_) inside the want_derivatives block before the blend.

## src/eos/scv/eqscve.f90:148 and 202 -- Upward-neighbour smoothing weight is inverted (creates two discontinuities instead of removing one)
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited -- F77 eqscve.f lines 71-80 (`FST = (DTABOVE+TOLST)/(2.0D0*TOLST)` for ISMT=-1, `FST = 0.5D0*DTBELOW/TOLST` for ISMT=+1) and 119-127 for pressure; application at eqscve.f ~ lines 175-185 identical to eqscve.f90:297-309, 355-367.
- detail: The stencil-blending near cell edges is meant to make the piecewise 4-point spline continuous. For a point within tol of the LOWER edge (direction -1) the weight of the own-cell result is (d+tol)/(2tol): 0.5 at the edge, 1 at distance tol -- correct. For a point within tol of the UPPER edge (direction +1) the neighbour weight is 0.5*d/tol applied as own + w*(shift-own): 0 at the edge and 0.5 at distance tol -- the reverse of what is needed (should be 0.5*(1 - d/tol)). Result: at each cell boundary the left side uses 100% of cell n while the right side uses 50/50 of cells n and n+1 (jump of half the stencil difference), and a second jump of the same size appears at d = tol where smoothing switches off. Applies in both logT (tol 0.032 dex, cell ~0.08) and logPgas (tol 0.08 dex, cell 0.2), so 40% of every cell is affected. Observable: small step discontinuities in rho, cp, grad_ad along the profile of every SCV-EOS run (low-mass/brown-dwarf models), i.e. in Gamma1/Brunt output and Henyey residual noise; the same code is copied into the "both directions" branch at lines 419-447.

## src/eos/scv/eqscve.f90:182 and src/eos/scv/eqscvg.f90:136 -- Upward pressure search can set idp = 0 (out-of-bounds tablenv/tablex index) on the first lookup
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- F77 eqscve.f line 109 `IDP = MIN(NPTSX(IDT)-3,JJ)` and eqscvg.f line 85, versus the downward branch which has `IDP = MAX(1,JJ)`.
- detail: The search hints idtt/idp are module variables that start at 0. On the first SCV lookup jjj = min(nptsx-3, 0) = 0 and the upward loop starts at j = 2; if 4.0 <= logPgas < tablenv(idtt,2,1) = 4.2 it exits with jj = 0 and idp = min(nptsx-3, 0) = 0 (no max(1,.) clamp, unlike line 165). Every subsequent stencil then reads tablenv(idtt, 0:3, *) -- index 0, which with column-major storage addresses memory before the array (silent garbage; abort under -fcheck=bounds). The same hazard exists in eqscvg.f90:136. It only fires if the very first SCV evaluation of the process lands in that one 0.2 dex pressure cell (after any call idp >= 1 and build_scv_envelope_table leaves idp = 74), which is why it has not been seen; the temperature search has the clamp (`idtt = max(1,ii)`) and this one should too.

## src/eos/eqstat.f90:549-551 and 578-580 -- Ion-fraction blend across the Saha ramp is a no-op (Saha values overwritten before blending)
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- F77 eqstat2.f lines 230-232 (`FXION(1)=1, FXION(2)=0, FXION(3)=1` before EQRELV) and 253-255 (`FXION(1) = 1.0D0 + WT0*(FXION(1)-1.0D0)` ...).
- detail: Unlike density and the derivatives, the Saha ion fractions are not saved before the fully-ionized values (1,0,1) are stored into `ion_fraction`; the "blend" at 578-580 is then 1 + f*(1-1) = 1, 0, 1 for every f. The intent (a smooth transition of FXION between eqsaha and the fully ionized limit) is defeated: ion_fraction jumps from the Saha values to (1,0,1) at logT = cutoff-0.5 instead of ramping. Consumers are condopacpint (conductive-opacity ion weights) and the atm/env ion-fraction output arrays. Physically small in the default cutoff region (logT 5.5-6.0, where H and He are essentially fully ionized anyway) but it is a copy-paste omission that grows if saha_log10t_cutoff is lowered.

## src/eos/scv/scv_envelope_table.f90:287-291 -- Ragged-edge clamp for the du/dT temperature derivative is dead (overwritten on the next line)
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- F77 setscv.f lines 212-217 (`jj = min(NPTSX(ii),j) / QR(K) = TABLENV(ii,jj,6) / QR(K) = TABLENV(IDT+K-1,J,6)`); somebody added the clamp and left the original assignment after it.
- detail: `jj = min(nptsx(ii), p_idx); interp_x(k_idx) = tablenv(ii,jj,6)` is immediately overwritten by `interp_x(k_idx) = tablenv(idtt+k_idx-1, p_idx, 6)`, so rows idtt+1, idtt+2 are read at p_idx even where p_idx > nptsx(row) -- never-written (zero) entries feed the derivative. The same unclamped pattern is used for dsmix_dlnt (lines 202-204) and dqdt_dlnt (270-272). Rows are filled in increasing logT and nptsx is non-decreasing, so for the FORWARD stencil this only hits entries above nptsx(idtt) which are themselves never used by eqscve; consequence today is confined to tablenv columns 7-12, which eqscve does not read. Reportable as intent/code mismatch and as a trap for anyone who starts using those columns.

## src/eos/scv/eqscvg.f90:307-310 -- log10 of the interpolated entropy of mixing, which is exactly 0 for X=0 or Y=0 envelope compositions
- class: numerical
- severity: low
- confidence: medium
- provenance: inherited -- F77 eqscvg.f lines 214-216 `QSMT = FTD(1)*LOG10(TEMP(1,3)) + ...`.
- detail: build_scv_envelope_table sets smix(t,p) = 0 whenever the H or He number density of the ENVELOPE mixture is zero (scv_envelope_table.f90:101-104), i.e. for Xenv = 0 or Xenv = 1 (Y computed as 1-Xenv). eqscvg then forms dlnsmix_dlnt from log10 of the pressure-interpolated smix -> -Inf, and 0*(-Inf) = NaN propagates into cp and grad_ad. Only an issue for pure-H or H-free envelopes run with the SCV EOS and a zone composition different from the envelope's (the eqscvg branch), so low impact; but no guard exists.

## src/eos/mhd/mhdpx1.f90:101-107 -- X-interpolation reads table_vars(1:3, 21:25), which the ZAMS tables never fill (uninitialized stack in the modern code)
- class: numerical
- severity: low
- confidence: medium
- provenance: modernization (weakly) -- F77 mhdpx1.f has `SAVE` so VAROUT was static (zero-filled); mhdpx1.f90:34 declares `table_vars` as a plain automatic local.
- detail: ZAMS tables carry ivarc = 20 variables (mhdpx2 fills table_vars(1:3, 1:20)), but the quint loop runs iv = 1..ivarx = 25, so mhd_output(21:25) is interpolated from garbage. meqos does not read elements 21-25, so the result is only wasted work -- unless the compiler traps on signalling NaN/denormal garbage (-ffpe-trap) in quint. Bounding the loop at ivarc for the ZAMS branch removes the read.

## Weak/uncertain observations

- eqstat.f90:473 defines beta14 = 4(1-beta)/beta in the Saha/relv path but the OPAL overlays (632, 648, 699, 715, 773, 789) and meqos (84) define it as 1-beta. Inherited (eqstat2.f 135 vs 288). No consumer outside eos reads i_beta14 (only envint_kernel copies it), so no observable effect today.
- eqstat.f90:466-472: beta < 0 (Prad > P, possible for trial (T,P) during envelope/atmosphere iteration) only prints a warning, then 1/beta and, on the SCV path, log10(beta*pressure) (eqscve.f90:97) -> NaN. Inherited. Not linked to the He-exhaustion MLT NaN symptom by evidence in this partition.
- eqstat.f90:560-598: the ramp derivative term d(ramp_factor)/dlnT is omitted (F77 line 239 `C WT1 = ...` commented out), so the blended cp_dt/grada_dt are the blend of derivatives, not the derivative of the blend. Inherited, deliberate-looking approximation.
- eqstat.f90:569-572: when want_derivatives is false the six saha_*_dt/_dp locals are subtracted while never assigned (F77 had SAVE'd stale values). Results are discarded; only matters under FPE trapping on uninitialized NaN.
- mhdpx1.f90:62-82 vs 82: table selection uses `< zams_centre_boundary` but X-interpolation uses `<= zams_centre_boundary`; at logT exactly equal to centre_log10t(1) the quint branch runs on table_vars(1:3) that the centre-table branch did not fill. Inherited (LT/LE in mhdpx1.f); measure zero.
- eqscve.f90:257-262 (direction -1) reads row idtt-1 at idp..idp+3 without a ragged clamp; rows with fewer pressure points (e.g. nptsx 38 -> 49 at logT 3.54/3.62) contribute never-written entries at 50% weight. Physically unreachable region (logPgas > ~11.4 at logT ~3.6). Inherited.
- scv_envelope_table.f90:82: helium fraction for the mixing entropy is 1 - Xenv (Z ignored). Inherited approximation.
- eqscvg.f90: the entropy of mixing is tabulated for the ENVELOPE composition yet eqscvg is used precisely when the zone composition differs from it; and eqscvg has no edge smoothing at all while eqscve does. Inherited design.
- eqstat.f90:494-506: no blend between SCV and Saha at the SCV table boundary (logPgas = 4.0, logT edges); the EOS switches abruptly (valid_table_point is not consulted by eqstat2 -- intended fallback to Saha). Inherited design.
- mhdst1.f90:145-160 compares log10t2 with log10t_down/up by exact floating equality (`.ne.`); fine as long as the tables are written from identical grids. Inherited.
- eos_lib.f90 OPAL path (documented) falls through with the previous point's opal eos_output_06 on an esac06 alternate return -- preserved behaviour, out of this partition's remit to judge.

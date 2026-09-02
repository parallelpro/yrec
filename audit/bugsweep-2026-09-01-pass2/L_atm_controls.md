# Pass-2 sweep L: atmosphere boundary conditions/tables, control-file parsing, logging

Date: 2026-09-01. READ-ONLY review; nothing built or run. Earlier-pass reports were not consulted.

## Files read in full

- /Applications/YREC/src/atm/atm_lib.f90 (229 lines)
- /Applications/YREC/src/atm/atm_table_lib.f90 (92)
- /Applications/YREC/src/atm/atmstruct_lib.f90 (38)
- /Applications/YREC/src/atm/envstruct_lib.f90 (41)
- /Applications/YREC/src/atm/ttau_lib.f90 (131)
- /Applications/YREC/src/atm/tables/alfilein.f90 (416)
- /Applications/YREC/src/atm/tables/alsurfp.f90 (302)
- /Applications/YREC/src/atm/tables/altabinit.f90 (147)
- /Applications/YREC/src/atm/tables/kcsurfp.f90 (187)
- /Applications/YREC/src/atm/tables/surfp.f90 (183)
- /Applications/YREC/src/io/read_controls.f90 (2508)
- /Applications/YREC/src/io/controls_lib.f90 (810)
- /Applications/YREC/src/io/check_controls.f90 (80)
- /Applications/YREC/src/io/run_log.f90 (100)
- /Applications/YREC/src/io/luout_lib.f90 (28)
- /Applications/YREC/src/io/print_allard_tables.f90 (66)

Read for context (not assigned, not fully audited): src/io/inlist_new_read.inc, src/io/inlist_new_decl.inc, src/core/atmosphere_derivs.f90, src/core/envint_kernel.f90 (lines 160-300), src/core/yrec_capi.f90, src/core/yrec_reset.f90, src/wind/massloss.f90 (lines 235-275), src/numerics/numerics_lib.f90 (inter3/kspline/ksplint/locate/polint), src/defaults/controls_registry.tsv (spot checks). F77 ancestors consulted at 6cd5673: surfp.f, kcsurfp.f, alsurfp.f, altabinit.f, alfilein.f, parmin.f (there is no atmos.f/ttau.f/readin.f at that revision; the T(tau) statement functions lived in the callers, and parmin.f is the read_controls ancestor).

## Verified clean (re-derived, found correct)

- T(tau) relations (ttau_lib.f90): Eddington log T = log Teff + 0.25 log10(0.75 (tau + 2/3)); the constant -0.031235 = 0.25 log10(0.75) to 6 digits. Krishna-Swamy 1966 T^4 = 0.75 Teff^4 [tau + 1.39 - 0.815 e^(-2.54 tau) - 0.025 e^(-30 tau)] with the tau=0 start value 0.550 = 1.39 - 0.815 - 0.025. HSRA polynomial is offset by harvard_t_tau(2/3) so T = Teff at tau = 2/3. Photosphere x_limit values: Eddington log10(2/3) = -0.176091; KS root of the bracket = 4/3 at tau = 0.312 (log10 = -0.5056) checked numerically.
- atmosphere_derivs.f90 (context): d log10 P / d log10 tau = g tau / (kappa P) is the correct log-log form of dP/dtau = g/kappa (ln10 factors cancel).
- surfp.f90 / kcsurfp.f90: goto-elimination vs surfp.f/kcsurfp.f is equivalent (Teff row search fall-through to row=nt; the `if (k .lt. gmin)` "first four points" guard reproduces the F77 loop fall-through; zero-trip loop semantics identical). kcsurfp differs from surfp only by table names/sizes (ntc=76 vs nt=57, ngc=ng=11). Sentinel test `.le. 0` is safe: no legitimate log10 P <= 0 exists in any shipped Kurucz/Castelli table (checked all input/atmos/kurucz/*.tab and CastelliKurucz/*.tab); -999 sentinels occur only on the low-gravity side of each row.
- atm_init (atm_lib.f90) Somers gmin/gmax index scan matches the F77 logic; `ng` vs `ngc` at line 182 is harmless (both 11).
- alsurfp.f90 vs alsurfp.f: 4x4 Lagrange (locate + polint) faithful; asymmetric handling (Teff above max -> lookup_failed, gray fallback; Teff below min -> fatal) is by design and identical to F77. Call sites (envint_kernel.f90:235/245/254, massloss via atm_get_surface_pt, test_atm) match dummies.
- altabinit.f90 / alfilein.f90 / print_allard_tables.f90 argument lists and the min/max/limit bookkeeping match the F77.
- No analytic dP/dTeff or dP/dg exists in the assigned files: the surface-pressure routines return log10 P only, and the solver's boundary derivatives come from the envelope triangle (surfbc.f90:151-154, finite differences over tri_delta_teffl / tri_delta_logl). Nothing to check there in this partition.
- read_controls.f90 adopt_canonical_names (1504-1751), the post-override re-syncs (1867-1871, 1900-1903, 2098-2105, 2250-2264) and echo_settings: every legacy-name -> canonical-name copy was checked pairwise (about 240 assignments); no crossed pairs found. The K97 branch reproduces the old wind law exactly (exr = 2 + 5m - 4mc with m = alfa/3, c = 2 gives 2 - alfa; exm = -alfa/3; omega exponent 3 = 1 + 4 alfa/3 + 2(1 - 2 alfa/3)).
- Units: calsolage / endage / setdt are in yr everywhere (stop_conditions compares target_end_age with dage*1e9; setup_solar_calibration sets 1.0d8 and target_solar_age); tdisk is Gyr and compared with disk_gate_age_gyr; clsun/crsun cgs; wmax/wmax_sun/pmmsolw rad/s. No log/linear or Gyr/yr crossing found.
- Legacy .nml1/.nml2 vs new-style &star_job/&controls: the generated inlist_new_read.inc seeds every new-style name from the legacy local, reads, and maps back to the same legacy locals before adopt_canonical_names runs, so both paths flow through identical post-processing. Every legacy namelist member absent from the new-style groups (hpttol, ies, imu, fstch, sstandard) is documented DEAD in controls_registry.tsv and is referenced nowhere.
- controls_check_lib rules match the gating code they cite; run_log_lib throttling/header logic is correct (final line dedup via last_printed_model).

---

## atm/atm_lib.f90:220 -- atm_get_surface_pt never propagates alsurfp's error to ierr
- class: interface
- severity: medium
- confidence: high
- provenance: new-code (facade added 2026; F77 massloss called ALSURFP directly, whose failure path was a STOP)
- detail: `ierr = 0; jerr = 0; call alsurfp(..., jerr); if (jerr == 0) then; return; end if` -- there is no assignment `ierr = jerr` in the non-zero branch, so ierr is always 0. When alsurfp hits its fatal branch (log Teff <= allard_teffl_min, alsurfp.f90:161-181) it prints "Program Terminated" and returns ierr=1, but the caller wind/massloss.f90:247 sees ierr=0, skips its `if (ierr /= 0) return`, and continues with whatever stale atm_table%atm_log10_pressure/temperature the last envelope integration left. The same caller also ignores the lookup_failed flag entirely, so the Teff-above-table / logg-out-of-range case (alsurfp returns without touching atm_table) likewise proceeds with stale surface P,T. Fires only in the accretion/mass-loss entropy path with atm_choice = 4 (Allard), but there it silently corrupts the accretion entropy term.

## atm/tables/surfp.f90:74-96 and kcsurfp.f90:73-97 -- hot-Teff gravity branches ignore the -999 sentinel and interpolate garbage
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited (surfp.f:51-64 has the identical `IF(ATMTL(J).GT.4.55D0)` linear / 3-point branches using ATMGL(NG-1),ATMGL(NG) and ATMPL(J,NG-2..NG) with no -999 check; the Somers 5/15 gmin/gmax guard was only added to the general branch)
- detail: For any table row with log Teff > 4.55 the code linearly interpolates between columns ng-1 and ng, and for 4.5 < log Teff <= 4.55 it uses a 3-point Lagrangian on columns ng-2..ng, without consulting kurucz_gmin_index/castelli_gmin_index. Checked against the shipped tables: in EVERY CastelliKurucz atmk2004*.tab the last row (log Teff 4.699, 50 kK) has P(10) = -999 with only P(11) valid, so the linear branch produces log P ~ -999*(1-fx) + 4.4*fx; in kurucz atmk1990p05/p075/p10 rows 56-57 the linear branch and in m20/p03/p05/p075/p10 rows 52-55 the 3-point branch consume -999 as well. The Teff spline stencil row_base = min(nt-3, ...) includes the last row for any log Teff above kurucz_castelli_teff_table(74) (~4.69), so with atm_choice = 5 a model hotter than ~48-49 kK gets a silently absurd surface pressure (and the envelope integration starts from it) instead of an error. Cold-side hot stars (massive MS, post-AGB) are the affected runs.

## atm/tables/alfilein.f90:81,132 -- allard_use_tau100 guard tests an uninitialized local; the intended flag is never checked
- class: logical
- severity: medium
- confidence: high
- provenance: inherited typo (alfilein.f:104 `If(LATMTPTau100)` vs the COMMON flag `LALTPTau100` at alfilein.f:47/53), made worse by modernization (the header comment at lines 74-80 claims the local is "always-default-valued", but `logical :: latmtptau100` has neither an initializer nor SAVE, so it is stack garbage)
- detail: The old-format Nextgen.p00 branch is supposed to abort when the user asks for tau=100 P,T (star%ctrl%allard_use_tau100) because those columns do not exist in the 1999 files. The test reads a distinct, never-assigned local instead. Two consequences: (1) with a compiler that does not zero the stack the abort can fire spuriously ("Invalid old Allard input file for requested PT,TL at Tau=100") on a perfectly valid run; (2) when allard_use_tau100 = .true. really is set with the old-format file, the guard never fires and alsurfp.f90:243-257 interpolates the -999-filled allard_log10_pressure_tau100/temp_tau100 arrays, yielding log P = log T = -999 at the surface. The fix is a one-token rename to star%ctrl%allard_use_tau100 (which is what the comment says the author meant).

## io/read_controls.f90:1082 (and every other DATA / initialized legacy local) -- re-entrant runs inherit the previous run's namelist values, contradicting the "pristine defaults" contract
- class: logical
- severity: medium
- confidence: high
- provenance: new-code (the four-step pristine-defaults design at read_controls.f90:9-19 and the re-entrancy contract at core/yrec_capi.f90:17-19 / core/yrec_reset.f90:13-14 are 2026 additions; the F77 parmin ran once per process)
- detail: read_controls resets star%ctrl and re-seeds the controls_lib BUFFER before every read, but the legacy-spelled namelist locals inside read_input are DATA-initialized (`data kttau/0/`, `data clsun,crsun/...`, `data wmax,wmax_sun/...`, `data awind,pmma.../...`, `data tdisk,pdisk,ldisk/...`, `data atmerr.../...`, `data zxmix,frac_c.../...`, `data toll,tolr,tolz,lcals,...,calsolage,calsolzx/...`, the s0_* cross sections, etc.) or declared with initializers (alphac, alphae, dtdif, wnew, niter1-4, lthoul, lsemic, lmonte, ...), which gives them the SAVE attribute. A second in-process yrec_run whose inlist omits, say, atm_choice keeps the first run's kttau (the legacy read only overwrites what the file provides, and the new-style path seeds `atm_choice = kttau` from the stale local at inlist_new_read.inc:176 before reading). Every such control then flows into star%ctrl via adopt_canonical_names. The re-entry acceptance test (core/test/test_reentry.f90) runs the SAME inlist twice, so it cannot detect this. Observable: a pyyrec grid that varies inlists between calls silently carries over any control the later inlist does not restate (atmosphere choice, solar constants, wind law, overshoot, tolerances...).

## atm/tables/altabinit.f90:88 -- low-gravity extrapolation limit is four column widths, not "one column's width" as documented
- class: physical
- severity: low
- confidence: high
- provenance: inherited (altabinit.f:78 `GLmin(i) = GLs(j1) - 4D0*(GLs(j1+1) - GLs(j1))`, with the same "one column's width" comment at altabinit.f:14/70)
- detail: allard_gl_row_max is one column width above the last valid column, but allard_gl_row_min is set 4*(dgl) BELOW the first valid column (2 dex for the 0.5-dex NextGen grids: a table starting at log g = 3.5 accepts log g = 1.5). alsurfp.f90:138 gates only on the overall min/max, so a pre-MS / giant model with log g well below the table is silently handled by 4-point Lagrange extrapolation two dex outside the data instead of triggering lookup_failed and the gray fallback. The comment/code mismatch makes the asymmetry look unintended; if the 4-column margin is deliberate the header should say so.

## atm/tables/surfp.f90:46 and kcsurfp.f90:45 -- no upper Teff / log g range guard; spline extrapolates above the tables
- class: numerical
- severity: low
- confidence: high
- provenance: inherited (surfp.f:25 `IF(TEFFL.LT.3.5D0 .OR. GL.LT.-0.5D0)` only)
- detail: Only the low side is guarded. For log Teff above the last table row the Teff row search falls through to row = nt and the 4-point cubic spline in Teff extrapolates; for log g above kurucz_logg_table(ng) = 5.0 the `log10_gravity .ge. logg(ng-1)` branch extrapolates the natural spline through the top four columns (and for the hot rows the linear branch extrapolates too). No warning is written. Combined with the previous surfp finding this is how a hot model gets an unphysical boundary pressure with no diagnostic; a symmetric `.gt.` guard (or at least a run-log warning) would match the Allard routine's behaviour.

## atm/tables/alsurfp.f90:280-293 -- two run-log lines on every envelope call when printing is NOT requested
- class: logical
- severity: low
- confidence: high
- provenance: inherited (alsurfp.f:150ff writes the same block to ISHORT under `.not. LPRT`), but it now contradicts the 2026 run-log design (run_log.f90 header: the run log is meant to be small and byte-pinnable)
- detail: envint_kernel.f90:254 always passes print_flag = .false., and the `if (.not. print_to_files)` block writes formats 74/75 to run_log_unit unconditionally. With atm_choice = 4 that is two lines per envelope integration (three integrations per triangle, several per Henyey iteration), i.e. tens of thousands of lines per run in what the log redesign calls the compact run log. It should be gated on solver_diagnostics() like the other envelope forensics, or dropped.

## atm/tables/surfp.f90:123-158 (kcsurfp.f90:126-162) -- general-case gravity search leaves a gap below the top valid column and reads out of bounds when a row is all -999
- class: logical
- severity: low
- confidence: medium
- provenance: inherited (surfp.f general-case loop after the Somers 5/15 IMINMAX edit; the "first 4 points" fallback indexes ATMGL(K+IMIN-1))
- detail: The loop `do k = gmax-3, gmin, -1` brackets log g in [logg(k+1), logg(k+2)), so its highest bracket ends at logg(gmax-1); the separate top branch (line 99) only starts at logg(ng-1). If gmax < ng (valid data missing at the HIGH-g end) any log g in [logg(gmax-1), logg(ng-1)) matches no bracket, the loop completes with k = gmin-1, and the `k .lt. gmin` fallback interpolates with columns gmin..gmin+3 -- the wrong end of the row. Separately, atm_init sets gmin = ng when a row is entirely -999 (atm_lib.f90:182) "so the code should break", but the fallback then reads kurucz_logg_table(k+ng-1) for k = 1..4, i.e. indices 11..14 of an 11-element array (silent out-of-bounds read, garbage log P) rather than an error; gmax_index is also left unassigned for such a row. Neither case is reachable with the shipped tables (all sentinels are on the low-g side, no row is all -999), so this is a latent hazard for user-supplied tables.

## io/read_controls.f90:2069-2070 vs 2027-2045 -- exr sign contradicts the header derivation (code is right for K97, comment is wrong)
- class: physical
- severity: low
- confidence: medium
- provenance: unclear (the derivation comment and the "JvS 09/25 REMOVED TYPO" edit are post-6cd5673; parmin.f:970-997 carried the K97/V13 block with a different EXR expression)
- detail: The header says JDOT ~ R^(2+5m+4mbc) (and, ignoring GM/R, R^(2+5m+4mbc)); the code sets `exr = 2 + 5m - 4mc` (b dropped, sign flipped). For awind='K97' the code value reproduces the historical Chaboyer/Krishnamurthi law exr = 2 - alfa exactly, so the code is the trusted side for K97; for 'V13' c = 0 so the term vanishes; only 'CUS' users with pmmc /= 0 and pmmb /= 1 see a formula that differs from the documented one. Reported as a comment/intent mismatch rather than a numerical bug -- worth having the wind-law owner confirm which sign is meant for the CUS case.

## Weak/uncertain observations

- read_controls.f90:2049,2057,2077: `pmmm = alfa/3.0`, `pmmm = 0.22`, `c_2 = 0.0506` -- single-precision literals into double (0.22 -> 0.21999999880); also DATA `pmmjd 1.32e30, pmmmd 1.27e12, pmmsolw 2.836e-6` (line 1197) and `frac_o 0.468195` (line 1031) without d exponents. Relative error ~1e-8 in wind-torque normalisation and the GS98 O fraction; inherited from parmin.f. Not in the burn/net/turnover set already fixed.
- read_controls.f90:2407: for kind=1 (evolve) cards rescale_params(3,nkind) keeps its default 0, which satisfies `.ge. 0`, so star%envelope_metal_fraction is set to 0 and eos_set_mixture is called with Z=0 during control parsing; inherited (parmin.f:1378 `IF(RESCAL(3,NKIND).GE.0.0D0) ZENV=...`), harmless because read_starting_model.f90:1143/1175 overwrites both before use, but the format text "0 or negative = keep current" does not describe what this line does.
- read_controls.f90:1824-1831: last_slash_idx is uninitialized if the log path has no '/' (a deck setting log_output_file = 'run.log'); the mkdir string then uses garbage. New-code (2026 output-dir logic). The shipped default '{YREC_OUTPUT}/run.log' always contains a slash.
- controls_lib.f90:745: exmd/extau/exr/exm/exl/expr/constfactor/excen/c_2/wind_law_omega_exponent have no initializers and are assigned only inside `if(lrot)` in echo_settings; a non-rotating run stores module-static zeros into star%ctrl, fine unless any wind-torque path runs with rotation off.
- alfilein.f90: FeH/alpha count scans are bounded by allard_num_gl (5) rather than by a FeH/alpha dimension (inherited; no effect for shipped files with <= 5 values); exact float compare `dabs(alpha_value - target_alpha) .gt. 0d0` where the first pass used 1d-5 and the second 1d-6 (inherited).
- alsurfp.f90:243-257: polint is called with an array element as the xa actual and `xa(20)` dummy (sequence association over the trailing 4 elements) -- legal, but fragile if the grid arrays are ever made allocatable/non-contiguous.
- print_allard_tables.f90 recomputes allard_teffl_min/max from the grid instead of using the values altabinit stored (harmless, identical result).
- read_controls.f90:1358-1363 vs 1453-1458: when control_nml_override is set but physics_nml_override is blank (single-file pyyrec call), getarg(2) of the HOST process is consulted for the physics file -- contradicts the "must not read the host's command line" note at yrec_capi.f90:13-14; only matters for paired-style inlists driven from a host with argv.
- nml_file_has_group (read_controls.f90:2429) matches `&controls` by prefix, so a legacy file whose first non-blank token on a line begins with `&controlsx` would be misdetected; no such file exists.

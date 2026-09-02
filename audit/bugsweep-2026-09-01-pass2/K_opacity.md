# Bug sweep pass 2 -- K: opacity (src/kap/)

Date: 2026-09-01. Read-only review; no file edited, nothing built or run.
F77 ancestors taken from git rev 6cd5673 (src/*.f).

## Files read in full

- /Applications/YREC/src/kap/kap_lib.f90
- /Applications/YREC/src/kap/setupopac.f90
- /Applications/YREC/src/kap/surfopac.f90
- /Applications/YREC/src/kap/opacity_table_lib.f90
- /Applications/YREC/src/kap/conductive_table_lib.f90
- /Applications/YREC/src/kap/conductive/condopacp.f90
- /Applications/YREC/src/kap/conductive/condopacpint.f90
- /Applications/YREC/src/kap/alex06/readalex06.f90, alex06tab.f90, getalex06.f90
- /Applications/YREC/src/kap/alex94/read_alex94_tables.f90, alex94_surface_table.f90, alex94_fixed_z_table.f90, alex94_interp3d.f90
- /Applications/YREC/src/kap/kurucz90/read_kurucz_tables.f90, build_kurucz_splines.f90, kurucz.f90, kurucz2.f90
- /Applications/YREC/src/kap/laol89/rdlaol.f90, rdzlaol.f90, sulaol.f90, zsulaol.f90, gtlaol.f90, gtlaol2.f90, gtpurz.f90
- /Applications/YREC/src/kap/opal92/opal92_read_tables.f90, opal92_surface_table.f90, opal92_fixed_z_table.f90, opal92_interp3d.f90, opal92_interp3d_z2.f90, opal92_interp2d.f90, opal92_interp2d_z2.f90
- /Applications/YREC/src/kap/opal95/ll95tbl.f90, getopal95.f90, opal95_fixed_z_table.f90, opal95_surface_table.f90, opal95_interp2d.f90, opal95_interp3d.f90, opal95_interp4d.f90
- Supporting (read for the algorithms the kap code depends on): /Applications/YREC/src/numerics/numerics_lib.f90 (findex, interp, intrp2, cspline, splint, locate, ysplin, intpol, hunt), /Applications/YREC/src/util/xrng4.f90
- F77 ancestors consulted: getopac.f, getopal95.f, ll95tbl.f, op952d/3d/4d.f, op95xtab.f, op95ztab.f, yllo3d.f, yllo3d2.f, yllo2d.f, yllo2d2.f, readalex06.f, alex06tab.f, getalex06.f, alxtbl.f, alxztab.f, alx8th.f, yalo3d.f, condopacp.f, condopacpint.f, gtlaol.f, gtpurz.f, sulaol.f, rdlaol.f, ll95tbl.f

Call-site check (grep -rn "call <name>" over /Applications/YREC/src, all hits read) done for every bare external in the assignment: getalex06, alex94_interp3d, alex94_surface_table, alex94_fixed_z_table, kurucz, kurucz2, gtpurz, gtlaol, gtlaol2, getopal95, opal95_surface_table, opal95_fixed_z_table, opal92_interp3d, opal92_interp3d_z2, opal92_surface_table, opal92_fixed_z_table, condopacpint, condopacp, setupopac, surfopac, ll95tbl, rdlaol, rdzlaol, sulaol, zsulaol, readalex06, read_alex94_tables, read_kurucz_tables, build_kurucz_splines. Actual argument lists match the dummy lists in count, order and type (including the `*100` alternate return of kurucz/kurucz2, whose label-100 path correctly leaves got_atmosphere_opacity=.false.).

## Verified clean (re-derived, found correct)

- 4-point Lagrange weights and their derivatives in numerics_lib interp/intrp2; all OPAL95 2D/3D/4D assemblies (weights applied to log10 kappa in logT and logR; dlnkap/dlnT|rho = dlnkap/dlnT|R - 3 dlnkap/dlnR applied exactly once, at the end).
- OPAL95 low-R and high-R linear extrapolation: derivative at the edge node times (logR - logR_edge), per T row, then interpolated in T.
- OPAL95 X-stencil selection in getopal95 (tables 9 and 10 dropped correctly for Z>0.04 / Z>=0.1; node 4 abscissa = 1-Z when table 10 is used) and Z-stencil clamp `min(index_z, num_z-3)`.
- OPAL95 slot arithmetic for X=1-Z tables: producer (ll95tbl) and consumers (interp4d, fixed_z_table, surface_table) use the same `start(iz)+ix` convention, so X=1-Z lookups land on the right table (see finding 1 for the X=0 end).
- OPAL92 natural cubic spline ysplin (boundary rows s1 = 1.5 f2 - s2/2, coefficient recovery c4 = d3/h^2, c3 = (f - s_i - d3)/h), intpol (subset spline + bisection), findex; yllo3d/yllo3d2 and yllo2d/yllo2d2 mirrors are byte-equivalent in logic.
- LAOL89: cspline/splint (NR), locate, xrng4 4-point window; gtlaol/gtlaol2/gtpurz are structurally identical copies.
- Kurucz90: table read (T ascending, 56 rows), log10 taken once in build_kurucz_splines, density-row search, alternate-return-on-out-of-table semantics (mirrors getopac.f).
- ALEX94: 105-table slot ii = ix + 7*(iz-1), fixed-Z table and surface-X slot 8, alex94_interp3d index ranges.
- ALEX06: 143-table storage jj = (i-1)*16+ii and jj = 128+ii for the X=1-Z set; the kk4=128 special case in alex06tab for (index_x=6, iz=13); readalex06 header check.
- Conductive/radiative blend in kap_eval: kappa = kr*kc/(kr+kc); dlnkappa = r' + c' - (kr r' + kc c')/(kr+kc) (both derivatives) -- correct.
- Hubbard-Lampe fallback in kap_eval (formula and both derivatives) and its selection region.
- condopacp: cinterp3 Hermite interpolation, hunt bracketing, the (-3.5194 = log10(16 sigma/3)) constant, kappa_cond = 16 sigma T^3 /(3 rho sigma_c).
- Atmosphere/interior ramp in kap_eval: weight w = (T_max - logT)/(T_max - T_min) applied consistently to kappa (see weak list for the derivative approximation).

# Findings

## /Applications/YREC/src/kap/opal95/ll95tbl.f90:80 -- OPAL95 X=0 tables at Z=0.08 and Z=0.10 are overwritten by X=1-Z tables (slot collision)
- class: logical
- severity: high
- confidence: high
- provenance: inherited -- ll95tbl.f line 60 `N = NZ(IZ)+IX`, line 33 `DATA NZ/0,10,...,100,109,118/`, line 35 `DATA NUMXM/10x10,9,9,8/`
- detail: Slot is `n = opal95_table_start_index(iz) + ix` with start = [0,10,...,100,109,118] and ix = 1..10 (ix=9 read only for Z<=0.04, ix=10 only for Z<=0.08; the file is X-outer, Z-inner so ix=10 is read last). For Z=0.06 (iz=11, start 100) ix=10 lands on slot 110, which is also (Z=0.08, ix=1, X=0); for Z=0.08 (start 109) ix=10 lands on slot 119, which is also (Z=0.10, ix=1, X=0). Because ix=10 is read last, the X=0 tables at Z=0.08 and Z=0.10 are destroyed and replaced by the X=0.94/Z=0.06 and X=0.92/Z=0.08 tables, while slots 109 and 118 are never written (contain whatever the array was initialised to). Consumers (opal95_interp4d `xz_table_index = opal95_index_x(k,j) + start(index_z+k-1)`, opal95_fixed_z_table, opal95_surface_table) use the same convention, so X=1-Z lookups are right but every X=0 lookup at Z=0.08 or 0.10 reads a hydrogen-rich table. Fires in the 4D path for X<0.2 (X stencil 0,0.1,0.2,0.35) whenever the Z stencil contains 0.08 (Z>0.04; weight ~-0.05 at Z=0.05, growing toward Z=0.08/0.10), and with weight 1 in the HB branch of kap_eval (Z>0.1, logT>7.0) which calls getopal95 with `table_metal_fraction = 0.1d0` exactly: there kappa(X=0,Z=0.1) is kappa(X=0.92,Z=0.08) -- at electron-scattering temperatures 0.38 instead of 0.20 cm^2/g, i.e. ~0.28 dex too high, blended linearly toward the pure-Z table as Z->1 (still ~0.15 dex at Z=0.5). Every core-He-burning model computed with OPAL95 is affected (core opacity, hence convective-core size, HB luminosity/lifetime), and the discontinuous, overestimated core opacity at core He exhaustion is a plausible contributor to the run_from_zahb_to_tahb GN93 NaN-in-MLT symptom near model 950. Solar-Z fixed-Z tables (nodes 0.01..0.04) are unaffected, so MS/RGB calibrations do not see it.

## /Applications/YREC/src/kap/conductive/condopacpint.f90:114-120 -- mixture conductive-opacity derivatives are wrong (multiplied by log10 conductivity; dlnkap/dlnT uses the rho derivatives)
- class: physical
- severity: medium
- confidence: high
- provenance: inherited -- condopacpint.f lines 61 and 65 `QODC = CONDL*(WtH1*QODCX+WtHe4*QODCY+WtOx*QODCZ)`, `QOTC = CONDL*(WtH1*QODCX+...)`
- detail: The mixture conductivity is 1/sigma = sum w_i/sigma_i (line 109, correct), so the correct derivative is dln sigma/dln rho = sum_i (w_i/sigma_i) (dln sigma_i/dln rho) / sum_i (w_i/sigma_i). The code instead multiplies the plain weighted sum of the per-species derivatives by mix_log10_cond (the log10 of the conductivity, a number of order 5-20), and for the temperature derivative (line 119) it reuses dlnkap_dlnrho_* instead of dlnkap_dlnt_*. Both are then folded into the final derivatives at lines 125-126 (`-1 - deriv`, `3 - deriv`), and into the radiative/conductive blend derivatives in kap_eval. The value of kappa_cond is right, only its derivatives are wrong, so the effect is on Henyey Jacobian accuracy and on any output that uses dlnkap/dlnT, dlnkap/dlnrho (e.g. the derivative terms in the MLT/convection setup and pulsation columns) wherever conduction matters: degenerate He cores on the RGB, white-dwarf-like cores, HB cores. The commented-out lines 116-117/121-122 show the intended (still not quite right, but derivative-weighted) formula.

## /Applications/YREC/src/kap/conductive/condopacpint.f90:130-158 -- low-density extrapolation branch: sign of log conductivity flipped and extrapolation applied in the wrong direction
- class: physical
- severity: low
- confidence: high
- provenance: inherited -- condopacpint.f lines 77-92 (`CONDL = LOG10(...)` without the minus sign, `OCL = OCL + QODC*(-6.0D0-DL)`)
- detail: In the in-range branch mix_log10_cond = -log10(sum w_i/sigma_i) = log10 sigma_mix, and log kappa = -3.5194 + 3 logT - log rho - log sigma_mix. In the extrapolation branch (log rho < -6 and log R >= 0) line 146 drops the minus sign, so line 151 adds +log sigma_mix instead of subtracting it (kappa_cond proportional to sigma instead of 1/sigma). Line 154 then extrapolates from the edge log rho = -6 to the actual density with delta = (-6 - log rho) > 0, but the correct displacement is (log rho - (-6)) < 0. The derivatives returned are those of pure hydrogen only (lines 132-133; He and O derivatives discarded into unused_deriv). Fires only at log rho < -6 with log T <= 4, where kappa_cond is huge and the blend returns ~kappa_rad, so the practical effect is small; but any diagnostic that prints kappa_cond there is wrong, and the mislabeled derivatives leak into the blend derivatives.

## /Applications/YREC/src/kap/opal95/getopal95.f90:136-183 -- `density_shifted` read uninitialised when the low-X/low-T flag is set (was SAVEd in F77)
- class: numerical
- severity: medium
- confidence: medium
- provenance: modernization -- getopal95.f line 53 `SAVE` (LSHIFT persisted across calls), lines 141-144 set LSHIFT only in the ELSE branch, identical control flow
- detail: `density_shifted` is a plain local; it is assigned only inside the `else` of `if (hydrogen_fraction.lt.0.2d0 .and. opal95_index_t.le.5)` (lines 136-180), yet it is read unconditionally at line 183 to choose between one shared R stencil and four per-row R stencils (with opal95_extrap_hi_row flags from the previous call). In the F77 the routine was under blanket SAVE, so LSHIFT kept the previous call's value; now the value is undefined (in practice stack garbage, frequently whatever the last call left there, but not guaranteed). Fires for X<0.2 at log T below about 4.0 (index_t<=5) -- cool He-rich or metal-rich gas, e.g. the outer layers of a He-burning model reached by the 4D path, or low-X test tables. When garbage is .true., the per-row indices opal95_index_rho(2..4) and extrap flags from an unrelated earlier point are used, giving a wrong but plausible opacity; with -finit-local-zero it is always .false., which is the historically most common outcome. Note also that `low_regime_flag` is set and never used (same in F77: LLOW set, never read), so the "postponed check" the comment at line 134 promises never happens; the 4x4 stencil for X<0.2, low T can straddle the empty (9.999) corner of the table.

## /Applications/YREC/src/kap/alex06/alex06tab.f90:128-132 -- fourth X node abscissa keyed on the Z index instead of the X index
- class: logical
- severity: medium
- confidence: high
- provenance: inherited -- alex06tab.f lines 100-101 `IF(IZ.EQ.NUMZ-3) QR(4) = 1.0D0-ZE`
- detail: The X stencil uses tables index_x..index_x+3; only when index_x = num_x-3 is the fourth table the X=1-Z set (storage 128+iz). The code sets interp_nodes(4) = 1-Z whenever iz = num_z-3 (Z stencil reaching the top, Z>0.06), regardless of index_x, and otherwise uses grid_x(index_x+3). Two mis-labellings result: (a) Z>0.06 with X<=0.8 (index_x<=5): node 4 is a real grid table (X=0.35..0.9) but is labelled 1-Z (e.g. X=0.6, Z=0.08: nodes (0.35,0.5,0.7,0.92) instead of (0.35,0.5,0.7,0.8); Lagrange weights change by ~0.1-0.15 and the interpolant is no longer through the correct points); (b) Z<=0.06 with X>0.8 (index_x = num_x-3): node 4 is the 1-Z table but labelled grid_x(9)=1.0 (small error, ~0.02-0.06 in the abscissa). Case (a) fires for any model using ALEX06 low-T opacities with Z>0.06 in the cool layers (metal-rich stars, or He-burning/AGB envelopes where the 4D path is used); the surface opacity error is a few percent in kappa, enough to shift Teff of cool metal-rich models.

## /Applications/YREC/src/kap/opal95/opal95_fixed_z_table.f90:34-47 -- fixed-Z table rejects 0.08 < Z <= 0.10 although the comment (and the table) allow Z up to 0.10
- class: logical
- severity: low
- confidence: medium
- provenance: inherited -- op95ztab.f lines 30-37 `DO I = 3,NUMZ-1 ... IF(ZZG(I).GE.Z)` then `DESIRED Z > 0.1D0; STOP`
- detail: The loop looks for the first i in 3..num_z-1 with grid_z(i) >= Z; num_z-1 = 12 is Z=0.08, so any Z in (0.08, 0.10] falls through to the "DESIRED Z > 0.1D0" error and returns ierr even though the Z=0.10 tables exist and the later stencil code (lines 98-116) explicitly handles a stencil containing Z=0.10. The 4D path (getopal95) handles the same Z correctly, so this only affects setupopac/rezone calls that build the envelope fixed-Z table for a very metal-rich envelope (Z between 0.08 and 0.10); the run dies with an out-of-range message instead of proceeding. Comment/code mismatch.

## /Applications/YREC/src/kap/laol89/sulaol.f90:58 -- second LAOL table's T grid is logged with the first table's row count
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- sulaol.f line 64 `DO IT=1, NUMT` followed by line 68 `DO IT=1, NT2`
- detail: For the second (two-Z) LAOL table the loop `do it=1, numt ... ot2(it)=log10(ot2(it))` uses the first table's numt, while the immediately following loop and every consumer (gtlaol2's locate on ot2) use opacity_table%nt2. If nt2 > numt the top nt2-numt entries of ot2 stay linear T (not monotonic with the logged rows below, so locate brackets the wrong row); if nt2 < numt, log10 of zero-padded entries is taken (-Inf, harmless if never indexed). No observable effect as long as both LAOL tables ship the same T grid, which is the case for the distributed files; it is a copy-paste latent bug that would fire silently with a differently gridded second table.

## /Applications/YREC/src/kap/kap_lib.f90:197-236 -- HB branch has no else-trap: if none of OPAL95/OPAL92/LAOL is selected, `table_metal_fraction` and the interior opacity are used uninitialised
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- getopac.f lines 81-91 same IF/ELSE IF chain without ELSE; the non-HB branch (getopac.f line 141, kap_lib.f90 line 322) does have the trap
- detail: In the Z>0.1, logT>7 branch the three table options are an if/else-if chain with no final else, whereas the parallel non-HB branch at line 322 stops with "no interior opacity" if no table is selected. With use_pure_z_table set but no high-T table chosen (a possible input combination since the flags are independent), `table_metal_fraction`, `opacity`, `log10_opacity` and the derivatives are undefined and are then blended with the pure-Z result at lines 225-235 (division by `table_metal_fraction - 1`). Only fires on a mis-specified input; the observable is a silent garbage opacity instead of the clean stop the other branch gives.

## /Applications/YREC/src/kap/laol89/gtlaol.f90:47 (also gtlaol2.f90:45, gtpurz.f90:43) -- extrapolation tolerance uses the natural log of tollaol on a log10 density grid
- class: physical
- severity: low
- confidence: medium
- provenance: inherited -- gtlaol.f line 28 `TOLL = LOG(TOLLAOL)`, line 63 `RDL.GT.TR(1)-TOLL`
- detail: The comment says tollaol "permits some extrapolation beyond table edge" and the default is tollaol = 10 (opacity_table_lib), suggesting the intended meaning is a factor of 10 in density, i.e. 1 dex on the log10 rho grid. `log(tollaol)` is ln(10) = 2.303, so linear extrapolation is allowed 2.3 dex beyond either density edge of the LAOL table, more than twice the evident intent. Only matters for LAOL89 runs (rarely used now), where opacities up to 2.3 dex outside the table are silently secant-extrapolated rather than flagged.

## Weak/uncertain observations

- kap_lib.f90:334-351: the atmosphere/interior ramp blends kappa linearly in T but blends the *logarithmic derivatives* linearly with the same weights, omitting the (kappa_i - kappa_a) dw/dlnT term and the kappa-weighting; inherited (getopac.f lines 160-170). Derivative-only error confined to the ramp zone.
- kap_lib.f90:180: HB-branch temperature threshold is logT>7.0 (F77: 7.7); documented as a deliberate JCZ 2021 change, noted only because it widens the reach of finding 1.
- condopacp.f90: hunt initial guesses are plain locals, so on every call the search starts from an undefined guess (F77 SAVEd them); hunt is robust to a bad guess, so only a speed cost, but the value is formally undefined.
- condopacpint.f90:92-93: the in-range test is `logT.le.9.0 .and. logrho.le.9.75` while condopacp treats the exact top edge as out of range (ierr), so a point exactly at logT=9.0 or logrho=9.75 errors instead of interpolating; inherited.
- condopacpint.f90:85-88: when no conductive opacity is available the "derivatives" are set to the opacity value itself (1e10); only consumed if a caller ignores got_conductive_opacity (none found); inherited.
- getopal95.f90:166-169: the high-rho row search tests opal95_fixed_z_opacity for the corner but walks opal95_full_opacity for the row, with an X index from the fixed-Z stencil; consistent only because both are 9.999-padded identically; inherited.
- opal95_surface_table.f90: X node abscissa for table 10 is grid_x(10)=1.0 rather than 1-Z, and for Z exactly 0.04 with X>0.9 the (existing but not built) row 9 is used; small abscissa error; inherited.
- opal92_interp3d.f90 / opal92_interp3d_z2.f90: `rho_search_index` / `dens_index_z2` are intent(inout) starting guesses to findex that are never initialised (F77 SAVEd); findex tolerates any guess; modernization, harmless.
- read_alex94_tables.f90: the file read clobbers opacity_table%alex95_index_t with the file's row counter; guarded by the later search, harmless.
- getalex06.f90 and alex94 (yalo3d) extrapolate_linear branches keep the cubic edge derivative while linearly extrapolating the value; inherited inconsistency, small.
- alex94_surface_table.f90: on a Z change the surface (fixed-X) table is rebuilt at the *current point's* X rather than the envelope X; functionally masked because the interior 4-X path is used when X differs; inherited.
- alex94_interp3d.f90: extrapolate_linear handling exists only in the surface-X branch, not the interior 4-X branch; inherited asymmetry.
- kurucz2.f90 lacks the -3/4.1 pre-check that kurucz.f90 has (documented as intentional); relies on the alternate-return path.
- read_kurucz_tables.f90 stores *linear* kappa in the array named kurucz_log10_opacity until build_kurucz_splines logs it in place; naming only, but a future reader of the array before build would get the wrong quantity.
- laol89 gtlaol/gtlaol2/gtpurz return secant-slope derivatives (dlnkap/dlnrho from neighbouring rows, not the spline derivative); inherited approximation.
- kap_lib.f90 kap_update_surface_tables: jerr from a callee is captured but not propagated in one path (no functional effect found).

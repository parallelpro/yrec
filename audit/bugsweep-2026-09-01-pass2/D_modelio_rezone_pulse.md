# Sweep D (pass 2): starting-model ingest, rezoning, rotation stability setup, stitched pulse model, pulse/profile/history writers, model file I/O

Date: 2026-09-01. Read-only review; no source edited, nothing built or run.

## Files read in full

- src/core/read_starting_model.f90 (1234 lines; ancestor starin.f)
- src/setup/rezone.f90 (999 lines; ancestor hpoint.f -- "points.f" does not exist at 6cd5673)
- src/setup/rotation_stability_setup.f90 (502 lines; ancestor setupv.f)
- src/core/stitched_model.f90 (646 lines; ancestor stitch.f for the envelope stitch; the pulse/profile/seismic column layer is new code)
- src/io/write_fgong_pulse.f90 (122), src/io/write_gyre_pulse.f90 (119), src/io/write_gsm_pulse.f90 (114)
- src/io/profile_output.f90 (164), src/io/history_output.f90 (305), src/io/output_columns.f90 (112), src/io/write_output_headers.f90 (49), src/io/yrec_output.f90 (260)
- src/io/read_yrec7.f90 (300), src/io/read_model2.f90 (137), src/io/read_mod_model.f90 (132), src/io/write_mod_model.f90 (83)
- F77 at 6cd5673 consulted: starin.f, hpoint.f, setupv.f, stitch.f, getyrec7.f, getmodel2.f, putyrec7.f, putmodel2.f, midmod.f, getw.f, coefft.f, wrtmod.f, wrtout.f, wrthead.f, rotgrid.f, getrot.f, setups.f
- Supporting modern files consulted for consumers/producers: core/henyey_coefficients.f90 (pulse_* fills, eps_channels), rotation/mid_timestep_model.f90, rotation/rotation_scratch_lib.f90, core/evolve_step.f90 (rezone call), envelope/envelope_derivs.f90 and envelope/envint_kernel.f90 (gradient slot order), eos/eqstat.f90 + eos/eos_lib.f90 (result slot semantics), numerics/numerics_lib.f90 (osplin signature), state/star_info_lib.f90 (index constants)

## Verified clean (re-derived and found correct)

- rotation_stability_setup.f90 vs setupv.f: line-by-line faithful. First/last interface blocks, the 4-point Lagrangian interface interpolation, delmi minimum, theta_mean, ES/GSF factor loop, LDIFAD block, `transport_zone_begin-1` access (same as IMIN-1 in F77). Call site secular_transport.f90:170 matches the dummy list.
- rezone.f90: envelope-T range check, flag scan and dedupe, overshoot/fine-zone insertion loops, chi/spline point placement, X- and Z-gradient index arithmetic (matches hpoint.f), locate_new_cz_edges, osplin argument order (xval,yval,xtab,ytab,n,k) for the structure arrays (logP/logT/logR/L/logRho/composition), the rotation remap (j_rot, omega, fp, ft, eta2, mean_radius moved in lockstep with the structure arrays), shell-mass/dm setup after transfer, tho/qwrst setup. Call site evolve_step.f90:297 matches the 7-argument dummy list; rezone always reaches the TRANSFER block (only early returns are ierr paths), so finding 1 below fires every step.
- read_starting_model.f90: build_shell_masses, initialize_rotation_state (J and KE sums), update_surface_mixture, YMOD branch (core_cz_top_index0/envelope_cz_bottom_index0/trial_sign_flag are unset there but never read before being reset at line 443), CNO/isotope block, envelope-refit `j` guard. Non-first-call path is a faithful transcription.
- stitched_model.f90 stitching: interior 1..nz copied verbatim; envelope points are env_struct in fitting-point-to-photosphere order with the `radius <= logR(nz)` skip; atmosphere is atmo_struct outermost-first with height accumulated from atmo_delta_depth; gradient slot orders differ between envelope (rad, ad, actual; envelope_derivs.f90:105-107) and atmosphere (rad, actual, ad; envint_kernel.f90:495-497) and the code indexes each correctly. Spotted-Teff correction matches stitch.f:88. The N2 second pass uses centred dlnP/dr and dlnrho/dr.
- FGONG (write_fgong_pulse.f90): var(1..40) mapping checked against the FGONG definition -- var(2)=ln(m/M), var(15)=N2 r/g (=A4), species var(21..25)=he3,c12,c13,n14,o16, var(29)=h2, var(30)=he4, var(31)=li7, var(33)=n15, var(34)=o17, var(35)=o18, surface-to-centre order, ivers=1300 wide layout with `nn iconst ivar ivers` record. eps (var 9) is nuclear-only because eps_total sums components 1..5 only (henyey_coefficients.f90:388-390); eps_grav in var(19). Except var(14) -- see finding 3.
- GYRE ext writer (stitched_model + write_gyre_pulse_ext path): 18-column schema 101 layout, dlnkap/dlnT and dlnkap/dlnrho convention (already-fixed item, not re-examined beyond confirming call), eps_T/eps_rho as eps*dlneps/dln{T,rho}.
- GSM (write_gsm_pulse.f90): version attribute 110, dataset names/units match the GYRE GSM v1.10 spec.
- history_output.f90: species/luminosity index mapping (vals 43-50 = xa(4..11,1); 51-62 = xa(4..15,nz); 23-27 = lum(1..5); 28 he_c; 29 grav; 30 neu; 33-42 neutrino fluxes) verified against star_info_lib indices.
- yrec_output.f90 isochrone he_core_mass = m(h_shell_zone_begin-1) matches wrtout.f:606.
- read_yrec7.f90 vs getyrec7.f: formats 10/30/40/70/100/200/300/500 identical; extended-composition pair loop parity identical to the putyrec7.f writer.
- read_model2.f90: format 40 `10X,1P7E17.9` differs from getmodel2.f `12X` but 10X matches the putmodel2 writer (`'TLUMX',5X`); the F77 12X was the bug (fixed in commit 956a5ed). Not a defect.
- read_mod_model.f90 / write_mod_model.f90: consistent YMOD 1 layout; luminosity_breakdown zeroed on read (cosmetic).

---

## setup/rezone.f90:845 -- Start-of-step eps/esum interpolated against a mass table that was already overwritten with the NEW grid
- class: logical
- severity: medium
- confidence: high
- provenance: inherited -- hpoint.f 650-657 (`DO 1000 J=1,MNEW: HS(J)=HSS(J)`) followed by 670/674 `CALL OSPLIN(HSS,ESUMO,HS,SESUM,NTAB,NTOT)` / `CALL OSPLIN(HSS,EPSO,HS,YTAB,NTAB,NTOT)`; the modern code reproduces it exactly.
- detail: osplin(xval,yval,xtab,ytab,n,k) (numerics_lib.f90:782) interpolates the table (xtab(1:n),ytab(1:n)) onto xval(1:k). At lines 824-831 the TRANSFER loop copies the new grid into `star%log_mass(1:new_num_zones)`. Lines 845-852 then call osplin with xtab = `star%log_mass` and ytab = `star%eps_total` / `spline_y` (which are still on the OLD grid, n = old_point_count), so the abscissa and ordinate of the table belong to different grids. When points were inserted, xtab(1:M) is the new grid truncated to M entries and each old eps value is attributed to the wrong mass; when points were deleted, xtab(MNEW+1:M) still holds stale old-grid masses that are smaller than xtab(MNEW), so the spline table is non-monotonic and the fitted slopes near the surface are garbage. rot_scr%old_esum / old_eps are consumed ONLY by rotation/mid_timestep_model.f90:165,174 (`esumm = old_esum + f*(eps_total - old_esum)`), which feeds rotation_stability_setup (Eddington-Sweet velocity through epsilm) for every rotating run and every overshoot+extended-composition run, on every timestep (rezone has no early "nothing changed" return). Fix: perform the two osplin calls before the TRANSFER loop, or pass a saved copy of the old mass grid as xtab.

## core/stitched_model.f90:448 -- Profile column 52 'mu' is the specific gas constant R/mu, not mu
- class: physical
- severity: medium
- confidence: high
- provenance: new-code (the pulse_* arrays and profile column layer have no F77 ancestor; the source slot is henyey_coefficients.f90:464 `pulse_mean_molecular_weight(im) = eos_res(i_gas_constant)`, which is the F77 coefft.f:419 `PRMU = RMU`)
- detail: eos_res(i_gas_constant) is eqstat's `specific_gas_constant = gas_constant*(1/mu_ion + 1/mu_e)` (eqstat.f90:479,547,575), i.e. R/mu in erg/g/K (~1.4e8 for solar composition). The array name `pulse_mean_molecular_weight` and the profile_output.f90:99 column name 'mu' both claim it is the dimensionless mean molecular weight (~0.6). Anyone reading the profile column gets a number 8 orders of magnitude off, and the seismic-column builder in this same file trusts the name (finding 4). Observable in every profile file that selects column 52.

## core/stitched_model.f90:449-455, 512-514 -- 'mu_e_inv' profile column and FGONG var(14) are mu_e, not 1/mu_e (inverted)
- class: physical
- severity: medium
- confidence: high
- provenance: new-code (henyey_coefficients.f90:465-466 stores `pulse_electron_mean_molecular_weight = eos_res(i_mu_e_inv)`; eqstat.f90:477,545 `electron_mean_weight_inverse = X*0.9921 + (1-X)*0.4995` is 1/mu_e; F77 coefft.f PEMU=EMU was used by wrtmod.f:216 as an inverse weight, `PELPF = CGAS*T*rho*PEMU`)
- detail: `pulse_electron_mean_molecular_weight` already holds 1/mu_e (~0.85 for X=0.7), but both the profile case (53) and the FGONG builder take `1.0d0/pulse_electron_mean_molecular_weight`, so column 53 'mu_e_inv' and FGONG var(14) (write_fgong_pulse.f90:92; the FGONG definition of var 14 is 1/mu_e) contain mu_e (~1.17) instead. Any FGONG consumer that uses var(14) (e.g. to reconstruct electron pressure or the Ledoux/thermohaline terms) gets the reciprocal. Silent -- values are of order unity either way.

## core/stitched_model.f90:252,281-282 -- Ledoux gradient (profile col 56 'gradL') has the mu-gradient term with the WRONG SIGN, plus a spurious spike at the interior/envelope junction
- class: physical
- severity: medium
- confidence: high
- provenance: new-code (compute_seismic_columns has no F77 ancestor)
- detail: The comment (line 229) and formula intend `gradL = grada + (1/delta) dln(mu)/dlnP`. lnmu is taken from `ext_profile_value(52, j)`, which by finding 2 is R/mu, so dln(R/mu)/dlnP = -dln(mu)/dlnP and the composition term enters with the opposite sign; in a region of inward-increasing mu (H-burning shell, He core boundary) gradL is reported BELOW grada instead of above it, i.e. the column would classify a Ledoux-stable region as more unstable. Independently, ext_profile_value(52,·) returns 0 for envelope and atmosphere points (the `k > nz` branches return 0 for pulse-only columns), so lnmu = log(1e-30) = -69 there, and the centred difference at j = nz and nz+1 straddles the jump from ~+18.8 to -69, producing an O(1e2/dlnP) spike in gradL at the fitting point. Observable in every profile that selects gradL.

## core/read_starting_model.f90:616-627,636 -- Core extension shifts every structure array except logRho
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- starin.f 390-401 shifts HS,HR,HL,HP,HT,LC,HCOMP,OMEGA but not HD; 409 `FACD = HP(M1)-HD(M1)-HT(M1)` mixes the shifted and unshifted arrays exactly as here.
- detail: The shift loop moves log_mass/logR/L/logP/logT/convective_flag/xa/omega up by num_core_shells_added but leaves `star%logRho` in place. `density_estimate_offset` (line 636) therefore combines P and T of original shell 1 with rho of original shell (MCORE+1), and the new core shells get densities estimated with that mismatched offset; shells first_original_shell..nz remain misaligned by MCORE shells until the first EOS pass recomputes rho. Consequence is confined to the first Henyey iteration's starting guess and any diagnostic printed before it (larger MCORE and steeper central rho gradients make the first-iteration guess worse); the run recovers, so severity low.

## core/read_starting_model.f90:604-613 -- JSON overflow in core extension writes "RUN TERMINATED" but does not terminate
- class: logical
- severity: medium
- confidence: high
- provenance: inherited -- starin.f 378-387 has the same three WRITEs and no STOP.
- detail: If `nz + num_core_shells_added > json` the code prints "STARIN: ***** RUN TERMINATED *****" and then falls straight into the shift loop `do i = nz,1,-1 ... star%log_mass(i+num_core_shells_added) = ...`, writing beyond the json bound of every structure array. The message/behaviour mismatch is itself reportable. The routine has an ierr argument path available (used elsewhere in this file); this branch should set ierr and return.

## core/read_starting_model.f90:956-961 -- Envelope refit interpolates X and Z with the sign reversed relative to rho/P/R/T
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- starin.f 668-669 `HCOMP(1,J)=HCOMP(1,M)+F*(HCOMP(1,M)-ENVX(K))`.
- detail: The four structure variables use `val(nz) + f*(env - val(nz))`, but X and Z use `xa(nz) + f*(xa(nz) - env)`, i.e. they extrapolate away from the envelope value. Harmless only because the envelope integrator normally carries the surface X,Z so the difference is zero; whenever env_hydrogen_fraction differs from xa(i_h1,nz) (e.g. a starting model whose surface composition was changed by the run's update_surface_mixture before the refit) the refitted shells get composition moving in the wrong direction.

## core/stitched_model.f90:414-419 with io/profile_output.f90:66-70 -- Profile columns 19-23 'eps_ppI..eps_he3' are fractions, and 'eps_he3' is triple-alpha
- class: interface
- severity: low
- confidence: high
- provenance: new-code (column naming); the storage as fractions is inherited (coefft.f stores SEG(1..5)=component/SESUM)
- detail: henyey_coefficients.f90:399-401 stores eps_channels(1:5) = component/eps_total (dimensionless fractions) while eps_channels(6,7) (neu, grav) are absolute rates. The profile names 'eps_ppI' etc. sit beside 'eps_nuc' and suggest erg/g/s. Further, energy_gen_component(5) is triple_alpha_gen and star_info_lib.f90:64 correctly calls the matching luminosity slot i_lum_3alpha, but the eps slot is named i_eps_he3 and the profile column 'eps_he3' -- a He-burning model's profile labels its 3-alpha fraction as an He3 channel. Documentation/interface bug; numbers are right, labels are wrong.

---

## Weak/uncertain observations

- read_starting_model.f90:422-425: `lexcp0` is printed in format 1040 on the mixing-length mismatch error path; on the first call it has not been read from the model header (uninitialised print, cosmetic).
- read_starting_model.f90:330 and 690: on non-first calls rescale_model can be invoked twice for the same run (inherited oddity from starin.f; second call is a no-op if the first converged).
- read_starting_model.f90:632: `actual_gradient = (logT(2)-logT(1))/(logP(2)-logP(1))` is evaluated after the shift; indices 1..MCORE still hold the stale originals, so for num_core_shells_added == 1 entries 1 and 2 are both original shell 1 and the quotient is 0/0 = NaN, which then propagates into the extrapolated core logT. Inherited from starin.f 406; needs MCORE == 1, which the chi_grid_scale(2) rounding rarely produces.
- rezone.f90:542-549,597-603: X- and Z-gradient point insertion has no json bound check (inherited hpoint.f 393-400); the flag-point overflow guard does not cover it.
- rezone.f90:103,258: pmax1..pmax5 are read from chi_grid_scale(11..15) and never used (dead controls, inherited).
- rezone.f90 flag dedupe: after a deletion pass the duplicate check is not re-run, so a later insertion can coincide with a kept flag (inherited; effect is one redundant point).
- rotation_stability_setup.f90:479: ht_temp_scale2 assigned, never read (same as HTP2 in setupv.f).
- write_gyre_pulse.f90 (interval writer): interior-only model with M = m(nz) and R at the fitting point, documented as byte-pinned historical; argument `log_luminosity` actually receives linear luminosity_lsun (naming only, arithmetic consistent).
- write_fgong_pulse.f90:5-6 vs 38,66-67: header comment says ivers 300 / 1P5E16.9; code writes ivers 1300 / E26.18E3 (comment lies; code matches MESA's wide layout).
- write_fgong_pulse.f90:49-50: glob(4)/glob(5) use surface Z/X of the outermost point; FGONG defines them as initial Z/X -- differs after diffusion.
- write_fgong_pulse.f90: no explicit r=0 centre point is emitted (innermost point is shell 1 at finite r); most FGONG readers tolerate this.
- output_columns.f90:85-112: append_column does not check nsel < max_cols(128) before `sel(nsel)=...`; would need >128 user-selected columns to fire.
- stitched_model.f90: for envelope/atmosphere points all pulse-only columns (cp, delta, kap derivatives, eps derivatives, mu, mu_e) are 0, so profile/ext-GYRE files carry zeros for these above the fitting point rather than a "not available" marker; the N2 and gradL centred differences (finding 4) are the only columns that mix them with interior values.

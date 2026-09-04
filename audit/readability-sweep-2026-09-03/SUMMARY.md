# Readability / maintainability sweep, 2026-09-03 -- consolidated summary

Seven read-only reviews, one per domain, each reading every file in
its area in full (about 230 files). Per-area reports are alongside
this file (core, burn, eos, kapatm, rotation, mixwind, io);
INSTRUCTIONS.md is the shared brief. This sweep was about clarity and
maintainability, not bugs -- the two prior bug sweeps
(bugsweep-2026-08-31, bugsweep-2026-09-01-pass2) and ROADMAP sections
10/11 were treated as known and not re-reported. Bugs that surfaced
anyway are listed in section 1.

Every proposal is tagged (A) byte-safe or (B) changes numbers. Nearly
everything here is class A. The few B items are named explicitly so
they can be batched into one reseed.

The eleven bug claims in 1.1 were checked against the source by hand
before this summary was written; the rest of the reports are the
reviewers' words.

## 1. Bugs found on the way (not in ROADMAP 10/11)

### 1.1 Verified against the source

| # | Where | What | Reach | Fix class |
|---|-------|------|-------|-----------|
| 1 | core/read_starting_model.f90:829 | `log10_radius = 0.5*(log_L + star%solar_luminosity_cgs - ...)` uses the LINEAR Lsun (3.8e33) where the F77 (`CLSUNL`) and every sibling site (rebuild_envelope:129, surfbc:185, stitched_model:147) use `log10_solar_luminosity`. Modernization regression from 6e31b68. | Only the "requested envelope mass shallower than the model's" branch of rescale_and_refit_envelope; no Stage-0 deck runs it. | A for every pinned case; B on that branch. |
| 2 | core/envint_kernel.f90:801-802 | Envelope inversion swaps `env_luminosity` one way only (missing `env_luminosity(i2) = swap_temp`); the other 19 triplets swap both ways. | Latent: env_luminosity is a constant in every entry today. | A. |
| 3 | eos/opal/rhoofp01.f90:73 | Table-priming call passes `rad_flag` with `ideriv_dbg=1`; radsub01's no-radiation block then divides by an uninterpolated cv (0/0). Same defect Batch 3 fixed in rhoofp06:87 (pass 0). | Only with use_opal2001_eos, which no test enables; a trap under -ffpe-trap=invalid. | A (output discarded). |
| 4 | core/check_solar_calibration.f90:72 | `dlum_dalpha = 0.0139d0 ! empirical result: +0.139` -- value and comment differ by 10x. Inherited from chkcal.f:43. | Solar-calibration Newton step only. | Needs the author: which is right? B either way. |
| 5 | kap/opal95/ll95tbl.f90:105-140,184 | Every body-row `read` has `end=9999`; label 9999 is `continue; return` with ierr still 0 and opal95_fixed_z_table never called. A truncated table returns success. | Error path only. | A. |
| 6 | io/equal_to_model.f90:81-84 | `k0 = min(max(j-2,1), num_equal_points-3)`; the "JVS fix for NPT=3" catches k0=0 but not k0=-1 (num_equal_points=2, reachable when two model intervals are exactly equal). Out-of-bounds read. | Rare; reachable. | A for n>=4 (clamp only). |
| 7 | mixing/temperature_gradients.f90:260 | `convective_velocity = 1.0d0-11` is -10, not 1e-11. Branch is unreachable (same sign test already taken at :202). | Dead branch. | A (delete branch). |
| 8 | core/henyey_iterate.f90:160 | `if (start_new_triangle.or.reset_triangle .and.iteration_level.eq.2)` parses as A .or. (B .and. C). Matches the F77 so the numbers are "right"; the intent is unclear. | Live. | A to parenthesize as-is; B if the intent was (A .or. B) .and. C. |
| 9 | mixing/mix.f90:466, rotmix.f90:220 | Settling substep ceiling uses `mod(max_settling_dt, timestep)` -- arguments reversed. Works by accident except when max_settling_dt is an exact multiple of timestep (mix: divide by zero; rotmix: one extra substep). | Live, rare. | B. |
| 10 | mixing/temperature_gradients.f90:126 | `star%alfmlt/phmlt/cmxmlt` are set to 0 here and nowhere else; henyey_coefficients:474 copies them to the output arrays valfmlt/vphmlt/vcmxmlt, which are therefore always zero. Either dead output or a lost assignment in the MLT block. | Output only. | A (delete) or needs the author. |
| 11 | setup/rescale_model.f90:96-103 | Loop variables swapped: `species_idx` runs over zones, `zone_idx` over species 5..11. Correct numerically, reads as wrong. | Cosmetic. | A. |

### 1.2 Reported, not independently checked (see the area reports)

- core: check_solar_calibration:113-117 reads uninitialised `log_zx_mismatch`; run_yrec:144-149 comment vs `cycle`; henyey_coefficients:245 vs 379 `.gt.`/`.lt.` equality edge; shell_physics:77 `saha_state` never zeroed; envelope_derivs:80 integrand writes star%iovim.
- burn/net: `en` assigned never read (burn_lib 1122,1157,1756); net_lib:2171 `o16_gamma_frac` unused; stale "LOG TO BASE 10" comments on absolute derivatives.
- eos: eqbound/eqbound01 compute `t_fraction` and never use it (06 uses it); mhd/rabu:31 bounds check after the read; `eos_get_gamma1` has no production caller; dead writes between `return` and `end if` in esac01/esac06.
- kap/atm: read_kurucz_tables silently drops malformed rows; readalex06 validates 16 of 19 logR columns; opal92_table_prep addresses one row with two strides; opal95_fixed_z_table "OMIT TABLE 10" comment vs unconditional code; alex94_interp3d honours `extrapolate_linear` in one branch only.
- rotation: viscos Endal-Sofia block computed and never read; am_advection_diffusion_coeffs residual_check/damping_factor dead; theta_term_n/p permanently zero; secular_transport:505 accepts an unconverged iterate silently; implicit_diffusion_coeffs half-cell offset (masked by a zeroed coefficient); zone_moments_of_inertia header vs code.
- mixing/wind/setup: `rot_scr%difad_shear_coeff2` only ever 0 so the diffusive-shear terms in circulation_velocities are identically zero; kawaler_wind vs wind_spindown disagree on the saturation gate; mdot remixes species 1..11 only under extended composition and stores grams in `accreted_mass_fraction`; massloss entropy path dead (`accretion_efficiency` = 1 constant); find_convection_zones:78 writes `convective_flag(num_zones+1)`; semiconvection:218 writes the EOS density back into star%logRho.
- io/numerics: interpret_kind_cards sets the EOS mixture to Z=0 transiently on every kind card; tridia:1318 no-op rescale; intpt silent return with outputs unset; model_to_equal dead derivative block and mixed index spaces; yrec_output 'GSM' vs 'gsm'; read_mod_model first read unguarded; write_gsm_pulse ignores herr.

## 2. Cross-domain themes

The same six patterns account for most of the reading cost in every area.

1. **Clones by suffix.** Second-Z / second-vintage copies of whole files: kap (kurucz2, gtlaol2, opal92_*_z2, kcsurfp + mirrored halves in 8 setup routines, ~1100 lines), eos (three OPAL vintages x 9 files, three 60-line blend blocks in eqstat2, 11-way MHD ladder), burn (13-reaction kernel in engeb and net_lib rates; liburn/liburn2 share ~250 lines), core (envelope refit written twice: read_starting_model 813-1034 vs rebuild_envelope 113-323), rotation (am_transport_grid vs composition_grid; microdiff vs microdiff_etm species blocks), mixing (He-settling substep in mix vs rotmix), numerics (cspline/kspline, splinc/splinj, splint/splintd2, interp/intrp2). Proposed everywhere: one implementation over a derived type or an argument, keeping operand order.
2. **Hand-rolled stencils and small formulas repeated.** 4-point Lagrange search (13 copies in kap, 3 in circulation_velocities, 36 unrolled in rotation_stability_setup, 7 in eqscvg, 3 in model_to_equal); log10 R from (L, Teff) (6 sites in core, 6 in wind); spot-adjusted Teff (5); shell masses from log mass (3); CZ-extent search (4 in seculr, 3 in mixing); disk-locking predicate (3); dt-factor clamp (3). All byte-safe if the helper keeps each caller's operand order and literal kinds.
3. **Packed arrays addressed by bare integers.** `mhd_output(25)`, OPAL `eos_output` slots, `reaction_rate_1..13`, `atime(1..14)`, `chi_grid_scale(1..12)`, `rescale_params(1..4)`, profile columns in stitched_model, species columns 1..15 in mixing/rotation, `laol_work_array` permutation. Proposed: `integer, parameter` index blocks in the owning lib module.
4. **Hidden channels through `star%` / scratch.** rezone uses `logP_start/logT_start/logRho_start/old_shell_mass` as spline scratch; rotmix smuggles `mix_scr%delm` through `star%gradT`; temperature_gradients reads `star%iovim` left by the last Henyey pass; atm_get configured by save/override/restore of `star%job%env_step_*` at three sites; kernel integrands set `star%iovim`; loaders write foreign state (rdlaol -> yale_eos_lib debye_huckel_z, setupopac -> star%use_two_z_tables, print_allard_tables mutates its table). Proposed: explicit arguments / local buffers.
5. **F77 control-flow residue.** Alternate returns (`*999`, `return 1`) in eos and kap; `-999` and `9.999` sentinels; `12345678` init flags; `continue` before `return`; goto-emulation blocks with unreachable tails; `if (i > n)` after loop completion; sign-encoded flags (`star%dt = -abs(dt)`, `fcorr0` sign as mode, negated stop thresholds); the intent(out) flag relay in henyey (four constants set at the bottom and threaded up to unread SAVEd driver locals); dead kenv/katm counters through six files.
6. **Stale headers.** Roughly 100 comments cite files, COMMON blocks or SAVE semantics that no longer exist (`core/read_input.f90`, `const_lib`, `setup/remap.f90`, "common/a/", "SAVE preserved", "stub stops", a botched-sed `star%*, star%*` line in observables_lib). Also headers that describe the opposite of the code (net_lib merge direction, gravitational_settling "metals later", zone_moments_of_inertia, massloss "bug not fixed").

Dead code worth a single purge commit (all A): burn_settle_mix + burn_mix_extrapolated and the BUR-ST branches in evolve_angular_momentum (switch hardwired false); mid_timestep_model 197-352; viscos Endal-Sofia; am_advection residual/damping; banded_solver 114-152; wcz; observables_lib locate_core_cz; read_controls 960-1138 DATA archaeology; splinj/splnr (no callers); homogenize_convection_zones 118-161; model_to_equal derivative block. Roughly 1500 lines.

## 3. Suggested batches (in order)

Each batch is one commit series on yrec-modern, Tier 1 gate per commit, with a byte-compare of every pin at the end of the batch. Class A batches need no reseed; the one B batch gets one reseed.

**R1 -- Comments and dead code (A, S-M).** Theme 6 plus the dead-code purge above, plus the eleven cosmetic/naming bugs in 1.1 that do not change numbers (#2, #5, #7, #11) and the unused-local lists in every report. No numerics touched. Add a gate1 warning that greps comments for `src/...f90` / `common/` paths that do not resolve (io recommendation 5).

**R2 -- Named indices and constants (A, M).** Theme 3: index parameter blocks for mhd_output, OPAL eos_output, reaction rates (`star%reaction_rate_1..13` -> array), atime / chi_grid_scale / rescale_params, profile columns, species columns; kap_table_dims module (one dimension block per family -- the ll95tbl 126/130 class of drift); named sentinels. Pure renames and same-value parameters.

**R3 -- Small shared helpers (A, M).** Theme 2: `log10_radius_from_l_teff`, `spot_adjusted_log_teff`, `shell_masses_from_log_mass`, `stencil4_locate`, `lagrange4`, `cz_extent`, `disk_locked`, `clamp_dt_factor`, envelope reversal by array section (fixes #2 by construction). Each helper must preserve the caller's operand order and literal kinds; byte-compare after each.

**R4 -- Explicit data flow (A, M).** Theme 4 and theme 5's flag relay: henyey flags to locals and delete kenv/katm; atm_get optional `envint_config` argument replacing the save/override/restore blocks and dummy arrays; rezone local spline buffers; rotmix gradient argument; temperature_gradients `zone_index` argument; microdiff `*_bl` local unit copies; loaders return what they learn. Signature changes at ~30 call sites, arithmetic untouched.

**R5 -- De-duplication by derived type (A, L).** Theme 1, one domain per commit series: kap table-set types (kurucz/opal92/laol/surface_p) deleting the `2`/`_z2`/`c` files; eos `opal_eos_vintage` type + `blend_opal_result`; MHD composition-indexed state; one envelope refit (read_starting_model calls rebuild_envelope); numerics spline family; am_transport_grid/composition_grid split; burn kernel (engeb vs rates) -- note the engeb rate block 1143-1149 has a measured bit-drift history, so that one extraction is B and belongs in R6.

**R6 -- Numbers change, one reseed (B, S).** Batch together: #1 (log10 Lsun in read_starting_model, no pinned case moves), #3 (rhoofp01 priming flag), #4 dlum_dalpha (after the author decides), #9 mod argument order, kawaler/wind_spindown saturation gate, `pow(1d1,x)` -> `exp10`, single-precision literal d0 suffixes (map_user_inputs 9.4E0, boole scalex/scaley, saha_eos DATA), tridia no-op rescale, engeb kernel extraction. Measure drift per the Batch 3 ritual, record it in the commit message.

Open decisions for the author before R6: dlum_dalpha 0.0139 vs 0.139; whether alfmlt/phmlt/cmxmlt output was ever meant to be live; whether BS extrapolation (burn_settle_mix) is ever coming back or can be deleted in R1; henyey_iterate:160 intended precedence.

## 4. Files reported fine as they are

core: henyey_eliminate, henyey_solve, envint_lib, point_scratch_lib, star_setup, yrec_capi, yrec_reset, neutrino_flux_table, stop_conditions, monte_carlo, shell_physics. eos: eos_mixture_lib, mu, yale_eos_lib, mhd_eos_lib, rtab/mhdtbl/mhdpx/mhdpx1, gmass01/06, t6rinteos01, radsub01, fully_ionized_eos, saha_eos, test_eos. kap/atm: atmstruct_lib, envstruct_lib, ttau_lib, altabinit, surfopac, conductive_table_lib, rdzlaol, zsulaol, alex94_fixed_z_table, alex94_surface_table, getalex06, opal95_interp{2,3,4}d, test_kap. rotation: am_convective_regions, am_diffusion_coeffs, composition_diffusion_coeffs, equal_grid_to_model, check_composition, zahn_coupling_factor, thoul_diffusion, implicit_diffusion_coeffs, equipotential_integrand, shell_inertia_integral, omega_from_j. mixing/wind/util/setup: version, xrng4, zero, timestep_limit_structure, locate_shell_boundaries, matt_structure_factor, wind_spindown, wcz (dead), overshoot_boundaries, setup_solar_calibration, setup_star_calibration. io/numerics/state: check_controls, history_output, luout_lib, run_log, write_mod_model, write_output_headers, profile_output, read_mod_model, math_lib, controls_sync_lib (generated), phys_const_lib, intpar_lib, star_info_lib. burn/net: see burn.md.

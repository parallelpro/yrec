# eos/ review (2026-10-29, read-only evaluation; paths under src/)

Scores vs rotation/: Readability 3, Organization 3, Maintainability 2, Extensibility 2, Docs accuracy 3.

In good shape: the single public facade (enforced), one settings record, one result record with units,
the derivative contract and its finite-difference check, Saha and fully ionized solvers as named steps,
2001/2006 OPAL merged in one module, named columns for SCV and MHD.
Main obstacle: implicit contracts and hidden state -- constants set by core at run time, in_atmosphere
leaving result fields stale, the Saha warm start in module variables (saha_state is only an iteration
counter), an SCV table not tied to its mixture, the Debye-Hueckel mixture depending on the opacity
choice -- and a golden-file-only test, which is how the MHD reader lost all its tables (newunit move).

## Possible bugs
1. MHD tables are never read -- eos/mhd/mhd_eos.f90:131-146,223,240,257. VERIFIED by the integrator.
   newunit= returns negative units, but mhdst reads a table only `if (unit > 0)`. With use_mhd_eos
   nothing is read and ierr stays 0; lines 256/258 index zams_lower_log10t/zams_upper_log10t at 0
   (num_t = 0); compute_mhd_eos then interpolates zero tables (1/chi_rho = 1/0). Introduced by c00a9a20.
   The opens have no iostat, so a missing file aborts, although eos_init's header (eos_lib.f90:94) says
   ierr = 1. Fix: drop the unit guards (all 8 files are opened unconditionally), add iostat. The MHD
   tables are not shipped (YREC_MHD_TABLES env for test_eos), so add a reader test on a synthetic
   unformatted file.
2. Debye-Hueckel metals depend on the opacity choice -- eos/eos_lib.f90:231-235 with kap/kap_lib.f90:537-541.
   debye_huckel_z(4:18) is filled only when the LAOL89 opacity tables are read; with DH on and other
   opacities, almix (eos/yale/yale_eos.f90:731-734) counts only C, N, O. Fix: eos builds the 18-element
   mixture from its own mixture, or the star layer always passes it. (Results change for DH decks
   without LAOL89, if any -- measure.)
3. in_atmosphere leaves first derivatives, cp and grad_ad stale on the Saha path -- eos/yale/yale_eos.f90:208.
   Below cutoff-0.5 and outside SCV/OPAL coverage, dlnrho_dlnt, dlnrho_dlnp, cp, grad_ad keep the caller's
   input values; in the ramp they are blended in (eos/blended_eos.f90:356-359). README sec.3 and
   eos/eos_results.f90:15 say every field is overwritten. Callers that pass .true. read them:
   convection/semiconvection.f90:146,179,186 (grad_ad at CZ edges), core/envelope/envelope_integration_kernel.f90:282.
   Fix: always compute first derivatives (class B; measure).
4. eos_get_gamma1 uses derivatives it asks not to compute -- eos/eos_lib.f90:188-189 passes
   in_atmosphere=.true. then builds gamma1 from dlnrho_dlnp, dlnrho_dlnt, cp of a fresh record (1/0 on
   the Saha path). The test comment (eos/test/test_eos.f90:249-255) blames lderiv=.false. instead.
   Passing .false. leaves the pinned points unchanged.
5. rhoofp06's bracketing fix not ported -- eos/opal/opal_eos.f90:366-421 (rhoofp01) and
   eos/opal/opal95_eos.f90:183-239 (rhoofp): shrink density once by x0.2, regula falsi from a possibly
   unbracketed pair, give up after 10 tries, return -999 silently. Batch 2's bracket loop, cap 30 and
   run-log message went into rhoofp06 only (opal_eos.f90:1511-1573). Also the priming twin at
   opal_eos.f90:308-309 (known, ROADMAP sec.12). Class B.
6. use_diffusion_z changes the hydrogen constant -- eos/blended_eos.f90:255-256 vs eos/eos_mixture.f90:60,70:
   the diffusion-of-Z branch divides by 1.008 (and 4.004); the default multiplies by 0.9921 (0.24975);
   1/1.008 = 0.992063, a 4e-5 shift of the H term of 1/mu_ion even for uniform Z. One constant (B).
7. The MHD "mean weights" check can never fail -- mhd_eos.f90:93-98 with eos_mixture.f90:103-107
   (compute_mean_weights back-solves 1/mu_e from beta P/(rho T), so its R equals P_gas/(rho T)
   identically); the "cross-check" header (eos_mixture.f90:80-82) misleads; compute_opal_eos discards the
   same output (opal_eos.f90:162-164).
8. MHD reports stale results as success -- mhd_eos.f90:569-572,783-785: out-of-range T or P above a table
   row returns ierr = 0 with the previous point's mhd_output.
9. Negative beta only warns -- eos/blended_eos.f90:317-323 (continues with negative beta_inverse;
   divides by zero at beta = 0). Low priority.
Known (ROADMAP), not re-reported: close(2) at opal_eos.f90:1011, stale var_idx column in log10_ne_grid,
SCV idp = 0, silent non-convergence in the Yale solvers, the 1995 t6rinterp dix.

## Findings (ranked)
1. Run-time constants -- state/physical_constants.f90:26-27 (gas_constant, radiation_constant_over_3, cmkh
   are variables) set by core/run/run_setup.f90:93-116; eos is wrong unless core setup ran, which is why
   test_eos.f90:47-49,203 boots a star. Make them parameters (gas_constant is a literal; the others need
   the compile-time-rounding check of Code style rule 3). A, S.
2. Per-vintage OPAL layout -- opal/opal_eos.f90:494-504,550-579,656-661,716-761,783-823,852-903,1019-1145,
   compute_opal_eos 67-68 and 90-151, blended_eos.f90:211-216,448-476: ~17 `if (vintage == iv_opal01)`
   branches. Replace by a parameter array indexed by vintage (nr, nt, nvar, x grid, slot numbers (0 when
   absent), message prefixes, (T,rho) window, radiation constants). Lines 852-872 can collapse today
   (i_opal01_p == i_opal_p, i_opal01_e == i_opal_e). A, M.
3. Merge radsub01 and radsub06 -- opal_eos.f90:1165-1249 and 1597-1682: with rad_flag = 1 (the only value
   passed, 876-883) they do the same operations in the same order; they differ in slot numbers,
   single/double literals, the dE/drho correction (2006) and the gamma3 difference (2001). One routine
   driven by finding 2's table with the literals as arguments. A, S.
4. blend_eos clean-up (A, S-M): fix the misindented nest at 336-419; compute_saha_eos is called identically
   in both branches (340-352) -- once, then SCV if selected; the fully ionized setup is duplicated
   (328-335, 371-378) -- one helper; derive saha_ramp_scale (2) from saha_ramp_width (0.5) (184); stop
   reusing beta14 as Prad/P before reassigning (314-324); num_species duplicates n_mix_species (181,289,298);
   the four _dt/_dp differences that run on undefined values when want_derivatives is false (389-392)
   go inside the if; stale comments ("defined in STARIN" 242, "or it has failed" 348, the comment above
   opal_selected 443-447), ierr declared after locals (232).
5. evaluate_shifted -- blended_eos.f90:119-153: lines 137-146 copy ten fields that the final blend_eos
   overwrites before reading; first/unused round-trip unread. Keep only the density chain and the four
   first derivatives. A (B only in the in_atmosphere corner), S.
6. saha_state -- yale/yale_eos.f90:44,315, eos_lib.f90:59, README:44, blended_eos.f90:50-51,75: an iteration
   counter nothing reads (s%central_saha_state only keeps it); the real state is saha_e_guess/x1/x2 (module
   variables), which couples callers incl. the shifted evaluations. Rename and fix docs (A); a
   caller-owned guess record is B, M.
7. eos_init's 14 positional path arguments -- eos_lib.f90:95-99, core/run/run_setup.f90:139-143,
   test_eos.f90:283,394. Use an eos_table_paths_t record built next to equation_of_state_settings. A, S.
8. OPAL vintage selection -- eos_configuration.f90:16-17: three logicals, "at most one" unenforced. One
   opal_vintage integer mapped from the controls (setup/microphysics_settings.f90) with an error on
   conflict. A, S.
9. Debye-Hueckel and species slots (A, S): name the bare indices in eos_set_debye_huckel_composition
   (eos_lib.f90:250-255) with i_h1, i_he4...; species_ion_fraction(7) (yale_eos.f90:205) -> ix_h; state or
   assert that ionization_temp's order (34-36) matches composition_layout's ix_* order; pass the DH
   composition as an argument instead of module variables (yale_eos_tables.f90:29-31).
10. SCV coupling -- scv/scv_eos.f90:907-912,999-1019: build_scv_envelope_table uses the lookup's warm-start
   caches idtt/idp as loop scratch (masks the known idp = 0 bug); compute_scv_envelope_eos (96-97) tests the
   composition against eos_mix, not the X, Z the table was built for, so eos_set_mixture and
   eos_update_envelope_tables must be called as a pair (core/model/model_rescaling.f90:158 and
   io/run_cards.f90:122 call only the first). Record the table's composition (A); local scratch (B). M.
11. Dead state and identity indirection (A, S): opal_eos_var_order and eos_index_inverse are identity maps
   (opal_eos_tables.f90:110, 78 references); write-only members amu_grid, log10_ne_grid,
   species_fraction_header, moles_per_gram_table, log10_r_value, mean_molecular_weight_header, xxh (keep
   the 1995 member layout); yale_eos.f90:244-246; opal_eos.f90:743,749,1200-1201; eos_get_gamma1's input
   copies (eos_lib.f90:156-172) with their stale comment.
12. No message policy -- kernels print to run_log_unit, to * (opal_eos_tables.f90:401,416), or both
   (mhd_eos.f90:99-106,599-602); warnings next to fatal errors. Use rotation_messages' pattern
   (an eos_messages module). A, S-M.
13. Test gaps -- test/test_eos.f90: eos_get never called with in_atmosphere=.true. (224,296,325,420-422,547);
   no eos_update_envelope_tables after a composition change, no eos_set_debye_huckel_z, no edge-ramp
   continuity, no MHD reader test; no ok/FAIL anywhere. Add physics-limit checks: ideal gas
   (rho = mu P_gas/(R T), grad_ad -> 0.4, cp -> 5R/2mu), radiation-dominated grad_ad -> 1/4, continuity
   of rho and grad_ad at cutoff-0.5, at the cutoff and at the OPAL edge ramp. A, M.
14. MHD readability -- table_selector encoding (-ic / +ic / 3+ic) gives one flag three meanings
   (mhd_eos.f90:649-705): pass family and index; arrays for the eight paths and units (120-158); name cnvs (53).
15. Stale comments: test_eos.f90:12-31 (unit numbers, iopale, env_comp) and 313; opal_eos.f90:203-207,
   1254-1258 (compute_opal2001/2006_eos) and 1012 ("opened by read_controls"); opal95_eos.f90:934;
   mhd_eos.f90:117-119 ("star%ctrl" unit numbers); eos_lib.f90:151 (eos_eval); opal_eos_tables.f90:421
   ("linear interpolation" but takes next row's edge); eos_mixture.f90:59 ("spare slot" is He's Z/A, used at
   blended_eos.f90:329); tools/check_boundaries.py:58-61 (eqstat, mu, eos_eval).

## Structural recommendations
1. eos owns its constants and gets a star-free test: constants as parameters; a test that calls eos_init
   directly with limit, continuity and in_atmosphere checks as ok/FAIL, mutation-checked, plus an MHD reader
   test on a synthetic unformatted file; drop the star boot. M.
2. A table-source descriptor in blend_eos: each source supplies selected, window(T,rho), evaluate (to an
   eos_result plus a covers flag) and edge weight; the three OPAL vintages become rows (findings 2, 3, 8);
   adopt_or_blend_opal_result becomes a generic blend_table_result(r, table_r, w). M-L; A vintage by vintage.
3. Explicit-state pass: a README "state" section listing each module-level state (Saha warm start, SCV
   idtt/idp and the table's composition, OPAL load flags, t_row_index, quad caches and the esac06 stencil
   find_opal2006_ramp reuses, DH composition, mhd_output); then reduce it: DH composition as argument (A),
   SCV table records its composition (A), OPAL reloads when its table path changes (yrec_reset.f90:14-19
   currently assumes "same configuration") (A), Saha guesses in a caller-owned record (B). M.

## Newcomer trace "add an EOS source" (today ~10-11 files vs ~5 for rotation's add-a-mechanism)
eos/<name>/<name>_eos.f90; eos/<name>_eos_tables.f90; eos_configuration flag (exclusivity unenforced);
eos_lib (15th positional eos_init argument, read call, branch in eos_get); for a blended region
blended_eos.f90 (OPAL-only loop, coverage decided in three places); setup/microphysics_settings.f90;
defaults/controls_registry.tsv + generators; core/run/run_setup.f90; test_eos.f90 (configure_eos, both
eos_init sites, regenerate expected); README. The README recipe (README.md:54) omits several. Stuck points:
constants only set at run time; beta14 (two meanings); in_atmosphere; MHD as template is broken.

Fine as is: eos_results.f90, eos_configuration.f90 (bar finding 8), yale_eos_tables.f90, scv_eos_tables.f90,
mhd_eos_tables.f90, eos_mixture.f90 (fix two comments).

# YREC modernization roadmap -- remaining work

This file lists what is NOT done yet. Completed campaigns are
documented in the git history of this file and in the commit
messages; the one-line summary:

- Phase 1: every COMMON block -> modules/derived types.
- Phase 2: per-domain facades (eos_get/kap_get/atm_get), domain
  folder reorganization, misplaced-file relocations.
- Phase 3 stages 1-4: public/private boundary closed + boundary
  checker; per-domain libraries + standalone byte-pinned tests;
  ierr-not-stop across the libraries; named-index result arrays;
  facades merged to one public query each (eos_get/kap_get).
- Phase 4/5: the star layer (one type(star_info) `star`), run_yrec/
  evolve_step drivers, star_job + evolve_state, re-entrancy
  (yrec_reset; second in-process run byte-identical), CI.
- Phase 6: state taxonomy finished; physics purity (eos/kap/net
  star-blind); wind reclassified star-layer.
- Legacy campaigns: 504 numbered DOs -> block DO; 249 blanket saves
  -> 10 annotated INTENTIONAL; 698 gotos -> 0. All byte-gated.
- libyrec + pyyrec: embeddable engine (yrec_capi, make lib), ctypes
  binding, CLI-oracle acceptance test; `make yrec_libs` layout check.
- MESA conventions: star%job/star%ctrl nested, everything else flat
  on MESA vocabulary; controls -> star% campaign complete
  (generated controls_state/sync from the reader-internal buffer).
- Inlist revamp: &star_job/&controls with readable names (registry-
  driven, defaults/*.defaults generated), converter + byte-identity
  tests. ALL 1163 legacy deck pairs migrated to single inlists
  (the 4 solar pairs remain as legacy-reader fixtures); stray
  fort.* writers retired.
- Readability: double negatives, named indices, dead F77 deleted,
  giant procedures decomposed; core/ driver program restructured.
- Stitched-model restructure: one authoritative full-star assembly
  per step (stitched_model_lib); writers are pure readers; turnover
  is an observable; profile gradT/grada swap bug fixed.
- Build correctness: generated deps.mk with real .mod prerequisites
  (stale-ABI segfault class eliminated); make clean purges *.mod.
- MESA-style output: ONE output path (use_legacy_output retired);
  history/profiles/pulse (GYRE/FGONG/GSM behind USE_HDF5=1) over the
  full extended model; io/yrec_output split into coordinator +
  output_columns + history_output + profile_output;
  defaults/{history,profile}_columns.list are the AUTHORITY on the
  default column selection (compiled in via gen_default_columns.py,
  make regenerates on edit).
- Seismic observables: nu_max + delta_nu_rho (scaling relations, with
  nu_max_sun/delta_nu_sun/Teff_sun controls), delta_nu (sound-travel
  integral) + delta_Pg (Brunt integral over the innermost cavity) --
  opt-in history columns; MESA-style structure-limit stops
  (log_L/Teff/log_g/nu_max upper+lower, per kind card, &star_job).
- Numbering/constants: fresh start-model loads reset the model
  counter (set_initial_model_number, default 1; <=0 keeps the stored
  counter); solar_mass_cgs + G_cgs are controls (G sentinel keeps the
  historical log10 G bit-for-bit); registry-default doc bugs fixed
  (lfirst is .true., niter3 is 2).
- Reproducibility (2026-08-30): math/math_lib.f90 with intrinsic and
  crmath backends (USE_CRMATH=1; MESA SDK is a prerequisite),
  -ffp-contract=off everywhere, all real-exponent ** -> pow()/exp10(),
  contract enforced by check_boundaries; ~35% cost, drift 1e-11.
  src/install script (prereq checks + build + env report).

Standing verification: the tiered policy in GUIDELINES.md
("Verification tiers") -- gate1 per commit, owning suite for the
touched subsystem, full battery only for type-layout/output-format/
cross-cutting changes, with the deliberate reseed ritual for
intentional output changes. Never run two batteries/builds
concurrently in the tree.

---

## 1. ierr-not-stop -- DONE 2026-08-30

Every optional-ierr fallback is gone: ierr is REQUIRED on all 16
library entries (module procedures, so the compiler enforces every
caller), the eos_get dispatch collapsed, and the integrand-callback
protocol (bsstep/mmid/atmosphere_derivs/envelope_derivs) carries
ierr -- resolving envelope_derivs' documented residual. New error
paths reach the drivers through compute_scale_height ->
overshoot_boundaries -> find_convection_zones, massloss/mdot,
burn_mix_extrapolated, turnover, and the opal92 interpolation chain
(whose EXTRAPOLATION/X-GRID stops are now ierr).

Remaining stops, all deliberate: core/main.f90 (the CLI exit),
numerics_lib x5 (ludcmp/tridia/rational-interp singular-matrix
stops -- a separate hardening if ever needed), condopacp x3
(conductive-table setup; init-time).

## 2. envint purity split -- DONE 2026-08-30

Phase A: the per-point side channel (star%pulse's q* scratch +
star%current_*) moved out of star_info into core/point_scratch_lib
(pt_scr; yrec_reset snapshots it) -- clearing the long-standing
"pulse q* residue" note. Phase B: core/envint_kernel.f90 hosts the
star-blind integrate_envelope_atmosphere (surface boundary, gray
atmosphere, envelope integration; atmo_struct/env_struct are its
product), configured by an envint_config the driver builds;
core/envint_lib's atm_get is now a 98-line driver applying the
star-state consequences (pphot, Allard gray fallback). Enforced by
check_boundaries' STAR_BLIND_FILES rule. atm_get's remaining 21
arguments are all genuinely per-call (each caller computes its own
L/g/Teff; the triangle-vertex protocol; the fixed integrand-callback
tail) -- no de-tramp forced.

## 3. De-tramp queue

By argument count: henyey_eliminate (39), secular_transport (37),
mhdst1 (36), engeb (35), microdiff_mte (33), temperature_gradients
(33), mid_timestep_model (29). CONSTRAINT for
temperature_gradients: it is on the star-blind kernel's call path
(envelope_derivs -> temperature_gradients_r), so its de-tramp must
NOT absorb star% reads -- explicit-argument consolidation only. The
eos engines (eos_eval/eqstat, 30-31) are documented internals --
lowest value.

## 4. Log verbosity + terminal output

CASE.log still carries full solver iteration traces; add a
verbosity control, MESA-style terminal output behind a control, and
quiet the .short-style config echo (which still prints the legacy
namelist group for new-style runs) in the same stroke.

## 5. Reproducibility follow-ups

- Cross-platform acceptance: same inlist under USE_CRMATH=1 on
  macOS-arm64 and Linux-x86_64, byte-compare (user runs on the HPC;
  needs the MESA SDK there and the same gfortran major -- ideally
  the SDK's own gfortran on both ends).
- Decide the default backend (intrinsic = fast, crmath = provenance;
  flipping the default is one Makefile line + one rebaseline).
- Optional: a CORE-MATH backend under math_lib (correctly rounded at
  near-libm speed) to shrink the ~35% cost; contained to one file +
  link flags.
- Algorithmic (helps both backends): cache linear+log forms in the
  EOS/opacity interpolation hot paths to cut exp10/log round-trips.

## 6. Blocked / small / optional

- MHD table data: no tables ship with the repo; obtaining them is
  the only way to pin the MHD path (starin bug fix and the
  hsubp/sconvec/massloss MHD extension remain unverified).
- FGONG glob(6) mixing-length alpha and the central d2P/drho2
  globals are zero -- fill if a consumer needs them.
- other/ hook system (MESA-style null-default procedure pointers):
  lowest value; only if a concrete A/B use case appears.
- Domain test binaries write test_X.short files in cwd (legacy
  naming leftover).

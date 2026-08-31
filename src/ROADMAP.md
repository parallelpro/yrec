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

## 3. De-tramp queue -- DONE 2026-08-30 (mhdst1 deferred)

All byte-gated (gate1 + 19-case config matrix per commit):
- henyey_eliminate 39 -> 3: star% arrays read directly; the 30
  per-shell equation terms are two henyey_shell_terms records
  (prev/cur, "prev = cur" replaces the 15-line copy block); module
  procedure, relocated setup/ -> core/.
- secular_transport 37 -> 15 and mid_timestep_model 29 -> 16: the 13
  midpoint-in-time structure arrays moved into rot_scr (they are the
  seculr/midmod pipeline workspace; yrec_reset already snapshots
  it); seculr's 10 star-model args read from star directly; both are
  module procedures; midmod relocated setup/ -> rotation/.
- engeb 35 -> 20: the 15 trailing reaction-rate arrays (vestigial
  COMMON relay) written to star%reaction_rate_* directly; knock-on:
  timestep_limit_hburn keeps them as locals (real rates->eqburn flow
  at one zone), hburn/heburn/compute_timestep drop the relay.
- microdiff pipeline: mte 33 -> 15, run 27 -> 13, coefficients
  13 -> 9 via a microdiff_grid record (eq/eq_mid instances); pinned
  by the cm_settle_he* matrix cases.
- temperature_gradients 33 -> 19: merged with its _r result-array
  wrapper (named-index values unpacked into locals at entry); all
  callers now use the plain name. Per the recorded constraint, no
  star% absorption -- it is on the star-blind kernel's integrand
  path, and the module header now says so.

DEFERRED: mhdst1 (36). Analysis done: its 13 trailing args
(log10t_down..mass_fraction_up) are the same caller scratch at all 8
call sites and would become locals (36 -> 23); but the MHD path has
ZERO runtime coverage (no tables ship -- item 6), so a byte-gate
passes trivially and verifies nothing. Do it together with the MHD
table acquisition in item 6.

The eos engines (eos_eval/eqstat, 30-31) are documented internals --
judged not worth touching.

## 4. Log verbosity + terminal output -- DONE 2026-08-31

Most of this had already landed with the MESA-style run log
(b38ec9c: compact progress line on terminal + run.log,
terminal_interval, report_solver_diagnostics, run.log default
location, inlist_used; the .short config echo was retired with the
one-output-path change). The 2026-08-31 close-out:

- Verbosity sweep: the settling ITERATION trace, the per-rezone
  TOTAL J / K.E. bookkeeping, and the per-model settling-suspension
  messages were still unconditional -- the GS_rot solar log was 8586
  lines, now 39. Forensics behind report_solver_diagnostics;
  suspension messages print once per suspension (reset-covered
  latch star%settling_suspended_reported) + SETTLING RESUMED.
- End-of-run summary: 'run finished: <reason>' + final model/age to
  terminal and run log (star%termination_reason, set by
  init_stop_conditions default 'model budget exhausted' and
  overridden by the age / central-abundance / structure-limit stops
  and the calibration verdicts); wall-clock line terminal-only so
  the run log stays byte-pinnable.
- Banner: single-inlist runs print 'inlist : <file>' instead of the
  phantom legacy 'PHYSICS namelist : yrec8.nml2'.

NOT done (judged out of scope): MESA-style retry/backup lines on
timestep cuts; a header-reprint cadence beyond the existing
per-card + every-10-lines reprints.

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
- PRE-EXISTING RUN FAILURE (found 2026-08-31 during the verbosity
  rebaseline, reproduced identically on 249b60d): run_from_zahb_to_
  tahb/Solar_m1p0feh+0p0_GN93_MESA_TAHB dies with a NaN in the MLT
  cubic (tpgrad sqrt-of-negative -> ierr stop) near model 950. Not
  part of the standing battery; needs its own investigation.

---

# 2026-08-31 full-codebase audit -- next roadmap

Folder-by-folder sweep (all 235 .f90, cross-cutting greps, -Wall
harvest) after the pulse-physics fixes (exact Brunt N2, kap/eps
derivative factors, GSM version, atmosphere radii). Items ordered
bugs -> design -> validation; each carries its verification tier.

## 7. Correctness bugs found

- **rezone flag_point overflow** (setup/rezone.f90): the discontinuity
  scan guards flag_count at 100, but the THREE appends after the loop
  (overshoot_base_zone, fine_zone_base, star%nz) have no bound check --
  a model with >~97 flagged points writes past flag_point(100).
  Fix: size by parameter, guard every append. Also
  radiative_zone_bounds(13,2)/convective_zone_bounds(12,2) have no
  overflow guards. Tier 1.
- **Library stops missed by the ierr campaign**: kap/conductive/
  condopacp.f90 (3 raw stops on table range) and numerics_lib
  (ludcmp 'Singular matrix', tridia x2, polint den<1e-20). Library
  code must return ierr; a singular Henyey matrix should surface as
  numerics_termination, not kill a batch job. Tier 1 + owning suite.
- **Single-precision literals in double expressions**: burn_lib (18:
  reaction-rate and neutrino-flux constants like 2.79e-8, 1.017677e-4),
  net_lib (3), turnover_timescale (1.0e20 guard). Silent truncation to
  ~7 digits; also off-message for the crmath reproducibility story.
  Convert to d-literals. BYTE-CHANGING: deliberate rebaseline; check
  the solar-pin drift is at rounding level. Tier 3 (output change).
- **rezone ceiling idiom** `mod(dp,dp).ne.0d0` (twice): float-equality
  as a ceiling test -- effectively always true, so it over-counts by
  one point when the division is exact-in-reals. Replace with
  ceiling(); byte-gate (expected: identical except pathological
  spacings). Tier 1.
- **eqstat dead+broken SCV derivative branch** (~line 525): reads
  specific_heat_cp_2/adiabatic_gradient_2 that nothing assigns;
  unreachable (do_scv_derivatives hardcoded .false.). Delete the
  branch. Byte-safe. Tier 0.
- **read_starting_model core extension** uses an ideal-gas
  logRho = logP - logT - offset ("MHP 4/12 replaced broken eqstat
  call"). eos_get exists now -- use it for consistent extended-core
  densities. Affects only loads that extend the model inward. Tier 2.
- **Float-equality guards on physics values**: burn_lib
  hydrogen_fraction.eq.0.0, microdiff_coefficients species .eq.0.0.
  Each is probably benign (exact-zero sentinels) -- audit and either
  document or convert to explicit sentinels. Tier 1.

## 8. Design debt

- **Silently-static local arrays** (24 x -Wsurprising): gfortran moved
  large locals (light_element_save, compute_seismic_columns' work
  arrays, henyey solve buffers, ...) to static storage -- implicit
  shared state that undercuts the yrec_reset re-entrancy guarantee and
  forbids threading libyrec. Convert to allocatable locals or owned
  module workspace; then the warning class becomes a boundary-checker
  assertion (zero tolerated). Tier 2 (reentry suite owns it).
- **Blanket save + manual reset lists** (evolve_step, run_yrec):
  adding a local silently leaks state across in-process runs unless
  the author remembers the reset block. Replace with a derived-type
  state (default initializers; reset = default-init assignment) so
  forgetting is impossible. Tier 2 (test_reentry).
- **-finit-local-zero masks uninitialized-variable bugs** (it made
  them deterministic zeros instead of crashes). Add a CI/dev lane
  building with -Wmaybe-uninitialized (+ GNU_DEV_FFLAGS run of the
  solar case), fix what it finds, then decide whether prod keeps the
  flag as belt-and-braces. Tier 2.
- **Numerical Recipes provenance** in numerics_lib (SPLINE, LUDCMP,
  LUBKSB, POLINT, MMID, bsstep lineage -- self-documented in the
  headers): NR's license does not permit source redistribution, and
  this repo is public. Replace: LAPACK (already an SDK dependency)
  for LU; independently-licensed spline + Bulirsch-Stoer, or
  clean-room rewrites. Byte-rebaseline per replacement. The single
  highest-liability item in the repo. Tier 3.
- **Pulse column magic numbers**: build_pulse_points fills pts(35,*)
  by bare index; every consumer (FGONG/GSM/GYRE-ext writers) indexes
  numerically. This audit had to reverse-engineer the map to fix
  N2/kap/eps. Add named indices (ipul_r=1, ..., ipul_N2=8, ...) in
  stitched_model_lib and use them everywhere. Tier 0/1.
- **Duplicate physics site**: io/write_gyre_pulse.f90 (the
  pulse_gyre_interval path) recomputes its columns independently of
  build_pulse_points -- that duplication is exactly why the N2 bug had
  to be fixed twice. Either make it consume the interior slice of
  stx_pulse, or retire the interval path outright (the profile-coupled
  pulse stream has superseded it). While there: its dummy
  `log_luminosity` actually receives LINEAR L/Lsun -- rename. Tier 2.
- **Controls sanity checker**: the template shipped
  overshoot_alpha_envelope = 0.5 with envelope_overshoot_active unset
  -- silently inert physics for years. Add a startup pass that warns
  on inconsistent combinations (alpha>0 with its active flag off,
  lovmax with betac=0 [caps overshoot to zero], rotation_active with
  omega=0 [exists], diffusion+rotation exclusive flags, ...). The
  registry gives the list a natural home. Tier 1.
- **Unused variables** (54, hotspot read_starting_model with 24) and
  unused dummies (7): delete; gets the -Wall build to (near) silence
  so new warnings are visible. Tier 0.
- **Rotation naming unification**: the alias table in
  secular_transport's header (QWRMAX/HRU/COD2/... each with 2-3 names
  across callees) is debt the de-tramp left; one name per quantity.
  Tier 1.
- **json = 5000 static sizing** (69 arrays in star_info alone):
  acceptable footprint (~75 MB RSS measured), but it is a hard zone
  cap and static state. Long-term: allocatable star arrays. Low
  urgency; pairs naturally with the threading item. Tier 3.
- **dlog/dexp/dabs archaic intrinsics** (~330 across 27 files):
  mechanical modernization sweep, byte-identical. Tier 1.

## 9. Physics validation tier (new)

The N2 bug survived every byte-pin because pins freeze outputs, not
truth. Add a small validation layer that checks physics against
independent references:

- **Pulse-file contract test** (pytest): read a produced GSM file and
  independently re-derive N2 from rho/P/Gamma_1/r, check kap/eps
  derivative column magnitudes, monotonic r, version attr, column
  set. Directly encodes this week's four bugs against regression.
  Cheap; run in CI. Tier: owning suite.
- **Cross-code benchmark**: one (M, Z, alpha) track YREC vs MESA;
  compare Teff/L/R at matched central-hydrogen points and GYRE
  frequencies at matched Delta_nu. Documented tolerances, manual
  script under testsuite/. Would have caught the N2 bug on day one.
- **Rotation benchmark**: reproduce a published YREC solar spin-down
  or open-cluster omega evolution; pin the omega(r) profile at
  landmark ages. The rotation folder is the largest body of physics
  with no truth-level test.
- **Seismic-observable spot checks**: delta_nu (sound-travel) and
  delta_Pg (Brunt integral) audited clean this pass -- add them to the
  contract test against GYRE-derived values from the same model so
  they stay clean.

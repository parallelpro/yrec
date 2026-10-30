# YREC modernization roadmap -- remaining work

> File paths in this document are relative to `src/` unless they start with `docs/` (moved here from `src/` in 2026-10).

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
- Inlist revamp: &control/&physics (named &star_job/&controls from
  2026-08 to 2026-10) with readable names (registry-driven, defaults/*.defaults generated), converter + byte-identity
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
  (log_L/Teff/log_g/nu_max upper+lower, per kind card, &control).
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

## 7. Correctness bugs found -- items 1-5 DONE 2026-08-31

Done (commits fb28dac, 212a652, d853a29): rezone flag_point guarded
appends + ceiling(); the inherited radiative_zone_bounds(1,2)=M typo
(convec.f) fixed alongside; ierr on ludcmp/polint/tridiag_gs/ctridi/
tridia + the condopacp/cinterp3 range stops, threaded through the
microdiff, seculr, alsurfp and rotation_shape_factors chains;
eqstat's dead SCV-derivative branch deleted; the 19 e-notation
single-precision constants in burn/net/turnover converted to
d-notation (solar drift at the 9th significant digit; net baseline
reseeded). Also done from section 8: the ipul_* named pulse-column
indices (008b12e) and the controls sanity checker (b57c1a3, warns
terminal-only on inert combinations). Float-equality audit and
read_starting_model's ideal-gas core extension remain below.

## 7-remaining. Correctness bugs still open

(2026-09-03 audit of this list: the rezone flag_point overflow, the
library stops in condopacp/ludcmp/tridia/polint, the rezone
mod(dp,dp) ceiling idiom and the eqstat dead SCV derivative branch
were all already fixed in earlier passes and are dropped from here;
the single-precision literal item is DONE in Batch 3.)

- **read_starting_model core extension** uses an ideal-gas
  logRho = logP - logT - offset ("MHP 4/12 replaced broken eqstat
  call") for the shells it adds inside the first original shell.
  Reviewed 2026-09-03: it is an initial guess that the first Henyey
  solve relaxes, calibrated on the first original shell, so the
  error is second order over a few shells -- not a correctness bug.
  eos_get could replace it for consistency (Tier 2); low value.
- **Float-equality guards on physics values**: burn_lib
  hydrogen_fraction.eq.0.0d0 (skip the pp/CNO block in a hydrogen-
  free zone), microdiff_coefficients species .eq.0.0d0 (skip Thoul
  where the species is absent in this and the next zone). Reviewed
  2026-09-03: both are exact-zero "absent" sentinels; a tiny nonzero
  value takes the full path to the same answer. Leave as is.

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
- **Numerical Recipes provenance** in numerics_lib (SPLINE, SPLINT,
  LUDCMP, LUBKSB, POLINT, HUNT, LOCATE, MMID, bsstep lineage --
  self-documented in the headers): NR's license does not permit source redistribution, and
  this repo is public. Replace: LAPACK (already an SDK dependency)
  for LU; independently-licensed spline + Bulirsch-Stoer, or
  clean-room rewrites. Byte-rebaseline per replacement. The single
  highest-liability item in the repo. Tier 3.
- **Pulse column magic numbers**: build_pulse_points fills pts(35,*)
  by bare index; every consumer (FGONG/GSM/GYRE-ext writers) indexes
  numerically. This audit had to reverse-engineer the map to fix
  N2/kap/eps. Add named indices (ipul_r=1, ..., ipul_N2=8, ...) in
  stitched_model_lib and use them everywhere. Tier 0/1.
- **Duplicate physics site -- DONE (86b3fa63, 2026-10-29: the
  interval file is the stitched model's interior points, written by the
  one GYRE text writer; bug_sweep M78)**: io/write_gyre_pulse.f90 (the
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

## 10. File-by-file bug sweep -- 2026-08-31

Eleven parallel domain reviews covered all 231 non-test .f90 files
(~58k lines); every finding carries F77 provenance against 6cd5673.
FULL PER-DOMAIN REPORTS: audit/bugsweep-2026-08-31/ (committed).
Summary below; "VERIFIED" = independently re-derived in the main
session, others are agent findings at the stated confidence.

### Fixed during the sweep

- rotation_shape_factors called with 10 args (missing the new
  required ierr) from henyey_iterate + rezone -- a 212a652
  regression (caller grep truncated by head); UB on every rotating
  run. Fixed + pushed (955a20d). Lesson: caller sweeps never through
  head; bare-external signature changes get a multiline-aware audit.

### VERIFIED high -- physics-affecting, fix candidates in order

1. FIXED afff605 (impact measured: p modes +0.10-0.72 uHz depending
   on nu_max, g-dominated modes untouched -- see the commit message).
   stitched_model.f90:529-530 (new-code): the PULSE builder swaps
   grad/grad_ad in the ENVELOPE region -- env_gradients order is
   (rad, ad, actual) per its own profile-writer comment, but the
   pulse block reads (2) as grad and (3) as grad_ad. Thermal N^2
   sign flips across the convective envelope of every GYRE/FGONG/GSM
   file. The profile-side twin of this exact swap was fixed; the
   pulse site was missed. FIX FIRST -- feeds the live science.
2. core/burn_lib.f90:1196 + :645 (inherited, engeb.f/deutrate.f):
   d(p,gamma)3He uses T9^(+2/3) where CF88 has T9^(-2/3) -- the
   code's own log-derivative (cc13*(-2)) proves the intent.
   Deuterium burning off by ~T9^(4/3): ~1e4 at 1e6 K. Affects
   pre-MS/birthline deuterium physics.
3. net/net_lib.f90:343 (inherited, sneut.f): tfac2 has
   `cvp*cvp - cap-cap` (parses as cvp^2 - 2*cap) where the Itoh/
   Timmes form is cvp^2 - cap^2 -- the classic sneut typo MESA
   fixed. Wrong flavor coefficient in pair/photo/brem neutrino
   losses (matters for late evolution / cores).
4. kap/conductive/condopacpint.f90:119 (inherited -- the F77's live
   line has the same copy-paste, correct formula commented out
   beside it): conductive dlnkap_dlnT is assembled from the RHO
   derivatives (dlnkap_dlnrho_*); the T derivatives are computed and
   discarded. Corrupts the Jacobian whenever Potekhin conduction is
   on.
5. io/model_to_equal.f90:125/232 (inherited, HCOMP(8,...)): the
   metal-diffusion ("MASS FRACTION OF METALS") run is built from
   composition(8,...) = N15 (~1e-8), and equal_to_model applies the
   result to slot 3 = Z. use_diffusion_z metal settling in the
   rotation path has been wrong since the COMMON era.
6. eos/opal/esac.f90:117-123 (modernization): the mass-fraction
   guard's GO TO 61 became a single-statement IF guarding only the
   first write; ierr=1+return run unconditionally on the first call
   -- every use_opal95_eos run dies at its first OPAL-regime EOS
   point (esac01/06 converted the same construct correctly).
7. numerics_lib intpt (:1864, modernization): the table-scan's
   GOTO-exit became `return` -- returns with outputs unset for any
   in-range temperature; breaks all 11 mhdpx2 call sites (entire MHD
   EOS path; zero coverage). Fix is `exit`.
8. util/timestep_limit_heburn + compute_timestep (modernization-
   aggravated; PROPOSED ROOT CAUSE of the known TAHB NaN): helium_dt
   is intent(inout) read-before-write on the "core Y below atime(1)"
   branch (documented in the file header); the caller's local lost
   the F77 SAVE carry, so it reads 0 under -finit-local-zero ->
   helium_dt = const*0 = 0 -> timestep collapses -> Inf/NaN -> the
   observed DELDEL=NaN in tpgrad at core-He exhaustion. Needs a
   reproduction run to confirm before fixing.
9. core_physics agent, VERIFIED reasoning, needs care:
   stitched_model's eps_eps_T/rho columns are DOUBLE-multiplied by
   eps -- engeb's accumulators are already eps*dlneps/dlnT, so the
   2026 "absolute derivative" pulse fix over-corrected (the kap half
   of that fix is right, the eps half wrong).

### High-confidence agent findings (spot-check before fixing)

- rebuild_envelope.f90:266 X/Z fitting-point interpolation sign
  reversed (inherited getnewenv.f; parallel branch correct);
  rebuild_envelope omega() read-before-set (comment claims a SAVE
  that does not exist).
- envelope_derivs.f90:94 dydx(3) multiplies by fp where dydx(1)
  divides -- rotating-envelope structure error ~fp^2 (inherited
  qenv.f).
- mid_timestep_model deuterium-rate arrays lost their F77 SAVE:
  mid-substep D-burning rates ~0 after the first substep
  (modernization; rotating runs).
- burn_mix_extrapolated lost cross-call SAVE state: BS extrapolation
  silently no-ops at order>=2 (modernization; rotation path).
- diffuse_composition missing `else i1=zone_end`: He4 rebalance over
  an undefined range when the unstable region's top is radiative
  (inherited mixcom.f).
- am_advection_diffusion_coeffs: shear/GSF eq-grid coefficients
  never assigned -> those mechanisms contribute ZERO in LDIFAD
  advection-diffusion mode (inherited).
- compute_quadrupole drho/dr uses rho(i) where rho(i-1) belongs
  (inherited copy-paste).
- trapzd interpolates rho/omega^2/eta2 at a loop-invariant offset
  (+slope*del) instead of at y (+slope*(y-b1)) -- the rotation <g>
  shape integrals converge to the wrong value (inherited).
- bsstep's step-underflow guard tests the pass-through hydrogen
  fraction instead of the independent variable (inherited NR-form
  deviation).
- map_user_inputs: qs0e/qqs0ee derivative scales for reactions 2-7
  divide by the PP S0 instead of each reaction's own (~1e-22
  factor); Seff derivative corrections effectively zeroed under
  use_new_nuclear_rates (inherited).
- Noerdlinger dlnLambda/dX settling correction added with the wrong
  sign (gravitational_settling_setup:379, inherited); grsett
  midpoint off-by-one + stale element 1; convergence tests missing
  abs() in both settling solvers (inherited).
- engeb Itoh branch: PET = PEP + DSNUDT copy-paste and log-vs-
  absolute derivative mixing (inherited).
- liburn/liburn2 radiative_frac = intended-1 (weights in [-1,0]) for
  zones leaving a retreating CZ -- pre-MS Li depletion (inherited).
- OPAL95 opacity ll95tbl slot collision (110/119 double-booked;
  X=0 tables clobbered for Z>0.04); getopal95 density_shifted read
  uninitialized on the low-X/low-T shortcut (modernization);
  alex06tab X-node(4)=1-Z regardless of index_x for Z>0.06.
- mhdpx1 out-of-range T returns ierr=0 with stale/zero output
  (inherited; MHD); eqstat ionization-cutoff blend adds raw Saha
  QDTT/QDTP instead of the differenced form (inherited, both eras);
  SCV ragged-edge unclamped reads (inherited);
  fully_ionized_eos Newton non-convergence returns unset outputs.
- read_controls:2317 lfirst(1)=.true. safety force hits a dead local
  (modernization; fixed 2026-10-29 in fb46cfe8, bug_sweep M47: card 1
  of ten shipped decks had evolved an empty model); output_columns
  append_column unbounded vs
  max_cols; read_starting_model core-extension overflow prints "RUN
  TERMINATED" but continues (inherited).
- atm: atm_get_surface_pt drops alsurfp's fatal ierr (modernization);
  surfp/kcsurfp high-gravity branch ignores the gmax index
  (splines over -999 sentinels; inherited); turnover_timescale
  spline_taucz_done never set on the innermost-point branch
  (goto-elimination regression) + core_cz_top_index read
  uninitialized for sun-like structures.
- matt_wind PMM path ignores wind_loss_active; check_angular_momentum
  reversal threshold 1.0 rad/s (dead guard) and cut-once logic;
  wind_spindown_matt gl-vs-cgl fcen ~1.5e7x (sharpened known item).
- rezone gradient-insertion has no json bound check (the audit fixed
  only flag_point).

### Medium/low + weak observations

See audit/bugsweep-2026-08-31/*.md -- ~35 medium/low findings and
~120 one-line weak observations, each with provenance. Also
recorded there: verified-clean lists (Henyey algebra, seismic
integrals, Thoul solver, FGONG slot map, mod-file symmetry,
controls adoption copies, conductive combination formula).

### Suggested fix order

Batch 1 (science-first, output-changing, one reseed): pulse
grad/grad_ad swap + eps double-multiplication + condopacpint
T-derivative + N15-vs-Z metal diffusion.
Batch 2 (physics constants, reseed): deuterium exponent, sneut
tfac2, engeb Itoh branch, liburn radiative_frac, Noerdlinger sign.
Batch 3 (crash/UB class, byte-safe or uncovered): esac guard, intpt
exit, heburn dt=0 (after reproducing the TAHB crash), dropped-SAVE
family (rebuild_envelope/mid_timestep_model/burn_mix_extrapolated/
massloss), read_starting_model overflow stop.
Every batch: verify the specific claim first (agents are good but
not infallible), fix, byte-gate or deliberately reseed per tiers.

## 11. Second independent bug sweep -- 2026-09-01

Twelve parallel reviews (A-L) on a partition cut along different
seams from section 10, with the section-10 findings deliberately
withheld from the reviewers so that overlap is corroboration rather
than echo. Same F77-provenance rule (6cd5673). One reviewer (C)
accidentally saw ~15 lines of section 10 mid-sweep and disclosed it;
its helium_dt finding is therefore counted as weakly corroborated
only. FULL REPORTS: audit/bugsweep-2026-09-01-pass2/ (committed;
INSTRUCTIONS.md records the protocol).

Totals: 10 high / 42 medium / 55 low + ~150 one-line weak items.
"VERIFIED" = re-derived in the main session against the code and
the F77.

### Concordance with section 10

Corroborated independently (pass 2 reached the same defect from the
code alone): sec-10 verified highs 2 (deuterium T9 exponent, B), 3
(sneut cap-cap, B), 4 (condopacpint dlnkap/dlnT from rho
derivatives, K), 5 (N15-as-Z metal diffusion, G -- rated latent),
6 (esac OPAL95 guard, I), 7 (intpt return-vs-exit, H AND J; J adds
the mechanism: zeros -> meqos 1/0 -> every lmhd run dies at its
first EOS call), 8 (helium_dt read-before-write, C, weak); and the
agent-level items rebuild_envelope omega read-before-set (A),
rebuild_envelope X/Z refit sign (D), mid_timestep_model D-rate SAVE
loss (E), burn_mix_extrapolated SAVE loss (E), diffuse_composition
missing else-i1 (F, now rated HIGH), shear/GSF coefficients never
assigned (F), compute_quadrupole rho(i) (F), trapzd offset (F, H),
bsstep guard (H), map_user_inputs S0 scaling (C), settling
convergence tests without abs (G), engeb Itoh copy-paste (B),
liburn radiative_frac (B), ll95tbl slot collision (K), getopal95
density_shifted (K), alex06tab X node (K), eqstat ramp blend (J),
mhdpx1 stale return (J), read_starting_model overflow-continues
(D), atm_get_surface_pt swallowed ierr (L), surfp/kcsurfp -999
sentinels (L -- now checked against the shipped tables: every
Castelli-Kurucz table's last row carries one, so atm_choice=5 above
~48 kK is silently garbage), matt_wind PMM gating (C),
check_angular_momentum 1 rad/s dead guard (F).

NOT corroborated -- verify before fixing: sec-10 item 9 (pulse
eps_eps_T/rho double-multiplied by eps; D read the same columns and
did not flag it); Noerdlinger settling sign (G instead found a
Coulomb-log double count at the same site, see below).

Disagreement resolved for section 10: envelope_derivs.f90:94
dydx(3)*fp vs dydx(1)/fp -- A rated the envelope ODEs clean, but
with x = log P the chain rule gives dlog r/dlog P = -rP/(G m rho
f_P), so both rows must carry 1/f_P whatever the f_P convention;
the inconsistency is real (inherited qenv.f:51/53), O(omega^2)
small for slow rotators.

### NEW in pass 2 -- VERIFIED

1. core/structure/structure_equations.f90:464 + stitched_model.f90:448-455,
   512-514 (new-code): pulse_mean_molecular_weight is filled from
   eos_res(i_gas_constant) = R/mu, and
   pulse_electron_mean_molecular_weight already holds 1/mu_e, so
   profile column 'mu' is ~1e8, 'mu_e_inv' and FGONG var(14) are
   mu_e (inverted), and the Ledoux column gradL takes dln(R/mu) =
   -dln mu (composition term wrong-signed) with a spurious spike at
   the interior/envelope junction where column 52 is 0 -> ln(1e-30).
   GYRE ignores var(14) for adiabatic modes, so NO frequency impact;
   every profile-file analysis of mu/gradL is wrong.
2. eos/opal/esac06.f90:406-411 + rhoofp06.f90:176 (inherited
   esac06.f:477; tolerance change 2025-10-10): on the deriv_order=1
   trial calls from rhoofp06 only eos(1) is re-interpolated, but the
   tail still rescales the cv slot (x moles*R/tmass ~1e2) from its
   stale value on every call, so cv compounds until the next
   full-order call resets it; harmless in production, but it is the
   "eos(5) growing without bound" that the 2025-10-10 loosening of
   the rho(P,T) inversion tolerance from 0.5e-7 to 1e-5 was papering
   over. Consequence today: every OPAL06 run carries ~1e-5 relative,
   non-smooth solver noise in rho(P,T), 200x worse than the 1995/
   2001 paths, with oeqos06's P-consistency check commented out.
   Fix order: guard the tail scalings by deriv_order (byte-safe),
   THEN restore 0.5e-7 (output-changing).
3. eos/scv/eqscve.f90:148,202 (inherited eqscve.f:78/126; eqscvg
   has no smoothing at all -- the "eqscvg same" first written here
   was wrong): the upward-neighbour smoothing weight is 0.5*d/tol
   (weight on the shifted stencil: 0 at the cell boundary, 0.5 at
   d=tol, then a jump to 0) where continuity requires 0.5*(1-d/tol)
   to mirror the downward branch's (d+tol)/(2tol). Instead of
   removing the stencil discontinuity it halves it and adds a second
   jump at d=tol, in both logT and logP, over 40% of every cell.
   Step noise in rho, cp, grad_ad for every SCV-EOS run (the
   run_from_start_to_dbl and giant_differential_rotation OPALSCV
   decks). Reseed on fix.
4. core/model/rezoning.f90:845-852 (inherited hpoint.f "MHP 6/00" block):
   the start-of-step osplin of eps/esum onto rot_scr%old_esum/
   old_eps passes star%log_mass as the table abscissa AFTER the
   transfer loop overwrote it with the new grid (old count, new
   abscissae) -- effectively "old eps at the same INDEX", so it is
   mis-registered wherever points were inserted/deleted below.
   Feeds mid_timestep_model:165,174 -> Eddington-Sweet velocities in
   rotation_stability_setup, every rotating timestep. Fix: move the
   two osplin calls above the transfer loop (they belong with the
   omega/j/fp/ft/r0/eta2 splines).
5. core/burn_lib.f90:1504,1519,1537 (inherited engeb.f:823/836/852):
   dlnrate_dlnT for C13(a,n), C12(a,g), N14(a,g) is dscreen +
   rho/rate*(dS/dlnT) with rate = 1.157e22*rho*exp(screen)*S, so the
   analytic part is divided by ~1e22 -- the Jacobian sees only the
   screening derivative (~0) where the true value is ~20-30. pp/CNO
   use the analytic form and are fine. Hits every He-burning core
   (Henyey convergence in exactly the regime of the TAHB NaN).
6. kap/opal95/ll95tbl.f90:80 sharpened (inherited, DATA
   NZ/...,100,109,118/): with the X=1-Z tables read last (ix=10),
   slot 110 (Z=0.08, X=0) and slot 119 (Z=0.10, X=0) are overwritten
   by the X=0.94/0.92 tables; kap_lib's HB branch (Z>0.1, logT>7)
   calls getopal95 at Z=0.1 exactly, so a He-burning core (X=0)
   gets kappa(X=0.92, Z=0.08) -- ~0.28 dex too high at electron-
   scattering temperatures. Third independent candidate for the
   TAHB NaN besides helium_dt=0 and item 5.
7. io/read_controls.f90 (new-code): 72 DATA statements plus
   initialized declarations make the legacy locals (kttau, clsun,
   wmax, awind, tdisk, tolerances, s0_*, alphac, ...) implicitly
   SAVEd, so a second in-process yrec_run inherits whatever the
   previous inlist set for any control the new inlist omits.
   Breaks the yrec_capi re-entrancy contract; test_reentry runs the
   same inlist twice and cannot see it. Fix: explicit reset block at
   entry (or move them into a derived type reset by yrec_reset).
8. rotation/seculr/zahn_coupling_factor.f90:62-66 (inherited
   getfc.f:44-45): the header defines alpha = 1/2 dln(r^2 Omega)/
   dln r and the variable is even named half_dlnj_dlnr, but the code
   omits the 1/2 -- the alpha*U term in Zahn's f_c is double-
   weighted.

### NEW in pass 2 -- high-confidence agent findings (spot-check)

- atmosphere_derivs.f90:61 (modernization): T(tau) integrand reads
  live star%job%atm_choice instead of the kernel's cfg; on an Allard
  table failure envint_kernel flips only cfg, so that integration
  mixes an HSRA T(tau) with an Eddington start point.
- henyey_iterate.f90:133 -> henyey_coefficients.f90:300
  (modernization): envelope_zone_index is a plain local set only by
  mix at level>2 but passed at every level; the accretion entropy
  correction is applied to every zone at levels 1-2.
- am_advection_diffusion_coeffs.f90:552-559 (modernization): the
  "CORRECTIONS TOO LARGE" GOTO 950 (timestep-cut block) became an
  `exit` of the innermost loop only; the rejected substep closes
  normally and a later converged substep resets the flag. LDIFAD.
- microdiff_coefficients.f90:113-118/184-189 (modernization):
  locals read after a `cycle` that skipped setting them (F77 SAVE);
  rescued only by -finit-local-zero. Li/Be diffusion every step.
- gravitational_settling_setup.f90:245,316-321,350 (inherited):
  variable Coulomb log both fed into every Thoul coulomb_log(i,j)
  and divided into the prefactor while the time unit already
  carries the 2.2 -- rate biased by 2.2/lnLambda (Thoul, non-fit,
  coulomb_log_choice 2|3).
- equal_to_model.f90 asymmetric X floor/cap (inherited): flagged as
  a candidate for the 0.8 Msun "UNABLE TO SOLVE FOR NEW ABUNDANCES
  IN SHELL 1" abort (medium confidence; that message prints
  zone_begin, i.e. the range starting at zone 1).
- burn_lib lirate88 mode 2 ignores its arguments, so "end-of-step"
  Li rates are never computed (both slots hold start-of-step; comment
  at evolve_step.f90:541 disagrees with the code); dburn/dburnm read
  total_shell_mass uninitialized on the single-zone accretion path
  (SAVE lost); Itoh branch also has a 0/0 NaN if TCUT(5) < 7.
- solve_composition.f90:60: composition is intent(out) but read for
  the Newton initial guess (works only because gfortran does not
  exploit undefined-on-entry).
- rescale_model.f90:372 (inherited rscale.f:343): the Z-ramp block
  scales CNO by (Zold-Zc)/Zold (the REMOVED fraction) where the main
  block uses Zc/Zold; negative CNO if Zc > Zold.
- mdot.f90:267 reads zone_mass_grams(0) for fully convective
  accretors; secular_transport.f90:270,289 apply the Matt torque
  without the use_wind_torque gate.
- SCV pressure search lacks max(1,jj): a first lookup at 4.0 <=
  logPgas < 4.2 reads index 0 (eqscve+eqscvg, inherited); SCV ion-
  fraction ramp blend is a no-op (Saha values overwritten first).
- esac.f90 (OPAL95) has no k3==nt guard (t6_grid(nt+1) at the
  coolest rows; 2001/2006 added it); 1995 radsub applies the unit
  revision to P and S but not E/cv; t6rinterp dix (comment claims a
  SAVE) is a plain local.
- condopacpint.f90:130-158: the log rho < -6 extrapolation branch
  flips the sign of log sigma and extrapolates the wrong way.
- alfilein.f90:81/132 tests an uninitialized local latmtptau100
  instead of allard_use_tau100 (inherited typo; the tau=100 guard
  for old NextGen files never works).
- shape.f90 Radau central seed 6 vs 6/7; solid_body_omega "dI/domega"
  is really a*dI/da; max_diffusion_iters unbounded vs 16/50-element
  history arrays; setups.f90:118 Avogadro digit transposition
  6.0222137e23 (only the Debye-Hueckel coefficient, ~2e-5); lir
  returns result_y unassigned on the degenerate-table path.

### Verified-clean (negative evidence, pass 2)

Henyey elimination/solve and all coefficient partials; surfbc
triangle; envelope+atmosphere ODEs (modulo the fp row above); MLT
cubic; turnover walker; seismic integrals; Thoul matrix + LU;
Lax-Wendroff fluxes/BCs; Saha and fully-ionized derivatives (non-DH);
composition Jacobian (kemcom); mass-weighted homogenization; all
convection-zone modes; Kawaler/Matt/Reimers mappings; all three
T(tau) relations; ~240 legacy->canonical control copies and the
legacy-vs-inlist equivalence; controls_sync (364 members each way);
ludcmp/lubksb/polint/mmid/ratext/qgauss/splines vs NR; Lagrange/
spline opacity machinery and the rad/cond blend; OPAL95
extrapolation; all bare-external call sites in every partition
(argument lists match everywhere).

### Revised fix order (supersedes section 10's)

Batch 0 (byte-safe or uncovered, do first): intpt exit; esac guard;
esac06 tail-scaling guard; rezone osplin ordering (rotating runs
only -- check the matrix); read_controls reset block; diffuse_
composition else-i1; microdiff_coefficients locals; heburn dt=0
(reproduce TAHB first); dropped-SAVE family.
Batch 1 (pulse/profile, output-changing, reseed mesa baselines):
mu/mu_e_inv/gradL columns; then adjudicate sec-10 item 9.
Batch 2 (physics, reseed): alpha-capture Jacobian; ll95tbl slots;
condopacpint derivatives; rhoofp06 tolerance restore; SCV weight;
deuterium exponent; sneut; engeb Itoh; liburn; N15-vs-Z; Zahn 1/2;
quadrupole; trapzd; envelope fp row.
Then re-run the TAHB case: items helium_dt / alpha-Jacobian /
ll95tbl are three independent candidates for that NaN; fix batch 0
first and see which survives.

**Batch 0 status -- 2026-09-02, DONE (byte-safe part).** Landed:
intpt exit (numerics_lib); esac first-call guard; esac/esac01/esac06
tail scaling gated on `eos_index_inverse(slot) <= deriv_order`;
read_controls DATA -> entry-time reset block (all 72 DATA statements
and 27 initialised declarations in read_input; values unchanged);
microdiff_coefficients explicit zeros at skipped zones; burn_mix_
extrapolated explicit SAVE for the seven step-state locals; massloss
seeds. Gates: gate1 IDENTICAL; config matrix (16) + 4 solar cases
byte-identical to standard/; MESA pins, test_eos pass. test_reentry
gained a two-inlist form (single-inlist style; A sets s0_pp, B omits
it; B must equal a fresh B) -- verified to FAIL on the pre-fix read_controls and pass
after. Deferred to a rotation-path reseed batch (all output-changing
for rotating cases -- difrot, 3Msolar_GS_rot, Test_solar_GS_rot, the
cm_rot_* matrix cases): rezone osplin ordering, diffuse_composition
else-i1, mid_timestep_model deuterium arrays, rebuild_envelope omega
read-before-set. heburn dt=0 still awaits the TAHB reproduction.

**Rotation-path batch -- 2026-09-02, DONE (0747603).** rezone osplin
ordering, diffuse_composition else-i1, mid_timestep_model deuterium
arrays -> rot_scr, rebuild_envelope -> star%omega. Only the three
testsuite difrotmix cases changed (first difference at model 2 from
the envelope omega seed; final calibrated models within 2e-5 in log
Teff, surface X within 0.2%); reseeded. Example GS_rot case and all
non-rotating pins byte-identical.

**sec-10 item 9 adjudicated -- NOT a bug. REVERSED 2026-10-10: it
was a bug (section 13).** GYRE's mesa_file_m reads eps_eps_T and
divides by eps to get eps_T = dln eps/dln T, and MESA writes
d_epsnuc_dlnT (= eps * dln eps/dln T) into that column. The premise
"pulse_dlneps_dlnt is the log-log derivative used in ql_dt" was
wrong: ql = cccql*eps, so ql_dt takes the ABSOLUTE d eps/dlnT, which
is what engeb returns and what the field held. eps_total*field was
eps times the right value.

**TAHB NaN -- 2026-09-02, ROOT CAUSE = heburn dt=0 (2be7c30).**
Reproduced on the post-batch-0 build: Solar_m1p0feh+0p0_GN93_MESA_TAHB
dies at model 977, the first model with Y_c < atime(1) = 0.002. With
the He-shell branch computing a fresh limit (time_dy_total /
time_dy_shell, previously dead controls) the case runs to its stop
condition at model 998. The alpha-capture Jacobian and ll95tbl slot
items stay in Batch 2 as physics fixes, not crash candidates.

**Batch 1 -- 2026-09-02, DONE.** pulse_mean_molecular_weight now
holds mu = R/eos_res(i_gas_constant); the electron slot was renamed
pulse_electron_mean_weight_inverse and is written straight to profile
column 53 and FGONG var(14); envelope/atmosphere points supply mu from
the gas-pressure identity (R rho T/(beta P), same as eqstat's
specific gas constant) and the envelope now also fills column 16
(beta); gradL drops its composition term where mu is unavailable
instead of taking log(1e-30). Checked on a 50-model solar MESA-mode
run with mu/mu_e_inv/beta/gradL/brunt_N2 enabled: mu 0.62 (centre)
to 1.28 (neutral atmosphere), continuous across the interior/envelope
junction; 1/mu_e 0.84 at the centre; gradL-grada bounded (max 0.046
in the H ionisation zone). No reseed needed: the default profile
column set and the GYRE pulse format do not carry these columns, so
the MESA pins are byte-identical.

**Batch 2 -- 2026-09-02, DONE.** All 14 physics items plus one
modernization regression found during the pre-fix review
(audit/batch2-review-2026-09-02.md has the derivations): alpha-
capture Jacobians (reactions 8/10/11: missing 1/S factor, missing
T9**-3/2 in the resonant terms, dr1 index), deuterium T9**-2/3 (both
deutrate and engeb), engeb Itoh branch (absolute derivatives, and
the T line no longer built from the rho line), liburn/liburn2
radiative_frac orientation, sneut cap*cap, eqscve upward smoothing
weight, rhoofp06 tolerance back to 0.5d-7 (+ bracket check, cap 30,
non-convergence written to the run log), condopacpint derivatives
(conductivity-weighted mean, T inputs), ll95tbl unpacked slot layout
(n_opal95_xz 130), secular_transport GETFC gate back to lvfc (94c7f45
had switched it to ldifad), Zahn alpha/2, quadrupole rho stencil
(zone-1), trapzd (y-b1) offsets, envelope dr/dlogP 1/f_P,
model_to_equal Z slot. 13 of 15 are inherited from the F77 source;
rhoofp06 (2025-10) and the lvfc gate (2026) are not. Reseeded every
pin. Drift: calibrated solar models end within 1.2e-5 in log Teff,
5.7e-5 in log L, 2.5e-4 in X_c, with ~25% fewer models (the
deuterium phase is shorter); the 12/22-model pre-MS matrix cases
move by 0.02-0.31 dex in log L at model 2 already (deuterium
luminosity), the 0.3 Msun start->D_bl case reaches its stop at
1.06e-3 instead of 5.6e-3 Gyr; dbl->ZAMS byte-identical.
expected_test_net.out (deutrate x5400, sneut pair/photo/brems tens
of %) and expected_test_eos.out (rho(P,T) at the 1e-7 level)
reseeded; test_kap and test_atm unchanged. Left open: reaction 11's
N14(a,g) polynomial coefficients (0.177/3.94 vs CF88), whether
rotation_stability_setup:482 wants a logarithmic dlneps/dlnT (resolved
2026-10-10: no -- it needs the absolute one, see section 13),
compute_quadrupole's 1/R^4 header vs 1/R^3 code, radsub06's 0/0 on
the priming call.

**Batch 3 -- 2026-09-03, DONE.** The Batch 2 leftovers plus the
single-precision literal item from section 7-remaining. N14(a,g)F18
linear coefficient 0.117 (both copies carried 0.177, a digit
transposition inherited from the F77; the 3.94 is correct -- it is
2 x 1.97 in the derivative). Every unsuffixed REAL literal in
burn_lib and net_lib (362 inexact, 570 total) converted to d0/d-n
form, the two `**(2./3.)` exponents included. compute_quadrupole
header fixed to say 1/R**3 (the code always used R**3). rhoofp06's
table-priming esac06 call passes radiation flag 0 like every other
call in the file (radsub06 was dividing by a zero cv on that call
-- harmless because discarded, but an -ffpe-trap=invalid hit).
New NaN guard: converged_model_is_nan (stop_conditions) checks
log_L, log_Teff, age/dt and the zone-by-zone logT/logRho/logP/logR/L
of every converged structure before it is written; a hit is a
positive-ierr error (exit 1) with the slot named in the terminal and
run log. Before it, a run that went non-finite kept "converging"
(every comparison against NaN is false) to the end of its model
budget and exited 0. Reseeded every pin. Drift: the 12 testsuite
solar runs end within 5e-6 in log Teff, 3.3e-5 in log L, 8.7e-5 in
X_c (one case 853->854 models); config matrix, standard-solar and
run_from_* cases byte-identical or unchanged in the printed
summary; the ZAHB->TAHB solar case ends at the same model 1108 with
the He-burning age changed by 2e-8 relative. expected_test_net.out
reseeded (5.5e-6 max relative change -- exponent coefficients
amplify the ~6e-8 literal rounding); test_eos/kap/atm unchanged.
Reviewed and closed without a change: rotation_stability_setup:482
is dimensionally consistent as written (engeb's derivatives are
absolute, d eps/d lnT, so f*deps/dlnT + eps*(1-f-chi_T) is fine);
the read_starting_model ideal-gas core extension and the exact-zero
"absent species" guards (see section 7-remaining). Still open and
by design: a negative-ierr numerics_termination (solution diverged)
exits 0; the envelope/atmosphere pulse rows of stitched_model leave
the mu_e_inv/kap_T/eps columns unfilled.


## 12. Readability / maintainability sweep -- 2026-09-03

Seven read-only per-domain reviews (every file read in full),
asking two questions: does the logic make sense on reading, and is
there a clearer way to write the same code. Reports and a
consolidated summary with cross-domain themes and a six-batch plan
(R1 comments/dead code, R2 named indices, R3 shared helpers, R4
explicit data flow, R5 derived-type de-duplication, R6 the
number-changing items with one reseed) are in
audit/readability-sweep-2026-09-03/ (SUMMARY.md first). Eleven bugs
found on the way and verified against the source are in its
section 1.1 -- among them a modernization regression
(read_starting_model:829 linear Lsun where log10 is needed, on a
branch no deck runs), the rhoofp01 twin of the Batch 3 rhoofp06
priming fix, and a one-way swap in envint_kernel's envelope
reversal.

R1-R5 were executed 2026-09-03..04 (7ef705c -> eece94b1, 75
commits, 225 files, +8309/-13186) in three waves of parallel
worktree agents; every commit class A, the 37-pin selection
byte-identical after every wave, aux battery 22/22, ZAHB->TAHB spot
check unchanged. Highlights: named species/rate/table indices
throughout, the henyey flag relays and kenv/katm gone, atm_get
takes an explicit step config, callees report decisions and the
drivers flip star%job%, one kap table-set type per family with two
instances (five clone files deleted), one OPAL-EOS vintage type
(quadeos01/06, t6rinteos*, eqbound*, gmass06 deleted), one envelope
refit helper, one rates vector, engeb returns a record. The
execution record with per-batch detail, the eleven-item list of
author decisions for R6, and the agents' corrections to the audit
is audit/readability-sweep-2026-09-03/EXECUTION.md; the seventeen
agent reports and their briefs are under its reports/. R6 is not
started.


## 13. Architecture audit -- 2026-10-10

Four newcomer tasks (an extra energy term, a new opacity/EOS option, a
rate update / new isotope, a new angular-momentum transport process)
were traced end to end; the findings and the proposed target
organization are docs/TARGET_ARCHITECTURE.md (Appendix B lists the
bugs). docs/GLOSSARY.md is the naming authority; docs/HOW_YREC_WORKS.md
now uses current names. Bugs found (origin noted per item):

- **ldifad dead path -- DONE.** The advective transport option reads
  equal-grid coefficients nothing computes (nor did F77); it is inert
  except in rotating runs with instability transport, where
  secular_transport now stops with ierr = 1 (test_error_exit.py
  test_advective_transport_rejected). No pinned output changes.
- **Pulse eps_eps_T / eps_eps_rho wrong by a factor eps -- DONE.** Not
  inherited: 2ed398ca (2026-08-30) multiplied engeb's absolute
  d eps/dln derivative by eps a second time, misled by the field names
  pulse_dlneps_* (section 11's sec-10 item 9 ruling, now reversed).
  The writers now copy the field; zone_dlnepsilon_*, pulse_dlneps_*
  and rot_scr%dlnepsilon_* are renamed *deps_*. The rotation use
  (rotation_stability_setup) was always right: there f*deps/dlnT is
  added to epsilm*(1-f-chi_T), consistent only for the absolute
  derivative -- this closes the open question in Batch 2's notes.
- **star%am_transport_convective_flag never assigned -- DONE.**
  Inherited: F77 main.f's LCZ(JSON) was a local passed only to
  CRRECT -> GETROT and never set, so the Omega(j) inversion that
  henyey_iterate runs on every Newton iteration (rotating runs) treated
  every zone as radiative instead of enforcing solid-body rotation over
  the AM-convective zones, as getw does. rezone computed the right
  flags but stored them in the write-only cz_flag_start;
  read_starting_model computed and discarded them. Both now write
  star%am_transport_convective_flag; cz_flag_start removed. Drift, all
  rotating cases only (byte-level: 9 legacy logs + the GS_rot MESA
  baseline; 26 non-rotating logs byte-identical): calibrated solar
  models at 4.568 Gyr |dlog L| <= 3e-5, |dlog R| <= 1.2e-5,
  |dX_c| <= 7e-5, log Teff unchanged to 6 digits; along the tracks
  (mostly the fully convective pre-MS) |dlog L| <= 0.0019,
  |dlog R| <= 9e-4. Reseeded.
- **Two d(p,g) rate copies differ by 0.68% -- DONE.** Inherited
  (deutrate.f vs engeb.f): deutrate turned CF88's N_A<sigma v> into a
  per-gram rate with 3.0115d23 = N_A/2, engeb with N_A/2.013553 (the
  deuteron mass, as in its isotope table) -- in the rate and in the
  convective-overturn cap. Now one burn_lib function, dpg_rate, plus
  the constants avogadro_dpg / deuteron_mass_amu, used by both; its
  arithmetic is engeb's token for token (bitwise-equal to the old
  expression at 16,441 (T, rho) points, rate and cap), so only the
  deuterium depletion rate moves (x 0.99327). Drift (32 of 35 legacy
  logs + both MESA baselines; byte-identical: cm_basic_composition and
  the m0030 seed->start and dbl->ZAMS runs): every run ends on the same
  stop; calibrated solar models at 4.568 Gyr |dlog L| <= 4.4e-5,
  |dlog R| <= 1.5e-5, |dlog Teff| <= 1e-5, |dX_c| <= 1.2e-4, models to
  solar age +-4; pre-MS config-matrix cases at model 12 (D burning)
  dlog L = +3.8e-4, dlog R = +2.0e-4; m0030 start->D_bl |dlog L| 1e-5.
  expected_test_net.out (deutrate line, x 0.99327) and all those pins
  reseeded; the GS_rot MESA case now writes 8 profiles instead of 9.

Not yet checked (reported by the audit agents): beta14 defined
differently on the Yale and OPAL/MHD paths; the opacity ramp's
derivative blending; max_delta_j_by_iter(16) overflow; the envelope
overshoot distance gated on envelope_overshoot_active in evolve_step
but on lovstm in evolve_angular_momentum.

### Phase 1 (network) on yrec-next -- see docs/TARGET_ARCHITECTURE.md

- **Table-driven implicit burning -- DONE (6846ee78).** net/network_lib.f90:
  the network as species + reaction rows (rate row, rate law, change
  list); one implicit_burn_step builds the equations and Jacobian from
  the table. kemcom's hand-coded 7-species system is gone
  (solve_composition 423 -> 104 lines). Agrees with kemcom to 3e-14 on
  20,000 random cases.
- **Table self-check -- DONE (9c216c31).** emitted_mass per row; every
  row conserves baryon number; the table is printed and pinned in
  test_net.
- **One rate evaluation -- DONE (35db7fa5).** net_lib's
  nuclear_reaction_rates feeds both rates (composition) and engeb
  (energy, derivatives, fluxes); the duplicated block is gone.
  Bit-identical to both former copies on a 6,150-point grid.
- **S-factor derivative scales normalised per reaction -- DONE.**
  Inherited from F77 remap.f: the S'/S and S''/S scales of reactions
  2-7 were divided by pp's Solar Fusion I S0 (4e-22 keV b), making
  them ~1e-22 and deleting the q2..q5 terms of He3+He3, He3+He4,
  C12+p, C13+p, N14+p and O16+p (the deck S' values were ignored).
  Each now uses its own S0 (5.4e3, 0.53, 1.34, 7.6, 3.5, 9.4 keV b,
  from the F77 data line). With the deck defaults the scales become
  1.24, 1.14, 1.00, 1.004, 0.54, 2.0. Rate changes near the solar
  centre (log T 7.1): He3+He3 -1.9%, He3+He4 -1.3%, C12+p +5.9%,
  C13+p +0.1%, N14+p -4.9%, O16+p -14% (pp unchanged). Solar models
  at 4.568 Gyr: |dlog L| <= 1.4e-4, |dX_c| <= 3.8e-4, |dlog Teff|
  <= 6e-5, |dlog R| <= 2e-5; 16 solar-age pins reseeded, 19 short
  pre-MS/config pins byte-identical.
- **Exact T-derivative of reactions 1-7 -- DONE (89c1fa42).** The
  Jacobian term now carries the S'/S'' scales; test_net checks every
  slot's analytic dln(rate)/dlnT and dln(rate)/dln rho against centred
  differences (pinned). N14+p was off by 5.4e-4, O16+p by 7.7e-3.
- **rates takes the composition vector -- DONE (fb6ddfd6).**
- **Deuterium in the composition path's screening -- DONE (b3f5a322).**
- **LNEWNUC = .false. restored -- DONE (fff3deb5).** A 2026 regression
  (rename campaign dropped the SSTANDARD copy): the flag left every
  pp/CNO rate zero. New test_physics::test_legacy_cross_sections.
- **neutrino_flux_table cutoff guard -- DONE (9be2a129).**
- Deliberately left as special operators (documented in network_lib):
  deuterium (exponential sub-stepping with the overturn cap), Li/Be
  (liburn), eqburn's equilibrium mode, the Be7/N15 branchings. Making
  D an explicit network species changes the pp-chain model itself
  (p+p -> d no longer instantaneous), a physics decision, not a
  restructuring. The species count from the table (the literal 15 in
  ~45 declarations) moves to Phase 5 (state), where the arrays become
  allocatable.

### Network reorganization (yrec-next) -- burning code into net/

Rate rows into net_lib (d3b3d440); nuclear energy generation as the
side-effect-free net routine nuclear_energy with engeb as its adapter
(b8fb1ff7); eqburn into network_lib (5153a25e); deuterium (deutrate ->
deuterium_rate, dburn, dburnm) into net/light_elements.f90 with the
model state as arguments. All byte-identical. Bugs seen on the way,
left for a decision (physics changes):
- **dburnm adds grams to a dimensionless 1 -- FIXED (mass-weighted).** star%accreted_mass_fraction
  is the mass accreted this step in g (wind/mdot.f90); dburn mixes it
  mass-weighted (X*M + X_acc*M_acc)/(M + M_acc), dburnm as
  (X + X_acc*f*M_acc)/(1 + f*M_acc). Active only in rotating runs with
  mass accretion.
- **dburn's single-shell accretion branch uses total_shell_mass unset
  -- FIXED (both integrators set it from the shell's mass).**
  (zero only through -finit-local-zero), so a lone radiative surface
  shell under accretion takes the accreted deuterium and ignores its
  own.
- dburn and dburnm are otherwise one algorithm (units, rate source and
  the 1e-11 / 1e-14 thresholds differ): merge once the two points above
  are settled.
- **Li/Be burning never sees the end-of-step model -- FIXED (flags
  swapped at the three lirate88 calls).** lirate88 was
  called with flag 1 at the start of the step (current model = the
  start-of-step model, stored in star%rate_*) and flag 2 at the end
  (recomputed from star%*_start into star%rate_*_start), so liburn's
  start/end rate interpolation uses start-of-step rates throughout.
  Inherited from F77 (the flags were swapped). Now flag 2 at the start
  (start-of-step model -> star%rate_*_start) and flag 1 at the end (the
  converged model, or the mid-sub-step model in the rotation loop ->
  star%rate_*). Structure unchanged; surface Li7 at 4.568 Gyr -0.088 dex
  (noGS_norot) and -0.113 dex (GS_rot rotating), Be9 < 0.001 dex. A
  calibration that tunes mixing to the solar Li abundance shifts.

### Information flow (yrec-next) -- results through arguments, stages in order

- **compute_zone_energy returns everything -- DONE (f1585425).** The
  neutrino loss, alpha-capture and He3 energies, fluxes and Be7 came
  back through star% scalars; per-shell reaction rates and branchings
  were stored and never read (removed).
- **mix -> core/step/composition_step.f90 -- DONE (86fd404d).**
- **The composition step is called by evolve_step, not by the
  structure solver -- DONE.** henyey_iterate (crrect) refreshed the
  surface boundary condition, re-ran the composition step at levels
  >= 3, then iterated. It is now two calls, refresh_surface_bc and
  henyey_iterate, and solve_level runs the three stages in that order
  (surfbc reads the surface composition, so it must stay first).
- **Open: the accretion term of the entropy equation covers every zone
  at levels 1-2.** henyey_coefficients adds the accreted layer's T and P
  changes (rot_scr%delta_log_*) to zones at or above
  envelope_zone_index. That index is the composition step's envelope
  edge at levels >= 3, but at levels 1-2 it was an unset local of
  crrect, 0 through -finit-local-zero, so every zone gets the term
  there. Preserved explicitly in solve_level. Only runs with
  use_mass_accretion and mass_accretion_rate > 0; levels 1-2 are
  predictors, so the converged model is affected only through the
  starting point of levels 3-4. Likely unintended; fixing it changes
  accretion runs.
- **core/step/convection_geometry.f90 -- DONE.** The CZ-base interpolation
  (three copies: observables_lib, turnover_timescale, burn_step's
  locate_cz_base), the pressure scale height behind the envelope
  overshoot depth (evolve_step, evolve_angular_momentum) and the
  overshoot-edge search (twice in locate_cz_base) are now
  cz_base_fraction, pressure_scale_height and first_zone_above;
  locate_cz_base moved there as locate_surface_mixing_base. Each caller
  keeps its own gradients (plain or rotation-corrected) and its own
  overshoot switch (envelope_overshoot_active, or lovstm in the rotation
  sub-steps -- the open gating question above). Byte-identical.
- **Open: the overshoot edge's zone width mixes two models.**
  locate_surface_mixing_base finds the end-of-step overshoot edge in the
  current radii but takes the width of the zone above it from the
  start-of-step radius (logR_start), so cz_base_frac mixes two models.
  Likely meant radius(cz_base_zone); preserved for now.
- **One trace-element integrator for deuterium and Li/Be -- DONE.**
  net/light_elements.f90: burn_exposure, surviving_fraction and
  accreted_survival (a species destroyed at a rate interpolated over the
  step, linearly in time for D or in log for Li/Be; closed form, plus a
  Simpson rule for steadily accreted gas). They replace dburn's
  sub-steps (which summed to the same integral without accretion) and
  liburn's sub-stepping plus rational extrapolation (which converged on
  it to 1e-6; checked against a replica of the old algorithm to 4e-8).
  average_li_be_over_cz no longer writes the adjusted CZ-base rate into
  star%rate_* (the rotation sub-steps carried it into the next
  sub-step's start rates) and lost its branch on star%log_rate_*_prev >
  0, which a natural log of a rate per second never is (the members are
  gone). A leftover write(*,...) that printed on every refinement under
  accretion is gone.
- **Open: surface Li7 has a noise floor of ~0.005 dex.** Found while
  measuring the integrator change: the deuterium closed form alone
  (mathematically the old sum, different only in roundoff) moves surface
  Li7 at 4.57 Gyr by +0.0053 dex in noGS_norot, because a marginal
  convective zone at the start-of-step CZ base flips in 46 of 348 models
  and the start-of-step CZ rate average uses that integer zone index (the
  end-of-step one is interpolated within its zone by cz_base_frac).
  Interpolating the start depth the same way would remove most of it.
- **rotmix burns on the rotation sub-step's own rates -- DONE.** The
  network burn in each rotation sub-step (mid_timestep_model -> rotmix)
  took its rates from rot_scr%reaction_rate_by_zone, a copy the last
  composition step left: the rates of the end-of-step structure from the
  previous iterate. Deuterium and Li/Be in the same sub-step already used
  the midpoint model; the network now does too (network_rates_of_model,
  passed to rotmix as an argument). The relay is gone, with rezone's
  regridding of it, which covered 7 of its 15 rows. Measured on
  GS_rot (the only shipped case with rotational mixing): |dlog L| 8.5e-5,
  |dlog Teff| 6e-6, |dX_c| 2.2e-4 (X_c lower), surface Li7 +1e-4 dex.
- **One name per index; Q-values in the reaction table -- DONE.** The
  rate vector's rows 1-13 had a second set of names (rr_*) for the r_*
  slots; the cross-section slots shared the s_ prefix of the network
  species (now xs_*). Each network-table row now carries its Q-value and
  chain, He3(He4,g) is split by Be7 branch, and nuclear_energy sums the
  energy and the chain shares from the table (reaction_energies) instead
  of a hand-written list of Q-values beside it. Adding a reaction is now
  its rate formula plus one table row. Byte-identical at print precision
  in solar noGS and the config matrix (roundoff only: f + (1-f) for the
  branches, and association order).

### Module naming (yrec-next, 2026-10) -- docs/GUIDELINES.md naming rules

- **Sweep -- DONE (cfd4b460, 166b30d2, 8d4bdcbb, 4766b8c0, bd754c07).**
  Every file is named after its module; _lib remains only on the domain
  facades (atm_lib, eos_lib, kap_lib, math_lib, net_lib, numerics_lib,
  wind_lib); every other module is a noun. Ten modules had carried _lib
  only because they contained a procedure of their own name; those
  procedures are now verbs (compute_nuclear_energy, advance_composition,
  advance_rotation_substep, apply_secular_transport,
  compute_temperature_gradients, run_microdiff, ...). In net/ every
  procedure follows the verb rule, and net_lib is a pure facade.
- **Open: the state modules.** star_info_lib, rotation_scratch_lib and
  point_scratch_lib keep their names until the Phase 5 state split
  (docs/TARGET_ARCHITECTURE.md section 4.2) replaces them.
- **Open: facade-only use for eos, kap, atm, numerics.** (Status
  2026-10-29: FACADE_ONLY_DOMAINS lists all seven star-free domains,
  and integration_parameters is gone since 97c00d5a; the text below is
  as written in 2026-10.) Only net/ is in
  check_boundaries' FACADE_ONLY_DOMAINS. The others are still entered
  through internal modules, mostly table state: read_controls and
  read_starting_model set yale_eos_tables, scv_eos_tables,
  opacity_tables and atm_tables directly; henyey_iteration and
  rebuild_envelope read yale_eos_tables (Debye-Hueckel); the envelope
  integrator and stitched_model use atm_tables, atmosphere_structure,
  envelope_structure and t_tau_relations; envelope_integration_kernel
  and read_controls use integration_parameters. Each needs its facade
  to expose an interface (setters for the table paths and flags, the
  structure profiles as types) before its domain joins the list.
- **Open: procedure names outside net/.** The verb rule was applied
  outside net/ only where a module clashed with its routine; the rest
  follow it as code is touched.

### net/ hygiene (yrec-next, 2026-10) -- the evaluation's six items

- **DONE, bitwise identical:** compute_nuclear_energy in explicit steps
  on one composition array (35e5f699); compute_reaction_rates in named
  steps (8f619439); the compiled-out Li/Be midpoint variant in its own
  file (9274053d); burn_step's five region loops as one burn_regions
  list, burn_step private, explicit imports (71a9bfb7).
- **DONE: test_net covers the solver, trace integrator, energy and
  table.** Closed forms for burn_implicitly, the integrator against
  closed forms and quadrature, eps at three reference points with
  derivative checks, and each table row's Q against AME2020 mass
  excesses.
- **FIXED (664c2fbc; see the 2026-10-11 follow-ups below): N14(a,g) ->
  O18 deposited too much.** q_mev = 5.815 implies a mean neutrino energy of 0.256 MeV.
  The F18 beta+ decay (endpoint 0.634 MeV, mean positron energy
  0.250 MeV) gives about 0.38 MeV, about 0.42 with the 3% electron
  capture branch. The inherited value looks like the mean positron
  energy was taken for the neutrino's. The fix is q_mev of about
  5.65 MeV; it changes helium-burning runs, so it needs a measured
  commit. Every other beta+ row (N13, O15, F17, Be7, B8) matches the
  literature to within a few keV.
- **Open: store_li_be_rates zeroes both rate sets past the first cool
  shell,** whichever set it is filling. Making each call clear only
  its own set changes surface Li slightly (measured earlier), so it
  needs a measured commit.

### net/ re-evaluation follow-ups (yrec-next, 2026-10-11)

All bitwise identical (test_net pins, gate1, config matrix, full run).
- **DONE: the mathematical constants are parameters** (529d3fb3). ln10,
  cc13, cc23, cpi and the rest were assigned by setups at run time, so
  net/ silently computed with zeros unless setups had run. Each
  compile-time value is bit-identical to the run-time one.
- **DONE: net_lib exports 35 names, not 62** (ad723f2a): only what the
  rest of the code uses. The individual rates, the table and the
  thermal-neutrino routines stay internal.
- **DONE: test_net pins the remaining pieces** (deaf3c55): the
  equilibrium burn, burn_deuterium_region (accretion included),
  compute_li_be_rates, the deuterium cap and the solar neutrino fluxes.
- **DONE: the equilibrium burn (eqburn) rewritten**, 201 lines to 90
  (6101fe95), and burn_deuterium_region in named steps (772c8220).
- **DONE: comments in current terms** (ff9b569b): no ALL-CAPS or
  "logic unchanged" boilerplate left in net/'s code, and one comparison
  style throughout.
- **DONE (ccd1918c, below): the level-1 equilibrium burn is dropped.**
- **Left as is: Timmes' sneut** (1,357 lines of external code with
  saved tables). It is why compute_nuclear_energy cannot be pure.
  (Pure since 78f58053, 2026-10-29: the saved tables are parameter
  arrays.)

### net/ to 9 (yrec-next, 2026-10-11)

- **DONE, bitwise identical: simeqc and burn_implicitly are pure**
  (ef4c50f1); in net/ only Timmes' neutrino routines and their callers
  are not (they are since 78f58053, 2026-10-29).
- **DONE, bitwise identical: net/ uses no star_info_lib** (d06db92b). The
  slot constants live in the constants-only modules composition_layout
  and neutrino_flux_layout (state/), which star_info_lib re-exports;
  check_boundaries forbids star_info_lib in net/. They cannot live in
  net/: numerics uses star_info_lib.
- **DONE, bitwise identical: each reaction's rate lives in its table
  row** (c2955bdd): rate law, S-factor fit, scale slot, screening
  constants, burning regime. Adding a reaction is one row plus its slot
  name. Two departures from TARGET_ARCHITECTURE 4.5(b): rate laws are
  integers dispatched in compute_screened_rate, not procedure pointers
  (the table stays a compile-time constant and everything stays pure);
  and each law applies the screening exponent inside its own formula,
  not the driver after a bare fit, because exp(a + b) and
  exp(a)*exp(b) differ in the last bit. Rate scale factors are still
  referred to by slot (scale_slot), not by name (by name since
  77dfc919, 2026-10-29: each row's scales are at its own reaction
  slot).
- **DONE, measured: the level-1 equilibrium burn is dropped** (ccd1918c).
  The predictor pass ran the implicit network and then overwrote X (and
  Z, C12, O16 in helium burning) with an explicit equilibrium-He3/CN
  estimate; both passes now burn the same way and the second copy of
  the burning equations is gone. The hydrogen-shell timestep limiter
  reads dX/dt from the table (compute_abundance_rates).
  logT_cutoffs(2) and (3) are now read and unused.
  Drift, full pin set: pre-MS, ZAMS and config-matrix cases
  byte-identical; the 16 solar-age cases move, most by 1e-6 to 2e-5 in
  log L at the solar age, but up to |dlog L| 3.8e-3 on the early main
  sequence (Test_solarGS98_dif) and |dX_c| 1.1e-3 at the solar age
  (Test_solar_GS_norot_gray). Root cause: time discretization. The
  composition step runs once per Newton level, on a structure the
  predictor influenced, so the result carries an operator-splitting
  error of the timestep. With Test_solarGS98_dif at a 4x finer timestep
  the two predictors agree to 8.5e-5 in log L (2.5e-7 at the solar age),
  while both standard-timestep runs sit 3.6-3.7e-3 from that solution,
  equally far: the predictor moves the answer inside the standard
  timestep's own error. (A 2% change in the convergence tolerance moves
  it by 1e-6, so this is not solver noise.) Reseeded: the 16 legacy logs
  and the two MESA baselines (same file sets); backup in the session's
  backup_pre_series4.
- **DONE, measured: N14(a,g) deposits 5.646 MeV** (664c2fbc; was 5.815,
  which implied the mean positron energy as the neutrino loss). No
  pinned case reaches helium burning; a horizontal-branch run moves by
  |dlog L| <= 4.2e-8 (its seed has already burned its N14).


### Phase 3: kap/ (yrec-next, 2026-10-11)

All bitwise identical (test_kap, gate1, test_net, config matrix).
- **DONE: kap is reached only through kap_lib** (ce886f00): the outside
  uses of opacity_tables were dead, as was the llaol switch.
- **DONE: test_kap pins every family** (a19e49a0): the production
  combination, OPAL92 and LAOL89 with two Z tables, Kurucz 1990 and
  Alexander 1994/95 (no deck uses the legacy families).
- **DONE: kap reads no star** (6e979f91): kap_settings_t, filled by the
  star layer's opacity_settings (setup/microphysics_settings.f90) and
  installed by kap_init; the readers open their files with newunit;
  the conductive table is read at kap_init from
  star%job%conductive_table_path (read_controls no longer opens it);
  the surface-mixture slots move to composition_layout.
- **DONE: one module per family** (fe63ec5c): 35 files, 32 of them bare
  routines, become 9 modules; routines named as verbs (legacy names in
  TARGET_ARCHITECTURE appendix A); the Kurucz lookup's alternate return
  becomes a logical in_table.
- **DONE: hygiene and comments** (f7107efd, f3b505c6): real(dp), private
  by default, explicit imports, one comparison style, comments in
  current terms; kap/README.md.
- **Open:** the readers still branch on end= labels; util's xrng4 is a
  bare routine kap calls without an interface; the unit-number controls
  (opal95_table_unit, icondopacp, ...) are still set by read_controls
  but read by nothing.

### Phase 3: eos/ (yrec-next, 2026-10-11)

All bitwise identical (test_eos, test_kap, test_atm, test_net, gate1,
config matrix) except where noted.
- **DONE: each OPAL vintage has its own table-loaded flag** (85f74924):
  the three readers shared one flag, so loading one vintage kept the
  next from reading its own table in the same process (a latent bug;
  no run loads two vintages).
- **DONE: test_eos pins every configuration** (8cd2d244): OPAL 1995 and
  2001, SCV only, Yale only, numerical derivatives, diffusion of Z and
  the Debye-Hueckel correction, besides the production SCV + OPAL 2006.
- **DONE: eos reads no star and is reached only through eos_lib**
  (c00a9a20): eos_settings_t, filled by the star layer's
  equation_of_state_settings and installed by eos_init; the OPAL table
  path comes from star%job%opal_eos_table_path and eos opens it; the
  Debye-Hueckel composition and the SCV envelope table are set through
  eos_lib.
- **DONE: alternate returns become in_table** (afbec93f): esac* and
  oeqos* report an OPAL miss through a logical, with identical control
  flow.
- **DONE: one module per family** (1ab2dc75): 42 bare routines become
  opal95_eos, opal2001_eos, opal2006_eos, mhd_eos, scv_eos, yale_eos and
  blended_eos; OPAL's and Daeppen's routines keep their names (legacy
  names in TARGET_ARCHITECTURE appendix A).
- **DONE: hygiene** (02edee58, a15d1dc7): real(dp), private by default,
  explicit imports, one comparison style, assumed-shape arrays (not in
  the vendored MHD readers), constants as parameters, seven procedures
  pure.
- **DONE: eos_init is a list of steps** (5cbfb7a0): each family reads its
  own tables; the Debye-Hueckel constants are parameters (bit-identical
  to the run-time values); a missing Fermi or SCV table returns ierr
  instead of crashing.
- **DONE: comments in current terms** (778208e2); eos/README.md.
- **Open:** the MHD path has no pin (its tables are not shipped).
  rhoofp01's density-edge table differs from opal01_density_index_edge
  at rows 78-80 (113 against 115), as in OPAL's code; making them one
  table would be a measured change. The three OPAL vintages are near
  copies of each other (esac*, rhoofp*, radsub*). beta14 still means
  different things on different paths (TARGET_ARCHITECTURE appendix B).
  The unit-number controls iopale, fermi_unit, iscvh/iscvhe/iscvz and
  unit_zams_*/unit_centre* are set by read_controls but read by nothing.
- **Fixed after the full battery** (841ff500): the module merge had
  renamed eqstat2 inside a message string too, so the 0.03 Msun
  dbl-to-ZAMS pin's log read "blend_eos: BETA is negative" where it
  read "EQSTAT2: ..." (numbers identical). The message keeps its text;
  merge_family no longer renames inside strings. The kap merge
  (fe63ec5c) renamed seven 2026-written error messages the same way;
  no pin reaches them and they now name the routine that prints them.

### Phase 3: atm/ (yrec-next, 2026-10-11)

All bitwise identical (test_atm, test_eos, test_kap, test_net, gate1,
config matrix, the testsuite's Allard and Kurucz solar models).
- **DONE: test_atm pins every path** (ecea0ce0): besides the Kurucz,
  Castelli/Kurucz and Allard NextGen 2 lookups, the Allard tables at
  tau = 100, the old NextGen 1 layout, a lookup outside the Allard
  tables, and the three T-tau relations. read_allard_tables now closes
  its unit, so the tables can be read twice in one process.
- **DONE: atm reads no star and is reached only through atm_lib**
  (85575e1b): atm_settings_t, built by the star layer's
  atmosphere_settings; tables opened with newunit; atm_get_surface_pt
  dispatches the tables itself and returns P and T as arguments, with
  the caller's current atm_choice (surfbc switches it); the atmosphere
  integration point and the envelope/atmosphere structure profiles
  move to core, which computes them.
- **DONE: one module per table family** (d9909c9f): kurucz_atmosphere
  and allard_atmosphere replace five bare-procedure files; atm_init is
  three steps.
- **DONE: hygiene and comments** (249b7a03, 3e5479cb): real(dp),
  private by default, explicit imports, one comparison style, eight
  procedures pure, comments in current terms; atm/README.md.
- **Open:** read_allard_tables's check that an old NextGen 1 file cannot
  give the tau = 100 point tests latmtptau100, a never-set misspelling
  of the allard_use_tau100 flag, so that misconfiguration reads unset
  tau = 100 tables instead of stopping. Fixing it changes no pinned
  output but is a behaviour change. The accretion branch of
  change_mass (formerly wind/massloss) reads the shared atmosphere point
  after an Allard lookup that may have failed (kept as it was). The atm_table_file_unit
  and allard_table_unit controls are no longer read.

### eos/ follow-up: interface, derivatives, solvers, OPAL (yrec-next, 2026-10-11)

From an evaluation of eos/ (readability ~6, maintainability ~5.5,
hygiene ~8.5 out of 10). All bitwise identical (domain tests, gate1,
config matrix, the testsuite's Allard and Kurucz solar models, the
Debye-Hueckel 0.03 Msun case) except where noted.
- **DONE: the derivative contract** (1c177c25, b96442be): test_eos
  checks every derivative against centred finite differences at six
  points per configuration; eos/README.md states what each one is. It
  found:
  - dlnrho_dlnp_dt was d(dln rho/dln T)/dln P (renamed dlnrho_dlnt_dp);
  - with an OPAL table, the six second derivatives are the Yale
    solution's while density and the first derivatives are OPAL's
    (off by 8-180%, one sign flip);
  - in the Saha-to-fully-ionized ramp the blend omits the ramp weight's
    T dependence (eqstat2.f's WT1 line is commented out) and blends
    QDTT/QDTP as fully + w*saha (it never forms the differences);
  - of eos_eval's inout arguments only the density (the fully ionized
    iteration's starting guess) and saha_state carry information in.
  use_numerical_derivatives agrees with finite differences to ~1e-6.
  All of these are recorded, not fixed: correcting them is a physics
  change (they feed the Newton Jacobian, so they should affect
  convergence rather than converged models, which is untested).
- **DONE: eos_result** (5d39f836): one named record replaces the 24
  slots and the 31-argument eos_eval. The families take it as r, so
  20-30 argument lists become 3-9. log T and log P are intent(in).
- **DONE: the Yale solvers as named steps** (cdda0898, 9e0169dd):
  compute_saha_eos (447 lines) and compute_fully_ionized_eos (288)
  become drivers over 4 and 6 steps, with records for the state the
  steps share. The Saha solver's hidden SAVE guesses are now module
  variables. Prather's symbols stay, with a notation key, rather than
  long names that would bloat the algebra.
- **DONE: the OPAL fold** (522fdf58): one compute_opal_eos for all
  vintages, one loop in blend_eos, and the 2001/2006 code merged into
  opal_eos with esac/readcoeos/t6rinteos folded, branching on the
  vintage only where they differ. The OPAL code drops from 3238 to
  2688 lines; rhoofp01/06 and radsub01/06 differ in algorithm and stay
  two, and the 1995 code keeps its own state layout.
- **Open:** the derivative findings above; the MHD path is still
  unpinned (it compiles against the new record); beta14's two
  definitions.

### wind/ dissolved: braking into rotation, mass change into core (yrec-next, 2026-10-11)

From an evaluation of wind/ (~4.5/10). Four files repeated one
sub-stepped loop and two solid-body forms. The structure factor
travelled through star%job. massloss and mdot took 21 and 25 star
fields as arguments, and the two halves of the domain shared nothing.
All bitwise identical: test_braking, domain tests, gate1, the config
matrix (cm_rot_kawaler and the two mass-loss decks among them), two
accreting variants of cm_massloss_const against the previous commit,
and the full battery.
- **DONE: test_braking pins every braking path** (8ac25953). It covers
  the convection-zone form, with:
  - the Matt law with and without Rossby or B-field scaling, saturated
    and sub-stepped;
  - the Kawaler law with and without tau_cz scaling, inactive, and on
    its error path.

  It also covers the solid-body form for both laws.
- **DONE: magnetic_braking** (be9f634d). rotation/magnetic_braking.f90
  replaces matt_wind, kawaler_wind, wind_spindown_matt, wind_spindown,
  matt_structure_factor and wind_lib with two entry points:
  - brake_convection_zone, one loop over a law-dependent domega;
  - compute_solid_body_braking, where the Kawaler form is the Matt one
    with f_cen = 1.

  The structure factor is a function result, so the
  star%job%structfactor control is gone. log10_radius_from_l_teff
  moves to physical_constants and also replaces observables's
  log_r_surface_cm, which was the same expression.
- **DONE: mass_change** (5e46b707). core/step/mass_change.f90's change_mass
  and apply_mass_change (formerly massloss and mdot) read the model
  from star, so they take 3 and 8 arguments. mdot's mu output, always
  overwritten before use, is a local.
- **DONE: the wind domain is gone.** wind/ and libwind.a are removed,
  and check_boundaries drops the domain.
- **Open (kept as they were, now each in one place):**
  - the solid-body Matt f_sun lacks 1/G (the gl item above);
  - the Matt law ignores wind_loss_active;
  - with B-field scaling but no Rossby scaling, the Matt convection-zone
    iteration scales Omega by a tau_cz that is only set with Rossby
    scaling. The torque is then zero, and test_braking shows the step
    removes no angular momentum. The 2026-08-31 sweep noted the read;
    the original read a stale SAVEd value instead. The solid-body form
    ignores B-field scaling without Rossby scaling;
  - change_mass resets the envelope entropy to zero, so accretion's
    R/P/T correction is a zero shift. It still rounds the radii through
    exp10/log10. The convection zone's new mean thermal energy is
    computed and never used.

### numerics/ revamp (yrec-next, 2026-10-11)

From an evaluation of numerics/ (~4.5/10). It was one 3012-line
numerics_lib with 38 procedures. About 1,300 of its lines were comments,
largely one repeated header. Four procedures carried one caller's
physics, eight read star_info_lib's json, six kept blanket SAVEs, and
there were duplicate families: three tridiagonal solvers and five
spline pairs. It had no test of its own. All bitwise identical (test_numerics,
test_rotation, test_eos, the domain tests, gate1, the config matrix,
the testsuite's Allard and Kurucz solar models, the Debye-Hueckel 0.03
Msun case, and the full battery) except three error messages.
- **DONE: test_numerics pins every routine** (94536a04), error paths
  included. Under gcov only unreachable lines go unexecuted.
- **DONE: single-caller routines move to their domains** (18f2b895,
  3d809025):
  - trapzd and qgauss (with the integrand func) become
    rotation/shape/shell_quadrature;
  - tridia becomes secular_transport's solve_angular_momentum_transport;
  - intpt becomes mhd_eos's interpolate_mhd_table.

  test_braking became test_rotation to pin the rotation pair, and
  test_eos pins intpt.
- **DONE: five topic modules behind numerics_lib** (b71ea887):
  cubic_spline, quadratic_spline, table_interpolation, linear_solvers,
  bulirsch_stoer. numerics is star-free and facade-only.
  - kspline/splinc became natural_spline; ksplint/splintd2 fold into
    splint; ctridi and tridiag_gs became solve_tridiagonal.
  - The saved state of bsstep, ratext, lir and intpol went (ratext's
    tableau is a record bsstep holds for one step), along with lir's
    unused resume-the-search option.
  - integration_parameters moved to state/.
- **DONE: bsstep integrates an ode_system** (278fa171). The envelope's
  twelve pass-through values became the kernel's envelope_equations,
  per the "physics never inside the numerics" principle (TARGET_ARCHITECTURE
  P2).
- **DONE: hygiene and a README** (8193b431).
- **Found, kept as they were (recorded, not fixed):**
  - bsstep's step-underflow test compares h with the hydrogen fraction
    (the F77 X+H.EQ.X, where X was the hydrogen-fraction argument, not
    the independent variable X0). bsstep now takes the reference as an
    argument, and the envelope passes X.
  - A step's error is the largest non-NaN scaled error, NaN only when
    every component is. This came from how the compiler handled
    max(0, |e|) and a NaN-skipping max, which Fortran leaves
    processor-dependent. Resizing bsstep's arrays changed the
    instructions chosen and turned the 0.03 Msun deck's clean stop
    (its atmosphere integration starts from NaN) into STOP 1. That
    behaviour is now written out explicitly and pinned. The gate1
    and config-matrix decks did not catch it; the Debye-Hueckel deck
    did.
  - The quadratic spline's knot cases 3 and 4, and its sign-flip
    guards, cannot be reached through osplin: slopes keeps every slope
    below twice both adjacent secants.
- **Open:** the Numerical Recipes routines (cspline, splint, polint,
  ludcmp/lubksb, hunt, locate, bsstep with mmid and ratext) remain a
  licensing liability
  (MODERNIZATION_CHANGES). Replacing them changes results at round-off.

### mixing/ split: convection/ and rotation's burn-settle-mix (yrec-next, 2026-10-11)

mixing/ held two unrelated things: convection physics, used by core
and rotation, and the burn-settle-mix sub-stepping, used only by
rotation. All bitwise identical (test_convection, test_rotation,
test_numerics, the domain tests, gate1, the config matrix with its
overshoot and semiconvection decks, the Allard, Kurucz and
Debye-Hueckel decks, and the full battery).
- **DONE: test_convection pins the convection physics** (1639147b):
  the gradient on every path, the zones with every kind of overshoot
  (Yale EOS scale heights), homogenization.
- **DONE: burn-settle-mix moves to rotation** (960e717e): rotmix (live)
  and the compiled-out Bulirsch-Stoer path (burn_settle_mix,
  burn_mix_extrapolated, burn_li_be_midpoint).
- **DONE: convection/ is a star-free domain behind convection_lib**
  (238a7be2):
  - modules temperature_gradients, convection_zones, semiconvection
    and convection_configuration;
  - settings installed every run card (mixing length, Lsun) plus a
    setter for zenvm;
  - the overshoot region (star%iov1/iov2) returned through arguments,
    and the gradient takes in_overshoot_region;
  - star%iovim and the always-zero star%alfmlt/phmlt/cmxmlt removed;
  - core and rotation now call convection through an explicit
    interface, which exposed am_convective_regions' intent(in)
    convective flag (element num_zones+1 is written).
- **Open:**
  - DONE 2026-10-11 (the author chose to comment it out): the
    Bulirsch-Stoer burn-plus-mix path (about 770 lines in rotation/) is
    now plain comments, its two call sites in evolve_angular_momentum
    commented out and their always-false switch folded away. Its
    burn_li_be_midpoint was that path's own variant of core's
    burn_li_be, not a live duplicate.
  - rotmix's name is misleading: it burns and settles helium, it does
    not mix. To be addressed in the rotation series.
  - Semiconvection has no standalone test (it needs the opacity
    tables); cm_semiconvection pins it.

### io/ revamp: the controls reader, the model files, the pulsation writers (yrec-next, 2026-10-11)

io/ held a 2,477-line controls reader (one routine, parmin's descendant)
and eight external procedures called through implicit interfaces. All
bitwise identical: test_controls compares every deck's star%ctrl and
star%job, and the quick check, test_reentry, test_inlist_convert,
test_mesa_output and the full battery compare the runs.
- **DONE: test_controls** (89b04418): reads all 313 decks in examples/
  and testsuite/ through read_controls and dumps every element of
  star%ctrl and star%job (generated dump code) against
  io/test/expected_controls/.
- **DONE: the generators run and are checked** (2b32c077, 384b3d73):
  gen_inlist had broken on a moved file; all five generators take
  --check, and test_generated_files runs them. sstandard is documented
  as a legacy-only control in the registry instead of by hand.
- **DONE: the reader split by job** (fe347c28, b10d37ff): controls_reader
  (read_controls, read_input), legacy_controls (the namelist variables),
  run_files, derived_controls and run_cards. The settling grid mapping
  moved to rotation/microdiff (4ce14a01).
- **DONE: the defaults in one routine** (a56aab9f): every legacy variable
  was zeroed and then defaulted on each read; set_legacy_defaults now
  assigns each once. The three binary-model units, never read, are gone.
- **DONE: derivations on the canonical names** (9364934f): the derived
  controls were computed on the legacy names after they had been copied,
  then copied again in three re-sync blocks (21 controls). They now work
  on the canonical names, and read_input reads as its sequence of steps.
- **DONE: model_files and pulsation_output modules** (f7d50b4c): the .mod
  writer and reader, the YREC7 and MODEL2 readers (26 arguments each),
  and the FGONG, GSM and two GYRE writers now have explicit interfaces.
- **DONE: controls.defaults documents the real defaults** (c3f40675): 188
  documented defaults were wrong (mostly "0.0 (implicit)"); gen_inlist
  now takes them from the code, and 69 registry descriptions that were
  migration notes describe their controls.
- **DONE: io/README.md and comments** (7ebca0e3, 064ca119): the COMMON-block
  history in controls_buffer and legacy_controls (about 900 lines) gave
  way to short headings.
- **DONE: I/O errors stop the run cleanly** (2026-10-11, after the
  entries above): every open and namelist read in io/ and every output
  writer's open checks its status; check_io reports "cannot <action>
  <path>: <reason>" and returns ierr. A missing starting model or a
  misspelled control used to end in a Fortran runtime error.
- **Open:**
  - DONE 2026-10-11 (the author chose "code is the source"): each
    control's default is written once. 50 legacy defaults that repeated
    their canonical twin's controls_buffer initializer now mirror it
    (alphac = overshoot_alpha_core), the registry's default column is
    gone, and gen_inlist stops if a live control has no default in the
    code. The legacy namelists are frozen: new controls are new-format
    only. Doing this found tscut's lost F77 default of 6.0 (restored).
    (Not yet true for 180 legacy names, whose literals in
    set_legacy_defaults disagreed with their canonical initializers;
    6085de4c, 2026-10-29, made the initializers the one source.)
  - read_yrec7 and read_model2 return 26 header values, most of which the
    caller discards; a header record would shorten both.
  - monte_carlo_file1_path and monte_carlo_file2_path are NUL-filled when
    a deck leaves them unset (pinned by test_controls).
  - Ten decks cannot start on their own: six class exercises name a
    relative startmodels/ directory, and four testsuite decks start from
    another run's output. They now stop with "cannot open the starting
    model"; test_controls records that.

### rotation/ series: element diffusion split out, modules, README (yrec-next, 2026-10-12)

All bitwise identical (the new tests, test_reentry, the quick check,
the solar GS decks rotating and not, the three testsuite
rotating-mixing decks, and the full battery).
- **DONE: pins first** (491adc36, 94c33ff4): test_diffusion (now
  test_element_diffusion) runs both settling routines on a synthetic
  solar-like model. No deck used the original gravitational_settling,
  which a deck gets whenever it turns on settling without LNEWDIF.
  test_rotation gains the rotating structure (moments of inertia,
  omega from j, the shape factors, viscos).
- **DONE: element_diffusion/ is its own domain** (f2d05f1e, bc215a7b,
  fb0aa9d6, ce7648e2): settling is not rotation physics (core settles
  non-rotating models through it). Topic modules behind
  element_diffusion_lib, whose diffuse_elements makes the
  original-or-newer choice both callers used to make; its own scratch,
  reset through the facade; a README.
- **DONE: rotation's procedures are module procedures** (6c08a206,
  b4ea1bd1, d380df08, 5e650c72): its 41 external-procedure files became
  topic modules (rotating_shape, moments_of_inertia, omega_profile,
  am_regions, viscosity, and in seculr am_transport_coefficients,
  circulation, rotational_mixing, am_grid, transport_checks; the driver
  is evolve_rotation). rotmix, which burned and settled but did not mix,
  is burn_and_settle in mid_timestep_model. Rotation stays facade-free
  (the author's decision, confirmed 2026-10-12); check_boundaries lists
  the six procedures core and setup call.
- **DONE: evolve_angular_momentum reads as its steps** (e7d66a5b).
- **DONE: comments and README** (rotation/README.md): the migration
  history, F77 banners and stale file references gave way to
  descriptions; rotation_scratch_lib's 32 "former common" groups have
  topical headings.
- Found: use_thoul_fit and use_thoul_diffusion act only in the original
  settling routines. rotmix declared the enclosed mass intent(in)
  although settling converts it in place (now inout).

Round two (2026-10-12), after a rating of rotation/ at 5.5/10; all
bitwise identical (the rotating solar deck, the three testsuite
rotating-mixing decks, the quick check with the config matrix at every
commit, and the full battery at 43f852c1):
- **DONE: the secular transport is pinned** (19f8ba45): test_rotation
  runs am_transport_grid -> am_diffusion_coeffs -> the solve -> the map
  back -> the composition mixing on a synthetic model, through the
  first and second iterations, the secular and dynamical shear bands, a
  mu-gradient band, the GSF modes, lvfc, constant background diffusion
  and no transport in the core.
- **DONE: the long routines read as their steps** (db16932a, f8efd576,
  2bcdc80a, 4f939ca0): apply_secular_transport, circulation_velocities,
  advance_rotation_substep, burn_and_settle, check_angular_momentum and
  am_transport_grid are hosts whose internal procedures are named steps;
  no two rotation procedures share a name.
- **DONE: the scratch holds only what rotation shares** (27d6e05a,
  58702a54): what nobody read went (the max_domega_dr loop, the
  always-zero mu-gradient velocity, fact1-6, fgsfj, the potential
  terms, two derivative copies), core's mass-change state moved to
  star, and the spline tables are the routines' locals. Names say what
  they hold: the interface structure (density_interface, ... for dm,
  pm, tm, hs3, delami, delmi, epsilm, qdtmi), mid_scr with the star%
  names it interpolates, the chi grid (eq_chi, num_eq_points).
- **DONE: layering** (65768aed): rotation_stability_setup left setup
  for rotation/seculr/stability_coefficients, its header giving the
  coefficients' equations; core imports rotation's scratch only where
  it uses it. shell_physics stays in core: every run calls it at the
  start, and rezone for envelope overshoot.
- **DONE: a settings record** (43f852c1): rotation_settings
  (rotation_configuration) holds the 54 controls rotation reads, named
  for what they do (walpcz is cz_rotation_law_exponent, constfactor
  wind_constant, ...); setup/rotation_settings_setup builds it and
  run_yrec installs it each run card. check_boundaries'
  CONTROLS_FREE_DIRS keeps rotation off star%ctrl/star%job.
- **DONE: comments** (7dabe3f8): the ALL-CAPS F77 comments became
  descriptions with equations and references; comments that
  contradicted the code now describe it (the mixing grid is in chi;
  convergence is on j; di_domega is (omega/2) dI/d omega; the "_mid"
  structure is at the sub-step's end, not its middle).
- Found:
  - mu_gradient_scale (fmu) has had no effect since the always-zero
    mu-gradient velocity went; the record leaves it out.
  - The reversal safeguard in check_angular_momentum enters only where
    omega exceeds 1 rad/s (as F77 checkj.f), so it practically never
    acts. checkj.f's authors made the first negative j stop the run
    (IREDO.GT.0 replacing .GT.3); the message's "3 attempts" is stale.
- **Decided (2026-10-12):**
  - Kept and documented (docs/HOW_YREC_WORKS.md section 4.6, the
    README, inhibit_by_mu_gradient): the mu gradient is counted twice
    in the Eddington-Sweet inhibition (inhibit_by_mu_gradient adds
    d ln mu/d ln P to a margin that already includes it), so mu
    gradients cut U_ES by at most a half. Inherited from F77 (setupv.f,
    vcirc.f); changing it would move every rotating model.
  - Fixed: apply_secular_transport returned without setting redo_flag
    when the model was stable on the first iteration, so after a redo
    the stale flag would rewind the clock again; it now starts false.
  - Commented out, like the Bulirsch-Stoer path: the advection-diffusion
    transport (LDIFAD: am_advection_diffusion_coeffs, banded_solver, the
    advection branches of the live routines, the theta term, some 40
    scratch arrays and henyey's derivative copies; a rotating run that
    sets the flag still stops), and mass_change's rotation-scaled solar
    wind (and its four star members).
- **Open (behaviour; each needs a decision and a measured commit):**
  - At a transport region's ends, the equal-grid zone sums take in a
    convection zone below (above) whenever the neighbouring zone is
    convective, but the change goes back into it only when the end
    zone itself is (am_grid and rotational_mixing alike).
  - check_composition checks species 5-11 on the last iteration,
    before secular_transport mixes them; 12-15 are never checked.
  - max_delta_j_by_iter has 16 slots; nothing caps max_diffusion_iters.
- **Open (structure):**
  - Settling converts star%m to Bahcall & Loeb units and back in
    place; the round trip can change its last bit every step. Working
    on a copy would fix it but changes results at round-off; it needs
    a measured commit.
  - Rotation calls core's shell_physics and the light-element burn;
    core's convection_geometry reads mid_scr.
  - burn_and_settle and core's composition_step each size the settling
    sub-steps from the helium settling timescale with near-identical
    code; element_diffusion could own it.

### core/ Stage A: modules, layout, named steps, comments (yrec-next, 2026-10-12)

Readability without changing results: every commit byte-identical on the
quick check with the config matrix, the rotating solar and testsuite
decks, and three "core extras" no pinned deck covers (a TAHB deck that
rescales its model and burns helium, a solar deck writing the
neutrino-flux table, the accreting cm_massloss_const variant), compared
with a frozen 66d7788f binary; the step splits also on pulse output
(GYRE and FGONG, 84 files).
- **DONE: yrec_libs links again** (23779e3b): the per-domain archives
  missed math/ and element_diffusion/.
- **DONE: util/ dissolved** (91011457): core/step/timestep_control,
  core/run/yrec_version; xrng4 into kap's LAOL89, zero into the MHD EOS.
- **DONE: every procedure a module procedure** (8c9a69d3): the 23
  external procedures of core/ and setup/ in 16 topic modules; noun
  routine names became verbs (compute_surface_boundary for surfbc,
  compute_henyey_coefficients, prepare_run, set_up_constants_and_tables,
  ...). The compiler now checks every call; none needed a change.
- **DONE: core/ by stage** (f7789ad4): run/, step/, structure/,
  envelope/, model/ (with setup's rescaling, rezoning, shell
  boundaries), diagnostics/; setup/ keeps the settings builders.
- **DONE: the long routines in named steps** (442ea9b2, fef25296,
  aab69ed9, 5d7fbd0e, 102e7713): compute_henyey_coefficients,
  henyey_iterate, compute_surface_boundary, the envelope integrations,
  rebuild_envelope, read_starting_model's steps, rescale_model,
  map_user_inputs, rezone's steps, compute_timestep, the mass change,
  the turnover walker, the pulse points.
- **DONE: the math-lib contract per scope** (2ae7ee15): check_boundaries
  walks the scopes; refresh_surface_bc lacked math_lib (crmath builds
  only).
- **DONE: comments and core/README.md** (9337e932): 688 ALL-CAPS lines
  became descriptions with equations; henyey_solve's rows derived.
- **FIXED (behaviour, one path): the shallower envelope refit**
  (4f34e95f): read_starting_model's log R added Lsun in erg/s for log
  Lsun since the 2026-08 modernization (6e31b68d), so moving the
  fitting point outward failed in the envelope integration. No pinned
  deck takes the path.
- **Spot-checked (2026-10-12): all real.** The 0.03 Msun ZAMS pin
  (run_from_dbl_to_zams Test_m0030_feh0_noGS_ZAMS) was re-pinned for the
  refit fix: it used to stop in the envelope integration, now it refits
  and stops at its first evolution model (MODEL FAILED TO CONVERGE,
  OPAL06 density inversion in the degenerate interior) -- as it does
  with the refit switched off; this deck cannot evolve in the current
  code, a separate question.
  - **FIXED, group A (0d66503f): guards** that act only where the code
    misbehaved (overruns in rezone's X/Z insertion,
    extend_core_toward_center, convection_geometry, the turnover
    stencils; divisions by zero for a zero mass-change rate and a zero
    settling sub-step count; rescale_z = 0 zeroing Z; the senv = -1e-12
    reversal; the Allard-fallback integrand; the envelope shifts'
    convective flag; the burn locals below the cutoff in the pulse
    output; the mass change's debug print; the calibration's 14-cycle
    message).
  - **FIXED, group B (bc874267..3d065901, one measured commit each;
    all inherited, the HR-motion read from the 2025 ENTIME fix):**
    update_total_mass's m(0) (dm(1) from the center, as the rezone);
    the eps_grav partials' 1/chi (Newton path; every deck moves at the
    solver tolerance, |dlog L| peaks 1.5-1.8e-4 near the ZAMS as a
    ~1e4 yr timing shift, finals within 2e-5); the end-age floor, 1e-3
    of the age (only final dt's change); the HR-motion limit's start L,
    now star%surface_luminosity_lsun_start (nz_start goes; one
    testsuite deck moves within its own timestep sensitivity); the
    atmosphere's first slope with f_P (rotating decks, <= 1.9e-4 near
    the ZAMS); each triangle vertex's own Teff (spot decks, <= 3.2e-5).
  - **FIXED, group C (91a6b567..67d711a7; all inherited, each shown on
    a scratch deck, the 21 decks of the check byte-identical):** the
    wind removes dJ once (it removed dJ + <j>|dM|); the convection
    zone's law is omega ~ r^walpcz in mass_change and envelope_rebuild
    too; the Z ramp's CNO factor Z_new/Z_old (CNO/Z was 0.23 in the
    ramped core against 0.69); a reused model is rescaled once, and the
    relative species target no longer overwrites its control (Li7/X
    was 1e-12 after a second card, now 10^(value-12)); switching back
    from gray to a table starts a new triangle.
  - **Found while checking, open:** with -fcheck=all the decks tried
    (cm_base and the two mass-loss decks) stop at the start, in
    stitched_model's atmosphere_radius_cm, which reads
    env_log10_radius(num_env_points) with num_env_points = 0 when the
    starting model's stitched model is built (before any envelope is
    stored): the first stitched model's atmosphere radii come from
    memory outside the array. Output only. Also: three testsuite
    Allard decks (0.03, 0.05 and 0.1 Msun) stop with MODEL FAILED TO
    CONVERGE after 50-100 models, the same before this series.
- **The 0.03 Msun deck (e16b304a): a stop it can reach, and why it
  fails.** Its central-H stop never triggers in a brown dwarf, and
  every YREC version fails to converge (a Newton two-cycle) later in
  the contraction: the F77 original after 210 models (0.047 Gyr),
  370ae06d..b4673635 after 179, 0304e269 after 83, today after
  77. The unified deck's coarser mesh, tolerances, iteration levels,
  timestep limits and triangle made it fail on model 1; the deck now
  has the original's values and stops at 50 models (re-pinned).
  - Drift against the F77 original at equal age (original settings):
    dlog L -2.3e-3 at 0.006 Gyr, -2.0e-2 at 0.008, -8.4e-2 at 0.012;
    up to b4673635 it was below 2e-7. Bisected to Batch 2 (0304e269),
    and within it to item 2, the CF88 d(p,g)3He rate's T9**(-2/3)
    (was +2/3, about 5000x slower at 1e6 K): restoring the old power
    alone brings the drift to -9.4e-5. Deliberate physics; no fix.
  - The failure point is not one bug: on b4673635, Batch 2's eqscve
    blend fix alone moves it from 179 to 62 models (drift 1e-4),
    its rhoofp06 fix alone to 210 (drift 1e-9), the deuterium rate
    alone to 83; the other files leave it at 179. The solve is
    marginal here, and small changes move where it breaks. Open if
    the deck should evolve further: the two-cycle itself.

### core/ Stage B: kernels with unit tests, the star as an argument (yrec-next, 2026-10-13)

Every commit byte-identical on the check's 21 decks and the core extras
(and, where they apply, the Kurucz, Allard and Debye-Hueckel decks);
the full battery runs on the series' end.

- **DONE: the rule and the ratchet** (49cea924). GUIDELINES' star rule
  6: equations and algorithms live in star-blind kernels with strict
  intent and a unit test (STAR_BLIND_FILES); drivers gather from the
  star and store back. check_boundaries' STAR_ACCESS_BUDGET may only go
  down per directory.
- **DONE: four kernels with tests.**
  - The structure equations at one shell, `evaluate_shell_equations`
    (14c2215a): `test_structure` checks every partial against finite
    differences on a TAMS model (OPAL, tolerance 3e-3) and a birthline
    model (Yale EOS, 1e-4); with 18869257's 1/chi removed it fails.
  - The mass change's masses, wind angular momentum and accreted
    composition (a7349330): `test_mass_change`, conservation to
    round-off; fails with bc874267 or 91a6b567 reverted.
  - The rezone's regridding (2f0eb509): `test_regrid` pins the
    regridding's non-conservation (60-point model: J -1.3e-3, X +1.7e-4,
    Y -3.3e-4, falling as h^2).
  - The envelope triangle's geometry and fit (7d904caf):
    `test_triangle_fit`, exact planes, the inside test, the reflections.
  - Also star-free now: the envelope integrands (87348acc), and
    circulation_velocities' formulas (e9c00d7f; apply_secular_transport
    in named steps, fbc153a8).
- **DONE: state documented and tidied.** state/README.md and
  setup/README.md (4bea0c8b): every star member's meaning, unit,
  writers, readers and lifetime. 14 misleading members renamed, 4 dead
  ones removed, 1 made a local (b3a6675b).
- **DONE: the star as an argument** (802b3cdb). 81 core procedures take
  `type(star_info), intent(inout) :: s`; only run_yrec, yrec_reset, the
  C API and the tests name the global (core: 124 references, 2
  importers, from 3166 and 29). io/, rotation/, element_diffusion/ and
  setup/ still use the global.
- **Open, with proposed fixes (to decide later).** Found by the tests
  and the reviews of this series. Each entry says what the fix would
  change.
  1. **Triangle** (core/structure/triangle_fit.f90): reflect_until_inside
     never ends for a zero or negative tri_delta_teffl/tri_delta_logl, or
     a model file whose stored orientation disagrees with its vertices;
     the fit has no pivoting, so collinear fitting points, or two
     vertices at the same fitting-point log T, give Inf/NaN silently.
     The triangle is refitted only between iteration levels, so the
     Newton iteration itself always sees one linear boundary row.
     Proposed: check the two controls (> 0) at the read and the stored
     orientation at the model read; cap the walk (e.g. 1000 steps) with
     a clean error; stop with a message when the fitting points are
     near-collinear. Keep the arithmetic (pivoting would move every deck
     at round-off). Byte-identical for valid decks.
     **DONE** (d2ff6886): the walk gives up after 1000 reflections (the
     battery's decks need at most 11), the fit refuses a relative pivot
     below 1e-10 (the decks never go below 6e-3, on the 0.05 and 0.1
     Msun Allard decks), and the read refuses non-positive sizes. The
     model-read orientation check was left out: with the walk capped, a
     wrong orientation ends in the same clean error, and a stored
     triangle that the run discards would have been refused for
     nothing.
  2. **Rotation** (rotation/seculr): cz_moment_of_inertia is intent(out)
     but read before it is set (works because the caller set it); the
     dynamical-shear limit wmin in compute_shear_velocities comes from
     the plain difference quotient while its trigger uses the centred
     gradient, so it can overshoot omega(i-1) or go <= 0 (log NaN); the
     U_GSF correction picks its branch with wmin instead of omega(i-1).
     Proposed: intent(inout) (byte-identical); wmin from the same
     gradient as the trigger and the branch from omega(i-1), after
     measuring which decks reach that branch; one measured commit.
     **intent(inout) DONE** (2d62df76); the wmin fix remains.
  3. **Rezone** (core/model/rezoning.f90): the deuterium gate reads
     xa(h2, nz) with the pre-rezone nz; when it is false,
     deuterium_burning_rate_start stays in old-mesh order (read by
     mid_timestep_model). Proposed: read the new surface; regrid the
     rate unconditionally. Measured; expected drift at round-off or none.
  4. **Accretion mixing** (core/step/mass_change_kernel.f90): mixes as if
     the zone were uniform (X(jenv) M_cz); the step uses two definitions
     of M_cz; with jenv = 1 the zone compresses about m(1). Proposed: mix
     with the zone's true content, the sum of X dm, so species are
     conserved exactly; one definition of M_cz. Changes accreting runs
     only (no pinned deck accretes).
  5. **Jacobian** (core/structure/shell_equations.f90, eos/): the
     convective gradient's L dependence is not linearized; its P and T
     derivatives near the adiabat rest on the EOS's grad_ad derivatives
     (off by up to 72% of the entry, < 1% of qt); the OPAL tables' dln
     rho derivatives agree with finite differences to 2e-3. These change
     convergence, not converged results. Proposed: first run the 0.03
     Msun deck with numerical derivatives for these entries (item 8);
     fix them only if that cures the two-cycle.
  6. **State** (state/README.md section 5): about 19 members are written
     but never read; envelope_hydrogen_fraction is read by kap_init and
     the kind cards before the first card sets it; envelope_radius is
     written twice per model; check_star_contracts does not split
     statements at ';'. Proposed: remove the dead stores (the diagnostics
     stored but never written -- fxion_zone, scp, the two He3
     luminosities, pulse_dlnrho_dlnp -- and the always-zero MLT pulse
     arrays are removed, the author's decision of 2026-10-13); check whether
     kap_init uses the unset envelope X (likely overwritten before use)
     and set it first if it does. **Dead stores DONE** (c4f359c8,
     452101af), except dlnrho_dlnp and dlnrho_dlnt, which reach the next
     coefficient pass through an intent(out) dummy; the envelope X,
     envelope_radius and the ';' split remain.
  7. **Out-of-range read** (core/diagnostics/stitched_model.f90):
     atmosphere_radius_cm reads env_log10_radius(0) when the starting
     model's stitched model is built before any envelope is stored (it
     stops -fcheck builds at startup). Proposed: take the outer shell's
     radius as the photosphere when the envelope is empty (changes only
     the starting model's atmosphere radii in profile output); then add
     a -fcheck=all run of the config matrix to the battery.
  8. **Low-mass convergence** (one investigation): the 0.03 Msun Newton
     two-cycle and the 0.05 and 0.1 Msun Allard testsuite decks that
     stop converging. Proposed, on the 0.03 Msun deck: log the triangle
     fit's conditioning up to the failure (item 1), try numerical
     derivatives (item 5), and find where the corrections alternate
     (surface rows, convective interior, EOS table edges); a fix only
     once the cause is known.
  Proposed order: the byte-identical safety work (items 1 and 7's
  -fcheck run, 2's intent); the small measured fixes (3, 4, 2's limit,
  7's radius); 6 after the output decision; then the investigation (8,
  with 5).

### Scratch with a one-step lifetime; inlist group names (yrec-next, 2026-10-13)

Every commit byte-identical on the check's 21 decks, the core extras,
the Kurucz, Allard and Debye-Hueckel decks and test_reentry. The test
for scratch: poison a record with NaN at every step start; if every
deck stays byte-identical, nothing in it outlives a step.

- **DONE: the envelope integration's results are the star's.** The
  envelope and atmosphere points (s%env_points, s%atm_points,
  s%atm_point) go to the kernel as arguments (b7c5ef69); the
  integrands' point physics is the kernel's own variable (50bd4313).
- **DONE: the stitched model is the star's record, s%stitched**
  (1cc08918).
- **DONE: every surface EOS call passes its composition**; the
  Debye-Hueckel side channel (eos_set_debye_huckel_composition) is gone
  (fbdc6b97).
- **DONE: the central EOS call's Saha state is the star's** (9837d1a8).
- **DONE: inlist groups &control and &physics** (01931c03), the legacy
  YREC names, in place of &star_job and &controls; an inlist with the
  old names stops with the command that renames them
  (tools/upgrade_inlist.py DECK.inlist).
- **DONE: rotation** (0f2d50a9). moment_of_inertia_tolerance (the F77
  TOLERI, never set) is an optional argument of omega_from_j,
  enforce_rotation_profile and solid_body_omega; omega_mid, the
  solid-body iteration's first guess at the next step, is the star's;
  yrec_reset no longer restores rot_scr, mid_scr and circ_scr.
- **DONE: element diffusion** (a82ada96). The run log's "settling
  suspended" latch is the star's; reset_element_diffusion is gone.
- **Open, with proposed fixes.**
  1. **TOLERI is zero**, as in the F77, so solid_body_omega always takes
     its 20 Newton corrections except in the transport's loosened
     iterations. Proposed: a real tolerance (e.g. 1e-12), measured; it
     moves omega at round-off and saves time.
  2. **DONE (2d62df76): check_angular_momentum loosened
     rotation_settings%shape_tolerance** by writing the settings record
     and restoring it; the shape takes it as an argument now, as
     inertia_tolerance is passed.
  3. **Core writes rotation's scratch within a step**: the rezone stores
     the previous model's physics on the new mesh in rot_scr%old_* (and
     del_grad_diff_interface), and convection_geometry reads mid_scr's
     gradients. Proposed: the step driver owns that record and hands it
     to the rezone and the rotation driver as an argument, so
     rotation_scratch_lib is rotation's alone.

### rotation/ Stage B; element diffusion, setup and io take the star (yrec-next, 2026-10-13)

Every commit byte-identical on the check's 21 decks, the core extras,
the Kurucz, Allard and Debye-Hueckel decks, test_reentry and the unit
tests; the full battery runs on the series' end.

- **DONE: omega_profile is star-free** (25f04291): the disk-locking test
  takes the model's age (disk_gate_age_gyr) as an argument.
- **DONE: rotation's leaves are star-free** (ae031e5b): circulation's
  velocities and the Zahn factor are arguments; magnetic braking reads
  the solar units from the settings record and the step's tau_cz and
  P_phot from a wind_conditions record; viscos takes kappa and c_P and
  returns nu and K; the grids take L_sun from the settings record.
- **DONE: the drivers take s** (5241a1db); rotation names the global
  nowhere (budget (0, 0)).
- **DONE: core no longer writes rotation's scratch** (56aa2e6f). The
  start-of-step physics on the step's mesh is the star's step snapshot
  (s%old_gradT, ...); the rezone's store into
  rot_scr%del_grad_diff_interface was dead (poisoned after the rezone,
  every rotating deck byte-identical) and is gone; burn_li_be takes the
  gradients it locates the mixed region's base with. check_boundaries:
  rotation_scratch_lib and element_diffusion_scratch are their domains'
  alone.
- **DONE: eight rotation kernels are star-blind, with physics tests**
  (b98eb3bf). test_rotation_kernels: the spherical limits of I, fp and
  ft and their growth with omega (1 - fp as omega^2, 1 - ft as
  omega^4); omega_from_j gives each radiative zone and the convective
  envelope back their J; the AM diffusion conserves J (round-off), takes
  exactly the wind's surface term, keeps a solid body; the composition
  mixing conserves each species; the Kawaler exponents. The AM solve
  moved next to its system (am_transport_coefficients).
- **DONE: element diffusion is a star-free domain** (926f34ef): a
  settings record (diffusion_settings) installed per run card; the
  helium switch and the run log's latch come in through
  diffuse_elements. It joins STAR_FREE_DIRS.
- **DONE: setup takes s** (3c6b54fe): the settings builders read it
  intent(in); map_user_inputs writes it.
- **DONE: io takes s** (02c388b2): the 30 readers and writers that touch
  the star. solver_diagnostics() stays star-free (the run log's copy of
  the control, installed by the controls read), so its star-free callers
  are star-free in fact.
- Now only run_yrec, yrec_reset, controls_sync and the tests name the
  global star.
- **Open, with proposed fixes.** (Items 1-3 done 2026-10-28, apart
  from item 1's kernels and tests; see the next section.)
  1. **Rotation's module scratch.** secular_transport,
     mid_timestep_model, circulation, stability_coefficients, am_grid,
     equal_grid and rotational_mixing still share rot_scr, mid_scr and
     circ_scr, and so cannot be star-blind kernels yet. Proposed: per
     sub-step records passed explicitly, kernel by kernel, each with a
     check in test_rotation_kernels (the circulation and instability
     velocities against their limits: no shear velocity for a solid
     body, U_ES ~ omega^2 at fixed structure).
  2. **The controls read and the global.** read_controls(s) is handed
     the global star by run_yrec, while controls_sync's
     store_controls_to_star writes the global directly: the same object
     reached two ways inside one call. Proposed: controls_sync takes s
     (tools/gen_controls_state.py).
  3. **The global itself.** With every routine taking s, the star can
     become run_yrec's local (or a handle table), which the C API and a
     binary or multi-star driver would need (docs/TARGET_ARCHITECTURE.md
     section 4.2).

### No module-level state in the star layer; extra energy sources (yrec-next, 2026-10-28)

An overnight series. Every refactor commit is byte-identical on a
40-deck set (the config matrix, an accretion variant of
cm_massloss_const, the 12 testsuite pins, the four solar decks, the
m0030 chain, the two light-element decks, a ZAHB to TAHB run and the
3 Msun rotating TAMS deck) against a snapshot of bd56b718, and passes
the unit tests, test_reentry and (for the global star) test_pyyrec.

- **DONE: controls_sync takes s** (4dd2a637; the previous section's
  open item 2). The generator emits seed_controls_buffer(s) and
  store_controls_to_star(s); state/ is at (0, 0).
- **DONE: rotation's scratch records are arguments** (174a54a4; open
  item 1, the records part). evolve_angular_momentum allocates rot_scr,
  mid_scr and circ_scr for the step and passes them down;
  rotation_scratch_lib declares the types only, which the new
  SCRATCH_TYPES_ONLY rule in check_boundaries.py enforces. Filling the
  records with 0xff bytes, then zero bytes, at every allocation left
  every rotating deck byte-identical. Still open from item 1: the
  routines that take the records (circulation, stability_coefficients,
  am_grid, equal_grid, rotational_mixing) are not yet star-blind kernels
  with physics checks in test_rotation_kernels (no shear velocity for a
  solid body, U_ES ~ omega^2 at fixed structure); they still take json
  from star_info_lib.
- **DONE: element diffusion's scratch record is an argument**
  (7219c426), the same way: diffuse_elements allocates diff_scr.
- **DONE: no global star** (c35778aa; open item 3). run_yrec owns the
  run's star as an allocatable local; yrec_run_prologue copies a
  never-written static template into it (zero plus the declaration
  defaults, what the module variable was at process start), replacing
  the first-entry snapshot. The test programs declare their own star,
  with SAVE (8 MB overflows the stack otherwise; the first try
  segfaulted). Every directory's STAR_ACCESS_BUDGET is (0, 0).
- **DONE: extra energy sources** (6fc1a3f8; TARGET_ARCHITECTURE 4.5a,
  the Phase 2 exit criterion). core/structure/energy_sources.f90, a
  star-blind registry a program fills with add_energy_source before
  run_yrec; form_extra_energy_terms in the shell kernel; the
  i_lum_extra channel (n_lum_channels 9), counted by both
  renormalizations; history column lum_extra, profile column eps_extra;
  examples/other_energy (yrec_extra_heat, `make other_energy`);
  test_structure's 'extra' case and test_other_energy. Not yet: the
  built-in terms as list entries, the (lnP, lnT, lnR, L) partials of
  the design, and a warning for a source above the fitting point.
- **DONE: O23, the two dead EOS derivatives** (6252b537). dlnrho_dlnp and
  dlnrho_dlnt are locals of compute_henyey_coefficients; the star
  members and the arguments of henyey_iterate and read_starting_model
  are gone. Seeding the first shell's EOS with NaN left all 40 decks
  byte-identical, which confirmed that the seed was dead.
- **DONE: the aging flag replaces the negative dt** (88b32a90; a Phase 2
  item). solve_structure no longer negates s%dt for a model that is not
  aged: the first guess, the rates of change, the composition rerun,
  the light-element burn read evolve_model_flag, and the Henyey pass
  gets a zero step. The rotation step kept the step negated, so as to
  keep O24 (below) byte-identical until its own commit (2bdbfbad), which gives
  it a zero step too.
- **DONE: the binary's version line follows the commit** (31bbd69b).
  core/run/yrec_version.o depends on src/.git_hash_stamp, which make
  rewrites only when the hash changes; before, the main tree's binary
  reported e669801a long after that commit (the pinned logs carry it;
  every comparator skips the line).
- **Extra check.** 31 more decks, run with the bd56b718 binary and with
  88b32a90's, are byte-identical (outputs and terminal text): the 28
  ZAHB to TAHB decks (rescale-only cards, so the un-aged path) and the 3
  giant differential-rotation decks. Five of them fail the same way with
  both binaries, not investigated: the three giant decks stop at once
  ("LROT set to TRUE, but OMEGA(1) = 0": their shipped starting models
  carry no rotation), Solar_m1p0feh+0p0_PARSEC_TAHB stops while setting
  up its second card, and ZAHB_m0p7feh+0p0_TAHB goes NaN in the
  convection routine after about 190 s.
- **Found and FIXED: O24** (2bdbfbad; docs/bug_sweep.md). An un-aged
  rotating model without instability transport passed its negative step
  to the wind brake, which spun the surface convection zone up.
  Inherited from GETW. The rotation call now gets a zero step for a
  model not aged, and the no-transport branch applies the wind only over
  a positive step. No pinned deck reaches it (40 decks byte-identical);
  on a deck built to reach it (a rescale-only card at 0.1 Gyr, solid
  body, transport off) J no longer rises by 3.5% over the two un-aged
  models, and five models later J and Omega are 3% lower than before.
- **Method.** scratchpad check.py builds a worktree, runs the deck set
  in it in parallel, and byte-compares every output against the base
  snapshot and the main tree's pins (ignoring the version line);
  tests.sh builds the test programs (`make tests`; without them the
  domain tests skip) and runs every pytest file the tree has. The poison
  tests fill an allocated record through c_loc with 0xff or 0x00 bytes
  under YREC_POISON. A new deck under examples/ needs its controls hash
  (YREC_REGEN_CONTROLS=1 python3 -m pytest test_controls.py; check that
  only its line is added), and a new module or use statement needs
  `python3 tools/gen_makefile_deps.py`.

### rotation/ Stage C: tested kernels, a mechanism list, records by stage (yrec-next, 2026-10-28)

From the 2026-10-28 evaluation of rotation/ (readability 3.5, organization
2.5, maintainability 2.5, extensibility 2, out of 5): the drivers read as
named steps and the comments are good, but the three mechanisms were
hard-wired through about ten places, one 70-field scratch record coupled
seven modules, the velocities had no physics tests, and names misled.
Every commit is byte-identical on 56 decks: the 40-deck set, 15 variants
of `Test_solarGS98_difrotmixper18p4d1p0` stopped at 0.3 Gyr that switch
on the options no shipped deck uses (the Zahn factor, GSF inhibition
modes 0 and 2, solid-body modes 1 to 3 with a negative cz exponent,
walpcz = -0.5, no transport in the core, no background diffusion, the
basic composition, one iteration, no disk, the Kawaler law, no
transport), and a rescale-only rotating deck. The unit tests and the
MESA baselines pass on each.

- **DONE: hygiene** (effa1b10): the species ranges named
  (`last_mu_species`, `num_tracked_species()`); `solid_body_everywhere()`
  for the six copies of the solid-body test; the chi scales named
  settings (rotation no longer imports io's controls_buffer); the
  double-meaning `unstable_zone_found` split; `dlnomega_dlnr`, which held
  dω/d ln r, renamed; 17 never-written inout dummies made intent(in); the
  registry's descriptions of walpcz (an exponent, not a weight) and lvfc
  corrected.
- **DONE: tested kernels** (8c1a7f33): `instability_velocities` (the
  four velocity kernels, out of the 1043-line circulation.f90) and
  `transport_regions` (the region scan, which two routines had copied),
  both star-blind, with physics checks in `test_rotation_kernels`: ω²
  scalings, no shear term or GSF for a solid body, the μ factor between
  1/2 and 1, the shear's quadratic onset and its stabilization. Each
  check fails under a mutation of its kernel.
- **DONE: a mechanism list** (09d9c730): `state/transport_mechanism_layout`
  slots, `star%mechanism_velocity(n, json)`, per-mechanism
  `velocity_scale(:)` and `mixing_scale(:)`; the bookkeeping loops over
  the slots, and `mechanism_sum` reproduces both former sums bit for bit.
  Adding a mechanism is now five steps (rotation/README.md §4).
- **DONE: records by stage, one name per quantity** (f0e8b4b6): `rot_scr`
  split into `sub_model`, `stab` and `chi_grid` (with `sub_physics`, the
  former `mid_scr`, and `circ_scr`), each routine taking the records it
  uses; the `_mid` suffixes gone from the sub-step model; seven setup-only
  interface fields made locals; the secular_transport alias table
  emptied (`total_velocity`, `eq_*`, `diffusion_converged`).
- **DONE: tests of the transport grid and the map back** (d793d630):
  equally spaced chi points, grid masses that telescope to the region's,
  J conserved by the map back in a redistribution and a loss taken to
  its documented 1e-5, the convection zone below taking the region's
  first change.
- **DONE: the region's grid quantities in chi_grid** (78cb3b43): the six
  `eq_*` arrays and `single_interface` live in the grid record, so
  `am_transport_grid` (28 arguments to 20), the composition routines and
  `apply_secular_transport` stop passing them; `grid_spacing`, a copy of
  `eq_chi_spacing`, is gone.
- **Found: O25** (docs/bug_sweep.md): walpcz ≤ -2 (constant j) in a fully
  convective star gives ω ∝ r⁻² at the centre and NaN at the first model.
  Inherited from WCZ; needs the authors' choice of fix.
- **Seen: O13** in the grid test: with a radiative first zone the grid
  counts a convection zone below it and the map back leaves it alone.
- **Still open** (items 1, 4 and 5 done in Stage D, below):
  1. `apply_secular_transport` (558 lines): 16 internal procedures share
     about 20 host locals. Proposed: the per-region transport (grid,
     solve, map back, mixing) as a module procedure with explicit
     records, and the surface wind terms beside magnetic_braking.
  2. The sub-step re-burns, settles and burns deuterium on the
     interpolated model, calling back into core (rotation and core
     depend on each other). Proposed: a core service "advance the
     composition over an interval on a given model", which both the
     step and the sub-steps call.
  3. Argument lists of 15 to 26 remain where the drivers pass model
     arrays one by one to routines that tests call with synthetic
     arrays.
  4. `star%omega_mid` keeps the misleading `_mid` (it is the sub-step
     end's omega).
  5. Physics checks for `rotation_stability_setup`'s coefficients and
     `check_angular_momentum`.
- **Method.** scratchpad `decks/rv_*` (the variants) and `check.py --set
  rot|all`; `tighten.py` (every inout dummy flipped to intent(in),
  compiled, and the ones gfortran names as written put back, per
  procedure); each new check confirmed by a mutation of the code it
  covers.

### rotation/ Stage D: validated, restructured, the open rotation bugs fixed (yrec-next, 2026-10-28)

Asked for rotation/ to be "a really good module" (the re-evaluation
put it near 7 of 10), in six phases. Refactors were checked
byte-identical (the 56-deck set, or the 26 rotating decks with every
output file compared, and the test outputs); each results-changing fix
is its own commit, measured against the commit before on the rotating
decks (each run in its own directory, history.data compared column by
column; scratchpad measure.py, histcmp.py, cmp_runs.py) and, for the
settling fix, on the 50 self-contained decks of the set.

- **Phase 1, a safety net** (4e764df2): analytic and law checks in
  test_rotation_kernels (the stability coefficients against their
  formulas; the wind's spin-down against the saturated and unsaturated
  solutions, second order in the step; constant-D diffusion against its
  eigenmodes, h^2 and dt; the checks' contracts), each confirmed by a
  mutation; CI runs both rotation test programs; src/tools/coverage.py
  (gcov; rotation/ at 92.7%).
- **Phase 2, validation** (908b5501): docs/ROTATION_VALIDATION.md;
  sample_cases/rotation_vs_MESA (a solid body under J conservation
  against MESA 26.04.1: the main-sequence omega within 0.8% slow, 0.4%
  at omega/omega_crit = 0.25); src/tools/substep_convergence.py (first
  order in time_dw_mix). Found O26, O27, O28.
- **Phase 3, apply_secular_transport's steps as module procedures**
  (6dbf17b0): explicit arguments and intents for every step; the
  region's angular momentum diffusion (region_transport) and the outer
  point's wind terms (magnetic_braking's outer_point_wind_terms) as
  star-blind kernels with tests; the interface quantities in stab. O28
  fixed with it.
- **Phase 4, a mechanism interface** (8716c7e7):
  transport_mechanisms' compute_mechanism_velocity, slot by slot;
  the layout declares each mechanism's correction (the shear corrects
  the GSF) and controls; tests check the slots written and that the
  controls are the registry's and are copied.
- **Phase 5, the open rotation bugs**, one measured commit each:
  O26 J and E_rot of the model the step ends with (1add4edd; history
  columns only), O27 the surface region's implicit wind term
  (0ceb4632; one variant moves), O12 the Newton tolerance 1e-14
  (f1cbabae; all move at round-off, amplified to 2e-4 dex), O11 the
  dynamical-shear limit (42bd6b41; the five decks with dynamical
  shear), O17 the braking's 1/G, wind switch and tau_cz (920e6022,
  63b8325c, c5520bc0; the 14 Matt-law decks with surface transport for
  the first), O13 a region's end cell and write-back (ed8f7af6; no deck
  reaches it), O14 the late species checked, at most 16 iterations
  (a62b9741), O15 the safeguard removed (6574fd95), O25 constant j at
  the centre (8fb3c67b), max_diffusion_iters defaults to 10 (94755a59),
  O16 settling on copies and the MOD (c66901c9; the 23 settling decks,
  3e-5 dex in log L at 4.57 Gyr). The author questions on O13, O15 and
  O25 were answered with the bug list's proposed defaults; the bug list
  records them as revisable. Found O29 (open).
- **Phase 6, cleanup** (d51f2ea7, 13f9a1ba, 1ee8e6ed, 21d88db0,
  2f4a8f78): the commented-out code (830 lines) in rotation/retired;
  the thresholds as named parameters; the shear kernel pure (the
  dispatcher reports dynamical shear); star%omega_mid is
  star%omega_substep; rotation_messages (rotation_report to the terminal
  and the run log, rotation_note to the run log too with solver
  diagnostics).
- **Learned:** build with the default target (`make`, not `make yrec`):
  only the default adds -O3, and make does not rebuild objects when the
  flags change, so an unoptimized test build looked like a
  code-generation sensitivity until the flags were found. git's stash
  list is shared by every worktree. The machine-local MESA-format
  baseline of the rotating solar model pins history.data: reseed it
  after a fix that moves it.
- **Still open:**
  1. O29: the dynamical-shear limit takes the mean gravity's sign (only
     beyond break-up), and that run ends with exit status 0 after a
     BSSTEP error.
  2. The sub-step re-burns, settles and burns deuterium on the
     interpolated model through core (Stage C's item 2).
  3. Argument lists of 15 to 20 remain where drivers pass model arrays
     one by one (Stage C's item 3).
  4. The constant-j law (walpcz <= -2) spins a convection zone's base
     beyond break-up once a small radiative core forms (O25's variant
     stops at 1.67 Myr): physics for the authors.

### rotation/ Stage E: the Solberg-Hoiland instability (yrec-next, 2026-10-28)

Asked for the full list of hydrodynamic mechanisms against
Pinsonneault's thesis (section II.E: convection, Solberg-Hoiland and
the dynamical shear; Eddington-Sweet, GSF and the secular shear), the
2025 code lacked Solberg-Hoiland. Added as the fourth transport slot,
following the rotation README's recipe: the kernel
`compute_solberg_hoiland_velocity`, its dispatcher case, the controls
`solberg_hoiland_velocity_scale` (default 0, off) and
`solberg_hoiland_mixing_scale` (1), and the opt-in profile column v_sh.
With the shear's own N, its threshold is never below the
dynamical-shear limit, so it acts only where the dynamical shear does
(docs/ROTATION_VALIDATION.md section 5): off, the 50 self-contained
decks are byte-identical; on, 23 of the 26 rotating decks are, and the
three with the most dynamical shear move by at most 3.1e-4 dex in log
L. Seven kernel checks, each confirmed by a mutation (eight mutations,
all caught). The thesis's other exclusions stay out: ABCD and the triply
diffusive instability (no diffusion coefficient; YREC's GSF note equates
ABCD with the uninhibited GSF, its default), Schatzman's turbulent
diffusion, magnetic fields. The thesis's mu currents (eqs. 14-16) were
replaced by the Maeder & Zahn (1998) inhibition in 2002 (MHP 06/02).

### Module evaluation follow-up: eos, net, io, numerics, math, state, setup (yrec-next, 2026-10-29)

**Why.** After rotation/'s Stage D the owner asked for the other
modules to reach the standard rotation/ sets, on five criteria:
readability, organization, maintainability, extensibility (a newcomer
can add physics) and documentation that matches the code, and for the
bugs found on the way to be fixed. "Concise" meant the logic, never
fewer features or controls. Five read-only reviews scored the modules
first (908b983c, audit/module-evaluation-2026-10-29/SUMMARY.md; out of
5 relative to rotation/, for readability, organization,
maintainability, extensibility and docs): net 4/4/3/3/3, setup
4/3/3.5/3/3, math 4/3.5/3/3.5/2, numerics 3/4/2.5/3/3, eos 3/3/2/2/3,
io 2.5/2.5/2/2/3, state 2.5/2/2.5/2.5/2.5. The obstacles they named:
net's hand-kept index layouts; numerics' unchecked preconditions and
an ierr that reported success; io's controls pipeline (two default
sources that disagreed for 107 controls, the ctrl/job split kept in a
Python list); star_info in COMMON-block order with hand-kept member
tables; physical constants assigned at run time; tests that pinned
numbers and checked no physics; 5-10 false claims in every README.
Six bugs were confirmed by hand before any work began.

**How.** Six work packages, each in its own git worktree from
2c88d896 with the files it owned (NET: net/ and its builder; STATE:
star_info's layout and the state/setup READMEs; OUT: io's output, model
files, run cards and run files; NUM: numerics/, math/ and the
duplicated numerics in rotation/ and kap/; CTL: the controls pipeline;
EOS: eos/ and the physical constants). Each package kept a hand-off
note (its commits, bugs in this list's format, decisions with
revisable defaults, leftovers); bugs went into docs/bug_sweep.md
section 7 as M1-M34. Every commit was built, checked with
check_boundaries and the unit tests, and run on a deck set whose every
output file was byte-compared with a snapshot of 2c88d896 and whose
run logs were compared with the pins (scratchpad check.py): the quick
set (the configuration matrix, an accretion variant, a solar deck; 18
decks), the full set (40, with the testsuite and the pre-MS, ZAMS, ZAHB
and rotating TAMS decks) or all 56 (with the rotating variants). Each
results-changing fix is its own commit, measured against its parent
(history.data column by column, scratchpad histcmp.py). The packages
were merged onto one integration branch in order: the audit commit,
NET's 22 commits, STATE's 5, OUT's 17 with a regeneration of
state/README.md's tables (fb83611d), NUM's 23, CTL's 15, this
documentation, EOS's 27, and last STATE's member renames (a568e084). CTL's
s%job%-to-s%ctrl% rename was re-run as a script on the merged tree
rather than merged as a diff, and check_controls takes NUM's named
bsstep limits. On the merged tree: after NET and STATE all 56 decks
byte-identical; after OUT and NUM only the run logs' "stop conditions"
line differs (55 logs, M19); after CTL also the 23 rotating Matt-law
decks (9b58f6f6, below); after EOS every deck (M37-M39, below), and
the renames byte-identical. Unit tests at the end: 63 passed, 1
skipped (HDF5). The pins in standard/ were reseeded once, after the
series: pins reseeded at 2f84ae0e, the series' last code commit (the
previous pins are 2c88d896's outputs, which a build of it reproduces).

#### net/ (NET, 22 commits)

- **Layout** (2cba8030, 77dfc919, fa8cdec2): the chain energies are
  one array, eps_by_chain, which core loops over through
  luminosity_channel_of_chain and eps_channel_of_chain. One reaction
  layout in the new net_configuration: 11 reaction slots (the dead
  slots 9 and 13 gone, M6), the composition rows after them, the extra
  scales by name (xs_pep, xs_be7e, xs_be7p, xs_hep), per-slot
  derivative scales; star_info's scale arrays take their sizes from
  net_lib, and Monte Carlo and map_user_inputs address them by name.
  The builder is setup/net_settings_setup's pure (declared so since
  c96e6109; a function of an intent(in) star before)
  net_settings_from_star, installed by run_yrec with
  set_net_settings(..., ierr) after the table's self-check (every slot
  computed by exactly one row; the rate law, chain, branch, species
  known). net_settings_t's declared defaults are the controls'.
- **One job per routine** (78f58053, 61bbbd64, b5d5c721, f9d923c4,
  4ad0fce3, 3e9132bb, a1cf682b, a605b827): the thermal-neutrino chain
  is pure (its coefficients are parameter arrays), so GUIDELINES' net/
  purity exception went; nuclear_rates is filled in one place, the
  d(p,g) fit and the deuterium cap are written once; the energy sum
  only sums; every threshold and prefactor is a named parameter with
  its unit and source (true since a902a116, round 2: the deuterium
  threshold and the reference S0 values had copies); the screening is
  its own module; the network
  species are one table; a chain without a luminosity channel is
  skipped. One Sources text per file (ee5c6122), its contradictions
  resolved against the code, the F77 comments kept verbatim in
  net/retired/f77_rate_comments.txt.
- **Tests** (104f3b80, f2c6ebe6, bdad7653): test_net's ok/FAIL checks
  (the rates' and energy's derivatives against centred differences,
  both neutrino fits, deuterium burning with its cap, strong screening
  against Graboske's formula, the burn's non-convergence, the table's
  invariants), each confirmed by a mutation: 26 of 26 caught (two
  were not at first; the two later commits tightened those checks).
- **Bugs:** M1-M7. **Docs:** README to rotation's shape with recipes
  traced by adding a reaction (90d041e0), the glossary rows (7b0764af),
  coverage 90.6% of 1582 lines (4eed10cb).

#### state/ and setup/ (STATE, 5 commits)

- 49aa2d6b: star_info's 239 members in thirteen lifetime classes,
  longest-lived first, each class with a banner (what it holds, when
  it is valid) and each member on its own line under one comment in a
  fixed form (meaning, unit, writers, notes); history in a Sources
  section. 8ddbf6db: tools/gen_star_docs.py builds state/README.md's
  tables from those comments and the code; test_generated_files runs
  its --check, which also fails on a member without a comment, a wrong
  writer or a member nothing reads.
- f2558123: setup/README.md (each builder, where and when it is
  installed, the rules, recipes); 368078df: the solar calibration's
  four fixed derivatives are named parameters; 61e44415: setup/'s
  coverage, 95.2% of 252 lines.
- Prepared as a script: 62 member renames (one convention for earlier
  values, readable legacy names, the unit-suffix rule) and the stotal
  duplicate merged into log_total_mass, byte-identical over the 56
  decks on a scratch branch; applied after all merges (a568e084,
  below).

#### io/: output, model files, run cards, run files (OUT, 17 commits)

- **Bugs:** M18-M26, each its own commit (49e417d7, 34fc2988,
  8c80883a, b1e7ae2c with 30fae4e5, ba392477, e2335b37, 466b46ff,
  226fabae, 770f6ce2).
- **Structure, byte-identical:** the history and profile columns as
  named tables of (name, kind) with an index per column, checked by
  gen_default_columns.py and a sync test (770f6ce2); the output's state
  one record, output_state, that run_yrec owns (226fabae); run_cards
  with named card kinds and one echo (ccc50cca); model_files with each
  format's layout and one slot list (0103db35); declarations, messages
  through report_error, real(dp), use-only (b1f61cbc, 7bc6d4a6);
  stitched_model's comments and declarations (d8ae6c7b, 4d9b76d1).
- **Docs:** io/README.md with the recipes "add a history column", "add
  a profile column" and "add a field to the .mod model", traced on a
  copy; coverage 86.5% of 2065 lines (ed229bdb).

#### numerics/ and math/ (NUM, 23 commits)

- **Bugs:** M8-M17, each its own commit.
- **One error contract** (ecce23f0): numerics_status names the codes
  (singular, coincident points, bad input, step underflow, unsorted
  points, derivatives failed); the routines print nothing, and the
  callers that can meet a failure write the same messages to the same
  places; check_boundaries checks that numerics/ and math/ hold no I/O
  statement; numerics depends on nothing outside itself; the envelope's
  status is the kernel's envelope_termination (numerics_termination is
  kept for now).
- **One routine per job:** one natural cubic spline (4eaaaece); the
  tridiagonal solve in two halves that rotation's angular momentum
  solve reuses bit for bit (9c1c03b1); one 4-point Lagrange formula
  for numerics and rotation's three copies (f402d470); kap's hunt in
  numerics and LAOL's natural splines (a25ee9ca).
- **Readability:** quadratic_spline, lir and simeqc out of their F77
  shape (7ed8722a, 90475cd2, f1108365), bsstep, mmid and findex
  indented (61b0dfcb), the comments that named wrong callers
  (10fe4a65); every use of numerics_lib has an only-list, enforced
  (38ac9180).
- **Tests:** test_numerics keeps its pins and gains ok/FAIL property
  checks (exactness on polynomials, convergence orders, the error
  codes), each confirmed by numerics/test/mutate_numerics.py (55 of 55
  caught; 54bf5e79); test_math checks each math backend's contract
  (c2638e4b). math_lib's README, safe_divide_by_exp removed (023d5080).
  numerics/README.md with "which routine for what" and recipes; coverage
  91.1% of 827 lines (90619ab5); math 100%.

#### io/: the controls pipeline (CTL, 15 commits)

- **Bugs:** M27 (a0b4f434), M28 (af946c0c), M29 and M30 (db99ec9a).
- **One default source** (6085de4c): each control's default is its
  controls_buffer initializer (180 literals moved from
  set_legacy_defaults); test_controls checks that a deck setting
  nothing reads exactly controls_state().
- **One registry** (5f3835a3, 88b6c96f, 8398d475): the registry says
  each control's canonical member, its home (ctrl or job) and its use
  (run, read, none), with a row for every buffer member;
  gen_controls_state and gen_inlist check the declared types against
  it; the copies between the legacy and canonical names
  (adopt_canonical.inc, legacy_mirror.inc) are generated from it;
  controls_sync is generated into io/, private.
- **Homes** (97c00d5a, 937ada41, 96b1d4ab): the Bulirsch-Stoer
  controls are controls and envint_config fields
  (state/integration_parameters.f90 deleted), checked at read time
  against bsstep's limits (with nuse <= imax, 396bbb7c); the wind law's
  coefficients come from magnetic_braking's pure derive_wind_law,
  called by rotation_settings_from_star (ten pseudo-controls gone);
  star%job holds exactly the 36 controls that code outside the read
  assigns (34 moved to star%ctrl), and check_boundaries fails when the
  registry's home or use disagrees with the code.
- **Hygiene:** zsi and alex_mixture_x/z deleted, mu_gradient_scale
  retired (465bfc60); one validation module, controls_check, with
  check_controls, the warnings and report_control (db99ec9a); the
  modules' declarations, and io/read_controls.f90 renamed
  io/controls_reader.f90 after its module (7f111d37).
- **Docs:** the controls parts of io/, setup/ and state/README.md, with
  an "add a control" recipe traced end to end on a throwaway branch.

#### eos/ and the physical constants (EOS, 27 commits)

From 43e9b717 to 72de2124. From the first fix that moved every deck
(15af2ef3) on, deck identity was checked on 50 decks against the
previous commit's binary (every output file, and histcmp of
history.data), the quick set against a private baseline.

- **Bugs:** M35-M46, each its own commit. Results move with three:
  the Debye-Hueckel metals (15af2ef3, M37: every deck), the
  diffusion-of-Z atomic weights (3db309e9, M38: the 31 decks with
  diffusion of Z) and the Saha ramp's blend (f8114dca, M39: 49 of 50
  decks); the rest are byte-identical (MHD, which no deck uses; latent
  paths; tests).
- **Constants** (d6736788, d2c578ea, e78f6898): the physical constants
  of state/physical_constants (gas_constant, radiation_constant_over_3,
  cmkh, csig, csigl, cdelrl, cmixl2, cmixl3) are compile-time
  parameters, from CODATA 1986 with their sources, each checked
  bit-identical in hex against the run-time value; ca3l and cmkhn are
  gone. cgl (log10 G) depends on the control G_cgs, so it stays the one
  run-time value, set by run_setup from the pure log10_newton_g;
  convection reads it from its settings record
  (log10_gravitational_constant). test_eos and test_convection no
  longer boot a star. The two rotation-test commits before it made
  test_rotation_kernels and test_rotation set log10 sigma and a/3,
  which they had left at zero (the wind laws' radii had been 1/130 of
  the true ones, review_state_setup bug 3), so that the constants
  commit is byte-identical (the all set and the unit tests).
- **Tests** (c71c3ba0, fc15b3fd): test_eos builds the SCV envelope
  table as a run does (M46), and gains ok/FAIL checks of the physics
  (the ideal-gas limits with the first-order Fermi-Dirac correction,
  gas plus radiation, grad_ad -> 1/4, continuity at both ends of the
  Saha ramp and of the OPAL 2006 edge ramp, in_atmosphere, the
  envelope table after a composition change, the Debye-Hueckel
  mixture), each confirmed by a mutation; a FAIL exits 1.
- **Structure, byte-identical** (33d0dc45 to e87752bc): one table,
  opal_layout, of what tells the OPAL vintages apart (sizes, grids,
  slots, radiation constants, formats, windows), with one add_radiation
  for radsub01/06; MHD's table selector a family and an index;
  eos_init's 14 path arguments one eos_table_paths_t record, built in
  setup/ (equation_of_state_table_paths); the numerical derivatives'
  shifted evaluations trimmed; SCV's 4-point sums numerics' lagrange4;
  one opal_vintage setting (two selected is an error at eos_init);
  blend_eos as named steps (zone_mixture, the radiation, the gas, the
  selected table, through the generic adopt_table_result and
  blend_table_result); saha_state renamed saha_iterations, the counter
  it always was; eos_messages (eos_log to the run log, eos_report to
  both) for every message; beta14 one meaning, 4(1 - beta)/beta.
- **Docs** (6927b7ce): eos/README.md to the rotation standard (what
  lives where, the table of which source covers which log T, the two
  blends, conventions, every module-level state with its lifetime, the
  tests, recipes for a control, an EOS source and a blended region);
  the eos solvers in STAR_BLIND_FILES; coverage 85.1% of 3365 lines.
- **On the integration branch,** for decks the package did not run:
  the 0.03 Msun dbl-to-ZAMS deck, a marginal one, stops at a different
  age (final dlog L +6.9e-3 after 15af2ef3, +2.0e-3 after f8114dca),
  the 1 Msun ZAHB deck ends at dlog L -5.1e-4 and +5.9e-4 after the
  same two, and o24_rescale_rotating ends on its model budget at a
  different age; the refactors after f8114dca are byte-identical on
  all 56 decks.

#### The member renames (a568e084)

STATE's rename map applied after every package was merged, so that it
covers their code: 62 star members renamed (one convention for earlier
values, `<member>_start`, `<member>_prev`, `<member>_iteration_start`;
readable names for the legacy ones, jcz, iov1, senv, pphot, ...; the
unit suffix for the members not in CGS) and stotal merged into
log_total_mass. Only star accesses change, so output column names keep
their MESA spelling. Byte-identical: the 56 decks against the merged
tree's outputs before the renames, every per-deck file, and the unit
tests. The comments and documents that still named the old members
followed in a documentation commit (comments only: each changed source
compiles to the same object).

#### Results-changing commits

| Commit | Change | Moves |
|---|---|---|
| ee268284 | M1: eps_total includes the α captures | GYRE's eps in helium-burning cores only; no deck output of the set |
| 0ca5ff3e | M2: the capped deuterium rate's derivatives are zero | the Jacobian only; no deck |
| 85af8bef | M3: sneut's recombination derivative | no deck of the set; the 1 Msun ZAHB deck with the Itoh losses moves late in helium burning (final dlog L -2.9e-3, the runs parting near core-helium exhaustion) |
| 0c4fc520 | M4: triple α's z53 | no deck (none screens triple α strongly) |
| 4eaaaece | M11: one natural cubic spline | test_kap's OPAL92 and Kurucz 1990 pins only (<= 6.3e-16 relative; dlnkap/dlnT -0.32606 -> -0.35187 at OPAL92's top T node) |
| 49e417d7 | M18: the Monte Carlo model file keeps every run's model | only a calibrating Monte Carlo deck; none shipped |
| 34fc2988 | M19: the central-D stop label | the run logs' stop line (55 logs, 34 pinned) |
| 8c80883a | M20: a kind-3 card's stops | the run-plan line of 7 unpinned seed decks |
| af946c0c | M28: the settling latch per card | multi-card decks only; on a four-card deck built to reach it, at 1.5 Gyr dlog L +5.4e-5, dX_c -1.7e-4 |
| 9b58f6f6 | the controls' inexact single-precision defaults (sstandard, the CUS O fraction, the Matt law's solar constants, V13's m and c_2) are d0 | the 23 rotating Matt-law decks: at the logs' last line, \|dlog L\| <= 6.0e-5 (rv_sbm3), \|dlog Teff\| <= 2.7e-5, \|dX_c\| <= 5.6e-5 (rv_ref), model counts within 4; Test_solar_GS_rot final dlog L -1.2e-7 |
| d6736788, d2c578ea | test_rotation_kernels and test_rotation set log10 sigma and a/3, which they had left at zero | those tests' pins only (3 and 22 lines) |
| c71c3ba0 | M46: test_eos builds the SCV envelope table | test_eos's pins only (46 lines) |
| 15af2ef3 | M37: the Debye-Hueckel metals from eos's own mixture | every deck: the OPAL and SCV solar decks final \|dlog L\| <= 2.7e-5; the Yale-EOS solar deck +1.09e-3 (dX_c -7.9e-4); the 0.03 Msun dbl-to-ZAMS deck, a marginal one, +6.9e-3, stopping at another age; the ZAHB deck -5.1e-4 |
| 3db309e9 | M38: the diffusion-of-Z atomic weights | the 31 decks with diffusion of Z: final \|dlog L\| <= 3.0e-5, age-matched up to 3.5e-4 at 35-45 Myr |
| f8114dca | M39: the Saha ramp's blend | 49 of 50 decks: final \|dlog L\| <= 3.0e-5, age-matched up to 3.1e-4; the 0.03 Msun deck +2.0e-3, the ZAHB deck +5.9e-4 |

#### Still open

1. **M31:** use_structure_dt_limits is switched off for the rest of the
   run, and use_mass_accretion is the same kind of latch. Proposed: a
   star latch armed from the control at every card start, as M28.
2. **M32, M33, M34** and the read-time envelope Z: docs/bug_sweep.md
   section 7 and its "Decisions needed". The read-time Z (the run plan
   sets the envelope Z to the last card's positive rescale Z, which
   kap_init reads before the starting model sets it again; PARMIN did
   the same): proposed, move that assignment into prepare_run just
   before kap_init, so io stops reaching into eos (byte-identical).
   **Round 2:** the read-time Z and the EOS mixture moved to
   prepare_run (498cb8c7, byte-identical; which card's Z is meant is
   still the authors' question); M32 has four such controls (fo too),
   and M33 cannot be reached from a deck (derive_controls turns
   use_thoul_diffusion on with use_diffusion_z, as PARMIN did).
3. **net/:** the Caughlan & Fowler rates take no cross-section scale
   (M5; decided, revisable); sneut's bremsstrahlung derivative is good
   to 2e-6 of the total and was not investigated; the CF88 rates'
   resonant terms (above T9 ~ 1) have no check; the screening
   constants are literals per row (proposed: the charges Z1, Z2 in
   screening_t, so that the self-check can compute z_p); adding an
   isotope still means the literal 15 in about 40 declarations and the
   model format (Phase 5). **Round 2:** the screening constants follow
   from each row's two charges (68d01eec, measured: M63); the rest
   stands.
4. **state/ and setup/:** constants and copies on the star, each with
   a byte-identical plan in the STATE notes: solar_bolometric_magnitude
   (a parameter in physical_constants), hsra_log_t_offset
   (envelope_integration calls hsra_t_tau_offset),
   envelope_store_log_temperature and use_two_z_tables (computed where
   they are read). timestep_yr against dt/seconds_per_year is not
   byte-identical (the .mod header), and the cp, mu, dlnrho_dlnt pulse
   pairs are evaluations at different times, not duplicates: decisions
   with measured commits. The builders under one
   name, <domain>_settings_from_star, installed from one module
   (setup/settings_install.f90: per run, per Monte Carlo run, per card),
   with check_boundaries allowing the set_* calls only there.
   **Round 2:** done but for envelope_store_log_temperature and
   use_two_z_tables, which stay on the star: the two constants
   (08a876fd), the builders' one name and setup/settings_install.f90
   with its check (a51dc9b3).
5. **io/:** model_files' readers return 26 header values (a header
   record next); the card-kind constants are run_cards', while core
   still tests bare kinds (step_driver, calibration, starting_model;
   proposed, a constants-only module in state/); the history columns
   neut_diag_1/2 keep their opaque names (a rename is an output
   change); log_output_file could become a canonical control, so that
   the output takes its directory from star%ctrl; the namelist groups
   still read 99 controls into locals named after their inlist names
   and copy them through the legacy names (inlist_new_read.inc;
   binding the canonical names directly is a generator change touching
   every row). **Round 2:** the readers return one model_file_header
   record (641ce4c4); the rest stands.
6. **numerics/:** pure for splint, intpol, polint, ludcmp and
   solve_tridiagonal, now that they print nothing; ysplin, intpol,
   intrp2 and numerics_termination removed once their callers in kap/
   and element_diffusion/ move; descriptive names for the YREC-written
   routines (renames were outside this wave); lir onto interp's
   weights (results-changing: lir corrects the closest weight so the
   four sum to one). Rotation's solve_angular_momentum_transport still
   prints its own "singular matrix" line, against rotation's message
   rule; eos/scv's two splint calls no longer print splint's message
   (the EOS package's choice). **Round 2:** all done: pure (93977508;
   intpol is gone), the four removed (dadfb0e1), descriptive aliases
   beside the old names (38dbdeb9), lir on interp's weights (9b05b334;
   the weight it replaces is the first point's beyond the target, not
   the closest: be48e7e6), rotation's solve messages through
   rotation_report (2de142b0), scv's splint failures reported again
   (8ea4f7d2).
7. **eos/ and the constants:** cgl is still a run-time variable for
   core/, rotation/ and io/ (proposed: a log10_gravitational_constant
   field in rotation_settings and in the envelope kernel's cfg, from
   log10_newton_g, then cgl core-private). The Debye-Hueckel
   composition reaches the solver through yale_eos_tables' module
   variables (proposed: a dh_composition record built in eos_get and
   passed down). The Saha warm start is the previous call's solution,
   whoever made it (proposed: a caller-owned guess record; not
   byte-identical, so measured), and build_scv_envelope_table uses
   SCV's warm-start cell as loop scratch (proposed: locals; measured).
   A table-source descriptor for blend_eos (adopt_table_result and
   blend_table_result are its start); the OPAL tables are read once per
   process, so a rerun with another path keeps the first (proposed:
   record the path and re-read); rhoofp01 and rhoofp06 now differ only
   in data and could be one routine; in_atmosphere's one remaining
   effect, skipping the Saha second derivatives, costs 13 call sites
   outside eos to remove. Found, not fixed: R/mu jumps by 7e-4 at the
   Saha cutoff (the ramp's effective beta P/(rho T) against the
   composition's R(1/mu_ion + 1/mu_e); bug_sweep's "Decisions
   needed"), a point outside SCV's table falls back to Saha's result
   without a word, and a vintage switched on by set_eos_settings
   without eos_init would read unit 0 (test code only). **Round 2:**
   log10 G is a field of rotation's settings and of the envelope
   kernel's and the shell equations' records, and check_boundaries keeps
   the domains, rotation and the star-blind kernels off cgl (454014c0;
   core and io's pulsation writer still read it); the Debye-Hueckel
   composition is a record (8bcfb584); the Saha start is the caller's
   (be874f6c, measured: M76); the SCV builder no longer moves the
   lookups' cell (d8a034d8); one density search for 2001 and 2006
   (7bc47cd6); the unit-0 read is an error (771557ec, M75). Open: the
   table-source descriptor, the OPAL tables read once per process,
   in_atmosphere's call sites, the R/mu step and SCV's silent table
   edge (M88).

#### Learned

- check.py copies its extra decks (the accretion variant, rv_*, o24)
  into the worktree's examples/ while it runs, and test_controls reads
  every deck there: the unit tests run during a deck run fail
  test_controls falsely. Run one after the other.
- Load: five packages shared 12 cores, at load averages of 150-330; a
  full 56-deck run took 26-40 minutes. The packages grouped their
  verification (one rotating-set run covering several rotation
  refactors, the final full set on the last code commit, quick sets on
  copies of the tree while the worktree ran the long set), and each
  commit message says which run covers it.
- Generated files at merge: deps.mk, state/README.md (gen_star_docs),
  the controls' includes and dumps and the star contracts all follow
  the code, and a merge has to regenerate them. A test commit that
  did not regenerate deps.mk broke test_generated_files for every
  later commit of its package until a rebase fixed each one. A
  mechanical rename that other packages' files also need (s%job% to
  s%ctrl%) goes in as a script re-run on the merged tree, not as a
  diff; renames other packages would conflict with wait as a script.
- A mutation that a new check does not catch means the check is too
  weak, not that the mutation is harmless: NET's two misses (the
  strong-screening regime, a 5e-3 slip under a 5e-3 tolerance) each
  led to a sharper check.
- A commit message's "In the original F77?" is a claim to check
  against the import: M2's says the cap is not in the F77, and ENGEB
  has it.
- test_element_diffusion's stale metal coefficients (M33) make the
  unit tests fail now and then (3 in 100 runs at the base); a failure
  there alone, passing on a rerun, is that.

### Round 2: the re-evaluation's items (yrec-next, 2026-10-29/30)

**Why.** The five reviews of the module evaluation were repeated at
7aa138b2, the end of the module series above, on the same scale
(d4862175, audit/module-reevaluation-2026-10-29/SUMMARY.md; out of 5
relative to rotation/, for readability, organization, maintainability,
extensibility and docs): math 4.5/4.5/4/4.5/4.5 (average 4.4), net
4.5/4.5/4/3.5/3.5 (4.0), numerics 4/4/3.5/4/3.5 (3.8), setup
4/3.5/3.5/3.5/4 (3.7), state 3.5/3/3.5/3.5/3.5 (3.4), eos
3.5/3.5/3/3/3.5 (3.3), io 3.5/3/3/3.5/3.5 (3.3). Every module had
risen (io and state from 2.4, eos from 2.6). The owner asked for about
4.5 everywhere: every issue the reports list fixed, their "what would
raise each score by one point" items done, and their possible bugs
settled. Each report marked the first evaluation's findings fixed,
partly fixed or open, listed what the series had introduced (stale
comments after the merges and renames, a half-done rename in
star_info's Sources notes), and ranked new issues. Six claims were
checked by hand before any work began: card 1 of ten shipped decks
evolves an empty model (io/run_cards.f90:39 forced the legacy
variable after it had been copied into the control the run reads; a
one-line note in section 10 since August); meval's arcs after a point
on node n-1; eos/README's recipe step after OPAL, which blend_eos's
early returns skipped; the damaged Sources notes; O16(p,g)'s
single-precision reference S0; and log10 G without an initializer.

**How.** Six work packages again, each in its own git worktree from
7aa138b2 with the files it owned: NUM2 (numerics/, math/, and the
callers' handling of numerics' status in core/, rotation/, kap/, atm/
and element_diffusion/), STATE2 (state/, setup/, a new
setup/settings_install.f90, log10 G in the settings records, the star
contracts and the two star tools), NET2 (net/, its builder, the
channel sums and the burning gate in core/), CTL2 (the controls
pipeline in io/, defaults/, the generators and the registry, the
card-1 bug, the table paths end to end, map_user_inputs moved to io/),
OUT2 (io's writers, model files and run files, one io/README.md, an io
unit test, the stitched model's pulse points) and EOS2 (eos/, the
Debye-Hueckel leftovers in kap/, the read-time EOS mixture). Each kept
a hand-off note with the same sections as last time. Renames that
touch other packages' files were prepared as scripts and tested on a
scratch branch, not committed: STATE2's unit suffixes and `_prev`
rule, NET2's reaction-scale members, CTL2's 34 legacy-spelled control
code names, EOS2's Saha counter at eos_get's callers. Verification as
in round 1: every commit built and checked (check_boundaries, the star
contracts and docs, the generated files), the unit tests, and the
scratchpad check.py, which byte-compares every per-deck output file of
the 56-deck set with a snapshot of 7aa138b2 and the run logs with the
pins (which were 7aa138b2's outputs). Several decks of one directory
share history.data and inlist_used, a race under parallel runs, so
only per-deck files are compared. A results-changing fix is its own
commit, measured against its parent (history.data, scratchpad
histcmp.py); decks the sets do not hold were run on their own before
and after (the card-1 decks, a ZAHB-to-TAHB deck with GYRE files, a
solar deck with the default screening threshold, decks reading a .mod
or a YREC7 starting model).

The packages were merged onto one branch (mods2) in order: the audit
commit d4862175; NUM2's 22 commits; STATE2's 13; NET2's 23; CTL2's 17;
OUT2's 9 and 685f0d70 (CTL2's ready-to-paste controls text in OUT2's
merged io/README.md, and the other READMEs' pointers to its section
numbers); EOS2's 34; then the four rename scripts, each its own commit:
0197827c (STATE2's), 6882c0eb (NET2's), 692bc883 (CTL2's), 810a552c
(EOS2's). Checked on the merged tree after each segment, per-deck
files: after NUM2 and STATE2 all 56 decks identical to the snapshot;
after NET2 every per-deck file equal to NET2's own final outputs (its
first two fixes move every deck, below); after CTL2 only the run logs
differ, each gaining the controls warnings' lines (55 logs; one deck
identical); after OUT2 identical; EOS2 merged in segments, so that
each results-changing commit was measured on the merged code: 4a06cd08
moves two decks by at most 1e-8, the Debye-Hueckel metals (d3a4e8c1)
and its mass unit (bc1bb502) every deck, by the orders of magnitude
EOS2 measured (the Yale-EOS solar deck +1.73e-4 in final log L, the
0.03 Msun deck 3.9e-4 age-matched with a later ZAMS stop), the Saha
warm-start record (be874f6c) 12 decks (the settling solar decks by up
to 2.4e-5), every other EOS2 commit identical; each rename commit
byte-identical on all 56 decks. At the merge, test_structure's Yale pin
(NET2's two fixes and EOS2's metals both moved it) was regenerated on
the merged code; the state2 script missed one file (test_io declares
its star at module level, which its scope test did not take for the
star) and the integrator's copy covers it. Unit tests at 810a552c: 74
passed, 1 skipped (HDF5); 63 before the round. The card-1 fix was
checked on the merged binary on the nine decks it reaches, to a cap of
300 models. Bugs went into docs/bug_sweep.md section 8 as M47-M88. The
pins in standard/ were reseeded once, after the round: pins reseeded
at 810a552c, the round's last code commit (for NET2's two fixes, CTL2's run-log warnings and EOS2's
four moving fixes, below).

#### numerics/ and math/ (NUM2, 22 commits)

- **Bugs** (M53-M57): meval's unset tail arcs (6721d0b3: meval in three
  steps around one "arcs of interval i" record, checked exhaustively,
  every sorted subset of 16 candidate points against the points one at
  a time, bit for bit); osplin's NaN from a zero-width end interval
  (d367ecd0); simeqc's rank-deficient success (def22898, ludcmp's
  relative test; the smallest scaled pivot on five decks is 6.2e-7);
  the conductive opacity's undefined hunt guesses (4a814b28); bsstep's
  failure state (acb65d9d: y restored; the X + h == X test at X = 0
  documented, not changed).
- **Status codes reach a message:** every osplin call checks its ierr
  (rotation's nine through a new rotation_numerics_report, c7d8beff;
  the regridding through rezone, ed87d288); every caller that reports
  a numerics failure writes "<routine>: <numerics_status_text>", and
  the rotating shape's polint failure, silent since rotation's kernels
  stopped printing, is reported (d4bcaf27); compute_quadrupole's
  private Thomas solve is solve_tridiagonal, bit for bit, and
  rotation's three solve messages go through rotation_report
  (2de142b0).
- **One routine per job:** ysplin, intpol, intrp2 and
  numerics_termination removed, intrp2's 15 callers calling
  interp(..., x_eval=) bit for bit (dadfb0e1); findex is hunt with
  YCK's codes outside the grid (c6e3baf8); lir takes interp's weights
  and lagrange4's sum (9b05b334; unit-test pins move at round-off, no
  deck reaches lir).
- **Contracts:** cspline and natural_spline check n >= 2 and monotone
  abscissae through an optional ierr (75a98edc); splint, polint, ludcmp
  and solve_tridiagonal pure (93977508); the singular contract stated
  once (dense solvers relative, 1e-14 of the row; Thomas an exact
  zero); check_boundaries reads math_lib's shadowed functions from its
  public lines (hypot, expm1, log1p and log2 were missing from the
  hand-kept list) and rejects `stop` in numerics/ and math/
  (f062c956).
- **Readability:** hunt converted (4a814b28); polint, ludcmp and lubksb
  indented, their steps named (03976ee3); descriptive aliases beside the
  old names (lagrange3_weights, lagrange4_weights,
  equidistant_quadratic, interpolate_tabulated_functions,
  bracket_index; 38dbdeb9).
- **Tests:** simeqc on a pivoting matrix, a round-off pivot, a row
  scaled by 1e-30 and a NaN column; bsstep's large steps from their own
  derivatives, through every acceptance branch (fd8b0e1e); math's
  contract stated once, a real single-precision check and
  math/test/mutate_math.py (42048e86); the mutation scripts check their
  texts, run by the pytests (c184bfc0): 71 numerics and 17 math
  mutations, all caught.
- **Docs:** the comments that disagreed with the code, GUIDELINES'
  numerics rule and HOW_YREC_WORKS' failure rule (be48e7e6); both
  READMEs' recipes re-traced by adding a module on a throwaway branch,
  with the coverage of now (2cd53247; numerics 90.6% of 790 lines).

#### state/ and setup/ (STATE2, 13 commits)

- **log10 G** (454014c0, M81): cgl starts at its default; rotation's
  settings, the envelope kernel's configuration and the shell
  equations' structure_settings carry log10_gravitational_constant,
  which their kernels read; check_boundaries keeps the physics domains,
  rotation and the star-blind kernels off cgl; the rotation tests'
  hand-set -7.17571 went. 56225890: physical_constants is private and
  real(dp), and the version strings live in core/run/yrec_version.f90
  (a new commit hash rebuilds one file).
- **One install module** (a51dc9b3): setup/settings_install.f90 holds
  every star-layer installer call, one procedure per lifetime
  (install_run_settings, install_monte_carlo_settings,
  install_card_settings, install_surface_mixture_settings);
  check_boundaries' INSTALL_FILES rule allows the set_* calls only
  there; one builder name, `<domain>_settings_from_star`; the controls'
  slot names in state/control_slots.f90, so rotation's builder no longer
  imports io's buffer.
- **The star's contents:** be7_mass_fraction_zone and
  neutrino_flux_zone are locals of the neutrino table (M82),
  zahn_coupling rotation's scratch, envelope_specific_entropy a local of
  change_mass (33760aa0); solar_bolometric_magnitude a parameter and
  hsra_log_t_offset atm's pure function at the envelope kernel's
  configuration (08a876fd).
- **The checkers:** check_star_contracts keyed by module and host,
  resolving calls by scope and use, a read statement's target a write
  (f5a88324, M83); contracts on every core/run and core/step driver
  (10d42185); gen_star_docs' class rules (earlier values in the
  snapshot class and nothing else there; members one module uses in
  "Private to one module", which replaces "Scratch"; 737c00a1);
  check_settings_fields: every field of a domain's settings record is
  set by its builder (4f99f42c: the seven builders set all 143 fields).
- **Docs:** state's prose after the merges (eps_total since M1, the
  damaged Sources entries, README section 4 a table that gen_star_docs
  checks against the live members; b96337ee); the READMEs' three
  recipes traced by doing them on a throwaway branch (a star field, a
  settings field, a physics domain), and what they missed added
  (e17c33a7, de59b0da); setup/'s coverage 95.3% of 316 lines.
- **Renames** (0197827c, the script after the merges): the unit-suffix
  rule (h_shell_*_radius_fraction, fit_point_logP/logT/logR,
  neutrino_flux_total_e10_per_cm2_s), `_prev` only for the previous
  model (cz_base_radius_last_call, envelope_cz_base_zone_last_call),
  i_eps_3alpha, and the names at the star's edge (the envelope kernel's
  configuration, rotation's sub-step physics, two dummies) taking the
  members' names.

#### net/ (NET2, 23 commits)

- **Results-changing, each measured against its parent:** O16(p,g)'s
  reference S0 in double precision (efefc632, M58); one Avogadro number,
  N_A = 1/m_u, for every rate (dbb5f570, M59); the Itoh et al. losses'
  mixture, every nucleus plus the rest of Z as Ne20 (e1d89c2f, M60); the
  screening constants computed from each reaction's two charges
  (68d01eec, M63). Moves: the table below.
- **Organization:** network.f90 is the reaction and species tables and
  their readers; the self-check and set_net_settings are network_checks
  (which also checks that a slot's branch factors add up to 1;
  2ca3ce94); the per-reaction slots are state/reaction_layout, a
  constants-only module, so star_info uses no physics facade (ec2eaaf6,
  2094bf02); the energy path's rate floor goes by rate law, not slot
  range (fb771cfa); core's shell equations gate on the burn result's
  `burning` instead of their own copy of the cutoff (a8b59825); the hep
  factor is computed where the fluxes read it (23a134b7); the facade
  exports only what code outside net/ uses (2c485f8a);
  net_settings_from_star is pure (c96e6109).
- **One definition each** (a902a116, d5199a28): erg_per_mev, the
  deuteron mass, the deuterium threshold (three copies), and the
  reference S0 values that map_user_inputs repeated: a table row's
  reference S0 and reaction_rates' pep, Be7, hep and N15 values, beside
  the named multipliers they set.
- **Extensibility:** the luminosity channels summed by one loop that
  both renormalizations share (b2f0ece2); burn_network_region reads the
  network species at their slots in the full composition (3d4ce49b,
  M62), its arguments sized n_species_extended (72b28bfc); the species
  table names each species' screening isotope (10580145); a reaction's
  screening is its two charges.
- **Errors:** burn_implicitly names its two failures, and burn_step
  reports each through report_error (360944cf, M61).
- **Tests:** test_net's 25 checks (the builder and the install
  component by component, the Be7 branch fractions, every slot's rate
  positive at a hot point, the chains' sum against eps, the singular
  burn); compositions by slot name; a FAIL exits 1 (a4b42588,
  5d055e98); 14 more mutations, all caught (40 in all).
- **Docs:** net/README.md against the code, every user outside net/,
  the four recipes rewritten from a trace on a throwaway branch (a
  reaction with a new rate law, an S-factor reaction with its controls,
  a chain, Ne20 as a species: the trace stopped at the star layer's
  literal 15s and the 15-slot model formats); coverage 90.5% of 1611
  lines (90d39ba2).
- **Renames** (6882c0eb, the script): the star's qs0e_scale,
  qqs0ee_scale, be7p_qs0e_scale, be7p_qqs0ee_scale, o16_gamma_scale and
  c12_alpha_scale take the settings' names (s_prime_scale,
  s_double_prime_scale, be7p_s_prime_scale, be7p_s_double_prime_scale,
  n15_pg_scale, n15_pa_scale: the last two scale N15 + p's branches,
  not the reactions their old names said), and map_user_inputs'
  Solar Fusion I references follow.

#### io/: the controls pipeline (CTL2, 17 commits)

- **Bugs:** card 1 starts from the starting model whatever
  first_call_flag(1) says, as PARMIN forced (fb46cfe8, M47); num_runs
  outside 0 to max_runs stops the read (6730c3ae, M48); the copy of
  the deck goes to the run's output directory (34997232, M49); the
  Monte Carlo paths expanded under their canonical names (be059f10,
  M50); every read resets star%job too (957ec3da, M51); four table
  paths' registry descriptions (745af296, M52).
- **The read is one pipeline:** it opens no file of the run; run_yrec's
  open_run_files makes the output directory, copies the deck and opens
  the run's files after it, and the read's messages reach the run log
  in order when that opens, through io_units' pending lines (c5532dfe).
  A controls error now stops the run before the run log exists, so its
  message is on the terminal only.
- **Homes:** every table path is a control with a registry row (home
  ctrl) and reaches its domain in a paths record built in setup/ (kap,
  eos, atm); the hand-declared star_job paths, the table-path side
  channel and its load/store are gone, and the surface-mixture weights
  and kap's LAOL header mixture are locals (745af296). star%ctrl holds
  the 255 controls the run reads; the 76 input-only and unused ones
  stay in the read's buffer, registry home `buffer` (88b19e32). The 20
  dead unit-number controls are deleted, and the run's files get
  newunit units in io_units (6682d60c). The solar units are three
  controls copied to the star the same way (99ab5565); the OPAL EOS
  vintage is one derived control, which the eos builder copies, with a
  warning when two are selected (12faddeb).
- **The registry:** one status vocabulary (live, legacy, unused,
  derived) and an aliases column: six misleading inlist names replaced,
  the old ones still read (new_zone_omega_weight, alex95_table_z_list,
  kurucz_atm_table_file, kurucz_castelli_atm_table_file, ll_data_file,
  opal_eos_table_file); surface_composition_file accepted and ignored
  (1cd89324). Each limit and reference value written once, and the
  diffuse_lithium check that could not fire removed (e22cf809).
- **Messages:** one API in io_units (report_error, report_warning,
  report_note, log_line, beside check_io); a controls warning now
  reaches the run log too, as rotation's reports do (d694adf9: the run
  logs change, below).
- **Code style:** map_user_inputs moved to io/user_inputs.f90 with git
  mv (9cef61be), its parameters named and its derivative scales one
  function (ba210753); derived_controls and run_files import what they
  take from the buffer by only-list (7833983b).
- **Docs:** io/README.md's controls parts (pasted into OUT2's merged
  README by 685f0d70), with recipes for a control, a control core reads
  directly and a table path, traced with three dummies on a throwaway
  branch.
- **Renames** (692bc883, the script): the 34 controls whose code name
  was a legacy spelling (npenv, rsclzc, acfpft, atime, alpham, htoler,
  fcorr, fo, ies, imu, lvfc, tgcut, walpcz, ...) take readable code
  names, most of them the inlist name the decks already use; the decks
  read exactly as before.

#### io/: output, model files, io's own tests (OUT2, 9 commits)

- **Bugs:** GYRE's eps derivative columns are eps_nuc's, as its eps
  column is (9168a576, M77); the interior GYRE file is the stitched
  model's interior points, so N^2 has one formula (86b3fa63, M78); an
  empty starting model is a clean error (55cf35e5, M79); the profile's
  eps_ppI .. eps_he3 where the α captures burn, decided and documented
  (aca05304, M80).
- **Structure, byte-identical** (the quick set and eight paired decks
  that reach what the sets do not: GYRE, FGONG and isochrone files, a
  .mod and a YREC7 starting model, the two-run calibrated Monte Carlo
  model file): model_files names its layouts and record formats, and
  the three readers return one model_file_header record instead of up
  to 26 values (641ce4c4); the writers read s%stitched from the star
  (star rule 3) and name their formats (97d5370b); build_pulse_points'
  regions and species count named (64228e14).
- **Tests:** test_io, an ok/FAIL program for io's own procedures
  (expand_value, directory_of, pulse_format_code, parse_columns,
  nml_file_has_group, a bit-exact .mod round trip, synthetic YREC7 and
  MODEL2 files good and damaged, the column tables): 65 checks, each
  confirmed by a mutation (d4b2cb38).
- **Docs:** io/README.md one README in rotation's structure (74fda0d5,
  then CTL2's controls text in 685f0d70); io/'s coverage 94.0% of 2275
  lines (86.5% of 2065 in round 1).

#### eos/ (EOS2, 34 commits)

- **Bugs** (M64-M76): SCV's general path searched an unbuilt column
  (4a06cd08); under SCV, R/mu and 1/mu_e were Saha's (26751395); the
  Debye-Hueckel metals did not add up to the zone's Z (d3a4e8c1); the
  diffusion of Z divided by Z_env = 0 (c6efd005); the OPAL 1995 and 2001
  trial calls rewrote every slot (12c42467); read_scv_tables' stale
  envelope table, leaked unit and unchecked reads (835a923b);
  eos_set_mixture and eos_update_envelope_tables as a pair (debdb797);
  the Debye-Hueckel coefficient's transposed Avogadro digits
  (bc1bb502); eos_settings_t's defaults (65ca207c); errors that reached
  the run log only, and lir's and splint's dropped statuses (8ea4f7d2);
  blend_eos's early return after OPAL (2b671c3a); the Fermi-Dirac and
  OPAL readers without iostat, and a read of unit 0 (771557ec); the
  Saha solver's start from the last solve anywhere (be874f6c). Four
  more were found and not fixed (M85-M88).
- **io's bug 10:** the run plan's envelope Z and the EOS mixture moved
  from the controls read to prepare_run, unchanged (498cb8c7): io no
  longer reaches into eos.
- **Structure:** blend_eos's four steps are subroutines (zone_mixture,
  radiation_share, gas_eos, opal_step; 2b671c3a); one density search
  and one dispatch for OPAL 2001 and 2006, their instances an array,
  the edge ramps by layout (7bc47cd6), the 1995 code on its layout row
  (5f6a9cf8); one SCV cell search and one ionization step (4f70b867),
  and the builder no longer uses the lookups' warm-start cell as
  scratch (d8a034d8); the Saha solver's nests split (charge_balance,
  eliminate_helium) with named matrix terms (c7adcf24); the
  Debye-Hueckel composition a record eos_get builds and passes down
  (8bcfb584), and kap's LAOL reader keeps its header mixture to itself
  (dd7349be); one table of the composition terms' 1/A (e732612b); the
  constants from physical_constants (58c2525d); the envelope mixture's
  members named (inverse_atomic_weight, number_fractions; 3096b145);
  eos_lib exports eos_init, eos_get and eos_set_mixture only
  (64500ba9).
- **Tests:** test_eos's ok/FAIL checks from 27 to 49: continuity of
  every blended field, under SCV too, for any source (c2863a17); the
  derivative check judged where agreement is expected (da05c29f); the
  SCV, Fermi-Dirac and OPAL readers on missing and truncated tables;
  the settings' defaults and refusals; R/mu under SCV; the
  Debye-Hueckel record's metals adding up to Z (1fb891ca); the OPAL
  trial calls; lir's status. Each is confirmed by a mutation.
- **Docs:** eos/README.md rewritten against the code, with recipes for
  a source, an OPAL vintage, a blended region and a control; the source
  recipe traced by adding a throwaway source (9ce8e797); coverage 85.5%
  of 3262 lines.
- **Rename** (810a552c, the script): eos_get's Saha counter is
  saha_iterations at all 16 callers outside eos/ (it was saha_state,
  ksaha or saha_flag).

#### Results-changing commits

| Commit | Change | Moves |
|---|---|---|
| fb46cfe8 | M47: card 1 starts from the starting model | the ten card-1 decks, none in the sets: the nine that find their starting model stopped at model 1 and now run (to a 300-model cap: the 3 Msun decks to 0.34 Gyr, the TAHB decks through card 2; uncapped, RGp2diffrotI0 reaches central Y < 0.5 at model 2514 and the PARSEC TAHB deck central Y < 0.001) |
| 34997232 | M49: inlist_used in the output directory | the copy's location for the decks whose run log is not in {YREC_OUTPUT} (six class_exercises decks and about 190 sample_cases decks) |
| d694adf9 | controls warnings reach the run log (CTL2's D8) | the run logs of the decks that trip one (55 of 56: the envelope-overshoot warning in 54, the core-overshoot one in 25, the solar-wind one in 7); no other file |
| efefc632 | M58: O16(p,g)'s reference S0 is 9.4d0 | 35 of 56 decks; Test_solarGS98_base final dlog L +4.1e-6, dX_c +7.0e-5 (age-matched 1.0e-4 near the ZAMS, that deck's noise floor) |
| dbb5f570 | M59: one Avogadro number | all 56; the solar deck final dlog L -2.0e-5, dX_c +5.2e-5; the 3 Msun TAMS deck (testsuite, outside the set) -3.4e-4 in log L at fixed X_c late on the main sequence, a linear response to the pre-MS deuterium burning |
| e1d89c2f | M60: the Itoh et al. losses' mixture | no deck of the set; the 1 Msun ZAHB deck with the Itoh losses within 1e-6 at fixed Y_c down to Y_c 0.2, then parting at the breathing pulses (final dlog L +2.9e-3) |
| 68d01eec | M63: screening constants from the charges | no shipped deck (all screen weakly); with the default threshold 0.03 the solar deck final dlog L +4.2e-6 |
| 9b05b334 | lir on interp's weights | unit-test pins only (test_eos one line, <= 9e-16; two test_numerics pins, 1 ulp); no deck calls lir |
| 9168a576 | M77: GYRE's eps derivatives | GYRE and GSM files only, above log T 7.5 (no shipped deck writes them) |
| 86b3fa63 | M78: the interior GYRE file | its N2 column only (round-off, and its last point; no shipped deck writes it) |
| 55cf35e5 | M79: an empty starting model | an error path: exit 1 with a message instead of a runtime error |
| 4a06cd08 | M64: SCV's general path | 2 of 56 decks, by at most 1e-8 |
| d3a4e8c1 | M66: the Debye-Hueckel metals beyond O are Z less C, N and O | every deck: the configuration matrix 4e-9; the solar, rotation and pre-MS decks up to 2e-5 (age-matched up to 3e-4); the Yale-EOS solar deck +1.8e-4; the ZAHB deck +7e-5; o24_rescale_rotating +9e-4 (its last card ends 7% later); the 0.03 Msun deck 3.9e-4 age-matched (its ZAMS stop at 1.242e7 instead of 1.216e7 yr) |
| bc1bb502 | M71: the Debye-Hueckel mass unit | every deck: the configuration matrix and pre-MS decks 5e-13 to 6e-11; the solar and rotation decks 1e-9 to 4e-5 (age-matched up to 3e-4); the 0.03 Msun deck 4.5e-4 and o24_rescale_rotating 3.9e-4 age-matched |
| be874f6c | M76: the Saha solver's start | 14 of 56 decks in EOS2's measurement (12 on the merged tree): the settling solar decks by up to 1.5e-5 (2.4e-5 on the merged tree), rv_sbm2 8e-6, the rest below 3e-9 |

#### Still open

1. **numerics/:** rotation's log_lagrange4 still takes an assumed-size
   a(*) (a(:) changes the calling convention: a rotating-deck run);
   cspline's end slopes are not optional arguments; kap and atm pass
   numerics' codes up as their own ierr (revisable: kap-own codes);
   splint's, osplin's and lir's searches keep their own node
   conventions beside locate and hunt (folding them in is not
   byte-identical in principle); bsstep's underflow reference stays the
   hydrogen fraction, as in the F77 (revisable: the independent
   variable would end a failing X = 0 envelope step sooner). CI builds
   and runs test_math but neither test_numerics nor test_io (proposed:
   add both programs to the build line and their .py files to the
   pytest line).
2. **state/ and setup/:** check_settings_fields checks that every field
   is assigned, not that the right control feeds it (a registry column,
   control to settings field); star_job's members have no comments and
   gen_star_docs does not cover them; dlum_dx ... drad_dalpha stay on
   the star (M34's fix, one record) and star_mass and log_total_mass
   stay a pair (they advance by different arithmetic: not
   byte-identical). The registry's description of tolerance_fraction
   still names senv (log_fit_point_mass_fraction since a568e084).
3. **net/:** adding an isotope still means the literal 15 in about 40
   declarations outside net/ and the 15-slot model formats (Phase 5);
   the Caughlan & Fowler rates take no cross-section scale (M5,
   decided); test_net has no finite-difference check of the new
   eps_nuc derivatives that M77 writes (proposed: against
   total_energy_gen_rate - eps_neutrino_loss).
4. **io/, controls:** controls_buffer and legacy_controls are still
   public modules, imported whole by the read's own files (proposed:
   gen_inlist emits the only-list include); the namelist groups still
   read into locals named after the inlist names and copy through the
   legacy names (binding the canonical names directly changes every
   row's generated code; the aliases would carry over); the diffusion
   implications and the molecular-table precedence stay in
   derive_controls (core reads two of the switches they set); M31 and
   the use_mass_accretion latch; M84 (a sample deck's missing starting
   model).
5. **io/, output:** FGONG's var(9) is eps_total, where MESA writes
   eps_nuc + eps_grav (the owner's call); an α-capture profile column
   (M80) needs an eps_channels row; the history columns neut_diag_1/2
   keep their opaque names; open_run_files' `print *` lines; core still
   tests bare card kinds (step_driver, calibration, starting_model).
6. **eos/:** the R/mu step at the Saha cutoff (bug_sweep's "Decisions
   needed"); M85 and M86 (the EOS mixture after a later card's Z
   rescale, and before the first card's update_surface_mixture; both
   core/model, results-changing), M87 (a sharper Debye-Hueckel check),
   M88 (SCV's table edge a step); on a second eos_init an OPAL vintage
   keeps the table read first (proposed: open_opal_eos_table clears the
   vintage's loaded and first-call flags); SCV's search start (idtt,
   idp) is still module state (proposed: in eos_result, as the Saha
   start); a table-source descriptor for blend_eos (worth it at a third
   table source); OPAL's report_* keep ten per-vintage message branches;
   set_eos_settings stays public in eos_configuration (test_eos installs
   settings without eos_init); in_atmosphere's call sites outside eos
   were not audited.
7. **Comments and docs the renames left:** comments still name old
   code names or old homes (setup/net_settings_setup.f90's "run_yrec
   installs it", io/io_units.f90's "the files the controls read opens",
   state/control_slots.f90's `atime` and `lnewvars`, `atime(itime_*)` in
   core/step/timestep_control.f90 and core/model/shell_boundaries.f90),
   and docs/GLOSSARY.md names about 25 rotation files that no longer
   exist and HOW_YREC_WORKS section 4.1 the F77-era EOS file names
   (both stale before the round). A comments-only commit each, objects
   compared.

#### Learned

- Merging a results-changing package in segments, each with its own
  deck run on the merged code, keeps every move attributable to one
  commit: EOS2's three moving fixes were each measured again after the
  merge, and the commits between them shown identical.
- A rename script tested on its own branch can still miss code another
  package adds: STATE2's scope test did not take test_io's module-level
  star for the star. Running each script after all merges, as the last
  round did, caught it at the build.
- A deck's noise floor is worth measuring before a drift is read as a
  response: a 1e-10 relative change of N_A moves the solar deck by
  1e-4 in log L age-matched near the ZAMS, the size of NET2's first
  fix, while N_A x (1 + 1e-6) showed the 3 Msun deck's move to be a
  linear response, 140 times larger than that perturbation's.
- A fix measured before anything is built on it can be rewritten
  cheaply: EOS2's first Debye-Hueckel rule (Ne..Fe scaled by the zone's
  Z) moved the ZAHB deck by 6e-3, because a helium-burning core's C and
  O became Ne..Fe; the commit was rewritten before any later
  measurement.
- A log change can be proved to be nothing else: CTL2 compared each run
  log with the base and accepted only added " WARNING: " lines.
- Paired decks cover what the sets do not reach: OUT2 ran eight (GYRE,
  FGONG and isochrone files, .mod and YREC7 starting models, the
  two-run calibration) before and after each refactor of the writers
  and readers.
- Tracing a recipe by doing it finds what reading it does not: STATE2's
  traces found the math-lib rule, the Makefile's second list of test
  binaries and a gen_star_docs run missing from the steps; NET2's found
  the star layer's literal 15s.
- A package's "In the original F77?" is a claim to check against the
  import, again: the documentation pass changed five verdicts (meval's
  arcs and the conductive opacity's guesses are 2026 defects, not the
  F77's; GYRE's eps pair repeats the F77 pulse file's; the GYRE writer
  is 2026, not 2025; the F77's Debye-Hueckel metals were zone-blind
  too), and two figures (N_A m_u is 1.000143, not 1.000135; LowM_HB's
  missing model is all 28 decks', not one's).
- A mutation that loops forever hangs a mutation script, so both
  numerics scripts now time out at 120 s; an equivalent mutant (one the
  code overwrites later) says nothing about a check's strength.

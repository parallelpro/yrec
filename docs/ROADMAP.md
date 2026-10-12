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
  (modernization); output_columns append_column unbounded vs
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

1. core/henyey_coefficients.f90:464 + stitched_model.f90:448-455,
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
4. setup/rezone.f90:845-852 (inherited hpoint.f "MHP 6/00" block):
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
- **mix -> core/composition_step.f90 -- DONE (86fd404d).**
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
- **core/convection_geometry.f90 -- DONE.** The CZ-base interpolation
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
- **Open: facade-only use for eos, kap, atm, numerics.** Only net/ is in
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

### net/ to 9 (yrec-next, 2026-10-11)

- **DONE, bitwise identical: simeqc and burn_implicitly are pure**
  (ef4c50f1); in net/ only Timmes' neutrino routines and their callers
  are not.
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
  referred to by slot (scale_slot), not by name.
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
- **DONE: mass_change** (5e46b707). core/mass_change.f90's change_mass
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
  ludcmp/lubksb, bsstep) remain a licensing liability
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
- **Open:**
  - Settling converts star%m to Bahcall & Loeb units and back in
    place; the round trip can change its last bit every step. Working
    on a copy would fix it but changes results at round-off; it needs
    a measured commit.
  - Rotation calls back into core (shell_physics, the light-element
    burn) and setup (rotation_stability_setup): the layering runs both
    ways.
  - rotation_scratch_lib holds core's mass-change state (masschg2/3)
    and the commented-out path's tables; they belong elsewhere.
  - burn_and_settle and core's composition_step each size the settling
    sub-steps from the helium settling timescale with near-identical
    code; element_diffusion could own it.

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

**sec-10 item 9 adjudicated -- NOT a bug.** GYRE's mesa_file_m reads
eps_eps_T and divides by eps to get eps_T = dln eps/dln T, and MESA
writes d_epsnuc_dlnT (= eps * dln eps/dln T) into that column. YREC's
eps_total*pulse_dlneps_dlnt is that same absolute derivative
(pulse_dlneps_dlnt is the log-log derivative used in ql_dt). Leave.

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
rotation_stability_setup:482 wants a logarithmic dlneps/dlnT,
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
reversal. Nothing has been edited yet; awaiting the choice of
batch.

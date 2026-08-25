# YREC modernization roadmap -- remaining work

This file lists what is NOT done yet. The completed campaigns
(2026-08-20 .. 2026-08-23) are documented in the git history of this
file and in the commit messages; the one-line summary:

- Phase 1: every COMMON block -> modules/derived types.
- Phase 2: per-domain facades (eos_get/kap_get/atm_get), domain
  folder reorganization, misplaced-file relocations.
- Phase 3 stages 1-3: public/private boundary closed + boundary
  checker; per-domain libraries + standalone byte-pinned tests
  (test_eos/kap/atm/net); ierr-not-stop across the libraries.
- Phase 4/5: the star layer (one type(star_info) `star`), run_yrec/
  evolve_step drivers, star_job + evolve_state, re-entrancy
  (yrec_reset; second in-process run byte-identical), CI.
- Phase 6: state taxonomy finished; physics purity (eos/kap/net
  star-blind); wind reclassified star-layer.
- Legacy campaigns: 504 numbered DOs -> block DO; 249 blanket saves
  -> 10 annotated INTENTIONAL; 698 gotos -> 0. All byte-gated.
- libyrec + pyyrec: embeddable engine (yrec_capi, make lib), ctypes
  binding, CLI-oracle acceptance test.
- MESA conventions: star%job/star%evo nested; star_info members on
  MESA vocabulary (nz, xa, logT/logRho/logP/logR, m, dm, j_rot,
  i_rot, fp_rot/ft_rot, dt, *_start); state/ consolidated to one
  file, domain table state moved into eos//kap//atm/.
- Inlist revamp: &star_job/&controls with readable names (registry-
  driven, defaults/*.defaults generated), legacy decks dispatch to
  the untouched byte-pinned path, converter + byte-identity tests.

Standing verification for everything below: full clean build,
gate2 (dual solar cases, byte-diff vs standard), the standalone
test suite, and for anything near the drivers test_reentry +
test_pyyrec. make clean after any module-TYPE or signature change.

---

## Named-index result arrays (phase-3 stage 4)

`eos_get` has 27 positional arguments; MESA's `eosDT_get` returns
one `res(:)` array indexed by named constants (`i_lnPgas`, `i_Cp`,
...). Adopt the same:

- Index-constant block (`i_pressure`, `i_grad_ad`, `i_cp`, ...,
  `num_eos_results`) in eos_lib; facade packs/unpacks around the
  unchanged eqstat/eqstat2/meqos.
- Migrate the ~10 caller files; each site's local variable soup
  collapses to one array + named indexing.
- Then kap_get (4 outputs) if worth it.
- Byte-identical verifiable throughout (packing the same doubles
  changes no arithmetic). This is also the natural basis for a
  pyyrec in-memory results API (profile/history accessors over
  star%), which is the reason to do it sooner rather than later.

## Numerics-gate ierr opt-in

The numerics_lib procedures carry OPTIONAL ierr gates, but their
callers do not yet pass ierr -- a bsstep/splint failure still stops
the process. Thread ierr from the gates up through the callers to
evolve_step's step_status. Payoffs: embedding safety (a failed
model becomes YrecError in pyyrec instead of killing the worker),
m0030-style configurations become re-enterable, and test_reentry
can use a numerics-terminated case.

## io-writer stops

putstore/wrtmod-family write-path stops -> ierr, same pattern as
stage 3. Small, mechanical.

## controls -> star% campaign (ACTIVE 2026-08-24; supersedes the old
## "const umbrella dissolution" section)

User-approved target shape: nested `star%ctrl` (namelist controls,
immutable after read, structural reset via `star%ctrl =
controls_state()`) + nested `star%job` (paths, run-list arrays, MC
config, nk); everything the code COMPUTES flattens to direct star%
members. phys_const_lib (relocated to state/, Phase 0 DONE) stays
the one legitimately blanket-use module. Pinned rule: mutables never
live in star%ctrl -- a control that seeds a working value (cmixl,
atm_choice) gets a star% working copy.

The namelist wall: Fortran forbids namelist reads into derived-type
components, so controls_lib survives as parmin's read BUFFER. Read
sequence: (1) star%ctrl = controls_state() [pristine defaults];
(2) seed buffer from star%ctrl; (3) parmin reads over it (namelist
only overwrites what the file provides -- semantics unchanged);
(4) store buffer -> star%ctrl. The seed/store copies and the type
body are GENERATED from controls_lib's declarations (single source
of defaults; members lacking an initializer get the explicit
zero/blank/.false. that static storage gave them). This replaces
controls_reset_lib's snapshot-restore outright.

Phases (each byte-gated + test_net + test_reentry):
- Phase 0 (DONE): phys_const_lib -> state/.
- Phase A (DONE 2026-08-24, 6 batches): every non-namelist straggler
  evicted -- atm trio + tenv, cross_section_scale family, iolaol2/
  ioopal2 -> luout_lib, cmixl -> star%mixing_length_alpha (cmixl2/3
  stay: MLT formula constants despite the names), the solar octet,
  nk -> star%job%nk (DO-variable cannot be a component: local
  kind_card drives the loop, post-loop fix-up preserves DO
  semantics). controls_lib is now a pure namelist target set.
- Phase B (DONE 2026-08-24): tools/gen_controls_state.py generates
  controls_state's 418 default-initialized components
  (state/controls_state_def.inc) + the seed/store copies
  (state/controls_sync_lib.f90) from the buffer's declarations;
  star%ctrl added; read_controls runs reset/seed/read/store;
  controls_reset_lib DELETED (test_reentry byte-identical without
  it). Consumers still read the buffer -- star%ctrl is written, not
  yet read.
- Phase C controls stream (DONE 2026-08-24, batches 1-3): all 418
  buffer members classified and migrated -- 346 immutable ->
  star%ctrl (batch 1 mega-sweep); 68 mutable namelist/card members
  -> star%job (batches 2-3: card arrays, calibration protocol,
  model-restore set, driver toggles, physics-adjusted config); 23
  working/diagnostic members -> flat star% (vfc, calcad outputs,
  chkscal bookkeeping, ...). The buffer now has NO production
  readers/writers outside the read path (parmin/remap/sync).
  Read-path invariants learned: remap is read-path (reads+writes
  the buffer); parmin stores buffer->star after remap so
  output_init_mesa sees real values; model-restore readers keep
  their control-named intent(out) dummies bare.
- Phase C flattening stream (DONE 2026-08-24): evidence-driven via
  defaults/flatten_rename_map.tsv (writer/reader analysis per
  member). 16 dead members deleted (user-approved); MC sample
  arrays -> star%job; prev/diag/turnover/light_burn/engeb/flux/
  env_comp/thermo/run/evo dissolved to flat star% (158 members,
  names kept -- MESA micro-renames like mean_molecular_weight->mu
  and the diag cryptic-name modernization deferred to a later
  scriptable pass); pulse's 9 pulse_* physics arrays flattened.
  STAY NESTED as documented solver workspace: rot, mix_phys, circ,
  pulse's q* print scratch -- flattening them awaits the solver
  cleanup (and rot's dm/pm/tm collide with model members).
- Phase D (DONE 2026-08-24): const_lib umbrella deleted (89 files ->
  phys_const_lib, 71 dead imports removed, read path -> controls_lib);
  controls_lib relocated to io/ and declared parmin-private;
  luout_lib -> io/, intpar_lib -> numerics/; const/ folder gone;
  idt/idd dead stores deleted. THE CAMPAIGN IS COMPLETE. Follow-ups
  queued separately: MESA micro-rename pass (map QUERY rows), solver
  cleanup (rot/mix_phys/circ/pulse-residue out of star_info),
  registry-driven MESA renames for ctrl members.

## Library-based yrec link

`make libs` builds 15 per-domain archives but yrec still links the
flat object list. Link yrec (and the test programs) from the
libraries so the boundary is enforced at link time, not just by the
checker script.

## Inlist registry follow-ups

- 8 controls still status=todo in defaults/controls_registry.tsv
  (provisionally named, semantics unverified). Verify and freeze.
- Consider converting the examples/ decks to new-style copies
  (keeping the legacy originals while the pinned baselines read
  them) and documenting the new format in README.
- The .short echo still prints the legacy namelist group for
  new-style runs; acceptable, revisit if user-facing.

## star_info flattening (decided 2026-08-24; folded into Phase C)

User decision: only star%ctrl and star%job stay nested (the two
input bundles -- MESA's actual principle: s% job/s% pg are nested,
physics state is flat). Every other sub-struct (prev/run/evo/
turnover/flux/diag/rot/thermo/circ/mix_phys/engeb/env_comp/
light_burn/pulse) dissolves to direct star% members, fused with the
MESA-vocabulary rename inside the controls-campaign Phase C batches
so each member is touched once.

## atm envelope-integration purity split

core/envint_lib's envelope integration half still reads ~90 star%
references; the tables half of atm/ is already pure. A driver/
kernel split like burn_lib's would finish domain purity. Largest
remaining architecture item, low urgency (cleanly quarantined).

## MHD table data half

No MHD tables ship with the repo; test_eos's MHD section is
skip-unless-YREC_MHD_TABLES infrastructure. Obtaining tables from
outside the repo is the only way to pin the MHD path (including the
starin bug fix and the hsubp/sconvec/massloss MHD extension).

## other/ hook system (optional, unchanged)

MESA-style null-default procedure pointers (other_eos_get, ...)
checked before standard dispatch. Lowest value; only if a concrete
A/B use case appears.

## MESA-style output (in progress; user direction 2026-08-23)

use_legacy_output landed (path-dependent default: legacy decks
.true., new-style inlists .false.) with the first MESA-mode stream:
CASE.history in MESA's history.data layout, .track v0's 83 columns
under MESA vocabulary, numerically interchangeable with the track
(test_mesa_output.py). Remaining, in order:

- DONE (2026-08-23): profiles (profile{N}.data, zone 1 = surface,
  every profile_interval models; the history profile_number column
  replaces profiles.index per user design), column selection
  (history_columns_file/profile_columns_file, one name per line,
  unknown names fatal with the valid list logged),
  star_history_name (default history.data in the output dir), and
  the FGONG pulse writer (io/write_fgong_pulse.f90, ivers 300,
  ICONST 15 / IVAR 40, layout ported from MESA pulse_fgong.f90;
  selected via pulse_format = 'FGONG', GYRE remains the default).
  All in io/yrec_output.f90's registry of names + pure readers over
  star_info; test_mesa_output.py pins the whole contract.
- MESA-style terminal output behind a control; quiet the .short
  config echo / log verbosity control at the same time.
- DONE (2026-08-23): GYRE-HDF5 (GSM) pulse writer behind an optional
  HDF5 dependency (make USE_HDF5=1, MESA SDK static libs); profiles
  AND pulse files now cover the full extended model -- interior +
  re-integrated envelope + atmosphere (build_extended in
  io/yrec_output.f90, stitch's recipe), photospheric M/R/L globals;
  fixed the dead-code stray RETURN that zeroed all 12 MESA history
  star%run% sources (log_R, log_g, inertias, rotation, snu, h-shell).
- Pulse follow-ups: FGONG glob(6) mixing-length alpha and the
  central d2P/d2rho globals are currently zero -- fill if a consumer
  needs them.
- DONE (2026-08-23): output I/O centralized in io/yrec_output.f90;
  MESA mode produces exactly CASE.history + CASE.log (short_file_unit
  retargeted to the log; no legacy opens, no stubs). Special rules:
  LMONTE forces legacy; helium-flash keeps ilast (restore file);
  calibration-summary/rewind chain is legacy-only, so MESA-mode
  history accumulates calibration cycles.
- Log verbosity: CASE.log currently carries the full solver
  iteration traces; add a verbosity control with the terminal step.

## Readability (DONE 2026-08-23)

Four-item program, all landed with byte-identical gates:
double negatives flipped (88 sites, 7 kept where the .not.(...)
conjunction is the semantic condition); named indices for the
positional arrays (488 sites: xa species / del_grad / seg /
luminosity_breakdown / neutrino flux, parameters in star_info_lib);
commented-out F77 code deleted (813 lines, lineage notes kept);
giant procedures decomposed into contains-based named phases --
parmin (6), starin (3), hpoint (5), atm_get (4), engeb (3; its
reaction-rate/screening core stays inline because extraction
perturbs FP scheduling by 1 ulp -- documented in the code). sneut
untouched by design.

## core/ driver program (DONE 2026-08-24)

Phases 1-5 + calibration, all gated byte-identical (three pinned
cases incl. Test_solarcal for the calibration work): inherited
stop-disarm bug fixed; LNUTAB resurrected as compute_neutrino_fluxes;
stop_conditions module (protocol constants, age predicates, D/X/Y
stop table); crrect ladder via solve_level; both drivers decomposed
into named phases (run_yrec 629->432 incl. new docs); MC/.snu writer
-> io/write_run_summaries; 114 duplicate json decls deleted; 20 live
placeholder flags named from verified usage (10 dead carriers keep
the suffix as documentation); calibration protocols (solar triples /
star pairs) stated once at run_yrec's verdict call site with named
cycle constants, end_of_card_calibration extracted, chkcal's stale
pair-protocol header fixed.

## Science wiring (sg-rotation side)

Wire pyyrec into the fit_*_mcmc.py workflow: a YREC rotational-
evolution track producer (multiprocessing, one engine per worker
process), feeding the same fitting machinery that currently
consumes MESA tracks. Depends on nothing above, benefits from the
numerics-gate opt-in (worker survival) and stage 4 (in-memory
reads).

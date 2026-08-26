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

## Named-index result arrays (phase-3 stage 4) -- DONE 2026-08-25

eos_lib carries the index-constant block (i_temperature ...
i_cp_dp, num_eos_results = 24) and `eos_get_r`, which packs the 24
outputs into one intent(inout) res(:) around the unchanged eos_get
(inout slots -- i_log10_density, i_beta, the ionization fractions,
the gradient/cp guesses -- keep their historical carry through the
array). kap_lib gained the matching `kap_get_r` (i_kap,
i_log10_kap, i_dlnkap_dlnrho, i_dlnkap_dlnt; ion_fraction stays an
explicit inout arg), and temperature_gradients_r unpacks both
arrays around the unchanged temperature_gradients.

Migrated (each gate2 byte-identical): observables_lib's central
conditions, compute_scale_height, massloss's two accretion-entropy
sites, then the nine formerly-deferred plumbing sites --
shell_physics, henyey_coefficients, semiconvection x3,
read_starting_model's convective-flag test, atmosphere_derivs,
envelope_derivs, envint's atmosphere start (this last one keeps
its host-associated scalars and wraps the call with a symmetric
prepack/unpack, since the scalars are shared across envint_lib's
contained routines). Every production call site now goes through
the _r facades; the scalar eos_get/kap_get remain only for the
standalone domain tests, and scalar temperature_gradients is
called only by its _r wrapper.

LOAD-BEARING RULE (gate-proven twice): eqstat updates its density
argument IN PLACE, and downstream code -- kap_get, the tail of
henyey_coefficients, semiconvection's next call -- must see the
UPDATED value, not the pre-call one. Every migrated site either
passes eos_res(i_log10_density) onward or writes it back to the
local immediately after the call. The indices are the basis for
the pyyrec in-memory results API.

## Numerics-gate ierr opt-in -- DONE 2026-08-25

numerics_termination (-2) protocol: bsstep failures flow
bsstep -> envint -> surfbc -> henyey_iterate -> evolve_step ->
run_yrec -> main (clean exit 0, matching the legacy stop) and to
pyyrec as a distinct negative status; ksplint/splint/splintd2/
intpol gates wired into their hosts' stage-3 ierr channels.
Residual (documented): opal92 interp chain and the calcad/
tauintnew output-path gates keep absent-ierr stops.

## io-writer stops -- DONE 2026-08-25

parse_columns config stops -> ierr through output_init_mesa ->
read_input; GSM-without-HDF5 rejected at configuration time via
gsm_supported(); residual grid-overflow stops in rebuild_envelope/
read_starting_model converted. Remaining stops by design: main's
exit-1, envint's absent-ierr funnels, qenv's tpgrad residual.

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
  Solver cleanup DONE 2026-08-25: rot/mix_phys/circ moved OUT of
  star_info into rotation/rotation_scratch_lib.f90 (instances
  rot_scr/mix_scr/circ_scr, private-by-default; yrec_reset
  snapshots them alongside star0); their 14 PROPERTY members
  (bl_* scales, MLT alpha vectors, metal_abundance_change, the
  es/ss/gsf circulation velocities) flattened onto star% first.
  Micro-renames DONE same day: del_grad -> gradr/gradT/grada,
  sesum/seg/sbeta/seta/svel/so/locons/sfxion modernized,
  mean_molecular_weight -> mu; scp deliberately kept (fill-time
  distinct from cp, documented). Only pulse's q* print scratch
  remains nested in star_info.
- Phase D (DONE 2026-08-24): const_lib umbrella deleted (89 files ->
  phys_const_lib, 71 dead imports removed, read path -> controls_lib);
  controls_lib relocated to io/ and declared parmin-private;
  luout_lib -> io/, intpar_lib -> numerics/; const/ folder gone;
  idt/idd dead stores deleted. THE CAMPAIGN IS COMPLETE. Follow-ups
  queued separately: MESA micro-rename pass (map QUERY rows), solver
  cleanup (rot/mix_phys/circ/pulse-residue out of star_info),
  registry-driven MESA renames for ctrl members.

## Library-based yrec link -- DONE 2026-08-25

`make yrec_libs` links core/main.o against the 13 per-domain
archives (list repeated for single-pass resolution); an object
missing from its domain library surfaces as an undefined symbol,
which the flat link would paper over. Verified byte-identical to
the flat binary on the solar case; CI builds it as the layout
check. The default `yrec` stays flat-linked (byte-stable object
order).

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

# YREC phase-three roadmap: adopting MESA's module architecture

Written 2026-08-21, immediately after the phase-two domain sweep
finished (see GUIDELINES.md for phases one and two). The reference
architecture studied is MESA 26.04.1 (local copy at $MESA_DIR),
whose every physics module follows one identical skeleton:
`public/` (exactly `<mod>_def.f90` + `<mod>_lib.f90`, the only
surface other modules may touch), `private/` (all implementation),
`defaults/` (namelist controls), `other/` (user hook templates),
`test/` (standalone per-module test program), and `make/` (each
module builds its own static library). MESA even ships a
`package_template/` directory for stamping out new modules.

Phase two gave YREC the *semantic* half of this: one public
`<domain>_lib` facade per domain, domain folders with
pipeline/component subfolders, and derived-type state modules. This
roadmap covers the *mechanical* half -- what actually stops a file
from bypassing the facade, what makes a domain testable on its own,
and what makes errors survivable. Each stage below is independently
shippable and ends with the standing verification (full clean build +
Stage-0 byte-identical regression), plus, once stage 2 lands, the new
per-module tests.

**Explicit non-goal: handles / multi-instance support.** MESA's
`alloc_eos_handle` machinery exists so one process can run multiple
differently-configured stars (binaries). YREC runs one star per
process; the phase-one decision to use single module-level instances
(GUIDELINES.md, `prev_model` discussion) stands. Do not add handles.

---

## Stage 1 -- finish the public/private boundary (fix the bypasses)

The facades exist but nothing routes around-the-side calls through
them. Current inventory of bypasses (verified by grep, 2026-08-21):

1. `atm/turnover/calcad.f90` calls `eqstat2` and `esac06` directly
   (documented exception since the eos sweep). MESA's answer is
   `eosDT_get_component` -- a public accessor for *one named
   component's* results, bypassing the blend but not the boundary.
   Fix: add an `eos_get_component`-style entry to `eos_lib`
   (log T, log P in; the OPAL-2006 or eqstat2 result set out) and
   migrate calcad to it.
2. `wind/massloss.f90` calls `atm/tables/alsurfp.f90` directly.
   Fix: add a small `atm_lib` public entry (surface P/T lookup at
   given Teff/g -- the thing `alsurfp` does) and migrate.
3. `setup/setups.f90` performs three domains' table loads itself
   (`setupopac` for kap, `mhdtbl` for eos/MHD, `alfilein` for
   atm/Allard). MESA's answer is the `<mod>_init` lifecycle entry.
   Fix: add `kap_init`, extend `eos_lib` with an init that loads MHD
   tables when enabled, add `atm_init`; `setups.f90` then calls three
   facade inits instead of three internals. (`setupopac.f90` and
   friends move from de-facto-public to genuinely internal.)
4. `core/starin.f90` and `setup/hpoint.f90` call `kap/surfopac.f90`
   (refresh cached surface-composition table slices). Fix: make this
   a public `kap_lib` entry (`kap_update_surface_tables` or similar)
   -- it is a legitimate lifecycle operation, just unnamed as such.
5. `setup/grsett.f90` / `setup_grsett.f90` call
   `rotation/microdiff/` internals (`lax_wendrof1/2`,
   `get_imp_diffco`, `thdiff`). This is the legacy
   backwards-compatibility diffusion path. RESOLVED (2026-08-21, user
   choice): a third option beat the two anticipated here -- the pair
   is settling physics that lived in `setup/` only by name, so it was
   relocated into `rotation/microdiff/` itself, dissolving the bypass
   (the calls became intra-folder, symmetric with `microdiff`'s own
   use of the same kernels).
6. `core/main.f90` does `use opacity_table_lib` directly -- check
   what it touches and either route through `kap_lib` or drop the
   use-statement if vestigial.
7. State-boundary oddity: the OPAL-2006 EOS files
   (`eos/opal/*06.f90`) store their working arrays in
   `atm_table_lib` (historical accident of the original COMMON
   layout). Move those members to `opal_eos_lib` so eos state lives
   in eos's state module. Pure member relocation, byte-identical
   verifiable.

Deliverable: zero cross-domain calls that do not go through a
`<domain>_lib` facade, verified by a grep-based checker script (see
stage 2 for making that enforcement automatic).

Found during stage-1 execution (2026-08-21): `setup/setups.f90` also
read the **Fermi-Dirac table** (the degenerate-electron EOS table,
former `common/ccr/`) inline, into members that lived in
`atm_table_lib` but are consumed only by `eos/yale/eqrelv.f90` --
the same two-part pattern as the OPAL-2006 oddity. RESOLVED same day:
the `fermi_table_*` members moved to a new `state/yale_eos_lib.f90`
(the Yale/Prather EOS's own state module) and the load into
`eos_init`; the Kurucz/Castelli/Allard surface-table loads (the
inline reads noted under item 3) moved into a full `atm_init` at the
same time; and the SCV EOS table reads (a third inline setups.f90
block, spotted last) moved into `eos_init` as well -- order-safe
since the SCV read shares no file units or state with the atm reads
it now precedes, and directly regression-verified since the Stage-0
OPALSCV cases exercise the SCV tables. **Stage 1 is complete**:
every inventoried bypass and all three execution-time discoveries
are resolved, each verified byte-identical. `setup/setups.f90` now
performs no domain table I/O at all -- it computes constants and
calls `kap_init`/`eos_init`/`atm_init`.

## Stage 2 -- per-module standalone builds and tests

MESA compiles each module into its own static library and ships a
standalone test program per module; the library boundary *is* the
enforcement (private objects simply are not exported).

- Makefile: give each domain a library target (`libeos.a`,
  `libkap.a`, `libatm.a`, `libnuclear.a`, ...), link `yrec` from the
  libraries. The existing recursive-find single-build keeps working
  during the transition; the library split can go domain by domain.
- Add `test/` programs per domain, starting with eos: a small main
  that loads tables, evaluates `eos_get` (and `eos_get_component`)
  over a checked-in grid of (rho, T, X, Z) points, and byte-compares
  against stored expected output. Same pattern for kap
  (`kap_get` over a T/rho/composition grid) and atm (`atm_get` for a
  few (Teff, L) vertices).
- **This is what finally closes the LMHD coverage gap** -- with one
  caveat found during execution (2026-08-21): **no MHD table files
  ship with the repository** (`find` turns up only source/objects
  under `eos/mhd/`), so the gap splits into an infrastructure half
  and a data half. `test_eos` implements the infrastructure half: its
  MHD section runs `eos_get` under `use_mhd_eos` when the
  `YREC_MHD_TABLES` environment variable points at the 8 table files,
  and self-reports SKIPPED otherwise. The pinned baseline covers the
  shipped-data configuration only; actually pinning the MHD path
  needs tables obtained from outside the repo. No Stage-0 case sets
  LMHD, so until then the MHD path of `eos_get` (including the
  starin.f90 bug fix and the hsubp/sconvec/massloss MHD extension
  from the eos sweep) remains verifiable only in structure, not in
  numbers.
- Add the boundary-checker script from stage 1 to the test suite so
  facade bypasses fail CI rather than accumulating again.

STATUS (2026-08-21): stage 2 substantially complete. Delivered:
per-domain library targets (`make libs`, 15 archives; the yrec
executable still links the flat object list -- the library-based
yrec link remains open, transitional by design); the boundary checker
enforced via `pytest test_boundaries.py`; and three standalone test
programs with pinned baselines wired into pytest --
`eos/test/test_eos.f90` (eos_get over an 8-point solar-track grid
with derivatives, both eos_get_gamma1 branches including the
Yale/SCV branch no Stage-0 configuration reaches, MHD as
skip-unless-tables infrastructure), `kap/test/test_kap.f90` (kap_get
over 7 points on the OPAL95/GS98 table), and
`atm/test/test_atm.f90` (all three tabulated atmosphere options --
Kurucz via surfp, Castelli/Kurucz via kcsurfp, Allard via the
atm_get_surface_pt facade -- whose mutual ~0.01-0.05 dex agreement at
matching (Teff, g) doubles as a cross-validation). Open items carried
forward: the MHD data half (above), the library-based yrec link, and
a full `atm_get` envelope-integration test (needs the entire
eos+kap+solver state booted -- effectively a full-model concern,
currently covered by Stage-0; revisit if the envelope integrator
ever needs pointwise pinning).

## Stage 3 -- the ierr-not-stop error discipline

Inventory (2026-08-21): 90 `stop` statements in library domains --
eos 39, kap 29, atm 7, rotation 5, numerics 5, wind 2, mixing 1,
misc 2, nuclear 0. In MESA no library routine ever stops; everything
returns `ierr` and the *application* decides.

- Add an `intent(out) :: ierr` argument to each facade
  (`eos_get`, `kap_get`, `atm_get`, and the new stage-1 entries),
  defaulting to the current behavior at the top-level caller: on
  nonzero ierr, `core/main.f90` prints the context and stops -- so
  the observable behavior of a failing run is unchanged, but the
  library itself becomes embeddable.
- Convert leaf `stop`s bottom-up per domain (kap's table-miss stops,
  eos's out-of-range stops, ...), threading ierr through the
  intermediate calls. Do this per domain, one verified commit each.
- Sequencing note: this stage deliberately comes *after* stage 2
  because Stage-0 regression runs never trigger any of these stops
  (a passing run proves nothing about the error paths). The
  standalone per-module tests can deliberately feed out-of-range
  inputs and assert `ierr /= 0` with no crash -- the first time
  YREC's error paths become testable at all.

STATUS (2026-08-21): **kap converted** -- all 29 library stops became
required-ierr returns in the leaves (each failure's diagnostic still
prints at the point of failure), threaded up to OPTIONAL `ierr`
arguments on `kap_get` and `kap_init`: caller passes ierr -> error
surfaces MESA-style with no stop; caller omits it -> identical
historical behavior, with the stop relocated to a single labeled
funnel per facade (2 stops now stand where 29 did, removable when the
last caller opts in). `kap_update_surface_tables`'s chain has no
error paths and needed nothing. test_kap gained the first error-path
assertions in YREC's history (out-of-table and no-table-chosen both
return ierr=1 without crashing) and its baseline was regenerated.
Execution surfaced a build-system footgun now recorded in
GUIDELINES.md: no `.mod` dependency tracking means any module-
procedure signature change requires `make clean`, or stale callers
segfault. Remaining after kap: eos, atm, and the facade-less
domains (rotation 5, numerics 5, wind 2, mixing 1, misc 2).

STATUS (2026-08-21): **eos converted** -- all 39 library stops (37
live plus 2 that were already dead code behind an unconditional
return, meqos/mhdpx's commented-out F77 999-handling, converted in
place for uniformity) became required-ierr returns in the leaves,
threaded through the full chain (eqstat/eqstat2 -> oeqos*/eqbound* ->
esac* -> readco*/t6rint*, the rhoofp* function trio, and the MHD
chain mhdtbl -> mhdst -> mhdst1 plus meqos -> mhdpx -> mhdpx1) up to
OPTIONAL `ierr` on all three facades: eos_get, eos_init (whose own
inline Fermi-glitch stop joined the funnel), and eos_get_gamma1.
3 stops now stand where 40 did -- one deliberate funnel per facade.
The domain's pre-existing *recoverable* error channels were preserved
untouched and are now clearly distinguishable from fatal errors: the
esac* out-of-table alternate returns (labels 62/65), oeqos*'s *998
fall-back-to-Yale exits, and rhoofp*'s -999 sentinel are recoverable
by design; ierr covers only the corrupt-data/misconfiguration class
that used to stop. Because of that design, healthy tables cannot
drive the facade to ierr /= 0, so test_eos asserts both halves:
the success path threads ierr = 0 end to end through the facade, and
a white-box call into esac06 with an invalid rad_flag returns ierr=1
with no crash. Remaining: atm (7), then the facade-less domains.

STATUS (2026-08-21): **atm converted** -- all 7 stops: the table
leaves surfp/kcsurfp (out-of-table below logTeff 3.5 / logG -0.5,
historically always fatal here, unlike the recoverable channels
elsewhere), alsurfp's 9999 fatal exit (its lookup_failed flag stays
the recoverable channel), and the Allard load chain
alfilein -> altabinit gained required ierr; atm_get's own two
integration-failure stops (atmosphere and envelope MAXSTP exhaustion)
became jerr -> funnel. All three facades -- atm_get, atm_init,
atm_get_surface_pt -- carry the OPTIONAL ierr contract; 3 funnel
stops stand where 7 did. atm_get's internal calls into eos_get /
kap_get deliberately do NOT yet pass ierr: each domain's funnel
preserves the historical stop for callers that haven't opted in, and
threading the cross-domain calls is the natural follow-on once the
solver itself can consume ierr. test_atm asserts surfp's
out-of-table error returns ierr=1 without crashing and that the
facade success path threads ierr=0 (its own iowr now points at the
scratch .short file so the diagnostics stay off the byte-compared
stdout). Remaining: the facade-less domains (rotation 5, numerics 5,
wind 2, mixing 1, misc 2), where the surfacing decision is per
public entry rather than per facade.

STATUS (2026-08-21): **facade-less domains converted -- STAGE 3
COMPLETE.** Three surfacing patterns, chosen per the call-graph
reality:

1. *numerics* (5: ksplint, bsstep, intpol, splint, splintd2): the
   module's explicit interfaces allow the optional-ierr trick without
   a facade, so each procedure carries its own gate at the failure
   point -- ierr present: error return; absent: the historical stop
   stands. Zero caller churn; callers opt in per call site.
2. *the rotation/wind/mixing/misc cluster* (10 counted + 2 found
   during execution): every one of these plain-external routines
   funnels through just two public entries -- getw (rotation) and mix
   (mixing) -- into four driver files. So required ierr threads the
   whole graph (checkc/checkj/bandw/dadcoeft -> seculr;
   wind -> mwind, cowind -> mcowind; rotmix <- bursmix/midmod;
   simeqc <- kemcom <- rotmix/mix, where kemcom's own uncounted stop
   surfaced during execution; tpgrad <- physic/coefft/sconvec;
   getw and mix on top), and the four driver call-site files --
   core/main, core/crrect, core/starin, setup/hpoint -- preserve the
   historical stop via `if (jerr /= 0) stop`. This is the full MESA
   discipline for these domains: the stop is now driver policy, not
   library behavior. One documented residual: atm/qenv.f90's tpgrad
   call keeps a local stop because qenv's signature is fixed by the
   bsstep integrand-callback protocol (same class as qgauss's
   hard-coded call into rotation func).
3. *eos stragglers found during execution*: setup/rtab.f90 (2 stops)
   and util/rabu.f90 (1 stop, uncounted) are the MHD table-read
   helpers called only from eos/mhd/mhdst1 -- they joined the eos
   chain (required ierr into eos_init's funnel). By the
   misplaced-file test both belong under eos/mhd/; relocation left
   as a follow-up candidate, not done here.

Final stop inventory in library code: 13 = the 8 facade funnels
(eos 3, kap 2, atm 3) + numerics' 5 per-procedure gates, every one
behind an opt-in ierr; plus qenv's documented callback residual.
Everything else that stops is driver code (core/, io/, setup/'s
hpoint/rscale/parmin), where stopping is policy, which stage 3
deliberately does not touch.

## Stage 4 -- named-index result arrays

`eos_get` currently has 27 positional arguments; MESA's `eosDT_get`
returns one `res(:)` array indexed by named constants (`i_lnPgas`,
`i_Cp`, `i_gamma1`, ... in `eos_def`). Adopt the same:

- Add the index-constant block (`i_pressure`, `i_density`,
  `i_grad_ad`, `i_cp`, ... plus a `num_eos_results` count) to
  `eos_lib` (YREC's analogue of `eos_def` -- a separate `_def` file
  is optional; the constants matter, the file split does not).
- New facade signature: inputs stay explicit (T, P, composition),
  the ~20 outputs collapse into `res(num_eos_results)` (+
  `d_dlnT(:)`/`d_dlnP(:)` arrays for the derivative outputs).
  `eqstat`/`eqstat2`/`meqos` keep their historical signatures; the
  facade packs/unpacks -- same thin-wrapper principle as phase two.
- Migrate the ~10 caller files; each site's local variable soup
  (`beta, beta_inverse, beta14, ...` declared in every caller)
  collapses to one array + named indexing, which is the real
  readability win.
- Then the same for `kap_get` (smaller: 4 outputs) if it proves
  worth it there.
- This is the churniest stage (every call site rewritten), which is
  why it comes after stage 2's tests exist. It remains byte-identical
  verifiable throughout: packing/unpacking the same doubles changes
  no arithmetic.

## Stage 5 (optional) -- the other/ hook system

MESA's `other/` modules hold null-default procedure pointers
(`s% eos_rq % other_eos_component => my_eos`) so users override
physics without editing the library. For YREC this is the lowest-
value item: the user base edits the source directly, and without
handles the pointers would live as module-level state. Only build it
if a concrete use case appears (e.g. an experimental EOS someone
wants to A/B without forking eqstat2). If built: one procedure
pointer per facade (`other_eos_get`, `other_kap_get`), checked
before the standard dispatch, null by default -- behavior-identical
when unset.

---

## Sequencing rationale and cross-cutting notes

Order is 1 -> 2 -> 3 -> 4 (-> 5) because each stage makes the next
one safer: 1 gives every cross-domain interaction a named facade
entry, so 2's library boundaries have nothing to strand; 2's
standalone tests are what make 3's error paths and 4's signature
churn verifiable beyond Stage-0's happy path; 4 is the largest
mechanical churn and should not go first on that ground alone. 5 is
optional and independent.

The deferred blanket-`SAVE` cleanup (GUIDELINES.md, "Deferred:
unneeded blanket SAVE statements") is independent of all five stages
and can interleave anywhere; the natural moment for each file is
whenever a stage already has that file open for edits.

Verification discipline is unchanged throughout: full clean build +
Stage-0 byte-identical diff for every commit, per-module tests
additionally from stage 2 on, and the boundary-checker keeps stage
1's win locked in. As with phase two, investigate before assuming:
each stage's first commit should re-verify this roadmap's inventory
against the then-current source rather than trusting these counts.

---

# Phase four: the star layer (star_info)

Added 2026-08-21, after stage 3 completed, from a measured
investigation of the non-physics folders (core, io, misc, mixing,
nuclear, numerics, rotation, setup, util). Question asked: are they
modular or intertwined? Answer: the physics domains that phases 1-3
worked on are now genuinely modular -- entered only through facades,
calling only downward, owning their state, standalone-tested. The
star-level folders are not, and the entanglement has one root cause:
**the stellar model has no single owner or representation.**

## Evidence (2026-08-21 dependency scan)

Cross-domain call counts (rows call into columns; facade calls
included -- the point is which edges exist at all):

```
          atm  core eos  io   kap  misc mix  nuc  num  rot  setup util wind
atm       -    .    4    1    3    1    .    .    42   .    .     .    .
core      3    -    1    9    2    6    4    4    .    8    10    3    1
eos       .    .    -    .    .    .    .    .    33   .    5     4    .
io        4    .    1    -    .    2    .    .    3    1    .     .    .
kap       .    .    .    .    -    .    .    .    62   .    .     3    .
misc      1    .    2    .    2    -    .    1    1    2    1     .    .
mixing    .    .    4    .    3    3    -    7    .    5    .     4    .
nuclear   .    .    .    .    .    .    .    -    2    .    .     .    .
numerics  .    .    .    .    .    10   .    .    -    2    .     .    .
rotation  1    .    .    2    .    14   6    3    30   -    2     .    5
setup     1    .    1    .    2    2    1    2    31   7    -     .    .
wind      1    .    2    .    .    3    .    .    .    .    .     .    -
```

The four concrete mechanisms of entanglement:

1. **The model is smeared across ~40 arrays declared in
   `program main`** (`composition(15,json)`, `log_mass`, `omega`,
   `moment_of_inertia`, ..., `json=5000`) **and threaded through
   giant positional argument lists**: `crrect` 60 arguments,
   `starin` 50, `midmod`/`wrtout` 43, `seculr` 37, `hpoint` 35,
   `engeb` 35, `coefft` 33, `getw` 29. Every signature is a
   hand-maintained slice of the model vector; every intent bug
   phases 2-3 surfaced (3 of them) was an argument-list mismatch.

2. **The other half of the model lives in ~10 shared mutable
   state modules** (former COMMONs in `state/`): `oldmod_lib` IS the
   previous model (`old_pressure`, `old_temperature`,
   `old_composition`, ...); `rotdiff_lib`, `turnover_lib`,
   `run_diag_lib`, `scrtch_lib` are each read by 6-9 domains. Two
   representations, no owner.

3. **`misc/` is a grab-bag, not a domain**, manufacturing fake
   dependency edges (rotation->misc 14, numerics->misc 10):
   `tpgrad` is convection/atm physics, `coefft` is the Henyey
   matrix builder (core), `physic` is the all-zones thermo driver
   (core), `solid`/`shape` are rotation geometry,
   `spline`/`splinnr`/`slopes`/`search`/`choose`/`cases`/`simeqc`
   are numerics, `chkcal`/`chkscal`/`pdist` are core diagnostics,
   `stitch` is a core/model operation.

4. **Cycles and layering violations**: rotation<->mixing is a true
   cycle (`getw` -> `bursmix`/`mixcz`/`convec`; `rotmix`/`mix` ->
   `microdiff`/`grsett`/`ndifcom`) -- rotational mixing genuinely
   couples them. `io` COMPUTES at write time (`wrtout` ->
   `calcad`/`gettau`, `wrtmod`/`getnewenv` -> `atm_get`), and
   rotation calls io's `model_to_equal`/`equal_to_model` to copy
   model state. `numerics` is not a leaf (`qgauss`'s hard-coded
   call into rotation's `func`; `meval` -> misc's
   `choose`/`cases`/`search`). `setup/hpoint` and `setup/midmod`
   are solver orchestration living under a "setup" name.

## Target architecture

Three layers (MESA's shape, adapted to the single-star,
no-handles decision):

- **Drivers** (`main`, and the stage-3 driver-side stop sites):
  policy -- when to step, when to stop.
- **Star layer** (core solver, setup's stepping routines, io,
  mixing, rotation): owns `type(star_info) :: star` -- the model
  arrays, one instance, modified in place through the run. The
  rotation<->mixing cycle stops being an architectural problem
  here: both are star-layer routines operating on the same struct,
  which is exactly how MESA treats rotational mixing.
- **Physics services** (eos/kap/atm/nuclear/wind facades +
  numerics as a pure leaf): take plain arguments, return results,
  NEVER see star_info. This is MESA's boundary too (star_def is
  unknown to eos/kap/net), and it preserves the standalone
  testability stage 2 built. The phase-one single-instance
  decision stands: one module-level `star`, no handles;
  `oldmod_lib` becomes simply a second instance (`prev`).

## Steps (each byte-identical-verifiable, in order)

1. **Dissolve `misc/`** -- pure `git mv` relocations per the
   phase-two misplaced-file test: `tpgrad` -> atm (or a convection
   home decided on investigation), `coefft`/`physic`/`stitch`/
   `chkcal`/`chkscal`/`pdist` -> core, `solid`/`shape` -> rotation,
   `spline`/`splinnr`/`slopes`/`search`/`choose`/`cases`/`simeqc`
   -> numerics. Also `rtab`/`rabu` -> eos/mhd (queued since stage
   3). Removes the fake edges so the real graph is visible.
2. **Make numerics a pure leaf**: the misc spline helpers move in
   (step 1); `qgauss`'s hard-coded call to rotation `func` becomes
   a procedure dummy argument (F2003 procedure interface), passed
   by rotation at the call site. After this, numerics' column in
   the matrix has no outgoing edges.
3. **Introduce `star_info`** (`state/star_info_lib.f90` or a new
   `star/` home): fields = exactly `program main`'s model arrays,
   same names, same shapes. `main` declares `star` and passes it;
   convert the worst signatures top-down -- `crrect` (60),
   `starin` (50), `midmod` (43), `wrtout` (43), `hpoint` (35),
   `getw` (29), `mix` (21) -- each routine replacing its model-
   slice arguments with `star`, one routine per commit,
   byte-identical each time. Non-model arguments (tolerances,
   flags, per-call scalars) stay as arguments.
4. **Fold the `state/` model modules into star_info**:
   `oldmod_lib` -> `type(star_info) :: prev` (model save/restore
   becomes a struct copy); `scrtch_lib`/`run_diag_lib` ->
   `star%diag...`; `rotdiff_lib` -> `star%rot...`; decided
   member-by-member with the same relocation discipline as the
   OPAL-2006/Fermi state moves in stage 1.
5. **Make io a pure reader**: the quantities `wrtout`/`wrtmod`
   recompute at write time get computed in the star layer, stored
   in `star`, and io only formats. Kills the io -> atm edges and
   the rotation -> io state-copy calls.

Relation to the phase-three stages: stage 4 (named-index result
arrays) concerns the physics facades' result vectors and is
orthogonal -- it can interleave; stage 5 (other/ hooks) remains
optional and independent. The blanket-SAVE cleanup interleaves as
before -- note `program main`'s arrays and the star-layer SAVEs are
exactly the files steps 3-4 will have open.

Verification discipline unchanged: full clean build + Stage-0
byte-identical diff per commit, standalone tests, boundary checker;
re-verify this section's inventory against the then-current source
before each step's first commit.

STATUS (2026-08-21): **steps 1-4 COMPLETE, step 5 partial.**

- Step 1 done: misc/ dissolved (17 relocations incl. spline/splinnr,
  dead but kept per the standing decision), kemcom -> mixing,
  rtab/rabu -> eos/mhd; tpgrad and solid added to the boundary
  checker's mixing/rotation allowlists as deliberate public entries.
- Step 2 done: qgauss takes its integrand as a procedure dummy
  (fpft passes func); numerics' matrix row now has zero outgoing
  edges.
- Step 3 done: star_info introduced; main owns no model arrays;
  crrect 60->21 args, starin 50->23, midmod 43->31, wrtout 43->19,
  hpoint 35->12, getw 29->10, mix 21->10. Two separate-storage traps
  found and preserved as arguments: mix's
  mixed_zone_bounds_no_overshoot (crrect passes its own local) and
  getw's local radiative_zone_bounds at the midmod call. The first
  of these was initially missed and caught ONLY by the Stage-0
  byte-diff (a 1e-9 drift in the solar cases; the short m0030 run
  stayed identical) -- the standing lesson: argument-count audits at
  every call site are necessary, and the byte-diff is the last line
  of defense that actually caught it.
- Step 4 done: prev_model -> star%prev, shell_diag -> star%diag,
  run_diag -> star%run, rot_diff -> star%rot (61 files). Physics
  internals (atm, nuclear, wind) that read these former COMMONs now
  reference star% visibly -- the remaining physics-domain
  star-coupling is grep-able, deliberately.
- Step 5 SUBSTANTIALLY COMPLETE (2026-08-21, second pass): the
  state-computing blocks moved out of wrtout into
  core/update_output_diagnostics.f90, called by main immediately
  before wrtout: the luminosity-breakdown renormalization (which
  MUTATES the model), the core-CZ mass (with the preserved FX/FX2
  stale-SAVE bug carried along intact), the central-conditions
  eos_get evaluation, and the surface-CZ base interpolation -- all
  results now stored in star%run% (8 new fields) and only READ by
  the writers. Three documented residuals remain in io: (1) gettau
  stays at its original wrtout position because its atm-side
  integration prints progress lines into the .short stream --
  hoisting it reorders the stream (values identical, layout not), a
  fact discovered by the byte-diff, so it is blocked by print
  interleaving, not data flow; (2) wrtmod/putstore's print-mode
  envelope integration (the atm profile printer -- io orchestrates,
  atm computes-and-prints); (3) wrtout's own display-only
  derivations (SNU sums, H-shell locations, moments of inertia,
  rotation period), which are history-column evaluations in the
  MESA sense and legitimately live with the writer. Also from the
  first pass: getnewenv (model construction -- moves the outer
  fitting point under mass loss -- misfiled in io/) relocated to
  core/. Original deferral analysis, kept for the record: wrtout is not a
  formatter with two stray calls -- it has a compute PROLOGUE (a
  central-point eos_get evaluation plus the CZ-base interpolation
  that fills star%run%envelope_mass/envelope_radius/central_* and
  the envelope_cz_* track columns) whose values feed gettau and
  calcad mid-routine, so the pure-reader refactor means extracting
  that prologue into a star-layer update_output_diagnostics that
  stores its results in star (star%run has most fields; the
  envelope_cz_* column values need new ones) with wrtout reading
  only. wrtmod/putstore's atm_get calls are the atm profile printer
  running in print mode -- io orchestrates, atm computes-and-prints;
  arguably correct as is. Neither was attempted mechanically after
  the step-3 lesson; both are bounded, described here, and next.


FOLLOW-ONS COMPLETED SAME DAY (2026-08-21/22):

- **Model scalars absorbed into star_info** (the deferral in
  star_info_lib's original header, resolved by the per-call-site
  audit it called for): num_zones, model_number, core_cz_top_index,
  envelope_cz_bottom_index, log_total_mass, total_mass_msun,
  log_teff, log10_luminosity. Kept as arguments (separate storage at
  some call site): mix's core_cz_edge/envelope_cz_edge (crrect
  passes its own locals, mix writes them) and
  mixed_zone_bounds_no_overshoot. Signatures now: starin 15, crrect
  17, wrtout 10, mix 7, hpoint 6, getw 5 (from 50-60 at phase-four
  start); update_output_diagnostics takes no arguments. Call-site
  slimming is now positional with per-site argument-count
  assertions (the step-3 lesson, made structural).
- **CI** (.github/workflows/ci.yml): build + the three standalone
  table tests + boundary checker on every push/PR; sparse blob-less
  checkout pulls only src/ + ~17 MB of tables. First runs surfaced
  and fixed two real portability gaps: gfortran 13 hard-errors on
  >132-column free-form lines (now -ffree-line-length-none), and
  the macOS/arm64-pinned baselines differ from x86_64 in the last
  printed digit -- amplified to ~1e-11 relative by the eos grid's
  differenced second derivatives -- handled by a CI-only 1e-9
  relative tolerance (conftest.compare_output; local runs remain
  byte-exact, structure and non-numeric tokens exact everywhere).
  The Stage-0 byte regression remains a local/manual gate (its
  standards are largely untracked).

Next milestone (per the phase-four discussion): the embeddable
evolve-one-star entry -- run_star(controls, star, ierr) -- built on
stage 3's ierr plumbing and this phase's star_info; the driver-side
`if (jerr /= 0) stop` sites are the seams. Folds in the
library-based yrec link and the SAVE cleanup (module state must be
resettable for re-entrancy).

STATUS (2026-08-22): **phase A (first carve) done.** The entire body
of `program main` -- controls read (parmin), setup, and the
Monte-Carlo/run loop with the full evolution driver -- moved verbatim
into core/run_yrec.f90 as `subroutine run_yrec(ierr)`; main is now a
ten-line CLI wrapper (`call run_yrec(ierr); if (ierr /= 0) stop 1`).
The blanket `save` in run_yrec is load-bearing and documented:
program-unit variables were implicitly static, so the moved locals
(data statements, rescale bookkeeping) keep exactly their old
semantics -- and, equally, run_yrec is NOT yet re-entrant; that is
phase C's job. ierr is wired but always 0 in phase A: the body's own
stops and the stage-3 driver seams remain, and converting them into
error returns is phase B. The finer split (read_controls /
star_setup / evolve_step) deliberately waits for phase B: the setup
section's locals (table paths, mixture weights) and the evolution
loop's label structure (15/100/110/200/250) want a controls
structure ("star_job") to carve against, not another blind
mechanical cut -- design that belongs with the ierr work.

STATUS (2026-08-22): **phase B done.** Every driver-layer fatal
condition returns through run_yrec's ierr, and the CLI wrapper exits
1 -- the first time a pipeline can detect YREC failure by exit
status (bare `stop` exits 0). Converted: parmin's config check,
rscale (3), getyrec7's model-glitch check, starin (2 stops + the two
formerly-silent `stop 9998`s, which now print), hpoint (2), the
stage-3 jerr seams in starin/crrect/hpoint/run_yrec, and run_yrec's
own rotation-config stops; its end-of-job stop became a clean
return. update_output_diagnostics threads ierr into its eos_get
call. New end-to-end test (test_error_exit.py, local -- needs the
full input/ tree): a config error exits 1 cleanly with the
diagnostic still printed. Remaining stops below run_yrec: the
physics facades' opt-in funnels (callers that pass ierr already
bypass them -- the driver chain now does), numerics' gates, qenv's
callback residual, and io's writer stops (putstore/wrtmod family),
which are next if full no-stop operation is wanted before phase C.
The star_job controls structure and the evolve_step split remain
the open design work of this phase.

STATUS (2026-08-22, later): **star_job + the setup split done.**
state/star_job_lib.f90 holds the job-level configuration -- the ~26
table/file paths parmin fills and setups consumes (formerly run_yrec
locals threaded through both calls), the LAOL mixture-weight array,
and the Monte-Carlo run-range bounds -- as one module-level `job`
instance (the star_info single-instance pattern). run_yrec's setup
section became two star-layer routines: read_controls(ierr) (parmin
into const_lib targets + job paths) and star_setup(ierr) (setups +
the MC per-run parameter read). setups itself gained required ierr,
threading the kap_init/eos_init/atm_init facades' optional ierr --
a missing or corrupt table at startup is now an error return through
run_yrec rather than a stop in a facade funnel; the three standalone
test programs pass and check it too. run_yrec now reads: declare ->
read_controls -> star_setup -> the run loop. Remaining in this
phase: the evolve_step split (wants an evolution-state inventory of
the loop's ~dozens of saved locals -- same method as the star_info
scalars audit), then phase C (reset/re-entrancy).

STATUS (2026-08-22, later still): **evolve_state + evolve_step
done.** state/evolve_state_lib.f90 holds the 26 scalars the model
loop genuinely shares with the run/Monte-Carlo levels (instance
`evo`), completing the state taxonomy: star_info (model) / star_job
(job config) / evolve_state (driver working state) / const_lib
(physics controls). Measured twice: the first audit counted
declaration lines as outside-uses and swept 28 loop temporaries into
the struct; the Stage-0 byte-diff caught the resulting 1e-16 codegen
drift in the per-zone prediction arithmetic, a 4-build variable-set
bisection isolated the scratch delta_*_step temporaries, and the
corrected executable-lines-only audit gave the clean 26/28 split.
core/evolve_step.f90 then took the model-loop body verbatim (the
divergence/redo cycle, corrector, mixing, rotation, rezoning, output,
timestep update, stop criteria) with a step_status contract
(0=continue, 1=run finished, 2=leave run loop; historical goto
110/200 exits retarget local labels) and the 28 region-only locals
under a load-bearing save. run_yrec's model loop is now four lines.
nk stayed const_lib module state (former common/zramp/), not an
argument. Phase C (reset/re-entrancy) is next: its reset surface is
now enumerable -- evolve_state, evolve_step's saved locals, star,
and the domains' table/lazy-init state.
---

# Phase six: one state taxonomy (finishing the MESA shape)

Added 2026-08-22 after the phase-five evolve_step extraction, from a
survey of the eleven remaining shared state/ modules plus const_lib.
Design question (user): group run-lifetime state under state/ and
move structure/current-step variables into star_info? Sharpened by
the survey into ownership + lifetime, with per-module verdicts:

| module (instance)             | contents                            | verdict |
|-------------------------------|-------------------------------------|---------|
| envelope_comp (env_comp)      | surface X/Z, mixture, stotal/senv   | star%env_comp |
| turnover (turnover)           | convective turnover time, pphot     | star%turnover |
| light_burn (light_burn)       | Li/Be/D burn rates + prev trackers  | star%light_burn |
| fluxes (flux_diag)            | neutrino fluxes, SNU rates          | star%flux |
| engeb_diag (engeb_diag)       | nuclear-energy diagnostics          | star%engeb |
| mdphy (mix_phys)              | per-zone mixing/diffusion physics   | star%mix_phys |
| temp (shell_temp)             | per-zone cp/mu/qdt/thdif/visc       | star%thermo |
| temp2 (circ_vel)              | per-zone circulation velocities     | star%circ |
| pulse_diag (pulse_diag)       | pulsation-output diagnostics        | star%pulse |
| envstruct (env_struct)        | envelope-integration profile        | stays: atm-owned (like atm_table) |
| atmstruct (atmo_struct)       | atmosphere-integration profile      | stays: atm-owned |

A correction from the survey: temp/temp2 looked like dissolvable
scratch by name, but their arrays are cross-routine carriers (physic
fills them, rotation/microdiff read them) -- per-step derived model
state, so they fold into star_info rather than into locals. The
genuinely-scratch lesson (delta_*_step, phase five) stands: scratch
lives in locals, never in structs.

Steps:
1. Fold the nine model-state modules into star_info (step-4
   mechanics: instance moves in as a component, type stays; users'
   `inst%` -> `star%comp%`). Physics-domain files that read these
   (eos's env_comp use is the deepest) gain visible star% references
   -- the same deliberate, grep-able coupling as the step-4 folds.
2. envstruct/atmstruct stay standalone modules, reclassified as
   atm-domain state alongside atm_table_lib.
3. Split const_lib -- the real remaining mess, three things in one
   name: physical constants (immutable), the ~200 namelist targets
   (controls, read-only after parmin), and stray mutable working
   state (debye_huckel_*, the zramp rescale state including nk, the
   SCV EOS tables, which belong in an scv_eos state module). Target:
   phys_const_lib (parameters) / controls as part of star_job (the
   MESA inlist->s% pattern) / strays to their owners. This is what
   makes "read-only after setup" checkable, and it defines the
   reset boundary for phase C (controls survive a reset, working
   state does not).

End state: everything about a run reachable from four roots -- star,
job, evo, and true constants -- physics domains opaque behind their
facades, and phase C's reset implementable by enumerating those
roots. Verification discipline unchanged (the phase-five codegen
lesson applies doubly here: eqstat2's hot loops read env_comp).

STATUS (2026-08-22): **PHASE 6 COMPLETE.**
- Step 1 done (commit 7594654): nine model-state modules folded into
  star_info; 57 files; module-library Makefile prerequisites updated
  to match their actual star_info use (kap_lib raced star_info.mod
  under -j4 until then).
- Step 2 done by classification: envstruct/atmstruct recorded as
  atm-domain state, peers of atm_table_lib; no move needed.
- Step 3 done: const_lib is now a compatibility umbrella
  (use-and-reexport of the two new homes; every existing
  `use const_lib` unchanged):
    * const/phys_const_lib.f90 -- const1/2/3 physical constants and
      the version string (cmixl's per-kind-card caveat documented);
    * const/controls_lib.f90 -- every namelist-target block (~60
      former COMMONs) with original commentary, documented as
      read-only after read_controls, with the mixed-in stragglers
      (zramp's nk run index, computed members, working scalars)
      listed as the remaining member-level cleanup;
    * evictions to owners, NOT re-exported (users reference the
      domain modules directly): Kurucz/Castelli surface-pressure
      tables -> atm_table_lib; the LAOL pure-Z table block ->
      opacity_table_lib; the SCV EOS tables -> new
      state/scv_eos_lib.f90; the Debye-Huckel working members ->
      yale_eos_lib. Control flags that rode along in those blocks
      (use_pure_z_table, use_scv_eos, use_debye_huckel_correction)
      currently live with their tables -- functional, documented
      stragglers for the member-level pass.
  Execution facts: the model-file I/O routines (getyrec7/getmodel2/
  putmodel2) have DUMMY ARGUMENTS named after control flags
  (use_pure_z_table...); a blanket `use` of the new homes collides
  (gfortran rejects a dummy redeclaring a use-associated name), so
  those files keep only-clause uses -- the reason they had
  only-clause const uses all along. An only-clause through the
  umbrella (use const_lib, only: x) still works, so re-export is
  fully transparent.
End state reached: star / job / evo / phys_const as the four roots,
controls_lib as the job's control surface behind the const_lib
umbrella, physics tables in domain state. Remaining (recorded, not
phase-6): member-level eviction of the documented stragglers, and
per-file migration off the umbrella.


---

STATUS (2026-08-22, night): **phase-5 step C COMPLETE + nuclear test.**
Re-entrancy achieved (commit d5798fd): run_yrec is re-enterable via
pristine snapshot (core/yrec_reset.f90 + the GENERATED
const/controls_reset_lib.f90 mirroring every control AND phys_const
member), with flag-gated SAVEd-local re-initialization in evolve_step
/ update_output_diagnostics and a unit-closing sweep. Acceptance test
(test_reentry + test_reentry.py): the second of two in-process calls
is byte-identical to a fresh process, on the solar case. Two findings
the test forced: (1) cases that historically TERMINATE via opt-in
numerics-gate stops (m0030's BSSTEP ending) kill the process mid-call
and are not re-enterable until their callers pass ierr -- documented
residual; (2) cmixl (the const3 caveat member) leaked run-1's
calibrated mixing length into run-2's startup echo until phys_const
joined the snapshot. Nuclear (commit 3e8b544): module map added and
the fourth per-module standalone test (rates/sneut/azbar/deutrate +
remap's SFII cross-section scales) pinned, in `make tests` and CI.
Next milestone unchanged: libyrec + the Python binding.

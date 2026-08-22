# YREC module-modernization guidelines

This document records the ground rules developed while converting YREC's
COMMON-block-based Fortran into modules. It exists so that future work
(by a person or an agent) on any not-yet-converted folder follows the same
reasoning instead of re-deriving it, and so a not-yet-converted file's
treatment is decided by its actual coupling, not by which folder it
happens to live in.

## The two kinds of data hiding behind COMMON, and their two treatments

Before converting any file, classify what its COMMON block(s) actually
carry. There are two genuinely different cases, and they call for
different fixes -- do not default to one for everything.

**1. Global constants / configuration** -- set once (typically at
startup or namelist-read time), read broadly, never really an
"input" or "output" of any one call. Examples: physical constants
(`common/const1/`, `common/const2/`), I/O logical unit numbers
(`common/luout/`), broadly-shared timestep-limit parameters
(`common/ctlim/`).

  -> Convert to a **module holding module-level variables**, named
  `<domain>_lib` (see naming, below). Every file that declared the old
  COMMON block switches to `use <domain>_lib` instead. This requires
  **no argument-list or call-site changes anywhere** -- it is a pure
  declaration swap, mechanically safe regardless of how many files are
  touched. This matches MESA's own practice (`const_def`, `chem_def`,
  etc. are modules with shared module-level values, not threaded
  through argument lists).

  Do not be deterred by a wide closure (a block used by dozens or
  nearly 100 files) -- widen the batch, not the technique. Threading a
  "constant" through every intermediate argument list between where
  it's set and where it's used (including subroutines that don't use
  it themselves, just relay it) is *more* invasive than a module swap
  and buys no extra safety.

  **1b. Genuinely evolving global state (not constant, but still
  module territory).** Some blocks change every model and are written
  from several distinct places rather than once at startup (example:
  `common/oldmod/`, the previous-model snapshot, written from
  `core/starin.f90`, `core/main.f90`, and `setup/hpoint.f90`). This is
  *not* case 2 -- nothing here flows cleanly from one caller to one
  callee within a single call sequence, it's read and written from
  many distant, unrelated points in the call graph, exactly like a
  constant is, just mutable. The mechanical conversion is exactly as
  safe as case 1 (a pure storage-association swap, still verified
  byte-identical) -- don't let "this is real physics state, not a
  constant" talk you out of the module treatment; that's about the
  data's *nature*, not the conversion's *risk*. The one added cost is
  verification thoroughness: trace *every* write site (not just one,
  the way a namelist-read case has exactly one writer) before
  converting.

  For this sub-case, consider a **real derived type** instead of bare
  module variables (MESA's `prev_star_info` pattern) rather than
  defaulting to the bare-variable style: `type :: prev_model_state`
  with the block's members as components, plus one module-level
  `type(prev_model_state) :: prev_model` instance. Two reasons this is
  worth the extra work here specifically: the `prev_model%` prefix
  makes it unambiguous at every use site that a field is *previous*-
  model state rather than some current-model quantity that happens to
  share a name, and bundling the fields into one type makes later
  iteration/serialization (profile/history/pulse output) much easier
  than chasing bare module variables one at a time. This does not
  imply multi-instance support the way MESA's real `star_info` needs
  (YREC runs one star at a time) -- it's a single module-level
  instance, same usage pattern as the COMMON block was, just through
  named/typed fields instead of bare positional slots. The real cost
  is that every read/write site needs a `%` inserted, not just the
  declaration -- budget for rewriting the whole file body, not a
  one-line swap. Bare module variables remain the default for
  genuinely constant/config data (case 1 proper); reach for the type
  only when the data is itself evolving state naturally suited to
  being bundled and later iterated.

**2. Real per-call data flow currently smuggled through COMMON** --
varies call to call, and is genuinely an input or output of a specific
computation, just hidden in a side channel instead of the argument
list. Example: `common/tridi/` (the caller fills in a tridiagonal
matrix's diagonals, calls the solver, reads the solution back out --
all via COMMON rather than arguments).

  -> Convert to **explicit subroutine arguments** (`intent(in)` /
  `intent(out)` as appropriate). This makes the real data dependency
  visible at the call site instead of implicit. Only do this where the
  closure (the set of files sharing the block) is narrow -- a handful
  of files. A wide closure of this kind would have the same blast-radius
  problem as case 1 without the "it's just a declaration swap" safety
  margin; treat it as a module instead, or handle case-by-case.

**Closure size is a hint, not the test -- always trace the actual data
flow before picking a treatment.** `common/intpar/` looked like a
textbook case 2 candidate at first glance: shared by exactly two
files (`numerics/bsstep.f90`, `atm/atm_lib.f90`), same narrow-closure
shape as `common/tridi/`. But tracing what each file actually reads
showed they use *disjoint* members -- `atm_lib.f90` only reads
`tolerance_fraction`, `bsstep.f90` only reads `max_stage_index`/
`extrap_order` -- and all three are NAMELIST /physics/ values set once
at startup. Nothing computes a value in one file for the other to
consume; they just happen to share a block that a third file
(`core/parmin.f90`) initializes. That's case 1, not case 2. The tell:
in a genuine case 2, the *same* member is written by one file and read
by another *within one call sequence* (`common/tridi/`'s sub_diag
etc., filled by `ccoeft`/`dcoeft`, read moments later by `ctridi`/
`tridia`). If each file's actual usage doesn't overlap, or nothing
about the value changes call to call, it's case 1 regardless of how
few files share the block.

## Deciding whether a file's COMMON usage is even real

Before doing either conversion, check whether the block is actually
read or written in the file's body, or just declared and never
touched. This codebase has many files that declare a COMMON block
purely as boilerplate inherited from the original F77 source (a
"declared only to preserve layout" placeholder, sometimes documented
as such in a comment from the earlier readability-refactor pass, but
not always -- verify by grep, don't trust the comment alone). If
nothing in the body references any member, the block can simply be
deleted from that file's declarations -- nothing to convert, no module
or argument needed for that file.

## Folder names describe function, not coupling -- verify before trusting

Two folders so far (`util/`, `numerics/`) turned out to contain files
whose *behavior* matched the folder name but whose *coupling* did not:
`util/entime.f90`/`htimer.f90` etc. are a real timestep-control
subsystem sharing genuine model state, not trivial utilities; and
`numerics/findsh.f90`, `lax_wendrof1.f90`/`lax_wendrof2.f90` are
physics routines (shell-boundary finding, element-diffusion PDE
solving) that happen to be implemented with numerical techniques, not
generic numerics.

Before converting any not-yet-touched file, check its actual COMMON
usage (`grep "^\s*common" <file>`, then check who else shares each
block, then check whether the file's own body actually reads/writes
each member) rather than assuming from its folder or its name. If a
file's real coupling shows it belongs to a different physics domain
than its current folder suggests, relocate it there rather than
force-converting it in place as if it were a generic utility.

## The categorization recipe (apply to any not-yet-converted set of files)

For each file:
1. Does it use COMMON at all? If not: it's a zero-risk conversion
   target (module, body-included, same as the `numerics_lib` core 20).
2. If it does, is each declared block actually used in the body, or
   dead boilerplate? Dead blocks: just delete the declaration.
3. For blocks that are genuinely used: is this global constant/config
   (case 1 above) or real per-call data flow (case 2 above)?
4. For case 2: is the closure (files sharing the block) narrow (a
   handful -> argument-convert now) or wide (treat as case 1, or defer
   to a case-by-case decision)?
5. Does the file's actual behavior match the folder it's currently
   filed under? If not, relocate rather than convert in place.

## Naming

- `_lib` is reserved for exactly one module per domain: the small,
  public, callable-from-anywhere surface for that domain (`numerics_lib`,
  `const_lib`, `luout_lib`, and later `eos_lib`, `kap_lib`, etc.) --
  matching MESA's own convention (`eos_lib`, `kap_lib`, `num_lib`).
- Internal implementation detail (where a domain is large/entangled
  enough to need a facade over existing legacy code rather than a
  direct merge) does not need a uniform suffix -- name it for what it
  does, the way MESA's private implementation modules are (`eosDT_eval`,
  not `eosDT_mod`). Only the one public facade per domain gets `_lib`.
- For a domain small and self-contained enough to convert directly
  (no facade needed, as with `numerics_lib`), merge routines into the
  one `_lib` file rather than keeping one module per routine -- Fortran
  can't split a single `module` across files (only submodules can, a
  heavier mechanism nothing else here uses), and one module per routine
  would mean either breaking the one-routine-per-file convention
  everywhere else in this codebase's history, or accepting an
  aggregator-of-aggregators. Merge when the domain is a real "library"
  of leaf routines (as `numerics_lib` is); keep the facade-over-legacy-
  code pattern for domains still entangled with the solver (`eos`,
  `kap`, etc. -- see below).

## Physics domains still entangled with the solver (eos, kap, nuclear,
## atm, wind, mixing, rotation)

As of 2026-08-21 every domain, including these, is COMMON-free -- all
COMMON blocks repo-wide are now module/derived-type state (see
`state/*.f90`, `const_lib.f90`). That mechanical conversion (phase
one) is done. What these domains still have is architectural
entanglement: external files reach directly into a domain's internal
module state (`use opal_eos_lib` and touch `opal_eos%...` fields
directly, etc.) instead of going through a small, explicit-interface
entry point, and that reach-through is woven directly into the Henyey
coefficient-building code (`misc/coefft.f90`) or otherwise central to
the solver. This is phase two, deliberately sequenced after phase one
finished. Sequencing it this way isn't just cleanup ordering -- designing each
facade's return arguments (e.g. `eos_get`'s `p, gamma1, chit, chirho,
cp, ...`) is much easier once the domain's per-model state already
lives in an explicit derived type (oldmod_lib/scrtch_lib/
turnover_lib/light_burn_lib/engeb_diag_lib-style) instead of being
scattered across raw COMMON; doing the facade design first would mean
redoing it once the state conversion catches up anyway. When that
phase starts:

1. Design one small, explicit-interface entry point per domain (e.g.
   `eos_get(rho, t, ... -> p, gamma1, chit, chirho, cp, ...)`).
2. Make that entry point's body a thin facade that calls the existing,
   unchanged legacy COMMON-based code and copies results into the new
   return arguments. This is a pure refactor -- verify it produces
   byte-identical output by construction, since no numerics changed.
3. Migrate call sites to the new facade one at a time, each its own
   verified commit.
4. Only after every external caller has migrated, and no file outside
   the domain still reaches into its internal COMMON blocks, consider
   (as an optional, separate, later step) removing COMMON from the
   domain's own internals. The swap-value goal is already achieved
   after step 3; step 4 is not urgent.

**First facade, `eos_get` (2026-08-21)**: `eos/eos_lib.f90`, `module
eos_lib` / `contains subroutine eos_get(...)`, is the first domain
facade built under this plan -- a second precedent (alongside
`numerics_lib`) for "a module hosting real callable subroutines, not
just derived-type state." It has the same 27-arg signature as
`eqstat.f90` plus one new trailing `optional` argument,
`composition_at_zone(15)`; Fortran requires the module-procedure form
(not a bare external subroutine) for the optional argument to be
legal. Internally it dispatches to `meqos` (MHD path) or `eqstat`
(non-MHD path) based on `use_mhd_eos`, both called unchanged -- and,
on the non-MHD path, centralizes the Debye-Hückel composition setup
(`debye_huckel_x/y/z_total/z(3)`) that 6 of its 10 callers used to
duplicate verbatim, gated on `present(composition_at_zone) .and.
use_debye_huckel_correction`. Migrated 10 files / 13 call sites
(`atm/atm_lib.f90`, `atm/qatm.f90`, `atm/qenv.f90`, `misc/coefft.f90`,
`io/wrtout.f90`, `misc/physic.f90`, `core/starin.f90`,
`mixing/hsubp.f90`, `mixing/sconvec.f90` x3, `wind/massloss.f90` x2)
to call `eos_get` instead of duplicating the `use_mhd_eos` if/else at
each site. `atm/turnover/calcad.f90` was deliberately left alone -- confirmed
against the original F77 source, it never participated in the
`eqstat`/`meqos` dispatch (calls `esac06`/`eqstat2` directly under its
own `use_opal2006_eos` check) and isn't part of this pattern. The
migration also fixed a real bug in `core/starin.f90` (a missing
`ELSE`, confirmed via git archaeology against the pre-modernization
source, that made it call both `meqos` and `eqstat` when MHD was on)
and extended real MHD coverage to `hsubp.f90`/`sconvec.f90`/
`massloss.f90`, none of which ever checked `use_mhd_eos` even in the
original decades-old source -- both changes were deliberate,
user-approved parts of the migration, not accidents. Caveat: no
Stage-0 regression case sets `LMHD`, so the byte-identical verification
covers the non-MHD path only; the `starin.f90` fix and the 3
newly-MHD-capable files have no reference output to check against.

**`eos/` reorganized by EOS type (2026-08-21)**: per the YREC public
release paper, the code offers three user-facing EOS choices --
**OPAL** (recommended, tabulated, three table vintages), **SCV**
(Saumon-Chabrier-Van Horn, tabulated, covers the cool/dense phase
space OPAL's tables don't), and **Yale** (the Prather 1976 Saha solve
plus a fully-ionized/relativistic regime, analytic, always available)
-- plus a fourth, undocumented-in-the-paper legacy option, **MHD**
(Mihalas-Hummer-Dappen, tabulated), which the tie-breaker in
`eqstat2` never touches (`meqos`/`eos_get` route around it entirely
via `use_mhd_eos`). `eos/` now has one subfolder per type:
`eos/opal/` (27 files: `eqbound*`, `esac*`, `gmass*`, `oeqos*`,
`quad*`, `radsub*`, `readco*`, `rhoofp*`, `t6rinterp*`/`t6rinteos*`),
`eos/mhd/` (`mhdpx`, `mhdpx1`, `mhdpx2`, `mhdst`, `mhdst1`, `mhdtbl`,
and `meqos` -- the MHD top-level entry point, grouped here the same
way `oeqos*` sit inside `opal/` as OPAL's top-level entries),
`eos/scv/` (`eqscve`, `eqscvg`), and `eos/yale/` (`eqsaha`, `eqrelv`).
`eos/eos_lib.f90` (the cross-type facade), `eos/eqstat.f90` (see
below), and `eos/mu.f90` (a mean-molecular-weight helper genuinely
shared by both MHD and all three OPAL vintages, so it doesn't belong
to any one subfolder) stay at `eos/` root. The corresponding
derived-type module state, `state/opal_eos_lib.f90` and
`state/mhd_eos_lib.f90`, deliberately stayed in `state/` rather than
moving alongside their routines -- keeps the established convention
that all per-domain module state lives together in `state/*_lib.f90`,
by user's explicit choice. Before moving anything, checked for
cross-family calls (none: MHD/OPAL/SCV/Yale don't call into each
other, only into `eos/mu.f90` or shared `numerics_lib.f90`/`util/`
helpers) and fixed every stale `eos/<file>.f90` path reference in
comments across the repo (`state/atm_table_lib.f90`,
`state/mhd_eos_lib.f90`, `state/opal_eos_lib.f90`,
`const/const_lib.f90`, `numerics/numerics_lib.f90`, `eos/eos_lib.f90`
itself) to point at the new subfolder paths.

**`eqstat`/`eqstat2` co-located in one file, same reorg**: `eqstat` is
a pure numerical-differentiation wrapper around `eqstat2` and the two
have always been a matched pair, so `eqstat2.f90`'s subroutine moved
into `eqstat.f90` (both remain plain external subroutines, not module
procedures -- callers, `eos_lib.f90`'s `eos_get` for `eqstat` and
`atm/turnover/calcad.f90` calling `eqstat2` directly, needed no changes).
Co-locating them let gfortran see both signatures in the same
compilation unit for the first time, which surfaced a genuine (though
harmless) pre-existing interface violation: `eqstat` was passing its
own `intent(in) metal_fraction` into `eqstat2`'s `intent(inout)`
dummy of the same name -- illegal by the Fortran standard, previously
invisible since the two were separate translation units with no
shared interface. Traced rather than silently "fixed": the only write
to `metal_fraction` inside `eqstat2`'s body was a save-then-restore
pair around the whole routine, and none of the callees it's passed to
(`opal/oeqos*.f90`, `scv/eqscve.f90`) mutate their own `metal_fraction`
dummy (all declare it `intent(in)`) -- so the restore was provably a
no-op. Removed the dead save/restore, changed `eqstat2`'s dummy to
`intent(in)`; zero behavior change, verified byte-identical.

**Second facade, `kap_get` (2026-08-21)**: unlike `eos`, investigation
found `kap/getopac.f90` was *already* a single, clean, explicit-
interface entry point -- all 8 external call sites did a plain
`call getopac(...)`, no duplicated dispatch logic anywhere (`eos`'s
core problem). So this wasn't a new-facade build the way `eos_get`
was: `getopac.f90` was renamed to `kap/kap_lib.f90` (`module kap_lib`
/ `contains subroutine kap_get(...)`), body and argument list
unchanged, purely to give `kap` the same public-facade shape as
`eos_lib` (matching this doc's own naming rule, which already
anticipated `kap_lib` by name) -- a third precedent (alongside
`numerics_lib`, `eos_lib`) for "a module hosting a real callable
subroutine." The 8 callers got `use kap_lib` added and their calls
renamed to `kap_get`; no Makefile change beyond adding
`kap/kap_lib.f90` to `MODULE_SRCS`.

Investigation also found two misplaced-domain file groups inside
`kap/`, same pattern as `meval.f90` (single external caller elsewhere,
zero domain content) -- moved out:
- `alsurfp.f90`/`alfilein.f90`/`altabinit.f90` -> `atm/`. These
  interpolate Allard NextGen atmosphere tables for Log(P)/Log(T) at
  fixed Teff/GL/tau=100 -- a boundary condition, not opacity -- and
  sit right alongside `atm/tables/kcsurfp.f90` (the Kurucz/Castelli
  equivalent, already correctly filed) at their only call sites
  (`atm/atm_lib.f90`, `wind/massloss.f90`).
- `ifermi12.f90`/`zfermim12.f90` -> `nuclear/`. Fermi-Dirac degeneracy
  integrals called only from `nuclear/nuclear_lib.f90`, nowhere in
  `kap/` itself.

The remaining 34 files split cleanly into 7 physical groups, matching
`kap_get`'s own flags exactly: `kap/alex06/` (Ferguson et al. 2005
low-T opacities, 3 files), `kap/alex94/` (Alexander 1994 low-T
opacities -- note the code's own flag is `use_alex95_tables` despite
the tables being 1994-dated; not renamed, just flagged), `kap/kurucz90/`
(4 files), `kap/opal95/` (Iglesias & Rogers 1996 H/He/metal mixtures,
7 files), `kap/opal92/` (the Livermore "LL"-named tables, 7 files),
`kap/laol89/` (Los Alamos, including its pure-Z variant used for the
He-burning correction, 7 files), `kap/conductive/` (Potekhin, per
Cassisi et al. 2007, 2 files). `kap_lib.f90` and `setupopac.f90` (the
startup-time table loader, called once from `setup/setups.f90`) stay
at `kap/` root. Per the YREC public release paper, this maps onto
three *user-facing* opacity families (conductive, atomic, molecular);
the paper's third atomic family (OP/Badnell 2005) doesn't have an
obvious match among `opal95`/`opal92`/`laol89` -- left unresolved,
not blocking, since the reorg reflects the code's actual structure as
it exists today rather than the paper's summary framing.

Verification: full clean build + Stage-0 byte-identical regression,
each step (facade rename, misplaced-file moves, subfolder split)
checked independently before combining into the commit.

**Deferred: unneeded blanket `SAVE` statements.** Most F77-era
routines in this codebase open with a bare `save` (saves every local
in that scoping unit), inherited caution from the original source
rather than a real cross-call dependency in most cases -- distinct
from genuine caches like `eos/opal/quad.f90`'s `cache_slot`-indexed
coefficient arrays or `kap/opal92/yllo3d.f90`'s `abund_index`/
`temp_index`/`dens_index`, which really do need a value written on
one call to survive and be read on the next. `kap/kap_lib.f90`'s
`kap_get` was traced by hand (2026-08-21) and found to have no such
dependency -- every local is written before it's ever read, on every
path through the routine, within a single call -- so its `save` looks
safe to drop. Decided to hold off on doing this (or auditing any other
file for the same thing) until the rest of the phase-two domain sweep
(`nuclear`, `atm`, `wind`, `mixing`, `rotation`) is finished, so this
kind of cleanup doesn't get interleaved with the facade/reorg work.
When it's picked back up: trace each candidate file's locals by hand
(read-before-write on every path, including early returns) before
removing its `save`, same rigor as the `metal_fraction` intent fix
above, and verify with the standard full-build + Stage-0 byte-diff
per file or small batch.

**Third facade, `atm_get` (2026-08-21)**: per the YREC public release
paper's section 2.3.6 ("Boundary Conditions"), `atm/` computes the
(P, T, R) envelope solution at the model's fitting point for a given
(Teff, L). Investigation found the domain has two callable entry
points at genuinely different grains, not one -- different from both
`eos` and `kap`: `envint.f90` (renamed `atm_lib.f90`'s `atm_get`) is
the generic "solve one envelope for this exact (Teff, L)" primitive,
already called uniformly (no duplicated dispatch) from 7 sites; and
`surfbc.f90` is a specialized, single-caller (`core/crrect.f90`, the
solver's Newton-Raphson corrector) wrapper that adds a cached
(Teff, L) triangle for cheap derivative interpolation plus the paper's
"hot edge -> fall back to gray atmosphere" logic. Only `envint`/
`atm_get` was renamed to match `eos_lib`/`kap_lib`'s shape (a module
hosting the domain's generic primitive); `surfbc.f90` stays as its own
file, already correctly using `atm_get` as a plain, non-duplicated
call.

Renaming `envint` surfaced a second instance of the `metal_fraction`-
shaped bug: `atm_get`'s own header already documented
`hydrogen_fraction`/`metal_fraction`/`pressure_rotation_factor`/
`temperature_rotation_factor` as legitimately `intent(inout)` (they
flow through to `bsstep`/`mmid`'s callback machinery), but
`surfbc.f90`, which relays them straight through, had declared its
own same-named dummies `intent(in)` -- too narrow, invisible while
`envint` was a bare external subroutine with no interface to check
against, caught the moment `atm_lib.f90` gave it one. Fixed by
widening `surfbc.f90`'s dummies to `intent(inout)`, matching
`core/crrect.f90` (its only caller), which already passes real local
variables there -- zero behavior change.

Investigation also found one misplaced file (the reverse move of the
`kap` sweep's `alsurfp`/`alfilein`/`altabinit` discovery):
`surfopac.f90` had zero atm-domain content -- it only refreshes
cached table slices in `kap/opal95/`, `kap/opal92/`, `kap/alex94/` --
and its callers (`core/starin.f90`, `setup/hpoint.f90`) aren't in
`atm/` either. Moved to `kap/`, alongside `setupopac.f90`.

The remaining 9 files split into 2 groups matching the paper's own
structure: `atm/tables/` (the 3 table-lookup boundary-condition
sources the paper names -- Kurucz 1993 via `surfp.f90`, Castelli &
Kurucz 2003 via `kcsurfp.f90`, Allard & Hauschildt 1995 via
`alsurfp.f90`/`alfilein.f90`/`altabinit.f90`) and `atm/turnover/`
(`gettau.f90`/`taucal.f90`/`tauint.f90`/`tauintnew.f90` -- the
convective-overturn-timescale calculation, a distinct sub-concern
built on `atm_get`'s machinery but not itself boundary-condition
physics; consumed cross-domain by `core/`, `io/`, and
`rotation/getw.f90`, not misplaced, just multi-consumer). `atm_lib.f90`,
`qatm.f90`, `qenv.f90`, and `surfbc.f90` -- the core boundary-condition
machinery -- stay at `atm/` root.

Verification: full clean build + Stage-0 byte-identical regression,
checked after the facade rename, after the `surfbc.f90` intent fix,
and after the subfolder split, before combining into one commit.

**`nuclear` completed its already-started library merge (2026-08-21)**:
unlike `eos`/`kap`/`atm`, `nuclear_lib.f90` was already a `numerics_lib`-
style module aggregating 6 routines (`neutrino`, `nulosses`, `azbar`,
`sneut`, `rates`, `eqburn`) -- this domain was partway through the
"small library, merge directly" treatment (`GUIDELINES.md`'s
alternative to the facade pattern), not the facade pattern itself.
Investigation found no `eos`-style duplicated dispatch (the domain's
only 2-choice call site, `rotation/getw.f90`'s `liburn` vs `liburn2`,
is single-sited, not repeated) -- the actual gap was 7 remaining files
(`dburn.f90`, `dburnm.f90`, `deutrate.f90`, `engeb.f90`, `liburn.f90`
-- which itself bundled a second small subroutine, `safedivexp` --
`liburn2.f90`, `lirate88.f90`) plus 2 private Fermi-Dirac integral
functions (`ifermi12.f90`/`zfermim12.f90`, only called from this
module's own `rates`) still sitting as separate plain-external files.
`dburn.f90`/`dburnm.f90` were investigated for a possible
`eos_get`-style optional-arg collapse and explicitly kept separate:
they differ in rate-data source, timestep units (seconds vs. Gyr),
convergence threshold (1e-11 vs. 1e-14, plus an accretion-aware
pre-check `dburn` has and `dburnm` doesn't), and the accretion-
weighting formula structure -- genuine algorithmic differences for two
different call-site contexts, not textually-identical duplication.
All 9 files moved into `nuclear_lib.f90` unchanged; the 6 caller files
that didn't already `use nuclear_lib` (`setup/midmod.f90`,
`misc/coefft.f90`, `core/main.f90`, `util/ytime.f90`,
`rotation/getw.f90`, `mixing/bursmix.f90`) got it added.

This merge surfaced three mechanical issues worth remembering:

1. `engeb.f90`'s own `use nuclear_lib` (it called `neutrino`) had to be
   deleted -- a module cannot `use` itself; host association already
   gives every contained procedure access to every sibling procedure.
2. `sneut`'s existing declaration `double precision ifermi12,
   zfermim12, ...` (an external-function type declaration, needed
   when those two were separate files) had to be stripped down to just
   its genuine locals once `ifermi12`/`zfermim12` moved into the same
   module -- otherwise the plain local-type declaration shadows the
   host-associated module function, producing an `_ifermi12_` linker
   symbol the build can no longer provide (the function is now
   `___nuclear_lib_MOD_ifermi12`). Watch for this shape whenever a
   function that used to be called externally moves into the same
   module as one of its existing callers.
3. **A third instance of the `eqstat`/`eqstat2`-shaped intent bug, but
   mirror-image this time.** `rotation/getw.f90` calls `liburn`/
   `liburn2` passing its own `shell_mass`, which it had declared
   `intent(in)`. `liburn`/`liburn2` need `intent(inout)` for
   `shell_mass` -- unlike `eqstat2`'s `metal_fraction` (a dead
   save/restore that could just be deleted), here the CZ-base mass is
   genuinely perturbed and the perturbed value is read by the
   following rate/abundance sums before being restored, so the write
   is real and the callee's `intent(inout)` was already correct; the
   fix was to widen the *caller*'s declaration to match, not narrow
   the callee. Net effect on `getw.f90`'s own array is still zero
   (every path that writes also restores before return) -- that's
   exactly why its `intent(in)` had gone unnoticed for so long. Lesson
   generalized: when this class of error appears, don't assume the
   narrow-the-callee fix from the first two instances applies --
   trace whether the callee's write is truly dead (delete + narrow) or
   genuinely used before being restored (keep + widen the caller).
   3 for 3 so far on "co-locating/adding an explicit interface
   surfaces a latent intent bug" -- treat it as expected, not
   surprising, whenever this kind of rename happens.

The Makefile's `MODULE_SRCS` list also needed a genuine fix, not just
reordering: a comment already on the existing `nuclear/nuclear_lib.o`
order-only-prerequisite line states that list *position* in
`MODULE_SRCS` does not guarantee build order between module objects
under parallel (`-j`) builds -- only an explicit order-only
prerequisite does. `atm_lib.f90`'s earlier `MODULE_SRCS` reorder (this
same day) had worked by luck, not by guarantee. Replaced both ad hoc
reorders with explicit `nuclear/nuclear_lib.o: | ...` and
`atm/atm_lib.o: | ...` prerequisite lines naming every module object
each actually `use`s, matching the established precedent's own stated
correct pattern.

Verification: full clean build (through 3 rounds of intent/linker
fixes) + Stage-0 byte-identical regression.

**`wind` needed no facade, only relocation (2026-08-21)**: investigation
found no `eos`-style duplicated dispatch -- `mwind.f90` and
`mcowind.f90` each already internally dispatch Kawaler-type
(`wind.f90`/`cowind.f90`) vs. Matt-type angular-momentum-loss physics,
and where `rotation/seculr.f90` calls both `mwind` and `mcowind`, it's
for genuinely different rotation-model branches (differential vs.
solid-body, different `fix_omega_at_surface` cases), not repeated
duplication of the same choice. `amcalc.f90`, `cowind.f90`/
`mcowind.f90`, `mdot.f90`/`massloss.f90`, `mwind.f90`/`wind.f90` are
all genuine wind-domain physics (magnetic-wind angular-momentum loss,
mass loss/accretion), correctly filed even though their consumers are
in `rotation/`/`core/` -- physics feeding a solver isn't misplacement.

Three files were, though, found and relocated by the same
zero-domain-content/callers-elsewhere test used for `meval.f90`/
`alsurfp.f90`/`surfopac.f90`:
- `wczimp.f90` -> `rotation/`. Enforces a user-specified rotation
  profile (solid body / constant specific angular momentum / power
  law) in a convective zone -- rotation-state control, not wind-driven
  angular-momentum loss. Called only from `rotation/getw.f90`/
  `rotation/getrot.f90`.
- `viscos.f90` -> `rotation/`. Computes microscopic viscosity/thermal
  diffusivity "for use by the rotational-mixing/instability diffusion
  routines" (its own header). Zero wind content. Called only from
  `misc/physic.f90`.
- `calcad.f90` -> `atm/turnover/`. The acoustic-depth diagnostic
  (already familiar from the `eos`/`kap` sweeps as the deliberately-
  unmigrated `eqstat2`/`esac06` caller) integrates sound-speed
  profiles for asteroseismology output -- no wind content, and the
  same category of work as `atm/turnover/`'s other structure-
  integration diagnostics (`gettau`/`taucal`/`tauint`/`tauintnew`).
  Called only from `io/wrtout.f90`.

`wcz.f90` -- confirmed zero callers anywhere, its own header already
says every caller switched to `wczimp.f90` -- was deliberately **left
in place** rather than deleted, per explicit user choice: it's
documented dead code, not an accident, and the project isn't treating
"unreachable" as automatic grounds for removal here.

Verification: full clean build + Stage-0 byte-identical regression
after the relocations.

## Build mechanics

- Any file introducing `module ... contains` must be added to the
  Makefile's `MODULE_SRCS` list so it's compiled before anything that
  `use`s it (GNU Make has no built-in Fortran module dependency
  scanning; this is an order-only prerequisite, safe under `-j`
  parallel builds).
- `.mod` files are build artifacts (like `.o`), gitignored, never
  committed.
- Once a subroutine gains an explicit interface (by moving into a
  module), gfortran starts checking argument count/type/rank between
  caller and callee for the first time -- calls that were silently
  accepted under the old implicit-interface F77-style convention can
  newly fail to compile. This is virtually always a genuine
  pre-existing bug being surfaced, not something the conversion
  introduced. Investigate the actual call site and the real established
  calling convention elsewhere in the codebase (grep for how the same
  routine is called correctly elsewhere) rather than guessing a fix or
  changing the routine's own dummy-argument contract to paper over a
  mismatch.

## Verification discipline (non-negotiable for every conversion)

1. Full clean build (`make clean && make -j4`) succeeds with no new
   errors. Pre-existing warnings unrelated to the change are fine.
2. Byte-identical diff, **excluding only the git-hash version-stamp
   line** (`# YREC v...`), against the trusted Stage 0 reference
   (`examples/run_standard_solar_model/standard/`, all 4 cases: `.short`,
   `.track`, `.store`). This is a strict byte-diff, not the looser
   fractional-tolerance check `test_all.py`'s own default uses --
   floating-point-precision regressions (e.g. an "improved" literal
   losing precision) have historically been invisible under a loose
   tolerance and only caught this way.
3. When a change touches a code path Stage 0 doesn't exercise (no
   example namelist enables the relevant flag), say so explicitly
   rather than silently claiming full coverage, and prefer running at
   least one additional real test case that does exercise it if one is
   available in the repo (e.g. `examples/run_from_dbl_to_zams`'s
   low-mass cases for anything ZAMS-adjacent). If no such case exists,
   rely on direct code-inspection verification (checking the change is
   behavior-preserving by construction) and disclose that limitation.
4. Only commit after verification passes. One conversion (or one
   tightly related batch, like a domain's full module conversion) per
   commit, with a message that states what changed, why, and what was
   verified -- never batch multiple unverified changes into one commit.
5. Never run the full multi-directory regression suite (`test_all.py`
   across every `examples/` directory) as a matter of course -- it has
   at least one known runaway case (`run_from_dbl_to_zams`'s 14 Msun
   ZAMS model never converges and will fill disk if left running).
   Stage 0 plus a targeted extra case is the right default scope;
   widen deliberately, not by habit.

## Regression-testing gotchas learned the hard way

- `test_all.py`'s own comparison auto-creates a reference from a run's
  own output if none exists yet (`if not os.path.isfile(ref):
  shutil.copyfile(...)`). A directory with no pre-existing `standard/`
  reference gives you a "does it run" smoke test, not a real regression
  comparison, even though it reports as a pass. Check the reference's
  timestamp/git history before trusting a green result from a directory
  you haven't deliberately baselined.
- Precision-losing "corrections" are the single most common way an
  otherwise-mechanical change silently breaks output: e.g. an
  unsuffixed literal like `83.14511` rewritten to `83.14511d0` is
  numerically different (single-precision-then-widened vs. direct
  double), and can cascade into thousands of differing lines through
  an iterative solver. Never "improve" a literal's precision suffix
  while doing an unrelated mechanical conversion.

## Re-triage after every conversion, not just at the start

A COMMON-to-module conversion can leave *other* files newly eligible
for a treatment they weren't eligible for before, as a side effect --
and nothing flags this automatically. Concretely: converting
`common/const1/`/`common/luout/`/`common/intpar/`/`common/tridi/` left
9 files in `numerics/` (`ctridi`, `tridia`, `bsstep`, `intpol`,
`splint`, `splintd2`, `trapzd`, `qgauss`, `intpt`) fully COMMON-free,
but each of those commits was correctly scoped to just its own COMMON
block, so none of them checked whether the change also unlocked a
merge-into-`numerics_lib`-style conversion elsewhere. They sat
unconverted for several commits before anyone (the user, in this
case) asked why. After converting a COMMON block, re-run the
COMMON-usage census (`grep "^\s*common" <file>`) on every file that
*used to* share it, not just the files directly touched by that one
conversion -- a file can cross from "still coupled" to "fully clean"
without being the direct target of the change that did it.

## Known interface-violation shapes

Every module conversion so far has surfaced at least one pre-existing
bug via the new explicit interface (per the build-mechanics note
above). Recognizing which of these three shapes you're looking at
saves re-deriving the fix from scratch:

- **Array too small.** Caller declares a local array shorter than the
  callee's dummy (e.g. `atm/tables/alsurfp.f90` passed 4-element arrays where
  `polint`'s dummy is `xa(20)`/`ya(20)`). Tell: "Actual argument
  contains too few elements." Fix: widen the caller's declared array
  size to match the callee's real contract (check other correct
  callers, e.g. `rotation/fpft.f90`, to confirm the right size) --
  never shrink the callee's dummy to match the caller.
- **Scalar vs. length-1 array.** Caller passes a `(1)`-dimensioned
  array where the callee's dummy is a plain scalar, or vice versa
  (`atm/turnover/calcad.f90` did this twice, against `numerics/boole.f90` and
  against `splint`). Tell: "Rank mismatch ... (scalar and rank-1)."
  Both sides are legal, standalone F77 idioms (a length-1 array and a
  scalar share identical memory layout, so this always worked via
  sequence association under the old implicit interface) -- fix by
  passing the array *element* (`x(1)`) instead of the whole array at
  the call site, not by changing either declaration.
- **Intent mismatch.** Caller declares its own dummy argument
  `intent(in)`, then passes it as the actual argument to a callee
  whose corresponding dummy is `intent(inout)`/`intent(out)`
  (`atm/atm_lib.f90` did this passing values through `bsstep`/`mmid` to
  an arbitrary `deriv` callback that might modify them). Tell: "Dummy
  argument ... with INTENT(IN) in variable definition context." Fix:
  first determine which side is actually right by checking whether the
  callee (or something *it* calls) ever assigns to the value -- if so,
  the callee's wider intent reflects real behavior and the caller's
  declaration should be widened to match (safe whenever every call
  site already passes a real variable, never a literal/expression, in
  that position -- check this before widening). Don't narrow the
  callee's intent to make the error go away without checking first.
- **Dummy argument collides with a module variable.** A file's own
  subroutine signature declares a dummy argument whose name happens to
  match a `const_lib` (or other module) variable added by an *earlier,
  unrelated* conversion -- surfaces only when this file's own COMMON
  block finally gets converted and needs its first `use const_lib`.
  Real example: `io/getyrec7.f90`/`getmodel2.f90`/`putyrec7.f90`/
  `putmodel2.f90` read/write old model-file fields named
  `rotation_active`, `envelope_overshoot_active`, `lovstc`, etc as
  `intent(in)`/`intent(out)` dummy arguments -- same names as unrelated
  runtime-config variables already in `const_lib`. Tell: "ambiguous
  reference to ... from current program unit" or "has no IMPLICIT
  type" errors pointing at the dummy-argument declaration line, right
  after adding a `use` statement. Fix: `use const_lib, only: <the one
  or two members actually needed>` instead of the blanket `use
  const_lib` -- never rename the dummy arguments themselves when
  they're part of a shared calling convention across sibling files
  (check for siblings with the same signature pattern before deciding
  a fix is file-local).

## Two correctness traps specific to the "keep local + copy-assign" pattern

The `stolr0`/`tscut`/`lovste`-style pattern (local NAMELIST-spelled
variable, copy-assigned into a differently-spelled const_lib canonical
name right after the namelist read) is safe in the common case, but
two real bugs have come from it -- both silent (wrong numbers, not
build failures), both only surfacing in test cases that actually
exercise the affected code path:

- **The namelist member isn't really the source of truth.** Some
  NAMELIST-spelled locals are themselves overwritten later by other
  logic before their "real" value is used elsewhere -- most commonly
  `setup/remap.f90` recomputing a value from other, more-user-friendly
  namelist inputs (`common/ctol/`'s `hpttol` -> `chi_grid_scale`: once
  `lnewvars` is set, remap.f90 derives `chi_grid_scale` from
  `tol_dm_min`/`tol_dm_max`/etc, not from `hpttol` at all). Copying the
  raw namelist-read local into const_lib right after the namelist read
  captures the *wrong* value if anything downstream (especially
  `call remap`) recomputes it under the canonical name directly.
  Before wiring the copy-assignment, check whether the local name is
  read anywhere *else* in core/parmin.f90 after the namelist read (not
  just written) -- if the only remaining readers are core/parmin.f90's
  own diagnostic writes, point those writes at the canonical
  const_lib name instead of the stale local, rather than copying.
- **Local values can be overridden after the first copy-assignment.**
  core/parmin.f90 has several post-namelist-read consistency-
  enforcement blocks (e.g. "DBG 12/95 ENSURE CORRECT PARAMETERS FOR Z
  DIFFUSION": `if (ldifz) ldify=.true.`) that mutate a NAMELIST-spelled
  local *after* it may have already been copy-assigned into its
  const_lib canonical name elsewhere in the same subroutine. The
  const_lib copy silently keeps the pre-override value unless a second
  copy-assignment runs after the override block too. After adding any
  copy-assignment, grep the rest of core/parmin.f90 for further plain
  `name = ...` assignments to that same local (not just its
  declaration and namelist appearance) and re-copy-assign after each
  one found. This is easy to miss in Stage 0 testing specifically
  because Stage 0's own namelist value for the affected flag may
  already match the override's forced value (masking the bug) --
  another reason the extra `run_from_dbl_to_zams` case, or any test
  case exercising the *opposite* namelist setting, earns its place in
  the verification set.

## Git hygiene when interleaving sub-tasks in one session

`git mv`/`git add` stage changes that persist across unrelated later
work. Doing a `git mv` for one sub-task (e.g. relocating misfiled
files) and then, several steps later, running `git commit -m "..."`
without an explicit pathspec for an unrelated sub-task will silently
sweep the still-staged renames into that commit -- the result isn't
wrong (everything staged did get verified together), but the commit
message won't describe what's actually in it. Re-run `git status`
immediately before every commit, not just after making the edits you
were thinking about, and scope `git commit` with an explicit file
pathspec whenever more than one sub-task's changes might be staged at
once.

## Scripting a COMMON-to-module conversion: three parser traps

All three of these produced silent, wrong results (not crashes) until
caught by manually inspecting a sample conversion before applying
broadly -- always do that inspection, don't trust a script that
reports 100% success without having read at least one real diff.

- **macOS/BSD `sed -E` silently ignores `\s`.** It's a GNU-sed/Perl
  extension, not POSIX; on macOS a pattern like
  `s/^\s*(double precision|...) *:: *//` simply fails to match (no
  error) and the input passes through unchanged. Downstream steps in
  the same pipe (e.g. a later `tr -d ' '`) then silently glue the
  unstripped prefix onto the next token (`logical ::
  use_extended_composition` became `logical::use_extended_composition`
  as one string), corrupting a symbol-list extraction without any
  error message. Use `[[:space:]]`/`[ \t]` in shell `sed`, or do this
  kind of parsing in Python instead (as most of this project's
  conversion scripts already do) rather than a quick ad hoc `sed`
  one-liner.

- **A naive `str.split(",")` breaks on 2D array bounds.** A
  declaration like `old_composition(15,json)` has a comma *inside*
  its dimension spec; splitting the whole declaration list on every
  comma turns this into two bogus tokens (`old_composition(15` and
  `json)`), silently corrupting the parsed member list. Use a
  paren-depth-aware splitter (only split on commas at depth 0) instead
  of a plain `.split(",")` whenever parsing a Fortran declaration or
  argument list that might contain arrays.
- **Declaration order does not have to match COMMON order.** Fortran
  groups type declarations by type (all the `double precision`
  members together, then all the `logical` members, then `integer`,
  etc.), which is often a different order than the COMMON statement's
  own positional listing -- this is legal and common, not a bug in the
  source. Comparing the declared-name list against the COMMON-member
  list must be done as an unordered set comparison, not a sequence
  comparison, or a correct file will be wrongly rejected as
  "mismatched."

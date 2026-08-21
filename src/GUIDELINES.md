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
files (`numerics/bsstep.f90`, `atm/envint.f90`), same narrow-closure
shape as `common/tridi/`. But tracing what each file actually reads
showed they use *disjoint* members -- `envint.f90` only reads
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

These are not COMMON-free and their COMMON is woven directly into the
Henyey coefficient-building code (`misc/coefft.f90`) or otherwise
central to the solver. Do not attempt full COMMON elimination on the
first pass for these. This is deliberately phase two: finish
converting every domain's COMMON blocks to modules/derived types
first (the mechanical work this whole document otherwise describes),
then come back and do this disentangling as a separate later pass.
Sequencing it this way isn't just cleanup ordering -- designing each
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
  callee's dummy (e.g. `kap/alsurfp.f90` passed 4-element arrays where
  `polint`'s dummy is `xa(20)`/`ya(20)`). Tell: "Actual argument
  contains too few elements." Fix: widen the caller's declared array
  size to match the callee's real contract (check other correct
  callers, e.g. `rotation/fpft.f90`, to confirm the right size) --
  never shrink the callee's dummy to match the caller.
- **Scalar vs. length-1 array.** Caller passes a `(1)`-dimensioned
  array where the callee's dummy is a plain scalar, or vice versa
  (`wind/calcad.f90` did this twice, against `numerics/boole.f90` and
  against `splint`). Tell: "Rank mismatch ... (scalar and rank-1)."
  Both sides are legal, standalone F77 idioms (a length-1 array and a
  scalar share identical memory layout, so this always worked via
  sequence association under the old implicit interface) -- fix by
  passing the array *element* (`x(1)`) instead of the whole array at
  the call site, not by changing either declaration.
- **Intent mismatch.** Caller declares its own dummy argument
  `intent(in)`, then passes it as the actual argument to a callee
  whose corresponding dummy is `intent(inout)`/`intent(out)`
  (`atm/envint.f90` did this passing values through `bsstep`/`mmid` to
  an arbitrary `deriv` callback that might modify them). Tell: "Dummy
  argument ... with INTENT(IN) in variable definition context." Fix:
  first determine which side is actually right by checking whether the
  callee (or something *it* calls) ever assigns to the value -- if so,
  the callee's wider intent reflects real behavior and the caller's
  declaration should be widened to match (safe whenever every call
  site already passes a real variable, never a literal/expression, in
  that position -- check this before widening). Don't narrow the
  callee's intent to make the error go away without checking first.

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

## Scripting a COMMON-to-module conversion: two parser traps

Both of these produced silent, wrong results (not crashes) until
caught by manually inspecting a sample conversion before applying
broadly -- always do that inspection, don't trust a script that
reports 100% success without having read at least one real diff.

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

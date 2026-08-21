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
first pass for these. Instead:

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

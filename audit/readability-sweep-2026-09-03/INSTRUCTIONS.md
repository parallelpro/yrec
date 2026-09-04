You are doing a READ-ONLY readability/maintainability review of part of YREC, a Fortran stellar-evolution code at /Applications/YREC (git branch yrec-modern). Do NOT edit, create, move or delete any file in /Applications/YREC, do not run make, do not run any test or the yrec binary. Read files with the Read tool (whole files -- this is a review, not a search) and use grep/Bash only for read-only queries. Do not pipe exhaustive sweeps through `head`.

Background you need:
- The code was modernized in 2026 from F77 COMMON blocks to modules. Conventions: one model in `star` (type star_info, src/state/star_info_lib.f90), job/controls in `star%job` / `star%ctrl` (src/io/controls_lib.f90, src/state/controls_sync_lib.f90), physical constants in phys_const_lib, errors via `ierr` (positive = error, exit 1; negative = numerics termination), `stop` is driver policy only. Domains: core (solver/driver), net (nuclear rates), core/burn_lib (burn drivers), eos, kap, atm, mixing, rotation (seculr = angular-momentum transport, microdiff = element settling, shape), wind, setup, io, util (timestep), numerics (Numerical-Recipes-lineage kernels), math.
- src/GUIDELINES.md holds the ground rules; src/ROADMAP.md section 8 "Design debt" already lists: silently-static local arrays, blanket save + manual reset lists, -finit-local-zero masking, Numerical Recipes provenance, pulse column magic numbers, the write_gyre_pulse duplicate, a controls sanity checker, unused variables, rotation naming aliases, json=5000 static sizing, dlog/dexp archaic intrinsics. Sections 10 and 11 hold two prior file-by-file BUG sweeps. Do not re-report those items; this review is about CLARITY and MAINTAINABILITY, and about logic that does not make sense on reading.
- Every output-changing edit costs a full ~45-minute reseed of byte-pinned regression outputs, so classify each suggestion as (A) byte-safe restructuring (pure renames, comment fixes, extracting a subroutine with identical arithmetic, replacing goto with structured control flow with identical evaluation order, dead-code removal, named constants replacing magic numbers with the same value) or (B) changes the numbers (reordering floating-point operations, replacing an algorithm, changing tolerances).

What to look for, per file (read every file assigned to you, all of it):
1. Logic that does not make sense on reading: branches that cannot be taken, conditions that contradict their comment, variables computed and never used, values overwritten before use, suspicious index offsets, control flow that only makes sense as an F77 artefact (goto ladders, computed gotos, alternate returns, ENTRY-style flags, sentinel values like -999 threaded across many callers), comments that describe a different code than what is there. If something looks like an actual BUG, say so explicitly and separately -- but check sections 10/11 of ROADMAP.md first so you do not re-report a known one.
2. Readability: misleading or stale names, functions doing several unrelated jobs, deep nesting, duplicated code blocks that should be one routine (the three OPAL EOS vintages and the opal92/95 table families are KNOWN near-duplicates -- only report if you have a concrete, safe unification), magic numbers, flags that encode several meanings, argument lists where the caller cannot tell what is in/out, header comments that no longer match.
3. Maintainability: the top 3-5 structural changes in your area that would make new physics easier to add (e.g. a table of reactions instead of unrolled code, a registry, a named-index scheme, an explicit state type) -- concrete, with the files and the mechanism, and honest about cost.

Deliverable (your final message; be specific, cite file:line, keep prose tight, no filler):

## <area> -- summary (5-10 lines: what is in good shape, what is the main obstacle to understanding)

## Possible bugs / logic that does not make sense
- file:line -- what, why it looks wrong, what you would expect. (Say "not in ROADMAP 10/11" or "known, skip".)

## Readability findings (ranked, most valuable first, max ~25)
- file:line -- issue -- proposed change -- class A or B -- effort (S/M/L)

## Structural recommendations (3-5)
- what, where, mechanism, what it enables, cost, class A or B

## Files that are fine as they are
- list

Assigned files:

# YREC bug-identification sweep, PASS 2 — shared instructions

You are reviewing part of the YREC stellar-evolution code at
/Applications/YREC/src (branch yrec-modern). It was modernized from
F77 (COMMON blocks, gotos) to modern Fortran over many byte-gated
refactors. Your job: find REAL BUGS in your assigned files — not
style, not refactoring opportunities.

This is a SECOND, INDEPENDENT pass. A first sweep was done on
2026-08-31 with a different file partition; its findings are being
deliberately withheld from you so that agreement or disagreement
between the two passes is meaningful. Do NOT go looking for the
earlier reports (audit/bugsweep-2026-08-31/ or ROADMAP section 10) —
work from the code alone.

## What to look for

**Physical**: wrong formula (compare against the comment/paper it cites
or standard stellar-structure expressions), sign errors, wrong
constant values (check against CODATA/astro constants), unit
inconsistencies, wrong variable used in a formula (e.g. T where P
belongs), asymmetries between parallel branches that should mirror
each other (core vs envelope, up vs down, X vs Y vs Z, the three
OPAL-generation copies of a routine).

**Numerical**: float equality tests on computed values, catastrophic
cancellation, mixed single/double precision, division by possibly-zero
quantities, uninitialized variables read before set (trace the actual
flow — local variables that were SAVEd or in COMMON in the F77 and
are now plain locals are a recurring modernization hazard),
convergence tests that can't terminate, accumulation in the wrong
precision, log/exp of quantities that can be <= 0.

**Logical**: off-by-one in loop bounds or array indexing, array bounds
violations, index transposition (i,j vs j,i), copy-paste errors
between similar blocks (a classic in this codebase), inverted
conditions, dead code that was clearly MEANT to run, variables
assigned but the wrong one read back, save/state leaking between
calls, branches that can never be taken, GOTO-elimination mistakes
(an exit that became a return, a loop that lost its fall-through
label, an ELSE that vanished).

**Interface**: many routines are still bare external subroutines
(not module procedures), so a call site with the wrong number,
order, or type of actual arguments compiles silently. For every bare
external you review, check at least the call sites in your own
files: `grep -rn "call <name>" /Applications/YREC/src` (never pipe
that through head — read all of them) and compare the actual list
against the dummy list, including intent(out) arguments that the
caller then never reads and intent(in) arguments passed as
expressions or literals that the callee writes to.

## Provenance (required for every finding)

Most modernized files carry a header naming the original F77 file
("unchanged from the original seculr.f"). The pre-modernization
source is at git revision 6cd5673 (files src/<name>.f). For every
finding, run e.g.:
    cd /Applications/YREC && git show 6cd5673:src/seculr.f | grep -n ...
and classify:
- **inherited**: same defect present in the F77 original (cite it)
- **modernization**: introduced during conversion (the original did
  it differently — cite both)
- **new-code**: file/feature has no F77 ancestor
- **unclear**: could not locate the original site

If `git show 6cd5673:src/<name>.f` fails, list candidates with
`git show 6cd5673:src | tr ' ' '\n' | grep -i <part-of-name>`.

## Already fixed — do NOT re-report

Brunt N2 Ledoux term in pulse output; kap derivative factor
convention in GYRE columns; starin missing-ELSE MHD dispatch;
radiative_zone_bounds(1,2)=M typo; rezone flag_point overflow +
mod(dp,dp) ceiling idiom; eqstat dead SCV-derivative branch;
single-precision e-notation literals in burn/net/turnover;
condopacp + numerics_lib library stops (now ierr); intent bugs in
eqstat2/surfbc/liburn; registry defaults lfirst/niter3; profile
gradT/grada column swap AND the pulse-builder grad/grad_ad envelope
swap in stitched_model (both fixed); rotation_shape_factors callers
missing the new ierr argument (fixed); -Wsurprising static locals;
unused variables; NR-provenance licensing.

Known-but-undiagnosed SYMPTOMS (root-causing these IS valuable): the
run_from_zahb_to_tahb Solar_m1p0 GN93 case dies with NaN in the MLT
cubic (tpgrad sqrt-of-negative, DELDEL=NaN) near model 950, at core
helium exhaustion; a 0.8 Msun run died with "UNABLE TO SOLVE FOR NEW
ABUNDANCES IN SHELL 1 / RUN STOPPED AFTER 50 ATTEMPTS".

## Discipline

- READ-ONLY. Do not edit any file. Do not build. Do not run yrec.
  git commands must be read-only (show/log/grep).
- Read every assigned file IN FULL. For big files read in chunks.
- Be adversarial with yourself: before reporting, try to explain the
  code as CORRECT (historical convention, guarded upstream, value
  provably in range). Only report what survives. Quality over volume;
  cap at your ~15 strongest findings, mention weaker ones in one line
  each at the end.
- Comments lie. Trust the code, use comments as intent evidence.
- A mismatch between a comment's stated intent and the code IS a
  reportable finding.
- Where you can, state the OBSERVABLE consequence (which runs, which
  output quantity, roughly how large) — that is what decides fix
  priority.

## Output

Write your findings to the output file named in your task, format:

    ## <file>:<line> -- <one-line title>
    - class: physical | numerical | logical | interface
    - severity: high | medium | low   (impact if it fires)
    - confidence: high | medium | low (that it is really a bug)
    - provenance: inherited <orig.f evidence> | modernization | new-code | unclear
    - detail: 2-6 sentences: what the code does, what it should do,
      when it fires, and the provenance evidence.

Start the file with the list of files you read in full and a short
"verified clean" list of the algorithms you re-derived and found
correct (that negative evidence is valuable too). End the file with a
"## Weak/uncertain observations" one-liner list. Also return a brief
summary (counts by severity, and your top 3) as your final reply.

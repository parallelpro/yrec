# YREC bug-identification sweep — shared instructions

You are reviewing part of the YREC stellar-evolution code at
/Applications/YREC/src (branch yrec-modern). It was modernized from
F77 (COMMON blocks, gotos) to modern Fortran over many byte-gated
refactors. Your job: find REAL BUGS in your assigned files — not
style, not refactoring opportunities.

## What to look for

**Physical**: wrong formula (compare against the comment/paper it cites
or standard stellar-structure expressions), sign errors, wrong
constant values (check against CODATA/astro constants), unit
inconsistencies, wrong variable used in a formula (e.g. T where P
belongs), asymmetries between parallel branches that should mirror
each other (core vs envelope, up vs down, X vs Y vs Z).

**Numerical**: float equality tests on computed values, catastrophic
cancellation, mixed single/double precision, division by possibly-zero
quantities, uninitialized variables read before set (trace the actual
flow), convergence tests that can't terminate, accumulation in the
wrong precision.

**Logical**: off-by-one in loop bounds or array indexing, array bounds
violations, index transposition (i,j vs j,i), copy-paste errors
between similar blocks (a classic in this codebase), inverted
conditions, dead code that was clearly MEANT to run, variables
assigned but the wrong one read back, save/state leaking between
calls, branches that can never be taken.

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

## Known findings — do NOT re-report

Already found/fixed: Brunt N2 Ledoux term in pulse output; kap/eps
derivative factor convention in GYRE columns; starin missing-ELSE MHD
dispatch; radiative_zone_bounds(1,2)=M typo; rezone flag_point
overflow + mod(dp,dp) ceiling idiom; eqstat dead SCV-derivative
branch; single-precision e-notation literals in burn/net/turnover;
condopacp + numerics_lib library stops (now ierr); intent bugs in
eqstat2/surfbc/liburn; registry defaults lfirst/niter3; profile
gradT/grada column swap; -Wsurprising static locals; unused
variables; NR-provenance licensing; esac.f90's mistranslated
mass-fraction guard (OPAL95 EOS first call always errors — found by
the eos/opal sweep already). Known-but-undiagnosed SYMPTOMS
(root-causing these IS valuable): the run_from_zahb_to_tahb
Solar_m1p0 GN93 case dies with NaN in the MLT cubic (tpgrad
sqrt-of-negative, DELDEL=NaN) near model 950; a 0.8 Msun run died
with "UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL 1 / RUN STOPPED
AFTER 50 ATTEMPTS". Known open audit items you may EXAMINE and
sharpen but must label as known: float-equality guards
(burn_lib hydrogen_fraction.eq.0.0; microdiff species .eq. 0.0);
read_starting_model's ideal-gas core extension.

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

## Output

Write your findings to the output file named in your task, format:

    ## <file>:<line> -- <one-line title>
    - class: physical | numerical | logical
    - severity: high | medium | low   (impact if it fires)
    - confidence: high | medium | low (that it is really a bug)
    - provenance: inherited <orig.f evidence> | modernization | new-code | unclear
    - detail: 2-6 sentences: what the code does, what it should do,
      when it fires, and the provenance evidence.

End the file with a "## Weak/uncertain observations" one-liner list.
Also return a brief summary (counts by severity) as your final reply.

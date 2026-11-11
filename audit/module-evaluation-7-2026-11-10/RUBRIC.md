# The evaluation rubric (as used for the seventh evaluation, 2026-11-10)

Copied from the session notes so that later evaluations use the same text. Paths and the commit named
below are the seventh evaluation's; a later evaluation substitutes its own commit and scratch directory.

---

# Seventh evaluation (2026-11-10): every package of src/, read-only, out of 10

SP = /private/tmp/claude-501/-Applications-yrec/dbaf25ef-4aaa-4542-8cc3-27775b685d6f/scratchpad
MAIN = /Applications/yrec, branch yrec-next at 823326a7 (pushed to GitHub today; CI's first GitHub run is going).

## Rules
- READ ONLY. Never edit, build, run, `make` or git-write anything in MAIN or in any worktree; read files there
  (and `git log`/`git show`/`git grep` are fine). Scratch scripts and counts go in SP/review7_<area>/ (create it).
  Run any Python over the tree with `python3 -I`. Do not run the YREC binary or the test suite (other reviewers read
  the same tree at the same time). Git needs `export DEVELOPER_DIR=/Library/Developer/CommandLineTools`.
- Write your report to SP/notes/review7_<area>.md yourself, as the last step of your task, then hand back a short
  summary (scores, the top findings). Paths in the report are under src/ unless they start at the repository root;
  cite file:line at 823326a7.
- Use they/them for the owner and anyone else.

## The bar (what 10 means)
Scores 1-10 per criterion: 10 = at or above the standard below with nothing material left; 9 = only minor items left;
8 = a few real gaps; 7 = several; below 7 = structural problems. Be calibrated and consistent with the earlier
reports' scale (audit/module-evaluation-4-2026-10-31/ and module-evaluation-6-2026-11-02/): no grade inflation for
work done, no deflation for nitpicks. Score what the code IS, then explain each point lost.

The standard (src/rotation/ set it; rotation is scored against it too, as written, not against itself):
- README.md per directory: a "what lives where" table that matches the code; the flow (what calls what, in order);
  conventions (units CGS unless the name has a unit suffix, log is log10, errors return ierr, where messages go);
  what each test checks; a coverage figure; "how to extend" recipes that name EVERY file a newcomer must touch
  (verify at least one recipe by tracing it).
- Domains read controls from a settings record built in setup/, never s%control/s%physics directly; kernels are
  star-free (arguments in, results out) and listed in tools/check_boundaries.py STAR_BLIND_FILES where they qualify.
- Kernels print nothing; a domain's messages go through one module (a report = terminal + run log; a note).
  Failures return ierr (positive = error); no silent success.
- Numerical thresholds are named parameters with unit and meaning; physical constants are named, sourced,
  compile-time parameters.
- Tests check behaviour, not just pinned numbers: ok/FAIL checks (limits, conservation laws, scalings, exactness,
  continuity, error contracts), each confirmed to fail under a mutation of the code it covers.
- docs/GUIDELINES.md "Naming", "The star object: six rules", "Declarations and interfaces" (private by default,
  `use ..., only:`, real(dp), assumed-shape arrays, pure where possible) and "Code style" (named steps, each formula
  once, constants are parameters, named slots, comments describe the code as it is now, no ALL-CAPS prose, no
  commented-out code, `<`/`==` not `.lt.`/`.eq.`, ok/FAIL tests). Read those sections first.
- Commented-out code the authors kept lives in <dir>/retired/*.f90.txt with a header (where from, how to revive).
- Star contracts (`! star reads:` / `! star writes:` lines, tools/check_star_contracts.py) where the directory is
  required to have them.

## Inputs
- Your package's previous report, if any (listed in your task); its items' status is part of your report.
- docs/bug_sweep.md: open defects of your package weigh in the scores but are NOT re-reported (cite their numbers);
  report only findings not already there (grep it before you report one).
- docs/ROADMAP.md sections since the previous report; docs/TARGET_ARCHITECTURE.md where it describes your package.
- tools/check_boundaries.py, tools/check_star_contracts.py, tools/gen_star_docs.py: what they enforce for you.

## Report format (follow the sixth reports, audit/module-evaluation-6-2026-11-02/io.md, as the model)
1. `# <area> -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)` and a Scope paragraph (what you
   read whole, what skimmed, what you counted and how).
2. `## Scores (1-10)`: per criterion (readability, organization, maintainability, extensibility, docs), the
   history where there is one (earlier /5 scores doubled) -> now, one sentence why, and "Left:" bullets. The average.
   For two directories in one review, score each directory separately.
3. `## Summary` (a short paragraph or bullets).
4. `## Previous report's items: status` (DONE / PARTLY / OPEN, with where) -- only if there is a previous report.
5. `## New findings (not in bug_sweep)`: numbered; each with file:line, **Issue** (concrete; for a bug, the input that
   triggers it and what goes wrong; whether the F77 had it, with `git show 6cd56739:src/<file>.f` evidence),
   **Fix**, **Kind** (A = no result changes; B = results change), **Size** (S/M/L). Check each by reading the code
   twice; drop what you cannot confirm. Bugs first, then docs/structure.
6. `## What would raise each score` (per criterion, concrete items).

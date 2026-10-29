# Module evaluation, 2026-10-29 -- summary

Five read-only reviews of eos/, net/, io/, numerics/ with math/, and
state/ with setup/, each reading every file in its area in full, scored
against the standard src/rotation/ set after its Stage D (a README with
a what-lives-where table and extension recipes, star-free kernels with
mutation-checked physics tests, one message module, named thresholds).
Known debt in docs/ROADMAP.md and docs/bug_sweep.md was not re-reported.
The per-area reports are alongside this file; they are the reviewers'
words, condensed, with file:line references at yrec-next 2c88d896.

## Scores (out of 5, relative to rotation/)

| Module | Lines | Readability | Organization | Maintainability | Extensibility | Docs | Main obstacle |
|---|---|---|---|---|---|---|---|
| net | 4.8k | 4 | 4 | 3 | 3 | 3 | hand-kept index layouts around the reaction table |
| setup | 0.6k | 4 | 3 | 3.5 | 3 | 3 | map_user_inputs and net's builder live elsewhere |
| math | 0.2k | 4 | 3.5 | 3 | 3.5 | 2 | no README, no test |
| numerics | 2.4k | 3 | 4 | 2.5 | 3 | 3 | unchecked preconditions, ierr that reports success |
| eos | 8.2k | 3 | 3 | 2 | 2 | 3 | hidden state, golden-file-only test |
| io | 4.6k | 2.5 | 2.5 | 2 | 2 | 3 | the controls pipeline |
| state | 1.6k | 2.5 | 2 | 2.5 | 2.5 | 2.5 | star_info is a grab bag in COMMON-block order |

## Bugs checked against the source by hand

1. The MHD EOS reads no tables: open(newunit=) returns negative units and
   mhdst reads a table only if its unit is positive (eos/mhd/mhd_eos.f90).
2. The profile and GYRE eps omit the alpha captures: eps_total sums
   energy components 1-5 only (core/structure/structure_equations.f90).
3. bsstep never checks its limits; imax and nuse are inlist controls, and
   imax = 12 overruns the 11-entry substep_sequence.
4. The settling switch, a user control used as a latch, is never re-armed
   for a later card or Monte Carlo run.
5. chi_grid_scale (the zoning tolerances) is zero unless
   use_new_rate_variables is true; its default is false (every shipped
   deck sets it true).
6. write_mod_model always rewinds, so the Monte Carlo model file keeps
   only the last run's model.

The reports list further suspected bugs that were not checked by hand
when this summary was written; the fixes made from this evaluation, and
their verification, are in docs/ROADMAP.md and docs/bug_sweep.md.

## Themes across modules

- Physical "constants" assigned at run time (state/physical_constants.f90
  from core/run/run_setup.f90), so star-free domains need a booted star
  in their tests.
- Controls without a single source of truth: two default sources that
  disagree for 107 controls, the ctrl/job split in a Python list, five
  hand edits and three generators to add a control.
- Tests that pin numbers but check no physics (eos, numerics), which
  pinned existing bugs and missed the MHD breakage.
- No error and message contract in the kernels: printing to three places,
  ierr = 0 on failure.
- Index layouts kept by hand: net's scale slots, history's vals(22+i),
  profile cases that write zeros for a forgotten region.
- Documentation drift: every README had 5-10 false claims; the star
  member tables were hand-maintained and had lost 12 members.

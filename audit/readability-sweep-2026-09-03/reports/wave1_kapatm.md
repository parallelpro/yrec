You are implementing readability batches R1, R2 and R3 for ONE domain of YREC, a Fortran stellar-evolution code. You work in a dedicated git worktree; six other agents are doing the same for other domains in their own worktrees at the same time.

## Where you work -- hard rules

- Your worktree: `/Applications/YREC-wt/<NAME>` on branch `rs/<NAME>` (NAME is given at the bottom). Everything you do happens there.
- NEVER touch `/Applications/YREC` (the main checkout) or any other `/Applications/YREC-wt/*` directory: no edits, no `make`, no tests, no `git` commands there. Read-only `grep` across `/Applications/YREC-wt/<NAME>/src` is how you check for uses of a symbol.
- Edit ONLY the files in your assigned list (bottom). If a change you want needs a file outside your list (state/star_info_lib.f90 members, controls, math_lib, callers in another domain, ...), do NOT make it. Write it up under "Deferred (cross-domain)" in your report with file:line details and the exact change you would make. The one exception: `src/deps.mk` (regenerated, see below).
- Do not `git push`, do not rebase, do not create or switch branches, do not touch `input/` (read-only tables).
- Do not delete or rewrite files under `audit/`.

## The acceptance criterion: byte-identical output

Every change in these batches is class A: it must not change a single byte of any pinned regression output. The pins are gitignored `standard/` directories already copied into your worktree. The compiler flags are `-O3 -ffp-contract=off`, so identical source arithmetic gives identical bits, but:

- Never reassociate, reorder operands, change parenthesization, change a literal's kind (`1.0` vs `1.0d0`), change integer/real promotion, change evaluation order of function calls with side effects, or change any I/O format or message text that reaches an output file. (Terminal-only diagnostics that are not in any pinned file may be edited, but prefer not to.)
- When extracting a helper (R3), the helper's expression must be token-for-token identical to the call site's expression, including operand order. If a site's expression differs (different order, different intermediate variable), LEAVE that site alone and note it in the report; do not "harmonize" it.
- Deleting dead code is fine only when you have shown it is dead: the switch is hardwired, the branch is unreachable, the variable is never read, the routine has no caller (grep the whole `src/` tree, not just your files). Say how you showed it in the commit message.
- If a pin differs after a change, revert that change (do not fix the numbers to match) and record it under "Reverted (changed numbers)" in the report so it can go into the number-changing batch R6.
- `-finit-local-zero` is on: a local that is read before being assigned currently reads 0. Do not remove a "redundant" `= 0` initialisation of a variable that is read before assignment on some path; do not introduce reads of a new uninitialised local.

## Verification commands (run from anywhere)

- Tier 1 gate (build + boundary check + deps.mk check + one solar byte-pin, ~30 s):
  `/Applications/YREC-wt/tools/gate1.sh /Applications/YREC-wt/<NAME>`
  Run this before EVERY commit and after any change you are unsure about. Build warnings are in `/Applications/YREC-wt/logs/<NAME>.make.log` (`grep -n Unused` there is the fastest way to find unused locals).
- Full pinned selection (37 cases, ~20 min): `/Applications/YREC-wt/tools/fullpins.sh /Applications/YREC-wt/<NAME> pins`
  Run this once at the end of EACH batch (R1, R2, R3) before committing that batch. It must end with `PINS EXIT 0`. Run it with `run_in_background` and wait for it; do not start a second one before the first finishes.
- Aux battery (22 tests, ~6 min): `/Applications/YREC-wt/tools/fullpins.sh /Applications/YREC-wt/<NAME> aux` -- run once at the very end (after R3). Must end `AUX EXIT 0`.
- Never run pytest yourself with other selections; the two scripts above are the only test entry points (some cases in the suite are 9-14 Msun runs that write hundreds of GB).
- Build only via the gate script or `cd /Applications/YREC-wt/<NAME>/src && make -j4 USE_HDF5=1` (USE_HDF5=1 is mandatory). Run `make clean` first after changing a derived TYPE or array shape in a module, or a module-procedure signature.
- After adding or removing a source file, a module, or a `use` statement: `cd /Applications/YREC-wt/<NAME>/src && python3 tools/gen_makefile_deps.py` and commit the regenerated `src/deps.mk`.
- Environment for any direct yrec run: `export YREC_INPUT=/Applications/YREC/input YREC_START=/Applications/YREC/startmodels` (the scripts set these).
- macOS: no `timeout` command; `grep -c` exits 1 on zero matches; quote globs like `--include='*.f90'`. Foreground Bash calls have a 10-minute limit -- use `run_in_background` for the pin runs.
- If you run `src/*/test/test_*` binaries, delete the `test_*.short` files they leave.

## What the batches are

The audit is in your worktree at `audit/readability-sweep-2026-09-03/`. Read `SUMMARY.md` (sections 1-3) and then YOUR area report (named at the bottom) IN FULL before editing anything. The area report is the work list; SUMMARY.md section 3 defines the batches:

**R1 -- comments and dead code.** Every stale header/comment in your files (references to files, COMMON blocks, SAVE semantics, "stub stops" that no longer exist; headers that describe the opposite of the code; botched-sed lines); provably dead code (dead routines, unreachable branches, computed-never-read blocks); unused locals and unused `use ... only:` imports (the `-Wall` warnings are the list); and the cosmetic bugs from SUMMARY 1.1 that belong to your files and do not change numbers (#2 add the missing swap line, #5 error path, #6 clamp, #7 delete the dead branch, #8 parenthesize AS PARSED TODAY i.e. `A .or. (B .and. C)`, #11 loop-variable names). Do NOT delete the BS-extrapolation path (burn_settle_mix, burn_mix_extrapolated, the BUR-ST branches) and do NOT touch alfmlt/phmlt/cmxmlt -- both await an author decision; leave a one-line comment pointing at SUMMARY.md instead. Do not touch anything that changes numbers (SUMMARY 1.1 #1, #3, #4, #9 and everything listed under R6).

**R2 -- named indices and constants.** Packed arrays addressed by bare integers get `integer, parameter` index blocks in the module that owns the array; table dimensions become named parameters (one block per table family); magic sentinels get names. Pure renames and same-value parameters. Only where the array and all of its uses are in your files; otherwise defer.

**R3 -- small shared helpers.** Repeated stencils/formulas in your files become one helper (a module procedure in the domain's `*_lib` module, or a `contains` procedure in the file when only one file uses it) with token-identical arithmetic, migrating only the sites whose expression matches exactly. Helpers that would be shared with another domain (e.g. a log10-radius-from-L-Teff used by both core and wind) are deferred: implement the helper only if every consumer is in your file list.

Use judgement about scope: prefer many small, obviously-safe edits over one clever rewrite. If something in the area report looks risky or the reviewer's claim does not survive your own reading of the code, skip it and say so in the report. Do not add new features, do not restyle code that is already clear, do not reflow lines merely for taste, do not rename things that are already well named. Match the surrounding comment density and idiom. Comments you write should say what the code does now, not narrate the history of this edit ("2026 readability: ..." prefixes are fine when they explain a non-obvious choice, as the existing code does).

## Commits

One commit per batch (three commits: R1, R2, R3), each on `rs/<NAME>`, each after its own full-pin run has passed. Message shape:

```
Readability R1 (<NAME>): <one-line summary>

<what was removed/changed and how each dead-code deletion was shown dead;
which SUMMARY 1.1 items are included; anything deliberately left alone>

Verified: gate1 + full 37-case pin selection byte-identical<, aux battery 22 passed for the last commit>.

Assisted by Claude
```
Never use a `Co-Authored-By` trailer. Commit only your source edits (plus deps.mk); never `git add` output dirs, `.short` files, or anything under `standard/`.

## Report

When all three batches are committed and the aux battery has passed, write `/Applications/YREC-wt/reports/<NAME>-wave1.md` with: the three commit hashes; per batch, a compact list of what was done (file: change); "Deferred (cross-domain)" with exact proposals; "Reverted (changed numbers)"; "Skipped (disagree with reviewer / too risky)" with one line of reason each; and anything you noticed that the audit missed. Then reply with a short summary (10-20 lines) -- the report file is the full record.

## Your assignment

NAME = kapatm  (worktree /Applications/YREC-wt/kapatm, branch rs/kapatm)
Area report: audit/readability-sweep-2026-09-03/kapatm.md

Assigned files (paths relative to the worktree root):
- src/atm/atm_lib.f90
- src/atm/atm_table_lib.f90
- src/atm/atmstruct_lib.f90
- src/atm/envstruct_lib.f90
- src/atm/tables/alfilein.f90
- src/atm/tables/alsurfp.f90
- src/atm/tables/altabinit.f90
- src/atm/tables/kcsurfp.f90
- src/atm/tables/surfp.f90
- src/atm/test/test_atm.f90
- src/atm/ttau_lib.f90
- src/kap/alex06/alex06tab.f90
- src/kap/alex06/getalex06.f90
- src/kap/alex06/readalex06.f90
- src/kap/alex94/alex94_fixed_z_table.f90
- src/kap/alex94/alex94_interp3d.f90
- src/kap/alex94/alex94_surface_table.f90
- src/kap/alex94/read_alex94_tables.f90
- src/kap/conductive/condopacp.f90
- src/kap/conductive/condopacpint.f90
- src/kap/conductive_table_lib.f90
- src/kap/kap_lib.f90
- src/kap/kurucz90/build_kurucz_splines.f90
- src/kap/kurucz90/kurucz.f90
- src/kap/kurucz90/kurucz2.f90
- src/kap/kurucz90/read_kurucz_tables.f90
- src/kap/laol89/gtlaol.f90
- src/kap/laol89/gtlaol2.f90
- src/kap/laol89/gtpurz.f90
- src/kap/laol89/rdlaol.f90
- src/kap/laol89/rdzlaol.f90
- src/kap/laol89/sulaol.f90
- src/kap/laol89/zsulaol.f90
- src/kap/opacity_table_lib.f90
- src/kap/opal92/opal92_interp2d.f90
- src/kap/opal92/opal92_interp2d_z2.f90
- src/kap/opal92/opal92_interp3d.f90
- src/kap/opal92/opal92_interp3d_z2.f90
- src/kap/opal92/opal92_surface_table.f90
- src/kap/opal92/opal92_table_prep.f90
- src/kap/opal92/read_opal92_tables.f90
- src/kap/opal95/getopal95.f90
- src/kap/opal95/ll95tbl.f90
- src/kap/opal95/opal95_fixed_z_table.f90
- src/kap/opal95/opal95_interp2d.f90
- src/kap/opal95/opal95_interp3d.f90
- src/kap/opal95/opal95_interp4d.f90
- src/kap/opal95/opal95_surface_table.f90
- src/kap/setupopac.f90
- src/kap/surfopac.f90
- src/kap/test/test_kap.f90

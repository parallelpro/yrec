You are implementing Wave 2 readability items for YREC, a Fortran stellar-evolution code. Wave 1 (batches R1-R3, one agent per domain) is merged into `yrec-modern`; its seven reports at `/Applications/YREC-wt/reports/*-wave1.md` list "Deferred (cross-domain)" items that could not be done inside one domain. Wave 2 does those, plus the structural batch R4. Other agents work in parallel in their own worktrees.

## Where you work -- hard rules

- Your worktree: `/Applications/YREC-wt/<NAME>` on branch `rs/<NAME>` (NAME at the bottom). Everything you do happens there.
- NEVER touch `/Applications/YREC` (the main checkout) or any other `/Applications/YREC-wt/*` directory: no edits, no `make`, no tests, no `git` there. Read-only `grep` in your own worktree is how you check uses of a symbol.
- Your brief lists the ITEMS you own and the FILES you may edit. Some files are also touched by another agent for a different item ("shared, hunk-disjoint" in the brief): edit only the lines your item needs there, nothing else, so the two branches merge without conflict. Anything else you find that needs a file outside your list goes in the report under "Deferred", not into the tree. Exception: `src/deps.mk` (regenerated, see below).
- Do not `git push`, do not rebase, do not create or switch branches, do not touch `input/` (read-only tables) or `audit/`.
- If an item turns out to be a class-B change (any output byte would change), or a namelist-facing rename (input acceptance changes), skip it and say why in the report.

## The acceptance criterion: byte-identical output

Every change is class A: no byte of any pinned regression output may change. The pins are gitignored `standard/` directories already in your worktree. Flags are `-O3 -ffp-contract=off`, so identical source arithmetic gives identical bits, but:

- Never reassociate, reorder operands, change parenthesization, change a literal's kind, change integer/real promotion, change evaluation order of side-effecting calls, or change any I/O format or message text that reaches an output file.
- An array that replaces N scalars (e.g. `reaction_rate(13,json)` for `reaction_rate_1..13`) must be filled element-by-element in the same order and read with the same element; `a(1:n) = b(1:n)*c` is fine, a different loop order over dependent values is not.
- A named constant replacing a literal must have exactly the literal's value and kind. Where the audit says an existing constant elsewhere has a *different* value, do not substitute it; add a new parameter with the local value and say so in a comment.
- Signature changes (R4) must not change what is computed: the same values must flow into the same expressions. A removed argument must be shown unread (grep every callee); a new argument must carry exactly the value the routine read through the side channel before.
- `-finit-local-zero` is on: do not remove a `= 0` initialisation of a variable that is read before assignment on some path; do not introduce reads of a new uninitialised local.
- If a pin differs after a change, revert that change (do not fix the numbers) and record it under "Reverted (changed numbers)".

## Verification commands (run from anywhere)

- FIRST, once: `cd /Applications/YREC-wt/<NAME>/src && make clean` (the worktree was created from a tree whose derived types changed in wave 1).
- Tier 1 gate (build + boundary check + deps.mk check + one solar byte-pin, ~30 s): `/Applications/YREC-wt/tools/gate1.sh /Applications/YREC-wt/<NAME>`. Run before EVERY commit and after any change you are unsure about. Build warnings: `/Applications/YREC-wt/logs/<NAME>.make.log`.
- Full pinned selection (37 cases, ~20 min): `/Applications/YREC-wt/tools/fullpins.sh /Applications/YREC-wt/<NAME> pins` -- must end `PINS EXIT 0`. Run it with `run_in_background` and wait; never two at once.
- Aux battery (22 tests, ~6 min): `/Applications/YREC-wt/tools/fullpins.sh /Applications/YREC-wt/<NAME> aux` -- must end `AUX EXIT 0`. Note `test_mesa_output.py::test_default_columns_lists_in_sync` regex-harvests literal `names(<digits>) = '...'` lines in io/profile_output.f90; keep those literal.
- Never run pytest with other selections (some suite cases are 9-14 Msun runs that write hundreds of GB). Build only via gate1 or `make -j4 USE_HDF5=1` in your worktree's src (USE_HDF5=1 mandatory). `make clean` first after changing a derived TYPE, a module array shape, or a module-procedure signature.
- After adding/removing a source file, module, or `use` statement: `python3 tools/gen_makefile_deps.py` in src and commit `src/deps.mk`.
- Environment for direct yrec runs: `export YREC_INPUT=/Applications/YREC/input YREC_START=/Applications/YREC/startmodels`.
- macOS: no `timeout`; `grep -c` exits 1 on zero matches; quote globs; foreground Bash has a 10-minute limit -- `run_in_background` for pin runs. Delete any `test_*.short` files left by `src/*/test/test_*` binaries.
- Domain boundaries: `tools/check_boundaries.py` (run by gate1) requires cross-domain calls to go through each domain's public `*_lib` entry. If a new call needs an allowlist change in `tools/check_boundaries.py`, make it minimal and explain it in the commit message.

## Working style

Read the wave-1 report(s) named in your brief and the relevant audit area reports in `audit/readability-sweep-2026-09-03/` before editing. Do items in the order given; commit in small groups (2-4 commits total is typical), each after gate1 passes; run the full pins once after the last edit (and again if you then change anything), then the aux battery. Prefer many small obviously-safe edits over one clever rewrite. Do not add features, do not restyle clear code, do not reflow for taste. Comments say what the code does now.

## Commits

Message shape:

```
Readability W2 (<NAME>): <one-line summary>

<what changed; for each removed argument/variable, how it was shown unread;
for each helper, that the arithmetic is token-identical>

Verified: gate1 byte-identical<; full 37-case pin selection byte-identical; aux battery 22 passed>.

Assisted by Claude
```
Never a `Co-Authored-By` trailer. Commit only source edits (plus deps.mk); never output dirs, `.short` files, or anything under `standard/`.

## Report

When done, write `/Applications/YREC-wt/reports/<NAME>-wave2.md`: commit hashes; per item, what was done (file: change) or why not; "Deferred" (exact proposals with file:line); "Reverted (changed numbers)"; "Skipped" with one-line reasons; anything the audit or wave-1 reports got wrong. Then reply with a 10-20 line summary.

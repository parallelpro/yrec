# Eighth evaluation: report findings to bug_sweep numbers (registered 2026-11-17)

One row per report finding (split findings one row per part), from the "New findings (not in
bug_sweep)" sections of audit/module-evaluation-8-2026-11-17/*.md. "M319".."M355" are the new entries
of docs/bug_sweep.md §19 (branch eval8-docs). Kind: A byte-identical, B changes results or pinned
text, as the report says, corrected where the check found otherwise.

Every finding was read again in the code at 1179afbf, and against the F77 at 6cd56739, before its
row was written; four were also run on a scratch build (the "checked" column). The dispositions:
- **M-number**: a new entry of its own;
- **grouped into M-number**: the package's entry for stale or inexact text (M333, M337, M344, M347,
  M351, M355), or core's messages (M332);
- **already there as**: an existing entry;
- **not registered (style or structure)**: by bug_sweep §10's rule, the reports' style and structure
  points are not defects and get no entry; SUMMARY.md's "What is left under 8" carries them, as the
  seventh evaluation's PLAN carried its "plan items";
- **dropped**: not confirmed in the code. **None was dropped**: every defect the six reports list is
  in the code as described, apart from the corrections in the second table.

| report | finding | disposition | Kind | checked | one-line title |
|---|---|---|---|---|---|
| core_a | 1 | M319 | B | read; deck count | An accreting, disk-locked model older than the disk-locking age takes a negative step |
| core_a | 2 | M320 | B | read | A non-positive L at the core edge gives a negative burning limit, taken by `abs` |
| core_a | 3 | M321 (with atm 2) | B | read | The Allard and gray switch flips at every model when log g leaves the tables |
| core_a | 4 | M322 | A | read | The envelope kernel reads locals nothing set (`-finit-local-zero`) |
| core_a | 5 | M323 | B | read; numbers recomputed | Accretion's dead entropy path can still stop the run |
| core_a | 6 | grouped into M332 | B | read (every place) | Seventeen messages not in the stated form |
| core_a | 7 (the README's sentences) | grouped into M333 | A | read; manifests recounted | "Each check ... confirmed to fail under a mutation"; the coverage date; two stale lines |
| core_a | 7 (the missing mutations, mutate.py's second argument list) | not registered (style or structure) | A | counts recounted | Mutations for the radiative `qt` partials and the gray-photosphere checks |
| core_a | 8 | grouped into M333 | A | read | "Add a timestep limit" leaves out four things; two recipes missing |
| core_a | 9 (stale comments) | grouped into M333 | A | read (each) | Eight stale comments in structure/, step/ and core_messages |
| core_a | 9 (unnamed values, unused names, loose intents) | not registered (style or structure) | A | sampled | Thresholds and slots by number; 30 unused `only:` names; `intent(inout)` by habit |
| core_a | 10 (the README's "kernels print nothing") | grouped into M333 | A | read | `evaluate_energy_sources` reports from the shell kernel's call tree |
| core_a | 10 (the checker lists, the nine-line wrapper) | not registered (style or structure) | A | read | energy_sources.f90 and henyey_elimination.f90 outside STAR_BLIND_FILES |
| core_a | 11 | not registered (style or structure) | A | sampled | Gamma1, the point store, the mixed mean, the interpolation, ambient Teff, atomic weights written more than once |
| core_b | 1 | M324 | B | read; F77 nesting | The turnover time's innermost-point branch is overwritten (2026 regression) |
| core_b | 2 | M325 | B | read | Each card's first step runs with the turnover time zeroed |
| core_b | 3 | M326 | A (part a), B (part b) | read; deck count | The core extension: 0/0 for one added point; the density not shifted |
| core_b | 4 | M327 | A/B | read | The stitch swallows an envelope-integration failure |
| core_b | 5 | M328 | A | read | The star calibration's false crossing when the model starts above R* |
| core_b | 6a | M329 | A | read | The turnover walk reads flags beyond its last point |
| core_b | 6b | M330 | A | read | The rezone's added outer point without a json test |
| core_b | 6c | M331 | A | read | A NULL path from a C host reaches `strlen` |
| core_b | 7 (the statements) | grouped into M333 | A | read | "Nothing it computes feeds the physics" |
| core_b | 7 (moving the renormalization) | noted in M333 as its option | B | -- | The luminosity renormalization into the step |
| core_b | 8 | not registered (style or structure) | A | read (the shell masses) | One formula in two places: shell masses, the jump insertion, ambient Teff |
| core_b | 9 | grouped into M332 | B | read (every place) | Ten messages not in the form; four events on the terminal only |
| core_b | 10 | grouped into M333 | A | read (each) | Stale or inexact text in model/, run/ and diagnostics/ |
| core_b | 11 | not registered (style or structure) | A | sampled | Unnamed slots, thresholds and unsourced tables |
| core_b | 12 | not registered (style or structure) | A | counts | diagnostics/ and model/'s drivers have no test of their own |
| eos | 1 | M334 | A | read | The MHD reader reads past its arrays before checking the count; one anonymous message |
| eos | 2 | M335 | A | read | SCV rows' counts and the first OPAL block's T6 grid unchecked |
| eos | 3 | M336 | A (statement), B (fix) | read | SCV's two interpolations differ at the composition handover |
| eos | 4 | grouped into M337 | A | read | Three recipe gaps |
| eos | 5 | grouped into M337 | A | read | The per-path table's 1/mu_ion row |
| eos | 6 | grouped into M337 | A | read (each) | Stale comments and references |
| eos | 7 | not registered (style or structure) | A | not re-counted | 22 unused imports; four texts that repeat a named value |
| kap | 1 | M338 | A | read; deck count | kap's envelope Z never follows the model |
| kap | 2 | M339 | A | read | OPAL 1995's surface table X nodes |
| kap | 3 | M340 | A | run | The conductive opacity NaN for a gas with no electrons |
| kap | 4 | M341 | B | run | The conductive derivative where a species' density is clamped |
| kap | 5 | M342 | A | read | The Alexander lookups extrapolate with ierr 0 |
| kap | 6 | M343 | A | read | What `kap_init` still accepts |
| kap | 7 | not registered (style or structure) | A | not re-counted | Behaviour without a check: `open_table`, the high-R shift, the precedence, five readers |
| kap | 8 | grouped into M344 | A | read | A second copy of the molecular precedence in io; the recipe's omissions |
| kap | 9 | grouped into M344 | A | read (each) | Stale or inexact text |
| kap | 10 | not registered (style or structure) | A | sampled | `pure`, a loose intent, unused imports, unsourced atomic weights, long routines |
| atm | 1 | M345 | A | run | The Kurucz lookup has no per-row log g limit |
| atm | 2 | already there as M321 (core_a 3; one entry for both) | B | read | The switch back to the Allard tables asks log Teff only |
| atm | 3 | M346 | A | read | Unguarded paths of the Allard reader |
| atm | 4 | grouped into M347 | A | read (a search for readers) | A dead public entry and six dead fields |
| atm | 5 | grouped into M347 | A (B for generated text) | read | The choices listed by hand in two messages |
| atm | 6 | grouped into M347 | A | read (the formats) | The reader's log Teff format is one of two, chosen by row count |
| atm | 7 | grouped into M347 | A | read | Three stale README passages |
| atm | 8 (text items) | grouped into M347 | A | read | Argument names, `atm_init`'s order, test titles, a broken comment |
| atm | 8 (the literal 0, the literals 999, a public parameter) | not registered (style or structure) | A | read | Unnamed values and one export |
| convection_diffusion | 1 | M348 | A | read | `convective_flag(num_zones+1)` written past a json-shell model's arrays |
| convection_diffusion | 2 | M349 | B | read | The merge is not the union when the upper region reaches below the lower one |
| convection_diffusion | 3 | M352 | A | read | `change_zone` divides by Z |
| convection_diffusion | 4 | M353 | B (test pins) | read | The original routines' derivative factors |
| convection_diffusion | 5 | M354 | A | run | A one-zone settling region returned as found |
| convection_diffusion | 6 | M350 | B (change), A (document) | read | Semiconvection's first shell joins whatever the timestep |
| convection_diffusion | 7 | grouped into M351 | A | read (a search for callers) | `pressure_scale_height` exported without an outside caller |
| convection_diffusion | 8 | grouped into M351 | A | read | The coverage figure predates round 1B |
| convection_diffusion | 9 (TARGET_ARCHITECTURE.md:508, :553) | no entry: corrected in this registration's docs commit | A | read | Two stale rows of the routine table |
| convection_diffusion | 9 (M308) | already there as M308, whose status is corrected (fixed, ebf4d507) | A | read; the pin searched | M308 listed open though fixed |
| convection_diffusion | 10 | grouped into M355 | A | read | The settings table's cell; `solve_diffusion_term`'s dead work |
| convection_diffusion | 11 | not registered (style or structure) | A | not re-counted | Long tails of `intent(out)` results |

## What the check corrected

No finding was dropped. These are the places where the code, or a run, did not bear a report out in
full; each entry gives the corrected account.

| report | finding | what the report says | what the check found |
|---|---|---|---|
| core_a | 2 | the F77's core division is at ytime.f:68 | It is at ytime.f:82 (and xtime.f:38 for hydrogen). |
| core_a | 4 | of the 22 scalars copied into `eos_res`, "`log10_density` is the density the EOS starts from" | Only the fully ionized solver starts from the record's density, and only when it is not the default 0 (eos/yale/yale_eos.f90:785-788); the Saha solver, which an atmosphere's start reaches, starts from the record's `saha_*` fields. The unset `prev_tau`, `prev_density` and `prev_opacity` are as described. |
| core_b | 3 | after the extension `logRho(k)` for k > N is old point k's | True for N < k <= nz; the top N points hold whatever the array held beyond the old model, which the report does not say. |
| eos | 1 | "the 17 `ierr = 1` sites" of the MHD reader | The file holds 22 `ierr = 1` statements; the one message for all of them is as described. |
| kap | 6 | (b) and (c) are in the F77, citing rdlaol.f and setkrz.f | True, with other evidence: the tolerance's log is taken in GTLAOL and GTPURZ at each lookup (gtpurz.f:31), not in RDLAOL; SETKRZ does have an `IPUT.LT.1` stop, but inside its row loop, so a file without rows still ends with `NUMTM = 0` (setkrz.f:51-67). |
| convection_diffusion | 5 | with a one-zone region "the grids are built from NaN" | The one-zone region is real (run: zone_begin = zone_end = 3, `region_found`). `size_equal_grid` then takes `int(0/0)`, which is undefined: on the build here (gfortran, arm64) it gives 0, so the grid is 2 points at a spacing of 0, not NaN; another platform gives another count, and a devel build traps. |
| convection_diffusion | 9 | M308's pin holds "no Saha line" | Confirmed (0 lines at 1179afbf). M308's entry counted 12 such lines in its pin; the fixing commit ebf4d507 counts seven removed. |
| eos | 3 | the two SCV paths differ "by up to 5e-4 dex" | 5e-4 is test_eos's tolerance, and 3e-4 at log T 4.13 its comment's figure; neither was measured for the registration, and the entry says so. |

Two slips of docs/bug_sweep.md itself, found while counting for the registration and corrected in
it: M308's status (above), and the "Where the defects came from" table, whose fixed row gave 118
modernization items beside a list of 119 (the total, 289, was right).

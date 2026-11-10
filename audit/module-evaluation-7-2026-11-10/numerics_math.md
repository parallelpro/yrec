# numerics_math -- seventh evaluation (2026-11-10, yrec-next 823326a7, read-only)

Scope: I read whole every .f90 in numerics/ and math/, both READMEs, both test programs (test_numerics.f90, test_math.f90), their expected outputs, both mutation scripts, and test_numerics.py and test_math.py at the repository root. I also read:
- check_math_lib, check_no_io and the domain lists in tools/check_boundaries.py;
- the numerics rows of docs/HOW_YREC_WORKS.md §5;
- the CI workflow;
- docs/GUIDELINES.md "Naming", "Declarations and interfaces" and "Code style";
- the bug_sweep entries for this area (M99-M104, M124, M125, M132-M135, M142, M145, M158-M160).

Neither directory has changed since the fourth evaluation: the last commit touching src/numerics or src/math is 90628915 (2026-10-30). Since e8c847ba, the only change around them is in callers' `use` lines: io_units and run_log became run_messages. The question was therefore which items are still open, plus anything the fourth review missed.

Counts came from scratch Python run with `python3 -I` over the tree. The scripts are in SP/review7_numerics_math/:
- `mutcheck.py`: dead mutation strings;
- `checkcov.py`: which check each mutation covers;
- `calls.py` and `calls.txt`: the 92 call sites outside numerics/.

Nothing was built or run, and nothing in /Applications/yrec was edited. Paths are under src/ unless they start at the repository root.

## Scores (1-10)

Earlier 1-5 scores are doubled. Each line reads first -> second -> third -> fourth -> now.

**numerics**
- **Readability 6 -> 8 -> 9 -> 9 -> 9.** The headers give conventions and failure states, osplin's steps have names, and the Numerical Recipes transliterations read cleanly. Left:
  - about 85 lines of knot-case 3/4 branches that osplin never reaches: quadratic_spline.f90:218-255, :340-358, :462-489;
  - 27 inline history lines ("formerly", "went in 2026-10", "KC 2025-05-31", "llp 8/19/08");
  - cspline's literal 0.99d30 (cubic_spline.f90:66, :80);
  - ludcmp's parameter `tiny`, which shadows the intrinsic (linear_solvers.f90:48).
- **Organization 8 -> 8 -> 9 -> 9 -> 9.** Five topic modules, a status module and the facade; no private spline, dense solve or tridiagonal solve outside the domain kernels the README names. Left:
  - the five aliases are used only by test_numerics. A grep finds no caller outside numerics/;
  - the searches kept in eos/opal and kap (`stencil4_locate`, opacity_tables.f90:334, :370) are named in neither README §1 nor §3 (M160 in part).
- **Maintainability 5 -> 7 -> 8 -> 7.5 -> 7.5.** The pins, 64 checks and 71 mutations are strong, and CI now lists test_numerics.py (.github/workflows/ci.yml:126-131). Left:
  - CI's first GitHub run of that list is still pending;
  - M99-M103 and M132-M135 are open: wrong values or lost codes with ierr 0, all latent;
  - 23 dead mutation strings;
  - two checks have no mutation (finding 1);
  - two checks assert behaviour lir's interface does not promise (finding 3).
- **Extensibility 6 -> 8 -> 9 -> 9 -> 9.** The scheme recipe was traced, and the Makefile, deps.mk and the per-directory library pick up a new numerics/<noun>.f90. `ode_system` carries the physics. Left:
  - the HOW_YREC_WORKS table the recipe points to has drifted (M160);
  - the recipes leave out README §6's check list and the header lists that name every routine (finding 5d);
  - osplin's and lir's optional `ierr` lets a new caller drop the status without a compile error (finding 4).
- **Docs 6 -> 7 -> 9 -> 8 -> 8.** The facade lists, purity, the no-I/O rule, the message paths, bsstep's limits and the retry arithmetic (16 per retry, about 270 retries, 84 000 calls) all check out against the code. Left:
  - M100, M103, M135, M158, M159 and M160 are open;
  - README §6's "each check was confirmed to fail under a mutation" is false for 2 of 64 checks (finding 1);
  - smaller inexact statements (finding 5).

Average **8.5** (6.2 -> 7.6 -> 8.8 -> 8.5 -> 8.5).

**math**
- **Readability 8 -> 9 -> 10 -> 10 -> 10.** 172 lines in two mirrored `#ifdef` branches, each explained. M124's wrong comment counts under docs. Nothing material is left. Optionally, MESA's bound `100.0d0` (math_lib.f90:72, :127) could get a name.
- **Organization 7 -> 9 -> 10 -> 10 -> 10.** One module, a dependency root; the backend is chosen at build time, and the checker reads the shadowed names from the module's own `public` lines.
- **Maintainability 6 -> 8 -> 9 -> 8.5 -> 8.5.** Both backends are tested, with 17 mutations, and CI now lists test_math.py. Left:
  - CI's first run is pending;
  - M145 (pow_i's zero-base branch) and M142 (the checker's `**` hole) are open;
  - a second MESA-definition boundary has no check: the |n| < 100 bound of repeated multiplication under crmath (finding 2).
- **Extensibility 7 -> 9 -> 10 -> 10 -> 10.** The traced recipe for a shadowed function still holds: the checker picks up a new name by itself, and the recipe covers a function crmath lacks.
- **Docs 4 -> 9 -> 9 -> 8 -> 8.** Left:
  - M124's two comment statements: math_lib.f90:106-108's "matches x**y exactly" under crmath, and README:26's "a few ulp of y ln x";
  - M142's three README §3 statements, unchanged;
  - README §2's "each is checked by `test_math`" is not true for two rows (finding 2).

Average **9.3** (6.4 -> 8.8 -> 9.6 -> 9.3 -> 9.3).

## Summary

- **No change in either directory since the fourth evaluation.** Of its ten findings:
  - nine became bug_sweep entries (M132-M135, M142, M145, M158-M160), and all nine are open;
  - the tenth, the dead mutation strings, was left out of bug_sweep on purpose (bug_sweep.md:3287) and is unchanged (23 strings).
  - Of its "What would raise" items, only CI moved. The workflow runs both pytest files, and the branch was pushed today, but no GitHub run has finished. M125 is therefore PARTLY done.
- **The scores stand at 8.5 and 9.3.** Nothing regressed, and nothing was fixed. The new findings are small and of the same kinds the fourth review weighed:
  - test claims a little wider than the mutation script backs;
  - a MESA-definition bound without a check;
  - a few inexact sentences.
- **numerics is structurally at rotation's bar:**
  - one facade and named status codes;
  - pure routines, no I/O and no state, all enforced by check_boundaries (STAR_FREE_DIRS, NO_IO_DIRS, FACADE_ONLY_DOMAINS, ONLY_LIST_MODULES);
  - traced recipes.
  What keeps it below is its error contract, which is thinner than README §3 says (M134, M135, M133, M99-M101).
- **Callers' status handling** (92 call sites outside numerics/, tests excluded):
  - Every call that returns a status is checked: splint 20, osplin 8, lir 3, polint 7, ludcmp 1, simeqc 1, solve_tridiagonal 4, the forward sweep 1 and bsstep 2.
  - Exceptions already recorded:
    - none of the 23 natural_spline calls passes `ierr` (M99);
    - the envelope drops bsstep's `derivatives_ierr` and exits 0 (M132, envelope_integration_kernel.f90:542-556, :859-872, unchanged);
    - the network maps every simeqc code to "singular" (M102).
  - Message routing: callers in 12 files write `'<routine>: '//numerics_status_text(ierr)` directly to `run_log_unit` or `*` (31 writes), outside a domain message module. This is M112's class (direct run-log writes outside io), not a numerics defect. rotation (`rotation_numerics_report`) and SCV (`eos_report`) route theirs.
- **math's code is still at the bar.** Its open items are tests and statements: M145, M142, M124, and finding 2.

## Previous report's items: status

**Findings**
1. splint on a decreasing table that cspline accepts: **OPEN** (M134). cubic_spline.f90:120-155 is unchanged.
2. The envelope drops `derivatives_ierr`, and the allard header gives the wrong ierr: **OPEN** (M132). envelope_integration_kernel.f90:542-556, :859-872 (only its `use` lines changed); allard_atmosphere.f90:57-58.
3. math README §3 against the checker: **OPEN** (M142). check_boundaries.py:1218-1292, with the `**` regex at :1282; math/README.md:37.
4. README §3's "cheap preconditions are checked": **OPEN** (M135). numerics/README.md:43. osplin's table order, simeqc's num_cols, ludcmp's n <= np, polint's n and lir's sizes are still unchecked.
5. pow_i's zero base has no check: **OPEN** (M145). math_lib.f90:95-98; test_math.f90:115-123.
6. bsstep accepts a step with a NaN in y: **OPEN** (M133). bulirsch_stoer.f90:144-158; README:72.
7. numerics_status leaves out cspline's cases: **OPEN** (M158). numerics_status.f90:18-24, :53.
8. hunt's header for decreasing tables, and the `interval_arcs` name: **OPEN** (M159). table_interpolation.f90:369-371; quadratic_spline.f90:9-10.
9. The HOW_YREC_WORKS table lacks polint and lists `search`; §3 names no other domain's searches: **OPEN** (M160). docs/HOW_YREC_WORKS.md:538-551.
10. 23 dead strings in mutate_numerics.py: **OPEN**, not in bug_sweep by decision (bug_sweep.md:3287). Recounted: 23 strings in the 20 list-form mutations.

**What would raise (numerics)**
- Readability:
  - remove or test the knot-case 3/4 branches: **OPEN**;
  - move inline history to a Sources section: **OPEN**;
  - name 0.99d30: **OPEN**.
- Organization:
  - name the other domains' searches: **OPEN**;
  - decide on the aliases: **OPEN**. Still used only by test_numerics.f90:83-84.
- Maintainability:
  - fix M99-M103 and findings 1, 4 and 6 with checks and mutations: **OPEN**;
  - CI: **PARTLY** (ci.yml:126-131 lists test_numerics.py; the first GitHub run is pending);
  - delete the dead strings: **OPEN**;
  - update the coverage figure: n/a, since nothing was removed.
- Extensibility:
  - the precondition rule in "Add a scheme": **OPEN**. README:95 still says only "check its cheap preconditions";
  - HOW_YREC_WORKS in step with the README: **OPEN** (M160).
- Docs:
  - §3's sentence: **OPEN**;
  - hunt's header, M103 and M100: **OPEN**;
  - numerics_status: **OPEN**;
  - the quadratic_spline header name: **OPEN**;
  - the envelope comment: **OPEN** (envelope_integration_kernel.f90:546-548).

**What would raise (math)**
- integer zero-base checks with a mutation: **OPEN** (M145);
- the checker parses the whole exponent and is tested on synthetic cases: **OPEN** (M142);
- CI runs test_math: **PARTLY** (listed; first run pending);
- M124's statements: **OPEN**. math_lib.f90:106-108 and README:26. The CI part waits on M125;
- §3 matching the checker: **OPEN** (M142).

## New findings (not in bug_sweep)

1. **numerics/test/mutate_numerics.py:32-104 against numerics/README.md:90 and numerics/test/test_numerics.f90:25-26**
   - **Issue:**
     - Two of the 64 checks have no mutation: "splint: the slope's error falls as h^3 (order to 0.3)" (test_numerics.f90:624) and "interp: weights sum to 1, derivatives to 0" (:819). No mutation label is a prefix of either; `checkcov.py` maps every check to the labels the script matches with `startswith`.
     - Both the README ("Each check was confirmed to fail under a mutation of the code it covers ... 71 mutations, all caught") and the test header claim every check. For these two, nothing the script runs backs that claim. A rewrite of interp's derivative weights or of splint's `dy` that kept the sums or the order wrong would not be reported by the script.
     - The other 62 checks each have at least one mutation.
     - Related, docstring only: mutate_numerics.py:9-10 says "a nonzero exit ... counts as caught", but :150-151 counts a nonzero exit only when the check's line is missing. mutate_math.py:12-13 says this correctly.
   - **In the F77?** No; the checks and the script are 2026.
   - **Fix:**
     - add a mutation for each check and confirm it prints FAIL (`--only`). For the sums check, interp's `dweight(4) = ...` plus 1.0d-9. For the order check, a change to splint's `dy` term that the run shows fails it;
     - correct the docstring.
   - **Kind:** A. **Size:** S.
2. **math/math_lib.f90:72 and :127 (the bound `abs(y) < 100.0d0`, `abs(x) < 100.0d0`) against math/test/test_math.f90:47-68, :108-111 and math/README.md:25, :28**
   - **Issue:** under crmath, README §2's table says pow(x, n) takes repeated multiplication for an integer-valued |n| < 100, and so does exp10(n) for |n| < 100. It adds "each is checked by `test_math`". The checks reach only:
     - |n| <= 12 for pow;
     - |n| <= 22 for exp10;
     - n = 100 and 101, which are beyond the bound.
     So no check sees n in 13-99 (pow) or 23-99 (exp10). Lowering pow_r's bound to 13, or exp10's to 23, would leave all ten crmath checks ok, by tracing each check:
     - integer paths: n -12..12 and 0..22;
     - real paths and the single-precision exponent: non-integer y;
     - (-2)**100 and (-2)**101: the composition either way, within the check's tolerance;
     - 3e9 and the zero-base checks: untouched.
     The one mutation of the bound (mutate_math.py:55-57) moves it up, to 1e10, and is caught only by the 3e9 check. exp10's bound has no mutation at all.
   - **Why it matters:** the bound is MESA's definition. Below it, results are repeated multiplication; above it, the composition. The two differ in the last bits, so the bit-for-bit reproducibility the backend exists for depends on the bound.
   - **Also:** four of the ten intrinsic-backend checks have no intrinsic mutation (the three negative-base checks and y = 3e9). They test libm through the `x**y` passthrough, which "the intrinsic backend is x**y, bit for bit" already pins, so README §5's "each check" holds per label only under crmath.
   - **In the F77?** No; math_lib is 2026.
   - **Fix:**
     - under crmath, check pow(x, 99.0d0) and pow(x, -99.0d0) against the 99-fold product (and its reciprocal), and exp10(±99.0d0) against 99 multiplications by 10;
     - add mutations lowering each bound to 23 (19 mutations then);
     - optionally, name the bound (`repeated_multiplication_limit = 100`).
   - **Kind:** A (test only). **Size:** S.
3. **numerics/table_interpolation.f90:206-208, :216 against numerics/test/test_numerics.f90:937-941 (and the pins at :414-421)**
   - **Issue:**
     - lir's header says a failure "leaves result_y unset", and `result_y` is `intent(out)`, so it is undefined on entry. Two checks nonetheless assert that it is unchanged after a failure (`all(res == -1.0d0)`, in "lir: one point is numerics_bad_input, result untouched" and "lir: a flat table ..."). Two pins print it ("lir one point (no change)", "lir flat table (no change)").
     - They pass because gfortran leaves the array alone today. bug_sweep.md:3458 (M135) records that gfortran 15 at -O3 already treats polint's `intent(out)` scalar `dy` as clobbered, and test_numerics is built at -O3 (Makefile:193, :362).
     - The checks thus pin a compiler behaviour, not lir's contract. osplin, which promises "untouched", declares `yval` `intent(inout)` (quadratic_spline.f90:66).
   - **In the F77?** No; lir's `ierr` and these checks are 2026.
   - **Fix:** either make `result_y` `intent(inout)` and say "untouched" in the header (all three MHD calls pass a defined array), or drop the `all(res == -1.0d0)` conditions and the two "(no change)" pins.
   - **Kind:** A (byte-identical; the test output changes only if the pins are dropped). **Size:** S.
4. **numerics/quadratic_spline.f90:68 and numerics/table_interpolation.f90:220 (with numerics/README.md:43)**
   - **Issue:**
     - osplin's and lir's `ierr` are optional ("where it was added to an existing routine"), yet all 8 osplin calls and all 3 lir calls in the tree pass it.
     - Because the status is optional, a new caller that leaves it out gets a silent success: osplin leaves `yval` as it was, and lir leaves `result_y` undefined with the flag 0. That is the failure mode GUIDELINES "Code style" rule 7 forbids ("never reports a success it did not have"), and the compiler cannot catch it.
     - cspline and natural_spline have the same optional argument and 23 callers that omit it (M99).
   - **In the F77?** No (2026).
   - **Fix:**
     - make `ierr` required on osplin and lir now. This is byte-identical, and every call already compiles;
     - do the same for cspline and natural_spline together with M99's fix;
     - drop the README's "optional where ...".
   - **Kind:** A. **Size:** S.
5. **Docs and style, small**
   - a. **cubic_spline.f90:164-166.** "Column n holds the value, slope and half the second derivative ... (c(4,n) = 0), for a linear extrapolation beyond it." That is linear only when y2(n) = 0, the natural end both callers use. After cspline with a clamped end, c(3,n) ≠ 0 and the column is a quadratic (test_numerics.f90:590-593 checks exactly that). Fix: "the end tangent for a natural end".
   - b. **test_numerics.f90:116.** "bisection gives a zero-width interval only between equal last points". The same happens between equal first points for x below the table: xa = (1, 1, 2, 3) at x = 0.5 bisects to klo 1 and khi 2, giving coincident points (cubic_spline.f90:134-148).
   - c. **mutate_numerics.py:9-10.** See finding 1.
   - d. **numerics/README.md:94-101.** "Add an interpolation scheme" names README §1 and §2 but not:
     - §6, which lists every check;
     - the routine lists in numerics_lib.f90:7-16 and in the module's own header (for example table_interpolation.f90:6-14).
     "Add a solver" gives a new status code its text, export and `codes` entry, but not test_numerics' `use numerics_lib, only:` list, README §3's list of codes or HOW_YREC_WORKS.md:554-555's.
   - e. **Style:**
     - ludcmp's `real(dp), parameter :: tiny` shadows the intrinsic TINY (linear_solvers.f90:48);
     - `one = 1.0d0` (bulirsch_stoer.f90:92);
     - `double precision` in test_numerics.f90:1162, against GUIDELINES "Declarations and interfaces" rule 3.
   - **In the F77?** a, b, d and e no; e's `tiny` is Numerical Recipes' TINY, kept.
   - **Fix:** text only.
   - **Kind:** A. **Size:** S.

## What would raise each score

**numerics**
- **Readability 9 -> 10:**
  - remove the knot-case 3/4 branches osplin never reaches (byte-identical, as rotation's O15 did), or test `choose` and `cases` directly;
  - move the 27 history lines into a Sources section per file;
  - name 0.99d30, and rename `tiny`.
- **Organization 9 -> 10:**
  - name the OPAL EOS bisections and kap's stencil walks beside the domain solvers in README §1 (M160);
  - use the five aliases at call sites, or drop them.
- **Maintainability 7.5 -> 9+:**
  - fix M99-M103 and M132-M135, each with a check and a mutation;
  - a green first CI run (M125);
  - findings 1 and 3;
  - delete the 23 dead mutation strings.
- **Extensibility 9 -> 10:**
  - make `ierr` required (finding 4);
  - add the precondition rule and the missing files and lists to both recipes (finding 5d);
  - keep HOW_YREC_WORKS §5 in step (M160).
- **Docs 8 -> 10:**
  - M135's sentence, M100's and M103's headers, M158 and M159;
  - the envelope comment (M132);
  - finding 1's claim and finding 5's statements.

**math**
- **Maintainability 8.5 -> 10:**
  - M145 and finding 2 (checks and mutations for the zero-base branch and both bounds);
  - M142's checker fix, tested on `x**(2-a)`, `x**((a))`, `x**n` and a parameter initializer;
  - a green first CI run.
- **Docs 8 -> 10:**
  - M124's two statements;
  - make §3 match the checker (M142);
  - make §2's "each is checked" true (finding 2).
- **Readability, organization and extensibility** are at 10.

# numerics/ + math/ -- fourth evaluation (2026-10-31, yrec-next e8c847ba, read-only)

Read whole: every .f90 in src/numerics/ and src/math/, their test/ dirs (test programs, expected outputs, mutate_*.py), both READMEs, test_numerics.py and test_math.py at the root. The standard I used: src/rotation/README.md and the GUIDELINES sections "Naming", "Declarations and interfaces" and "Code style". I traced every caller of numerics_lib outside numerics/ and grepped src/ for private searches, splines, solves and interpolation. Scratch work is under SP/review4_num/:
- `build/probe.f90` and `build/probe2.f90`: numerics edge cases;
- `mathcr/`: a crmath build with pow_i mutated;
- `chk/` and `chk2/`: a copy of check_boundaries.py run on synthetic files.

Nothing in /Applications/yrec was edited. One caveat: the first checker run imported check_boundaries.py from src/tools, which may have refreshed the gitignored `src/tools/__pycache__/check_boundaries.cpython-314.pyc` (mtime 12:49). It is a regenerable cache, and git status is unchanged.

## Scores (1-10)

The earlier scores (1-5 scale) are doubled; each line reads first -> second -> third -> now.

**numerics**
- Readability: 6 -> 8 -> 9 -> **9**. The headers give units, conventions and failure states, osplin's steps have readable names, and the Numerical Recipes transliterations are clean. Two things cost the last point: about 85 lines of ACM-574 branches that osplin can never reach (knot cases 3 and 4, 28% of quadratic_spline), and history kept inline ("formerly ...", "went in 2026-10").
- Organization: 8 -> 8 -> 9 -> **9**. There are five topic modules, a status module and the facade, and solvers that carry a domain's physics live with their domain. Two gaps: the searches kept in other domains (10 OPAL EOS bisections, kap's two stencil walks) are not cross-referenced, and five aliases double the public names.
- Maintainability: 5 -> 7 -> 8 -> **7.5**. The pins, the checks and the 71 mutations (all present) are strong. Against that: CI runs none of it (M125); M99-M103 are open; and this review finds more preconditions that are not checked and give wrong values with ierr 0 (splint on a decreasing table, osplin on a decreasing table, simeqc without a right-hand side, polint at n = 1). The mutation script also carries 23 dead strings.
- Extensibility: 6 -> 8 -> 9 -> **9**. There are traced recipes for a scheme and a solver, `ode_system` carries the physics, and adding a status code is spelled out. One sign of drift: the HOW_YREC_WORKS table that step 4 says to maintain is already missing polint.
- Docs accuracy: 6 -> 7 -> 9 -> **8**. Most claims check out: the mutation count, the facade lists, purity, no use of other modules, the bsstep limits in controls_check, the message paths and the removed names. But §3's "the cheap preconditions are checked" overstates five routines. hunt's header is wrong for decreasing tables, the bad-input text is misleading, and the HOW_YREC_WORKS table is out of date. M100 and M103 (wrong headers) are still open.

**math**
- Readability: 8 -> 9 -> 10 -> **10**. 172 lines in two mirrored `#ifdef` branches, with every branch explained. M124's wrong comment counts under docs.
- Organization: 7 -> 9 -> 10 -> **10**. One module that is a dependency root; the backend is chosen at build time; the shadowed names are read from its public lines.
- Maintainability: 6 -> 8 -> 9 -> **8.5**. Both backends are tested, and there are 17 mutations. But one MESA-definition branch (pow_i's zero base) has no check: a scratch mutant passed every check. The checker that guards the usage rule misses `x**(2-a)` and `x**((a))`. CI does not run (M125).
- Extensibility: 7 -> 9 -> 10 -> **10**. The recipe for adding a shadowed function was traced, the checker picks up new names by itself, and a function crmath lacks is covered.
- Docs accuracy: 4 -> 9 -> 9 -> **8**. M124's three statements are still open. §3 adds three more that the checker contradicts: `x**n` "stays as `**`", parameter initializers are "exempt", and integer exponents are "exact repeated multiplication" (§2's own table says power by squaring). That is six wrong statements in a 66-line README.

## Summary

- Neither directory changed after the third evaluation. Its items M99-M104, M124 and M125 are all still open, and they count against both scores.
- numerics is structurally at rotation's bar: a facade, named status codes, pure routines, no I/O, no state, traced recipes and mutation-checked checks.
- Outside the documented M99 and M102, every caller handles the status of splint, osplin, lir, polint, ludcmp, simeqc and the Thomas solves. The exception is the envelope. It drops bsstep's `derivatives_ierr` and turns an EOS or opacity failure at a midpoint substep into the exit-0 "solution diverged" stop (finding 2).
- What keeps numerics below the third evaluation's maintainability and docs marks: its error contract is thinner than §3 says. Several preconditions are not checked and give plausible wrong values with ierr 0, on top of M99-M101:
  - splint on a decreasing table that cspline accepts;
  - osplin on a decreasing table;
  - simeqc with no right-hand side;
  - polint at n = 1.
- math's code is still at the bar. Its README's usage rule (§3) describes a checker other than the one that runs: `x**n` is refused, `**` in a parameter is not exempt, and `x**(2-a)` is not caught. pow_i's zero-base branch has no check.
- No private splines or dense or tridiagonal solves exist outside numerics beyond the domain kernels the docs name. Private searches remain in eos/opal (10 bisections) and kap (2 stencil walks), and nothing links them to numerics' list of searches.

## New findings (not in bug_sweep), ranked

1. **src/numerics/cubic_spline.f90:34-35, :61, :113-141.**
   - Issue: cspline and natural_spline accept a strictly decreasing table, and test_numerics pins that ("a decreasing table, the increasing one's y2", test_numerics.f90:643-647). But splint, the only routine that evaluates those second derivatives, bisects as if the table increased (:134-141). It checks nothing, and its own header does not say it needs an increasing table. Scratch probe2: x = (4,3,2,1,0), y = x² gives 8.375 at 0.5 (the same spline on the increasing table: 0.339) and 20.375 at 3.5 (12.339), both with ierr 0.
   - Latent: every caller's table increases (rezoning's chi and log m, the kap and atm grids, the turnover time's r − r_bcz − H_p).
   - Fix: splint returns numerics_bad_input when xa(n) < xa(1) (one comparison), or bisects in either direction as locate does. Or cspline refuses a decreasing table. Add a check and a mutation.
   - A, S.
2. **src/core/envelope/envelope_integration_kernel.f90:544-557, :861-872** (a caller of numerics).
   - Issue: both bsstep calls leave out `derivatives_ierr`. Every failure other than bad input becomes envelope_termination (−2), which main.f90:20-23 exits 0 on as the historical "solution diverged". So an EOS or opacity failure inside a midpoint substep (envelope_integrands.f90:109-117, :224-232):
     - prints "bsstep: derivatives failed",
     - loses its code,
     - and ends the run with exit 0.
   - The same failure at the step's end (:881-884) returns its positive code and exits 1. The comment at :546-549 ("the solution diverged", "the historical message") describes neither path.
   - Related (atm): allard_atmosphere.f90:57-58 says ierr is 1 on a polint failure, but the routine returns polint's code, 2.
   - Fix: pass `derivatives_ierr` and return it when jerr == numerics_derivatives_failed. Keep envelope_termination for numerics_step_underflow, the real divergence.
   - A: converging runs are unchanged; a run that fails this way exits 1 instead of 0. S.
3. **src/math/README.md:37 vs src/tools/check_boundaries.py:903-913.** The usage rule as written is not the one enforced. Each point below was confirmed with a copy of the checker on scratch files.
   - "Integer exponents (`x**2`, `x**n`) ... stay as `**`": the checker refuses `x**n` for an integer variable, because any non-digit token after `**` counts as a real exponent.
   - "Parameter initializers are exempt": that exemption applies only to the elementary-function test (`code_kept`). The `**` test scans `code_lines`, so `real(dp), parameter :: c = 2.0d0**0.5d0` is refused.
   - The `**` test reads only the first token after `**(`, so a real exponent written `x**(2-a)` or `x**((a))` passes. That is a hidden libm pow that bypasses crmath. None exists today; bulirsch_stoer.f90:170's integer `2.0d0**((...)/2)` gets through the same way.
   - "Exact repeated multiplication" contradicts §2's table, which says the compiler's power by squaring.
   - Fix: say that integer variables go through `pow(x, n)` and that parameters are not exempt from the `**` rule (or exempt them in the checker). Make the `**` test scan the whole parenthesized exponent and accept only integer literals and names declared integer.
   - A, S/M.
4. **src/numerics/README.md:43.** "The cheap preconditions are checked: sizes, ..., the order of a spline's abscissae" holds for cspline, osplin's evaluation points, lir's point count and bsstep. It does not hold for these, each confirmed on a scratch build or by reading:
   - osplin's table order (quadratic_spline.f90:53-54, :63-91): x = (4,3,2,1,0), y = x² gives 2.875 at 1.5 and 6.375 at 2.5 (x²: 2.25 and 6.25), with ierr 0.
   - simeqc's num_cols > num_unknowns (linear_solvers.f90:194-203): a 2x2 call with num_cols = 2 returns ierr 0 and overwrites A's last column with a "solution".
   - ludcmp's n ≤ np (:47-51).
   - polint's n ≥ 1 (table_interpolation.f90:137-189): at n = 1, dy (intent(out)) is never set, and ierr is 0.
   - lir's sizes (:228): it compares y_stride, which it never otherwise uses, with size(table_y, 1). It does not compare num_y, which indexes table_y and result_y, or num_points with table_z and table_y.
   - All latent: the callers pass consistent sizes and increasing tables.
   - Fix: add each check (one comparison each) with a check and a mutation, or narrow the sentence to the checks that exist. Drop lir's y_stride.
   - A, S.
5. **src/math/math_lib.f90:95-98.**
   - Issue: under crmath, pow_i's x == 0 branch is what makes README §2's pow(0, 0) = 0 and pow(0, y < 0) = 0 hold for integer exponents, and no check covers it. test_math's y ≤ 0 cases use real exponents only (test_math.f90:117-123), and pow(0.0d0, 3) gives 0 even without the branch. On a scratch crmath build with the branch deleted, pow(0, 0) = 1 and pow(0, -1) = +Inf for integer exponents, and all ten crmath checks still print ok.
   - Fix: add pow(0.0d0, 0) and pow(0.0d0, -1) to each backend's zero-base check (crmath 0 and 0; intrinsic 1 and +Inf), and a mutation of the branch in mutate_math.py.
   - A, S.
6. **src/numerics/bulirsch_stoer.f90:144-158.**
   - Issue: a step is accepted when the components that are not NaN meet the tolerance, so bsstep returns ierr 0 with a NaN in y. This is pinned (expected_test_numerics.out:250, "NaN in component 3 of 3: ierr 0 ... y(3) is NaN T") and stated in README §4:72. It is a success it did not have (GUIDELINES "Code style" rule 7). The envelope catches it only because its next derivative call fails.
   - Fix: keep the NaN-skipping maximum (the pinned 0.03 Msun one-equation case depends on it, and it never accepts a step), but after acceptance return numerics_bad_input, with y restored, if any y is NaN.
   - A on the shipped decks if no accepted envelope step carries a NaN (not measured); otherwise it changes how those failing runs end. S.
7. **src/numerics/numerics_status.f90:18-24, :53.**
   - Issue: the comments on the codes leave out cspline's cases (coincident for equal neighbours; bad input for n < 2 and for abscissae out of order). The text for numerics_bad_input is "bad input (a size or limit out of range, too few points, or NaN)", so a caller's message for cspline's zigzag table names three causes, none of which applies.
   - Fix: list cspline in the comments and add "abscissae out of order" to the text; regenerate the expected output.
   - A (one message text), S.
8. **src/numerics/table_interpolation.f90:369-371; src/numerics/quadratic_spline.f90:9-10.**
   - Issue: hunt's "Outside the table it returns 0 (before the first point) or n (at or beyond the last)" is true only for an increasing table. For a decreasing one, x = table(1) gives 0 and x = table(n) gives n−1 (scratch probe; test_numerics.f90:1037-1042 encodes that). Separately, quadratic_spline's header calls the function `interval_arcs`; that is the type, and the function is `interval_arcs_of`.
   - Fix: state the decreasing case in hunt's header; correct the name.
   - A, S.
9. **docs/HOW_YREC_WORKS.md:538-551; src/numerics/README.md:45.**
   - Issue: the numerics table, which README step 4 tells authors to update, has no row for polint (used by the Allard atmosphere and by rotation's distortion-integral extrapolation, rotating_shape.f90:126). Its row at :549 lists the private `search` as a public routine. Separately, README §3's "One bisection (`locate`) and one warm-started search (`hunt`)" does not mention:
     - the ten private bisections in eos/opal (opal95_eos.f90:144, 156, 398, 423, 443; opal_eos.f90:295, 306, 538, 564, 586), each with its own 1e-7 slack and early exit on equality;
     - kap's two stencil walks (opacity_tables.f90:334, 370), which §1's list of what lives with its domain does not name either.
   - Fix: add polint and drop search from the table; list the searches kept in other domains (the authors' code) in §1.
   - A, S.
10. **src/numerics/test/mutate_numerics.py:46-91.**
    - Issue: 20 of the 71 mutations carry lists of "the forms the code has taken", and 23 of those strings no longer occur in the code (for example `system_matrix(ib)=...`, `if (err == 2 .and. present(ierr))`, `weight(closest)=...`). `--check` passes as long as one form matches, so a dead form is never noticed.
    - Fix: keep only the current form of each, so that `--check` guards every string.
    - A, S.

## What would raise each score (concrete)

**numerics**
- Readability 9 -> 10:
  - Remove the knot-case 3/4 branches of `choose`, `cases` and `spline` that osplin cannot reach (byte-identical; rotation did the same with its unreachable branches, O15). Or keep them and test `choose` and `cases` directly.
  - Move the inline history ("formerly ...", "went in 2026-10", "(2026-10: on a 2-D array ...)") into a Sources section at the end of each file.
  - Give cspline's 0.99d30 natural-end threshold a name, as natural_spline names its 1.0d30.
- Organization 9 -> 10:
  - In §1 and §3, name the searches kept in other domains (the OPAL EOS's bisections, kap's stencil walks, osplin's internal search) beside the solvers that live with their domains. Or move the OPAL bisections onto locate and hunt in a measured commit.
  - Decide on the five aliases: use them at call sites, or drop them.
- Maintainability 7.5 -> 9+:
  - Fix M99-M103 and findings 1, 4 and 6, each with a check and a mutation: splint's direction, osplin's table order, simeqc's num_cols, polint's n, lir's num_y and num_points (and drop y_stride).
  - Get CI to build test_numerics and run test_numerics.py (M125).
  - Delete the 23 dead mutation strings.
  - Bring the coverage figure up to date after the dead branches go.
- Extensibility 9 -> 10:
  - Add a rule to the "Add a scheme" recipe: check every precondition the routine relies on, or name the unchecked ones in its header.
  - Keep HOW_YREC_WORKS §5 in step with the README (add polint, drop search).
- Docs 8 -> 10:
  - Rewrite §3's preconditions sentence to match the code.
  - Fix hunt's header for decreasing tables, plus M103's NaN sentence and M100's "each point gets the value it has alone".
  - Complete numerics_status's comments and the bad-input text.
  - Correct the header name in quadratic_spline.
  - Fix the envelope comment that calls every bsstep failure divergence (finding 2).

**math**
- Maintainability 8.5 -> 10:
  - Add integer-exponent zero-base checks per backend, with a mutation of pow_i's x == 0 branch.
  - Make the checker's `**` test parse the whole exponent (finding 3), and test the checker on `x**(2-a)`, `x**((a))`, `x**n` and a parameter initializer.
  - Make CI run test_math (M125).
- Docs 8 -> 10:
  - Fix M124's three statements.
  - Make §3 match the checker, or the checker match §3: integer variables, parameter initializers, and "power by squaring" in place of "exact repeated multiplication".
- Readability, organization and extensibility are at 10; nothing material is left there.

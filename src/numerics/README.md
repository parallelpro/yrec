# numerics/ — the numerical library

This directory holds the general-purpose numerical methods the physics uses: splines, interpolation in tables, searches, linear solvers and an ODE integrator. Code outside numerics/ uses it only through `numerics_lib`, and names what it takes (`use numerics_lib, only: ...`); `tools/check_boundaries.py` enforces both. numerics reads nothing from `star` (`STAR_FREE_DIRS`), keeps no state from one call to the next, does no I/O, and knows no physics: a routine that needs a caller's physics, such as an integrand, receives it as an object. It uses no other YREC module (not even `math_lib`, since no routine here calls an elementary function).

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `numerics_lib` (`numerics_lib.f90`) | **The public interface.** It re-exports what the modules below make public, the status codes and the bsstep limits, and gives the YREC-written routines with short legacy names a descriptive alias (the same procedure under a second name): `lagrange3_weights` (`inter3`), `lagrange4_weights` (`interp`), `equidistant_quadratic` (`quint`), `interpolate_tabulated_functions` (`lir`), `bracket_index` (`findex`). |
| `numerics_status` | The status codes every routine returns in `ierr` (§3), and `numerics_status_text`, a short description of each for a caller's message. |
| `cubic_spline` | `cspline` and `splint` (Numerical Recipes): the second derivatives of the cubic spline for given end slopes, and its value at x, with its first derivative when asked (`dy=`). `natural_spline`: `cspline` with both ends natural. `compute_spline_coefficients`: the spline's cubic on each interval in powers of x − x(i), for tables that store it (OPAL92 and Kurucz 1990 opacities). |
| `quadratic_spline` | `osplin`: the shape-preserving osculatory quadratic spline of ACM Algorithm 574 (McAllister & Roulier 1981), which keeps the data's monotonicity and convexity. Its steps are private: `slopes`, the slope at each point; `meval`, the evaluation in three steps (below, inside and above the table) around one interval's arcs, `interval_arcs_of` (`choose`'s knot case and `cases`' knots), which `spline` evaluates; `search` and `safedivide`. |
| `table_interpolation` | Lagrange weights for 3 and 4 points with their derivatives (`inter3`, `interp`), applied by `lagrange4`; `quint` for three equidistant points; Neville's `polint` for any number; `lir`, cubic interpolation of several tabulated functions at once; and the searches `locate` (bisection) and `hunt` (warm-started, Numerical Recipes), with `findex`, hunt with YCK's codes outside the grid (−1, −n). |
| `linear_solvers` | `ludcmp`/`lubksb` (LU with partial pivoting, Numerical Recipes), `simeqc` (Gaussian elimination for one or more right-hand sides, IBM SSP's SIMQ), `solve_tridiagonal` (the Thomas algorithm), and its two halves `tridiagonal_forward_sweep` and `tridiagonal_back_substitution`. |
| `bulirsch_stoer` | `bsstep`: one adaptive Bulirsch–Stoer step (modified midpoint rule with rational extrapolation, Numerical Recipes) of an `ode_system`; its limits `bsstep_max_equations`, `bsstep_max_stages`, `bsstep_max_extrap_order`. |

Solvers that carry a domain's physics in their interface live with that domain:
- the shape-factor quadratures are in `rotation/shape/shell_quadrature.f90`;
- the angular-momentum solve, `solve_angular_momentum_transport` in `rotation/seculr/am_transport_coefficients.f90`, is the Thomas algorithm with a step between its halves (it calls `tridiagonal_forward_sweep` and `tridiagonal_back_substitution`);
- the MHD table lookup is `interpolate_mhd_table` in `eos/mhd/mhd_eos.f90`.

## 2. Which routine for what

| To ... | use | notes |
|---|---|---|
| interpolate a smooth function through tabulated points | `natural_spline` (or `cspline` with known end slopes), then `splint` | C², exact on cubics with the right end slopes; `splint(..., dy=)` gives the slope; extrapolates the end cubic |
| store a spline as polynomial coefficients | `compute_spline_coefficients` after `natural_spline` | c(1:4, i) on [x(i), x(i+1)]; column n holds the end value and slope |
| interpolate where overshoots must not appear (monotone or convex data) | `osplin` | C¹, keeps monotonicity and convexity; at many sorted points in one call |
| interpolate with 4 points (or 3) and want the derivative | `interp` = `lagrange4_weights` (`inter3` = `lagrange3_weights`), then `lagrange4` | exact on cubics (quadratics); the weights' derivatives with `dweight` |
| interpolate at three equidistant points | `quint` = `equidistant_quadratic` | exact on quadratics |
| interpolate (or extrapolate to zero) with n points and an error estimate | `polint` | Neville; exact on polynomials of degree n − 1 |
| interpolate many functions tabulated on the same points | `lir` = `interpolate_tabulated_functions` | cubic (4 points or more) or linear |
| find the interval holding x | `locate` (bisection), `hunt` or `findex` (start from a guess: for a sequence of nearby x) | the conventions differ at the nodes and outside: see each header |
| solve a dense linear system | `ludcmp` then `lubksb` (several right-hand sides: one `lubksb` each) | `simeqc` solves in place for a few right-hand sides at once |
| solve a tridiagonal system | `solve_tridiagonal` | `tridiagonal_forward_sweep` + `tridiagonal_back_substitution` when a step goes between the two |
| integrate ODEs | `bsstep` in your stepping loop, on an `ode_system` (§4) | smooth right-hand sides, tight tolerances |

## 3. Conventions

- Real kind `dp` (`real64`). Arrays are sized by the arguments (explicit shape with a size argument, as the Numerical Recipes routines have it, or assumed shape); nothing is dimensioned by the model's maximum number of zones.
- **The error contract.** A routine that can fail has an `ierr` argument (optional where it was added to an existing routine: `osplin`, `lir`, `cspline` and `natural_spline`) and returns one of `numerics_status`'s codes: `numerics_ok` (0), `numerics_singular` (for `ludcmp` and `simeqc` a matrix singular to working precision: a zero row, or no pivot candidate above `singular_pivot_tolerance`, 1e-14, of its row's largest element; for the Thomas algorithm, which does not pivot, a zero pivot), `numerics_coincident_points` (two abscissae equal where they must differ, `lir`'s flat table, `osplin`'s end interval of zero width under a point beyond it), `numerics_bad_input` (sizes or limits out of range, too few points, abscissae not monotone, a NaN pivot), `numerics_step_underflow` (bsstep), `numerics_unsorted_points` (osplin), `numerics_derivatives_failed` (bsstep; the derivatives' own code in `derivatives_ierr`). The cheap preconditions are checked: sizes, bsstep's limits, the first pivot, the order of a spline's abscissae (`cspline`: strictly monotone; `osplin`: sorted points). A routine never stops the run; what its outputs hold after a failure is in its header (`osplin` leaves `yval` untouched, `cspline` does not set `y2`, `bsstep` restores `y`).
- **No I/O.** numerics prints nothing and never stops the run (`tools/check_boundaries.py` checks that numerics/ and math/ hold no I/O statement and no `stop`): the caller reports a failure, where the failure means something, as `<routine>: <numerics_status_text(ierr)>` (`splint: coincident abscissae`), and passes the code up. The messages go where the routines' own used to: splint's to the run log (eos's, in SCV's envelope table, to the terminal and the run log through `eos_report`, as its lir failures in the MHD lookup are), ludcmp's, solve_tridiagonal's, polint's and bsstep's to the terminal; rotation's through `rotation_numerics_report` (the terminal and the run log), rezone's (osplin) to both.
- **Searches.** One bisection (`locate`, with its 1% slack beyond the ends) and one warm-started search (`hunt`; `findex` is hunt with −1 and −n for outside, which OPAL92 and Kurucz 1990 test). `locate`, `hunt`, and the bisections inside `splint` and `osplin`, differ at the nodes and outside the table; each header states its convention.
- Procedures without side effects are `pure`: all but `bsstep` and `mmid`, which call an `ode_system`'s derivatives (they may change their object).
- Numerical Recipes and IBM SSP routines keep their authors' names (`cspline`, `splint`, `polint`, `ludcmp`, `lubksb`, `hunt`, `locate`, `bsstep` with `mmid` and `ratext`; `simeqc`). The routines written for YREC since 2026 are named for what they do (`natural_spline`, `compute_spline_coefficients`, `lagrange4`, `solve_tridiagonal`); the older ones with short names (`inter3`, `interp`, `quint`, `lir`, `findex`) have a descriptive alias in `numerics_lib` (§1), which new code may use, and the existing callers keep the short names; `osplin` keeps ACM 574's name.

## 4. Integrating ODEs

`bsstep` integrates any system you describe as an extension of the abstract type `ode_system`. Its deferred procedure `derivatives(self, x, y, dydx, ierr)` returns dy/dx. The extension carries whatever the derivatives need. The envelope integration's `envelope_equations` (`core/envelope/envelope_integration_kernel.f90`) holds the luminosity, rotation factors, composition, the atmosphere's current point and the EOS's Saha iteration counter, and calls `compute_atmosphere_derivs` or `compute_envelope_derivs`. To integrate a new system:

```fortran
type, extends(ode_system) :: my_system
   real(dp) :: k                     ! whatever the derivatives need
contains
   procedure :: derivatives => my_derivatives
end type

subroutine my_derivatives(self, x, y, dydx, ierr)
   class(my_system), intent(inout) :: self
   real(dp), intent(in) :: x, y(:)
   real(dp), intent(out) :: dydx(:)
   integer, intent(out) :: ierr
   ierr = 0
   dydx(1) = -self%k*y(1)
end subroutine
```

and call `bsstep(system, y, dydx, num_eqs, x, h_try, tolerance, y_scale, h_did, h_next, step_err, max_stage_index, extrap_order, underflow_reference, ierr)` in your own stepping loop. On success x has advanced by `h_did`, and `h_next` is the suggested next step. The limits: 1 ≤ `num_eqs` ≤ `bsstep_max_equations` (15), 2 ≤ `extrap_order` ≤ `bsstep_max_extrap_order` (7) and `extrap_order` ≤ `max_stage_index` ≤ `bsstep_max_stages` (11); outside them bsstep returns `numerics_bad_input` without a step. The envelope integration's two, the controls `max_stage_index` and `extrap_order` (legacy `IMAX`, `NUSE`), are checked against these limits when the deck is read (`io/controls_check.f90`, from the same parameters), so a deck outside them stops at startup. A failure of the derivatives returns `numerics_derivatives_failed` (their code in the optional `derivatives_ierr`), and a step that vanishes (`underflow_reference + h == underflow_reference`) `numerics_step_underflow`; on either, `y` and `x` are as they came in (`h_did`, `h_next` and the error estimate are not meaningful). `underflow_reference` sets the scale against which h vanishes: the envelope passes its hydrogen fraction, as the F77 code did. A reference of zero makes the test h == 0, which a failing step reaches only after h has shrunk by 16 per retry to the smallest subnormal, about 270 retries (some 84 000 derivative calls with 11 stages); pass a value of the order of your independent variable.

The step's error is the largest scaled error among the components that are not NaN; it is NaN only when every component is. That reproduces how the original compiled loop behaved, and `test_numerics` pins it.

## 5. Removed names

`ysplin` and `intpol` (2026-10: wrappers of `natural_spline` with `compute_spline_coefficients`, and with `splint(..., dy=)`), `intrp2` (`interp` without `dweight`) and `numerics_termination` (the envelope's status is `envelope_termination`, `core/envelope/envelope_integration_kernel`) went once their callers had moved, bit for bit.

## 6. Tests

`test_numerics` (`make tests` in `src/`, then `pytest test_numerics.py` from the repository root) has two parts. The pins print what every public routine returns on small polynomial and rational data, error paths included, for a byte comparison with `numerics/test/expected_test_numerics.out`; they guard refactors. After a deliberate change, regenerate the file with `src/test_numerics > src/numerics/test/expected_test_numerics.out` (from the repository root, after `make tests`) and read the diff. The checks test behaviour, each printing its error and `ok` or `FAIL` (the program stops with status 1 on a FAIL; the pytest compares a check's line by label and status):
- **cubic splines:** the clamped spline is exact on a cubic, value and `dy`; `compute_spline_coefficients`' power form is `splint`'s, its last column too; natural ends have zero second derivatives and a line is reproduced, extrapolation included; the clamped spline's error falls as h⁴ (its slope's as h³) on 1/(1+x²); the natural spline through two points is their line, and through (0,0), (1,1), (3,0) has the slope 1/2 at the middle; `cspline`'s codes for n < 2, a repeated and a non-monotone abscissa, and a decreasing table giving the increasing one's second derivatives;
- **osplin:** monotone data give a monotone spline and convex data a convex one; it interpolates the nodes and reproduces a line; every nonempty sorted subset of 16 candidate points (below the table, on and between nodes, node n−1 and node n twice, above: 65535 runs) equals its points one at a time, bit for bit; n < 3 and decreasing points give their codes and leave the output untouched, and so does a point beyond an end interval of zero width;
- **Lagrange:** `inter3` exact on a quadratic and `interp`/`lagrange4` on a cubic, with their derivative weights; `interp`'s weights sum to 1 and their derivatives to 0, and they are the same bit for bit without `dweight`; `quint` exact on a quadratic; `polint` exact on a quartic, free of scale (abscissae of 1e-25), and coincident abscissae give their code; the descriptive aliases give their legacy routines' results;
- **lir:** exact on a cubic and a line, a decreasing table the same as the increasing one (on a rational function, where the stencil matters), linear below four points, the nodes' values on a node, outside the table flag 0, and its error codes;
- **searches:** `locate` brackets x in increasing and decreasing tables, gives 1 and n−1 on the end nodes, the end interval within 1% of the span beyond an end and 0 or n beyond that; `findex` and `hunt` bracket x from any starting guess, nodes included, with their codes outside;
- **linear algebra:** an LU solve's residual; rank-deficient matrices (a round-off pivot, a zero pivot) and a zero row (the matrix untouched) are singular, a NaN column is bad input; `simeqc` gives LU's solution for a matrix that needs interchanges, flags a singular matrix and a round-off pivot, solves a matrix with a row scaled by 1e-30, and a NaN column is bad input; `solve_tridiagonal`, and its two halves, give a dense solve's solution, and a zero first or later pivot is singular;
- **bsstep:** sin x, cos x and exp(−0.3x) to x = 5 within 1e-9, x advancing by `h_did`; limits out of range, failed derivatives and a vanishing step give their codes, and leave y and x as they came in (the pins show each branch of a step's acceptance: the sixth stage, the seventh, a later one, and after a retry);
- **status:** every code has its own text.

Each check was confirmed to fail under a mutation of the code it covers: `python3 src/numerics/test/mutate_numerics.py [--only REGEX]` (from the repository root, after `make tests`) applies each mutation, rebuilds `test_numerics`, expects the check to print FAIL and restores the file whatever happens (71 mutations, all caught, on 2026-10-29). `--check` only verifies that each mutation's text is still in its file; `test_numerics.py` runs it, so a rewrite of a mutated routine that strands a mutation fails the unit tests.

`python3 src/tools/coverage.py --work DIR --dir numerics --test test_numerics --deck examples/run_config_matrix/cm_base.inlist` measures the line coverage of numerics/ (a -O0 --coverage build in DIR): 90.6% of 790 lines (2026-10-29, after the round-2 changes). The other files are at 98.7–100%, and `quadratic_spline` at 71.7%, whose unexecuted lines are the branches that the slopes from `slopes` never reach (the module header says why): knot cases 3 and 4 in `choose`, `cases` and `spline`, an end slope opposing the secant, a zero secant, and the chord for a knot on an interval's end; and `compute_spline_coefficients`' one-point column. math/ is at 100% under `test_math` and `test_math_crmath` (src/math/README.md §5).

**Add an interpolation scheme** (traced 2026-10-29 by adding a module `numerics/two_point_weights.f90` with a `linear_weights` and following these steps on a throwaway branch: the build, the checker, the check, its mutation and the pytest all behaved as written):
1. Write it as a `pure` module procedure with explicit arguments (assumed-shape arrays, `real(dp)`) in the module of its kind: `table_interpolation.f90` for weights, searches and table lookups, `cubic_spline.f90` or `quadratic_spline.f90` for splines, or a new module `numerics/<noun>.f90` named after its module (the Makefile finds it; run `python3 tools/gen_makefile_deps.py` in `src/` after adding a module or a `use`, and commit `deps.mk`, which `test_generated_files.py` checks). If it can fail, give it an `ierr` from `numerics_status` (`use numerics_status, only: ...`), check its cheap preconditions, and print nothing (the checker refuses I/O and `stop` in numerics/; a new module is internal to numerics/ at once, so other directories reach it only through `numerics_lib`).
2. Add it to its module's `public` list and to `numerics_lib.f90`'s `use <module>, only:` and `public` lines.
3. In `numerics/test/test_numerics.f90`, add it to the program's `use numerics_lib, only:` list and add a check that a mistake would break: exactness on the polynomials of its degree, its measured order on a smooth function, its error codes. Add a mutation of it to `numerics/test/mutate_numerics.py` (a tuple: the file under `src/`, a text that occurs once, its replacement, the check's label) and confirm that the check fails (`python3 src/numerics/test/mutate_numerics.py --only '<label>'`). Regenerate the expected output: `make -C src tests && src/test_numerics > src/numerics/test/expected_test_numerics.out` from the repository root, and look at the diff (only your lines should be new; the pins are compared byte for byte, the checks' lines by label and status).
4. Add a row to §2's table and the routine to §1's, and to the numerics table of `docs/HOW_YREC_WORKS.md` §5 (routine, what it computes, its callers).
5. Callers take it with `use numerics_lib, only: <name>` (the checker refuses a bare `use numerics_lib`) and report its failures as `<routine>: <numerics_status_text(ierr)>` through their domain's messages (§3).

**Add a solver** (a linear or ODE solver): the same steps, in `linear_solvers.f90` or `bulirsch_stoer.f90` (or a new module), with a check against a known solution (a residual, a dense solve, an analytic ODE solution) and one per error path. A failure the existing codes do not name gets a new code in `numerics_status.f90`, with its text in `numerics_status_text`, its export in `numerics_lib.f90`, and its place in the `codes` list of test_numerics' `status_checks` (which checks that codes and texts are distinct; the loop takes the list's size). A solver whose interface carries a domain's physics belongs with that domain (§1).

Several routines are adapted from *Numerical Recipes* (Press et al.): `cspline`, `splint`, `polint`, `ludcmp`, `lubksb`, `hunt`, `locate`, `bsstep` (with `mmid` and `ratext`). `docs/MODERNIZATION_CHANGES.md` flags this provenance as a licensing liability; replacing them changes results at round-off and needs a rebaseline.

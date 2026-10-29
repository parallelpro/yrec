# numerics/ — the numerical library

This directory holds the general-purpose numerical methods the physics uses: splines, interpolation in tables, searches, linear solvers and an ODE integrator. Code outside numerics/ uses it only through `numerics_lib`, and names what it takes (`use numerics_lib, only: ...`); `tools/check_boundaries.py` enforces both. numerics reads nothing from `star` (`STAR_FREE_DIRS`), keeps no state from one call to the next, does no I/O, and knows no physics: a routine that needs a caller's physics, such as an integrand, receives it as an object. It uses no other YREC module (not even `math_lib`, since no routine here calls an elementary function).

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `numerics_lib` (`numerics_lib.f90`) | **The public interface.** It re-exports what the modules below make public, the status codes and the bsstep limits. |
| `numerics_status` | The status codes every routine returns in `ierr` (§3), and `numerics_status_text`, a short description of each for a caller's message. |
| `cubic_spline` | `cspline` and `splint` (Numerical Recipes): the second derivatives of the cubic spline for given end slopes, and its value at x, with its first derivative when asked (`dy=`). `natural_spline`: `cspline` with both ends natural. `compute_spline_coefficients`: the spline's cubic on each interval in powers of x − x(i), for tables that store it (OPAL92 and Kurucz 1990 opacities). `ysplin`, `intpol`: kept for now (§5). |
| `quadratic_spline` | `osplin`: the shape-preserving osculatory quadratic spline of ACM Algorithm 574 (McAllister & Roulier 1981), which keeps the data's monotonicity and convexity. Its steps (`slopes`, `meval`, `choose`, `cases`, `search`, `spline`, `safedivide`) are private. |
| `table_interpolation` | Lagrange weights for 3 and 4 points with their derivatives (`inter3`, `interp`; `intrp2` kept for now), applied by `lagrange4`; `quint` for three equidistant points; Neville's `polint` for any number; `lir`, cubic interpolation of several tabulated functions at once; and the searches `locate` (bisection), `hunt` (warm-started, Numerical Recipes) and `findex` (warm-started, YCK). |
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
| interpolate with 4 points (or 3) and want the derivative | `interp` (`inter3`), then `lagrange4` | exact on cubics (quadratics); the weights' derivatives with `dweight` |
| interpolate at three equidistant points | `quint` | exact on quadratics |
| interpolate (or extrapolate to zero) with n points and an error estimate | `polint` | Neville; exact on polynomials of degree n − 1 |
| interpolate many functions tabulated on the same points | `lir` | cubic (4 points or more) or linear |
| find the interval holding x | `locate` (bisection), `hunt` or `findex` (start from a guess: for a sequence of nearby x) | the conventions differ at the nodes and outside: see each header |
| solve a dense linear system | `ludcmp` then `lubksb` (several right-hand sides: one `lubksb` each) | `simeqc` solves in place for a few right-hand sides at once |
| solve a tridiagonal system | `solve_tridiagonal` | `tridiagonal_forward_sweep` + `tridiagonal_back_substitution` when a step goes between the two |
| integrate ODEs | `bsstep` in your stepping loop, on an `ode_system` (§4) | smooth right-hand sides, tight tolerances |

## 3. Conventions

- Real kind `dp` (`real64`). Arrays are sized by the arguments (explicit shape with a size argument, as the Numerical Recipes routines have it, or assumed shape); nothing is dimensioned by the model's maximum number of zones.
- **The error contract.** A routine that can fail has an `ierr` argument (optional where it was added to an existing routine, `osplin` and `lir`) and returns one of `numerics_status`'s codes: `numerics_ok` (0), `numerics_singular` (a zero row or pivot), `numerics_coincident_points` (two abscissae equal where they must differ), `numerics_bad_input` (sizes or limits out of range, too few points), `numerics_step_underflow` (bsstep), `numerics_unsorted_points` (osplin), `numerics_derivatives_failed` (bsstep; the derivatives' own code in `derivatives_ierr`). The cheap preconditions are checked: sizes, bsstep's limits, the first pivot. A routine never stops the run.
- **No I/O.** numerics prints nothing (`tools/check_boundaries.py` checks that numerics/ and math/ hold no I/O statement): the caller reports a failure, where the failure means something. The callers that report keep the messages the routines used to print (`ERROR IN SPLINT ROUTINE.` to the run log; `ludcmp: singular matrix`, `solve_tridiagonal: singular matrix`, `polint: repeated abscissa (den ~ 0)` and `ERROR IN BSSTEP` to the terminal); `numerics_status_text` gives a short description of each code for a new caller's message.
- **Searches.** `locate`, `hunt` and `findex`, and the bisections inside `splint` and `osplin`, differ at the nodes and outside the table; each header states its convention.
- Procedures without side effects are `pure`. `bsstep` and `mmid` are not (they call an `ode_system`'s derivatives, which may change their object); `splint`, `intpol`, `polint`, `ludcmp` and `solve_tridiagonal` are not yet, since they printed until 2026-10-29.
- Numerical Recipes and IBM SSP routines keep their authors' names (`cspline`, `splint`, `polint`, `ludcmp`, `lubksb`, `hunt`, `bsstep` with `mmid` and `ratext`; `simeqc`); routines written for YREC are named for what they do.

## 4. Integrating ODEs

`bsstep` integrates any system you describe as an extension of the abstract type `ode_system`. Its deferred procedure `derivatives(self, x, y, dydx, ierr)` returns dy/dx. The extension carries whatever the derivatives need. The envelope integration's `envelope_equations` (`core/envelope/envelope_integration_kernel.f90`) holds the luminosity, rotation factors, composition and EOS state, and calls `compute_atmosphere_derivs` or `compute_envelope_derivs`. To integrate a new system:

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

and call `bsstep(system, y, dydx, num_eqs, x, h_try, tolerance, y_scale, h_did, h_next, step_err, max_stage_index, extrap_order, underflow_reference, ierr)` in your own stepping loop. On success x has advanced by `h_did`, and `h_next` is the suggested next step. The limits: 1 ≤ `num_eqs` ≤ `bsstep_max_equations` (15), 2 ≤ `extrap_order` ≤ `bsstep_max_extrap_order` (7) and `extrap_order` ≤ `max_stage_index` ≤ `bsstep_max_stages` (11); outside them bsstep returns `numerics_bad_input` without a step. The envelope integration's two, the controls `max_stage_index` and `extrap_order` (legacy `IMAX`, `NUSE`), are checked against these limits when the deck is read (`io/controls_check.f90`, from the same parameters), so a deck outside them stops at startup. A failure of the derivatives returns `numerics_derivatives_failed` (their code in the optional `derivatives_ierr`), and a step that vanishes (`underflow_reference + h == underflow_reference`) `numerics_step_underflow`.

The step's error is the largest scaled error among the components that are not NaN; it is NaN only when every component is. That reproduces how the original compiled loop behaved, and `test_numerics` pins it.

## 5. Kept for now

These public names stay until the callers outside numerics/ move off them (no rename or removal in the 2026-10-29 wave); new code uses the replacement.

| Name | Replacement | Callers left |
|---|---|---|
| `ysplin` | `natural_spline` + `compute_spline_coefficients` (ysplin is now that) | none (test only) |
| `intpol` | `natural_spline` + `splint(..., dy=)` (intpol is now that) | none (test only) |
| `intrp2` | `interp` without `dweight` (intrp2 calls it) | kap (alex94, alex06, opal95), element_diffusion (settling_grid, microdiff_mapping) |
| `numerics_termination` | `envelope_termination` (`core/envelope/envelope_integration_kernel`), the envelope integration's status | none |

## 6. Tests

`test_numerics` (`make tests` in `src/`, then `pytest test_numerics.py` from the repository root) has two parts. The pins print what every public routine returns on small polynomial and rational data, error paths included, for a byte comparison with `numerics/test/expected_test_numerics.out`; they guard refactors. The checks test behaviour, each printing its error and `ok` or `FAIL` (the program stops with status 1 on a FAIL; the pytest compares a check's line by label and status):
- **cubic splines:** the clamped spline is exact on a cubic, value and `dy`; `compute_spline_coefficients`' power form is `splint`'s; natural ends have zero second derivatives and a line is reproduced, extrapolation included; the clamped spline's error falls as h⁴ (its slope's as h³) on 1/(1+x²); `ysplin`/`intpol` are the natural spline for n = 2, 3, 4 and 8, with the slope 1/2 at the middle of (0,0), (1,1), (3,0);
- **osplin:** monotone data give a monotone spline and convex data a convex one; it interpolates the nodes and reproduces a line; a run of points equals its points one at a time (below the table, on nodes, on the last node, above, an interior point then points above); n < 3 and decreasing points give their codes and leave the output untouched;
- **Lagrange:** `inter3` exact on a quadratic and `interp`/`lagrange4` on a cubic, with their derivative weights; `interp`'s weights sum to 1 and their derivatives to 0, and `intrp2`'s are `interp`'s bit for bit; `quint` exact on a quadratic; `polint` exact on a quartic, free of scale (abscissae of 1e-25), and coincident abscissae give their code;
- **lir:** exact on a cubic and a line, a decreasing table the same as the increasing one (on a rational function, where the stencil matters), linear below four points, the nodes' values on a node, outside the table flag 0, and its error codes;
- **searches:** `locate` brackets x in increasing and decreasing tables, gives 1 and n−1 on the end nodes, the end interval within 1% of the span beyond an end and 0 or n beyond that; `findex` and `hunt` bracket x from any starting guess, nodes included, with their codes outside;
- **linear algebra:** an LU solve's residual; rank-deficient matrices (a round-off pivot, a zero pivot) and a zero row (the matrix untouched) are singular, a NaN column is bad input; `simeqc` gives LU's solution and flags a singular matrix; `solve_tridiagonal`, and its two halves, give a dense solve's solution, and a zero first or later pivot is singular;
- **bsstep:** sin x, cos x and exp(−0.3x) to x = 5 within 1e-9, x advancing by `h_did`; limits out of range, failed derivatives and a vanishing step give their codes;
- **status:** every code has its own text.

Each check was confirmed to fail under a mutation of the code it covers: `python3 src/numerics/test/mutate_numerics.py [--only REGEX]` (from the repository root, after `make tests`) applies each mutation, rebuilds `test_numerics`, expects the check to print FAIL and restores the file whatever happens (55 mutations, all caught, on 2026-10-29).

`python3 src/tools/coverage.py --work DIR --dir numerics --test test_numerics --deck examples/run_config_matrix/cm_base.inlist` measures the line coverage of numerics/ (a -O0 --coverage build in DIR): 91.1% of 827 lines (2026-10-29). Every file is at 99–100% but `quadratic_spline` (74.5%), whose unexecuted lines are the branches that the slopes from `slopes` never reach (the module header says why): knot cases 3 and 4 in `choose`, `cases` and `spline`, an end slope opposing the secant, a zero secant, and the chord for a knot on an interval's end; and `compute_spline_coefficients`' one-point column. math/ is at 100% under `test_math` (src/math/README.md).

**Add an interpolation scheme** (traced 2026-10-29):
1. Write it as a `pure` module procedure with explicit arguments (assumed-shape arrays, `real(dp)`) in the module of its kind: `table_interpolation.f90` for weights, searches and table lookups, `cubic_spline.f90` or `quadratic_spline.f90` for splines, or a new module `numerics/<noun>.f90` named after its module (the Makefile finds it; run `python3 tools/gen_makefile_deps.py` in `src/` after adding a module or a `use`). If it can fail, give it an `ierr` from `numerics_status` (`use numerics_status, only: ...`), check its cheap preconditions, and print nothing.
2. Add it to its module's `public` list and to `numerics_lib.f90`'s `use <module>, only:` and `public` lines.
3. In `numerics/test/test_numerics.f90`, add it to the program's `use numerics_lib, only:` list and add a check that a mistake would break: exactness on the polynomials of its degree, its measured order on a smooth function, its error codes. Add a mutation of it to `numerics/test/mutate_numerics.py` and confirm that the check fails (`python3 src/numerics/test/mutate_numerics.py --only '<label>'`). Regenerate `expected_test_numerics.out` (the checks' lines are compared by label and status).
4. Add a row to §2's table and the routine to §1's.
5. Callers take it with `use numerics_lib, only: <name>` (the checker refuses a bare `use numerics_lib`).

**Add a solver** (a linear or ODE solver): the same steps, in `linear_solvers.f90` or `bulirsch_stoer.f90` (or a new module), with a check against a known solution (a residual, a dense solve, an analytic ODE solution) and one per error path. A failure the existing codes do not name gets a new code in `numerics_status.f90`, with its text in `numerics_status_text`, its export in `numerics_lib.f90`, and its place in test_numerics' `status_checks` list (which checks that codes and texts are distinct). A solver whose interface carries a domain's physics belongs with that domain (§1).

Several routines are adapted from *Numerical Recipes* (Press et al.): `cspline`, `splint`, `polint`, `ludcmp`, `lubksb`, `hunt`, `locate`, `bsstep` (with `mmid` and `ratext`). `docs/MODERNIZATION_CHANGES.md` flags this provenance as a licensing liability; replacing them changes results at round-off and needs a rebaseline.

# numerics/ — the numerical library

This directory holds the general-purpose numerical methods the physics uses: splines, interpolation in tables, linear solvers and an ODE integrator. Code outside numerics/ uses it only through `numerics_lib` (enforced by `tools/check_boundaries.py`). numerics reads nothing from `star`, keeps no state from one call to the next and knows no physics. A routine that needs a caller's physics, such as an integrand, receives it as an object.

---

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `numerics_lib` (`numerics_lib.f90`) | **The public interface.** It re-exports everything below, plus `numerics_termination`, the status a driver gives a failed envelope integration. |
| `cubic_spline` | `cspline` and `splint` (Numerical Recipes): the second derivatives for given end slopes (a slope of 10³⁰ or more means a natural end), and the value at x. `natural_spline` calls `cspline` with both ends natural. `ysplin` and `intpol` are YCK's natural spline in coefficient form, used for the Cox and Kurucz opacity tables. |
| `quadratic_spline` | `osplin`: the shape-preserving osculatory quadratic spline of ACM Algorithm 574 (McAllister & Roulier 1981), which keeps the data's monotonicity and convexity. Its helpers (`slopes`, `meval`, `choose`, `cases`, `search`, `spline`) are private. |
| `table_interpolation` | Lagrange weights for 3 and 4 points (`inter3`, `interp`, `intrp2`), with `lagrange4` to apply them; `quint` for equidistant points; Neville's `polint`; `lir`, cubic interpolation of several tabulated functions at once; and the bracketing searches `findex` and `locate`. |
| `linear_solvers` | `ludcmp`/`lubksb` (LU with partial pivoting), `simeqc` (Gauss-Jordan for one or more right-hand sides) and `solve_tridiagonal` (the Thomas algorithm). |
| `bulirsch_stoer` | `bsstep`: one adaptive Bulirsch-Stoer step (modified midpoint rule with rational extrapolation) of an `ode_system`. |

Solvers that carry a domain's physics in their interface live with that domain:
- the shape-factor quadratures are in `rotation/shape/shell_quadrature.f90`;
- the angular-momentum solve is `solve_angular_momentum_transport` in `rotation/seculr/secular_transport.f90`;
- the MHD table lookup is `interpolate_mhd_table` in `eos/mhd/mhd_eos.f90`.

## 2. Integrating ODEs

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

and call `bsstep(system, y, dydx, num_eqs, x, h_try, tolerance, y_scale, h_did, h_next, step_err, max_stage_index, extrap_order, underflow_reference, ierr)` in your own stepping loop. On success x has advanced by `h_did`, and `h_next` is the suggested next step. ierr is the derivatives' failure, or 1 when the step underflows (`underflow_reference + h == underflow_reference`).

The step's error is the largest scaled error among the components that are not NaN; it is NaN only when every component is. That reproduces how the original compiled loop behaved, and `test_numerics` pins it.

## 3. Conventions

- Real kind `dp` (double precision). Arrays are sized by the arguments; nothing is dimensioned by the model's maximum number of zones.
- Failure is reported through `ierr` (nonzero), never by stopping; a short diagnostic goes to standard output or the run log. A singular matrix or coincident abscissae give ierr = 1.
- Procedures without I/O are `pure`.

## 4. Tests and extending

`test_numerics` (`make tests` in `src/`, then `pytest test_numerics.py` from the repository root) calls every public routine on small polynomial and rational data, error paths included. It drives the quadratic spline through extrapolation, single points, nodes and unsorted points, and takes `bsstep` through each acceptance branch, a retry, an integrand failure, step underflow and NaN errors. Under gcov the only lines it leaves unexecuted are unreachable ones (see the `quadratic_spline` header). A deliberate change to its output means regenerating `numerics/test/expected_test_numerics.out`.

**Add a routine:** put it in the module of its topic (or a new module), export it from `numerics_lib`, and add a section to `test_numerics`. A routine with a single caller that knows that caller's physics belongs with that caller instead.

Several routines are adapted from *Numerical Recipes* (Press et al.): `cspline`, `splint`, `polint`, `ludcmp`, `lubksb`, `bsstep` (with `mmid` and `ratext`). `docs/MODERNIZATION_CHANGES.md` flags this provenance as a licensing liability; replacing them changes results at round-off and needs a rebaseline.

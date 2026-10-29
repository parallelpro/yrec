#!/usr/bin/env python3
"""Mutation check of test_numerics: each check must fail under a mutation
of the code it covers.

    python3 src/numerics/test/mutate_numerics.py [--only REGEX] [--check]

from the repository root, after `make tests` in src/. For each mutation
it edits one source file, rebuilds test_numerics, runs it and expects the
named check to print FAIL (a crash, a nonzero exit or a run longer than
120 s also counts as caught); the file is restored afterwards, whatever
happens. make uses
1 s mtimes, so the script sleeps a second around each write. Exit status
1 if any mutation survives. --check only verifies that every mutation's
text occurs exactly once in its file (run it after editing a mutated
routine: a rewrite can leave a mutation pointing at text that is gone).
"""
import argparse
import pathlib
import re
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[3]
SRC = REPO / "src"
BINARY = SRC / "test_numerics"

# (file under src/, the text replaced (it must occur exactly once), the
# replacement, the label of the check that must fail). Where a routine
# was rewritten, old and new are lists of the forms the code has taken,
# and the one present is used.
MUTATIONS = [
    ('numerics/cubic_spline.f90', '*((y(2)-y(1))/(x(2)-x(1))-yp1)', '*((y(2)-y(1))/(x(2)-x(1))+yp1)', 'cspline/splint: clamped spline exact on a cubic'),
    ('numerics/cubic_spline.f90', '- ((3.0d0*a**2-1.0d0)*y2a(lo)', '- ((3.0d0*a**2-2.0d0)*y2a(lo)', 'splint: its dy exact on the clamped cubic'),
    ('numerics/cubic_spline.f90', 'c(4,i) = (y2(i+1) - y2(i))/(6.0d0*h)', 'c(4,i) = (y2(i+1) - y2(i))/(3.0d0*h)', "compute_spline_coefficients: the power form is splint's"),
    ('numerics/cubic_spline.f90', '         y2(1) = 0.0d0\n         u(1) = 0.0d0', '         y2(1) = 0.0d0\n         u(1) = 1.0d0', 'natural_spline/splint: exact on a line'),
    ('numerics/cubic_spline.f90', '      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0)', '      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0d0) + 1.0d-3', 'natural_spline: zero second derivatives at the ends'),
    ('numerics/cubic_spline.f90', '         u(i) = (6.0d0*((y(i+1)-y(i))', '         u(i) = (6.1d0*((y(i+1)-y(i))', 'cspline: the error falls as h^4'),
    ('numerics/cubic_spline.f90', '+ h*(y2(n-1) + 2.0d0*y2(n))/6.0d0', '+ h*(y2(n-1) + 2.0d0*y2(n))/5.0d0', "compute_spline_coefficients: the power form is splint's, column n too"),
    ('numerics/cubic_spline.f90', '         u(1) = 0.0d0\n      else', '         u(1) = 1.0d-3\n      else', 'natural_spline: slope 1/2 at the middle'),
    ('numerics/cubic_spline.f90', '         qn = 0.0d0\n         un = 0.0d0', '         qn = 0.0d0\n         un = 1.0d-3', 'natural_spline: two points give their line'),
    ('numerics/cubic_spline.f90', '         if (x(i+1) == x(i)) then', '         if (.false.) then', 'cspline: n < 2, repeated and non-monotone'),
    ('numerics/cubic_spline.f90', '         if ((x(i+1) > x(i)) .neqv. (x(2) > x(1))) then', '         if (.false.) then', 'cspline: n < 2, repeated and non-monotone'),
    ('numerics/cubic_spline.f90', '      if (n < 2) then\n         if (present(ierr)) ierr = numerics_bad_input', '      if (n < 1) then\n         if (present(ierr)) ierr = numerics_bad_input', 'cspline: n < 2, repeated and non-monotone'),
    ('numerics/cubic_spline.f90', '         p = sig*y2(i-1)+2.0d0', '         p = sig*y2(i-1)+2.0d0+merge(1.0d-9, 0.0d0, x(2) < x(1))', "natural_spline: a decreasing table"),
    ('numerics/quadratic_spline.f90', ['         first_derivs(idx)=y_diff1/(table_x(idx) - x_hat)', '            first_derivs(idx) = y_diff1/(table_x(idx) - x_hat)'], ['         first_derivs(idx)=-y_diff1/(table_x(idx) - x_hat)', '            first_derivs(idx) = -y_diff1/(table_x(idx) - x_hat)'], 'osplin: monotone data give a monotone spline'),
    ('numerics/quadratic_spline.f90', ['         first_derivs(idx)=y_diff1/(table_x(idx) - x_hat)', '            first_derivs(idx) = y_diff1/(table_x(idx) - x_hat)'], ['         first_derivs(idx)=0.5d0*y_diff1/(table_x(idx) - x_hat)', '            first_derivs(idx) = 0.5d0*y_diff1/(table_x(idx) - x_hat)'], 'osplin: convex data give a convex spline'),
    ('numerics/quadratic_spline.f90', ['         eval_y(start_idx) = table_y(table_idx)\n      start_idx1', '               eval_y(start_idx) = table_y(table_idx)\n', '            eval_y(i) = table_y(j+1)\n'], ['         eval_y(start_idx) = table_y(table_idx) + 1.0d-12\n      start_idx1', '               eval_y(start_idx) = table_y(table_idx) + 1.0d-12\n', '            eval_y(i) = table_y(j+1) + 1.0d-12\n'], 'osplin: interpolates the nodes'),
    ('numerics/quadratic_spline.f90', ['         x_hat=(table_x(prev_idx) + x_bar)/2.d0', '            x_hat = (table_x(prev_idx) + x_bar)/2.d0'], ['         x_hat=(table_x(prev_idx) + x_bar)/2.2d0', '            x_hat = (table_x(prev_idx) + x_bar)/2.2d0'], 'osplin: exact on a line'),
    ('numerics/quadratic_spline.f90', '         call take_arcs(arcs, n - 1)\n', '         if (arcs%interval == 0) call take_arcs(arcs, n - 1)\n', 'osplin: every sorted subset'),
    ('numerics/quadratic_spline.f90', '         if (arcs%interval /= j) &', '         if (arcs%interval == 0) &', 'osplin: every sorted subset'),
    ('numerics/quadratic_spline.f90', '      if (n < 3) then', '      if (n < 2) then', 'osplin: n < 3 is numerics_bad_input'),
    ('numerics/quadratic_spline.f90', ['if (err == 2 .and. present(ierr))', '         if (xval(i+1) >= xval(i)) cycle'], ['if (err == 3 .and. present(ierr))', '         if (.true.) cycle'], 'osplin: decreasing points are numerics_unsorted_points'),
    ('numerics/quadratic_spline.f90', '          (xval(k) > xtab(n) .and. xtab(n) == xtab(n-1))) then', '          (xval(k) > xtab(n) .and. .false.)) then', 'osplin: no extrapolation from a zero-width end interval'),
    ('numerics/quadratic_spline.f90', '      if ((xval(1) < xtab(1) .and. xtab(2) == xtab(1)) .or. &', '      if ((xval(1) < xtab(1) .and. .false.) .or. &', 'osplin: no extrapolation from a zero-width end interval'),
    ('numerics/table_interpolation.f90', '      denom2 = -diff21*diff32\n', '      denom2 = diff21*diff32\n', 'inter3: exact on a quadratic'),
    ('numerics/table_interpolation.f90', '      dweight(3) = (dx1 + dx2)/denom3', '      dweight(3) = (dx1 - dx2)/denom3', 'inter3: its derivative weights exact'),
    ('numerics/table_interpolation.f90', ['      weight(4) = (dx1*dx2*dx3)/denom4\n      dweight(1)', '      weight(4) = (dx1*dx2*dx3)/denom4\n      if (present(dweight))'], ['      weight(4) = (dx1*dx2*dx4)/denom4\n      dweight(1)', '      weight(4) = (dx1*dx2*dx4)/denom4\n      if (present(dweight))'], 'interp/lagrange4: exact on a cubic'),
    ('numerics/table_interpolation.f90', ['      dweight(1) = (dx3*dx4 + dx2*dx4 + dx2*dx3)/denom1', '         dweight(1) = (dx3*dx4 + dx2*dx4 + dx2*dx3)/denom1'], ['      dweight(1) = (dx3*dx4 + dx2*dx4 - dx2*dx3)/denom1', '         dweight(1) = (dx3*dx4 + dx2*dx4 - dx2*dx3)/denom1'], 'interp: its derivative weights exact on a cubic'),
    ('numerics/table_interpolation.f90', '      if (present(dweight)) then\n', '      if (present(dweight)) then\n         weight(1) = weight(1)*(1.0d0+epsilon(1.0d0))\n', 'interp: the weights without dweight are those with it'),
    ('numerics/table_interpolation.f90', '      s = w(1)*y(1)+w(2)*y(2)+w(3)*y(3)+w(4)*y(4)', '      s = w(1)*y(1)+w(2)*y(2)+w(3)*y(3)+w(3)*y(4)', 'interp/lagrange4: exact on a cubic'),
    ('numerics/table_interpolation.f90', '      y = y0 + t*d1 + 0.5d0*t*(t-1.d0)*d2', '      y = y0 + t*d1 + 0.5d0*t*(t+1.d0)*d2', 'quint: exact on a quadratic'),
    ('numerics/table_interpolation.f90', '          den = w/den', '          den = w/(2*den)', 'polint: exact on a quartic'),
    ('numerics/table_interpolation.f90', 'abs(den) <= coincidence_tolerance*max(abs(ho), abs(hp))', 'abs(den) < 1.0d-20', 'polint: scale-free'),
    ('numerics/table_interpolation.f90', ['             ierr = numerics_coincident_points\n             return', '               ierr = numerics_coincident_points\n               return'], ['             ierr = numerics_bad_input\n             return', '               ierr = numerics_bad_input\n               return'], 'polint: coincident abscissae'),
    ('numerics/table_interpolation.f90', ['      weight(closest)=(1.d0+weight(closest))-diff', '         weight(closest) = (1.d0 + weight(closest)) - weight_sum', '         weight(corrected) = (1.d0 + weight(corrected)) - weight_sum'], ['      weight(closest)=(1.d0+weight(closest))-0.999d0*diff', '         weight(closest) = (1.d0 + weight(closest)) - 0.999d0*weight_sum', '         weight(corrected) = (1.d0 + weight(corrected)) - 0.999d0*weight_sum'], 'lir: exact on a cubic and a line'),
    ('numerics/table_interpolation.f90', ['            if (table_z(search_idx) < target_z) exit', '            if (table_z(above) < target_z) exit'], ['            if (table_z(search_idx) > target_z) exit', '            if (table_z(above) > target_z) exit'], 'lir: a decreasing table as the same table increasing'),
    ('numerics/table_interpolation.f90', ['      y2=1.0d0-y1', '      frac_upper = 1.0d0 - frac_lower'], ['      y2=1.001d0-y1', '      frac_upper = 1.001d0 - frac_lower'], 'lir: linear below four points'),
    ('numerics/table_interpolation.f90', ['               result_y(y_idx)=table_y(y_idx+base_idx)\n', '               result_y(i) = table_y(i, above)\n'], ['               result_y(y_idx)=table_y(y_idx+base_idx) + 1.0d-9\n', '               result_y(i) = table_y(i, above) + 1.0d-9\n'], 'lir: on a node, the tabulated values'),
    ('numerics/table_interpolation.f90', ['            interp_flag=0\n            exit', '            interp_flag = 0\n            exit'], ['            interp_flag=1\n            exit', '            interp_flag = 1\n            exit'], 'lir: outside the table, flag 0'),
    ('numerics/table_interpolation.f90', ['      if (num_points < 2) then', '      if (num_points < 2 .or.'], ['      if (num_points < 1) then', '      if (num_points < 1 .or.'], 'lir: one point is numerics_bad_input'),
    ('numerics/table_interpolation.f90', ['         if (diff == 0.0) then', '      if (span == 0.0d0) then'], ['         if (diff == 1.0d300) then', '      if (span == 1.0d300) then'], 'lir: a flat table is numerics_coincident_points'),
    ('numerics/table_interpolation.f90', '         if ((xx(n) > xx(1)).eqv.(x > xx(jm))) then', '         if ((xx(n) > xx(1)).eqv.(x >= xx(jm))) then', 'locate: brackets x in increasing and decreasing tables'),
    ('numerics/table_interpolation.f90', '      if (j == 0 .and. abs(x - xx(1)) <= slack) then\n         j = 1', '      if (j == 0 .and. abs(x - xx(1)) <= slack) then\n         j = 0', 'locate: the end nodes give 1 and n-1'),
    ('numerics/table_interpolation.f90', 'edge_slack = 0.01d0', 'edge_slack = 0.001d0', 'locate: within 1% of the span beyond an end'),
    ('numerics/table_interpolation.f90', 'edge_slack = 0.01d0', 'edge_slack = 0.1d0', 'locate: beyond the slack, 0 and n'),
    ('numerics/table_interpolation.f90', ['            if(grid_x(j+1) > x_eval)then', '            if (grid_x(j+1) > x_eval) then', '      if (index == 0) index = -1'], ['            if(grid_x(j+1) >= x_eval)then', '            if (grid_x(j+1) >= x_eval) then', '      if (index == 0) index = 1'], 'findex: brackets x from any guess'),
    ('numerics/table_interpolation.f90', '      if (index == n_grid) index = -n_grid', '      if (index == n_grid) index = n_grid - 1', 'findex: brackets x from any guess'),
    ('numerics/table_interpolation.f90', ['      if (target_value >= table(index_mid).eqv.ascending) then', '         if ((x >= table(index_mid)) .eqv. ascending) then'], ['      if (target_value > table(index_mid).eqv.ascending) then', '         if ((x > table(index_mid)) .eqv. ascending) then'], 'hunt: brackets x from any guess'),
    ('numerics/table_interpolation.f90', '         index = 0\n         index_hi = n + 1', '         index = 0\n         index_hi = n', 'hunt: brackets x from any guess'),
    ('numerics/linear_solvers.f90', 'b(i) = sum/a(i,i)', 'b(i) = sum/(a(i,i)*(1.0d0+1.0d-12))', 'ludcmp/lubksb: the residual of a 6x6 solve'),
    ('numerics/linear_solvers.f90', 'if (aamax <= singular_pivot_tolerance) ierr', 'if (aamax <= 0.0d0) ierr', 'ludcmp: a round-off pivot'),
    ('numerics/linear_solvers.f90', 'if (aamax <= singular_pivot_tolerance) ierr', 'if (aamax < 0.0d0) ierr', 'ludcmp: a zero pivot'),
    ('numerics/linear_solvers.f90', 'if (aamax == 0d0) then', 'if (aamax < 0d0) then', 'ludcmp: a zero row is numerics_singular, the matrix untouched'),
    ('numerics/linear_solvers.f90', 'if (ieee_is_nan(a(imax,j))) then', 'if (.false.) then', 'ludcmp: a NaN column is numerics_bad_input'),
    ('numerics/linear_solvers.f90', ['      system_matrix(ib)=system_matrix(ib)-system_matrix(ia)*system_matrix(ic)', '            a(row,num_cols) = a(row,num_cols) - a(row,k)*a(k,num_cols)'], ['      system_matrix(ib)=system_matrix(ib)+system_matrix(ia)*system_matrix(ic)', '            a(row,num_cols) = a(row,num_cols) + a(row,k)*a(k,num_cols)'], "simeqc: LU's solution for two right-hand sides"),
    ('numerics/linear_solvers.f90', ['      if (abs(biga) == 0.0d0) then', '         if (abs(pivot) == 0.0d0) then', '         if (largest_scaled <= singular_pivot_tolerance) then'], ['      if (abs(biga) < 0.0d0) then', '         if (abs(pivot) < 0.0d0) then', '         if (largest_scaled < 0.0d0) then'], 'simeqc: a singular matrix is numerics_singular'),
    ('numerics/linear_solvers.f90', '         if (largest_scaled <= singular_pivot_tolerance) then', '         if (largest_scaled <= 0.0d0) then', 'simeqc: a round-off pivot'),
    ('numerics/linear_solvers.f90', '                    largest_scaled = abs(a(i,j))/row_size(i)', '                    largest_scaled = abs(a(i,j))/maxval(row_size)', 'simeqc: a row scaled by 1e-30'),
    ('numerics/linear_solvers.f90', '         if (ieee_is_nan(pivot)) then', '         if (.false.) then', 'simeqc: a NaN column is numerics_bad_input'),
    ('numerics/linear_solvers.f90', ['         ex(j) = ex(j) - gama(j+1)*ex(j+1)', '         x(j) = x(j) - gama(j+1)*x(j+1)'], ['         ex(j) = ex(j) - gama(j+1)*ex(j+1)*(1.0d0+1.0d-12)', '         x(j) = x(j) - gama(j+1)*x(j+1)*(1.0d0+1.0d-12)'], "solve_tridiagonal: a dense solve's solution"),
    ('numerics/linear_solvers.f90', '         x(j) = (r(j) - a(j)*x(j-1))/bet', '         x(j) = (r(j) + a(j)*x(j-1))/bet', 'tridiagonal sweep + back substitution'),
    ('numerics/linear_solvers.f90', '      bet = b(1)\n      if (bet == 0) then', '      bet = b(1)\n      if (.false.) then', 'solve_tridiagonal: a zero first pivot'),
    ('numerics/linear_solvers.f90', '         bet = b(j) - a(j)*gama(j)\n         if (bet == 0) then', '         bet = b(j) - a(j)*gama(j)\n         if (.false.) then', 'solve_tridiagonal: a zero later pivot'),
    ('numerics/bulirsch_stoer.f90', '       y_out(i) = 0.5d0*(y_mid(i) + y_new(i) + h_sub*y_out(i))', '       y_out(i) = 0.5d0*(y_mid(i) + y_new(i) + 1.01d0*h_sub*y_out(i))', 'bsstep: sin x, cos x and exp(-0.3x)'),
    ('numerics/bulirsch_stoer.f90', '          indep_var = indep_var + h\n', '          indep_var = indep_var + 0.5d0*h\n', 'bsstep: x advances by h_did'),
    ('numerics/bulirsch_stoer.f90', '          max_stage_index < extrap_order .or. ', '', 'bsstep: limits out of range'),
    ('numerics/bulirsch_stoer.f90', 'if (present(derivatives_ierr)) derivatives_ierr = jerr_integrand', 'if (present(derivatives_ierr)) derivatives_ierr = 0', 'bsstep: failed derivatives'),
    ('numerics/bulirsch_stoer.f90', '       ierr = numerics_step_underflow', '       ierr = numerics_bad_input', 'bsstep: zero tolerance ends in numerics_step_underflow'),
    ('numerics/bulirsch_stoer.f90', ' derivatives_ierr = jerr_integrand\n               y = y_sav\n', ' derivatives_ierr = jerr_integrand\n', 'bsstep: on a failure y and x are as they came in'),
    ('numerics/bulirsch_stoer.f90', '            ierr = numerics_step_underflow\n            y = y_sav\n', '            ierr = numerics_step_underflow\n', 'bsstep: on a failure y and x are as they came in'),
    ('numerics/numerics_lib.f90', 'interpolate_tabulated_functions => lir, bracket_index => findex', 'interpolate_tabulated_functions => lir, bracket_index => hunt', 'aliases: each gives'),
    ('numerics/numerics_status.f90', "         text = 'step underflow'", "         text = 'singular matrix'", 'numerics_status: distinct codes'),
]


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        missing = 0
        for path, old, new, label in MUTATIONS:
            text = (SRC / path).read_text()
            olds = old if isinstance(old, list) else [old]
            if not any(text.count(o) == 1 for o in olds):
                print(f"MISSING  {label}  ({path})")
                missing += 1
        print(f"{len(MUTATIONS) - missing} of {len(MUTATIONS)} mutations apply")
        return 1 if missing else 0
    survived = 0
    for path, old, new, label in MUTATIONS:
        if a.only and not re.search(a.only, label):
            continue
        f = SRC / path
        original = f.read_text()
        olds, news = (old, new) if isinstance(old, list) else ([old], [new])
        present = [k for k, o in enumerate(olds) if original.count(o) == 1]
        assert present, (path, olds)
        old, new = olds[present[0]], news[present[0]]
        try:
            time.sleep(1)
            f.write_text(original.replace(old, new))
            b = run(["make", "test_numerics"], cwd=SRC)
            if b.returncode != 0:
                print(f"BUILD FAILED  {label}\n{b.stdout[-800:]}{b.stderr[-800:]}")
                survived += 1
                continue
            try:
                r = run([str(BINARY)], cwd=SRC, timeout=120)
            except subprocess.TimeoutExpired:
                print(f"caught    {label}  ({path}; the mutant did not finish in 120 s)")
                continue
            lines = [l for l in r.stdout.splitlines() if l.startswith(label)]
            caught = any(l.rstrip().endswith("FAIL") for l in lines) or \
                (not lines and r.returncode != 0)
            print(f"{'caught  ' if caught else 'SURVIVED'}  {label}  ({path})")
            survived += 0 if caught else 1
        finally:
            f.write_text(original)
            time.sleep(1)
            f.touch()
    run(["make", "test_numerics"], cwd=SRC)
    print(f"{survived} mutation(s) survived")
    return 1 if survived else 0


if __name__ == "__main__":
    sys.exit(main())

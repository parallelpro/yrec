# numerics/ + math/ review (2026-10-29, read-only evaluation; paths under src/)

Scores vs rotation/: numerics Readability 3, Organization 4, Maintainability 2.5, Extensibility 3, Docs 3;
math 4 / 3.5 / 3 / 3.5 / 2.

Good skeleton: numerics_lib facade enforced (FACADE_ONLY_DOMAINS), star-free (STAR_FREE_DIRS), no saved
state, Bulirsch-Stoer takes an ode_system. Weak edges: unchecked preconditions (OOB or wrong results in
osplin, ysplin, bsstep, ludcmp, solve_tridiagonal off the happy path; one reachable by user controls
(bsstep's imax/nuse), one via OPAL92 at the top T node); ierr that lies (README says singular systems give
ierr = 1; several give 0); printing kernels (three destinations). Duplicated algorithms inside: two natural
cubic splines, three 4-point Lagrange weight formulas, six interval searches. The test byte-pins outputs,
including two locate bugs. Outside numerics: four exact Lagrange-4 copies (class A to remove) and one
Numerical Recipes routine (hunt) missing from ROADMAP's licensing list. math_lib is well motivated and
should stay separate; it needs a README and a backend test; safe_divide_by_exp has no caller.
The reviewer confirmed each bug by reading or by compiling a copy with -fcheck=bounds.

## Possible bugs
1. bulirsch_stoer.f90:40, 82, 98, 211, 227 -- none of bsstep's limits is checked. VERIFIED by the
   integrator: num_eqs <= 15 (nmax), max_stage_index <= 11, extrap_order <= 7; imax and nuse are inlist
   controls (namelist /physics/, io/legacy_controls.f90:280/:288, defaults :409-410; copied in
   io/read_controls.f90:253-254 into state/integration_parameters). imax = 12 or nuse = 8 writes past
   substep_sequence, x_hist, fx, tableau; >= 16 equations silently writes into the next tableau column;
   nuse = 1 reads substep_sequence(0) (line 129).
2. quadratic_spline.f90:596, 694, 622 -- osplin reads out of bounds: >= 2 points all below the table
   reads eval_x(k+1) at 596; >= 2 points all above reads eval_x(0) at 694; a point exactly on the last node
   followed by points above reads table_x(n+1) at 622. osplin discards meval's err (37, 43); decreasing
   points leave yval undefined with no flag (the test pins zeros); n < 3 reads table_y(3) (73). Latent.
3. cubic_spline.f90:152-173 -- ysplin wrong for n = 3 (row 2 gets left-end coefficients, then right-end
   overwrite them; on (0,0),(1,1),(3,0) middle slope 0.9545 vs natural spline 0.5) and reads uninitialized
   values for n = 2. Reachable: kap/opal92/opal92_opacity.f90:363 accepts t6 equal to the top grid T; then
   lines 430-432 give jt = 3 and intpol at :500 returns a wrong d ln kappa/d ln T (the value is right).
   prepare_opal92_tables lets density rows with 2-3 points through (`if (jd <= 1) cycle`, :180).
4. linear_solvers.f90:98, 106 and :80-88 -- ludcmp ierr unreliable: a singular matrix without a zero row gets
   a 1e-20 pivot and ierr = 0 (confirmed on a rank-2 3x3), contradicting linear_solvers.f90:11 and README
   sec.3; a column whose candidate pivots are all NaN leaves imax unset (0 at j = 1) and the row swap at :88
   writes a(0,k). Caller element_diffusion/thoul_solver.f90:287.
5. linear_solvers.f90:251-252 -- solve_tridiagonal never checks the first pivot b(1): b(1) = 0 returns NaN
   with ierr = 0, against its header. rotation/seculr/am_transport_coefficients.f90:206-207 has the same gap.
6. table_interpolation.f90:415-420 -- locate's 1% slack is relative to the end value: for a decreasing
   table every out-of-range x gets clamped (expected_test_numerics.out:141); when xx(1) <= 0, an x exactly
   on xx(1) returns 0 (pinned at :132). No effect on today's callers (kap/laol89/laol_opacity.f90:267-272,
   select_four_point_window :626-631; atm/tables/allard_atmosphere.f90:206-208, 221-222 treat 0 like 1 and
   n like n-1): dead logic with a wrong contract.
7. table_interpolation.f90:244-246, 263-265 -- lir fails silently: < 2 points returns with intent(out)
   interp_flag undefined (the test pins the caller's -7 surviving, :119); a flat table returns result_y
   unset and interp_flag = 1. No ierr. Caller eos/mhd/mhd_eos.f90:637.
8. quadratic_spline.f90:74, 76, 111, 140 -- the 2025 safedivide guard protects only the last interval;
   other secants divide plainly; header (:55-56) overstates.
9. math/math_lib.f90:90-94, 113-116, 142 -- crmath and intrinsic backends disagree at the edges: crmath
   pow(0,y) returns 0 for y <= 0 (intrinsic 1 or Inf); negative x with integer-valued |y| >= 100 gives NaN
   under crmath; floor(y) overflows default integer for |y| >= 2^31; header (:136) calls exp10 "exact for
   integer x" (exact only up to 1e22). From MESA.
10. table_interpolation.f90:190 -- polint's coincidence test is an absolute 1e-20 on a dimensional
   difference (scale-dependent; harmless for the Allard grids).

## Findings (ranked)
1. bsstep limits (bulirsch_stoer.f90:63-144): check 1 <= num_eqs <= nmax, 2 <= extrap_order <=
   min(ncol, max_stage_index), max_stage_index <= size(substep_sequence), returning ierr; tie imax to
   size(substep_sequence); list the limits in README sec.2; drop the commented-out F77 line (:136);
   pow(2.0d0, int(...)) (:135) can be 2.0d0**k (removes numerics' use math_lib) -- check bit-identity. A, S.
2. osplin's error contract (quadratic_spline.f90:28-44, 519-726): ierr (meval's err = 2, n < 3); handle
   all-outside and last-node-then-above cases; add to the test. A for valid inputs, M.
3. Two natural cubic splines (cubic_spline.f90:126-264): ysplin/intpol duplicate natural_spline/splint in
   slope form; fixed np = 100 with explicit-shape c(4,np) and no n <= 100 check; a private Thomas solve
   (:174-185); constants reassigned every row (:152-173); stale "second order terms" comment (:151);
   meaningless "<1" (:227-228); unreachable error branch with an opacity message (:246-251). Proposed: an
   optional derivative for splint, drop intpol; build kap's coefficient tables from y2 (also fixes bug 3).
   laol's four cspline(...,1.0d30,1.0d30) calls (kap/laol89/laol_opacity.f90:208,345,479,584) are
   natural_spline calls (A). The sentinel is documented three ways (1e30 or more; 0.99d30 or more; README)
   while the code tests > 0.99d30; make end slopes optional arguments. B, M.
4. False successes (linear_solvers.f90:48, 98, 106, 251; table_interpolation.f90:244, 263): ludcmp flags a
   pivot below a relative tolerance instead of substituting 1e-20; solve_tridiagonal checks b(1); lir
   returns ierr. A unless a run hits the 1e-20 path. S.
5. Docs that disagree: README sec.3 and linear_solvers.f90:11 "singular matrix gives ierr = 1" (cspline,
   osplin, lir have no ierr at all); README:20 says solve_angular_momentum_transport is in
   secular_transport.f90 (it is rotation/seculr/am_transport_coefficients.f90:178); stale caller notes
   ("Cox" in README:12, cubic_spline.f90:12, :204, :247, table_interpolation.f90:348 -- users are OPAL92 and
   Kurucz; inter3 "Livermore opacity" (:25-26) -- callers atm/kurucz_atmosphere and eos/scv; interp
   "VandenBerg CAPPA" (:57) -- callers kap alex94, alex06, opal95); simeqc called "Gauss-Jordan"
   (linear_solvers.f90:7, :153) but it is Gaussian elimination with back substitution (IBM SSP SIMQ);
   test_numerics.f90:11-12 and :147 claim all four knot cases and "two knots per interval" (contradicts
   quadratic_spline.f90:11-13; a brute force reached neither case 3 nor 4); quadratic_spline.f90:547-549
   says out-of-range points get the end value but the code extrapolates (:571-585); README sec.1 helper list
   omits safedivide. A, S.
6. quadratic_spline.f90 still F77 inside: nested blocks unindented, end if closing blocks 50-100 lines up
   (meval :595, :640, :697; cases :180-238); ALL-CAPS comments naming dead variables (XVALS, XTAB, S, FIRST;
   :403-413, :448, :473); commented-out code (:112-113, :124-125); six "(Restructured 2026 ... unchanged)"
   boilerplate notes (:79, :364, :388, :442, :619, :657); cases and spline want select case; slopes computes
   the harmonic mean 2m1m2/(m1+m2) of the two secants -- say so. A, M.
7. test/test_numerics.f90 pins values, checks nothing: no ok/FAIL, no orders; it pins bugs (out:132, :141,
   :119). Add: cspline/splint exact on cubics with clamped ends and an interior O(h^4) check; inter3 and
   interp exact on quadratics/cubics; osplin keeps monotone and convex data monotone and convex; LU
   residuals and a rank-deficient matrix; tridiagonal vs a dense solve; bsstep vs sin/cos; ysplin vs
   natural_spline for n = 3, 4, 8; mutation checks. Drop use star_info_lib, only: json (:62). A, M.
8. Kernels print, to three destinations: stdout at linear_solvers.f90:49, :257, table_interpolation.f90:191,
   bulirsch_stoer.f90:138; the run log at cubic_spline.f90:116; both at :247-248; only for that numerics
   depends on io/io_units' mutable unit numbers (cubic_spline.f90:17). Return distinct codes and let callers
   report. A for results (terminal text changes -- check whether any pinned .log contains these lines), S.
9. Names: inter3, interp, intrp2, polint, lir, quint (six overlapping interpolators), findex/locate,
   cspline/natural_spline/ysplin/intpol, simeqc. Add a README "which routine for what" table and an "add an
   interpolation scheme" recipe (module, facade, a property test); descriptive facade names for YREC-authored
   routines (lagrange4_weights, bracket_index, ...); keep NR and IBM SSP code under its authors' names.
   A, S README, M renames.
10. Six interval searches: splint :102-111, intpol :236-245, locate :404-413, quadratic search :377-424,
   lir's linear scan :256-281, findex :352-385. One bisection plus one warm-start search (findex and kap's
   hunt). A/B, M.
11. table_interpolation.f90:58-122, 298-311 -- three 4-point Lagrange formulas: intrp2 is interp without
   dweight with identical expressions (make dweight optional, delete intrp2: A); lir's inline weights onto
   interp is B. S.
12. lir (table_interpolation.f90:222-343): assumed-size (*) args (:225-227); stride constant 1 so stride_m1 = 0
   and result_y(y_idx+stride_m1)=0 lines (:276, :327, :341) do nothing; y1..y4 hold abscissae, z1..z4
   differences (:298-305); diff means three things (:254, :313, :323); the diff == 0 test inside the search
   loop (:263); comments name dead NT, IL, M, K, A(K) (:243, :286, :312). A, S-M.
13. linear_solvers.f90:160-231 -- simeqc is flat-index SSP arithmetic, 16 cryptic integers, unindented; one
   caller net/network.f90:560 (n <= 8). Rewrite on a 2-D a(n, n+nrhs) keeping operation order (A, careful)
   or call LU (B). M.
14. Facade hygiene: numerics_termination (numerics_lib.f90:32) is the envelope integration's status (used in
   core/envelope, main, run_driver, yrec_capi) -> move to core/envelope or state/. 12 whole-module
   `use numerics_lib` break Declarations rule 2: core/step/shell_physics.f90:39 (uses nothing),
   core/envelope/envelope_integration_kernel.f90:156, core/diagnostics/turnover_timescale.f90:118,
   rotation/seculr/rotational_mixing.f90:169, am_grid.f90:60 and :316, equal_grid.f90:173,
   element_diffusion/original_settling.f90:75, thoul_solver.f90:62, settling_grid.f90:27 and :287,
   microdiff_run.f90:30. Have check_boundaries enforce only: (math_lib the exception). A, S.
15. math/math_lib.f90: delete safe_divide_by_exp (:187-202; no caller; magic 709.7827; .lt.); add a README;
   a small test_math under both backends (pow and exp10 integer paths, pow(0,.), negative bases); fix the
   exp10 "exact" claim; real(dp); drop ALL-CAPS prose. A, S.

## Duplicated numerics outside numerics/
- rotation/seculr/stability_coefficients.f90:226-250 lagrange_weights: same expressions, same order as intrp2. A.
- rotation/seculr/stability_coefficients.f90:257-264 interp4: = lagrange4 (left-to-right sum). A.
- rotation/seculr/circulation.f90:267-276 log_lagrange4: = exp(lagrange4(w, log(a(k:k+3)))); assumed-size a(*). A.
- eos/scv_eos_tables.f90:75-78 scv_weighted_sum4: = lagrange4. A. (EOS package's file -- leave to EOS.)
- rotation/seculr/am_transport_coefficients.f90:178-229 solve_angular_momentum_transport: a Thomas copy with a
  round trip on solution(n) between sweeps (:220-221); deliberate (ROADMAP numerics revamp); README path stale.
  If numerics split the tridiagonal solve into forward sweep + back substitution it could reuse them bit for
  bit (A); a plain solve_tridiagonal call is B.
- kap/conductive/conductive_opacity.f90:487-550 hunt: Numerical Recipes 1986; missing from ROADMAP sec.8's NR
  list; can move verbatim to table_interpolation. A.
- eos/opal/opal_eos.f90 and opal95_eos.f90: 12 inline bisections with +-1e-7 offsets (vendor code; leave
  unless each comparison is reproduced).
- eos/opal_eos_tables.f90:474-512 quad (vendor; leave).
- Fine where they are: kap/opacity_tables.f90:319, :355 (stencil4_locate; deliberately unmerged),
  kap/conductive/conductive_opacity.f90:430 (cinterp3, Potekhin), core/structure/triangle_fit.f90:196-230.
- No private spline or LU copies outside numerics.

## Structural recommendations
1. One error contract for numerics: named status codes (singular, bad input, out of range, step underflow);
   cheap precondition checks (sizes, monotonicity, bsstep limits); no I/O, so the io_units dependency goes;
   numerics_termination to the envelope.
2. One routine per job, chosen through a README table: one bracket search + one warm-start search (absorbing
   hunt); Lagrange-3/4 weights with optional derivatives plus general polint; one cubic spline
   (cspline/splint with optional derivative, replacing ysplin/intpol); osplin; LU; tridiagonal with factor
   and back-substitute; the ODE integrator. Property tests per scheme (exactness on polynomials of its
   degree, a measured order).
3. Keep math/ separate (dependency root below numerics; its own use rule; its concern is reproducibility):
   README (what lives there, two backends, usage contract, how to add a shadowed function), test_math,
   delete safe_divide_by_exp.

Fine as is: the facade, star-free directory, no SAVE; ratext's tableau as a per-step record; ode_system and
README sec.2; bsstep's NaN-skipping error maximum; the NR routines cspline, splint, lubksb, polint, mmid,
ratext are faithful copies (licensing known); inter3, interp, quint formulas and findex's contract;
ysplin correct for n >= 4; slopes keeps choose to cases 1 and 2; math_lib's two-backend design.

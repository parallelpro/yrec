# Validating rotation

How rotation/ is checked beyond the pinned decks: against analytic
solutions and physical laws (unit tests), in its sub-step (a
convergence study), against another code (MESA), and by line coverage;
and why the Solberg-Hoiland instability, added in 2026, is off by
default.
Measured on yrec-next, 2026-10-28.

## 1. Laws and analytic solutions

`test_rotation_kernels` (`make tests` in `src/`, then `pytest
test_rotation_kernels.py`) runs the kernels on synthetic inputs, with
each expected value formed in the test, not by the kernel. Each check
was confirmed to fail under a mutation of the code it covers. The
analytic ones:

| Kernel | Solution | Measured |
|---|---|---|
| `brake_convection_zone` (Kawaler), saturated | Omega_0 exp(-K Omega_sat^2 T) | Error (KT)^3/(12 N^2) to 1% in N calls: the implicit midpoint rule's leading term. Order 2.00. |
| the same, unsaturated (n = 3) | Omega_0 (1 + 2 K Omega_0^2 T)^(-1/2) | Order 2.00. |
| `am_diffusion_coeffs` + `solve_angular_momentum_transport`, constant D and I | The cosine modes of dOmega/dt = kappa d^2 Omega/dchi^2 with no flux at the ends | Each mode is an eigenvector of the implicit step to round-off. The decay rate converges to kappa (pi/L)^2 as h^2 (order 2.00); the amplitude to exp(-kappa (pi/L)^2 t) as dt (backward Euler, order 1.0). |
| `rotation_stability_setup` | Its documented coefficients, on cubics in r (exact under its 4-point interpolation) | 1e-11 or better. |
| `check_composition` | mu of a fully ionized gas: A/2 (H), A/3 (He4), 12/7 (C12) | Round-off. |

The others check conservation laws (J through the diffusion and the
map back, each species through the mixing), limits (the spherical
limit of the shape and the moments of inertia; no shear term or GSF
velocity for a solid body), scalings (U_ES and U_GSF ~ omega^2; 1 - fp
~ omega^2; the shear's quadratic onset) and the checks' contracts.
`src/rotation/README.md` section 4 lists them.

## 2. The sub-step

The rotation step is cut into sub-steps that each change omega by at
most `time_dw_mix`. `src/tools/substep_convergence.py` reruns a deck
with a sequence of values and compares the last models. On the
rotating reference variant (the testsuite's
`Test_solarGS98_difrotmixper18p4d1p0`, stopped at 300 Myr: transport,
the Matt wind, disk locking to 1 Myr), relative differences between
successive factors:

| time_dw_mix | 0.04 to 0.02 | 0.02 to 0.01 | 0.01 to 0.005 | 0.005 to 0.0025 | 0.0025 to 0.00125 |
|---|---|---|---|---|---|
| Surface omega | 1.1e-3 | 7.6e-4 | 5.3e-4 | -3.0e-4 | 1.3e-4 |
| log L | -5.0e-4 | -2.3e-4 | -6.1e-5 | -1.0e-4 | 3.8e-6 |

- **The sub-step error is small:** halving the factor from the
  default 0.01 moves the 300 Myr surface omega by 5e-4 and L by
  1e-4.
- **The order is not clean.** Each run chooses its own model
  timesteps from its own omega, so the runs differ by their step
  sequences too, at about 1e-4. At 200 Myr the total angular
  momentum's differences from the finest run halve with the factor
  (6.6e-3, 3.4e-3, 1.5e-3, 7.6e-4 and 2.7e-4 for 0.04 to 0.0025):
  first order, as the operator splitting of wind, transport and
  burning implies.
- **The J figures above predate the fix of O26**, when the history's
  J column lagged its model by a step (docs/bug_sweep.md O26), so its
  value at the stop depended on the last step's length; the age-matched
  figures at 200 Myr, where every run takes steps of about the same
  length, were affected far less.

## 3. Against MESA

`sample_cases/rotation_vs_MESA` runs a 1 Msun star as a solid body
with no wind, from the pre-main sequence to 1 Gyr, in both codes with
the same J. Under J conservation omega = J/I, so the comparison checks
the J-to-omega inversion and the moment of inertia, with and without
the rotational distortion.

- **Main sequence** (50 Myr to 1 Gyr): the surface omega agrees to
  0.8% at omega/omega_crit = 0.006 and to 0.4% at 0.25.
- **Pre-main sequence:** the differences follow the radius. The two
  stars contract on different schedules from different starting
  models.
- **J:** YREC conserves it to 2.3e-5 over the run, MESA to 3e-12.

The README there has the table and the steps.

## 4. Coverage

`src/tools/coverage.py --work DIR --deck ...` builds with -O0
--coverage in DIR and measures gcov line coverage of rotation/ under
the two test programs and the decks given. With the rotating
configuration-matrix decks, the rotating variants (50 Myr) and the
rescale-only rotating deck: 96.9% of 2168 executable lines (2026-10-29) (92.7% of
2203 before the checks' unreachable branches and the gradient-reversal
safeguard, which acted only above 1 rad/s, were removed: O15). Every
file is above 90%; the checks, the braking, the mechanisms' dispatcher
and the region's transport are at 100%.

## 5. Solberg-Hoiland

Pinsonneault (1988, section II.E) lists six mechanisms: three dynamical
(convection, the Solberg-Hoiland instability, the dynamical shear) and
three secular (Eddington-Sweet circulation, GSF, the secular shear).
The 2025 code had all but Solberg-Hoiland. It is now the fourth
transport slot, `compute_solberg_hoiland_velocity`, off unless
`solberg_hoiland_velocity_scale` is set (2026-10-28).

An interface is unstable where N^2 + kappa^2 < 0, with kappa^2 = 4
omega^2 + 2 omega d omega/d ln r and N^2 = N_T^2 + N_mu^2 (N_T is half
the dynamical-shear limit, N_mu^2 the secular shear's mu term). The
threshold, |d omega/d ln r| > (N^2 + 4 omega^2)/(2 |omega|), is never
below the dynamical-shear limit 2 N_T, since N^2 + 4 omega^2 >= 4 N_T
|omega|. So the instability can act only where the dynamical shear
already does, and its velocity there, (8/45) K (|d omega/d ln r|/q_SH)^2/r,
is at most the dynamical shear's. `test_rotation_kernels` checks the
criterion on single interfaces and the implication on a grid of 432
interfaces (N_T, N_mu, omega and d omega/d ln r).

- **Off (the default):** the 50 self-contained decks of the check set
  are byte-identical to the build before the slot, every output file.
- **On (scale 1):** 23 of the 26 rotating decks are byte-identical,
  among them two that report dynamical shear (3 and 4 times). The
  three that report it most move: rv_alpham05 (231 reports) by 1.7e-4
  dex in log L at most and 0.2% in the surface omega; rv_nocore (820)
  by 8.5e-5 dex and 0.3%; the solid-core rv_sbm2 (6591) by 3.1e-4 dex,
  4.4% in the final surface omega and 0.6% in J.

MESA's D_SH (D_SH_factor, default 0) is the counterpart; its source
notes that Langer's code said "DO NOT USE".

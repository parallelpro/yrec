# The physics validation tier

Tests that check YREC against something it did not produce: a law
re-derived from primitive columns, another code, a published number or
an observation of the Sun. The byte pins (docs/GUIDELINES.md,
"Verification tiers") freeze the output and catch any change; they do
not say whether the output is right. A wrong N^2 once passed every pin.
This tier is the complement (ROADMAP section 9).

Every check here has three properties:

- **its expected value comes from outside the code under test;**
- **its tolerance was measured**, is stated beside the measured value,
  and names what sets it (round-off, a difference scheme's order, a
  known difference between two codes);
- **it was seen to fail**, under a mutation of the quantity or of the
  code; each test file has a `can_fail` test that keeps this true.

Measured 2026-11-17 on valid-tier (yrec-next at 1179afbf), macOS/arm64,
gfortran 15.2, MESA 26.04.1, GYRE 9.1.1.

## Running it

| Test | Compares YREC with | CI | Needs | Time |
|---|---|---|---|---|
| `test_validation_contract.py` | the structure equations, from its own files | yes | the solar deck's tables | 45 s |
| `test_validation_gyre.py` | GYRE on its pulse files | no | `GYRE_DIR` | 80 s |
| `test_validation_solar.py` | the Sun | no | the full `input/`; `GYRE_DIR` for the modes | 3.7 min |
| `test_validation_mesa.py` | MESA (stored reference) | no | the solar deck's tables; `GYRE_DIR` for the modes | 100 s |
| `test_validation_rotation.py` | a published calibration, a budget, a law | no | the full `input/` | 3.2 min |

```
python3 -m pytest test_validation_contract.py                       # the CI part
GYRE_DIR=/path/to/gyre python3 -m pytest test_validation_*.py       # the whole tier, 11 minutes
python3 -m pytest -s test_validation_gyre.py                        # -s prints the measured values
```

- A test whose need is missing skips with the reason; nothing reaches
  for the network.
- `YREC_EXE` names another `yrec` to validate (a build with a change, an
  HDF5 build); the default is `src/yrec`.
- `YREC_TEST_HDF5=1` with an HDF5 build adds the GSM file to the contract
  and GYRE tests.
- `YREC_VALIDATION_REPORT=<file>` writes a test file's measured values.
- `python3 -I src/tools/validate_model.py FILE.GYRE|.FGONG|.GSM` audits
  one pulse file from the command line. The same module holds the
  readers and the re-derivations the tests use (standard library only).
- Reference data: `validation/` (the Sun, the rotation landmarks) and
  `sample_cases/solar_vs_MESA/mesa/` (MESA's side), 100 kB in all, each
  with a note on how it was made.

An expected failure (`xfail`, strict) marks a defect the tier found and
did not fix (section 6). Its check states the law at the tolerance a
correct code would meet. When the defect is fixed the test reports
XPASS, which fails the file until the mark is removed.

## 1. The structure and pulse-file contract

`test_validation_contract.py` runs the standard solar deck
(`examples/run_standard_solar_model`, `Test_solar_noGS_norot`) and
checks six of its models, from 2 Myr on the pre-main sequence to
3.4 Gyr, from the profile and the pulse files alone. Each row is the
largest value over the six models.

| Check | Measured | Tolerance | What sets it |
|---|---|---|---|
| Header: L = 4 pi R^2 sigma Teff^4, g = G M/R^2, L_star, R_star and the age against the history, a mesh point at R_star | 8e-15 | 1e-12 | round-off |
| r increases strictly, m never falls, P falls strictly, T never rises, all positive | holds | exact | |
| Interior, hydrostatic equilibrium and mass continuity in YREC's own scheme (trapezoid in log m) | 1.3e-6, 2.0e-6 | 6e-6 | the solver's convergence on this deck |
| The same one order higher (three-point rule) | 3.5e-4, 2.6e-4 | 1.1e-3, 8e-4 | the scheme's truncation, h^2/12 with h = 0.046 in ln P |
| Interior, dlnT/dlnP = nabla | 4.2e-4 | 1.3e-3 | nabla is centred on the pair in the scheme, at the points in the file |
| Envelope, the same three laws in ln P (three-point rule) | 5.7e-4, 3.8e-4, 6.2e-3 | 1.8e-3, 1.2e-3, 2e-2 | the rule's error on the output mesh (0.05 dex in P); the envelope itself is integrated to 1e-9 |
| Atmosphere, heights against hydrostatic equilibrium | 6.1e-3 | 1.8e-2 | plane-parallel gravity (2.5e-3 at the top) and the depth quadrature (V11) |
| Atmosphere, dlnT/dlnP = nabla in radiative layers; m = M and L = L_star exactly | 3.1e-3 | 1e-2 | the pair rule on the T(tau) mesh |
| Luminosity equation, (L - int (eps_nuc + eps_neu + eps_grav) dm)/L; L = eps m at the first point | 4.6e-4; 8.5e-4 | 1.5e-3; 2.5e-3 | this deck's solver tolerances (V1); a converged model reaches 4e-6 |
| N^2 against the rule the file documents (thermal form where negative, else the centred density form) | 1.5e-15 | 1e-9 | round-off of a small difference |
| N^2's two forms, dlnrho/dlnP - 1/Gamma_1 against delta (grad_ad - grad), in uniform radiative layers | 9.0e-4 | 3e-3 | the centred difference over one shell |
| N^2 > 0 in radiative layers, < 0 where nabla > nabla_ad | holds | exact | |
| Lamb frequency, sound speed, gravity, gradr/grada against their definitions; the profile's N^2 is the pulse file's | 1.2e-14 | 1e-12 | round-off |
| The same quantity in the profile and the GYRE file (r, m, P, T, rho, L, Gamma_1, the gradients, delta, kap, omega) | 7.4e-15 | 1e-12 | round-off of the profile's log columns |
| kap_kap_T and kap_kap_rho against the opacity's own change along the model | 1.4e-2 (interior), 0.41 (above) | 4e-2, 1 | centred differences; the left side reaches 8 in the hydrogen ionization zone |
| eps_eps_rho/eps and eps_eps_T/eps in burning shells; zero above the fitting point | 1.02 to 1.08; 4.9 to 13.8 | 0.95 to 1.2; 3 to 25 | two-body reactions; He3 out of equilibrium |
| Fitting point, the jumps in ln r, ln T, ln rho beyond the laws; the step in m/M | 5.1e-6, 2.0e-4, 3.3e-4; 2.3e-8 | 2e-5, 6e-4, 1e-3; 1e-7 | the surface boundary condition is interpolated in the envelope triangle (a triangle 10 times smaller gave jumps 100 times smaller) |
| The interior-only GYRE file is the full file's first nz points | exact | exact | |
| FGONG: the globals, each shared column, m from ln(m/M), N^2 from A4 | 3e-16 | 1e-12 | round-off |
| GSM (HDF5 build): version 110, the header and the 18 columns | exact | 1e-15 | |

What a check can see follows from its tolerance: an error of 1e-5 in
the interior's hydrostatic balance (GYRE's G for YREC's is 2.7e-4), of
2e-3 in the envelope's, of 2% in the atmosphere's.

`test_each_check_can_fail` applies 19 mutations to a copy of a model
and requires the named check to fail. Among them are states this code
has been in: N^2 from the thermal form only; the envelope's actual and
adiabatic gradients swapped; atmosphere heights not accumulated; the
opacity derivatives zero above the fitting point; the energy
derivatives multiplied by eps once more. In scratch, a build with the
first and third of these in the Fortran failed `check_brunt`,
`check_monotonic` and `check_atmosphere_equations`.

## 2. The seismic observables against GYRE

`test_validation_gyre.py`: the standard solar model at 4.568 Gyr, its
pulse files in GYRE (adiabatic), and the history's seismic columns.

| | YREC | GYRE | |
|---|---|---|---|
| Formats read | GYRE text, FGONG, GSM | all three; FGONG's p modes within 1.4e-12 uHz of the text file's, GSM's equal | tolerance 1e-9 uHz |
| `delta_nu` (sound travel to the photosphere) | 142.677 uHz | radial large separation, n = 15..25: 135.770 | ratio 1.051; band 1.03 to 1.07 (V7) |
| The same integral through the atmosphere | 137.832 uHz | n = 25..31: 138.070 | -0.17%; tolerance 0.6% |
| `delta_nu_rho` (scaling) | 135.115 uHz | 135.770 | -0.48%; tolerance 1.5% |
| `delta_Pg` (buoyancy sum) | 1641.50 s | l = 1 period spacing, n_g = 120..200, centre resolved: 1546.79 | +6.1% (V5, expected failure) |
| The buoyancy integral with the centre's part | 1547.11 s | 1546.79 | +2.0e-4; tolerance 3e-3 |
| GYRE on the file as written | | 1510.36 s | -2.4% (V6, expected failure) |

- `delta_nu` and `delta_Pg` are the integrals observables.f90
  documents, re-derived from the pulse file to 1e-9; `nu_max` and
  `delta_nu_rho` are their definitions to 1e-9.
- GYRE must be told YREC's G (`&constants G_GRAVITY = 6.6725218e-8`):
  the GYRE text and GSM files carry none, and GYRE's own, 6.67430e-8,
  is 2.7e-4 away. FGONG carries it.
- GYRE reads N^2 from the file, so it cannot tell whether the file's
  N^2 is right; the contract test's two-forms check and the MESA
  comparison do that. What GYRE checks is that YREC's scalars describe
  the file.
- `test_checks_can_fail` gives GYRE a file with Gamma_1 2% high and one
  with N^2 5% high; each comparison fails.

## 3. The Sun

`test_validation_solar.py`: the calibrated solar model with the GS98
mixture and diffusion (`testsuite/Test_solarGS98_dif`) at 4.568 Gyr.

| Quantity | Observed (source) | Model | Offset | Band | A standard model with GS98 (literature) |
|---|---|---|---|---|---|
| L | 3.828e33 erg/s (IAU 2015 B3) | | +6.5e-5 | 1e-3 | calibrated |
| R | 6.957e10 cm (IAU 2015 B3) | | -1.6e-4 | 1e-3 | calibrated |
| Teff | 5772 K (IAU 2015 B3) | 5772.5 K | +0.5 K | 1 K | |
| R_cz/R | 0.713 +- 0.001 (Basu & Antia 1997) | 0.7127 | -0.3 sigma | 3 sigma | 0.713 to 0.714 (Bahcall et al. 2005; Serenelli et al. 2009) |
| Surface Y | 0.2485 +- 0.0035 (Basu & Antia 2004) | 0.2434 | -1.5 sigma | 3 sigma | 0.242 to 0.243, -1.6 to -1.8 sigma (the same) |
| Surface Z/X | 0.0229 (Grevesse & Sauval 1998) | 0.022894 | -0.02% | 1% | calibrated |
| c(r), 0.1 to 0.95 R, 71 radii | MDI inversion (Rhodes et al. 1997) | | rms 8.9e-4, largest 2.7e-3 | 2e-3, 5e-3 | rms about 1e-3 (the same) |
| nu(l = 0), n = 6..10 | BiSON (Broomhall et al. 2009) | | +3.9e-4 to +9.4e-4 | 2e-3 | |
| Large separation, n = 15..25 | 134.59 uHz (the same fit to BiSON) | 135.84 uHz | +0.93% | 0 to 2% | the surface effect: +8.9 uHz at n = 21 |

The literature values are quoted from memory of those papers' tables
and should be checked against them before they are cited further; the
observed values and their sources are in the test's header.

The checks can fail: the standard solar deck without diffusion has the
same L, R and Teff and fails three (R_cz/R = 0.7267, 13.7 sigma; Y =
0.2638, +4.4 sigma; sound speed rms 4.2e-3).

## 4. Against MESA

`test_validation_mesa.py` and `sample_cases/solar_vs_MESA` (its README
has the set-up, the stored files and how to remake them): one 1 Msun
track with the standard solar deck's physics in both codes.

YREC minus MESA (dex; the age relative):

| At | log L | log Teff | log R | log T_c | log rho_c | age |
|---|---|---|---|---|---|---|
| X_c = 0.6 | -0.0094 | -0.0031 | +0.0016 | -0.0032 | -0.0046 | +2.6% |
| X_c = 0.35 | -0.0076 | -0.0030 | +0.0022 | -0.0026 | -0.0009 | +4.0% |
| X_c = 0.1 | +0.0111 | -0.0021 | +0.0098 | +0.0030 | +0.0336 | +7.7% |
| 4.568 Gyr | -0.0141 | -0.0034 | -0.0003 | -0.0058 | -0.0145 | (X_c +0.0148) |
| Band | 0.030 | 0.007 | 0.020 | 0.012 | 0.070 | 15% |

At 4.568 Gyr GYRE's frequencies of the two models differ by 1.6e-3 to
2.3e-3 (band 5e-3), the radial large separation by 3.1e-4 (band 2e-3),
and the codes' own `delta_nu` by 0.54% (band 1%).

The bands are about twice the measured differences: they bound how far
apart two independent codes are here, and YREC may move inside them.
What explains the differences, measured (the README has the numbers):

- **The opacity, more than half of the luminosity offset.** At equal
  temperature along the two solar-age models YREC's opacity is higher
  by 0.3% at log T = 7.0, rising to 3.1% below the convection zone; the
  two codes' OP tables differ by 0.7 to 1.6% at their common nodes
  between log T = 6.2 and 6.8 (V8). MESA with its opacity raised by
  2.5% there loses 0.0059 dex in log L at X_c = 0.6 (63% of the offset).
- **Not the equation of state:** P/(rho T) agrees to 3e-4.
- **Not separated:** the nuclear rates, the opacity above log T = 6.9,
  the mixing-length formulations and the difference schemes: 0.0035 dex
  in log L and 0.0019 in log Teff at X_c = 0.6.

This check is coarse by nature: it would see a 5% error in luminosity
or a 1% error in mean density, not the 1e-5 changes of most fixes.

## 5. Rotation

`test_validation_rotation.py`; docs/ROTATION_VALIDATION.md section 6
has the account. On the testsuite's rotating solar deck:

| Check | Measured | Tolerance |
|---|---|---|
| The paper's rotation period, 25.4 d at 4.568 Gyr | 25.412 d | 1% (the paper's) |
| The paper's L, R (1e-4) and Z/X (1% of 0.02292) | +7.0e-5, +3.8e-5, -0.01% | 2e-4, 1% |
| The paper's Li7 depletion, 2.3 dex within 1% | 2.086 dex | V9, expected failure |
| J(t) = J_0 - int torque dt, the torque from the published law and the history | 4.9e-3 of J_0 | 1.5e-2 |
| Per step on the main sequence, J lost over torque dt | 0.9993 (median) | 3e-3 |
| Spin-down exponent, 1 to 4.568 Gyr, against a solid body under the law | -0.677 against -0.689 | 0.03 |
| That body with frozen structure against the law's analytic solution | 1.5e-3 | 4e-3 |
| omega(m) at five ages against the stored landmarks | reproduced | 1% (a pin by tolerance) |
| Interior over surface rotation at the solar age | 1.04 to 1.06 (the Sun: 0.94) | 0.85 to 1.15 |

## 6. Defects and surprises the tier found

None is fixed here (this round changes no result). V-numbers are this
document's; M317 is docs/bug_sweep.md's.

| | What | Where | Size | The check that shows it |
|---|---|---|---|---|
| M317 (known, open) | The star's mass in grams is not star_mass times solar_mass_cgs; the profile's mass column jumps at the photosphere | the solar starting models; `stitched_model.f90` `atmosphere_value` | 5.03e-5 | contract: two expected failures |
| V1 | The standard solar deck's solver tolerances (3e-7 to 9e-7) leave the luminosity equation open, in proportion to the time step | `examples/run_standard_solar_model/Test_solar_noGS_norot` (HTOLER 6e-7; the other three solar examples and the testsuite have 6e-9) | 4.6e-4 of L at 3.4 Gyr; log L at 4.568 Gyr moves by +1.4e-4 dex with tolerances 1000 times tighter (the solar calibration's own tolerance is 3e-5); 8.9e-7 on `Test_solarGS98_dif`, whose tolerances are 100 times tighter | contract: `test_luminosity_equation_is_converged` |
| V2 | The stitched model keeps the envelope's repeat of the fitting point when its radius lands above the interior's: two points 1e-6 R apart whose density rises outward | `stitched_model.f90` `build_stitched_model` (the `<=` test on log R) | ln rho +1.1e-4 to +3.3e-4 over dr/r = 1e-6 to 5e-6, in 5 of 6 models; harmless here (the layer is convective, so N^2 is the thermal form); where the fitting point is radiative the centred N^2 is differenced across the pair (not measured) | contract: `test_no_density_inversion_at_the_fit` |
| V3 | In atmosphere layers that are convectively unstable the files' nabla is grad_ad, while the temperature follows the T(tau) relation | `envelope_integration_kernel.f90` (the atmosphere's `pulse_gradient`) | 15% in nabla just above the photosphere of the 2 Myr model (Teff = 4436 K); none in the solar-age model | contract: `test_atmosphere_gradient_everywhere` |
| V4 | Where YREC sets nabla = nabla_ad exactly (a convective core), the thermal N^2 is exactly 0, the rule's `>= 0` test takes the centred density difference, and N^2 is that difference's error: positive through the outer half of the core | `stitched_model.f90` `difference_brunt_n2` | 80 points of the 32 Myr model's convective core, up to 6e-5 of g^2 rho/P (N^2 = +5e-11 s^-2) | contract: `test_brunt_is_not_positive_in_convection` |
| V5 | `delta_Pg` leaves out the part of the buoyancy integral inside the innermost mesh point (N(r_1), since N ~ r there) | `observables.f90` `compute_seismic_observables` (the sum starts at point 2) | +6.1% on the standard solar model (1641.5 s for 1547 s); +2.7% on `Test_solarGS98_dif`, whose `extend_core_inward` puts the innermost point at 0.0043 R (1568.5 s for 1527 s) | gyre: `test_delta_pg_against_gyre` |
| V6 | The pulse files' innermost point is at r = 0.0097 R (0.0043 R with `extend_core_inward`), so GYRE's g-mode periods depend on how it fills the centre | the mesh's inner fitting point, 1e-4 M | GYRE's l = 1 period spacing 1510.4 s as written, 1546.8 s with the centre resolved (-2.4%); -1.2% with `extend_core_inward` (1509.2 s for 1526.8 s); p modes move by 0.002 uHz | gyre: `test_g_modes_do_not_depend_on_the_centre_fill` |
| V7 | `delta_nu` is 5.1% above the radial large separation one observes | `observables.f90`: the integral stops at the photosphere | 142.68 against 135.77 uHz; 137.83 through the atmosphere. MESA's `delta_nu` has the same definition (141.89 for its model) | gyre: the band 1.03 to 1.07 |
| V8 | YREC's OP opacity table is not MESA's `OP_gs98` | `input/opacity/atomic/op05/GS98.OP17` | equal at log T >= 7.0, 0.7 to 1.6% higher at log T = 6.2 to 6.8 at X = 0.7, Z = 0.02; 3.1% higher below the convection zone along the solar models; 0.0059 of the 0.0094 dex in log L against MESA | `sample_cases/solar_vs_MESA/README.md` |
| V9 | The rotating solar deck no longer gives the paper's lithium | `testsuite/Test_solarGS98_difrotmixper18p4d1p0` (`mixing_velocity_scale = 0.371`) | 2.086 dex for the paper's 2.3 dex within 1% | rotation: `test_published_lithium_depletion` |
| V10 | Two decks build the fixed-Z opacity table at a Z the model does not have | `opal95_single_table_z = 0.016232` with Z = 0.016492 in `Test_solar_noGS_norot` and `testsuite/Test_solarGS98_base` | none that matters: log L at 4.568 Gyr moves by 1.2e-6 dex between the two table Z (the benchmark deck, in scratch), so the table is at most a speed-up that never applies | measured in scratch |
| V11 | The atmosphere's heights depart from hydrostatic equilibrium | plane-parallel gravity; the depth step dtau/<rho kappa> | up to 0.6% of a step (6e-6 R over the atmosphere) | contract: the atmosphere tolerance |
| V12 | The GYRE text and GSM files carry no G, and YREC's is not GYRE's | `physical_constants.f90` (log G = -7.17571) | 2.7e-4 in G | section 2 |

## 7. What this tier would have caught

Result-changing bugs from docs/bug_sweep.md and the pulse writers'
history (`stitched_model.f90`, Sources), each tried as a mutation or
measured.

| Bug | Would a check here have seen it? | At what size |
|---|---|---|
| The N^2 bug: N^2 from the thermal form only, without the composition term (12% in a subgiant's period spacing) | Yes: `check_brunt` (the file against its rule, and the two forms). Tried as a mutation of the model and of the Fortran | the structure's N^2 is 7 times the thermal one at the solar centre; the tolerance is 1e-9 |
| The envelope's gradient slots swapped, flipping N^2 across the envelope (Sources 5) | Yes: `check_envelope_equations` (dlnT/dlnP = nabla) and `check_brunt`. Mutation | the two differ by up to 0.3 in the superadiabatic layer; the tolerance is 2e-2 |
| kap_kap_T and kap_kap_rho zero above the fitting point (Sources 3) | Yes: `check_opacity_and_energy_derivatives`. Mutation | up to 8 against 1 |
| The eps derivatives wrong by a factor eps (Sources 4) | Yes: the same check. Mutation | a factor 15 at the solar centre; the band is 3 to 25 |
| GSM stamped with the text format's version, which GYRE refuses | Yes, with an HDF5 build: `test_gsm_file` (version 110) and GYRE's reading of every format | exact |
| O22: the profile's N^2 column disagreed with the pulse files' | Yes: `check_derived_columns` requires them equal | sign flips in convection zones against exact equality |
| M317: 1.98900e33 g for a star_mass of 1 | Yes: found again by the header check, before the bug sweep's entry was known to this work | 5.03e-5 against 1e-9 |
| O26: the history's J one step behind its model | Yes: the rotation budget. Tried on the history | 2.4e-2 of J_0 against 1.5e-2; the median step ratio 1.021 against 3e-3 |
| O17: the solid-body Matt torque without 1/G | Only on a deck that takes that path with a wind (the budget would be off by orders of magnitude); the rotating solar deck does not | |
| M278: the conductive opacity's weights (the solar log L by -1.7e-5) | No | 1.7e-5 dex against bands of 1e-3 (the Sun) and 0.03 (MESA) |

The last row is the tier's limit: it sees wrong physics at the size of
a percent, or exactly where a law is exact. A change of 1e-5 in a solar
model remains the byte pins' and the measured commits' business.

## 8. What is not covered

- One mass, one composition: the contract's six models are 1 Msun.
  Radiative envelopes, convective cores on the main sequence (V4 was
  seen only in a transient one), giants and helium burning are not
  checked. `check_*` in the contract test take any model, so a second
  deck is a fixture away.
- Rotating models' structure equations (the f_P and f_T factors) are
  not in the contract.
- Nonadiabatic pulsation is not run; the derivative columns are checked
  only for consistency with the model.
- The cross-code reference is one MESA version with one set of options.
- No van Saders & Pinsonneault (2013) track: the repository does not
  hold that grid's parameters.

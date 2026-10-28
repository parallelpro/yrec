# Bug sweep on yrec-next, September to October 2026

## 1. About this report

This report lists every defect found on the `yrec-next` branch after it diverged from `main` at [0ff84d7e](https://github.com/parallelpro/yrec/commit/0ff84d7e) (2026-09-03), up to [88b32a90](https://github.com/parallelpro/yrec/commit/88b32a90) (2026-10-28). There are 64 items: 41 fixed and 23 open. Bugs fixed before the divergence (bug-sweep batches 1-3, such as the CF88 deuterium rate) are in `docs/MODERNIZATION_CHANGES.md` and are not repeated here.

Each item answers the same questions:

- **In the original F77?** Each item was checked against commit [6cd56739](https://github.com/parallelpro/yrec/tree/6cd56739/src) ("First commit", 2025-04-10), the import of the Fortran 77 code. "Yes" means the same logic is there, with a link and a quote. "No" names the commit that introduced the defect: either a 2025 YREC-team commit or a 2026 modernization commit.
- **What was (or is) wrong**, and the **expected behaviour, and why**.
- **What we did:** the fix commit and its measured effect. For an open item, what has been recorded or tested so far, and the proposed fix.
- **Now:** where the code is today.

**How the fixes were verified.**

- A pure refactor must be byte-identical on the test battery:
  - 37 pinned runs (the legacy decks' logs, and the MESA-format baselines' history and profile files);
  - a 21-deck quick check;
  - the Kurucz, Allard and 0.03 Msun decks;
  - the unit tests.
- A fix that changes results is its own commit. Its message records the drift, as matched-age differences in log L, log Teff, log R, X_c and surface Li, and the pins are reseeded in that commit.
- "Byte-identical" below means that no pinned file changed.
- Several fixes are protected by unit tests that fail if the fix is reverted: `test_structure` (every Henyey partial against finite differences), `test_mass_change`, `test_net`, `test_eos`, `test_triangle_fit` and `test_braking`.

**Conventions.**

- F77 links point at 6cd56739, and current-code links at e8549036.
- F77 routines are in capitals (CRRECT, SURFBC); current names are in code font (`compute_surface_boundary`).
- Numbers quoted from commit messages were measured at the time of the commit.

### Where the defects came from

| | In the F77 code | 2025 YREC-team commits | 2026 modernization | Total |
|---|---|---|---|---|
| Fixed | 31 | 2 (S5, L9) | 7 (N4, S13, L1, L4, L6, L7, L17) | 40 |
| Open | 20 | 1 (O4, with a 2026 part) | 2 (O5, O22) | 23 |

In addition, half of O10 (an overshoot gate) is a 2026 translation slip.

### Decisions needed from the authors

These open items depend on what the code was meant to do, not only on what it does:

- **O20:** is the calibration's dL/dα 0.0139 or 0.139?
- **O15:** is the gradient-reversal safeguard wanted? It has been inactive since the F77 code.
- **O12:** what tolerance should the solid-body Newton iteration stop at?
- **O13:** should the end cell of a transport region include an adjacent convection zone, both in its sums and in its write-back?
- **O1, O2, O3:** may the low-mass convergence work change the Jacobian? It would move every deck at the solver tolerance.

### Summary

| ID | Defect | In F77? | Status |
|---|---|---|---|
| **N** | **Nuclear network and burning (§2)** | | |
| N1 | Two copies of the d(p,γ) rate, 0.68% apart | Yes | Fixed, [05808dd7](https://github.com/parallelpro/yrec/commit/05808dd7) |
| N2 | S-factor derivative scales of reactions 2-7 divided by the pp S₀ | Yes | Fixed, [833b1a89](https://github.com/parallelpro/yrec/commit/833b1a89) |
| N3 | Rate temperature derivatives omit the S-factor scales | Yes | Fixed, [89c1fa42](https://github.com/parallelpro/yrec/commit/89c1fa42) |
| N4 | `LNEWNUC = .false.` gave zero pp/CNO rates | No (2026) | Fixed, [fff3deb5](https://github.com/parallelpro/yrec/commit/fff3deb5) |
| N5 | ¹⁴N(α,γ) Q value lost the positron energy, not the neutrino's | Yes | Fixed, [664c2fbc](https://github.com/parallelpro/yrec/commit/664c2fbc) |
| N6 | Li/Be burning never saw the end-of-step rates | Yes | Fixed, [0011e3a8](https://github.com/parallelpro/yrec/commit/0011e3a8) |
| N7 | Rotation sub-step burned with stale rates | Yes | Fixed, [4da5c0b4](https://github.com/parallelpro/yrec/commit/4da5c0b4) |
| N8 | DBURNM mixed grams with a dimensionless 1 | Yes | Fixed, [29048ffd](https://github.com/parallelpro/yrec/commit/29048ffd) |
| N9 | DBURN single-shell branch used an unset region mass | Yes | Fixed, [29048ffd](https://github.com/parallelpro/yrec/commit/29048ffd) |
| N10 | Mathematical constants zero until SETUPS ran | Yes | Fixed, [529d3fb3](https://github.com/parallelpro/yrec/commit/529d3fb3) |
| **S** | **Structure, timestep and boundary conditions (§3)** | | |
| S1 | Ω(j) inversion in the Newton loop had no convective flags | Yes | Fixed, [2ba4e32b](https://github.com/parallelpro/yrec/commit/2ba4e32b) |
| S2 | Mass change read m(0) for a fully convective star | Yes | Fixed, [bc874267](https://github.com/parallelpro/yrec/commit/bc874267) |
| S3 | ε_grav partials lacked 1/χ_T | Yes | Fixed, [18869257](https://github.com/parallelpro/yrec/commit/18869257) |
| S4 | End-age timestep floor in years instead of Gyr | Yes | Fixed, [38610071](https://github.com/parallelpro/yrec/commit/38610071) |
| S5 | HR-motion limit read start-of-step L at the wrong point | No (2025) | Fixed, [eee659ce](https://github.com/parallelpro/yrec/commit/eee659ce) |
| S6 | Atmosphere integration's first slope lacked f_P | Yes | Fixed, [e2c526aa](https://github.com/parallelpro/yrec/commit/e2c526aa) |
| S7 | Triangle vertex integration overwrote the model's Teff | Yes | Fixed, [3d065901](https://github.com/parallelpro/yrec/commit/3d065901) |
| S8 | Wind removed the convection zone's J twice | Yes | Fixed, [91a6b567](https://github.com/parallelpro/yrec/commit/91a6b567) |
| S9 | Wind and accretion extrapolated Ω ∝ r^(−walpcz) | Yes | Fixed, [d68c375e](https://github.com/parallelpro/yrec/commit/d68c375e) |
| S10 | Z ramp scaled CNO by the fraction of Z removed | Yes | Fixed, [1bcdee96](https://github.com/parallelpro/yrec/commit/1bcdee96) |
| S11 | Reused model rescaled twice; species target overwritten | Yes | Fixed, [df7ca6df](https://github.com/parallelpro/yrec/commit/df7ca6df) |
| S12 | Return from the gray atmosphere kept the old triangle | Yes | Fixed, [67d711a7](https://github.com/parallelpro/yrec/commit/67d711a7) |
| S13 | Envelope refit used L_sun where log L_sun is needed | No (2026) | Fixed, [4f34e95f](https://github.com/parallelpro/yrec/commit/4f34e95f) |
| **L** | **Latent bugs and error paths (§4)** | | |
| L1 | Pulse ε_T, ε_ρ columns multiplied by ε twice | No (2026) | Fixed, [1d664520](https://github.com/parallelpro/yrec/commit/1d664520) |
| L2 | Stable rotation sub-step kept a stale redo flag | Yes | Fixed, [66d7788f](https://github.com/parallelpro/yrec/commit/66d7788f) |
| L3 | Three OPAL EOS vintages shared one "loaded" flag | Yes | Fixed, [85f74924](https://github.com/parallelpro/yrec/commit/85f74924) |
| L4 | Saha cutoff default 6.0 lost | No (2026) | Fixed, [8a69dc72](https://github.com/parallelpro/yrec/commit/8a69dc72) |
| L5 | I/O errors ended in runtime errors | Yes | Fixed, [a9de5f99](https://github.com/parallelpro/yrec/commit/a9de5f99) |
| L6 | 188 documented control defaults wrong | No (2026) | Fixed, [c3f40675](https://github.com/parallelpro/yrec/commit/c3f40675) |
| L7 | Rename script changed a routine name inside a message | No (2026) | Fixed, [841ff500](https://github.com/parallelpro/yrec/commit/841ff500) |
| L8 | Thirteen overrun, divide-by-zero and edge-case guards | Yes (11 of 13) | Fixed, [0d66503f](https://github.com/parallelpro/yrec/commit/0d66503f) |
| L9 | Envelope reversal swapped the luminosity one way | No (2025) | Fixed, [4ce66e49](https://github.com/parallelpro/yrec/commit/4ce66e49) |
| L10 | Truncated OPAL95 table returned success | Yes | Fixed, [d62e5360](https://github.com/parallelpro/yrec/commit/d62e5360) |
| L11 | Dead branch set a velocity to `1.0D0-11` (= −10) | Yes | Fixed, [67d49e77](https://github.com/parallelpro/yrec/commit/67d49e77) |
| L12 | Ambiguous triangle-reset condition | Yes | Fixed, [4ce66e49](https://github.com/parallelpro/yrec/commit/4ce66e49) |
| L13 | LDIFAD transport never worked | Yes | Fixed (stops), [15d569e3](https://github.com/parallelpro/yrec/commit/15d569e3) |
| L14 | Envelope triangle could hang or go NaN | Yes | Fixed, [d2ff6886](https://github.com/parallelpro/yrec/commit/d2ff6886) |
| L15 | Diagnostics always zero or never written | Yes (6 of 8) | Removed, [c4f359c8](https://github.com/parallelpro/yrec/commit/c4f359c8) |
| L16 | Rezone's ∇_ad − ∇ store for rotation was dead | Yes | Removed, [56aa2e6f](https://github.com/parallelpro/yrec/commit/56aa2e6f) |
| L17 | `cz_moment_of_inertia` read before set (`intent(out)`) | No (2026) | Fixed, [2d62df76](https://github.com/parallelpro/yrec/commit/2d62df76) |
| **O** | **Open: structure, solver and composition (§5)** | | |
| O1 | Jacobian gaps: convective gradient, OPAL ρ derivatives | Yes | Open |
| O2 | EOS second derivatives with OPAL and in the Saha ramp | Yes | Open |
| O3 | Low-mass convergence (0.03, 0.05, 0.1 Msun) | Yes | Open |
| O4 | Starting model's stitched atmosphere reads out of bounds | No (2025/2026) | Open |
| O5 | Accretion entropy term at levels 1-2 (latent) | No (2026) | Open |
| O6 | Rezone's deuterium gate reads the old surface index | Yes | Open |
| O7 | Accretion mixing assumes a uniform convection zone | Yes | Open |
| O8 | Li/Be rate store zeroes both rate sets | Yes | Open |
| O9 | Surface ⁷Li noise floor of ~0.005 dex | Yes | Open |
| O10 | Overshoot edge width from two models; wrong overshoot gate | Width yes, gate no | Open |
| **O** | **Open: rotation, wind, microphysics and output (§6)** | | |
| O11 | Dynamical-shear limit can overshoot or go negative | Yes | Open |
| O12 | Solid-body Newton tolerance never set | Yes | Open |
| O13 | Transport region ends: sums and write-back disagree | Yes | Open |
| O14 | Species 5-15 unchecked after mixing; 16-slot history | Yes | Open |
| O15 | Reversal safeguard inactive; `mu_gradient_scale` inert | Yes | Open |
| O16 | Settling rewrites m in place; reversed `MOD` | Yes | Open |
| O17 | Magnetic braking: 1/G, wind switch, τ_cz | Yes | Open |
| O18 | `beta14`, opacity ramp derivatives, OPAL 2001 edge, MHD test | Yes | Open |
| O19 | Allard reader's misspelled flag; failed lookup used | Yes | Open |
| O20 | Calibration dL/dα 0.0139 against the recorded 0.139 | Yes | Open |
| O21 | Settling stencil out of bounds for a 2-point grid | Yes | Open |
| O22 | Profile N² column disagrees with the pulse files | No (2026) | Open |
| O23 | Two EOS derivatives threaded through but unused | Yes | Fixed, [6252b537](https://github.com/parallelpro/yrec/commit/6252b537) |
| O24 | Un-aged rotating model without transport spun up by the wind | Yes | Open |


## 2. Nuclear network and burning (fixed; results change)

Nine of these ten are in the original F77 code. N4 was introduced by the 2026 modernization.

### N1. Two copies of the d(p,γ)³He rate differed by 0.68%

- **Status:** fixed in [05808dd7](https://github.com/parallelpro/yrec/commit/05808dd7) (2026-10-10).
- **In the original F77?** Yes. DEUTRATE (the deuterium depletion) and ENGEB (the energy release) each coded the CF88 fit with a different divisor ([deutrate.f:34-40](https://github.com/parallelpro/yrec/blob/6cd56739/src/deutrate.f#L34-L40), [engeb.f:518-526](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L518-L526); `ANUC(3) = 2.013553`):
  ```fortran
  C FACTOR OF 3.0115D23 REFLECTS AVAGADROS NUMBER DIVIDED BY THE
  C MASS OF THE DEUTERON IN AMU
        RDEUT = RHO*2.240D3*T9P23*EXP(Z)*TFACDEUT*3.0115D23          ! deutrate.f
          RDEUT = RHO*2.240D3*T9P23*EXP(ZZ)*TFACDEUT*6.023D23/ANUC(3) ! engeb.f
  ```
- **Now:** `deuterium_proton_capture_rate`, [reaction_rates.f90:187-204](https://github.com/parallelpro/yrec/blob/e8549036/src/net/reaction_rates.f90#L187-L204), called by both the light-element burning and the energy generation.
- **What was wrong:** DEUTRATE divided N_A by the integer mass number (3.0115d23 = 6.023d23/2), although its own comment says "the mass of the deuteron". ENGEB divided by the nuclear mass, 2.013553 amu. The two differ by 2/2.013553 = 0.99327, in the rate and in the convective-overturn cap alike.
- **Expected behaviour, and why:** the deuterium destroyed and the energy released for it (Q = 5.494 MeV per reaction) must come from one rate. Otherwise the D-burning luminosity does not match the D actually consumed. The divisor is the nuclear mass in amu, the value in ENGEB's own isotope table.
- **What we did:** made one function, written token for token as ENGEB's expression. The energy is bitwise unchanged, and only the depletion rate moves (×0.99327).
  - Calibrated solar models at 4.568 Gyr: |Δlog L| ≤ 4.4e-5, |ΔX_c| ≤ 1.2e-4.
  - Pre-MS during D burning: Δlog L = +3.8e-4.
  - A related bug fixed before the divergence: both copies used T9^(+2/3) where CF88 has T9^(−2/3) (0304e269, 2026-09-02).

### N2. The S-factor derivative scales of reactions 2-7 were divided by the pp S₀ (REMAP)

- **Status:** fixed in [833b1a89](https://github.com/parallelpro/yrec/commit/833b1a89) (2026-10-10).
- **In the original F77?** Yes ([remap.f:52-53](https://github.com/parallelpro/yrec/blob/6cd56739/src/remap.f#L52-L53), [:126-192](https://github.com/parallelpro/yrec/blob/6cd56739/src/remap.f#L126-L192)):
  ```fortran
        DATA S0_1_1_A98,S0_3_3_A98,S0_3_4_A98,S0_1_12_A98,S0_1_13_A98/
       *     4.00D-22,5.4D3,5.3D-1,1.34D0,7.6D0/
              QS0E(2) = (S0P_3_3/S0_3_3)/(QS0E_3_3_A98/S0_1_1_A98)
              QS0E(7) = (S0P_1_16/S0_1_16)/(QS0E_1_16_A98/S0_1_1_A98)
              QS0E(8)=(S0P_1_BE7P/S0_1_BE7P)/(QS0E_1_BE7P_A98/S0_1_BE7P_A98)   ! correct
  ```
- **Now:** `set_first_derivative_scales` and `set_second_derivative_scales`, [user_inputs.f90:171-258](https://github.com/parallelpro/yrec/blob/e8549036/src/setup/user_inputs.f90#L171-L258).
- **What was wrong:** the rate fit is R = T9^(−2/3) + q1 T9^(−1/3) + s′(q2 + q3 T9^(1/3)) + s″(q4 T9^(2/3) + q5 T9). The scale s′ (s″) should be the deck's S′/S₀ (S″/S₀) relative to Solar Fusion I's. For reactions 2-7 the reference ratio was divided by the pp S₀ (4.00e-22 keV b) instead of the reaction's own. That made each scale too small by S₀(pp)/S₀(own), from 7e-26 (³He+³He) to 8e-22 (³He+⁴He). The q2…q5 terms of ³He+³He, ³He+⁴He, ¹²C+p, ¹³C+p, ¹⁴N+p and ¹⁶O+p were therefore effectively deleted, and the deck's S′ and S″ had no effect. pp and ⁷Be+p were right.
- **Expected behaviour, and why:** with S(E) = S₀ + S′E + S″E²/2, the effective S at the Gamow peak carries correction terms proportional to that reaction's S′/S₀ and S″/S₀ (Fowler, Caughlan & Zimmerman 1967, eq. 52; Bahcall 1989, eq. 3.14). Rescaling to a new compilation multiplies them by a number of order unity. With the deck defaults the scales become 1.24, 1.14, 1.00, 1.004, 0.54 and 2.0.
- **What we did:** each reaction is now normalized by its own S₀ (5.4e3, 0.53, 1.34, 7.6, 3.5 and 9.4 keV b, from the F77 DATA line).
  - Rates at log T = 7.1: ³He+³He −1.9%, ³He+⁴He −1.3%, ¹²C+p +5.9%, ¹⁴N+p −4.9%, ¹⁶O+p −14%. pp is unchanged.
  - Solar models at 4.568 Gyr: |Δlog L| ≤ 1.4e-4, |ΔX_c| ≤ 3.8e-4 (X_c higher, from less CN burning).

### N3. The rates' temperature derivatives omitted the S-factor scales (ENGEB)

- **Status:** fixed in [89c1fa42](https://github.com/parallelpro/yrec/commit/89c1fa42) (2026-10-10).
- **In the original F77?** Yes ([engeb.f:678-694](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L678-L694)):
  ```fortran
  C MHP 8/14 RATES CORRECTED TO PERMIT USER MODIFICATION OF REACTION
  C RATE DERIVATIVES
           R1=T9M23+Q1(I)*T9M13+QS0E(I)*(Q2(I)+Q3(I)*T9P13)+
       *      QQS0EE(I)*(Q4(I)*T9P23+Q5(I)*T9)
              DRATT(I)=DSCT(I)-(Q6(I)*T9M13+(2.*T9M23+Q1(I)*T9M13-Q3(I)*
       *      T9P13-2.*Q4(I)*T9P23-3.*Q5(I)*T9)/R1)/3.+2.*(Q8(I)*T9)**2
  ```
- **Now:** `compute_s_factor_rate`, [reaction_rates.f90:434-473](https://github.com/parallelpro/yrec/blob/e8549036/src/net/reaction_rates.f90#L434-L473).
- **What was wrong:** MHP 8/14 added the scales `QS0E` and `QQS0EE` to the rate but not to its derivative. In the derivative the q3 term lacks `QS0E(I)`, and the q4 and q5 terms lack `QQS0EE(I)`. `DRATT` feeds dln ε/dln T, which enters the Henyey matrix. In the F77 itself the scales were ~1e-22 (N2), so the rate had no q3…q5 terms while its derivative still carried them in full. After N2's fix the error is 5.4e-4 for ¹⁴N+p and 7.7e-3 for ¹⁶O+p.
- **Expected behaviour, and why:** the derivative must be that of the rate actually computed:
  dln r/dln T = dscr/dln T − (1/3)[q6 T9^(−1/3) + (2T9^(−2/3) + q1 T9^(−1/3) − s′q3 T9^(1/3) − s″(2q4 T9^(2/3) + 3q5 T9))/R] + 2(q8 T9)².
  A wrong Jacobian leaves the converged model unchanged but changes Newton's path. The same derivative also feeds the pulsation ε_T column.
- **What we did:** the scales are now carried in the derivative.
  - `test_net` checks every reaction's analytic dln r/dln T and dln r/dln ρ against centred differences; all are now at the floor set by the screening derivative.
  - At 4.568 Gyr: |Δlog L| ≤ 2.9e-5, |ΔX_c| ≤ 3.5e-5.

### N4. `LNEWNUC = .false.` left every pp and CNO rate at zero

- **Status:** fixed in [fff3deb5](https://github.com/parallelpro/yrec/commit/fff3deb5) (2026-10-10).
- **In the original F77?** No; the F77 path worked. `SSTANDARD` shared `COMMON /CROSS/` with the scales ENGEB multiplies the rates by ([parmin.f:369-370](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L369-L370), [remap.f:89-92](https://github.com/parallelpro/yrec/blob/6cd56739/src/remap.f#L89-L92), [engeb.f:685](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L685)):
  ```fortran
        DATA SSTANDARD/0.9828,1.0485,0.9815,0.9241,1.3818,1.0542,1.0,
        IF(LNEWNUC)THEN
           DO I = 1,17
              SSTANDARD(I) = 1.0D0
           RATE(I) = RATE(I)*SSTANDARD(I)
  ```
  The bug was introduced by [1439fb40](https://github.com/parallelpro/yrec/commit/1439fb40) (2026-08-20, "Replace common/cross/ and common/weak/ with const_lib module"). Its comment claimed that `sstandard`'s namelist value "is never referenced" because remap recomputes the scales. remap does that only when `LNEWNUC` is set.
- **Now:** [read_controls.f90:475-478](https://github.com/parallelpro/yrec/blob/e8549036/src/io/read_controls.f90#L475-L478).
- **What was wrong:** after 1439fb40 nothing copied `sstandard` into `cross_section_scale`. With `LNEWNUC = .false.`, the code default, the scales stayed at zero: every pp/CNO rate was zero and the ⁷Be branching ratio 0/0. All 509 inlists and 17 legacy decks set `LNEWNUC = .true.`, so no test caught it.
- **Expected behaviour, and why:** `LNEWNUC = .false.` scales the rates by the `SSTANDARD` ratios to Bahcall & Ulrich (1988). `.true.` derives them from the deck's S₀ inputs.
- **What we did:** with the flag off, `sstandard` is copied into `cross_section_scale` before the inputs are mapped. The new test `test_physics.py::test_legacy_cross_sections` runs the solar decks with `LNEWNUC = .FALSE.`. Shipped decks are byte-identical.

### N5. The ¹⁴N(α,γ)¹⁸O energy release subtracted the positron energy instead of the neutrino energy

- **Status:** fixed in [664c2fbc](https://github.com/parallelpro/yrec/commit/664c2fbc) (2026-10-11).
- **In the original F77?** Yes ([engeb.f:920-922](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L920-L922)):
  ```fortran
        DG(11)=EG(11)*5.815*CONVERT
  ```
- **Now:** [network.f90:206-210](https://github.com/parallelpro/yrec/blob/e8549036/src/net/network.f90#L206-L210), with the derivation at [:43-47](https://github.com/parallelpro/yrec/blob/e8549036/src/net/network.f90#L43-L47).
- **What was wrong:** the Q value of ¹⁴N(α,γ)¹⁸F(e⁺ν)¹⁸O from the masses is 6.071 MeV, so 5.815 MeV means 0.256 MeV is lost. That is the mean positron energy of ¹⁸F decay (0.250 MeV), not the mean neutrino energy. `test_net` compares every reaction's Q with AME2020 minus its neutrino loss. This was the only β⁺ reaction off by more than a few keV; ¹³N, ¹⁵O, ¹⁷F, ⁷Be and ⁸B all agree.
- **Expected behaviour, and why:**
  - From AME2020 mass excesses, Q = Δ(⁴He) + Δ(¹⁴N) − Δ(¹⁸O) = 6.071 MeV. Atomic masses already include the positron's annihilation, so only the neutrino energy is lost.
  - ¹⁸F decays 96.73% by β⁺ (⟨E_ν⟩ = 383.7 keV) and 3.27% by electron capture (E_ν = 1655.9 keV) (ENSDF; Tilley et al. 1995). So ⟨E_ν⟩ = 0.425 MeV, and the deposited energy is 6.071 − 0.425 = 5.646 MeV.
- **What we did:** set `q_mev = 5.646`. No pinned deck reaches helium burning, so all are byte-identical. A horizontal-branch run moves by |Δlog L| ≤ 4.2e-8, because its seed has already burned its ¹⁴N.

### N6. Li/Be burning never saw the end-of-step model (LIRATE88 flags swapped)

- **Status:** fixed in [0011e3a8](https://github.com/parallelpro/yrec/commit/0011e3a8) (2026-10-10). The flag was then replaced by named entry points in e6caa171.
- **In the original F77?** Yes. MAIN calls LIRATE88 with `J = 1` at the start of the step and `J = 2` at the end ([main.f:633-634](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L633-L634), [:842-843](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L842-L843)). GETW also calls `J = 2` ([getw.f:324-325](https://github.com/parallelpro/yrec/blob/6cd56739/src/getw.f#L324-L325)). Inside LIRATE88, however, `J = 1` evaluates the model passed in and stores the end-of-step slot, while `J = 2` evaluates the saved start-of-step model ([lirate88.f:31-37](https://github.com/parallelpro/yrec/blob/6cd56739/src/lirate88.f#L31-L37)):
  ```fortran
  C find burning rates at the beginning of the time step.
                 CALL LIRATE88(HCOMP,HD,HT,M,1)                 ! main.f:634
  C FIND BURNING RATES AT THE END OF THE TIME STEP.
                 CALL LIRATE88(HCOMP,HD,HT,M,2)                 ! main.f:843
           IF(J.EQ.1)THEN
              RHOX = EXP(CLN*HD(I))*HCOMP(1,I)                  ! lirate88.f: J=1 -> RLI6
           ELSE
              RHOX = EXP(CLN*HDO(I))*HCOMPP(1,I)                !             J=2 -> RLI60
  ```
- **Now:** `store_start_li_be_rates` and `store_end_li_be_rates`, [burn_step.f90:378-462](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/burn_step.f90#L378-L462).
- **What was wrong:** at the start of the step, `J = 1` filled the end-of-step slot from the start-of-step model. At the end, `J = 2` recomputed the start slot from the start-of-step model again. In GETW, `J = 2` ignored the midpoint model it was passed. So LIBURN interpolated between two copies of the start-of-step rates.
- **Expected behaviour, and why:** LIBURN integrates dX/dt = −λ(t)X with ln λ interpolated linearly in time between the start and end rates (liburn.f:155-163), so it needs the end-of-step rate. ⁷Li(p,α) scales roughly as T²⁰ near 2.5 MK. A contracting pre-MS star's convection-zone base is hotter at the end of the step, so start-only rates under-deplete Li.
- **What we did:** swapped the flags at the three calls. The structure is identical, since these are trace species. Surface ⁷Li at 4.568 Gyr: 4.52e-9 → 3.70e-9 (−0.088 dex) for the non-rotating solar deck, and 1.02e-11 → 7.83e-12 (−0.113 dex) for the rotating one. ⁹Be moves by less than 0.001 dex.

### N7. The rotation sub-step burned with stale end-of-step rates (ROTMIX)

- **Status:** fixed in [4da5c0b4](https://github.com/parallelpro/yrec/commit/4da5c0b4) (2026-10-10).
- **In the original F77?** Yes. MIDMOD calls ROTMIX with the midpoint temperature ([midmod.f:268-273](https://github.com/parallelpro/yrec/blob/6cd56739/src/midmod.f#L268-L273)), but ROTMIX takes its rates from `COMMON /BURN/` (rotmix.f:43-59). MIX fills that block during the last CRRECT iteration ([mix.f:230-246](https://github.com/parallelpro/yrec/blob/6cd56739/src/mix.f#L230-L246)):
  ```fortran
  C MHP 6/00 - FOR THERMODYNAMIC CONSISTENCY USE CURRENT T, NOT END OF
  C TIMESTEP T (REPLACED HT IN CALL WITH HTM)
        CALL ROTMIX(DT,HCOMP,HS2,HTM,5,M,MRZONE,MXZONE,
           HR1(II) = HCOMPM(1,II)                                 ! rotmix.f
           HCOMPM(1,I) = HR1(I)                                   ! mix.f
  ```
- **Now:** [mid_timestep_model.f90:102-112](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/mid_timestep_model.f90#L102-L112) calls `network_rates_of_model` ([burn_step.f90:616-647](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/burn_step.f90#L616-L647)).
- **What was wrong:** the rates came from the end-of-step structure of the previous Newton iterate. The midpoint temperature passed in was used only for the cutoff test, and HPOINT regridded only 7 of the 15 rate rows on a rezone. Deuterium and Li/Be in the same sub-step already used the sub-step's own model.
- **Expected behaviour, and why:** in an operator-split sub-step the network burns with the rates of the sub-step's own model, which is what MHP's 6/00 comment intends. pp and CNO rates are steep in T (about T⁴ and T¹⁸), so end-of-step rates applied at the midpoint bias the burn at first order in Δt.
- **What we did:** the rates are now computed on the midpoint model and passed as an argument; the COMMON relay and its partial regridding are gone. On the rotating solar deck (the only shipped deck with rotational mixing): |Δlog L| 8.5e-5, |ΔX_c| 2.2e-4, surface Li +1e-4 dex.

### N8. DBURNM added the accreted mass in grams to a dimensionless 1

- **Status:** fixed in [29048ffd](https://github.com/parallelpro/yrec/commit/29048ffd) (2026-10-10). DBURN and DBURNM were merged in c81d8ba5.
- **In the original F77?** Yes ([dburnm.f:98-103](https://github.com/parallelpro/yrec/blob/6cd56739/src/dburnm.f#L98-L103)). MDOT stores `FMASSACC` in grams ([mdot.f:202](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L202)):
  ```fortran
           FMASSA = FAC2*FMASSACC
           XH2N = (XH2NO+XH2NA*FMASSA)/(1.0D0+FMASSA)
        FMASSACC = DELM                                    ! mdot.f
  ```
- **Now:** `burn_deuterium_region`, [light_elements.f90:257-289](https://github.com/parallelpro/yrec/blob/e8549036/src/net/light_elements.f90#L257-L289).
- **What was wrong:** the comments define `FMASSACC` as ΔM/M_cz, a dimensionless number (dburnm.f:74-75), and DBURNM's blend is written for that. But MDOT stores ΔM in grams. DBURN and LIBURN correctly use it as a mass: (X M + X_acc ΔM)/(M + ΔM). In DBURNM, `1 + FAC2*ΔM` is dominated by grams, so any accretion replaced the surface region's D, H and ³He by the accreted values. This is active only in rotating runs with accretion.
- **Expected behaviour, and why:** each species' mass is conserved when the region takes in the sub-step's share of the accreted gas: X_new = (X M + X_acc f ΔM)/(M + f ΔM).
- **What we did:** the blend is now mass-weighted, as DBURN's is. On a scratch accreting deck (1e-8 Msun/yr of gas with D = 5e-5), surface D goes from 5.0e-5 (the accreted value) to 1.5e-5 (the star's own burning). No pinned deck accretes.

### N9. DBURN's single-shell accretion branch used a region mass it never set

- **Status:** fixed in [29048ffd](https://github.com/parallelpro/yrec/commit/29048ffd) (2026-10-10).
- **In the original F77?** Yes, though under `SAVE` the variable held a stale value rather than zero ([dburn.f:16-22](https://github.com/parallelpro/yrec/blob/6cd56739/src/dburn.f#L16-L22), [:99-101](https://github.com/parallelpro/yrec/blob/6cd56739/src/dburn.f#L99-L101)):
  ```fortran
        IF(IBEGIN.EQ.IEND)THEN
           X = HCOMPP(1,IBEGIN)
  ...
           XH2N = (XH2NO*SUMM+XH2NA*FMASSACC)/(SUMM+FMASSACC)
  ```
  `SUMM` is set only in the multi-shell branch, so the single-shell branch used the previous call's region mass. It became zero (through `-finit-local-zero`) when [1c510cc5](https://github.com/parallelpro/yrec/commit/1c510cc5) (2026-08-22) removed the blanket `SAVE`.
- **Now:** [light_elements.f90:226-234](https://github.com/parallelpro/yrec/blob/e8549036/src/net/light_elements.f90#L226-L234).
- **What was wrong:** for a single radiative shell under accretion, the blend weighted the shell's own D by a stale or zero mass. The shell took the accreted deuterium and dropped its own.
- **Expected behaviour, and why:** a single shell's region mass is its own mass, so the blend is (X m + X_acc ΔM)/(m + ΔM).
- **What we did:** `region_mass` is now set in both branches. No pinned deck reaches this branch.

### N10. The mathematical constants were zero until SETUPS had run

- **Status:** fixed in [529d3fb3](https://github.com/parallelpro/yrec/commit/529d3fb3) (2026-10-11).
- **In the original F77?** Yes, the same pattern, though harmless there. `CLN`, `CC13`, `CPI` and the rest live in `COMMON /CONST1/` and `/CONST3/` and are assigned in SETUPS, with no BLOCK DATA ([setups.f:75-95](https://github.com/parallelpro/yrec/blob/6cd56739/src/setups.f#L75-L95)):
  ```fortran
        CLNDP = DLOG(10.0D0)
        CLN = CLNDP
        CC13 = 1.0D0/3.0D0
        CPI = 3.1415926535898D0
  ```
  MAIN calls SETUPS before any physics, which hid the problem.
- **Now:** [physical_constants.f90:18-23](https://github.com/parallelpro/yrec/blob/e8549036/src/state/physical_constants.f90#L18-L23).
- **What was wrong:** once the network became a library, any caller that had not run the setup computed with zeros: exp(ln10 · log ρ) = 1 and T9^cc13 = 1. Test harnesses had to set the constants by hand.
- **Expected behaviour, and why:** mathematical constants do not depend on the run, so they should be compile-time parameters, valid from the first call.
- **What we did:** all ten are now `parameter`s, each bit-identical to the run-time value; byte-identical.

## 3. Structure, timestep and boundary conditions (fixed; results change)

Eleven of these thirteen are in the original F77 code. S5 comes from the 2025 ENTIME rewrite, and S13 from the 2026 modernization.

### S1. The Ω(j) inversion inside the Newton iteration never got its convective flags

- **Status:** fixed in [2ba4e32b](https://github.com/parallelpro/yrec/commit/2ba4e32b) (2026-10-10).
- **In the original F77?** Yes. `LCZ` is a local array of the main program. It is passed only to CRRECT, and nothing in `main.f` sets it ([main.f:188](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L188), [main.f:687-693](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L687-L693)). CRRECT hands it to GETROT ([crrect.f:238-239](https://github.com/parallelpro/yrec/blob/6cd56739/src/crrect.f#L238-L239)), whose only flag-setting call, OVROT, is commented out ([crrect.f:107-110](https://github.com/parallelpro/yrec/blob/6cd56739/src/crrect.f#L107-L110)). GETROT then treats every shell as radiative ([getrot.f:16-24](https://github.com/parallelpro/yrec/blob/6cd56739/src/getrot.f#L16-L24)):
  ```fortran
  C main.f
        DIMENSION LCZ(JSON)
       *           HG,HI,HJM, ITLVL,LCZ,MRZONE,MXZONE,NRZONE,NZONE,
  C crrect.f
  c      IF(LROT)THEN
  c         CALL OVROT(HCOMP,HD,HP,HR,HS,HT,LC,M,LCZ,MRZONE,MXZONE,
           IF(LROT)THEN
              CALL GETROT(HD,HJM,HR,HS,HS2,LCZ,M,ETA2,HI,OMEGA,QIW,R0)
  C getrot.f
              IF(.NOT.LCZ(I)) THEN
  C  RADIATIVE SHELL
                 JSTART = I
                 JEND = I
                 CALL SOLID(HD,HJM,HR,HS,HS2,JSTART,JEND,ETA2,HI,OMEGA,
  ```
  By contrast, GETW calls OVROT before GETROT ([getw.f:102, 119](https://github.com/parallelpro/yrec/blob/6cd56739/src/getw.f#L102-L119)). HPOINT computes the flags into a local and copies them into `LCZO`, which nothing reads ([hpoint.f:697-701](https://github.com/parallelpro/yrec/blob/6cd56739/src/hpoint.f#L697-L701), [:729](https://github.com/parallelpro/yrec/blob/6cd56739/src/hpoint.f#L729)). STARIN computes `LCZ` and discards it ([starin.f:902](https://github.com/parallelpro/yrec/blob/6cd56739/src/starin.f#L902)).
- **Now:** the flags are written by the rezone ([rezoning.f90:928-945](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/rezoning.f90#L928-L945)) and by the starting model ([starting_model.f90:220-225](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/starting_model.f90#L220-L225)). They are read in the Newton loop ([henyey_iteration.f90:425-432](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/henyey_iteration.f90#L425-L432)).
- **What was wrong:** with rotation on, CRRECT inverts j → Ω (GETROT) after every Newton correction. The flags it passed were never assigned, so they were all false (static storage is zero). Every shell was therefore inverted on its own (SOLID), and solid-body rotation (WCZ) was never imposed over the zones that are convective for angular-momentum transport. The 2026 tree kept the same gap. The rezone stored the right flags only in a write-only `cz_flag_start`, and `read_starting_model` discarded them.
- **Expected behaviour, and why:** zones that are convective for angular-momentum purposes (convection plus the OVROT overshoot) are coupled by turbulence. Inside them Ω follows the convection-zone law (solid body for `walpcz >= 0`), as GETW and HPOINT enforce. The structure solve must use the same Ω profile as the angular-momentum step, because that profile sets f_P, f_T and the rotational energy term.
- **What we did:** the rezone and `read_starting_model` now write `star%am_transport_convective_flag`, and `cz_flag_start` is gone. He-flash steps, which skip the rezone, keep the last flags on the same grid.
  - Only rotating runs change: 9 legacy logs and the GS_rot MESA baseline. The 26 non-rotating logs are byte-identical.
  - Calibrated solar models at 4.568 Gyr: |Δlog L| ≤ 3e-5, |Δlog R| ≤ 1.2e-5, |ΔX_c| ≤ 7e-5; log Teff is unchanged to 6 digits.
  - Along the tracks the largest change is on the fully convective pre-MS: |Δlog L| ≤ 1.9e-3, |Δlog R| ≤ 9e-4.

### S2. `update_total_mass` read m(0) when the convection zone reached the centre

- **Status:** fixed in [bc874267](https://github.com/parallelpro/yrec/commit/bc874267) (2026-10-12).
- **In the original F77?** Yes. In MDOT the shell masses are recomputed from `JENV` ([mdot.f:183-186](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L183-L186)). `JENV = 1` is a legal case that MDOT handles explicitly elsewhere ([mdot.f:52-56](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L52-L56)):
  ```fortran
  C CORRECT MASS CONTENTS OF INDIVIDUAL SHELLS (HS2, IN GM)
        DO I = JENV,M-1
           HS2(I) = 0.5D0*(HS1(I+1)-HS1(I-1))
        END DO
  ```
  HPOINT reflects the first point about the centre instead (`DS2 = - DS3`, [hpoint.f:680-687](https://github.com/parallelpro/yrec/blob/6cd56739/src/hpoint.f#L680-L687)).
- **Now:** [mass_change_kernel.f90:70-98](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change_kernel.f90#L70-L98); the fix is at lines 93-96.
- **What was wrong:** after mass loss or accretion, the convection zone's shell masses are set to dm(k) = (m(k+1) − m(k−1))/2 from k = jenv. In a fully convective star (jenv = 1) this reads `HS1(0)`, one word before the array, on every mass-changing step.
- **Expected behaviour, and why:** the innermost shell runs from the centre to the midpoint of m(1) and m(2), so dm(1) = (m(1) + m(2))/2. That is the m(0) = −m(1) convention that HPOINT and the rezone use. Otherwise Σdm no longer equals the total mass, and everything summed over dm is off: the zone's J, its mass and its mixed abundances.
- **What we did:** the loop starts at max(jenv, 2), and dm(1) = (m(1) + m(2))/2 when jenv = 1.
  - Over 12 fully convective pre-MS models: cm_massloss_const |Δlog L| ≤ 1e-6; cm_massloss_reimers |Δlog L| ≤ 2e-6.
  - The other 19 check decks are byte-identical.
  - `test_mass_change` fails if the fix is reverted.

### S3. The ε_grav partials lacked a 1/χ_T factor

- **Status:** fixed in [18869257](https://github.com/parallelpro/yrec/commit/18869257) (2026-10-13).
- **In the original F77?** Yes ([coefft.f:265-276](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L265-L276)):
  ```fortran
              ENTR1 = P*QDT/D
              ENTR2 = ENTR1/DELA
              ENTR = (ENTR2*HTTT - ENTR1*HPPP)*CLN
              ...
              QLP = QLP + DT*(ENTR*(1.0D0-QDP+QDTP)-ENTR1 - ENTR3*QAP)
              QLT = QLT + DT*(ENTR*(-QDT+QDTT) + ENTR2 - ENTR3*QAT)
  ```
  Elsewhere the F77 itself divides by `QDT`, for example `QAT = -1.0D0 - QDT - QCPT + QDTT/QDT` in EQSAHA ([eqsaha.f:301-302](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqsaha.f#L301-L302)). So `QDTT` and `QDTP` are dχ_T/dln T and dχ_T/dln P, not dln χ_T/dln T and dln χ_T/dln P.
- **Now:** [shell_equations.f90:298-315](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/shell_equations.f90#L298-L315).
- **What was wrong:** ε_grav = (P χ_T/ρ)(Δln T/∇_ad − Δln P)/Δt contains the factor χ_T = ∂ln ρ/∂ln T, so its log-derivatives need dln χ_T = dχ_T/χ_T. The partials used dχ_T without dividing by χ_T. The residual itself was right; only the Jacobian was wrong.
- **Expected behaviour, and why:** the Henyey Jacobian must hold the true partials of the residual. With a wrong Jacobian Newton still converges to the same solution, but along a different path. The accepted model then differs at the level of the solver tolerance, and the difference propagates through the timestep sequence.
- **What we did:** the partials now use `dlnrho_dlnt_dp/dlnrho_dlnt` and `dlnrho_dlnt_dt/dlnrho_dlnt`. Every deck moves at the solver-tolerance level:
  - Config matrix: |Δlog L| ≤ 3e-6.
  - Solar decks: 1e-6 (noGS) and 1.2e-5 (GS rot).
  - Rotating-mixing and neutrino-table decks: up to 1.5-1.8e-4 near the ZAMS. That phase arrives about 1e4 yr earlier; the final models agree within 2e-5 in log L.
  - The total Newton iteration count over the 21 check decks is unchanged (1042).
  - `test_structure`, a finite-difference check of every partial, fails if the factor is removed.

### S4. The end-age timestep floor was 10⁻³ of the age read in years instead of Gyr

- **Status:** fixed in [38610071](https://github.com/parallelpro/yrec/commit/38610071) (2026-10-13).
- **In the original F77?** Yes ([main.f:869-872](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L869-L872)). `DAGE` is in Gyr (`DAGE = DAGE + 1.0D-9*DDAGE`, [main.f:573](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L573)):
  ```fortran
           IF(LENDAG(NK) .AND. ENDAGE(NK).GT.0.0D0)THEN
              IF(ENDAGE(NK)-DAGE*1.0D9.LE.1.0D0)THEN
                 DELTS = MAX(DELTS0,1.0D-3*DAGE*CSECYR)
  ```
- **Now:** [step_driver.f90:197-212](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/step_driver.f90#L197-L212).
- **What was wrong:** at a fixed end age, the clipped (near-zero) last timestep is replaced in the saved model by the previous step, floored at a fraction of the age. `1.0D-3*DAGE*CSECYR` with `DAGE` in Gyr is 0.005 yr at 4.57 Gyr, so the floor never acted.
- **Expected behaviour, and why:** the floor gives a sequence restarted from the final model a sensible first timestep. It is meant to be 10⁻³ of the age, 10⁻³ × DAGE × 10⁹ yr.
- **What we did:** the floor is now `1.0d-3*s%age_gyr*1.0d9*seconds_per_year`. Only the final model's stored timestep changes: GS rot solar 4.06e6 → 4.568e6 yr, difrotmixper18p4d0p1 8.43e5 → 4.568e6 yr. Everything else is byte-identical.

### S5. The HR-motion timestep limit read the start-of-step L at the wrong mesh point

- **Status:** fixed in [eee659ce](https://github.com/parallelpro/yrec/commit/eee659ce) (2026-10-13).
- **In the original F77?** No, not in this form. It was introduced by the YREC team's 2025 ENTIME rewrite [979a9c26](https://github.com/parallelpro/yrec/commit/979a9c26) (2025-05-12; its call in HTIMER was made consistent in [956a5ede](https://github.com/parallelpro/yrec/commit/956a5ede)). HPOINT saves the old point count `MO` before the rezone ([hpoint.f:276-277 at 979a9c26](https://github.com/parallelpro/yrec/blob/979a9c26/src/hpoint.f#L276-L277)) and then re-splines `HLO` onto the new points. ENTIME reads `HLO(MO)` ([entime.f:52 at 979a9c26](https://github.com/parallelpro/yrec/blob/979a9c26/src/entime.f#L52)):
  ```fortran
  C  JVS 05/25 Added model number to list of saved values
        MO = M
  ...
        TESTL = ABS(DLOG10(HLO(MO))-DLOG10(HL(M)))
  ```
  The ENTIME at 6cd56739 had a different flaw on the same read: it took `MO` as the last non-zero `HLO` and compared linear L ([entime.f:60-80](https://github.com/parallelpro/yrec/blob/6cd56739/src/entime.f#L60-L80)).
- **Now:** the limit is at [timestep_control.f90:621-630](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/timestep_control.f90#L621-L630) and the save at [rezoning.f90:431-434](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/rezoning.f90#L431-L434).
- **What was wrong:** the limit compared the surface log L with `luminosity_lsun_start(nz_start)`, where `nz_start` is the point count before the rezone. But the rezone had re-splined `luminosity_lsun_start` onto the new mesh. Whenever the count changed, the read hit an interior point (more points) or a stale entry past the new count (fewer points).
- **Expected behaviour, and why:** the limit keeps a step's motion in (log Teff, log L) within the envelope triangle's size, so that the triangle interpolation of the surface boundary condition stays valid. That needs the previous model's surface L, saved before the mesh changes, as log Teff already is.
- **What we did:** a new field, `star%surface_luminosity_lsun_start`, replaces `nz_start`, which nothing else read.
  - In Test_m1p4feh+0p0GS98_dif_TAMS the point count changed in 430 of 558 limited steps; the value read was up to 1.9e-2 dex off.
  - That deck moves by |Δlog L| ≤ 2.7e-3 while X_c > 0.05, and core-H exhaustion comes 7.6 Myr (0.25%) earlier. This is within the deck's own timestep sensitivity: changing tri_delta_logl from 0.040 to 0.0399 alone moves it by 3e-2.
  - The other four testsuite decks that use the limit, and all 21 check decks, are byte-identical.

### S6. The atmosphere integration's first slope lacked f_P

- **Status:** fixed in [e2c526aa](https://github.com/parallelpro/yrec/commit/e2c526aa) (2026-10-13).
- **In the original F77?** Yes. ENVINT sets the starting derivative with the plain g ([envint.f:200-202](https://github.com/parallelpro/yrec/blob/6cd56739/src/envint.f#L200-L202)), while the derivative routine QATM uses g·f_P ([qatm.f:44, 67](https://github.com/parallelpro/yrec/blob/6cd56739/src/qatm.f#L44-L67)):
  ```fortran
  C envint.f
        X0 = PL - GL + DLOG10(O)
        DYDX(1) = DEXP(CLN*(GL+X0-OL-PL))
  C qatm.f
        G = DEXP(CLN*GL)*FPL
        DYDX(1) = G*TAU/(P*O)
  ```
- **Now:** [envelope_integration_kernel.f90:446-450](https://github.com/parallelpro/yrec/blob/e8549036/src/core/envelope/envelope_integration_kernel.f90#L446-L450).
- **What was wrong:** every atmosphere step after the first takes dlog P/dlog τ = g_eff τ/(κP) from QATM, with g_eff = g f_P. The start derivative handed to the first BSSTEP used plain g, so it equals 1 rather than f_P. A rotating star's first Bulirsch-Stoer step therefore began from the wrong slope.
- **Expected behaviour, and why:** with rotation, hydrostatic equilibrium in the atmosphere is dP/dτ = g_eff/κ, so at the starting point dlog P/dlog τ = f_P. The modified-midpoint step uses the supplied start derivative directly, so an inconsistent value biases that step.
- **What we did:** the start derivative is now multiplied by `pressure_rotation_factor`. Only rotating decks move: GS rot solar |Δlog L| ≤ 2.4e-6; rotating-mixing decks ≤ 8.3e-5 (1.9e-4 for the solid-body one near the ZAMS); the final models agree within 2.4e-5 in log L.

### S7. Integrating a triangle vertex overwrote the model's log Teff

- **Status:** fixed in [3d065901](https://github.com/parallelpro/yrec/commit/3d065901) (2026-10-13).
- **In the original F77?** Yes. SURFBC integrates each new vertex by writing its Teff into its `TEFFL` argument, which is the model's log Teff ([surfbc.f:129-136](https://github.com/parallelpro/yrec/blob/6cd56739/src/surfbc.f#L129-L136)). CRRECT recomputes `TEFFL` only after the corrections ([crrect.f:90-123](https://github.com/parallelpro/yrec/blob/6cd56739/src/crrect.f#L90-L123)), and TPGRAD's spot contrast at depth reads it in between ([tpgrad.f:90](https://github.com/parallelpro/yrec/blob/6cd56739/src/tpgrad.f#L90)):
  ```fortran
  C surfbc.f
           IF(.NOT.LTRI(I)) THEN
  C NEW ENVELOPE NEEDED
              LFLP = .TRUE.
              TEFFL = TRIT(I)
  C tpgrad.f
              ATEFFL = TEFFL - 0.25*LOG10(SPOTF * SPOTX**4.0 + 1.0 - SPOTF)
  ```
- **Now:** [surface_boundary.f90:224-245](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/surface_boundary.f90#L224-L245); `log10_teff` is `intent(in)`.
- **What was wrong:** after a new or reset triangle, the model's log Teff was the last integrated vertex's until CRRECT recomputed it. The first Newton pass's depth-dependent spot contrast (`LSDEPTH`) used that value.
- **Expected behaviour, and why:** the vertex integrations are auxiliary envelopes at the triangle's corners and must not change the model's state. The spot contrast should depend on the model's own Teff, not on whichever corner happened to be integrated last.
- **What we did:** each vertex now uses a local `vertex_log10_teff`. Only spot decks with depth-dependent spots change: 1.0 Msun |Δlog L| ≤ 3.2e-5, 0.3 Msun 1.1e-7, 0.1 Msun 3e-9. The 21 check decks are byte-identical.

### S8. The wind removed the convection zone's angular momentum twice

- **Status:** fixed in [91a6b567](https://github.com/parallelpro/yrec/commit/91a6b567) (2026-10-12).
- **In the original F77?** Yes. MDOT first scales j in the convection zone by (J_cz + ΔJ)/J_cz ([mdot.f:146-157](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L146-L157)). It then shrinks the zone's masses at that j ([mdot.f:162-186](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L162-L186)):
  ```fortran
           CJNEW = (CJTOT+DELJ)
           RATIO = CJNEW/CJTOT
           DO J = JENV,M
              HJM(J) = HJM(J)*RATIO
           END DO
        ...
              HS1(I) = HS1(JENV)+FX*(HS1(I)-HS1(JENV))
  ```
- **Now:** the kernel is at [mass_change_kernel.f90:133-179](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change_kernel.f90#L133-L179), and the driver's order of operations at [mass_change.f90:524-568](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change.f90#L524-L568).
- **What was wrong:** with ΔJ = j_s ΔM, shrinking the masses after the j scaling removed a further Σ j Δ(dm) ≈ ⟨j⟩|ΔM|. The zone lost ΔJ + ⟨j⟩|ΔM|.
- **Expected behaviour, and why:** the lost matter carries away only its own specific angular momentum j_s. The remaining zone must hold exactly J_cz + ΔJ on its new masses.
- **What we did:** J_cz is now measured before the masses change, and j is scaled after them by (J_cz + ΔJ)/Σ(j dm_new).
  - Scratch deck (cm_rot_base at 1e-9 Msun/yr, fully convective): 1.08e43 extra was removed per step, against ΔJ = −3.80e43. Now exactly ΔJ is removed.
  - No shipped deck both rotates and loses mass; all are byte-identical.
  - `test_mass_change` fails if the fix is reverted.

### S9. The wind and the accretion first guess extrapolated Ω ∝ r^(−walpcz)

- **Status:** fixed in [d68c375e](https://github.com/parallelpro/yrec/commit/d68c375e) (2026-10-12).
- **In the original F77?** Yes, in MDOT ([mdot.f:143](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L136-L145)) and GETNEWENV ([getnewenv.f:252-254](https://github.com/parallelpro/yrec/blob/6cd56739/src/getnewenv.f#L250-L257)). Both disagree with the law that WCZ imposes ([wcz.f:38-53](https://github.com/parallelpro/yrec/blob/6cd56739/src/wcz.f#L38-L53)):
  ```fortran
  C mdot.f
              WSURF = OMEGA(M)*10.0D0**(HR(M)*WALPCZ)/RTOT**WALPCZ
  C getnewenv.f
              W0 = OMEGA(MM)*10.0D0**(HR(MM)*WALPCZ)
                 OMEGA(J) = W0/10.0D0**(HR(J)*WALPCZ)
  C wcz.f
  C  GENERAL LAW FOR OMEGA IN C.Z.: OMEGA = C*R**WALPCZ,WHERE C IS A CONSTANT
              OMEGA(J) = W0*DEXP(CLN*WALPCZ*HR(J))
  ```
- **Now:** the wind is at [mass_change_kernel.f90:103-131](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change_kernel.f90#L103-L131), and the accretion guess at [envelope_rebuild.f90:267-297](https://github.com/parallelpro/yrec/blob/e8549036/src/core/envelope/envelope_rebuild.f90#L267-L297).
- **What was wrong:** the two places that extrapolate the convection-zone law beyond the last shell used Ω r^walpcz = const. These are the wind's surface Ω and the first guess for shells added by accretion. Only −2 < walpcz < 0 is affected, because walpcz ≤ −2 means constant j and walpcz ≥ 0 means solid body.
- **Expected behaviour, and why:** the extrapolation must follow the law WCZ imposes inside the zone. For walpcz = −1 the surface should spin slower than the fitting point, not faster.
- **What we did:** both now use Ω = Ω_old (r/r_old)^walpcz. In a scratch deck (walpcz = −1, 1e-9 Msun/yr) the first step's ΔJ went from 1.0032 to 0.9969 times the solid-body value. All shipped decks have walpcz = 0 and are byte-identical.

### S10. The Z ramp scaled CNO by the fraction of Z removed

- **Status:** fixed in [1bcdee96](https://github.com/parallelpro/yrec/commit/1bcdee96) (2026-10-12).
- **In the original F77?** Yes ([rscale.f:342-356](https://github.com/parallelpro/yrec/blob/6cd56739/src/rscale.f#L342-L356)). Compare the uniform rescale, which uses Z_new/Z_old ([rscale.f:74](https://github.com/parallelpro/yrec/blob/6cd56739/src/rscale.f#L74-L78)):
  ```fortran
                 IF (FS .LT. RSCLZM1(NK)) THEN
                    ZFACT = (HCOMP(3,I)-RSCLZC(NK))/HCOMP(3,I)
                 ...
                    ZFACT = (HCOMP(3,I)-ZP)/HCOMP(3,I)
  C uniform rescale:
                 ZFACT=DMAX1(0.D0,HCOMP(3,J)+DELZ)/(HCOMP(3,J)+1.D-30)
  ```
- **Now:** [model_rescaling.f90:404-439](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/model_rescaling.f90#L404-L439).
- **What was wrong:** the `LZRAMP` rescale set the interior Z to the ramp value, but multiplied C12…O18 by (Z_old − Z_new)/Z_old. That is right only when Z is halved. It gives zero when Z is unchanged and a negative value when Z rises.
- **Expected behaviour, and why:** CNO is part of Z. Changing Z at a fixed mixture scales each isotope by Z_new/Z_old, so CNO/Z keeps the mixture's value.
- **What we did:** the factor is now Z_new/(Z_old + 1e-30). In a scratch deck (Z = 0.0015 inside q = 0.2, ramping to 0.002 at q = 0.5), CNO/Z was 0.230 at the centre and 0.099 at q = 0.35, against the surface's 0.691; it is now 0.691 throughout. No shipped deck uses the ramp.

### S11. A reused model was rescaled twice, and the relative species target overwrote its control

- **Status:** fixed in [df7ca6df](https://github.com/parallelpro/yrec/commit/df7ca6df) (2026-10-12).
- **In the original F77?** Yes. On a card that reuses the model in memory, STARIN calls RSCALE, jumps to 3000 and calls it again ([starin.f:153-159](https://github.com/parallelpro/yrec/blob/6cd56739/src/starin.f#L153-L159), [:458-460](https://github.com/parallelpro/yrec/blob/6cd56739/src/starin.f#L458-L460)). RSCALE converts the relative target in place, in `COMMON /NEWCMP/` ([rscale.f:107](https://github.com/parallelpro/yrec/blob/6cd56739/src/rscale.f#L101-L110)):
  ```fortran
        IF(.NOT.LFIRST(NK)) THEN
           IF(IRESCA(NK).NE.1) CALL RSCALE(HL,HCOMP,HS,HSTOT,M,NK,
       *      SMASS,LC)
           GOTO 3000
  ...
        IF(IRESCA(NK).NE.1) CALL RSCALE(HL,HCOMP,HS,HSTOT,M,NK,
  C rscale.f
              IF(LREL) XNEWCP = DEXP(CLN*(XNEWCP-12.0D0))*HCOMP(1,M)
  ```
- **Now:** [starting_model.f90:254-260](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/starting_model.f90#L254-L260) and [model_rescaling.f90:176-201](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/model_rescaling.f90#L176-L201).
- **What was wrong:** a model read from a file was rescaled once, but a reused one twice. For the absolute targets (X, Z, mass) the second pass only repeats the first, up to round-off. The relative species target (`LREL`), however, overwrote its control `XNEWCP` with the converted mass fraction, and the second pass converted that again: 10^(~0 − 12) X ≈ 1e-12 X.
- **Expected behaviour, and why:** a rescale card is applied once, whatever the model's source. A user control must not be overwritten by its derived value, so that every card converts the requested log abundance afresh.
- **What we did:** the reused-model call is removed, and the conversion goes into a local `species_target`.
  - A card-2 Li7 rescale to 3.3 gave Li7/X = 1.0e-12; it now gives 1.9953e-9 = 10^(3.3−12).
  - Two TAHB decks that rescale a reused model change below the printed digits; the 21 check decks are byte-identical.

### S12. Switching back from the gray atmosphere to a table did not start a new triangle

- **Status:** fixed in [67d711a7](https://github.com/parallelpro/yrec/commit/67d711a7) (2026-10-12).
- **In the original F77?** Yes. In SURFBC the switch to gray sets `LNEW`, but the switch back does not ([surfbc.f:42-82](https://github.com/parallelpro/yrec/blob/6cd56739/src/surfbc.f#L42-L82)):
  ```fortran
              KTTAU = 0
              LNEW = .TRUE.
        ...
        IF(LTTAU)THEN
           IF(KTTAU0.EQ.3.AND.TEFFL.LT.3.95D0)THEN
              KTTAU = KTTAU0
              LTTAU = .FALSE.
  ```
- **Now:** [surface_boundary.f90:116-167](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/surface_boundary.f90#L116-L167).
- **What was wrong:** Kurucz (`KTTAU` = 3) and Allard (4) atmospheres fall back to gray above their Teff range and return once Teff is back in range. On the return, the old triangle, whose three vertex envelopes had been integrated with the gray T(τ), was kept. The surface boundary condition stayed partly gray until each vertex was reflected out.
- **Expected behaviour, and why:** the triangle interpolates the fitting-point (P, T, R) between three envelope integrations. All three must use the atmosphere now in force, so either switch must rebuild them.
- **What we did:** both switches now set `start_new_triangle`. In a scratch deck (Kurucz, 2.4 Msun pre-MS; gray from model ~810, back at ~1460) the track differs from model 1500 on, by |Δlog Teff| ≤ 1.9e-5 at model 1750. No shipped deck switches.

### S13. The shallower envelope refit used L_sun in erg/s where log L_sun is needed

- **Status:** fixed in [4f34e95f](https://github.com/parallelpro/yrec/commit/4f34e95f) (2026-10-12).
- **In the original F77?** No. The F77 is correct ([starin.f:582](https://github.com/parallelpro/yrec/blob/6cd56739/src/starin.f#L581-L583)). The bug was introduced by [6e31b68d](https://github.com/parallelpro/yrec/commit/6e31b68d) (2026-08-20, "Modernize starin.f"; [starin.f90:1122](https://github.com/parallelpro/yrec/blob/6e31b68d/src/starin.f90#L1121-L1123)):
  ```fortran
  C 6cd56739
              RL = 0.5D0*(BL + CLSUNL - 4.0D0*TEFFL - C4PIL - CSIGL)
  ! 6e31b68d
            log10_radius = 0.5d0*(log10_luminosity + solar_luminosity_cgs - &
                 4.0d0*log_teff - c4pil - csigl)
  ```
- **Now:** [starting_model.f90:828-836](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/starting_model.f90#L828-L836).
- **What was wrong:** when a card moves the fitting point outward, `read_starting_model` integrates a new envelope from L and Teff. Using L_sun = 3.8e33 where log L_sun belongs gave log R ≈ 2e33. The envelope integration then failed with a Saha failure and "ERROR IN BSSTEP", and the run ended without refitting.
- **Expected behaviour, and why:** from L = 4πR²σTeff⁴, log R = (log L/L_sun + log L_sun − 4 log Teff − log 4π − log σ)/2, exactly the F77 line.
- **What we did:** the formula uses `log10_solar_luminosity` again.
  - cm_base with senv0_array(1) = −1e-7 now refits and evolves its ten models.
  - The 0.03 Msun ZAMS pin was re-pinned: it now refits and stops at its first evolution model, instead of in the envelope integration.
  - No other pinned deck takes this path; all are byte-identical.

## 4. Latent bugs and error paths (fixed; pinned output unchanged)

These bugs did not change any pinned deck: either no shipped deck reaches the path, or the wrong value never reached the output. Each was still a real defect for a user who sets a different control, or for the next person to edit the code.

### L1. The pulse files' ε_T and ε_ρ columns were multiplied by ε a second time

- **Status:** fixed in [1d664520](https://github.com/parallelpro/yrec/commit/1d664520) (2026-10-10).
- **In the original F77?** No. It was introduced by [2ed398ca](https://github.com/parallelpro/yrec/commit/2ed398ca) (2026-08-30); the F77 had no GYRE or GSM writer. The F77 convention was already the absolute derivative: in COEFFT `QL` is ε itself, so `QET` and `QED` are dε/dln T and dε/dln ρ ([coefft.f:245-247](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L245-L247)):
  ```fortran
              QL = E
              QLT = QLT + QET + QED*QDT
              QLP = QLP + QED*QDP
  ```
- **Now:** [stitched_model.f90:584-597](https://github.com/parallelpro/yrec/blob/e8549036/src/core/diagnostics/stitched_model.f90#L584-L597) and [pulsation_output.f90:229-242](https://github.com/parallelpro/yrec/blob/e8549036/src/io/pulsation_output.f90#L229-L242).
- **What was wrong:** the stored fields were named `pulse_dlneps_*`. Misled by that name, 2ed398ca multiplied them by ε in both writers. That was right for κ, whose fields really are log-derivatives, but wrong for ε. At the solar centre the columns implied dln ε/dln T = 171 instead of 11.
- **Expected behaviour, and why:** the GSM/GYRE columns are ε · dln ε/dln T = dε/dln T (MESA's `d_epsnuc_dlnT`). The stored field already holds exactly that.
- **What we did:** both writers copy the field unchanged, and the fields are renamed `pulse_deps_*` so the names say what the values are. Only the 15 MESA-baseline `.GYRE` profiles changed (old = ε × new exactly). Only nonadiabatic GYRE runs read these columns.

### L2. A stable rotation sub-step left the redo flag from the previous attempt (SECULR)

- **Status:** fixed in [66d7788f](https://github.com/parallelpro/yrec/commit/66d7788f) (2026-10-12).
- **In the original F77?** Yes. SECULR returns early when no instability occurs ([seculr.f:196-197](https://github.com/parallelpro/yrec/blob/6cd56739/src/seculr.f#L196-L197)), before CHECKJ or CHECKC can write `LREDO`. GETW keeps `LREDO` in a SAVEd local and acts on it ([getw.f:283-293](https://github.com/parallelpro/yrec/blob/6cd56739/src/getw.f#L283-L293)):
  ```fortran
  C  IF LDO=F,NO INSTABILITIES OCCUR (STABLE AGAINST ALL MECHANISMS).
           IF(.NOT.LDO) GOTO 9999
  ...
           IF(LREDO)THEN
              FX = 2.0D0*DT/DELTS
              SUMDT = SUMDT - 2.0D0*DT
  ```
- **Now:** [secular_transport.f90:156-159](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/secular_transport.f90#L156-L159).
- **What was wrong:** after a redo, if the retried sub-step was stable against every mechanism, the flag was still `.TRUE.` from the failed attempt. GETW then rewound the clock again and restored J, μ and the composition, although nothing had failed. In the f90 the dummy is `intent(out)`, so its value on entry was undefined anyway.
- **Expected behaviour, and why:** a sub-step that needs no transport needs no redo. The flag must describe the current call only.
- **What we did:** `redo_flag = .false.` on entry. Byte-identical on every rotating deck.

### L3. The three OPAL EOS vintages shared one "tables loaded" flag

- **Status:** fixed in [85f74924](https://github.com/parallelpro/yrec/commit/85f74924) (2026-10-11).
- **In the original F77?** Yes. ESAC, ESAC01, ESAC06 and the three RHOOFP readers all test the same `COMMON /LREADCO/ ITIME` ([esac.f:89-90](https://github.com/parallelpro/yrec/blob/6cd56739/src/esac.f#L89-L90)):
  ```fortran
        COMMON/LREADCO/ITIME
        IF(ITIME .NE. 12345678) THEN
          ITIME=12345678
  ```
- **Now:** [opal_eos_tables.f90:244-246](https://github.com/parallelpro/yrec/blob/e8549036/src/eos/opal_eos_tables.f90#L244-L246).
- **What was wrong:** once one vintage (1995, 2001 or 2006) had loaded its table, a second vintage used later in the same process skipped its own read and interpolated in an empty table.
- **What we did:** the flag is now an array indexed by vintage. No production run uses two vintages. The new `test_eos` loads all three in one process, and before the fix its OPAL95 and OPAL2001 results came out identical.

### L4. The Saha cutoff's default of log T = 6.0 was lost

- **Status:** fixed in [8a69dc72](https://github.com/parallelpro/yrec/commit/8a69dc72) (2026-10-11).
- **In the original F77?** No; the F77 default is right ([parmin.f:397-398](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L397-L398), `TSCUT = 6.0`). It was lost in two steps on 2026-08-20:
  - [9acfb574](https://github.com/parallelpro/yrec/commit/9acfb574) moved the default to a renamed variable, which broke reading `TSCUT`.
  - [54965bac](https://github.com/parallelpro/yrec/commit/54965bac) restored the name as a new uninitialized local and copied it over the 6.0 after the read.
- **Now:** [legacy_controls.f90:340-349](https://github.com/parallelpro/yrec/blob/e8549036/src/io/legacy_controls.f90#L340-L349).
- **What was wrong:** a deck that did not set `TSCUT` got a Saha cutoff of log T = 0, and the documentation said 0.
- **What we did:** the default is 6.0 again. Every shipped deck sets `TSCUT`, so no run changes.

### L5. I/O errors ended in a Fortran runtime error

- **Status:** fixed in [a9de5f99](https://github.com/parallelpro/yrec/commit/a9de5f99) (2026-10-11).
- **In the original F77?** Yes. None of the opens or namelist reads has `IOSTAT` or `ERR` ([parmin.f:741-744](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L741-L744)):
  ```fortran
        OPEN(UNIT=ISTAND, FILE=YREC1, STATUS='OLD')
        OPEN(UNIT=IRUN, FILE=YREC2, STATUS='OLD')
        READ(UNIT=ISTAND, NML=CONTROL)
        READ(UNIT=IRUN, NML=PHYSICS)
  ```
  Before the fix, only 1 of the 18 opens and none of the 5 namelist reads in the I/O code checked their status.
- **Now:** `check_io`, [io_units.f90:19-40](https://github.com/parallelpro/yrec/blob/e8549036/src/io/io_units.f90#L19-L40).
- **What was wrong:** a missing starting model, a misspelled control or an unwritable output directory ended in a runtime error and a backtrace, with no message in YREC's terms.
- **Expected behaviour, and why:** the run should stop with a message that names the file and the reason, and exit with a nonzero status.
- **What we did:** every open and namelist read goes through `check_io`, which prints `ERROR: cannot <action> <path>: <reason>` and returns `ierr = 1`; the program exits with status 1. Successful runs are byte-identical.

### L6. 188 documented control defaults were wrong

- **Status:** fixed in [c3f40675](https://github.com/parallelpro/yrec/commit/c3f40675) (2026-10-11).
- **In the original F77?** No; the F77 defaults are right (for example `ATMERR = 3.0D-4`, [parmin.f:399-400](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L399-L400)). The wrong documentation came from the hand-extracted controls registry ([4b6b6470](https://github.com/parallelpro/yrec/commit/4b6b6470), 2026-08-22), from which the `.defaults` files were generated.
- **What was wrong:** 188 entries were wrong, 184 of them shown as "0.0 (implicit)". For example, the code really sets `atm_error_tol = 3.0d-4`, `absolute_tolerance = 1.0d-5`, `wind_law_index = 1.5` and the solar luminosity to 3.8515d33.
- **What we did:** `gen_inlist.py` now reads every default from the code and stops if a control has none. This was a documentation fix; the reader is unchanged.

### L7. A rename script changed a routine name inside an error message

- **Status:** fixed in [841ff500](https://github.com/parallelpro/yrec/commit/841ff500) (2026-10-11).
- **In the original F77?** No; the F77 text is right ([eqstat2.f:129-132](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L129-L132)). It was introduced by [1ab2dc75](https://github.com/parallelpro/yrec/commit/1ab2dc75) (2026-10-11), whose merge script renamed `eqstat2` inside a string literal.
- **What was wrong:** the 0.03 Msun deck's log read "blend_eos: BETA is negative" instead of "EQSTAT2: BETA is negative", with identical numbers. The full battery caught it.
- **What we did:** restored the text, and the merge script no longer renames inside strings.

### L8. Thirteen guards against overruns, divisions by zero and edge cases

- **Status:** fixed in [0d66503f](https://github.com/parallelpro/yrec/commit/0d66503f) (2026-10-12).
- **In the original F77?** Yes for 11 of the 13; numbers 9 and 13 were introduced in 2026.
- **Expected behaviour, and why:** each guard acts only where the code overran an array, divided by zero, hit an exact-equality edge or followed a failed lookup. Every case that already ran correctly gives the same numbers.

| # | Defect | F77 | Now |
|---|---|---|---|
| 1 | The rezone's X and Z point insertion could run past `JSON`; only the χ placement checked the limit. | [hpoint.f:386-448](https://github.com/parallelpro/yrec/blob/6cd56739/src/hpoint.f#L386-L448) | [rezoning.f90:577-584](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/rezoning.f90#L577-L584): stops with `ierr` |
| 2 | STARIN printed "RUN TERMINATED" when it could not extend the core, then went on shifting points past `JSON`. | [starin.f:377-391](https://github.com/parallelpro/yrec/blob/6cd56739/src/starin.f#L377-L391) | [starting_model.f90:553-567](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/starting_model.f90#L553-L567): returns `ierr = 1` |
| 3 | LIBURN read `HRO(0)` when the star was fully convective at the start of the step (only `JNENV`, not `JNENV0`, is tested `> 1`). | [liburn.f:89-91](https://github.com/parallelpro/yrec/blob/6cd56739/src/liburn.f#L89-L91) | [convection_geometry.f90:136-151](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/convection_geometry.f90#L136-L151) |
| 4 | The convective turnover search read an unset core index (SAVEd, so stale or 0), and its spline stencils could reach `M+1`. | [tauintnew.f:45-53](https://github.com/parallelpro/yrec/blob/6cd56739/src/tauintnew.f#L45-L53) | [turnover_timescale.f90:159-161](https://github.com/parallelpro/yrec/blob/e8549036/src/core/diagnostics/turnover_timescale.f90#L159-L161) |
| 5 | A zero mass-change rate divided by zero in the timestep limit. | [mdot.f:59-64](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L59-L64) | [mass_change.f90:438-455](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change.f90#L438-L455): no limit for a zero rate |
| 6 | The settling sub-step count could be zero, giving a division by zero and no settling. | [mix.f:413-415](https://github.com/parallelpro/yrec/blob/6cd56739/src/mix.f#L413-L415) | [composition_step.f90:289-303](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/composition_step.f90#L289-L303): at least one sub-step (the reversed `MOD` is open, O16) |
| 7 | A rescale card that left Z unset (`RSCLZ` defaults to 0) set Z = 0. | [rscale.f:67-75](https://github.com/parallelpro/yrec/blob/6cd56739/src/rscale.f#L67-L75) | [model_rescaling.f90:134-141](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/model_rescaling.f90#L134-L141): applied only when > 0 |
| 8 | An envelope mass of exactly −1e-12 was integrated but not reversed (`.GT.` versus `.LT.`). | [envint.f:412, 683](https://github.com/parallelpro/yrec/blob/6cd56739/src/envint.f#L683) | [envelope_integration_kernel.f90:928-931](https://github.com/parallelpro/yrec/blob/e8549036/src/core/envelope/envelope_integration_kernel.f90#L928-L931) |
| 9 | After the Allard tables failed, the fallback integration used HSRA instead of the gray T(τ): the kernel switched its own copy of the atmosphere choice, but the integrand read the run's. | No; introduced by [f570e41e](https://github.com/parallelpro/yrec/commit/f570e41e) (2026-08-30). The F77 switched the shared `KTTAU` ([envint.f:146-153](https://github.com/parallelpro/yrec/blob/6cd56739/src/envint.f#L146-L153)). | [envelope_integrands.f90:76-81](https://github.com/parallelpro/yrec/blob/e8549036/src/core/envelope/envelope_integrands.f90#L76-L81) |
| 10 | The envelope shifts moved every envelope array except the convective flag, which was then read one point off. | [getnewenv.f:130-143](https://github.com/parallelpro/yrec/blob/6cd56739/src/getnewenv.f#L130-L143), starin.f:625-652 | [envelope_rebuild.f90:223-243](https://github.com/parallelpro/yrec/blob/e8549036/src/core/envelope/envelope_rebuild.f90#L223-L243) |
| 11 | Below the burning cutoff, the shell kept the last burning shell's ε derivatives, which reached the pulse arrays. | [coefft.f:213-218](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L213-L218), 412-415 | [shell_equations.f90:244-256](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/shell_equations.f90#L244-L256): zeroed |
| 12 | An unlabelled debug print in the mass change (`WRITE(*,*)TAUCZ/CSECYR,RTOT/CRSUN`). | [massloss.f:71](https://github.com/parallelpro/yrec/blob/6cd56739/src/massloss.f#L71) | [mass_change.f90:155-165](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change.f90#L155-L165): labelled, diagnostics only |
| 13 | The calibration's give-up message said 15 cycles; the loop gives up after 14 (`ICONV` starts at 1). | No; introduced by [cc1b9ba7](https://github.com/parallelpro/yrec/commit/cc1b9ba7) (2026-08-31). The F77 gave up silently ([main.f:961-964](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L961-L964)). | [run_driver.f90:265-273](https://github.com/parallelpro/yrec/blob/e8549036/src/core/run/run_driver.f90#L265-L273) |

- **What we did:** byte-identical on the quick check, the rotating decks and the extra decks. Visible changes: the accretion run loses the debug print (12); the pulse files' two ε-derivative columns change below the burning cutoff (11); and the controls dumps were regenerated for the rescale card (7).

### L9. The envelope reversal swapped the luminosity array one way only

- **Status:** fixed in [4ce66e49](https://github.com/parallelpro/yrec/commit/4ce66e49) (2026-09-04). The same day, [17d8852a](https://github.com/parallelpro/yrec/commit/17d8852a) turned all 21 swap triplets into whole-array reversals.
- **In the original F77?** No; the F77 envelope has no luminosity array. It was introduced in 2025:
  - `ENVL` was added with a correct swap in de5a940f (2025-08-27).
  - In 812669e0 (2025-09-08), `ENVL(I2) = DUM1` became `EQDT(I2) = DUM1`.
  - 081a7361 (2025-09-17) removed the stray line, which left `ENVL` swapped one way ([envint.f:781-785 at 081a7361](https://github.com/parallelpro/yrec/blob/081a7361/src/envint.f#L781-L785)):
  ```fortran
              DUM1 = ENVL(I1)
              ENVL(I1) = ENVL(I2)
              DUM1 = EQDT(I1)
  ```
- **Now:** [envelope_integration_kernel.f90:927-954](https://github.com/parallelpro/yrec/blob/e8549036/src/core/envelope/envelope_integration_kernel.f90#L927-L954).
- **What was wrong:** after the reversal the array's upper half held the lower half's values, and the lower half was unchanged.
- **What we did:** every array is now reversed as a whole, so no single swap line can be dropped again. No output changes, because the envelope luminosity is the same constant at every envelope point today.

### L10. A truncated OPAL95 opacity table returned success (LL95TBL)

- **Status:** fixed in [d62e5360](https://github.com/parallelpro/yrec/commit/d62e5360) (2026-09-04); one gap remains (below).
- **In the original F77?** Yes. Every body-row read has `END=9999`, and label 9999 is a bare `RETURN` ([ll95tbl.f:72-148](https://github.com/parallelpro/yrec/blob/6cd56739/src/ll95tbl.f#L72-L148)):
  ```fortran
        READ(ILIV95,20,END=1000) (R0GR(I),I=1,NUMD)
           READ(ILIV95,30,END=9999)TGR(I),(CAPPA2(N,I,J),J=1,NUMD)
  C ...
        CALL OP95ZTAB(Z)
   9999 CONTINUE
        RETURN
  ```
- **Now:** `read_opal95_tables`, [opal95_opacity.f90:88-101](https://github.com/parallelpro/yrec/blob/e8549036/src/kap/opal95/opal95_opacity.f90#L88-L101) (end of file in a header) and [:185-188](https://github.com/parallelpro/yrec/blob/e8549036/src/kap/opal95/opal95_opacity.f90#L185-L188) (end of file in the body).
- **What was wrong:** a table cut off in mid-body returned silently. The unit was not closed and OP95ZTAB never ran, so the fixed-Z table was never built. An end of file at a table header (`END=1000`) was also silent, and OP95ZTAB then built the fixed-Z table from a partly filled `CAPPA2`. The 2026 conversion added an `ierr` argument but kept the bare return, so `ierr` stayed 0.
- **Expected behaviour, and why:** a short or truncated table set must stop the run. Otherwise the opacities come from unfilled table memory, with no warning.
- **What we did:** both paths now print a message, close the unit and return `ierr = 1`. Every shipped OPAL95 file has all 126 tables, so the result is byte-identical.
- **Still open:**
  - A file that ends exactly after a complete table still stops with a Fortran runtime error, because the next table's title read has no `end=` ([opal95_opacity.f90:169](https://github.com/parallelpro/yrec/blob/e8549036/src/kap/opal95/opal95_opacity.f90#L169)). This is loud, but not the clean path.
  - A file with no TABLE line prints an undefined table number.

### L11. A dead branch set the convective velocity to `1.0D0-11`, which is −10 (TPGRAD)

- **Status:** fixed in [67d49e77](https://github.com/parallelpro/yrec/commit/67d49e77) (2026-09-04).
- **In the original F77?** Yes ([tpgrad.f:149-154](https://github.com/parallelpro/yrec/blob/6cd56739/src/tpgrad.f#L149-L154)):
  ```fortran
        TEST = G*(-QDT)*PRESHT*DELDEL
        IF(TEST.GT.0.0D0)THEN
           VEL = V*DSQRT(G*(-QDT)*(CMIXL**2)*PRESHT*DELDEL/8.0D0)
        ELSE
           VEL = 1.0D0-11
        ENDIF
  ```
- **Now:** [temperature_gradients.f90:223-226](https://github.com/parallelpro/yrec/blob/e8549036/src/convection/temperature_gradients.f90#L223-L226).
- **What was wrong:** `1.0D0-11` is the expression 1 − 11 = −10, not the 1.0D-11 placeholder used on the other failure paths (tpgrad.f:60, 138). The branch can never run, though. This TEST is the earlier TEST at line 106 times PRESHT² > 0, and a non-positive earlier TEST already stops at line 117 (`STOP 912`, now `ierr = 1`).
- **Expected behaviour, and why:** no second sign test is needed. If one were kept, it should read `1.0D-11`.
- **What we did:** removed the branch and left a comment explaining why the radicand is positive there. No output changes.

### L12. The triangle-reset condition read ambiguously (CRRECT)

- **Status:** fixed (parenthesized as the compiler parses it) in [4ce66e49](https://github.com/parallelpro/yrec/commit/4ce66e49) (2026-09-04).
- **In the original F77?** Yes ([crrect.f:67](https://github.com/parallelpro/yrec/blob/6cd56739/src/crrect.f#L67)); the level driver is at [main.f:680-699](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L680-L699):
  ```fortran
        IF(LNEW.OR.LRESET .AND.ITLVL.EQ.2) LSBC = .TRUE.
  ```
- **Now:** [henyey_iteration.f90:82-83](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/henyey_iteration.f90#L82-L83).
- **What was wrong:** with Fortran precedence this is `LNEW .or. (LRESET .and. ITLVL==2)`, which a reader can easily take for `(LNEW .or. LRESET) .and. ITLVL==2`. The second term has no effect, because MAIN already sets `LSBC = .TRUE.` for level 2 ([main.f:697](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L697)).
- **Expected behaviour, and why:** the parsed reading is the intended one. A new triangle must be built at level 1 ("FOR THE FIRST MODEL OF A RUN, THE TRIANGLE IS GENERATED HERE", main.f:681). The other reading would make the statement a no-op and leave a new triangle unbuilt at level 1.
- **What we did:** parenthesized it as parsed; byte-identical. The redundant second term could also be dropped.

### L13. LDIFAD (advection-diffusion transport of angular momentum) never worked

- **Status:** fixed in [15d569e3](https://github.com/parallelpro/yrec/commit/15d569e3) (2026-10-10): a rotating run that sets it now stops. The path itself was commented out in [66d7788f](https://github.com/parallelpro/yrec/commit/66d7788f) (2026-10-12).
- **In the original F77?** Yes. SECULR calls DADCOEFT when `LDIFAD` is set ([seculr.f:296-307](https://github.com/parallelpro/yrec/blob/6cd56739/src/seculr.f#L296-L307)). DADCOEFT reads equal-grid coefficients that no routine assigns ([dadcoeft.f:533-551](https://github.com/parallelpro/yrec/blob/6cd56739/src/dadcoeft.f#L533-L551)):
  ```fortran
           IF(.NOT.LDIFAD)THEN
              CALL DCOEFT(ECOD,DR,DT,EI,EJ,EW,NTOT,WIND1,WIND2,LFIX)
           ELSE
              CALL DADCOEFT(DR,DT,EI,EJ,EW,NTOT,WIND1,WIND2,DJ,
  C dadcoeft.f
           ECOD3(I) = 0.2D0*V1*FGEOM(I)
           V0 = ECOD3(I)/FGEOM(I)
  ```
- **Now:** [secular_transport.f90:213-230](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/secular_transport.f90#L213-L230) (`refuse_advection_diffusion`).
- **What was wrong:** the coefficients DADCOEFT reads are never assigned anywhere in the F77. They are `FGEOM`, `ECOD5`, `ECOD6` and `EV0`…`EV2B` in `/DIFAD3/`, `QCHIRE` in `/EGRIDCHI/`, and `FGEOMIX` and `REQ` in `/DIFAD4/`. So `ECOD3 = 0`, `V0 = 0/0`, and the composition mixing is zero.
- **Expected behaviour, and why:** a run that asks for a transport the code cannot compute must stop with a message, not run on zero or NaN coefficients.
- **What we did:** `apply_secular_transport` returns `ierr = 1` with a message when the flag is set; `test_error_exit.py::test_advective_transport_rejected` covers it. 16 legacy decks and 316 inlists set `LDIFAD = .TRUE.`, but none of them rotate. The check therefore sits where the flag acts, not at startup. No pinned output changes.

### L14. The envelope triangle could hang, go NaN, or accept a non-positive size

- **Status:** fixed in [d2ff6886](https://github.com/parallelpro/yrec/commit/d2ff6886) (2026-10-14).
- **In the original F77?** Yes, all three parts, in SURFBC (the only copy of the walk and the fit):
  - The reflection walk has no limit ([surfbc.f:110-125](https://github.com/parallelpro/yrec/blob/6cd56739/src/surfbc.f#L110-L125)):
    ```fortran
     10      CONTINUE
              IF(TEMP.LT.-TRIERR) THEN
                 FTRI = -FTRI
                 TRIT(I1) = TRIT(I2) + TRIT(I3) - TRIT(I1)
                 TRIL(I1) = TRIL(I2) + TRIL(I3) - TRIL(I1)
                 GOTO 10
    ```
  - The plane fit has no pivot check ([surfbc.f:166-179](https://github.com/parallelpro/yrec/blob/6cd56739/src/surfbc.f#L166-L179)):
    ```fortran
           TEMP = 1.0D0/(TS(2)-TS(1))
           TEMP1 = (TS(3)-TS(1))*TEMP
           TEMP2 = 1.0D0/(PS(3)-PS(1)-TEMP1*(PS(2)-PS(1)))
    ```
  - `TRIDT` and `TRIDL` are read through the PHYSICS namelist with no check (defaults at [parmin.f:435](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L435)).
- **Now:**
  - Walk and fit: [triangle_fit.f90:154-179](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/triangle_fit.f90#L154-L179) and [:201-233](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/triangle_fit.f90#L201-L233); the limits are at lines 35-43.
  - Error reports: [surface_boundary.f90:204-214](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/surface_boundary.f90#L204-L214).
  - Size check: [derived_controls.f90:44-53](https://github.com/parallelpro/yrec/blob/e8549036/src/io/derived_controls.f90#L44-L53).
- **What was wrong:**
  1. The walk never ends for a zero-area triangle (`TRIDT = 0`), or for one whose orientation flag disagrees with its vertex order (`TRIDT < 0` gives a clockwise start with `FTRI = +1`).
  2. Two vertices at the same fitting-point log T, or fitting points collinear in (log P, log T), give Inf/NaN boundary coefficients silently.
  3. Non-positive sizes were accepted.
- **Expected behaviour, and why:** a degenerate triangle must end in a clean error, not a hang or NaN coefficients. Pivoting was not added, because it would move every deck at round-off.
- **What we did:**
  - The walk gives up after 1000 reflections with `ierr = 1`; the decks in the battery need at most 11.
  - The fit refuses a pivot smaller than 1e-10 of the differences it divides; the decks never go below 6e-3.
  - `derive_controls` refuses a non-positive `tri_delta_teffl` or `tri_delta_logl`.
  - Byte-identical on all decks. The unit test `test_triangle_fit` and `test_error_exit.py::test_nonpositive_triangle_size_refused` cover the new paths.

### L15. Diagnostics that were always zero or never written

- **Status:** removed in [c4f359c8](https://github.com/parallelpro/yrec/commit/c4f359c8) (2026-10-13).
- **In the original F77?** For six of the eight arrays, yes:
  - **`valfmlt`, `vphmlt`, `vcmxmlt`:** yes. TPGRAD only ever zeroes the scalars ([tpgrad.f:38-40](https://github.com/parallelpro/yrec/blob/6cd56739/src/tpgrad.f#L38-L40)). COEFFT copies them (coefft.f:422-424), and WRTMOD writes them only in the `IPVER = 3` pulse format (wrtmod.f:226-231; the default is `IPVER = 1`).
    ```fortran
          ALFMLT=0.0D0
          PHMLT=0.0D0
          CMXMLT=0.0D0
    ```
  - **`fxion_zone` (`SFXION`):** yes. It is stored at [coefft.f:371-373](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L371-L373); PUTSTORE writes zeros in its columns, and its `SFXION` write is commented out (putstore.f:221-226).
  - **He3 luminosities (`UMHE33`, `UMHEAL` in `/GRAB/`):** yes. They are set at [coefft.f:240-241](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L240-L241) and never read.
  - **`pulse_dlnrho_dlnp` (`PQDP`):** no. WRTMOD wrote it ([wrtmod.f:215-231](https://github.com/parallelpro/yrec/blob/6cd56739/src/wrtmod.f#L215-L231)). It lost its last reader when [6135144e](https://github.com/parallelpro/yrec/commit/6135144e) (2026-08-26) retired the .pmod/.penv/.patm files.
  - **`scp`:** no. It was added after 6cd56739 (812669e0 and 081a7361, 2025-09) and lost its reader when [49f9d57e](https://github.com/parallelpro/yrec/commit/49f9d57e) (2026-08-28) retired the .store file.
- **Now:** the arrays and their stores are gone.
- **What was wrong:** three arrays were always zero. The other five were computed on every Newton iteration and written by no output routine.
- **Expected behaviour, and why:** a stored per-zone quantity should be emitted by some writer, or not exist. Otherwise it costs time and misleads whoever reads the code.
- **What we did:** removed the eight members, their stores, and the network's `eps_he3_he3` and `eps_he3_burning`, which fed only the He3 luminosities. Byte-identical on all decks.

### L16. The rezone's store of ∇_ad − ∇ for the rotation step was dead (HPOINT)

- **Status:** removed in [56aa2e6f](https://github.com/parallelpro/yrec/commit/56aa2e6f) (2026-10-13).
- **In the original F77?** Yes. HPOINT stores `DDELM` for every zone ([hpoint.f:758-759](https://github.com/parallelpro/yrec/blob/6cd56739/src/hpoint.f#L758-L759)). SETUPV overwrites it ([setupv.f:281-286](https://github.com/parallelpro/yrec/blob/6cd56739/src/setupv.f#L281-L286)) before VCIRC reads it, over the same zone range and in the same SECULR call:
  ```fortran
  C MHP 06/02                                                       (HPOINT)
              DDELM(IM) = DELAO(IM) - DELO(IM)
           DDELM(I) = MAX(DELAMI(I)-DELMI(I)+DELMU,DMIN)          ! SETUPV
            DDEL=DDELM(I)+ QMU/QP                                 ! VCIRC
  ```
- **Now:** a comment remains at [rezoning.f90:1008-1023](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/rezoning.f90#L1008-L1023). The writer that matters is [stability_coefficients.f90:184-185](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/stability_coefficients.f90#L184-L185).
- **What was wrong:** HPOINT's zone-centred value (no μ term, no floor) was never read; it was a different quantity in any case. The store also made the core write into the rotation step's scratch.
- **What we did:** removed the store. As a check, the array was filled with NaN after the rezone, and every rotating deck stayed byte-identical.

### L17. `cz_moment_of_inertia` was `intent(out)` but read before it was set

- **Status:** fixed in [2d62df76](https://github.com/parallelpro/yrec/commit/2d62df76) (2026-10-14).
- **In the original F77?** No. It was introduced by [94c7f45c](https://github.com/parallelpro/yrec/commit/94c7f45c) (2026-08-20, "Modernize seculr.f"). In the F77, `HICZ` is an ordinary argument of SECULR. Its solid-body branch passes the value MIDMOD set, the whole star's I, to MWIND without changing it ([seculr.f:170-188](https://github.com/parallelpro/yrec/blob/6cd56739/src/seculr.f#L170-L188)):
  ```fortran
           ELSE IF(LCZSUR .AND. LSOLID)THEN
              LJDOT = LJDOT0
              HSBOTJ = 0.0D0
              HSTOPJ = EXP(CLN*HSTOT)
              IMAXJ = 1
              CALL MWIND(BL,DT,HSBOTJ,HSTOPJ,IMAXJ,M,LJDOT,OMEGAS,
  ```
- **Now:** [secular_transport.f90:99-103](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/secular_transport.f90#L99-L103) (the declaration) and [:369-380](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/secular_transport.f90#L369-L380) (the read).
- **What was wrong:** the solid-body spin-down branch passes the argument to the braking before assigning it. The value of an `intent(out)` dummy on entry is undefined; it worked only because compilers keep the caller's scalar.
- **What we did:** made it `intent(inout)` and documented the contract; byte-identical.

## 5. Open: structure, solver and composition

These are found, recorded and left as they are, because each fix changes results and needs its own measured commit. Where a fix is proposed, its expected effect is given.

### O1. Gaps in the Henyey Jacobian: the convective gradient and the OPAL density derivatives

- **Status:** open.
- **In the original F77?** Yes. In a convective shell COEFFT sets the L derivative of the temperature equation to zero ([coefft.f:203-208](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L203-L208)). Near the adiabat TPGRAD takes the gradient's P and T derivatives from the EOS's ∇_ad derivatives ([tpgrad.f:139-143](https://github.com/parallelpro/yrec/blob/6cd56739/src/tpgrad.f#L139-L143)). With an OPAL table, OPAL's ρ derivatives are adopted as they are ([eqstat2.f:292-293](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L292-L293)).
  ```fortran
  C TEMPERATURE GRADIENT IS CONVECTIVE
           QTL = 0.0D0
           QTP = QT*(-1.0D0 + QACP)
           QTT = QT*QACT
  ```
- **Now:** [shell_equations.f90:224-234](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/shell_equations.f90#L224-L234) ("not linearized"), [temperature_gradients.f90:128-133](https://github.com/parallelpro/yrec/blob/e8549036/src/convection/temperature_gradients.f90#L128-L133) and [blended_eos.f90:488-489](https://github.com/parallelpro/yrec/blob/e8549036/src/eos/blended_eos.f90#L488-L489).
- **What is wrong:** there are three gaps.
  - In the superadiabatic layers ∇ depends on L through ∇_rad, but that column of the Jacobian is zero.
  - The (T, ln P) and (T, ln T) entries rest on the EOS's ∇_ad derivatives. `test_structure` finds these off by up to 72% of the entry, which is under 1% of the residual.
  - OPAL's dln ρ/dln P and dln ρ/dln T agree with finite differences of the tabulated ρ(P, T) only to about 2e-3.
- **Expected behaviour, and why:** Newton converges quadratically only when the Jacobian matches the residuals it differentiates. These gaps leave the converged model unchanged but slow convergence and make it less robust. That matters most in thin, superadiabatic, partly degenerate envelopes: low-mass stars and brown dwarfs (see O3).
- **What we did / proposed fix:** `test_structure` (14c2215a) checks every partial against finite differences and marks these three as known. Proposed: run the 0.03 Msun deck with numerical derivatives for these entries first, and fix the analytic ones only if that cures the deck's two-cycle (O3).

### O2. The EOS second derivatives with an OPAL table, and in the Saha ramp

- **Status:** open. A mislabelled name was fixed; the numbers are unchanged.
- **In the original F77?** Yes. The ramp weight's derivative is commented out, and `QDTT`/`QDTP` are blended without forming the difference ([eqstat2.f:236-266](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L236-L266)). The OPAL branches replace only ρ and the first derivatives ([eqstat2.f:281-310](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L281-L310)):
  ```fortran
        WT0 = WTS2*WTS*(6.0D0*WTS2 - 15.0D0*WTS + 10.0D0)
  C      WT1 =-SCALES*30.0D0*WTS2*(1.0D0-WTS)**2 ! Derivative of WT0
  ...
           QDTT = QDTT + WT0*QDTT0
           QDTP = QDTP + WT0*QDTP0
  ```
  The mislabel was not in the F77. It was introduced by f0640110 (2026-08-19), which named dχ_T/dln P `dlnrho_dlnp_dt`; 5d39f836 renamed it `dlnrho_dlnt_dp`.
- **Now:** [blended_eos.f90:380-417](https://github.com/parallelpro/yrec/blob/e8549036/src/eos/blended_eos.f90#L380-L417) (the ramp) and [:478-514](https://github.com/parallelpro/yrec/blob/e8549036/src/eos/blended_eos.f90#L478-L514) (OPAL adoption).
- **What is wrong:**
  - Inside an OPAL table, ρ and its first derivatives are OPAL's, but the six second derivatives (dχ_T/dln T, dχ_T/dln P, d∇_ad/dln T and dln P, dc_P/dln T and dln P) are left from the Yale solution. They are 8-180% off, with one sign flip.
  - In the Saha-to-fully-ionized ramp the weight w(log T) is applied without dw/dlog T.
  - `QDTT` and `QDTP` are blended as f_full + w·f_saha instead of f_full + w(f_saha − f_full).
- **Expected behaviour, and why:** the second derivatives must belong to the EOS actually returned. For a blend f = f_full + w(T)(f_saha − f_full), d/dln T must include (dw/dln T)(f_saha − f_full). These derivatives feed the ε_grav and convective-gradient partials (S3, O1), so they affect convergence, not converged models.
- **What we did / proposed fix:** `test_eos` (1c177c25, b96442be) checks every derivative against centred finite differences and pins the numbers; `eos/README.md` documents the OPAL caveat. `use_numerical_derivatives` agrees with finite differences to ~1e-6. Proposed, tied to O3: use numerical derivatives inside OPAL, restore the dw terms, and keep the change only if it measurably helps convergence.

### O3. Low-mass convergence: the 0.03 Msun two-cycle and the 0.05 and 0.1 Msun Allard decks

- **Status:** open.
- **In the original F77?** Yes for the 0.03 Msun deck: the original code also hits a Newton two-cycle, after 210 models (0.047 Gyr). The Allard testsuite decks are newer than 6cd56739, so no F77 run exists. Neither code cuts the timestep and retries after a failed solve ([main.f:775-778](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L775-L778)).
- **Now:** [step_driver.f90:140-159](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/step_driver.f90#L140-L159) (a diverged step is not retried) and [starting_model.f90:175-186](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/starting_model.f90#L175-L186) ("MODEL FAILED TO CONVERGE"). The decks are `examples/run_from_dbl_to_zams/Test_m0030_feh0_noGS_ZAMS.inlist`, `testsuite/Test_m0p05feh+0p0GS98_allard_15Gyr.inlist` and `testsuite/Test_m0p1feh+0p0GS98_allard_ZAMS.inlist`.
- **What is wrong:**
  - The 0.03 Msun deck's Newton corrections alternate between two states late in the contraction, and the step fails. The failure point moves with any small change: F77 at 210 models, b4673635 at 179, 0304e269 at 83, today at 77. On b4673635, the eqscve blend fix alone moves it to 62, the rhoofp06 fix alone to 210, and the D rate alone to 83.
  - The 0.05 and 0.1 Msun Allard decks stop with "MODEL FAILED TO CONVERGE" after 50-100 models.
  - A new lead from the triangle work (L14): these decks have by far the thinnest envelope triangles. The relative pivot of the plane fit is ~6e-3, against ≥ 0.15 on every other deck. The fitting points are nearly collinear, so the boundary rows' coefficients are amplified roughly 25 times relative to the other decks.
- **Expected behaviour, and why:** with a consistent Jacobian and a smooth residual, Newton converges on each step. A two-cycle points to one of three things: an inaccurate Jacobian (O1, O2); a residual that jumps (a triangle refit between levels, an EOS table edge, a convective-flag flip); or an ill-conditioned surface boundary row.
- **What we did / proposed fix:** bisected the 0.03 Msun drift and failure point, and gave that deck a stop it can reach (e16b304a). L14 added the triangle guards. Proposed, on the 0.03 Msun deck:
  1. Log the triangle fit's conditioning up to the failure.
  2. Try numerical derivatives (O1).
  3. Find where the corrections alternate.
  4. Fix only once the cause is known.

  A step retry with a smaller timestep would also make the code more robust, independent of the cause.

### O4. The starting model's stitched atmosphere reads the envelope radius out of bounds

- **Status:** open (output only).
- **In the original F77?** No. 6cd56739's STITCH has no atmosphere section. The read was added by 081a7361 (J. van Saders, 2025-09-17, "Atm stiched onto envelope"; [stitch.f:160](https://github.com/parallelpro/yrec/blob/081a7361/src/stitch.f#L160)):
  ```fortran
              RAD = DLOG10(DEXP(CLN*ENVR(NUMENV)) + ATMOR(I))
  ```
  The starting-model call that reaches it with `NUMENV = 0` came with b116224d (2026-08-25).
- **Now:** [stitched_model.f90:426-433](https://github.com/parallelpro/yrec/blob/e8549036/src/core/diagnostics/stitched_model.f90#L426-L433).
- **What is wrong:** the starting model's stitched model is built before any envelope point is stored, so the code reads `env_log10_radius(0)`. The first model's atmosphere radii in the profile and pulse output are garbage. Under `-fcheck=all`, cm_base and the two mass-loss decks stop there.
- **Expected behaviour, and why:** the atmosphere sits on the photosphere. With no envelope points, its base radius should come from a defined value, the outer interior shell's radius.
- **Proposed fix:** use the outer shell's radius when the envelope is empty, which changes only the starting model's atmosphere radii in the output. Then add a `-fcheck=all` run of the config matrix to the test battery.

### O5. The accretion term of the entropy equation reaches every zone at iteration levels 1-2 (latent)

- **Status:** open, but latent: the shifts it gates are currently exactly zero.
- **In the original F77?** No, not in this form. COEFFT applies the shift from `JENV` outward ([coefft.f:252-259](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L252-L259)), and CRRECT has `SAVE` ([crrect.f:58](https://github.com/parallelpro/yrec/blob/6cd56739/src/crrect.f#L58)). So at levels 1-2 `JENV` kept the previous model's value, and only a run's very first solve saw zero. The defect came with [1c510cc5](https://github.com/parallelpro/yrec/commit/1c510cc5) (2026-08-22), which removed that `save`. `envelope_zone_index` became a plain local, zero (through `-finit-local-zero`) at levels 1-2.
  ```fortran
              IF(LMDOT.AND.DMDT0.GT.0.0D0)THEN
                 IF(IM.GE.JENV)THEN
                    HTTT = HTT(IM)+DLOGT
                    HPPP = HPP(IM)+DLOGP
  ```
- **Now:** [step_driver.f90:615-641](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/step_driver.f90#L615-L641) and [structure_equations.f90:138-145](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/structure_equations.f90#L138-L145).
- **What is wrong:** with accretion, the accreted layer's log T and log P shifts are added to every shell at levels 1-2, not only to the surface convection zone. The shifts are, however, exactly zero in both codes.
  - The F77 sets `SCEN = 0` (massloss.f:169, 172) before MDOT forms them with `DLNM = 0` (mdot.f:101-106).
  - The current code zeroes `envelope_specific_entropy` ([mass_change.f90:117](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change.f90#L117)) before it is used.
- **Expected behaviour, and why:** the accreted entropy changes only the convective envelope it mixes into. The term should start at the envelope edge at every level.
- **Proposed fix:** pass the predictor's envelope edge at levels 1-2. This is byte-identical today and matters only if the entropy correction is switched back on.

### O6. The rezone's deuterium gate reads the pre-rezone surface index

- **Status:** open.
- **In the original F77?** Yes. HPOINT tests `HCOMP(12,M)` after regridding, but before `M` becomes the new point count (that happens at hpoint.f:678; [hpoint.f:591-602](https://github.com/parallelpro/yrec/blob/6cd56739/src/hpoint.f#L591-L602)):
  ```fortran
        IF(LEXCOM .AND. HCOMP(12,M).GE.1.0D-14)THEN
           DO J = 1,M
              HPO(J) = DRATE0(J)
           END DO
           CALL OSPLIN(HSS,HTO,HS,HPO,NTAB,NTOT)
  ```
- **Now:** [rezoning.f90:787-794](https://github.com/parallelpro/yrec/blob/e8549036/src/core/model/rezoning.f90#L787-L794). The reader is [mid_timestep_model.f90:228-243](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/mid_timestep_model.f90#L228-L243).
- **What is wrong:** the gate does not read the new surface. If the mesh grew, it reads an interior point; if it shrank, a stale value past the new count. When the gate is false, the start-of-step deuterium burning rate is not regridded and stays in old-mesh zone order. The rotation sub-steps then read it zone by zone on the new mesh.
- **Expected behaviour, and why:** every per-zone array that survives the rezone must be moved onto the new mesh.
- **Proposed fix:** test the new surface and regrid the rate unconditionally. The expected drift is at round-off or none.

### O7. Accretion mixing assumes a uniform convection zone, and the step uses two definitions of M_cz

- **Status:** open.
- **In the original F77?** Yes ([mdot.f:45-56](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L45-L56), [:162-168](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L162-L168), [:195-201](https://github.com/parallelpro/yrec/blob/6cd56739/src/mdot.f#L195-L201)):
  ```fortran
        DMCZT = SUMDM
           DMCZ = HS1(M) - HS1(JENV)
              ABUND=(HCOMP(I,JENV)*DMCZT+COMPACC(I)*DELM)/(DELM+DMCZT)
  ```
- **Now:** [mass_change_kernel.f90:44-98](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change_kernel.f90#L44-L98) and [:195-210](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change_kernel.f90#L195-L210).
- **What is wrong:**
  - The new abundance is (X(jenv) M_cz + X_acc ΔM)/(M_cz + ΔM). It uses the base shell's X for the whole zone, so species are conserved only if the zone is already exactly uniform.
  - M_cz is defined two ways. One is m(nz) − m(jenv), the difference of zone-centre masses, used in the timestep limit, the entropy shift and the mass-loss compression. The other is Σdm(jenv…nz), used in the mixing and the thermal-energy means.
  - With jenv = 1, the zone-centre masses are compressed about m(1), which is inconsistent with either definition.
- **Expected behaviour, and why:** mixing should conserve each species exactly, X_new = (ΣX dm + X_acc ΔM)/(Σdm + ΔM). The zone mass should have one definition throughout the step.
- **What we did / proposed fix:** extracted the kernel with `test_mass_change` (a7349330), and fixed S2 and S8 in the same routine. Proposed: mix with ΣX dm and use one M_cz. This changes accreting runs only; no pinned deck accretes.

### O8. `store_li_be_rates` zeroes both rate sets past the first cool shell (LIRATE88)

- **Status:** open.
- **In the original F77?** Yes ([lirate88.f:105-113](https://github.com/parallelpro/yrec/blob/6cd56739/src/lirate88.f#L105-L113)):
  ```fortran
    110 CONTINUE
        DO 120 K = I,M
           RLI6(K)=0.0D0
  ...
           RLI60(K)=0.0D0
  ```
- **Now:** [burn_step.f90:409-435](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/burn_step.f90#L409-L435).
- **What is wrong:** each call zeroes the start-of-step and end-of-step rates from the first shell cooler than log T = 6 outward, whichever set it is filling. Since N6's fix the end-of-step and rotation sub-step calls really do fill the end set, so they can zero start rates stored earlier, and the reverse. A cutoff that differs between calls changes the rates LIBURN interpolates between.
- **Expected behaviour, and why:** a call should write, and clear, only the set it computes.
- **Proposed fix:** clear only the call's own set. A trial measured earlier changes surface Li slightly, so it needs its own measured commit.

### O9. Surface ⁷Li has a noise floor of about 0.005 dex

- **Status:** open.
- **In the original F77?** Yes. The start-of-step convection-zone average starts at the integer zone `JENV0`, while the end-of-step one is corrected with the fractional edge `FX` ([liburn.f:231-274](https://github.com/parallelpro/yrec/blob/6cd56739/src/liburn.f#L231-L274)):
  ```fortran
           DO 68 I = JENV0,M
              BLI60 = BLI60 + RLI60(I)*HS2(I)
  ...
        IF(JENV.GT.1.AND.JENV.LT.M)THEN
           HS2(JENV) = HS2(JENV)+FX*(HS1(JENV)-HS1(JENV-1))
  ```
- **Now:** [burn_step.f90:163-198](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/burn_step.f90#L163-L198).
- **What is wrong:** when a marginal shell at the start-of-step base flips between convective and radiative, the start average jumps by a whole shell. A change at the level of round-off (the deuterium closed form) moved surface ⁷Li at 4.57 Gyr by +0.0053 dex in the non-rotating solar deck; 46 of 348 models flipped.
- **Expected behaviour, and why:** Li depletion depends steeply on the base temperature. Both ends of the rate interpolation should locate the base the same way, continuously in its position, so that round-off cannot cause discrete jumps.
- **Proposed fix:** store a start-of-step edge fraction with the start zone and interpolate the start depth the same way. This changes every deck with extended composition at the ~0.005 dex level.

### O10. The overshoot edge's zone width mixes two models, and the rotation sub-steps gate envelope overshoot on the wrong switch

- **Status:** open.
- **In the original F77?**
  - **The width: yes.** LIBURN divides by a width whose outer radius is the start-of-step `HRO` and whose inner radius is the current `HR` ([liburn.f:111-118](https://github.com/parallelpro/yrec/blob/6cd56739/src/liburn.f#L111-L118)):
    ```fortran
                 IF(R.LT.ROS)THEN
                    JENV = J + 1
                    DR = EXP(CLN*HRO(J+1))-R
                    FX = 0.5D0-((ROS-R)/DR)
    ```
  - **The gate: no.** The F77 gates the envelope overshoot depth on `LOVSTE` in all three places (main.f:627, main.f:836, [getw.f:318](https://github.com/parallelpro/yrec/blob/6cd56739/src/getw.f#L316-L320): `IF(LOVSTE .AND. LCM(M))THEN`). The conversion of GETW in [180849b5](https://github.com/parallelpro/yrec/commit/180849b5) (2026-08-20) wrote `if(lovstm .and. convective_flag_mid(num_zones))`. `LOVSTM` is the intermediate-zone overshoot switch.
- **Now:** the width is at [convection_geometry.f90:162-173](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/convection_geometry.f90#L162-L173). The gate is at [evolve_rotation.f90:367-375](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/evolve_rotation.f90#L367-L375), against [step_driver.f90:651-662](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/step_driver.f90#L651-L662).
- **What is wrong:**
  - The edge fraction used in the end-of-step Li/Be average divides by a zone width taken from two models.
  - With `LOVSTE` on and `LOVSTM` off, a rotating model's sub-steps take the overshoot branch (which is gated on `envelope_overshoot_active`) with zero depth, so their Li/Be burning ignores the overshoot region.
- **Expected behaviour, and why:** the width should be r(k) − r(k−1) of one model. The same physical switch should set the envelope overshoot depth in both places, as in the F77.
- **Proposed fix:**
  - Restore the `envelope_overshoot_active` gate in the rotation sub-steps. No shipped deck has both envelope overshoot and rotation on, so this should change no pinned output.
  - The width fix changes decks with envelope overshoot and extended composition (for example cm_env_overshoot) and needs a measured commit.

## 6. Open: rotation, wind, microphysics and output

### O11. The dynamical-shear limit can overshoot the neighbouring Ω or go negative (VCIRC)

- **Status:** open.
- **In the original F77?** Yes. The trigger tests the centred gradient `QWLNR` (seculr.f:92-94), but `WMIN` is built from the plain difference quotient ([vcirc.f:378-407](https://github.com/parallelpro/yrec/blob/6cd56739/src/vcirc.f#L378-L407)):
  ```fortran
           IF(ABS(QWLNR(I)).GT.QWRMAX(I))THEN
  ...
              IF(OMEGA(I).LT.OMEGA(I-1))THEN
                 WMIN = OMEGA(I) + QWRMAX(I)*CLN*(HR(I)-HR(I-1))
              ELSE
                 WMIN = OMEGA(I) - QWRMAX(I)*CLN*(HR(I)-HR(I-1))
  ```
- **Now:** [circulation.f90:588-620](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/circulation.f90#L588-L620).
- **What is wrong:**
  - The centred gradient (ΔΩ/4)(1/a + 1/b), with a and b the two half-widths in ln r, is ≥ |ΔΩ|/Δln r, with equality only when the centre is midway in ln r. So the trigger can fire while |Ω_i − Ω_{i−1}| is below limit · Δln r.
  - `wmin` then lands beyond Ω_{i−1}: the "limited" profile is steeper than the real one, and the Goldreich-Schubert-Fricke velocity is scaled up instead of down.
  - For Ω_i > Ω_{i−1}, `wmin` can be ≤ 0, and its logarithm is NaN.
  - The branch for that velocity is chosen with `wmin`, while its normalization uses Ω_{i−1}.
- **Expected behaviour, and why:** `wmin` should be the Ω_{i−1} at which the trigger's own centred gradient equals the limit, wmin = Ω_i ± 4 · limit/DR1. Whenever the trigger fires, that value lies strictly between Ω_i and Ω_{i−1}, so it is positive and the correction reduces the velocity, as it should once the dynamical instability has removed the excess shear. The branch should follow the actual profile.
- **Proposed fix:** first measure which decks reach this branch (they print "DYNAMICAL SHEAR-SHELL"), then compute `wmin` from the trigger's gradient in one measured commit. Only rotating runs that hit dynamical shear change, and the NaN path goes away.

### O12. The solid-body Newton iteration's tolerance (`TOLERI`) is never set

- **Status:** open. The plumbing changed; the behaviour did not.
- **In the original F77?** Yes. `COMMON /ERRMOM/ TOLERI` appears only in CHECKJ and SOLID, and nothing sets it. CHECKJ only loosens it temporarily ([checkj.f:135-142](https://github.com/parallelpro/yrec/blob/6cd56739/src/checkj.f#L135-L142)); SOLID reads it ([solid.f:56-65](https://github.com/parallelpro/yrec/blob/6cd56739/src/solid.f#L56-L65)):
  ```fortran
        IF(DABS(DELJ/CZJ).GT.TOLERI) THEN
           IF(ICOUNT.LT.20) THEN
              ICOUNT = ICOUNT + 1
  ```
- **Now:** `default_inertia_tolerance = 0.0d0` at [omega_profile.f90:15-20](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/omega_profile.f90#L15-L20), used at [:211-222](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/omega_profile.f90#L211-L222).
- **What is wrong:** with a tolerance of zero, |ΔJ/J| > 0 practically always holds, so every call takes all 20 Newton corrections. Each correction is a shape pass and a moment-of-inertia pass, for every radiative zone, in every Newton iteration of a rotating run. Only the transport's last iteration passes a loosened tolerance. After the few corrections that reach round-off, the rest move Ω only in its last bits. The 2026 conversion first made it a scratch field that carried its value from one step to the next.
- **Expected behaviour, and why:** a Newton iteration on J = ΩI(Ω) should stop at a real relative tolerance near round-off. It converges quadratically, so Ω would be unchanged to round-off and most of the passes saved.
- **What we did / proposed fix:** [0f2d50a9](https://github.com/parallelpro/yrec/commit/0f2d50a9) made the tolerance an explicit optional argument with default 0, which is byte-identical. Proposed: give the default a real value, for example 1e-12, in a measured commit. Ω moves at round-off and rotating runs get faster.

### O13. At a transport region's ends, the cell sums and the write-back disagree about a neighbouring convection zone

- **Status:** open.
- **In the original F77?** Yes. The sums are in ROTGRID and MIXGRID ([rotgrid.f:104-145](https://github.com/parallelpro/yrec/blob/6cd56739/src/rotgrid.f#L104-L145), [mixgrid.f:55-88](https://github.com/parallelpro/yrec/blob/6cd56739/src/mixgrid.f#L55-L88)), and the write-back in EQ2MOD and MIXCOM ([eq2mod.f:46-69](https://github.com/parallelpro/yrec/blob/6cd56739/src/eq2mod.f#L46-L69), [mixcom.f:97-119](https://github.com/parallelpro/yrec/blob/6cd56739/src/mixcom.f#L97-L119)):
  ```fortran
        IF(IBEG.GT.1)THEN
           DO II = IBEG-1,1,-1
              IF(.NOT.LCZ(II))THEN          ! sums: take in the CZ below
  ...
        IF(.NOT.LCZ(IBEG) .OR. IBEG.EQ.1)THEN
           I0 = IBEG                        ! write-back: only if IBEG itself is convective
  ```
- **Now:** [am_grid.f90:178-211](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/am_grid.f90#L178-L211) (sums) and [:351-378](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/am_grid.f90#L351-L378) (write-back); for the composition, [rotational_mixing.f90:459-485](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/rotational_mixing.f90#L459-L485) and [:267-294](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/rotational_mixing.f90#L267-L294).
- **What is wrong:** the end cell takes in a whole convection zone just beyond the region whenever the neighbouring zone is convective: its mass, I and J, or its mass for the composition. The cell's change is written back into that convection zone only when the end zone itself is convective. When the end zone is radiative but its neighbour is convective, a change computed for the end zone plus the convection zone is applied to the end zone alone. For J, the conservation rescale then spreads the missing share over the region. For the composition there is no rescale (it is commented out in the F77 MIXCOM), so species are not conserved at such ends.
- **Expected behaviour, and why:** the zones summed into an end cell must be the zones that receive its change. Then J and every species are conserved by construction, and the rescale only removes interpolation error.
- **Proposed fix:** widen the write-back range with the same test the sums use, in a measured commit. Only rotating runs with instability transport whose regions end next to a convection zone change.

### O14. The rotational mixing checks miss species 5-15, and an iteration history has 16 slots with no cap (CHECKC, CHECKJ)

- **Status:** open.
- **In the original F77?** Yes. CHECKC checks species 1-11 only on the last iteration ([checkc.f:48-54](https://github.com/parallelpro/yrec/blob/6cd56739/src/checkc.f#L48-L54)), inside the loop. SECULR mixes species 5-11 (5-15 with extended composition) after the loop ([seculr.f:340-376](https://github.com/parallelpro/yrec/blob/6cd56739/src/seculr.f#L340-L376)). CHECKJ declares `DELJI(16)` and `IMAXI(16)` and indexes them by iteration ([checkj.f:27-30](https://github.com/parallelpro/yrec/blob/6cd56739/src/checkj.f#L27-L30)); the iteration count is never bounded.
  ```fortran
        IF(IT.EQ.ITDIF2)THEN
           JEND = 11
        ELSE
           JEND = 4
        ENDIF
  ```
- **Now:** [transport_checks.f90:431-437](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/transport_checks.f90#L431-L437) and [secular_transport.f90:179-209](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/secular_transport.f90#L179-L209).
- **What is wrong:**
  - Species 5-11 are mixed once, after the loop, and nothing checks them afterwards; the last-iteration test sees them before their mixing.
  - If the transport converges early, species 5-11 are never tested. Species 12-15 never are.
  - A negative abundance from that mixing is neither repaired nor answered with a sub-step cut.
  - A deck with more than 16 transport iterations writes past the two history arrays.
- **Expected behaviour, and why:** every mixed species should be checked after its own mixing, with the same repair or cut. The iteration count must fit the history arrays.
- **Proposed fix:**
  - Check species 5 to the last one after they are mixed. A cut there must restore the full composition, since the restart currently restores only slots 1-4.
  - Reject more than 16 iterations at read time. The cap is byte-identical for valid decks.

### O15. The gradient-reversal safeguard acts only above 1 rad/s, and `mu_gradient_scale` does nothing (CHECKJ)

- **Status:** open.
- **In the original F77?** Yes. The entry test uses 1 rad/s, and the routine's other tests use 1e-15 rad/s ([checkj.f:146-160](https://github.com/parallelpro/yrec/blob/6cd56739/src/checkj.f#L146-L160)). `FMU` appears only in commented-out formulas, and the μ-gradient velocity it scaled is zeroed ([vcirc.f:87-93](https://github.com/parallelpro/yrec/blob/6cd56739/src/vcirc.f#L87-L93)).
  ```fortran
        IF(OMEGA(I)-OMEGA(I-1).GT.1.0d0)THEN
  C  IF PREVIOUS GRADIENT WAS POSITIVE, LEAVE ALONE.
           IF(WOLD(I)-WOLD(I-1).GT.1.0D-15)THEN
  ```
- **Now:** [transport_checks.f90:232-245](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/seculr/transport_checks.f90#L232-L245).
- **What is wrong:** the safeguard enters only where Ω rises outward by more than 1 rad/s. The Sun rotates at 2.9e-6 rad/s, so it never fires, and reversals created by the transport go unnoticed. The 1e-15 thresholds elsewhere in the routine suggest that `1.0d0` was meant to be `1.0d-15`. The control `mu_gradient_scale` (`FMU`) is read and documented, but nothing uses it.
- **Expected behaviour, and why:** a safeguard should act on what it was written for, and a control that does nothing should not be offered as a physics knob.
- **Proposed fix:**
  - The safeguard needs the authors' decision. If it is wanted, use 1e-15 in a measured commit; it will then start forcing solid-body pairs. If not, delete it.
  - Drop `mu_gradient_scale` from the controls, or warn when it is set. This is byte-identical.

### O16. Settling rewrites the mass coordinate in place, and its sub-step count uses a reversed `MOD` (MIX)

- **Status:** open. The zero-count half of the `MOD` problem was guarded in L8.
- **In the original F77?** Yes. The settling setup multiplies the model's `HS1` by the unit conversion and divides it back afterwards (microdiff_setup.f:116-121, microdiff_etm.f:143-147; [microdiff_setup.f:113-121](https://github.com/parallelpro/yrec/blob/6cd56739/src/microdiff_setup.f#L113-L121)). The sub-step count is in MIX ([mix.f:412-415](https://github.com/parallelpro/yrec/blob/6cd56739/src/mix.f#L412-L415)):
  ```fortran
           NSTEP = INT(DELTS/DTMAX)
           IF(MOD(DTMAX,DELTS).NE.0.0D0)NSTEP=NSTEP+1
           DT = DELTS/DFLOAT(NSTEP)
  ```
- **Now:** [composition_step.f90:289-310](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/composition_step.f90#L289-L310), and the conversion at [microdiff_driver.f90:395-402](https://github.com/parallelpro/yrec/blob/e8549036/src/element_diffusion/microdiff_driver.f90#L395-L402) and [microdiff_mapping.f90:405-410](https://github.com/parallelpro/yrec/blob/e8549036/src/element_diffusion/microdiff_mapping.f90#L405-L410).
- **What is wrong:**
  - x·s/s need not equal x in floating point, so the structure's mass coordinate can change in its last bit on every settling sub-step.
  - The count tests `MOD(DTMAX, DELTS)` instead of `MOD(DELTS, DTMAX)`. When the timestep is an exact multiple of the settling limit, this takes one sub-step too many.
- **Expected behaviour, and why:** settling must not alter the structure's mass coordinate, so it should convert a copy. The count should be max(1, ⌈Δt/Δt_settle⌉).
- **Proposed fix:** convert local copies, and correct the `MOD` at both call sites (here and in the rotation sub-step), ideally in one helper. The effect is round-off drift, plus one fewer sub-step when the ratio is an exact integer.

### O17. Magnetic braking: the solid-body Matt factor lacks 1/G, the Matt law ignores the wind switch, and B-field scaling without Rossby scaling gives zero torque

- **Status:** open; `test_braking` pins the current behaviour.
- **In the original F77?** Yes, all three:
  - MCOWIND divides by `EXP(CLN*GL)`, where `GL` is an undeclared, never-set local (zero, so the factor is 1). MWIND uses `CGL` ([mcowind.f:107](https://github.com/parallelpro/yrec/blob/6cd56739/src/mcowind.f#L107), [mwind.f:86](https://github.com/parallelpro/yrec/blob/6cd56739/src/mwind.f#L86)).
  - MWIND's Matt branch never tests `LJDOT`, and SECULR calls it whatever `LJDOT0` is ([seculr.f:170-188](https://github.com/parallelpro/yrec/blob/6cd56739/src/seculr.f#L170-L188)).
  - `TCZ` is set only under `LROSSBY` but used under `LBSCALE` ([mwind.f:117-122](https://github.com/parallelpro/yrec/blob/6cd56739/src/mwind.f#L117-L122)).
  ```fortran
        FSUN = 0.5*PMMSOLW**2*CRSUN**3/EXP(CLN*GL)/CMSUN       ! mcowind.f
           IF(LBSCALE)THEN
              WNOW = W*TCZ/PMMSOLTAU                           ! mwind.f
  ```
- **Now:** [magnetic_braking.f90:198-206](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/magnetic_braking.f90#L198-L206) against [:250](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/magnetic_braking.f90#L250); [:100-107](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/magnetic_braking.f90#L100-L107); [:297-310](https://github.com/parallelpro/yrec/blob/e8549036/src/rotation/magnetic_braking.f90#L297-L310).
- **What is wrong:**
  1. The solid-body f_sun = ½Ω_sun²R_sun³/M_sun is about 7e-13 instead of 1.07e-5. The Sun's centrifugal factor is then c₂²/(c₂² + f) ≈ 1/1.004 instead of 1. By our estimate the solid-body Matt torque is low by about 1.004^m, with m the law's exponent.
  2. With the Matt law and the wind switched off, a decoupled surface convection zone or a solid-body star still loses J in the secular transport. The Kawaler branch returns correctly.
  3. With B-field scaling but no Rossby scaling, τ_cz stays 0, so the iterated torque is zero; `test_braking` confirms that no J is removed. The solid-body form ignores B-field scaling in this case.
- **Expected behaviour, and why:**
  1. f_sun = ½Ω_sun²R_sun³/(GM_sun), as in the other form, so that the factor is 1 for the Sun.
  2. With the wind off, neither law should remove J.
  3. B-field scaling needs τ_cz: compute it whenever the flag is set, or reject the combination at read time.
- **Proposed fix:** one measured commit each. Each affects only decks with that option: solid-body Matt, Matt with the wind off, and B-field scaling without Rossby scaling.

### O18. EOS and opacity details: two meanings of `beta14`, the opacity ramp's derivatives, OPAL 2001's density edge, and the untested MHD path

- **Status:** open.
- **In the original F77?** Yes for the first three; the missing test is a 2026 matter, since the F77 had no tests.
  - `BETA14` is 4(1−β)/β on the Yale/SCV path but 1−β after OPAL or MHD ([eqstat2.f:126-135](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L126-L135) against [:281-302](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L281-L302)):
    ```fortran
          BETA14 = 4.0D0*BETAI*DMIN1(1.0D0,BETA14)
    ...
                BETA14=1.0D0-BETA
    ```
  - The opacity ramp's derivatives are at [getopac.f:151-159](https://github.com/parallelpro/yrec/blob/6cd56739/src/getopac.f#L151-L159).
  - RHOOFP01's density-edge table differs from the loader's ([rhoofp01.f:23-26](https://github.com/parallelpro/yrec/blob/6cd56739/src/rhoofp01.f#L23-L26) `6*115,8*113` against [readcoeos01.f:33-36](https://github.com/parallelpro/yrec/blob/6cd56739/src/readcoeos01.f#L33-L36) `9*115,5*113`). OPAL's own code has the same mismatch.
- **Now:** [blended_eos.f90:314-324](https://github.com/parallelpro/yrec/blob/e8549036/src/eos/blended_eos.f90#L314-L324), [kap_lib.f90:326-338](https://github.com/parallelpro/yrec/blob/e8549036/src/kap/kap_lib.f90#L326-L338), [opal_eos.f90:276-290](https://github.com/parallelpro/yrec/blob/e8549036/src/eos/opal/opal_eos.f90#L276-L290) and [test_eos.f90:33-45](https://github.com/parallelpro/yrec/blob/e8549036/src/eos/test/test_eos.f90#L33-L45).
- **What is wrong:**
  - **`beta14`:** nothing outside the EOS uses it today, so its two meanings are a latent trap.
  - **The ramp:** κ = wκ_int + (1−w)κ_atm, but the dln κ derivatives are blended with w and 1−w instead of wκ_int/κ and (1−w)κ_atm/κ, and the (dw/dln T)(κ_int − κ_atm)/κ term is dropped. This affects only the Jacobian and any derivative output.
  - **The OPAL 2001 edge:** for the T6 rows 78-80 the density search stops at column 113, where the table is filled to 115. Pressures in between are reported as off the table (OPAL 2001 only).
  - **MHD:** its tables are not shipped, so `test_eos` runs the path only when `YREC_MHD_TABLES` is set.
- **Proposed fix:**
  - Give `beta14` one meaning, which is byte-identical outside the EOS.
  - Blend κ-weighted derivatives and add the dw term; this affects convergence only.
  - Use the loader's edge table, measured on OPAL 2001 decks.
  - Obtain the MHD tables and add a pin.

### O19. The Allard reader tests a misspelled flag, and the accretion surface uses a possibly failed Allard lookup (ALFILEIN, MASSLOSS)

- **Status:** open.
- **In the original F77?** Yes. ALFILEIN tests `LATMTPTau100`, which is never set; the real flag is `LALTPTau100`. MASSLOSS reads `AT` and `AP` after ALSURFP without testing `LAlFail` ([alfilein.f:102-112](https://github.com/parallelpro/yrec/blob/6cd56739/src/alfilein.f#L102-L112), [massloss.f:133-140](https://github.com/parallelpro/yrec/blob/6cd56739/src/massloss.f#L133-L140)):
  ```fortran
  C Ensure that we are not requesting PL,TL at Tau-100. These are not
  C present in the old Nextgen 1 files. If requested, fail
        If(LATMTPTau100) then
  ...
           CALL ALSURFP(TEFFL,GL,LPRT,LAlFail)
           T = 10.0D0**AT
  ```
- **Now:** [allard_atmosphere.f90:353-357](https://github.com/parallelpro/yrec/blob/e8549036/src/atm/tables/allard_atmosphere.f90#L353-L357) and [mass_change.f90:247-258](https://github.com/parallelpro/yrec/blob/e8549036/src/core/step/mass_change.f90#L247-L258).
- **What is wrong:** with τ = 100 requested and an old NextGen 1 file, the reader carries on with the τ = 100 tables still holding the −999 "invalid" marker, instead of stopping as its comment says. After a failed Allard lookup, the accretion surface state uses the previous atmosphere point.
- **Proposed fix:** test the real flag and return an error; check the failure flag and keep the lookup in locals. Neither changes a pinned deck.

### O20. The calibration's dL/dα is 0.0139, while its own comment records the fit as +0.139 (CHKCAL)

- **Status:** open; this needs the authors' decision.
- **In the original F77?** Yes ([chkcal.f:40-45](https://github.com/parallelpro/yrec/blob/6cd56739/src/chkcal.f#L40-L45)):
  ```fortran
        DLDX = -3.78D0      ! empirical result:  -3.783    RMS error .070
        DRDX = -0.89D0      ! empirical result:  -0.890    RMS error .048
        DLDA = 0.0139D0     ! empirical result:  +0.139    RMS error .0022
        DRDA = -0.050D0     ! empirical result:  -0.0504   RMS error .0059
  ```
- **Now:** [calibration.f90:72-93](https://github.com/parallelpro/yrec/blob/e8549036/src/core/run/calibration.f90#L72-L93).
- **What is wrong:** the value and the recorded fit differ by a factor of 10. The value enters the fixed 2×2 Newton Jacobian of the solar calibration. By our estimate the determinant is −0.0533 with 0.0139 and −0.0827 with 0.139, so each correction to α (and to X) differs by about 1.55×. The calibrated model does not depend on it, because the fixed point is log L = log R = 0; the iteration path and the number of iterations do.
- **Expected behaviour, and why:** the Jacobian should hold the fitted derivative. A closer Jacobian means fewer and more robust calibration iterations.
- **Proposed fix:** the authors say which number is right. Failing that, run the calibration deck with both values and compare iteration counts and the final α and X, which should agree within tolerance. Then correct whichever of the value or the comment is wrong.

### O21. The settling grid's interpolation stencil can read out of bounds for a 2-point grid

- **Status:** open; one of the two routines is guarded.
- **In the original F77?** Yes. The "JVS fix" catches only k0 = 0 ([microdiff_etm.f:55-68](https://github.com/parallelpro/yrec/blob/6cd56739/src/microdiff_etm.f#L55-L68), [equal_to_model.f:47-61](https://github.com/parallelpro/yrec/blob/6cd56739/src/equal_to_model.f#L47-L61)):
  ```fortran
               K0 = MAX(J-2,1)
  C  ENSURE THAT LAST INTERP. POINT NO GREATER THAN LAST EQUALLY SPACED POINT.
               K0 = MIN(K0,NPT-3)
  C JVS fix for NPT = 3?
               IF (K0 .EQ. 0) K0=1
  ```
- **Now:** the original routine is guarded (`k0 .lt. 1`, [settling_grid.f90:344-369](https://github.com/parallelpro/yrec/blob/e8549036/src/element_diffusion/settling_grid.f90#L344-L369), since 60788234). The newer routine, used by every deck with the new settling, still has the F77 test ([microdiff_mapping.f90:321-343](https://github.com/parallelpro/yrec/blob/e8549036/src/element_diffusion/microdiff_mapping.f90#L321-L343)).
- **What is wrong:** the grid size is ⌈Δr/Δr_min⌉, raised from 1 to 2, so a region of two equal intervals gets 2 points. Then k0 = −1, and the 4-point stencil reads `eq_radius(-1:2)`. Even with k0 = 1, a 2- or 3-point grid makes the stencil read points that this call never set.
- **Expected behaviour, and why:** interpolate only within the grid's own points (linear for 2, quadratic for 3), or guarantee at least 4.
- **Proposed fix:** enforce at least 4 points where both grids are built. This is byte-identical for regions that already have 4 or more, and changes only the small regions that misread now. Add a bounds-checked run to the battery (see O4).

### O22. Profile column 54 (N²) disagrees with the N² written to the pulsation files

- **Status:** open.
- **In the original F77?** No; the F77 wrote no N² or A4. WRTMOD's pulse model left N² to the pulsation code ([wrtmod.f:215-231](https://github.com/parallelpro/yrec/blob/6cd56739/src/wrtmod.f#L215-L231)). Column 54 was added by 0a75b9a8 (2026-08-24) as a raw centred difference. 2ed398ca and 67b36131 later changed only the pulse writers.
- **Now:** column 54 is at [stitched_model.f90:262-271](https://github.com/parallelpro/yrec/blob/e8549036/src/core/diagnostics/stitched_model.f90#L262-L271). The writers' N² is at [stitched_model.f90:693-731](https://github.com/parallelpro/yrec/blob/e8549036/src/core/diagnostics/stitched_model.f90#L693-L731) and [pulsation_output.f90:181-213](https://github.com/parallelpro/yrec/blob/e8549036/src/io/pulsation_output.f90#L181-L213).
- **What is wrong:** column 54 is g[(1/Γ₁) dln P/dr − dln ρ/dr] by centred differences everywhere. In convection zones that is cancellation noise of random sign: 67b36131 found 26 of 81 convective points spuriously positive on a subgiant. The pulse writers instead keep the thermal form g²ρδ(∇_ad − ∇)/P wherever it is negative. Column 54's header nevertheless claims the same content as FGONG's A4.
- **Expected behaviour, and why:** one N² everywhere: the thermal form where it is negative (exact for homogeneous mixed matter, and smooth), and the centred difference elsewhere, where the Ledoux μ term matters. Then column 54 shows what the pulsation codes receive. MESA takes a different route: it writes N² = (g²ρ/P)(χ_T/χ_ρ)(B − (∇ − ∇_ad)), with B the Brassard composition term formed from a difference at fixed T and ρ, which avoids the cancellation altogether. That is an option worth considering for both outputs.
- **Proposed fix:** apply the hybrid in the profile column through one function shared with the two writers. Only column 54 changes, and only in convection zones.

### O23. Two EOS derivatives are threaded through the Newton driver but carry nothing

- **Status:** fixed (2026-10-28), [6252b537](https://github.com/parallelpro/yrec/commit/6252b537): the derivatives are locals of `compute_henyey_coefficients`, and the arguments and star members are gone. Seeding the first shell's EOS with NaN left all 40 checked decks byte-identical, so the seed was dead, as stated below.
- **In the original F77?** Yes, structurally. MAIN passes `QDT` and `QDP` through CRRECT to COEFFT ([crrect.f:121-123](https://github.com/parallelpro/yrec/blob/6cd56739/src/crrect.f#L121-L123)), and STARIN seeds them as −1 and 1 ([starin.f:145-146](https://github.com/parallelpro/yrec/blob/6cd56739/src/starin.f#L145-L146)). EQSTAT lists them only as outputs ([eqstat.f:13-22](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat.f#L13-L22)).
- **Now:** [henyey_iteration.f90:134-157](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/henyey_iteration.f90#L134-L157) (`intent(out)` at line 155) and [structure_equations.f90:33-79](https://github.com/parallelpro/yrec/blob/e8549036/src/core/structure/structure_equations.f90#L33-L79).
- **What is wrong:** the coefficient pass reads an `intent(out)` dummy on entry, whose value is formally undefined. In practice nothing depends on it: every EOS family assigns these derivatives before using them. So `s%dlnrho_dlnp`, `s%dlnrho_dlnt` and the starting model's seeds are dead state.
- **Proposed fix:** make them locals of the coefficient pass and drop the arguments and the two star members. This is expected to be byte-identical.

### O24. An un-aged rotating model without instability transport is spun up by the wind

- **Status:** open (latent; found 2026-10-28 while tracing the negative `dt`). No pinned deck reaches it, and it has not been measured.
- **In the original F77?** Yes. MAIN marks a model that is solved but not aged with a negative timestep, `IF(.NOT.LAGE) DELTS = -DABS(DELTS)` ([main.f:644](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L644)), and calls GETW with it ([main.f:816](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L816)). GETW's branch without instability transport calls MWIND with `DELTS` whatever its sign ([getw.f:107-161](https://github.com/parallelpro/yrec/blob/6cd56739/src/getw.f#L107-L161)); only the transport branch skips a non-positive step ([getw.f:166-167](https://github.com/parallelpro/yrec/blob/6cd56739/src/getw.f#L166-L167)):
  ```fortran
            CALL MWIND(BL,DELTS,HSBOT,HSTOP,IMAX,M,LJDOT,OMEGAS,      ! getw.f:151, LINSTB = F
      IF(LSKIP.OR.DELTS.LE.0.0D0) GOTO 9999                         ! getw.f:167, LINSTB = T
  ```
- **Now:** since [88b32a90](https://github.com/parallelpro/yrec/commit/88b32a90) the aging flag replaces the negative `s%dt`, but the rotation call keeps the negated step for a model not aged, explicitly, so that the fix is a commit of its own ([step_driver.f90:538](https://github.com/parallelpro/yrec/blob/88b32a90/src/core/step/step_driver.f90#L538)); `spin_down_without_transport` hands it to `brake_convection_zone` ([evolve_rotation.f90:206-251](https://github.com/parallelpro/yrec/blob/4dd2a637/src/rotation/evolve_rotation.f90#L206-L251)), whose change of Ω is linear in the step ([magnetic_braking.f90:111-128](https://github.com/parallelpro/yrec/blob/4dd2a637/src/rotation/magnetic_braking.f90#L111-L128)).
- **What is wrong:** with a negative step the wind's ΔΩ changes sign, so the surface convection zone gains the angular momentum one step of wind would remove. It happens when all of these hold: rotation on, instability transport off, the wind torque on, no disk locking, a convective surface, and a model that is not aged (a rescale-only card, `run_kind = 2`, once log T_c ≥ 6.6, or a negative model number).
- **Expected behaviour, and why:** a model that is not aged spans no time, so the wind should remove nothing, as the transport branch already ensures.
- **Proposed fix:** pass `merge(s%dt, 0.0d0, evolve_model_flag)` to `evolve_angular_momentum`, as the Henyey pass already gets, so that the brake removes nothing over a zero step (or skip the wind in `spin_down_without_transport` for a non-positive step). One measured commit; only decks with the conditions above can move.

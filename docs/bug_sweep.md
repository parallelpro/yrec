# Bug sweep on yrec-next, September to October 2026

## 1. About this report

This report lists every defect found on the `yrec-next` branch after it diverged from `main` at [0ff84d7e](https://github.com/parallelpro/yrec/commit/0ff84d7e) (2026-09-03), up to [72de2124](https://github.com/parallelpro/yrec/commit/72de2124) (2026-10-29). There are 115 items: 95 fixed and 20 open. Bugs fixed before the divergence (bug-sweep batches 1-3, such as the CF88 deuterium rate) are in `docs/MODERNIZATION_CHANGES.md` and are not repeated here.

Each item answers the same questions:

- **In the original F77?** Each item was checked against commit [6cd56739](https://github.com/parallelpro/yrec/tree/6cd56739/src) ("First commit", 2025-04-10), the import of the Fortran 77 code. "Yes" means the same logic is there, with a link and a quote. "No" names the commit that introduced the defect: either a YREC-team commit (2025 or early 2026) or a 2026 modernization commit.
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
- Section 7's items, from the module evaluation of 2026-10-29, were verified on a 56-deck set against a snapshot; its introduction says how.
- Several fixes are protected by unit tests that fail if the fix is reverted: `test_structure` (every Henyey partial against finite differences), `test_mass_change`, `test_net`, `test_eos`, `test_triangle_fit` and `test_braking`.

**Conventions.**

- F77 links point at 6cd56739, and current-code links at e8549036.
- F77 routines are in capitals (CRRECT, SURFBC); current names are in code font (`compute_surface_boundary`).
- Numbers quoted from commit messages were measured at the time of the commit.

### Where the defects came from

| | In the F77 code | 2025-26 YREC-team commits | 2026 modernization | Total |
|---|---|---|---|---|
| Fixed | 72 | 4 (S5, L9, M21, M24) | 19 (N4, S13, L1, L4, L6, L7, L17, M7, M17, M18, M19, M22, M25, M26, M27, M35, M44, M45, M46) | 95 |
| Open | 16 | 1 (O4, with a 2026 part) | 3 (O5, O22, M33) | 20 |

In addition, half of O10 (an overshoot gate) is a 2026 translation slip. M10's unguarded divisions are the F77's, and its guard on the last interval alone a 2025 YREC-team commit; M33's branch is the F77's, and its stale values the port's.

### Decisions needed from the authors

These open items depend on what the code was meant to do, not only on what it does:

- **O20:** is the calibration's dL/dα 0.0139 or 0.139?
- **O1, O2, O3:** may the low-mass convergence work change the Jacobian? It would move every deck at the solver tolerance.
- **M32:** should `target_metal_fraction`, `envelope_generation_flag` and `pmm_norm_mdot`, which are read and used nowhere, be wired to what they name or retired?
- **M33:** should `use_diffusion_z` without Thoul's coefficients settle no metals (zero coefficients, as the F77 did) or be rejected at read time?
- **The EOS's R/μ at the Saha cutoff** (not a defect entry; `docs/ROADMAP.md`, module evaluation follow-up): below the cutoff the ramp sets R/μ = βP/(ρT), the effective value with the degeneracy correction, and at and above it R(1/μ_ion + 1/μ_e) from the composition; at log T = 6, log P = 13 they differ by 7e-4 (ρ, ∇_ad and c_P are continuous). Semiconvection's mean weights and the profile's μ read it. Which definition is meant?
- **The read-time envelope Z** (not a defect entry; `docs/ROADMAP.md`, module evaluation follow-up): the run plan sets the envelope Z, which the opacity setup reads before the starting model sets it again, to the last card's positive rescale Z, as PARMIN did ([parmin.f:1378](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L1378)). Is the last card's, rather than the first's, intended?

Answered with the bug list's proposed defaults on 2026-10-28, revisable: O12 (1e-14), O13 (the end cell includes an adjacent convection zone in its sums and its write-back), O15 (the safeguard removed), O25 (a convection zone reaching the centre is a solid body under walpcz ≤ -2).

Decided by default on 2026-10-29, revisable: M5 (the scales on the Caughlan & Fowler rates stay without effect, with a warning; making them apply would be byte-identical for the shipped decks), M12 (`ludcmp`'s singular tolerance, 1e-14 on the scaled pivot), M14 (`locate`'s slack, 1% of the table's span), M16 (`polint`'s coincidence, 4 ε of the larger distance), M17 (crmath's pow(0, y) = 0 for every y, MESA's definition), M22 (an unknown `pulse_format` stops the run only when pulse output is on), M27 (`hpttol`'s inlist name, `zoning_tolerances`), M28 (the settling latch re-armed at every run card, so at every Monte Carlo run too), after O15, `mu_gradient_scale` retired (legacy decks may still set `FMU`; the inlist format no longer offers it), M35 and M36 (the MHD EOS needs all eight tables, and a point outside them is an error), M37 (the Debye-Hückel metals from eos's own envelope mixture, the metals YREC's mixture lacks left out) and M40 (`in_atmosphere` skips only the Saha solver's second derivatives).

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
| O11 | Dynamical-shear limit can overshoot or go negative | Yes | Fixed |
| O12 | Solid-body Newton tolerance never set | Yes | Fixed (1e-14) |
| O13 | Transport region ends: sums and write-back disagree | Yes | Fixed |
| O14 | Species 5-15 unchecked after mixing; 16-slot history | Yes | Fixed |
| O15 | Reversal safeguard inactive; `mu_gradient_scale` inert | Yes | Fixed (removed; documented) |
| O16 | Settling rewrites m in place; reversed `MOD` | Yes | Fixed |
| O17 | Magnetic braking: 1/G, wind switch, τ_cz | Yes | Fixed |
| O18 | `beta14`, opacity ramp derivatives, OPAL 2001 edge, MHD test | Yes | Open |
| O19 | Allard reader's misspelled flag; failed lookup used | Yes | Open |
| O20 | Calibration dL/dα 0.0139 against the recorded 0.139 | Yes | Open |
| O21 | Settling stencil out of bounds for a 2-point grid | Yes | Open |
| O22 | Profile N² column disagrees with the pulse files | No (2026) | Open |
| O23 | Two EOS derivatives threaded through but unused | Yes | Fixed, [6252b537](https://github.com/parallelpro/yrec/commit/6252b537) |
| O24 | Un-aged rotating model without transport spun up by the wind | Yes | Fixed, [2bdbfbad](https://github.com/parallelpro/yrec/commit/2bdbfbad) |
| O25 | Constant j in a fully convective star diverges at the centre | Yes | Fixed |
| O26 | The history's total J and rotational KE lag their model by a step | Yes | Fixed (history columns only) |
| O27 | An inner region zeroes the wind term the surface region averages | Yes | Fixed |
| O28 | The transport check reads an unset convergence flag | Yes | Fixed (byte-identical), with the restructure of `apply_secular_transport` |
| O29 | The dynamical-shear limit takes the sign of the mean gravity | Yes | Open (beyond break-up only) |
| **M** | **Found in the module evaluation (§7)** | | |
| M1 | Total ε (profile eps_nuc, GYRE eps) left out the α captures | Yes | Fixed, [ee268284](https://github.com/parallelpro/yrec/commit/ee268284) |
| M2 | Capped deuterium rate kept the uncapped rate's derivatives | Yes | Fixed, [0ca5ff3e](https://github.com/parallelpro/yrec/commit/0ca5ff3e) |
| M3 | sneut's recombination derivatives missed dc₀₀/dν | Yes | Fixed, [85af8bef](https://github.com/parallelpro/yrec/commit/85af8bef) |
| M4 | Triple α's strong-screening z₅₃ was 9.104, not 9.014 | Yes | Fixed, [0c4fc520](https://github.com/parallelpro/yrec/commit/0c4fc520) |
| M5 | Scale entries that map to no reaction were silently inert | Yes | Warned, [77dfc919](https://github.com/parallelpro/yrec/commit/77dfc919) |
| M6 | Dead reaction slots 9 and 13; scale arrays a new reaction would alias | Yes | Fixed, [77dfc919](https://github.com/parallelpro/yrec/commit/77dfc919) |
| M7 | A `test_net` pin read the star's default initialization | No (2026) | Fixed, [104f3b80](https://github.com/parallelpro/yrec/commit/104f3b80) |
| M8 | bsstep checked none of its limits (`imax`, `nuse`) | Yes | Fixed, [33f6ca61](https://github.com/parallelpro/yrec/commit/33f6ca61) |
| M9 | osplin read out of bounds; meval's error dropped | Yes | Fixed, [c08569b0](https://github.com/parallelpro/yrec/commit/c08569b0) |
| M10 | osplin's divide-by-zero guard on the last interval only | Yes (guard 2025) | Fixed, [d4135f2d](https://github.com/parallelpro/yrec/commit/d4135f2d) |
| M11 | ysplin wrong for 3 points; OPAL92's dln κ/dln T at its top node | Yes | Fixed, [4eaaaece](https://github.com/parallelpro/yrec/commit/4eaaaece) |
| M12 | ludcmp passed singular matrices; a NaN column wrote a(0,k) | Yes | Fixed, [37296d25](https://github.com/parallelpro/yrec/commit/37296d25) |
| M13 | The tridiagonal solves never checked the first pivot | Yes | Fixed, [bcdcd649](https://github.com/parallelpro/yrec/commit/bcdcd649) |
| M14 | locate's 1% slack relative to the end value | Yes | Fixed, [a85bfe5e](https://github.com/parallelpro/yrec/commit/a85bfe5e) |
| M15 | lir failed silently | Yes | Fixed, [df8a8487](https://github.com/parallelpro/yrec/commit/df8a8487) |
| M16 | polint's coincidence test absolute | Yes | Fixed, [1c67bae5](https://github.com/parallelpro/yrec/commit/1c67bae5) |
| M17 | crmath's pow and exp10 at their edges | No (2026) | Fixed, [c2638e4b](https://github.com/parallelpro/yrec/commit/c2638e4b) |
| M18 | Monte Carlo model file kept only the last run's model | No (2026) | Fixed, [49e417d7](https://github.com/parallelpro/yrec/commit/49e417d7) |
| M19 | Run log labelled the central-D stop "log rho" | No (2026) | Fixed, [34fc2988](https://github.com/parallelpro/yrec/commit/34fc2988) (pins to reseed) |
| M20 | A rescale-and-evolve card ignored its stops | Yes | Fixed, [8c80883a](https://github.com/parallelpro/yrec/commit/8c80883a) |
| M21 | Most path controls' placeholders never expanded | No (2026, YREC team) | Fixed, [b1e7ae2c](https://github.com/parallelpro/yrec/commit/b1e7ae2c) |
| M22 | `pulse_format` case-sensitive at startup; unknown meant GYRE | No (2026) | Fixed, [ba392477](https://github.com/parallelpro/yrec/commit/ba392477) |
| M23 | The YREC7 and MODEL2 readers trusted the file | Yes | Fixed, [e2335b37](https://github.com/parallelpro/yrec/commit/e2335b37) |
| M24 | Output directory made by an unquoted shell command | No (2026, YREC team) | Fixed, [466b46ff](https://github.com/parallelpro/yrec/commit/466b46ff) |
| M25 | Output state in module variables; history file left open | No (2026) | Fixed, [226fabae](https://github.com/parallelpro/yrec/commit/226fabae) |
| M26 | A profile column missing in a region was a silent zero | No (2026) | Fixed, [770f6ce2](https://github.com/parallelpro/yrec/commit/770f6ce2) |
| M27 | Zoning tolerances zero without `use_new_rate_variables` | No (2026) | Fixed, [a0b4f434](https://github.com/parallelpro/yrec/commit/a0b4f434) |
| M28 | Settling switch overwritten as a run-long latch | Yes | Fixed, [af946c0c](https://github.com/parallelpro/yrec/commit/af946c0c) |
| M29 | Semiconvection error: run log only, flags mislabelled | Yes | Fixed, [db99ec9a](https://github.com/parallelpro/yrec/commit/db99ec9a) |
| M30 | Light-element-sum warning never printed its end | Yes | Fixed, [db99ec9a](https://github.com/parallelpro/yrec/commit/db99ec9a) |
| M35 | The MHD EOS read no tables | No (2026) | Fixed, [43e9b717](https://github.com/parallelpro/yrec/commit/43e9b717) |
| M36 | MHD returned stale results as success outside its tables | Yes | Fixed, [43e9b717](https://github.com/parallelpro/yrec/commit/43e9b717) |
| M37 | Debye-Hückel metals beyond C, N, O from the opacity tables | Yes | Fixed, [15af2ef3](https://github.com/parallelpro/yrec/commit/15af2ef3) |
| M38 | Diffusion of Z changed the EOS's H and He atomic weights | Yes | Fixed, [3db309e9](https://github.com/parallelpro/yrec/commit/3db309e9) |
| M39 | Saha ramp: two second derivatives summed, ionization fractions not blended | Yes | Fixed, [f8114dca](https://github.com/parallelpro/yrec/commit/f8114dca) |
| M40 | `in_atmosphere` left the Saha first derivatives stale | Yes | Fixed, [cd3b96ad](https://github.com/parallelpro/yrec/commit/cd3b96ad) |
| M41 | OPAL 2001/1995 density searches could start unbracketed | Yes | Fixed, [86d0b206](https://github.com/parallelpro/yrec/commit/86d0b206) |
| M42 | β ≤ 0 only warned | Yes | Fixed, [20d05da3](https://github.com/parallelpro/yrec/commit/20d05da3) |
| M43 | SCV envelope table used at a composition it was not built for | Yes | Fixed (latent), [72de2124](https://github.com/parallelpro/yrec/commit/72de2124) |
| M44 | A second `eos_init` asked about a freed unit number | No (2026) | Fixed, [33d0dc45](https://github.com/parallelpro/yrec/commit/33d0dc45) |
| M45 | `eos_get_gamma1` used derivatives it had asked not to compute | No (2026) | Fixed, [341d222a](https://github.com/parallelpro/yrec/commit/341d222a) |
| M46 | `test_eos`'s SCV values were Saha's (no envelope table) | No (2026) | Fixed, [c71c3ba0](https://github.com/parallelpro/yrec/commit/c71c3ba0) |
| M31 | `use_structure_dt_limits` off for the rest of the run | Yes | Open |
| M32 | Three controls stored and used nowhere | Yes | Open (documented) |
| M33 | `test_element_diffusion` reads stale metal coefficients | Branch yes; stale values no (2026) | Open |
| M34 | Monte Carlo summary prints zero derivatives after a first-check convergence | Yes | Open (cosmetic) |


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

- **Status:** fixed (2026-10-28). `compute_shear_velocities` takes wmin = Ω_i − (Ω_i − Ω_{i−1}) · limit/|dΩ/d ln r|, the Ω_{i−1} at which the trigger's own centred gradient is at the limit (strictly between Ω_i and Ω_{i−1} whenever the trigger fires), and scales U_GSF by the gradient term with wmin over the term with Ω_{i−1}, each in the form its own profile selects (`gsf_gradient_term`). `test_rotation_kernels` checks that on zones uneven in ln r, with the limit between the plain quotient and the centred gradient, the correction lowers U_GSF at every interface; the old wmin raised it at all 11. Measured on the 26 rotating decks: the five that report dynamical shear move, 21 are byte-identical; the testsuite's `Test_solarGS98_difrotmixper18p4d0p1` by 5e-5 dex in log L at most (1e-7 at 4.57 Gyr), the variants by up to 1e-4, and the solid-core variant rv_sbm2, which reports it 6438 times, by 3.6e-4 in log L and 4% in the surface omega.
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

- **Status:** fixed (2026-10-28): `default_inertia_tolerance` is 1e-14, near round-off (the bug's example was 1e-12; the smaller value keeps J within 1e-14 of the input, since the iteration ends by resetting j to I omega/dm). Measured on the 26 rotating decks: every one moves. In the fully convective configuration-matrix decks, which run few calls, by 2e-11 in log L and 3e-13 in J: the direct effect. Elsewhere the round-off differences are amplified, as by the marginal zones of O9, to at most 2e-4 dex in log L and 1e-3 in the surface omega (2-3% in the solid-core variant rv_sbm2), and the rescale-only deck, which stops at a model count, ends at another age. rv_ref runs 4% faster.
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

- **Status:** fixed (2026-10-28), by the proposal below, as the default the authors can revisit: the end cell includes an adjacent convection zone in both its sums and its write-back, as the sums (the F77 code's fuller treatment) always did. `equal_grid_to_model` and `diffuse_composition` widen the zones changed through the convection zones directly beyond the region, whether or not the end zone is convective itself. `test_rotation_kernels` checks a region whose first zone is radiative above a convection zone: the convection zone takes the region's first change of j and of each abundance (both caught under the old rule). No rotating deck reaches this case: the 26 are byte-identical; `test_rotation`'s synthetic transport case, whose region ends below the convective envelope at a radiative zone, changes (the envelope now takes the end zone's change). The composition's map back remains a spline without a rescale, so species are conserved only approximately there (2% in the test's coarse region).
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

- **Status:** fixed (2026-10-28). `apply_secular_transport` checks the species mixed after the iterations (5 to 11, or 15) with `check_mixed_species`, the repair-or-cut of slots 1-4 (`check_species_range`, now shared; a cut there is redone by the caller, which restores every tracked species); a rotating run's `max_diffusion_iters` outside 1..16 stops at read time (`derive_controls`, exit 1, `test_error_exit`). `test_rotation_kernels` checks the late species' repair and cut. The 26 rotating decks are byte-identical (none sets more than 16 iterations, and no late species goes anomalous in them). check_composition's last-iteration check of slots 5-11 before their mixing is kept, and harmless.
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

- **Status:** fixed (2026-10-28), as the default the authors can revisit: the safeguard is removed (git history keeps it; `check_angular_momentum`'s header says where it was), since it never acted and at the 1e-15 rad/s of its other tests it would act on the map back's interpolation noise. With it went reject_negative_j's unreachable cut-and-redo branch and `check_angular_momentum`'s unused `dt` argument. `mu_gradient_scale` stays readable (decks that set it still parse; legacy decks get FMU = 1 by default, so a warning would fire on every one) and its registry description says it has no effect. Byte-identical on the 26 rotating decks.
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

- **Status:** fixed (2026-10-28). The zero-count half of the `MOD` problem was guarded in L8.
- **What we did:** `microdiff` converts copies of the timestep, dlnP/dr, the mass coordinate and the total mass to Bahcall-Loeb units and discards them, so the model's own are never rescaled; `microdiff_etm` no longer converts back. Both settling sub-step counts (`settle` in the composition step, `burn_and_settle` in the rotation sub-step) take mod(timestep, max_settling_dt). Measured on the 50 self-contained decks of the checked set, each run alone: the 23 with settling move, 27 are byte-identical, as are the chained m0030, ZAHB and nuclear-rate decks. The copies alone account for all of it (the count changes only at an exact multiple, which no deck met): the testsuite's settling solar decks by 3e-5 dex in log L at 4.57 Gyr at most (1.6e-4 at any age), the rotating variants by up to 3.7e-4: round-off in the mass coordinate, amplified.
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

- **Status:** fixed (2026-10-28), one commit per item. `test_braking` pins the behaviour.
- **Item 3, what we did:** τ_cz is interpolated whenever B-field scaling is on, with or without Rossby scaling, in both forms; the solid-body form scales Ω by τ_cz/τ_cz,⊙ in that case too, with a fixed saturation. `test_rotation_kernels` checks that both forms brake and agree with B-field scaling alone (two mutations, one per form, are caught). No rotating deck uses B-field scaling without Rossby scaling: the 26 are byte-identical; `test_rotation`'s case now removes angular momentum.
- **Item 2, what we did:** `brake_convection_zone` returns without braking when the wind is off, for the Matt law as for Kawaler's. `test_rotation_kernels` checks both laws (a mutation that drops the test is caught). No rotating deck runs the Matt law with the wind off: the 26 are byte-identical; `test_rotation`'s wind-off case now keeps j.
- **Item 1, what we did:** `compute_solid_body_braking`'s solar centrifugal term divides by G, as `brake_convection_zone`'s does. `test_rotation_kernels` checks that the two forms give one star the same drop of omega over a short step (to 3e-6; 1.7e-3 apart before). Measured on the 26 rotating decks: the 14 Matt-law decks whose transport reaches the surface move (the Kawaler and solid-body decks do not use this form); the solid-body form's drop of omega rises by 0.17%, and J at 4.57 Gyr falls by 4e-4 in the testsuite's rotating solar decks (surface omega by 4e-4, log L by 2e-5 at most), by up to 1.6e-3 in the variants.
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

- **Status:** fixed (2026-10-28), [2bdbfbad](https://github.com/parallelpro/yrec/commit/2bdbfbad). Found the same day while tracing the negative `dt`. No pinned deck reaches it: the 40 checked decks are byte-identical before and after.
- **What we did:** a model that is not aged now gets a zero step in the rotation call, as the Henyey pass does, and `spin_down_without_transport` applies the wind only over a positive step. Measured on a deck built to reach the bug (the testsuite's solid-body solar deck with instability transport off and a rescale-only card at 0.1 Gyr, two un-aged models): before the fix the total angular momentum rose by 3.5% over the two un-aged models (2.042 to 2.114 × 10⁴⁹ g cm² s⁻¹) and the surface Ω by 3.5%; after it, J is unchanged to all printed digits. Five models later J is 3.2% lower and Ω 3.1% lower than before the fix, log L differs by 3 × 10⁻⁵ dex and log T_eff by 9 × 10⁻⁶.
- **In the original F77?** Yes. MAIN marks a model that is solved but not aged with a negative timestep, `IF(.NOT.LAGE) DELTS = -DABS(DELTS)` ([main.f:644](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L644)), and calls GETW with it ([main.f:816](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L816)). GETW's branch without instability transport calls MWIND with `DELTS` whatever its sign ([getw.f:107-161](https://github.com/parallelpro/yrec/blob/6cd56739/src/getw.f#L107-L161)); only the transport branch skips a non-positive step ([getw.f:166-167](https://github.com/parallelpro/yrec/blob/6cd56739/src/getw.f#L166-L167)):
  ```fortran
            CALL MWIND(BL,DELTS,HSBOT,HSTOP,IMAX,M,LJDOT,OMEGAS,      ! getw.f:151, LINSTB = F
      IF(LSKIP.OR.DELTS.LE.0.0D0) GOTO 9999                         ! getw.f:167, LINSTB = T
  ```
- **Now:** since [88b32a90](https://github.com/parallelpro/yrec/commit/88b32a90) the aging flag replaces the negative `s%dt`, but the rotation call keeps the negated step for a model not aged, explicitly, so that the fix is a commit of its own ([step_driver.f90:538](https://github.com/parallelpro/yrec/blob/88b32a90/src/core/step/step_driver.f90#L538)); `spin_down_without_transport` hands it to `brake_convection_zone` ([evolve_rotation.f90:206-251](https://github.com/parallelpro/yrec/blob/4dd2a637/src/rotation/evolve_rotation.f90#L206-L251)), whose change of Ω is linear in the step ([magnetic_braking.f90:111-128](https://github.com/parallelpro/yrec/blob/4dd2a637/src/rotation/magnetic_braking.f90#L111-L128)).
- **What is wrong:** with a negative step the wind's ΔΩ changes sign, so the surface convection zone gains the angular momentum one step of wind would remove. It happens when all of these hold: rotation on, instability transport off, the wind torque on, no disk locking, a convective surface, and a model that is not aged (a rescale-only card, `run_kind = 2`, once log T_c ≥ 6.6, or a negative model number).
- **Expected behaviour, and why:** a model that is not aged spans no time, so the wind should remove nothing, as the transport branch already ensures.
- **Proposed fix:** pass `merge(s%dt, 0.0d0, evolve_model_flag)` to `evolve_angular_momentum`, as the Henyey pass already gets, so that the brake removes nothing over a zero step (or skip the wind in `spin_down_without_transport` for a non-positive step). One measured commit; only decks with the conditions above can move.

### O25. Constant specific angular momentum in a fully convective star diverges at the centre

- **Status:** fixed (2026-10-28), with the first proposal below as the default the authors can revisit: under walpcz ≤ -2, a convection zone that reaches the centre rotates as a solid body, in every mode (`enforce_rotation_profile`). `test_rotation_kernels` checks a fully convective star under constant j: one omega, and its J kept (a mutation dropping the test is caught). The variant deck that stopped at its first model now runs; it stops at 1.67 Myr instead, once a small radiative core has formed and the constant-j surface convection zone spins its base beyond break-up (omega about 1e-3 rad/s, where the mean gravity turns negative; see O29). No shipped deck sets the law: the 26 rotating decks are byte-identical.
- **In the original F77?** Yes. WCZ imposes a uniform j on the convection zone when `WALPCZ .LE. -2` ([wcz.f:18-30](https://github.com/parallelpro/yrec/blob/6cd56739/src/wcz.f#L18-L30)):
  ```fortran
        ELSE IF(WALPCZ.LE.-2.0D0)THEN
  C  CONSTANT SPECIFIC ANGULAR MOMENTUM PER UNIT MASS IN THE C.Z.
  ```
- **Now:** [omega_profile.f90:401](https://github.com/parallelpro/yrec/blob/f0e8b4b6/src/rotation/omega_profile.f90#L401), `enforce_rotation_profile`'s constant-j branch.
- **What is wrong:** with uniform j, ω = j/(I/dm) ∝ r⁻². When the surface convection zone reaches the centre (a fully convective pre-main-sequence star), the innermost zone's ω is huge, the shape factors go non-physical, and the temperature gradient takes a square root of a negative number: the run stops at its first model with NaN. Measured: the testsuite's `Test_solarGS98_difrotmixper18p4d1p0` deck with `new_zone_omega_weight = -2.0` (walpcz) stops at once.
- **Expected behaviour, and why:** a rotation law should give a finite ω everywhere. Uniform j cannot hold down to r = 0.
- **Proposed fix (authors' choice):** treat a convection zone that reaches the centre as a solid body under this law (as zones below the surface already are), or reject walpcz ≤ -2 at read time with a message. Either is one measured commit; no shipped deck moves.

### O26. The history's total angular momentum and rotational kinetic energy lag their model by a step

- **Status:** fixed (2026-10-28). Found the same day by the sub-step convergence study (`src/tools/substep_convergence.py`).
- **What we did:** `evolve_angular_momentum` re-totals J and E_rot (`total_rotation_state`) for the j and omega it ends with, at each of its exits. Measured on the 26 rotating decks: only the history columns `total_angular_momentum` and `total_rotational_ke` change; every other column, and every pinned `.log` and `.mod`, is byte-identical. The history's J at the final model moves by -2.5% typically (one step's wind loss; -3e-4 to -2% for the testsuite's solar decks at 4.57 Gyr), and by up to 3-5% at any model. The MESA-format baseline of the rotating solar model, which pins history.data, is reseeded.
- **In the original F77?** Yes. HPOINT, the rezone at the start of the step, sets `SJTOT = SUMJ` ([hpoint.f:715-719](https://github.com/parallelpro/yrec/blob/6cd56739/src/hpoint.f#L715-L719)); GETW then moves and removes angular momentum, and nothing re-totals it before WRTOUT writes `SJTOT` to the track file ([wrtout.f:579](https://github.com/parallelpro/yrec/blob/6cd56739/src/wrtout.f#L579)).
- **Now:** `retally_angular_momentum` in the rezone sets `s%total_angular_momentum` and `s%total_rotational_ke` ([rezoning.f90:956-972](https://github.com/parallelpro/yrec/blob/e214d9cf/src/core/model/rezoning.f90#L956-L972)); `evolve_step` then runs the structure and the rotation step and writes the model ([step_driver.f90:150-226](https://github.com/parallelpro/yrec/blob/e214d9cf/src/core/step/step_driver.f90#L150-L226)) without re-totalling; the history columns `total_angular_momentum` and `total_rotational_ke` read the stale values ([history_output.f90:206-207](https://github.com/parallelpro/yrec/blob/e214d9cf/src/io/history_output.f90#L206-L207)).
- **What is wrong:** each history row's J is the start of its step's (after the rezone), the end of the previous model's, while omega and the rest are the end of the step's. Seen in the rotating reference variant: the last step to the 300 Myr stop, cut to a third of the previous step, loses the full previous step's 3.3% of J in the history.
- **Expected behaviour, and why:** a row describes one model; its J should be the sum of j dm of the model it prints.
- **Proposed fix:** re-total J and E_rot after the rotation step, before the output. Only the two history columns (and the rezone's diagnostic line, whose "previous" becomes the end of the last step) change; nothing computes with the totals.

### O27. An inner transport region zeroes the wind term that the surface region averages

- **Status:** fixed (2026-10-28). Found the same day while restructuring `apply_secular_transport`. Instrumented, the rotating decks (the rotating reference variants, the testsuite's rotating solar decks, the 3 Msun TAMS deck) never had an inner region before the surface region in an iteration after the first.
- **What we did:** the sub-step's `surface_wind_state` keeps the surface region's implicit term (`surface_implicit`, with `surface_implicit_set`) apart from the per-region terms the solve gets; `outer_point_wind_terms` averages with it only when an earlier iteration of the sub-step set it, and otherwise starts afresh (it passes iteration 1 to `compute_solid_body_braking`). Measured on the 26 rotating decks: 25 byte-identical; the GSF-inhibition-mode-2 variant moves (log L by 4.4e-7 at most, surface omega by 0.5% at most and 5e-4 at 300 Myr, J by 4e-4): there the surface region first appears in a later iteration of some sub-steps, whose implicit term was averaged with zero.
- **In the original F77?** Yes. SECULR sets `WIND1 = WIND2 = 0` for every region that does not reach the surface ([seculr.f:290-293](https://github.com/parallelpro/yrec/blob/6cd56739/src/seculr.f#L290-L293)); MCOWIND, for the surface region, averages `WIND2` with its value on entry from the second iteration on ([mcowind.f:121-125](https://github.com/parallelpro/yrec/blob/6cd56739/src/mcowind.f#L121-L125)).
- **Now:** `set_surface_wind_terms` in [secular_transport.f90](https://github.com/parallelpro/yrec/blob/e214d9cf/src/rotation/seculr/secular_transport.f90), and `compute_solid_body_braking`'s `domega_end` average ([magnetic_braking.f90](https://github.com/parallelpro/yrec/blob/e214d9cf/src/rotation/magnetic_braking.f90)).
- **What is wrong:** the implicit wind term is meant to be averaged with the previous iteration's (MHP's damping of the iteration). When another region lies below the surface region in the same iteration, it has just set the term to zero, so the average is with zero: the surface's implicit loss is halved in that iteration. The same holds when the first iteration had no surface region: the average is then with whatever the term held.
- **Expected behaviour, and why:** the surface region's implicit term should average with the surface region's previous iterate, or start afresh when there is none.
- **Proposed fix:** keep the surface region's previous implicit term apart from the per-region terms the solve gets, and average only with a value from an earlier iteration of the same sub-step. Byte-identical on the shipped decks.

### O28. The angular momentum check reads an unset convergence flag

- **Status:** fixed (2026-10-28) by the restructure of `apply_secular_transport` into module procedures, which passes false (`check_iterate`). Byte-identical on the 56 checked decks: the build's `-finit-local-zero` had made the unset flag false.
- **In the original F77?** Yes. SECULR sets `LOKAD = .FALSE.` only inside a region's angular momentum solve ([seculr.f:295](https://github.com/parallelpro/yrec/blob/6cd56739/src/seculr.f#L295)) and passes it to CHECKJ after the regions ([seculr.f:342](https://github.com/parallelpro/yrec/blob/6cd56739/src/seculr.f#L342)).
- **Now:** `diffusion_converged`, set only in `transport_angular_momentum_in_region` and read by `check_iterate` ([secular_transport.f90](https://github.com/parallelpro/yrec/blob/e214d9cf/src/rotation/seculr/secular_transport.f90)).
- **What is wrong:** under solid-body rotation (`force_solid_body_rotation` or `solid_body_mode_flag = 1`), and in an iteration whose regions are all single interfaces, no region solves for angular momentum, so the flag is never set; `check_angular_momentum` reads it, and a true value would mark the iteration converged and suppress its report. The solid-body decks print the report every sub-step, so the stack has held false.
- **Expected behaviour, and why:** the flag belongs to the commented-out advection-diffusion solve, which could report its own convergence; without that solve it is always false.
- **Proposed fix:** pass false. Byte-identical.

### O29. The dynamical-shear limit takes the sign of the mean gravity

- **Status:** open (found 2026-10-28 by the walpcz = -2 variant deck after O25's fix). Only a model rotating beyond break-up reaches it.
- **In the original F77?** Yes. VCIRC's QWRMAX is 2 g sqrt((rho/P) max(grad_ad - grad, 1e-3)) with the signed mean gravity.
- **Now:** `rotation_stability_setup`'s `stab%dynamical_shear_omega_limit` ([stability_coefficients.f90](https://github.com/parallelpro/yrec/blob/e214d9cf/src/rotation/seculr/stability_coefficients.f90)).
- **What is wrong:** the Richardson number holds g^2, so the largest stable |d omega/d ln r| is 2 |g| sqrt(...). Where the centrifugal term exceeds gravity the mean gravity is negative, the limit is negative, every interface reads as dynamically unstable, and the shear's correction of U_GSF takes an inner omega outside the two zones' values (seen: -2.9e-2 rad/s). The same run then ends with an integration error (BSSTEP) and exit status 0, where the error contract asks for 1.
- **Expected behaviour, and why:** a model beyond break-up should stop with a message; the limit, where it is used, should take |g|.
- **Proposed fix:** take abs(g) in the limit, and stop the run (ierr) when the mean gravity of a rotating model turns negative. One measured commit; only models beyond break-up change.

## 7. Found in the module evaluation (2026-10-29)

On 2026-10-29 five read-only reviews scored eos/, net/, io/, numerics/ with math/, and state/ with setup/ against the standard rotation/ set after its Stage D ([audit/module-evaluation-2026-10-29](https://github.com/parallelpro/yrec/blob/908b983c/audit/module-evaluation-2026-10-29/SUMMARY.md), [908b983c](https://github.com/parallelpro/yrec/commit/908b983c)). The work that followed, one package per area, fixed the bugs the reviews reported and found more while writing the new tests; it is described in `docs/ROADMAP.md` ("Module evaluation follow-up"). These items come from it.

- M1-M7 are the network's, M8-M17 the numerics' and math's, M18-M26 the output's, model files' and run files', M27-M30 the controls', and M35-M46 the EOS's (placed after M30). M31-M34 are open.
- **How they were verified.** Every commit was checked by a byte comparison of every output file of a 56-deck set (the configuration matrix, the solar and testsuite decks, the pre-main-sequence chain, the rotating variants) against a snapshot of [2c88d896](https://github.com/parallelpro/yrec/commit/2c88d896), the run logs also against the pins, and by the unit tests. A fix that changes results is its own commit, and its message records which decks move and by how much; where no deck of the set reaches the defect, a deck built to reach it was measured. Once an EOS fix had moved every deck, the EOS commits were compared on 50 decks with the previous commit's binary instead (every output file, and the drift from history.data). The pinned files in `standard/` were reseeded once, after the series, at [2f84ae0e](https://github.com/parallelpro/yrec/commit/2f84ae0e): for M19's run-log line, the precision change of the rotating decks' defaults ([9b58f6f6](https://github.com/parallelpro/yrec/commit/9b58f6f6), in the ROADMAP) and M37-M39.
- Current-code links point at [396bbb7c](https://github.com/parallelpro/yrec/tree/396bbb7c/src), and the EOS items' at [a568e084](https://github.com/parallelpro/yrec/tree/a568e084/src), the end of the series.

### M1. The total ε (the profile's eps_nuc and GYRE's eps) left out the α captures

- **Status:** fixed in [ee268284](https://github.com/parallelpro/yrec/commit/ee268284) (2026-10-29).
- **In the original F77?** Yes. ENGEB returns the five chain energies EPP1, EPP2, EPP3, ECN and E3AL as arguments and keeps the α captures in EALPCA ([engeb.f:1043-1051](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L1043-L1051)); COEFFT's per-shell total SESUM sums the five ([coefft.f:215-225](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L215-L225), [:347](https://github.com/parallelpro/yrec/blob/6cd56739/src/coefft.f#L347)), while ENGEB's derivatives sum every reaction ([engeb.f:942-954](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L942-L954)):
  ```fortran
  C ENERGY FROM ALPHA CAPTURE REACTIONS.
        EALPCA=DG(8)+DG(10)+DG(11)                        ! engeb.f
              SESUM(IM) = EG(1)+EG(2)+EG(3)+EG(4)+EG(5)   ! coefft.f
  ```
- **Now:** `store_zone_energy_channels`, [structure_equations.f90:240-265](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/structure/structure_equations.f90#L240-L265), sums `eps_by_chain` over every chain.
- **What was wrong:** `eps_total`, from which the profile's eps_nuc, GYRE's eps and rotation's interpolated ε (for the Eddington-Sweet velocities) are taken, summed the pp I-III, CNO and triple-α chains. ¹²C(α,γ), ¹³C(α,n) and ¹⁴N(α,γ) went only into the luminosity breakdown (`i_lum_he_c`). GYRE's eps_eps_T is the Henyey dε/dln T of the full ε, so the two columns disagreed wherever helium burns.
- **Expected behaviour, and why:** ε is the energy released per gram by every reaction. A pulse file's ε and dε/dln T must describe the same ε, since GYRE's nonadiabatic equations use both.
- **What we did:** `eps_total` is now the sum over every chain; the ε channels stay the five chains' fractions of it. The chain energies had become one array, `eps_by_chain`, in the commit before ([2cba8030](https://github.com/parallelpro/yrec/commit/2cba8030)).
  - The full set is byte-identical in every output file: hydrogen-burning zones carry no α captures, and no rotating deck reaches helium burning.
  - The 1 Msun ZAHB deck with profiles and GYRE files on: history, track, model and every default profile column byte-identical. GYRE's central ε is unchanged at the ZAHB (7144 erg/g/s), goes from 3361 to 5776 erg/g/s in mid core-helium burning and from 31.6 to 4217 near its end, where ¹²C(α,γ) dominates; there eps_eps_T/eps goes from 2230 to 16.7.

### M2. The capped deuterium rate kept the uncapped rate's derivatives

- **Status:** fixed in [0ca5ff3e](https://github.com/parallelpro/yrec/commit/0ca5ff3e) (2026-10-29).
- **In the original F77?** Yes. ENGEB caps the d(p,γ) rate in the surface convection zone (G. Somers 3/17) and then adds the deuterium energy to both derivative sums as if uncapped ([engeb.f:523-535](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L523-L535), [:960-962](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L960-L962)). The fix's commit message says the cap is not in the F77; it is.
  ```fortran
          IF(IU.GE.JCZ .AND. TAUCZ.GT.1.0D0)THEN
             RDEUTMAX = 6.023D23/ANUC(3)/TAUCZ
  ...
                   RDEUT = RDEUTMAX/X
        SUM1 = SUM1 + DGDEUT
        SUM2 = SUM2 + DGDEUT
        SUM3 = SUM3 + DGDEUT*QRTDEUT
  ```
- **Now:** `compute_deuterium_burning`, [nuclear_energy.f90:143-167](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/nuclear_energy.f90#L143-L167), with the cap in `cap_deuterium_rate`, [reaction_rates.f90:156-171](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/reaction_rates.f90#L156-L171).
- **What was wrong:** where the cap binds, deuterium burns no faster than the convective overturn: ε_D = (N_A/m_d/τ) X_D Q, which depends on neither ρ nor T. The energy sum still added ε_D to dε/dln ρ and ε_D times the fit's dln r/dln T to dε/dln T.
- **Expected behaviour, and why:** the derivatives must be those of the rate actually used, zero for a capped rate. Only the Henyey Jacobian was wrong, which slows the Newton iteration but does not move the converged model.
- **What we did:** `compute_deuterium_burning` returns dln ε/dln ρ as well: 1 and the fit's dln r/dln T uncapped, 0 and 0 where the cap binds. `test_net` checks the capped point against centred differences.
  - The three pre-MS decks (0.3 Msun) are byte-identical in every output file, history.data included: their deuterium burns near log T 6 in the core, four orders of magnitude below the cap. The full set is byte-identical too.

### M3. sneut's recombination-neutrino derivatives missed dc₀₀/dν

- **Status:** fixed in [85af8bef](https://github.com/parallelpro/yrec/commit/85af8bef) (2026-10-29). Found by the new `test_net` check of the Itoh losses.
- **In the original F77?** Yes, in Timmes' sneut ([sneut.f:1247-1251](https://github.com/parallelpro/yrec/blob/6cd56739/src/sneut.f#L1247-L1251)):
  ```fortran
         c00    = 1.0d0/(1.0d0 + f1*nu + f2*nu2 + f3*nu3)
         c01    = f1 + f2*2.0d0*nu + f3*3.0d0*nu2
         dum    = zeta*c00
         dumdt  = zetadt*c00 + zeta*c01*nudt
  ```
- **Now:** [thermal_neutrinos.f90:1339-1345](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/thermal_neutrinos.f90#L1339-L1345).
- **What was wrong:** in the recombination term (Itoh et al. 1996, eqs. 6.5-6.14) the derivatives used `c01`, the derivative of the polynomial, as if it were dc₀₀/dν; it lacks the factor −c₀₀². For ν ≥ 0 (the f's are zero below) the recombination term's own derivatives were off by up to 28%: at log T 8.6, log ρ 5 the total dε/dln ρ was 1.2e-3 off and dε/dln T 2.5e-5. The 2026 bug sweep had fixed another of sneut's typos (`tfac2`).
- **Expected behaviour, and why:** d(ζ c₀₀)/dν = ζ dc₀₀/dν with dc₀₀/dν = −c₀₀² (f₁ + 2f₂ν + 3f₃ν²).
- **What we did:** `c01 = -c00*c00*(f1 + f2*2.0d0*nu + f3*3.0d0*nu2)`; the check's tolerance is 1e-5 (the bremsstrahlung derivative is good to 2e-6 there). Only the Henyey Jacobian, and the pulse output's ε derivatives, of runs with `use_itoh_neutrino_loss` see it.
  - The full set is byte-identical: none of its decks uses the Itoh losses past the main sequence.
  - The 1 Msun ZAHB deck with the Itoh losses on: the age-matched drift in log L is 3.6e-9 over the first quarter of the run, 5.7e-7 to half and 2.8e-4 to three quarters (Y_c 6.6e-4); near core-helium exhaustion the runs part (1575 → 1656 models, final age 2.321e8 → 2.378e8 yr, final Δlog L −2.9e-3, Δlog Teff +1.5e-3). Only the Jacobian changes, so each converged model moves within the Newton tolerance, and the late core-helium-burning phase (convective-core growth, breathing pulses) amplifies such differences.

### M4. Triple α's strong-screening constant z₅₃ was 9.104, not 9.014

- **Status:** fixed in [0c4fc520](https://github.com/parallelpro/yrec/commit/0c4fc520) (2026-10-29). Found while naming the screening constants.
- **In the original F77?** Yes; the twelfth entry of `Z53` is triple α's ([engeb.f:393-395](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L393-L395)):
  ```fortran
        DATA ZPRD/1.,4.,4.,6.,6.,7.,8.,12.,16.,12.,14.,12.,36./,Z53/1.175,
       1   3.73,3.73,4.804,4.804,5.385,5.941,9.014,11.24,9.014,10.15,
       2   9.104,23.28/,Z43/0.52,1.31,1.31,1.488,1.488,1.61,1.721,2.577,
  ```
- **Now:** the triple-α row of the reaction table, [network.f90:242-246](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/network.f90#L242-L246).
- **What was wrong:** Graboske et al.'s (1973, table 4) constants are z_p = (Z₁+Z₂)^p − Z₁^p − Z₂^p for p = 1.86, 5/3, 4/3 and 2/3. Every row agrees with that to the table's rounding except triple α's z₅₃: Z = 2 and 6 give 9.014, the value ¹²C(α,γ) and ¹³C(α,n), with the same charges, carry. A digit transposition.
- **Expected behaviour, and why:** 9.014. z₅₃ enters only the strong form, so only strongly screened triple α moves: in pure helium at T9 = 0.08 the strong form applies from ρ ~ 2.5e6 g/cm³, and there the old value made triple α about 2.5% too fast (4% at 1e7).
- **What we did:** 9.014; `test_net` compares triple α's strong exponent with Graboske et al.'s formula at log T 7.8, log ρ 8, which catches the old value. The full set and the 1 Msun ZAHB deck with the Itoh losses are byte-identical: no deck screens triple α strongly. The degenerate helium cores at the flash of a 1 Msun star (ρ_c ~ 1e6) sit at the edge of that regime.

### M5. The scale-vector entries that map to no reaction were silently inert

- **Status:** warned in [77dfc919](https://github.com/parallelpro/yrec/commit/77dfc919) (2026-10-29); the behaviour is unchanged, by decision (below).
- **In the original F77?** Yes. ENGEB multiplies by `SSTANDARD(I)` only the seven S-factor rates and the pep, ⁷Be and hep terms ([engeb.f:678-685](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L678-L685), [:723-725](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L723-L725)); entries 8-13 are never read, and its header says why ([engeb.f:117-120](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L117-L120)):
  ```fortran
  C  REACTIONS J = 8, 13 ARE NOT RELEVANT FOR THE SOLAR INTERIOR; THEY ARE
  C   HOLDOVERS FROM THE EARLIER YALE CODE.
        DO 40 I=1,7
           RATE(I) = RATE(I)*SSTANDARD(I)
  ```
- **Now:** `cross_section_scales_from_sstandard`, [net_configuration.f90:108-130](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/net_configuration.f90#L108-L130), and `set_net_settings`, [network.f90:264-295](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/network.f90#L264-L295).
- **What was wrong:** with `use_new_nuclear_rates = .false.`, `sstandard(8..13)` went into `cross_section_scale` and nothing used them: every Caughlan & Fowler rate had no scale slot, and entries 9 (¹⁶O+α) and 13 (¹²C+¹²C) have no reaction at all. A deck setting `sstandard(10)` to scale ¹²C(α,γ) got no effect and no message.
- **Expected behaviour, and why:** a control that does nothing should say so.
- **What we did:** the legacy vector goes through `cross_section_scales_from_sstandard`, which warns (net's message, terminal and run log) when entry 9 or 13 is not 1; `set_net_settings` warns at configure time when a scale is set on a reaction whose rate law takes none (`inert_scale_count`). Making the scales apply to the Caughlan & Fowler rates is one factor in each rate law and byte-identical for every shipped deck (all have 1 there); it is listed under "Decisions needed". Byte-identical: the full set.

### M6. Two dead reaction slots, and scale arrays that a new reaction would alias

- **Status:** fixed in [77dfc919](https://github.com/parallelpro/yrec/commit/77dfc919) and [fa8cdec2](https://github.com/parallelpro/yrec/commit/fa8cdec2) (2026-10-29).
- **In the original F77?** Yes. ENGEB carries reactions 9 and 13 and zeroes them after the energy sums are formed ([engeb.f:931-938](https://github.com/parallelpro/yrec/blob/6cd56739/src/engeb.f#L931-L938)); the derivative scales `QS0E(8)` and `QQS0EE(8)` hold the seven S-factor reactions and, in slot 8, ⁷Be+p ([remap.f:25](https://github.com/parallelpro/yrec/blob/6cd56739/src/remap.f#L25), [:163](https://github.com/parallelpro/yrec/blob/6cd56739/src/remap.f#L163)):
  ```fortran
  C SET TO ZERO O16+ALPHA AND C12+C12 RATES.
        DRATRO(9) = 0.0
        DRATT(9) = 0.0
        DG(9) = 0.0
  ```
- **Now:** the layout in [net_configuration.f90:28-30](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/net_configuration.f90#L28-L30) (`r_*`, `num_reactions = 11`) and the table's self-check, `network_table_errors`, [network.f90:418-441](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/network.f90#L418-L441).
- **What was wrong:** slots 9 and 13 were carried, zeroed in the energy derivatives and printed, so a reaction placed in either would have lost its Jacobian silently. In the 2026 layout `num_reactions = 14` would have aliased the pep scale (slot 14 of the same vector), and a new S-factor reaction had no derivative-scale slot (slot 8 was ⁷Be+p's).
- **Expected behaviour, and why:** each reaction has its own slot in each scale array, and every slot is computed by exactly one table row.
- **What we did:** one layout in `net_configuration` (11 reaction slots, the extra scales `xs_pep`, `xs_be7e`, `xs_be7p`, `xs_hep` by name, per-slot derivative scales with ⁷Be+p's in their own fields), and a self-check of the table, run by `set_net_settings` at configure time with `ierr`: every slot is computed by exactly one row, a shared rate has a row that computes it, and the rate law, chain, branch, burning group and species exist. Dropping the dead slots removes summed zeros, an exact change: the full set and `test_net`'s values are byte-identical.

### M7. A `test_net` pin depended on the star's default initialization

- **Status:** fixed in [104f3b80](https://github.com/parallelpro/yrec/commit/104f3b80) (2026-10-29).
- **In the original F77?** No; the test came with [06b40060](https://github.com/parallelpro/yrec/commit/06b40060) (2026-10-10).
- **What was wrong:** the pinned `deuterium_rate(1.5d0, 6.2d0, 0.70d0, 1.ge.star%jcz, star%convective_turnover_timescale)` took the zone flag and the turnover time from whatever `star_info`'s default initialization left there.
- **What we did:** the call passes explicit outside-the-zone arguments, which give the same value ([test_net.f90:249](https://github.com/parallelpro/yrec/blob/396bbb7c/src/net/test/test_net.f90#L249)).

### M8. bsstep checked none of its limits; the inlist controls imax and nuse could write out of bounds

- **Status:** fixed in [33f6ca61](https://github.com/parallelpro/yrec/commit/33f6ca61) (2026-10-29), with the read-time checks in [97c00d5a](https://github.com/parallelpro/yrec/commit/97c00d5a) and [396bbb7c](https://github.com/parallelpro/yrec/commit/396bbb7c).
- **In the original F77?** Yes. BSSTEP's substep sequence has 11 entries and RATEXT's tableau 15 rows and 7 columns, while `IMAX` and `NUSE` come from the PHYSICS namelist unchecked ([bsstep.f:13-17](https://github.com/parallelpro/yrec/blob/6cd56739/src/bsstep.f#L13-L17), [:39-49](https://github.com/parallelpro/yrec/blob/6cd56739/src/bsstep.f#L39-L49), [ratext.f:7-8](https://github.com/parallelpro/yrec/blob/6cd56739/src/ratext.f#L7-L8), [parmin.f:402](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L402)):
  ```fortran
        COMMON/INTPAR/STOLR0,IMAX,NUSE
       *          YSEQ(3),NSEQ(11),ERR(3)
              HNEXT = H*DFLOAT(NSEQ(NUSE-1))/DFLOAT(NSEQ(I))
        H = 0.25D0*H/2.0D0**INT((IMAX-NUSE)/2)
  ```
- **Now:** `bsstep`, [bulirsch_stoer.f90:81-168](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/bulirsch_stoer.f90#L81-L168); `check_controls`, [controls_check.f90:75-118](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/controls_check.f90#L75-L118).
- **What was wrong:** `imax = 12` or `nuse = 8` wrote past the substep sequence, the abscissae and the tableau; 16 or more equations wrote into the next tableau column; `nuse = 1` read `substep_sequence(0)`; and `nuse > imax + 1` made the retry factor 0.25/2^((imax−nuse)/2) keep or grow h, so the step never underflowed: an endless loop.
- **Expected behaviour, and why:** a solver must refuse arguments its arrays cannot hold, and a step-size retry must shrink the step.
- **What we did:** `bsstep` returns `numerics_bad_input`, before any derivative call, unless 1 ≤ num_eqs ≤ 15, 2 ≤ extrap_order ≤ 7 and extrap_order ≤ max_stage_index ≤ 11 (the named `bsstep_max_equations`, `bsstep_max_extrap_order`, `bsstep_max_stages`); the envelope reports it as a configuration error. The controls read rejects the same limits (`imax` 1 to 11, `nuse` 2 to 7, `nuse` ≤ `imax`, a positive `stolr0`) with exit 1 and a message naming them, from the same parameters. Every shipped deck has 11 and 7: the full set is byte-identical.

### M9. osplin read out of bounds off the table, and discarded meval's error

- **Status:** fixed in [c08569b0](https://github.com/parallelpro/yrec/commit/c08569b0) (2026-10-29).
- **In the original F77?** Yes, in ACM 574's MEVAL (the YREC tail fix of LLP 8/19/08 covered one case). With every point below the table, `START` runs to `NUME+1` and the search reads `XVAL(START)`; OSPLIN never reads `ERR` ([meval.f:99-102](https://github.com/parallelpro/yrec/blob/6cd56739/src/meval.f#L99-L102), [:136](https://github.com/parallelpro/yrec/blob/6cd56739/src/meval.f#L136), [osplin.f:41](https://github.com/parallelpro/yrec/blob/6cd56739/src/osplin.f#L41)):
  ```fortran
    20  DO 30 I=1,NUME
            IF (XVAL(I) .GE. XTAB(1)) GO TO 40
            START=I+1
        CALL SEARCH(XTAB,NUM,XVAL(START),LCN,FND)
        CALL MEVAL(XVAL,YVAL,XTAB,YTAB,M,N,K,EPS,ERR)
  ```
- **Now:** `osplin` and `meval`, [quadratic_spline.f90:44-63](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/quadratic_spline.f90#L44-L63), [:466-627](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/quadratic_spline.f90#L466-L627).
- **What was wrong:** two or more points all below the table read `eval_x(k+1)` (and could write `eval_y(k+1)`); two or more all above read `eval_x(0)`; a point on the last node followed by points above read `table_x(n+1)` and on, as did a run that reached the last node from an interior one; n < 3 read `table_y(3)`. meval's err = 2 (decreasing points) was dropped, leaving the values as they were with no flag.
- **What we did:** the in-range part runs only when some point is in range, and the loop that advances over the nodes stops after the last point in range. osplin has an optional `ierr`: `numerics_bad_input` for n < 3, `numerics_unsorted_points` for decreasing points, either leaving `yval` as it was (`intent(inout)`, so the declaration says so). Every case that stayed in bounds takes the same branches with the same arithmetic: the full set is byte-identical. `test_numerics` checks that a run of points equals its points taken one at a time, across these cases.

### M10. osplin's divide-by-zero guard covered only the last interval

- **Status:** fixed in [d4135f2d](https://github.com/parallelpro/yrec/commit/d4135f2d) (2026-10-29).
- **In the original F77?** Yes: SLOPES divides by every interval's width unguarded ([slopes.f:49-53](https://github.com/parallelpro/yrec/blob/6cd56739/src/slopes.f#L49-L53), [:83-86](https://github.com/parallelpro/yrec/blob/6cd56739/src/slopes.f#L83-L86)). The YREC-team commit [a093427f](https://github.com/parallelpro/yrec/commit/a093427f) (2025-05-31, "avoided floating point exception in slopes.f") guarded the last secant and the last end slope.
  ```fortran
        M1=YDIF1/(XTAB(2) - XTAB(1))
        M2=YDIF2/(XTAB(3)-XTAB(2))
  ```
- **Now:** `slopes` and `safedivide`, [quadratic_spline.f90:80-157](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/quadratic_spline.f90#L80-L157).
- **What was wrong:** a repeated abscissa inside the table gave infinite or NaN slopes, which spread to the neighbouring intervals, while the header said that the divisions by zero were guarded.
- **What we did:** every secant and both end slopes go through `safedivide` (zero for a zero-width interval), which returns the plain quotient for distinct abscissae. `test_numerics` pins a table with a repeated abscissa: the spline at 2.5 was NaN and is 5.7708. Byte-identical: the quick set and the 29 rotating decks.

### M11. ysplin was wrong for 3 points and read unset values for 2; OPAL92's dln κ/dln T was wrong at its top temperature node

- **Status:** fixed in [4eaaaece](https://github.com/parallelpro/yrec/commit/4eaaaece) (2026-10-29).
- **In the original F77?** Yes (YCK 3/91). The first and last interior rows are both row 2 when N = 3, so the last overwrites the first; for N = 2 row 2 reads `H(3)` and `F(3)`, never set ([ysplin.f:26-47](https://github.com/parallelpro/yrec/blob/6cd56739/src/ysplin.f#L26-L47)):
  ```fortran
        I=2
        C(2,I)=CONST1*H(I)*F(I+1)+CONST2*H(I+1)*F(I)
  ...
        I=N-1
        C(2,I)=CONST1*H(I)*F(I+1)+CONST2*H(I+1)*F(I)
  ```
- **Now:** `ysplin` and `intpol` are wrappers of `natural_spline`, `splint` (with an optional derivative) and `compute_spline_coefficients`, [cubic_spline.f90:146-210](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/cubic_spline.f90#L146-L210); OPAL92's interpolation in T, [opal92_opacity.f90:509-520](https://github.com/parallelpro/yrec/blob/396bbb7c/src/kap/opal92/opal92_opacity.f90#L509-L520).
- **What was wrong:** on (0,0), (1,1), (3,0) the middle slope was 0.9545 instead of the natural spline's 1/2. OPAL92 reaches n = 3 at the top grid temperature (T6 = 100, log T = 8), where its six temperature rows clip to three, so intpol's derivative was wrong there (its value was right).
- **Expected behaviour, and why:** the natural cubic spline through the points, for any n ≥ 2.
- **What we did:** OPAL92 and Kurucz 1990 build their coefficient tables from `natural_spline`'s second derivatives and interpolate in T with `natural_spline` and `splint`'s derivative; the slope-form solve is gone.
  - At X = 0.7, Z = 0.019, log T = 8, log ρ = 3 (a new `test_kap` point), dln κ/dln T goes from −0.32606 to −0.35187; κ and dln κ/dln ρ are unchanged. Elsewhere `test_kap`'s OPAL92 and Kurucz 1990 lines move by at most 6.3e-16 relative.
  - No deck uses the OPAL92 or Kurucz 1990 opacities (no inlist sets `lopal92` or `lkur90`): the full set is byte-identical.

### M12. ludcmp passed a singular matrix as solved, and a NaN column wrote a(0,k)

- **Status:** fixed in [37296d25](https://github.com/parallelpro/yrec/commit/37296d25) (2026-10-29).
- **In the original F77?** Yes. LUDCMP is Numerical Recipes': a zero pivot becomes `TINY`, and the pivot row is set only when a candidate's scaled size is `.GE.` the running maximum, which a NaN never is ([ludcmp.f:9](https://github.com/parallelpro/yrec/blob/6cd56739/src/ludcmp.f#L9), [:63-67](https://github.com/parallelpro/yrec/blob/6cd56739/src/ludcmp.f#L63-L67), [:82](https://github.com/parallelpro/yrec/blob/6cd56739/src/ludcmp.f#L82)):
  ```fortran
        PARAMETER (NMAX=100,TINY=1.0D-20)
                IF (DUM.GE.AAMAX) THEN
                    IMAX = I
                IF (A(J,J).EQ.0D0) A(J,J) = TINY
  ```
- **Now:** `ludcmp`, [linear_solvers.f90:39-127](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/linear_solvers.f90#L39-L127).
- **What was wrong:** YREC's `ierr` reported only a zero row, so a rank-deficient matrix passed with `ierr = 0`: rows (1,2,3), (4,5,6), (7,8,9) give a pivot of 1.9e-17 of its row's scale, and row 3 = row 1 + row 2 an exact zero. A column whose candidate pivots were all NaN left the pivot row unset (0 under `-finit-local-zero`), and the row swap wrote `a(0,k)`.
- **What we did:** `numerics_singular` when |pivot| / max|row| ≤ `singular_pivot_tolerance` = 1e-14; the factorization still completes as Numerical Recipes' does, so a caller that ignores `ierr` behaves as before. The pivot row starts at j, and a NaN pivot returns `numerics_bad_input`. On the decks that run the Thoul solver, the only caller, the smallest scaled pivot is 0.54. Byte-identical: the quick set and the diffusion decks.

### M13. solve_tridiagonal never checked its first pivot, nor did rotation's copy

- **Status:** fixed in [bcdcd649](https://github.com/parallelpro/yrec/commit/bcdcd649) (2026-10-29); rotation's angular momentum solve has the check since it calls the same two halves ([9c1c03b1](https://github.com/parallelpro/yrec/commit/9c1c03b1)).
- **In the original F77?** Yes, in TRIDIAG_GS and in TRIDIA, which test every later pivot ([tridiag_gs.f:10-15](https://github.com/parallelpro/yrec/blob/6cd56739/src/tridiag_gs.f#L10-L15), [tridia.f:14-19](https://github.com/parallelpro/yrec/blob/6cd56739/src/tridia.f#L14-L19)):
  ```fortran
        BET = B(1)
        EX(1) = EX_PRIME(1)/BET
        DO 1 J=2,NPT
           GAMA(J) = C(J-1)/BET
           BET = B(J) - A(J)*GAMA(J)
           IF( BET.EQ.0 ) STOP '#TRIDIA:SINGULAR MATRIX'
  ```
- **Now:** `tridiagonal_forward_sweep` and `tridiagonal_back_substitution`, [linear_solvers.f90:244-309](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/linear_solvers.f90#L244-L309); `solve_angular_momentum_transport`, [am_transport_coefficients.f90:180-226](https://github.com/parallelpro/yrec/blob/396bbb7c/src/rotation/seculr/am_transport_coefficients.f90#L180-L226).
- **What was wrong:** b(1) = 0 divided by zero and returned NaN with `ierr = 0`, against the header ("a zero pivot returns ierr = 1").
- **What we did:** the first pivot is checked before the division. Rotation's solve now calls numerics' two halves around its wind step, in the same order of operations, so its results are the same bit for bit. Byte-identical: the quick set, the rotating set and the unit tests.

### M14. locate's 1% slack was relative to the end value

- **Status:** fixed in [a85bfe5e](https://github.com/parallelpro/yrec/commit/a85bfe5e) (2026-10-29).
- **In the original F77?** Yes; the slack is YREC's (Numerical Recipes' locate has none) ([locate.f:25-30](https://github.com/parallelpro/yrec/blob/6cd56739/src/locate.f#L25-L30)):
  ```fortran
        IF((J .EQ. 0) .AND. (X .GT. 0.99D0*XX(1)))THEN
           J = 1
        END IF
        IF((J .EQ. N) .AND. (X .LT. 1.01D0*XX(N)))THEN
           J = N-1
  ```
- **Now:** `locate`, [table_interpolation.f90:364-392](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/table_interpolation.f90#L364-L392).
- **What was wrong:** the factors 0.99 and 1.01 are a slack only for an increasing table of positive values. For a decreasing table every point outside it was clamped to the end interval, and for xx(1) ≤ 0 a point exactly on the first node returned 0 (out of range). The old test pinned both.
- **What we did:** the slack is 1% of the table's span |xx(n) − xx(1)| beyond either end, in either direction, and a point on the first node gives 1. The callers treat 0 like 1 and n like n−1 (LAOL89's X index, the Allard rows and columns), so their results do not move: the quick set and the Allard deck are byte-identical.

### M15. lir failed silently

- **Status:** fixed in [df8a8487](https://github.com/parallelpro/yrec/commit/df8a8487) (2026-10-29).
- **In the original F77?** Yes. LIR's two diagnostics have their formats but no write: fewer than two points returns before `INTER` is set, and equal end abscissae return with `Y` unset and `INTER = 1` ([lir.f:39-42](https://github.com/parallelpro/yrec/blob/6cd56739/src/lir.f#L39-L42), [:53](https://github.com/parallelpro/yrec/blob/6cd56739/src/lir.f#L53), [:124-132](https://github.com/parallelpro/yrec/blob/6cd56739/src/lir.f#L124-L132)):
  ```fortran
        IF(NT.LT.2) GO TO 101
        INTER=1
        IF(DIFF) 4,102,3
    101 CONTINUE
        RETURN
    102 CONTINUE
        RETURN
  ```
- **Now:** `lir`, [table_interpolation.f90:221-315](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/table_interpolation.f90#L221-L315).
- **What we did:** both cases set the interpolation flag to 0 and return `numerics_bad_input` or `numerics_coincident_points` through a new optional `ierr` (the MHD EOS's three calls compile unchanged). The old test pinned the caller's −7 in the flag. Byte-identical: the quick set.

### M16. polint's coincidence test was absolute

- **Status:** fixed in [1c67bae5](https://github.com/parallelpro/yrec/commit/1c67bae5) (2026-10-29).
- **In the original F77?** Yes, with a `STOP` ([polint.f:31-32](https://github.com/parallelpro/yrec/blob/6cd56739/src/polint.f#L31-L32)):
  ```fortran
            DEN = HO - HP
            IF(dabs(DEN).LT.1.0D-20)STOP
  ```
- **Now:** `polint`, [table_interpolation.f90:150-200](https://github.com/parallelpro/yrec/blob/396bbb7c/src/numerics/table_interpolation.f90#L150-L200).
- **What was wrong:** |x_i − x_{i+m}| < 1e-20 is an absolute test on a dimensional difference: abscissae 1e-25 apart were "coincident", nearly equal large ones were not. Rotation's shape quadrature extrapolates over steps 0.25^(j−1), which fall below 1e-20 at j = 34; a deck that refined that far would have stopped.
- **What we did:** |den| ≤ 4 ε max(|x_i − x|, |x_{i+m} − x|) (`coincidence_tolerance`), which an exact coincidence always meets. `test_numerics` checks polint on abscissae scaled by 1e-25. Byte-identical: the quick set, the Allard deck and `test_rotation`.

### M17. crmath's pow and exp10 at their edges

- **Status:** fixed in [c2638e4b](https://github.com/parallelpro/yrec/commit/c2638e4b) (2026-10-29), where cheap; one difference between the backends is kept by decision.
- **In the original F77?** No; `math_lib` and its crmath backend came with [71d1447b](https://github.com/parallelpro/yrec/commit/71d1447b) (2026-08-30).
- **Now:** `pow_r` and `exp10`, [math_lib.f90:87-162](https://github.com/parallelpro/yrec/blob/396bbb7c/src/math/math_lib.f90#L87-L162).
- **What was wrong:** under crmath, `pow(x, y)` and `exp10(x)` took `floor(y)` before testing |y| < 100, which overflows the default integer for |y| ≥ 2³¹ (undefined behaviour), and `pow` of a negative base with an integer-valued |y| ≥ 100 gave NaN where x**y is finite. The headers called exp10 exact for every integer.
- **What we did:** |y| is tested first, and a negative base with an integer-valued y takes exp(log(−x) y) with the sign of (−1)^y. exp10 is documented as exact for 0 ≤ x ≤ 22 and correctly rounded for −22 ≤ x < 0. Kept, and documented: crmath's pow(0, y) = 0 for every y (MESA's definition), where the intrinsic backend gives 1 for y = 0 and +Inf for y < 0. The new `test_math` checks each backend's contract (`test_math_crmath` is built whenever the MESA SDK is there). The intrinsic backend, which the decks use, is untouched.

### M18. The Monte Carlo model file kept only the last run's model

- **Status:** fixed in [49e417d7](https://github.com/parallelpro/yrec/commit/49e417d7) (2026-10-29).
- **In the original F77?** No. WRTLST rewound only the final-model unit, so WRTMONTE appended one model per Monte Carlo run ([wrtlst.f:141](https://github.com/parallelpro/yrec/blob/6cd56739/src/wrtlst.f#L141), [wrtmonte.f:25-27](https://github.com/parallelpro/yrec/blob/6cd56739/src/wrtmonte.f#L25-L27)):
  ```fortran
        IF(IWRITE.EQ.ILAST) REWIND ILAST
  ```
  The unconditional rewind came with the .mod format, [0066e7af](https://github.com/parallelpro/yrec/commit/0066e7af) (2026-08-26).
- **Now:** `write_mod_model(s, unit, replace)`, [model_files.f90:56-98](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/model_files.f90#L56-L98).
- **What was wrong:** `write_mod_model` rewound its unit before every model, and `write_monte_carlo_model` writes each run's model to `monte_carlo_file2_path`, so the file held one model, the last run's.
- **What we did:** `replace` is true for the final-model file (rewound: one model, the latest) and false for the Monte Carlo file (appended). On a two-run Monte Carlo deck with the solar calibration on (`cm_monte_carlo` with `calibrate_solar_model`, run 2 with S11 × 1.01) the .mc2 file now holds run 1's model followed by run 2's, which is byte-identical to the old file; the .mc1 summary, .mod, .log and history.data are byte-identical. The shipped `cm_monte_carlo` does not calibrate and writes no Monte Carlo model: byte-identical, as is the full set.

### M19. The run log labelled the central deuterium stop "central: log rho"

- **Status:** fixed in [34fc2988](https://github.com/parallelpro/yrec/commit/34fc2988) (2026-10-29). Changes pinned logs, reseeded at [2f84ae0e](https://github.com/parallelpro/yrec/commit/2f84ae0e).
- **In the original F77?** No; PARMIN printed "CENTRAL D" ([parmin.f:1321-1323](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L1321-L1323)). It was relabelled by the MESA-style run log, [b38ec9c3](https://github.com/parallelpro/yrec/commit/b38ec9c3) (2026-08-28).
  ```fortran
    370          FORMAT(1X,'EVOLVE TO AGE ',L1,' SET DELT ',
       *              L1,' FINAL AGE ', E9.2,' FIXED TSTEP ',E9.2,
       *  ' CENTRAL D ',E10.4,' CENTRAL X ',E12.4,' CENTRAL Y ',E12.4)
  ```
- **Now:** format 370, [run_cards.f90:127](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/run_cards.f90#L127).
- **What was wrong:** the run plan's stop line printed `central_deuterium_stop` (the central D mass fraction at which the card stops) as "central: log rho =". The registry called it MESA's `log_center_density_limit` too.
- **What we did:** "central: D =", next to "X =" and "Y ="; the registry row says what the control is (it has no MESA counterpart). On the full set every output file is byte-identical except the .log of the 55 decks whose cards have a stop, in that one line; 34 of them are pinned logs (the configuration matrix, the solar decks, the dbl and ZAMS decks, the testsuite), reseeded at 2f84ae0e.

### M20. A rescale-and-evolve card ignored its age, timestep and central stops

- **Status:** fixed in [8c80883a](https://github.com/parallelpro/yrec/commit/8c80883a) (2026-10-29).
- **In the original F77?** Yes. PARMIN sets `LENDAG` and `LSETDT` in the `KINDRN = 1` branch only; its `KINDRN = 3` branch, when it starts from the previous card's model, writes the rescale card's format 451, its own 551 being defined but unused ([parmin.f:1288-1298](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L1288-L1298), [:1355-1373](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L1355-L1373)):
  ```fortran
         IF(KINDRN(NKIND).EQ.1) THEN
  ...
                 LENDAG(NKIND)=.TRUE.
  ...
           ELSE IF(KINDRN(NKIND).EQ.3) THEN
  ...
                 WRITE(IOWR,451) NKIND
  ```
- **Now:** `set_stop_flags`, [run_cards.f90:68-73](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/run_cards.f90#L68-L73), called for both evolving kinds.
- **What was wrong:** a kind-3 card's `max_age`, `timestep_override` and `central_*_stop` were ignored (the flags kept their default, false), and its run-plan line said "rescale the previous card's model".
- **What we did:** the flags are derived for kinds 1 and 3; a kind-3 card echoes its stop line, and its run-plan line reads "rescale and evolve the previous card's model". No shipped deck sets a stop on a kind-3 card (314 decks scanned), so no output moves; the run-plan line changes in seven unpinned seed decks.

### M21. Most path controls never had their {YREC_...} placeholders expanded

- **Status:** fixed in [b1e7ae2c](https://github.com/parallelpro/yrec/commit/b1e7ae2c) and [30fae4e5](https://github.com/parallelpro/yrec/commit/30fae4e5) (2026-10-29).
- **In the original F77?** No; the placeholders came with the YREC-team commit [e5c35e5d](https://github.com/parallelpro/yrec/commit/e5c35e5d) (2026-02-11, #30), which expanded a fixed list of path controls.
- **Now:** `expand_paths`, [run_files.f90:69-114](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/run_files.f90#L69-L114).
- **What was wrong:** 15 of the path controls (and `inlist_used_file`) were expanded. A {YREC_INPUT}, {YREC_START} or {YREC_OUTPUT} in any other reached the open literally: the OPAL 1995 and 2001 EOS tables, the Kurucz, LAOL and OPAL 1992 opacity tables, the eight MHD tables, the seven Alexander 1995 tables, the isochrone, Monte Carlo and dynamics files, and the two columns files.
- **What we did:** every string control that names a file is expanded, from one list; `test_generated_files.py` checks the list against the registry (every string control but names, codes and labels) and fails when a line is removed. The columns files are expanded since the output reads them after the controls read (M25). No shipped deck uses a placeholder in the newly expanded controls: byte-identical.

### M22. pulse_format: a case-sensitive startup check, and an unknown name meant GYRE

- **Status:** fixed in [ba392477](https://github.com/parallelpro/yrec/commit/ba392477) (2026-10-29).
- **In the original F77?** No; the pulse output and its format control are 2026 ([eb624a69](https://github.com/parallelpro/yrec/commit/eb624a69), [ec77ad33](https://github.com/parallelpro/yrec/commit/ec77ad33), 2026-08-22/23).
- **Now:** `pulse_format_code`, [yrec_output.f90:112-128](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/yrec_output.f90#L112-L128), used by `output_open` and the writer.
- **What was wrong:** the startup check rejected only an exact 'GSM' on a build without HDF5, while the writer took 'gsm' as GSM too; so 'gsm' on such a build passed startup and then warned at every profile, writing nothing. Any other name, a typo or 'Fgong', silently meant GYRE.
- **What we did:** one function reads the name in any case for both; with `write_pulse_flag` on, an unknown name stops the run at startup with exit 1 ("pulse_format = 'GYR' is not GYRE, FGONG or GSM"). `test_error_exit` checks the unknown name. Byte-identical: the quick set.

### M23. The YREC7 and MODEL2 readers trusted the file

- **Status:** fixed in [e2335b37](https://github.com/parallelpro/yrec/commit/e2335b37) (2026-10-29).
- **In the original F77?** Yes. GETMODEL2 reads the zone count and indexes the mass array with it in the same statement, and neither GETMODEL2 nor GETY7 has an `IOSTAT` or a check against `JSON` ([getmodel2.f:37-39](https://github.com/parallelpro/yrec/blob/6cd56739/src/getmodel2.f#L37-L39)):
  ```fortran
        READ(IREAD,10) ATEMP,MODEL,M,SMASS,TEFFL,BL,HSTOT,DAGE,
       *   DDAGE,HS(1),HS(M)
  ```
- **Now:** `read_yrec7` and `read_model2`, [model_files.f90:206-567](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/model_files.f90#L206-L567).
- **What was wrong:** a header with more than `json` zones wrote past the star's arrays; a malformed or truncated model ended in a Fortran runtime error; `read_model2` had no `ierr`, and `read_mod_model`'s tag and named-global reads had no `iostat`.
- **What we did:** the zone count must be 1 to `json` (read_model2 reads its header into locals first), YREC7's convection-zone edges, which index its composition records, must be 1 to nz, and every read has an `iostat`: a failure reports "cannot read <part> of the starting model <file>: <reason>" and returns `ierr = 1`, which `read_starting_model` passes on. `test_error_exit` covers a MODEL2 model with 99999 zones and one cut inside its shells. Byte-identical: the quick set, whose decks read MODEL2 models, and a run from a YREC7 seed model.

### M24. The output directory was made by an unquoted shell command

- **Status:** fixed in [466b46ff](https://github.com/parallelpro/yrec/commit/466b46ff) (2026-10-29).
- **In the original F77?** No; the YREC-team commits [4d315ebd](https://github.com/parallelpro/yrec/commit/4d315ebd) (2026-02-20, #31) and [e66afb17](https://github.com/parallelpro/yrec/commit/e66afb17) (2026-02-24, #40) added `CALL system('mkdir -p ...')` with the path cut into a 256-character buffer.
- **Now:** `directory_of` and `make_directory`, [run_files.f90:201-237](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/run_files.f90#L201-L237).
- **What was wrong:** with no '/' in `log_output_file` mkdir ran with no operand (a usage error on the terminal); the slice arithmetic overran the buffer for a directory longer than about 247 characters; and the path went to the shell unquoted, so a directory with a space was never created and the run stopped at its first open. The directory was also computed a second way by the output.
- **What we did:** the directory is found once (up to the log path's last '/') and created with its parents through the C library's `mkdir` (iso_c_binding; mode 777 less the umask, as `mkdir -p`). `test_mesa_output` runs into a nested directory with a space in its name (the previous binary exits 1). Byte-identical: the quick set.

### M25. The output's state lived in module variables, and the history file stayed open across runs

- **Status:** fixed in [226fabae](https://github.com/parallelpro/yrec/commit/226fabae) (2026-10-29).
- **In the original F77?** No; the MESA-format writers are 2026 ([eb624a69](https://github.com/parallelpro/yrec/commit/eb624a69), split in [dad3c4de](https://github.com/parallelpro/yrec/commit/dad3c4de), 2026-08-29).
- **Now:** `output_state`, owned by `run_yrec`, and `output_open`, [yrec_output.f90:49-100](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/yrec_output.f90#L49-L100).
- **What was wrong:** the output directory, the profile counter, the history path and selection, the profile selection and the progress-line counters were module variables. The history file's unit (a `newunit`) was found again by `inquire(file=)` at every row, so `yrec_reset`'s sweep of units 7-99 never closed it, and a second in-process run into another directory left the first file open.
- **What we did:** one record, `output_state`, that `run_yrec` owns and passes down: `output_open` after the controls read (it first closes a history file a previous call left open), `output_write_model` each model, `output_close` at the end. The writers take the star `intent(in)`. Byte-identical: the quick set, `test_reentry` and `test_mesa_output`.

### M26. A profile column without a value in a region was a silent zero

- **Status:** fixed in [770f6ce2](https://github.com/parallelpro/yrec/commit/770f6ce2) (2026-10-29).
- **In the original F77?** No; the stitched model's profile values are 2026 ([065c064f](https://github.com/parallelpro/yrec/commit/065c064f), 2026-08-25).
- **Now:** `interior_value`, `envelope_value` and `atmosphere_value`, [stitched_model.f90:478-557](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/diagnostics/stitched_model.f90#L478-L557) and after; the error at [:896](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/diagnostics/stitched_model.f90#L896).
- **What was wrong:** the profile values were a numbered `select` whose envelope and atmosphere cases fell back to zero for any column they did not list, so a new column forgotten in a region printed zeros there.
- **What we did:** each region's value routine lists every column (those the envelope and atmosphere do not compute are an explicit zero group); a column without a case stops the run at the first model ("the profile column beta has no value in the envelope"), and `test_default_columns_lists_in_sync` checks the case lists statically. Byte-identical: the quick set and `test_mesa_baselines`.

### M27. The zoning tolerances were zero unless use_new_rate_variables was set

- **Status:** fixed in [a0b4f434](https://github.com/parallelpro/yrec/commit/a0b4f434) (2026-10-29).
- **In the original F77?** No. `HPTTOL` was the working vector with its DATA defaults, overwritten from the named tolerances only under `LNEWVARS` ([parmin.f:406-407](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L406-L407), [remap.f:194-198](https://github.com/parallelpro/yrec/blob/6cd56739/src/remap.f#L194-L198)):
  ```fortran
        DATA HPTTOL/1.0D-8,8.0D-2,5.0D-2,5.0D-2,1.0D0,1.0D0,0.0D0,5.0D-2,
       *            2.0D-2,5.0D-2,5.0D-2,1.0D-1/
        IF(LNEWVARS)THEN
  C MINIMUM LOG MASS SPACING BETWEEN SHELLS
           HPTTOL(1) = TOL_DM_MIN
  ```
  The copy was dropped by [a576ec06](https://github.com/parallelpro/yrec/commit/a576ec06) (2026-08-20, sweep batch 4, "lnewvars is true in every current example"), the pattern of N4.
- **Now:** `chi_grid_scale = hpttol` in the generated `adopt_canonical.inc`, included by `adopt_canonical_names`, [controls_reader.f90:243-259](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/controls_reader.f90#L243-L259).
- **What was wrong:** `chi_grid_scale(12)` (the rezoning tolerances, slots `ichi_*`) had no initializer and was filled only by `map_user_inputs` under `lnewvars`, whose default is false. A deck without the switch zoned with zero tolerances, and the core extension divided by them. The registry called `hpttol` unused.
- **Expected behaviour, and why:** without the switch the zoning uses `hpttol` and its defaults, as the F77 did.
- **What we did:** the copy is back; `map_user_inputs` still overwrites it under `lnewvars`. `hpttol` is a live control with the inlist name `zoning_tolerances`, documented slot by slot, and the `use_new_rate_variables` row says what the switch does (the zoning and timestep tolerances from the named controls; nothing to do with reaction rates). `test_controls.py` reads `cm_base` with the switch off and finds `chi_grid_scale` and `atime` at their documented defaults. All 539 shipped deck files set the switch: byte-identical.

### M28. The settling switch, a control, was overwritten as a run-long latch

- **Status:** fixed in [af946c0c](https://github.com/parallelpro/yrec/commit/af946c0c) (2026-10-29).
- **In the original F77?** Yes. MIX and ROTMIX switch `LDIFY` off once the central X falls below `XMIN` ([mix.f:377-381](https://github.com/parallelpro/yrec/blob/6cd56739/src/mix.f#L377-L381), [rotmix.f:98-101](https://github.com/parallelpro/yrec/blob/6cd56739/src/rotmix.f#L98-L101)), and only PARMIN, called once per run, sets it:
  ```fortran
        IF(LDIFY)THEN
           IF(HCOMP(1,1).LT.XMIN)THEN
              LDIFY=.FALSE.
  ```
- **Now:** `star%settling_active`, armed in `begin_kind_card` ([run_driver.f90:331](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/run/run_driver.f90#L331)) and switched off in [composition_step.f90:256-264](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/step/composition_step.f90#L256-L264) and [mid_timestep_model.f90:460-467](https://github.com/parallelpro/yrec/blob/396bbb7c/src/rotation/mid_timestep_model.f90#L460-L467).
- **What was wrong:** `diffuse_helium_active` is a user control, but the two settling routines overwrote it as a "stop for good" latch. A later run card, or a later Monte Carlo run, that starts from a hydrogen-richer model ran without settling.
- **Expected behaviour, and why:** each card starts with the deck's switch; the latch describes the current model sequence only.
- **What we did:** the control never changes after the read; the latch is the star's `settling_active`, armed from the control at the start of every card (so every Monte Carlo run too) and switched off by the settling routines as before. rotation/ reads no control any more. No shipped deck moves (the full set is byte-identical). On a four-card variant of `cm_settle_he` built to reach it (rescale and evolve, evolve to 1 Gyr, rescale X to 0.74 so that X_c goes from 0.647 to 0.659, evolve to 1.5 Gyr; `hydrogen_diffusion_floor` 0.655), card 4 now settles until X_c falls below the floor two models later, where before it did not settle at all; at 1.5 Gyr Δlog L = +5.4e-5, Δlog Teff = −4.0e-5, ΔX_c = −1.7e-4, ΔY_surf = −3.6e-4.

### M29. The semiconvection-with-overshoot error reached the run log only, with its flags mislabelled

- **Status:** fixed in [db99ec9a](https://github.com/parallelpro/yrec/commit/db99ec9a) (2026-10-29). Failure path only.
- **In the original F77?** Yes, both ([parmin.f:896-902](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L896-L902)):
  ```fortran
              WRITE(ISHORT,2)LSEMIC,LOVSTE,LOVSTC,LOVSTM
      2       FORMAT(1X,'ERROR IN SUBROUTINE PARMIN'/'SEMI-CONVECTION',
       * ' AND OVERSHOOT FLAGS BOTH TURNED ON'/'FLAGS LSEMIC',L2,
       * ' OVERSHOOT - CORE,ENVELOPE,INTERMEDIATE-',3L2/'RUN STOPPED')
  ```
- **Now:** `check_controls`, [controls_check.f90:75-118](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/controls_check.f90#L75-L118).
- **What was wrong:** the error went to the run log only, while the two errors beside it went to the terminal too, and it printed the envelope flag under the label "core".
- **What we did:** `check_controls` reports it on the terminal and in the run log through `report_control`, with the flags in the labelled order and "SEMI-CONVECTION" kept in the text (`test_error_exit` looks for it). No pinned log changes.

### M30. The light-element-sum warning never printed its end

- **Status:** fixed in [db99ec9a](https://github.com/parallelpro/yrec/commit/db99ec9a) (2026-10-29). Failure path only.
- **In the original F77?** Yes; format 595 has eleven descriptors for seven values, copied from format 596 above it ([parmin.f:1105-1110](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L1105-L1110)):
  ```fortran
                 WRITE(*,595)XH2_INI,XHE3_INI,XLI6_INI,XLI7_INI,XBE9_INI,
       * XB10_INI,XB11_INI
   595           FORMAT('SUM OF LIGHT ELEMENT MASS FRACTIONS EXCEEDS 1',
       * 11E12.4,' MIX NOT MODIFIED')
  ```
- **Now:** [derived_controls.f90:218](https://github.com/parallelpro/yrec/blob/396bbb7c/src/io/derived_controls.f90#L218).
- **What was wrong:** the output stopped at the eighth descriptor, so "); mixture not modified" never printed.
- **What we did:** `7es12.4`.

### M35. The MHD EOS read no tables

- **Status:** fixed in [43e9b717](https://github.com/parallelpro/yrec/commit/43e9b717) (2026-10-29).
- **In the original F77?** No. MHDTBL opens all eight tables on the controls' unit numbers, and MHDST reads a table only when its unit is positive, the F77's way of saying "no table" ([mhdtbl.f:29-32](https://github.com/parallelpro/yrec/blob/6cd56739/src/mhdtbl.f#L29-L32), [mhdst.f:59-62](https://github.com/parallelpro/yrec/blob/6cd56739/src/mhdst.f#L59-L62)):
  ```fortran
  C     USES 3 ZAMS-TYPE TABLES AND OPTIONALLY 5 CENTER-TYPE TABLES.
        OPEN(UNIT=IOMHD1, FILE=FMHD1, STATUS='OLD',
  C     READ ZAMS TABLES
        IF(IR1.GT.0) THEN
  ```
  It was introduced by [c00a9a20](https://github.com/parallelpro/yrec/commit/c00a9a20) (2026-10-11), which opened the tables with `open(newunit=)`, whose units are negative.
- **Now:** `read_mhd_tables` and `compute_mhd_eos`, [mhd_eos.f90](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/mhd/mhd_eos.f90).
- **What was wrong:** with `use_mhd_eos` nothing was read and `ierr` stayed 0: the temperature limits indexed row 0, and `compute_mhd_eos` interpolated zero tables (1/χ_ρ = 1/0). The opens had no `iostat`, so a missing file ended in a runtime error, although `eos_init` promises `ierr = 1`.
- **What we did:** all eight tables are opened with `iostat` (a missing one returns `ierr = 1` with a run-log message) and read; every read has an `iostat`, and a family with fewer than 4 rows or densities is an error. `test_eos` gains an MHD reader test on eight synthetic tables of an ideal gas, where the interpolation is exact (the density in each family and at the boundary, the gradients, both out-of-range errors, a missing and a truncated file), each check confirmed by a mutation. The MHD tables are not shipped and no deck uses them: the quick set and the unit tests are byte-identical.

### M36. The MHD EOS returned the previous point's result as success outside its tables

- **Status:** fixed in [43e9b717](https://github.com/parallelpro/yrec/commit/43e9b717) (2026-10-29).
- **In the original F77?** Yes. MHDPX1 leaves for its return with `IERR` untouched when log T is outside the tables ([mhdpx1.f:35-38](https://github.com/parallelpro/yrec/blob/6cd56739/src/mhdpx1.f#L35-L38)); INTPT returns without output above a row's highest pressure; and at log T equal to the first centre row MHDPX1 fills the centre rows but interpolates in X across the ZAMS rows, which are SAVEd from an earlier call:
  ```fortran
        IF(TL.LT.TMINI .OR. TL.GT.TMAXI) THEN
  C         IRANGE = 0
            GO TO 999
        END IF
  ```
- **Now:** `mhdpx1`, [mhd_eos.f90:561](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/mhd/mhd_eos.f90#L561), and `compute_mhd_eos`.
- **What was wrong:** a log T outside the tables returned `ierr = 0` with the previous point's result; a log P above a row left that row's values at zero (the mean-weight check then usually, but not always, raised `ierr`).
- **Expected behaviour, and why:** a point the tables do not cover is an error: the MHD path has no other EOS to fall back on.
- **What we did:** both out-of-range cases return `ierr = 1`, and `compute_mhd_eos` writes the point to the run log; the X interpolation follows the family the temperature selected, and the rows start at zero. No deck uses MHD.

### M37. The Debye-Hückel metals beyond C, N and O came from the opacity tables

- **Status:** fixed in [15af2ef3](https://github.com/parallelpro/yrec/commit/15af2ef3) (2026-10-29). Results change in every deck.
- **In the original F77?** Yes. The correction's metal mixture `ZDH` is read from the LAOL89 opacity table's header, and EQRELV sums it into `ALMIX` ([rdlaol.f:52](https://github.com/parallelpro/yrec/blob/6cd56739/src/rdlaol.f#L52), [eqrelv.f:31-32](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqrelv.f#L31-L32), [:66-72](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqrelv.f#L66-L72)):
  ```fortran
        READ(IOLAOL,130)ZHIT,ZDH                                ! rdlaol.f
  C     ZDH IS ARRAY OF RELATIVE MASS FRACTIONS OF LAOL METAL MIXTURE
              ALMIX = ALMIX + DHNUE(II)*ZDH(II)                 ! eqrelv.f
           CMFDH0 = 2.0D0*XXDH+1.5D0*YYDH+0.5D0*(ZZDH+ALMIX)
  ```
- **Now:** `eos_set_mixture`, [eos_lib.f90:213](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/eos_lib.f90#L213), which fills the correction's named slots `idh_*`.
- **What was wrong:** without the LAOL89 opacities, as in every shipped deck, the slots of the metals heavier than O were zero, so the charge-weighted sum `almix` counted only C, N and O, about 60% of it.
- **Expected behaviour, and why:** the correction is that of the gas being modelled; it should not depend on which opacity tables are read.
- **What we did:** `eos_set_mixture` fills Ne, Na, Mg, Al, Si, Ar and Fe from eos's envelope mixture, as mass fractions X_i = fx_i m_u A_i; P, S, Cl, Ca, Ti, Cr, Mn and Ni, which YREC's 12-species mixture lacks, stay zero (about 6% of `almix`). kap no longer hands eos the LAOL mixture, and `eos_set_debye_huckel_z` is gone. `test_eos` checks that the correction raises ρ, and more when the heavier metals are all Fe.
  - All 56 decks change (every one turns the correction on). Solar decks with OPAL 2006 and SCV: final |Δlog L| ≤ 2.7e-5, |ΔX_c| ≤ 3.3e-5; age-matched up to 2.9e-4 in log L before the main sequence.
  - The Yale-EOS solar deck `Test_solarGS98_yaleeos`: final Δlog L +1.09e-3, Δlog Teff +1.19e-4, ΔX_c −7.9e-4. Where OPAL covers a point it replaces the Yale result, so only the Yale-only decks feel the correction throughout the star; there `almix` grows by about two thirds and the density at fixed P by about 2%.
  - Measured on the integration branch, for decks the package did not run: the 0.03 Msun dbl-to-ZAMS deck, a marginal one, stops at a different age (1.2544e7 → 1.2255e7 yr; final Δlog L +6.9e-3); the 1 Msun ZAHB deck, final Δlog L −5.1e-4; `o24_rescale_rotating` ends on its model budget at a different age (1.645e8 → 1.345e8 yr).

### M38. With diffusion of Z, the EOS's hydrogen and helium terms used other atomic weights

- **Status:** fixed in [3db309e9](https://github.com/parallelpro/yrec/commit/3db309e9) (2026-10-29). Results change in the decks with diffusion of Z.
- **In the original F77?** Yes. EQSTAT2's diffusion-of-Z branch divides the zone's H and He excesses by `ATOMW2(7)` and `ATOMW2(12)` (1.008, 4.004), while its default branch multiplies by `ATOMWT` (0.9921, 0.24975) ([eqstat2.f:48](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L48), [:62-80](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L62-L80), [:106-108](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L106-L108)):
  ```fortran
        DATA ATOMWT/0.9921D0,0.24975D0,0.08322D0,0.4995D0/
           DXX = (X - XENV)/ATOMW2(7)
           DFX1 = DFX1*ATOMWT(1)
  ```
- **Now:** `zone_mixture`, [blended_eos.f90:262](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/blended_eos.f90#L262).
- **What was wrong:** 1/1.008 = 0.992063 against 0.9921, so a zone with X ≠ X_env got 4e-5 of its hydrogen term's excess wrong, even for a uniform Z.
- **Expected behaviour, and why:** one composition gives one 1/μ_ion, whichever branch computes it.
- **What we did:** the diffusion branch uses the default branch's constants; with Z = Z_env the two give the same 1/μ_ion bit for bit, and `test_eos` checks that. 31 of 50 decks change, those with diffusion of Z: final |Δlog L| ≤ 3.0e-5, |ΔX_c| ≤ 4.9e-5; age-matched Δlog L up to 3.5e-4 (rv_sbm2) at 35-45 Myr, the approach to the main sequence, where a small shift in timing shows at a matched age.

### M39. The Saha ramp blended two second derivatives as a sum, and not the ionization fractions

- **Status:** fixed in [f8114dca](https://github.com/parallelpro/yrec/commit/f8114dca) (2026-10-29). Found by the new continuity checks of `test_eos`. Results change in most decks.
- **In the original F77?** Yes. EQSTAT2 saves the Saha values, differences most of them against the fully ionized ones and blends r = r_FI + w (r_Saha − r_FI); `QDTT0` and `QDTP0` are never differenced, and `FXION` is overwritten with the fully ionized values before the blend ([eqstat2.f:216-222](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L216-L222), [:228-232](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L228-L232), [:253-262](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L253-L262)):
  ```fortran
        FXION(1) = 1.0D0
  ...
        FXION(1) = 1.0D0 + WT0*(FXION(1) - 1.0D0)
  ...
           QDTT = QDTT + WT0*QDTT0
  ```
- **Now:** `blend_saha_into`, [blended_eos.f90:360](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/blended_eos.f90#L360).
- **What was wrong:** within 0.5 dex below the Saha cutoff, `dlnrho_dlnt_dt` and `dlnrho_dlnt_dp` were r_FI + w r_Saha, the sum of the two solutions' at the ramp's foot; the ionization fractions were blended with themselves, so the ramp was fully ionized and the fractions jumped at its foot (the conductive opacity reads them).
- **What we did:** both are blended like every other quantity; the continuity checks at the ramp's foot and at the cutoff compare them, and fail with the old blend. 49 of 50 decks change: final |Δlog L| ≤ 3.0e-5, |ΔX_c| ≤ 8.6e-5; age-matched Δlog L up to 3.1e-4 at 26-46 Myr. On the integration branch, the 0.03 Msun dbl-to-ZAMS deck stops at 1.2163e7 instead of 1.2255e7 yr (final Δlog L +2.0e-3) and the 1 Msun ZAHB deck ends at Δlog L +5.9e-4.

### M40. With in_atmosphere, the Saha path left its first derivatives stale

- **Status:** fixed in [cd3b96ad](https://github.com/parallelpro/yrec/commit/cd3b96ad) (2026-10-29).
- **In the original F77?** Yes; EQSAHA jumps past its derivatives when `LATMO` is set ([eqsaha.f:173-174](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqsaha.f#L173-L174)):
  ```fortran
        IF(LATMO) GO TO 60
  C COMPUTE FIRST DERIVATIVES
  ```
- **Now:** `compute_saha_eos`, [yale_eos.f90:184](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/yale/yale_eos.f90#L184).
- **What was wrong:** dln ρ/dln T, dln ρ/dln P, c_P and ∇_ad kept whatever the caller's record held (and in the ramp that was blended in), although `eos_results` says every field is overwritten. The envelope integration, semiconvection's zone-edge ∇_ad and the starting model's surface pass `.true.` and read them.
- **What we did:** the first derivatives are always computed; `in_atmosphere` skips only the second ones. `test_eos` compares them with a full evaluation. All 56 decks are byte-identical: no shipped deck read the stale values.

### M41. The OPAL 2001 and 1995 density searches could start from an unbracketed pair

- **Status:** fixed in [86d0b206](https://github.com/parallelpro/yrec/commit/86d0b206) (2026-10-29).
- **In the original F77?** Yes. RHOOFP01 and RHOOFP shrink the lower trial density once, run regula falsi for 10 tries and return −999 silently ([rhoofp01.f:86-93](https://github.com/parallelpro/yrec/blob/6cd56739/src/rhoofp01.f#L86-L93), [:130-131](https://github.com/parallelpro/yrec/blob/6cd56739/src/rhoofp01.f#L130-L131); [rhoofp.f:75-81](https://github.com/parallelpro/yrec/blob/6cd56739/src/rhoofp.f#L75-L81)):
  ```fortran
          if(p1 .gt. pnr) then
            p2=p1
            rhog2=rhog1
            rhog1=0.2D0*rhog1
  ```
  The 2006 search was fixed in the 2026 bug sweep (batch 2).
- **Now:** `rhoofp01`, [opal_eos.f90:243](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/opal/opal_eos.f90#L243); `rhoofp`, [opal95_eos.f90:85](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/opal/opal95_eos.f90#L85).
- **What was wrong:** a pair that did not bracket the pressure made the search give up, and the blend fell back to the Yale/SCV result with `ierr = 0` and no message.
- **What we did:** 2006's bracket loop (at most 20 shrinks, down to the table's 1e-14), up to 30 regula falsi steps, and a run-log message when it falls back. No shipped deck uses the OPAL 1995 or 2001 EOS: 50 decks byte-identical.

### M42. beta ≤ 0 only warned

- **Status:** fixed in [20d05da3](https://github.com/parallelpro/yrec/commit/20d05da3) (2026-10-29).
- **In the original F77?** Yes ([eqstat2.f:126-134](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqstat2.f#L126-L134)):
  ```fortran
        if(BETA .LT. 0d0) then
          write(ishort,*) 'EQSTAT2: BETA is negative: TL,PL,BETA,1-BETA',
       *        TL,PL,BETA,1D0-BETA
        endif
        BETAI = 1.0D0/BETA
  ```
- **Now:** `blend_eos`, [blended_eos.f90:192](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/blended_eos.f90#L192).
- **What was wrong:** where the radiation pressure reaches P, the solve went on with a negative (or, at β = 0, infinite) 1/β.
- **What we did:** the point goes to the run log and `blend_eos` returns `ierr = 1`: there is no gas pressure to solve for. `test_eos` checks it. No deck reaches it: 50 decks byte-identical.

### M43. The SCV envelope table was used at a composition it was not built for

- **Status:** fixed in [72de2124](https://github.com/parallelpro/yrec/commit/72de2124) (2026-10-29). Latent.
- **In the original F77?** Yes; EQSCVE tests the zone against the envelope composition of `COMMON/COMP`, not against the composition the table was built for ([eqscve.f:27-29](https://github.com/parallelpro/yrec/blob/6cd56739/src/eqscve.f#L27-L29)):
  ```fortran
        IF(ABS(X-XENV).GT.1.0D-5 .OR. ABS(Z-ZENV).GT.1.0D-5)THEN
           CALL EQSCVG(TL,T,PL,P,DL,D,X,Z,BETA,BETAI,BETA14,FXION,RMU,
  ```
- **Now:** `compute_scv_envelope_eos` and `build_scv_envelope_table`, [scv_eos.f90](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/scv/scv_eos.f90).
- **What was wrong:** `eos_set_mixture` changes the envelope composition without rebuilding the table, and before the first build the table is all zeros, so a zone matching the new composition would read a table for another one.
- **What we did:** the table records its X and Z (−1 until built), and the lookup tests the zone against those; a zone at another composition takes SCV's general path. No shipped path reached it (every run rebuilds the table whenever it sets the mixture): 56 decks byte-identical. `test_eos`'s new check fails without the change.

### M44. A second eos_init could ask about a freed unit number

- **Status:** fixed in [33d0dc45](https://github.com/parallelpro/yrec/commit/33d0dc45) (2026-10-29).
- **In the original F77?** No; the OPAL table's unit is a `newunit` since [c00a9a20](https://github.com/parallelpro/yrec/commit/c00a9a20) (2026-10-11).
- **What was wrong:** the vintage's reader closes the table after reading it, but `opal_shared%table_unit` kept the number, and the next `eos_init` asked `inquire` about it: a freed `newunit` number that has not been reused is a runtime error. `test_eos` met it once its messages stopped opening a unit.
- **What we did:** the readers reset `table_unit` to 0 after the close ([opal_eos.f90:957](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/opal/opal_eos.f90#L957), [opal95_eos.f90:1057](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/opal/opal95_eos.f90#L1057)). Byte-identical.

### M45. eos_get_gamma1 built Γ₁ from derivatives it had asked not to compute

- **Status:** fixed in [341d222a](https://github.com/parallelpro/yrec/commit/341d222a) (2026-10-29). Only `test_eos` calls it.
- **In the original F77?** No; the routine came with [4b54a7c0](https://github.com/parallelpro/yrec/commit/4b54a7c0) (2026-08-21).
- **Now:** [eos_lib.f90:145](https://github.com/parallelpro/yrec/blob/a568e084/src/eos/eos_lib.f90#L145).
- **What was wrong:** its Yale/SCV branch passed `in_atmosphere = .true.` and then used dln ρ/dln P, dln ρ/dln T and c_P, which that flag skipped (1/0 on a fresh record; `test_eos` had left its two cool points out); outside the OPAL 2006 table its OPAL branch returned what the previous point left.
- **What we did:** it passes `.false.`, and outside the OPAL table it takes the Yale/SCV branch; `test_eos` pins all four points and checks the fallback.

### M46. test_eos never built the SCV envelope table, so its SCV values were Saha's

- **Status:** fixed in [c71c3ba0](https://github.com/parallelpro/yrec/commit/c71c3ba0) (2026-10-29). Test only.
- **In the original F77?** No; the pins came with [8cd2d244](https://github.com/parallelpro/yrec/commit/8cd2d244) (2026-10-11).
- **What was wrong:** `eos_init` reads the SCV tables but does not build the envelope-mixture table, which a run builds when it reads the starting model (`eos_update_envelope_tables`). `test_eos` never built it, so every SCV lookup at the envelope composition met an all-zero table, reported the point outside it, and the blend kept the Saha result: the pinned "SCV" values were Saha's.
- **What we did:** `test_eos` builds the table after `eos_init`; 46 pinned lines change. No code outside the test changed.

### M31. use_structure_dt_limits is switched off for the rest of the run

- **Status:** open (found 2026-10-29 while building M28's deck).
- **In the original F77?** Yes. HTIMER turns `LPTIME` off for good above log T_c 7.1, and MAIN saves and restores it only for the solar calibration ([htimer.f:65-71](https://github.com/parallelpro/yrec/blob/6cd56739/src/htimer.f#L65-L71), [main.f:261](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L261), [:955](https://github.com/parallelpro/yrec/blob/6cd56739/src/main.f#L955)):
  ```fortran
        IF(LPTIME) THEN
           IF(HT(1).GT.7.1D0 .AND. TLUMX(7).LT. 0.0D0) THEN
              WRITE(*, 100) HT(1), TLUMX(7)
   100        FORMAT('LPTIME SET FALSE - TC ',F7.4,' EGRAV ',E10.2)
              LPTIME = .FALSE.
  ```
- **Now:** `s%job%use_structure_dt_limits` is switched off in `begin_kind_card` ([run_driver.f90:399](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/run/run_driver.f90#L399)) and `evolve_step` ([step_driver.f90:198](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/step/step_driver.f90#L198)), and restored only per calibration cycle ([run_driver.f90:238](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/run/run_driver.f90#L238), [:271](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/run/run_driver.f90#L271)).
- **What is wrong:** the same class as M28. A later card that re-reads the pre-main-sequence starting model runs without the structure-based timestep limits. Seen: the four-card variant of `cm_settle_he` with card 3 re-reading the starting model after 1 Gyr diverges at the start of card 4 (NaN gradient in the temperature-gradient routine, then "LOG T OF 12.56 OUT OF OPAL 95 TABLE RANGE"), with and without M28's fix. `use_mass_accretion` is the same kind of latch, documented as "for the rest of the run" ([mass_change.f90:125](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/step/mass_change.f90#L125)).
- **Expected behaviour, and why:** each card starts with the deck's switch.
- **Proposed fix:** as for M28, a star latch armed from the control at every card start. One measured commit; only decks that re-read a pre-MS model after an evolve card move, and none is shipped besides the calibration, which already restores it.

### M32. Three controls are stored and used nowhere

- **Status:** open, documented: the registry's `used` column says `none` for each, and `tools/check_boundaries.py` keeps that true.
- **In the original F77?** Yes. `ZXMIX` is set from the mixture tables and read by no computation ([parmin.f:1065](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L1065)); `LENVG` is only printed ([parmin.f:1164](https://github.com/parallelpro/yrec/blob/6cd56739/src/parmin.f#L1164)); `PMMMD` appears only in COMMON statements.
- **Now:** the rows `zxmix` (`target_metal_fraction`), `lenvg` (`envelope_generation_flag`) and `pmmmd` (`pmm_norm_mdot`) of [controls_registry.tsv](https://github.com/parallelpro/yrec/blob/396bbb7c/src/defaults/controls_registry.tsv).
- **What is wrong:** a deck can set the CUS mixture's Z/X, the envelope generation flag and the Matt law's solar mass-loss rate, and none of them changes anything: the mixture change keeps the starting model's Z, the outer layers are built by the stitched model, and the mass-loss rate's use was the commented-out rotation-scaled solar wind.
- **Proposed fix (the authors' choice):** wire each to what it names, or retire it as `mu_gradient_scale` was (status `unused`: the legacy decks may still set it, the inlist format no longer offers it). Either is byte-identical for the shipped decks unless wired.

### M33. test_element_diffusion reads stale metal coefficients

- **Status:** open (found 2026-10-29: the unit tests failed twice in 14 runs, in `test_element_diffusion` only). Already present at 2c88d896: its test binary differs from the pin in 3 of 100 runs.
- **In the original F77?** The branch, yes; the stale values, no. SETUP_GRSETT sets the metal coefficients only under `LTHOUL` ([setup_grsett.f:248-272](https://github.com/parallelpro/yrec/blob/6cd56739/src/setup_grsett.f#L248-L272)), but keeps them in `COMMON/GRAVSZ`, which nothing else writes, so an F77 run without Thoul's coefficients settles no metals, the same way every time:
  ```fortran
           IF(LDIFZ)THEN
              FAC=FGRZ*HRU(I)**2*HTU(I)**2.5D0/LN_LAMBDA
              IF(LTHOUL)THEN
  ...
                 COD1Z(I) = CSTDIFFMIX*CSTMIXING*FAC*HQPR(I)*ZZ*(AP+AT)
  ```
- **Now:** `gravitational_settling_setup`, [original_settling.f90:616-640](https://github.com/parallelpro/yrec/blob/396bbb7c/src/element_diffusion/original_settling.f90#L616-L640).
- **What is wrong:** with `use_diffusion_z` on and `use_thoul_diffusion` off, the metals settle with whatever `diff_scr`'s metal coefficients hold: the previous call's values, or memory not written in this run. The unit test's case "old: without the Thoul et al. coefficients" is that combination, so its dY/dZ split varies from run to run and its pin holds one draw. No shipped deck sets `use_diffusion_z` without `use_thoul_diffusion`.
- **Expected behaviour, and why:** the same input gives the same output; the Bahcall & Loeb scheme has no metal coefficients of its own.
- **Proposed fix (element_diffusion's, and the authors' choice):** zero the metal coefficients outside the Thoul branch (the F77's behaviour), or reject or warn on the combination in `controls_check`; then repin the test case. One measured commit; no shipped deck moves.

### M34. The Monte Carlo summary prints zero calibration derivatives when a run converges at its first check

- **Status:** open (found 2026-10-29). Cosmetic: one record of the Monte Carlo summary file.
- **In the original F77?** Yes. CHKCAL returns on convergence before it sets the fixed derivatives, and WRTMONTE prints them from `COMMON/CALSUN` ([chkcal.f:25-44](https://github.com/parallelpro/yrec/blob/6cd56739/src/chkcal.f#L25-L44), [wrtmonte.f:39](https://github.com/parallelpro/yrec/blob/6cd56739/src/wrtmonte.f#L39)):
  ```fortran
        IF(ABS(BL).LT.TOLL .AND. ABS(RL).LT.TOLR)THEN
  ...
              LSOL = .TRUE.
              GOTO 9999
  ...
        DLDX = -3.78D0		! empirical result:  -3.783    RMS error .070
  ```
- **Now:** `check_solar_calibration`, [calibration.f90:53-87](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/run/calibration.f90#L53-L87), which reads its own parameters `dlogl_dx_fixed` ... since [368078df](https://github.com/parallelpro/yrec/commit/368078df); record 2 of the summary, [monte_carlo.f90:195](https://github.com/parallelpro/yrec/blob/396bbb7c/src/core/run/monte_carlo.f90#L195), prints the star members `dlum_dx` ... `drad_dalpha`.
- **What is wrong:** in the first Monte Carlo run of a process, a calibration that converges at its first check prints 0 for the four derivatives (later runs print the previous run's, the fixed values). The values used are always the fixed ones.
- **Proposed fix:** print the calibration's parameters in record 2 and delete the four members. That record changes only in this case; no deck exercises it.

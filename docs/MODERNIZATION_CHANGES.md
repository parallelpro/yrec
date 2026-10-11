# YREC modernization: what changed numerically, and what changed functionally

Baseline: the original fixed-form F77 import, commit `6cd56739` (2025-04-10).
Head: branch `yrec-modern` at `c994fa1d` (2026-09-04). 728 commits between them.

Sources: the full commit history (subjects and bodies), `docs/ROADMAP.md`,
`docs/GUIDELINES.md`, `audit/*`, the in-source "preserved" / "fixed here"
annotations, and a file-level diff of the two trees. Commit hashes are given so
every claim can be checked with `git show <hash>`.

Two headline facts:

1. **Almost all of the 437 modernization commits (2026-08-19 onward) are
   byte-identical refactors**, verified on every commit against a strict
   byte-diff of a legacy solar run. Roughly 390 of 437 are pure refactors.
2. **The numbers did change**, in a small number of deliberate, documented
   batches: 18 inherited physics bugs fixed in the 2026-09-02/03 bug-sweep
   batches, five output-only pulse/profile physics fixes, one reproducibility
   change to compiler flags, and a handful of 2025 fixes by other authors.
   Every one of these reseeded the regression pins and recorded the measured
   drift in its commit message.

---

## 1. Timeline

| Phase | Dates | Commits | Who | What |
|---|---|---|---|---|
| Import and warning cleanup | 2025-04 to 05 | 43 | Cao, Rendina, van Saders, Pinsonneault | F77 import, `run_tests.sh`, multi-compiler Makefile, ENTIME timestepper fixes, gfortran warning sweep |
| Test infrastructure | 2025-06 to 08 | 22 | Rendina, Cao | pytest runner with reference capture, ASan crash fixes, ifx support |
| Stitch and IMPJMOD | 2025-09 to 10 | 32 | van Saders, Tayar, Pinsonneault | Stitched `.store`, IMPJMOD rotation option, braking-law typo, `rhoofp06` tolerance |
| ZAHB/HB and guards | 2025-11 to 12 | 40 | Rendina, Byrom, Zinn, Tayar | LROT/OMEGA guards, ZAHB decks, HB opacity gate, semiconvection toggle |
| Docs and codespace | 2026-01 | 9 | Rendina, van Saders | Install docs, devcontainer, FGRY/FGRZ |
| Release-paper prep | 2026-02 to 05 | 145 | Patton, Rendina, others | Env-var paths, version stamp, namelist standardization, new opacity and atmosphere tables, manuscript |
| Free-form transliteration | 2026-08-19 to 20 | 74 | Li | All 253 `.f` files to free-form `.f90` with readable names, byte-identical on four solar cases |
| COMMON elimination and architecture | 2026-08-21 | 63 | Li | Zero COMMON blocks, `star` derived type, eos/kap/atm facades, ierr, CI |
| Legacy campaigns and MESA conventions | 2026-08-22 | 42 | Li | 504 DO loops, 249 SAVEs, 698 gotos removed, libyrec/pyyrec, controls registry, MESA history |
| Output and controls campaigns | 2026-08-23 to 24 | 53 | Li | Profiles, FGONG/GSM, extended model, `star%ctrl`/`star%job`, 115-file rename |
| Regression net and legacy retirement | 2026-08-25 to 28 | 44 | Li | 16-config matrix, mass-loss fix, stitched-model restructure, `.mod` format, old outputs retired, run log |
| Seismic, stops, reproducibility | 2026-08-29 to 30 | 62 | Li | Structure-limit stops, nu_max/delta_nu, `-ffp-contract=off`, math_lib, envint purity |
| Audits and bug sweeps | 2026-08-31 to 09-03 | 24 | Li | Two independent sweeps, Batches 0 to 3 of physics fixes, TAHB NaN |
| Readability waves | 2026-09-04 | 75 | Li | R1 to R3, W2, W3 per domain, all byte-identical |

---

## 2. Numerical changes

Legend for the **Origin** column: *F77* = bug inherited from the original
code; *Modern* = regression introduced during modernization and fixed;
*Deliberate* = intentional change to defaults, precision, or flags.

### 2.1 Cumulative effect on the standard cases

The best single measure of how far today's numbers sit from the original code,
for the cases the project pins:

| Case family | Where the numbers moved | Size |
|---|---|---|
| Calibrated solar models (12 testsuite runs) | Batch 2 (dominated by the d(p,g) fix) | log Teff <= 1.2e-5, log L <= 5.7e-5, X_c <= 2.5e-4, about 25 % fewer models to solar age |
| Same, Batch 3 on top | N14(a,g) coefficient and literal sweep | log Teff <= 5e-6, log L <= 3.3e-5, X_c <= 8.7e-5 |
| Pre-main-sequence starts (12 of 22 config-matrix cases) | Batch 2 d(p,g) fix | 0.02 to 0.31 dex in log L at model 2, D-burning onset moves from ~1.7e6 K to ~1e6 K |
| 0.3 Msun start to D-burning example | Batch 2 | Stops at 1.06e-3 Gyr instead of 5.6e-3 Gyr |
| dbl to ZAMS runs | all batches | byte-identical |
| Rotating solar case | Turnover-timescale semantics (b116224d) | final log Teff ~2e-4 relative |
| Three difrotmix rotation cases | Rotation-path batch (07476036) | first difference at model 2, final log Teff within 2e-5, surface X within 0.2 % |
| ZAHB to TAHB solar case | TAHB NaN fix plus Batch 2 | previously died with NaN at model 977; now runs to its stop at model 1108, 0.119 Gyr |
| Every case on FMA hardware | `-ffp-contract=off` (89c25e60) | last-bit ULP drift |
| Every case | e-notation literal sweep (d853a29c) | solar age column moves in the 9th significant digit |

Not covered by any pin, so unmeasured: OPAL95 EOS runs, MHD EOS runs,
`lvfc` rotating runs, Itoh neutrino runs above logT 7.5, He-burning cores
using the alpha-capture Jacobian, conductive-opacity-dominated degenerate
cores, mass-loss runs (which never ran at all before).

### 2.2 Nuclear network and energy generation

| Commit | Date | Origin | What was wrong | Impact |
|---|---|---|---|---|
| 0304e269 (Batch 2, item 2) | 2026-09-02 | F77 | d(p,g)3He rate used T9^(+2/3) instead of CF88's T9^(-2/3), in both `engeb` and `deutrate` | Rate x1e4 at 1e6 K, x260 at the solar centre. D lifetime at logT 6.2 goes from 6000 yr (coded) to 0.6 yr (correct). Dominant driver of every pre-MS shift above. |
| 0304e269 (item 1) | 2026-09-02 | F77 | Alpha-capture T-derivatives (C13(a,n), C12(a,g), N14(a,g)) had the wrong prefactor (nuclear term ~1e-22 too small), missing T9 powers | dln(rate)/dlnT was ~-0.1 to -0.5 where 15 to 40 is correct. Jacobian only, plus `rot_scr%dlnepsilon_dlnt` and the pulse eps_T column, above logT 7.73. Convergence in He-burning regime. |
| 0304e269 (item 3) | 2026-09-02 | F77 | Itoh neutrino branch: `dlnepsilon_dlnt` was a copy-paste of the rho line, and mixed log-derivatives into absolute sums | Jacobian and rotation/pulse consumers with Itoh losses above logT 7.5. Removes a 0/0 when snu = 0. |
| 0304e269 (item 4) | 2026-09-02 | F77 | `liburn`/`liburn2` radiative fraction ran over (-1, 0] instead of [0, 1] for zones behind a retreating convection zone | Li6/Li7/Be9 depletion exponent wrong by factors of several in pre-MS runs with extended composition. Every such pinned case moves. |
| 0304e269 (item 5) | 2026-09-02 | F77 | `sneut`: `cvp*cvp - cap - cap` where `cvp**2 - cap**2` is intended (the classic sneut typo MESA also fixed) | Pair, photo and bremsstrahlung neutrino losses change by tens of percent above logT 7.5 (RGB and He cores, massive MS). |
| 0ff84d7e (Batch 3) | 2026-09-03 | F77 | N14(a,g)F18 linear coefficient 0.177 was a digit transposition of CF88's 0.117 | ~0.5 % of that rate at T9 = 0.15. |
| 0ff84d7e (Batch 3) | 2026-09-03 | Deliberate | 362 unsuffixed single-precision REAL literals in `burn_lib`/`net_lib` promoted to double, including two `**(2./3.)` exponents | ~6e-8 relative in the constants; test_net moves up to 5.5e-6. |
| d853a29c | 2026-08-31 | Deliberate | 19 e-notation constants (rate coefficients, 1e7 T floor, sneut fit constants) promoted to d-notation | Solar age in the 9th digit; test_net at 1e-11. |
| 955a50c7 | 2026-08-23 | F77 | Central-abundance stop disarm negated the deuterium stop instead of the hydrogen stop when the start model was already below the X threshold | Runs with abundance stops halted after one model with a corrupted D stop. No pin affected. |
| b16c853c, d02bfee1 | 2025-05-31 | Deliberate (Cao) | `liburn` exp-overflow guard; one commit used an undefined variable, fixed next commit | Only when the exponent exceeds 709; abundance set to zero. |

### 2.3 Equation of state

| Commit | Date | Origin | What was wrong | Impact |
|---|---|---|---|---|
| 1bc94ae4 | 2025-10-10 | Deliberate (Rendina) | `rhoofp06` rho(P,T) inversion tolerance loosened 0.5e-7 to 1e-5 to stop cp growing without bound | Every OPAL06 run carried ~1e-5 relative, non-smooth noise in rho(P,T) from this date until Batch 2. |
| 40aa3197 (Batch 0) | 2026-09-02 | F77 | `esac`/`esac01`/`esac06` rescaled the E and cv slots by T6 and R/mu on every trial call, compounding stale values | Output-unchanged in itself (consumers re-call at full order), but this was the real cause of the cp blow-up that 1bc94ae4 papered over. |
| 0304e269 (item 7) | 2026-09-02 | Modern (2025) | `rhoofp06` tolerance restored to 0.5e-7, bracket verified, refinement cap 11 to 30 | rho(P,T) shifts up to 1e-5 relative in all 311 OPAL06 inlists; test_eos reseeded at 1e-7. |
| 0304e269 (item 6) | 2026-09-02 | F77 | `eqscve` SCV upward-neighbour smoothing weight was `0.5*d/tol` instead of `0.5*(1-d/tol)`, in logT and logP | Step noise in rho, cp, grad_ad for every SCV-EOS run: up to 5 % in dlnrho/dlnT and 2.7 % in grad_ad in M-dwarf and brown-dwarf photospheres. Solar pins unchanged. |
| 40aa3197 (Batch 0) | 2026-09-02 | Modern | `esac` OPAL95 first-call guard collapsed to a single-statement IF, so every OPAL95-EOS run died at its first OPAL point | Restores the OPAL95 EOS path. |
| 40aa3197 (Batch 0) | 2026-09-02 | Modern | `intpt` table scan's goto-exit had become `return`, leaving outputs unset; killed the whole MHD EOS path at the first evaluation | Restores the MHD path. Zero runtime coverage, so unverified beyond compiling. |
| 31dc110b | 2026-08-21 | F77 | `starin` called both `meqos` and `eqstat` when MHD was on (missing ELSE); MHD dispatch extended to `hsubp`, `sconvec`, `massloss` which never checked the flag | MHD-on runs only. No reference output exists. |
| 0ff84d7e (Batch 3) | 2026-09-03 | F77 | `rhoofp06` priming call passed a radiation flag that caused a 0/0 in a discarded slot | No numeric change; removes an FPE trap hit. The `rhoofp01` twin is still open. |

### 2.4 Opacity

| Commit | Date | Origin | What was wrong | Impact |
|---|---|---|---|---|
| 0304e269 (item 9) | 2026-09-02 | F77 | OPAL95 table slot layout assumed packed while readers used unpacked indices, so the (X=0, Z=0.08) and (X=0, Z=0.10) tables were overwritten by X=0.94/0.92 tables | A He-burning core got kappa about 0.28 dex too high at electron-scattering temperatures. All TAHB/ZAHB/He-ignition cases move physically. Solar, MS and RGB untouched. |
| 0304e269 (item 8) | 2026-09-02 | F77 | Conductive opacity: dlnkap/dlnT was built from the rho derivatives, and both derivatives were multiplied by log conductivity (~12 to 17) instead of conductivity-weighted | kappa unchanged; Jacobian corrupted whenever Potekhin conduction is on. Weight ~5e-4 in solar and PMS (pins unchanged) but dominant in degenerate cores. Also feeds rotating-RGB diffusion coefficients and the pulse kap_T column. |
| 5b1a5ed6 | 2025-12-17 | Deliberate (Zinn) | High-Z opacity ramp gate lowered from logT 7.7 to 7.0 | Z > 0.1 zones use OPAL below 1e7 K instead of 5e7 K; HB models with semiconvection plus overshoot. |

### 2.5 Mixing and convection

| Commit | Date | Origin | What was wrong | Impact |
|---|---|---|---|---|
| fb28dacf | 2026-08-31 | F77 | `find_convection_zones` wrote the top of radiative region 1 with the last zone whenever the outermost CZ was detached from the surface | Byte-identical on all pins; affects only detached-outer-CZ structures. |
| fb28dacf | 2026-08-31 | F77 | `rezone` flag-point appends could write past the array on heavily structured models | Byte-identical on pins. |
| 0304e269 (item 14) | 2026-09-02 | F77 | `model_to_equal` built the metal-settling run from the N15 slot instead of Z | Latent on today's paths; hygiene fix, no pin moved. |
| a093427f | 2025-05-31 | Deliberate (Cao) | `slopes` safe-divide returns 0 when numerator or denominator is 0 | Seven 2025 test cases changed (m10000/m16000/m6000 ZAMS, m25000 pre-MS, two D-burning). Old code occasionally gave zero core radius and turnover times off by orders of magnitude. |
| a8e60674 | 2025-12-17 | Deliberate (Zinn) | Semiconvection restricted to the core zone only | Semiconvection runs only. |
| 5bcbc1fe, 7ae2d8c5 | 2025-05/06 | F77 | Index guards and a missing DIMENSION in `starin` that let `convec` write 44 bytes past an array | Memory-corruption fixes; the pre-fix effect on outputs was never quantified. |

### 2.6 Rotation, angular momentum, diffusion

| Commit | Date | Origin | What was wrong | Impact |
|---|---|---|---|---|
| 0304e269 (item 9b) | 2026-09-02 | Modern (94c7f45c, 2026-08-20) | `zahn_coupling_factor` gate switched from `lvfc` to `ldifad`, silently zeroing all rotational mixing when `velocity_form_correction` is on and advection-diffusion is off | Restored. lvfc-mode rotating runs; no pin sets lvfc. Live for 13 days. |
| 0304e269 (item 10) | 2026-09-02 | F77 | Zahn's alpha = 1/2 dln(r^2 Omega)/dln r was coded without the 1/2 | Order-unity in the velocity-form correction for lvfc runs. Any historically tuned FC absorbed this. |
| 0304e269 (item 11) | 2026-09-02 | F77 | `compute_quadrupole` drho/dr stencil used zone i where i-1 belongs, collapsing to half a one-sided difference | Quadrupole moment 10 to 20 % in the solar radiative interior; Eddington-Sweet correction a few to 10 %; surface omega and Li at percent level. All rotating pins move. |
| 0304e269 (item 12) | 2026-09-02 | F77 | `trapzd` interpolated rho, omega^2, eta^2 at a loop-invariant offset instead of at y | Rotation shape integrals biased 1 to 3 %; f_P by 1e-8 (Sun) to 1e-6 (rapid rotators). Bit-level on rotating pins. |
| 0304e269 (item 13) | 2026-09-02 | F77 | `envelope_derivs` divided dm/dlogP by f_P but multiplied dr/dlogP by f_P, violating mass continuity by f_P^-2 | Non-rotating bit-identical (f_P = 1). Rotating: R and Teff ~1e-6 for the Sun, ~1e-4 for pre-MS and rapid rotators. |
| 07476036 | 2026-09-02 | F77 | `rezone` interpolated old eps/esum onto the new grid after the grid had already been replaced (mis-registered wherever points were inserted) | Feeds Eddington-Sweet velocities every rotating timestep. |
| 07476036 | 2026-09-02 | F77 | `diffuse_composition` upper bound never set when the diffused region's top is radiative or the surface, so sums ran over 1..0 | Rotating runs. |
| 07476036 | 2026-09-02 | Modern | Mid-step deuterium rate arrays relied on the F77 blanket SAVE; after the SAVE campaign they were ~0 after the first substep | Rotating pre-MS D burning. |
| 07476036 | 2026-09-02 | Modern | `rebuild_envelope` read a stale local omega and dropped the `omega_from_j` result | Rotating envelope seed. First difference at model 2 in the difrotmix cases. |
| b116224d | 2026-08-25 | Deliberate | Turnover timescale computed once per step from the stitched model; wind saturation and D limiter read the step-start value instead of a mid-step recompute | Non-rotating byte-identical; rotating solar final log Teff ~2e-4 relative. |
| 955a20da | 2026-08-31 | Modern (212a652c, same day) | Two callers of `rotation_shape_factors` missed the new `ierr` argument, so the callee wrote through a garbage stack slot on every Newton iteration of rotating runs | Undefined behaviour; happened to leave pins byte-identical on this compiler. |
| a90c0cad | 2025-09-30 | F77 (van Saders) | Braking-law radius exponent typo: an extra `*PMMB` factor | All rotating runs with PMM wind loss. |
| d5cd30e7 | 2026-02-23 | F77 (van Saders) | Surface-CZ search in `tauintnew` snagged on models whose last points are radiative | Turnover timescale and wind saturation. |
| dbdb7b89 | 2025-09-03 | Deliberate (Tayar) | Moment of inertia recomputed before rotation in `hpoint` | Rotating runs; author noted uncertainty about necessity. |
| d91c3a81 and following | 2025-09 | New physics (Tayar, Morales) | IMPJMOD option: solid body, constant j, or power-law omega in the CZ | New option, off by default. |
| eeeb2420 then c95697e1 | 2025-12 / 2026-02 | Reverted | `walpcz` gate flipped `.GE.` to `.GT.` and back | Net zero; rotating runs between those dates differed. |
| 224c9f0d | 2026-01-26 | New physics (van Saders) | FGRY/FGRZ multipliers in the new diffusion routines | Only when set to values other than 1. |

### 2.7 Wind and mass loss

| Commit | Date | Origin | What was wrong | Impact |
|---|---|---|---|---|
| 50944f10 | 2026-08-25 | F77 | `massloss` passed 24 actuals into 23 dummies; LMDOT runs dereferenced a scalar as the composition array and segfaulted. Mass loss had never worked. Behind that, the Reimers rate was computed in CGS and converted to CGS again (~1e25 inflation). | Mass-loss runs now run, with the standard Reimers law `Mdot = creim (L/Lsun)(R/Rsun)/(M/Msun)` Msun/yr. Also documented that PDISK was always rad/s and TDISK Gyr, contrary to 70 deck comments. |

### 2.8 Atmosphere, envelope, and the stitched model

| Commit | Date | Origin | What was wrong | Impact |
|---|---|---|---|---|
| 2be7c301 | 2026-09-02 | Modern (SAVE loss) plus F77 | `timestep_limit_heburn` read `helium_dt` before writing it; the F77 relied on the caller's blanket SAVE, the modern code got a zeroed local, so dt = 0 and NaN at core-He exhaustion. The shell branch now uses `time_dy_total`/`time_dy_shell`, which were mapped but never consumed. | Root cause of the TAHB NaN at model 977. Nothing else fires this branch; all pins byte-identical. |
| 5e98ee6a | 2026-08-26 | F77 | The `.store` stitch's own atmosphere integration clobbered the solver's adaptive step sizes from the first write onward | Legacy LSTCH runs only: first model identical, later models ~1e-7. |
| eccf7133 | 2026-08-25 | Deliberate | Stitched build made side-effect-free (pphot saved and restored) | History pphot on profile steps no longer depends on `profile_interval`. Pins byte-identical. |
| 979a9c26, 956a5ede | 2025-05 | Deliberate (van Saders) | ENTIME envelope-triangle timestepper rewritten | Only decks with `LTRIST = T` and `LPTIME = F`; magnitude unstated. |

### 2.9 Output-file-only physics fixes (evolution unchanged)

These change what profile and pulse files contain, and therefore seismic
frequencies computed from them, but not the evolution itself.

| Commit | Date | Origin | What was wrong | Impact |
|---|---|---|---|---|
| 2ed398ca | 2026-08-30 | F77 | Pulse N^2 used the composition-homogeneous thermal form, omitting the Ledoux mu-gradient term | On a 1.2 Msun subgiant the g-cavity N was low by ~30 %, the period spacing high by 12 % (581 s vs 510 s), l=1 mixed modes off by up to ~16 uHz. Also `kap_T`, `kap_rho`, `eps_T`, `eps_rho` were bare log-derivatives instead of absolute ones. |
| 67b36131 | 2026-08-30 | Follow-up | Hybrid N^2: thermal form kept inside convective zones, centred difference elsewhere | 26 of 81 convective points had spurious positive N^2; ~0.4 uHz bias even in radial modes. |
| afff605a | 2026-08-31 | F77 | Pulse builder swapped grad and grad_ad across the envelope region | p-mode shifts: median 0.43 / max 0.72 uHz (MS), 0.26 / 0.69 (subgiant), 0.10 / 0.16 (RGB). |
| b116224d | 2026-08-25 | F77 | Profile envelope columns 13/14 (gradT/grada) swapped | Every profile file's envelope region. |
| a23a4846 | 2026-08-29 | F77 (stitch.f) | Atmosphere radii added the raw per-step length instead of accumulating it; heights ~1e-5 R* instead of ~1e-3 and non-monotonic | GYRE rejected the files. Frequency effect up to ~0.5 uHz at the high-frequency end. |
| b4673635 (Batch 1) | 2026-09-02 | Modern | Profile `mu` column held R/mu; `mu_e_inv` was inverted twice; gradL had a spike at the envelope junction | Non-default columns; no frequency impact. |
| 1694dbc4 | 2026-08-23 | Modern | A stray RETURN made the history sources block dead: log_R, log_g, moments of inertia, rotation period, snu rates were all zero in `history.data` | Fixed. |
| 59e8b788 | 2026-08-29 | Modern | GSM writer stamped version 101 instead of 110 | GYRE rejected GSM files. |

### 2.10 Reproducibility, precision, and defaults

| Commit | Date | What | Impact |
|---|---|---|---|
| 89c25e60 | 2026-08-30 | `-ffp-contract=off` on gfortran, `-fp-model precise` on Intel | Forbids FMA contraction. Last-bit drift on FMA hardware including arm64; all pins reseeded; 314a4153 reseeded the four committed domain baselines it missed. |
| 71d1447b and following | 2026-08-30 | `math_lib`: all transcendental calls and real-exponent `**` routed through one module; optional `USE_CRMATH=1` correctly-rounded backend | Default backend byte-identical. crmath drifts ~1e-11, ~35 % slower, needed for cross-platform bit identity. |
| 031be858 | 2026-08-29 | `set_initial_model_number` (default 1) | First converged model is numbered 1 instead of the stored counter (~64 from the birthline). Changes model numbers everywhere and the profile cadence phase. All pins reseeded. |
| 64693ebc | 2026-08-29 | `solar_mass_cgs` and `G_cgs` become controls | Defaults reproduce the old constants bit-for-bit (`G_cgs = -1` keeps the historical log10 G = -7.17571). Overriding G to CODATA shifts log g by +3.5e-4. |
| c55710d3 to de95bfd6 | 2026-08-21 | CI comparator tolerance 5e-13 to 1e-9 relative | CI only (Linux x86 vs macOS arm64 ULP noise, worst observed 9.6e-12). Local discipline stays byte-exact. |
| 0ff84d7e | 2026-09-03 | NaN guard on every converged model | A run that goes non-finite now ends with exit 1 at the first NaN instead of running to its budget and exiting 0. |
| f0a430c4 to d8e864f7 | 2025-11/12 | `LROT` with `OMEGA(1) = 0` | First silently disabled rotation, then replaced by a STOP. |
| d5207b56, 2c0bfe91 | 2026-08-29 | Documentation corrections | `max_iter_level3` default is 2, not 0; `first_call_flag` default is true. Both mis-documented defaults had caused users to shift results at 1e-5 or restart cards unintentionally. |

### 2.11 In-process re-entrancy

| Commit | Date | What | Impact |
|---|---|---|---|
| 40aa3197 (Batch 0) | 2026-09-02 | `read_controls` DATA statements replaced by an entry-time reset | A second in-process run inherited any control the new inlist omitted. Fresh-process results unchanged. |

---

## 3. Functional changes

### 3.1 Source language and layout

- **Then:** one flat `src/` with 253 fixed-form `.f` files (42,862 lines),
  1,490 COMMON statements over 236 block names, `IMPLICIT REAL*8` everywhere,
  an unnamed 1070-line main program, argument tramps of up to 60 arguments,
  and a 32.8 MB compiled binary committed to git.
- **Now:** 234 free-form `.f90` files (54,387 lines), 55 modules, zero COMMON
  blocks, `implicit none` in every file, organized by domain
  (`atm`, `core`, `eos`, `io`, `kap`, `math`, `mixing`, `net`, `numerics`,
  `rotation`, `setup`, `state`, `util`) plus `defaults/` and `tools/`; the
  `wind` directory of the first reorganization was dissolved into
  `rotation` and `core` in 2026-10.
- **One derived type** `star` in `src/state/star_info_lib.f90`, with
  `star%job` (per-card state) and `star%ctrl` (controls, immutable after
  read). Named index constants for every species, energy channel, and
  gradient slot.
- **Star-blind kernels:** `envint_kernel` and `burn_lib`'s kernel take a
  config record and never touch `star`. Enforced by `tools/check_boundaries.py`.
- **Domain facades:** `eos_get`, `kap_get`, `atm_get` are the only public
  entry points, returning named-index result arrays in MESA's shape.
- **Generated code:** `deps.mk` (module dependency order), the `star%ctrl`
  definition and sync module, the inlist declaration and read includes, and
  the default output column tables are all generated by scripts under
  `src/tools/` and checked for freshness in CI.
- **Legacy constructs removed:** 504 numbered DO loops, 698 gotos, 249
  blanket SAVEs, 813 lines of commented-out F77, 720 unused labels, all
  byte-gated.

### 3.2 Input and controls

- **Single inlist** with `&star_job` and `&controls`, readable names, run as
  `yrec CASE.inlist`. The reader auto-detects style, so legacy `.nml1`/`.nml2`
  pairs still run unchanged.
- **Controls registry** `src/defaults/controls_registry.tsv` (384 rows) is the
  source of truth, mapping every legacy name to a canonical one with type,
  default, status, and MESA equivalent. Generated `star_job.defaults` and
  `controls.defaults` document every control.
- **Representative renames:** `cmixla` to `mixing_length_alpha`, `endage` to
  `max_age`, `nmodls` to `max_model_number`, `lrot` to `rotation_flag`,
  `lcals` to `calibrate_solar_model`, `kttau` to `atm_choice`,
  `alphac` to `overshoot_alpha_core`, `lsemic` to `use_semiconvection`,
  `creim` to `reimers_scaling_factor`, `dmdt0` to `mass_change`,
  `end_xcen` to `central_hydrogen_stop`, `htoler` to `solver_tolerances`.
- **Converter** `tools/upgrade_inlist.py` turns a legacy pair into an inlist,
  preserving values and comments; `test_inlist_convert.py` proves both forms
  give byte-identical output. All 1,163 legacy decks in the repo were
  converted; four `.nml1` pairs are kept as legacy-reader fixtures.
- **Path placeholders** `{YREC_INPUT}`, `{YREC_START}`, `{YREC_OUTPUT}` from
  environment variables replace hard-coded relative paths.
- **Controls added:** `history_columns_file`, `profile_columns_file`,
  `profile_interval`, `write_profile_flag`, `write_pulse_flag`,
  `pulse_format`, `pulse_gyre_interval`, `terminal_interval`,
  `report_solver_diagnostics`, `inlist_used_file`, `set_initial_model_number`,
  `compute_neutrino_fluxes`, `solar_mass_cgs`, `G_cgs`, `Teff_sun`,
  `nu_max_sun`, `delta_nu_sun`, eight structure-limit stops
  (`log_L`, `Teff`, `log_g`, `nu_max` upper and lower), `impjmod`.
- **Controls removed:** the whole legacy output-selector family (`ldebug`,
  `lcorr`, `lmilne`, `ltrack`, `lstore`, `lstpch`, `lscrib`, `lstatm`,
  `lstenv`, `lstmod`, `lstphys`, `lstrot`, `lpulse`, `ipver`, `itrver`,
  `npoint`, `nprtmod`, `nprtpt`, `poa/pob/poc/pomax/lpout`) and their file
  paths; `lnewtcz` (legacy turnover path); `lacout`/`ageout` (acoustic-depth
  machinery); the transitional `use_legacy_output`.
- **Startup sanity checker** warns about inert or self-cancelling
  combinations (overshoot alpha set with the flag off, both calibration modes
  on, and so on).

### 3.3 Output

- **`history.data`** in MESA layout, 88 available columns, 30 by default,
  selectable through `history_columns_file`; flushed per row.
- **`profile<N>.data`** every `profile_interval` models, 60 available columns,
  12 by default, surface first, over the stitched interior plus envelope plus
  atmosphere point set.
- **Pulse files** `profile<N>.data.GYRE` (schema 101), `.FGONG` (ivers 1300,
  ported from MESA 26.04.1), `.GSM` (HDF5, needs `USE_HDF5=1`), all over the
  stitched model with seismic columns `brunt_N2`, `lamb_S2`, `gradL`,
  `csound`.
- **`.mod` model file** replaces `.last`/MODEL2: named header fields, full
  double precision (MODEL2 truncated to 9 digits), one model per file.
  Legacy MODEL2 and YREC7 start models still load.
- **Run log** in MESA style: a compact progress line every `terminal_interval`
  models, solver diagnostics behind a flag, an end-of-run summary with the
  termination reason. A solar log went from 16,713 lines to about 30.
- **`inlist_used`** verbatim copy of the inputs replaces the namelist echo.
- **Version and git hash** baked in at build time and printed in the banner.
- **Dropped:** `.track`, `.short` (renamed `.log` and slimmed), `.store`,
  `.pmod`/`.penv`/`.patm`, `.milne`, `.modpt`, `.debug`, `.snu`, `.excomp`,
  `.dyn`, `.calcad`, `.full`, and all stray `fort.*` writers.

### 3.4 Error handling and robustness

- **Then:** about 96 bare `STOP` sites, all exiting with status 0.
- **Now:** `ierr` is a required argument on all 16 library entries and on the
  integrator callbacks. Eleven `stop` statements remain, all in standalone
  test programs. Positive `ierr` exits 1; the historical "solution diverged"
  case returns negative `ierr` and still exits 0 by design.
- **NaN guard** on every converged model before it is written.
- **Stop conditions** consolidated into one module with a termination reason
  string; MESA-style structure limits added.
- **Re-entrancy:** `run_yrec` can be called repeatedly in one process with
  byte-identical results, verified by `test_reentry` and `test_pyyrec`.

### 3.5 Build

- **Then:** one `gfortran -w -O3 *.f` command, warnings suppressed, no
  incremental build, a duplicate Makefile for ifort, a bash Makefile
  generator.
- **Now:** one Makefile with `FC=gfortran|ifort|ifx`, per-file objects,
  generated dependency order, targets `default`/`release`, `devel`
  (bounds and FPE traps), `debug` (sanitizers), `prof`, `libs`, `lib`,
  per-domain static libraries, `USE_CRMATH=1`, `USE_HDF5=1`, an `install`
  script, a `Makefile.legacy` escape hatch, and a devcontainer.

### 3.6 Testing and CI

- **Then:** nothing at the first commit; a shell runner arrived two weeks later.
- **Now:** a 44-case `testsuite/`, a pytest harness with 14 test files, a
  16-config regression matrix, byte-exact MESA-format pins for two solar
  cases, standalone per-domain programs (`test_eos`, `test_kap`, `test_atm`,
  `test_net`) with committed baselines, contract tests for inlist conversion,
  re-entrancy, pyyrec, error exit, stop conditions, MESA output layout, and a
  loose physics sanity test. GitHub Actions builds, runs the domain tests,
  checks generated files, and runs the solar case twice to prove determinism.
- **Verification tiers** are codified in `docs/GUIDELINES.md`: Tier 1 (every
  commit) is a strict byte-diff of a legacy solar run; Tier 3 (full battery)
  for layout or format changes, with a documented reseed ritual.

### 3.7 Library and Python

- `libyrec.dylib`/`.so` with a C-interoperable `yrec_run` entry.
- `pyyrec` ctypes package: `pyyrec.run(control_nml, physics_nml)`.
- `wrappers/modelgrid_tools/` grid runners, track loaders, Slurm scripts.

### 3.8 Physics options

- **No physics prescription was dropped.** Every EOS, opacity set,
  atmosphere choice, mixture, wind law, diffusion, overshoot, semiconvection,
  spot, disk-locking, and S-factor control survives under a new name.
- **Added:** seismic observables `nu_max`, `delta_nu_rho`, `delta_nu`,
  `delta_Pg`; structure-limit stops; controls for G and Msun; the opt-in
  per-zone neutrino table; exact Ledoux `brunt_N2` and other seismic profile
  columns; the IMPJMOD rotation option; FGRY/FGRZ; a correctly-rounded math
  backend; new OP05 (AAG21, Magg22, Lodders25) opacity tables, Castelli-Kurucz
  atmosphere tables, and an OPAL EOS table generator.
- **Removed:** the legacy `taucal` turnover path, the acoustic-depth
  machinery (superseded by the `csound` column), and the HR-travel pulse-dump
  trigger.

### 3.9 Documentation

- New: `docs/GUIDELINES.md`, `docs/ROADMAP.md`, the controls registry and
  defaults files, column lists, `docs/HOW_YREC_WORKS.md`, 60 audit reports
  under `audit/`, per-case READMEs, wrapper READMEs, a 43-line README with
  citations, and the release paper.
- **Gap:** `docs/source/quickstart.md` and the Sphinx docs still describe the
  F77 workflow (`yrec file.nml1 file.nml2`, `.track`/`.last`/`.short`).
  They have not been updated for inlists, `history.data`, `.mod`, path
  placeholders, or the new build targets.

---

## 4. Quirks deliberately preserved byte-for-byte

These are known defects or oddities kept on purpose so the regression pins
stay valid. Each is annotated in the source.

- `rotation/magnetic_braking.f90` (formerly `wind_spindown_matt.f90`): the
  solid-body centrifugal normalization omits 1/G, where the
  convection-zone form divides by it. The original divided by
  exp(ln10*gl) with `gl` a never-assigned local (reads 0, so by exactly
  1). Sharpened by the sweep as a ~1.5e7 factor.
- `opal_eos_lib.f90`, `esac01/06`, `radsub01`, `rhoofp*`,
  `microdiff_coefficients.f90`: single-precision literals (`83.14511`, the
  OPAL x-grid, `0.5`, `2.2`) parsed as single then widened, kept verbatim.
- `burn_lib.f90` reaction-rate block kept inline because extracting it moves
  the ppII luminosity by 1 ulp.
- `shell_physics.f90` passes the EOS-returned density to `kap_get` rather
  than the model's, because the other choice drifts rotating runs.
- `esac.f90`: dE/dRho slot deliberately not rescaled (commented out in the
  original); `rhoofp`: radiation pressure deliberately not subtracted.
- `alfilein.f90`: the `latmtptau100` typo means the tau=100 guard for old
  NextGen files never works.
- `atm_lib.f90`: `alsurfp`'s error code deliberately not propagated.
- `observables_lib.f90`: the FX/FX2 typo block (a ratio that is identically 1).
- `envint_kernel.f90`: the first atmosphere step inherits the previous
  integration's tau/density/opacity, affecting only the outermost point.
- `stitched_model.f90`: seismic columns not recomputed on the interior-only
  branch.
- `find_convection_zones.f90` silently discards single-shell convection zones.
- `G_cgs = -1` sentinel keeps the historical rounded log10 G.

---

## 5. Open numerical issues (not yet fixed)

The bug sweeps left a triaged backlog in `docs/ROADMAP.md` sections 10 and 11
and in `audit/readability-sweep-2026-09-03/`. The ones most likely to matter
to users:

- **Atmosphere tables above ~48 kK:** `surfp`/`kcsurfp` spline over the -999
  sentinels every Castelli-Kurucz table's last row carries, so `atm_choice=5`
  at high gravity and temperature is silently garbage.
- **Advection-diffusion rotation mode:** shear and GSF equal-grid
  coefficients are never assigned, so those mechanisms contribute zero under
  `ldifad`; the "corrections too large" timestep cut became an inner-loop exit.
- **Gravitational settling:** the Coulomb log is applied twice (Thoul,
  `coulomb_log_choice` 2 or 3), a midpoint off-by-one, and convergence tests
  missing `abs()`. The Noerdlinger sign question is unresolved.
- **`rebuild_envelope`:** X/Z fitting-point interpolation has a reversed sign
  on one branch.
- **`map_user_inputs`:** S-factor derivative scales for reactions 2 to 7
  divide by the pp S0 instead of their own, zeroing the Seff derivative
  corrections under new nuclear rates.
- **`rescale_model`:** the Z-ramp block gives negative CNO if the target Z
  exceeds the old Z.
- **`check_solar_calibration`:** `dlum_dalpha = 0.0139` with a comment saying
  0.139; which is right needs the author.
- **Conductive opacity** low-density extrapolation branch flips the sign of
  log sigma.
- **OPAL95 EOS** (`esac.f90`) has no top-of-table guard at the coolest rows;
  its `radsub` applies the unit revision to P and S but not E and cv.
- **`eqstat`** ionization-cutoff blend adds raw Saha derivatives instead of
  the differenced form.
- **The Kawaler law's two forms** (`magnetic_braking`'s `set_kawaler_law`)
  disagree on the saturation gate.
- **Settling substep ceiling** in `mix` and `rotmix` has `mod` arguments
  reversed (divide by zero when `max_settling_dt` is an exact multiple).
- **`setups.f90`:** Avogadro digit transposition 6.0222137e23, affecting only
  the Debye-Hueckel coefficient at ~2e-5.
- **Numerical Recipes provenance** in `numerics_lib` (spline, LU, polint,
  Bulirsch-Stoer) is flagged as the single highest-liability item; each
  replacement will need a rebaseline.
- **No truth-level validation tier exists.** The N^2 bug survived every byte
  pin because pins freeze outputs, not physics. A cross-code MESA benchmark
  and a rotation benchmark are proposed but not built.

Unverified paths with zero regression coverage: MHD EOS (no tables ship),
OPAL95 EOS, `lvfc` rotation, mass loss, Itoh neutrinos, He-burning Jacobians.

---

## 6. Documentation inconsistencies to resolve

- `docs/HOW_YREC_WORKS.md` states the physics and numerics are unchanged and
  byte-identical throughout. That was true until 2026-08-29 and is false now.
- The single-precision literal policy is domain-dependent and undocumented:
  burn/net literals were promoted in Batch 3, while `microdiff_coefficients`
  and `opal_eos_lib` still say in capitals not to fix theirs.
- `docs/GUIDELINES.md` still defines the Stage-0 byte-diff in terms of
  `.short`, `.track`, and `.store`, all of which were retired.
- ROADMAP section 6 still lists the TAHB failure as needing investigation;
  section 11 records it fixed. Section 8 still lists the pulse column indices
  and the controls sanity checker as open; section 7 records them done.
- ROADMAP lists `atm_get_surface_pt` dropping `alsurfp`'s error as an open
  bug; the code documents the same behaviour as deliberate.
- `COMMENTS.md` (gitignored) reports as live the stop-disarm bug fixed in
  955a50c7.
- `wcz.f90` was recorded as deliberately kept dead code, then deleted in the
  readability sweep.

---

## 7. How byte-identity was verified

- Every change is classed A (byte-safe) or B (changes numbers) up front.
  Class B work is batched so one reseed covers many fixes, and the measured
  drift is recorded in the commit message.
- The per-commit gate is a strict byte-diff of a legacy solar run excluding
  only the git-hash line. The looser 2 % tolerance of the legacy pytest runner
  is explicitly not the gate, because it hides precision regressions.
- Full-model pins live under gitignored `standard/` directories and are
  machine-local. The only committed numeric baselines are the four domain
  test outputs. CI proves determinism by running the solar case twice rather
  than comparing to a committed baseline.
- Known trap recorded in the guidelines: the legacy runner seeds a reference
  from the run's own output if none exists, so a directory with no prior
  `standard/` is a smoke test reporting as a pass.

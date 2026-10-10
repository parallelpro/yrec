# How YREC Works

*A technical report on the equations, numerical schemes, and physics of the
Yale Rotating Evolution Code — written against the 2026 modernized source
(branch `modernize-readable-refactor`). File references use the current
layout; the physics and numerics are unchanged from the classic code
(Pinsonneault 1988; Guenther et al. 1992; van Saders & Pinsonneault 2012;
and the YREC release paper lineage), verified byte-identical throughout the
modernization.*

---

## 1. The problem YREC solves

YREC evolves a single, spherically averaged star through time as a sequence
of quasi-static models. The star is represented by ~1000–2000 concentric
mass shells (`star%nz`, capacity `json = 5000`). Each timestep the code:

1. advances the **composition** (nuclear burning, convective mixing,
   microscopic diffusion) at fixed structure,
2. relaxes the **structure equations** to a new hydrostatic equilibrium at
   the new composition and time,
3. transports **angular momentum** and applies rotational mixing (if
   rotation is on), iterating structure ⇄ rotation,
4. writes diagnostics and picks the next timestep.

This is the classic operator-split Henyey approach: composition and
structure are advanced in alternating half-steps rather than solved as one
coupled system, with a fourth iteration level available to converge the
splitting for high-precision (solar) work.

### 1.1 The structure equations

With mass `m` as the independent (Lagrangian) coordinate, YREC solves the
four standard stellar-structure ODEs, in log₁₀ form internally
(`logP`, `logT`, `logR`, and linear `L` as dependent variables, `log m` as
the coordinate; all cgs):

| # | Equation | YREC form (per shell, `core/coefft.f90`) |
|---|----------|------------------------------------------|
| 1 | Mass conservation | dr/dm = 1/(4πr²ρ) — residual `eq_r_val = exp(ln10·(log m − log ρ − log 4πr³))` |
| 2 | Hydrostatic equilibrium | dP/dm = −Gm/(4πr⁴) · **f_P** — residual `qp = −exp(ln10·(log G + 2 log m − log P − log 4πr³ − log r))·f_P` |
| 3 | Energy transport | dlnT/dm = ∇ · dlnP/dm — residual `qt = ∇ · qp`, with ∇ = ∇_rad (radiative) or ∇_conv (MLT) |
| 4 | Energy conservation | dL/dm = ε_nuc − ε_ν + ε_grav — assembled from `engeb` (burning + neutrinos) plus the entropy ("gravothermal") term −T ∂S/∂t evaluated with the stored per-shell rates of change from the previous step |

**Rotation enters the structure equations** through the Kippenhahn &
Thomas (1970) / Endal & Sofia (1976) form factors `f_P` (`star%fp_rot`) and
`f_T` (`star%ft_rot`): surface-averaged corrections computed by
`rotation/shape/fpft.f90` from ω(m) and the shell shape (oblateness via the
η₂ Radau parameter, `star%eta_squared`), multiplying the hydrostatic and
radiative-transport terms respectively. This is the standard 1.5-D scheme:
the star remains a 1-D object, but each shell is an oblate equipotential.

The **entropy term** deserves note: rather than differencing S directly,
YREC stores per-shell rates of change of (log T, log P, L, log R) from the
last converged step (`star%run%*_entropy_term`) and uses them both in the
energy equation and to extrapolate the first guess of the next model
(`improved_first_guess_flag` / LNEWS).

### 1.2 Composition equations

15 species are tracked per shell (`star%xa(1:15,k)`): H1, He4, Z (bulk
metals), He3, C12, C13, N14, N15, O16, O17, O18, H2, Li6, Li7, Be9.
Composition evolves by:

- **nuclear burning** — an implicit (backward-Euler) network solve per
  region, Newton iteration on the coupled H/He3/CNO abundances
  (`mixing/kemcom.f90`: up to 50 Newton iterations on a 7-species implicit
  system with an analytic Jacobian, solved by `numerics/simeqc`), simpler
  explicit updates for slow species (`core/burn_lib.f90: dburn/dburnm`);
- **instantaneous convective mixing** — convective zones (plus overshoot
  extensions) are homogenized every step (`mixing/mix.f90`, `mixcz`);
- **light elements** (Li6, Li7, Be9, H2) — burned with sub-stepped
  Bulirsch–Stoer-extrapolated updates that account for the time-dependent
  depth of the surface convection zone (`burn_lib: liburn/liburn2`,
  `deutrate`; the CZ-depth history enters via the pressure-scale-height
  bounds stored in `star%light_burn`);
- **microscopic diffusion** (optional) — gravitational settling and
  thermal diffusion of He and metals in the Bahcall & Loeb / Michaud
  formulation with Thoul, Bahcall & Loeb (1994) coefficients
  (`rotation/microdiff/`, legacy path `grsett`), sub-stepped so no
  settling step exceeds a set fraction of the local settling timescale.

### 1.3 Angular momentum equation (rotation on)

Specific angular momentum j(m) (`star%j_rot`) evolves by:

- **magnetized wind torque** at the surface CZ (Kawaler 1988 power-law
  wind or the Matt-style "PMM" law; saturation at ω_sat; optional disk
  locking on the pre-MS holds ω fixed for `tdisk`),
- **solid-body rotation enforced in convective zones** (`rotation/solid`),
- **diffusive transport in radiative zones**:
  ∂ω/∂t = (1/ρr⁴) ∂/∂r [ρ r⁴ D(r) ∂ω/∂r], with the total diffusion
  coefficient D built from the active hydrodynamic instabilities.

## 2. Numerical scheme for the structure solve

### 2.1 The fitting-point (shooting-splicing) method

YREC does not integrate to the photosphere with the Henyey grid. The
domain is split at a **fitting point** near mass fraction 1 − 10⁻¹²:

- **Interior** (center → fitting point): implicit Henyey relaxation.
- **Envelope** (fitting point → photosphere): direct inward ODE
  integration at fixed (Teff, L).
- **Atmosphere** (τ ≈ 0 → photosphere): direct ODE integration in optical
  depth.

The envelope+atmosphere integrations provide the interior's outer boundary
condition as functions P_fit(Teff, L), T_fit(Teff, L), R_fit(Teff, L).

### 2.2 The envelope triangle

Because each envelope integration is expensive, `core/surfbc.f90`
maintains a **triangle of three (log Teff, log L) vertices** at which full
envelope integrations have been performed, and interpolates the
fitting-point conditions linearly inside it. As the Henyey iteration moves
the surface, the triangle "walks": when the trial (Teff, L) exits the
triangle, the far vertex is reflected through the near edge and one new
envelope integration replaces it. Only when the triangle degenerates is it
rebuilt from scratch (`getnewenv`). This amortizes typically one envelope
integration (or none) per Henyey iteration instead of three.

### 2.3 Linearization and the block elimination

At each interior shell `core/coefft.f90` evaluates the four residuals and
their analytic partial derivatives with respect to the corrections
(δlogP, δlogT, δlogR, δL) at the shell and its neighbor — the EOS
(`eos_get`), opacity (`kap_get`), gradient (`tpgrad`), and energy
generation (`engeb`) derivatives feed the chain rule, e.g. the radiative
branch of equation 3 carries
∂qt/∂logT = qt·(−4 + ∂lnκ/∂lnT + ∂lnκ/∂lnρ·∂lnρ/∂lnT).

The 4-equation, 4-unknown per-shell-pair system is **reduced on the fly**:
`reduce` (called inside `coefft`'s shell loop) builds the trapezoidal
residuals, Gauss-Jordan-eliminates the local 4×4 (pivoting the L, R, T, P
rows in turn) while carrying the coupling to the previous shell, and so
fuses the forward elimination into the coefficient pass itself. What is
stored (`star%elim_coeff(4,2,:)` / `star%elim_rhs(4,:)`) is not a matrix
to be factored later but an already-eliminated **bidiagonal recurrence in
(δlogP, δlogT)**, with rows 3 and 4 retaining the recipe to recover
(δlogR, δL). `core/hsolve.f90` then only closes the recurrence at the
surface — a 2×2 solve against the 6-element linearized surface boundary
condition (`star%surface_bc`) delivered by the envelope triangle — and
back-substitutes inward, recovering δlogR and δL from the stored rows.
Corrections are applied with under-relaxation factors, and convergence is
judged on the maximum residuals/corrections (`star%max_residual`,
`max_correction_index` track the worst shell for diagnostics).

Central boundary conditions come from series expansions at m → 0 (the
central shell uses the expansion values of r, L ∝ m^{1/3}, m rather than
the singular ODE forms).

### 2.4 The four iteration levels

`core/crrect.f90` is called up to four times per timestep by
`evolve_step`, with different boundary-condition policies
(`max_iter_level1..4` control the iteration budget):

1. **Level 1** — iterate with a *fresh* surface BC (envelope triangle
   updated as the solution moves),
2. **Level 2** — BC recomputation forced on,
3. **Level 3** — BC frozen; converges the interior at fixed envelope,
4. **Level 4** — the high-precision pass inside the structure⇄rotation
   loop, repeating after composition/rotation updates (used for solar
   calibration; harmful where semiconvection matters, per the code's own
   comments).

Failure at any level sets the divergence flag; `evolve_step`'s retry loop
(`retry_step`) then calls `starin` to **cut the timestep** and rebuild the
trial model, repeating until convergence or the cut limit.

### 2.5 Envelope/atmosphere integrators

`core/envint_lib.f90` (`atm_get`) integrates:

- **Atmosphere**: from τ near 0 to the photospheric pressure, in the
  variable x = log(κP/g) (one ODE for log P vs. τ via dP/dτ = g/κ), with
  T(τ) supplied by a **gray Eddington relation, the Krishna-Swamy (1966)
  solar-calibrated T(τ)**, or (options 3/4/5) tabulated **model-atmosphere
  boundary conditions** (Kurucz 1990, Castelli–Kurucz, Allard NextGen)
  looked up in (Teff, log g) instead of integrated. A retry loop divides
  the starting density by 10 when the starting point overshoots the
  endpoint.
- **Envelope**: from the photosphere to the fitting mass, three ODEs
  (mass, log T, log R as functions of log P) with full EOS/opacity/MLT
  physics at every substep.

Both integrations use the **Bulirsch–Stoer method** —
`numerics/numerics_lib.f90: bsstep` (modified-midpoint `mmid` sub-stepping
+ `ratext` diagonal rational-function extrapolation to zero step size,
Numerical Recipes RZEXTR), with adaptive step-size memory carried between
calls (one of the six deliberately retained `save` states) and per-step
error control against `y_scale`. Up to 2000 steps; step-halving retry on
error-budget failure.

### 2.6 Timestep control

`util/htimer.f90` selects Δt for the next model from the most restrictive
of a set of governors (all namelist-controlled `time_*` / `tol_*`
factors): fractional changes in L, R, P, T per shell, central and total
hydrogen/helium depletion, CZ-base migration, ω changes (global and in the
mixing region), the H-shell mass scale for shell-burning stages
(`hydrogen_dt`), phase-dependent factors (`timestep_factors` array keyed to
evolutionary phase), and a hard `time_max_dt_frac` growth limit. A
structure-based limiter (`use_structure_dt_limits`) switches itself off
above a central-temperature threshold (with `LPTIME SET FALSE` in the log).
Timestep override (`timestep_override`) forces fixed Δt per run.

### 2.7 Rezoning

`setup/hpoint.f90` re-meshes after each converged model: shells are
inserted/deleted so that a weighted structure function χ (weights
`chi_grid_scale(1:12)` over ΔlogP, ΔlogT, Δlog r, ΔL, composition change,
…) stays within tolerance per shell; extra resolution is forced around the
H-burning shell (three flag points at its base/mid/top), the convective
boundaries, and — with envelope overshoot on — a finely zoned region below
the CZ base. Interpolation of all quantities onto the new grid uses
osculatory (shape-preserving quadratic) splines (`numerics/meval` +
McAllister–Roulier `choose/cases/spline` machinery) and cubic splines
(`kspline/ksplint`).

## 3. Order of operations within one timestep

`core/evolve_step.f90` (the `retry_step` loop; parenthesized routines):

```
 0. On re-entry after divergence: rebuild trial model, cut Δt   (starin)
 1. Mass loss / accretion; re-fit envelope mass                 (massloss, wind/mdot)
 2. Save start-of-step composition -> star%prev%xa_start
 3. COMPOSITION HALF-STEP                                       (mix)
      convec: locate CZs (Schwarzschild), overshoot extensions  (convec, oversh, sconvec)
      homogenize mixed zones                                    (mixcz)
      nuclear burning per radiative region + CZs                (kemcom, dburn/dburnm via burn_lib)
      light elements with CZ-depth history                      (liburn/liburn2, deutrate)
      microscopic diffusion (He, Z settling), sub-stepped       (microdiff/ or grsett)
      renormalize abundances
 4. REZONING                                                    (hpoint)  [skipped during He-flash]
 5. STRUCTURE SOLVE, levels 1-3                                 (crrect -> coefft/reduce/hsolve/surfbc)
      divergence at any level -> cycle retry_step (cut Δt)
 6. Store new entropy terms (rates of change) for next model
 7. STRUCTURE <-> ROTATION LOOP (itrot = 1..num_rotation_structure_iters)
      level-4 structure pass                                    (crrect)
      angular momentum evolution over sub-steps                 (getw)
          wind torque (Kawaler/PMM; disk locking)               (mwind/mcowind)
          solid-body CZs                                        (solid/wczimp)
          instability diffusion + rotational mixing             (seculr -> codiff, rotgrid,
                                                                 dcoeft/tridia or dadcoeft,
                                                                 ndifcom/mixcom, checkj, checkc)
      new fp/ft from new ω                                      (fpft)
 8. Locate H shell and CZ boundaries                            (findsh)
 9. Light-element burning bookkeeping (non-rotating path)       (lirate88/liburn)
10. Output quantities computed into star_info                   (update_output_diagnostics, gettau)
11. Output written (legacy wrtout/.track/.store or MESA-mode
    history/profiles/pulse via yrec_output)                     (output_write_model)
12. Next Δt                                                     (htimer/wtime/ytime)
```

Convergence-of-splitting: because composition moved before structure, the
level-4 loop (`itdif1` iterations) optionally re-does mixing + structure
until self-consistent — the recommended configuration for solar models.

## 4. Microphysics, in dispatch order

### 4.1 Equation of state (`eos/`, facade `eos_get`)

Dispatch inside `eqstat` per thermodynamic regime and namelist choices:

1. **OPAL EOS tables** (1995, 2001, or 2006 vintages; `esac/esac01/
   esac06` interpolation with `rhoofp*` secant inversion from (P,T) to ρ)
   in their validity range (T ≳ 5000 K bands, with a blend/ramp to the
   analytic EOS at the table edges via stored Saha-side derivatives);
2. **SCV** (Saumon–Chabrier–van Horn 1995) tables for cool/dense regimes
   (H/He/Z tables with numerical second derivatives; documented historical
   quirks preserved);
3. **Yale analytic EOS** elsewhere: Saha ionization solve
   (`eos/yale/eqsaha`: H, He, He⁺ plus a set of metals in a
   Newton-iterated 3×3 reduced system, with degenerate-electron pressure
   from tabulated Fermi–Dirac integrals in `eqrelv`, density found by a
   damped Newton iteration in log ρ), optional Debye–Hückel Coulomb
   correction;
4. **fully ionized** branch below/above the Saha range (`eqrelv` alone);
5. **MHD EOS** (Mihalas–Hummer–Däppen) tables optionally replacing the
   whole dispatch (no tables ship with the repo).

Outputs: ρ, ∂lnρ/∂lnT, ∂lnρ/∂lnP, c_P, ∇_ad, Γ₁, ionization fractions,
degeneracy η, plus the T- and P-derivatives of these when the solver needs
Jacobians (`want_derivatives`).

### 4.2 Opacity (`kap/`, facade `kap_get`)

Rosseland-mean κ blended from:

- **conductive opacities** (Potekhin tables, `condopacp`, added in
  parallel: 1/κ = 1/κ_rad + 1/κ_cond);
- **radiative tables**, chosen by regime and namelist: OPAL 1995 (GS98
  mixture as shipped; 4-D interpolation in (Z, X, T, log R ≡ ρ/T₆³),
  including the X-, Z-direction table stacking and edge extrapolation
  logic in `getopal95`), OPAL 1992, LAOL 1989; **low-temperature
  molecular opacities** from Ferguson/Alexander 2006 or Alexander 1995,
  or Kurucz 1990, blended across a T_mol transition window;
- optional pure-Z tables supporting metal diffusion.

Derivatives ∂lnκ/∂lnρ and ∂lnκ/∂lnT are returned from the interpolants
(quadratic/cubic interpolation with stored weights per table family).

### 4.3 Convection (`mixing/tpgrad.f90`)

Schwarzschild criterion (∇_rad vs ∇_ad; Ledoux/semiconvection is a
separate flag treatment, mutually exclusive with core overshoot).
Convective gradient from standard **Böhm-Vitense (1958) mixing-length
theory**: the cubic MLT equation is solved for the convective efficiency,
giving ∇_conv between ∇_ad and ∇_rad, with α = `mixing_length_alpha`
(per-run; the solar-calibrated value is the free parameter of the solar
calibration). Analytic ∂∇/∂(lnP, lnT, ln r) derivatives feed the Henyey
Jacobian. Starspot-modified variants (Somers: spot filling factor/flux
contrast `spotf/spotx`) adjust the effective ∇_rad. **Overshoot**: the
mixed region is extended a fixed number of pressure scale heights
(`overshoot_alpha_core/envelope`, `oversh.f90`) — mixing only ("step"
overshoot), optionally with an adiabatic temperature gradient in the
overshoot region (`adiabatic_overshoot_active`); `sconvec` implements the
rate-limited "semiconvective" overshoot variant with a propagation-speed
cap.

### 4.4 Nuclear network (`net/net_lib.f90` kernels + `core/burn_lib.f90` drivers)

- 13 principal reactions: pp, He3+He3, He3+He4, p-captures on C12, C13,
  N14, O16, Be7(e⁻,ν), Be7(p,γ), plus He-burning (3α, C12+α, C13+α,
  N14+α) — rates in the Caughlan & Fowler (1988) / Bahcall parameterized
  forms with the N15 branching computed explicitly; an optional
  newer-rate set with **Solar Fusion II** S-factors and per-reaction
  namelist scale factors (`s0_*` controls, `setup/remap.f90`);
- **weak and intermediate screening** (Salpeter/Graboske) selected by a
  threshold; the code's own comment notes strong screening is untested;
- **neutrino energy losses** from each pp/CNO chain branch, plus
  plasma/photo/pair/bremsstrahlung losses from **Itoh et al. (1996)**
  (`sneut`) for advanced stages; the solar-neutrino diagnostic computes
  Cl37 and Ga71 SNU rates from the flux spectrum;
- energy generation `engeb` assembles ε per shell and its T, ρ
  derivatives for the Henyey Jacobian, and stores the per-reaction rate
  arrays used by the composition solve.

### 4.5 Atmospheres (`atm/`)

T(τ) relations (gray Eddington; Krishna-Swamy) integrated by the
atmosphere ODE, or tabulated surface boundary conditions
P_surf(Teff, log g) from Kurucz 1990, Castelli–Kurucz, or Allard NextGen
grids (with τ=100 P, T columns for the pulse output path), cubic-spline
interpolated in both directions with edge-validity bookkeeping
(gmin/gmax per Teff row). Out-of-range Allard lookups fall back to gray.

### 4.6 Rotation physics (`rotation/`)

- **Structural**: fp/ft factors and shell oblateness (η₂ from the Radau
  equation integrated over the model; `shape/`), surface-average gravity
  corrections.
- **Instability diffusion coefficients** (`seculr/codiff`, from the
  Pinsonneault et al. 1989 formulation): dynamical instabilities
  (checked/enforced instantaneously) and secular ones providing D(r) —
  **Eddington–Sweet meridional circulation**, **Goldreich–Schubert–
  Fricke**, **secular shear**, with mean-molecular-weight-gradient
  inhibition (f_μ weighting) and per-mechanism efficiency factors
  (`fes/fgsf/fss/fc` — the f_c "coupling" factor relates composition
  mixing to angular momentum transport, D_mix = f_c · D_AM);
  optionally a constant additional D.
- **Solvers**: each unstable radiative region is mapped to an
  equally-spaced grid (`rotgrid`), the AM diffusion equation is solved
  implicitly — tridiagonal (`dcoeft`+`tridia`) for the pure-diffusion
  form, or a banded advective-diffusive variant (`dadcoeft`, with its own
  timestep-cutting loop) — then mapped back (`eq2mod`). `checkj` guards
  against unphysical ω-gradient reversals (enforcing local solid-body
  where diffusion overshoots); composition mixing uses the same grid
  (`ndifcom/mixcom`, species-by-species with conservation-checked
  updates); `checkc` validates and can force timestep cuts.
- **Sub-stepping**: AM evolution runs on sub-timesteps (`dtdif` fraction
  of the diffusion time), with the wind torque applied per sub-step.

## 5. The numerics toolbox (`numerics/numerics_lib.f90`)

| Kernel | Method | Used by |
|---|---|---|
| `bsstep`/`mmid`/`ratext` | Bulirsch–Stoer with modified midpoint + rational extrapolation | atmosphere & envelope integration, Li/Be burning extrapolation |
| `hsolve` (in core/) | Henyey surface closure + back-substitution (elimination lives in `reduce`) | structure solve |
| `simeqc` | Gauss-Jordan with pivoting | kemcom's Newton systems, calibration |
| `tridia`/`ctridi` | tridiagonal (Thomas) solvers | AM diffusion, splines |
| `kspline`/`ksplint`, `spline`(osculatory), `meval`/`choose`/`cases` | cubic and shape-preserving quadratic splines | rezoning, table interpolation, turnover timescale |
| `lir`/`lir1` | 4-point Lagrange / linear table interpolation with index memory | EOS/opacity table lookups |
| `hunt`, `locate`, `search`, `findex` | bisection / hunting searches | all table families |
| `qgauss` | 5-point Gauss–Legendre quadrature | shape integrals |
| `intpol` | spline interpolation with derivative | opacity (Cox) path |

Six of these deliberately retain internal `save` state (algorithm memory:
step-size history, refinement accumulators, search-index memory) — the
only remaining `save`s in the code outside the three driver files.

## 6. Calibration modes

- **Solar calibration** (`calibrate_solar_model`, kind-card triples): runs
  the model to the solar age, compares (log L, log R) — and optionally
  Z/X — to solar targets, and updates (X₀, α_MLT) (and Z₀ for Z/X) with
  stored empirical partial derivatives (∂logL/∂X = −3.78 etc., hardwired
  from experiment in `chkcal`), iterating whole evolutionary runs to
  tolerance (`toll/tolr` ≈ 3×10⁻⁵). Yields the classic calibrated solar
  model in a handful of runs.
- **Star calibration** (`lcalst`/`chkscal`): same machinery against
  arbitrary (Teff, L) targets.
- **Monte-Carlo mode** (`lmonte`): repeated runs drawing physics
  parameters (rate S-factors, diffusion factors, Z/X…) per realization —
  the classic solar-model uncertainty machinery (requires legacy output).

## 7. Program architecture (2026)

```
main (CLI, 10 lines)                      pyyrec (ctypes) ── libyrec
  └─ run_yrec(ierr)                            └─ yrec_capi: yrec_run()
       ├─ yrec_run_prologue     (re-entrancy: snapshot/restore all state)
       ├─ read_controls         (io/parmin: legacy $CONTROL/$PHYSICS or
       │                         MESA-style &star_job/&controls inlists)
       ├─ star_setup            (eos_init/kap_init/atm_init table loads)
       └─ run loop (kind cards) ── model loop ── evolve_step
                                                    │
     star layer (owns type(star_info) star): starin, crrect, hsolve,
       coefft, envint_lib, surfbc, burn_lib, mix, hpoint, getw, wind
                                                    │
     physics services (never see star_info): eos_get, kap_get, net_lib
       kernels, atm tables, numerics_lib
                                                    │
     io: yrec_output (legacy .short/.track/.store OR MESA-style
       history.data / profile{N}.data / pulse files), yrec_reset
```

All mutable state hangs off one root — `star` (`state/star_info_lib.f90`),
with `star%prev` (start-of-step model), `star%job` (run configuration),
`star%evo` (driver scalars), nested MESA-style; physics-domain table state
lives inside each domain. Two `run_yrec` calls in one process are
byte-identical to two fresh processes (the re-entrancy contract), which is
what makes the Python embedding sound.

## 8. Key lineage references

- Henyey, Forbes & Gould 1964 (relaxation method); Larson & Demarque 1964
  (the Yale code ancestry).
- Pinsonneault 1988 (PhD thesis: the rotating code); Pinsonneault et al.
  1989, ApJ 338, 424 (rotational instabilities & solar spin-down).
- Kippenhahn & Thomas 1970; Endal & Sofia 1976, 1978 (fp/ft rotating
  structure).
- Böhm-Vitense 1958 (MLT); Krishna-Swamy 1966 (T–τ).
- Guenther, Demarque, Kim & Pinsonneault 1992, ApJ 387, 372 (standard
  solar model with YREC).
- Bahcall & Loeb 1990; Thoul, Bahcall & Loeb 1994 (diffusion).
- Kawaler 1988 (wind law); Matt et al. 2012-era laws (the PMM option);
  Reimers 1975 (red-giant mass loss).
- Rogers & Iglesias / Iglesias & Rogers (OPAL opacities & EOS); Ferguson
  et al. 2005 / Alexander & Ferguson 1994 (low-T opacities); Saumon,
  Chabrier & van Horn 1995 (SCV EOS); Potekhin (conduction).
- Caughlan & Fowler 1988; Adelberger et al. 2011 (Solar Fusion II);
  Itoh et al. 1996 (neutrino losses); Grevesse & Sauval 1998 (mixture).
- van Saders & Pinsonneault 2012; Somers & Pinsonneault 2015-16 (spots,
  modern YREC applications).

---

*Caveats: hsolve's algebra is dense and was renamed structurally, not
rederived, during modernization — the description in §2.3 reflects the
code's own documentation plus the coefft equation forms. The MHD EOS path
is unexercised (no tables ship). The Monte-Carlo and helium-flash paths
are legacy-output-only.*

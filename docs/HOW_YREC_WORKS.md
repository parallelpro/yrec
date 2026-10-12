# How YREC Works

*A technical report on the equations, numerical schemes, and physics of the
Yale Rotating Evolution Code — written against the 2026 modernized source
(branch `yrec-modern`; file, routine and field references refreshed
against `src/` on 2026-10-10). File references use the current layout,
relative to `src/`. Routines are given by their current name, with the
legacy F77 name in parentheses at the first mention in each section; the
full map is Appendix A of `TARGET_ARCHITECTURE.md` (this directory). The
physics and numerics are unchanged from the classic code
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
4. picks the next timestep and writes diagnostics.

This is the classic operator-split Henyey approach: composition and
structure are advanced in alternating half-steps rather than solved as one
coupled system, with a fourth iteration level available to converge the
splitting for high-precision (solar) work.

### 1.1 The structure equations

With mass `m` as the independent (Lagrangian) coordinate, YREC solves the
four standard stellar-structure ODEs, in log₁₀ form internally
(`logP`, `logT`, `logR`, and linear `L` as dependent variables, `log m` as
the coordinate; all cgs):

| # | Equation | YREC form (per shell, `henyey_coefficients` (`coefft`), `core/henyey_coefficients.f90`) |
|---|----------|------------------------------------------|
| 1 | Mass conservation | dr/dm = 1/(4πr²ρ) — residual `cur%qr = exp(ln10·(log m − log ρ − log 4πr³))` |
| 2 | Hydrostatic equilibrium | dP/dm = −Gm/(4πr⁴) · **f_P** — residual `qp = −exp(ln10·(log G + 2 log m − log P − log 4πr³ − log r))·f_P` (held in `pt_scr%qp`, mirrored into `cur%qp`) |
| 3 | Energy transport | dlnT/dm = ∇ · dlnP/dm — residual `qt = ∇ · qp`, with ∇ = ∇_rad (radiative) or ∇_conv (MLT) |
| 4 | Energy conservation | dL/dm = ε_nuc − ε_ν + ε_grav — assembled from `compute_zone_energy` → `nuclear_energy` (burning + neutrinos) plus the entropy ("gravothermal") term −T ∂S/∂t, evaluated from the change in (log T, log P) since the start of the step divided by Δt |

**Rotation enters the structure equations** through the Kippenhahn &
Thomas (1970) / Endal & Sofia (1976) form factors `f_P` (`star%fp_rot`) and
`f_T` (`star%ft_rot`): surface-averaged corrections computed by
`rotation_shape_factors` (`fpft`,
`rotation/shape/rotation_shape_factors.f90`) from ω(m) and the shell shape
(oblateness via the η₂ Radau parameter, `star%eta_squared`), multiplying
the hydrostatic and radiative-transport terms respectively. This is the
standard 1.5-D scheme: the star remains a 1-D object, but each shell is an
oblate equipotential.

The **entropy term** deserves note: rather than differencing S directly,
YREC evaluates it from the per-shell changes of log T and log P since the
start of the step (`star%log_temperature_delta`, `star%log_pressure_delta`,
accumulated with every Newton correction) and the local EOS derivatives.
Separately, after the second iteration level it stores per-shell rates of
change of (log T, log P, L, log R) (`star%temperature_entropy_term`,
`star%pressure_entropy_term`, `star%luminosity_entropy_term`,
`star%radius_entropy_term`) and uses them to extrapolate the first guess of
the next model (`star%ctrl%improved_first_guess_flag`, namelist `lnews`); the
extrapolated change also seeds the log T / log P accumulators.

### 1.2 Composition equations

15 species are tracked per shell (`star%xa(1:15,k)`): H1, He4, Z (bulk
metals), He3, C12, C13, N14, N15, O16, O17, O18, H2, Li6, Li7, Be9.
Composition evolves by:

- **nuclear burning** — an implicit (backward-Euler) network solve per
  region, Newton iteration on the coupled H/He3/CNO abundances
  (`burn_implicitly`, `net/network.f90`, driven per region by
  `burn_network_region` in `core/burn_step.f90`: the equations and
  Jacobian come from a reaction table, solved by `simeqc` in
  `numerics/numerics_lib.f90`), plus explicit updates: equilibrium-He3/CN
  hydrogen and helium burning rates on the first iteration level
  (`compute_equilibrium_burn_rates`) and sub-stepped deuterium burning
  (`burn_deuterium_region`, `net/light_elements.f90`). `src/net/README.md`
  walks through these routines and their call order;
- **instantaneous convective mixing** — convective zones (plus overshoot
  extensions) are homogenized every step (inside `advance_composition`,
  `core/composition_step.f90`, and again after the structure converges by
  `homogenize_convection_zones` (`mixcz`));
- **light elements** (Li6, Li7, Be9, H2) — burned in the trace
  approximation by exponential decay with the rate interpolated over the
  step (closed form, `surviving_fraction` in `net/light_elements.f90`),
  accounting for the time-dependent depth of the surface convection zone
  (`burn_li_be` in `core/burn_step.f90`, rates from `compute_li_be_rates`, the CZ
  base from `core/convection_geometry.f90`; the CZ-depth history enters via the pressure-scale-height
  bounds `star%pressure_scale_height_start`/`_end` and the previous CZ-base
  radius `star%cz_base_radius_prev`);
- **microscopic diffusion** (optional) — gravitational settling and
  thermal diffusion of He and metals in the Bahcall & Loeb / Michaud
  formulation with Thoul, Bahcall & Loeb (1994) coefficients
  (`rotation/microdiff/`: `microdiff`, or the older
  `gravitational_settling` (`grsett`)), sub-stepped so no
  settling step exceeds a set fraction of the local settling timescale.

### 1.3 Angular momentum equation (rotation on)

Specific angular momentum j(m) (`star%j_rot`) evolves by:

- **magnetized wind torque** at the surface CZ (Kawaler 1988 power-law
  wind or the Matt-style "PMM" law; saturation at ω_sat; optional disk
  locking on the pre-MS holds ω fixed for `tdisk`),
- **solid-body rotation enforced in convective zones** (`solid_body_omega`
  (`solid`), `rotation/solid_body_omega.f90`),
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

Because each envelope integration is expensive, `surfbc`
(`core/surfbc.f90`) maintains a **triangle of three (log Teff, log L)
vertices** at which full envelope integrations have been performed, and
interpolates the fitting-point conditions linearly inside it
(`star%envelope_fit_coeffs`). As the Henyey solution moves the surface,
the triangle "walks": when the trial (Teff, L) exits the triangle, the far
vertex is reflected through the near edge and one new envelope
integration replaces it. All three vertices are recomputed only when a new
triangle is requested: the first model of a run, the `lnew0` control, a
switch of atmosphere table, or a moved fitting point (the outermost
interior point added or dropped by `rezone` (`hpoint`), or a mass change
from mass loss/accretion). This amortizes typically one envelope
integration (or none) per boundary-condition refresh (at most one per
call to `henyey_iterate` (`crrect`)) instead of three.

### 2.3 Linearization and the block elimination

At each interior shell `henyey_coefficients` (`coefft`,
`core/henyey_coefficients.f90`) evaluates the four residuals and
their analytic partial derivatives with respect to the corrections
(δlogP, δlogT, δlogR, δL) at the shell and its neighbor — the EOS
(`eos_get`), opacity (`kap_get`), gradient (`compute_temperature_gradients`
(`tpgrad`)), and energy generation (`compute_zone_energy`) derivatives feed the chain
rule, e.g. the radiative
branch of equation 3 carries
∂qt/∂logT = qt·(−4 + ∂lnκ/∂lnT + ∂lnκ/∂lnρ·∂lnρ/∂lnT).

The 4-equation, 4-unknown per-shell-pair system is **reduced on the fly**:
`henyey_eliminate` (`reduce`, `core/henyey_elimination.f90`; called inside
`henyey_coefficients`' shell loop) builds the trapezoidal
residuals, Gauss-Jordan-eliminates the local 4×4 (pivoting the L, R, T, P
rows in turn) while carrying the coupling to the previous shell, and so
fuses the forward elimination into the coefficient pass itself. What is
stored (`star%elim_coeff(4,2,:)` / `star%elim_rhs(4,:)`) is not a matrix
to be factored later but an already-eliminated **bidiagonal recurrence in
(δlogP, δlogT)**, with rows 3 and 4 retaining the recipe to recover
(δlogR, δL). `henyey_solve` (`hsolve`, `core/henyey_solve.f90`) then only
closes the recurrence at the surface — a 2×2 solve against the 6-element
linearized surface boundary condition (`star%surface_bc`, built by
`henyey_iterate` from the envelope triangle's `star%envelope_fit_coeffs`)
— and back-substitutes inward, recovering δlogR and δL from the stored rows.
Corrections are applied with under-relaxation factors, and convergence is
judged on the maximum residuals/corrections (`star%max_residual`,
`star%max_correction_index` track the worst shell for diagnostics).

Central boundary conditions come from series expansions at m → 0 (the
central shell uses the expansion values of r, L ∝ m^{1/3}, m rather than
the singular ODE forms).

### 2.4 The four iteration levels

`henyey_iterate` (`core/henyey_iteration.f90`) is called up to four times
per timestep by `evolve_step` (through its `solve_level` wrapper), with
different boundary-condition policies (`star%ctrl%max_iter_level1..4`
control the iteration budget). The surface BC is refreshed (`surfbc`) at
most once per call, before the Newton loop:

1. **Level 1** — iterate on the *cached* envelope fit of the previous
   model (the triangle is built here only when a new one is requested,
   e.g. for the first model of a run),
2. **Level 2** — BC recomputation forced on: the triangle is re-checked
   (and walked) at the level-1 solution before iterating,
3. **Level 3** — BC frozen; the composition is re-mixed from the
   start-of-step abundances at the current structure (`advance_composition`,
   the corrector pass, called from `solve_level`), then the interior is converged at fixed envelope,
4. **Level 4** — the high-precision pass inside the structure⇄rotation
   loop, also re-mixing first, repeating after composition/rotation
   updates (used for solar calibration; harmful where semiconvection
   matters, per the code's own comments).

Corrections beyond the divergence tolerances (`star%ctrl%htoler(:,2)`) at
any level, or a model still unconverged after level 4, set the divergence
flag (`star%model_diverged_flag`). `evolve_step`'s retry loop (`retry_step`)
then calls `read_starting_model` (`starin`), which **stops the run** with
an error ("MODEL FAILED TO CONVERGE", `ierr = 1`). Despite the loop's
shape and that routine's header comment, the current code does not cut
the timestep and retry.

### 2.5 Envelope/atmosphere integrators

`atm_get` (`envint`; driver in `core/envelope_integration.f90`, star-blind kernel in
`core/envelope_integration_kernel.f90`, integrands in `core/atmosphere_derivs.f90` and
`core/envelope_derivs.f90`) integrates:

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

`compute_timestep` (`htimer`, `util/compute_timestep.f90`, with the
per-criterion limiters `timestep_limit_structure`/`_omega`/`_hburn`/
`_heburn`/`_hr` (`ptime`/`wtime`/`xtime`/`ytime`/`entime`) in
`util/timestep_limit_*.f90`) selects Δt for the next model from the most
restrictive of a set of governors (all namelist-controlled `time_*`
factors): fractional changes in L, R, P, T per shell, central and total
hydrogen/helium depletion, CZ-base migration, ω changes (global and in the
mixing region), the H-shell mass scale for shell-burning stages
(`hydrogen_dt`), and a hard `time_max_dt_frac` growth limit. The `time_*`
controls are the named slots of the 14-element per-phase factor array
`star%ctrl%atime` (inlist name `timestep_factors`). A structure-based
limiter (`star%job%use_structure_dt_limits`, legacy `LPTIME`) switches
itself off once log T_c > 7.1 with L_grav < 0 (the terminal reports
"structure-based limits disabled"). Timestep override
(`star%job%timestep_override`) forces fixed Δt per run.

### 2.7 Rezoning

`rezone` (`hpoint`, `setup/rezone.f90`) re-meshes after each converged
model: shells are inserted/deleted so that a structure function
χ = Δlog m/δm + ΔL/δL − ΔlogP/δP, with separate ΔX and ΔZ limits, stays
within tolerance per shell (the scales come from
`star%ctrl%chi_grid_scale(1:12)`, whose named slots are the `tol_*`/`flag_*`
controls; the log P scale differs between the envelope, the CZ base and
the radiative interior); extra resolution is forced around the
H-burning shell (three flag points at its base/mid/top), the convective
boundaries, and — with envelope overshoot on — a finely zoned region below
the CZ base. Interpolation of all quantities onto the new grid uses
osculatory (shape-preserving quadratic) splines (`osplin` → `meval` +
McAllister–Roulier `choose/cases/spline` machinery,
`numerics/quadratic_spline.f90`) and natural cubic splines
(`natural_spline`/`splint`, `numerics/cubic_spline.f90`).

## 3. Order of operations within one timestep

`evolve_step` (`core/evolve_step.f90`; stages 0-4 run inside the
`retry_step` loop, then 5-12 once; current names, legacy in parentheses):

```
 0. reload_model_if_diverged: if the last attempt
      diverged, stop the run (ierr = 1); no Δt cut   read_starting_model (starin)
 1. advance_composition_and_age (model being aged):
      mass loss / accretion (flags step 7)           change_mass (massloss, mdot)
      save start-of-step composition -> star%xa_start
      COMPOSITION HALF-STEP (predictor pass)         advance_composition (mix)
        locate CZs (Schwarzschild), overshoot        find_convection_zones (convec)
                                                       -> extend_by_overshoot (oversh)
        rates of every shell                         compute_burn_rates
        implicit burn per radiative shell            burn_network -> burn_network_region
          and per mixed zone (homogenized as unit)     -> burn_implicitly (kemcom)
        explicit equilibrium burn, level 1 only      burn_equilibrium (eqburn)
        semiconvection (optional, levels > 1)        semiconvection (sconvec)
        microscopic diffusion, sub-stepped           microdiff or gravitational_settling
                                                       (grsett)
        renormalize abundances
        deuterium burning                            burn_deuterium (dburn)
      advance the age
 2. rezone_or_snapshot: REZONING [not in He flash]   rezone (hpoint)
      start-of-step Li/Be rates                      store_start_li_be_rates (lirate88)
 3. solve_structure: STRUCTURE SOLVE, levels 1-3     solve_level, per level:
      surface boundary condition                     refresh_surface_bc -> surfbc
      levels >= 3: composition step again from
        star%xa_start (corrector pass)               advance_composition
      Newton loop                                    henyey_iterate (crrect) ->
                                                       henyey_coefficients (coefft),
                                                       henyey_eliminate (reduce),
                                                       henyey_solve (hsolve)
      store new entropy terms (rates) after level 2
      divergence at any level -> cycle retry_step (stage 0)
 4. converge_with_rotation (itrot = 1..num_rotation_structure_iters; 1 if no rotation)
      level-4 structure pass                         solve_level (as in 3)
      re-homogenize convection zones                 homogenize_convection_zones (mixcz)
      AM evolution over sub-steps (rotation on)      evolve_angular_momentum (getw)
        wind torque (Kawaler/PMM; disk locking)      brake_convection_zone,
                                                       compute_solid_body_braking
                                                       (mwind, mcowind)
        rotation law / solid body in CZs             enforce_rotation_profile (wczimp),
                                                       solid_body_omega (solid)
        instability diffusion + rotational mixing    apply_secular_transport (seculr) and the
                                                       rotation/seculr/ solvers (see 4.6)
        Li/Be burning (with rotational mixing)       store_end_li_be_rates, burn_li_be
      new fp/ft from new ω                           rotation_shape_factors (fpft)
 5. Locate H shell and CZ boundaries                 locate_shell_boundaries (findsh)
 6. Light-element burning (non-rotating path)        burn_light_elements: store_end_li_be_rates,
                                                       burn_li_be (liburn)
 7. After a mass change: move the fitting point      rebuild_envelope (getnewenv),
                                                       rotation_shape_factors
 8. Next Δt                                          compute_timestep (htimer)
                                                       -> timestep_limit_*
 9. Assemble the stitched full-star model            build_stitched_model
10. Observables computed into star_info              compute_observables
                                                       (update_output_diagnostics), incl.
                                                       compute_turnover_timescale (gettau)
11. Output: history / profiles / pulse files         output_write_model (yrec_output;
                                                       legacy wrtout retired)
12. Stop conditions (age, abundances, radius)        check_stop_conditions
```

Convergence-of-splitting: because composition moved before structure,
levels 3 and 4 re-do the mixing from the start-of-step composition at the
current structure before iterating (`advance_composition`, the corrector pass, run by `solve_level` between the surface boundary condition and `henyey_iterate`), and the
level-4 loop runs `num_rotation_structure_iters` (`itdif1`) passes — the
recommended configuration for solar models.

## 4. Microphysics, in dispatch order

### 4.1 Equation of state (`eos/`, facade `eos_get`)

`eos_get` (`eos/eos_lib.f90`) calls `eqstat`, which dispatches (with
`eqstat2`, both in `eos/eqstat.f90`) per thermodynamic regime and namelist
choices:

1. **OPAL EOS tables** (1995, 2001, or 2006 vintages; `esac/esac01/
   esac06` interpolation with `rhoofp*` secant inversion from (P,T) to ρ)
   in their validity range (T ≳ 5000 K bands, with a blend/ramp to the
   analytic EOS at the table edges via stored Saha-side derivatives);
2. **SCV** (Saumon–Chabrier–van Horn 1995) tables for cool/dense regimes
   (H/He/Z tables with numerical second derivatives; documented historical
   quirks preserved);
3. **Yale analytic EOS** elsewhere: Saha ionization solve
   (`saha_eos` (`eqsaha`), `eos/yale/saha_eos.f90`: H, He, He⁺ plus a
   set of metals in a Newton-iterated 3×3 reduced system, with
   degenerate-electron pressure from tabulated Fermi–Dirac integrals in
   `fully_ionized_eos` (`eqrelv`, `eos/yale/fully_ionized_eos.f90`),
   density found by a damped Newton iteration in log ρ), optional
   Debye–Hückel Coulomb correction;
4. **fully ionized** branch below/above the Saha range (`fully_ionized_eos`
   alone);
5. **MHD EOS** (Mihalas–Hummer–Däppen) tables optionally replacing the
   whole dispatch (`meqos`, `eos/mhd/meqos.f90`, called directly by
   `eos_get`; no tables ship with the repo).

Outputs: ρ, ∂lnρ/∂lnT, ∂lnρ/∂lnP, c_P, ∇_ad, Γ₁, ionization fractions,
degeneracy η, plus the T- and P-derivatives of these when the solver needs
Jacobians (`want_derivatives`).

### 4.2 Opacity (`kap/`, facade `kap_get` over `kap_eval` (`getopac`), `kap/kap_lib.f90`)

Rosseland-mean κ blended from:

- **conductive opacities** (Potekhin tables,
  `kap/conductive/conductive_opacity.f90`, added in parallel:
  1/κ = 1/κ_rad + 1/κ_cond);
- **radiative tables**, chosen by regime and namelist: OPAL 1995 (GS98
  mixture as shipped; 4-D interpolation in (Z, X, T, log R ≡ ρ/T₆³),
  including the X-, Z-direction table stacking and edge extrapolation
  logic in `kap/opal95/opal95_opacity.f90`), OPAL 1992, LAOL 1989;
  **low-temperature molecular opacities** from Ferguson/Alexander 2006 or Alexander 1995,
  or Kurucz 1990, blended across a T_mol transition window;
- optional pure-Z tables supporting metal diffusion.

Derivatives ∂lnκ/∂lnρ and ∂lnκ/∂lnT are returned from the interpolants
(quadratic/cubic interpolation with stored weights per table family).

### 4.3 Convection (`convection/`, through `convection_lib`; `compute_temperature_gradients` (`tpgrad`))

Schwarzschild criterion (∇_rad vs ∇_ad; Ledoux/semiconvection is a
separate flag treatment, mutually exclusive with core overshoot).
Convective gradient from standard **Böhm-Vitense (1958) mixing-length
theory**: the cubic MLT equation is solved for the convective efficiency,
giving ∇_conv between ∇_ad and ∇_rad, with α = `star%mixing_length_alpha`
(per-run; the solar-calibrated value is the free parameter of the solar
calibration). Analytic ∂∇/∂(lnP, lnT, ln r) derivatives feed the Henyey
Jacobian. Starspot-modified variants (Somers: spot filling factor/flux
contrast `star%ctrl%spot_filling_factor`/`spot_temp_contrast`, namelist
`spotf/spotx`) adjust the effective ∇_rad. **Overshoot**: the
mixed region is extended a fixed number of pressure scale heights
(`star%ctrl%overshoot_alpha_core`/`_envelope`; `extend_by_overshoot`
(`oversh`), `convection/convection_zones.f90`, called from
`find_convection_zones` (`convec`)) — mixing only ("step"
overshoot), optionally with an adiabatic temperature gradient in the
overshoot region (`star%ctrl%ladov`, inlist name
`adiabatic_overshoot_active`); `extend_by_semiconvection` (`sconvec`,
`convection/semiconvection.f90`) implements the
rate-limited "semiconvective" overshoot variant with a propagation-speed
cap.

### 4.4 Nuclear network (`net/` physics + `core/burn_step.f90` drivers)

The routines, their call order and how to extend them are described in
`src/net/README.md`. In brief:

- 13 reaction slots: pp, He3+He3, He3+He4, p-captures on C12, C13,
  N14 and O16, plus He burning (3α and α-captures on C12, C13 and N14);
  two slots (O16+α, C12+C12) are reserved and never computed, and the
  Be7(e⁻,ν) and Be7(p,γ) branches of He3+He4 enter as a branching
  fraction — rates in the Caughlan & Fowler (1988) / Bahcall parameterized
  forms with the N15 branching computed explicitly; an optional
  newer-rate set with **Solar Fusion II** S-factors and per-reaction
  namelist scale factors (`s0_*` controls, turned into
  `star%cross_section_scale` by `map_user_inputs` (`remap`) in
  `setup/map_user_inputs.f90`);
- **weak and intermediate screening** (Salpeter/Graboske) selected by a
  threshold; the code's own comment notes strong screening is untested;
- **neutrino energy losses** from each pp/CNO chain branch, plus
  plasma/photo/pair/bremsstrahlung losses from **Itoh et al. (1996)**
  (`sneut`, `net/thermal_neutrinos.f90`) for advanced stages; the solar-neutrino
  diagnostic computes Cl37 and Ga71 SNU rates from the flux spectrum;
- one rate evaluation, `compute_reaction_rates`, feeds both consumers:
  `nuclear_energy` (`net/nuclear_energy.f90`) assembles ε per shell and
  its T, ρ derivatives for the Henyey Jacobian, and `compute_zone_energy`
  returns everything it computes in its result record; the composition
  solve takes the same rates repacked per s per amu by
  `compute_network_rates` (called from `compute_burn_rates` in `advance_composition`).

### 4.5 Atmospheres (`atm/`)

T(τ) relations (gray Eddington; Krishna-Swamy) integrated by the
atmosphere ODE, or tabulated surface boundary conditions
P_surf(Teff, log g) from Kurucz 1990, Castelli–Kurucz, or Allard NextGen
grids (with τ=100 P, T columns for the pulse output path), cubic-spline
interpolated in both directions with edge-validity bookkeeping
(gmin/gmax per Teff row). Out-of-range Allard lookups fall back to gray.

### 4.6 Rotation physics (`rotation/`)

- **Structural**: fp/ft factors and shell oblateness (η₂ from the Radau
  equation integrated over the model; `rotation/shape/`, entry point
  `rotation_shape_factors` (`fpft`)), surface-average gravity
  corrections.
- **Instability diffusion coefficients** (`apply_secular_transport` (`seculr`)
  and `diffusion_velocity_scales` (`codiff`) in `rotation/seculr/`, from the
  Pinsonneault et al. 1989 formulation): dynamical instabilities
  (checked/enforced instantaneously) and secular ones providing D(r) —
  **Eddington–Sweet meridional circulation**, **Goldreich–Schubert–
  Fricke**, **secular shear**, with mean-molecular-weight-gradient
  inhibition (f_μ weighting) and per-mechanism efficiency factors
  (`fes/fgsf/fss/fc` — the f_c "coupling" factor relates composition
  mixing to angular momentum transport, D_mix = f_c · D_AM);
  optionally a constant additional D.
- **Solvers**: each unstable radiative region is mapped to an
  equally-spaced grid (`am_transport_grid` (`rotgrid`)), the AM diffusion
  equation is solved implicitly — tridiagonal (`am_diffusion_coeffs` (`dcoeft`) + `solve_angular_momentum_transport` (`tridia`))
  for the pure-diffusion form, or a banded advective-diffusive variant
  (`am_advection_diffusion_coeffs` (`dadcoeft`) + `banded_solver`, with its
  own timestep-cutting loop) — then mapped back (`equal_grid_to_model`
  (`eq2mod`)). `check_angular_momentum` (`checkj`) guards
  against unphysical ω-gradient reversals (enforcing local solid-body
  where diffusion overshoots); composition mixing uses the same grid
  (`diffuse_composition_driver`/`diffuse_composition` (`ndifcom`/`mixcom`),
  species-by-species with conservation-checked updates);
  `check_composition` (`checkc`) validates and can force timestep cuts.
- **Sub-stepping**: AM evolution (`evolve_angular_momentum` (`getw`)) runs
  on sub-timesteps (`dtdif` fraction of the diffusion time), with the wind
  torque (`rotation/magnetic_braking.f90`: `brake_convection_zone`/
  `compute_solid_body_braking` (`mwind`/`mcowind`)) applied per sub-step.

## 5. The numerics toolbox (`numerics/`, through `numerics_lib`)

| Kernel | Method | Used by |
|---|---|---|
| `bsstep` (with `mmid`, `ratext`) | Bulirsch–Stoer with modified midpoint + rational extrapolation, for an `ode_system` | atmosphere & envelope integration (`envelope_equations`) |
| `henyey_solve` (`hsolve`; in `core/henyey_solve.f90`) | Henyey surface closure + back-substitution (elimination lives in `henyey_eliminate` (`reduce`)) | structure solve |
| `simeqc` | Gauss-Jordan with pivoting | `burn_implicitly`'s Newton systems |
| `ludcmp`/`lubksb` | LU with partial pivoting | Thoul diffusion coefficients |
| `solve_tridiagonal` (`tridiag_gs`, `ctridi`) | tridiagonal (Thomas) solver | gravitational settling, rotational composition diffusion |
| `solve_angular_momentum_transport` (`tridia`; in `rotation/seculr/secular_transport.f90`) | Thomas solver with the angular-momentum bookkeeping | AM diffusion |
| `cspline`/`splint`, `natural_spline` (`kspline`, `splinc`), `osplin`→`meval`/`choose`/`cases`/`spline` | cubic and shape-preserving quadratic splines | rezoning, table interpolation, turnover timescale |
| `lir` | 4-point Lagrange / linear table interpolation (`lir1` was removed as dead code) | MHD EOS table lookups |
| `hunt`, `locate`, `search`, `findex` | bisection / hunting searches (`hunt` lives in `kap/conductive/conductive_opacity.f90`) | all table families |
| `average_over_equipotential` (`qgauss`; in `rotation/shape/shell_quadrature.f90`) | 10-point Gauss–Legendre quadrature | shape integrals |
| `intpol` | spline interpolation with derivative | OPAL 1992 and Kurucz 1990 opacity interpolation |

None of the numerics keeps `save` state (2026-10: the six that did were
shown not to need it, or given explicit state, such as `ratext`'s
tableau). The only routine-local `save`s are in the two driver files
(`core/run_yrec.f90`, `core/evolve_step.f90`) and
`rotation/burn_mix_extrapolated.f90` (commented out). `numerics/README.md` describes the
modules.

## 6. Calibration modes

- **Solar calibration** (`star%ctrl%calibrate_solar_model`, kind-card
  triples; cards prepared by `setup_solar_calibration` (`setcal`)): runs
  the model to the solar age, compares (log L, log R) — and optionally
  Z/X — to solar targets, and updates (X₀, α_MLT) (and Z₀ for Z/X) with
  stored empirical partial derivatives (∂logL/∂X = −3.78 etc., hardwired
  from experiment in `check_solar_calibration` (`chkcal`,
  `core/check_solar_calibration.f90`)), iterating whole evolutionary runs
  to tolerance (`toll/tolr` ≈ 3×10⁻⁵). Yields the classic calibrated solar
  model in a handful of runs.
- **Star calibration** (`star%ctrl%calibrate_star_flag`, namelist `lcalst`;
  `setup_star_calibration`/`check_star_calibration` (`setscal`/`chkscal`)):
  same machinery against arbitrary (Teff, L) targets.
- **Monte-Carlo mode** (`star%ctrl%monte_carlo_active`, namelist `lmonte`;
  `core/monte_carlo.f90`): repeated runs drawing physics
  parameters (rate S-factors, diffusion factors, Z/X…) per realization —
  the classic solar-model uncertainty machinery (compile-verified only: no
  regression deck enables it).

## 7. Program architecture (2026)

```
main (CLI, ~10 lines)                     pyyrec (ctypes) ── libyrec
  └─ run_yrec(ierr)                            └─ yrec_capi: yrec_run()
       ├─ yrec_run_prologue     (re-entrancy: snapshot/restore all state;
       │                         core/yrec_reset.f90)
       ├─ read_controls         (io/read_controls.f90; read_input (parmin):
       │                         legacy $CONTROL/$PHYSICS or MESA-style
       │                         &star_job/&controls inlists)
       ├─ star_setup            (eos_init/kap_init/atm_init table loads)
       └─ run loop (kind cards) ── model loop ── evolve_step
                                                    │
     star layer (owns type(star_info) star): read_starting_model (starin),
       henyey_iterate (crrect), henyey_coefficients (coefft),
       henyey_eliminate (reduce), henyey_solve (hsolve), envelope_integration,
       surfbc, burn_step, composition_step, rezone (hpoint), evolve_angular_momentum
       (getw), wind
                                                    │
     physics services: eos_get, kap_get, net_lib kernels, atm tables,
       numerics_lib (per-point inputs explicit; configuration still read
       from star%ctrl / star%job, see below)
                                                    │
     io: yrec_output (MESA-style history.data / profile{N}.data /
       GYRE/FGONG pulse files; legacy .short/.track/.store retired)
```

All model state hangs off one root — `star` (`state/star_info_lib.f90`),
flat MESA-style: the start-of-step snapshot is the `*_start` family
(`star%xa_start`, `star%logP_start`, …) and driver-step scalars are plain
members (`star%dt`, `star%model_diverged_flag`, …). Only the two input
bundles are nested: `star%job` (run list and calibration cards) and
`star%ctrl` (namelist controls, fixed after `read_controls`). Solver
scratch lives in `rotation/rotation_scratch_lib.f90` and
`core/point_scratch_lib.f90`, table state inside each physics domain. Two
`run_yrec` calls in one process are byte-identical to two fresh processes
(the re-entrancy contract), which is what makes the Python embedding sound.

**The physics boundary today.** Per-point evaluation takes explicit
inputs: `eos_get`, `kap_get` and `net_lib`'s `compute_reaction_rates` receive T, P or ρ and
the composition as arguments, and the envelope kernel
`core/envelope_integration_kernel.f90` is star-blind (enforced by
`tools/check_boundaries.py`). Configuration is still read from `star`: 26
files under `eos/`, `kap/`, `atm/` and `net/` `use star_info_lib, only:
star` and read about 60 distinct members — 54 `star%ctrl`/`star%job`
controls (table choices, paths, units, thresholds, the mixture) plus
derived configuration scalars (`star%cross_section_scale` and the other
rate scales, the envelope X/Z used to prepare tables,
`star%use_two_z_tables`), but no per-zone model state. The target, in
which each library receives its controls once in `*_init` and never
imports `star`, is `TARGET_ARCHITECTURE.md` §4.1 and §4.6.

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

*Caveats: `henyey_solve`'s (`hsolve`) algebra is dense and was renamed
structurally, not rederived, during modernization — the description in
§2.3 reflects the code's own documentation plus the
`henyey_coefficients` (`coefft`) equation forms. The MHD EOS path
is unexercised (no tables ship). Legacy output (`wrtout`) was retired on
2026-08-28, so the Monte-Carlo and helium-flash paths now write through
`io/yrec_output.f90` like every other run; the Monte-Carlo path is
compile-verified only.*

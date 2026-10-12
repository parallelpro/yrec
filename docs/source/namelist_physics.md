Physics Options
===============

YREC is generally called by two namelists: the "control" namelist, or `nml1`; and the "physics" namelist, or `nml2`. This page is a reference that will document various current and deprecated possible commands in the namelist files.

## Physics (`.nml2`)

### Constants

YREC does computations in CGS but outputs in solar units.

| Parameter     | Description    |
| ------------ | -------------- |
| `CLSUN`      | Definition of $L_\odot$ used for tracks and output |
| `CRSUN`      | Definition of $R_\odot$ used for tracks and output |

### Core fitting options

| Parameter | Description                                                                              |
| --------- | ---------------------------------------------------------------------------------------- |
| `LCORE`   | If `.TRUE.`, move the central fitting point inward using a Taylor expansion approximation. |
| `FCORE`   | Factor in mass by which the central fitting point is moved inward.                       |

### Gravitational settling & microscopic diffusion

| Parameter  | Description                                                                                     |
| ---------- | ----------------------------------------------------------------------------------------------- |
| `LDIFY`    | Enable gravitational settling of helium.                                                        |
| `LDIFZ`    | Enable gravitational settling of metals.                                                        |
| `LDIFLI`   | Enable gravitational settling of light species.                                                 |
| `FGRY`     | Scaling factor for helium-settling diffusion coefficients (e.g., 0.75 if no rotational mixing). |
| `FGRZ`     | Scaling factor for metal-settling diffusion coefficients.                                       |
| `DT_GS`    | Settling solved in multiple timesteps if global timestep < `DT_GS * τ(CZ base)`.                       |
| `XMIN`     | Disable diffusion once central hydrogen `X < XMIN`.                                             |
| `YMIN`     | Disable diffusion once surface helium `Y < YMIN`.                                               |
| `GRTOL`    | Numerical tolerance for diffusion abundance changes.                                            |
| `NITER_GS` | Maximum iterations for diffusion abundance solver.                                              |
| `ILAMBDA`  | Formula for Coulomb logarithm, used in collision integrals; 4 (recommended) is a constant, 2.2.                                                    |

### Numerical convergence options

| Parameter     | Description                                                                                             |
| ------------- | ------------------------------------------------------------------------------------------------------- |
| `LNEWS`       | Extrapolate changes from prior step for next step in initial entropy change estimate; can improve convergence.                  |
| `NITER1`      | Iterations before SBC (surface boundary condition) is checked.                                          |
| `NITER2`      | Iterations after SBC check, assumes start of T nuclear burning rates.                                               |
| `NITER3`      | Iterate between model solution and nuclear burning rates (set to 0 if you encounter convergence difficulties).                          |
| `NITER4`      | Iterate betweeen model solution, rotation, and rotationally induced mixing.  Set to 0 for models without rotation.                           |
| `ITDIF1`      | Iterations between structure equations and rotational mixing (1=none). This is the same as NITER3 unless >1. Recommend 2 for rotational mixing models.           |
| `FCORR0`      | The code predicts corrections to the structure variables, which can initially be overestimated. This is the iteration 1 scale factor by which they are multipled.                                      |
| `FCORRI`      | This is the increment to the factor for each successive interation, for example it reaches the max of 1 in 8 steps if (0.6,0.05). larger factor = faster convergence, less stable.                                                      |
| `HTOLER(1,1)` | Convergence tol for log P. If the maximum correction in a given stage anywhere is below this threshold (for all 4 state variables), the run is converged. |
| `HTOLER(2,1)` | Convergence tol for log T. |
| `HTOLER(3,1)` | Convergence tol for log R. |
| `HTOLER(4,1)` | Convergence tol for L/Ltot. Note: all are local |
| `HTOLER(5,1)` | If the difference between analytic and numerical derivatives is better than this, iteration is skipped. (5,1) is for P+R |
| `HTOLER(5,2)` | If the difference between analytic and numerical derivatives is better than this, iteration is skipped. (5,2) is for T+L |

### Convective mixing & overshoot

| Parameter | Description                                                          |
| --------- | -------------------------------------------------------------------- |
| `LSEMIC`  | Enable semi-convection.                                              |
| `LOVSTC`  | Enable core overshoot mixing. Note: unless `LOVAD`=`.TRUE.`, this is overmixing only, and does not impact the thermal structure                                        |
| `ALPHAC`  | Core overshoot in pressure scale heights.                            |
| `LOVMAX`  | Alternative test (only checked if `LOVSTC` = `.TRUE.`) to avoid pathologies in overshoot of small cores.                 |
| `BETAC`   | Alternative test - maximum overshoot in fraction of physical size of core                           |
| `LOVSTE`  | Enable envelope overshoot.                                           |
| `ALPHAE`  | Envelope overshoot in pressure scale heights.                        |
| `LOVSTM`  | Enable overshoot in intermediate convection zones (recommend not on, needs to be checked). |
| `ALPHAM`  | Intermediate overshoot (above and below) in pressure scale heights.                               |
| `TGCUT`   | Adiabatic convection assumed above this log T.                 |

### Equation of state (EOS)

| Parameter   | Description                                                     |
| ----------- | --------------------------------------------------------------- |
| `LSCV`      | Use SCV EOS in domain where it is valid and OPAL is not.                                 |
| `LOPALE06`  | Use OPAL EOS. Note - a simplified Boltzmann/Saha solver, fully ionized at high T, is used when this is disabled, or outside the table.                                                   |
| `LNumDeriv` | Use numerical derivatives instead of analytic ones (very slow). |
| `LDH`       | Use Debye–Hückel correction for the Yale EOS.                                    |
| `ETADH0`    | Lower Debye–Hückel parameter limit.                             |
| `ETADH1`    | Upper Debye–Hückel parameter limit.                             |
| `TSCUT`     | Threshold above which the default EOS assumes full ionization.            |

### Opacity options

| Parameter    | Description                                                     |
| ------------ | --------------------------------------------------------------- |
| `LCondOpacP` | Use updated conductive opacity tables.                          |
| `OPTOL`      | The code generates a 2D surface abundance opacity table, useful for outer CZs and envelope integration. It switches to a 3D/4D solution when x, z differ from it by more than this. |
| `TMOLMIN`    | Below this `log T`, use only molecular opacity (if available).       |
| `TMOLMAX`    | Above this `log T`, use only atomic opacity. Ramp used between `TMOLMIN` and `TMOLMAX`.                         |

### Nuclear reaction & burning options

```{admonition} Numerics
:class: note
The numerics of abundance change solutions and thresholds for computing reaction rates. Only `TCUT1` is adjusted from default. Do not change without good reason. 
```

| Parameter       | Description                                           |
| --------------- | ----------------------------------------------------- |
| `WEAKSCREENING` | Threshold for using weak-screening in reaction rates. |
| `LSNU`          | Compute solar neutrino fluxes.                        |
| `TCUT(1)`       | Temperature threshold for computing $\epsilon _{\mathrm{nuc}}$. Others are thresholds for progressively more advanced burning stages (Default 6.5, bad for D burning in preMS)     |
| `CMIN`          | Minimum abundance tracked in solver.                  |
| `ABSTOL`        | Absolute convergence criterion for abundance changes. |
| `RELTOL`        | Relative convergence criterion for abundance changes. |
| `KEMMAX`        | Maximum substeps for abundance integration.           |

### Light-element burning cross-sections

| Parameter | Description                                 |
| --------- | ------------------------------------------- |
| `LXLI6`   | Use manually specified p+⁶Li cross-section. |
| `LXLI7`   | Use manually specified p+⁷Li cross-section. |
| `XSLI6`   | Value of p+⁶Li cross-section (keV barns).   |
| `XSLI7`   | Value of p+⁷Li cross-section (keV barns).   |
| `LXBE91`  | Toggle p+⁹Be → γ channel.                   |
| `LXBE92`  | Toggle p+⁶Li → d channel.                   |
| `LXBE93`  | Toggle p+⁶Li → α channel.                   |
| `XSBE91`  | Cross-section for p(⁹Be,γ).                 |
| `XSBE92`  | Cross-section for p(⁹Be,d).                 |
| `XSBE93`  | Cross-section for p(⁹Be,α).                 |

### Surface boundary conditions

| Parameter                           | Description                                                           |
| ----------------------------------- | --------------------------------------------------------------------- |
| `KTTAU`                             | SBC model selector (0=Eddington, 1=Krishnaswamy, 3=Kurucz, 4=Allard). |
| `TRIDT`                    | The code solves for a SBC between 3 models. This is the spacing between 2 of them in Log Teff. The third is a different L at the average Teff.                          |
| `TRIDL`                    | Spacing between the lower pair and upper point in Log L.                          |

### Atmospheric integration

| Parameter | Description                                                                                             |
| --------- | ------------------------------------------------------------------------------------------------------- |
| `ATMERR`  | Numerical tolerance for atmosphere integration.                                                         |
| `ATMMAX`  | Maximum step size in log τ for atmosphere integrator.                                                   |
| `ATMMIN`  | Minimum step size in log τ for atmosphere integrator.                                                   |
| `ATMBEG`  | Initial step size in log τ for atmosphere integrator.                                                   |
| `ATMD0`   | log10 of starting density at τ = 0 in the atmosphere integrator (may need adjusting for massive stars). |

### Envelope integrator

| Parameter | Description                                                                                  |
| --------- | -------------------------------------------------------------------------------------------- |
| `ENVERR`  | Numerical tolerance for envelope integration.                                                |
| `ENVMAX`  | Maximum step size in log P for envelope integrator.                                          |
| `ENVMIN`  | Minimum step size in log P for envelope integrator.                                          |
| `ENVBEG`  | Initial step size in log P for envelope integrator.                                          |
| `STOLR0`  | Required fractional agreement between envelope mass and mass at fitting point.               |
| `IMAX`    | Integrators start with a large step, divide to smaller, and extrapolate to zero. Maximum number of loops to infer integration step size (envelope + atmosphere). |
| `NUSE`    | Number of prior steps used for extrapolation to zero step size (atm+env).                                     |

### Spatial tolerances

| Parameter               | Description                                    |
| ----------------------- | ---------------------------------------------- |
| `TOL_DM_MIN`            | Minimum mass spacing (most important for thin outer fitting points).                          |
| `TOL_DM_MAX`            | Maximum log mass spacing (most important for stellar cores).                      |
| `TOL_DL_MAX`            | Maximum spacing in L/Lsurf (most important for energy generating regions).                    |
| `TOL_DP_ENV_MAX`        | Maximum log P spacing in outer CZ (most important for outer layers).             |
| `TOL_DP_CZBASE_MAX`     | Maximum log P spacing at SCZ base.              |
| `TOL_DP_CORE_MAX`       | Maximum log P spacing below SCZ (important between burning, SCZ base).     |
| `FLAG_DX`               | Minimum threshold for preserving discontinuities in X. |
| `FLAG_DZ`               | Minimum threshold for preserving discontinuities in Z. |
| `TOL_CZBASE_FINE_WIDTH` | Width of fine-zoned region at CZ base (0.4 for red giant models).         |
| `FLAG_DW`               | Minimum threshold for preserving Ω discontinuities.    |

### Timestep criteria

| Parameter           | Description                                           |
| ------------------- | ----------------------------------------------------- |
| `LTRIST`            | Restrict timestep during Hertzsprung gap.             |
| `TIME_CORE_MIN`     | Minimum X,Y for using core vs shell burning criteria. |
| `TIME_DX_CORE_TOT`  | Max absolute ΔX in core H-burning.                    |
| `TIME_DX_CORE_FRAC` | Max fractional ΔX in core H-burning.                  |
| `TIME_DY_CORE_TOT`  | Max absolute ΔY/Z in core He-burning (red clump/HB).                 |
| `TIME_DY_CORE_FRAC` | Max fractional ΔY/Z in core He-burning (red clump/HB).               |
| `TIME_DX_TOTAL`     | Max absolute H burned per timestep (post-MS).         |
| `TIME_DX_SHELL`     | Max fractional ΔX in shell H-burning (post-MS).                 |
| `TIME_DY_TOTAL`     | Max absolute He burned (AGB).                         |
| `TIME_DY_SHELL`     | Max fractional ΔY in He shell.                        |
| `LSTRUCT_TIME`      | Restrict timestep based on structural changes - disabled once Lgrav<0 (minimum radius, close to ZAMS).        |
| `TIME_DT`           | Max local ΔT (pre-MS).                                |
| `TIME_DP`           | Max local ΔP (pre-MS).                                |
| `TIME_DR`           | Max local ΔR (pre-MS).                                |
| `TIME_DL`           | Max local ΔL (pre-MS).                                |
| `TIME_MAX_DT_FRAC`  | Max fractional timestep change step-to-step.          |
| `TIME_DW_MIX`       | Max ΔΩ per mixing substep (can be more than one per full model step).                            |
| `TIME_DW_GLOBAL`    | Max ΔΩ per model.                                     |

### Rotation

| Parameter | Description                                     |
| --------- | ----------------------------------------------- |
| `LROT`    | Enable rotation.                                |
| `LNEW0`   | Force new SBC each model when rotation enabled. Required for models with rotation, as rotation is held fixed in the solution. `LROT` = `.TRUE.` sets this `.TRUE.` |

### Rotation & mixing timesteps

| Parameter        | Description                                               |
| ---------------- | --------------------------------------------------------- |
| `TIME_DW_GLOBAL` | Maximum ΔΩ per model step (replaces DTWIND).                                |
| `TIME_DW_MIX`    | Maximum ΔΩ per mixing step (replaces DTDIF).                               |
| `ITDIF2`         | Iterations between nuclear burning and rotational mixing (1=none). |
| `DJOK`           | Numerical tolerance for rotational mixing.                           |

### Structural effects of rotation

| Parameter | Description                                                      |
| --------- | ---------------------------------------------------------------- |
| `ACFPFT`  | Fractional accuracy for departures from spherical symmetry (cubed).       |
| `ITFP1`   | Iterations between equipotential shape and effective radius solutions (structural effects of rotation). |
| `ITFP2`   | Iterations for measuring effective radius given shape parameter (structural effects of rotation).                        |

### Initializing rotation

| Parameter | Description                                 |
| --------- | ------------------------------------------- |
| `LWNEW`   | Convert non-rotating model to rotating.     |
| `WNEW`    | Uniform angular velocity if `LWNEW = TRUE`. If star-disk coupling is being used, should be set to the same value as the disk period. Check your units! |
| `LDISK`   | Enforce disk locking.                       |
| `TDISK`   | Disk lifetime in Gyr. 1e5yr = 1e-4 Gyr is the D burning BL. 1e-4 to 1e-2 is the recommended range                       |
| `PDISK`   | Disk rotation period (rad/s). 8 days (9.0902e-6) is the observed median period in solar-mass disked stars in star forming regions.               |

### Angular momentum loss & winds

| Parameter   | Description                                             |
| ----------- | ------------------------------------------------------- |
| `LJDOT0`    | Enable angular momentum loss.                           |
| `LMWIND`    | Use Matt et al. 2012 wind model.                        |
| `AWIND`     | Wind law selector (`V13` (van Saders & Pinsonneault 2013 model), `K97` (legacy Krishnamurthi et al. 1997), `CUS` is custom. `PMMA`, `PMMB`, `PMMC`, `PMMD`, `PMMM` are only used if this is set to CUS.).          |
| `FK`        | Scale factor for angular momentum loss. Typically used as a free parameter to reproduce the solar rotation at the solar age.      |
| `LROSSBY`   | Scale the loss rate relative to the Rossby number. Should be used with a saturation threshold (`WMAX` or `WMAXSUN` below).                          |
| `WMAX`      | Maximum angular velocity (rad/s) for scaling the mass loss rate and magnetic fields in the scaled solar wind. This is over-written if WMAX_SUN is below 1000; the latter expresses WMAX relative to the solar angular velocity.                    |
| `WMAX_SUN`  | WMAX expressed in solar Ω units.                        |
| `PMMSOLP`   | log P at Teff in solar model (for B-field calibration). |
| `PMMSOLTAU` | Solar convective overturn timescale.                    |

### Angular momentum transport

| Parameter | Description                                                        |
| --------- | ------------------------------------------------------------------ |
| `LINSTB`  | Include rotationally induced mixing and hydrodynamic AM transport. |
| `LSOLID`  | Enforce solid-body rotation everywhere. Set true for strongly coupled magnetic models.                            |
| `WALPCZ`  | Angular velocity in CZ scales as R^WALPCZ. 0 = solid body -2 = uniform specific angular momentum                          |
| `IMPJMOD` | Global angular momentum coupling settings. 0 = no effect 1 = solid body everywhere 2 = solid body core, decoupled cz 3 = solid body core, coupled at base of CZ rate.                                         |
| `LCODM`   | Include magnetic diffusion term.                                   |
| `CODM`    | Magnetic diffusion term for angular momentum transport. 0 = Pure hydro. Nonzero=extra diffusive mechanism. 9e4 matches the solar core rotation constraints.                                    |
| `FW`      | Global factor for scaling the angular momentum transport diffusion coefficients.                 |
| `FC`      | Scale factor for mixing diffusion coefficeints relative to angular momentum. Typicaly calibrated to match Li depletion in a benchmark (cluster stars or the Sun).                              |
| `IGSF`    | Toggle for the stability condition used to evaluate the GSF and ABCD instabilties. 1 is recommended.           |

### Starspots

| Parameter | Description                       |
| --------- | --------------------------------- |
| `LSDEPTH` | Control that descreases spot contrast with depth, assuming a constant B. This does not have a dramatic effect on the results.  |
| `SPOTF`   | Starspot filling factor. 0.25 is the active star bound in main sequence stars.              |
| `SPOTX`   | Starspot temperature contrast; only alter if SPOTF > 0. 0.85 is a typical sunspot value. |

### Rarely used options

| Parameter     | Description                                 |
| ------------- | ------------------------------------------- |
| `LTHOUL`      | Use Thoul formalism for settling.           |
| `LTHOULFIT`   | Use published fit instead of full solution (not recommended). |
| `HTOLER(1,2)` | Abort threshold for log P divergence greater than this.       |
| `HTOLER(2,2)` | Abort threshold for log T divergence greater than this.       |
| `HTOLER(3,2)` | Abort threshold for log R divergence greater than this.       |
| `HTOLER(4,2)` | Abort threshold for L/Ltot divergence greater than this.      |
| `TENV0`       | Minimum T for envelope integration. Code stops below this.         |
| `TENV1`       | Maximum T for envelope integration. Code stops above this.         |
| `ATMSTP`      | Max step for atmosphere integrator (`LENVG`). |
| `ENVSTP`      | Max step for envelope integrator (`LENVG`).   |

### Legacy options

| Parameter    | Description                                       |
| ------------ | ------------------------------------------------- |
| `LEXCOM`     | Extend tracked isotopes to include D, Li6/7, Be9. |
| `LNEWDIF`    | Enable updated settling routines.                 |
| `LNEWNUC`    | Use new nuclear cross-section format.             |
| `LNULOS1`    | Use Itoh et al. (1996) neutrino loss rates.       |
| `LNEWVARS`   | Use updated variable names for tolerances.        |
| `TOL_DX_MAX` | Legacy mesh spacing in X.                         |
| `TOL_DZ_MAX` | Legacy mesh spacing in Z.                         |

### Experimental options

| Parameter  | Description                                              |
| ---------- | -------------------------------------------------------- |
| `DPENV`    | Full mixing from center to fractional mass depth.        |
| `LADOV`    | Experimental adiabatic overshoot.                        |
| `LVFC`     | Zahn 1991 variable-FC model.                             |
| `LDIFAD`   | Not available: the diffusion + advection treatment of the circulation never worked and is commented out. A rotating run with `LDIFAD = .TRUE.` stops. |
| `LMDOT`    | Enable experimental mass-loss module.                    |
| `LSOLWIND` | Not used: the solar-wind mass loss scaled by rotation it would select is commented out. |

### Unused parameters

| Parameter | Description                                          |
| --------- | ---------------------------------------------------- |
| `TCUT(2)` | Not used. Former pp-chain equilibrium threshold.               |
| `TCUT(3)` | Not used. Former CNO-cycle equilibrium threshold.              |
| `MCORE`   | Not used. Obsolete core mesh allocation parameter.             |
| `LENVG`   | Not used. Outer-layer structure calculation after convergence. |
| `FMU`     | Not used. Scaled the μ-gradient velocity, which was always zero and was removed (2026-10). |
| `ALFA`    | Not used. Legacy exponent for Kawaler (1988) wind law.         |

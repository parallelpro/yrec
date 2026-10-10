# YREC glossary

This file is the authority on names in YREC (see `TARGET_ARCHITECTURE.md` section 4.7). Each row is one physical quantity. It gives the symbol, the units, the name or names the code uses now, the legacy F77 name or names, the file that owns the quantity, and a one-line definition. Review checks new code against it.

**Ways to use it**

- **You know the F77 name** (`HL`, `HCOMP`, `DELTS`, `SESUM`, ...): search this file for it. Legacy names appear in upper case.
- **You know the current name**: search for it, e.g. `luminosity_lsun` or `rot_scr%dm`.
- **You are adding physics**: read section 1 (conventions) and section 3 (naming traps) first.

**Columns**

| Column | Meaning |
|---|---|
| Symbol | Textbook symbol. Subscript *i* is a species and *k* is a zone. |
| Units | cgs unless stated. "dex" means a log₁₀ quantity. "—" means dimensionless or a flag. "unconfirmed" means the units could not be settled from the code. |
| Current name(s) | Name in the code today. Several names separated by commas are copies or aliases of the same quantity. "also" lists the name a different routine uses for it. |
| Legacy | F77 name before the 2026 modernization: a COMMON member, a `main.f` array, or a namelist variable. An arrow (`A → B`) records an intermediate 2026 name. "—" means the quantity is new since 2026. |
| Owner | File that computes (writes) the value. Every `star%` member is *declared* in `state/star_info_lib.f90`. `rot_scr`, `mix_scr` and `circ_scr` are declared in `rotation/rotation_scratch_lib.f90`, and `pt_scr` in `core/point_scratch_lib.f90`. Controls are declared in the generated `state/controls_state_def.inc` (`ctrl`) and `state/job_controls_def.inc` (`job`). |

**Naming rules** (from `TARGET_ARCHITECTURE.md` 4.7, plus what the code already practises):

1. `_start` means the value at the start of the step. Section 3 lists the places where that is not yet true.
2. Derivatives are named `d<f>_dln<x>`, and the declaration comment says whether the derivative is **absolute** (d*f*/d ln *x*) or **logarithmic** (d ln *f*/d ln *x*). Several current names break this rule (section 3, traps 1–3).
3. No name is reused for a different quantity.
4. A `log` prefix means log₁₀ (`logP`, `log_mass`, `log_L`). `ln` appears only in derivative names, where it means a natural log.
5. A unit suffix is used where units are not cgs: `_lsun`, `_msun`, `_rsun`, `_gyr`, `_yr`, `_kms`, `_days`, `_cgs`.
6. Legacy names stay in a comment at each declaration. `TARGET_ARCHITECTURE.md` Appendix A maps legacy routines to current files.

Line numbers cited here are for the working tree of 2026-10-10. `src/` is being edited, so if a line has moved, search for the quoted name.

---

## 1. Conventions

### 1.1 Units and logarithms

| Quantity | Convention |
|---|---|
| General | cgs: g, cm, s, K, erg, dyn cm⁻². |
| Luminosity in the structure | `luminosity_lsun` is **linear** L/L☉. The interface copies used by rotation (`rot_scr%interface_luminosity`) and the pulse files are in erg s⁻¹. |
| Mass | `m`, `dm` in g; `log_mass`, `log_total_mass` in log₁₀ g. Observables (`star_mass`, `core_cz_mass`, `envelope_mass`, `h_shell_*_mass`) are in M☉. |
| Radius | `logR` in log₁₀ cm. `log_R_surface` and `envelope_cz_log_radius` are log₁₀(*r*/R☉). `envelope_radius` is in R☉, and `h_shell_*_radius` is *r*/R. |
| Time | `star%dt` in s, `timestep_yr` in yr, the age `dage` in Gyr. `seconds_per_year` = 3.1558×10⁷ (`setup/setups.f90:57`). |
| Stored logs | All log₁₀. ln 10 is `ln10` (also `clndp`, legacy `CLN`); 1/ln 10 is `clni` (`CLNI`). |
| Specific energy rates | erg g⁻¹ s⁻¹. Exception: `eps_channels(1:5,:)` holds fractions (trap 4). |
| Solar units | Set once by `setup/setups.f90` from the namelist values `clsun`/`crsun` and the control `solar_mass_cgs`. The luminosity values are rescaled per Monte-Carlo run. |

### 1.2 Grid and zone indexing

- Zone 1 is the innermost mesh point and zone `nz` is the **fitting point**, the outermost interior point. Above `nz` the envelope and atmosphere are separate integrations (`env_struct`, `atmo_struct`). The stitched model (`core/stitched_model.f90`) joins all three for output.
- `dm(k)` is a point-centred mass **weight**, not the mass between two points: `dm(k) = ½(m(k+1) − m(k−1))`. `dm(nz)` also contains the entire envelope mass above the fitting point, so Σ`dm` = M (`setup/rezone.f90:723-735`).
- The envelope integration (`core/envelope_derivs.f90`) uses log₁₀ P as its independent variable and carries `y(1)` = log₁₀(m/M), `y(2)` = log₁₀ T and `y(3)` = log₁₀ r, at constant L.

### 1.3 The Henyey difference equations

The authoritative statement is the header of `core/henyey_eliminate.f90`. In summary:

- **Dependent variables:** log₁₀ P, log₁₀ T, log₁₀ R and *L* in L☉ (linear). **Coordinate:** log₁₀ m.
- **Equation terms:** `q*` = d(variable)/d log₁₀ m.
- **Residual for equation *k* across shell *i*:** `elim_rhs(k,i) = ½ Δlog₁₀m · (q_k(i−1) + q_k(i)) − (y_k(i) − y_k(i−1))`.
- **Partials:** `q*_dp`, `q*_dt`, `q*_dr` are **natural-log** derivatives ∂q/∂ln P, ∂q/∂ln T, ∂q/∂ln R. The eliminator multiplies ½Δlog₁₀m by ln 10 once, which turns them into derivatives with respect to the log₁₀ unknowns.
- **Luminosity partial:** `q*_dl` is ∂q/∂L **divided by ln 10**, for example `qt_dl = qt/(L ln10)`, so that the same factor leaves it as ∂q/∂L.
- **Luminosity term:** `ql = (ln10 · m / L☉) · ε`, with ε in erg g⁻¹ s⁻¹ (`cccql = ln10_over_lsun · m`). Therefore `ql_dt` and `ql_dp` take the **absolute** derivatives dε/d ln T and dε/d ln ρ, which is what `nuclear_energy` returns and `compute_zone_energy` passes on (`burn%deps_dlnt`, `burn%deps_dlnrho`).
- **Missing partials:** `henyey_shell_terms` has no `ql_dr` or `ql_dl`. An energy term that depends on R or L cannot be linearized without extending the type and the pivots.
- **Boundaries:** the centre condition is set at zone 1 in `henyey_coefficients`. The surface condition comes from the envelope fit (`surface_bc`, built from `envelope_fit_coeffs` in `henyey_iterate`).

### 1.4 Timestep sign convention

- **A negative dt means "do not age this model".** `star%dt` (legacy `DELTS`) is set to `-dabs(star%dt)` when the model is not to be aged, i.e. on rescale-only kind cards (`core/evolve_step.f90:325-331`).
- **Code that reads the sign:**
  - `henyey_coefficients` includes the entropy (ε_grav) terms only for `delta_time > 0`;
  - the improved-first-guess test in `evolve_step`;
  - the rotation step (`evolve_angular_momentum`), which also does arithmetic on dt;
  - `compute_timestep`, whose header says it exits on a negative dt (`evolve_step` makes dt positive before its call).
- **Reset:** dt is made positive again before `compute_timestep` (`evolve_step.f90:141`).
- **Target:** the flag `s%step%aging` replaces the sign (`TARGET_ARCHITECTURE.md` 4.2).

### 1.5 Species index constants

Declared in `state/star_info_lib.f90:40-55`. The composition array is `star%xa(n_species_extended, json)`. The **basic** layout uses slots 1–11 (`n_species_basic = 11`). The **extended** layout (`use_extended_composition`, legacy `LEXCOM`) adds slots 12–15 (`n_species_extended = 15`). Arrays are always sized for 15.

| Constant | Slot | Species | Note |
|---|---|---|---|
| `i_h1` | 1 | ¹H (X) | |
| `i_he4` | 2 | ⁴He (Y) | |
| `i_metals` | 3 | Z (total metals) | YREC evolves Z as one species. Slots 5–11 and 13–15 are tracers carried *inside* Z, not additions to it. |
| `i_he3` | 4 | ³He | |
| `i_c12` | 5 | ¹²C | |
| `i_c13` | 6 | ¹³C | |
| `i_n14` | 7 | ¹⁴N | |
| `i_n15` | 8 | ¹⁵N | |
| `i_o16` | 9 | ¹⁶O | |
| `i_o17` | 10 | ¹⁷O | |
| `i_o18` | 11 | ¹⁸O | last slot of the basic set |
| `i_h2` | 12 | ²H (D) | extended only |
| `i_li6` | 13 | ⁶Li | extended only |
| `i_li7` | 14 | ⁷Li | extended only |
| `i_be9` | 15 | ⁹Be | extended only |

The 12-element **surface mixture** has its own slot constants `ix_na`=1, `ix_al`, `ix_mg`, `ix_fe`, `ix_si`, `ix_c`, `ix_h`=7, `ix_o`, `ix_n`, `ix_ar`, `ix_ne`, `ix_he`=12 (`n_mix_species = 12`; `star_info_lib.f90:61-64`). It is used by `job%mixture_weights`, `ctrl%mixture_weights_seed` and `fxenv`.

### 1.6 Snapshot families

| Pattern | Meaning today |
|---|---|
| `*_start` (structure: `logP_start`, ...) | The previous converged model, written by `rezone` onto the new grid (or by `evolve_step`'s `rezone_or_snapshot` when rezoning is skipped). |
| `xa_start` | Composition at the start of the step, stored after mass loss (`evolve_step.f90:245`). The rotation sub-steps overwrite it (trap 9). |
| `old_*` (rotation: `old_omega`, ...) | Previous model's rotation state on the new grid, stored by `rezone` (`setup/rezone.f90:789-797`). Legacy `/OLDROT/`. |
| `orig_*` | Start-of-step copies kept across the structure⇄rotation iterations (`evolve_step.f90:414-431`). Legacy `/ORIGSTART/`. |
| `*_prev` | Previous iteration, sub-step or model, depending on the routine. Read the declaration. |
| `rot_scr%*_mid` | Structure interpolated to the middle of the current rotation sub-step (`rotation/mid_timestep_model.f90`). |
| `eq_*`, `*_eqgrid` | Values on the equally spaced χ grid of the rotation solver. |

---


### Units rule (2026-10)

Internal variables are CGS unless their name carries a unit suffix:
`mass` g, `mass_msun` M☉; `radius` cm, `radius_rsun` R☉; `luminosity`
erg/s, `luminosity_lsun` L☉; `age` s, `age_gyr` Gyr; `rate` per s,
`rate_per_gyr` per Gyr. `log_x` is log₁₀ of x in the unit its name
implies. Output columns keep their MESA names. Many current names
predate the rule (the solar-unit masses and luminosity channels, `dage`,
the per-Gyr rate arrays); the tables
below record their real units, and renames follow as code is touched.

## 2. Quantities by domain

### 2.1 Grid and structure (per zone unless noted)

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| N | — | `star%nz` | `M` → `num_zones` | setup/rezone.f90 | Number of interior mesh points; zone `nz` is the fitting point. |
| — | — | `json` (parameter, 5000) | `JSON` | state/star_info_lib.f90 | Fixed capacity of every per-zone array. |
| log₁₀ m | dex (g) | `star%log_mass` | `HS` | setup/rezone.f90 | Enclosed mass at the point; the Henyey coordinate. |
| m | g | `star%m` | `HS1` → `enclosed_mass` | setup/rezone.f90 | Enclosed mass, 10^`log_mass`. |
| Δm | g | `star%dm` | `HS2` → `shell_mass` | setup/rezone.f90 | Point-centred mass weight ½(m(k+1)−m(k−1)); `dm(nz)` includes the whole envelope (section 1.2). |
| log₁₀ P | dex (dyn cm⁻²) | `star%logP` | `HP` → `log_pressure` | core/henyey_iterate.f90 | Total pressure; Newton corrections are applied here. |
| log₁₀ T | dex (K) | `star%logT` | `HT` → `log_temperature` | core/henyey_iterate.f90 | Temperature. |
| log₁₀ r | dex (cm) | `star%logR` | `HR` → `log_radius` | core/henyey_iterate.f90 | Radius. With rotation, the volume-equivalent radius r_ψ of the isobar (`rotation/shape/shape.f90`). |
| log₁₀ ρ | dex (g cm⁻³) | `star%logRho` | `HD` → `log_density` | core/henyey_coefficients.f90 | Density returned by the EOS at (P, T). |
| L | L☉ (linear) | `star%luminosity_lsun` | `HL` | core/henyey_iterate.f90 | Luminosity through the point. |
| — | logical | `star%convective_flag` | `LC` | core/henyey_coefficients.f90, core/shell_physics.f90 | Point is convective according to `temperature_gradients`. |
| log₁₀ M | dex (g) | `star%log_total_mass`, `star%stotal` | `HSTOT`, `STOTAL` | core/read_starting_model.f90:312, wind/mdot.f90:249-250 | Total stellar mass. Two copies, kept equal by hand (trap 17). |
| M | M☉ | `star%star_mass` | `SMASS` → `total_mass_msun` | core/read_starting_model.f90 (model read), setup/rescale_model.f90, wind/mdot.f90:242 | Total stellar mass in solar units. |
| log₁₀ L | dex (L☉) | `star%log_L` | `BL` → `log10_luminosity` | core/henyey_iterate.f90 | log₁₀ `luminosity_lsun(nz)`, the surface luminosity. |
| log₁₀ T_eff | dex (K) | `star%log_Teff` | `TEFFL` → `log_teff` | core/henyey_iterate.f90 | Effective temperature, from `envelope_fit_coeffs(7:9)`. |
| log₁₀(m_fit/M) | dex | `star%senv` | `SENV` | core/henyey_iterate.f90:86 | `log_mass(nz) − log_total_mass`: mass fraction inside the fitting point. |
| log₁₀ P etc. at step start | as above | `star%logP_start`, `logT_start`, `logR_start`, `logRho_start`, `luminosity_lsun_start` | `HPO`, `HTO`, `HRO`, `HDO`, `HLO` (`/OLDMOD/`) → `old_pressure` ... | setup/rezone.f90, core/evolve_step.f90 | Previous converged model on the current grid (section 1.6). |
| log₁₀ T_eff, N at step start | dex, — | `star%log_Teff_start`, `star%nz_start` | `TEFFLO`, `MO` | setup/rezone.f90, core/evolve_step.f90 | Used by the HR-track timestep limiter (`util/timestep_limit_hr.f90`). |
| — | logical | `star%convective_flag_start` | `LCO` | setup/rezone.f90 | Previous model's `convective_flag`. (`cz_flag_start`, legacy `LCZO`, was write-only and was removed in 2026-10.) |
| — | — | `star%old_shell_mass` | `HSS` | none | **Unused.** Historically held log-mass coordinates, not shell masses (trap 6). |

### 2.2 Henyey solver: equation terms and work arrays

The `cur%`/`prev%` terms are components of `henyey_shell_terms` (`core/henyey_eliminate.f90`). The legacy names for the shell below have a `0` suffix (`QP0`, `QPR0`, ...).

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| q_P | — | `cur%qp`, `pt_scr%qp` | `QP` | core/henyey_coefficients.f90 | dlog₁₀P/dlog₁₀m = −G m² f_P/(4π r⁴ P). |
| ∂q_P/∂ln r | — | `cur%qp_dr` | `QPR` | core/henyey_coefficients.f90 | = −4 q_P. |
| ∂q_P/∂ln P | — | `cur%qp_dp` | `QPP` | core/henyey_coefficients.f90 | = −q_P. |
| q_T | — | `cur%qt`, `pt_scr%qt` | `QT` | core/henyey_coefficients.f90 | dlog₁₀T/dlog₁₀m = ∇ · q_P. |
| ∂q_T/∂ln r | — | `cur%qt_dr` | `QTR` | core/henyey_coefficients.f90 | −4 q_T, plus q_T · ∂ln∇/∂ln R in convective zones. |
| ∂q_T/∂L ÷ ln10 | L☉⁻¹ | `cur%qt_dl`, `pt_scr%qtl` | `QTL` | core/henyey_coefficients.f90 | q_T/(L ln10) where radiative; 0 where convective. |
| ∂q_T/∂ln P | — | `cur%qt_dp` | `QTP` | core/henyey_coefficients.f90 | Radiative: q_T κ_ρ (∂lnρ/∂lnP). Convective: q_T(−1 + ∂ln∇/∂lnP). |
| ∂q_T/∂ln T | — | `cur%qt_dt` | `QTT` | core/henyey_coefficients.f90 | Radiative: q_T(−4 + κ_T + κ_ρ ∂lnρ/∂lnT). Convective: q_T ∂ln∇/∂lnT. |
| q_r | — | `cur%qr` | `QR` | core/henyey_coefficients.f90 | dlog₁₀r/dlog₁₀m = m/(4π r³ ρ). |
| ∂q_r/∂ln r | — | `cur%qr_dr` | `QRR` | core/henyey_coefficients.f90 | = −3 q_r. |
| ∂q_r/∂ln P | — | `cur%qr_dp` | `QRP` | core/henyey_coefficients.f90 | = −q_r ∂lnρ/∂lnP. |
| ∂q_r/∂ln T | — | `cur%qr_dt` | `QRT` | core/henyey_coefficients.f90 | = −q_r ∂lnρ/∂lnT. |
| q_L | L☉ per dex of m | `cur%ql` | `QL` | core/henyey_coefficients.f90 | dL/dlog₁₀m = (ln10 m/L☉)(ε + ε_grav − ε_rot), with ε = `burn%total_energy_gen_rate`. |
| ∂q_L/∂ln P | L☉ | `cur%ql_dp` | `QLP` | core/henyey_coefficients.f90 | `cccql`·(dε/dlnρ · ∂lnρ/∂lnP + entropy-term partials); absolute dε. |
| ∂q_L/∂ln T | L☉ | `cur%ql_dt` | `QLT` | core/henyey_coefficients.f90 | `cccql`·(dε/dlnT + dε/dlnρ · ∂lnρ/∂lnT + entropy-term partials). |
| ln10 m/L☉ | g s erg⁻¹ | `cccql` (local) | `CCCQL` | core/henyey_coefficients.f90:308 | `ln10_over_lsun · m`; converts ε to q_L. |
| (∂lnρ/∂lnP) at fixed T | — | `dlnrho_dlnp` (henyey argument; also `star%dlnrho_dlnp`, `eos_res(i_dlnrho_dlnp)`) | `QDP` | eos/eos_lib.f90 (via `eos_get`) | Logarithmic density derivative; 1/χ_ρ. |
| (∂lnρ/∂lnT) at fixed P | — | `dlnrho_dlnt` (also `star%dlnrho_dlnt`, `eos_res(i_dlnrho_dlnt)`) | `QDT` | eos/eos_lib.f90 | Logarithmic density derivative; −δ. |
| ∂ln∇/∂ln{T,P,R} | — | `dgrad_dt_component`, `dgrad_dp_component`, `dgrad_dr_component` | `QACT`, `QACP`, `QACR` | mixing/temperature_gradients.f90 | Logarithmic derivatives of the convective gradient. |
| — | mixed | `star%elim_coeff(4,2,json)` | `HHA` | core/henyey_eliminate.f90 | Recursion coefficients expressing (P,T)-row unknowns of shell *i* in terms of (R,L). |
| — | dex / L☉ | `star%elim_rhs(4,json)` | `HHB` | core/henyey_eliminate.f90, core/henyey_solve.f90 | Residuals during elimination; overwritten in place with the corrections (ΔlogP, ΔlogT, ΔlogR, ΔL) by `henyey_solve` (trap 11). |
| — | mixed | `star%surface_bc(6)` | `HSRF` | core/henyey_iterate.f90:202-212 | Outer-boundary rows built from `envelope_fit_coeffs`. |
| — | dex / L☉ | `star%max_residual(4)` | `HMAX` | core/henyey_eliminate.f90, core/henyey_iterate.f90 | Holds three different things within one iteration (trap 10). |
| — | — | `star%max_correction_index(4)` | `NHMAX` | core/henyey_iterate.f90 | Zone of the largest correction per variable. |
| Δlog₁₀P over step | dex | `star%log_pressure_delta` | `HPP` | core/evolve_step.f90, core/henyey_iterate.f90 | Change of logP since the step began (first guess plus accumulated corrections); enters ε_grav. |
| Δlog₁₀T over step | dex | `star%log_temperature_delta` | `HTT` | core/evolve_step.f90, core/henyey_iterate.f90 | Same for logT. |
| dlog₁₀T/dt etc. | dex s⁻¹ | `star%temperature_entropy_term`, `pressure_entropy_term`, `radius_entropy_term`; `luminosity_entropy_term` (≈ dlnL/dt, s⁻¹) | `HTT0`, `HPP0`, `HRR0`, `HLL0` (`/ENTROP/`) | core/evolve_step.f90:388-391 | Rates of change from the last step, used to extrapolate the next first guess. Not entropy (trap 12). |
| f_corr | — | `star%job%fcorr` | `FCORR` | core/evolve_step.f90, core/henyey_iterate.f90 | Newton under-relaxation factor applied to the corrections. |

### 2.3 Thermodynamics, opacity and transport coefficients (per zone)

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| ∇_rad | — | `star%gradr` | `SDEL(1,:)` → `DEL_GRAD(1,:)` | core/henyey_coefficients.f90, core/shell_physics.f90 | Radiative temperature gradient. |
| ∇ | — | `star%gradT` | `SDEL(2,:)` → `DEL_GRAD(2,:)` | as above | Actual temperature gradient dlnT/dlnP. |
| ∇_ad | — | `star%grada` | `SDEL(3,:)` → `DEL_GRAD(3,:)` | as above | Adiabatic gradient (from the EOS). |
| κ | cm² g⁻¹ | `star%opacity_zone` | `SO` | as above | Rosseland opacity. Not ¹⁶O (an earlier rename guessed `o16_zone`). |
| β | — | `star%beta` | `SBETA` | core/henyey_coefficients.f90 | Gas-pressure fraction from the EOS. |
| η | — | `star%eta` | `SETA` | core/henyey_coefficients.f90 | Electron degeneracy parameter. |
| — | — | `star%fxion_zone(3,:)` | `SFXION` | core/henyey_coefficients.f90 | Three ionization fractions from the EOS; the species order per slot is unconfirmed. |
| v_conv | cm s⁻¹ | `star%conv_vel` | `SVEL` | core/henyey_coefficients.f90, core/shell_physics.f90 | MLT convective velocity (10⁻¹¹ placeholder where radiative). |
| c_P | erg g⁻¹ K⁻¹ | `star%scp` | `SCP` | core/henyey_coefficients.f90:366 | Specific heat from the Henyey pass; read by the profile writers. |
| c_P | erg g⁻¹ K⁻¹ | `star%cp` | `CP` (`/TEMP/`) | core/shell_physics.f90:127 | Same EOS quantity, filled by `shell_physics` for the rotation pipeline (trap 18). |
| μ | — | `star%mu` | `HAMU` → `mean_molecular_weight` | core/shell_physics.f90:125 | Mean molecular weight assuming full ionization. Differs from `pulse_mean_molecular_weight` (trap 18). |
| (∂lnρ/∂lnT) at fixed P | — | `star%qdt` | `SQDT` | core/shell_physics.f90:128 | Per-zone copy for the rotation pipeline. |
| Γ₁ | — | `star%adiabatic_index_gamma1` | `GAM1` (`/SOUND/`) | core/henyey_coefficients.f90 | First adiabatic exponent χ_ρ c_P/c_V (Guenther 1995). |
| ν | cm² s⁻¹ | `star%visc` | `VISC` | rotation/viscos.f90 | Kinematic viscosity, radiative plus molecular. |
| K | cm² s⁻¹ | `star%thdif` | `THDIF` | rotation/viscos.f90 | Radiative thermal diffusivity 16σT³/(3κρ²c_P). |
| EOS outputs | mixed | `eos_res(i_*)`: `i_log10_density`, `i_beta`, `i_fxion`(3), `i_gas_constant` (R/μ), `i_mu_ion_inv`, `i_mu_e_inv`, `i_eta`, `i_dlnrho_dlnt`, `i_dlnrho_dlnp`, `i_cp`, `i_grada`, and the second derivatives `i_*_dt`, `i_*_dp` | positional `eqstat` arguments | eos/eos_lib.f90:71-80 | Named slots of the 24-element `eos_get` result. `i_gas_constant` is R/μ, not μ. |
| opacity outputs | mixed | `kap_res(i_kap)`, `i_log10_kap`, `i_dlnkap_dlnrho`, `i_dlnkap_dlnt` | `O`, `OL`, `QOD`, `QOT` | kap/kap_lib.f90:36-38 | κ, log₁₀ κ and the **logarithmic** derivatives κ_ρ and κ_T. |
| mid-step copies | as base | `mix_scr%amum`, `cpm`, `delm`, `del_adiabatic_mix`, `del_radiative_mix`, `esumm`, `om`, `qdtm`, `thdifm`, `velm`, `viscm`, `epsm` | `AMUM`, `CPM`, `DELM`, `DELAM`, `DELRM`, `ESUMM`, `OM`, `QDTM`, `THDIFM`, `VELM`, `VISCM`, `EPSM` (`/MDPHY/`) | rotation/mid_timestep_model.f90 | μ, c_P, ∇, ∇_ad, ∇_rad, ε_nuc, κ, ∂lnρ/∂lnT, K, v_conv, ν interpolated to mid-sub-step. `epsm` is ε_nuc+ε_ν+ε_grav; `esumm` is ε_nuc only. |

### 2.4 Composition

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| X_i | mass fraction | `star%xa(i,k)` | `HCOMP` → `composition` | mixing/mix.f90, rotation/*, core/burn_step.f90 | Abundance of species slot *i* (section 1.5) in zone *k*. |
| — | — | `n_species_basic` (11), `n_species_extended` (15) | `NSPEC` (local) | state/star_info_lib.f90:54-55 | Species counts of the two layouts. |
| X_i at step start | mass fraction | `star%xa_start` | `HCOMPP` (`/OLDMOD/`) | core/evolve_step.f90:245 | Start-of-step composition. Burning and mixing start from it, and Newton levels > 2 reset `xa` to it (`henyey_iterate.f90:136`). Also a rotation sub-step checkpoint (trap 9). |
| X_i at step start (kept) | mass fraction | `star%orig_composition` | `HCOMPP0` (`/ORIGSTART/`) | core/evolve_step.f90:420 | Restores `xa_start` on structure⇄rotation passes after the first. |
| X_i snapshot | mass fraction | `rot_scr%composition_snapshot` | `HCOMPS` (`/OLDAB/`) | rotation/* | Composition saved before a rotational-diffusion iteration; restored between iterations. |
| X_env | mass fraction | `star%envelope_hydrogen_fraction` | `XENV` | core/read_starting_model.f90 | Surface X used for the envelope and atmosphere. |
| Y_env | mass fraction | `star%envelope_helium_fraction` | `YENV` (`/COMP2/`) | core/read_starting_model.f90 | 1 − X_env − Z_env − X(³He). |
| Z_env | mass fraction | `star%envelope_metal_fraction` | `ZENV` | core/read_starting_model.f90 | Surface Z. |
| X(³He)_env | mass fraction | `star%envelope_he3_fraction` | `Y3ENV` | core/read_starting_model.f90 | Surface ³He. |
| X, Z for the surface BC | mass fraction | `star%xnew`, `star%znew` | `XNEW`, `ZNEW` | core/read_starting_model.f90, setup/rescale_model.f90 | Envelope X and Z passed to `surfbc`. |
| Z − Z_CNO | mass fraction | `star%zenvm` | `ZENVM` | core/read_starting_model.f90:901 | Surface Z excluding C, N and O. |
| 1/μ_ion | — | `star%amuenv` | `AMUENV` | core/read_starting_model.f90:915 | Σ X_i/A_i of the 12-element surface mixture (the inverse mean ion weight; a nearby comment calls it "mean atomic weight"). |
| — | number fraction | `star%fxenv(12)` | `FXENV` | core/read_starting_model.f90 | Number fractions of the surface-mixture elements (`ix_*` order) for the EOS. |
| — | — | `star%job%mixture_weights(12)`, `star%ctrl%mixture_weights_seed(12)` | `V` | setup/setups.f90, core/read_starting_model.f90 | Relative abundances of the 12-element surface mixture. |

### 2.5 Energy generation and nuclear network

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| ε_nuc | erg g⁻¹ s⁻¹ | `star%eps_total` | `SESUM` | core/henyey_coefficients.f90:342 | Sum of the five chain rates (pp I incl. d(p,γ), pp II, pp III, CNO, 3α). Excludes ν losses and α captures. |
| ε_chain/ε_nuc | fraction | `star%eps_channels(1:5,k)` (`i_eps_pp1`, `i_eps_pp2`, `i_eps_pp3`, `i_eps_cno`, `i_eps_he3`) | `SEG(1:5)` | core/henyey_coefficients.f90:346-354 | **Fractions** of `eps_total`, not erg g⁻¹ s⁻¹ (trap 4). Slot 5 holds the 3α share despite its `he3` name (trap 5). |
| ε_ν | erg g⁻¹ s⁻¹ | `star%eps_channels(i_eps_neu,k)`; `burn%eps_neutrino_loss` | `SEG(6)`; `ENU` (`/NEWEPS/`) | core/henyey_coefficients.f90:345; net/nuclear_energy.f90:547, 587 (`thermal_neutrino_loss`) | Thermal neutrino loss rate (stored negative). |
| ε_grav | erg g⁻¹ s⁻¹ | `star%eps_channels(i_eps_grav,k)`, `star%gravitational_luminosity(k)` | `SEG(7)`, `HHC` | core/henyey_coefficients.f90:289-290 | −T ds/dt from `log_*_delta`. Uses 1/yr instead of 1/dt when X > 0.01 and dt < 1 yr (line 280-281). `gravitational_luminosity` is a specific rate, not a luminosity (trap 7). |
| ε | erg g⁻¹ s⁻¹ | `burn%total_energy_gen_rate` (local `energy_gen_rate`) | `SUM1` (`engeb`), `E` (`coefft`) | net/nuclear_energy.f90 (`sum_energy_generation`) | Total nuclear rate including deuterium burning and (negative) ν losses; the ε in q_L. |
| ε_ppI, ε_ppII, ε_ppIII, ε_CNO, ε_3α | erg g⁻¹ s⁻¹ | `burn%pp_chain_energy_gen`, `he3he4_be7_electron_energy_gen`, `he3he4_be7_proton_energy_gen`, `cno_cycle_energy_gen`, `triple_alpha_energy_gen` | `EPP1`, `EPP2`, `EPP3`, `ECN`, `E3AL` | net/nuclear_energy.f90:423-437 | Chain energy rates. ppI includes d(p,γ). ppII is the ⁷Be electron-capture branch, ppIII the ⁷Be proton-capture branch. |
| dε/dln ρ | erg g⁻¹ s⁻¹ | `burn%deps_dlnrho`; also `zone_deps_dlnrho`, `star%pulse_deps_dlnrho`, `rot_scr%deps_dlnrho` | `PEP`; `QED`, `PQED`, `SQED` | net/nuclear_energy.f90:356 | **Absolute** Σ εⱼ dln rⱼ/dlnρ (minus ν-loss term). The copies were named `*dlneps*`/`*dlnepsilon*` until 2026-10 (trap 1). |
| dε/dln T | erg g⁻¹ s⁻¹ | `burn%deps_dlnt`; also `zone_deps_dlnt`, `star%pulse_deps_dlnt`, `rot_scr%deps_dlnt` | `PET`; `QET`, `PQET`, `SQET` | net/nuclear_energy.f90:357 | **Absolute** dε/dlnT (trap 1). |
| ε_αcap | erg g⁻¹ s⁻¹ | `burn%eps_alpha_capture` | `EALPCA` (`/NEWEPS/`) | net/nuclear_energy.f90 (`chain_energy_shares`) | Energy from ¹³C, ¹²C and ¹⁴N α captures. |
| L_channel | L☉ | `star%luminosity_breakdown(8)` (`i_lum_pp1`, `i_lum_pp2`, `i_lum_pp3`, `i_lum_cno`, `i_lum_3alpha`, `i_lum_neu`, `i_lum_grav`, `i_lum_he_c`) | `TLUMX` | core/henyey_coefficients.f90; renormalized in core/henyey_iterate.f90:173-182 | Integrated luminosity per channel. Rescaled each iteration so the sum equals L(nz) (trap 15). |
| ε_rot | erg g⁻¹ s⁻¹ | `rot_scr%rotational_energy_term` | `DEROT` (`/ROTEN/`) | core/henyey_coefficients.f90:302 | (E_rot − E_rot,old)/(Δt Δm), subtracted from ε in q_L. |
| E_rot | erg (per point) | `star%kinetic_energy_rot` | `HKEROT` | core/henyey_iterate.f90:320, rotation/* | ½ Ω j Δm of the point. |
| E_rot at step start | erg | `star%kinetic_energy_rot_old` | `HKEROTO` | core/evolve_step.f90:297 | Start-of-step copy (an `_old`, not `_start`, name). |
| r_j (mixing) | as above | `rot_scr%reaction_rate_by_zone(15,k)` (rows `rr_*`) | `HCOMPM` (`/BURN/`) | core/burn_step.f90 (`compute_burn_rates`, from `network_rate_vector`) | Rates and branchings packed for the network burn. Rows 9 and 13 are zero placeholders. |
| — | Gyr⁻¹ amu⁻¹ | `gyr_amu_per_sec_gram` (parameter) | `C21` | state/phys_const_lib.f90 | Converts a rate per s per g into the units of `reaction_rate`. |
| Φ_ν | 10¹⁰ per g; totals 10¹⁰ (unconfirmed) | `burn%neutrino_flux(10)` (slots `i_nu_*`); `neutrino_flux_total(10)`; `neutrino_flux_zone(10,k)` | `FLUX`, `FLUXTOT` (`/FLUXES/`); `FLUXI0` | net/nuclear_energy.f90 (`solar_neutrino_fluxes`); totals in core/henyey_coefficients.f90 | Solar-neutrino source rates: one shell's record, Σ·Δm over zones, and per zone for the table. Slots 9–10 are the ³He+³He and ³He+⁴He rates, not neutrinos. |
| — | SNU | `star%cl37_snu_rate`, `star%ga71_snu_rate` | `CLSNU`, `GASNU` | core/observables_lib.f90 | ³⁷Cl and ⁷¹Ga capture rates from the flux totals. |
| X(⁷Be) | mass fraction | `burn%be7_mass_fraction`, `star%be7_mass_fraction_zone` | `BE7MASSFRACTION`, `BE7M` | net/nuclear_energy.f90 (`solar_neutrino_fluxes`); per zone core/neutrino_flux_table.f90 | Equilibrium ⁷Be abundance (one shell's record, and per zone; 0 where X ≤ 10⁻⁶). |
| r_D | per Gyr per amu | `star%deuterium_burning_rate`, `deuterium_burning_rate_start` | `DRATE`, `DRATE0` (`/DEUTER/`) | core/burn_step.f90:1346 (`compute_burn_rates`, from `deuterium_rate`) | Deuterium burning rate, end and start of step. |
| ΔM_acc | g | `star%accreted_mass` | `FMASSACC` | wind/mdot.f90:272; read by `burn_deuterium_region` | Mass accreted during the step; the accreted composition is blended into the burned surface region by mass. Named `accreted_mass_fraction` until 2026-10. |
| r_Li6, r_Li7, r_Be9 | s⁻¹ (N_A⟨σv⟩·ρX_H, taking A_H = 1) | `star%rate_li6`, `rate_li7`, `rate_be9`; `*_start` | `RLI6`, `RLI7`, `RBE9`; `RLI60`, `RLI70`, `RBE90` | core/burn_step.f90 (`store_li_be_rates`) | Light-element burning rates at end and start of step; `burn_li_be` multiplies them by the step in seconds. |
| — | — | `star%cross_section_scale(17)`; `qs0e_scale(8)`, `qqs0ee_scale(8)`, `o16_gamma_scale`, `c12_alpha_scale` | `SSTANDARD`; `QS0E`, `QQS0EE`, `FO16`, `FC12` (`/CROSS/`) | setup/map_user_inputs.f90, core/monte_carlo.f90 (`apply_monte_carlo_parameters`) | Multipliers on the S-factors (Monte-Carlo runs overwrite 1–3 and 16), their first- and second-derivative scales, and the ¹⁶O(α,γ)/¹²C(α,γ) multipliers. |

### 2.6 Envelope, atmosphere and surface fit

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| — | mixed | `star%envelope_fit_coeffs(9)` | `CFENV` | core/surfbc.f90:215-223 | Linear fits in (logP, logT) at the fitting point: (1–3) logR, (4–6) log L, (7–9) log T_eff. |
| log₁₀ T_eff (vertices) | dex (K) | `star%trial_log_temperature(3)`; also `tri_teffl` | `TRIT` | core/surfbc.f90 | log T_eff of the three envelope-triangle vertices. Not an interior temperature (trap 16). |
| log₁₀ L (vertices) | dex (L☉) | `star%trial_log_luminosity(3)`; also `tri_logl` | `TRIL` | core/surfbc.f90 | log L of the three vertices. |
| log₁₀ P, T, r at fit point | dex | `star%fit_point_pressure(3)`, `fit_point_temperature(3)`, `fit_point_radius(3)`; also `vtx_logp`, `vtx_logt`, `vtx_logr` | `PS`, `TS`, `RS` | core/envint_kernel.f90 | Envelope solution at the fitting mass for each vertex. log₁₀ values despite the names (trap 16). |
| — | — | `star%trial_sign_flag`; also `tri_orientation` | `FTRI` | core/surfbc.f90 | Orientation (sign) of the envelope triangle; a real number, not a logical. |
| — | dex | `star%stored_envelope_state(4)` | `ESTORE` | core/envint_kernel.f90:702-705 | (logP, logT, logR, log m/M) where the envelope integration first reaches T ≥ `tenv`. |
| — | — | `star%istore_flag`; also `stored_vertex_index` | `ISTORE` | core/envint_kernel.f90 | Vertex whose state is in `stored_envelope_state`. |
| — | logical | `star%recompute_envelope_triangle`; also `start_new_triangle` | `LNEW` | core/evolve_step.f90, core/surfbc.f90 | Build a fresh triangle around (T_eff, L). |
| — | logical | `star%reset_triangle` | `LRESET` | core/surfbc.f90, setup/rezone.f90 | Redo all three envelopes and retriangulate if needed. |
| log₁₀ T | dex (K) | `star%tenv` | `TENV` (`/CTLIM/`) | io/read_controls.f90 | ½(`tenv0` + `tenv1`); T threshold for storing envelope state. |
| log₁₀ T(τ=2/3) | dex (K) | `star%atm_hras` | `HRAS` (`/ATMOS/`) | setup/setups.f90:124 | HSRA T–τ fit at τ = 2/3, the offset for `atm_choice` = 2. The `star_info` comment says Krishna-Swamy; the code uses `hsra_t_tau_offset`. |
| — | logical | `star%use_ttau_relation` | `LTTAU` | core/surfbc.f90, core/envint_lib.f90 | Runtime switch to a T–τ relation, e.g. the gray fallback outside table range. |
| — | mixed | `env_struct%env_log10_*`, `env_gradients(3,:)`, ... | `/ENVSTRUCT/`-style arrays | core/envint_kernel.f90 | Envelope points (fitting point to photosphere). Gradient rows are (rad, **ad**, actual) (trap 19). |
| — | mixed | `atmo_struct%atmo_log10_*`, `atmo_gradients(3,:)`, ... | `/ATMSTRUCT/`-style arrays | core/envint_kernel.f90 | Atmosphere points. Gradient rows are (rad, **actual**, ad) (trap 19). |

### 2.7 Rotation state (per zone unless noted)

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| Ω | rad s⁻¹ | `star%omega` | `OMEGA` | rotation/omega_from_j.f90, rotation/solid_body_omega.f90 | Angular velocity; solid body within each convective region. |
| j | cm² s⁻¹ | `star%j_rot` | `HJM` → `specific_angular_momentum` | rotation/evolve_angular_momentum.f90, rotation/solid_body_omega.f90 | Specific angular momentum I Ω/Δm. |
| I | g cm² (per point) | `star%i_rot` | `HI` → `moment_of_inertia` | rotation/shape/zone_moments_of_inertia.f90 | Moment of inertia of the point's shell. |
| dI/dΩ | g cm² s | `star%qiw`; also `di_domega` | `QIW` | rotation/shape/zone_moments_of_inertia.f90 | Used in the Newton iteration for Ω(j). |
| r₀ | cm | `star%mean_radius`; also `r0` | `R0` | rotation/shape/shape.f90 | Characteristic radius of the distorted isobar, r = r₀(1 − A P₂). |
| η₂ | — | `star%eta_squared`; also `eta2` | `ETA2` | rotation/shape/shape.f90 | Radau distortion parameter from Radau's equation (not η squared). |
| ⟨g⟩ | cm s⁻² | `star%mean_gravity`; also `hg` | `HG` | rotation/shape/rotation_shape_factors.f90:111 | Effective gravity averaged over the isobar. |
| f_P | — | `star%fp_rot` | `FP` → `pressure_rotation_factor` | rotation/shape/rotation_shape_factors.f90 | Rotational correction in hydrostatic equilibrium (Endal & Sofia). |
| f_T | — | `star%ft_rot` | `FT` → `temperature_rotation_factor` | rotation/shape/rotation_shape_factors.f90 | Rotational correction in radiative transfer. |
| — | logical | `star%am_transport_convective_flag` | `LCZ` | setup/rezone.f90, core/read_starting_model.f90 (`am_convective_regions`) | T where a zone is convective for AM purposes (overshoot included). Read by `henyey_iterate`'s Ω(j) inversion. Never assigned before 2026-10 (inherited from F77; Appendix B #3). |
| Ω_old etc. | as base | `star%old_omega`, `old_specific_angular_momentum`, `old_moment_of_inertia`, `old_hg`, `old_mean_radius`, `old_eta_squared` | `WOLD`, `HJX`, `HIO`, `HGO`, `R0X`, `ETA2X` (`/OLDROT/`) | setup/rezone.f90:789-797 | Previous model's rotation state on the new grid. |
| j at step start | cm² s⁻¹ | `star%orig_specific_angular_momentum` | `HJM0` (`/ORIGSTART/`) | core/evolve_step.f90:418 | Restored on structure⇄rotation passes after the first. |
| J | g cm² s⁻¹ | `star%total_angular_momentum` | `SJTOT` | setup/rezone.f90, core/read_starting_model.f90 | Total angular momentum. |
| E_rot,tot | erg | `star%total_rotational_ke` | `SKEROT` | setup/rezone.f90, core/read_starting_model.f90 | Total rotational kinetic energy. |
| max ΔΩ/Ω | — | `star%max_domega_frac` | `DWMAX` | util/timestep_limit_omega.f90 | Largest fractional change of Ω since the previous model; drives the rotation timestep. |
| — | Gyr | `star%disk_gate_age_gyr` | `SAGE` (`/DISK/`) | core/evolve_step.f90:457 | Age compared with `disk_locking_age_gyr` by `disk_locking_engaged()`. |
| Ω at sub-step start | rad s⁻¹ | `rot_scr%omega_substep_start` | `WMST` (`/OLDROT2/`) | rotation/evolve_angular_momentum.f90 | Ω at the start of the current diffusion sub-step; read at the surface for the wind torque. |
| mid-sub-step structure | as base | `rot_scr%log_density_mid`, `log_pressure_mid`, `log_temperature_mid`, `log_radius_mid`, `luminosity_lsun_mid` (linear), `omega_mid`, `hg_mid`, `eta_squared_mid`, `moment_of_inertia_mid`, `mean_radius_mid`, `qiw_mid`, `convective_flag_mid`, `am_transport_convective_flag_mid` | locals of `getw`/`midmod` | rotation/mid_timestep_model.f90 | Structure interpolated to mid-sub-step for `secular_transport`. |

### 2.8 Rotational transport (`secular_transport` pipeline)

Alternate names are those listed in the alias table at the top of `rotation/seculr/secular_transport.f90`.

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| U_ES | cm s⁻¹ | `star%es_circulation_velocity`; `circ_scr%es_circulation_velocity_prev` | `VES`; `VES0` (`/TEMP2/`) | rotation/seculr/circulation_velocities.f90 | Eddington–Sweet circulation velocity (current and previous iteration). |
| U_SS | cm s⁻¹ | `star%secular_shear_velocity`; `..._prev` | `VSS`; `VSS0` | rotation/seculr/circulation_velocities.f90 | Secular-shear velocity (set to the dynamical-shear value where that operates). |
| U_GSF | cm s⁻¹ | `star%gsf_circulation_velocity`; `..._prev` | `VGSF`; `VGSF0` | rotation/seculr/circulation_velocities.f90 | Goldreich–Schubert–Fricke velocity. |
| U_μ | cm s⁻¹ | `circ_scr%mu_gradient_velocity` | `VMU` | rotation/seculr/circulation_velocities.f90 | μ-gradient (inhibiting) velocity. |
| U | cm s⁻¹ | `diffusion_velocity` (`secular_transport` argument); also `total_circulation_velocity` | `HV` | rotation/seculr/circulation_velocities.f90 | Total circulation velocity at interfaces; transport regions are where U > 0. |
| D_Ω | cm² s⁻¹ | `star%am_diffusion_coeff` (profile `D_omega`); local `am_diffusion_coeff` | `COD` | rotation/seculr/diffusion_velocity_scales.f90; copied at secular_transport.f90:295 | AM diffusion coefficient on the model grid; `star` keeps the last sub-step. |
| D_mix | cm² s⁻¹ | `star%mixing_diffusion_coeff` (profile `D_mix_rot`); argument `mixing_diffusion_coeff`; also `diffusion_coeff` | `COD2` | as above | Rotational mixing coefficient for composition. |
| f_c | — | `star%vfc` | `VFC` (`/VARFC/`) | rotation/seculr/zahn_coupling_factor.f90:112 | Zahn (1992) coupling factor min(1, ...), applied when `lvfc` is set. |
| (dΩ/dln r)_max | rad s⁻¹ | `rot_scr%max_domega_dr` (shell_physics); local `dynamical_shear_omega_limit` (secular_transport); also `dlnomega_dlnr_max` | `QWRMAX` (`/DWMAX/`); `QWRMX0` → `max_domega_dr_old` | core/shell_physics.f90:182, setup/rotation_stability_setup.f90:201 | Largest absolute dΩ/dln r stable against dynamical shear (Ri = ¼). Absolute and per ln r, despite `dr`/`dln` in the names (trap 3). |
| dΩ/dln r | rad s⁻¹ | `dlnomega_dlnr` (local) | `QWLNR` | rotation/seculr/secular_transport.f90:163 | **Absolute** dΩ/dln r at interfaces (trap 3). |
| r | cm | `radius_unlogged`; also `radius`, `radius_mid_prev` | `HRU` | setup/rotation_stability_setup.f90 | Unlogged radius at zone centres. |
| equal-grid set | as base | `eq_mixing_diffusion_coeff`, `eq_angular_momentum`, `eq_am_diffusion_coeff`, `eq_delta_angular_momentum`, `eq_mass`, `eq_moment_of_inertia`, `eq_omega`; aliases `equally_spaced_diffusion_coeff`, `angular_momentum`, `diffusion_coeff`, `delta_angular_momentum`/`dj`, `equally_spaced_mass` | `ECOD2`, `EJ`, `ECOD`, `DJ`, `EM`, `EI`, `EW` | rotation/seculr/am_transport_grid.f90, equal_grid_to_model.f90 | The same quantities on the equally spaced χ grid. |
| m on χ grid | g | `rot_scr%eq_enclosed_mass` | `ES1` (`/EGRID/`) | rotation/equal_spaced_grid.f90 | Enclosed mass at the χ points. |
| P, T, ρ, m, r at interface | cgs | `rot_scr%pm`, `tm`, **`dm`**, `hs3`, `interface_radius` | `PM`, `TM`, `DM`, `HS3`, `RM` (`/INTVAR/`) | setup/rotation_stability_setup.f90:57, 311-325 | Interface (cell-edge) values. `rot_scr%dm` is a **density**, not a mass (trap 8). |
| L, ε, g at interface | erg s⁻¹; erg g⁻¹ s⁻¹; cgs | `rot_scr%interface_luminosity`, `epsilm`, `interface_gravity_factor` | `ALM`, `EPSILM`, `HGM` | setup/rotation_stability_setup.f90 | Interface luminosity (cgs), ε_nuc (from `esumm`) and local gravity. |
| ∇, ∇_ad, ∂lnρ/∂lnT at interface | — | `rot_scr%delmi`, `delami`, `qdtmi` | `DELMI`, `DELAMI`, `QDTMI` | setup/rotation_stability_setup.f90:91, 319-323 | `delmi` holds interpolated ∇_rad until line 91 replaces it with min(∇_rad, ∇_ad). |
| μ, K, ν, Ω at interface | —, cm² s⁻¹, cm² s⁻¹, rad s⁻¹ | `rot_scr%mean_molecular_weight_interface`, `thermal_diffusivity_interface`, `kinematic_viscosity_interface`, `omega_interface` | `AMUMI`, `THDIFMI`, `VISCMI`, `WM` (`/INTVR2/`) | setup/rotation_stability_setup.f90 | Interface transport properties. |
| κ_ρ, κ_T | — | `rot_scr%dlnkappa_dlnrho`, `dlnkappa_dlnt` | `SQOD`, `SQOT` (`/ROTDER/`) | core/henyey_coefficients.f90 | Logarithmic opacity derivatives (genuinely logarithmic). |
| — | — | `rot_scr%neutrino_loss_fraction` | `FRACNUC` | core/henyey_coefficients.f90:382 | (ε_nuc − ε_ν − ε_grav)/ε_nuc with the stored signs; not a neutrino fraction (trap 14). |

### 2.9 Mixing, convection and zone bookkeeping

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| α_MLT | — | `star%mixing_length_alpha` | `CMIXL` (`/CONST3/`) | core/run_yrec.f90:270, core/check_solar_calibration.f90 | Working mixing length in pressure scale heights, copied from `job%mixing_length_array(nk)` each kind card. |
| — | zone index | `star%mixed_zone_bounds(12,2)` | `MXZONE` | mixing/find_convection_zones.f90 (called from mixing/mix.f90, core/evolve_step.f90:496) | First and last zone of each mixed (convective plus overshoot) region. |
| — | zone index | `star%mixed_zone_bounds_no_overshoot(12,2)` | `MXZON0` | mixing/mix.f90 | Same without overshoot. |
| — | zone index | `star%radiative_zone_bounds(13,2)` | `MRZONE` | mixing/find_convection_zones.f90 (called from core/evolve_step.f90:496) | First and last zone of each radiative region. |
| — | zone index | `star%core_cz_top_index` | `JCORE` | setup/rezone.f90:558, mixing/mix.f90 | Top zone of the convective core (1 if none). |
| — | zone index | `star%envelope_cz_bottom_index` | `JENV` | setup/rezone.f90:567, mixing/mix.f90 | Bottom zone of the surface convection zone (`nz` if none). |
| — | zone index | `star%iov1`, `star%iov2`, `star%iovim` | `IOV1`, `IOV2`, `IOVIM` (`/DPMIX/`) | mixing/overshoot_boundaries.f90 | Overshoot zone range, and the zone index passed to `temperature_gradients` for the adiabatic-overshoot test. |
| α_ov H_P | cm | `star%pressure_scale_height_start`, `pressure_scale_height_end` | `PSCAL0`, `PSCAL1` (`/LIOV/`) | core/evolve_step.f90:537-545, rotation/evolve_angular_momentum.f90:343-348 | Envelope **overshoot distance** below the CZ base, not H_P itself (trap 13). |
| r_CZ,prev | cm | `star%cz_base_radius_prev` | `RBASE0` (`/PREVCZ/`) | core/burn_step.f90:318, 336 | CZ-base radius at the end of the previous step, for light-element burning. |
| — | zone index | `star%envelope_cz_base_zone_prev`, `star%jcz` | `JENVOS`, `JCZ` | core/burn_step.f90:355, core/evolve_step.f90:280 | CZ-base zone of the previous step; CZ-base zone used by deuterium burning. |
| τ_cz | s | `star%convective_turnover_timescale`, `..._old` | `TAUCZ`, `TAUCZ0` (`/OVRTRN/`) | core/turnover_timescale.f90 | Convective turnover time, H_P/v_conv (or CZ width/v_conv) at the CZ base. |
| — | — | `star%alfmlt`, `phmlt`, `cmxmlt`; `valfmlt`, `vphmlt`, `vcmxmlt` | `ALFMLT`, `PHMLT`, `CMXMLT`, `VALFMLT`, ... (`/PUALPHA/`) | mixing/temperature_gradients.f90 | MLT pulse-output scalars and vectors; currently always zero (henyey_coefficients.f90 comment). |

### 2.10 Timestep and driver state

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| Δt | s | `star%dt` | `DELTS` → `delta_time` | core/evolve_step.f90, util/compute_timestep.f90 | Timestep. Negative means "do not age" (section 1.4). |
| Δt (saved) | s | `star%dt_saved` | `DELTS0` → `delta_time_saved` | core/evolve_step.f90:142-160 | dt before the end-age adjustment, written to the model when stopping at a fixed age. |
| Δt (chosen) | s | `star%chosen_dt` | `DELTSH` → `hydrogen_dt` | util/compute_timestep.f90 | Timestep `compute_timestep` adopted (minimum over all limiters). Reused as \|dt\| in `read_starting_model` (trap 20). |
| Δt | yr | `star%timestep_yr` | `DDAGE` | core/evolve_step.f90:261 | dt in years. |
| t | Gyr | `star%dage` | `DAGE` (`/THEAGE/`) | core/evolve_step.f90:262 | Stellar age. |
| — | — | `star%model_number` | `MODEL` | core/evolve_step.f90:162 | Model counter (not advanced on rescale-only cards). |
| — | — | `star%job%nk` | `NK` (`/ZRAMP/`) | core/run_yrec.f90 | Index of the current kind card (run). |
| — | logical | `star%model_diverged_flag` | `LARGE` | core/henyey_iterate.f90 (`corrections_too_large`) | Newton corrections exceeded `htoler(:,2)`. No timestep is cut: `reload_model_if_diverged` calls `read_starting_model`, which stops the run with ierr = 1 (`core/evolve_step.f90:94-99`, `:200-206`). |
| — | zone index | `star%h_shell_zone_begin`, `h_shell_midpoint_zone`, `h_shell_end_index` | `JXBEG`, `JXMID`, `JXEND` | setup/locate_shell_boundaries.f90 | Bottom, middle and top zones of the H-burning shell. |
| — | logical | `star%has_h_shell` | `LSHELL` | setup/locate_shell_boundaries.f90 | A hydrogen-burning shell exists. |
| — | logical | `star%solar_calibration_active`; `delta_x`, `delta_alpha`; `log_l_prev`, `log_r_prev` | `LSOL`; `DX`, `DA`; `BLP`, `RLP` (`/CALSUN/`) | core/check_solar_calibration.f90 | Solar calibration in progress; the next X and α corrections; the previous run's surface log L and log R. |
| ∂logL/∂X, ∂logR/∂X, ∂logL/∂α, ∂logR/∂α | dex | `star%dlum_dx`, `drad_dx`, `dlum_dalpha`, `drad_dalpha` | `DLDX`, `DRDX`, `DLDA`, `DRDA` | core/check_solar_calibration.f90:49-51 | Calibration Jacobian (empirical constants). |

### 2.11 Diagnostics and observables (filled after convergence by `core/observables_lib.f90`)

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| log₁₀ R | dex (R☉) | `star%log_R_surface` | — | core/observables_lib.f90:282-285 | Photospheric radius from L and T_eff. |
| log₁₀ g | dex (cgs) | `star%log_g_surface` | — | core/observables_lib.f90:283 | Surface gravity. |
| log₁₀ T_c, P_c, ρ_c | dex (cgs) | `star%central_log10_temperature`, `central_log10_pressure`, `central_log10_density` | `TCENTER`, `PCENTER`, `DCENTER` (`/CENT/`) | core/observables_lib.f90:190-192 | Centre values extrapolated from zone 1. |
| β_c, η_c | — | `star%central_beta`, `central_degeneracy_eta` | — | core/observables_lib.f90:193-194 | Central gas-pressure fraction and degeneracy. |
| M_cc | M☉ | `star%core_cz_mass` | — | core/observables_lib.f90:147 | Mass of the convective core. |
| M_CZ | M☉ | `star%envelope_mass` | `ENVLM` | core/observables_lib.f90:217 | Mass of the **surface convection zone** (above its base), not of the envelope above the fitting point (trap 21). |
| r_CZ | R☉ | `star%envelope_radius`; `star%envelope_cz_log_radius` (log₁₀ R☉) | `ENVR` | core/observables_lib.f90:219-221 | Radius of the surface-CZ base. |
| T, ρ, P, κ at CZ base | K, g cm⁻³, dyn cm⁻², cm² g⁻¹ | `star%envelope_cz_temperature`, `envelope_cz_density`, `envelope_cz_pressure`, `envelope_cz_opacity` | — | core/observables_lib.f90:222-232 | Conditions at the surface-CZ base (centre values if fully convective). |
| I_tot, I_CZ | g cm² | `star%total_moment_of_inertia`, `cz_moment_of_inertia` | — | core/observables_lib.f90:368-415 | Total moment of inertia, and that of the surface CZ. |
| P_rot | d | `star%rotation_period_days` | — | core/observables_lib.f90 | 2π/Ω(nz), capped at 9999. |
| v_eq | km s⁻¹ | `star%surf_velocity_kms` | — | core/observables_lib.f90 | Ω(nz) R. |
| ν_max | μHz | `star%nu_max` | — | core/observables_lib.f90:306 | Scaling-relation ν_max from (log g, T_eff). |
| Δν_ρ | μHz | `star%delta_nu_rho` | — | core/observables_lib.f90:310 | Scaling-relation Δν from the mean density. |
| Δν | μHz | `star%delta_nu` | — | core/observables_lib.f90:327 | Asymptotic large separation from the sound-travel integral. |
| ΔΠ₁ | s | `star%delta_Pg` | — | core/observables_lib.f90:355 | Asymptotic l=1 g-mode period spacing. |
| m_H-shell | M☉ | `star%h_shell_bot_mass`, `h_shell_mid_mass`, `h_shell_top_mass` | — | core/observables_lib.f90:442-447 | Masses at the H-shell boundaries. |
| r_H-shell/R | fraction | `star%h_shell_bot_radius`, `h_shell_mid_radius`, `h_shell_top_radius` | — | core/observables_lib.f90:448-456 | Radii of the H-shell boundaries as fractions of R. |

### 2.12 Pulse output

The per-zone `star%pulse_*` arrays are filled by `core/henyey_coefficients.f90` for every converged model. The column set `stx_pulse(35, n_ext)` is assembled over the full stitched model by `build_pulse_points` in `core/stitched_model.f90`. The FGONG, GYRE and GSM writers only read these arrays.

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| (∂lnρ/∂lnP) at fixed T | — | `star%pulse_dlnrho_dlnp` | `PQDP` (`/PULSE1/`) | core/henyey_coefficients.f90 | Logarithmic (correct). |
| (∂lnρ/∂lnT) at fixed P | — | `star%pulse_dlnrho_dlnt` | `PQDT` | core/henyey_coefficients.f90 | Logarithmic (correct); δ = −this. |
| dε/dlnρ | erg g⁻¹ s⁻¹ | `star%pulse_deps_dlnrho` | `PQED` | core/henyey_coefficients.f90 | **Absolute**; written to `eps_eps_rho` as is (trap 1). |
| dε/dlnT | erg g⁻¹ s⁻¹ | `star%pulse_deps_dlnt` | `PQET` | core/henyey_coefficients.f90 | **Absolute**; written to `eps_eps_T` as is (trap 1). |
| κ_ρ, κ_T | — | `star%pulse_dlnkap_dlnrho`, `pulse_dlnkap_dlnt` | `PQOD`, `PQOT` | core/henyey_coefficients.f90 | Logarithmic (correct). |
| c_P | erg g⁻¹ K⁻¹ | `star%pulse_specific_heat` | `PQCP` | core/henyey_coefficients.f90 | Specific heat. |
| μ | — | `star%pulse_mean_molecular_weight` | `PRMU` | core/henyey_coefficients.f90:412 | R/(R/μ) from the EOS; the legacy slot held R/μ. |
| 1/μ_e | — | `star%pulse_electron_mean_weight_inverse` | `PEMU` | core/henyey_coefficients.f90 | Inverse electron mean weight (already inverted). |
| pulse columns 1–18 | cgs | `stx_pulse(ipul_r ... ipul_omega, :)`: r, m, L (erg s⁻¹), P, T, ρ, ∇, N², Γ₁, ∇_ad, δ, κ, κ·κ_T, κ·κ_ρ, ε, ε·ε_T, ε·ε_ρ, Ω | — | core/stitched_model.f90:47-55 | GYRE schema-101 set. Columns 16–17 multiply `eps_total` by the already-absolute derivative (trap 1). |
| pulse columns 19–34 | mixed | `ipul_cp`, `ipul_mu_e_inv`, `ipul_h1`, `ipul_z`, species at `ipul_species_base+k` (k=1..11: ³He, ¹²C, ¹³C, ¹⁴N, ¹⁶O, ²H, ⁴He, ⁷Li, ¹⁵N, ¹⁷O, ¹⁸O), `ipul_eps_grav` | — | core/stitched_model.f90:700-706 | FGONG extras. Column 35 is never filled. |
| q* point scratch | mixed | `pt_scr%qqdp`, `qqdt`, `qqod`, `qqot`, `qdel`, `qdela`, `qqcp`, `qrmu`, `qemu`, `qd`, `qdl`, `qo`, `qol`, `qp`, `qpl`, `qt`, `qtl`, `current_*` | same names (`/PULSE2/`), `current_*` → `/ENVPRT/` | core/envelope_derivs.f90, core/atmosphere_derivs.f90, core/henyey_coefficients.f90 | Physics of the last evaluated point for the envelope/atmosphere recorder. `qp`, `qt`, `qtl` mean different things per writer (trap 2). |

### 2.13 Physical constants and unit scales

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| ln 10 | — | `ln10`, `clndp` | `CLN`, `CLNDP` | setup/setups.f90:43-44 | Natural log of 10 (two copies). |
| 1/ln 10 | — | `clni` | `CLNI` | setup/setups.f90:45 | |
| 4π, log₁₀ 4π, log₁₀(4π/3) | — | `c4pi`, `c4pil`, `c4pi3l` | `C4PI`, `C4PIL`, `C4PI3L` | setup/setups.f90 | |
| log₁₀ G | dex (cgs) | `cgl` | `CGL` | setup/setups.f90:78-80 | From control `G_cgs` if > 0, else −7.17571. |
| R_gas | erg mol⁻¹ K⁻¹ | `gas_constant` | `CGAS` | setup/setups.f90:73 | 8.314510×10⁷. |
| a/3 | erg cm⁻³ K⁻⁴ | `radiation_constant_over_3`, `ca3l` (log₁₀) | `CA3`, `CA3L` | setup/setups.f90:70-71 | Radiation constant over 3 (4σ/3c). |
| σ | erg cm⁻² s⁻¹ K⁻⁴ | `csig`, `csigl`; `sigma_sb_cgs_legacy` (5.669×10⁻⁵, viscos only) | `CSIG`, `CSIGL` | setup/setups.f90:67-68, state/phys_const_lib.f90 | Stefan–Boltzmann constant (two values coexist). |
| — | — | `cmixl2`, `cmixl3`, `cdelrl` | `CMIXL2`, `CMIXL3`, `CDELRL` | setup/setups.f90:91-93 | MLT constants (1/3, 16√2 σ) and the ∇_rad constant. Not powers of the mixing length. |
| yr | s | `seconds_per_year` | `CSECYR` | setup/setups.f90:57 | 3.1558×10⁷. |
| L☉ | erg s⁻¹ | `star%solar_luminosity_cgs`, `log10_solar_luminosity` | `CLSUN`, `CLSUNL` (`/CONST/`) | setup/setups.f90 | Solar luminosity (namelist `clsun`; rescaled per Monte-Carlo run). |
| ln10/L☉ | s erg⁻¹ | `star%ln10_over_lsun` | `CLNSUN` | setup/setups.f90:48 | ln 10 divided by L☉. **Not a logarithm** (trap 22). |
| M☉ | g | `star%solar_mass_cgs` (control `ctrl%solar_mass_cgs`), `log10_solar_mass` | `CMSUN`, `CMSUNL` | setup/setups.f90 | Solar mass. |
| R☉ | cm | `star%solar_radius_cgs`, `log10_solar_radius`; `rsun_cgs_legacy` (microdiff) | `CRSUN`, `CRSUNL` | setup/setups.f90 | Solar radius (namelist `crsun`). |

### 2.14 Key controls

- **Inlist name:** the new-style `&controls` name, or `&star_job` where marked (sj).
- **Legacy:** the F77 namelist name, read by the legacy reader (`io/read_controls.f90`).
- **Code field:** where the value lives at run time. `ctrl` means `star%ctrl%`. `job` means `star%job%`, which holds run-list and calibration state and may change during a run.
- **Full list:** the ~380 controls are in `defaults/controls_registry.tsv` and the generated `defaults/*.defaults`.

| Symbol | Units | Current name(s) | Legacy | Owner | Definition |
|---|---|---|---|---|---|
| α_MLT | — | inlist `mixing_length_alpha` (sj, per card); code `job%mixing_length_array(50)` | `CMIXLA` | job | Mixing length per kind card; copied to `star%mixing_length_alpha`. |
| α_ov,core | H_P | `overshoot_alpha_core`; code `ctrl%overshoot_alpha_core` | `ALPHAC` | ctrl | Core overshoot extent. |
| α_ov,env | H_P | `overshoot_alpha_envelope`; code `ctrl%overshoot_alpha_envelope` | `ALPHAE` | ctrl | Envelope (undershoot) extent below the surface CZ. |
| — | logical | `core_overshoot_active`, `envelope_overshoot_active`, `intermediate_overshoot_active`; code `job%core_overshoot_active`, `job%envelope_overshoot_active`, `job%lovstm`; `adiabatic_overshoot_active` → `ctrl%ladov` | `LOVSTC`, `LOVSTE`, `LOVSTM`; `LADOV` | job / ctrl | Enable overshoot at the core, at the surface-CZ base, and at intermediate CZs; treat overshoot regions as adiabatic. |
| — | logical | `use_semiconvection`; code `job%use_semiconvection` | `LSEMIC` | job | Semiconvection treatment (conflicts with core overshoot). |
| — | logical | `use_diffusion_z`, `diffuse_helium_active`, `diffuse_lithium`, `use_thoul_diffusion`; code `job%use_diffusion_z`, `job%diffuse_helium_active`, `ctrl%diffuse_lithium`, `ctrl%use_thoul_diffusion` | `LDIFZ`, `LDIFY`, `LDIFLI`, `LTHOUL` | job / ctrl | Gravitational settling of Z, of He and of Li/Be; use the Thoul et al. coefficients. |
| f_Y, f_Z | — | `settling_y_factor`, `settling_z_factor`; code `ctrl%fgry`, `job%fgrz` | `FGRY`, `FGRZ` | ctrl / job | Multipliers on the He and metal settling velocities. |
| — | logical | `rotation_flag`; code `job%rotation_active` | `LROT` | job | Rotation on (structure corrections, Ω, j). |
| — | logical | `instability_transport_active`; code `job%instability_transport_active` | `LINSTB` | job | Rotational AM and composition transport by instabilities. |
| Ω₀ | rad s⁻¹ | `initial_omega` (with `set_initial_omega`); code `ctrl%initial_omega`, `ctrl%set_initial_omega` | `WNEW`, `LWNEW` | ctrl | Initial solid-body Ω assigned to every zone (`read_starting_model.f90:355-359`). |
| — | — | `num_rotation_structure_iters`; code `job%num_rotation_structure_iters` | `ITDIF1` | job | Structure⇄rotation passes per step (the level-4 loop). |
| — | — | `max_diffusion_iters`; code `ctrl%max_diffusion_iters` | `ITDIF2` | ctrl | Iteration limit of the AM diffusion solve. |
| — | — | `convergence_tolerance`; code `ctrl%convergence_tolerance` | `DJOK` | ctrl | Relative Δj (and ΔΩ) tolerance of the AM diffusion iteration. |
| max ΔΩ/Ω per sub-step | — | `diffusion_timestep_factor`; code `job%diffusion_timestep_factor` | `DTDIF` | job | Largest fractional ΔΩ allowed per rotational sub-step; sets the sub-step count (`evolve_angular_momentum.f90:191`). Overwritten by `time_dw_mix` when `use_new_rate_variables` is set (trap 23). |
| f_ES, f_GSF, f_SS, f_μ, f_c | — | `es_velocity_scale`, `gsf_velocity_scale`, `secular_shear_velocity_scale`, `mu_gradient_scale`, `mixing_velocity_scale`; code `ctrl%` same names | `FES`, `FGSF`, `FSS`, `FMU`, `FC` (`/VMULT/`) | ctrl | Efficiency multipliers of the transport velocities; `fc` also caps the Zahn factor. |
| — | — | `fpft_average_weight`; code `job%acfpft` | `ACFPFT` | job | Used as the **relative tolerance** of the r₀/η₂ iterations (`shape.f90:65,112`), not as an averaging weight (trap 24). |
| — | — | `solver_tolerances`; code `ctrl%htoler(5,2)` | `HTOLER` | ctrl | Rows 1–4: per-variable convergence (col 1) and too-large (col 2) correction limits; row 5: residual test. |
| — | — | `max_iter_level1` ... `max_iter_level4`; code `ctrl%max_iter_level1..4` | `NITER1`..`NITER4` | ctrl | Newton iteration limits of the four solve levels (level 3 defaults to 2). |
| — | — | `solver_underrelax_start`, `solver_underrelax_increment`; code `ctrl%fcorr0`, `ctrl%fcorri` | `FCORR0`, `FCORRI` | ctrl | Start value of `fcorr` (\|fcorr0\| − fcorri) and its per-iteration growth; the sign of `fcorr0` enables the ramp. |
| — | — | `atm_choice`; code `job%atm_choice` | `KTTAU` | job | Surface BC: 0 Eddington gray, 1 Krishna-Swamy, 2 HSRA, 3–5 tabulated (Kurucz, Allard, Kurucz/Castelli). May be reset to 0 at run time. |
| — | — | `env_error_tol`, `atm_error_tol`; code `ctrl%env_error_tol`, `ctrl%atm_error_tol` | `ENVERR`, `ATMERR` | ctrl | Error tolerances of the envelope and atmosphere integrations. |
| Δlog T_eff, Δlog L | dex | `tri_delta_teffl`, `tri_delta_logl`; code `ctrl%tri_delta_teffl`, `ctrl%tri_delta_logl` | `TRIDT`, `TRIDL` | ctrl | Size of the envelope-fit triangle. |
| — | logical | `use_opal95_eos`, `use_opal2001_eos`, `use_opal2006_eos`, `use_mhd_eos`, `use_scv_eos`; code `ctrl%use_opal95_eos`, ... `ctrl%use_mhd_eos`; `use_scv_eos` lives in `eos/scv_eos_lib` | `LOPALE`, `LOPALE01`, `LOPALE06`, `LMHD`, `LSCV` | ctrl | EOS table selection on top of the Yale/Saha base. |
| — | logical | `use_opal95_tables` (sj), `use_alex06_tables` (sj), `use_conductive_opacity`; code `ctrl%` same names | `LOPAL95`, `LALEX06`, `LCONDOPACP` | ctrl | Opacity table selection. |
| — | logical | `use_extended_composition`; code `job%use_extended_composition` | `LEXCOM` | job | Evolve the 15-species layout (adds ²H, ⁶Li, ⁷Li, ⁹Be). |
| log₁₀ T | dex (K) | `nuclear_logT_cutoffs(5)`; code `ctrl%nuclear_logT_cutoffs` | `TCUT` | ctrl | Temperatures below which burning groups are off; (1) gates `compute_zone_energy` in `henyey_coefficients`. |
| — | logical | `helium_flash_active`; code `ctrl%helium_flash_active` | `LKUTHE` | ctrl | He-flash mode: no rezoning, no luminosity renormalization, adaptive `fcorr`. |
| — | — | `timestep_factors`; code `ctrl%atime(14)` (slots `itime_*` in `io/controls_lib.f90:113-126`); `use_structure_dt_limits` → `job%use_structure_dt_limits` | `ATIME`; `LPTIME` | ctrl / job | Timestep limiters; (13) is the largest dt growth factor, and (8–11) apply only while `use_structure_dt_limits` is on (switched off at run time after log T_c > 7.1). Overwritten from `time_*` controls when `use_new_rate_variables` is set (trap 23). |
| — | logical | `improved_first_guess_flag`; code `ctrl%improved_first_guess_flag` | `LNEWS` | ctrl | Extrapolate the next first guess from `*_entropy_term` rates. |
| — | logical | `use_new_rate_variables`; code `ctrl%lnewvars` | `LNEWVARS` | ctrl | Remaps `tol_*`, `flag_*` and `time_*` into the rezoning tolerances `ctrl%chi_grid_scale(12)` (legacy `HPTTOL`, slots `ichi_*`) and `atime`. Not a nuclear-rate switch (trap 23). |
| t_end | yr | `max_age` (sj, per card); code `job%target_end_age(50)` | `ENDAGE` | job | Stop age per kind card. |
| X₀, Z₀ | mass fraction | `initial_x_array`, `initial_z_array` (sj); code `job%initial_x_array(50)`, `job%initial_z_array(50)` | `XENV0A`, `ZENV0A` | job | Initial X and Z per kind card. |
| — | — | `run_kind` (sj); code `job%rescale_kind(50)` | `KINDRN` | job | 1 evolve, 2 rescale (no aging), 3 rescale and evolve (pre-MS). |
| — | logical; dex | `calibrate_solar_model` (sj), `luminosity_tolerance`, `radius_tolerance` (sj); code `ctrl%` same names | `LCALS`, `TOLL`, `TOLR` | ctrl | Automatic solar calibration and its L and R tolerances. |
| — | logical; Gyr; rad s⁻¹ | `disk_locking_active`, `disk_locking_age_gyr`, `disk_omega_rad_s`; code `job%` same names | `LDISK`, `TDISK`, `PDISK` | job | Disk locking: surface CZ held at `disk_omega_rad_s` while age ≤ `disk_locking_age_gyr`. Old deck comments give wrong units. |
| Ω_sat | rad s⁻¹ or Ω_sat·τ_cz | `wind_saturation_omega`; code `job%wind_saturation_omega` | `WMAX` | job | Wind saturation threshold. Its meaning depends on its magnitude (trap 25). |
| η_R | M☉ yr⁻¹ per (L R/M) in solar units | `use_reimers_wind`, `reimers_scaling_factor`; code `ctrl%use_reimers_wind`, `ctrl%reimers_scaling_factor` | `LREIMER`, `CREIM` | ctrl | Reimers mass loss Ṁ = η_R (L/L☉)(R/R☉)/(M/M☉) (`wind/massloss.f90:121-128`). |

---

## 3. Known naming traps

Every item was confirmed by reading the code. Lines are for the 2026-10-10 working tree.

1. **`dln` names that held ABSOLUTE ε derivatives (renamed 2026-10).**
   - `compute_zone_energy` (formerly `engeb`) returns `burn%deps_dlnrho` and `burn%deps_dlnt` as absolute dε/d ln ρ and dε/d ln T in erg g⁻¹ s⁻¹ (a `nuclear_energy_result`, computed in `net/nuclear_energy.f90`).
   - Three copies carried logarithmic-looking names until 2026-10: `zone_dlnepsilon_*` (now `zone_deps_*`, `core/henyey_coefficients.f90`), `star%pulse_dlneps_*` (now `star%pulse_deps_*`), `rot_scr%dlnepsilon_*` (now `rot_scr%deps_*`).
   - The Jacobian always used them correctly, and so does the rotation code (`setup/rotation_stability_setup.f90`): there f·dε/dlnT is added to ε·(1 − f − χ_T), which is consistent only for the absolute derivative.
   - Misled by the old names, commit 2ed398ca (2026-08-30) multiplied them by ε a second time in both pulse writers. The `eps_eps_T`/`eps_eps_rho` columns were wrong by a factor of ε until the 2026-10 fix (`TARGET_ARCHITECTURE.md` Appendix B #1).
   - The `engeb` comment used to say "the callers divide by eps themselves". None does; the comment was corrected in 2026-10.
2. **`pt_scr%qp`, `qt`, `qtl` mean different things depending on who wrote them last.**
   - In `henyey_coefficients` they are Henyey terms: q_P, q_T = ∇ q_P, and q_T/(L ln10) (`core/henyey_coefficients.f90:189-205`).
   - In the envelope and atmosphere integrands they are P (dyn cm⁻²), T (K) and log₁₀ T (`core/envelope_derivs.f90:127-130`, `core/atmosphere_derivs.f90:84-87`).
   - `qrmu` is R/μ, not μ.
3. **`dlnomega_dlnr` is the absolute dΩ/d ln r, and `max_domega_dr` is per ln r.**
   - `secular_transport.f90:163` forms ¼ΔΩ·(Σ 1/Δln r), which has units of rad s⁻¹.
   - `rot_scr%max_domega_dr`, the local `dynamical_shear_omega_limit` and `dlnomega_dlnr_max` (all legacy `QWRMAX`) are the limiting |dΩ/d ln r|, not dΩ/dr (`core/shell_physics.f90:182`, `setup/rotation_stability_setup.f90:198-203`).
4. **`eps_channels(1:5,:)` are fractions, but slots 6 and 7 are rates.**
   - The `star_info` comment says "[erg/g/s]" (`state/star_info_lib.f90:89`).
   - `henyey_coefficients` stores `energy_gen_component(j)/eps_total` in slots 1–5 (`:346-354`), and erg g⁻¹ s⁻¹ in slots 6 (ν) and 7 (grav).
   - The profile columns `eps_ppI`…`eps_he3` (`core/stitched_model.f90:483-487`) are therefore fractions of `eps_nuc`, while `eps_neu` and `eps_grav` are rates.
5. **`i_eps_he3` (slot 5) and the profile column `eps_he3` hold the triple-α share.**
   - Slot 5 is filled from `burn%triple_alpha_energy_gen` (`core/henyey_coefficients.f90:227`).
   - The luminosity twin is correctly named `i_lum_3alpha`.
6. **`old_shell_mass` is not a mass and is unused.**
   - The `star_info` comment (`state/star_info_lib.f90:339-342`) says it holds log-mass coordinates of the pre-rezoning grid. `setup/rezone.f90:67-74` says `HSS` was the *new* log-mass grid, and that this use moved to locals.
   - Nothing reads or writes `star%old_shell_mass` today.
   - `converged_zone` (`LOCONS`) is likewise unused.
7. **`gravitational_luminosity` is a specific rate (erg g⁻¹ s⁻¹), not a luminosity.**
   - It is ε_grav (`core/henyey_coefficients.f90:289-290`), legacy `HHC`.
   - The luminosity is `luminosity_breakdown(i_lum_grav)`.
8. **`rot_scr%dm` is the interface DENSITY, not a shell mass.**
   - `setup/rotation_stability_setup.f90:315` interpolates log ρ to the interface and unlogs it.
   - `star%dm` is the point mass weight in g.
   - Its siblings `rot_scr%pm`, `tm` and `hs3` are P, T and m at the interface.
9. **`xa_start` has two meanings.**
   - It is the start-of-step composition (`core/evolve_step.f90:245`). Burning starts from it, and Newton levels > 2 reset `xa` to it (`core/henyey_iterate.f90:136`).
   - Inside the rotation step it is also the **sub-step checkpoint**. It is restored on a redo (`rotation/evolve_angular_momentum.f90:310`) and overwritten with the current composition after each accepted sub-step (`:324`, `:388`, `:394`).
   - `orig_composition` exists only to undo that overwriting between structure⇄rotation passes (`evolve_step.f90:420`, `:430`).
10. **`max_residual(4)` holds three quantities within one Newton iteration.**
    - Zeroed at `core/henyey_iterate.f90:155-158`.
    - Max |residual| per equation, accumulated in `henyey_eliminate.f90:90-99` and tested against `htoler(5,*)` (`henyey_iterate.f90:186-189`).
    - Max |correction| after `henyey_solve`, where the L correction is the smaller of absolute and relative (`:216-241`). This is tested against `htoler(1:4,*)`.
    - The **signed** correction at the location of the maximum, for the log line (`:255-258`).
11. **`elim_rhs` becomes the corrections.** `henyey_solve` overwrites the residuals in place with ΔlogP, ΔlogT, ΔlogR and ΔL (`core/henyey_solve.f90`), which `henyey_iterate.f90:294-304` then applies.
12. **`*_entropy_term` arrays are rates of change, not entropy.**
    - `temperature_entropy_term` = Δlog T/Δt, `pressure_entropy_term` = Δlog P/Δt, `radius_entropy_term` = Δlog R/Δt, and `luminosity_entropy_term` ≈ Δln L/Δt (`core/evolve_step.f90:388-391`).
    - They exist only to extrapolate the next first guess.
    - The entropy term itself is ε_grav, formed from `log_*_delta`.
13. **`pressure_scale_height_start/end` hold α_ov·H_P, and the two writers gate it differently.**
    - The value is the envelope overshoot distance, and zero when overshoot is off (`core/evolve_step.f90:537-545`).
    - `evolve_step` gates it on `envelope_overshoot_active` (`LOVSTE`). `rotation/evolve_angular_momentum.f90:343-348` gates it on `job%lovstm` (intermediate overshoot), which looks like a bug.
14. **`neutrino_loss_fraction` is not a neutrino fraction.**
    - It is (ε_nuc − ε_ν − ε_grav)/ε_nuc with ε_ν stored negative (`core/henyey_coefficients.f90:378-385`), legacy `FRACNUC`.
15. **`luminosity_breakdown` is renormalized every Newton iteration.**
    - `henyey_iterate.f90:173-182` scales all eight channels so that they sum to L(nz).
    - Unknown heating therefore leaks into the nuclear channels. This is the trap named in `TARGET_ARCHITECTURE.md` 4.5a.
16. **The envelope-triangle arrays are named after the wrong quantity or omit "log".**
    - `trial_log_temperature` holds log **T_eff** of the vertices (`surfbc`'s own name is `tri_teffl`).
    - `fit_point_pressure/temperature/radius` hold **log₁₀** P, T, r (`vtx_logp/logt/logr`).
17. **`stotal` duplicates `log_total_mass`.**
    - Both are log₁₀ M in g.
    - Mass loss updates the two separately (`wind/mdot.f90:249-250`). The envelope code reads `stotal`, and the structure code reads `log_total_mass`.
18. **Several quantities have two copies with different fill times or definitions.**
    - c_P: `scp` (Henyey pass) and `cp` (`shell_physics`).
    - μ: `star%mu` (full-ionization estimate, `core/shell_physics.f90:125`), `pulse_mean_molecular_weight` (R/(R/μ) from the EOS, `henyey_coefficients.f90:412`), and the EOS slot `i_gas_constant`, which is R/μ.
19. **The 3-row gradient carriers use two row orders, and the named constants document the wrong one.**
    - `i_grad_rad=1, i_grad_actual=2, i_grad_ad=3` (`state/star_info_lib.f90:83-87`) give the order (rad, actual, ad), the order of the former `DEL_GRAD`. Their comment says they apply to `current_gradients` and `env_gradients`, but no code uses them.
    - `pt_scr%current_gradients` and `env_struct%env_gradients` are in fact (rad, **ad**, actual) (`core/envelope_derivs.f90:108-110`).
    - `atmo_struct%atmo_gradients` is (rad, **actual**, ad) (`core/envint_kernel.f90:459-461`).
    - The output code has already mixed up the two orders once (`core/stitched_model.f90:598-607`).
20. **`chosen_dt` is reused transiently.**
    - It is the timestep adopted by `compute_timestep` (legacy `DELTSH`, later named `hydrogen_dt` although it is not a hydrogen-burning governor).
    - `read_starting_model` is passed `star%chosen_dt` as `delta_time_abs` and sets it to |dt| (`core/read_starting_model.f90:311`; call sites `core/evolve_step.f90:206`, `core/run_yrec.f90:280`).
21. **`envelope_mass` and `envelope_radius` refer to the surface convection zone.**
    - They are the mass above the surface-CZ base in M☉, and that base's radius in R☉ (`core/observables_lib.f90:217-221`).
    - They do not describe the envelope above the fitting point.
    - `envelope_cz_log_radius` is log₁₀(r/R☉).
22. **`ln10_over_lsun` is not a log.**
    - It is ln 10/L☉ in s erg⁻¹ (`setup/setups.f90:48`; comment at `state/star_info_lib.f90:264-267`).
    - Multiplied by m, it gives `cccql` (`core/henyey_coefficients.f90:308`).
23. **`use_new_rate_variables` (`LNEWVARS`) has nothing to do with reaction rates.**
    - It overwrites the rezoning tolerances `chi_grid_scale` (legacy `HPTTOL`) and the timestep limiters `atime` from the `tol_*`, `flag_*` and `time_*` controls (`setup/map_user_inputs.f90:196-258`).
    - It also overwrites `diffusion_timestep_factor`, `max_domega_global` and `use_structure_dt_limits`.
    - The registry description "use the new nuclear-rate variable set" is wrong.
24. **`fpft_average_weight` (`acfpft`) is a convergence tolerance.**
    - The registry describes it as an averaging weight. The code uses it as the relative tolerance of the r₀ Newton iteration and, as its cube root, of the η₂ iteration (`rotation/shape/shape.f90:65`, `:98`, `:112`).
    - `check_angular_momentum.f90:168-178` loosens and restores it at run time, which is why it lives in `job`.
25. **`wind_saturation_omega` (`WMAX`) changes meaning with its value.**
    - With instability transport on and a value > 1, it is read as Ω_sat·τ_cz(☉), and Ω_sat = value/τ_cz.
    - Otherwise it is Ω_sat in rad s⁻¹ (`wind/kawaler_wind.f90:51-62`).
26. **The structure `_start` arrays are written by `rezone`, not at a single "start of step" point.**
    - `setup/rezone.f90` stores `logP_start`, `luminosity_lsun_start`, `log_Teff_start` and `nz_start` while regridding (`:345-355`), and `convective_flag_start` after regridding.
    - `xa_start` is stored later, after mass loss.
    - What all of these have in common is "previous converged model, remapped to the current grid".

---

## 4. Not confirmed

**Marked "unconfirmed" in the tables.** A reviewer who knows the answer should fill these in:

- the slot order of `fxion_zone(3,:)` and `eos_res(i_fxion:i_fxion+2)`;
- the physical units of `neutrino_flux`, `neutrino_flux_total` and `neutrino_flux_zone` (the code comments say "10¹⁰ ergs per gm");

**Left out of the tables** because their meaning or units could not be settled from the code:

- `star%fracstep` (`FRACSTEP`, initialized to 0.5);
- the scalar `star%convective_velocity` (`VEL`);
- `star%log_rate_li6/li7/be9_prev` (`BLI6OLD`, ...; log base unknown);
- `rot_scr%quadrupole_moment` and `local_gravity` (`QUAD`, `GG`).

**Discrepancy:** the value and comment of `dlum_dalpha` differ by a factor of 10 (`core/check_solar_calibration.f90:51`: value 0.0139, comment "+0.139").

### Network (2026-10)

| Symbol | Units | Current name | Legacy name | Owning file | Definition |
|---|---|---|---|---|---|
| R_j | per Gyr per amu | `reaction_rate(j, ...)` in network_lib | — | net/network_lib.f90 | Rate of table reaction j: b_j g_j prod X(law_j). |
| g_j | per Gyr per amu per (X_i X_k) | `rate_vec(rr_*)`, `rate_by_zone(rr_*, k)` | `HR1..HR13` | net/net_lib.f90 (`network_rate_vector`) | Composition rate rows; abundance factors omitted. |
| λ, dlnλ/dlnT, dlnλ/dlnρ | g⁻¹ s⁻¹ per (X_i X_k); — | `nuclear_rates%rate`, `%dlnrate_dlnt`, `%dlnrate_dlnrho` | `RATE`, `PET`-terms | net/net_lib.f90 (`nuclear_reaction_rates`) | The one per-zone rate evaluation (screened, scaled), with logarithmic derivatives. |
| f(N15→C12), f(Be7 e-capture) | — | `nuclear_rates%frac_n15_pa`, `%frac_be7_ecap` | `F3`, `F1` | net/net_lib.f90 | Branching fractions. |

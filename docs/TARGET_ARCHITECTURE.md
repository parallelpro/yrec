# YREC target architecture

**Status:** draft for discussion in the YREC group, 2026-10-10.
**Applies to:** branch `yrec-modern`.
**Companion documents:** `docs/HOW_YREC_WORKS.md` (physics and numerics, currently using legacy routine names), `docs/MODERNIZATION_CHANGES.md` (what the 2025–26 modernization changed), `docs/ROADMAP.md` (open work items).

---

## Summary

The 2026 modernization changed how YREC is *written*: free-form source, no COMMON blocks, no gotos, readable names, MESA-style inlists and output, a regression net. It changed much less of how YREC is *organized*. The COMMON blocks became one global variable, `star`, that 145 of 229 source files import, with about 265 state fields and 399 controls. Physics is still evaluated inside the solver loops. Most conventions are still undocumented.

As a result, adding physics remains hard for anyone who does not already know the code (section 1). This document proposes a target organization that a newcomer can extend, and a path to reach it without stopping the code from producing science.

The proposal **does not change the method**. Henyey relaxation with an envelope fitting point, operator-split composition, the four iteration levels and shellular rotation with Endal–Sofia shape factors all stay. It changes three things:

1. **State is passed explicitly.** Each stage of a timestep declares what it reads and writes.
2. **Physics is kept out of the numerics.** Solvers see residuals, Jacobians and diffusion coefficients, never table lookups.
3. **Each physics domain has one extension point.** Adding an energy source, a reaction, a table source or a transport mechanism means writing one module and registering it.

The work proceeds one subsystem at a time inside the running code, in the order network → structure → microphysics → rotation. Each replacement is checked against the code it replaces and against independent physics references.

---

## 1. Why the organization needs to change

### 1.1 What a newcomer faces today

In October 2026, four common tasks were traced through the code from start to finish. Difficulty is rated 1 (easy) to 5 (very hard) for someone new to the code.

| Task | Difficulty | Files touched | Main obstacle |
|---|---|---|---|
| Add an extra energy term (interior) | 3; 5 if it depends on r or L, or sits above the fitting point | ~5 plus generated files | Derivative conventions are undocumented. A luminosity renormalization absorbs the new term into the nuclear channels. The envelope integration drops it silently. |
| Add an opacity source | 3 | ~18 | Threading through init signatures, table paths and controls. The low-T table precedence is coded in three places. |
| Add an EOS option | 4 | ~16 | Two incompatible integration patterns; the 24 output slots have no written contract. |
| Update a reaction rate | 3–4 | both rate copies plus Monte Carlo | Every rate is coded twice (`rates` in `net_lib`, `engeb` in `burn_lib`), with derivatives worked out by hand. |
| Track a new isotope | 5 | 34+ | The species count 15 is a literal in 45 declarations. The composition solve has a hand-coded 7×7 Jacobian. |
| Add an angular-momentum transport process | 4 | ~12 | A process can only enter as a "velocity". Zones where the summed velocity is zero are skipped, so a new process tested on its own silently does nothing. |

### 1.2 Four root causes

1. **Nowhere to plug new physics in.** Each physics choice is a hard-coded if-chain. Adding an option means editing every place the existing options pass through.
2. **Physics runs inside the numerics.**
   - `henyey_coefficients` (legacy `coefft`) calls `eos_get`, `kap_get` and `engeb` inside the loop that feeds the eliminator. It also writes about 25 diagnostic arrays on every Newton iteration.
   - `henyey_iterate` (legacy `crrect`) runs envelope integrations, composition mixing (at levels ≥ 3) and the Ω(j) inversion inside the Newton driver.
   - In rotation, `secular_transport` (legacy `seculr`) mixes physics, wind torques, grid mapping, linear solves and convergence logic in one loop. The convergence checks halve the timestep and change controls as a side effect.
3. **Contracts that live only in people's heads.**
   - The structure solver's own header says its partial derivatives are "inferred from their use, not from surviving documentation".
   - Several variables mean different things at different times: `xa_start` (start of step, or a sub-step checkpoint), `max_residual` (three meanings in one iteration), a negative `dt` meaning "do not age", and `rot_scr%dm` (an interface density) versus `star%dm` (a shell mass).
   - Scratch modules carry data between routines, so results depend on call order.
4. **Duplication frozen in place.** Byte-identical verification was the right tool for the line-by-line conversion. It also forbade merging duplicates, since merging changes the last bits; `burn_lib.f90:1166` says so explicitly. The duplicates have since drifted apart: the two d(p,g) copies differ by 0.68%.

### 1.3 Bugs that came from the organization

These were found while tracing the tasks above. Each traces to one of the root causes. Details are in Appendix B.

- The GYRE pulse-file columns `eps_eps_T` and `eps_eps_rho` are wrong by a factor of ε. A field named as a log derivative holds the absolute derivative and is multiplied by ε a second time. *(Undocumented contract.)*
- Deuterium depletion and deuterium energy generation use rates that differ by 0.68%. *(Duplication.)*
- `star%am_transport_convective_flag` is read by the Ω(j) inversion on every Newton iteration but is never assigned. Every assignment goes to a local array with the same name. *(Global state plus name shadowing.)*
- The advective rotational-transport option (`ldifad`) cannot work, because its coefficient arrays are never filled. *(Scratch-module data flow.)*

---

## 2. Principles and non-goals

### 2.1 Principles

- **P1. Explicit state.** No routine below the driver imports `star`. Each routine receives the parts of the model it needs as arguments, with `intent(in)`, `intent(out)` or `intent(inout)`. The compiler then enforces the read/write contract.
- **P2. Physics never inside the numerics.** A solver receives either a callback that returns residuals and Jacobian entries, or precomputed coefficient arrays. Physics routines never see solver state.
- **P3. One extension point per domain.** The extension points are energy sources, reactions, opacity and EOS table sources, and transport mechanisms. Each is a list of entries sharing a single interface, the same idea as MESA's `other_*` hooks.
- **P4. Contracts written down.** Every interface states the units, the log base and the derivative convention of each argument.
- **P5. One name per quantity, one quantity per name.** A glossary (`docs/GLOSSARY.md`) is the authority on names.
- **P6. Validation against physics, not only against previous output.** Byte pins guard refactors. Physics references guard correctness.
- **P7. The timestep is a visible sequence of stages.** Every iteration that couples two stages, such as the mixing-and-structure splitting iteration or the structure-and-rotation iteration, is written in the step driver, not hidden inside a callee.

### 2.2 Non-goals

- **No change of method.** Henyey with an envelope fitting point, operator splitting, the four iteration levels, shellular rotation and the envelope "triangle" all stay. Changing any of them would be a new code, and a separate decision.
- **No change of language.** Modern Fortran (2008 or later) supports everything here. It keeps the code familiar to the YREC and MESA communities.
- **No new physics within this programme.** New physics comes after, through the extension points.
- **No silent changes to results.** Every change to output is deliberate, measured and recorded (section 5).

---

## 3. The timestep today

The order of operations in `core/evolve_step.f90` is sound: it is the physics. This section records it with current names (legacy names in parentheses) so that section 4 can say precisely what changes.

```
evolve_step
 retry_step loop (cut dt and retry on divergence)
  reload_model_if_diverged      read_starting_model (starin)
  advance_composition_and_age   massloss; xa_start := xa; mix (mix) at level 1; age += dt
  rezone_or_snapshot            rezone (hpoint), or keep the old grid during the He flash
  solve_structure               levels 1-3, each a call to henyey_iterate (crrect):
                                  surfbc (surfbc)          envelope triangle, when asked
                                  mix (mix)                again, at levels >= 3
                                  henyey_coefficients (coefft) per shell:
                                      eos_get, kap_get, engeb, gradients, diagnostics
                                      henyey_eliminate (reduce)
                                  henyey_solve (hsolve)
                                  omega_from_j + rotation_shape_factors, when rotating
  converge_with_rotation        num_rotation_structure_iters times:
                                  level 4; homogenize_convection_zones (mixcz)
                                  evolve_angular_momentum (getw)
                                      sub-steps: mid_timestep_model (midmod)
                                                 secular_transport (seculr) ...
                                  rotation_shape_factors (fpft)
 locate_shell_boundaries (findsh)
 burn_light_elements            lirate88, liburn
 rebuild_envelope               after mass loss
 compute_timestep (htimer)
 build_stitched_model           interior + envelope + atmosphere
 compute_observables
 output_write_model
 check_stop_conditions
```

The order is fine. The problems lie in how these stages communicate:

- **Through the global object.** The stages exchange data through roughly 265 fields of `star`, plus the scratch modules `rot_scr`, `mix_scr` and `pt_scr`. No stage states its inputs or outputs.
- **Stages nest inside one another.** Mixing, envelope integration and rotation shape updates run inside the Newton driver. The splitting iteration between composition and structure exists, but only as "`mix` is called again when `iteration_level > 2`" inside `henyey_iterate`. (Fixed for the composition step on `yrec-next`, 2026-10: `solve_level` in `evolve_step` now runs the surface boundary condition, the composition step and the Newton loop as three calls.)
- **Flags are encoded in data.** "Do not age this model" is a negative `star%dt`, and four routines read its sign.
- **Snapshots are scattered.** There are three separate snapshot families in `star`: `*_start` (start of step), `orig_*` (for the structure⇄rotation iteration) and `old_*` (rotation). A fourth, `rot_scr%old_*`, holds the rezone snapshot. Each is copied with hand-written loops.

---

## 4. Target architecture

### 4.1 Layers

```
  run driver          run_yrec: controls, tables, kind cards, calibration loops
       |
  step driver         evolve_step: the stage sequence and the coupling iterations
       |
  stages              composition, rezone, structure, rotation, post-step,
       |              diagnostics, output, timestep
       |
       +--> numerics          Henyey eliminate/solve, network solver, diffusion
       |                      solver, ODE integrators, interpolation
       |
       +--> physics providers EOS, opacity, atmosphere, reactions, energy sources,
                              transport mechanisms, mixing criteria
                                   |
                              tables and their readers
```

The rules are:

- **Only stages touch the model state.** They read it and write it.
- **Numerics never calls a physics provider directly.** It receives a callback or coefficient arrays from the stage.
- **Physics providers never see the model state or solver state.** They receive point values and their own controls.

`check_boundaries` already enforces some import rules. It gains one per rule above.

### 4.2 State

The single `star_info` splits into typed parts, and the model itself becomes a type that can be copied:

```fortran
integer, parameter :: dp = kind(1.0d0)

type :: model_t                      ! everything that defines the star at one instant
   integer :: nz = 0
   real(dp), allocatable :: log_mass(:), m(:), dm(:)               ! grid
   real(dp), allocatable :: logP(:), logT(:), logR(:), logRho(:), &
                            luminosity_lsun(:)                      ! structure
   logical,  allocatable :: convective_flag(:)
   real(dp), allocatable :: xa(:,:)                                 ! composition (n_species, nz)
   real(dp), allocatable :: omega(:), j_rot(:), i_rot(:), &
                            eta_squared(:), fp_rot(:), ft_rot(:)    ! rotation
   type(envelope_fit_t)  :: envelope                                ! fit coefficients, triangle, Teff
end type

type :: star_t
   type(model_t)       :: now, start        ! start-of-step snapshot: s%start = s%now
   type(step_state_t)  :: step              ! dt, age, model number, aging flag, divergence
   type(diagnostics_t) :: diag              ! eps channels, luminosity breakdown, pulse columns,
                                            ! observables; written only after convergence
   type(star_job)      :: job
   type(controls_t)    :: ctrl              ! per-domain parts, section 4.6
end type
```

What this changes:

- **Snapshots become one assignment.** `s%start = s%now` replaces the hand-copied `*_start`, `orig_*` and `old_*` families. The structure⇄rotation iteration keeps its own copy the same way.
- **No fixed shell cap.** Arrays are allocatable, sized when a model is loaded or rezoned. This removes the `json = 5000` limit and the static storage that blocks running several stars, or threads, in one process.
- **Flags become fields.** `s%step%aging` replaces the negative `dt`.
- **Diagnostics are separate from state.** They live in `s%diag` and are filled by a pass after convergence, not during Newton iterations.
- **The global goes away.** `star` stops being a module variable. `run_yrec` owns a `type(star_t)` and passes it down; the C API (`yrec_capi`) and `pyyrec` hold a handle to it.

### 4.3 Stage contracts

Each stage is a procedure whose arguments are the parts of the state it may touch. This table is the contract. The intents in the code must match it.

| # | Stage (legacy) | Reads | Writes | After the stage |
|---|---|---|---|---|
| 1 | `begin_step` | `now`, `step` | `start`, `step%aging` | `start` equals `now` |
| 2 | `mass_change` (massloss, mdot) | `now`, `ctrl%wind` | `now` grid and total mass, envelope target | mass budget closed |
| 3 | `composition_step` (mix) | `now` structure, `start%xa`, `ctrl%mixing`, `ctrl%net` | `now%xa`, convective bounds | abundances sum to 1 |
| 4 | `rezone` (hpoint) | `now` | `now` (all per-zone arrays) | the only stage allowed to change `nz` |
| 5 | `structure_solve` (crrect levels 1–3) | `now%xa`, `now` rotation (held fixed), `now%envelope`, `ctrl%solver` | `now` structure, `now%envelope` | residuals below tolerance |
| 6 | `structure_rotation_iteration` (level 4, getw, fpft) | `now`, `start` | `now` rotation, `now%xa`, `now` structure | splitting converged |
| 7 | `post_step` (findsh, liburn) | `now` | `now%xa` light elements, shell bounds | |
| 8 | `diagnostics` (stitched model, observables) | `now`, `start` | `diag` only | |
| 9 | `output` | `diag`, `now` | files only | |
| 10 | `next_timestep` (htimer) and stops | `now`, `diag`, `step` | `step` | |

Coupling iterations are written in the step driver. Two cases need care to keep today's numerics exactly:

- **The splitting iteration at levels ≥ 3.** The structure solver is split at the point where it used to call the composition step: `refresh_surface_bc` (which reads the surface composition, so it must come first), then the composition step, then `henyey_iterate`. The step driver makes the three calls, so the order of operations is unchanged and is visible in `evolve_step`. (Done on `yrec-next`; a plain split reads more easily than the `before_iteration` procedure argument first proposed here.)
- **The Ω(j) inversion during Newton iterations.** It is passed in the same way, as an `after_correction` procedure argument.

### 4.4 Structure: separating physics from the Henyey numerics

**Today.** `henyey_coefficients` evaluates the microphysics and forms the per-shell terms `q*` and their partials in one loop, interleaved with calls to `henyey_eliminate`. The per-shell type `henyey_shell_terms` carries only the partials that happen to be nonzero for today's physics: there is no `ql_dr` or `ql_dl`.

**Target.**

- **A pure function `shell_terms(point, providers) -> terms`.** It returns the four equation terms and all sixteen partials, plus the side information the diagnostics pass needs.
- **The eliminator** (`henyey_eliminate`, `henyey_solve`) takes arrays of `terms` and knows nothing about physics. Its pivots become a generic 4×4 block elimination rather than hand-derived algebra.
- **The energy equation's source is a sum over registered energy sources** (section 4.5a).
- **The diagnostics pass** calls `shell_terms` once on the converged model and fills `s%diag`.

**Conventions.** These must go into the code header and the glossary on day one. They hold today and stay:

- **Variables and coordinate.** The dependent variables are log₁₀ P, log₁₀ T, log₁₀ R and L in L☉. The coordinate is log₁₀ m.
- **The difference equation.** For equation k at shell i:
  `rhs_k = ½ Δlog₁₀m · (q_k(i−1) + q_k(i)) − (y_k(i) − y_k(i−1))`
  (`henyey_elimination.f90`, top of the routine).
- **Partials.** `q*_dX` is ∂q/∂ln X: a natural-log derivative for P, T and R. The eliminator multiplies ½Δlog₁₀m by ln 10 to convert to the log₁₀ variables. The luminosity partial is therefore pre-divided by ln 10 (`qt_dl = qt/(L ln 10)`, `henyey_coefficients.f90`).
- **The luminosity term.** `ql = (ln 10 · m / L☉) · ε`, where ε is the total specific energy rate in erg g⁻¹ s⁻¹. So `ql_dT` and `ql_dP` take the **absolute** derivatives dε/dln T and dε/dln ρ, which is what `engeb` returns.
- **The envelope.** The envelope integration (`envelope_derivs`) carries m, T and R against log P at constant L. Any energy source above the fitting mass is excluded. This must be stated in the energy-source interface, and checked: a source that is nonzero above the fitting point should cause a warning.

### 4.5 The four extension points

All four use the same pattern: a derived type with a name, its parameters and a procedure-pointer component, collected in a list that is built at startup from the controls. A user-supplied entry is added in exactly the same way as a built-in one.

#### (a) Energy sources

```fortran
type :: energy_source
   character(len=32) :: name = ''
   integer :: lum_channel = 0         ! slot in diag%luminosity_breakdown it reports to
   procedure(energy_source_eval), pointer, nopass :: eval => null()
end type

abstract interface
   subroutine energy_source_eval(pt, eps, deps, ierr)
      import :: shell_point, dp
      type(shell_point), intent(in) :: pt   ! m, r, P, T, rho, xa, dt, start-of-step values
      real(dp), intent(out) :: eps          ! erg/g/s
      real(dp), intent(out) :: deps(4)      ! d eps / d(lnP, lnT, lnR, L): absolute, not logarithmic
      integer,  intent(out) :: ierr
   end subroutine
end interface
```

- **Built-in entries:** nuclear, neutrino losses, gravitational (the entropy terms), rotational kinetic energy.
- **User hook:** `other_energy`, for tidal or magnetic heating and similar.
- **Luminosity channels:** each source reports to its own channel. This removes the trap where the luminosity renormalization in `henyey_iterate` folds unknown heating into the nuclear channels.

#### (b) Reactions and the network

```fortran
type :: reaction
   character(len=16) :: name = ''               ! e.g. 'd_pg_he3'
   integer  :: in_species(3) = 0, in_count(3) = 0
   integer  :: out_species(3) = 0, out_count(3) = 0
   real(dp) :: q_mev = 0                        ! energy deposited
   real(dp) :: q_nu_mev = 0                     ! energy lost to neutrinos
   procedure(rate_eval), pointer, nopass :: rate => null()
end type

abstract interface
   subroutine rate_eval(t9, rho, lam, dlnlam_dlnT, dlnlam_dlnrho)
      import :: dp
      real(dp), intent(in)  :: t9, rho
      real(dp), intent(out) :: lam, dlnlam_dlnT, dlnlam_dlnrho   ! unscreened N_A<sigma v>
   end subroutine
end interface
```

- **The species table** (names, A, Z, the index constants) is the single source of the species count. Arrays are sized from it, never from a literal.
- **Screening** is applied by the network driver, using the charges in the species table, so rate functions stay bare fits.
- **One rate evaluation per zone feeds both consumers.** Composition uses the rates and the stoichiometry. The nuclear energy source uses the same rates, the Q-values and the derivatives. The duplication between `rates` and `engeb` disappears.
- **A generic backward-Euler Newton solver** builds its Jacobian from the reaction table. It replaces the hand-coded 7×7 system in `solve_composition` (legacy `kemcom`). Light elements become ordinary species; species can be flagged as trace if the group wants to keep their separate treatment.
- **Rate scale factors** (`cross_section_scale`, the Monte Carlo perturbations) become one multiplier per reaction entry, referred to by name rather than by slot number.

#### (c) Microphysics table sources

**Opacity:**

```fortran
type :: kap_source
   character(len=32) :: name = ''
   real(dp) :: logT_min, logT_max, blend_dex    ! validity window and blend width
   procedure(kap_load), pointer, nopass :: load => null()
   procedure(kap_eval), pointer, nopass :: eval => null()   ! (logRho, logT, X, Z) -> kap, dlnkap/dlnRho, dlnkap/dlnT
end type
```

- The ordered list of sources replaces the if-chain in `kap_eval`.
- One generic driver does the temperature and Z blending. It blends derivatives consistently, including the weight-derivative term that is missing today.
- Opening and reading tables moves out of `io/read_controls.f90` into each source's `load`.

**EOS:**

- **A named result type** (`eos_result`): log ρ, ∂ln ρ/∂ln P at fixed T, ∂ln ρ/∂ln T at fixed P, c_P, ∇_ad, Γ₁, β, η, the ionization fractions and the second derivatives the solver uses. Each field states its units and its definition. It replaces the 24 positional slots. (Done 2026-10-11: `eos/eos_results.f90`; Γ₁ is still derived by the caller.)
- **One integration pattern instead of two:** an ordered list of sources with validity windows and a blend. The Yale/Saha base, SCV, OPAL 95/01/06 and MHD all become entries.

#### (d) Rotational transport mechanisms

```fortran
type :: transport_mechanism
   character(len=32) :: name = ''
   real(dp) :: am_efficiency = 1, mix_efficiency = 1
   procedure(transport_eval), pointer, nopass :: eval => null()
end type

abstract interface
   subroutine transport_eval(iface, nu_am, d_mix, ierr)
      import :: interface_state, dp
      type(interface_state), intent(in) :: iface   ! r, rho, P, g, N_T^2, N_mu^2, K, nu, Omega, q
      real(dp), intent(out) :: nu_am(:), d_mix(:)  ! cm^2/s at each cell interface
      integer,  intent(out) :: ierr
   end subroutine
end interface
```

- **The interface state is built once per Picard iteration.** Today it is assembled across `mid_timestep_model`, `rotation_stability_setup`, `rezone` and `mix`.
- **Each mechanism returns diffusion coefficients directly.** The existing mechanisms (Eddington–Sweet, GSF, secular shear, dynamical shear) each convert their velocity to D inside their own routine, using today's formula from `diffusion_velocity_scales`.
- **The solver takes only the summed `nu_am` and `d_mix` arrays.** Transport regions are found from D > 0, not from the summed velocity.
- **The control logic lives in the rotation stage driver,** not inside the physics or the solver:
  - sub-stepping, with a criterion that includes D·Δt/Δr² as well as ΔΩ;
  - Picard averaging;
  - convergence checks, which return a verdict instead of halving the timestep themselves.
- **Advective transport.** If the group wants it (Zahn 1992; Maeder & Zahn 1998), the interface gains an optional advective velocity U(r). Today's `ldifad` path, which attempted this, is dead code and would be removed.

### 4.6 Controls

- **Each control gets a domain.** The registry (`defaults/controls_registry.tsv`) gains a `domain` column (`solver`, `mixing`, `net`, `eos`, `kap`, `atm`, `rotation`, `wind`, `output`, `job`).
- **Each domain gets a generated controls type.** The generator writes a type per domain, and `s%ctrl` contains one of each.
- **Libraries receive their own controls at startup.** Each one gets its type once, in `*_init`, and never imports `star`. Today eos, kap, atm and net read about 60 `star%ctrl` fields between them.
- **Startup consistency checks** (`warn_inconsistent_controls`) gain one rule per new registry entry where it applies.

### 4.7 Names and conventions

- **The glossary** (`docs/GLOSSARY.md`) has one row per quantity: symbol, units, current name, legacy name or names, and the file that owns it. It is the authority on names, and review checks new code against it.
- **Naming rules:**
  - `_start` means start of step.
  - Derivatives are named `d<f>_dln<x>` and state whether they are absolute or logarithmic.
  - No name is reused for a different quantity.
- **Units: CGS unless the name says otherwise** (adopted 2026-10). Every internal variable is in CGS by default: a mass in g, a radius in cm, a luminosity in erg/s, a time in s, a rate per s. A variable in any other unit carries the unit as a suffix: `mass_msun`, `radius_rsun`, `luminosity_lsun`, `age_gyr`, `dt_yr`, `nu_max_uhz`, `rate_per_gyr`. A `log` prefix means log₁₀ of the quantity in the unit its name implies (`log_mass` is log₁₀ of grams; `log_l_lsun` is log₁₀ L/L☉). Dimensionless quantities need no suffix. When restructuring code that works in another unit system (the network's per-Gyr-per-amu rates, the solver's L☉ terms), prefer converting to CGS; otherwise rename with the suffix. The history and profile output columns keep their MESA names (`star_mass` in M☉, `log_L` in L☉); the output layer maps internal names to them.
- **Legacy names** stay in comments at each declaration (as the code does now). Appendix A gives the routine-level map.

### 4.8 What a newcomer task looks like afterwards

| Task | Today | Target |
|---|---|---|
| Extra energy term | ~5 files, undocumented conventions, renormalization trap | one module with an `energy_source`, one registration line, one control row |
| Update a reaction rate | two copies, derivatives by hand | edit one `rate` function |
| Track a new isotope | 34+ files | one species-table entry plus its reactions |
| Opacity source | ~18 files | one module with a `kap_source`, one registration line, one control row |
| EOS option | ~16 files, two patterns | one module with an EOS source, one registration line |
| AM transport process | ~12 files, hidden velocity gating | one module with a `transport_mechanism`, one registration line |

---

## 5. Verification policy

Three tiers, each matched to the kind of change.

### Tier R: refactor (today's practice)

- **Applies to:** renaming, explicit argument passing, splitting routines, generated controls types.
- **Requirement:** byte-identical on the regression net. This continues unchanged.

### Tier S: subsystem replacement

- **Applies to:** replacing a subsystem (the network, the structure terms, opacity and EOS dispatch, rotational transport).
- **Side by side during transition.** The old and new implementations coexist behind a control (for example `net_implementation = 'legacy' | 'table'`) until the new one is accepted. The group can then compare them on any case.
- **Measured on the full configuration matrix.** For each case, record the difference in log L, log Teff, log R, X_c and age at fixed landmarks (ZAMS, X_c = 0.3, TAMS, RGB base, solar age), and in GYRE frequencies for the solar model.
- **Tolerance principle.** Differences must be small compared with the effect of any physics choice the group would care about. Any larger difference must be explained, usually by a bug in one of the two implementations, and the explanation recorded.
- **Acceptance.** The pins are reseeded once, deliberately, and the old implementation is deleted.

### Tier P: physics validation (new; a prerequisite for tier S)

- **Pulse-file contract test:** re-derive N² from ρ, P, Γ₁ and r. Check the units of the derivative columns (this would have caught the ε factor).
- **Network unit tests:** each rate and its derivatives on a T₉ grid against REACLIB or NACRE II values, plus Q-value bookkeeping (energy released equals the mass defect).
- **EOS and opacity spot checks** against table values and thermodynamic consistency relations.
- **Cross-code track:** one (M, Z, α) track compared with MESA at matched central-hydrogen points, and GYRE frequencies at matched Δν. Tolerances to be documented.
- **Rotation benchmark:** reproduce a published YREC solar spin-down or open-cluster Ω(t) result, and pin Ω(r) at landmark ages. Rotation has no truth-level test today.

---

## 6. Migration plan

Each phase leaves the code working and producing science. Every phase also does its share of the state work: it converts the parts of `star` it owns into explicit arguments.

### Phase 0: groundwork (no change to results, except the bug fixes)

- **Agree on this document.**
- **Glossary and conventions.** Write `docs/GLOSSARY.md` and put the Henyey conventions of section 4.4 into the code headers.
- **Refresh `docs/HOW_YREC_WORKS.md`** to current names; ten of its file references point to files since renamed. Add a developer guide, "how to add X", for today's code.
- **Build the tier P tests.**
- **Add the `domain` column** to the controls registry.
- **Fix the four bugs in Appendix B,** each with a recorded reseed.

### Phase 1: network

- **Status: complete (2026-10-11, branch `yrec-next`; ROADMAP section 13).** The reaction table and its generic implicit solver; the single rate evaluation behind composition and energy; the exact Jacobian terms with a finite-difference test; each row's Q-value, chain and rate (rate law, S-factor fit, scale slot, screening), so a reaction is added as one row; the equilibrium burn removed. net/ is star-free and reached only through net_lib. Rate laws are dispatched by an integer rather than a procedure pointer (the table stays a compile-time constant). Deuterium and the light elements stay special operators; the species count moves to Phase 5.
- **Work:**
  - species table and reaction table;
  - `rates` returns derivatives, and `engeb` consumes it, one deliberate reseed (the d(p,g) pair was already merged in Phase 0);
  - generic backward-Euler solver;
  - species count taken from the table.
- **Exit criteria:** network unit tests pass; drift on the configuration matrix recorded; one copy of every rate.
- **Why first:** it is self-contained, it already has a standalone test (`net/test`), and its duplication is already producing inconsistencies.

### Phase 2: structure

- **Work:**
  - explicit state passing through `core/`;
  - `shell_terms` with the full Jacobian;
  - energy-source list;
  - diagnostics pass after convergence;
  - coupling iterations moved into the step driver (section 4.3);
  - the `aging` flag replaces the negative `dt`.
- **Exit criteria:** tier S on the configuration matrix; an `other_energy` example in `examples/`.

### Phase 3: microphysics

- **Work:**
  - domain controls types;
  - `eos_result`;
  - opacity and EOS source lists with the generic blend drivers;
  - table loading moved into the sources;
  - eos, kap, atm and net import nothing from the star layer, enforced by `check_boundaries` (done 2026-10-11: each installs a settings record at init).
- **Exit criteria:** tier S, including the MHD EOS path, which gets a pin first (test_eos pins every other EOS configuration; it runs MHD only when `YREC_MHD_TABLES` points at the tables, which the repository does not ship).

### Phase 4: rotation

- **Prerequisite:** the rotation benchmark from tier P.
- **Work:**
  - interface state;
  - mechanism list;
  - solver driven by D arrays;
  - sub-stepping, Picard averaging and convergence moved into the stage driver;
  - one name per quantity (the alias table in the `secular_transport` header goes away);
  - `ldifad` removed or rebuilt as an advective mechanism.
- **Why last:** rotation is the largest body of physics and has no truth-level test today.

### Phase 5: finish the state

- **Work:**
  - `star_t` as in section 4.2;
  - allocatable arrays;
  - scratch modules and the SAVEd driver state removed;
  - `star` no longer a module variable.
- **Result:** several stars per process, and threads within `libyrec`, become possible.

---


## Appendix A: legacy routine names and current files

Resolved from the git rename history of `src/`. Routines merged into a module are listed under that module.

| Legacy | Current file |
|---|---|
| `starin` | `core/read_starting_model.f90` |
| `crrect` | `core/henyey_iteration.f90` (`refresh_surface_bc`, `henyey_iterate`) |
| `coefft` | `core/henyey_coefficients.f90` |
| `reduce` | `core/henyey_elimination.f90` |
| `hsolve` | `core/henyey_solve.f90` |
| `surfbc` | `core/surfbc.f90` |
| `envint` | `core/envelope_integration.f90`, `core/envelope_integration_kernel.f90` |
| `qatm`, `qenv` | `core/atmosphere_derivs.f90`, `core/envelope_derivs.f90` |
| `physic` | `core/shell_physics.f90` |
| `gettau` | `core/turnover_timescale.f90` |
| `chkcal`, `setcal` | `core/check_solar_calibration.f90`, `setup/setup_solar_calibration.f90` |
| `chkscal`, `setscal` | `core/check_star_calibration.f90`, `setup/setup_star_calibration.f90` |
| `engeb` | `net/nuclear_energy.f90` (`nuclear_energy`), adapter `compute_zone_energy` in `core/burn_step.f90` |
| `deutrate`, `dburn`, `dburnm` | `net/light_elements.f90` (`deuterium_rate`, `burn_deuterium_region`) |
| `liburn`, `lirate88` | `core/burn_step.f90` (`burn_li_be`, `store_start_li_be_rates`/`store_end_li_be_rates`) |
| `nrate` (rates) | `net/net_lib.f90` (`compute_network_rates`, `compute_reaction_rates`) |
| `mix` | `core/composition_step.f90` (`advance_composition`) |
| `convec`, `oversh`, `hsubp`, `mixcz` | `convection/convection_zones.f90` (`find_convection_zones`, `extend_by_overshoot`, `compute_scale_height`, `homogenize_convection_zones`) |
| `sconvec` | `convection/semiconvection.f90` (`extend_by_semiconvection`) |
| `kemcom` | `net/network.f90` (`burn_implicitly`) |
| `tpgrad` | `convection/temperature_gradients.f90` (`compute_temperature_gradients`) |
| `rotmix`, `bursmix`, `bsrotmix` | `rotation/rotmix.f90`, `rotation/burn_settle_mix.f90`, `rotation/burn_mix_extrapolated.f90` (the last two commented out) |
| `hpoint` | `setup/rezone.f90` |
| `findsh` | `setup/locate_shell_boundaries.f90` |
| `setups` | `setup/setups.f90` |
| `htimer` | `util/compute_timestep.f90` |
| `getw` | `rotation/evolve_angular_momentum.f90` |
| `midmod` | `rotation/mid_timestep_model.f90` |
| `solid` | `rotation/solid_body_omega.f90` |
| `viscos` | `rotation/viscos.f90` |
| `fpft` | `rotation/shape/rotation_shape_factors.f90` |
| `seculr` | `rotation/seculr/secular_transport.f90` |
| `codiff` | `rotation/seculr/diffusion_velocity_scales.f90` |
| `rotgrid` | `rotation/seculr/am_transport_grid.f90` |
| `dcoeft`, `dadcoeft` | `rotation/seculr/am_diffusion_coeffs.f90`, `am_advection_diffusion_coeffs.f90` |
| `checkj`, `checkc` | `rotation/seculr/check_angular_momentum.f90`, `check_composition.f90` |
| `mwind`, `wind`, `mcowind`, `cowind`, `amcalc` | `rotation/magnetic_braking.f90` (`brake_convection_zone` for mwind/wind, `compute_solid_body_braking` for mcowind/cowind) |
| `massloss`, `mdot` | `core/mass_change.f90` (`change_mass`, `apply_mass_change`) |
| `eqstat`, `eqstat2` | `eos/blended_eos.f90` (`compute_blended_eos`, `blend_eos`) |
| `meqos`, `mhdtbl` | `eos/mhd/mhd_eos.f90` (`compute_mhd_eos`, `read_mhd_tables`) |
| `oeqos`, `oeqos01`, `oeqos06`; `eqbound`, `eqbound01`, `eqbound06` | `eos/opal/opal_eos.f90` (`compute_opal_eos` for all three; `find_opal_edge_ramp` over `find_opal95_ramp` (`opal95_eos.f90`), `find_opal2001_ramp`, `find_opal2006_ramp`) |
| `esac01`/`esac06`, `readcoeos01`/`06`, `t6rinteos01`/`06` | `eos/opal/opal_eos.f90` (`interpolate_opal_eos`, `readcoeos`, `t6rinteos`, branching on the vintage; `esac01`/`esac06` remain as entry points) |
| `eqscve`, `eqscvg`, `setscv` | `eos/scv/scv_eos.f90` (`compute_scv_envelope_eos`, `compute_scv_eos`, `build_scv_envelope_table`) |
| `eqsaha`, `eqrelv` | `eos/yale/yale_eos.f90` (`compute_saha_eos`, `compute_fully_ionized_eos`) |
| `mu` | `eos/eos_mixture.f90` (`compute_mean_weights`) |
| `surfp`, `kcsurfp` (`surface_p_interp`) | `atm/tables/kurucz_atmosphere.f90` (`interpolate_surface_pressure`; the reader is `read_surface_pressure_table`) |
| `alsurfp`, `alfilein`, `altabinit` | `atm/tables/allard_atmosphere.f90` (`interpolate_allard_surface`, `read_allard_tables`, `find_allard_table_bounds`) |
| `getopac` (`kap_eval`), `setupopac` (`load_opacity_tables`), `surfopac` (`refresh_surface_tables`) | `kap/kap_lib.f90` (all three) |
| `ll95tbl`, `getopal95` | `kap/opal95/opal95_opacity.f90` (`read_opal95_tables`, `compute_opal95_opacity`) |
| `rdlaol`, `gtlaol`, `rdzlaol`, `gtpurz` | `kap/laol89/laol_opacity.f90` (`read_laol_tables`, `compute_laol_opacity`, `read_pure_z_table`, `compute_pure_z_opacity`) |
| `readalex06`, `getalex06` | `kap/alex06/alex06_opacity.f90` (`read_alex06_tables`, `compute_alex06_opacity`) |
| `kurucz`, `condopacpint`, `condopacp` | `kap/kurucz90/kurucz_opacity.f90` (`compute_kurucz_opacity`), `kap/conductive/conductive_opacity.f90` (`compute_conductive_opacity`, `condopacp`) |
| `trapzd`, `qgauss`, `func` | `rotation/shape/shell_quadrature.f90` (`refine_distortion_integral`, `average_over_equipotential`, `equipotential_integrand`) |
| `tridia` | `rotation/seculr/secular_transport.f90` (`solve_angular_momentum_transport`) |
| `intpt` | `eos/mhd/mhd_eos.f90` (`interpolate_mhd_table`) |
| `tridiag_gs`, `ctridi` | `numerics/linear_solvers.f90` (`solve_tridiagonal`) |
| `kspline`, `splinc`; `ksplint`, `splintd2` | `numerics/cubic_spline.f90` (`natural_spline`; `splint`) |
| `bsstep`, `mmid`, `ratext` | `numerics/bulirsch_stoer.f90` (`bsstep` integrates an `ode_system`) |
| `parmin` | `io/read_controls.f90` (`read_controls`), with its steps in `io/legacy_controls.f90`, `run_files.f90`, `derived_controls.f90` and `run_cards.f90` |
| `getyrec7`, `getmodel2`, `wrtlst` | `io/model_files.f90` (`read_yrec7`, `read_model2`; the `.mod` format's `write_mod_model` and `read_mod_model` replaced `wrtlst`) |
| `wrthead` | `io/yrec_output.f90` (`output_run_header`) |
| `wrtout` | retired; replaced by `io/yrec_output.f90` |

## Appendix B: issues found in the October 2026 audit

The first four were confirmed by reading the code; none has been run or fixed yet. The last three were reported during the audit but not yet checked.

1. **Pulse ε-derivative columns wrong by a factor of ε.**
   - `engeb` returns the absolute dε/dln T = Σ εᵢ·dln rᵢ/dln T (`burn_lib.f90:1772`).
   - The structure Jacobian uses it correctly as that (`henyey_coefficients.f90:253`).
   - The same value is stored in `star%pulse_dlneps_dlnt`, a log-derivative name.
   - Both pulse writers then multiply it by ε again (`stitched_model.f90:579–580`, `write_gyre_pulse.f90:112–113`).
   - `ROADMAP.md` closed this as "not a bug" on the premise that the field holds a log derivative.
   - Affects non-adiabatic GYRE runs only.
   - Origin: not F77. Commit 2ed398ca (2026-08-30) multiplied the field by ε, misled by its name.
   - Fixed (2026-10): the writers copy the field, and the field and its two copies are renamed `*deps_*`.
2. **Two d(p,g) rates.**
   - `deutrate` divides by A = 2 (`3.0115d23`, `burn_lib.f90:703`).
   - `engeb` divides by the deuteron mass 2.013553 (`burn_lib.f90:1221`).
   - The two differ by 0.68%.
   - Fixed (2026-10): one shared `dpg_rate` in `burn_lib`, using the deuteron mass. Its arithmetic is `engeb`'s, bit for bit (checked at 16,441 (T, ρ) points), so only deuterium depletion moves: at the solar age |Δ log L| ≤ 4.4×10⁻⁵ and |ΔX_c| ≤ 1.2×10⁻⁴; during pre-main-sequence deuterium burning |Δ log L| ≈ 4×10⁻⁴.
3. **`star%am_transport_convective_flag` is never assigned.**
   - It is read at `henyey_iteration.f90:311` by `omega_from_j`, so with rotation on, every zone is treated as radiative in the Ω(j) inversion during Newton iterations.
   - Every assignment in the code goes to a local array of the same name.
   - Origin: inherited. In the F77 code, `main.f`'s `LCZ(JSON)` was a local passed only to `CRRECT` → `GETROT` and never set. `rezone` computed the right flags but stored them in the write-only `star%cz_flag_start`; `read_starting_model` computed them and discarded them.
   - Fixed (2026-10): both now write `star%am_transport_convective_flag`, and `cz_flag_start` is removed. Only rotating runs change: at the solar age |Δ log L| ≤ 3×10⁻⁵, |Δ log R| ≤ 1.2×10⁻⁵, |ΔX_c| ≤ 7×10⁻⁵; along the pre-main-sequence tracks |Δ log L| ≤ 0.002.
4. **`ldifad` cannot work.**
   - The advective path in `am_advection_diffusion_coeffs` reads 11 `rot_scr` arrays that nothing writes (`geometric_factor`, `dchi_dr_edge`, `am_2nd/3rd_deriv_coeff`, `eq_velocity_coeff*`, `mixing_geometric_factor`, `equatorial_radius`, ...), so it divides by a zero geometric factor.
   - It is inherited: in the original F77 only `dadcoeft.f` referenced `/DIFAD4/`, `/EGRIDCHI/`, `FGEOM`, `ECOD5/6` and the `EV*` arrays.
   - The flag takes effect only in rotating runs with instability transport. 16 legacy decks and 316 inlists set it, but none of them rotate, so for them it is inert.
   - Fixed (2026-10): `secular_transport` now stops with an error when the flag would take effect.
5. **`beta14` is defined differently by EOS path** (*reported, not yet checked*). It means 4(1−β)/β on the Yale path but 1−β after OPAL or MHD (`blended_eos.f90:324` vs `:484`).
6. **Opacity ramp derivatives are inconsistent** (*reported, not yet checked*). The ramp blends κ linearly but averages the dln κ derivatives without weighting them and drops the dw/dln T term (`kap_lib.f90:312–317`).
7. **`max_delta_j_by_iter(16)` can overflow** (*reported, not yet checked*). It overflows when `max_diffusion_iters > 16` (`check_angular_momentum.f90:91`).

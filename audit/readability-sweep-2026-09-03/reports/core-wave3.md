# core -- readability wave 3 report

Worktree `/Applications/YREC-wt/core`, branch `rs/core`, based on
yrec-modern `6449462`. Started with `make clean`. Every commit passed
gate1 (build + pinned case byte-identical); the full 37-case pin
selection and the 22-case aux battery were run once after the last
commit (results at the end).

## Commits (in order)

| hash | item | summary |
|---|---|---|
| `ceac51e` | 1 | share the envelope-refit shell fill (new `core/envelope_refit_lib.f90`) |
| `fdc2177` | 2 | one kernel per repeated profile/pulse column expression in stitched_model |
| `6b50553` | 3 | drivers apply compute_timestep's structure-limit switch-off to star%job |
| `789adce` | 4 | surfbc's luminosity dummy is log10 L -- `log_luminosity` |
| `2bc14b6` | 5 | name the surface-geometry expressions the envelope drivers share |
| `5d75705` | 6a | read_starting_model passes `zone_index=-1` instead of writing `star%iovim` |
| `1530b2b` | 6b | atm_get / integrate_envelope_atmosphere `print_flag` is intent(in) |

All commits: source edits only (plus `src/deps.mk` in `ceac51e` and
`2bc14b6`); no output dirs, no `.short` files, nothing under
`standard/`.

## Per item

### 1. One envelope refit -- extracted the identical inner part only

Not mergeable as one routine: after normalising names the two refits
differ in arithmetic and in side effects (list under Deferred 1). The
per-shell fill loop ("ASSIGN NEW POINTS" through the photosphere
boundary guard) is token-identical in both apart from one comparison
(`.lt.` vs `.le.` against `star%senv`).

- `src/core/envelope_refit_lib.f90` (new): module `envelope_refit_lib`,
  public `append_envelope_points(accept_point_at_fit_mass,
  species_end_index, composition, log_density, log_luminosity,
  log_pressure, log_radius, log_mass, log_temperature,
  convective_flag, num_zones, ierr)`. Body is the former loop text
  under the renames `j -> zone_index`, `star%<arr> -> dummy of the same
  array`, `num_species -> species_end_index`, `lower/target/upper ->
  mass_interp_x0/x1/x2`, `envelope_interp_fraction ->
  interp_fraction`, `do k = 4,... -> do k = i_he3,...` (`i_he3 = 4`).
  The one differing comparison is chosen by
  `accept_point_at_fit_mass` (`.le.` when true, `.lt.` when false;
  operands identical). Degenerate mass interval returns `ierr = 1`
  without printing.
- `src/core/rebuild_envelope.f90`: loop replaced by the call with
  `.false.`; on `ierr /= 0` prints its former `write(*,*)` message and
  returns (was `stop 9998` -> ierr path, unchanged text). Locals
  `env_point_index, k, mass_interp_x0/x1/x2, interp_fraction` removed.
- `src/core/read_starting_model.f90` (`rescale_and_refit_envelope`):
  loop replaced by the call with `.true.`; on `ierr /= 0` prints its
  former `write(run_log_unit,*)` message and returns. Locals
  `env_point_index, lower/target/upper_mass_coord,
  envelope_interp_fraction` removed.
- `src/deps.mk` regenerated.

### 2. Profile-column registry, core side

The brief's premise was slightly off: `compute_seismic_columns` already
reads through `ext_profile_value`; the actual duplication is between
`profile_value`/`ext_profile_value` and `build_pulse_points` (and
`compute_seismic_columns`' gravity). Five private module functions in
`src/core/stitched_model.f90` (placed before `ideal_gas_mu`, each body
token-for-token a former site, each site still reading the array it
read before):

- `sound_speed_cgs(gamma1, log10_p, log10_rho)` -- profile column 58,
  interior / envelope / atmosphere branches.
- `local_gravity_cgs(mass_g, radius_cm) = exp(ln10*cgl)*mass_g/(radius_cm*radius_cm)`
  -- `compute_seismic_columns` gravity and both `build_pulse_points` passes.
- `convective_indicator(is_convective)` -- column 9, interior and envelope.
- `envelope_mass_g(i)` -- ext column 2 (then `/star%solar_mass_cgs` as
  before) and the pulse mass of envelope points.
- `atmosphere_radius_cm(i)` -- ext column 3 (inside `log10` as before)
  and the pulse radius of atmosphere points.

Not unified (operand order differs, named in the module comment):
profile column 18 interior gravity vs `local_gravity_cgs`;
`compute_seismic_columns`' N^2 (derivatives divided by dr first) vs
`build_pulse_points` second pass (bracket divided by dr last). io/
duplicates are listed under Deferred 2; io/ untouched.

### 3. Driver-side `star%job%` flips (sidechan Deferred 4)

Applied exactly as proposed, in BOTH drivers (the brief named only
evolve_step; the user message added run_yrec:332):

- `src/core/evolve_step.f90`: new local `logical :: disable_structure_dt_limits`
  (line 54); line 140 `disable_structure_dt_limits = .false.`; the
  `call compute_timestep(...)` at 141-143 gained the trailing actual
  `disable_structure_dt_limits`; line 144
  `if (disable_structure_dt_limits) star%job%use_structure_dt_limits = .false.`.
- `src/core/run_yrec.f90`: same pattern -- local at line 53, pre-set at
  336, call 337-341 with the trailing actual, flip at 342.

For mixburn (callee side, `util/compute_timestep.f90`): add
`logical, intent(out) :: disable_structure_dt_limits` as the new last
dummy (after `log_teff`), set `.false.` on entry and `.true.` where
`structure_limits_active` is cleared, delete the
`star%job%use_structure_dt_limits = .false.` write. `compute_timestep`
is an external procedure with no interface block, so the extra
trailing actual is ignored until the dummy exists; the driver flip is
a no-op on the explicit `.false.`. Merge order: core alone is safe;
mixburn's callee change must merge together with or after core.

The `massloss` half of Deferred 4 (call at `core/evolve_step.f90:233`)
was NOT applied: `wind/massloss.f90`, where the `disk_exhausted` flip
now lives, has no wave-3 owner (mixburn owns `wind/mdot.f90` only).
See Deferred 3.

### 4. `surfbc` dummy `luminosity_linear -> log_luminosity`

`src/core/surfbc.f90` only (dummy declaration, the two uses at the
`tri_logl` assignments, and the header comment). No caller change.

### 5. Geometry helpers per site

Six `log10 R(L, Teff)` sites in core/ come in two operand orders (A and
B below) plus one outlier; they are NOT unified. Where two or more
sites were token-identical they now share one function; helper
headers point at the other order so the R6 decision is visible.

- `src/core/envint_lib.f90`: new `surface_log10_radius_cm(log_luminosity_lsun, log_teff)`
  (order A) and `ambient_log10_teff(log_teff)` (spot formula, real*4
  literals). Used at `rebuild_envelope.f90:103,131`,
  `surfbc.f90:179,186`, `stitched_model.f90:142` (radius only),
  `read_starting_model.f90:739` (spot only).
- `src/core/observables_lib.f90`: new public `log_r_surface_cm(...)`
  (order B), used in `compute_surface_globals` (:282) and
  `run_yrec.f90:224` (new `use observables_lib, only: log_r_surface_cm`
  at :36; deps.mk regenerated).

Sites left as they were: `read_starting_model.f90:718` (outlier,
uses `star%solar_luminosity_cgs`) and `stitched_model.f90:147-148`
(spot formula with d0 literals). Full list under Deferred 4.

### 6. Wave-2 Deferred leftovers (core-only, class A)

- 6a (`5d75705`, sidechan Deferred 2) `src/core/read_starting_model.f90`: deleted
  `star%iovim = -1` (was :688) and appended `zone_index=-1` to the
  `call temperature_gradients(...)` below it. The only remaining
  reader of `star%iovim` is the optional-argument fallback at
  `mixing/temperature_gradients.f90:104`, which now has no caller
  that relies on it.
- 6b (`1530b2b`, henyey Deferred 3) `print_flag` `intent(inout) -> intent(in)` in
  `src/core/envint_lib.f90` (atm_get) and `src/core/envint_kernel.f90:88`
  (integrate_envelope_atmosphere); it is only read, at
  `envint_kernel.f90:689`. `make clean` before gate1 (signature change).
- henyey Deferred 1 (`star%converged_zone` member in state/) and 2
  (dt sign -> flag, cross-domain) are not core-only; not touched.

## Deferred (exact proposals)

1. **Envelope refit, remaining differences** between
   `core/rebuild_envelope.f90` and `core/read_starting_model.f90`
   (`rescale_and_refit_envelope`) -- all outside the shared loop:
   - surface radius: `read_starting_model.f90:718`
     `0.5d0*(star%log_L + star%solar_luminosity_cgs - 4.0d0*star%log_Teff - c4pil - csigl)`
     uses the LINEAR solar luminosity (3.8e33) where every other site
     uses `star%log10_solar_luminosity` (audit core.md bug 1). Fixing
     it changes numbers -> R6.
   - env_struct shift: `rebuild_envelope.f90:167` drops points closer
     than `1.0D-10` in log mass to the last interior shell;
     read_starting_model shifts unconditionally.
   - spot-Teff gate: `convective_flag(num_zones)` (rebuild_envelope:129)
     vs `star%envelope_cz_bottom_index.eq.star%nz`
     (read_starting_model:735).
   - `.lt.` (rebuild_envelope) vs `.le.` (read_starting_model) for an
     envelope point exactly at the fitting mass -- now the
     `accept_point_at_fit_mass` argument; choosing one is a decision.
   - after the loop rebuild_envelope rebuilds shell masses and runs the
     three `walpcz` rotation branches + `omega_from_j`;
     read_starting_model sets solid-body `j_rot` and leaves shell
     masses to `build_shell_masses`.
   - the degenerate-interval and summary prints go to different units
     (`*` vs `run_log_unit`) with different text.

2. **Profile-column expressions duplicated in io/** (io untouched; the
   core kernels in `stitched_model.f90` could be exported through
   `stitched_model_lib` once io agrees):
   - `io/write_gyre_pulse.f90:66` `grav_const_cgs = exp(ln10*cgl)` and
     `:67-70` gravity `grav_const_cgs*mass/(radius_cm*radius_cm)`
     = `local_gravity_cgs` (core `stitched_model.f90`).
   - `io/write_fgong_pulse.f90:46` same `exp(ln10*cgl)` and the same
     gravity at `:96-97`.
   - `io/write_gyre_pulse.f90:49` `exp(ln10*log_radius(num_shells))`,
     `:68,:93` `exp(ln10*log_radius(i))`, `:96-98` the same
     `exp(ln10*log_*)` unit conversions for P, T, rho that
     `build_pulse_points` does.
   - `io/write_gyre_pulse.f90:74` `exp(ln10*(log_density(i) - log_pressure(i)))`
     -- reciprocal of the `sound_speed_cgs` argument (operand order
     reversed; not identical).
   - `io/yrec_output.f90:142,218` `exp(ln10*star%log_L)*star%solar_luminosity_cgs`
     (once via `exp10`, once via `exp(ln10*...)` -- not identical to
     each other), `:216` `exp(ln10*star%log_total_mass)`, `:217`
     `exp(ln10*(star%log_R_surface + star%log10_solar_radius))`.

3. **massloss driver half of sidechan Deferred 4.** When
   `wind/massloss.f90` gets an owner (its :283
   `if(disk_exhausted) star%job%use_mass_accretion = .false.` is the
   callee write): add a `logical :: disk_exhausted` local in
   `core/evolve_step.f90`, set `.false.` before the `call massloss(...)`
   at `:233`, pass it as the trailing-but-one actual (before `ierr`,
   as sidechan proposed; massloss is external with no interface block,
   so unlike the trailing compute_timestep actual this one must land
   TOGETHER with the callee change), and after the `if (ierr /= 0)
   return` add `if (disk_exhausted) star%job%use_mass_accretion = .false.`.

4. **log10 R(L, Teff) sites (R6 decision: pick one operand order).**
   All evaluate `L = 4 pi R^2 sigma Teff^4`:
   - order A `0.5d0*(logL + star%log10_solar_luminosity - 4.0d0*logTeff - c4pil - csigl)`
     -- now `envint_lib::surface_log10_radius_cm`, sites
     `core/rebuild_envelope.f90:103` (`log_luminosity_lsun, log_teff`),
     `core/surfbc.f90:179` (`tri_logl(i), log10_teff`),
     `core/stitched_model.f90:142` (`star%log_L, star%log_Teff`).
   - order B `0.5d0*(logL + star%log10_solar_luminosity - c4pil - csigl - 4.0d0*logTeff)`
     -- now `observables_lib::log_r_surface_cm`, sites
     `core/observables_lib.f90:282` (`compute_surface_globals`) and
     `core/run_yrec.f90:224` (then `- star%log10_solar_radius`).
     `wind/wind_lib.f90:22 log10_radius_from_l_teff` is this same
     order (`4.d0`); once the order is decided, one of the two could
     go.
   - outlier `core/read_starting_model.f90:718`
     `0.5d0*(star%log_L + star%solar_luminosity_cgs - 4.0d0*star%log_Teff - c4pil - csigl)`
     (bug; see Deferred 1).
   Spot ("ambient") Teff sites:
   - `log_teff - 0.25*log10(star%ctrl%spot_filling_factor * pow(star%ctrl%spot_temp_contrast, 4.0) + 1.0 - star%ctrl%spot_filling_factor)`
     (real*4 literals) -- now `envint_lib::ambient_log10_teff`, sites
     `core/rebuild_envelope.f90:131`, `core/surfbc.f90:186`,
     `core/read_starting_model.f90:739`;
     `mixing/temperature_gradients.f90:181-182` is the same text (not
     mine).
   - `core/stitched_model.f90:147-148` same formula with `0.25d0`,
     `4.0d0`, `1.0d0` -- different literal kinds (the real*4 `0.25`,
     `1.0` are exact in double, `pow(x, 4.0)` goes through
     `pow_r_sp -> pow_r(x, dble(4.0))`, so the values are probably
     equal, but that is a numbers claim, not a token identity):
     making it use `ambient_log10_teff` is an R6 decision.

5. **Retire `star%iovim`** (after 6a no caller writes it): delete the
   member at `state/star_info_lib.f90:308` (and the comment at :288),
   make `zone_index` non-optional in
   `mixing/temperature_gradients.f90:49,95` and drop the fallback at
   `:96-105` (plus the W2 comment at :33-39). state/ and mixing/ are
   not mine.

6. From henyey-wave2, still open and not core-only: Deferred 1
   (`star%converged_zone` deletion in state/) and Deferred 2 (dt sign
   -> `age_this_model` flag; crosses rotation/, mixing/, wind/).

## Reverted (changed numbers)

None. Every gate1 run was IDENTICAL on the first build.

## Skipped

- Unifying the two log10-radius operand orders, or the d0-literal spot
  site -- numbers-changing (R6).
- `read_starting_model.f90:718` `solar_luminosity_cgs` bug -- R6.
- massloss driver edit -- no wave-3 owner for `wind/massloss.f90`.
- io/ profile duplicates -- not my files (listed above).
- `star%iovim` member / `temperature_gradients` fallback removal --
  state/ and mixing/ are not mine.
- henyey Deferred 1 and 2 -- not core-only.

## What the audit / earlier reports got wrong or understated

- Brief item 3 named only `core/evolve_step.f90`; `core/run_yrec.f90`
  has an identical `compute_timestep` call (:337) that needed the same
  edit (the user message caught it). Sidechan's Deferred 4 line
  numbers (`evolve_step.f90:141/235`) are now `:141` and `:233`.
- Sidechan's Deferred 4 assumes a wave-3 owner for `wind/massloss.f90`;
  there is none (mixburn has `wind/mdot.f90` only), so the massloss
  half cannot be completed this wave.
- Brief item 2: `compute_seismic_columns` does not duplicate
  `profile_value` -- it already reads through `ext_profile_value`;
  the duplicated expressions are in `build_pulse_points`.
- core-wave1 described the two refits as differing in "two" places;
  they differ in six (Deferred 1). Only the fill loop is shared.
- Brief item 5 said "six sites differ in operand order"; in fact they
  form two token-identical groups (3 + 2) plus one outlier, so two
  helpers were possible without changing any order.

## Verification

- gate1: IDENTICAL for every commit (`/Applications/YREC-wt/logs/core.make.log`
  shows only the pre-existing stitched_model `-fmax-stack-var-size`
  and macOS linker warnings).
- Full 37-case pin selection after `1530b2b`: `37 passed, 272 deselected in 1598.93s`, `PINS EXIT 0` (`/Applications/YREC-wt/logs/core.pins.w3.log`).
- Aux battery: `22 passed in 575.25s`, `AUX EXIT 0` (`/Applications/YREC-wt/logs/core.aux.w3.log`).

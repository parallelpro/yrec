# Readability wave 2 -- henyey (core solver / envelope integration)

Worktree `/Applications/YREC-wt/henyey`, branch `rs/henyey` (based on
yrec-modern after the wave-2a merge). Three commits:

| hash | summary |
|---|---|
| `20e0871` | drop the four unread henyey flag relays (item 1) |
| `88dedaf` | remove kenv/katm; atm_get optional vertex group and step config (items 2, 3) |
| `35838b0` | name the fcorr0 sign mode; comment the dt sign, the run-log rewind, the two luminosity renormalizations (items 4, 5, 6) |

Verification: gate1 byte-identical before each commit; full 37-case pin
selection: 37 passed, PINS EXIT 0 (25 min); aux battery: 22 passed, AUX EXIT 0.
`make clean` was run before the first gate and again after the
module-procedure signature changes of items 2/3. `src/deps.mk` needed no
regeneration (no source files, modules or `use` lines added or removed;
gate1's `--check` reports it up to date).

## Items

### 1. Henyey flag relay -> locals (core 3) -- done (`20e0871`)

Shown unread: a whole-tree grep for `in_atmosphere`, `want_derivatives`,
`mixing_active`, `conductive_opacity_flag` restricted to the henyey chain
finds only the four assignments in henyey_coefficients, the two dummy
lists / two call lines (henyey_coefficients <- henyey_iterate <-
evolve_step), and the evolve_step reset block. `mixing_active` had no
reader at all. `star%converged_zone(im) = conductive_opacity_flag`
(always `.true.`) was the only store into `converged_zone`, and
`converged_zone` is read nowhere in `src/`, `tools/`, `tests/` or
`pyyrec/` -- so the store was removed rather than kept.

- `core/henyey_coefficients.f90`: `want_derivatives = .true.` and
  `in_atmosphere = .false.` are now `logical, parameter` locals (same
  values, same kind) fed to the same `intent(in)` positions of
  `eos_get`/`temperature_gradients`; the other two flags, the
  `converged_zone` store and the four dummies are gone. The
  `call temperature_gradients(...)` lines and `star%iovim = im` were not
  touched (sidechan's).
- `core/henyey_iterate.f90`: four `intent(out)` dummies dropped from the
  signature; call to henyey_coefficients shortened.
- `core/evolve_step.f90`: the four SAVEd locals and their
  `evolve_step_reset_pending` reset lines removed; henyey_iterate call
  shortened.
- `core/read_starting_model.f90`: the write-only local
  `local_conductive_opacity_flag` removed (same family; it was assigned
  once and never read).

### 2. Delete kenv/katm (core 4) -- done (`88dedaf`)

Traced: `kenv`/`katm` (dummy names `env_call_count`/`atm_call_count`,
`call_count` in numerics) were zeroed in henyey_iterate,
rebuild_envelope, read_starting_model and build_stitched_model,
incremented once per call in `atmosphere_derivs`/`envelope_derivs`, and
otherwise only passed through surfbc -> atm_get ->
integrate_envelope_atmosphere -> bsstep -> mmid -> deriv. No read
anywhere (the only former consumer, a format in the old envint, is
gone). Removed from all six signatures and the four drivers:

- `core/atmosphere_derivs.f90`, `core/envelope_derivs.f90`: dummy,
  declaration and `x = x + 1` removed.
- `numerics/numerics_lib.f90` (this item only): `mmid` and `bsstep`
  dummy lists / declarations and their `deriv`/`mmid` calls; the two
  stale header pointers "atm/atm_lib.f90" now say
  `core/envint_kernel.f90`. Nothing else in numerics touched.
- `core/envint_kernel.f90`: dummy/declaration and the six integrand /
  bsstep calls; no arithmetic touched.
- `core/envint_lib.f90`, `core/surfbc.f90`, `core/henyey_iterate.f90`,
  `core/rebuild_envelope.f90`, `core/read_starting_model.f90`,
  `core/stitched_model.f90`: counters and their zeroing removed.

Note: `core/atmosphere_derivs.f90` is not in the brief's explicit file
list but is unavoidable for this item (it is the fourth integrand
signature); sidechan's brief assigns all of `src/core/*` to henyey, so
it was edited here.

### 3. `atm_get` signature (core 5) -- done (`88dedaf`)

`envint_config` already existed in `envint_kernel` (the phase-B split);
what the drivers still did was overwrite `star%job%{atm,env}_step_*`
around the call so atm_get would snapshot their values into it.

- `core/envint_lib.f90`: new public `type envint_step_config`
  (`step_begin, step_min, step_max`) and `pure function
  fixed_envint_step(h)` (all three = h, the only pattern used). atm_get
  keeps its mandatory arguments positional, then `ierr`, then OPTIONAL
  keyword arguments: the fit-vertex group (`vertex_index`,
  `stored_envelope_state`, `stored_vertex_index`, `vtx_logp/logr/logt`)
  and `atm_steps`, `env_steps`. `cfg%*_step_*` come from the optional
  when present, else from `star%job` as before. When the vertex group is
  absent, zero-initialised scratch is passed to the kernel via a
  contained `run_kernel` wrapper (so the kernel call is written once).
- `core/surfbc.f90`: passes the vertex group by keyword (all six).
- `core/rebuild_envelope.f90`, `core/read_starting_model.f90`: the
  save/override/restore of `star%job%env_step_*`, the four dummy arrays,
  `ixx_flag`/`vertex_index = 0`, `atm_get_unused_flag = 0` are gone;
  `env_steps = fixed_envint_step(star%ctrl%chi_grid_scale(ichi_dp_env_max))`
  is passed instead.
- `core/stitched_model.f90`: same for both triples (`atm_step_size`,
  `envelope_step_size`); its six saved copies and four dummies gone.

Value fidelity: the double-precision step values reach
`cfg%env_step_*` through the type unchanged -- exactly the values the
override used to write into `star%job` and atm_get used to copy out.
Vertex scratch: the three non-surfbc drivers pass
`save_boundary_flag = .false.`; the kernel writes `vtx_*` /
`stored_*` only under that flag, and otherwise only executes
`if (stored_vertex_index.eq.vertex_index) stored_vertex_index = 0`
(0 vs 0, a no-op) -- identical to what the old `ixx = 0`, `idum = 0`
and dummy arrays produced.

Error-return paths: rebuild_envelope and read_starting_model returned
before their restore, leaving `star%job%env_step_*` overridden. With
the config approach `star%job` is never written. Shown harmless:
`star%job%env_step_*` is read only by atm_get itself (which builds a
fresh cfg every call) and by `seed_controls_buffer`
(`state/controls_sync_lib.f90`), which runs only inside `read_controls`
-- not reached after those error returns (the run terminates), and on a
multi-run re-entry `yrec_reset` restores `star = star0` first.
stitched_model restored before its error check, so nothing changes
there either.

### 4. Sign-encoded flags (core 14) -- half done (`35838b0`)

- **fcorr0 sign (done).** `core/henyey_iterate.f90`: new local
  `underrelax_ramp_active = star%ctrl%fcorr0.gt.0.0d0` set once on
  entry (`fcorr0` is a `/physics/` NAMELIST control; nothing writes it
  during the solve -- grep: only `read_controls`/`controls_sync_lib`),
  used in the ramp test `if (underrelax_ramp_active) star%job%fcorr =
  dmin1(...)`. The comparison expression is token-identical. The
  magnitude use `dabs(star%ctrl%fcorr0)` in evolve_step stays, with a
  comment naming the encoding. Not a namelist-facing rename.
- **dt sign (deferred, see below).** A named logical cannot replace it
  byte-identically: the negative `star%dt` is passed as
  `full_timestep` to `evolve_angular_momentum`
  (`rotation/evolve_angular_momentum.f90:164` feeds it to `matt_wind`
  before the `.gt.0` test at :187 and the sub-timestep arithmetic at
  :196-224), to `mix`, and to henyey_iterate's `delta_time` where
  `compute_entropy_term = delta_time.gt.0.0d0` and the
  `star%dt.le.0` first-guess test read the sign. A comment at the
  assignment (`core/evolve_step.f90:316-322`) now states the encoding
  and its readers.

### 5. Small core items (core 15) -- done (`35838b0`)

- `use_correct_gradients = .true.`: already deleted in wave 1 (R1,
  `reports/core-wave1.md` "Dead code"); nothing left to do.
- Retired-machinery comment in `update_output_flags_for_step`
  (`core/evolve_step.f90`): rewritten to say what the routine does now
  (rewind `run_log_unit` -- log_output_file / run.log -- when the lrwsh
  control `rewind_short_file` is set).
- Duplicated pressure-scale-height formula: already one contained
  function `envelope_overshoot_depth()` since wave 1 (R3); nothing left
  to do.
- Extra, same spirit: the two `do j = 1,15` copies of the
  start-of-step composition (`orig_composition`/`xa_start`, both
  declared `(n_species_extended,json)`) now loop to
  `n_species_extended` (= 15).

### 6. `luminosity_breakdown` renormalisation (core-wave1 Deferred 9) -- done (`35838b0`)

Not merged. One cross-reference comment at each site:
`core/henyey_iterate.f90` (per-iteration, guards `total > 0`) and
`core/observables_lib.f90::renormalize_luminosity_breakdown`
(post-convergence, no guard). They are numerically different routines
(the guard) run at different times, so they stay separate. (Note
`core/observables_lib.f90` is not in the explicit list but is core/*,
assigned to henyey by sidechan's brief; only the comment was added.)

## Deferred

1. **`star%converged_zone(json)` member deletion** --
   `src/state/star_info_lib.f90:355` (`logical :: converged_zone(json)`)
   and its mention at :346 (`converged_zone (LOCONS)`). After item 1 it
   is neither written nor read anywhere. `state/` is not mine; delete
   the member and the comment word.

2. **dt sign -> explicit flag.** Exact proposal (all Class-A if done
   together, but crosses into rotation/, mixing/, and the model file):
   add `logical :: age_this_model` to `star_info` (job section), set in
   `core/evolve_step.f90:322` (`age_this_model = evolve_model_flag`)
   instead of `star%dt = -dabs(star%dt)`, and replace the sign readers:
   `core/evolve_step.f90:330` (`star%dt.le.0.0D0` ->
   `.not.age_this_model`), `:377` (`star%dt.gt.0.0D0` ->
   `age_this_model`), `:134` (`star%dt = dabs(star%dt)` -> drop),
   `core/henyey_coefficients.f90` `compute_entropy_term =
   delta_time.gt.0.0d0` (-> a new logical dummy), `rotation/
   evolve_angular_momentum.f90:187` (`full_timestep.gt.0.0D0`) and the
   `matt_wind` call at :164 which today receives the NEGATIVE
   full_timestep when not aging -- the wind torque is multiplied by it,
   so that call would have to be gated (or keep receiving a signed
   value) to stay byte-identical; likewise `mix` (:244) and `massloss`
   (:227) receive `star%dt` and would need auditing for sign use. Do it
   as one cross-domain change with the pins, not piecemeal.

3. **`print_flag` intent(inout) in atm_get/integrate_envelope_atmosphere**
   (`core/envint_lib.f90`, `core/envint_kernel.f90:88`): it is only
   read; could be `intent(in)` (would let stitched_model pass a
   literal). Not done -- not in the brief.

## Reverted (changed numbers)

None.

## Skipped

None of the six items skipped.

## Notes on the audit / wave-1 reports

- core.md item 15 lists `use_correct_gradients` and the pressure-scale-
  height duplication as open; both were already closed by wave 1
  (`reports/core-wave1.md`, R1 dead code and R3 `envelope_overshoot_depth`).
  Only the `update_output_flags_for_step` comment was still open.
- core.md item 5 describes four `atm_get_dummy*` arrays and the
  `env_step_*` smuggling; build_stitched_model additionally smuggled
  the `atm_step_*` triple (via `star%ctrl%atm_step_size`) -- covered by
  the same `atm_steps` optional.
- core.md item 3 says the flags are relayed into "SAVEd driver locals
  nobody reads"; correct, and `star%converged_zone` was additionally
  write-only across the whole tree (tools/tests/pyyrec included).

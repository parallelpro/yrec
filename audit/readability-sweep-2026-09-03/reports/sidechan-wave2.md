# Readability wave 2 -- sidechan (batch R4, explicit data flow)

Worktree `/Applications/YREC-wt/sidechan`, branch `rs/sidechan`, base `2ff54e0`
(yrec-modern after wave 2a).

## Commits

| hash | summary |
|---|---|
| `76b08ec` | Readability W2 (sidechan): rezone keeps its spline scratch in locals (item 1) |
| `2854d8d` | Readability W2 (sidechan): temperature_gradients takes zone_index (item 2) |
| `1e7d1eb` | Readability W2 (sidechan): pass the DEL gradient into microdiff explicitly (item 3) |
| `525f75c` | Readability W2 (sidechan): set the Bahcall-Loeb scales in one place; drop dead un-scaling (item 4, byte-safe part) |
| `41408f1` | Readability W2 (sidechan): report job-flag flips from mdot / compute_timestep to the caller (item 5) |
| `0400aed` | Readability W2 (sidechan): share the ion-mean-weight excess block between mu and eqstat (item 6) |
| `f9b3869` | Readability W2 (sidechan): am_advection_diffusion_coeffs takes timestep and wind loss as inputs (item 7) |

Verification: gate1 IDENTICAL before every commit (`make clean` first for
`2854d8d`, `1e7d1eb`, `525f75c`: module-procedure signature / new module).
On the final tree (`f9b3869`): the full 37-case pin selection ended `PINS EXIT 0` (37 passed, 25:40,
`logs/sidechan.pins.log`) and the aux battery `AUX EXIT 0` (22 passed,
`logs/sidechan.aux.log`). Only source files and
`src/deps.mk` committed; no `test_*.short` left behind.

## Item 1 -- `rezone` local spline buffers (`76b08ec`)

`setup/rezone.f90`. The four `star%` members used as "dummy arrays" became
host-level locals of `rezone`, shared by its contained phases:

| was | now | role |
|---|---|---|
| `star%old_shell_mass` | `new_log_mass` | the NEW run of log-mass points (HSS in hpoint.f); the member name was a misnomer |
| `star%logRho_start` | `compact_log_mass` | buffer of the too-close-point deletion pass |
| `star%logP_start` | `regrid_in` | y-table copy in `regrid_in_place` |
| `star%logT_start` | `regrid_out` | osplin output in `regrid_in_place` |

What was checked:

- `old_shell_mass` has no reader or writer outside rezone (whole-tree grep);
  the `star_info_lib` member is now dead (state's file: Deferred 1).
- The osplin calls that write `star%logP_start / logR_start /
  luminosity_lsun_start / logT_start / logRho_start` for the P, R, L, T, rho
  run are NOT scratch and are unchanged: after the transfer loop they hold
  the start-of-step model on the new grid (what `evolve_step`'s no-rezone
  branch stores explicitly at `core/evolve_step.f90:270-308`) and are read
  by `util/timestep_limit_structure.f90:47-72`,
  `rotation/mid_timestep_model.f90:87-96`, `core/burn_lib.f90:2059-2083,
  2876-2882` over `1..star%nz`. Those calls overwrite elements `1..new nz`
  after the last scratch use, so every element a later reader touches is
  the same as before; only elements above the new nz differ (stale scratch
  before, pre-rezone values now) and nothing reads them.
- Every element of each new local is written before it is read on every
  path (`new_log_mass(1)` and `2..new_num_zones` in `assign_new_points`; the
  compaction buffer `1..j`; `regrid_in 1..old nz`; `regrid_out 1..new nz`
  by osplin), so `-finit-local-zero` zeros are never observed.
- `star%logP_start(1) = star%logP(1)` (~line 347) untouched: it is the
  start-of-step save, not scratch.

## Item 2 -- `temperature_gradients` `zone_index` argument (`2854d8d`)

`mixing/temperature_gradients.f90`: new trailing `integer, intent(in),
optional :: zone_index`; the ladov adiabatic-overshoot gate now compares
`overshoot_zone_index`, which is `zone_index` when present and `star%iovim`
otherwise. Gate expression unchanged apart from the operand name.

Call sites carry exactly what `star%iovim` held at each call, and the three
`star%iovim = ...` writes whose only reader was the following call are
deleted (only those lines changed in the shared core files):

| call site | passes | was |
|---|---|---|
| `core/henyey_coefficients.f90` | `zone_index=im` | `star%iovim = im` |
| `core/shell_physics.f90` | `zone_index=im` | `star%iovim = im` |
| `core/envelope_derivs.f90` | `zone_index=-1` | `star%iovim = -1` (kernel path stays explicit-args) |
| `mixing/semiconvection.f90` (3 calls) | `zone_index=-1` | nothing written (stale value read) |
| `core/read_starting_model.f90:697` | (absent -> reads `star%iovim`) | `star%iovim = -1` still written at :696 |

Deviation from the brief, semiconvection: the brief said "pass `star%iovim`
itself so today's value flows". That would not have preserved today's value
anyway, because deleting the three core writes changes what `star%iovim`
holds when semiconvection runs. Instead: semiconvection reads only
`radiative_gradient` from the call (`actual_gradient`, `is_convective`,
`convective_velocity`, `dgrad_*` are declared and never read -- grep), and
the gate only rewrites `actual_gradient`, so the index cannot reach any
output. Passed `-1` with a comment at the first call saying exactly this.
Byte-identical confirmed by gate1 + full pins.

Why `optional`: `core/read_starting_model.f90:696-697` (henyey's file) also
calls `temperature_gradients`, which the brief did not list; a mandatory
argument would not have built. The fallback reads `star%iovim` exactly as
before at that one site. Removal is Deferred 2.

## Item 3 -- `rotmix` gradient argument (`1e7d1eb`)

Added `del_grad(json)` (intent(in), placed after `log_temperature`) to
`microdiff`, `microdiff_mte`, `gravitational_settling`,
`gravitational_settling_setup` and replaced every `star%gradT` read under
`rotation/microdiff/` with it: `microdiff_mte.f90` 5 sites (the two
`lagrange4(facinterp, star%gradT(k0:k0+3))` calls, the linear
interpolation at the first equal-grid point, and the `zone_begin`/`zone_end`
end values), `gravitational_settling_setup.f90` 8 sites. Callers:

- `mixing/rotmix.f90`: passes `mix_scr%delm`; the save/overwrite loop
  (`del_grad2_save(i) = star%gradT(i); star%gradT(i) = mix_scr%delm(i)`),
  the restore loop and the `del_grad2_save` local are gone.
- `mixing/mix.f90`: passes `star%gradT` (what those readers saw before).

Checked that nothing else on the settling path reads `star%gradT`
(`implicit_diffusion_coeffs`, `tridiag_gs`, `lax_wendroff_step*`,
`model_to_equal`, `equal_to_model`, `thoul_diffusion`, `microdiff_run`,
`microdiff_coefficients`: none). `check_boundaries.py` already allowlists
`microdiff` / `gravitational_settling` as mixing -> rotation entries.

## Item 4 -- microdiff `*_bl` local copies (`525f75c`, partial; main part Deferred 3)

Checked the restore first, as instructed. It is NOT exact:

```
microdiff_setup.f90:161-165                       microdiff_etm.f90:168-174 / equal_to_model.f90:148-154
enclosed_mass(i)=enclosed_mass(i)*bl_mass_scale   enclosed_mass(i)=enclosed_mass(i)/bl_mass_scale
dlnp_dr(i)=dlnp_dr(i)/bl_radius_scale             dlnp_dr(i)=dlnp_dr(i)*bl_radius_scale
timestep=timestep/bl_time_scale                   timestep=timestep*bl_time_scale
total_mass=total_mass*bl_mass_scale               total_mass=total_mass/bl_mass_scale
```

with `bl_radius_scale = 1/6.9598d10`, `bl_mass_scale = 1d-2*bl_radius_scale**3`,
`bl_time_scale = 2.7d13*3.1558d7` -- none a power of two. A 1e5-sample
double-precision check with those constants: `m*s/s /= m` for 9.3 % of
values, `t/s*s /= t` for 13.0 %, `d/s*s /= d` for 0.8 %. And the
post-restore values ARE read: `mix.f90` passes `star%m` as `enclosed_mass`
(the model's mass array, perturbed by up to 1 ulp on ~9 % of zones every
settling step and carried forward), and `rotmix.f90` reuses
`enclosed_mass`, `dlnp_dr_settling`, `settling_dt`, `total_mass` across
`num_settling_substeps`. Keeping the callers' arrays intent(in) would
therefore change numbers -> class B, deferred with this evidence (Deferred 3).

Done, byte-safe:

- New module `rotation/microdiff/bahcall_loeb_units.f90`
  (`bahcall_loeb_units_lib`, `set_bahcall_loeb_scales`) holds the four
  `star%bl_*_scale` assignments once; `microdiff_setup` and
  `gravitational_settling_setup` call it and lose their private copies of
  the constants (identical values: `rsun_cgs_legacy`, `3.1558d7`). Same
  expressions, same order. `src/deps.mk` regenerated.
- The `radius_bl`/`temperature_bl` "restores" in `microdiff_etm` and
  `equal_to_model` were dead: those are locals of `microdiff` /
  `gravitational_settling`, which `return` right after the call. Deleted;
  `radius_bl` is now intent(in) in both and the unused `temperature_bl`
  dummy is dropped from both (and from the two call lines). This is the
  "item 4 restore code only" scope in `io/equal_to_model.f90`.

## Item 5 -- physics routines flipping `star%job%` flags (`41408f1`, prepare-only)

- `wind/mdot.f90`: no longer writes `star%job%use_mass_accretion`. New
  trailing-but-one dummy `logical, intent(out) :: disk_exhausted`
  (initialised `.false.` next to `ierr = 0`, so the `ierr = 1` early
  return leaves it defined; set `.true.` where the local
  `disk_exhausted_flag` was). `wind/massloss.f90` -- its only caller,
  which returns immediately after -- does
  `if (disk_exhausted) star%job%use_mass_accretion = .false.`. Nothing in
  `mdot` read the flag after the old write, so this is a pure move.
- `util/compute_timestep.f90`: the decision is now a local
  `structure_limits_active` (initialised from
  `star%job%use_structure_dt_limits`, cleared in the log Tc > 7.1 /
  L_grav < 0 block) and the two later reads (`timestep_limit_structure`
  gate, envelope-triangle gate) use the local. The
  `star%job%use_structure_dt_limits = .false.` write stays in the block, as
  the brief asked, because both callers are core (`evolve_step.f90:141`,
  `run_yrec.f90:332`). Driver edit: Deferred 4.

## Item 6 -- eos `ion_mean_weight_excess` (`0400aed`)

Compared the two blocks statement for statement: `mu.f90:35-46` and
`eqstat.f90:422-439` (pre-edit numbering) are identical (`dfx1`, `dfx12`,
threshold test, else `dfx1*w(1)`, `dfx12*w(3)`, `dfx4`, `mu = amuenv + dfx1
+ dfx4 + dfx12`), with the same `inverse_atomic_weights` DATA in both files.
eqstat's Saha setup sits AFTER the arithmetic in each branch (`saha =
fxenv` in the envelope branch; `amu_inverse = 1/mu` then the fractions in
the other), so hoisting the arithmetic into a call reorders nothing.

Added `ion_mean_weight_excess(x, z, mu_ion_inv, dfx1, dfx12, dfx4,
use_envelope)` to `eos/eos_mixture_lib.f90` (module now has a `contains`;
`eos_mix` stays public, procedure exported explicitly). `mu.f90` calls it
and reads only `mu_ion_inv` (its `inverse_atomic_weights` copy is gone);
`eqstat2` calls it and branches on `use_envelope_mixture`, keeping only the
Saha statements. `dfx4` is defined only in the non-envelope branch, which
is the only place either caller reads it. `eqstat`'s
`inverse_atomic_weights` table stays: it is also used for the electron
mean weight (`eqstat.f90:467-468, 535-536`) -- I first deleted it and the
build caught it.

## Item 7 -- `am_advection_diffusion_coeffs` (`f9b3869`)

- `json=5000, nmax=8000` shadowing: already gone before this wave (`json`
  from `star_info_lib`, `band_nmax = 8000` from `rotation_scratch_lib`).
  Nothing to do.
- Unused computed locals: a read/write scan of every local found exactly
  four that are written and never read: `omega_prev_medium_iter_avg`
  (3 writes), `omega_mid_prev` (2), `domega_dr_prev` (2),
  `max_omega_change_medium_iter_zone` (2). Deleted with their assignment
  statements; the `theta_iter_idx.le.2` if/else around the
  `omega_prev_medium_iter` update then had identical branches and was
  merged. `total_velocity` (named in the brief) IS read -- it is printed
  in the format-911 diagnostic -- and stays.
- `timestep` mutated-then-restored: now `intent(in)`, with the local
  `substep_timestep = timestep/dfloat(num_substeps)` used at the six
  former read sites. Byte-safe because the restore was `timestep =
  full_timestep`, an exact copy of the entry value (not a multiply/divide
  round trip). The `ierr = 1` early return skipped the restore, but
  `secular_transport`, `evolve_angular_momentum` and `converge_with_rotation`
  all `return` on nonzero ierr and the run stops; `sub_timestep` is a
  local of `evolve_angular_momentum` and is never printed on that path.
- `wind_loss_implicit` leaking: now `intent(in)`, with the local
  `wind_loss_implicit_iter` reset from it at every timestep-cut restart
  and rescaled per theta iteration. Byte-safe because
  `rotation/seculr/secular_transport.f90:339-364` assigns
  `wind_loss_implicit` on every path before each use and nothing reads it
  after the call (it is a local of `secular_transport`).

## Deferred (exact proposals)

1. **state, `state/star_info_lib.f90:337`**: delete
   `double precision :: old_shell_mass(json)` (and its mention in the
   nearby comment). No reader or writer remains (grep: only the historical
   note in `setup/rezone.f90:72`).
2. **henyey/state, `star%iovim` retirement**:
   `core/read_starting_model.f90:696` delete `star%iovim = -1`; at :697-700
   append `, zone_index=-1` to the `call temperature_gradients(...)`
   argument list. Then (mixing, any batch) in
   `mixing/temperature_gradients.f90` drop `optional` from `zone_index`,
   replace the `if (present(zone_index)) ... else overshoot_zone_index =
   star%iovim` block by `overshoot_zone_index = zone_index`, and delete
   the header paragraph about the fallback; then (state) delete `iovim`
   from `state/star_info_lib.f90:308` and the comment at :288, and from
   `io/controls_lib.f90:146`'s comment. Output-neutral: the only remaining
   reader is the fallback, and read_starting_model passes what it wrote.
3. **Item 4 main part -- class B**: making `dlnp_dr`, `enclosed_mass`,
   `timestep`, `total_mass` intent(in) in `microdiff` /
   `gravitational_settling` (local `*_bl` copies in the setups, restore
   loops in `microdiff_etm.f90:167-172` and `io/equal_to_model.f90:146-152`
   deleted) changes numbers: today `star%m` (via `mix.f90:443/449`) and
   rotmix's substep arrays come back perturbed by the inexact `a*s/s`
   round trip (evidence in Item 4 above). Proposal when a class-B window
   opens: exactly that change; expected effect is <= 1 ulp in
   `star%m`, `dlnp_dr`, `settling_dt`, `total_mass` after each settling
   call, which the pins will show as last-digit drift in long runs.
4. **henyey, item 5 driver edits**:
   - `core/evolve_step.f90:235-237` (`massloss` call inside `if
     (evolve_model_flag)`): add a `logical :: disk_exhausted` local,
     pass it as a new trailing-but-one actual (`..., new_atmosphere_fit_needed,
     disk_exhausted, ierr)`), and after the `if (ierr /= 0) return` add
     `if (disk_exhausted) star%job%use_mass_accretion = .false.`. Then in
     `wind/massloss.f90` (mine) turn the local `disk_exhausted` into
     `logical, intent(out) :: disk_exhausted` placed before `ierr` and
     delete the flip at the end of `massloss`.
   - `core/evolve_step.f90:141-143` and `core/run_yrec.f90:332-336`
     (`compute_timestep` calls): add a `logical ::
     disable_structure_dt_limits` local at each site, pass it as a new
     trailing actual, and after each call add `if
     (disable_structure_dt_limits) star%job%use_structure_dt_limits =
     .false.`. Then in `util/compute_timestep.f90` (mine) add `logical,
     intent(out) :: disable_structure_dt_limits`, set it `.false.` at
     entry and `.true.` where `structure_limits_active = .false.` is set,
     and delete the `star%job%use_structure_dt_limits = .false.` line.
     Output-neutral: `compute_timestep` already reads only the local after
     the decision, and nothing between the routine's return and the
     driver's next statement reads the flag (`run_yrec.f90:192/224` only
     save/restore it around calibration cycles).

## Reverted (changed numbers)

None. Every gate1 run reported IDENTICAL on the first build except item 6,
where the first build failed (missing `inverse_atomic_weights` in eqstat,
fixed before any pin ran).

## Skipped

- Item 4 main part: skipped as class B (see Deferred 3), per the brief's own
  instruction to check the restore first.
- Item 2 "pass `star%iovim` itself" for semiconvection: not done as
  written (see Item 2 for why it could not preserve today's value and why
  `-1` is output-neutral).

## What the audit / brief / wave 1 got wrong

- Item 1 (mixwind audit): `star%old_shell_mass` is not the OLD shell mass
  and not the scratch of a spline; it is the NEW run of log-mass points
  that rezone builds (HSS in hpoint.f). And of the four members listed as
  scratch, `logP_start/logT_start/logRho_start` are scratch only in two
  places; their osplin outputs for the P/T/rho run are the real
  start-of-step model on the new grid, read by three later routines.
  Replacing those would have broken timestep control.
- Item 2 (brief): `core/read_starting_model.f90:697` is a fourth caller of
  `temperature_gradients` with its own `star%iovim = -1`; the brief listed
  it only as a writer. It is why `zone_index` had to be optional.
- Item 4 (rotation audit item 9 / Structural 5): the restore is inexact,
  so the "callers' arrays intent(in)" part is not byte-safe -- the brief
  anticipated this. Additionally, `temperature_bl` was never read by
  `microdiff_etm` / `equal_to_model` at all; it existed only to be
  un-scaled.
- Item 7 (rotation audit item 17): the `json`/`nmax` shadowing was already
  fixed in wave 1; `total_velocity` is not unused (it is printed); the
  truly unused locals are the four named above, none of which the audit
  listed.

# Readability wave 1 -- rotation domain

Worktree `/Applications/YREC-wt/rotation`, branch `rs/rotation`, base `0d5c11a`.
Area report: `audit/readability-sweep-2026-09-03/rotation.md`.

| Batch | Commit | Verification |
|-------|--------|--------------|
| R1 comments and dead code | `a4fccb1` | gate1 IDENTICAL; full 37-case pin selection 37 passed, `PINS EXIT 0` |
| R2 named constants / sentinels / renames | `c660ee6` | gate1 IDENTICAL; 37 passed, `PINS EXIT 0` |
| R3 small shared helpers | `f3e841b` | gate1 IDENTICAL; 37 passed, `PINS EXIT 0`; aux battery 22 passed, `AUX EXIT 0` |

Flags throughout: `-O3 -ffp-contract=off -finit-local-zero -Wall`, `USE_HDF5=1`.
No file outside the assigned list was edited except the regenerated `src/deps.mk`.
Nothing under `audit/` or `input/` touched; no push, rebase or branch change.

Note on R1: the R1 commit was made once with foreign edits accidentally applied
to the tree by another agent's script sharing the scratchpad; it was undone
(`reset --soft`), the foreign edits were inverted exactly, and a rebuild showed
every rotation object file byte-identical to the build that had been pinned,
before `a4fccb1` was committed. The pin run covers exactly the committed tree.

## R1 -- comments and dead code (`a4fccb1`, 37 files, +226/-1267)

Dead code was shown dead by whole-`src/` grep (no reader / no caller) or by
control flow (flag never set, branch after a `.false.` hardwire, loop bounds).

- `viscos.f90`: Endal-Sofia/Thomas radiative and molecular viscosity comparison
  values (`viscosity_*_endal_sofia*`) deleted -- assignments were the only
  references. One-line comment left.
- `mid_timestep_model.f90`: 155-line `new_cz_detected` block deleted (flag
  initialised `.false.`, never set). `time_fraction` intent moved to R2.
- `evolve_angular_momentum.f90`: `skip_diffusion_flag` deleted (computed only on
  the `instability_transport_active` else-branch where it is always `.false.`);
  write-only `mix_grads_flag` and its MIXCZ/IMIX comment blocks deleted;
  commented F77 (`CALL PHYSIC`, `CALL SECULR`, `CALL LIBURN*`, `DO 20/WOLD`)
  deleted; stray `continue` deleted; `ierr` declaration moved next to the
  dummies; header text (main.f / "SECULR not yet converted") corrected.
  **BUR-ST branches kept**; the `LBURS = .false.` hardwire now carries a
  comment pointing at `SUMMARY.md`.
- `am_advection_diffusion_coeffs.f90`: `rhs_orig`, `coeff_matrix_orig`,
  `residual_check`, `damping_factor`, `omega_mid_init`, `domega_dr_init`,
  `omega_mid_it`, `domega_dr_it`, `total_delta_angular_momentum_alt`,
  `total_velocity` deleted (write-only); dead `lrossby` branch deleted
  (`lrossby` never set) -- `wind_saturation_threshold` now assigned directly
  from `star%job%wind_saturation_omega` (the live value).
- `check_angular_momentum.f90`: unused dummies `diffusion_velocity`, `zone_min`,
  `zone_max`, `omega_start`, `print_zone_id`, `print_zone_count` removed; sole
  caller `secular_transport` updated. `cut_count.gt.0` branch left (see Skipped).
- `secular_transport.f90`: dead `constant_diffusion_coeff_flag` branch removed
  (flag never `.true.`); the post-loop `diffusion_solve_ok = .true.` now has an
  explicit comment stating that an unconverged final iterate is accepted
  (audit "Possible bugs" #4 -- behaviour unchanged, comment only).
- `diffuse_composition.f90`: `dcomp`, `dcomp2`, `sum_species_orig/updated`
  accumulators deleted (never read); the resulting unused `shell_mass` dummy
  (flagged by `-Wunused-dummy-argument`) removed from the signature; both call
  sites in `diffuse_composition_driver.f90` updated (only caller).
- `diffusion_velocity_scales.f90`: writes to `circ_scr%hle` deleted, member
  `hle` and the never-read `rot_scr%del_grad_diff_prev`, `del_grad_diff_new`,
  `radius_prev` removed from `rotation_scratch_lib.f90` (no readers in `src/`);
  its `use rotation_scratch_lib` dropped (deps.mk regenerated).
- `banded_solver.f90`: 40-line commented F77 solver deleted; unused `json`
  parameter deleted.
- `circulation_velocities.f90`: write-only `v2`, `wmid`; ~100 lines of
  commented F77 deleted.
- `am_transport_grid`, `composition_grid`, `equal_grid_to_model`,
  `compute_quadrupole`, `zahn_coupling_factor`, `check_composition`,
  `am_diffusion_coeffs`, `rotation_shape_factors`, `zone_moments_of_inertia`,
  `shape`, `solid_body_omega`, `equal_spaced_grid`, `enforce_rotation_profile`:
  write-only locals (e.g. `prev_log_*`, `spherical_moment_of_inertia`, `aintt`),
  duplicate `use star_info_lib` lines, blank COMMON-residue blocks, stray
  `continue`, stale header text (qgauss/func, common/quadd/, "1/R**4").
- microdiff: `microdiff.f90` loop-exit index made explicit
  (`species_fraction(3,num_eq_points)`); `microdiff_mte.f90` unreachable
  `print*,'mte line 47'` (loop starts at 2); `microdiff_coefficients.f90`
  write-only `bl_radius_scale_local`, `DATA FGRLI` residue;
  `gravitational_settling.f90` write-only `z_change_first/last`;
  `microdiff_setup.f90` / `gravitational_settling_setup.f90` stray `continue`,
  `(1)` literals, commented old-version lines, "HYDROGEN-FREE" comment on the
  helium check corrected; `lax_wendroff_step2`, `implicit_diffusion_coeffs`,
  `thoul_diffusion` commented FORMAT statements; `microdiff_run`,
  `microdiff_coefficients`, `microdiff_etm`: headers now name the actual
  routines and say what `use_generic_diffusion_vectors` gates.

Post-loop `if (idx > max)` guards after loops that `exit` were kept: they are
semantic (the loop-exhausted case), not goto residue.

## R2 -- named constants, sentinels, renames (`c660ee6`, 14 files, +83/-72)

- `rotation_scratch_lib.f90`: `integer, parameter, public :: band_nmax = 8000`
  replaces the two local `nmax = 8000` in `banded_solver.f90` and
  `am_advection_diffusion_coeffs.f90` (the latter also drops its local
  `json = 5000` for `star_info_lib`'s).
- `am_advection_diffusion_coeffs.f90`: `iter_history_max = 50` for
  `max_omega_change_history` / `max_omega_change_zone_history`;
  `tiny = 1.0d-30` made a parameter.
- `microdiff_mte.f90`: local `half_json = 5000` replaced by `json`.
- `microdiff_setup.f90` (`crsun_bah`, `csecyr_bah`),
  `gravitational_settling_setup.f90` (`solar_radius_bl`, `seconds_per_year_bl`),
  `microdiff_coefficients.f90` (`bl_temp_scale_local`): parameters.
- Rename `fully_convective_flag` -> `settling_skipped_flag` in `microdiff`,
  `microdiff_setup`, `gravitational_settling`, `gravitational_settling_setup`
  (all 20 uses; set for fully-convective, H-exhausted and He-exhausted models).
- `gravitational_settling_setup.f90`: `metal_fraction_total` -> `z_plus_he3_fraction`,
  `iron_fraction` -> `metal_fraction`; `ac_scratch` split into
  `mass_weighted_conc` and `helium_ah_coeff`.
- `mid_timestep_model.f90`: `time_fraction` `intent(inout)` -> `intent(in)`.
- Comment fixes: `zone_moments_of_inertia` header (thin-shell 2/3 for every
  zone under `walpcz /= 0`, then 0.4 solid-sphere for zone 1),
  `gravitational_settling` header (metal settling exists),
  `enforce_rotation_profile` header (no `wcz.f90`; base modes are in
  `omega_from_j`), `implicit_diffusion_coeffs` (correction is a no-op from
  `microdiff_run`, live from `gravitational_settling`).

## R3 -- small shared helpers (`f3e841b`, 8 files)

All helper bodies are token-for-token the former inline expressions.

- `shape.f90`: contained `pure function radau_rk4_step(...)` replaces the two
  RK4 Radau steps. `rho_avg`/`rho_bar_avg` (independent assignments that sat
  between `deta1` and `deta2`) are computed before the call. The r0 Newton
  loops were **not** merged (`.le.` at one site, `.lt.` at the other).
- `microdiff_mte_lib`: `pure function lagrange4(fac, a, k0)` for the ten 1-D
  sites in `microdiff_mte` and two in `microdiff_etm` (new
  `use microdiff_mte_lib, only: lagrange4`; deps.mk regenerated). The 2-D
  `composition(row,k0..)` / `light` sites stay inline (a strided section would
  be needed).
- `circulation_velocities.f90`: contained `log_lagrange4(a, k, i)` for the six
  `exp(log(viscm/thdifm)*weights)` sites.
- `rotation_scratch_lib.f90`: `public logical function disk_locking_engaged()`
  for the disk-locking predicate in `evolve_angular_momentum`,
  `secular_transport`, `solid_body_omega` (module now imports `star`).

## Deferred (cross-domain)

1. `chi_grid_scale(2)/(9)/(11)` at `src/rotation/equal_spaced_grid.f90:39-41`
   and `src/rotation/seculr/am_transport_grid.f90:230-233`. The mapping to
   `tol_dm_max`/`tol_dl_max`/`tol_dp_core_max` lives in
   `src/setup/map_user_inputs.f90:203-209` and is applied only under
   `lnewvars`, so the named controls cannot be read directly. Proposal: in
   `src/io/controls_lib.f90` (owner of `chi_grid_scale`) add
   `integer, parameter, public :: ichi_dm_max = 2, ichi_dl_max = 9, ichi_dp_core_max = 11`
   next to the array declaration, use them in `map_user_inputs.f90:205/207/209`
   and at the two rotation sites above.
2. `log_luminosity*` arrays hold linear L/Lsun (audit #16). Rename touches
   `src/rotation/rotation_scratch_lib.f90` members and
   `src/mixing/burn_settle_mix.f90`, `src/mixing/semiconvection.f90`,
   `src/mixing/temperature_gradients.f90` besides five rotation signatures --
   one commit across rotation + mixing: `log_luminosity` -> `luminosity_lsun`,
   `log_luminosity_start` -> `luminosity_lsun_start` (pure rename).
3. `no_am_transport_in_core` (`src/io/controls_lib.f90:679`, original LNOJ):
   the audit's `skip_am_transport` rename would change a control name owned by
   io/ and its reader in `state/controls_sync_lib.f90`; the three rotation
   uses (`diffusion_velocity_scales`, `circulation_velocities`,
   `secular_transport`) follow automatically once renamed there.
4. Physical constants (audit #21): `1.6726d-24` (proton mass) at
   `src/rotation/microdiff/gravitational_settling_setup.f90:270` and
   `src/rotation/microdiff/microdiff_coefficients.f90:126`; `6.9598d10` (Rsun)
   in `microdiff_setup.f90:65` / `gravitational_settling_setup.f90:71`;
   `1.6605655d-24` (amu), `5.669d-5` (sigma), `6.7282653d-26` in
   `src/rotation/viscos.f90:25,66,101`; `5.486d-4` (m_e in amu) in
   `gravitational_settling_setup.f90:66` / `microdiff_coefficients.f90:74`.
   `src/state/phys_const_lib.f90` is owned by state/. Proposal: add
   `m_proton_cgs = 1.6726d-24`, `rsun_cgs = 6.9598d10`, `amu_cgs = 1.6605655d-24`,
   `m_electron_amu = 5.486d-4` there **with exactly these values** (they differ
   from the CODATA values used elsewhere, so an existing constant must not be
   substituted), then reference them from the rotation sites.
5. `rotation_scratch_lib.f90` legacy member names (`facd2/3`, `vesd2/3`, `hs3`,
   `pm`, `tm`, `tho`, `qwrst`, `wmst`, `fact1-6`, `fgsfj`, `es1`, `dm`,
   `delami`, `masschg2/3`): several are read by `src/mixing/burn_settle_mix.f90`
   and `src/mixing/*`; rename in a joint rotation+mixing commit. Moving
   `masschg2/3` to `star%` needs `state/star_info_lib.f90`.
6. `check_angular_momentum.f90:278` `write(*,120)` and
   `circulation_velocities.f90` `write(6,9911)` per-shell output: gating on
   `solver_diagnostics()` changes terminal output -- left for an author decision
   (not strictly cross-domain, but excluded by the "prefer not to" rule).

## Reverted (changed numbers)

None. All three pin runs and the aux battery passed on the first attempt.

## Skipped (disagree with reviewer / too risky)

- Audit #3 `am_transport_grid`/`composition_grid` extraction, #4 `cz_extent`
  helper, #5 `secular_transport` first/subsequent-iteration fold, #6/#20
  species-descriptor tables and `interp_to_equal`, #12 `apply_power_law_profile`,
  structural recommendations 1-5: restructuring, not token-identical helper
  extraction; the four CZ-extent copies differ in their gating (the i1 fix is
  present in only one), so a shared helper would change behaviour in three of
  them.
- Audit #13 `solve_r0` helper: the two Newton loops differ (`.le.` vs `.lt.`).
- Audit #7 (`species_fraction` returns a delta -> separate output argument),
  #9 (Bahcall-Loeb unit copies), #10 (drop `use_generic_diffusion_vectors`,
  call the LW steps twice): interface/data-flow changes beyond wave 1; comments
  now state the current behaviour.
- Audit "Possible bugs" #3 `theta_term_n/p` permanently zero: left as is (no
  comment change beyond what was there) -- deleting them would touch the
  coefficient expressions (`+theta_term_n`, a `+0.0d0` that -O3 folds, but the
  expression text would change); an author should decide whether the scheme is
  wanted.
- `check_angular_momentum` `cut_count.gt.0` branch: reachable only if a caller
  ever passes `cut_count > 0`; no caller does today, but it is a documented
  time-step-cut hook -- comment added, code kept.
- Band column layout parameters in `am_advection_diffusion_coeffs` (audit #17):
  the loops index columns with `10+imj`-style arithmetic; naming would either
  restyle every loop or be misleading.
- `zone_moments_of_inertia` `0.4d0`: single use inside an expression; naming it
  adds an indirection without clarifying.
- `compute_quadrupole` inline Thomas solver -> `tridia` (audit #22): pivot-check
  parity not established; comment fixed only.
- BUR-ST extrapolation branches, `burn_settle_mix`, `burn_mix_extrapolated`,
  `alfmlt/phmlt/cmxmlt`: excluded by the assignment.

## Things the audit missed

- `band_nmax = 8000` gives room for 2000 equal-grid points, but `rot_scr%ntot`
  (set in `equal_spaced_grid.f90:51-54` to the number of model zones in the
  transport range) is only bounded by `json = 5000`, and
  `am_advection_diffusion_coeffs` builds `num_equations = 4*num_eq_points-2`
  rows with no check against `band_nmax`. A rotating model with more than
  2000 zones in the transport region would overrun `coeff_matrix`/`rhs`
  silently. Recorded in the `band_nmax` comment; a bounds check (class B) is
  an author decision.
- `diffuse_composition`'s `shell_mass` dummy was dead (only its deleted
  conservation-sum read it); the audit listed the accumulators but not the
  argument.
- `evolve_angular_momentum` `skip_diffusion_flag` is not just "always
  `.false.`" by data: it is computed inside the `.not.instability_transport_active`
  else-branch and tested only on the `instability_transport_active` branch, so
  it could never be true by construction.

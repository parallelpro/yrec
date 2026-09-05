# mixwind wave 1 report (readability batches R1-R3)

Worktree `/Applications/YREC-wt/mixwind`, branch `rs/mixwind`, base `0d5c11a`.

## Commits

| Batch | Hash | Verification |
|-------|------|--------------|
| R1 comments / dead code | `9b61dbc` | gate1 IDENTICAL; full 37-case pins `PINS EXIT 0` (logs/mixwind.pins.R1.log) |
| R2 named constants / argument names | `4298d7b` | gate1 IDENTICAL; full 37-case pins `PINS EXIT 0` (logs/mixwind.pins.R2.log, 37 passed in 1147 s) |
| R3 shared helpers | `6a35c83` | gate1 IDENTICAL after `make clean`; full 37-case pins `PINS EXIT 0` (logs/mixwind.pins.R3.log); aux battery `AUX EXIT 0` (logs/mixwind.aux.log) |

Note on R2: the session was cut by a rate limit while the R2 pin run was in
flight. The run itself completed normally (`37 passed, 272 deselected`,
`PINS EXIT 0` in the log) against the binary built by the R2 gate1 at
00:20; no `make` ran between that build and the commit, and the R3 edits
were only applied to the working tree (never built) while the pins ran.
The commit was made on that evidence rather than re-running the same
selection on the unchanged binary. (The R3 pins and the aux battery
then ran on a `make clean` rebuild that includes every R2 change, so
R2 is additionally covered by those two green runs.)

Note on R3's commit message: the first `git commit` picked up another
domain's message file (shared scratchpad name collision); the message
was fixed with `git commit --amend` before anything else happened. The
tree of 6a35c83 is the one the R3 pins and aux ran on.

All flags as shipped (`-O3 -ffp-contract=off -finit-local-zero -Wall`).
No `make.log` warning was introduced in any batch; the only warning in
mixwind files is the pre-existing `age_gyr` unused dummy in
`wind/massloss.f90` (deferred, see below).

## R1 (9b61dbc): comments and dead code

- wind/wcz.f90: deleted (no `call wcz` anywhere in src/; not in the wind public surface). deps.mk regenerated.
- setup/rotation_stability_setup.f90: dead opacity_interface, grav_const_sq, ht_temp_scale2, dhtscale_dr and the ZM98 "Q variables" (qc1r, qqc1rr, qdr_local, dr_local, qqchitr, v0_local, f2_local, f3_local), dead ff_factor/ht_temp_scale_prev updates; duplicate header, stale commented code.
- util/timestep_limit_heburn.f90: dead locals; energy_gen_terms(6) out-array (never read by the caller) replaced by a scalar and dropped from the signature; util/compute_timestep.f90 (only caller) updated.
- util/timestep_limit_hburn.f90, timestep_limit_hr.f90, timestep_limit_omega.f90, timestep_limit_structure.f90, compute_timestep.f90: dead locals, duplicate `use`, commented debug writes.
- wind/massloss.f90: always-true apply_mass_change flag and dead branch removed, dead locals, header contradiction fixed.
- wind/mdot.f90: wrong "call does not match dummy list" header, dead local.
- wind/kawaler_wind.f90, matt_wind.f90, wind_spindown.f90, wind_spindown_matt.f90, matt_structure_factor.f90: F77 header fragments, commented F77 blocks, COMMON wording, stray `continue`s.
- mixing/*: stale COMMON/SAVE references and pre-modernization argument-list headers rewritten; temperature_gradients SUMMARY 1.1 #7 (dead duplicate `test` branch).
- setup/rezone.f90, setups.f90, map_user_inputs.f90, locate_shell_boundaries.f90, setup_solar_calibration.f90, setup_star_calibration.f90: comment cleanup.
- setup/rescale_model.f90: SUMMARY 1.1 #11 loop variable names; the invalid-species error write (run_log_unit, format 1004; an error path no pinned case reaches) now prints the value actually tested (star%ctrl%new_species_index) instead of the loop counter.

## R2 (4298d7b): named indices / constants, pure renames

- mixing/compute_scale_height.f90: dummies renamed to log_density/log_pressure/log_radius/log_mass/log_temperature.
- mixing/semiconvection.f90: log_luminosity -> luminosity_lsun, log_luminosity_zone -> luminosity_lsun_zone (it is L/Lsun, not a log); only_check_core -> logical parameter.
- mixing/temperature_gradients.f90: header note on unlogged luminosity; max_cubic_iterations = 25, superadiabatic_tol = 1.0d-6.
- mixing/solve_composition.f90: rhs_column_idx -> diagonal_idx (+ stride-8 comment).
- setup/locate_shell_boundaries.f90: luminosity_change_tol, hydrogen_surface_tol as parameters (were data).
- setup/setup_solar_calibration.f90: num_calibration_runs = 48.
- setup/rezone.f90: reaction_rate_species_index(7) parameter array (H1, He4, He3, C12, C13, N14, O16).
- util/compute_timestep.f90: intent(out) dummy hydrogen_dt -> chosen_dt (it stores the chosen step); local hydrogen_dt for the H-burn limit; dt_unlimited = 1.0d20.
- util/timestep_limit_hburn.f90: h_burn_energy_per_gram = 6.00d18. util/timestep_limit_heburn.f90: he_burn_energy_per_gram = 5.85d17, dt_unlimited. util/timestep_limit_omega.f90: dt_unlimited.
- wind/massloss.f90: accretion_efficiency = 1.0d0 parameter.

## R3 (6a35c83): small shared helpers (bodies token-identical to the call sites)

- wind/wind_lib.f90 (new): `log10_radius_from_l_teff(logL, logTeff)` and `matt_centrifugal_factor(omega, fsun, log10R, M_msun)`.
- wind/kawaler_wind.f90, matt_wind.f90, wind_spindown.f90, wind_spindown_matt.f90, matt_structure_factor.f90, massloss.f90: use `log10_radius_from_l_teff` (massloss's `4.0d0` vs the others' `4.d0` is the same double).
- wind/matt_wind.f90 (2 sites), wind/wind_spindown_matt.f90 (fcorr1/fcen1, fcorr2/fcen2): use `matt_centrifugal_factor`; the fsun line with the uninitialised `gl` untouched (PRESERVED BUG comment updated).
- wind/kawaler_wind.f90, wind/matt_wind.f90: contained `wind_domega(dt, omega[, omega_scaled, fcen])` for the torque expression shared by the first guess and the fixed-point loop.
- setup/rezone.f90: eight copy/osplin/copy-back triplets -> sibling internal `regrid_in_place(field)` (same scratch arrays, bounds, osplin call).
- setup/rotation_stability_setup.f90: three Lagrange-weight blocks and three 11-variable interface blocks -> contained `lagrange_weights`, `interp4`, `interpolate_to_interface`; the interior two-pass order (all weights, then all interpolations) kept; dr*/lag_* locals moved into the helper.
- src/deps.mk regenerated.

## Deferred (cross-domain)

1. **core/evolve_step.f90:261 -- drop the unused `age_gyr` argument of massloss.** `wind/massloss.f90` no longer reads `age_gyr` (R1 removed the dead age_seconds). Proposed: remove `star%dage` (2nd actual) from the call at core/evolve_step.f90:261 and then remove the `age_gyr` dummy/decl from massloss.f90:28,45 (I will do the wind side once the core side lands). Removes the only -Wall warning in the domain.
2. **state/star_info_lib.f90:445 / core/main.f90:82-86 -- `hydrogen_dt` is a misnomer.** compute_timestep's 2nd argument is intent(out) and receives the *chosen* timestep (min over all limiters), not the H-burn governor; nothing reads `star%hydrogen_dt` inside compute_timestep or its callees (grep). The comment at core/main.f90:82-86 ("read back by the next model's compute_timestep.f90 call") is wrong. Proposed: rename the member `hydrogen_dt` -> `chosen_dt` in state/star_info_lib.f90:445 and at core/evolve_step.f90:145,223, core/run_yrec.f90:280,348; rewrite the main.f90 note to "DELTSH: chosen timestep written by compute_timestep; read_starting_model reuses the slot transiently as |delta_time|".
3. **io/controls_lib.f90:48-55 (and the other `(50)` run arrays) -- name the run-array dimension.** Proposed: `integer, parameter, public :: max_runs = 50` in controls_lib, dimension the run arrays with it, then setup/setup_star_calibration.f90:43,44,62 and setup/setup_solar_calibration.f90 can use `max_runs` / `max_runs-2` instead of literals 50 / 48.
4. **rotation/rotation_scratch_lib.f90:31 -- named row indices for `reaction_rate_by_zone(15,json)`.** setup/rezone.f90's `reaction_rate_species_index = [1,2,4,5,6,7,9]` should be built from named species-row parameters owned by the rotation domain (e.g. `rr_h1 = 1, rr_he4 = 2, rr_he3 = 4, rr_c12 = 5, rr_c13 = 6, rr_n14 = 7, rr_o16 = 9`). Same value, pure rename.
5. **io/controls_lib.f90:94 (`atime(14)`) and :307 (`chi_grid_scale(12)`)** -- used from setup/map_user_inputs.f90 and setup/rezone.f90 under their old cryptic names; a rename needs the io domain first.
6. **rotation/mid_timestep_model.f90:363 -- `rotmix` rename.** If the mixing domain ever renames `rotmix` (SUMMARY suggests `rotational_mixing_step`), this caller (and mixing/burn_settle_mix.f90:117, mine) must follow; left as is.
7. **state/star_info_lib.f90:274 -- `radiative_zone_bounds(13,2)` / `convective_zone_bounds(12,2)`** dimensions have no named parameters; mixing/find_convection_zones.f90 and setup/locate_shell_boundaries.f90 use the literals 12/13. Proposed: `max_convective_zones = 12`, `max_radiative_zones = 13` parameters in star_info_lib.
8. **state/star_info_lib.f90:226 -- `ln_solar_luminosity` is a misnomer.** It is set to `ln10/solar_luminosity_cgs` (setup/setups.f90:48, core/monte_carlo.f90:100) and used at core/henyey_coefficients.f90:341. Proposed: rename to `ln10_over_lsun` (member + 3 sites). Pure rename.
9. **core/surfbc.f90:185, core/stitched_model.f90:147, core/observables_lib.f90:294-295 -- log10-radius expressions.** Each spells `0.5d0*(logL + log10 Lsun - 4 logTeff - c4pil - csigl)` in a *different operand order* from wind_lib's `log10_radius_from_l_teff` (which follows the wind files' order), so they cannot share the helper without changing rounding. Not proposed for wave 1; a wave-2 numerics-changing batch could unify them.
10. **tools/check_boundaries.py:89 -- stale comment** still says `solid` is "legitimately called from setup/midmod and wind/wcz"; wcz.f90 was deleted in R1. Proposed: drop "and wind/wcz".

## Reverted (changed numbers)

None. Every gate1 was IDENTICAL on the first build of each batch and every full-pin run passed.

## Skipped (disagree with reviewer / too risky)

- Single-precision literals (`0.5*`, `3.0*` in matt_wind/wind_spindown_matt, `1.` in matt_structure_factor, the spot block in temperature_gradients): changing kind would change numbers; kept verbatim (also inside wind_lib's helper).
- `gl` uninitialised in wind_spindown_matt.f90 (fsun line): fixing it changes numbers; documented in the PRESERVED BUG comment.
- `mod` argument order and the `9.4E0` literal flagged in the audit: same-value risk not zero; left.
- semiconvection density write-back and the rezone X/Z gradient blocks: the audit suggested merging them, but the two blocks are not token-identical (different expressions), so no helper.
- solve_composition format-1000 text and rescale_model's unconditional `write(*,*)` blocks: terminal-only, but the assignment says prefer not to touch.
- A util `clamp_dt_factor` helper for the 3-line clamp repeated in timestep_limit_omega/structure/hr: would need a new util module for three lines; not worth it.
- The BS-extrapolation path (burn_settle_mix / burn_mix_extrapolated) and alfmlt/phmlt/cmxmlt: untouched per the assignment.

## Things the audit missed (fixed in R1 unless noted)

- setup/rezone.f90: dead `log10_omega` loop and the `pmax` block (assigned, never read).
- setup/rotation_stability_setup.f90: `opacity_interface` and the ZM98 Q-variables were dead at all three interpolation sites (the audit only listed the header duplication).
- wind/massloss.f90: header contradicted the code, and `apply_mass_change` was always true after the early return.
- setup/rescale_model.f90: `icomp` loop counter printed in the error message instead of the tested value.
- util/timestep_limit_heburn.f90: the `energy_gen_terms` out-array was never read by its only caller.
- core/main.f90:82-86: wrong claim that compute_timestep reads `hydrogen_dt` back (deferred item 2, not fixed -- core file).
- tools/check_boundaries.py:89: stale wcz reference (deferred item 10).

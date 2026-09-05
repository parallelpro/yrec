# Wave 2 readability report: ioctrl (io/ controls, kap/, atm/, eos facade, setup/util read-side)

Worktree `/Applications/YREC-wt/ioctrl`, branch `rs/ioctrl` (base c764c46).

Commits (in order):

- `2759189` — items 1-2: `laol_table_unit` -> `opal92_table_unit`, `ikur2` -> `kurucz_table2_unit`; LAOL / pure-Z LAOL table units via `open(newunit=...)`
- `952aa62` — item 3: `tollaol`, `use_pure_z_table` become `star%ctrl` members (generated state/ files regenerated)
- `93d6d6e` — items 4-5: `rdlaol` returns the metal mixture, `eos_set_debye_huckel_z` setter; `star%use_two_z_tables` decided in core/run_yrec.f90
- `e50494d` — item 6: `alex95_*` -> `alex94_*` identifiers (namelist names untouched)
- `8ae9f96` — item 7: `max_runs = 50` parameter dimensions the per-run arrays
- `1c76b25` — item 8: `ichi_*` / `itime_*` slot names for `chi_grid_scale(12)` / `atime(14)`
- `ccb1867` — item 9: net_lib reaction-slot constants in setup/map_user_inputs.f90; constants-only `use controls_lib` exemption in check_boundaries
- `4478245` — items 10-13: comment/doc fixes

Every commit passed gate1 (build with USE_HDF5=1, domain tests,
check_boundaries, deps.mk check, solar byte-pin) before being
committed. The full 37-case pin selection and the aux battery were
run once on the final commit `4478245` -- see "Verification" at the
end. No pin ever differed, so nothing was reverted.

## Item 1 (2759189): unit-variable renames

- io/controls_lib.f90: `laol_table_unit` -> `opal92_table_unit` (the unit is opened on the OPAL92 table in kap/opal92/read_opal92_tables.f90, never on a LAOL file), `ikur2` -> `kurucz_table2_unit` (kap/kurucz90/read_kurucz_tables.f90). Neither is a namelist name (both are assigned fixed numbers in io/read_controls.f90), so input acceptance is unchanged.
- io/read_controls.f90, kap/opal92/read_opal92_tables.f90, kap/kurucz90/read_kurucz_tables.f90: use sites.
- state/controls_state_def.inc, state/controls_sync_lib.f90 regenerated (member names only).

## Item 2 (2759189): LAOL units via newunit

- kap/laol89/rdlaol.f90: `open(newunit=laol_unit, ...)` / `laol2_unit` locals replace the fixed `iolaol = 61`; kap/laol89/rdzlaol.f90: `pure_z_unit` replaces `iopurez = 62`. Both units were opened and closed within the one routine, so the fixed numbers were never observed elsewhere. `iolaol`/`iopurez` removed from kap/opacity_table_lib.f90 and the assignments removed from io/read_controls.f90 (a comment there records the move).

## Item 3 (952aa62): tollaol / use_pure_z_table into star%ctrl

- io/controls_lib.f90: `tollaol = 10.0d0` and `use_pure_z_table = .false.` declared (former common/nwlaol/ block from kap/opacity_table_lib.f90). `tollaol` is itself the `/physics/` namelist name (canonical `laol_z_tolerance` in the new-style inlist) and `use_pure_z_table` is the new-style canonical name of `lpurez`, so both buffer names are unchanged and namelist acceptance is unchanged.
- Readers: kap/laol89/gtlaol.f90, gtlaol2.f90, gtpurz.f90 (`star%ctrl%tollaol`), kap/setupopac.f90 and kap/kap_lib.f90 (`star%ctrl%use_pure_z_table`). Tests kap/atm/eos/net set `star%ctrl%use_pure_z_table = .false.` (test_net's `use opacity_table_lib, only: use_pure_z_table` removed).
- `llaol` NOT moved: it is not a namelist input and nothing assigns it anywhere (always `.false.`); the only read is core/read_starting_model.f90:1011 `if (.not.llaol)`. Left in kap/opacity_table_lib.f90 with a comment saying so -- see Deferred.
- io/read_yrec7 / read_model2's `use_pure_z_table` intent(out) dummies receive core/read_starting_model.f90's local `use_pure_z_table0`, so no model read writes the control.
- state/ files regenerated (366 members).

## Item 4 (93d6d6e): rdlaol no longer writes debye_huckel_z

- kap/laol89/rdlaol.f90: `subroutine rdlaol(laol_work_array, laol_debye_huckel_z, laol_table_path, laol_table2_path, ierr)`; the 18-element mixture read from the table goes to the new `intent(out)` argument; `use yale_eos_lib` dropped.
- eos/eos_lib.f90: new public `eos_set_debye_huckel_z(metal_mixture)` (plain array copy into yale_eos_lib's `debye_huckel_z`).
- kap/setupopac.f90: local `laol_debye_huckel_z(18)`, `call rdlaol(...)`, `if (ierr /= 0) return`, `call eos_set_debye_huckel_z(laol_debye_huckel_z)`, then `call sulaol` as before. Nothing reads `debye_huckel_z` between the table READ and the old in-place store (its only reader is eos/yale/fully_ionized_eos.f90, called much later), and the ierr early return leaves it untouched exactly as before.
- tools/check_boundaries.py: `eos_set_debye_huckel_z` added to the eos allowlist (comment explains).

## Item 5 (93d6d6e): use_two_z_tables decided in core

- core/run_yrec.f90: after `call read_controls(ierr)` and before `call star_setup(ierr)`: `star%use_two_z_tables = star%ctrl%use_z_ramp .or. star%job%use_diffusion_z`. Equivalent placement: `setups` (the only caller of `kap_init` -> `setupopac`) is called only from core/star_setup.f90, which run_yrec calls right after read_controls; `use_z_ramp` / `use_diffusion_z` are written only by read_controls / inlist_new_read; the four domain tests set `star%use_two_z_tables = .false.` themselves before calling setups.
- kap/setupopac.f90: the assignment removed (comment points at run_yrec).

## Item 6 (e50494d): alex95 -> alex94

- kap/opacity_table_lib.f90 (`n_alex94_*`, `alex94_*` arrays), kap/surfopac.f90, kap/alex94/*.f90, kap/kap_lib.f90 comment, io/controls_lib.f90 (`alex95_table_unit` -> `alex94_table_unit`) and io/read_controls.f90 use site. No quoted strings, formats or namelist names changed: `lalex95`, `use_alex95_tables`, `alex95_table_z_list` keep "95" (comment in controls_lib says so). `state/star_info_lib.f90`'s `alex95_table_paths` member is not mine -- see Deferred.

## Item 7 (8ae9f96): max_runs

- io/controls_lib.f90: `integer, parameter, public :: max_runs = 50`; all 27 per-run buffer arrays (the eight structure-limit stop arrays, target_end_age / timestep_override / central_*_stop / *_active, rescale_params(4,max_runs), rescale_kind, num_models, first_call_flag, initial_x/z_array, mixing_length_array, senv0_array, has_senv0_array, rsclzc/rsclzm1/rsclzm2) dimensioned `(max_runs)`.
- tools/gen_controls_state.py: collects `integer, parameter :: name = <digits>` lines and substitutes the literal into emitted dims, because state/star_info_lib.f90 includes the generated type bodies without `use controls_lib`. The regenerated state/ files are byte-identical (so state/ is untouched by this commit).
- setup/setup_star_calibration.f90 (`50`/`49` -> `max_runs`/`max_runs-1`) and setup/setup_solar_calibration.f90 (`num_calibration_runs = max_runs-2`) `use controls_lib, only: max_runs`.
- tools/check_boundaries.py: a `use controls_lib, only:` list drawn entirely from named constants (`CONTROLS_LIB_CONSTANT_RE`, initially `max_runs`) is not a buffer import. This is the one allowlist change of the wave beyond the eos entry; every other `use controls_lib` outside CONTROLS_BUFFER_IMPORTERS still fails. An alternative that needs no exemption -- declaring `max_runs` in state/star_info_lib.f90 -- is the state owner's call (Deferred).

## Item 8 (1c76b25): chi_grid_scale / atime slot names

- io/controls_lib.f90: `ichi_dm_min=1, ichi_dm_max=2, ichi_flag_dx=3, ichi_flag_dz=4, ichi_dx_max=5, ichi_dz_max=6, ichi_czbase_fine_width=7, ichi_dp_env_max=8, ichi_dl_max=9, ichi_dp_czbase_max=10, ichi_dp_core_max=11, ichi_flag_dw=12` and `itime_core_min=1, itime_dx_core_tot=2, itime_dx_core_frac=3, itime_dy_core_tot=4, itime_dy_core_frac=5, itime_dx_total=6, itime_dx_shell=7, itime_dt=8, itime_dp=9, itime_dr=10, itime_dl=11, itime_dy_shell=12, itime_max_dt_frac=13, itime_dy_total=14`. The brief's `ichi_dm_max = 2, ichi_dl_max = 9, ichi_dp_core_max = 11` verified against setup/map_user_inputs.f90:198-255; every one of the 26 slots is explicit there, so all are named.
- Readers subscript by name (numeric -> same-valued constant; one `use controls_lib, only: ...` line after `use star_info_lib`): core/read_starting_model.f90, core/rebuild_envelope.f90, rotation/equal_spaced_grid.f90, rotation/seculr/am_transport_grid.f90, rotation/seculr/composition_grid.f90, setup/rezone.f90, util/timestep_limit_{structure,omega,hburn,heburn,hr}.f90, util/compute_timestep.f90, setup/locate_shell_boundaries.f90, setup/map_user_inputs.f90. The shared files got subscripts plus the use line only; the wholly-owned util/timestep_limit_* files also had their header comments' `atime(N)` renamed. The `atime(N)` mentions left in the comments of util/compute_timestep.f90:12,153, setup/locate_shell_boundaries.f90:24,41 and `chi_grid_scale(7,8,10)` at setup/rezone.f90:220 are untouched (shared files) -- see Deferred.

## Item 9 (ccb1867): net_lib reaction constants in map_user_inputs

- setup/map_user_inputs.f90: `use net_lib, only: r_pp, r_he3he3, r_he3he4, r_pc12, r_pc13, r_pn14, r_po16, s_pep, s_be7e, s_be7p, s_hep, iq_be7p`; `cross_section_scale(14..17)` -> `s_pep/s_be7e/s_be7p/s_hep`, `qs0e_scale(8)`/`qqs0ee_scale(8)` -> `iq_be7p`, and (beyond the brief, same hunk) `cross_section_scale(1..7)`, `qs0e_scale(1..7)`, `qqs0ee_scale(4,5,7)` -> `r_*` (net_lib's `rates`/core `engeb` index these arrays with the reaction loop index, so slots 1-7 are the reaction slots). io/read_controls.f90 has no such subscripts. net is all-public and setup is app-side, so no allowlist change was needed.
- tools/check_boundaries.py: the constants exemption extended to `ichi_*`/`itime_*` and made to follow `&` continuation lines.

## Items 10-13 (4478245; the check_boundaries hunks landed in ccb1867)

- 10: tools/check_boundaries.py rotation comment: "solid ... legitimately called from setup/midmod and wind/wcz" -> "rotation_shape_factors (the former solid) ... called from core and setup/rezone" (its callers today: core/read_starting_model, core/evolve_step, core/henyey_iterate, setup/rezone; no wind/ file calls it).
- 11: eos/eos_lib.f90 public-line comment marks `eos_get_gamma1` as exported for eos/test/test_eos.f90 only; the check_boundaries eos allowlist comment says the same (the entry stays so the test build passes). Not deleted.
- 12: io/read_controls.f90: the eight fixed MHD units (`unit_zams_a/b/c` = 40/41/42, `unit_centre1-5` = 43-47) stay. They are provably NOT equal to the declared defaults: the star%ctrl members have no declaration-time default (generator seeds 0) and eos/mhd/mhdst.f90 treats a unit `.le. 0` as "table absent" (mhdst.f90:63,100,...,163), so dropping the assignments would change use_mhd_eos runs. A comment at the assignment block records this.
- 13: GUIDELINES.md "Regression-testing gotchas": one bullet on the aux test's `names(<digits>) = '...'` regex over io/profile_output.f90.

## Deferred (not mine or outside the brief)

1. state/star_info_lib.f90: `star%job%alex95_table_paths` -> `alex94_table_paths` (use sites: io/read_controls.f90:1401,1428 `opecalex = star%job%alex95_table_paths` / store, plus any kap/ readers) to finish item 6. Safe to rename: the namelist names are `opecalex` (legacy) and `alex95_table_z_list` (new-style), both separate strings.
2. state/star_info_lib.f90 ~248 comment "use_two_z_tables is derived by setupopac" is stale after item 5: -> "decided by core/run_yrec.f90 right after the controls read (use_z_ramp .or. use_diffusion_z)".
3. eos/yale_eos_lib.f90:49-52: comment says `debye_huckel_z` is "recomputed by mixing/compute_scale_height.f90"; it is not -- it is set only by `eos_set_debye_huckel_z` (kap/setupopac.f90 after `rdlaol`). One-line fix for the eos owner.
4. kap/opacity_table_lib.f90 `llaol`: never assigned anywhere (always `.false.`), read only at core/read_starting_model.f90:1011 `if (.not.llaol) then`. Author decision: delete the variable and the guard (the guarded branch is the only branch ever taken), or wire it to an input. Class-B if it were ever meant to be true; today deleting is a no-op.
5. Item 7 alternative: a `max_runs` parameter in state/star_info_lib.f90 (the star%job arrays' home) would let the two setup files use it without any `use controls_lib`; then CONTROLS_LIB_CONSTANT_RE could shrink back. State owner's call.
6. Item 8 leftovers in shared files (subscripts-only rule): comments at util/compute_timestep.f90:12 and :153 (`atime(13)` -> `atime(itime_max_dt_frac)`), setup/locate_shell_boundaries.f90:24 and :41 (`atime(1)` -> `atime(itime_core_min)`), setup/rezone.f90:220 (`chi_grid_scale(7,8,10)` -> `(ichi_czbase_fine_width, ichi_dp_env_max, ichi_dp_czbase_max)`), rotation/seculr/composition_grid.f90 / am_transport_grid.f90 / equal_spaced_grid.f90 have no numeric mentions left.
7. kapatm wave-1 deferred 9 (names for the `laol_work_array` index slots) untouched -- not in this brief.
8. io/read_controls.f90 ~245-251 comment block still describes `chi_grid_scale` as "const_lib's"; harmless but stale (const_lib is gone). Not touched because that hunk was not in any item.

## Reverted (changed numbers)

None. Every gate1 was IDENTICAL; the full pin selection and aux battery passed on the final commit.

## Skipped

- `llaol` move (item 3): not a namelist input and never written -- nothing to move; see Deferred 4.
- Dropping the MHD unit assignments (item 12): would change behaviour (units would seed as 0 = "absent"); documented instead.
- `alex95_table_paths` (item 6) and the star_info_lib comment: file not mine.

## What the audit / wave-1 reports got wrong

- kapatm wave-1 (Deferred 3) called `tollaol`, `llaol` and `use_pure_z_table` "three namelist inputs living in kap/opacity_table_lib". Only `tollaol` and `lpurez -> use_pure_z_table` are namelist inputs; `llaol` is neither read from any namelist nor assigned anywhere.
- The `iolaol`/`iopurez` units were described as "shared" fixed units; each is opened and closed within one routine (rdlaol / rdzlaol), so `newunit` was a pure local change.
- The mixwind wave-1 note behind item 7 assumed the setup files could see a controls_lib constant directly; they cannot without the check_boundaries exemption (or a state-side parameter), which the report above records.

## Verification

- gate1 before every commit: IDENTICAL (8 commits).
- Full 37-case pin selection on 4478245: `37 passed, 272 deselected in 1583.34s`, `PINS EXIT 0` (machine shared with two other worktrees' pin runs, hence the wall time).
- Aux battery on 4478245: `22 passed in 570.85s`, `AUX EXIT 0`.
- `make clean` was done at the start of the wave and again after the item-3 type change; the pins/aux run above used the gate1 build of 4478245.

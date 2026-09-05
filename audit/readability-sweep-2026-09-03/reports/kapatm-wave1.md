# Wave 1 readability report: kap + atm domain

Worktree `/Applications/YREC-wt/kapatm`, branch `rs/kapatm` (base 0d5c11a).

Commits:

- R1 `fed1d0c` — stale comments, dead code, unused locals/params, ll95tbl EOF error path
- R2 `53e7724` — one named-dimension block per table family, named sentinels and thresholds
- R3 `2959555` — shared stencil-search and Z-blend helpers

Every commit passed gate1 (build, tests, check_boundaries, deps.mk
check, solar byte-pin) and the full 37-case pin selection
(`PINS EXIT 0`) before being committed; the aux battery
(`AUX EXIT 0`) was run on the final commit. No pin ever differed, so
nothing was reverted.

## R1 (fed1d0c): comments and dead code

Comments/headers:

- kap/opal95/getopal95.f90: header no longer lists the F77 OP952D/... entry names; `if (.true.) then` wrapper removed; two `continue` before `return` removed; unused `use star_info_lib` (x2) and `num_xz` removed.
- kap/opal92/opal92_interp3d.f90: header no longer names YLLO3D/SETLLO; "LN" corrected to log10.
- kap/kurucz90/kurucz.f90, kurucz2.f90: M1/M2/M3/GETD/'TEM' vocabulary replaced by the current names; `local_logt/local_logrho/t_row_index` aliases replaced by the intent(in) dummies they copied.
- atm/tables/alsurfp.f90: header argument list rewritten (log_teff, log_g, print_to_files, lookup_failed, ierr, results in atm_table%...); COMMON /ATMPRT/ / KTTAU / LPRT text removed; orphan "End if inputs are valid" removed; write-only `bad_point = .true.` at the extended-range exit removed.
- kap/alex06/alex06tab.f90: common/alot06/ text removed; X weights go to their own `weight_x(4)` instead of overwriting `weight_z` (not read after that point).
- kap/alex94/alex94_surface_table.f90: "common/alot/" text removed.
- atm/atm_lib.f90: header no longer describes atm_get/envint (now in core/envint_lib.f90).
- atm/atmstruct_lib.f90, atm/envstruct_lib.f90: "computed by atm/atm_lib.f90" and a reader list naming non-existent files replaced.
- atm/tables/surfp.f90, kcsurfp.f90: "4-POINT LAGRANGIAN" corrected (it is kspline/ksplint); write-only atm_table members no longer assigned.
- atm/tables/alfilein.f90: fixed-form essay trimmed to a one-line note on the no-op EXTERNAL; ierr declaration moved; dead `continue` removed from sort_shell.
- atm/tables/altabinit.f90: ierr declaration moved after `implicit none`.
- kap/kap_lib.f90: "OPTIONAL ierr"/"single failure funnel"/"<.15" remarks removed; write-only `jerr` in kap_update_surface_tables removed.
- kap/test/test_kap.f90, atm/test/test_atm.f90: "optional ierr" comments corrected.
- kap/opal92/opal92_surface_table.f90: "CONVERSION FROM THE DERIVATIVE" comment above plain assignments removed; write-only `opacity_final` removed.
- kap/opal95/ll95tbl.f90: /newopac/ block comment removed; SUMMARY 1.1 #5: a short table set (EOF at the next header) closes the unit and returns ierr=1 instead of silently succeeding; EOF inside a table body reports and returns ierr=1. Dead in every pinned run (every OPAL95 file under input/opacity has exactly 126 tables).
- kap/conductive/condopacp.f90: stale comment block replaced.
- Banners still carrying ykoeff/setkrz/setllo/ylloc/ll4th/op95* names updated across kap/.
- Runs of blank lines trimmed in the alex94/alex06/atm table files.

Dead code (how shown dead is in the commit message):

- kap/laol89/gtlaol.f90, gtlaol2.f90: unreachable "OUTSIDE OPACITY TABLE(#2)" else-branch and its `if (num_valid_x .ge. 2)` wrapper removed; `local_x/local_logt/local_logrho` aliases replaced.
- kap/laol89/gtpurz.f90: `local_logt/local_logrho` aliases replaced.
- kap/kurucz90/build_kurucz_splines.f90: two unreachable `if (jd.le.1) cycle`; `index1/ids/idf` aliases replaced.
- kap/opal92/opal92_table_prep.f90: `index1/ids/idf` aliases replaced.
- kap/opal92/opal92_interp2d.f90, opal92_interp2d_z2.f90: write-only `lmore`.
- kap/opal92/read_opal92_tables.f90: write-only `local_grid_y`, `const_unused`.
- kap/kurucz90/read_kurucz_tables.f90: write-only `x_table_count`.
- atm/atm_table_lib.f90: write-only members teff_interp_start_index, gravity_interp_indices(4), imax1_placeholder/imax2_placeholder/ljvs_placeholder.
- kap/opacity_table_lib.f90: never-referenced kurucz_ix_x / kurucz2_ix_x.
- 30 files: unreferenced per-file `integer, parameter` copies of table dimensions.

## R2 (53e7724): named indices and constants

- kap/opacity_table_lib.f90: new parameters `n_opal95_tables = 126`, `opal95_missing_opacity = 9.999d0`, `opal95_missing_test = 9.9d0`, `opal92_missing_opacity = -9.999d0`, `n_laol_x/rho/t = 12/104/52`, `laol_opacity_cap = 1.0d35`, `alex_composition_tol = 1.0d-8`, `opal92_x_match_tol = 1.0d-5`, `opal95_composition_tol = 1.0d-4`; members `abund_index/temp_index/dens_index(_z2)` renamed `opal92_index_x/t/rho(_z2)`; bare LAOL module variables `olaol/oxa/ot/orho/numofxyz/numrho/numt` and the `*2` twins moved into the type as `laol_opacity/laol_grid_x/laol_grid_t/laol_grid_rho/laol_num_x/rho/t` and `laol2_*`.
- kap/conductive_table_lib.f90: `cond_table_loaded_marker = 12345`.
- atm/atm_table_lib.f90: `atm_table_ng = 11`; arrays declared with `atm_table_nt/ntc x atm_table_ng`.
- kap/opal95/{getopal95,ll95tbl,opal95_fixed_z_table,opal95_interp2d,opal95_interp3d,opal95_interp4d,opal95_surface_table}.f90, kap/opal92/*.f90 (7), kap/alex94/*.f90 (4), kap/alex06/*.f90 (3), kap/kurucz90/*.f90 (4): per-file `num_t/num_d/num_x/num_z` copies deleted; uses name the family parameter.
- kap/conductive/condopacp.f90: `n_temp_grid/n_rho_grid/n_z_grid` -> `cond_n_temp/rho/z`.
- kap/laol89/{gtlaol,gtlaol2,gtpurz,sulaol,zsulaol,rdlaol,rdzlaol}.f90: work arrays and size checks use `n_laol_*` (rdlaol's `11` kept as `n_laol_x-1`); LAOL member renames.
- kap/opal95/opal95_surface_table.f90: `x_stencil_6_7_8_10 = 11`; consumer test `x_table_index.lt.10` -> `.ne.x_stencil_6_7_8_10` (values are 1..7, 6 or 11).
- kap/kurucz90/kurucz.f90: `kurucz_logrho_max = -3.0d0`, `kurucz_logt_max = 4.1d0`.
- atm/tables/surfp.f90, kcsurfp.f90: `table_logteff_min = 3.5d0`, `table_logg_min = -0.5d0`, `logteff_three_point = 4.5d0`, `logteff_two_point = 4.55d0`; local nt/ng/ntc/ngc copies dropped in favour of atm_table_*.
- atm/atm_lib.f90, atm/tables/alfilein.f90: use atm_table_* parameters.

## R3 (2959555): small shared helpers

- kap/kap_lib.f90: module procedure `blend_in_z(value_z1, value_z2, metal_fraction, table_z1, table_z2)` replaces the four copies (x3 quantities) of the slope-then-blend Z interpolation in kap_eval (Kurucz, pure-Z with `table_z2 = 1.0d0`, OPAL, LAOL); local `slope` removed.
- kap/opacity_table_lib.f90: contains section with `stencil4_locate` (ALEX94/ALEX06 warm-start T and R search), `stencil4_locate_opal95` (OPAL95 variant) and `alex_x_stencil_start` (X ladder over alex95_grid_x(3:5)); all integer/comparison only.
- kap/alex94/alex94_interp3d.f90, kap/alex06/getalex06.f90: T and R searches -> `stencil4_locate`; X ladder -> `alex_x_stencil_start` (also in kap/alex94/alex94_surface_table.f90).
- kap/opal95/getopal95.f90: T and rho searches -> `stencil4_locate_opal95`; the X-stencil block duplicated between the 3D (k = 1) and 4D (k loop) paths -> contained `opal95_x_stencil(k, z_here, X)`; unused `x_shift_base` removed.
- kap/opal95/opal95_surface_table.f90: identical `X <= 0.8` and `fixed Z <= 0.04` branches merged with `.or.`.

## Deferred (cross-domain)

1. io/controls_lib.f90:285 — `laol_table_unit` is used only for the OPAL92 file (kap/opal92/read_opal92_tables.f90:35-47). Proposed: rename the member to `opal92_table_unit` in controls_lib and read_controls, then in read_opal92_tables.f90.
2. io/controls_lib.f90:646 — `ikur2` (second Kurucz table unit, kap/kurucz90/read_kurucz_tables.f90:74-80). Proposed: rename to `kurucz_table2_unit` alongside the existing `kurucz_table_unit` (io/controls_lib.f90:286).
3. io/read_controls.f90:1442,1444 — `iolaol = 61`, `iopurez = 62` are hard-coded unit numbers stored in bare opacity_table_lib variables and used in kap/laol89/rdlaol.f90 and rdzlaol.f90:24-43. Proposed: `open(newunit=...)` in rdlaol/rdzlaol and delete the two assignments and the bare variables; or, minimally, move them into `star%ctrl` next to the other table units.
4. kap/opacity_table_lib.f90: `tollaol, iolaol, iopurez, llaol, use_pure_z_table` stay bare module variables because they are read/written from io/read_controls.f90:914, io/read_yrec7.f90:24-109, io/read_model2.f90:24-82 and core/. Proposed: move `tollaol/llaol/use_pure_z_table` into `star%ctrl` (they are namelist inputs) and change those io/ sites.
5. kap/laol89/rdlaol.f90:59 — `read(iolaol,130) zhit, debye_huckel_z` writes a yale_eos_lib variable from the opacity reader. Proposed: read into a local and have the caller (kap/setupopac.f90) or the eos setup assign `debye_huckel_z`; needs an eos_lib entry.
6. kap/setupopac.f90:35 — `star%use_two_z_tables = star%ctrl%use_z_ramp .or. star%job%use_diffusion_z` is a star_info_lib member being decided inside the opacity setup. Proposed: compute it in core/ (where `use_z_ramp`/`use_diffusion_z` are resolved) and let setupopac only read it.
7. kap/opacity_table_lib.f90 `n_alex95_*` parameters and `alex95_*` type members describe the Alexander 1994 tables (files kap/alex94/*). Rename to `n_alex94_*`/`alex94_*` touches io/read_controls.f90:605-870 (`lalex95`, `alex95_table_unit`) and io/controls_lib.f90; deferred.
8. numerics_lib: the reviewer suggested a domain-neutral stencil search helper. `stencil4_locate`/`stencil4_locate_opal95` were kept in opacity_table_lib (the only users are kap/). If another domain grows the same pattern, move them to numerics_lib unchanged.
9. kap/laol89/rdlaol.f90:50-51 — `laol_work_array(6), (9), (8), (11), (1), ...` element meanings are defined by core/read_starting_model.f90:135-1148 (`species_mix_weights(12)`, "V(7)=H, V(12)=HE"). Proposed: an `integer, parameter` index block (`ix_h = 7, ix_he = 12, ...`) in the state/ or core/ module that owns the 12-element mixture layout, then used by rdlaol.

## Reverted (changed numbers)

None. Every batch pin run was byte-identical on the first attempt.

## Skipped (disagree with reviewer / too risky)

- SUMMARY 1.1 #1 alternate-return removal: out of scope by assignment.
- kap/kurucz90/read_kurucz_tables.f90 `cycle` -> error return on a short table: changes behaviour on a path no pin reaches but that a user table could; not a readability edit.
- kap/alex06/readalex06.f90 header/row reader helper: the two loops (X blocks, then the last-X Z blocks) differ in their index arithmetic (`jj`) and column checks; not token-identical.
- kap/opal92/opal92_table_prep.f90 two strides: not provably equal; left as is.
- kap/kurucz90 and kap/opal92 `_z2` clone merges: M-sized restructurings, not R1-R3 helpers.
- kap/opal95/ll95tbl.f90 row-band reads (formats 30-35 over 57/1/2/4/5/1-row bands): a loop needs a per-band width table; the Fortran rule that a shorter I/O list stops at the first unused edit descriptor makes it byte-safe, but the gain is small.
- BS-extrapolation path and alfmlt/phmlt/cmxmlt: untouched by instruction.

## Things the audit missed

- kap/laol89/sulaol.f90:58 — the second-table log10 loop runs over `opacity_table%laol_num_t` (first table's count) while indexing `laol2_grid_t`; the loop just below uses `laol2_num_t`. Harmless when both tables share a T grid (all shipped LAOL pairs), otherwise it under- or over-converts. Behaviour preserved; flagged for review.
- atm/tables/altabinit.f90:106-108 — format 910 `(a,i3)` is used with two character items and an integer; the second string and the integer are lost to format reversion ('Less than 4 rows: nTeff = ' prints, the count does not). Bad-table path only, run-log; preserved.
- kap/alex06/readalex06.f90:65 — the row-R check validates only the first 16 of 19 columns (`do kk=1,16`). Left as is.
- atm/atmstruct_lib.f90 header listed writer files that do not exist (fixed in R1).
- kap/laol89/rdlaol.f90 size check compares the X count against 11 while the first array extent is 12 (`n_laol_x`); preserved as `n_laol_x-1`.

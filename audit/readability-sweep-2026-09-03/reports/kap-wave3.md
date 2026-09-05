# Readability wave 3 (R5) -- kap/atm

Worktree `/Applications/YREC-wt/kap`, branch `rs/kap`, based on yrec-modern
`6449462`. Seven commits, all class-A (no output byte changed).

| commit    | item | summary |
|-----------|------|---------|
| `6b29141` | 1    | merge the two Kurucz90 table sets into `type kurucz_table_set`; `kurucz2.f90` deleted |
| `fd65557` | 1    | merge the two LAOL89 table sets into `type laol_table_set`; `gtlaol2.f90` deleted |
| `53177df` | 1    | merge the two OPAL92 table sets into `type opal92_table_set`; `opal92_interp2d_z2.f90`, `opal92_interp3d_z2.f90` deleted |
| `0f860da` | 1    | merge the Kurucz and Kurucz/Castelli surface-pressure tables into `type surface_p_table`; `kcsurfp.f90` deleted |
| `64d5967` | 2    | name the remaining table-dimension literals in `kap/**` |
| `7b03580` | 3    | use `numerics_lib`'s `lagrange4` for the contiguous 4-point stencil sums |
| `3c36ef5` | 4    | pass the OPAL95 stencil to `opal95_interp2d/3d/4d` as explicit arguments |

## Verification

- gate1 (`tools/gate1.sh`) IDENTICAL before every commit; `make clean`
  before each build that changed a derived type, module array or
  procedure signature.
- Full 37-case pin selection (`fullpins.sh ... pins`) after the last
  commit (`3c36ef5`): 37 passed, 272 deselected, `PINS EXIT 0`.
- Aux battery (`fullpins.sh ... aux`): 22 passed, `AUX EXIT 0`.
- The pins exercise only OPAL95 + ALEX06 (+ pure-Z tables, atm choices
  0/1/3/4). Every family this wave touched that no pin covers was
  cross-checked by running a variant of `Test_solar_noGS_norot` with
  the `6449462` binary and the new one and `diff -r`-ing the output
  directories (`.log`, `.mod`, `history.data` and everything else,
  ignoring only the version/wall-clock banner lines):
  - Kurucz90 two-Z (`LKUR90`, `kapp00`/`kapm03`, `LZRAMP`) after
    `6b29141` and again after `fd65557`: byte-identical.
  - LAOL89 two-Z (`LLAOL89`, `Z0188AG89`/`Z0190AG91` `.DBGLAOL`,
    `LZRAMP`) after `fd65557`: byte-identical.
  - OPAL92 two-Z (`LOPAL92`, `Z0188.OPALN`/`Z0200.OPALN`, `LZRAMP`)
    after `53177df`: byte-identical.
  - Castelli/Kurucz atmosphere (`KTTAU=5`, `atmk2004p00.tab`) after
    `0f860da`: byte-identical.
  - ALEX94 (`LALEX95`, `OPACALEXANDER.X00..X08`) after `7b03580`:
    byte-identical.
  Each variant's output was also confirmed to differ from the
  OPAL95/ALEX06 reference so the alternative path really ran.

## Item 1 -- table-set derived types

### Kurucz90 (`6b29141`)

`type kurucz_table_set` in `opacity_table_lib.f90`: `grid_logt`,
`log10_opacity`, `log10_rho`, `num_temps`, `ix_t = 1`, `ix_rho = 1`,
`spline_coeffs`, `density_start_index`, `density_count` (the former
`kurucz_*` members, prefix dropped, same shapes) plus
`logical :: check_range = .true.`. `opacity_table%kurucz(2)` replaces
the `kurucz_*`/`kurucz2_*` pairs; module variable, so `save`
semantics unchanged.

`diff -w kurucz.f90 kurucz2.f90`: exactly one non-rename difference,
the early-out `if (.not.(log10_density.le.kurucz_logrho_max .and.
log10_temperature.le.kurucz_logt_max)) return 1` that only
`kurucz.f90` had. It is kept in place under a plain
`if (tbl%check_range) then`; `read_kurucz_tables` clears
`opacity_table%kurucz(2)%check_range` before reading. All other
statements token-identical after the rename (one `grdnt` line was
merely continued differently). `build_kurucz_splines(tbl, ierr)` now
builds one set (its halves were identical except the stdout strings
"kurucz table"/"second kurucz table", now the generic first form);
`read_kurucz_tables` keeps its signature and calls an internal
`read_kurucz_table(unit, path, tbl, ierr)` twice with the FORMAT/OPEN/
READ/CLOSE statements moved verbatim, working scalars host-associated
so the carry-over between tables is as before, order read 1, read 2,
splines 1, splines 2 preserved. `kap_eval` passes `kurucz(1)`/`kurucz(2)`.

### LAOL89 (`fd65557`)

`type laol_table_set`: `opacity`, `grid_x`, `grid_t`, `grid_rho`,
`num_x/num_rho/num_t`, `slaol_opacity`, `slaol_log_rho`,
`slaol_d2opacity`, `slaol_num_points`, plus `integer :: table_number = 1`;
`opacity_table%laol(2)`. The pure-Z `zlaol_*`/`zslaol_*` members were
not touched.

`diff -w gtlaol.f90 gtlaol2.f90`: member renames, comments, and the two
run-log FORMATs 120/121 whose text in `gtlaol2` reads
`' OUTSIDE OPACITY TABLE #2, IN DENSITY.  '` /
`' OUTSIDE OPACITY TABLE #2, IN TEMPERATURE.  '`. These reach the run
log, so both wordings are kept verbatim (labels 220/221 for the #2
form) under `if (tbl%table_number == 2)`. `rdlaol` keeps its
signature; its two read sequences became an internal
`read_laol_table(path, tbl, work_array, debye_huckel_z, ierr)`. The
only textual difference between the halves was the header-mixture read
list -- `work_array2(6),(9),(8),(11),(1),(3),(2),(5),(10),(4)` versus
`laol_work_array(ix_c, ix_n, ix_o, ix_ne, ix_na, ix_mg, ix_al, ix_si,
ix_ar, ix_fe)` -- and with `ix_na=1, ix_al=2, ix_mg=3, ix_fe=4, ix_si=5,
ix_c=6, ix_o=8, ix_n=9, ix_ar=10, ix_ne=11` these are the same
positions, so one read list serves both (table 2 still reads into the
discarded `work_array2`/`zdh2`). `' SECOND OPACITY ARRAY TOO LARGE.'`
kept verbatim for table 2. `sulaol`'s spline loop became the internal
`spline_laol_table(tbl)`; the log10 conversion of the T grid stays in
`sulaol` because the second half loops to the FIRST table's `num_t`
(see R6 note below).

### OPAL92 (`53177df`)

`type opal92_table_set`: `grid_logt`, `grid_x`, `grid_logr`,
`log10_opacity`, `num_x`, `num_temps`, `spline_coeffs`,
`density_start_index`, `density_count`, `surface_x`, `surface_z`,
`surface_spline_coeffs`, `surface_x_index`, `index_x/index_t/index_rho
= 1`; `opacity_table%opal92(2)`. The set-2 cursors were never
DATA-initialised in the original (`kipmll2`) and so started at 0; they
now share set 1's `= 1` default. Every read of a cursor goes through
`findex`, which resets a value outside `1..n` to 1 before use, so 0 and
1 were indistinguishable (header comment in the type says so).

`diff -w` of each interpolator pair: member suffix, comment style, the
routine name in stdout-only messages (`'opal92_interp2d_z2: CHECK
NDSS'`, `'extrapolation fails'`, `'opal92_interp3d_z2: error in X2
grid'`, `'T out of table'`) and one `grdnt` assignment continued over
two lines. Both `opal92_interp3d` versions contain the X-grid `findex`
block -- see "Wrong in the brief/audit" below. `read_opal92_tables` /
`opal92_table_prep` / `opal92_surface_table`: the second-Z halves
became internal helpers (`read_opal92_table`, `prep_table`,
`build_surface_table`) called once per set in the original order; both
OPEN statements stay in `read_opal92_tables` (fixed unit vs NEWUNIT),
READ/FORMAT/CLOSE moved verbatim; `rho_search_index` stays
host-associated so it carries over between the two surface-table
builds as before; the surface Z is passed in (`star%ctrl%opal_table_z1`
/ `opal_table_z2`). The "second"/"second-Z" wording in stdout-only
diagnostics is dropped. `kap_eval` passes `opal92(1)`/`opal92(2)` at
the three call sites.

### Surface pressure (`0f860da`)

`type surface_p_table` in `atm_table_lib.f90`: `num_teff`, `teff`,
`logg`, `log10_pressure`, `gmin_index`, `gmax_index`; two named
instances `atm_table%kurucz_surface_p` / `castelli_surface_p` (named
rather than an array because the two are selected by `KTTAU`, not by
an index). The two tables have different row counts (`atm_table_nt=57`
vs `atm_table_ntc=76`), so the components are `allocatable`:
`atm_init` allocates each instance once (guarded by `allocated()`, so
`test_atm`'s repeated reads reuse the storage as the static arrays did)
with exactly the former extents and zero-fills it, which reproduces the
static-zero initial state -- an all `-999` row never sets `gmax_index`
and `surfp` then reads 0, as before. This is the one place where "same
shape, bounds and element order" needed allocatables to hold; each
instance keeps its former shape exactly. `kurucz_table_z` and
`atm_table_file_unit` stay bare module variables.

`diff -w surfp.f90 kcsurfp.f90` differed only in the table names,
`atm_table_nt` vs `atm_table_ntc` (now `tbl%num_teff`: loop bound,
clamp, `min(tbl%num_teff-3,row_base)`), header comments and line
continuation; the FORMATs 911/70/71 were already identical. The body is
now `surface_p_interp(tbl, ...)`; `surfp` and `kcsurfp` remain as the
entries `core/envint_kernel` and `test_atm` call, each passing its
instance. In `atm_init` the READ statements are unchanged apart from
names; the two identical "G Somers 5/15" edge scans became the internal
`find_gravity_edges(tbl)`.

## Item 2 -- dimension parameters (`64d5967`)

Wave 1 had already put a per-family `integer, parameter` block at the
top of `opacity_table_lib.f90` (opal95, opal92, alex94, alex06, kurucz,
laol) and `atm_table_lib.f90` (`atm_table_nt/ntc/ng/nta/nga`), and the
arrays are declared from them; no new parameters were needed.
Remaining literal copies replaced where provably that dimension:

- `kap/laol89/rdzlaol.f90:32`: header guard `.gt.104`/`.gt.52` ->
  `n_laol_rho`/`n_laol_t` (extents of `zlaol_opacity` and the grids;
  `rdlaol`'s identical guard already used the names).
- `kap/opal95/ll95tbl.f90`: missing-corner fill `(n,58,19)`, `do j =
  18,19`/`17,19`/`16,19`/`15,19`, `i = 70` -> `n_opal95_d`,
  `n_opal95_d-1..-4 .. n_opal95_d`, `n_opal95_t` (the READs on the
  same rows already used `n_opal95_d-k`).
- `kap/alex94/read_alex94_tables.f90`: `row_density_count.ne.17` ->
  `n_alex94_d` (the row's read list is `k=1,n_alex94_d`).

Left as literals (physical layout, thresholds or scratch sizes):
`ll95tbl`'s corner-row offsets 57/58/59..70, `build_kurucz_splines`'
`.lt.25` minimum, `rdlaol`/`setupopac`'s metal-mix arrays `(12)`/`(18)`,
`rdzlaol`'s `dummy(104)` scratch read target (11 elements used),
`alsurfp`'s `pressure_row(20)`.

126 vs 130 in `ll95tbl`: both numbers are real and both were already
named by wave 1 -- `n_opal95_tables = 126` (tables per file) and
`n_opal95_xz = 130` (slot count in the `(x,z)` table axis, with the
`opal95_table_start_index` offsets leaving unused slots). Nothing was
"fixed".

## Item 3 -- shared stencil kernels (`7b03580`)

`lagrange4` was used where the inline sum is
`w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4)` in that order with
contiguous unit-stride operands (or a `y(i0:i0+3)` slice along the
first dimension):

- `opal95_interp2d`: 3 T sums; `opal95_interp3d`: 3 T sums in the X
  loop (and, after `3c36ef5` made the X weights contiguous, its 3 X
  sums); `opal95_interp4d`: 3 Z sums.
- `getalex06`: 3 T sums; `alex94_interp3d`: 6 T sums + 3 X sums.
- `alex06tab`: the k=1..3 Z blend (tables `kk..kk+3` consecutive ->
  `alex06_full_opacity(kk:kk+3,i,j)`; the loop's dead `kk2/kk3/kk4`
  assignments dropped -- they are recomputed before their next use) and
  the final X blend over `opacity_by_x(1:4,i,j)`.

Left inline, each with a comment saying why:
- `opal95_interp4d` X sums: weights are the strided row
  `weight_x(k,1:4)`.
- `opal95_interp4d` T sums: the original accumulates three terms,
  stores, then adds the fourth (different association).
- `getalex06` / `alex94_interp3d` R sums: `y` is
  `alex06_opacity(ii,r:r+3)` / `alex94_opacity(8,ii,r:r+3)`, strided.
- `alex06tab` 4th-X Z blend: `kk4` is not `kk+3` in the Z=0.10 case.
- `opal95_fixed_z_table`: the four Z tables are `n_opal95_x` apart.

`stencil4_locate`: no new helper. Wave 1 already extracted
`stencil4_locate` (4 sites: `alex94_interp3d`, `getalex06` T and R) and
`stencil4_locate_opal95` (`getopal95` T and rho). The remaining
searches are pairwise different (see Deferred).

## Item 4 -- OPAL95 hidden-argument chain (`3c36ef5`)

`opal95_interp2d/3d/4d` now take every per-call value that `getopal95`
used to leave in `opacity_table` as an explicit `intent(in)` dummy of
the same shape: `index_t`, `index_rho(4)`, `weight_t(4)`,
`dweight_t(4)`, `weight_rho(4,4)`, `dweight_rho(4,4)`, `logr`,
`logr_lo_edge`, `logr_hi_edge(4)`, `extrap_lo`, `extrap_hi`,
`extrap_hi_row(4)`; `opal95_interp3d` additionally `index_x(4)`,
`weight_x(4)` (it only ever read row 1 of the X stencil, so
`getopal95` passes `opal95_index_x(1,:)`/`opal95_weight_x(1,:)`);
`opal95_interp4d` additionally `index_z`, `weight_z(4)`,
`index_x(4,4)`, `weight_x(4,4)`. The extrapolation flags/edges are
included because they are the same kind of per-call state (set fresh
in `getopal95` each call) even though the brief lists only
`logr/weight_*/index_*`. The `opacity_table%opal95_*` members remain
(`getopal95` still warm-starts its searches from the cursors); the
opacity tables and `opal95_table_start_index` are still read from
`opacity_table`. Inside the interpolators the change is a pure rename
of each read; the only expression whose text changed is the three X
sums in `opal95_interp3d`, now `lagrange4` (left-to-right, contiguous,
so identical). Alternate-return/error conventions and `kap_eval`'s
dispatch untouched.

## Deferred

Stencil searches not merged into a helper (none is token-identical with
any other after renaming; each line is the site's comparison / scan
range / fallback):

- `kap/alex94/alex94_fixed_z_table.f90:23-31`: `do i=3,n-1`,
  `x .lt. grid(i)`, `idz = i-2`; fallback `n-3`.
- `kap/opal95/opal95_surface_table.f90:37-45`: `do i=3,n-1`,
  `grid(i) .ge. x`, `i-2`; fallback `n-3`.
- `kap/opal95/opal95_fixed_z_table.f90:29-43`: `do i=3,n-1`,
  `grid(i) .ge. x`, `i-2`; fallback is an error return with a run-log
  message (no clamp).
- `kap/alex06/alex06tab.f90:47-55` (Z): `do i=3,n-2`, `x .le. grid(i)`,
  `i-2`; fallback `n-3`.
- `kap/alex06/alex06tab.f90:57-70` (X): `do i=3,n-2`, `x .le. grid(i)`,
  `i-2`; fallback `n-4` if Z >= 0.10 else `n-3`.
- `kap/opal95/getopal95.f90:185-209` (Z warm start): decision
  `x .gt. grid(cursor+2)`, up-scan `do i=cursor+3,n-1` with `x .lt.
  grid(i)` -> `i-2` (fallback `n-3`), down-scan `do i=cursor+1,2,-1`
  with `x .ge. grid(i)` -> `i-1` (fallback 1), then `min(cursor,n-3)`.
- `kap/opal95/getopal95.f90:258-274` (X warm start, in
  `opal95_x_stencil`): same shape as the Z warm start but without the
  final `min`, followed by the high-X missing-table overrides.
  A shared warm-start helper for these two would need the `min`
  failsafe hoisted out; it is only two sites, so not done.

Structural 2 beyond the class-A slice (the alternate returns in
`kurucz`/`gtlaol`/`opal92_*`, `kap_eval`'s dispatch) -- R6, not touched.

`sulaol.f90:36-43`: the log10 conversion of the SECOND LAOL table's T
grid loops to `laol(1)%num_t` (the first table's extent), as the
original did. Kept as is and flagged in a comment; whether it is a
latent bug when the two tables differ in T extent is an R6 author
decision (`sulaol.f90:33-36`).

## Reverted (changed numbers)

None. No pin or functional cross-check ever differed.

## Skipped

- No `kap_table_dims` block was created: the parameter blocks already
  exist per family at the top of the two `*_table_lib.f90` files from
  wave 1 and adding a second home for them would duplicate names.
- `check_boundaries.py` allowlist: no edit needed; every new call goes
  through the existing `*_lib` entries.
- Tests under `src/kap/test` / `src/atm/test`: unchanged; `test_atm`
  still calls `surfp`/`kcsurfp` (kept as wrappers) and `test_kap` does
  not touch the merged families.

## Wrong in the brief / audit / earlier reports

- The brief's "Known" list (and the wave-1 report it came from) says
  `opal92_interp3d.f90` has an X-grid `findex` block that
  `opal92_interp3d_z2.f90` lacks. Both files had it; the pair differed
  only in member names, message strings and a line continuation. No
  gate was needed.
- The brief says `surfp`/`kcsurfp` "differ in table length parameters":
  correct (`atm_table_nt` vs `atm_table_ntc`), and that is what forced
  the allocatable components.
- The audit's "13 hand-rolled 4-point window searches": after wave 1's
  helpers, seven sites remained and none pair up token-identically
  (list above), so the "shared `stencil4_locate`" recommendation could
  not be taken further without changing scan bounds or fallbacks.

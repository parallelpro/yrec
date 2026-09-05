# eos -- wave 3 (R5 structural) report

Worktree `/Applications/YREC-wt/eos`, branch `rs/eos`, base `yrec-modern`
6449462. Eight commits, all class A, each gate1 IDENTICAL before commit.

| # | hash | one line |
|---|------|----------|
| C1 | `1c9b1c29` | per-vintage OPAL EOS state instances opal95/opal01/opal06 |
| C2 | `0427922e` | merge quadeos01/quadeos06 into quad |
| C3 | `5cb3f143` | share the 2001/2006 T6-density interpolation (t6rint_core) |
| C4 | `bd9b7975` | share the 1995/2001 eqbound row search (eqbound_core) |
| C5 | `e56b0bef` | drop gmass06 (identical to gmass01) |
| C6 | `66a6478b` | name the per-vintage OPAL gate constants in eqstat2 |
| C7 | `5b243eb8` | composition-indexed MHD table state |
| C8 | `03c0291d` | named indices for MHD extension slots, SCV columns, mixture species |

Verification: gate1 IDENTICAL for every commit; full 37-case pin
selection after C8: 37 passed, PINS EXIT 0; aux battery: 22 passed, AUX EXIT 0.
The pins and aux exercise only the OPAL 2006 path, so the 1995/2001
OPAL vintages and the MHD path (no MHD tables ship with YREC;
`test_eos`'s MHD section runs only with `YREC_MHD_TABLES` set) were
covered by two private probe programs (scratch, not committed) compiled
against the pre-change and post-change trees and compared with `cmp`:

- OPAL probe: `oeqos`/`oeqos01`/`oeqos06` + `eqbound*` over
  6 X x 91 log T x 25 log P points per vintage (648 of the 1995 points
  in the stencil-overrun regime described under item 1), every
  `eos_output` slot and the run log printed in hex -- bit-identical for
  all three vintages after every commit C1-C8.
- MHD probe: synthetic unformatted tables in the `rtab`/`rabu` format
  for ZAMS a/b/c (lower + upper) and centre 1..5 (three composition
  passes each), `mhdst(40..47)`, then `mhdpx` over 5 X x 41 log T x 21
  log P points spanning the lower-ZAMS, upper-ZAMS and centre regions
  (480 / 1355 / 830 non-zero points), `mhd_output(1:25)` in hex --
  bit-identical after C7 and C8.

## Item 1 -- `opal_eos_vintage` derived type (C1-C5)

`src/eos/opal_eos_lib.f90`: `opal_eos_state`'s three parallel member
sets (unsuffixed / `_01` / `_06`) became per-vintage instances. Two
state types, not one, for a numerical reason the audit did not know:

- `type opal_eos_vintage` at the 2006 dimensions (`n_eos_nr_max =
  n_eos06_nr = 169`, `n_eos_nt_max = n_eos06_nt = 197`), instances
  `opal01`, `opal06`. Checked before choosing it: every esac01 /
  t6rinteos01 loop and index runs over the declared `n_eos01_nr/nt`
  parameters (169/191), never the physical bound; the 2001/2006
  stencils are clamped (`t6_order_hi`/`density_order_hi`); no
  `size`/`ubound`/`shape` or whole-array operation reads any state
  array. So the larger physical arrays change nothing (probe-verified).
- `type opal_eos_vintage95` at the exact 1995 dimensions (77 x 56) and
  in the previous member order, instance `opal95`. `esac.f90`'s
  boundary-sum stencils overrun the table (`t6_index_1+3 = nt+1`,
  `x_index_lo+3 = 6`) whenever T6 lies in the lowest table interval,
  and `t6rinterp` then reads `t6_grid(nt+1)` = the member that follows
  `t6_grid` (`x_interp_result(1,1)`) and `t6_grid_spacing_inv(nt+1)`
  (= `density_grid_spacing_inv(1)`). Re-dimensioning or re-ordering
  would change 1995 numbers. The header comment explains this; it is
  an R6 item for the owner (below).
- Explicit `nr`/`nt`/`nvar` count components were NOT added: every
  loop uses the file-local `integer, parameter :: nt = n_eosXX_nt`
  already, and adding a component that nothing reads would be a
  half-done change. Listed under Deferred.
- The read-only tables that were declaration-time initialisers
  (`x_grid`, `eos_var_order`, `t6_index_lo`, `density_index_edge*`)
  are named constants of the same shape and literals
  (`opal95_x_grid`, `opal01_x_grid`, `opal06_x_grid` keep the
  single-precision literals; `opal_eos_var_order` is one shared identity
  permutation; `opalXX_t6_index_lo`, `opalXX_density_index_edge`).
  `readco_init_flag` and `table_metal_fraction` are 3-element arrays of
  `opal_eos` indexed by `iv_opal95/01/06 = 1/2/3`, keeping their `= 0`
  initialisers. The vintage instances are module siblings of `opal_eos`
  rather than components so they stay initialiser-free (bss): gfortran
  puts a whole module object into `__DATA` if any component has a
  default initialiser; `yrec` shrank from 112.8 MB to 45.2 MB of
  initialised data. (Deviation from the brief's "instances inside
  `opal_eos_state`"; Deferred.)
- `i_opal_*` vs `i_opal01_*` index blocks stay, documented as the layout
  descriptor of each vintage.

Instance plumbing: `esac/readco/rhoofp/radsub/t6rinterp` (+01/+06) take
`type(opal_eos_vintageXX), intent(inout) :: v` first; `eqbound*`/`oeqos*`
(called by name from eqstat2) and `eos_lib`/`test_eos` bind
`opal95/opal01/opal06` directly. `xxh` is referenced nowhere (kept,
commented).

Clone families, each `diff`ed with names normalised (full diffs in the
commit messages):

- **quad / quadeos01 / quadeos06** (C2): only the cache member names and
  a redundant pair of parentheses around the single variable
  `x1_plus_x2(cache_slot)` in the 2006 copy. One `quad(cache, ...)` with
  `type(opal_quad_cache), intent(inout) :: cache` (`opalXX%quad`); the
  two clones `git rm`'d. Wave-1's "needs the cache derived type first"
  was right; done here.
- **t6rinteos01 / t6rinteos06** (C3): differ only in the RECURSIVE
  attribute and the out-of-range message text. New
  `opal/t6rint_core.f90` holds the 64-line kernel; each wrapper keeps its
  attribute, message and `ierr = 1`. **t6rinterp (1995) not merged**: its
  esactq2 block is ordered differently (and it is the routine that
  performs the overrun read above).
- **eqbound / eqbound01** (C4): the row search / edge / ramp tail (71
  lines) is identical; the density pre-check differs (`.gt.5.0d0` vs
  `.gt.7.0d0`) and stays in each wrapper; `opal/eqbound_core.f90` takes
  explicit-shape `t6_grid(nt)`, `density_edge_at_t(nt)`,
  `density_index_edge(nt)`, `density_grid(nr)`. Sequence association:
  eqbound01 passes `opal01%t6_grid` (physical length 197) to a dummy of
  length 191 -- a contiguous leading slice, the only rows the 2001 code
  addressed. **eqbound06 not merged**: no `ierr`, and different
  structure.
- **gmass01 / gmass06** (C5): identical apart from comments; gmass06
  `git rm`'d, esac06 calls gmass01. **gmass (1995) stays** (different
  ionisation energies, electron-mass term).
- **rhoofp / rhoofp01 / rhoofp06, radsub / radsub01 / radsub06**: NOT
  merged, instance argument only. rhoofp01 has the private
  `density_index_edge` DATA table differing from the module's (wave 1
  documented it); rhoofp06 omits the radiation term from
  `pressure_max/min` and passes `0` instead of `rad_flag` to esac06;
  rhoofp (1995) uses a different table layout. radsub06 has the
  amu/log10ne column difference the audit noted; radsub/radsub01 differ
  in slot layout (`i_opal_*` vs `i_opal01_*`).
- **esac* / readco* / oeqos***: thin vintage wrappers, instance
  argument only, as the brief asked.

## Item 2 -- `blend_opal_result` out of `eqstat2` (C6)

Not extracted: after substituting vintage names the three blocks still
differ in the gate expression itself (1995: `log10_temperature.le.8.0d0`,
`log10_density.le.5.0d0`; 2001: `temperature.le.100d6`, `.le.7.0d0`;
2006: `temperature.le.200d6`, `.le.7.0d0`) and in `eqbound06` taking no
`ierr`, so they are not token-identical. The identical remainder was
already the wave-1 `adopt_or_blend_opal_result()` contained procedure.
Done instead: the nine gate literals are local `double precision,
parameter`s with exactly the former values and kinds
(`opal95_t_min = 5.0d3, opal95_log10_t_max = 8.0d0,
opal95_log10_rho_max = 5.0d0; opal01_t_min = 2.0d3, opal01_t_max =
100d6, opal01_log10_rho_max = 7.0d0; opal06_t_min = 1.870d3,
opal06_t_max = 200d6, opal06_log10_rho_max = 7.0d0`) with a comment
saying the three gates are different tests. beta14 untouched.

## Item 3 -- composition-indexed MHD state (C7)

`src/eos/mhd_eos_lib.f90`: the per-composition members became one array
each with a trailing composition index -- `zams_lower_table(nt1m,nr1m,
ivarc,3)`, `zams_upper_table(nt2m,nr2m,ivarc,3)`, `zams_atomic_weight /
_number_abundance / _mass_fraction(nchem0,3)`,
`zams_mean_molecular_weight(3)`, `centre_table(ntxm,nrxm,ivarx,5)`,
`centre_*(nchem0,5)`, `centre_mean_molecular_weight(5)` -- with
`mhd_n_zams = 3`, `mhd_n_centre = 5`, `imhd_zams_a/b/c = 1/2/3`,
`imhd_centre_1..5 = 1..5`. Each slice `(:,:,:,ic)` / `(:,ic)` / `(ic)`
has the shape and element order of the member it replaced. (The brief's
sketch `imhd_zams_a_lo = 1, imhd_zams_a_hi = 2, ...` would put lower and
upper ZAMS tables, which have different shapes, in one array; the
lower/upper split is kept as two arrays and the index is the
composition.)

`src/eos/mhd/mhdst.f90`: the eight unit dummies are copied
element-by-element into `unit_zams(3)`/`unit_centre(5)` and the eight
`mhdst1` calls are two loops, `do ic = 1, mhd_n_zams` then `do ic = 1,
mhd_n_centre`. (a) Loop order = former statement order (a, b, c, centre
1..5). (b) Every actual is a contiguous section of the same length as
the explicit-shape dummy (full leading ranges; `nm` on `mhdst.o` and
`mhdpx2.o` shows no `_gfortran_internal_pack`), the scalar
`mean_molecular_weight(ic)` is an element as before. The `if
(unit.gt.0)`, `table_index = 0/1` and `if (ierr /= 0) return` inside the
loop are token-identical; the `mhdst` interface still takes the same
eight units (io/read_controls.f90's 40-47); the TEMPERATURE LIMITS block
is unchanged.

`src/eos/mhd/mhdpx2.f90`: the 11-way ladder is three range branches:
selector in [-3,-1] -> lower ZAMS composition `-selector`, row
`-selector`; [1,3] -> upper ZAMS composition `selector`, row `selector`;
[4,8] -> centre composition `selector-3`, row `selector`. (c) This is
exactly the table/row/mass-fraction triple each ladder branch passed;
any other selector still falls through and does nothing.
`interpolate_and_store` and its `intpt` call are unchanged.

## Item 4 -- named indices (C8)

- `mhd_output(25)`: wave 1 had already named the slots meqos reads
  (`i_mhd_*`), which the brief did not reflect. Added the ones provable
  from the only writer, `mhdst1`: `i_mhd_log10_u = 3` (the source of
  slot 22's derivative, per mhdst1's comment "log10(P)/DX, log10(U)/DX,
  DDELAD/DX, log10(CP)/DX"), `i_mhd_dlog10p_dx = 21`, `i_mhd_dlog10u_dx
  = 22`, `i_mhd_dgrad_ad_dx = 23`, `i_mhd_dlog10cp_dx = 24`,
  `i_mhd_placeholder = 25`; mhdst1's five literal writes and four literal
  source slots use them. Slots 6, 7, 17, 19 have no reader or writer that
  fixes their meaning; left literal (nothing addresses them).
- SCV columns (`src/eos/scv_eos_lib.f90`): `iscv_*` for the H/He tables
  (1 log10 P, 2 neutral-molecule/neutral-He fraction, 3 neutral-H/He+
  fraction, 4 log10 rho, 5 log10 S, 7 dlnrho/dlnT, 8 dlnrho/dlnP, 9
  dlnS/dlnT, 12 du/dT built by scv_envelope_table), `iscvz_*` for the Z
  table (4, 7, 10, 13) and `iscvenv_*` for all 12 envelope-table
  columns (writer scv_envelope_table; reader eqscve's `interp_result(1..5)`
  = columns 2..6). Every literal third subscript of
  `tablex/tabley/tablez/tablenv` in `eqscvg.f90`, `eqscve.f90`,
  `scv_envelope_table.f90` now uses them (eqscve's `j+1` loop over
  columns 2..6 left as a loop). H/He columns 6, 10, 11 and Z columns
  1-3, 5, 6, 8, 9, 11, 12 are read by nothing and stay unnamed.
- 12-species mixture (`src/eos/eos_mixture_lib.f90`): `ix_na = 1 ..
  ix_he = 12, n_mix_species = 12` -- the same names and values as the
  state agent's block in `state/star_info_lib.f90` (which the star layer
  uses to fill `fxenv`). Duplicated rather than `use`d so
  `eos_mixture_lib` and `yale/saha_eos.f90` stay free of `star_info_lib`
  (the module header's physics-purity contract; same precedent as
  net_lib's `n_reactions`). Used in eqstat2 (`atomic_weights_full(ix_h)`,
  `(ix_he)`, `saha_mass_fractions(ix_c/ix_h/ix_he)`, the metal loops
  `ix_na..ix_c` and `ix_o..ix_ne`, array dimensions), `eos_lib`
  (`species_fractions(n_mix_species)`) and `saha_eos`
  (`saha_mass_fractions(ix_he)`). Listed under Deferred as the
  duplication the brief anticipated.

C8 evidence beyond gate1: `otool -tv` disassembly of every changed
object is identical to the pre-change build except for the source line
numbers gfortran embeds in I/O statements (eqstat.o, saha_eos.o,
mhdst1.o gained 2-3 lines of `use`/comment); eqscvg.o, eqscve.o,
scv_envelope_table.o, eos_lib.o, mu.o are md5-identical.

## Deferred (exact proposals)

- `opal_eos_lib.f90:` `nr`/`nt`/`nvar` count components on
  `opal_eos_vintage` (brief item 1). Proposal: add `integer :: nr, nt`
  set once in readcoeos01/06 from `n_eos01/06_nr/nt`, and switch the
  `integer, parameter :: nt = n_eosXX_nt` in esac01/06, t6rinteos*,
  rhoofp01/06, eqbound01/06 to `v%nt` -- a loop-bound variable in place
  of a parameter (same values; gfortran may unroll differently, so
  verify with the OPAL probe).
- `opal_eos_lib.f90:` instances as siblings (`opal95`, `opal01`,
  `opal06`) instead of components of `opal_eos_state`. Making them
  components requires dropping the `= 0` initialisers of
  `readco_init_flag`/`table_metal_fraction` (moving them to an init
  call) or accepting the 112.8 MB `__DATA` image; author decision.
- `opal/t6rinterp.f90` vs `t6rint_core.f90`: the 1995 esactq2 block
  ordering; merging changes 1995 numbers (R6, see below).
- `opal/rhoofp01.f90:53` private `density_index_edge` DATA table vs
  `opal01_density_index_edge` in `opal_eos_lib.f90`: unchanged (wave 1
  skipped it, still B).
- `opal/eqbound06.f90`: no `ierr`; giving it one and sharing
  `eqbound_core` changes eqstat2's call structure (wave-1 item 16).
- `eqstat.f90:` the three OPAL vintage blocks share only
  `adopt_or_blend_opal_result`; a `blend_opal_result(...)` taking the
  gate as an argument would need the gate expressions harmonised
  (log10 T vs T) -- B.
- `eos_mixture_lib.f90:` `ix_*`/`n_mix_species` duplicate
  `state/star_info_lib.f90:56-59`. To de-duplicate with no numeric
  effect: replace the parameter block with `use star_info_lib, only:
  ix_na, ..., n_mix_species` plus `public :: ix_na, ...` (names already
  identical), at the cost of a `star_info_lib` dependency in
  `eos_mixture_lib` and `saha_eos`.
- `mhd/mhdpx1.f90`, `mhd/mhdst1.f90`: `saha_eos.f90`'s
  `data saha_weight_term/13 values/` and `nz/11/` (ionisation stages per
  species 1..11 + He twice) are indexed by the same species order as
  `ix_*`, but the 13th slot's species is not provable from the code;
  left literal.
- `scv/scv_envelope_table.f90:` `tablenv` columns 7-12 (the numerical
  derivatives) and `tablex/tabley` column 12 (du/dT) are written and read
  by nothing in `src/`. Candidate dead code (B to remove because the
  computation writes nothing else; safe to delete only after an owner
  confirms no external reader).
- `scv_eos_lib.f90:` `tablenv` is dimensioned `scv_nvar = 12` and
  `tablex/tabley` too, though the files supply 11 columns; the 12th is
  the derived du/dT. Not changed.
- `mhd/mhdst1.f90:159` comment "SPACE-HOLDER VARIABLE (LIKE VAR(20))"
  says slot 20 is a placeholder, whereas wave 1 named slot 20
  `i_mhd_log10_pgas` from meqos's use (`beta = exp10(out(20) -
  out(2))`). Either the comment is stale (MHD table format later filled
  slot 20) or meqos reads a placeholder; owner should check against the
  MHD table documentation. Not changed.

## Reverted (changed numbers)

None. (One private-probe segfault after C1 was an argument-count
mismatch in the three `oeqos*` `rhoofpXX(` calls missed by the rewrite
script; fixed before the commit, never committed.)

## Skipped

- Brief item 1, a single `opal_eos_vintage` type for all three vintages:
  the 1995 overrun read makes 1995's dimensions/order load-bearing (see
  Item 1); two types instead.
- Merging rhoofp*/radsub*: real differences (listed under Item 1).
- Brief item 2 extraction: blocks not token-identical (Item 2).
- Item 4, `mhd_output` slots 6, 7, 17, 19 and the unread SCV columns:
  not provable.

## What the audit / wave-1 report got wrong or missed

- Audit Structural 1 assumes the three OPAL state sets can share one
  type at the largest dimensions. For 1995 they cannot: `esac.f90`'s
  boundary sums (`esac.f90:205-231`) run `t6_scan_idx` to
  `t6_index_1+3 = nt+1` and `x_loop_index` to `x_index_lo+3 = 6`
  (`n_eos_mx = 5`; the X overrun lands on `eos_table(1,2,...)`, inside
  the array but the wrong slot), and `t6rinterp.f90` then reads
  `t6_grid(nt+1)` / `t6_grid_spacing_inv(nt+1)`, i.e. the first element
  of the member that follows each, whenever T6 is in the lowest table
  interval (5000 K <= T < 5500 K). Out-of-bounds reads that happen to
  land on neighbouring state; a pre-existing R6 item (fixing it changes
  1995 numbers in that regime). The 1995 path is not pinned; the private
  probe shows 648 of its 13 650 points are in that regime.
- Audit Structural 5 / brief item 4 list `mhd_output(25)` as unnamed;
  wave 1 (`i_mhd_*` in `mhd_eos_lib.f90`) had already named every slot
  meqos reads. Only the mhdst1-side slots were left.
- Wave-1 report "Skipped": "mhdst1's slot-21..25 fills: index variables,
  not bare literals" -- they were bare literals (21..25 and source 2, 3,
  8, 9); named in C8.
- Brief item 3's index sketch (`imhd_zams_a_lo = 1, imhd_zams_a_hi = 2`)
  mixes two differently-shaped table families in one index space; the
  composition index with separate lower/upper arrays is what the code
  structure supports without copies.
- Binary size: with any component default-initialised, gfortran emits
  the whole module object into `__DATA`; the old `opal_eos_state`'s
  DATA-initialised tables cost 67 MB of on-disk image. Not an audit
  finding; fixed as a side effect of C1.

## Verification

Run once after the last edit (C8 `03c0291d`), so they cover C1-C8
together; the tree at that point is byte-identical to the committed
head (`git status --short` empty).

- `tools/gate1.sh /Applications/YREC-wt/eos`: IDENTICAL before each of
  C1-C8 (make warnings in `logs/eos.make.log` unchanged from the base:
  the three pre-existing `-Wsurprising` stack-to-static notes at
  `mhdst.f90`).
- `tools/fullpins.sh /Applications/YREC-wt/eos pins`:
  `37 passed, 272 deselected in 1555.61s (0:25:55)`, `PINS EXIT 0`.
- `tools/fullpins.sh /Applications/YREC-wt/eos aux`:
  `22 passed in 546.41s (0:09:06)`, `AUX EXIT 0`.
- No `test_*.short` files remain; `git status --short` is empty (no
  output dirs, `.short` files or `standard/` committed).

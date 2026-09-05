# eos -- wave 1 readability (R1-R3)

Worktree `/Applications/YREC-wt/eos`, branch `rs/eos`, now based on
`yrec-modern` at `272d1e2` (the branch was rebased there by the
coordinator after the three commits were made on base `0d5c11a`; the
original hashes were a830f4d / ae5adf4 / c2ff5ec, picked cleanly with
identical file stats). Verification below refers to the pre-rebase
builds; the rebased tree has not been re-pinned.

| batch | commit | verification |
|---|---|---|
| R1 comments & dead code | `227953e` | gate1 IDENTICAL; 37-case pins `PINS EXIT 0` (`logs/eos.pins.R1.out`) |
| R2 named indices/constants | `4fba3c8` | gate1 IDENTICAL; 37-case pins `PINS EXIT 0` (`logs/eos.pins.R2.out`) |
| R3 shared helpers | `c764c46` | gate1 IDENTICAL; 37-case pins `PINS EXIT 0` (`logs/eos.pins.R3.out`); aux battery 22 passed `AUX EXIT 0` (`logs/eos.aux.out`) |

Every build went through gate1 (`-O3 -ffp-contract=off`, USE_HDF5=1);
`make clean` preceded the R3 gate run because scv_eos_lib gained a module
procedure. `src/deps.mk` was regenerated in R1 (removed `use`
statements); R2 and R3 added no `use` statements. The only compiler
warnings from eos files are the four pre-existing `-fmax-stack-var-size`
notes in mhdst.f90. Nothing under `/Applications/YREC` or any other
worktree was touched; `input/` and `audit/` are untouched.

(The R1 commit message carries "Assisted by Claude" twice -- a paste slip
in the message file; not amended because R2/R3 sit on top of it.)

## R1 -- comments and dead code (`227953e`, 45 files + deps.mk)

Headers and notes rewritten to describe the code as it is now (no COMMON
names, SAVE claims, deleted files or stub STOPs): eos_lib, eqstat, mhdst,
mhdst1, mhdtbl, readco/readcoeos01/readcoeos06, oeqos/oeqos01/oeqos06,
quad, rhoofp/rhoofp01/rhoofp06, gmass01/gmass06, t6rinterp,
fully_ionized_eos, scv_eos_lib, opal_eos_lib, mhd_eos_lib. Commented-out
F77 residue removed across mhd/, opal/, scv/, yale/ (READ/WRITE/FORMAT/
STOP/PAUSE/DATA/DIMENSION lines, "KC 2025-05-3x fixed ..." tags, F77 label
comments). Deliberately preserved behaviours are now documented where
they live: the no-D-suffix DATA literals in opal_eos_lib x_grid_01/06 and
saha_eos (area item 24, option A), radiation pressure not subtracted in
rhoofp*, the dE/dRho slot not rescaled in esac, the non-SAVE'd `dix` in
t6rinterp, rhoofp01's private edge table differing from the module's.

Dead code, each shown dead by grep over the whole `src/` tree or by
reading the control flow (details per item in the commit message):

- esac/esac01/esac06: write blocks between `return` and `end if`
  (goto-61/62/65/66 emulation) are unreachable; duplicated
  `x_index_3 = 1` deduplicated.
- eqbound/eqbound01: `t_fraction` computed, never read (only
  density_fraction feeds ramp_factor); dead "error exit" text after the
  final return. eqbound06: dead tail after return.
- mhdst1: two ierr=1/return blocks after an unconditional return; mhdst:
  empty `if (unit_centre1.gt.0) then / end if`.
- gmass01/gmass06: `electron_mole_excess` assigned, never read.
- meqos: dead tail after return; `ier_flag` never read.
- eqscvg: `total_entropy` assigned, never read.
- scv_envelope_table: first `log_rho_mix` assignment overwritten before
  any read; `dlnp_dlnt_gas` and the `dlncp_dlnp` loop+sum never read
  (only dlncp_dlnt reaches tablenv(..,10)).
- saha_eos: `rmub`; readco: `t6_count_used`.
- Unused locals in oeqos*/esac*/readcoeos01/rhoofp*/eqscve (blank_line,
  hydrogen_fraction_copy, density_copy, fill_idx, ivarx, cnvs, zero,
  metal_fraction_table, helium_fraction, gas_pressure, nps) and the
  unused `mx/mv/nr/nt` local parameters (area item 25).
- Unused imports: `use star_info_lib` (oeqos, oeqos01, oeqos06, eqscvg,
  meqos), `use luout_lib` / `use scv_eos_lib` / `use eos_mixture_lib` in
  eqstat's wrapper, `use phys_const_lib` / `use scv_eos_lib` in
  eos_eval/eos_get_gamma1, `use luout_lib` in mhdpx. Every `use math_lib`
  kept.
- `continue`/blank no-op runs before returns collapsed (rhoofp*,
  eqbound*, mhdst1, saha_eos, eqstat).
- `intent(in)` for log10_temperature/log10_pressure in saha_eos and
  fully_ionized_eos (area item 8; never assigned in either).

SUMMARY 1.1 items: none apply to eos files. #3 (rhoofp01 priming flag)
is number-changing and stays for R6.

## R2 -- named indices and constants (`4fba3c8`, 35 files)

- opal_eos_lib: `i_opal_p, i_opal_e, i_opal_s, i_opal_dedrho, i_opal_cv,
  i_opal_chi_rho, i_opal_chi_t, i_opal_gamma1, i_opal_gamma2_ratio,
  i_opal_gamma3m1` (1-10; slots of eos_output/eos_output_06 and the
  arguments of eos_index_inverse/_06) and `i_opal01_*` (1-9, the 2001
  layout has no dE/dRho); `opal_rho_not_found = -999.0d0`;
  `opal_flag_set = 12345678` (area items 10, 11 naming half, 12 naming
  half).
- oeqos/oeqos01/oeqos06, rhoofp*, eos_lib (eos_get_gamma1): every
  `eos_output*(N)` uses the names; esac*/radsub*: every
  `eos_index_inverse*(N)`; rhoofp*: all 26 `-999.0d0` returns ->
  `opal_rho_not_found`; esac*/readco*/rhoofp*: `12345678` ->
  `opal_flag_set`. No literal slot/sentinel remains in code.
- eqbound*, esac*, readco*, rhoofp*: the local `mx/mv/nr/nt` parameters
  are defined from `n_eos_mx/n_eos_mv/n_eos95_*/n_eos01_*/n_eos06_*`
  (same names at the use sites, one source of truth).
- mhd_eos_lib: `i_mhd_log10_rho=1, i_mhd_log10_p=2, i_mhd_chi_rho=4,
  i_mhd_chi_t=5, i_mhd_grad_ad=8, i_mhd_log10_cp=9,
  i_mhd_dgrad_ad_dlog10_rho=10, i_mhd_dgrad_ad_dlog10_t=11,
  i_mhd_dlogcp_dlogrho=12, i_mhd_dlogcp_dlogt=13, i_mhd_ion_frac_1=14,
  i_mhd_eta=18, i_mhd_log10_pgas=20`, with a note that the meanings are
  inferred from meqos's use (the MHD table documentation is not in the
  repository); applied in meqos and mhdpx (area item 9). mhdst, mhdpx1,
  mhdpx2: local `ivarc/ivarx/nchem0/nt*m/nr*m` defined from the `mhd_*`
  module parameters.
- scv_eos_lib: `scv_nt = 63, scv_np = 76, scv_nvar = 12, scv_nvar_z = 13`
  in all six table declarations; the `nts = 63` locals in eqscve, eqscvg,
  scv_envelope_table and eos_init are `nts = scv_nt` (area item 21).
- Renames: `atomic_weights -> inverse_atomic_weights` in eqstat and mu
  (`atomic_weights_full` unchanged; item 5); `species_mass_fraction ->
  species_number_fraction` in gmass*/esac* plus a slot comment (item 6);
  scv_envelope_table `helium_fraction_local` split into
  `helium_fraction_hhe` (1 - X, entropy-of-mixing loop) and
  `helium_fraction_local` (1 - X - Z) (item 7); `interp_x` split into
  `interp_x` (inter3 abscissae) and `interp_f` (values combined with the
  inter3 weights) (item 22, rename half).

## R3 -- shared helpers (`c764c46`, 7 files)

Every helper body is token-for-token the call-site text; no operand
order, parenthesisation or literal changed. Host-associated `contains`
procedures inherit the host's `use math_lib`, so `exp10` in
adopt_or_blend_opal_result resolves as before.

- eqstat.f90: the three identical adopt-or-ramp blocks of the 1995/2001/
  2006 branches in eqstat2 are one contains procedure
  `adopt_or_blend_opal_result()` (area item 1); each branch keeps its
  gate, oeqos*/eqbound* calls and comments.
- eos_lib.f90: eos_get makes one eos_eval call, passing the optional
  `composition_at_zone` through whether or not it is present (item 19);
  the 3-wide `i_fxion` slot is documented at the index block.
- mhd/mhdpx2.f90: contains procedure `interpolate_and_store(table, dims,
  nvar, log10t, num_t, num_r, row, hfrac)` replaces the eleven
  intpt-call + copy-loop + hfrac-store bodies; the selector if-chain and
  per-table arguments are unchanged (structural item 3, the part that
  needs no state change).
- scv/eqscve.f90: contains procedure `spline_cell(t_shift, p_shift,
  result)` replaces the four pressure-then-temperature spline blocks
  (home, T-shifted, P-shifted, both-shifted; item 3). The blocks differed
  only in integer index offsets, now the two arguments.
- scv_eos_lib.f90: module function `scv_weighted_sum4(w, y)` =
  `w(1)*y(1) + w(2)*y(2) + w(3)*y(3) + w(4)*y(4)`; replaces 46 inline
  4-point Lagrange sums (38 in eqscvg, 8 in eqscve; item 4). Call sites
  pass the table column as an array section (`tablex(ii,idp:idp+3,c)`,
  `temp_work(1:4,m)`). The log10-weighted `dlnsmix_dlnt` sum in eqscvg
  has a different form and stays inline.
- eos_mixture_lib.f90: comment giving the fxenv(12) species order (Na,
  Al, Mg, Fe, Si, C, H, O, N, Ar, Ne, He), read off eqstat2's
  atomic_weights_full and core/read_starting_model's setup (structural
  item 5, mixture part; comment only -- no i_* names added because the
  only numeric uses are eqstat2's 6/7/12 and loop bounds).

## Deferred (cross-domain)

1. `mhd/mhdst.f90` optional-centre-table branches are dead because
   `core/read_controls.f90:1431-1438` hard-codes units 40-47. Proposal:
   in read_controls drop the eight hard-coded unit assignments in favour
   of the namelist values (or document that they are always set), after
   which mhdst's `if (unit_centreN.gt.0)` guards can go. Not touched:
   removing the guards alone would change behaviour if a control ever
   set a unit to 0.
2. `eos_lib.f90` `eos_get_gamma1` has no production caller (audit); only
   `src/eos/test/test_eos.f90:233,254` use it. Proposal: either delete
   it together with the test sections, or move it out of the
   `check_boundaries.py` public allowlist and mark it test-only in the
   `public ::` line. Left as is because the test file's expectations and
   the allowlist (`tools/check_boundaries.py`) are outside the batch
   scope; its slot reads now use `i_opal_gamma1/i_opal_gamma2_ratio`.
3. `mhd/rabu.f90:31-33` bounds check after the implied-do read
   (overrun before the check fires). Proposal: read the count first,
   test it, then read the array -- a two-statement split of the READ.
   This changes the read sequence on malformed tables only, so it is a
   B-class edit; not done here.
4. `eqstat.f90` / `mu.f90` ion-mean-weight block (item 20): a shared
   function in eos_mixture_lib would need to return dfx1, dfx12, dfx4 and
   the envelope-branch flag because eqstat's copy interleaves the
   saha_mass_fractions setup in both branches. Proposal for R4:
   `subroutine ion_mean_weight_excess(x, z, mu_ion_inv, dfx1, dfx12,
   dfx4, use_envelope)` in eos_mixture_lib with the exact statement
   order of mu.f90; eqstat then reads the returned dfx values.

## Reverted (changed numbers)

None. Every gate1 run reported IDENTICAL on the first build and all
three full pin runs and the aux battery passed on the first attempt.

## Skipped (disagree with reviewer / too risky)

- item 2 (quad/quadeos01/quadeos06 dedupe): the three differ in the
  cache member names (`q/h/s/...` triples with `_01`/`_06`), so a shared
  routine needs the cache derived type of structural item 1 first.
- item 11 (alternate returns / `-999` sentinel -> status codes): the
  sentinel now has a name; replacing `*999`/`return 1` in esac*/rhoofp*/
  oeqos*/eqstat2 changes every caller's control flow -- R4 material.
- item 12 (readco 12345678 guard "redundant"): kept, now named
  `opal_flag_set`; the guard also protects the direct readco* call in
  test_eos.
- item 13 (zero-T6 row: 1995 exits the loop, 01/06 raise ierr):
  harmonising changes the 1995 behaviour on a malformed table (B).
- item 14 (T6 validity gates): the oeqos* gates are in T6 and the
  eqstat2 gates in K / log10 T with different constants, so there is no
  token-identical single copy; naming both sets would add six parameters
  for three comparisons each.
- item 15 (1995 pressure consistency check kept, 01/06 commented out):
  enabling or removing changes behaviour (A/B decision for the owner).
- item 16 (eqbound06 has no ierr; eqbound writes to `*`): aligning the
  signature touches eqstat2's call and the unit change alters where a
  message lands.
- item 17 (mhdst1/rabu/rtab silent ierr paths) and item 18
  (esac01/06 `deriv_order.gt.9` warn-only): adding messages or an ierr
  return changes output/behaviour on those paths.
- item 23 (radsub06 `rad_flag.ne.0` always true): left; removing the
  branch is safe but the duplicated comment block is the only
  documentation of the two formulas.
- SUMMARY 1.1 #3 (rhoofp01 priming flag) and the rhoofp01 edge-table
  discrepancy vs the module's `density_index_edge_at_t_01`: both change
  numbers (R6/B).
- scv_envelope_table :290-291 duplicate `interp_f(k_idx)` overwrite in
  the DU/DT loop: audit "known, skip" (a possible bug -- the first
  assignment clamps jj to the row length, the second does not); left
  verbatim.
- scv_envelope_table `t_interp_weight`/`p_interp_weight` (item 22,
  "delete dead weights"): they are inter3 output arguments, not
  deletable without changing inter3's interface (numerics domain).
- mhdpx1 `mhd_output(iv)` loop and mhdst1's slot-21..25 fills: index
  variables, not bare literals; left as is.

## Things the audit missed (fixed in R1 unless noted)

- eos_lib: eos_eval/eos_get_gamma1 headers still described the old
  error-funnel and COMMON-era table structure.
- eqbound/eqbound01: `t_fraction` dead (the audit lists it in the
  summary but not in the ranked items).
- scv_envelope_table: `dlnp_dlnt_gas`, `dlncp_dlnp` (loop + sum) and the
  first `log_rho_mix` assignment dead.
- mhdst: empty `if` block; mhdst1: two unreachable ierr blocks.
- t6rinterp and readcoeos06: comments claimed locals were SAVE'd; they
  are not (they hold `-finit-local-zero` values).
- saha_eos `rmub`, meqos `ier_flag`, gmass01/06 `electron_mole_excess`,
  eqscvg `total_entropy`: write-only locals.
- rhoofp01's private edge table differs from opal_eos_lib's (documented,
  not changed).
- gmass*: `species_number_fraction(7)` is never assigned (documented).
- eos_mixture_lib: fxenv species order was undocumented anywhere in the
  eos domain (comment added in R3).

# Readability wave 3 (R5, structural) -- rotation domain

Worktree `/Applications/YREC-wt/rotation`, branch `rs/rotation`, base
`6449462` (yrec-modern with waves 1, 2a, 2b merged).
Brief: `tools/wave3_rotation.md`; rules: `tools/wave3_prompt.md`.
Area report: `audit/readability-sweep-2026-09-03/rotation.md`.

| Item | Commit | Verification |
|------|--------|--------------|
| 1 equal-grid toolkit | `2e92a80` | gate1 IDENTICAL; cm_rot_base, cm_rot_disk, cm_rot_kawaler, cm_rot_settle_he, cm_rot_solid byte-identical |
| 2 dead-code purge (leftovers) | `4c9509d` | gate1 IDENTICAL; cm_rot_base, cm_rot_settle_he, cm_settle_he_z byte-identical |
| 3 species tables | `f8bced8` | gate1 IDENTICAL; cm_settle_he, cm_settle_he_z, cm_rot_settle_he, cm_rot_base byte-identical |
| 4 am_advection_diffusion_coeffs | (no edit needed; see item 4) | -- |
| 5 solve_r0 comments | `1830868` | gate1 IDENTICAL; cm_rot_base, cm_rot_solid byte-identical |

After the last edit (`1830868`): full 37-case pin selection 37 passed,
`PINS EXIT 0` (1582 s); aux battery 22 passed, `AUX EXIT 0` (569 s).

Flags throughout: `-O3 -ffp-contract=off -finit-local-zero -Wall`, `USE_HDF5=1`.
`make clean` was run at the start and again before the item-3 build (new
module).  Only `src/rotation/**` and the regenerated `src/deps.mk` were
edited; `check_boundaries.py` needed no allowlist change (all new modules
are rotation-internal, and `use controls_lib, only: ichi_*` is already
allowed).  Nothing under `audit/` or `input/` touched; no push, rebase or
branch change.  18 files, +404/-261 over the four commits.

## Item 1 -- equal-grid toolkit, identical part only (`2e92a80`)

`diff -w` of `seculr/am_transport_grid.f90` against `seculr/composition_grid.f90`
(after renaming) gave four token-identical blocks, now the four public
procedures of the new rotation-internal module
`src/rotation/seculr/equal_grid_lib.f90` (149 lines), each taking the
`rot_scr%xtab/ytab/xval/yval` scratch as explicit dummies:

| Helper | Former inline block |
|--------|---------------------|
| `edge_grid_abscissae(chi, echi, dchi, ntab, ntot, xtab, xval, ntabb)` | `xtab(1)=chi(1)`, midpoints `0.5d0*(chi(i)+chi(i-1))`, `ntabb=ntab+1`, `xtab(ntabb)=chi(ntab)`, `xval(1)=chi(1)`, `xval(i)=echi(i)-0.5d0*dchi` |
| `interp_edge_coeff(coeff, zone_begin, zone_end, ntab, ntabb, ntot, xval, xtab, ytab, eq_coeff)` | `ytab(1)=coeff(zone_begin+1)`, `ytab(i)=coeff(zone_begin+i-1)`, `ytab(ntabb)=coeff(zone_end)`, `osplin(xval,eq_coeff,xtab,ytab,ntabb,ntot)` |
| `dchi_dr_jacobian(log_density, log_radius, log_mass, log_pressure, enclosed_mass, epsm, surface_luminosity_lsun, zone_begin, ntab, ntot, chi, xval, xtab, ytab, yval)` | the mass/luminosity/pressure scales from `chi_grid_scale(ichi_*)`, the `four_pi_rho_r2` / `dchi_dr` loop writing `ytab(i)`, `osplin(xval,yval,xtab,ytab,ntab,ntot)` |
| `multiply_by_exp10(coeff, log_factor, ntot)` | `coeff(i) = coeff(i)*exp(ln10*log_factor(i))` |

Every statement inside the helpers is the former text; the callers call
each helper at the position the block occupied with the same inputs
(`am_transport_grid.f90:174-183, 214-234`; `composition_grid.f90:105-119`).
`luminosity_lsun(num_zones)` is passed as the scalar
`surface_luminosity_lsun`; `mix_scr%epsm` is passed explicitly.  The only
textual reordering: in `am_transport_grid` the xtab/xval construction was
interleaved with the first `ytab` fill for `am_diffusion_coeff`; the two
are independent assignments to distinct arrays, so filling xtab/xval first
then ytab is the same sequence of stores (commit message records this).
The three `exp10` multiply loops in `am_transport_grid` (mixing coeff, AM
coeff, and the advective+diffusive pair under
`use_diffusion_advection_transport`) became four `multiply_by_exp10` calls
in the same order.

Not extracted (non-identical between the two files) -- see Deferred 1-2.

## Item 2 -- dead-code purge with byte parity (`4c9509d`)

Each candidate in the brief was checked against the current tree:

| Candidate | Status |
|-----------|--------|
| `mid_timestep_model.f90` 197-352 block (`new_cz_detected`) | already deleted by wave 1 (`a4fccb1`); no trace |
| `evolve_angular_momentum.f90` BUR-ST branches (`burs_extrapolation_active`, lines 64, 96, 321, 356, 391, 402) | left -- author decision per the brief (BS-extrapolation path) |
| `viscos.f90` Endal-Sofia arithmetic | already deleted by wave 1 |
| `am_advection_diffusion_coeffs.f90` `residual_check/coeff_matrix_orig/damping_factor` | already deleted by wave 1 |
| `am_advection_diffusion_coeffs.f90` `theta_term_n/p` (lines 199-204) | left: they are read (`+theta_term_n`, `+ theta_term_p`) inside the coefficient expressions; deleting them changes the expression text (`x + 0.0d0` -> `x` is not byte-safe for signed zeros).  Deferred 4 |
| `banded_solver.f90` (audit 114-152) | already deleted by wave 1 (file is 111 lines) |
| `diffuse_composition.f90` `dcomp/sum_species` | already deleted by wave 1 |
| `microdiff_mte.f90` debug print | already deleted by wave 1 |
| `zone_moments_of_inertia.f90` `prev_log_*/spherical_moment_of_inertia` | already deleted by wave 1 |
| `rotation_shape_factors.f90` `aintt` | **still present** (wave-1 report listed it as removed): declaration and `aintt = 0.0d0` deleted now |
| `circulation_velocities.f90` commented F77 | already deleted by wave 1 |

A write-only scan of every rotation local (scratchpad script: locals with
assignments but no read in their unit, checked by hand) found three more
survivors, deleted with their assignment statements:

- `microdiff/gravitational_settling.f90`: `metal_x_prev_iter(json)` (two
  assignments inside the `use_diffusion_z` block, never read);
- `microdiff/gravitational_settling_setup.f90`: `hydrogen_fraction_cubed`
  (assigned, never read);
- `seculr/am_advection_diffusion_coeffs.f90`: `omega_mid_start(json)`
  (assigned in the domega_dr loop, never read).

No flag-guarded block was removed.  `use_generic_diffusion_vectors`
(`am_transport_grid`) and `use_diffusion_advection_transport` are live
namelist flags and their blocks stay.

## Item 3 -- species-descriptor tables (`f8bced8`)

Value comparison of the lists first:

| Routine | H | He | metal / diffused | e- / other |
|---------|---|----|------------------|------------|
| `gravitational_settling_setup` DATA | 1.008 / 1 | 4.004 / 2 | 55.86 / 26 (Fe) | `m_electron_amu` / -1 |
| `microdiff_coefficients` assignments | 1.008 / 1 | 4.004 / 2 | argument | `m_electron_amu` / -1 |
| `microdiff` (Fe, light elements) | -- | -- | 55.86 / 26; Li6 6.015/3, Li7 7.016/3, Be9 9.012/4, rows 13/14/15 | -- |
| `viscos.f90` 11-row table | 1.007825 / 1 | 4.0026 / 2 | Z 1.0 / 0 | He3 3.01603/2, C12..O18 |
| `check_composition.f90` 4-row table | 1.007825 | 4.002603 | 12.0 | He3 3.01603 |

The first three agree value-for-value and now read one table; the last two
differ (different A for H and He, and `viscos` treats Z as a single
species of weight 1/charge 0) and stay local, each with a comment saying
so.

New rotation-internal module `src/rotation/species_table_lib.f90`
(58 lines): `type species_props (name, weight, charge)`; the column
numbers `thoul_col_h=1, thoul_col_he=2, thoul_col_metal=3,
thoul_col_electron=4`; parameters `thoul_h1, thoul_he4, thoul_fe,
thoul_electron`; `num_light_diffused=3`, `light_diffused(3)` and
`light_diffused_row = [i_li6, i_li7, i_be9]` (= the former
`data light_element_id/13,14,15/`).  Readers:

- `gravitational_settling_setup.f90`: the DATA statements take
  `thoul_*%weight` / `%charge` (constant subobjects; compiles cleanly);
  `species_mass_fraction(1..3)` and `settling_ap/at(1)`, `(3)` become
  `thoul_col_*`.
- `microdiff_coefficients.f90`: the eight assignments read `thoul_*`;
  `mass_frac(1..3)`, `species_fraction(1..3,i)`, `concen(1)`,
  `conc_coeff(species_col,1)` and `species_col.eq.1/.eq.3` become
  `thoul_col_*`.
- `microdiff.f90`: the three local DATA arrays are gone;
  `num_light = num_light_diffused`; `light_diffused_row` is passed to
  `microdiff_mte`/`microdiff_etm` (both `intent(in)`); Fe literals and the
  light-element loop read the table; `species_col` and the
  `species_fraction(_mid)(1..3, ...)` subscripts become `thoul_col_*`.

Bare composition-row subscripts became `i_*` names from `star_info_lib`
in `gravitational_settling_setup`, `microdiff_setup`, `microdiff_mte`,
`microdiff_etm`, `check_composition`, `diffuse_composition` and `viscos`
(the `species_idx.eq.3` skips become `i_metals`).  Mapping proof: the
`composition(15,json)` dummy of every one of these routines is
`star%xa` passed down from `mix.f90` -> `secular_transport` /
`gravitational_settling` / `microdiff`, and `shell_physics` indexes the
same array with `i_h1`/`i_metals`; rows 1..4 are `i_h1, i_he4, i_metals,
i_he3` and rows 5..11 are `i_c12..i_o18` (`microdiff_etm`'s three
`do j = 5,11` loops that rescale the CNO isotopes with Z became
`do j = i_c12,i_o18`, same bounds).

The per-species diffusion descriptor loop was not built (Deferred 5).

## Item 4 -- `am_advection_diffusion_coeffs` cleanups (no commit)

Everything in the brief's list was already done by sidechan W2
(`f9b3869`, `reports/sidechan-wave2.md` item 7): the `json/nmax`
shadowing had gone before that wave (`json` from `star_info_lib`,
`band_nmax` from `rotation_scratch_lib`); the four write-only locals were
deleted; `timestep` and `wind_loss_implicit` are `intent(in)` with
`substep_timestep` / `wind_loss_implicit_iter` locals -- sidechan did
those as byte-safe, so nothing was left to "leave".  A fresh write-only
scan of the routine after `4c9509d` (which removed the one survivor,
`omega_mid_start`) finds nothing.  Recorded in the `1830868` commit
message.

## Item 5 -- `solve_r0` Newton loops (`1830868`)

Comment-only.  `shape/shape.f90` lines 53-55 (center loop, exits on
`.le. acfpft`) and 100-102 (per-zone loop, exits on `.lt. acfpft`) each
carry a one-line note pointing at the other and the operator difference,
both kept as in the original.

## Deferred (exact proposals)

1. **Zone-mass construction, `am_transport_grid.f90:107-170` vs
   `composition_grid.f90:57-100`.**  Same mass arithmetic, but
   `am_transport_grid` interleaves the moment of inertia and angular
   momentum with the mass in the same loops (`eq_moment_of_inertia(i) =
   eq_reduced_moment_of_inertia(i)*eq_mass(i)*...`, `eq_angular_momentum(i)
   = eq_angular_momentum(i)*eq_mass(i)` at 113-115, 132-134, 156-158, and
   the CZ additions at 139-142, 164-166), and `eq_mass(1)` is read for I
   and J before the convective shells are added to it.  A shared
   `edge_zone_masses` helper would have to be called before the I/J
   products and the CZ sums split into two loops, i.e. a reorder of
   dependent stores -- not token-identical, so not extracted.  Proposal
   (R6): share only the interior loop `eq_mass(i) =
   0.5d0*(eq_enclosed_mass(i+1)-eq_enclosed_mass(i-1))` plus the
   `em_top/em_bot` end-point values, returning them for the caller's I/J
   products; still needs a pin check because the CZ sums would move.
2. **Scaled edge interpolations, `am_transport_grid.f90:185-200`.**  The
   two `es_advective_velocity` / `es_diffusive_velocity` tables under
   `use_diffusion_advection_transport` are `interp_edge_coeff` with
   `scale_factor*` on every table entry.  Extending the helper with an
   optional scale would make the unscaled path's text `1.0d0*coeff(...)`
   or a branch -- either changes the expression; left inline.  Proposal:
   a second helper `interp_edge_coeff_scaled` with the `scale_factor*`
   text (token-identical to the two blocks); trivial, but adds a third
   copy of the fill loop, so left for the author.
3. **CZ-extent gating (R6, per the brief).**  `am_transport_grid.f90:124,
   136` and `composition_grid.f90:69, 77` gate on `zone_begin.gt.1` /
   `zone_end.lt.num_zones` with accumulate loops that `exit` on the first
   non-convective shell; `diffuse_composition.f90:136`
   (`convective_flag(zone_begin).and.zone_begin.gt.1`, `i1` fix in the
   else branch) and `equal_grid_to_model.f90:68`
   (`.not.convective_flag(zone_begin) .or. zone_begin.eq.1`) are the
   logical complement of each other and use a search that stops at the
   first convective shell.  Not touched.
4. **`theta_term_n/p`, `am_advection_diffusion_coeffs.f90:199-204`.**
   Both are `0.0d0` on every iteration and added into `diffusive_term2`
   and `advective_term1`.  Deleting them turns `x + 0.0d0` into `x`,
   which differs only in the sign of a zero result (`-0.0d0 + 0.0d0 =
   +0.0d0`); gfortran at `-O3` without `-ffast-math` keeps the addition,
   and the coefficients can be exactly zero in convective zones, so this
   is not provably byte-safe by inspection.  Proposal: delete the two
   locals and the `+theta_term_*` text, run the rotation pins; expected
   identical but must be pinned.
5. **Per-species diffusion descriptor loop (class B, per the brief).**
   `microdiff.f90:146-229` has three near-identical blocks (H, metals,
   light elements) that set `species_col`, the diffused (A, Z), fill
   `species_fraction(thoul_col_metal, :)` and call `microdiff_run`.  A
   loop over a descriptor list would reorder the per-species blocks'
   restore logic (the H block restores `species_fraction(thoul_col_h,:)`
   after its run, the metal block does not) and is not token-identical.
   `species_table_lib` now holds the descriptors it would need.
6. **`viscos.f90:28-34` and `check_composition.f90:50-51` tables.**  Differ
   value-for-value from the Thoul set (see the item-3 table); a shared
   table would have to carry two weight columns.  Author decision whether
   the 1.007825/4.0026 vs 1.008/4.004 difference is intended.
7. **`check_composition.f90:67-75`.**  `num_diffused_species` is `4`
   (or `11` on the last diffusion iteration) and the loop skips
   `i_metals`, so it walks rows `i_h1, i_he4, i_he3` (plus
   `i_c12..i_o18`) -- correct, but it depends on `i_he3 = 4` and
   `i_o18 = 11` being the last rows of each set; explicit row lists would
   make it self-describing.  Not done (changes loop structure).

## Reverted (changed numbers)

None.  No pin differed at any step.

## Skipped

- BS-extrapolation path (`evolve_angular_momentum.f90` 64-402
  `burs_extrapolation_active`, `rotation_stability_setup.f90` ~458-480,
  `difad_shear_coeff2` consumers): author decision per the brief.
- Item 4: nothing left to edit (see above).
- `equal_grid_lib` was not registered in `check_boundaries.py`: the
  script passed unchanged ("every cross-domain call goes through a public
  entry"), so no allowlist edit was needed.

## What the audit / earlier reports got wrong

- `reports/rotation-wave1.md` R1 lists `aintt` (`rotation_shape_factors`)
  among the write-only locals deleted; it was still declared and assigned
  (`aintt = 0.0d0`) at base `6449462`.  Removed in `4c9509d`.
- Audit Structural 1 describes the two grid files as sharing "the
  CZ-extent searches"; they share the accumulate-style search, but
  `diffuse_composition` / `equal_grid_to_model` use a different one
  (Deferred 3), and the zone-mass block differs by the I/J interleaving
  (Deferred 1), so the identical part is the four helpers above, not the
  whole grid construction.
- Audit items 6/20 say the three setup routines "each carry their own
  species list"; `microdiff_coefficients` takes the diffused species as an
  argument and `microdiff` carries the Fe and light-element descriptors,
  so the shared set is H/He/e- plus Fe and Li6/Li7/Be9, and `viscos` is
  not a member of the value-identical group (item-3 table).
- `reports/sidechan-wave2.md` item 7 covered every part of brief item 4,
  including the two the brief expected to be deferred; nothing was left
  as "rotation-only, safe".

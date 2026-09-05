# Readability wave 2 -- rotmix

Worktree `/Applications/YREC-wt/rotmix`, branch `rs/rotmix`, base `c764c46`.

## Commits

| hash | summary |
|---|---|
| `796e1ae` | Readability W2 (rotmix): luminosity_lsun renames, named reaction_rate_by_zone rows (items 1, 2) |
| `7fedff4` | Readability W2 (rotmix): rot_scr legacy names, one lagrange4 for microdiff (items 3, 4) |

Verification: gate1 IDENTICAL before each commit; on the final tree
(`7fedff4`) the full 37-case pin selection ended `PINS EXIT 0` (37 passed,
`logs/rotmix.pins.W2.log`) and the aux battery `AUX EXIT 0` (22 passed,
`logs/rotmix.aux.W2.log`). A first pin run started at 11:55 died with the
machine (empty log) and was re-run from scratch at 13:04. No `test_*.short`
files left behind, no orphaned yrec processes. Working tree clean; only
source files and `src/deps.mk` committed.

## Item 1 -- `luminosity_lsun*` renames

The brief named `log_luminosity`/`log_luminosity_start`; the actual state of
the tree was different: `star%luminosity_lsun` / `star%luminosity_lsun_start`
had already been renamed by earlier waves, and the surviving misnomer was the
scratch member `rot_scr%log_luminosity_mid` plus the `log_luminosity` dummies
that receive it. Renamed:

- `rot_scr%log_luminosity_mid` -> `rot_scr%luminosity_lsun_mid`
  (`rotation/rotation_scratch_lib.f90`, with a comment that it is linear
  L/Lsun unlike its `log_*` neighbours); users in
  `rotation/mid_timestep_model.f90`, `rotation/evolve_angular_momentum.f90`,
  `rotation/seculr/secular_transport.f90` (5 sites).
- dummy `log_luminosity` -> `luminosity_lsun` in
  `rotation/equal_spaced_grid.f90`, `rotation/seculr/am_transport_grid.f90`,
  `rotation/seculr/composition_grid.f90`,
  `rotation/seculr/diffuse_composition_driver.f90`,
  `mixing/burn_settle_mix.f90`.
- `mixing/semiconvection.f90` and `mixing/temperature_gradients.f90` were
  already correct (mixwind R2); nothing to do.

How linear vs. true-log were told apart (stated in the commit body): by data
flow, not by name. `luminosity_lsun_mid` is filled in `mid_timestep_model`
from `star%luminosity_lsun_start` + fraction * (`star%luminosity_lsun` -
start), and `secular_transport` multiplies it by `star%solar_luminosity_cgs`
to get erg/s -- a log would be exp10'd. `equal_spaced_grid` uses it as
`L/(LTOT*DL)` next to `LOG(M)/DM` and `LOG(P)/DP`. Conversely the wind
scalars `log_luminosity_lsun` (`wind/*.f90`, `core/rebuild_envelope.f90:40`)
are exp10'd before use and `star%trial_log_luminosity` feeds the Henyey
log variables; those are true logs and were left alone.

## Item 2 -- named rows of `reaction_rate_by_zone(15,json)`

Writer verified: `mixing/mix.f90:230-244` stores the 15 outputs of
`net_lib`'s `rates` (`net/net_lib.f90:1669`, originally HR1..HR13, HF1, HF2)
into rows 1..15 in dummy-argument order. Added to `rotation_scratch_lib`:

```
rr_pp=1, rr_he3_he3=2, rr_he3_he4=3, rr_c12_p=4, rr_c13_p=5, rr_n14_p=6,
rr_o16_p=7, rr_c13_alpha=8, rr_zero9=9, rr_c12_alpha=10, rr_n14_alpha=11,
rr_triple_alpha=12, rr_zero13=13, rr_frac_c12_alpha=14, rr_frac_be7_electron=15
```

Used in `mixing/rotmix.f90` (12 literal-row reads) and in
`setup/rezone.f90` (`use rotation_scratch_lib, only: rr_*`;
`reaction_rate_species_index(7) = [rr_pp, rr_he3_he3, rr_c12_p, rr_c13_p,
rr_n14_p, rr_o16_p, rr_zero9]`, plus the loop comment). Only the item-2 lines
of `rezone.f90` were touched (shared with ioctrl).

Finding: the mixwind wave-1 proposal described these rows as species
(H1/He4/He3/C12/C13/N14/O16). They are reaction rates. The original HCOMPM
regrid loop in `rezone` copied rows 1,2,4,5,6,7,9 as if they were the
composition slots, so rows 3 (He3+He4), 8, 10-15 are never regridded. That
is a historical behaviour, kept byte-for-byte (fixing it is class B); the
comment at the parameter now says so.

## Item 3 -- `rot_scr` legacy member names

Rule applied: rename or delete only where every reader and writer is in my
file list (rotation/, mixing/), so the change is self-contained and the
meaning is provable from the code; everything else is deferred with a
proposal.

Done:
- `es1(json)` -> `eq_enclosed_mass(json)`: set by `equal_spaced_grid`
  (osplin of log M on the equal-chi grid, then exp10 -> grams), differenced
  by `am_transport_grid` and `composition_grid`. Comment added.
- `wmst(json)` -> `omega_substep_start(json)`: written in
  `evolve_angular_momentum` from `omega` at the start of each sub-step, read
  by `secular_transport` as `omega_surface_start`. Comment added.
- `vesd2(json)`, `vesd3(json)` deleted: written in
  `seculr/circulation_velocities.f90` (`0.2d0*facd2*omega_interface`, then
  `*qqq` rescale) and never read anywhere in `src/` (whole-tree grep).
- `qwrmst(json)` deleted: two assignments in `evolve_angular_momentum`,
  no reader anywhere. The local `domega_dr` is kept because `theta_prev`
  still uses it.
- `facd2/facd3`: kept (written by `setup/rotation_stability_setup.f90`,
  outside my files) but now provably unread; comment says so. See Deferred.
- `qwrst`: kept (written by `setup/rezone.f90:832`, an ioctrl line) but now
  unread; comment says so. See Deferred.

Not renamed (writer in `setup/rotation_stability_setup.f90`, not mine):
`hs3`, `pm`, `tm`, `dm`, `delami`, `delmi`, `epsilm`, `qdtmi`, `fact1-6`,
`fgsfj`, `tho`. Proposals in Deferred. `masschg2/3` -> `star%`: deferred
with proposal.

## Item 4 -- one `lagrange4`

`numerics_lib` already had `pure function lagrange4(w, y)` with
`s = w(1)*y(1)+w(2)*y(2)+w(3)*y(3)+w(4)*y(4)` -- token-identical, same
order, to the microdiff copy `lagrange4(fac, a, k0)` which summed
`fac(1)*a(k0)+...+fac(4)*a(k0+3)`. Deleted the copy from
`rotation/microdiff/microdiff_mte.f90` (module `microdiff_mte_lib`); the 10
call sites there and 3 in `microdiff_etm.f90` now call
`lagrange4(facinterp, X(k0:k0+3))`. The 8 inline
`facinterp(1)*composition(row,k0)+...+facinterp(4)*composition(row,k0+3)`
sums in `microdiff_mte.f90` were token-identical in order and were folded
into `lagrange4(facinterp, composition(row,k0:k0+3))`. `microdiff_etm.f90`
no longer uses `microdiff_mte_lib`; `src/deps.mk` regenerated (one line).
Pins byte-identical.

## Item 5 -- literal 12/13 zone-bound sites (not changed; state owns the
constants)

Sites that should become `max_convective_zones` (12) /
`max_radiative_zones` (13) once state adds them:

- `rotation/am_convective_regions.f90:29-30, 34`
- `rotation/mid_timestep_model.f90:51-52`
- `rotation/evolve_angular_momentum.f90:57-58`
- `mixing/rotmix.f90:35-36`
- `mixing/burn_settle_mix.f90:52-53`
- `mixing/semiconvection.f90:57`
- `mixing/overshoot_boundaries.f90:26-27`
- `mixing/find_convection_zones.f90:63-65`
- `mixing/mix.f90:48, 62`
- `setup/rezone.f90:48`

(`mid_timestep_model.f90:156` `xa_start(12,...)` and
`evolve_angular_momentum.f90:378-393` `12,15` are species indices, not zone
bounds; `mix.f90:134` `accreted_composition(12)` likewise.)

## Deferred (exact proposals)

1. `mixing/mix.f90:230-244` (not in my list): the writer should use the
   `rr_*` names, e.g. `rot_scr%reaction_rate_by_zone(rr_pp,zone_idx) =
   rate_pp(zone_idx)` ... `(rr_frac_be7_electron,zone_idx) =
   frac_be7_electron(zone_idx)`.
2. `core/shell_physics.f90:14,31` and `core/rebuild_envelope.f90:12,30`:
   dummy `log_luminosity(json)` receives `star%luminosity_lsun` (linear;
   `shell_physics.f90:82` copies it straight into `luminosity_lsun`).
   Rename to `luminosity_lsun` (core domain). `rebuild_envelope.f90:40`
   `log_luminosity_lsun` is a true log (exp10'd at :109) -- keep.
3. `setup/rotation_stability_setup.f90` (setup domain) with matching
   `rotation_scratch_lib` members:
   - `fact1..fact4`, `fact6`, `fgsfj` are write-only (lines 97-102, 164,
     170, 175-176, 183, 196; no reader in `src/`): delete the assignments
     and the members. There is no `fact5` in the tree (the brief's
     `fact1-6` range is off by one).
   - `facd2`, `facd3` (lines 245-246): now unread after the `vesd2/3`
     deletion -- delete the assignments and the members.
   - interp4 outputs at lines 311-327: `hs3` -> `interface_mass`,
     `pm` -> `interface_pressure`, `tm` -> `interface_temperature`,
     `dm` -> `interface_density` (also read at
     `rotation/seculr/zahn_coupling_factor.f90:80,87,95` as
     `dm*interface_radius`), `delami` -> `interface_grad_ad`,
     `delmi` -> `interface_grad_rad`, `epsilm` -> `interface_eps`,
     `qdtmi` -> `interface_qdt`.
4. `setup/rezone.f90:831-832` (outside my item-2 hunk): `tho` -> e.g.
   `theta_rezoned` (read at `evolve_angular_momentum.f90:231,234` into
   `theta_prev`); `qwrst` is now unread -- delete line 832 and the member.
5. `masschg2/3` -> `star%`: members `accretion_specific_entropy`,
   `envelope_specific_entropy`, `updated_mass_msun`, `delta_log_pressure`,
   `delta_log_temperature` (used at `core/henyey_coefficients.f90:266,268`,
   `wind/massloss.f90:172-259`, `wind/mdot.f90:160,241`) and
   `solar_wind_mass_loss_rate_msun_yr`, `wind_reference_omega`,
   `wind_max_omega`, `use_rotation_scaled_solar_wind` (`wind/mdot.f90:86-89`).
   Proposal: move the first group to `star%` as
   `star%accretion_specific_entropy` etc. (they are per-model wind/accretion
   state, not rotation scratch) and the second group to `star%ctrl%` (they
   are namelist controls). Author decision needed first:
   `rot_scr%use_rotation_scaled_solar_wind` is never assigned anywhere in
   `src/` -- the inlist reads `lsolwind` into a `read_controls` local and
   `inlist_new_read.inc:265` sets a *local* `use_rotation_scaled_solar_wind`,
   neither of which reaches `rot_scr`; the three `wind_*` reals are likewise
   never assigned. So the `mdot.f90:86` branch is dead code today. Either
   wire the control through (class B, changes numbers for inlists with
   `lsolwind=.true.`) or delete the branch and the four members.
6. `rotation/seculr/composition_grid.f90:129-132` carries literal
   `chi_grid_scale(2)/(9)/(11)` subscripts that wave-1 deferred #1 did not
   list (it listed only `equal_spaced_grid.f90:40` and
   `am_transport_grid.f90:231`). Same treatment as those.
7. Merge note for the integrator: `equal_spaced_grid.f90:40` and
   `am_transport_grid.f90:231` (ioctrl's `chi_grid_scale(9)` lines) now also
   contain my renamed `luminosity_lsun(num_zones)`. If ioctrl edited the
   subscript on the same line, the merge conflicts textually; the resolution
   is `chi_grid_scale(<ioctrl name>)*luminosity_lsun(num_zones)`.

## Reverted (changed numbers)

None. Every pin stayed byte-identical.

## Skipped

- Namelist-facing names (`lsolwind`, `chi_grid_scale`, etc.): out of scope
  by the rules.
- Fixing the `rezone` regrid row selection (rows 3, 8, 10-15 not
  regridded): class B.
- `fact*`, `fgsfj`, `facd*`, `hs3`/`pm`/`tm`/`dm`/... renames: writers live
  in `setup/rotation_stability_setup.f90`, outside my file list (see
  Deferred 3).

## What the audit / wave-1 reports got wrong

- mixwind wave-1 #4: `reaction_rate_by_zone` rows are reaction rates
  (net_lib `rates` outputs), not species; the proposed
  `H1/He4/He3/C12/C13/N14/O16` labels for rows 1,2,4,5,6,7,9 were wrong and
  the old `rezone` comment saying the same was misleading (now replaced).
- The brief's item 1 names (`log_luminosity_start`) no longer existed;
  the live misnomers were `rot_scr%log_luminosity_mid` and the
  `log_luminosity` dummies.
- The brief said ioctrl edits `atime(` subscripts in `setup/rezone.f90`;
  there is no `atime` in that file. ioctrl's rezone hunks are the
  `chi_grid_scale(...)` subscripts (lines 202-367 and 522), none of which
  overlap my item-2 hunks (lines 18-25, 44-59, 604-619 in the new
  numbering).
- The stale HCOMPM comment in `rotation/mid_timestep_model.f90` described a
  composition reset that the code does (`star%xa` <- `star%xa_start`), not
  reaction rates; rewritten.

# Readability wave 2 -- state (star_info_lib / phys_const_lib and their consumers)

Worktree `/Applications/YREC-wt/state`, branch `rs/state`.

## Commits

| hash | summary |
|---|---|
| `3c5f090` | reaction_rate array, named slot counts, physical constants into phys_const_lib (items 1-6 declarations, item 2, 3, 4) |
| `18bd280` | use the named counts; ln10_over_lsun, chosen_dt; drop massloss age_gyr; ix_* mixture indices (items 5-10 consumers) |

`src/deps.mk` did not change (`gen_makefile_deps.py --check`: "deps.mk up to
date" -- no file gained or lost a `use` of a module it did not already use).

## Verification

- gate1 (tier 1): `GATE1: IDENTICAL` before each commit and on the final tree.
- Full 37-case pin selection: `37 passed, 272 deselected` -- `PINS EXIT 0` (byte-identical) on 18bd280. (A first run was killed by the harness while the machine slept, before producing any output; it was rerun from scratch.)
- Aux battery: `22 passed` -- `AUX EXIT 0` on 18bd280. No `test_*.short` files left behind.
- Only build warning in touched files is the pre-existing
  `core/burn_lib.f90:2272 light_element_save(3,json)` -Wsurprising.

## Per item

1. **reaction_rate_1..13 -> reaction_rate(n_reactions,json)** -- done.
   `state/star_info_lib.f90`: thirteen `reaction_rate_N(json)` members replaced by
   `reaction_rate(n_reactions,json)`, `n_reactions = 13` declared next to it
   (same value as net_lib's `num_reactions`; a copy because net_lib uses
   star_info_lib, so state/ cannot import it -- said so in the comment).
   `core/burn_lib.f90` (compute_neutrino_emission): the thirteen stores become
   `star%reaction_rate(r_pp,shell_index) = reaction_rate(r_pp)*gyr_amu_per_sec_gram`
   ... `star%reaction_rate(r_c12c12_unused,shell_index) = ...`, same order, same
   right-hand sides; the two header comments updated. grep of `src/` for
   `reaction_rate_[0-9]` finds no reader in io/ (or anywhere) -- no deferral needed.

2. **i_nu_he3he3 = 9, i_nu_he3he4 = 10** -- done.
   `star_info_lib.f90`: added to the i_nu_* block together with `n_nu_fluxes = 10`;
   the "spare slots" comment now says slots 9-10 hold the He3+He3 and He3+He4
   reaction rates that compute_neutrino_emission stores in flux units.
   `burn_lib.f90`: `star%neutrino_flux(9)`/`(10)` -> `(i_nu_he3he3)`/`(i_nu_he3he4)`;
   engeb's only-list extended.

3. **gyr_amu_per_sec_gram** -- done. `phys_const_lib.f90`: `double precision,
   parameter :: gyr_amu_per_sec_gram = 5.240358d-8`. `net/net_lib.f90`: the
   parameter deleted (a two-line comment points at phys_const_lib); rates,
   deutrate and engeb already `use phys_const_lib`.

4. **Rotation literal constants** -- done, literal lines only (hunk-disjoint with
   rotmix). `phys_const_lib.f90` gains `m_proton_cgs = 1.6726d-24`,
   `amu_cgs_legacy = 1.6605655d-24`, `m_electron_amu = 5.486d-4`,
   `rsun_cgs_legacy = 6.9598d10`, `sigma_sb_cgs_legacy = 5.669d-5`,
   `rad_viscosity_coeff_cgs = 6.7282653d-26`, all double precision, each with a
   comment. The two "identify" constants: 5.669d-5 is the older Stefan-Boltzmann
   value (differs from `csig = 5.67051d-5`, hence a separate legacy parameter);
   6.7282653d-26 = 4a/(15c) in cgs, the coefficient of the Ledoux (1958) radiative
   viscosity nu_rad = 4aT^4/(15 c kappa rho^2). `rsun_cgs_legacy` equals the
   namelist default at io/read_controls.f90:1131 but that is a variable, so a
   parameter is still needed.
   Sites: `rotation/viscos.f90` (`data amu/amu_cgs_legacy/`,
   `viscosity_radiative = rad_viscosity_coeff_cgs*...`,
   `1.6d1*cc13*sigma_sb_cgs_legacy*...`);
   `rotation/microdiff/gravitational_settling_setup.f90` (`data atomic_weight/...,m_electron_amu/`,
   `solar_radius_bl = rsun_cgs_legacy`, `rho_local/(m_proton_cgs*mass_weighted_conc)`);
   `rotation/microdiff/microdiff_coefficients.f90` (`atomic_weight(4) = m_electron_amu`,
   `ne = rho/(m_proton_cgs*ac)`); `rotation/microdiff/microdiff_setup.f90`
   (`crsun_bah = rsun_cgs_legacy`). All four already `use phys_const_lib`.

5. **Named counts** -- done. `star_info_lib.f90`: `n_species_basic = 11`,
   `n_species_extended = 15`, `n_lum_channels = 8`, `n_nu_fluxes = 10`,
   `max_mc_runs = 1000`; `xa`, `xa_start`, `orig_composition`,
   `luminosity_breakdown`, `neutrino_flux`, `neutrino_flux_total`,
   `neutrino_flux_zone` and the nine Monte-Carlo sample arrays declared with them.
   Consumers: `core/henyey_iterate.f90` (num_species, `do j = 1,n_lum_channels`),
   `core/evolve_step.f90` (num_species), `core/read_starting_model.f90`
   (num_species, `reference_composition(n_species_extended)`),
   `core/rebuild_envelope.f90` (`composition(n_species_extended,json)`,
   `species_end_index`), `core/henyey_coefficients.f90` (`n_lum_channels`,
   two `n_nu_fluxes` loops, `luminosity_breakdown(8)` -> `(i_lum_he_c)`),
   `core/observables_lib.f90` (renormalize_luminosity_breakdown loop),
   `core/neutrino_flux_table.f90` (three loops and the two implied-do write
   lists -- record layout unchanged), `core/monte_carlo.f90`
   (`min(star%job%mc_run_end,max_mc_runs)`).
   The three `do k = 4,species_end_index` in rebuild_envelope are now
   `do k = i_he3,species_end_index`: 4 is provably i_he3 -- the routine sets
   slots 1 (H) and 2 (Z) explicitly and computes He4 (slot 3) from slots 1, 3, 4
   after the copy, so the loop copies exactly the species from He3 up; stated in
   the commit message.

6. **max_convective_zones = 12, max_radiative_zones = 13** -- done for my files.
   `star_info_lib.f90`: declared; `mixed_zone_bounds`,
   `mixed_zone_bounds_no_overshoot` (both (max_convective_zones,2)) and
   `radiative_zone_bounds(max_radiative_zones,2)`.
   `mixing/find_convection_zones.f90`: three dummy declarations and
   `if (j_idx.lt.max_convective_zones) cycle`. `mixing/mix.f90`: dummy
   `mixed_zone_bounds_no_overshoot(max_convective_zones,2)` and local
   `radiative_zone_bounds(max_radiative_zones,2)`. `core/henyey_iterate.f90`
   and `core/read_starting_model.f90`: locals.
   The mix.f90:48 dummy was NOT changed to (13,2): both actuals are (12,2)
   (`star%mixed_zone_bounds_no_overshoot` from evolve_step, henyey_iterate's
   local), and the star_info_lib member is (12,2), so the audit's claim of a
   (13,2)/(12,2) mismatch is wrong; only the dimension is named.
   `setup/locate_shell_boundaries.f90` contains no 12/13 literal at all -- the
   (13,2)/(12,2) locals the brief had in mind are in `setup/rezone.f90:46`
   (not in my file list; see Deferred).

7. **ln_solar_luminosity -> ln10_over_lsun** -- done. Member renamed with a
   comment saying it is ln(10)/solar_luminosity_cgs, not a logarithm of Lsun;
   `setup/setups.f90`, `core/monte_carlo.f90`, `core/henyey_coefficients.f90`
   (`cccql = star%ln10_over_lsun*star%m(im)`). No io/ user; grep for the old
   name over src/ is empty.

8. **hydrogen_dt -> chosen_dt** -- done. Member renamed with the comment
   "chosen_dt (former DELTSH): the timestep compute_timestep adopted (the
   minimum over all its limiters, seconds); read_starting_model reuses the slot
   transiently as |delta_time|. Nothing reads it as an H-burning governor."
   Sites: `core/evolve_step.f90` (compute_timestep and read_starting_model
   calls), `core/run_yrec.f90` (two). `util/compute_timestep.f90`: comment only
   ("the caller stores it in star%chosen_dt") -- its local `hydrogen_dt` and the
   dummy of `util/timestep_limit_hburn.f90` are the real H-burn limiter and stay.
   Verified by reading compute_timestep: the actual passed for the member is the
   final `min` of all limiters, never fed back as an H-burn quantity.
   `core/main.f90:82-86`: main.f90 is 23 lines; the note was already removed by
   core R1. Nothing to rewrite.

9. **massloss age_gyr** -- done. `wind/massloss.f90`: `age_gyr` removed from the
   dummy list and its declaration; grep of the file shows the name occurred only
   on those two lines (no contained procedures), so it was unread.
   `core/evolve_step.f90`: the call drops `star%dage`; the remaining arguments are
   unchanged and in the same order.

10. **ix_* mixture index block** -- done. `star_info_lib.f90`:
    `ix_na = 1, ix_al = 2, ix_mg = 3, ix_fe = 4, ix_si = 5, ix_c = 6, ix_h = 7,
    ix_o = 8, ix_n = 9, ix_ar = 10, ix_ne = 11, ix_he = 12, n_mix_species = 12`;
    `mixture_weights` and `fxenv` declared `(n_mix_species)`. The layout was
    proven from read_starting_model's atomic-weight DATA
    (23.0, 26.99, 24.32, 55.86, 28.1, 12.015, 1.008, 16.0, 14.01, 39.96, 20.19,
    4.004) and the LAOL table header "CS mix: C,N,O,Ne,Na,Mg,Al,Si,Ar,Fe".
    `core/read_starting_model.f90`: the mixture block uses the names (sum of ten
    metals, zenvm minus ix_c/ix_o/ix_n, `species_mix_weights(ix_h)`, `(ix_he)`),
    `species_mix_weights`/`atomic_weight(n_mix_species)`, three `do i = 1,12` ->
    `n_mix_species`. `kap/laol89/rdlaol.f90:50-51` only: read list
    `laol_work_array(ix_c), (ix_n), (ix_o), (ix_ne), (ix_na), (ix_mg), (ix_al),
    (ix_si), (ix_ar), (ix_fe)` -- same positions as 6,9,8,11,1,3,2,5,10,4; the
    only-list gains the ten names.

## Deferred (outside my file list or the item's site list)

- Zone-bound dimension literals that should become
  `max_convective_zones`/`max_radiative_zones` (all `(12,2)`/`(13,2)`):
  `mixing/rotmix.f90:35-36`, `mixing/burn_settle_mix.f90:52-53`,
  `mixing/semiconvection.f90:57`, `mixing/overshoot_boundaries.f90:26-27`,
  `setup/rezone.f90:46`, `rotation/evolve_angular_momentum.f90:57-58`,
  `rotation/am_convective_regions.f90:29-30,34`,
  `rotation/mid_timestep_model.f90:51-52`. Mechanical, value-preserving; the
  names are public in star_info_lib.
- `num_species = 15`/`11` literals in `mixing/mix.f90:92,94`,
  `mixing/burn_settle_mix.f90:77,79`, `mixing/homogenize_convection_zones.f90:63-64`,
  `mixing/rotmix.f90:66,68` -> `n_species_extended`/`n_species_basic`. mix.f90 is
  in my list but these lines are not in item 5's site list and mixing is the
  rotmix agent's rename territory, so left alone.
- `core/observables_lib.f90:373 do i = 1, 8` runs over the eight real neutrino
  sources of clsnuf_diag/gasnuf_diag (not the 10 flux slots and not the 8
  luminosity channels); it needs its own name (e.g. `n_nu_sources = 8`), which
  the brief did not ask for.
- Seconds-per-year duplicates `3.1558d7`: `core/henyey_coefficients.f90:116`
  (`one_year_sec`), `setup/setups.f90:57` (`seconds_per_year`),
  `rotation/microdiff/gravitational_settling_setup.f90:72` (`seconds_per_year_bl`),
  `rotation/microdiff/microdiff_setup.f90:66` (`csecyr_bah`) -- one
  `phys_const_lib` parameter would replace all four; not in the brief.

## Reverted (changed numbers)

None.

## Skipped

None of the ten items was skipped.

## Corrections to the audit / wave-1 reports

- mix.f90:48: the dummy `mixed_zone_bounds_no_overshoot(12,2)` matches both
  actuals and the star_info_lib member (all (12,2)); there is no (13,2) mismatch.
- `setup/locate_shell_boundaries.f90` has no 12/13 literals; the zone-bound
  locals are in `setup/rezone.f90:46`.
- `core/main.f90:82-86` (the DELTSH/hydrogen_dt note) no longer exists after
  core R1; main.f90 is 23 lines.
- `ln_solar_luminosity` and `reaction_rate_1..13` have no io/ readers, so the
  conditional deferrals in items 1 and 7 did not apply.

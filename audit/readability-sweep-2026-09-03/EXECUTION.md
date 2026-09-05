# Readability sweep -- execution record (R1-R5, 2026-09-03 .. 2026-09-04)

Batches R1-R5 of SUMMARY.md section 3 are done on `yrec-modern`
(7ef705c -> eece94b1, 75 commits, 225 files, +8309/-13186 lines). R6
(the number-changing batch) is not started; see "Open for R6" below.

Every commit is class A: the 37-case pinned selection (run_config_matrix,
run_standard_solar_model incl. the MESA baselines, testsuite Test_solar,
the three m0030 run_from_* cases) is byte-identical after every wave, the
22-test aux battery passes, and the ZAHB->TAHB spot check reaches the same
final model (1108). No pin was ever reseeded. Per-commit gate: clean
build with USE_HDF5=1, unit tests, `tools/check_boundaries.py`,
`tools/gen_makefile_deps.py --check`, `tools/gen_controls_state.py` drift
check, one solar byte-pin (`GATE1: IDENTICAL`).

## How it was run

Three waves of parallel agents, each in its own git worktree
(`rs/<name>` branch off `yrec-modern`), with a shared rulebook and one
brief per agent (both copied to `reports/wave*_*.md`). Agents edited only
the files their brief listed; files two agents both needed were declared
"shared, hunk-disjoint" and merged by rebase (the only conflicts were
overlapping renames and the generated `deps.mk`). Each agent's own
report is in `reports/<name>-wave<N>.md`: per item what was done, exact
proposals for what was deferred, and corrections to the audit.

| Wave | Batches | Agents (branch) | Merged at |
|---|---|---|---|
| 1 | R1-R3 per domain | burn, mixwind, rotation, core, kapatm, io, eos | c764c46 |
| 2a | cross-domain deferred | state, rotmix, ioctrl | 2ff54e0 |
| 2b | R4 explicit data flow | henyey, sidechan | 6449462 |
| 3 | R5 derived-type de-dup | core, rotation, mixburn, kap, eos | eece94b1 |

## What changed, by batch

**R1 (comments, dead code).** Stale headers naming COMMON blocks, F77
files and `SAVE` semantics rewritten or deleted across every domain;
commented-out F77 removed; unused locals/imports/dummies deleted after
grep showed them unread; dead branches behind constant conditions
removed (never behind a namelist flag). `wcz.f90` and other zero-caller
routines deleted. `ll95tbl` silent EOF became an error path.

**R2 (named indices and constants).** `star%reaction_rate_1..13` ->
`reaction_rate(13,json)`; named rows for `rot_scr%reaction_rate_by_zone`
(`rr_*`); species slots `i_h1..i_be9` used instead of bare `1..15` at
~300 sites; `n_species_basic/extended`, `n_lum_channels`,
`n_nu_fluxes`, `max_convective_zones/max_radiative_zones`,
`max_mc_runs`, `max_runs`; `chi_grid_scale`/`atime` slot names
(`ichi_*`, `itime_*`); kap table dimensions and missing-value sentinels
per family; OPAL/MHD/SCV output-slot names; physical constants from the
rotation files into `phys_const_lib`; `alex95` -> `alex94`; opacity unit
numbers named or replaced by `newunit=`. Namelist variable names were
not changed anywhere.

**R3 (small shared helpers, token-identical bodies).** `lagrange4`,
`kspline`/`splinc` as wrappers over `cspline`, `splintd2` over `splint`;
`blend_in_z`, `stencil` search helpers in kap; `locate_cz_base` /
`average_cz_abundances` for liburn/liburn2; wind radius/centrifugal
helpers; MHD ladder and SCV `spline_cell` helpers; `luminosity_breakdown`
renormalisation shared by observables_lib/neutrino_flux_table.

**R4 (explicit data flow).** The four henyey flag relays deleted;
`kenv`/`katm` removed from six signatures; `atm_get` takes an optional
`envint_step_config` and optional vertex outputs (the
save/override/restore of `star%job%*_step_*` and the dummy arrays are
gone); `temperature_gradients(zone_index)` replaces `star%iovim`;
rezone keeps its spline scratch in locals; microdiff receives the DEL
gradient as an argument instead of rotmix swapping `star%gradT`;
`mdot`/`compute_timestep` report their decisions and the drivers flip the
`star%job%` flags; `rdlaol` returns `debye_huckel_z` to `setupopac`,
which hands it to eos through `eos_set_debye_huckel_z`;
`use_two_z_tables` is decided in `run_yrec` right after the controls
read; `tollaol`/`use_pure_z_table` live in `star%ctrl`; the fcorr0 sign
mode is a named logical.

**R5 (derived types, single implementations).** kap: `kurucz_table_set`,
`laol_table_set`, `opal92_table_set`, `surface_p_table` with two
instances each; `kurucz2.f90`, `gtlaol2.f90`, `opal92_interp2d_z2.f90`,
`opal92_interp3d_z2.f90`, `kcsurfp.f90` deleted; OPAL95 interpolators
take their stencil explicitly. eos: `opal_eos_vintage` (2001/2006) and
`opal_eos_vintage95` instances; `quadeos01/06` -> `quad`, `t6rint_core`,
`eqbound_core`, `gmass06` deleted; composition-indexed MHD tables
(`mhdst`/`mhdpx2` loops). core: `envelope_refit_lib::append_envelope_points`
shared by `rebuild_envelope`/`read_starting_model`; per-column kernels in
`stitched_model`; two named surface-geometry helpers (one per operand
order). rotation: `equal_grid_lib` (the token-identical grid blocks of
`am_transport_grid`/`composition_grid`), `species_table_lib`. burn/mixing:
`rates` writes one `rr_*`-indexed vector, `eqburn`/`solve_composition`
consume it; `engeb` returns a `burn_result` record; dburn/dburnm
difference table with named thresholds.

## Bugs from SUMMARY 1.1 -- status

Fixed (cosmetic/naming, no output change): the mixwind rows mislabelled
as species, mix.f90's `(12,2)` dummy vs `(13,2)` actual, unused
`massloss` age argument, several stale "recomputed by" comments. Not
fixed (numbers would change; R6): read_starting_model log10-vs-linear
Lsun, rhoofp01 priming flag, envint_kernel reversal, `mod` argument
order, kawaler/wind_spindown saturation gate, dlum_dalpha,
henyey_iterate precedence. Found on the way (all documented in place,
none fixed): the OPAL 1995 stencil overrun (`esac.f90` reads one past
`nt`; `t6rinterp` reads `t6_grid(nt+1)`, i.e. the next member -- this is
why the 1995 vintage keeps its own exact-dimension type); `sulaol`
log10-converts the second LAOL T grid over the first table's `num_t`;
read_starting_model's refit uses `star%solar_luminosity_cgs` where
rebuild_envelope uses the log form; `wind_spindown_matt` `gl`
uninitialised.

## Open for R6 (numbers change -- one reseed, author decisions first)

Author decisions still unanswered:
1. `dlum_dalpha` 0.0139 vs 0.139 (setup_solar_calibration).
2. `alfmlt/phmlt/cmxmlt` (always zero in temperature_gradients): make
   live or delete.
3. BS-extrapolation path (`burn_settle_mix`, `burn_mix_extrapolated`,
   `evolve_angular_momentum` BUR-ST branches, `difad_shear_coeff2`
   consumers): delete or keep; ~700 lines untouched pending this.
4. `henyey_iterate` (formerly :160) intended operator precedence.
5. `llaol` is never assigned (always `.false.`); delete it and the guard
   at `core/read_starting_model` (the guarded branch is dead).
6. `star%dt = -dabs(star%dt)` as "do not evolve": the negative value
   reaches `matt_wind` arithmetically, so a logical is not byte-safe
   (henyey-wave2 Deferred 2 has the proposal).
7. Microdiff Bahcall-Loeb unit copies: the in-place `a*s/s` restore is
   inexact (about 9 % of masses differ by 1 ulp today), so `intent(in)`
   caller arrays change numbers (sidechan-wave2 Deferred 3).
8. One operand order for the six `log10 R(L, Teff)` sites
   (core-wave3 Deferred 4 lists all six).
9. CZ-extent gating in the four rotation copies (the i1 fix is in only
   one; rotation-wave3 Deferred 3).
10. `theta_term_n/p` (`x + 0.0d0`) in am_advection_diffusion_coeffs.
11. `eqstat2` blend extraction (evaluation-order change) and the beta14
    convention discrepancy.

Class-A leftovers that need a state/io owner (exact proposals in the
reports): delete `star%converged_zone`, `star%old_shell_mass`,
`star%iovim`; `star%job%alex95_table_paths` -> `alex94_table_paths`;
`masschg2/3` into `star%`; `max_runs` into star_info_lib;
`gyr`/`3.1558d7` seconds-per-year duplicates; the io/ profile-column
expressions duplicated from core (core-wave3 Deferred 2); the
`ix_*`/`n_mix_species` duplication between eos_mixture_lib and
star_info_lib; `rr_*` indices moving from rotation_scratch_lib into
net_lib (mixburn-wave3 Deferred 4); the burn-region iterator and
reaction table proposals (mixburn-wave3 Deferred 2-3); the remaining
legacy `rot_scr` member names whose meaning could not be established.

## Corrections to the audit recorded by the agents

Collected in each report's last section. The recurring ones: line
numbers drifted after R1; several "duplicate" pairs were not
token-identical (envelope refits differ in six places, the four
CZ-extent searches differ in gating, `opal92_interp3d_z2` did have the
X-grid block, the q-tables are not 13 long); `aintt` was reported
removed in wave 1 but was not (removed in wave 3); `heburn` and
`burn_mix_extrapolated` do not call `rates`; `mhd_output` slots were
already named by wave 1.

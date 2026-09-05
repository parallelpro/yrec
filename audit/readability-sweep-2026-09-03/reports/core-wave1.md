# core -- wave 1 readability (R1-R3)

Worktree `/Applications/YREC-wt/core`, branch `rs/core`, base `0d5c11a`.

| batch | commit | verification |
|---|---|---|
| R1 comments & dead code | `8d6d9b6` | gate1 IDENTICAL; 37-case pins `PINS EXIT 0` (`logs/core.pins.R1.txt`) |
| R2 named indices/constants | `ec13358` | gate1 IDENTICAL; 37-case pins `PINS EXIT 0` (`logs/core.pins.R2.txt`) |
| R3 shared helpers | `02b9418` | gate1 IDENTICAL; 37-case pins `PINS EXIT 0` (`logs/core.pins.R3.txt`); aux battery 22 passed `AUX EXIT 0` (`logs/core.aux.txt`) |

Every build was `-O3 -ffp-contract=off` via gate1 (`make clean` before the
R2 and R3 gate runs because a module procedure / parameter list changed).
`src/deps.mk` was regenerated in R1, R2 and R3 (unused `use` removed,
`rebuild_envelope -> envint_kernel`, `read_starting_model` and
`neutrino_flux_table -> observables_lib`).

## R1 -- comments and dead code (`8d6d9b6`, 26 files)

Headers rewritten (what the routine does now; no label numbers, COMMON
names, or retired-routine names): main, run_yrec, evolve_step,
henyey_iterate ("formerly crrect"), henyey_coefficients ("formerly
coefft"), henyey_solve ("formerly hsolve"), read_starting_model,
rebuild_envelope ("formerly getnewenv"), check_solar_calibration,
check_star_calibration, surfbc, shell_physics, envint_lib, yrec_capi
(re-entrancy paragraph rewritten), observables_lib (botched-sed line
fixed), point_scratch_lib, stop_conditions, neutrino_flux_table.

Commented-out F77 removed: STARIN/HTIMER/LIBURN/CHKCAL/OVROT/QFPR/
TAUINTNEW/OPEN/WRITE call sites, "G Somers END"/"JVS end" markers, bare
CONTINUEs, LONG/LSHORT and IF(LEXCOM) wrappers in henyey_coefficients
(the four extended-species stores kept and re-indented), the 13-17
line blank blocks in rebuild_envelope/henyey_iterate/henyey_coefficients.

Dead code, each proven by grep over the whole `src/` tree:

- evolve_step: locals `nao`, `use_correct_gradients`.
- run_yrec / monte_carlo: `write_run_summaries` lost the intent(inout)
  `surface_z_over_x` and the write-only `monte_helium_diffusion_fraction`.
- envint_kernel: parameters nt/ng/ntc/ngc/nta/nga, the `xyz(22)` data
  array, write-only `unused_chdelj/unused_chdeld`; predicted next mass
  gets its own local `mass_next_step` instead of reusing `swap_temp`;
  `want_derivatives` set once.
- observables_lib: locate_core_cz's `core_boundary_fx2 = x/x` block and
  the two values nothing read (SUMMARY #18); `envelope_boundary_fx`
  became a local of locate_surface_cz_base (module SAVE and the reset
  line removed; `ksaha_center = 0` kept).
- turnover_timescale: third kspline/ksplint block (`log10_radius_interp`
  never read; ksplint's only side effect is an error stop on the same
  xa/x as the two kept calls).
- surfbc: `numenv` DATA counter, unused `use atm_lib`, duplicate `use
  star_info_lib`; `envelope_coeffs(i+i+i-3+j)` -> `envelope_coeffs(3*(i-1)+j)`.
- read_starting_model: `kk` and 23 physics scalars declared never read;
  nts/nps parameters; unused `use atm_lib`; `integer*4` -> `integer`.
- shell_physics: write-only `local_conductive_opacity_flag`; `saha_state`
  zeroed explicitly once before the loop (same value as the
  -finit-local-zero read).
- check_solar_calibration: write-only `log_zx_mismatch_prev`.
- henyey_coefficients: duplicate stores of `star%conv_vel(im)` and
  `star%pulse_dlnrho_dlnp/dlnt(im)` (overwritten with the same value
  before any read); unused imports i_grad_*/json.
- henyey_iterate: unused `json`; star_setup: unused `star`;
  stitched_model: unused i_grad_* (3 routines) and `star` in
  species_slot; neutrino_flux_table: unused i_metals/i_n15/i_o17;
  rebuild_envelope: unused `use atm_lib`, `species_index` renamed
  `env_point_index` (it indexes envelope points, 42 sites).

SUMMARY 1.1 items: #2 (envint_kernel swaps `env_luminosity` both ways;
constant along the envelope so output unchanged), #8
(`start_new_triangle .or. (reset_triangle .and. iteration_level.eq.2)`,
parenthesised as parsed today), #10 (comment only, at the
valfmlt/vphmlt/vcmxmlt copies).

## R2 -- named indices and constants (`ec13358`, 9 files + deps.mk)

- henyey_coefficients: `star%xa(1..15,im)` -> `star%xa(i_h1..i_be9,im)`
  (all 15 slots, incl. the `xa(i_h1,im).gt.0.01d0` test).
- shell_physics: `composition(1..4,im)` and `/atomic_weight(1..4)` ->
  i_h1/i_he4/i_metals/i_he3; comment on the atomic_weight data order.
- monte_carlo: `star%xa(3/1/2,..)` -> i_metals/i_h1/i_he4 in
  write_monte_carlo_model.
- rebuild_envelope: `composition(1..11,..)` -> i_* names; the thin-envelope
  clamp `min(target_envelope_mass, -1.0d-12)` -> `senv_thin_envelope`.
- envint_kernel: `double precision, parameter, public :: senv_thin_envelope
  = -1.0d-12` replacing the two inline tests (SUMMARY area item 6).
- read_starting_model: `reference_composition(5..11)` -> i_c12..i_o18;
  the CNO renormalisation loops `do i = i_c12,i_c13` / `i_n14,i_n15` /
  `i_o16,i_o18`; four implied-do write lists `k=i_c12,i_o18`,
  `k=i_he3,i_be9`.
- check_star_calibration: `log_r_rsun_current` -> `r_rsun_current` (holds
  R/Rsun, not its log; area item 11); `x_calibration_step = 0.01d0`.
- stop_conditions: `limit_active_below = 0.9d99` for the "limit not set"
  sentinel tests (area item 23, the rename half only).
- stitched_model: `max_ext = 3*json`; `ip_gamma1 = 10, ip_delta = 51,
  ip_mu = 52`, `ip_seismic_first/last`; compute_seismic_columns and
  ext_profile_value use the ip_* names instead of 2,3,4,5,6,10,12,14,51,
  52,53 (area item 13).

## R3 -- shared helpers (`02b9418`, 7 files + deps.mk)

All helper bodies are token-for-token the call-site text; `use math_lib`
is present in every helper that evaluates `exp` (math_lib shadows the
intrinsic with the crmath backend, so the resolution is unchanged).

- observables_lib: new public `shell_masses_from_log_mass(log_mass,
  log_total_mass, nz, m, dm)`; read_starting_model's build_shell_masses
  and neutrino_flux_table now call it (area item 21 / structural 5). The
  two sites were the same stencil under a consistent renaming of the
  three rolling scalars. rebuild_envelope's copy uses exp10 and a
  different operand order and is left alone.
- envint_kernel: the 21 explicit swap triplets of the envelope inversion
  replaced by `arr(1:n) = arr(n:1:-1)` per array (`(:,1:n)` for the two
  (3,json) arrays) inside `associate (n => env_struct%num_env_points)`;
  `swap_temp`, `swap_temp_logical`, `inversion_index1/2` deleted (area
  item 6; pure copies, no arithmetic).
- evolve_step: contained `double precision function
  envelope_overshoot_depth()` holds the if/else pressure-scale-height
  formula; rezone_or_snapshot and burn_light_elements assign
  `star%pressure_scale_height_start/_end = envelope_overshoot_depth()`
  (area item 15).
- stop_conditions: new public `check_rotation_initialised(ierr)` with the
  1611 format/print/ierr = 1; run_yrec and evolve_step replace their
  duplicate blocks (formats 1611 / 18) by `call
  check_rotation_initialised(ierr); if (ierr /= 0) return`. Print text,
  unit and order unchanged; both callers have ierr = 0 on entry (each
  is preceded by `if (ierr /= 0) return`) (area item 16).

## Deferred (cross-domain)

Exact proposals; none applied because the declaration lives outside
the assigned files.

1. `state/star_info_lib.f90:143` `xa(15,json)`: add
   `integer, parameter, public :: n_species_basic = 11,
   n_species_extended = 15` next to the i_* slots. Then replace
   `num_species = 11 / = 15` in `core/henyey_iterate.f90:130-131`,
   `core/evolve_step.f90:239-240`, `core/read_starting_model.f90:616-617`
   and `species_end_index = 15 / = 11` in `core/rebuild_envelope.f90:86-88`;
   also the `do k = 4,...` loops there could start at `i_he3`.
2. `state/star_info_lib.f90:156` `luminosity_breakdown(8)`: add
   `n_lum_channels = 8`; users `core/henyey_iterate.f90:175`,
   `core/henyey_coefficients.f90:118`, `core/observables_lib.f90`
   (renormalize_luminosity_breakdown `do i = 1,8`).
3. `state/star_info_lib.f90:367` `neutrino_flux(10), neutrino_flux_total(10)`
   and `neutrino_flux_zone(10,json)`: add `n_nu_fluxes = 10`; users
   `core/henyey_coefficients.f90:104,231`, `core/neutrino_flux_table.f90:59,91`.
4. `state/star_info_lib.f90:116-119` nine Monte-Carlo arrays `(1000)`:
   add `max_mc_runs = 1000`, use it there and in
   `core/monte_carlo.f90:44` `min(star%job%mc_run_end,1000)` (area item 24).
5. `mixing/mix.f90:48` dummy `mixed_zone_bounds_no_overshoot(12,2)` vs the
   actual `(13,2)` in `state/star_info_lib.f90:273`: declare the dummy
   `(13,2)` (or assumed-shape) so the shape matches (area item 17 tail).
6. kenv/katm counters (area item 4): zeroed and incremented in
   henyey_iterate, rebuild_envelope, read_starting_model, never read.
   Their signature chain runs surfbc -> atm_get (envint_lib) ->
   envint_kernel -> numerics bsstep -> derivs, so deleting them touches
   numerics; propose dropping the two arguments from all six signatures
   in one cross-domain commit.
7. Surface-radius-from-L-Teff (area item 8) and spot-adjusted Teff (area
   item 9) helpers: consumers include `atm/temperature_gradients` (spot
   Teff) and wind/mixing (log R). A shared `log10_radius_from_l_teff`
   would have to keep each caller's operand order to stay byte-identical
   (the six sites differ), so this is a B item unless done per-site.
8. henyey flag relay (area item 3) and the atm_get 22-argument signature
   (area item 5) are structural (R4+), not attempted here.
9. `core/henyey_iterate.f90:168-180` renormalises luminosity_breakdown
   with a `total .gt. 0` guard; `observables_lib.f90`
   renormalize_luminosity_breakdown has no guard. Not merged: the guard
   would change behaviour on a zero total (division vs skip). A merge
   needs a decision on which semantics is wanted.

## Reverted (changed numbers)

None. Every gate1 run reported IDENTICAL on the first build and all
three full pin runs passed on the first attempt.

## Skipped (disagree with reviewer / too risky)

- area item 22 (`if (.false.) print *` warning silencers in
  envelope_derivs/atmosphere_derivs): removing them re-introduces
  unused-dummy warnings; changing the signatures is cross-domain.
- area item 14 (sign-encoded flags `star%dt = -dabs(star%dt)`,
  `fcorr0` sign as mode): behaviour-preserving rewrite needs a new
  state member and touches henyey_iterate's contract; R4 material.
- area item 19 (turnover_timescale reading stx_prof in place, explicit
  ierr): restructure, not a small edit; only the dead third spline was
  removed (R1).
- area item 24 `convergence_iterations.ge.11` in monte_carlo: the
  meaning of 11 is not documented anywhere; naming it would be a guess.
- area item 1 (envelope refit duplication read_starting_model vs
  rebuild_envelope): the two copies differ in print units and in which
  step sizes they restore; a merge is a behaviour decision.
- area item 12 (check_solar_calibration empirical derivatives as
  parameters): they are stored in star% and read by the io domain.
- area item 10, `surfbc` dummy `luminosity_linear` (it is log10 L):
  noticed after R2 was committed; a dummy-only rename, safe to do in the
  next naming pass.
- evolve_step still carries `use math_lib` in rezone_or_snapshot and
  burn_light_elements although the only math there moved into the
  helper; harmless, left to avoid a rebuild after the R3 pin run.

## Things the audit missed

- henyey_coefficients stored `star%conv_vel(im)` and
  `star%pulse_dlnrho_dlnp/dlnt(im)` twice per zone with the same value
  (dead first store); removed in R1.
- observables_lib's `envelope_boundary_fx` was module SAVE state with a
  reset entry only because locate_core_cz once read it; after the
  dead-block removal it is a plain local.
- turnover_timescale's third spline evaluation was dead but its ksplint
  call could still error-stop; the kept two calls stop on the same
  condition, so removal is output-neutral.
- neutrino_flux_table imported i_metals/i_n15/i_o17 without using them
  (not flagged by -Wall, which does not report unused imported
  parameters or unused whole-module `use`).
- Concurrency hazard for the multi-agent sweep: the scratchpad directory
  is shared by all seven agents, and a generic script name
  (`scratchpad/r2.py`) was overwritten by another domain's script
  between my write and my run. Nothing in this worktree was affected
  (git status verified), but wave-2 assignments should mandate
  domain-prefixed scratch names or inline heredocs.

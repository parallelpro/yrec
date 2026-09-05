## NAME = core

Worktree `/Applications/YREC-wt/core`, branch `rs/core`. This is batch **R5 (structural)** for the core/atm chain. Read first: audit `core.md` (item 1 envelope refit; items 8-10; "Structural" 1, 4, 5), `reports/core-wave1.md` (Skipped: area item 1; Deferred 7; the `surfbc luminosity_linear` note), `reports/henyey-wave2.md` (what R4 already changed in `atm_get`/`envint_config`, and its Deferred list -- any core-only leftover there is yours), `reports/sidechan-wave2.md` (the `temperature_gradients` `zone_index` argument and the `star%job%` driver-side edit it deferred to `core/evolve_step.f90`), `SUMMARY.md` section 3 (R5).

You own `src/core/**` and `src/atm/atmstruct_lib.f90`/`envstruct_lib.f90` if an envelope item needs them (the `kap` agent owns the rest of `src/atm/**`; do not touch `atm/tables`, `atm/ttau_lib.f90`, `atm/atm_table_lib.f90`). `src/io/**` and `src/numerics/**` are not yours: the profile-column item stays inside core. `src/util/compute_timestep.f90` and `src/wind/mdot.f90` are the `mixburn` agent's for a different item; for item 3 you edit only `core/evolve_step.f90`.

### Items, in order

1. **One envelope refit** (core item 1 / Structural 1). `core/read_starting_model.f90` (the refit block around its `atm_get` call and the shell rebuild after it) and `core/rebuild_envelope.f90` implement the same envelope refit twice. Wave 1 recorded two differences: print units, and which `env_step_*` sizes each restores (the wave-2 `envint_config` argument may have removed the second one -- check henyey-wave2.md). `diff` the two blocks after normalising names. If the remaining differences are (a) the output unit and (b) nothing else, give `rebuild_envelope` a `log_unit` argument (and any other value that differed as an explicit argument) and make `read_starting_model` call it; every `write` keeps its format and its unit value; the shell arrays must be filled in the same order with the same expressions. If any arithmetic differs, do not merge: extract only the identical inner part (e.g. the per-shell fill loop) into a private helper both call, and put the exact remaining difference under Deferred with file:line -- that is then a behaviour decision for the author.

2. **Profile-column registry, core side only** (core Structural 4). `core/stitched_model.f90` `profile_value` (a `select case` over column names) and `compute_seismic_columns` compute the same quantities the io writers compute. Within core: collect the column-name -> expression mapping into one contained function used by both core sites, keeping each site's operand order (if the two core sites differ in an expression, keep both and name the difference). Do not touch io/; list under Deferred which io/ routine duplicates which core expression (file:line pairs) for the io owner.

3. **Driver-side `star%job%` flips** (sidechan-wave2 Deferred): if sidechan prepared `compute_timestep`/`mdot` to return a logical instead of writing `star%job%...` and wrote the exact `core/evolve_step.f90` edit under Deferred, apply that edit now -- only the driver line(s) it specified -- so the callee-side write can be removed by `mixburn` (tell them in your report which line you changed; they remove the callee write). If sidechan did not prepare it, skip.

4. **`surfbc` dummy `luminosity_linear`** (core-wave1: it is log10 L): rename the dummy (and the matching local names in `surfbc.f90` only) to `log_luminosity`; callers pass the same actual. Pure rename; grep that no other file declares a dummy of that name that would be confused.

5. **Geometry helpers, per site** (core Structural 5, wave-1 Deferred 7): the six `log10_radius_from_L_Teff` sites differ in operand order, so a shared helper is class B -- do NOT unify. Instead, where two or more sites inside core/ are token-identical, use one contained function for those; leave the others and list every site with its exact expression under Deferred so the author can decide the canonical order (R6).

6. **Leftovers from henyey-wave2**: anything in its Deferred list that is core-only and class A (e.g. a `read_starting_model` `star%iovim = -1` shown dead by sidechan, a sign-encoded flag proposal that it judged safe but ran out of time on). Do them last; each its own small commit.

### Files you may edit

`src/core/**`, `src/atm/atmstruct_lib.f90`, `src/atm/envstruct_lib.f90`, tests under `src/core/test` if any, `src/deps.mk`.

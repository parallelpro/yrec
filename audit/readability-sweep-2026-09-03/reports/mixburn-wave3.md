# Readability wave 3 -- mixburn (batch R5, structural: mixing + burn/net)

Worktree `/Applications/YREC-wt/mixburn`, branch `rs/mixburn`, base `6449462`
(yrec-modern with waves 1, 2a, 2b merged).

## Commits

| hash | item | summary |
|---|---|---|
| `796be57` | 1 | relay nuclear rates as one rr_*-indexed vector |
| `ac61ec9` | 2 | engeb returns its energy outputs as a burn_result record |
| `dfb592d` | 3 | dburn/dburnm interface table and named constants |

Verification: gate1 IDENTICAL before every commit (`make clean` at the start
and before item 2's build, which changes a module-procedure signature and
adds a derived type). On the final tree (`dfb592d`): full 37-case pin
selection -- see the "Final verification" section at the end; aux battery --
same. Only source files and `src/deps.mk` committed; no `test_*.short` left
behind.

Items 4, 5 and 6 produced no commits: 4 cannot land on this branch (the
callers that pass the new actual live on `rs/core` only), 5 is deferred with
findings, 6 is deferred with a mechanical proposal. Details below.

## Item 1 -- nuclear-rate relay as one vector (`796be57`)

Chosen form: `rates` writes one explicit-shape vector
`rate_vec(num_rate_rows)` for one zone; each caller passes the column
`rate_by_zone(:, zone_idx)` of an explicit-shape local
`rate_by_zone(num_rate_rows, json)` (sequence association: a whole column of
a column-major array is contiguous, so the callee's explicit-shape dummy
receives the same 15 doubles it used to write one at a time). Nothing is
copied on the way in or out.

- `src/net/net_lib.f90`: new `integer, parameter, public :: num_rate_rows = 15`
  next to `num_reactions = 13`. `rates(log_density, ..., o18_fraction, rate_vec)`
  replaces the fifteen `(json)` output arrays plus the `zone_idx` dummy.
  `zone_idx` was used only as the store subscript of those fifteen arrays
  (grep of the body: no other reference), so it is unread once the arrays
  are gone. The fifteen store statements at the end of the routine are the
  same expressions in the same order, each now assigned to
  `rate_vec(rr_<row>)`; rows 9 and 13 (`rr_zero9`, `rr_zero13`) are still
  zero-filled placeholders and say so in a comment. The `rr_*` names come
  from `rotation_scratch_lib` (`use rotation_scratch_lib, only: rr_...`),
  which was already the single index scheme for this vector; all fifteen
  names already existed, none was added.
- `src/mixing/solve_composition.f90`: `rate_by_zone(num_rate_rows, json)`
  dummy replaces twelve positional `(json)` rate dummies; 24 reads changed
  `rate_x(k)` -> `rate_by_zone(rr_x, k)`, same row for the same quantity.
- `src/core/burn_lib.f90` `eqburn`: `rate_by_zone` first dummy replaces nine
  positional rate dummies; 18 reads renamed by row. Zone-block reads are
  unchanged.
- `src/mixing/mix.f90`: local `rate_by_zone(num_rate_rows, json)` replaces
  the fifteen `(json)` locals; `call rates(..., rate_by_zone(:, zone_idx))`;
  the clear loop is one slice assignment `rate_by_zone(:, clear_idx) = 0.0d0`;
  `rot_scr%reaction_rate_by_zone(:, k) = rate_by_zone(:, k)` for
  `k = 1..star%nz` replaces the fifteen-line element copy (same row order,
  same columns); `solve_composition` and `eqburn` get the array.
- `src/mixing/rotmix.f90`: the twelve `(json)` unpack locals and the loop
  that filled them from `rot_scr%reaction_rate_by_zone` are gone;
  `solve_composition` receives `rot_scr%reaction_rate_by_zone` directly at
  both call sites. Safe because solve_composition reads only zones inside
  the radiative/convective bounds, all <= the `num_zones` columns that mix
  filled.
- `src/util/timestep_limit_hburn.f90`: same local-array pattern for the
  single shell-midpoint column; `eqburn` call updated.
- `src/net/test/test_net.f90`: `rate_vec(num_rate_rows)` replaces the fifteen
  `json = 5000` arrays; the three print statements print the same twelve
  rows in the same order with the same format. `./test_net` output is
  `cmp`-identical to `src/net/test/expected_test_net.out`.
- `src/rotation/rotation_scratch_lib.f90`: only the comment block above the
  `rr_*` parameters (the row layout is now "what rates() writes"); the
  `reaction_rate_by_zone(15,json)` declaration and the index values are
  untouched.
- `src/deps.mk` regenerated (net_lib.o now depends on
  rotation_scratch_lib.mod; solve_composition.o and test_net.o gained the
  same).

New `-Wsurprising` "moved from stack to static storage" warning for
`rate_by_zone` in mix.f90 and timestep_limit_hburn.f90 (a 15 x 5000 double
local exceeds `-fmax-stack-var-size`). Value-neutral: every element read is
written earlier in the same call (mix clears columns nz+1..json each call
before any read; timestep_limit_hburn writes and reads one column). The same
warning class pre-exists for many other `(json)`-sized locals in the tree.

## Item 2 -- `engeb` output record (`ac61ec9`)

`src/core/burn_lib.f90`: `type, public :: burn_result` with eight
components, all default-initialised to `0.0d0`:
`pp_chain_energy_gen, he3he4_be7_electron_energy_gen,
he3he4_be7_proton_energy_gen, cno_cycle_energy_gen,
triple_alpha_energy_gen, deps_dlnrho, deps_dlnt, total_energy_gen_rate`
(the former first eight positional dummies EPP1, EPP2, EPP3, ECN, E3AL, PEP,
PET, SUM1). `engeb(res, log_density, ...)` takes it as
`type(burn_result), intent(inout) :: res` -- inout, not out, so the
low-temperature early return keeps behaving exactly as before (it zeroes
only `deps_dlnrho`/`deps_dlnt` and leaves the energy components and the
total untouched). The default initialisation replaces the callers'
`-finit-local-zero` scalars: gfortran does not zero derived-type locals
without `-finit-derived`, so the record must carry its own zeros.

Inside engeb the 30 references (the early-return block and the two
internal subroutines `compute_energy_generation` /
`compute_neutrino_emission`) were renamed `x` -> `res%x`; no statement moved,
and the inline rate block (the bit-drift note at burn_lib.f90:1166-1169) is
untouched. The three callers:

- `src/core/henyey_coefficients.f90`: `type(burn_result) :: burn`; after the
  call, `zone_dlnepsilon_dlnrho = burn%deps_dlnrho`,
  `zone_dlnepsilon_dlnt = burn%deps_dlnt`,
  `energy_gen_rate = burn%total_energy_gen_rate`, and the five
  `energy_gen_component(1..5) = burn%...` copies -- the two
  `zone_dlnepsilon_*` locals survive because they are read further down
  outside my hunk (lines ~250-252, 380-381, 396-397).
- `src/core/neutrino_flux_table.f90`: eight scalar locals -> `burn`; nothing
  downstream reads them (they were discarded there before, too).
- `src/util/timestep_limit_heburn.f90`: `burn` plus
  `helium_energy_gen = burn%triple_alpha_energy_gen`, which the existing
  `if(helium_energy_gen.lt.1.d-22)` floor then reads as before.

The brief's `eps_nu` and `neutrino_flux(n_nu_fluxes)` components were NOT
added: those outputs already live on `star%` (`star%neutrino_loss_rate`,
`star%neutrino_flux(:)`, `star%alpha_capture_energy`), where the output
writers (outside my files) read them. Moving them would have touched io/.
The `compute_neutrino_emission` split into three internal subroutines was
not needed for this item and was not done.

## Item 3 -- `dburn`/`dburnm` interface note (`dfb592d`)

`src/core/burn_lib.f90` only.

- The dburn header's prose list became a side-by-side table (columns
  dburn / dburnm) of the four differences: timestep units (Gyr as passed /
  seconds converted by `timestep*gyr_per_year/seconds_per_year`),
  "nothing to burn" threshold (accretion-mixed `deuterium_fraction_test`
  `.lt. dburn_min_deuterium` = 1.0d-11 / plain CZ average
  `.lt. dburnm_min_deuterium` = 1.0d-14), accretion weighting (absolute
  masses `(D*M_cz + D_acc*f)/(M_cz + f)` / normalised
  `(D + D_acc*f_sub)/(1 + f_sub)` with `f_sub = step_fraction*star%accreted_mass_fraction`),
  rate source (`star%deuterium_burning_rate(_start)` / the dummies
  `deuterium_rate_start/_end`). The header says explicitly that the bodies
  are kept separate on purpose.
- Three private module parameters at the top of burn_lib:
  `dburn_min_deuterium = 1.0d-11`, `dburnm_min_deuterium = 1.0d-14`,
  `gyr_per_year = 1.0d-9`. Each replaces exactly one literal of the same
  value and kind in place (dburn's threshold test, dburnm's threshold test,
  dburnm's `timestep_gyr = timestep*gyr_per_year/seconds_per_year`);
  operand order unchanged. The other `1.0d-11` literals in burn_lib
  (deutrate's and liburn/liburn2's guards) are unrelated and untouched.
- Bodies not unified (class B).

The call sites confirm the units column: mix.f90 passes `dt_gyr` to dburn;
mid_timestep_model.f90 passes `sub_timestep` (seconds) to dburnm.

## Item 4 -- callee-side `star%job%` writes: NOT applied, reason

`reports/core-wave3.md` exists and confirms the driver edit: commit
`6b50553` on `rs/core` adds the local `disable_structure_dt_limits` in BOTH
`core/evolve_step.f90` and `core/run_yrec.f90`, pre-sets it `.false.`,
passes it as the new trailing actual of `compute_timestep`, and flips
`star%job%use_structure_dt_limits` when it comes back `.true.`.

But `6b50553` is on `rs/core` only (`git branch --contains` lists just
`rs/core`), and my branch is based on `6449462`, where both callers still
end their `compute_timestep` actual lists at `log_Teff`. `compute_timestep`
is an external procedure with no interface block, so if I add the trailing
dummy here and store `.false.`/`.true.` into it, the callee writes through a
stack slot the callers never filled -- undefined behaviour, and the pins in
this worktree could not be byte-identical (or would crash) with the
`star%job%` write removed. Per the shared rules (no rebase, no merge, do
not touch `core/evolve_step.f90`), the callee edit cannot be made
self-consistent on this branch. It is left in place, and the exact edit is
under Deferred 1 so it can be applied in the merge commit together with or
after `6b50553` (core's report says the same about merge order).

`wind/mdot.f90` needs nothing: it already returns `disk_exhausted`
(sidechan wave 2) and has no `star%job%` write. The remaining callee-side
write for that flag is `wind/massloss.f90:283`, which is not in my file
list (core's Deferred 3 covers the driver half).

## Item 5 -- burn-region iterator: deferred, with findings

Not implemented. The instrumented proof would have been feasible, but the
sites do not support the audit's rationale, and the mechanical replacement
has a pitfall that the audit's "once per call" description misses. What I
found:

1. The five burn loops are: mix.f90 lines 178-200 (`solve_composition`),
   222-281 (`eqburn`, explicit H/He burning), 442-463 (`dburn`), and
   rotmix.f90 lines 76-99 (`solve_composition`) and 182-195. In each, the
   radiative pass visits shells one at a time, per radiative region, with
   `exit` on `logT <= nuclear_logT_cutoffs(1)` (which ends the inner loop
   only; the next radiative region starts again), then the convective pass
   visits each mixed zone as a unit. An ordered `(begin, end, is_convective)`
   list built from `radiative_zone_bounds`, `mixed_zone_bounds` and the
   cutoff would reproduce this sequence exactly.
2. The list cannot be built "once per call" in mix: `find_convection_zones`
   is called again at mix.f90:428-436 (`dpenv < 1` and `iteration_level > 1`)
   between the eqburn loop and the dburn loop, rewriting both bounds arrays.
   The list would have to be rebuilt there. rotmix uses the mid-step
   `log_temperature` (rot_scr), not `star%logT`, so mix's list cannot be
   reused by rotmix either.
3. The eqburn loop's radiative and convective bodies differ (per-shell
   `star%xa(i_h1,inner_zone_idx) = hydrogen_fraction + dx_dt*dt_gyr` vs
   one averaged update written to every shell of the zone), so with a
   region list that loop still needs an `if (is_convective)` with both
   bodies. Only the four solve_composition/dburn loops collapse.
4. The audit's stated benefit -- "one place to encode the single-shell CZ
   is treated as radiative unless it is the centre rule" -- is not
   delivered by a burn iterator, because that rule is not in the burn
   loops. It is in `find_convection_zones.f90:91`
   (`if (zone_start.ne.zone_idx-1)`: a single-shell CZ is never recorded in
   `mixed_zone_bounds`, centre included, so it lands in the radiative
   bounds and burns as a radiative shell) and in the two homogenise sites
   `mix.f90:307` and `homogenize_convection_zones.f90:70`
   (`if (mix_start.ne.1 .and. mix_start.ge.mix_end) cycle`: skip single-shell
   zones unless they start at shell 1). The "unless it is the centre"
   exception exists only at the homogenise sites; find_convection_zones
   drops central single-shell zones too, so the exception can only fire if
   `semiconvection` or `overshoot_boundaries` later creates a
   single-shell zone starting at 1. The three sites are therefore
   consistent by accident, not by construction. Centralising the rule means
   changing how the bounds are built, which is a class-B question for
   find_convection_zones, not for mix.

Given 2-4, the iterator would add a module, a derived type and a
`json + max_convective_zones` scratch list to remove four short nested
loops, without reaching the rule it was meant to centralise. Proposal, if
the author still wants it, under Deferred 2.

## Item 6 -- reaction table: deferred, mechanical proposal

Not attempted (optional, and the finding below changes its shape). The
parallel tables are NOT all 13 long: `charge_product, z53, z43, z23, z86`
have `num_reactions = 13` entries (screening, one per reaction), but
`q1..q5` have `iq_be7p = 8` entries and `q6..q8` have `r_po16 = 7` (the
Bahcall Q-form fits for the seven proton-capture reactions plus the Be7+p
slot). A single `reaction_data(num_reactions)` array would need padding of
the Q-forms for reactions 8-13, which is a change of table content, and the
engeb copy of the rate block is class B in practice (burn_lib.f90:1166-1169).
Proposal under Deferred 3.

## Deferred

1. **Item 4 callee edit for `util/compute_timestep.f90`** (apply together with
   or after `rs/core` `6b50553`):
   - line 44: `log_teff)` -> `log_teff, disable_structure_dt_limits)`;
   - after line 68 (`double precision, intent(in) :: log_teff`): add
     `logical, intent(out) :: disable_structure_dt_limits` with a comment
     ".true. when the structure-based limits were switched off in this
     call; the driver flips star%job%use_structure_dt_limits";
   - first executable statement (before line 80's
     `if (previous_timestep.ge.0.0d0)`): `disable_structure_dt_limits = .false.`
     -- it must be set on every path, including the timestep-override
     branch and the negative-timestep branch;
   - line 95: replace `star%job%use_structure_dt_limits = .false.` with
     `disable_structure_dt_limits = .true.`;
   - update the comment at lines 71-75 (drop "the star%job% flag is still
     written here (2026 W2) until the driver takes over the flip").
   The driver already pre-sets `.false.` and applies the flip after the
   call, so the numbers are unchanged.

2. **Burn-region iterator** (item 5), if wanted despite the findings:
   new `src/mixing/burn_regions_lib.f90` with
   `type burn_region; integer :: zone_begin, zone_end; logical :: is_convective; end type`
   and
   `build_burn_regions(log_temperature, logT_cutoff, radiative_zone_bounds, num_radiative_zones, mixed_zone_bounds, num_mixed_zones, regions, num_regions)`
   filling `regions(1:num_regions)` in the order: for each radiative region
   in order, shells `b..e` one by one until the first with
   `log_temperature <= logT_cutoff` (then the next region); then each mixed
   zone `(b, e, .true.)`. Replace mix.f90:178-200 and 442-463 and the two
   rotmix loops by `do r = 1, num_regions ... end do`; rebuild the list
   after mix.f90:428-436; leave the eqburn loop (222-281) or give it the
   `is_convective` branch with both bodies verbatim. Proof: a debug print
   of `(zone_begin, zone_end)` before each solve_composition/dburn call on
   the solar pin and the run_config_matrix rotating cases, before and after.

3. **Reaction table** (item 6), net_lib `rates` only (engeb's copy stays
   inline): `type reaction_data` with components
   `charge_product, z53, z43, z23, z86` (13 rows) and a separate
   `type qform_data` with `q1..q5` (8 rows) and `q6..q8` (7 rows) -- or one
   type with the Q-form components padded by an explicit "unused" marker,
   which the author must approve since it changes the table's content.
   Build the parameter arrays with a scratch script that reads the existing
   `data` lines verbatim (same literal text, same order); add to
   `net/test/test_net.f90` a loop printing every old-array element next to
   the new component with `es24.16` and compare exactly; then change reads
   `charge_product(i)` -> `rxn(i)%charge_product` etc. inside the existing
   loops, element by element, no expression merged.

4. **Layering note** (item 1): `net_lib` now `use`s `rotation_scratch_lib`
   for the `rr_*` row names (the rotation agent's brief kept those
   parameters in rotation_scratch_lib and my brief said not to add names
   there). It is a one-way dependency (rotation_scratch_lib uses only
   star_info_lib), so no cycle, but the natural home for the row layout of
   what `rates()` writes is `net_lib`. Proposal: move the fifteen `rr_*`
   parameters and `num_rate_rows` into net_lib and re-export them from
   rotation_scratch_lib with `use net_lib, only: rr_...` plus `public`
   statements (rotation_scratch_lib.f90:17-18 is where the `use` line goes)
   -- values unchanged, so byte-safe, but it edits the rotation agent's
   hunk.

5. **`star%deuterium_burning_rate_start` / `_rate` vs dburnm dummies**: with
   the table in place, the remaining asymmetry that keeps dburn and dburnm
   apart is the rate source and the weighting formula (bugsweep
   core_physics.md:47, known). Nothing to do here without a decision.

## Reverted (changed numbers)

None. gate1 was IDENTICAL after every item.

## Skipped

- Item 2's `eps_nu`/`neutrino_flux` record components: already on `star%`,
  read by io/ writers outside my files (see item 2).
- Item 2's `compute_neutrino_emission` three-way split: not required for
  the record; not done.
- Item 4: see above (cannot be made self-consistent on this branch).
- Item 5, 6: deferred, see above.

## Corrections to the audit / earlier reports

- Brief item 1 says `timestep_limit_heburn.f90` and
  `burn_mix_extrapolated.f90` call `rates`. Neither does: the callers of
  `rates` are `mix.f90`, `timestep_limit_hburn.f90` and
  `net/test/test_net.f90`. `timestep_limit_heburn` calls `engeb` (item 2);
  `burn_mix_extrapolated` calls neither.
- Brief item 2 lists `eps_nu` and `neutrino_flux(n_nu_fluxes)` as engeb
  outputs to move into the record; they were never engeb arguments in this
  tree (they are `star%` members written by engeb's internal subroutine).
- Audit mixwind Structural 2: "build ... once per call" is wrong for mix
  (bounds are recomputed at mix.f90:428-436 mid-routine), and the
  single-shell-CZ rule it wants to centralise is in
  find_convection_zones/homogenise, not in the burn loops (item 5).
- Audit burn Structural 1 describes the q-tables as 13-reaction tables;
  they are 8- and 7-long (item 6).
- Audit burn Structural 5 says "move the timestep_gyr conversion (470)
  and thresholds (336/503) into two named parameters each" -- there are
  three literals in total (one conversion factor, two thresholds), now
  three parameters.

## Final verification

On the final tree (`dfb592d`, worktree clean):

- full 37-case pin selection: `37 passed, 272 deselected in 1577.18s
  (0:26:17)`, `PINS EXIT 0` (`logs/mixburn.fullpins.out`; slower than the
  nominal 19 min because three other worktrees' pin runs shared the
  machine);
- aux battery: `22 passed in 552.16s`, `AUX EXIT 0` (`logs/mixburn.aux.out`);
- no `test_*.short` left under `src/`; only source files and `src/deps.mk`
  committed.

The commit messages carry "gate1 byte-identical" only; the full pins and
aux were run once after the last edit as the rules require, and their
results are recorded here.

Coordinator note received after the aux run: item 4's callee edit is to be
applied by the coordinator on the merged tree from Deferred 1 above --
consistent with what was done (compute_timestep's signature and its
`star%job%` write are untouched on this branch).

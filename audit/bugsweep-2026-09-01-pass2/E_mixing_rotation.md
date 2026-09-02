# Pass 2, section E: convective mixing, overshoot, semiconvection, composition solve, rotation top level

Reviewer discipline: read-only; every file below read in full; F77 ancestors at
rev 6cd5673 (wczimp.f absent there, taken from c95697e) diffed by hand
against the modern text; all bare-external call sites grepped without `head`.

## Files read in full

Assigned:
- src/mixing/burn_mix_extrapolated.f90 (bsrotmix.f)
- src/mixing/burn_settle_mix.f90 (bursmix.f)
- src/mixing/find_convection_zones.f90 (convec.f)
- src/mixing/homogenize_convection_zones.f90 (mixcz.f)
- src/mixing/mix.f90 (mix.f)
- src/mixing/overshoot_boundaries.f90 (ovrsht.f)
- src/mixing/rotmix.f90 (rotmix.f)
- src/mixing/semiconvection.f90 (sconvec.f)
- src/mixing/solve_composition.f90 (kemcom.f)
- src/rotation/am_convective_regions.f90
- src/rotation/enforce_rotation_profile.f90 (wczimp.f / wcz.f)
- src/rotation/equal_spaced_grid.f90 (getgrid.f)
- src/rotation/evolve_angular_momentum.f90 (getw.f)
- src/rotation/mid_timestep_model.f90 (midmod.f)
- src/rotation/omega_from_j.f90 (getrot.f)
- src/rotation/rotation_scratch_lib.f90 (COMMON replacement module)
- src/rotation/solid_body_omega.f90 (solid.f)
- src/rotation/viscos.f90 (viscos.f)

Supporting (read for interface checks):
- src/core/burn_lib.f90 (eqburn, dburnm, liburn, liburn2, lirate88 dummy lists)
- src/mixing/compute_scale_height.f90
- src/rotation/seculr/check_angular_momentum.f90 (tolerance save/restore)
- src/rotation/seculr/am_transport_grid.f90, composition_grid.f90 (equal_spaced_grid callers)
- src/core/shell_physics.f90 (viscos caller), src/core/evolve_step.f90 (evolve_angular_momentum caller)
- src/wind/matt_wind.f90 (dummy list), src/numerics/numerics_lib.f90 (osplin, simeqc)
- src/rotation/mid_timestep_model.f90 callers and dburnm callers

## Verified clean (re-derived, matches F77 and/or physics)

- evolve_angular_momentum (getw.f): the retry/redo control flow (GOTO 40 -> `cycle retry`,
  GOTO 30 -> outer loop `exit`), diffusion sub-step count, disk-lock branch, restore of
  j_rot/amum/xa on redo, xa_start refresh for species 1-11, homogenize + find_convection_zones
  ordering, lirate88/liburn/liburn2 argument lists.
- omega_from_j (getrot.f): identical zone walk; enforce_rotation_profile call arguments.
- solid_body_omega (solid.f): Newton update delta/(I + omega dI/domega) and final
  j = I*omega/dm identical.
- enforce_rotation_profile (wczimp.f): all five branches (solid body; walpcz<=-2 & mode 0 constant j;
  power law modes 2/3 with interior solid body; omega = C R^walpcz) identical, including the
  mode-3 pooling of interior AM with the CZ.
- equal_spaced_grid (getgrid.f): identical; osplin(xval,yval,xtab,ytab,n,k) matches callers.
- viscos (viscos.f): radiative viscosity constant 6.7282653e-26 = 4a/(15c); thermal diffusivity
  16 sigma T^3/(3 kappa rho^2 cp); Coulomb log and molecular viscosity as in the original.
- mid_timestep_model (midmod.f): time interpolation of structure into rot_scr%*_mid,
  core/envelope boundary-zone determination, rotmix argument order, dburnm argument order.
- solve_composition (kemcom.f): 7x8 Jacobian entries (rows for H, He3, He4, C12, C13, N14, O16)
  and the mass-weighted averaging in eqburn; O18 update; Z re-closure when X < 5e-7.
- rotmix: De Morgan of the original `GOTO 170` guard (line 182); settling substep count with
  `.or. .eq.0` guard; gradT save/restore around the settling burn.
- overshoot_boundaries (ovrsht.f) and semiconvection (sconvec.f): boundary extension and zone
  merge logic identical to the originals; homogenize_convection_zones mass-weights by dm
  (no radius weighting anywhere in the mixing path, so no mass-vs-radius mismatch found).
- Interfaces: matt_wind (13 args incl. ierr), shell_physics, lirate88, liburn, liburn2, shape,
  zone_moments_of_inertia, dburnm, rotmix, burn_mix_extrapolated, diffuse_composition_driver,
  find_convection_zones, homogenize_convection_zones -- actual lists match dummy lists at every
  call site under src/.

Root cause of the two known symptoms was NOT located in these files. The
"UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL 1" message is emitted at
solve_composition.f90:398-401 (written with zone_begin, so "SHELL 1" means the core
convection zone starting at zone 1, not shell index 1 in isolation); nothing in kemcom's
Newton loop changed relative to the original, so the trigger is upstream (rates or timestep).

---

## src/rotation/mid_timestep_model.f90:86 -- deuterium burning rate accumulator lost its SAVE; sub-steps 2+ burn D with a zeroed/garbage rate

- class: numerical
- severity: high
- confidence: high
- provenance: modernization -- `git show 6cd5673:src/midmod.f | grep -n "DRATEM\|SAVE"` gives
  `53: DIMENSION DRATEM(JSON),DRATEM0(JSON)` immediately followed by `54: SAVE` (blanket SAVE of
  every local), and lines 284-290 are the same two-branch update as modern 380-390. Modern line
  86 declares `double precision :: deuterium_rate_mid(json), deuterium_rate_mid_start(json)` with
  no `save`, and `grep -n save src/rotation/mid_timestep_model.f90` finds none.
- detail: On the first diffusion sub-step (`first_call`) the routine sets
  `deuterium_rate_mid = drate0 + f*(drate-drate0)`; on every later sub-step it does
  `deuterium_rate_mid_start = deuterium_rate_mid; deuterium_rate_mid = deuterium_rate_mid + f*(...)`
  (lines 385-391), i.e. it depends on the array persisting between calls. As a plain automatic
  local it does not persist: with gfortran's `-finit-local-zero` (in GNU_BASE_FFLAGS) both
  arrays are zero on entry, so sub-step n>=2 passes rate_start = 0 and rate_end = f*(drate-drate0)
  to dburnm (burn_lib.f90:438), grossly under-burning deuterium; without that flag the values
  are stack garbage. Fires only when rotation with angular-momentum transport is on, the model
  still carries deuterium (xa_start(12,nz) > 1e-14), and the diffusion step is split into more
  than one sub-step -- i.e. rotating pre-main-sequence tracks. Observable: D depletion history,
  and the associated luminosity plateau, wrong in rotating pre-MS runs; non-rotating runs and
  post-D-burning runs are unaffected. Fix is to give the two arrays the `save` attribute (or
  move them into rot_scr).

## src/mixing/burn_mix_extrapolated.f90:32-44 -- blanket SAVE dropped; extrapolation state set at order 1 is gone at order >= 2

- class: logical
- severity: low (currently dead code) / high if ever enabled
- confidence: high
- provenance: modernization -- `git show 6cd5673:src/bsrotmix.f | grep -n SAVE` gives
  `12: SAVE` and `13: C SAVE X,D,JJ,LDO,NMAX,LCNO,LCNCHECK`; the modern file keeps only the
  commented line (44) and no `save` attribute anywhere.
- detail: `species_active`, `active_species_id`, `num_active_species`, `step_size_squared`,
  `use_cno_ratio_method`, `he3_extrapolate_log`, `use_cno_ratio_species` are assigned only
  inside `if (extrapolation_order.eq.1)` (lines 46-100) and read unconditionally on later
  calls (103-111, 120, 150, 186-264). With `-finit-local-zero` the order-2 call sees
  `num_active_species = 0`, extrapolates nothing, leaves `max_relative_error` at its 1e-30
  seed, and reports `converged = .true.`, so burn_settle_mix returns the raw 2-sub-step
  composition instead of the Richardson-extrapolated one. The tabulated
  `rot_scr%bs_extrapolation_table` is a module array and is fine; it is only the locals that
  are lost. Not reachable today because evolve_angular_momentum.f90:104 hardcodes
  `burs_extrapolation_active = .false.` (as did getw.f:84), but the routine is silently
  broken for anyone who flips that switch.

## src/mixing/solve_composition.f90:60 -- `composition` declared intent(out) but read for the initial guess and only partially written

- class: interface
- severity: medium
- confidence: high
- provenance: modernization -- `git show 6cd5673:src/kemcom.f | grep -n "HCOMP("` shows HCOMP as a
  plain dummy (line 18) read at 192-197 (`W(1) = HCOMP(1,IBEGIN-1)` ...) and written only for
  the burned species in IBEGIN..IEND; the F77 has no intent, so the read was legal.
- detail: Line 60 is `double precision, intent(out) :: composition(15,json)`. Lines 241-255
  read `composition(1..9, zone_begin-1)` (the previous shell) as the Newton starting point, and
  lines 423-433 write only species 1,2,3,4,5,6,7,9,11 in [zone_begin,zone_end]; species 8,10,
  12-15 and every zone outside the range are never assigned. With intent(out) the actual
  argument (star%xa at every caller) becomes undefined on entry per the standard, so any
  compiler is entitled to skip copy-in or poison the array; in practice the read at 241-255
  works today because gfortran does not exploit this, but `-finit-local-zero`/`-finit-derived`
  style options or a copy-in/copy-out temporary (e.g. if a caller ever passes a section) would
  zero the initial guess and the untouched species. Should be intent(inout).

## src/rotation/am_convective_regions.f90:26 -- intent(in) `convective_flag` passed to a routine that writes element num_zones+1

- class: interface
- severity: low
- confidence: high
- provenance: modernization -- am_convective_regions has no F77 ancestor as a separate file
  (`git show 6cd5673:src | tr ' ' '\n' | grep -i amcz` empty); the original convec.f line 55
  `LC(M+1) = .FALSE.` wrote to a COMMON array with no intent.
- detail: `logical, intent(in) :: convective_flag(json)` (line 26) is passed straight to
  `find_convection_zones`, whose dummy is `intent(inout)` and which executes
  `convective_flag(num_zones+1) = .false.` (find_convection_zones.f90:78). Modifying an
  intent(in) actual through a callee is non-conforming; the compiler cannot see it because
  find_convection_zones is a bare external. Harmless today because the actual is always a
  module array with json > num_zones, but it will bite if a caller ever passes an expression or a
  contiguous section of exactly num_zones elements.

## src/mixing/find_convection_zones.f90:78 and src/mixing/homogenize_convection_zones.f90:58 -- sentinel write at index num_zones+1

- class: logical
- severity: low
- confidence: high
- provenance: inherited -- convec.f:55 `LC(M+1) = .FALSE.`; mixcz.f:42 `LCZ(MP1) = .FALSE.`.
- detail: Both routines write `convective_flag(num_zones+1) = .false.` and then loop
  `do zone_idx = 1, num_zones+1` so that a CZ reaching the surface is closed. This is an
  out-of-bounds write whenever num_zones == json (the array is dimensioned json), and it
  silently clobbers whatever the caller stored at nz+1 (star%lc(nz+1) is also read by the
  rezone code after nz changes). Inherited, but worth an explicit `num_zones < json` guard or
  a local copy since bounds checking will trip it.

## src/mixing/mix.f90:465-467 -- settling sub-step count can be zero and the mod() operands are reversed

- class: numerical
- severity: low
- confidence: medium
- provenance: inherited -- mix.f:413-414 `NSTEP = INT(DELTS/DTMAX)` /
  `IF(MOD(DTMAX,DELTS).NE.0.0D0)NSTEP=NSTEP+1`; rotmix.f:146 has the corrected
  `IF(MOD(DTMAX,DELTS).NE.0.0D0.OR.NSTEP.EQ.0)NSTEP=NSTEP+1`.
- detail: `num_settling_steps = int(timestep/max_settling_dt)` followed by
  `if (mod(max_settling_dt,timestep).ne.0) num_settling_steps = num_settling_steps+1`. The
  test `mod(dtmax, dt)` is the wrong way round (the intent is "does dt divide evenly into
  dtmax-sized pieces", i.e. `mod(dt, dtmax)`); when timestep < max_settling_dt the int() gives
  0 and if max_settling_dt happens to be an exact multiple of timestep the increment is skipped,
  so `substep_dt = timestep/dfloat(0)` -> Inf. rotmix.f90:219-222 carries the `.or. .eq.0` guard
  that fixes exactly this; mix.f90 does not. Practically the exact-multiple condition is rare, so
  the observable consequence is usually just one extra sub-step, but the asymmetry between the
  two copies is a copy-paste divergence.

## src/mixing/semiconvection.f90:106 -- `only_check_core` inspects zone list entry 1, which is not necessarily the core CZ

- class: logical
- severity: low
- confidence: medium
- provenance: inherited -- introduced in F77 commit a8e6067 (sconvec.f `L_ONLY_CORE`,
  2025-12-17), carried over verbatim.
- detail: When `only_check_core` is set the loop upper bound becomes 1, i.e. the first entry of
  `convective_zone_bounds`. That entry is the innermost CZ in the list, which is the core CZ
  only if a convective core exists; for a star with a radiative core and an envelope CZ the
  routine applies its semiconvection test to the envelope. The comment says "core"; the code
  says "first CZ". Only fires when the option is enabled.

## src/rotation/solid_body_omega.f90:89-105 -- iteration tolerance is never initialised, and a non-converged iterate is written back as j (AM not conserved)

- class: numerical
- severity: low
- confidence: high
- provenance: inherited -- solid.f:9 `COMMON/ERRMOM/TOLERI`, 58-59 `IF(DABS(DELJ/CZJ).GT.TOLERI)
  THEN / IF(ICOUNT.LT.20)`, 67-68 `HJM(J) = HI(J)*OMEGA(J)/HS2(J)`; `git grep -n ERRMOM 6cd5673 -- src`
  shows TOLERI only in checkj.f and solid.f, never assigned.
- detail: `rot_scr%moment_of_inertia_tolerance` is declared in rotation_scratch_lib.f90:81 and
  only ever saved/relaxed/restored in check_angular_momentum.f90:187-197
  (`max(convergence_tolerance*1e-2, saved_tolerance)`), so its base value is the default 0.
  With tolerance 0 the test `dabs(delta/J) > 0` is always true and every call performs the full
  20 Newton iterations (cost only; the iterate is converged long before). The real defect is
  the fall-through: after 20 iterations the loop stops regardless of residual and line ~103
  writes `j = I*omega/dm` from the last iterate, so if the Newton step ever stalls (e.g. dI/domega
  large near break-up) the CZ's total angular momentum silently changes. Also NaN if J == 0
  (division by total_angular_momentum). Inherited, but a `stop`/ierr on non-convergence would be
  the right behaviour for a conservation-critical routine.

## src/rotation/evolve_angular_momentum.f90:155 -- envelope boundary read from convective_zone_bounds(num_convective_zones,1) with num_convective_zones possibly 0

- class: logical
- severity: low
- confidence: medium
- provenance: inherited -- getw.f:131 `IMAX = MXZONE(NZONE,1)`.
- detail: find_convection_zones does not record single-shell convective zones, so a model whose
  only convection is a one-zone surface layer has num_convective_zones = 0 and the index
  becomes convective_zone_bounds(0,1), an out-of-bounds read. The value feeds
  envelope_boundary_zone used later for the wind torque coupling. Rare (fully radiative
  rotating star), inherited.

## src/rotation/omega_from_j.f90:60 -- single-shell convection zone at the surface never gets enforce_rotation_profile applied

- class: logical
- severity: low
- confidence: medium
- provenance: inherited -- getrot.f loop `DO 20 I = 1,M` (line 17) with the same in-CZ/exit-CZ
  pattern.
- detail: The zone walk opens a CZ when the flag turns true and closes it (calling
  enforce_rotation_profile) when the flag turns false at i+1. If the flag is true at i = nz and
  false at nz-1 the CZ is opened at nz and the loop ends with no closing call, so omega(nz) is
  left at its old value instead of being recomputed from j. Inherited; consequence is one
  stale surface omega in the (rare) single-shell surface CZ case.

## src/rotation/mid_timestep_model.f90:198 -- comment says new-CZ detection redistributes j; the flag is cleared immediately so the block (203-353) is dead

- class: logical
- severity: low
- confidence: high
- provenance: inherited -- midmod.f:128 `LNEWCZ = .TRUE.` in the detection loop, then 138
  `LNEWCZ = .FALSE.` before 143 `IF(LNEWCZ)THEN`.
- detail: The detection loop sets `new_cz_detected = .true.` when a zone switches convective
  state, and line 198 resets it to `.false.` unconditionally before the `if (new_cz_detected)`
  block that homogenises j across the new CZ. The 150-line block is therefore unreachable in
  both the F77 and the modern code; the header comment ("SAVE CHANGES IN R, CZ DEPTH ... ")
  and the in-block comments describe behaviour that does not happen. Either the reset is a
  deliberate disable (then the block and the comment should go) or it is a leftover debugging
  short-circuit. Reported because comment/code mismatch is in scope.

---

## Weak/uncertain observations

- rotmix.f90:250/255: `if (ierr /= 0) return` inside the settling burn skips the gradT restore; harmless because the driver stops on ierr.
- mix.f90 / rotmix.f90 / mid_timestep_model.f90:398: several F77 `GOTO <outer-label>` became inner-loop `exit`; equivalent only because T is monotonic within a radiative zone list (checked, holds).
- semiconvection merge uses `.gt.` where overshoot_boundaries uses `.ge.` on the same zone-adjacency test; inherited asymmetry, effect is whether two zones touching at one shell merge.
- burn_settle_mix.f90:166-170 restores composition(1..3) from `composition_kept` but leaves star%xa_start(1..3,:) at the last sub-step value (inherited from bursmix.f); dead code (LBURS false).
- mid_timestep_model.f90:401/415: dburnm overwrites composition(1) and (4) from xa_start, discarding whatever H/He3 burning rotmix did in the same sub-step (inherited; second-order in dt).
- viscos.f90: metals counted with weight 1 in the ion number density and electron density taken as <Z> rho/amu; both only enter the Coulomb logarithm (inherited approximations).
- enforce_rotation_profile.f90: modes 2/3 are only reachable when walpcz < 0 and the CZ reaches the surface; for walpcz >= 0 every CZ is forced solid body regardless of mode (matches wczimp.f, but the namelist documentation implies otherwise).
- enforce_rotation_profile.f90:20 duplicated `use star_info_lib`; mid_timestep_model.f90:406-407 empty `if ... then / end if` -- no-ops.
- burn_mix_extrapolated.f90: `converged` is intent(inout) but is always assigned before being read; cosmetic.
- solve_composition.f90:398: the failure message prints zone_begin (comment at 396-397 notes the original printed an undefined `IU`); the "SHELL 1" in the known 0.8 Msun failure therefore means the zone range starting at 1 (core CZ or first radiative shell), not shell 1 alone.

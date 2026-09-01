# Bug sweep: src/rotation/*.f90 + src/rotation/seculr/*.f90 (+ wind torque files named in the task)

Files reviewed in full: am_convective_regions, enforce_rotation_profile, equal_spaced_grid,
evolve_angular_momentum, mid_timestep_model, omega_from_j, rotation_scratch_lib,
solid_body_omega, viscos (rotation/); am_advection_diffusion_coeffs, am_diffusion_coeffs,
am_transport_grid, banded_solver, check_angular_momentum, check_composition,
circulation_velocities, composition_diffusion_coeffs, composition_grid, compute_quadrupole,
diffuse_composition, diffuse_composition_driver, diffusion_velocity_scales,
equal_grid_to_model, secular_transport, zahn_coupling_factor (seculr/). Because the task
explicitly asks for the Kawaler/Matt wind formulas, src/wind/matt_wind.f90 and
src/wind/kawaler_wind.f90 (and the call sites of wind_spindown_matt) were also checked.

## mid_timestep_model.f90:86 -- deuterium-rate arrays lost their SAVE across substeps
- class: numerical (uninitialized read) / logical
- severity: medium
- confidence: high
- provenance: modernization. `git show 6cd5673:src/midmod.f` line 54 is a blanket `SAVE`;
  DRATEM/DRATEM0 are locals that persist between calls. The modernized
  `deuterium_rate_mid(json), deuterium_rate_mid_start(json)` (line 86) are plain automatic
  locals of a module procedure, and src/Makefile has no `-fno-automatic` (only
  `-finit-local-zero`).
- detail: In the `.not.first_call` branch (lines 386-391) the code reads
  `deuterium_rate_mid(i)` carried over from the *previous* call (previous diffusion
  substep) and increments it: original semantics is a running time-interpolation of the
  deuterium burning rate across substeps. With automatic zero-initialized locals, on every
  substep after the first the accumulated rate resets to 0 and becomes
  `step_fraction_ratio*(rate_end-rate_start)` (a tiny fraction of the true rate), and
  `deuterium_rate_mid_start` is 0. Fires whenever `xa_start(12,nz) > 1e-14` (pre-MS
  deuterium burning) with rotation on and more than one diffusion substep; dburnm then
  burns deuterium at a drastically wrong (near-zero, or start-anchored-wrong) rate.

## seculr/diffuse_composition.f90:158-166 -- missing `else i1 = zone_end`; He4 rebalance uses undefined i1
- class: logical
- severity: high (when it fires: mass fractions in the diffused region no longer sum to 1)
- confidence: high that i1 is read undefined; medium-high on practical impact
- provenance: inherited, with modernization-changed failure mode.
  `git show 6cd5673:src/mixcom.f` lines 109-119 has the identical
  `IF(LCZ(IEND).AND.IEND.LT.M)THEN ... ENDIF` with no ELSE, but under a blanket `SAVE`
  (line 12), so I1 held a stale value from a previous region/call; the modernized routine
  dropped the SAVE, so i1 is freshly zero (via -finit-local-zero) on every call.
- detail: The lower-CZ check has an ELSE setting `i0 = zone_begin`; the upper-CZ check has
  none. Whenever the unstable region's upper edge is a radiative zone (common: a region
  ending where the circulation velocity drops to zero mid-radiative-zone) or
  `zone_end == num_zones`, i1 is never assigned. The subsequent
  `do zone_idx = i0, i1` He4 adjustment (lines 211-216, `X(2)=1-X(1)-X(3)-X(4)`) then
  iterates over 1..0 and does nothing, so H/He3 are diffused but He4 is never re-derived
  and the abundances stop summing to unity for that region (in the F77 original the loop
  instead ran over whatever stale I0..I1 range the previous invocation left). The
  conservation-sum locals fed by i1 (sum_species_orig/updated, dcomp2) are dead
  diagnostics, so the He4 line is the only live consumer.

## seculr/am_advection_diffusion_coeffs.f90:249-252,280-283,683-686 -- shear/GSF eq-grid coefficients are never filled anywhere: permanently zero
- class: logical (dead code that was clearly meant to run) / physical
- severity: medium-high in the diffusion+advection (LDIFAD) transport mode
- confidence: high
- provenance: inherited. In the whole 6cd5673 tree, ESS/EGSF appear only in common/DIFAD3
  declarations and in dadcoeft.f lines 187-188/215/550 (in-place multiply and reads); no
  file ever assigns them. Same in the modernized tree: `shear_diffusion_coeff_eqgrid` /
  `gsf_diffusion_coeff_eqgrid` are assigned only by the in-place multiply at line 249.
- detail: circulation_velocities.f90 (lines 524-527) computes the model-grid DSS/DGSF
  (`shear_diffusion_coeff`, `gsf_diffusion_coeff`) explicitly "so we need to multiply the
  original velocity by W ... to cast it as an entry in the diffusion equation", but they
  are never interpolated onto the equal grid (am_transport_grid interpolates only
  VESA/VESD). Since common/module storage zero-initializes, the secular-shear and GSF
  instabilities contribute exactly nothing to the LDIFAD angular-momentum diffusive
  coefficient (line 280) or to the mixing coefficient DCMIX (lines 683-686), despite being
  computed and advertised. Additionally, even if seeded, the in-place update
  `ESS = ESS*(qwr_new/qwr_start)**2` inside the coefficient-iteration loop *compounds*
  the correction each iteration/substep instead of rescaling from a stored baseline.

## seculr/compute_quadrupole.f90:92-95 -- interior d(rho)/dr uses rho(i) where rho(i-1) belongs
- class: physical (copy-paste asymmetry between parallel expressions)
- severity: medium
- confidence: high
- provenance: inherited. `git show 6cd5673:src/getqua.f` lines 56-57:
  `DRHO = 0.5D0*(FPL*EXP(CLN*HD(I+1))+(FMI-FPL)*EXP(CLN*HD(I))-FMI*EXP(CLN*HD(I)))/DR`.
- detail: The weighted centered derivative used for the parallel quantity
  d(rho*omega^2)/dr (lines 101-103) is `0.5*(wp*f(i+1)+(wm-wp)*f(i)-wm*f(i-1))/dr`, but
  the density derivative's last term subtracts `rho(i)` instead of `rho(i-1)`, collapsing
  to `0.5*wp*(rho(i+1)-rho(i))/dr` -- a half-weighted one-sided difference. This enters
  the diagonal (`-4*pi*G*drho_dr/g`) of the tridiagonal quadrupole solve at every interior
  zone, so the Sweet/Zahn quadrupole moment -- and thence the Eddington-Sweet circulation
  correction ratio in circulation_velocities -- is systematically biased whenever
  rotational transport runs.

## seculr/check_angular_momentum.f90:207 -- omega-gradient-reversal detector threshold is 1.0 rad/s
- class: logical (comment/intent vs code mismatch; guard can essentially never fire)
- severity: medium
- confidence: high (that it disables the guard; medium that it was unintended rather than
  a deliberate quiet-disable)
- provenance: inherited. `git show 6cd5673:src/checkj.f` line 150:
  `IF(OMEGA(I)-OMEGA(I-1).GT.1.0d0)` (note lowercase `d0`, suggesting a late edit).
- detail: The primary scan for "POSITIVE OMEGA GRADIENT ENCOUNTERED" triggers only when
  omega increases outward by more than 1.0 rad/s between adjacent zones; every other
  comparison in the same machinery (old-omega checks, the boundary re-checks at lines
  242/257) uses 1.0d-15. Stellar omegas here are ~1e-6..1e-4 rad/s, so the entire
  "SEARCH FOR REVERSAL OF OMEGA GRADIENTS / ENFORCE SOLID-BODY ROTATION IN THE OFFENDING
  ZONES" block (lines 205-288) is effectively unreachable, contradicting the header's
  stated function #2.

## seculr/check_angular_momentum.f90:128-152 -- negative-J timestep cutting is dead; run always stops on first negative J/M
- class: logical
- severity: medium (turns a recoverable condition into a hard stop)
- confidence: high (as comment/code mismatch)
- provenance: inherited. `git show 6cd5673:src/checkj.f` lines 89-92:
  `IREDO = IREDO + 1` / `C IF(IREDO.GT.3)THEN` / `IF(IREDO.GT.0)THEN`.
- detail: cut_count is incremented and then immediately tested `.gt.0`, which is always
  true, so the error branch ("...3 ATTEMPTS AT CUTTING TIMESTEP FAILED / RUN STOPPED")
  fires on the *first* negative specific angular momentum; the else branch that halves dt
  and sets redo_flag (the whole redo/retry path wired through secular_transport and
  evolve_angular_momentum's `retry` loop, including the `elapsed -= 2*dt` rollback) is
  unreachable from this check. The original `IF(IREDO.GT.3)` is left commented above it,
  and the printed message misstates what happened.

## seculr/am_advection_diffusion_coeffs.f90:609-615 -- medium-iteration average averages a value with itself
- class: logical (variable assigned, then wrong one read back)
- severity: low-medium (intended NNN-level damping/averaging is a no-op)
- confidence: high
- provenance: inherited. `git show 6cd5673:src/dadcoeft.f` (block after label 900):
  `EWIT(K) = EWPREV(K)` then `EWIT2(K) = 0.5D0*EWIT(K)+0.5D0*EWPREV(K)`.
- detail: For theta_iter_idx > 2 the code sets
  `omega_prev_medium_iter(k) = omega_working(k)` and on the next line computes
  `omega_prev_medium_iter_avg(k) = 0.5*omega_prev_medium_iter(k) + 0.5*omega_working(k)`,
  i.e. exactly omega_working -- the previous medium iterate has already been overwritten,
  so the branch is identical to the theta_iter_idx <= 2 branch. The intended
  under-relaxation of the omega values used to build the advective/diffusive coefficients
  (omega_mid_it / domega_dr_it at lines 241-247) never happens.

## src/wind/matt_wind.f90:57-201 -- PMM (Matt) wind path ignores wind_loss_active
- class: logical / physical
- severity: medium (removes angular momentum when the wind torque is disabled)
- confidence: medium-high
- provenance: inherited. `git show 6cd5673:src/mwind.f`: the WIND dispatch passes LJDOT and
  wind.f gates its whole body with `IF(LJDOT)`, but the PMM body of mwind.f itself
  contains no LJDOT test (checked lines 30-160).
- detail: kawaler_wind.f90 wraps its spin-down in `if (wind_loss_active)`; matt_wind's PMM
  branch does not, and always subtracts `delta_j_per_mass` from the surface CZ.
  The exposed call site is secular_transport.f90:271-284: the isolated-surface-CZ branch is
  entered whenever `surface_cz_active .and. diffusion_velocity(zone_max)==0` (no gate on
  use_wind_torque); it sets `wind_loss_active = star%job%use_wind_torque` and calls
  matt_wind, so with `use_pmm_wind_law = T` and `use_wind_torque = F` the star is still
  spun down. (evolve_angular_momentum's call sites happen to be gated by use_wind_torque,
  so only the secular_transport path exposes it.)

## mid_timestep_model.f90:198 -- `new_cz_detected = .false.` immediately after the detection loop kills the CZ-boundary J bookkeeping
- class: logical (dead code that was clearly meant to run -- but disabled upstream)
- severity: medium if ever re-enabled; none as shipped
- confidence: high (inherited-deliberate)
- provenance: inherited. `git show 6cd5673:src/midmod.f` lines 128/138: LNEWCZ set .TRUE.
  in the loop, then unconditionally `LNEWCZ = .FALSE.` before `IF(LNEWCZ)THEN`.
- detail: The MHP 06/02 block (lines 203-353) that re-partitions specific angular momentum
  when shells drop out of / into a convection zone can never execute; local J conservation
  is silently assumed across CZ-boundary migration instead. Worth recording because
  (a) it contains raw `write(*,*)` debug spam and a guaranteed-1.0 diagnostic
  (`angular_momentum_check_ratio` computed *after* `star%j_rot(j)` is overwritten,
  line 287-289), (b) if re-enabled it clobbers the intent(inout) dummy `time_fraction`
  (fx) that callers still need (lines 260, 313), and (c) it indexes
  `del_grad_diff_new(change_region_start-1)` = index 0 when a changed region starts at
  zone 1 (line 222/228). All of these are faithful to the original.

## seculr/secular_transport.f90 + evolve_angular_momentum.f90 -- redo/rollback bookkeeping (verified correct, recording to prevent re-flagging)
- class: (non-finding after adversarial check)
- severity: n/a
- confidence: n/a
- provenance: matches 6cd5673 getw.f exactly
- detail: The `elapsed_substep_time - 2.0D0*sub_timestep` rollback in
  evolve_angular_momentum:337 looks wrong in isolation but is correct: check_angular_momentum
  / check_composition halve `dt` (an intent(inout) chain all the way up), so subtracting
  2x the *halved* step removes exactly one original substep, and `fx = 2*sub/full` is the
  pre-halving fraction. Noted here because it is the most tempting false positive in the file.

## Weak/uncertain observations (one line each)
- am_transport_grid.f90:151,182 and composition_grid.f90:82,106: `i0 = i + 1` / `i1 = i - 1`
  use the stale outer loop variable (loop variable is `ii`/`search_idx`) -- inherited
  verbatim from rotgrid.f (`I0 = I + 1` with loop var II); harmless today only because
  i0/i1 are never read in these two files (unlike diffuse_composition, where the same
  pattern's result *is* used).
- viscos.f90:69: `electron_number_density = mean_charge*density/amu` omits the
  sum(X_i/A_i) factor (n_e = rho/amu * sum(X Z/A)), overestimating n_e ~25% for solar mix;
  inherited (viscos.f `ELECN = ZF*RHO/AMU`); only enters the Coulomb logarithm, so effect
  on viscosity is at the percent level.
- check_angular_momentum.f90:104: `max_delta_j_by_iter(16)` indexed by iteration_number up
  to user-set max_diffusion_iters -- out-of-bounds if ITDIF2 > 16 (inherited dimension);
  similarly am_advection_diffusion_coeffs' history arrays are fixed at 50.
- matt_wind.f90:153: `scale_by_b_field` is honored inside the fixed-point loop even when
  `scale_by_rossby_number` is false, in which case `current_turnover_timescale` is read
  uninitialized (inherited shape; original read a SAVEd stale TAUCZP).
- am_advection_diffusion_coeffs.f90:191-204: lrossby/pmmsoltau read as uninitialized
  locals (never wired to the intended COMMON/PMMWIND) -- inherited and already
  self-documented in the file header; with -finit-local-zero the branch now
  deterministically takes the non-Rossby path.
- secular_transport.f90:334: `constant_diffusion_coeff_flag` (LCODM) is never assigned
  anywhere (self-documented dead branch, inherited), so the MHP 8/13 "treat entire domain
  as unstable when a constant D is added" intent is unwired even though
  use_constant_background_diffusion *does* add a constant D in diffusion_velocity_scales
  -- background diffusion therefore only acts inside velocity-detected unstable regions.
- evolve_angular_momentum.f90:201-205: the `int(a/b)` + `mod(a,b).ne.0.0` ceiling idiom on
  doubles (same class as the already-reported rezone finding) -- for generic real inputs
  the mod is essentially never exactly zero, so the step count is int(a/b)+1 nearly always;
  inherited from getw.f.
- mid_timestep_model.f90:457: "surface C.Z. DEEP ENOUGH for angular momentum loss" test is
  `(mass)/Msun .gt. 0.0` -- any nonempty CZ passes; inherited ((HSTOP-HSBOT)/CMSUN.GT.0.0).
- check_composition.f90:96: `if(zone_index.eq.1.or.zone_index.eq.num_zones)` is dead --
  the enclosing loop runs 2..num_zones-1; inherited.
- circulation_velocities.f90: `wmid` (lines 403-409) and `v2` (line 285) are computed and
  never used; the ZM98 mu-inhibition loop (268-291) does not skip convective interfaces
  (both loops that feed it do) -- all inherited verbatim from vcirc.f; harmless because
  VES is zero at convective interfaces.
- equal_spaced_grid.f90 / mid_timestep_model.f90 naming: `log_luminosity*` arrays actually
  hold *linear* L/Lsun (original HL; getw.f does `BL = LOG10(HL(M))`), and the chi
  luminosity term/scale (`luminosity_scale = scale*L_surf`) is correct for linear L --
  misnomer only, but a trap for future edits (a true log L here would flip signs and can
  cross zero).
- solid_body_omega.f90:89: `delta_j/total_angular_momentum` is 0/0 for a region with
  exactly zero total J (non-rotating start); guarded in practice by the LROT gate upstream.
- enforce_rotation_profile.f90:19-20: duplicated `use star_info_lib, only: star` (cosmetic);
  modes 2/3 (solid_body_mode_flag) faithfully match the 2025 F77 wczimp.f, including the
  r(istart)-weighted inertia sum for the sub-CZ region in mode 3 (intent-consistent, not a bug).

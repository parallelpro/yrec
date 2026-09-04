# Bug sweep: core solver / envelope integration / driver layer

Files reviewed in full: henyey_coefficients.f90, henyey_iterate.f90,
henyey_eliminate.f90, henyey_solve.f90, envint_kernel.f90, envint_lib.f90,
envelope_derivs.f90, atmosphere_derivs.f90, surfbc.f90, rebuild_envelope.f90,
shell_physics.f90, evolve_step.f90, run_yrec.f90, stop_conditions.f90
(all under /Applications/YREC/src/core/). Each was compared line-by-line
against its F77 ancestor at 6cd5673 (coefft.f, crrect.f, reduce.f, hsolve.f,
envint.f, qenv.f, qatm.f, surfbc.f, getnewenv.f, physic.f, main.f).
The four structure-equation residual/derivative blocks in
henyey_coefficients (qr/qp/qt/ql and all partials), the central boundary
rows, the full elimination algebra in henyey_eliminate, and the
back-substitution in henyey_solve were re-derived and match both the
original F77 and standard stellar-structure expressions; no defect found
there beyond the inherited quirks listed at the end.

## core/henyey_iterate.f90:391 -- rotation_shape_factors called without its non-optional ierr argument
- class: logical
- severity: high
- confidence: high
- provenance: modernization (commit 212a652 "ierr campaign close-out" added
  `integer, intent(out) :: ierr` as the 11th dummy of
  rotation/shape/rotation_shape_factors.f90:21-23; the call sites in
  core/read_starting_model.f90:1090 and core/evolve_step.f90:128,512 were
  updated, but core/henyey_iterate.f90:391 and setup/rezone.f90:898 still
  pass only 10 actual arguments)
- detail: rotation_shape_factors is a bare external subroutine (not in a
  module), so the 10-argument calls compile silently against the 11-dummy
  signature. The callee unconditionally executes `ierr = 0` (line 61) and
  polint writes ierr again on interpolation failure — every such call
  stores through whatever garbage address occupies the missing 11th
  argument slot on the stack: undefined behavior, silent memory
  corruption, or a crash. It fires on EVERY Newton iteration of every
  rotating run (the call is inside the `if (star%job%rotation_active)`
  block executed after each correction application), and additionally any
  real error detected inside rotation_shape_factors is lost. Fix: pass a
  jerr actual (and check it) at henyey_iterate.f90:391 and rezone.f90:898.

## core/rebuild_envelope.f90:266-269 -- X/Z interpolation sign inverted (extrapolates away from envelope value)
- class: physical
- severity: medium
- confidence: high
- provenance: inherited — getnewenv.f has the identical lines:
  `HCOMP(1,J) = HCOMP(1,M)+FX*(HCOMP(1,M)-ENVX(JM))` and the same for
  HCOMP(3,J), while the surrounding HD/HP/HR/HT lines all use
  `V(M)+FX*(ENV(JM)-V(M))`
- detail: In the species_index==1 branch (interpolating the new fitting
  point between the last interior point and the first envelope point) the
  structure variables interpolate correctly toward the envelope value,
  but X and Z use `comp + f*(comp - env)` — the difference is reversed, so
  the composition is extrapolated AWAY from the envelope value by the same
  fraction it should move toward it. The sibling branch at lines 302-305
  (between two envelope points) has the correct sign, proving the intent.
  Y then inherits the error via the 1-X-Z-He3 closure at line 273. Fires
  whenever mass loss moves the outer fitting point and the first envelope
  point already straddles the new SENV; the error is bounded by the
  interior-vs-envelope composition difference at the fit point (usually
  small, but nonzero with diffusion/accretion), and creates a spurious
  composition kink at the new outer zone.

## core/envelope_derivs.f90:94-95 -- dlogR/dlogP multiplies by fp where it should divide (rotating envelopes)
- class: physical
- severity: medium
- confidence: medium-high
- provenance: inherited — qenv.f:
  `DYDX(1) = -DEXP(CLN*(C4PIL+4.0D0*RL+PL-CGL-SL-SL))/FPL` and
  `DYDX(3) = -DEXP(CLN*(PL+RL-CGL-SL-DL))*FPL`
- detail: In the Kippenhahn-Thomas/Endal-Sofia formalism used by this
  code, dP/dm = -(Gm/4pi r^4)*fp (henyey_coefficients builds exactly that:
  qp is multiplied by fp_rot, and qr — dlogR/dlogm — carries no fp at
  all). With P as the independent variable, both dlogm/dlogP and
  dlogR/dlogP = (dlogR/dlogm)/(dlogP/dlogm) must therefore carry 1/fp.
  dydx(1) correctly divides by pressure_rotation_factor; dydx(3)
  multiplies by it. The two lines are mutually inconsistent whatever fp
  convention is assumed, and inconsistent with the interior equations in
  henyey_coefficients. Relative error in dR/dP is fp^2. Only fires in
  rotating runs (surfbc passes fp=1 otherwise), where fp<1 by up to a few
  percent for fast rotators — a systematic radius-structure error of the
  envelope, hence of the surface boundary condition fed back into the
  Henyey solve.

## core/rebuild_envelope.f90:103,352,358 -- omega() read before set; comment claims a blanket SAVE that does not exist
- class: logical / numerical
- severity: medium
- confidence: high
- provenance: modernization — getnewenv.f has a blanket `SAVE` (so OMEGA
  persisted across calls, however accidentally); the modernized file has
  NO save statement (grep confirms), yet its own comment at lines 92-103
  says "Because of the blanket SAVE below, its value nonetheless persists
  across calls ... Preserved exactly as such."
- detail: omega(json) is now an automatic local. In the solid-body branch
  (line 352, `omega(zone_index) = omega(old_num_zones)`) and the general
  case (line 358, `omega_ref = omega(old_num_zones)*...`) it is read
  before this call ever writes it. Under the project's gfortran flags
  (-finit-local-zero) that read is deterministically 0.0, so every new
  envelope point gets omega=0 and j_rot=0 — silently zero angular momentum
  injected at the surface of a rotating, mass-losing model (the original,
  for better or worse, reused the previous call's omega). Under the
  Makefile's ifort/ifx configurations, which have no local-zeroing flag,
  it is genuinely undefined data. The comment/code mismatch is itself a
  defect: either restore the save (and the documented behavior) or,
  better, seed omega from the caller's star%omega.

## core/envint_kernel.f90:316-337,479-483 -- blanket SAVE lost: eos guess chain and prev_tau/prev_opacity/prev_density read before set
- class: numerical
- severity: low-medium
- confidence: high (that the code deviates from its own documented
  contract and reads locals before set); medium (that it changes physics
  under gfortran)
- provenance: modernization — envint.f has a blanket `SAVE`; the
  modernized integrate_envelope_atmosphere has none, yet its comments
  claim the persistence: line 313-315 "the host scalars keep carrying the
  historical inout guesses across calls", line 476-479 "prev_tau/
  prev_density/prev_opacity ... values entering the FIRST step are
  whatever the previous integration left (a historical quirk ...
  preserved)"
- detail: (a) The prepack at 316-337 loads `pressure`, `log10_density`,
  `density`, `beta`, `ion_fraction`, `specific_heat_cp`, `grada`, etc.
  into eos_res before eos_get; none of these locals is assigned earlier
  in the call, and eos_lib.f90:56-59 documents that the historically-inout
  slots (i_log10_density, i_beta, fxion, i_mu_ion_inv, i_dlnrho_dlnt/p,
  i_cp, i_grada) "carry their previous-call values in" as live eqstat
  state/guesses. Without SAVE they carry stack contents: zeros under
  gfortran's -finit-local-zero (i.e. the original's first-call state on
  every call — different Saha-iteration seeding than the original's
  later calls), undefined garbage under the ifort/ifx flag sets.
  (b) prev_tau/prev_density/prev_opacity are read at line 483 in the
  first atmosphere step before their first assignment (499-500), so
  atmo_delta_depth(1) is computed from zeros/garbage every call, not from
  the previous integration as documented (output-only impact: the
  outermost atmosphere point's geometric depth). Restoring `save` on the
  relay scalars (or explicit seeding) would recover the documented
  behavior and remove the ifort UB.

## core/envint_kernel.f90:801-803 -- env_luminosity inversion is a half-swap (index2 never written)
- class: logical
- severity: low
- confidence: high
- provenance: new-code — the F77 envint.f inversion block (lines 682-731)
  swaps ENVD/ENVP/ENVR/ENVS/ENVT/ENVX/ENVZ/LCENV/EDELS/EVELS/EBETAS only;
  env_luminosity/env_opacity/env_dlnrho_dlnt/gamma1/cp/ion_fraction are
  modern additions
- detail: In the envelope-vector inversion, every other array does the
  three-line swap, but for env_luminosity the code saves
  `swap_temp = env_luminosity(i1)`, assigns
  `env_luminosity(i1) = env_luminosity(i2)`, and then immediately reuses
  swap_temp for env_dlnrho_dlnt without ever writing
  `env_luminosity(inversion_index2) = swap_temp`. Currently masked because
  every stored point receives the identical `luminosity_linear` (lines
  617, 718), so the array is uniform and the missing write changes
  nothing — but the moment env_luminosity ever varies along the envelope,
  the inverted array is corrupted (i2 keeps its pre-swap value). A latent
  copy-paste bug worth fixing while it is still free.

## core/henyey_iterate.f90:160 -- operator precedence: reset_triangle alone forces a surface-BC recompute on any level only via .and., LNEW on every level
- class: logical
- severity: low
- confidence: low (that it differs from intent)
- provenance: inherited — crrect.f:67 is the identical
  `IF(LNEW.OR.LRESET .AND.ITLVL.EQ.2) LSBC = .TRUE.`
- detail: `.and.` binds tighter than `.or.`, so this reads
  `start_new_triangle .or. (reset_triangle .and. iteration_level==2)`.
  The comment two blocks above ("SET UP SURFACE BOUNDARY CONDITIONS - 2ND
  AND 3RD LEVELS OF ITER ONLY") suggests the author may have intended
  `(LNEW .or. LRESET) .and. level==2`; as written, a pending
  start_new_triangle flips recompute_surface_bc to true even on levels 1,
  3, 4. Listed because it survived verbatim from the F77 (so behavior is
  unchanged), but the intent/precedence mismatch is real and worth a
  parenthesization with a decided semantic.

## Weak/uncertain observations
- henyey_iterate.f90:313: `elim_rhs(4,i)/star%luminosity_lsun(i)` divides by a zone luminosity that can be ~0 or cross zero (He-flash, deep interior); IEEE Inf is then discarded by dmin1, but the dev flag set (-ffpe-trap=zero) will trap. Inherited (crrect.f:178).
- surfbc.f90:159-176: the tri_scan re-triangulation loop has no iteration cap and resets tri_err to 0 on every flip; a degenerate/collapsed triangle could in principle cycle forever. Inherited (surfbc.f label 10 loop).
- surfbc.f90:223-225: 1/(vtx_logt(2)-vtx_logt(1)) and the temp2 denominator can be zero for degenerate envelope vertices; no guard. Inherited.
- henyey_coefficients.f90:245 vs 379: hot-zone gate uses `.gt.` cutoff, output gate uses `.lt.` — a zone with logT exactly equal to the cutoff would read energy_gen_component stale from the previous zone (measure-zero; inherited, coefft.f had the same pair plus a blanket SAVE).
- henyey_coefficients.f90:443-444: rot_scr%dlnepsilon_dlnrho/dlnt stored for every zone but only computed in hot zones — cold zones get the previous hot zone's values (first zones of a call: zeros). Inherited quirk of coefft.f's SAVE, now init-zero-dependent.
- henyey_coefficients.f90:405: `star%converged_zone(im) = conductive_opacity_flag` stores a flag that is unconditionally .true. (set at line 127); name is a double misnomer. Inherited (coefft.f LOCONS(IM)=LOCOND, LOCOND=.TRUE.).
- evolve_step.f90:411: luminosity_entropy_term divides by (L+L_start); zero if the zone luminosity crosses -L_start. Inherited (main.f:714 HLL0).
- evolve_step.f90:483: `use_correct_gradients = .true.` is assigned but never used — the original passed it to MIXCZ as IMIX ("so the correct grads are used"); the modern homogenize_convection_zones call takes no such flag, so that behavior toggle may have been silently dropped (needs a look at homogenize_convection_zones, outside this batch).
- evolve_step.f90:505: envelope_cz_zone_prev is passed to evolve_angular_momentum but only assigned when use_extended_composition is true (otherwise it is the SAVEd stale/zero value). Inherited (main.f JENV0 set only under LEXCOM, passed to GETW regardless).
- evolve_step.f90:224 / run_yrec.f90:298: `star%omega(1) .eq. 0` float-equality (vs integer 0) as the rotation-configuration guard; works only because omega comes verbatim from the model file. Inherited pattern.
- henyey_iterate.f90:229: `recompute_surface_bc .and. .not.envelope_recomputed_flag .and. converged` — Fortran does not guarantee short-circuit, and envelope_recomputed_flag is undefined when surfbc was not called (recompute_surface_bc false); value-irrelevant because of the .and., but formally an undefined-variable evaluation the original avoided via crrect's blanket SAVE.
- surfbc.f90:198-199 / rebuild_envelope.f90:163-166: spot branch compares double controls against single-precision literals (`.ne.0.0`, `.ne.1.0`) and passes `4.0`/mixed literals into pow/log10; benign for namelist-read values and inherited verbatim (KNOWN-class float-equality guard, listed for completeness).
- stop_conditions.f90:161: reached_end_age tests `(target-age) .le. 1.0` with no abs(), so it also reports "reached" for any age past the target; matches the original convention (main.f:889) and is presumably intended, but differs from the abs() form used at main.f:518/806.

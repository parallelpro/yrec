# Bug sweep: src/mixing/, src/wind/, src/util/

All files read in full; every finding checked against the F77 original at
6cd5673 (or, where the F77 evolved after that, the last pre-modernization
F77 revision a8e6067). Build context that matters below: `GNU_BASE_FFLAGS`
includes `-finit-local-zero` and does NOT include `-fno-automatic`, so
un-SAVEd locals are automatic and deterministically zeroed at every entry.

## util/compute_timestep.f90:67,129 + util/timestep_limit_heburn.f90:34,111 -- dropped blanket SAVE zeroes the He-burning timestep; likely ROOT CAUSE of the TAHB NaN crash
- class: logical / numerical
- severity: high
- confidence: high (that it is a real bug); medium-high (that it is the TAHB crash root cause)
- provenance: modernization. Original htimer.f line 28 has a blanket
  `SAVE` (and ytime.f line 22 too); modern compute_timestep.f90 declares
  `structure_dt, rotation_dt, helium_dt, ...` (line 67) with no SAVE.
- detail: timestep_limit_heburn's "not a helium flash, core Y below
  threshold" branch (line 111) computes
  `helium_dt = (5.85d17/Lsun)*helium_dt*(M/L)` -- it READS helium_dt on
  the RHS before this call ever assigns it (the file's own header comment
  at lines 29-33 documents the read-before-write and says it "relies on
  whatever value the caller's actual argument already held"). In the F77,
  HTIMER's blanket SAVE meant DELTSY carried the previous model's value
  (a large number, harmlessly discarded by the MIN). In the modern code
  the caller's `helium_dt` is an automatic local zeroed by
  `-finit-local-zero` on every entry, so the branch yields helium_dt = 0
  exactly. compute_timestep then takes
  `previous_timestep = min(..., helium_dt, ...)` = 0 (lines 157/162 --
  note there is a `<1.0d0 -> 1e20` guard for structure_dt and rotation_dt
  but NONE for helium_dt), so the model timestep collapses to 0 and stays
  0 on all subsequent models (the `previous_timestep.gt.1` test then
  routes to line 162, which again min's in the zero).
  When it fires: L_He > 1e-34 (heburn called), log rho_c < 5 (not the
  flash branch), and core_helium_fraction = 1 - Z(core) < atime(1). The
  TAHB inlists set time_core_min (= atime(1)) = 0.002, and on the
  ZAHB->TAHB run 1-Z(core) crosses 0.002 exactly at core-He exhaustion,
  i.e. near the TAHB -- matching "dies near model 950". A zero DELTS then
  produces Inf/NaN downstream in the DELTS divisions (e.g. the
  gravothermal HLL0 = 2(HL-HLO)/((HL+HLO)*DELTS) update in the main
  loop), NaN propagates into the shell luminosity B, and in
  temperature_gradients DELDEL = DELR - DELA = NaN fails the
  `deldel.le.1.0d-6` radiative test (NaN comparisons are false), falls
  into the convective branch, fails `test.gt.0`, and dies with the exact
  observed signature "TPGRAD TRIED TO TAKE SQUARE ROOT OF NEGATIVE
  NUMBER ... DELDEL = NaN". Suggested confirmation: log DELTS around the
  crash, or breakpoint the heburn else branch.

## mixing/burn_mix_extrapolated.f90:32-44 -- dropped blanket SAVE destroys all cross-call extrapolation state; Bulirsch-Stoer extrapolation silently no-ops
- class: logical
- severity: high (whenever the rotational burn/settle/mix path runs)
- confidence: high
- provenance: modernization. Original bsrotmix.f line 12 has blanket
  `SAVE` (the modern file even preserves the original's commented-out
  `SAVE X,D,JJ,LDO,NMAX,LCNO,LCNCHECK` line as a comment at line 44, but
  the active blanket SAVE was dropped).
- detail: `species_active`, `active_species_id`, `num_active_species`,
  `use_cno_ratio_species`, `use_cno_ratio_method`, `he3_extrapolate_log`
  and `step_size_squared` are set only when extrapolation_order == 1 and
  READ on every later call (order 2..11) from burn_settle_mix's loop.
  They are un-SAVEd automatic locals, so each later call sees
  num_active_species = 0 and all-false logicals (zeroed by
  -finit-local-zero; garbage without it). Consequences per order>=2 call:
  every extrapolation loop runs zero iterations, the error scan never
  runs, max_relative_error stays at its 1.0d-30 seed, `converged` is set
  .true. immediately at order 2, and the write-back loop (line 312) never
  copies the extrapolated composition into `composition`. Net effect:
  burn_settle_mix (called from rotation/evolve_angular_momentum.f90:446)
  performs the 2-substep pass, "converges" instantly with a bogus
  MAX ERR 1.0E-30 diagnostic, and returns the UNextrapolated 2-substep
  composition -- the entire Richardson extrapolation and the finer
  substep sequences (4..20) are silently disabled. Also `step_size_squared`
  elements from earlier orders read as 0 in the Neville weights (line 221),
  though that is masked by the empty species loop. The Stage 0 solar
  regression never exercises this path, which is how it slipped through.

## mixing/temperature_gradients.f90:260 -- convective velocity fallback `1.0d0-11` evaluates to -10.0, not 1e-11
- class: numerical (typo)
- severity: low (branch is currently unreachable)
- confidence: high
- provenance: inherited -- tpgrad.f has the identical `VEL = 1.0D0-11`.
- detail: The else branch of the convective-velocity guard assigns
  `convective_velocity = 1.0d0-11`, i.e. 1.0 - 11 = -10 cm/s, where every
  other fallback in the file uses 1.0d-11. Mitigating: the guard
  `test = g*(-dlnrho_dlnt)*presht*deldel .gt. 0` at line 257 is provably
  true whenever this point is reached (the earlier cubic guard at line
  202-203 already required deldel*g*(-dlnrho_dlnt)/presht > 0, and
  presht = P/(rho*g) > 0 always), so the else branch is dead code. Still
  a booby trap for anyone refactoring the guards; worth fixing to 1.0d-11.

## wind/massloss.f90:259-260 -- log10_density_local / beta_local read before ever being set, and the in-file comment falsely claims they are SAVEd
- class: logical (uninitialized read) + comment/code mismatch
- severity: medium-low
- confidence: high (uninitialized read); the practical impact is EOS-seed
  quality only
- provenance: modernization. Original massloss.f had blanket SAVE so DL/
  BETA carried the previous call's converged values (and were static-zero
  on the first call). The modern file has NO save statement (only the
  `data accretion_efficiency/1.0d0/`), yet the comment at lines 118-121
  claims the eos_res slots are "blanket-SAVEd like the locals it
  replaces, so the inout slots keep their historical cross-call carry".
  That claim is false as written.
- detail: At line 259 `eos_res(i_log10_density) = log10_density_local`
  and line 260 `eos_res(i_beta) = beta_local` seed the EOS call with
  locals whose first assignments are at lines 265-266, AFTER the call.
  With -finit-local-zero every call seeds log10(rho) = 0 (rho = 1 g/cc)
  for a photospheric-pressure EOS solve where the true density is
  ~1e-6..1e-8 -- a 6-8 dex bad initial guess instead of the historical
  previous-call carry. Whether that hurts depends on how eqstat uses the
  seed, but the standard-level read-of-undefined plus the false SAVE
  comment are both defects. Fires only on the accretion path
  (use_mass_accretion with mass_loss_rate > 0).

## mixing/mix.f90:465-467 -- settling substep count: mod() arguments transposed, and unlike rotmix there is no nstep==0 guard
- class: logical
- severity: low-medium (division by zero only on exact divisibility;
  otherwise a harmless extra substep)
- confidence: high (that the formula is not what was meant)
- provenance: inherited -- mix.f: `NSTEP = INT(DELTS/DTMAX);
  IF(MOD(DTMAX,DELTS).NE.0.0D0)NSTEP=NSTEP+1`. rotmix.f line 146 has the
  same transposed mod but WITH `.OR.NSTEP.EQ.0` (mirrored at modern
  rotmix.f90:220-222); mix lacks the guard in both eras.
- detail: The intended ceiling of timestep/max_settling_dt requires
  `mod(timestep, max_settling_dt)`; the code tests
  `mod(max_settling_dt, timestep)`. When max_settling_dt >= timestep and
  is an exact multiple of it (e.g. equal after a timestep clamp),
  num_settling_steps = 0 and line 468 divides by zero
  (settling_dt = timestep/0 -> Inf with default non-trapping, feeding
  microdiff/grsett an infinite dt). In the common cases it merely takes
  one more substep than needed. Same idiom family as the already-known
  rezone mod(dp,dp) finding, but this is a distinct live instance in the
  helium-diffusion path of mix. Fix is mod(timestep,max_settling_dt) plus
  the rotmix-style nstep==0 guard.

## util/timestep_limit_heburn.f90:111 -- the shell-burning branch formula is dimensionally incoherent even as originally written
- class: physical
- severity: medium (masked in F77; becomes finding #1 after modernization)
- confidence: medium-high
- provenance: inherited -- ytime.f:86
  `DELTSY = (5.85D17/CLSUN)*DELTSY*(HS1(JXBEG-1)/HL(JXBEG-1))`.
- detail: In the sibling core branch DELTSY is first a Delta-Y budget
  (min(atime(4), atime(5)*Y)) and is then converted to a time by
  (5.85e17/CLSUN)*(M/L). The else branch applies the same conversion
  factor to the PREVIOUS value of DELTSY -- which is a time in seconds,
  not a Delta-Y -- so its output is prev_timestep * O(1e15): dimensional
  nonsense that only "worked" because the resulting huge value was
  discarded by HTIMER's MIN. There is no shell-He Delta-Y criterion here
  at all despite the branch's evident intent. Also note
  `enclosed_mass(h_shell_zone_begin-1)` reads index 0 if
  h_shell_zone_begin == 1. Repairing finding #1 by just SAVE-ing
  helium_dt would restore the F77 behavior, which is itself broken; the
  branch needs an actual formula (e.g. atime-based Delta-Y at the He
  shell) or an explicit 1e20 sentinel.

## mixing/solve_composition.f90:60 vs 241-255 -- composition declared intent(out) but read before written (Newton initial guess)
- class: logical (interface/conformance)
- severity: low (benign with gfortran reference passing; real hazard with
  copy-in/copy-out or -finit treatment of undefined dummies)
- confidence: high
- provenance: modernization (original kemcom.f dummy HCOMP had no intent;
  it is read for the starting guess and partially rewritten).
- detail: `composition` is read at lines 241-255 (previous-shell or
  current-shell abundances as the Newton starting guess) and, on the
  early-return path at line 124-131, zones outside [zone_begin,zone_end]
  are never written at all -- both incompatible with intent(out), which
  leaves the dummy undefined on entry. Should be intent(inout). Same
  class as the already-fixed eqstat2/surfbc/liburn intent bugs, so this
  one appears to have been missed.

## mixing/semiconvection.f90:262-269 -- log_radius_zone never updated in the semiconvective search loop
- class: physical
- severity: low-medium
- confidence: high that the behavior exists; medium that it is a bug
  rather than tolerated approximation
- provenance: inherited -- sconvec.f's DO 32 loop refreshes B,SL,PL,TL,DL
  (and X,Z) per scanned zone but not RL; RL is set once from HR(JMR)
  (itself only added by "MHP 10/02 added definition of RL", before which
  RL was whatever was left over).
- detail: Inside the zone-by-zone stability scan, temperature_gradients
  is called with log_mass/log_pressure/log_temperature of the scanned
  zone but log_radius of the original CZ-adjacent zone, so the local
  gravity g = 10^(cgl+SL-2*RL) mixes the scanned zone's mass with the
  boundary zone's radius. DELR (the actual stability criterion) does not
  use RL, so the Schwarzschild test is unaffected; the inconsistent g
  only enters the MLT cubic/velocity when TL <= tgcut, and can make the
  cubic operate on inconsistent state (contributing spurious
  "sqrt of negative"/non-convergence diagnostics) for semiconvection
  scans reaching cool zones.

## mixing/find_convection_zones.f90:161 vs semiconvection.f90:321 -- merge criteria disagree (.ge. vs .gt.) between the two parallel merge scans
- class: logical (asymmetry between mirrored blocks)
- severity: low
- confidence: high (asymmetry is real and inherited); low that it causes
  damage
- provenance: inherited -- convec.f merges overshoot zones when
  `MXZONE(J,2).GE.MXZONE(J+1,1)`, sconvec.f merges semiconvective zones
  only when strictly `.GT.`.
- detail: After overshoot, two zones that merely touch (top == next
  bottom) are merged; after semiconvective extension the same touching
  configuration is left as two separate "zones" that share a boundary
  shell, and each is subsequently homogenized separately in mix (the
  shared shell gets the second zone's average). Also, both merge scans
  take the merged bounds as [bottom(k), top(k+1)] without max(), so a
  zone fully containing the next would get its top TRUNCATED to the
  contained zone's top -- geometrically impossible for well-ordered
  bounds but not enforced after semiconvective extension.

## wind/wind_spindown_matt.f90:135 -- fsun uses never-assigned `gl` instead of cgl (KNOWN/preserved; sharpened here)
- class: physical (uninitialized variable in formula)
- severity: medium when it fires (Matt-wind spindown path via mcowind)
- confidence: high
- provenance: inherited -- mcowind.f used bare GL where every neighboring
  line uses CGL; explicitly documented as a preserved defect in the
  modern file header, so listing only to quantify: with static-zero F77
  storage or -finit-local-zero, exp(ln10*gl) = 1, so
  fsun = 0.5*w_sun^2*R_sun^3/M_sun MISSES the 1/G factor and is ~1.5e7x
  too big. Since fcen1/fcen2 = ((c_2^2+fsun)/(c_2^2+fcorr))^excen, the
  centrifugal-reduction factors are enormously inflated whenever
  excen /= 0, distorting domega on the use_pmm_wind_law spindown path.
  matt_wind.f90:118 computes the same fsun correctly with cgl, so the
  torque applied by matt_wind and the domega reported by
  wind_spindown_matt are mutually inconsistent.

## mixing/mix.f90:119-133 -- deep-mixing (DPENV) comment says "outer part of the star is mixed"; the code force-mixes the INNER dpenv mass fraction
- class: comment/code mismatch (reportable per instructions)
- severity: low
- confidence: high
- provenance: inherited -- mix.f carries both the original 1980s comment
  ("DPENV = 0.7 means the outer .3 of the star is mixed") and the MHP
  1/95 change; the code (deep_mix_flag(zone) = .true. for
  m(zone) <= dpenv*Mtot, convective flag above) matches the newer inline
  comment "MIX FROM CENTER TO A FIXED MASS FRACTION" and not the header.
  Anyone setting dpenv from the header comment gets the complement of
  what they asked for.

## util/timestep_limit_structure.f90:56,77 -- luminosity "relative change" denominator mixes adjacent zones of the current model instead of old/new values
- class: numerical (inconsistent averaging)
- severity: low
- confidence: medium
- provenance: inherited (ptime.f identical).
- detail: The P/T/R criteria use |old(i) - new(i)|; the L criterion uses
  |old(i)-new(i)| * 2/(new(i)+new(i-1)) -- i.e. the normalization is the
  spatial average of the CURRENT model's zones i and i-1, not the
  temporal average 2/(old(i)+new(i)) used by the exactly analogous
  timestep_limit_omega. Near the center or near L sign changes
  (gravothermal), new(i)+new(i-1) can pass through ~0 and blow up the
  "relative change", spuriously crushing the structure timestep; the
  zero-sum guard only checks .gt.0, so a small negative sum flips the
  sign into a negative test_l that is then ignored by the max -- another
  quiet inconsistency.

## Weak/uncertain observations (one line each)

- find_convection_zones.f90:78 and homogenize_convection_zones.f90:58 write convective_flag(num_zones+1): out of bounds if a model ever reaches num_zones == json (inherited, LC(M+1) same).
- find_convection_zones.f90:91: single-shell convection zones are silently discarded (`zone_start.ne.zone_idx-1`) -- inherited design quirk, means a 1-zone CZ is treated as radiative everywhere downstream.
- semiconvection.f90's dummy `log_luminosity` actually carries LINEAR L/Lsun (HL is linear; verified via physic.f `B=HL(IM)` and main.f's linear HL updates) -- numerics correct, name actively misleading; temperature_gradients' header already documents the mirror-image confusion.
- semiconvection.f90:252-255 uses the eos mu of the PREVIOUS scanned zone with the current radius increment (offset-by-one accumulation in eq. 5' of Castellani et al.) -- inherited.
- burn_mix_extrapolated.f90:196: log(composition(4,i)) can be log(0) = -Inf if He3 hits exactly 0 during substeps after passing the order-1 threshold check -- inherited exposure (guard is evaluated only at order 1) but currently unreachable while finding #2 stands.
- timestep_limit_heburn.f90:75-79/85: h2_fraction (and li6/li7/be9) passed to engeb undefined when use_extended_composition = .false. -- F77 read a stale SAVE'd value, modern reads a zero; effectively benign in a He-flash zone.
- mdot.f90:267-269: shell-mass recompute reads zone_mass_grams(zone_idx-1) = index 0 when envelope_boundary_zone == 1 (fully convective star with mass change) -- inherited (HS1(I-1)).
- timestep_limit_omega.f90:32-36: division by omega(i)+old_omega(i) -- NaN if a zone has omega == 0 in both models (rotating runs normally never do).
- mix.f90:442 checks star%convective_flag(star%nz) for the "no surface CZ" diffusion bail-out while the zone search used deep_mix_flag -- inconsistent when dpenv deep-mixing is active (inherited).
- solve_composition.f90 (kemcom) 50-iteration failure at shell 1 on the 0.8 Msun run (known symptom): no coding defect found in the solver itself; two plausible physical triggers noted -- (a) compute_timestep applies NO He timestep limit until L_He > 1e-34, so pre-flash triple-alpha runaway can hand kemcom an unsolvable (rate*dt) system at the center, and (b) the min_abundance zeroing in the Newton loop (line 385-390) can oscillate: a species clamped to 0 each iteration keeps |correction| >= absolute_tolerance forever. Neither confirmed; instrument iteration_count/abundances at failure to discriminate.
- semiconvection.f90:106 `only_check_core = .true.` hardcodes checking only convection zone #1 (which need not be a core) -- inherited verbatim from the last F77 revision a8e6067 (L_ONLY_CORE), not a modernization change; noted because 6cd5673's sconvec.f still looped over all zones.

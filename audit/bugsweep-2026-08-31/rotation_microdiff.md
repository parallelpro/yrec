# Bug sweep: rotation/microdiff/ + rotation/shape/

Files reviewed in full: gravitational_settling.f90, gravitational_settling_setup.f90,
implicit_diffusion_coeffs.f90, lax_wendroff_step1.f90, lax_wendroff_step2.f90,
microdiff.f90, microdiff_coefficients.f90, microdiff_etm.f90, microdiff_mte.f90,
microdiff_run.f90, microdiff_setup.f90, thoul_diffusion.f90 (microdiff/);
equipotential_integrand.f90, rotation_shape_factors.f90, shape.f90,
shell_inertia_integral.f90, zone_moments_of_inertia.f90 (shape/).
Every file was diffed against its F77 original at 6cd5673 (setup_grsett.f, grsett.f,
get_imp_diffco.f, lax_wendrof1.f, lax_wendrof2.f, microdiff.f, microdiff_cod.f,
microdiff_etm.f, microdiff_mte.f, microdiff_run.f, microdiff_setup.f, thdiff.f,
func.f, fpft.f, shape.f, intmom.f, momi.f), including a mechanical numeric-literal
comparison (no transcription typos found). The thoul_diffusion.f90 Burgers-equation
assembly was checked line-by-line against thdiff.f (= Anne Thoul's published
diffusion.f) and against the Thoul, Bahcall & Loeb 1994 coefficient structure
(0.6/1.5/1.6/2.7/0.8 x_ij,y_ij,yy_ij couplings, z/A constraint rows, K0=2 p-scaling):
it is a faithful port; no assembly defect found there.

## rotation/microdiff/gravitational_settling_setup.f90:379-382 -- Noerdlinger dlnLambda/dX correction has the wrong sign
- class: physical
- severity: medium
- confidence: high
- provenance: inherited (6cd5673 setup_grsett.f:276-279, identical `QCOD1X(I)=QCOD1X(I) + COD1(I)*1.5D0/(LN_LAMBDA*(3.0D0*X+1.0D0))`)
- detail: For coulomb_log_choice==2, ln_lambda = -19.7 - 0.5*log(1/(0.75X+0.25)) + ...,
  so d(lnLambda)/dX = +1.5/(3X+1) (lnLambda increases with X; verified by hand from
  the code's own line 235). Since D1 is proportional to 1/lnLambda, the extra term in
  d(D1)/dX must be **-**D1*1.5/(lnLambda*(3X+1)); the code **adds** it with a plus
  sign, i.e. the correction is applied with inverted sign (error = 2x the term).
  Fires whenever Noerdlinger's Coulomb log is selected; corrupts the diffusion-
  coefficient derivative used to update D1 at midpoints in gravitational_settling.f90:177-179
  and hence the accuracy of the settled X profile (a few-percent-level term).

## rotation/microdiff/gravitational_settling.f90:222-226 -- Iteration midpoint deltas are shifted one cell w.r.t. the coefficient arrays (and element 1 is stale)
- class: logical
- severity: medium
- confidence: high
- provenance: inherited (6cd5673 grsett.f:159-161 `DO 30 I = 2,NPT / EX_H(I)=...` vs get_imp_diffco.f:12-14 reading `EX_H(1..NPT-1)`)
- detail: Throughout the pipeline "mid" index i means the midpoint between zones i and
  i+1 (lax_wendroff_step1 fills 1..npt-1 that way; implicit_diffusion_coeffs.f90:31-34
  corrects diffusion_coeff_mid(i), i=1..npt-1, which multiplies the x(i+1)-x(i) flux).
  But the implicit-iteration loop stores the change at the midpoint between i-1 and i
  at index i (`do eq_idx = 2,num_equal_points`), so the D2 correction
  dD2/dX*deltaX_mid is evaluated one cell inward of the coefficient it corrects.
  Worse, equal_hydrogen_fraction_mid(1) is never rewritten inside the loop: it still
  holds the Lax-Wendroff step-1 provisional change, and since implicit_diffusion_coeffs
  *accumulates* (`coeff_mid += deriv*delta`) with a delta that should shrink to zero as
  the iteration converges, the innermost midpoint coefficient is instead incremented by
  the same nonzero amount on every iteration -- a systematic drift of D2(1) proportional
  to the iteration count. Fires on every grsett settling step. (In microdiff_run.f90 the
  same shifted loop exists at lines 170-174 but is harmless there because
  diffusion_coeff2_deriv_mid is zeroed each iteration, lines 179-182.)

## rotation/microdiff/gravitational_settling.f90:249-257 -- Convergence test ignores negative abundance corrections (no abs)
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited (6cd5673 grsett.f:179-187 `DX = EX(I)-EX_P(I); IF(DX.GT.DXMAX)...`; same in microdiff_run.f:137-145)
- detail: max_delta_x is the maximum of the *signed* per-zone change; any zone whose
  hydrogen abundance decreased by a large amount contributes nothing, so the loop can
  declare convergence (max_delta_x < settling_tolerance) after one iteration when the
  corrections are predominantly negative, leaving the implicit D2 solve effectively
  un-iterated. Should be dabs(delta_x_local). Identical defect in
  microdiff_run.f90:199-207 for the new pipeline (both inherited verbatim). Impact:
  silently under-converged settling solutions rather than a crash, so it is easy to miss.

## rotation/microdiff/microdiff_coefficients.f90:164-169 -- fgrz is applied to light-element (Li/Be) diffusion, not just metals
- class: physical
- severity: medium
- confidence: medium
- provenance: new-code (comment "JvS 01/26"; 6cd5673 microdiff_cod.f:112-113 has plain `FAC=HRU_I**2*HTU_I**2.5D0/LN_LAMBDA` with the FGRLI variant commented out)
- detail: The 2026 addition scales fac by star%ctrl%fgry when species_col==1 and by
  star%job%fgrz when species_col==3. But microdiff.f90 passes species_col=3 for BOTH
  heavy-metal diffusion (line 185) and each light element Li6/Li7/Be9 (line 214), so a
  user who tunes the metal-settling factor fgrz away from 1 silently rescales lithium
  and beryllium settling by the same factor. The F77 original applied no factor to any
  microdiff species (a per-light-element FGRLI array existed but was commented out).
  If the intent was grsett-parity (FGRY/FGRZ scale He/metal settling), light elements
  should keep fac unscaled or get their own factor.

## rotation/microdiff/microdiff_coefficients.f90:113-118,184-189 -- Uninitialized coeff_scale/pressure_term/... read after the zero-species `cycle`
- class: numerical
- severity: medium
- confidence: high (uninitialized read is certain; bad consequence is probabilistic)
- provenance: modernization-degraded (6cd5673 microdiff_cod.f:68-73 has the same GOTO 5 skip, but all its locals are under `SAVE`, so the second loop read static, finite, stale values; the modern arrays are automatics containing stack garbage)
- detail: When species_fraction(species_col,i)==0.0 at i and i+1, the first loop sets
  only diffusion_term(i)=0 and cycles, leaving coeff_scale(i), pressure_term(i),
  temp_term(i), hydrogen_term(i) unset. The second loop then computes
  diffusion_coeff1(i) = fac*(dlnp_dr*(ap+at)+dlncdr*ah)*0.0 and
  diffusion_coeff2(i) = fac*0.0 from those unset values -- 0*garbage, which is 0
  unless the stack garbage happens to be Inf/NaN bit patterns, in which case NaNs
  enter the tridiagonal solve. Fires on every lithium/beryllium diffusion call for
  models with depleted zones (exactly the case the guard exists for). The `.eq.0.0`
  float-equality guard itself is a KNOWN open audit item ("microdiff species .eq. 0.0")
  and is not re-reported here; the uninitialized-read consequence is the new sharpening.
  Cheap fix: initialize the five stash arrays (or cycle in the second loop too).

## rotation/shape/shell_inertia_integral.f90:46-51 -- "dI/d(omega)" is dimensionally not dI/domega (missing 2/omega factor)
- class: physical
- severity: low
- confidence: medium
- provenance: inherited (6cd5673 intmom.f:15-21 identical SUM2; solid.f:61 uses `DELJ/(CZI + WGUESS*CZQIW)`)
- detail: I_per_mass = H1*(1+0.15*sum_j c_j*(j*eta+5)*a^j) with a proportional to
  omega^2, so dI/domega = H1*0.15*sum_j j*c_j*(j*eta+5)*a^j * (2/omega). The routine
  returns that sum WITHOUT the 2/omega factor, i.e. QIWM = (omega/2)*dI/domega, which
  has units of I, not I/omega. solid_body_omega.f90 then forms the Newton slope
  dJ/domega = I + omega*QIW, which should be I + omega*dI/domega = I + 2*QIW; the
  correction term is off by omega/2 (dimensionally inhomogeneous as written). Because
  the term is a small rotational correction and the omega iteration checks a
  J-residual, this only degrades Newton convergence rather than the converged answer --
  which is presumably why it survived since Law's thesis.

## rotation/shape/shape.f90:28,148-150 -- r0 declared intent(out) but read at r0(zone_start-1) before being set
- class: logical
- severity: low
- confidence: high (conformance violation certain; misbehavior compiler-dependent)
- provenance: modernization (F77 shape.f had no intents; R0(JSTART-1) legitimately carried the caller's previous value in)
- detail: When shape is called with zone_start>1 (solid_body_omega.f90:71 does this for
  every internal solid-body CZ), the eta2/r0 correction branch computes
  `dr = r0(i) - r0(i-1)` and `deta1 = .../r0(i-1)` with i = zone_start, reading
  r0(zone_start-1), which the routine has not written this call. With intent(out) that
  element is formally undefined on entry; it only works because gfortran does not
  scrub intent(out) explicit-shape arrays. eta2 was correctly declared intent(inout)
  for the same access pattern one line earlier -- r0 should be intent(inout) too.
  Latent divide-by-garbage if built with a compiler/sanitizer that poisons intent(out).

## rotation/microdiff/microdiff_etm.f90:101-104 -- Loop-exhaust fallback k0=num_eq_points-3 lacks the k0>=1 guard (OOB for npt<=3)
- class: logical
- severity: low
- confidence: high (defect), low (probability of firing)
- provenance: inherited (6cd5673 microdiff_etm.f:66-67 `K0 = NPT-3` unguarded; the JVS "fix for NPT = 3?" guard at line 61 exists only on the in-loop path)
- detail: The equally-spaced grid can legally have num_eq_points = 2 or 3 (the JVS
  NPT=1 trap in microdiff_mte.f90:80 makes 2 the minimum). If the inner search loop
  exhausts (radius rounding at the outer edge), k0 becomes -1 or 0 and
  tabler(k)=eq_radius(k0+k-1) indexes element 0 or -1 -- an array-bounds violation.
  The in-loop branch got a `if (k0 .eq. 0) k0=1` fix; the exhaust branch (a copy of
  the same logic) never did -- a classic asymmetric-fix copy-paste bug.

## rotation/microdiff/lax_wendroff_step1.f90:41-43 -- "Central" boundary condition assumes the diffusion region starts at m=0
- class: physical
- severity: low
- confidence: medium
- provenance: inherited (6cd5673 lax_wendrof1.f:19-21 `EMASS=EM(2); DEX_H=DT2*ECOD1(2)/EMASS`; lax_wendrof2.f:17-19 likewise)
- detail: The inner-boundary update uses diffusion_coeff1(2)/eq_mass(2), i.e. it
  implicitly sets COD1(1)=0 (zero flux) AND uses the *enclosed* mass eq_mass(2) as the
  zone mass. That is only right when the grid starts at the stellar center. When
  zone_begin>1 (diffusion region sits above a convective or H-exhausted core),
  eq_mass(1) is large and the zone mass should be eq_mass(2)-eq_mass(1); using
  eq_mass(2) suppresses the boundary update by orders of magnitude (and COD1(1) is
  not actually zero there). Same structure in lax_wendroff_step2.f90:39-41 with
  eq_mass_mid(1). Affects both grsett and microdiff pipelines for models with
  convective cores; the alpha(1) term of the implicit solve
  (gravitational_settling.f90:212, microdiff_run.f90:161) shares the assumption.

## rotation/microdiff/gravitational_settling_setup.f90:339-344 -- constant_mixing_coeff applied to D1/D2 but not to their X-derivatives
- class: physical
- severity: low
- confidence: high (asymmetry is certain; "bug vs. accepted sloppiness" is the question)
- provenance: inherited (6cd5673 setup_grsett.f:241-244: `COD1=CSTMIXING*FAC*...` but `QCOD1X(I) = FAC*HQPR(I)*(AP+AT+X*(QAPX+QATX))` without CSTMIXING; metal block 264-269 likewise applies CSTDIFFMIX*CSTMIXING to COD1Z but only CSTMIXING to QCOD1Z)
- detail: In the Thoul branch, diffusion_coeff1/diffusion_coeff2 are scaled by
  star%ctrl%constant_mixing_coeff (and the metal coeff1 additionally by
  constant_settling_reduction), but diffusion_coeff1_dx/diffusion_coeff2_dx (and
  src_grid_metal_diffusion_coeff1_dz) are not, so d(D)/dX is inconsistent with D
  whenever those mimic-mixing knobs differ from 1. The derivative arrays feed the
  midpoint-coefficient corrections in gravitational_settling.f90:177-187 and the
  implicit iteration, so the correction terms are over-weighted by 1/CSTMIXING.
  The "old ver" comments show the factors were bolted onto the coefficients (CFD 10/09)
  without touching the derivatives.

## rotation/microdiff/gravitational_settling_setup.f90:226 vs 256 -- Comment says DEL*6*(X-0.32), code implements (X+0.32)
- class: physical (comment/code mismatch)
- severity: low
- confidence: low (that the code half is the wrong half)
- provenance: inherited (6cd5673 setup_grsett.f:140 comment `(X-0.32)`, code line 168-169 `(X+0.32D0)`)
- detail: The header derivation of the Bahcall & Loeb D1 states the thermal-diffusion
  bracket as 5/4 + DEL*6*(X-0.32)/(5.4+6.3X-4.5X^2), but the code uses (X+0.32) --
  and the hand-verified analytic derivative on lines 259-264
  ((3.384+2.88X+4.5X^2)/FAC2^2) is consistent with the **+0.32** form, so the code is
  at least self-consistent. Either the comment or the code misquotes the B&L (1990)
  fit; if the paper's sign is the comment's, D1's thermal term is wrong for all
  non-Thoul (use_thoul_diffusion=false) settling runs. Flagged per the
  comment-vs-code rule; deciding it needs the B&L 1990 expression.

## rotation/numerics/numerics_lib.f90:1729-1732 (trapzd, drives rotation_shape_factors.f90) -- Midpoint interpolation evaluates rho/w2/eta22 at a fixed offset, not at y
- class: numerical
- severity: low
- confidence: high
- provenance: inherited (6cd5673 trapzd.f:38-41 identical `RHOT=RHOP+DRHO*DEL` inside the J loop)
- detail: (Out of assigned dirs but it is the integrand engine of the assigned
  fpft/rotation_shape_factors loop, so recorded here.) In the refinement branch,
  rhot/w2t/eta22t are computed as endpoint_prev + slope*DEL -- constant for all it
  midpoints -- instead of + slope*(y-b1); only smt tracks y. At the n=2 level (the
  only refinement fpft ever reaches, since jmax=2 -- inherited from fpft.f:24) the
  single midpoint at y=b1+dr/2 gets rho of the *outer endpoint* instead of the
  midpoint value. Biases the distortion integral AINT and thence fp/ft; bounded by
  the shell-to-shell variation of rho, so small per shell but systematic.

## rotation/microdiff/gravitational_settling.f90:280-288 -- z_change_first/z_change_last computed and never used
- class: logical (dead code that was meant to do something)
- severity: low
- confidence: medium
- provenance: inherited (6cd5673 grsett.f:206-207: ZZ1, ZZ2 assigned, never read)
- detail: Before converting the metal run to changes, the code computes the change at
  the first and last equally spaced points into z_change_first/z_change_last and then
  never reads them. In equal_to_model the boundary extrapolation of the metal change
  uses star%metal_abundance_change(1)/(npt) directly, so these look like leftovers of
  an intended boundary hand-off (the hydrogen path has no analogue). Harmless today,
  but it is the kind of "assigned but the wrong thing read back" residue the sweep
  asks to record; deletion or wiring-up should be a conscious choice.

## Weak/uncertain observations
- microdiff_etm.f90:63,123,153 and the etm original: zz2 = zz/composition(3,i) divides by Z with no guard; Z==0 anywhere gives NaN scaling of species 5-11 (inherited).
- microdiff_etm.f90:52-55 applies eq_delta_hydrogen unconditionally (no diffuse_helium_active guard); currently safe only because both call sites (mix.f90:423, rotmix.f90:156) gate the whole microdiff call on that flag -- defensive gap, not a live bug. Same for eq_delta_hydrogen/eq_delta_metal/eq_delta_light being automatics in microdiff.f90 (F77 originals were SAVE'd statics).
- microdiff_run.f90/microdiff_coefficients.f90: hydrogen_dlnc_dr is only computed on the species_col==1 call; a hypothetical metal-only invocation would read it uninitialized (currently unreachable, same gating as above).
- microdiff_mte.f90:78 `mod(drtot,drmin).ne.0.0d0` ceiling idiom on doubles -- same family as the KNOWN rezone mod(dp,dp) finding; also `half_json = 5000` is a hardcoded copy of json (original HALFJSON=JSON=5000), silently wrong if json ever changes.
- microdiff_mte.f90:144-147/239-242 loop-exhaust fallback sets k0 = num_zones-3 (outermost model points) rather than zone_end-3; if it ever fires with a deep surface CZ it extrapolates from points far above the diffusion region (inherited, float-edge only).
- microdiff_coefficients.f90:207 `dlncdr = 1.0` immediately overwritten in every branch (inherited dead assignment); single-precision literals 0.0/1.0/0.5/2.2 are byte-pinned intentionally per the file header (not re-reported).
- microdiff_coefficients.f90:213-219 one-sided stencil denominators both use dradi = r(i+1)-r(i); correct only because the grid is equally spaced (inherited).
- gravitational_settling_setup.f90:333-336: when the full Thoul solve is used (use_thoul_fit=false), the X-derivatives dap_dx/dat_dx/dac_dx still come from the *fit* polynomials -- an inherited approximation, not a transcription error (setup_grsett.f:235-238).
- gravitational_settling_setup.f90:82: atomic_weight He=4.004 (true 4.0026), Fe=55.86 (55.845) -- inherited, sub-0.1% effect.
- thoul_diffusion.f90:65,123 mutates the caller's mass_fraction(num_species) (electron mass fraction side effect) -- inherited from thdiff.f X(M)=A(M)/AC; callers pass scratch arrays so currently benign.
- rotation_shape_factors.f90:68 jmax=2 caps Richardson refinement at one extrapolation with no convergence fallback (dint test is vacuous) -- inherited design (fpft.f:24), flagged only because eps=1d-6 suggests a real convergence loop was intended; aintt (fpft AINTT) is assigned and never used (inherited dead).
- equipotential_integrand.f90:47-53: comments quote 12piG/5R^4 and 4piG/5R^3 with opposite signs and a factor G, while the code uses 4pi/3/r^4 and opposite signs -- inherited verbatim from func.f; consistent with how trapzd builds AINT (which absorbs G and the 3/5 bookkeeping), so treated as an inherited comment/code tension, not a code bug.
- zone_moments_of_inertia.f90:61-63 central-zone I uses 0.4*m*r^2 with r = geometric mean of the log radii of zones 1-2 (inherited momi.f:34-35); prev_log_mean_radius/prev_log_true_radius and spherical_moment_of_inertia are computed but never read (inherited dead, already noted in file comments).

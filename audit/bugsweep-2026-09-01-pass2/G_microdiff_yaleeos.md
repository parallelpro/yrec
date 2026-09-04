# Sweep G (pass 2): microscopic diffusion / gravitational settling, model<->equal-grid I/O, Yale analytic EOS

Read-only audit of branch yrec-modern, 2026-09-01. F77 ancestors taken from
`git show 6cd5673:src/<name>.f`.

## Files read in full

- rotation/microdiff/microdiff.f90 (248) -- ancestor microdiff.f
- rotation/microdiff/microdiff_setup.f90 (181) -- microdiff_setup.f
- rotation/microdiff/microdiff_mte.f90 (303) -- microdiff_mte.f
- rotation/microdiff/microdiff_etm.f90 (184) -- microdiff_etm.f
- rotation/microdiff/microdiff_run.f90 (233) -- microdiff_run.f
- rotation/microdiff/microdiff_coefficients.f90 (235) -- microdiff_cod.f
- rotation/microdiff/gravitational_settling.f90 (297) -- grsett.f
- rotation/microdiff/gravitational_settling_setup.f90 (385) -- setup_grsett.f
- rotation/microdiff/lax_wendroff_step1.f90 (72) -- lax_wendrof1.f
- rotation/microdiff/lax_wendroff_step2.f90 (71) -- lax_wendrof2.f
- rotation/microdiff/implicit_diffusion_coeffs.f90 (51) -- get_imp_diffco.f
- rotation/microdiff/thoul_diffusion.f90 (323) -- thdiff.f
- io/model_to_equal.f90 (333) -- model_to_equal.f
- io/equal_to_model.f90 (162) -- equal_to_model.f
- eos/yale/fully_ionized_eos.f90 (325) -- eqrelv.f
- eos/yale/saha_eos.f90 (483) -- eqsaha.f
- eos/yale_eos_lib.f90 (62) -- new (state module)
- eos/mu.f90 (61) -- mu.f
- eos/eos_mixture_lib.f90 (32) -- new (state module)

Also consulted for call-site / interface checks: mixing/mix.f90 (calls at
473/479), mixing/rotmix.f90 (247/252), numerics/numerics_lib.f90
(tridiag_gs, ludcmp, lubksb, interp, intrp2), eos/eqstat.f90 (callers of
fully_ionized_eos / saha_eos), src/makefile (compiler flags).

## Verified clean (re-derived and found correct)

- Thoul et al. (1994) Burgers-equation matrix in thoul_diffusion.f90:
  xx/y/yy/k definitions, gamma (concentration-gradient RHS) rows, the
  momentum rows (delta(i,j), heat-flux couplings 0.6*xx / -0.6*y), the
  energy rows (1.5*xx, -y*k*(1.6*xx+yy), -0.8*k(i,i), 2.7*k*xx*y), the
  charge-neutrality and zero-net-mass-flux constraint rows, LU solve for
  the three RHS families, and the final p/K_0 scaling ko*ac*cc -- all
  match Anne Thoul's published diffusion.f line for line and thdiff.f.
- Coulomb-log / Debye-length / inter-ionic-distance formulas in
  microdiff_coefficients.f90 and gravitational_settling_setup.f90
  (ne = rho/(m_u*ac), ao = (3/(4 pi n_i))^(1/3) = (0.23873/ni)^(1/3),
  lambda_D = 6.9010*sqrt(T/(ne*cz)), xij = 2.3939e3*T*lambda/|ZiZj|,
  cl = 0.81245*ln(1+0.18769*xij^1.2)) -- Thoul's form, coefficients
  consistent with CGS constants.
- Sign conventions ap=-A_p, at=-A_T*gradT, ad=-A_c(j,j) are used
  identically in microdiff_coefficients and gravitational_settling_setup,
  and agree with the sign of the hard-coded Thoul fits (H rises, Fe sinks;
  self-diffusion coefficient positive as the implicit scheme requires).
- Bahcall-Loeb unit conversions (radius, mass = 1e-2*rad^3, temperature
  1e-7, time 2.7e13 yr = 6e13 yr / 2.2) are applied in setup and exactly
  reversed in etm / equal_to_model; the constant-2.2 prefactor combined
  with the 6e13/2.2 time unit reproduces Thoul's tau_0 = 6e13 yr when the
  Thoul solver is given real Coulomb logs (coulomb_log_choice=4 and the
  microdiff path). See finding 1 for the case that does NOT cancel.
- Two-step Lax-Wendroff: step 1 (half step to midpoints from
  zone-centre fluxes) and step 2 (full step to centres from midpoint
  fluxes) have consistent centre (zero flux at r=0) and surface (zero
  flux at M) boundary treatment, and are identical to lax_wendrof1.f /
  lax_wendrof2.f. There is no CFL check, but the explicit term is an
  advective settling velocity whose crossing time of one grid cell is
  many orders of magnitude longer than any evolutionary timestep, so the
  explicit half-step is stable in practice; the second-derivative term is
  fully implicit (unconditionally stable).
- Tridiagonal matrix assembly in implicit_diffusion_coeffs.f90
  (sub/diag/super = -a*D(i-1), 1+a*(D(i-1)+D(i)), -a*D(i)) with alpha =
  4 pi dt/(dr*dm) is a correct conservative implicit discretization of
  4 pi d/dm (D dX/dr); tridiag_gs is a standard Thomas solve.
- fully_ionized_eos.f90: chi_T (qpt), chi_rho (qpd), the second
  derivatives qqptt/qqpdt/qqpdd, the density derivatives
  qdp=1/chi_rho, qdt=-chi_T/chi_rho and their T/P derivatives, the
  internal-energy derivative qudt = (P/rho u)(1-chi_T) (thermodynamic
  identity), c_P = (P delta/(rho T)) chi_T + (du/dT)_rho, grad_ad =
  P delta/(rho T c_P), and all the c_P / grad_ad log-derivatives were
  re-derived and are thermodynamically consistent for the NON-Debye-
  Huckel parts. The Fermi-table indexing (X = log rho - log mu_e -
  1.5 log T, y = log T on 0.2-dex nodes) and the d/dlnT = d/dy - 1.5 d/dx
  chain rule are correct. Byte-identical logic to eqrelv.f.
- saha_eos.f90: the K-loop/GOTO elimination (exit on K < -tol, 1e16 +
  nz1 bump on K > tol), the He I/He II three-way GOTO ladder (labels
  13/14/15), the block-elimination guards (labels 200/210/100/105/110/
  115) and the derivative sections were compared statement by statement
  with eqsaha.f and match. c33 < 0, c11 > 0, c22 > 0 always, so no
  division by zero; E is clamped to [1e-11, 1+Y]. The DATA-initialized
  iterates (E, He ion fractions) are implicitly SAVE, reproducing the
  F77 warm-start behaviour.
- mu.f90 identical to mu.f; call sites in oeqos/oeqos01/oeqos06/meqos
  pass 9 arguments matching the dummy list.
- Interface checks: microdiff / gravitational_settling call sites in
  mix.f90 and rotmix.f90 (11 actuals, matching dummies, star%m passed
  intent(inout) and restored); lax_wendroff_step1/2,
  implicit_diffusion_coeffs, tridiag_gs, thoul_diffusion, ludcmp/lubksb,
  model_to_equal, equal_to_model, fully_ionized_eos, saha_eos -- all
  actual lists match the dummy lists in count, order and type.

---

## rotation/microdiff/gravitational_settling_setup.f90:235-245,316-321,350 -- Coulomb logarithm counted twice when coulomb_log_choice = 2 or 3 is combined with the full Thoul solver
- class: physical
- severity: medium
- confidence: medium
- provenance: inherited -- setup_grsett.f lines 146-159 (LN_LAMBDA from Noerdlinger/Loeb, FAC = ...**2.5/LN_LAMBDA), 221-228 (TCL(II,JJ) = LN_LAMBDA for ILAMBDA /= 4), 250 (metal FAC /LN_LAMBDA)
- detail: The Thoul solver returns A_p, A_T, A_c that already scale as 1/lnLambda because every k(i,j) is proportional to coulomb_log(i,j). The code's unit system absorbs a *constant* 2.2 into the time unit (bl_time_scale = 6e13 yr / 2.2, comment "INCLUDES FACTOR OF 2.2 FROM LN LAMBDA") and then divides the prefactor by ln_lambda, so the physical settling rate is r^2 T^2.5 A * (2.2/ln_lambda) / 6e13 yr. That is correct when ln_lambda is the constant 2.2 (choice 1), when the Thoul fit is used (the fit is rescaled from Thoul's tabulated lnLambda), and when choice 4 feeds real Coulomb logs to Thoul while ln_lambda stays 2.2. But with use_thoul_diffusion=.true., use_thoul_fit=.false. and coulomb_log_choice=2 or 3, the Noerdlinger/Loeb lnLambda is passed into every coulomb_log(i,j) AND the prefactor is divided by the same lnLambda, giving a rate proportional to 1/lnLambda^2 * 2.2 -- an error factor of 2.2/lnLambda (roughly 0.9-1.1 in the solar interior, but 2.2/lnLambda ~ 0.5 or worse where lnLambda reaches 4-5 in hotter, less dense low-mass-star envelopes). The same double count is applied to the metal prefactor at line 350. Observable: the H and Z settling rate in that configuration is biased by 2.2/lnLambda relative to the intended Thoul-with-variable-lnLambda physics.

## rotation/microdiff/microdiff_coefficients.f90:113-118,184-189 -- coeff_scale/pressure_term/temp_term/hydrogen_term read uninitialized for depleted zones after the `cycle`
- class: numerical
- severity: medium
- confidence: high
- provenance: modernization -- microdiff_cod.f lines 25 (SAVE), 68-71 (GOTO 5 after setting only ADS(I)=0.0), 129-133 (FAC=FACS(I), AP=APS(I), ...); the F77 arrays FACS/APS/ATS/AHS were SAVEd statics, the f90 arrays are automatic locals
- detail: When the diffused species has zero abundance at point i and i+1, the first loop sets only diffusion_term(i)=0 and `cycle`s, leaving coeff_scale(i), pressure_term(i), temp_term(i), hydrogen_term(i) unassigned; the second loop then reads all four (fac = coeff_scale(i) etc.) and forms diffusion_coeff1(i) = fac*(...)*X and diffusion_coeff2(i) = fac*ad. In the F77 the arrays were SAVEd so the reads returned finite stale values from the previous species call and the exact-zero multipliers made the products 0. In the modernized module they are plain locals: the gfortran build is rescued only by -finit-local-zero in src/makefile line 85, while the ifort/ifx flag sets (makefile 162-171) have no such flag, so the reads are stack garbage; any NaN/Inf pattern makes fac*0 = NaN and poisons the diffusion coefficients. This path fires on every timestep for lithium/beryllium diffusion (diffuse_lithium), where the light element is exactly zero throughout the burned interior, so the exposure is large on Intel builds. Fix: assign all five arrays (zero) in the cycle branch.

## rotation/microdiff/gravitational_settling.f90:222-238 with implicit_diffusion_coeffs.f90:31-34 -- midpoint index misalignment and non-decaying correction at midpoint 1 in the implicit Newton update
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- grsett.f lines 159-161 (DO 30 I=2,NPT: EX_H(I)=0.5*(EX(I)+EX(I-1)-...)) and get_imp_diffco.f lines 12-14 (DO 10 I=1,NPT-1: ECOD2_H(I)=ECOD2_H(I)+EQCOD2X_H(I)*EX_H(I))
- detail: The midpoint arrays (equal_mass_mid, equal_diffusion_coeff2_mid, equal_diffusion_coeff2_dx_mid) are defined by model_to_equal with element i = midpoint between centre points i and i+1 (i=1..n-1). The iteration loop, however, stores the change at the midpoint between i-1 and i into element i (i=2..n). implicit_diffusion_coeffs then corrects D2_mid(i) with delta_abundance_mid(i) for i=1..n-1, i.e. with the abundance change one cell below the coefficient's own midpoint. Element 1 is never rewritten inside the loop, so it retains lax_wendroff_step1's provisional half-step change and, because the correction is cumulative (intent(inout) D2_mid += dD/dX * delta each iteration), that fixed increment is re-added on every iteration instead of decaying. Impact: the converged X profile of the old grsett path is slightly inconsistent at the base of the diffusing region and the coefficients are evaluated at a shifted composition; second order in dX per step, so small. microdiff_run has the same indexing (lines 172-176) but zeroes the derivative array so it is inert there.

## rotation/microdiff/microdiff_run.f90:204-213 and gravitational_settling.f90:250-258 -- convergence test ignores negative corrections
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited -- microdiff_run.f line 141 (IF(DX.GT.DXMAX)) and grsett.f lines 182-183 (DX = EX(I)-EX_P(I); IF(DX.GT.DXMAX))
- detail: max_abundance_change starts at 0 and is updated only when dx = X_new - X_prev is positive, so the loop declares convergence as soon as the largest POSITIVE correction is below settling_tolerance, regardless of how large the negative corrections are. In every settling step roughly half the grid (the region losing the species, e.g. the base of the diffusing region for hydrogen, the surface for metals) has negative corrections. The implicit solve is a linear tridiagonal system per iteration, so in practice the iterate changes are of similar magnitude with both signs and the test usually still works; but it can and will exit after one iteration in cases where the positive corrections happen to be small (e.g. a nearly uniform loss profile), leaving the coefficient correction unconverged. The printed DXMAX diagnostic is also misleading. Should be abs(dx).

## rotation/microdiff/gravitational_settling.f90 (whole implicit block, 222-262) and gravitational_settling_setup.f90:362-373 -- metal self-diffusion (second-derivative) term is computed, interpolated and never applied in the grsett path
- class: physical
- severity: low
- confidence: high
- provenance: inherited -- grsett.f lines 141-199 solve the implicit system for EX only; ECOD2Z(I) = ABS(FAC*AH) is set in setup_grsett.f (line ~262) and interpolated in model_to_equal.f but never used
- detail: For hydrogen the algorithm is Lax-Wendroff (pressure/temperature "settling" term) plus an implicit tridiagonal solve for the concentration-gradient term. For metals (use_diffusion_z in the grsett path) only the two Lax-Wendroff steps are executed on star%metal_abundance_change; rot_scr%src_grid_metal_diffusion_coeff2 / eq_metal_diffusion_coeff2_mid are filled and carried around but nothing reads metal_diffusion_coeff2_mid. The metal transport therefore has no gradient-limiting (diffusive) term and no hydrogen-gradient (A_H) term at all, so Z profiles produced by this path can steepen without bound at the base of the convection zone and at zone_begin; the newer microdiff path (species_col=3) does include both terms. Observable: differences in surface Z depletion between the two settling drivers and Z-profile spikes at zone boundaries in long grsett runs. (Also the name/comment "COD2Z = ABS(FAC*AH)" uses the hydrogen-gradient coefficient A_H^Fe as if it were the self-diffusion coefficient A_Fe,Fe; it is unused, but if someone wires it up it is the wrong coefficient.)

## eos/yale/fully_ionized_eos.f90:244-255,290-302 -- Debye-Huckel contributions to the second derivatives double-count the DH cross term
- class: physical
- severity: low
- confidence: medium
- provenance: inherited -- eqrelv.f lines 206-217 (QQPTT/QQPDT/QQPDD DH terms), 252-263 (QUDT/QQUTT/QQUDT DH terms), identical expressions
- detail: With P = N + P_dh, P_dh ~ rho^1.5 T^-0.5 and pdhp = P_dh/P, chi_T = N_T/P - 0.5 pdhp. Differentiating that exactly gives d chi_T/dlnT = (N_TT)/P - chi_T^2 + 0.25 pdhp (the -0.5 chi_T pdhp from -(N_T/P) chi_T cancels the +0.5 chi_T pdhp from d(pdhp)/dlnT). The code instead adds (0.5*qpt + 0.25)*pdhp on top of -qpt**2 where qpt already contains the DH part, leaving a spurious +0.5 chi_T pdhp; the same pattern gives spurious +0.5 chi_T pdhp in qqpdt and -1.5 chi_rho pdhp in qqpdd, and +0.5 qutd udhu in qqutt. Likewise qudt = pdu*(1-qpt) already contains the full DH part of (dlnu/dlnrho)_T through chi_T (the identity (du/dlnrho)_T = (P/rho)(1-chi_T) holds for the DH free energy, u_dh = 3 P_dh/rho), so the extra "+0.5*udhu" doubles it. First-order quantities (chi_T, chi_rho, c_P, grad_ad, dlnrho/dlnT, dlnrho/dlnP) are correct; only the derivative outputs (dlnrho_dlnt_dt, dlnrho_dlnp_dt, specific_heat_cp_dt/dp, adiabatic_gradient_dt/dp) are off by O(pdhp) ~ 1e-3 relative, i.e. Henyey Jacobian entries only, when use_debye_huckel_correction is on. Effect is on convergence rate, not the converged model.

## io/equal_to_model.f90:46-49 vs 96-106 vs 129-133 -- asymmetric clamping of the diffused hydrogen: no lower bound in the interior, no upper bound in the centre block
- class: logical
- severity: low
- confidence: medium
- provenance: inherited -- equal_to_model.f (same MAX in the centre loop, MIN(...,XMAX) in the interior and surface loops)
- detail: Zones 1..zone_begin get X = max(X + dX(1), 0) (floor, no cap); zones zone_begin+1..zone_end-1 and zone_end..num_zones get X = min(X + dX, 1 - Z - He3) (cap, no floor). Because hydrogen rises, the largest negative dX occurs just above zone_begin, in the interior loop, which has no floor, while the centre block (where X is by construction below hydrogen_diffusion_floor) gets the same dX(1) applied uniformly to every zone below zone_begin, including zones that are hydrogen-free. A negative X then feeds Y = 1 - X - Z - He3 > 1 and log(X) in the next microdiff_coefficients call (X is floored at 1e-24 there, but the model composition itself is left negative). This is the kind of composition corruption that can precede the "UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL 1" abort in the burning solver (the abort names shell 1, which is exactly the block that receives the uniform dX(1)); I could not prove that linkage from the code alone, so confidence is medium.

## io/model_to_equal.f90:125-127,233-236,255,316-319,329 -- "metal abundance" on the equal grid is taken from composition(8,:) (N15), not composition(3,:) (Z)
- class: logical
- severity: low
- confidence: high
- provenance: inherited -- model_to_equal.f lines 72, 153, 173, 224, 235 all use HCOMP(8,...)
- detail: star%metal_abundance_change (EZ) and rot_scr%metal_abundance_change_mid (EZ_H) are initialised from species column 8, which is i_n15 in the current composition layout, while the diffusion coefficients use composition(3,:) = Z. Today this is harmless because lax_wendroff_step1 overwrites EZ_H(1..n-1) with the provisional change and gravitational_settling only ever uses EZ - EZ_orig (the base value cancels); so it is a latent, not active, error. It is reported because it is one refactor away from mattering (any use of EZ as an absolute abundance, e.g. in a floor test) and because the code comment says "metal abundance".

## rotation/microdiff/microdiff_mte.f90:97-102 (and the analogous first-midpoint search in model_to_equal.f90) -- first-midpoint search bounded by the equal-grid count rather than the model grid
- class: logical
- severity: low
- confidence: medium
- provenance: inherited -- microdiff_mte.f lines 59-62 (DO IU=2,NPT over HRU with NPT = number of equal points)
- detail: `do iu=2,num_eq_points; if(radius_bl(iu).ge.eq_mid%radius(1)) exit` scans the MODEL radius array with the EQUAL-grid count as upper bound and starts at 2 rather than zone_begin+1. The target eq_mid%radius(1) = radius_bl(zone_begin) + dr/2 lies between model points zone_begin and zone_begin+1 (or a few above), so the scan needs iu to reach zone_begin+1. If num_eq_points < zone_begin+1 (a deep hydrogen-depleted core with many shells plus a coarse equal grid, num_eq_points is driven only by the 20 shells under the convection zone) the loop falls through, iu is clamped to num_eq_points and the linear interpolation is done between two unrelated deep-interior shells, silently corrupting mass/density/T/dlnP/dr/gradT/composition at the first midpoint. No abort, no message.

---

## Weak/uncertain observations

- microdiff_coefficients.f90:87 and diffusion_term(i)=0.0 / dlncdr=0.5*(...) literals are single precision; header says deliberate, not reported.
- gravitational_settling_setup.f90:335-343: diffusion_coeff1/coeff2 carry constant_mixing_coeff but diffusion_coeff1_dx/coeff2_dx do not (inherited); since the implicit update is cumulative in dD/dX this changes the converged answer, not just the rate, whenever constant_mixing_coeff /= 1.
- gravitational_settling_setup.f90:331-336: with the full Thoul solver only A_p and A_T come from the solver; A_c (self) and dA/dX still come from the Thoul fits (inherited partial implementation).
- microdiff_etm.f90 (zz2 = zz/composition(3,i)) and equal_to_model.f90:57,115,138 divide by Z; NaN for Z = 0 models with use_diffusion_z (inherited).
- microdiff_etm.f90: light elements in the envelope block receive no floor (can go negative), unlike X (inherited).
- microdiff_mte.f90:76 `mod(drtot,drmin).ne.0.0d0` float test -- effectively always true, only changes the count by one (inherited).
- microdiff_mte.f90:82 half_json hard-coded 5000 instead of json (fine while json = 5000).
- microdiff_mte.f90 and model_to_equal.f90 clamp k0 with num_zones-3 rather than zone_end-3, so the 4-point stencil may reach into the convective envelope (values exist, harmless).
- microdiff.f90:223 `species_fraction(3,i)` relies on the loop index value after the loop (i = num_eq_points); correct but fragile (same in microdiff.f line 180).
- model_to_equal.f90: dr1/dr2/fac1/fac2/delr computed and never used (inherited dead code).
- fully_ionized_eos.f90:204-209: on non-convergence (20 iterations) the routine writes a log line and returns with density and all derivative intent(out) arguments unassigned and log10_density unchanged; eqstat then continues with stale values, no ierr. Inherited silent failure.
- fully_ionized_eos.f90:90-91: jt1 clamped to 1 means log T < 5.2 is a 3-point Newton extrapolation of the Fermi table; only reachable if saha_log10t_cutoff is set below ~5.2.
- saha_eos.f90:311-315 and 387-391: for in_atmosphere or .not.want_derivatives the derivative intent(out) arguments are left unassigned; eqstat.f90 (Sweep J's file) subtracts saha_specific_heat_cp_dt etc. unconditionally at its interpolation block (~line 568), reading those undefined values (zeroed by -finit-local-zero on gfortran; results are unused when want_derivatives is false, so benign). Noted for Sweep J.
- thoul_diffusion.f90:123 writes mass_fraction(num_species) (electron mass fraction) into the caller's array -- both callers pass 4-element locals, so harmless, but intent(inout) on what is conceptually an input is a trap.
- lax_wendroff_step1 assumes the flux at the first equal-grid point (r = r(zone_begin), not r = 0 when the core is hydrogen-depleted) is zero; consistent no-flux BC but not documented.

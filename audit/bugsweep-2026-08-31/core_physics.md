# Bug sweep — core physics files (burn_lib, neutrino_flux_table, read_starting_model, turnover_timescale, observables_lib, stitched_model, monte_carlo, check_star_calibration, check_solar_calibration)

All files read in full; provenance checked against `git show 6cd5673:src/<name>.f`.

## core/burn_lib.f90:1196 (and :645) -- deuterium d(p,g)He3 rate uses T9^(+2/3) where CF88 has T9^(-2/3)
- class: physical
- severity: high
- confidence: high
- provenance: inherited (engeb.f:520 `RDEUT = RHO*2.240D3*T9P23*...`; deutrate.f:37 identical)
- detail: The nonresonant CF88 fit for H2(p,g)He3 is NA<sigma*v> = 2.24e3 * T9^(-2/3) * exp(-3.720/T9^(1/3)) * (1 + 0.112*T9^(1/3) + 3.38*T9^(2/3) + 2.65*T9) -- every other coefficient in the code (2.240d3, 3.72, 0.112, 3.38, 2.65) matches CF88 exactly, but the prefactor is coded as `t9_p23` = T9^(+2/3) (t9_p13 = pow(t9,cc13) is explicitly the PLUS power) in both engeb (line 1196) and deutrate (line 645). The routine's own temperature derivative proves the intent: `qrtdeut = cc13*((tfacdeut2/tfacdeut) - 2.0d0 - zz)` (line 1199) contains the -2/3 power-law term, i.e. the derivative is for a T9^(-2/3) rate while the rate itself uses T9^(+2/3). Net error factor T9^(4/3): at 1e6 K the D-burning rate is ~1e4 too SMALL, at solar center ~260x too small; pre-MS deuterium burning (dburn/dburnm consume deutrate's DRATE, engeb's dgdeut energy term uses the same rdeut) is drastically underestimated, and the internally inconsistent qrtdeut corrupts d(eps)/dlnT wherever D contributes.

## core/burn_lib.f90:1961 -- Itoh-branch dlnepsilon_dlnt built from dlnepsilon_dlnrho (copy-paste)
- class: logical
- severity: medium
- confidence: high
- provenance: inherited (engeb.f:1168-1169 `PEP = PEP + DSNUDD` / `PET = PEP + DSNUDT`)
- detail: In compute_neutrino_emission, when use_itoh_neutrino_loss is on and logT > nuclear_logT_cutoffs(5): `dlnepsilon_dlnrho = dlnepsilon_dlnrho + neutrino_dlnq_dlnd` then `dlnepsilon_dlnt = dlnepsilon_dlnrho + neutrino_dlnq_dlnt`. The second line reads the just-updated RHO derivative instead of dlnepsilon_dlnt, so the T-derivative of the energy generation returned to the Henyey solver is (rho-derivative + neutrino T-term) -- the nuclear sum3 T-derivative is discarded. Fires for every zone of an evolved star in the neutrino-cooled-core regime; degrades the luminosity-equation Jacobian (convergence, not converged-model values). Verbatim in the F77 original, so inherited.

## core/burn_lib.f90:1953-1961 -- Itoh branch adds bare log-derivatives to engeb's absolute-derivative accumulators
- class: physical
- severity: medium
- confidence: medium
- provenance: inherited (engeb.f:1160-1169)
- detail: engeb's returned "dlnepsilon" derivatives are actually ABSOLUTE derivatives d(eps)/dln(T,rho): sum3 = sum_i DG(i)*dlnrate_i = eps * dln eps/dlnT (erg/g/s), and both downstream Henyey use (henyey_coefficients.f90:291-294 scales ql and ql_dt by the same cccql factor; henyey_eliminate uses them with a ln10 factor as absolute d(eps)/dlogT) and the old (non-Itoh) neutrino branch (qetnx accumulates dln-terms multiplied by the absolute rates ex1..ex3 before subtracting, lines 2006-2021) are consistent with that. The Itoh branch instead converts the NEUTRINO derivative to a dimensionless log-derivative (`neutrino_dlnq_dlnt = -dsnudt*T/enu`, line 1953) and adds it to the absolute accumulator -- a units mismatch that mis-weights the neutrino contribution by a factor |eps_nu| relative to the nuclear terms. Same in the original.

## core/stitched_model.f90:529-530 -- pulse builder swaps grad/grad_ad for the envelope region
- class: logical (copy-paste / parallel-branch asymmetry)
- severity: high (for GYRE/FGONG output)
- confidence: high
- provenance: new-code (stitched restructure 065c064; the matching PROFILE-column swap was later fixed -- see the comment at lines 333-339 -- but the pulse builder site was missed)
- detail: env_gradients ordering is (1) radiative, (2) ADIABATIC, (3) ACTUAL (envelope_derivs.f90:105-107 fills current_gradients(2)=grada, (3)=actual_gradient). ext_profile_value's envelope case correctly maps col 13 (gradT) -> (3) and col 14 (grada) -> (2). But build_pulse_points region 2 sets `nab = env_gradients(2,i)` (adiabatic!) and `nab_ad = env_gradients(3,i)` (actual!). Consequences in every envelope-region pulse point: GYRE columns nabla and nabla_ad are swapped, and the first-pass thermal N^2 = g^2 (rho/P) delta (nab_ad - nab) has its SIGN flipped -- so convective envelope points get N^2 > 0 and are then overwritten by the centered-difference form (exactly the cancellation noise the hybrid scheme was designed to avoid there), while radiative envelope points keep a wrong-sign negative thermal N^2 (GYRE sees convection where there is none). The atmosphere region is unaffected because atmo_gradients uses a different ordering ((2)=actual, (3)=adiabatic, envint_kernel.f90:495-497) that happens to match the code. Note the two structs' differing row orders are the trap that caused this.

## core/stitched_model.f90:508-509 -- eps_eps_T / eps_eps_rho double-multiplied by eps (over-correction of a known finding)
- class: physical
- severity: medium (nonadiabatic GYRE runs only)
- confidence: medium-high
- provenance: new-code (commit 2ed398c "Fix pulse-output physics: ... absolute kap/eps derivatives"); sharpens the known "kap/eps derivative factor convention in GYRE columns" finding
- detail: The fix writes `pts(ipul_eps_eps_T) = eps_total * pulse_dlneps_dlnt` on the premise that pulse_dlneps_dlnt is a bare log-derivative. For KAP that premise holds (kap_res(i_dlnkap_dlnt) is dlnkap/dlnT), but for EPS it does not: pulse_dlneps_dlnt = engeb's sum3 = sum_i DG(i)*dlnrate_i = d(eps)/dlnT = eps*dlneps/dlnT ALREADY (see previous finding's evidence from henyey_coefficients/henyey_eliminate scaling). So the written GYRE column is eps^2*dlneps/dlnT -- too large by a factor eps. The pre-fix code, which wrote pulse_dlneps_dlnt directly, was already correct for eps (while genuinely wrong for kap). Only eps_eps_T/eps_eps_rho need reverting to the bare stored value; the kap half of the fix should stay.

## core/burn_lib.f90:328-332,393-403 -- dburn reads total_shell_mass uninitialized for a single-zone CZ with accretion
- class: numerical (uninitialized read)
- severity: medium
- confidence: high (that the read is undefined; low probability of the trigger)
- provenance: inherited defect, modernization-aggravated (dburn.f: SUMM only set in the IBEGIN.NE.IEND branch, but the blanket `SAVE` made the read stale-yet-finite; the modernized routine has no save, so the value is genuinely undefined)
- detail: When zone_begin == zone_end, the multi-zone branch that sums total_shell_mass is skipped, yet `deuterium_fraction_test = (D*total_shell_mass + Dacc*amf)/(total_shell_mass+amf)` (line 330) and the accretion-weighted update (lines 393-403) still read total_shell_mass whenever use_mass_accretion .and. zone_end == num_zones (a one-zone surface convection zone). In the F77 original the SAVE meant it silently reused the previous call's CZ mass (wrong but bounded); after modernization it is an uninitialized automatic variable.

## core/burn_lib.f90:393-403 vs 565-575 -- dburn/dburnm accretion weighting asymmetry (units mixing)
- class: physical
- severity: medium (accretion runs only)
- confidence: medium
- provenance: inherited (dburn.f:101 `XH2N=(XH2NO*SUMM+XH2NA*FMASSACC)/(SUMM+FMASSACC)` vs dburnm.f:101 `XH2N=(XH2NO+XH2NA*FMASSA)/(1.0D0+FMASSA)`)
- detail: star%accreted_mass_fraction is documented (in both routines) as DMDT*DT / ORIGINAL CZ MASS -- a dimensionless fraction. dburnm combines original and accreted deuterium as (D + Dacc*f)/(1+f), which is dimensionally consistent with that definition. dburn instead combines (D*M_cz + Dacc*f)/(M_cz + f), adding a dimensionless fraction to the absolute CZ mass (sum of shell_mass, in the model's mass units). The two near-duplicate routines cannot both be right; unless shell_mass happens to sum to 1, dburn essentially ignores the accreted deuterium (f << M_cz in cgs-like units). Same asymmetry in both F77 originals.

## core/burn_lib.f90:2553-2565 (liburn) and 2922-2934 (liburn2) -- radiative_frac is (intended - 1): weights fall in [-1,0] instead of [0,1]
- class: physical
- severity: medium
- confidence: high (mismatch with its own comment); medium-high (on the intended form)
- provenance: inherited (liburn.f:448-458: DMASS = HSBEG - HSEND; FRAD = (HS1(I)-HSBEG)/DMASS; identical in liburn2.f)
- detail: For the zones that "began convective and ended radiative" (CZ base retreated outward, so mass_coord_end > mass_coord_beg), the comment defines the recipe: a point x of the way (in mass) from the old to the new base spends x of the timestep in the CZ, 1-x radiative, so FRAD (the radiative-rate weight) should be 1 - x with x = (m - beg)/(end - beg). The code computes radiative_frac = (m - beg)/(beg - end) = -x: 0 at the old base (should be 1) and -1 at the new base (should be 0). The blend exp(frad*ln r_rad + (1-frad)*ln r_cz) is therefore evaluated with a negative weight on the radiative rate and a weight up to 2 on the CZ rate -- the per-zone depletion in the transition region is scaled by r_cz/r_rad relative to the intended value. Affects Li6/Li7/Be9 depletion whenever the surface CZ retreats (pre-MS spin-down/Li-depletion science). Verbatim in both originals.

## core/read_starting_model.f90:604-613 -- core-extension JSON overflow prints "RUN TERMINATED" but does not terminate
- class: logical
- severity: high (when it fires: out-of-bounds writes across every model array)
- confidence: high
- provenance: inherited (starin.f:378-387: same IF/WRITE block with no STOP, falls through into the shift loop)
- detail: extend_core_toward_center checks `num_shells_extended .gt. json`, writes "STARIN: ***** RUN TERMINATED *****" -- and then continues: no `ierr = 1`, no `return`. The very next loop shifts every model array to index i + num_core_shells_added up to num_shells_extended > json, writing past the declared bounds of star%log_mass/logR/xa/etc. The original had the identical dead warning (no STOP). The modernized routine has the ierr channel one line away, so this is also the cheapest possible fix.

## core/read_starting_model.f90:956-961 -- inverted interpolation sign for X and Z at the new envelope fitting point
- class: physical (sign error)
- severity: low-medium
- confidence: high
- provenance: inherited (starin.f:668-669 `HCOMP(1,J) = HCOMP(1,M)+FX*(HCOMP(1,M)-ENVX(JM))`, while the neighboring branch at 691-692 interpolates correctly)
- detail: When the requested envelope mass falls between the last interior point and the first envelope point, every structural quantity is interpolated as star + frac*(env - star) (lines 943-955), but X and Z are computed as star + frac*(star - env) -- extrapolated AWAY from the envelope value by exactly the amount they should move toward it. He4 then absorbs the error via the 1-X-Z-He3 closure. Harmless when the envelope composition is uniform (the usual case), wrong whenever a composition gradient reaches the fitting region (diffusion/accretion runs). The parallel two-envelope-point branch (lines 1005-1012) is correct, which is what marks this one as the copy-paste error.

## core/burn_lib.f90:2194 (liburn) / :2718 (liburn2) -- overshoot CZ-base search mixes current and previous-model radii
- class: logical (asymmetry between mirrored blocks)
- severity: low-medium
- confidence: medium-high
- provenance: inherited (liburn.f: ending-depth block uses `DR = EXP(CLN*HRO(J+1))-R` with R from HR(J))
- detail: In the ENDING CZ depth search (envelope overshoot active), shell_radius and search_radius are built from the current model's radius array, but the zone width used to place cz_base_frac is `delta_radius = exp(ln10*star%logR_start(zone_idx+1)) - shell_radius` -- the START-of-step radius of the zone above minus the current radius of the zone below. The matching STARTING-depth block a few lines up consistently uses logR_start throughout, so the ending block's mixture looks like a copy-paste from it. Between models the two grids differ by the structural change over the step, so cz_base_frac (and hence the burn-rate/mass adjustment of the CZ base zone) is systematically perturbed. Inherited from the original.

## core/check_solar_calibration.f90:72 -- dlum_dalpha value contradicts its own empirical comment by a factor of 10
- class: physical (comment-vs-code mismatch)
- severity: low
- confidence: medium (that one of the two numbers is wrong; unclear which)
- provenance: inherited (chkcal.f:43 `DLDA = 0.0139D0  ! empirical result: +0.139  RMS error .0022`)
- detail: The four fixed Newton partials each carry a comment with the measured value; three agree with the code (-3.78/-3.783, -0.89/-0.890, -0.050/-0.0504) but the fourth codes 0.0139 against a stated "+0.139". The 2x2 solve (lines 79-81) is algebraically correct, and because dlogL/dalpha is the smallest coupling the calibration still converges with either value (an approximate-Jacobian iteration), just at a different rate -- but one of the two numbers is a typo, and if it is the code's the alpha corrections are damped ~10x through the L-coupling term.

## core/burn_lib.f90:649/664 vs 1196/1203 -- deutrate and engeb disagree on the NA/m_d factor (and share none of it with rdeutmax)
- class: physical (constant inconsistency)
- severity: low
- confidence: high
- provenance: inherited (deutrate.f uses 3.0115D23 = NA/2.000 with a comment claiming "Avogadro / mass of deuteron in amu"; engeb.f uses 6.023D23/ANUC(3) = NA/2.013553 = 2.9912e23)
- detail: The same physical factor (Avogadro's number over the deuteron mass in amu) is 3.0115e23 in deutrate (divides by exactly 2, contradicting its own comment) and 2.9912e23 in engeb -- a 0.7% inconsistency between the abundance-update rate (DRATE, used by dburn/dburnm) and the energy-generation rate for the identical reaction. deutrate's overturn cap (rdeutmax) uses the same NA/2 value. Trivially fixable; noted mostly because the two sites are supposed to be the same rate.

## Known-item sharpening: core/read_starting_model.f90:616-673 -- extend_core never shifts logRho, and the ideal-gas offset mixes shells
- class: logical/physical
- severity: low-medium
- confidence: high
- provenance: inherited (starin.f: the shift loop moves HS/HR/HL/HP/HT/LC/HCOMP/OMEGA but not HD; RHOC = HD(1); FACD = HP(MP1)-HD(MP1)-HT(MP1)). LABELED KNOWN per instructions ("read_starting_model's ideal-gas core extension").
- detail: Sharpening the known audit item with the specific defects: (a) the array-shift loop omits star%logRho, so after shifting, all interior densities sit mcore indices too low until shell_physics recomputes them -- but central_log_density (= unshifted logRho(1), coincidentally the correct old-center value) and especially density_estimate_offset = logP(mp1) - logRho(mp1) - logT(mp1) are formed from MIXED shells: P and T are the shifted old shell 1, rho is the unshifted old shell mcore+1. The new core points' densities inherit the log(rho_1/rho_{mcore+1}) offset. (b) The new core radii use the same slightly-off central density. Self-limiting because shell_physics re-derives density from (P,T), but the geometric quantities (logR of the new points) keep the error.

## Weak/uncertain observations
- burn_lib.f90:1364-1372 (engeb): the rate<1e-30 branch zeroes dlnrate_dlnt but not dlnrate_dlnrho; the F77 SAVE made the subsequent 0*stale product harmless, the modernized automatic local is read undefined on the first call (0*undef is only safe if the garbage isn't NaN/Inf). Inherited pattern, modernization-degraded.
- burn_lib.f90:346-351 (dburn) / 513-518 (dburnm): comment promises max 0.5 per substep in ln, code caps 2*rate*dt at 0.2 (divisor 0.1) -- conservative, comment stale (inherited).
- burn_lib.f90:1328 (engeb): `hydrogen_fraction.eq.0.0` float-equality gate -- KNOWN audit item (burn_lib float-equality guards), not re-argued here.
- burn_lib.f90:1953-1954: Itoh branch divides by star%neutrino_loss_rate with no zero guard (NaN if snu underflows to 0 just above the T cutoff); inherited.
- burn_lib.f90:2212 vs 2736: liburn `exit`s the radiative loop on a below-threshold rate while liburn2 `cycle`s -- both faithful to the originals (GOTO 60 vs GOTO 50); the asymmetry is inherited and matters only if rates are non-monotonic across zones.
- burn_lib.f90:2528 (liburn/liburn2): mixed CZ abundances are written from max_zone (=old base when the CZ deepened) outward, leaving the newly engulfed zones [cz_base_zone, max_zone-1] with their radiative-loop abundances even though they were averaged into the CZ sums -- inherited (liburn.f DO 130 I=MAXJ,M).
- check_star_calibration.f90:63-64: `log_r_rsun_current` (and the ALR comment block) actually holds LINEAR R/Rsun -- faithful to chkscal.f (ALR=SQRT(...)/(TEFF^2*CRSUN)); the math is self-consistent in linear R, only names/comments lie.
- check_solar_calibration.f90:116: log_zx_mismatch_prev is assigned from log_zx_mismatch, which is unset on any call that misses the L/R tolerance gate; the original's blanket SAVE is gone, so this is an undefined read -- but the value is never consumed (dead store).
- observables_lib.f90:158-159: core_boundary_fx2 = (grada(i+1)-gradr(i))/(grada(i+1)-gradr(i)) is identically 1.0 AND is then not even used (the stale envelope_boundary_fx is) -- both quirks documented in-file as the preserved FX/FX2 bug; consumer is dead code.
- neutrino_flux_table.f90:131: "shell mass fraction" divides by the SOLAR mass, not the star's mass -- faithful to main.f's LNUTAB block (FM = HS2/1.9891D33); only mislabels the diagnostic for non-solar stars.
- monte_carlo.f90:202-203: rescale_params(:,run_index-2) underflows the array for run_index < 3 (inherited RESCAL(2,NK-2) in wrtmonte.f); MC protocol likely guarantees nk>=3 at this call site.
- turnover_timescale.f90:202-207: convective_flag(cz_base_index+1/+2) can index past num_points near the surface (in-bounds of json, but reads stale data from a previous, longer walk); inherited shape from tauintnew.f.
- burn_lib.f90:84 etc.: engeb's early-return leaves total_energy_gen_rate (intent(inout)) untouched -- documented in-file as a preserved property of the original; not re-reported.

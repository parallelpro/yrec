# Sweep B (pass 2): nuclear energy generation, rates, neutrino losses, burning

Date: 2026-09-01. Read-only review; no files edited, nothing built or run.

## Files read in full

- /Applications/YREC/src/core/burn_lib.f90 (3074 lines): eqburn, dburn, dburnm,
  deutrate, engeb (+ internal setup_abundances_and_composition,
  compute_energy_generation, compute_neutrino_emission), liburn, liburn2,
  lirate88.
- /Applications/YREC/src/net/net_lib.f90 (2416 lines): neutrino, nulosses,
  azbar, sneut, rates, safedivexp, ifermi12, zfermim12.

F77 ancestors consulted at rev 6cd5673 (git show only): engeb.f, deutrate.f,
sneut.f, dburn.f, dburnm.f, liburn.f, liburn2.f, lirate88.f, eqburn.f,
rates.f, main.f, getw.f.

Call sites checked with `grep -rn "call <name>" /Applications/YREC/src`
(no head): engeb (henyey_coefficients.f90:247, neutrino_flux_table.f90:89,
timestep_limit_heburn.f90:81); rates (mix.f90:182,
timestep_limit_hburn.f90:99, net/test/test_net.f90:150); eqburn
(mix.f90:299,331, timestep_limit_hburn.f90:105); dburn (mix.f90:525,538);
dburnm (mid_timestep_model.f90:401,415); deutrate (mix.f90:192,
test_net.f90:179); lirate88 (evolve_step.f90:336 mode 1, evolve_step.f90:541
mode 2, evolve_angular_momentum.f90:383 mode 2); liburn (evolve_step.f90:543,
evolve_angular_momentum.f90:391); liburn2 (evolve_angular_momentum.f90:404,
burn_settle_mix.f90:136); neutrino (burn_lib.f90:1947); safedivexp
(burn_lib.f90:2569,2572). All are module procedures (burn_lib / net_lib), so
argument count/type/order are compiler-checked; every actual list was compared
with the dummy list and matches.

## Verified clean (re-derived, found correct)

- Bahcall-form rates 1-7 (pp, He3+He3, He3+He4, C12+p, C13+p, N14+p, O16+p):
  r1 polynomial, exp(Q6 T9^-1/3 + Q7 + (Q8 T9)^2 + screening), and
  dlnrate/dlnT = dscreen - [Q6 T9^-1/3 + (2 T9^-2/3 + Q1 T9^-1/3 - Q3 T9^1/3
  - 2 Q4 T9^2/3 - 3 Q5 T9)/r1]/3 + 2 (Q8 T9)^2 -- all terms verified
  (burn_lib.f90:1355-1372), and the identical forms in net_lib rates
  (1920-2248).
- Alpha-capture RATE values (not derivatives) C13(a,n)O16, C12(a,g)O16,
  N14(a,g)F18 against CF88 tables; conversion factors N_A/(4*A) =
  1.157126e22 (A=13), 1.25388e22 (A=12), 1.07452e22 (A=14).
- Triple-alpha rate (2.79e-8 T9^-3 exp(-4.4027/T9) + 1.35e-8 T9^-3/2
  exp(-24.811/T9)) and both derivatives; 1.565315e21 = N_A^2/(3! * 4^3 ...)
  reproduces Kippenhahn's 5.09e11 rho^2 Y^3 T8^-3 form.
- Screening (Graboske/DeWitt): weak u = zeta z1z2 (rho/T^3)^1/2 with
  0.5u / -1.5u derivatives; intermediate 0.86 exponent with xtr = 1.58 =
  3*0.86-1, derivatives 0.43u / -1.29u; strong Lambda0^(2/3) form with
  derivatives; Lambda0 constant 5.9426e-6 = 1.88e8/10^13.5.
- Energy releases per reaction (6.664, 12.860, 1.586 + f1 17.394 + f2 11.499,
  3.457, 7.551, 9.054 + f3 4.966 + f4 12.128, 3.553, 2.216, 7.162, 5.815,
  7.275 MeV) and the Be7 e-/p-capture and N15 branching used to build them.
- Fermi-momentum constant 1.017677e-4 and 5.92986e9 K in the degenerate
  screening branch; the old (Beaudet-type) pair/photo/plasma neutrino fits and
  their rho^-1 / rho^0 / rho^2 scalings in the derivative bookkeeping.
- eqburn: He3 equilibrium quadratic, CN closure, dlnrate handling; c21 =
  5.240358e-8.
- dburn/dburnm: the 0.5 factor in the D(p,g) rate is compensated by the
  exp(-2 rate dt) update; sub-stepping and accretion bookkeeping (apart from
  finding 8).
- liburn refine/exit ladder (GOTO elimination faithful); timestep is passed in
  seconds and rates are per second (the "IN YEARS" comment is stale, code is
  consistent).
- lirate88 Li6/Li7/Be9 CF88 rate forms; neutrino -> nulosses -> azbar -> sneut
  plumbing; ifermi12/zfermim12 Antia rational-fit structure.
- Interfaces: all argument lists at all call sites listed above.

---

## core/burn_lib.f90:1196 -- D(p,gamma)He3 rate uses T9^(+2/3) instead of T9^(-2/3)
- class: physical
- severity: high
- confidence: high
- provenance: inherited -- `git show 6cd5673:src/deutrate.f | grep -n RDEUT` -> line 36 `RDEUT = RHO*2.240D3*T9P23*EXP(Z)*TFACDEUT*3.0115D23`; `git show 6cd5673:src/engeb.f` line 520 `RDEUT = RHO*2.240D3*T9P23*EXP(ZZ)*TFACDEUT*6.023D23/ANUC(3)`.
- detail: CF88 D(p,g)He3 is N_A<sv> = 2.24e3 T9^-2/3 exp(-3.720 T9^-1/3)(1 + 0.112 T9^1/3 + 3.38 T9^2/3 + 2.65 T9). Both the engeb copy (1196) and the deutrate copy (645) multiply by `t9_p23` (T9^+2/3). The derivative on the very next lines, `qrtdeut = cc13*((tfacdeut2/tfacdeut) - 2.0d0 - zz)` (1199), is the derivative of the T9^-2/3 form (the -2 = 3 * (-2/3)), so the code's own d ln r/d ln T contradicts its rate: intent evidence. At T = 1e6 K the coded rate is low by T9^-4/3 ~ 1e4; the D-burning threshold is pushed from ~1e6 K to ~2e6 K. Observable: pre-MS deuterium burning (dburn/dburnm, eps_D in engeb) occurs later/hotter than it should, altering the pre-MS luminosity plateau and the D depletion age in every run with use_extended_composition; also the deuterium rate derivative fed to the Henyey Jacobian is inconsistent with the rate.

## core/burn_lib.f90:1499-1540 -- alpha-capture dlnrate/dlnT (reactions 8, 10, 11) misses the 1e22*exp(screening) factor: derivative is effectively zero
- class: numerical
- severity: medium
- confidence: high
- provenance: inherited -- `git show 6cd5673:src/engeb.f | grep -n "RHO/RATE"` -> 823 `DRATT(8) = DSCT(8)+RHO/RATE(8)*(DR1*A1 + R1*DA1 +`, 836 `DRATT(10) = DSCT(10)+RHO/RATE(10)*(...`, 852 `DRATT(11) = DSCT(11)+RHO/RATE(11)*(...`.
- detail: reaction_rate(8) = 1.157126e22 * rho * exp(scr) * S(T), but dlnrate_dlnt(8) = dscreen + rho/rate * dS/dlnT, i.e. the ratio is (1/1.157e22/exp(scr)) * dlnS/dlnT ~ 1e-22 instead of dlnS/dlnT (which is ~20 for C12(a,g) at T8 ~ 1.5). Same for (10) and (11). Additionally, for (8) the a2..a5 terms in the derivative omit the t9_m32 factor present in the rate, and for (10) the r1/r2 terms `r1*cc13*0.1956` and `cc13*1.0616*r2` omit the T9^-2/3 factor (d/dlnT of 1/(1+c T9^-2/3) = (2/3) c T9^-2/3 r^2). Net effect: the Henyey Jacobian term ql_dt (henyey_coefficients.f90:291-294) gets only the screening derivative for C13(a,n), C12(a,g) and N14(a,g); the temperature sensitivity of these reactions is invisible to the Newton solver. Fires whenever log T >= nuclear_logT_cutoffs(4) and He-burning eps is significant. This is the regime of the known Solar_m1p0 GN93 ZAHB->TAHB NaN near core He exhaustion, where C12(a,g) becomes the dominant energy source (Y -> 0, 3-alpha ~ Y^3 dies first): a Jacobian that says d eps/dT ~ 0 while the true value is ~20 eps is a plausible contributor to the corrector overshooting into an unphysical (T, rho) where the MLT cubic goes negative; not proven, but the Jacobian is provably wrong there.

## net/net_lib.f90:343 -- Itoh tfac2: `cap-cap` should be `cap*cap` (C_V'^2 - C_A'^2)
- class: physical
- severity: high (only with use_itoh_neutrino_loss = .true.)
- confidence: high
- provenance: inherited -- `git show 6cd5673:src/sneut.f | grep -n tfac2` -> 124-125 `tfac2 = cv*cv - ca*ca + (xnufam-1.0d0)*(cvp*cvp - cap-cap)`.
- detail: Itoh et al. (1996) define the second coupling combination as C_V^2 - C_A^2 + (n-1)(C_V'^2 - C_A'^2). With cv = 0.5+2 sin^2(theta_W) = 0.9618, ca = 0.5, cvp = 1-cv = 0.0382, cap = 0.5, the correct tfac2 = 0.675 + 2*(0.00146 - 0.25) = 0.1781, but the coded `(cvp*cvp - cap - cap)` = 0.00146 - 1.0 = -0.9985 gives tfac2 = 0.675 + 2*(-0.9985) = -1.3220. Hence tfac3 = tfac2/tfac1 = -0.786 (should be +0.106) and tfac5 = -0.661 (should be +0.089). tfac3 enters the pair term `a1 = tfac4*(1+tfac3*qpair)` (net_lib.f90 ~574), the photo term (~1052), and bremsstrahlung (~1233, ~1389); tfac5 enters brem. Pair and photo neutrino losses are therefore wrong in sign/magnitude in the "(1 + tfac3*q)" correction, by tens of percent at high T where qpair, qphoto -> O(1). Observable: with LNULOS1 on, neutrino cooling of degenerate cores / late He-burning cores is wrong; default is off, so standard runs are unaffected.

## core/burn_lib.f90:1961 -- Itoh branch: dlnepsilon_dlnt built from dlnepsilon_dlnrho (copy-paste)
- class: logical
- severity: high (only with use_itoh_neutrino_loss = .true.)
- confidence: high
- provenance: inherited -- `git show 6cd5673:src/engeb.f | grep -n "PET = "` -> 1169 `PET = PEP + DSNUDT` (line 1168 is `PEP = PEP + DSNUDD`).
- detail: the code does `dlnepsilon_dlnrho = dlnepsilon_dlnrho + neutrino_dlnq_dlnd` and then `dlnepsilon_dlnt = dlnepsilon_dlnrho + neutrino_dlnq_dlnt`, discarding the nuclear sum3 (d eps/d ln T) accumulated at 1900-1930 and replacing it with the density derivative plus the neutrino T-derivative. The old-routine branch (1975-1980) correctly does `dlnepsilon_dlnt = dlnepsilon_dlnt - qetnx`. Whenever the Itoh path is taken the Henyey Jacobian's d eps/d ln T is garbage (equals d eps/d ln rho + neutrino term, typically ~1 instead of ~5-40). Nuclear eps itself is unaffected, so structure converges slowly/erratically rather than to a wrong answer -- but combined with finding 5 the two neutrino terms also have the wrong units.

## core/burn_lib.f90:1953-1961 -- Itoh branch adds dimensionless d ln q/d ln T to dimensional d eps/d ln T; and 0/0 when sneut returns zero
- class: physical / numerical
- severity: medium
- confidence: high
- provenance: inherited -- `git show 6cd5673:src/engeb.f` 1161 `DSNUDT = -DSNUDT*T/ENU`, 1168-1169; old branch 1227-1228 `PET = PET - QETNX`, `PEP = PEP - QEDNX`.
- detail: sum2/sum3 in engeb are Sum_i eg(i)*dlnrate_i, i.e. dimensional d eps/d ln X in erg/g/s (henyey_coefficients.f90:291-294 consumes them that way). The Beaudet branch subtracts qetnx/qednx, which are likewise dimensional (q * dlnq). The Itoh branch instead converts the sneut derivatives into logarithmic derivatives (`-dsnudt*T/snu`) and adds those dimensionless numbers to the dimensional sums: the neutrino contribution to the Jacobian is off by a factor |snu|. Secondly, sneut returns snu = 0 (early `if (temp .lt. 1.0d7) return`, net_lib.f90:397) while engeb calls it for any zone with log T >= nuclear_logT_cutoffs(5); the examples set TCUT(5)=7.5 but the registry default is 0.0, so a user who enables Itoh losses without setting TCUT(5) >= 7 gets `neutrino_dlnq_dlnt = -0*T/0` = NaN in every zone below 1e7 K, propagated into the Jacobian.

## core/burn_lib.f90:2558 (and 2927 in liburn2) -- radiative_frac for zones that changed from convective to radiative is in [-1, 0], not [1, 0]
- class: physical
- severity: medium
- confidence: high
- provenance: inherited -- `git show 6cd5673:src/liburn.f | grep -n "DMASS\|FRAD"` -> 449 `DMASS = HSBEG - HSEND`, 454 `FRAD = (HS1(I)-HSBEG)/DMASS`; liburn2.f 292, 297 identical.
- detail: For zones between the old CZ base (cz_base_zone_old) and the new, shallower base (cz_base_zone), the comment says a point 1/3 of the way from the old to the new base "SPENT 1/3 OF THE TIME IN THE CZ AND 2/3 OUT", so radiative_frac should run from 1 (at the old base, zone left the CZ immediately) to 0 (at the new base). With delta_mass = mass_coord_beg - mass_coord_end < 0 (mass coordinate increases outward, beg < end) and numerator (m - m_beg) >= 0, the code gives radiative_frac in [-1, 0]: exactly 1 less than intended. The depletion exponent radiative_frac*log(rate_rad) + (1-radiative_frac)*log(rate_cz) therefore weights the CZ-base rate by 1..2 and the local radiative rate by 0..-1. Since the local radiative rate at the old base is smaller than the CZ-average start rate, the mis-weighting over-depletes Li/Be in the just-decoupled shells (by exp of the log-ratio, easily a factor of several per step). Observable: pre-MS Li7 profile just below the retreating envelope CZ, and thus the surface Li that is dredged later; also the Li6/Be9 profiles. Fires every step in which the envelope CZ retreats (all of the pre-MS).

## core/burn_lib.f90:2985-3060 + core/evolve_step.f90:541 -- "end of the time step" Li/Be rates are never computed; both rate arrays hold start-of-step values
- class: logical
- severity: medium
- confidence: high
- provenance: inherited -- `git show 6cd5673:src/lirate88.f | grep -n "J.EQ\|HDO\|RLI60"` -> 31 `IF(J.EQ.1)THEN`, 35 `RHOX = EXP(CLN*HDO(I))*HCOMPP(1,I)`, 100 `RLI60(I) = FLI6*RHOX`; `git show 6cd5673:src/main.f` 596-640: `HDO(I) = HD(I)` (611) then `CALL LIRATE88(HCOMP,HD,HT,M,1)` (634) before Henyey, and 843 `CALL LIRATE88(HCOMP,HD,HT,M,2)` after; getw.f 325 same mode-2 call.
- detail: lirate88 with use_current_model /= 1 ignores its passed composition/log_density/log_temperature arguments and evaluates star%logRho_start, logT_start, xa_start into star%rate_*_start. evolve_step.f90:541 calls it in mode 2 with star%xa/logRho/logT under the comment "FIND BURNING RATES AT THE END OF THE TIME STEP"; the mode-1 call that fills star%rate_* (the "end" slot used by liburn) is at evolve_step.f90:336, before the Henyey solve, from the same start-of-step model. Result: rate_li6 == rate_li6_start (etc.) at every liburn call, the log-linear rate interpolation in liburn (2232-2241, 2745) collapses to the start-of-step rate, and the end-of-step model's temperature never enters the light-element burn: an explicit-Euler-in-rate scheme mislabelled as trapezoidal. The same holds in the rotation path (evolve_angular_momentum.f90:383, then `rate_*_start = rate_*` copies identical values). Consequence: timestep-dependent under-depletion of Li when the CZ base is heating (pre-MS), the opposite when cooling. The comment/code mismatch is itself the finding; the fix is either a mode-1 call after the Henyey solve or having mode 2 use the passed arrays.

## core/burn_lib.f90:296-330 (dburn) and dburnm -- total_shell_mass read uninitialized in single-zone accretion path
- class: numerical
- severity: medium
- confidence: medium
- provenance: modernization -- `git show 6cd5673:src/dburn.f | grep -n "SAVE\|SUMM"` -> 16 `SAVE`, SUMM set only inside the multi-zone branch (24-42) but read at 45 `XH2TEST = (XH2*SUMM+COMPACC(12)*FMASSACC)/(SUMM+FMASSACC)` and 99-101; dburnm.f 17 `SAVE`, same structure. The F90 total_shell_mass is a plain local.
- detail: When zone_begin == zone_end the first branch is taken and total_shell_mass is never assigned; the accretion block (`use_mass_accretion .and. zone_end.eq.num_zones .and. mass_accretion_rate > 0`) then divides by (total_shell_mass + accreted_mass_fraction) and uses it in the composition update (lines ~400-410 as well). mix.f90:525 calls dburn one zone at a time (zone_begin = zone_end) for radiative zones, so this fires for an accreting model whose surface zone is radiative. In the F77 the SAVE made SUMM a stale value from a previous multi-zone call (also wrong, but finite); in the F90 it is whatever is on the stack, so the surface D abundance under accretion is garbage or NaN. Same defect in dburnm (mid_timestep_model.f90:401/415).

## core/burn_lib.f90:1364-1367 -- dlnrate_dlnrho(i) left unset when a Bahcall-form rate underflows
- class: numerical
- severity: low
- confidence: medium
- provenance: modernization -- `git show 6cd5673:src/engeb.f` 688-691: the RATE<1e-30 branch sets `RATE(I)=0.`, `DRATT(I)=0.` but not DRATRO(I); F77 arrays were static (SAVE'd by default in that compiler model), the F90 `dlnrate_dlnrho(13)` is an automatic local declared at 885-890.
- detail: For reactions with reaction_rate(i) < 1e-30 only dlnrate_dlnt is zeroed. eg(i) = 0 so the later sum2 accumulation `eg(i)*dlnrate_dlnrho(i)` multiplies 0 by an uninitialized double; that is 0 unless the stack garbage is NaN/Inf, in which case the whole d eps/d ln rho becomes NaN. Deterministic garbage on a given build makes this hard to see, but it is a genuine read-before-write that valgrind/-finit-real=snan would flag, and the fix is a one-line `dlnrate_dlnrho(i)=0.` in the same branch.

## core/burn_lib.f90:2194 (and 2718 in liburn2) -- delta_radius mixes previous-model logR(zone+1) with current-model radius(zone)
- class: logical
- severity: low
- confidence: medium
- provenance: inherited -- `git show 6cd5673:src/liburn.f | grep -n "DR = "` -> 115 `DR = EXP(CLN*HRO(J+1))-R`; liburn2.f 99 same.
- detail: The pressure-scale-height / overshoot geometry for the CZ base compares exp(logR_start(zone+1)) (start-of-step radius of the next shell) with shell_radius = exp(ln10*radius(zone)) from the current model, whereas the sibling expression two lines above (2185) uses the current model for both. After a rezone or a step with appreciable contraction the two radii belong to different grids/epochs and delta_radius can be of either sign; it feeds cz_base_frac (consumed at 2356-2363 / 2823-2829) which weights the CZ-base rate. Small effect per step but of undefined sign; probably the intended line was `exp(ln10*radius(zone_idx+1)) - shell_radius`.

## net/net_lib.f90:1675-1680 -- rates: intent(out) json-length arrays, only element zone_idx assigned
- class: interface
- severity: low
- confidence: medium
- provenance: modernization (F77 rates.f wrote into COMMON arrays; no intent existed).
- detail: rate_pp(json) ... frac_be7_electron(json) are declared intent(out) but the routine assigns only index zone_idx. By the standard, intent(out) dummies become undefined on entry, so a compiler is entitled to leave the other json-1 elements undefined (allocatable/derived components would be deallocated). Callers (mix.f90:182 in a zone loop, timestep_limit_hburn.f90:99) pass the star%reaction_rate_* arrays and later read other elements. With gfortran this is benign today; with -finit-real or a different compiler it is a latent corruption of the rate arrays. Should be intent(inout).

## Weak/uncertain observations

- burn_lib.f90:2380-2384 (liburn): rate_li6/li7/be9(cz_base_zone) are modified in place (blend with zone-1) and never restored; in the rotation path they are then copied into rate_*_start for the next step (inherited liburn.f:271).
- burn_lib.f90:2585-2600 (liburn) when the CZ deepens (cz_base_zone < cz_base_zone_old): zones cz_base_zone..cz_base_zone_old-1 receive neither the burned CZ average nor a radiative burn (loops are 1..min_zone-1 and max_zone..num_zones); mix.f90 homogenises 15 species beforehand so only this step's burn is lost there (inherited `DO 130 I = MAXJ,M`).
- stitched_model.f90:508-509 / write_gyre_pulse.f90:111-112 (not my files): `star%eps_total(i)*star%pulse_dlneps_dlnt(i)` -- if pulse_dlneps_dlnt is filled from engeb's dimensional sum2/sum3 this is eps^2 * dlneps; if it is a true log-derivative it is fine. Cross-check for Sweep D.
- burn_lib.f90:1156/1190: `hydrogen_fraction.eq.0.0` float equality guarding rdeutmax/x; benign (X is exactly zero only when set so).
- burn_lib.f90 eqburn (~150-200): cn_ratio uses c12_p_rate in a denominator; net_lib rates zeroes any rate <= 1e-5 (2217) so a zone with log T just above nuclear_logT_cutoffs(3) but a tiny C12+p rate would divide by zero. Not reachable with the example cutoffs.
- burn_lib.f90:2483 (dburn): unconditional `write(*,913)` inside the accretion branch spams stdout every substep of every accreting model.
- burn_lib.f90:3062-3067 (lirate88): the tail loop after `exit` zeroes both rate_* and rate_*_start regardless of mode; harmless because the other array is zero there too, but it silently couples the two modes.
- burn_lib.f90 liburn: the argument is documented as "TIMESTEP IN YEARS" but callers pass star%dt in seconds and the rates are per second; code consistent, comment stale.
- burn_lib.f90:2313 (liburn): early return when the CZ falls below the burning threshold also skips the convective->radiative block that follows; matches the original GOTO but means Li in newly radiative shells is not burned on the step the CZ cools below threshold.
- The alpha-capture derivative defect (finding 2) cannot itself produce NaN: reaction_rate(8/10/11) cannot underflow to zero above realistic nuclear_logT_cutoffs(4), so the division is safe; the risk is a wrong Jacobian, not a FP exception.
- "UNABLE TO SOLVE FOR NEW ABUNDANCES IN SHELL 1" lives in src/mixing/solve_composition.f90 (outside this sweep); the inputs it receives from `rates` (with the <= 1e-5 -> 0 clamp at net_lib.f90:2217) are the only coupling seen from here. Nothing in rates/eqburn was found that would make the central-zone linear system singular.

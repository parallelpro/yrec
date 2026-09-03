# Batch 2 pre-fix review -- 2026-09-02

Read-only review of every Batch 2 item (ROADMAP sec. 11 "Revised fix
order") before any edit. Four independent reviewers (burn/net, eos,
kap, rotation/io/envelope) re-derived each claim from the code, the
F77 original (6cd5673) and the primary references; the main session
re-verified the three surprising findings (Itoh units, condopacpint
weighting, LVFC gate) directly. Nothing was edited, built or run.

Provenance: 13/14 inherited from the F77; only rhoofp06's tolerance
is modernization-era (1bc94ae, 2025-10-10). New this review: a
modernization REGRESSION in secular_transport's zahn_coupling_factor
gate (94c7f45) -- see item 9b.

Convention needed for items 1, 3: engeb's `dlnepsilon_dlnt/dlnrho`
are ABSOLUTE derivatives, sum3 = sum_i eps_i * dln r_i/dlnT (erg/g/s),
consumed as ql / ql_dt in henyey_coefficients:291-294. The
per-reaction `dlnrate_dlnt(i)` are dimensionless log-derivatives.

## 1. Alpha-capture Jacobian -- burn_lib.f90:1491-1540 (engeb.f:819-856)

Three defects, all inherited:
- prefactor: `density/reaction_rate(i)` where `1/S` is needed
  (S = the CF88 polynomial sum; rate = K rho e^U S, K ~ 1e22). The
  nuclear term is ~1e-22 x too small for reactions 8, 10, 11, so
  dlnrate_dlnt = dscreen_dlnt alone (~ -0.1..-0.5, wrong sign).
  True values at T9 = 0.1/0.2/0.3: 22.8/18.0/15.7 (C13an),
  21.5/16.6/14.3 (C12ag), 39.5/23.8/15.3 (N14ag).
- reaction 8: the resonant a2..a5 derivative terms omit `t9_m32`
  (present in the rate). <0.3% for T9 <= 0.3, grows above 0.5.
- reaction 10: `r1*cc13*0.1956` and `cc13*1.0616d0*r2` omit
  `t9_m23`. ~1.5% of the derivative at T9 = 0.1-0.3.
dlnrate_dlnrho (1 + dscreen) and triple-alpha are correct.
Verified by finite differences of the coded rate (scratchpad).

Fix (declare `s_sum`):
```fortran
! C13(ALPHA,N)O16   (replace 1499-1507)
      s_sum = a1*r1+t9_m32*(a2+a3+a4+a5)
      reaction_rate(8) = 1.157126d22*density*exp(screening_factor(8))*s_sum
      dlnrate_dlnrho(8)=1.0d0+dscreen_dlnrho(8)
      dr1 = cc13*(-2.0d0*t9_m23-0.0129d0*t9_m13+0.184d0*t9_p13)
      da1 = a1*(cc13*32.329d0*t9_m13 - 2.0d0*(t9/1.284d0)**2)
      dlnrate_dlnt(8) = dscreen_dlnt(8)+(dr1*a1 + r1*da1 + &
           t9_m32*(a2*(9.373*t9_m1-1.5d0)+a3*(11.873*t9_m1-1.5d0)+ &
           a4*(20.409*t9_m1-1.5d0)+a5*(29.283*t9_m1-1.5d0)))/s_sum
! C12(ALPHA,GAMMA)O16   (replace 1516-1523)
      s_sum = a1*(a2+a3)+a4+a5
      reaction_rate(10) = 1.25388d22*density*exp(screening_factor(10))*s_sum
      dlnrate_dlnrho(10) = 1.0d0+dscreen_dlnrho(10)
      dlnrate_dlnt(10) = dscreen_dlnt(10)+ &
           (a1*((cc13*32.120*t9_m13-2.0d0)*(a2+a3)+ &
           a2*(cc13*0.1956*t9_m23*r1-2.0d0*(t9/3.496)**2)+ &
           a3*(cc13*1.0616d0*t9_m23*r2))+ &
           a4*(27.499*t9_m1-1.5d0)+a5*(5.0d0+15.541*t9_m1))/s_sum
! N14(ALPHA,GAMMA)F18   (replace 1531-1540)
      s_sum = a1*r1+a2+a3+a4
      reaction_rate(11)= 1.07452d22*density*exp(screening_factor(11))*s_sum
      dlnrate_dlnrho(11)=1.+dscreen_dlnrho(11)
      dr1 = cc13*(-2.0d0*t9_m23-0.012d0*t9_m13+0.177d0*t9_p13+3.94d0*t9_p23)+0.406d0*t9
      da1 = a1*(cc13*36.031d0*t9_m13-2.0d0*(t9/0.881d0)**2)
      dlnrate_dlnt(11) = dscreen_dlnt(11)+(dr1*a1+r1*da1+ &
           a2*(2.798d0*t9_m1-1.5d0)+a3*(5.054d0*t9_m1-1.5d0)+ &
           a4*(12.310d0*t9_m1-cc23))/s_sum
```
Effect: logT >= cutoffs(4) = 7.73 only (core/shell He burning).
Jacobian only (converged models unchanged in principle) EXCEPT two
physical consumers: rot_scr%dlnepsilon_dlnt (rotation_stability_setup
:482) and the pulse/GYRE eps_T column. No pinned case reaches 7.73.
OPEN: reaction 11's r1 polynomial (0.177 T9^1/3, 3.94 T9^2/3) does not
obviously match CF88's N14(a,g) polynomial (0.012 T9^1/3, 1.45 T9^2/3,
0.117 T9, 1.97 T9^4/3, 0.406 T9^5/3) -- check the table before
touching r1/dr1 of reaction 11. Confidence: certain on the three
defects.

## 2. Deuterium d(p,g)3He exponent -- burn_lib.f90:1196 and :645

CF88: 2.24e3 T9^(-2/3) exp(-3.720 T9^-1/3)(1+0.112 T9^1/3+3.38 T9^2/3
+2.65 T9). Code uses `t9_p23` (+2/3). The code's own derivative
(qrtdeut, :1199) is exactly the -2/3 form, proving intent. Lifetime
check at logT 6.2, logrho 1.5: corrected 0.6 yr (first principles
0.54 yr), coded 6000 yr.
Fix: `t9_p23` -> `t9_m23` at :1196 (engeb), `t9p23` -> `t9m23` at
:645 (deutrate); both already computed.
Effect: rate x T9^-4/3 (x1e4 at 1e6 K, x260 at solar centre). Pre-MS
D-burning onset moves from ~1.7-2e6 K to ~1e6 K; every pinned case
starts from a *_Dbl.first model with D active, so ALL pins move
(early track; final solar calibration at round-off level);
expected_test_net.out deutrate line 1.57e5 -> ~8.5e8. Side note: the
two sites disagree by 0.7% in N_A/m_D (unrelated). Confidence: certain.

## 3. engeb Itoh branch -- burn_lib.f90:1952-1961

`dlnepsilon_dlnt = dlnepsilon_dlnrho + neutrino_dlnq_dlnt` is a
copy-paste (F77 `PET = PEP + DSNUDT`). But the intended one-liner is
still unit-inconsistent: lines 1953-1954 convert sneut's absolute
dsnu/dT, dsnu/drho into dimensionless log-derivatives, then add them
to the ABSOLUTE erg/g/s sums. The Beaudet branch (2020-2021) subtracts
absolute `qetnx = sum ex_k dln ex_k/dlnT`. Consistent form:
```fortran
          star%neutrino_loss_rate = -neutrino_loss_snu
          total_energy_gen_rate = total_energy_gen_rate + star%neutrino_loss_rate
! neutrino() returns absolute d snu/dT, d snu/d rho; dlnepsilon_* are
! absolute d eps/d ln T, d eps/d ln rho (same convention as the
! qetnx/qednx terms of the Beaudet branch).
          dlnepsilon_dlnrho = dlnepsilon_dlnrho - neutrino_density*neutrino_dlnq_dlnd
          dlnepsilon_dlnt   = dlnepsilon_dlnt   - neutrino_temp*neutrino_dlnq_dlnt
```
(rename the two locals; also removes the 0/0 when snu = 0.) Effect:
use_itoh_neutrino_loss + logT > cutoffs(5)=7.5 only; Jacobian plus
rot_scr%dlnepsilon_dlnt and pulse eps_T. Never fires in pinned solar
cases. Caveat: rotation_stability_setup:482-487 looks like it expects
a LOGARITHMIC dlneps/dlnT -- a units question for a later rotation
pass. Confidence: certain.

## 4. liburn/liburn2 radiative_frac -- burn_lib.f90:2547-2565, 2916-2934

Block runs only for cz_base_zone > cz_base_zone_old (retreating CZ),
so delta_mass = beg - end < 0 and radiative_frac in (-1, 0]. Intended
per the comment: FRAD = (end - m)/(end - beg), 1 at the old base, 0 at
the new. Coded exponent differs from the intended one by
ln(r_cz,start/r_rad,local) for every zone (a constant-factor error in
the local depletion rate, not a blend error). No clamp hides it.
Fix (both routines):
```fortran
      delta_mass = mass_coord_end - mass_coord_beg
      ...
         radiative_frac = (mass_coord_end - mass_coordinate(zone_idx))/delta_mass
```
Effect: Li6/Li7/Be9 in shells left behind by the retreating pre-MS CZ
(use_extended_composition) -- factors of several in the depletion
exponent; every pinned case with extended composition moves.
Confidence: certain.

## 5. sneut tfac2 -- net_lib.f90:342-343

`cvp*cvp - cap-cap` -> `cvp*cvp - cap*cap` (Itoh 1996 eq. 2.2-2.5;
Timmes sneut5 current; MESA neu/private/mod_neu.f90:39). tfac2
-1.31847 -> 0.18153, tfac3 -0.784 -> 0.108, tfac5 -0.659 -> 0.091.
Consumers: pair (a1 = tfac4(1+tfac3 qpair)), photo (a1 =
tfac4(1-tfac3 qphot): currently x1.52 instead of x0.93 at T <~ 1e8),
bremsstrahlung (z = tfac4 f - tfac5 g: sign of g flipped). Plasma and
recombination untouched. Effect: Itoh runs above logT 7.5 (RGB
degenerate cores, He cores, massive MS); expected_test_net.out sneut
lines 2-4 move by tens of %. Confidence: certain.

## 6. SCV smoothing weight -- eos/scv/eqscve.f90:148 and :202

Only eqscve has smoothing; eqscvg has NONE (ROADMAP:569 "eqscvg
same" is wrong -- correct it). Upward branch weight is
`0.5*d/tol` (0 at the node, 0.5 at d = tol) but the blend form
`own + w*(shift-own)` needs w = 0.5 at the node, 0 at d = tol:
```fortran
               temp_smooth_weight = 0.5d0*(1.0d0 - temp_dist_below/tol_temp_smooth)
               press_smooth_weight = 0.5d0*(1.0d0 - press_dist_below/tol_pressure_smooth)
```
Downward branch and the both-directions block need no change. The
spline stencils share knots so there is no jump AT the node; the
current jump is one step of 0.5|S_n+1 - S_n| at 60% of every cell.
Measured on h_tab_i.dat: T-direction up to 0.12 in dlnrho/dlnT (5%),
0.083 in dlnS/dlnT, 2.7% in grad_ad at logT 3.46-3.62, logP 4-4.4
(M-dwarf/BD photospheres); P-direction ~1e-3 dex.
Effect: use_scv_eos (all decks) wherever SCV is live: T < 1870 K,
the OPAL06 ramp band, or rhoofp06 fallbacks. Solar pins probably
unchanged; SCV-only and cool decks change. Confidence: 0.95.

## 7. rhoofp06 tolerance -- eos/opal/rhoofp06.f90:176

Restore 0.5d-7 (original, and rhoofp/rhoofp01). Mechanism verified:
trial calls pass deriv_order=1; esac06's tail rescaled the cv slot
(x ~1e2) from stale data on EVERY trial call; batch 0 gated the
tail on eos_index_inverse(slot) <= deriv_order, so nothing in the
inversion touches slot 5 now. The corruption only persisted across
inversions that hit the 11-refinement cap (silent -999 -> oeqos06
return 1 -> eqstat *998 -> Yale/SCV values, ierr=0, NO message);
the tight tolerance made the cap hit more often (plain regula falsi
is linear once an endpoint sticks; the table-edge clamp of the upper
bracket makes it hopeless), so loosening it hid the symptom.
Fix: tolerance line + comment; header lines 9-18 drop item (3);
rhoofp01.f90:11-13 cross-reference. Recommended in the same change:
raise the cap 11 -> ~30, re-enable the non-convergence write to the
run log, and verify the bracket after the x0.2/x5 start (audit note:
regula falsi may start unbracketed). Effect: every OPAL06 deck (311
inlists, all pins): rho(P,T) shifts up to ~1e-5 relative; all pins +
expected_test_eos.out reseed. Confidence: 0.9 that the crash cannot
return; 0.6 operational as a one-liner (hence the cap/diagnostic).
Side note: the once-per-process priming call takes radsub06 on
zero-initialised slots (0/0 at radsub06.f90:68) -- would trap under
-ffpe-trap=invalid; out of scope.

## 8. condopacpint derivatives -- kap/conductive/condopacpint.f90:114-120

Two defects, both inherited: (a) 119-120 feed dlnkap_dlnrho_* into
the T-derivative (dlnkap_dlnt_h1/he4/ox are computed and never read);
(b) BOTH derivatives are multiplied by mix_log10_cond (~12-17)
instead of divided by the summed conductivity. With
mix = -log10(sum w_i 10^-K_i), d mix/dx = sum w_i cond_i dK_i/dx /
sum w_i cond_i. At logT 7, logrho 3-5: correct ~ -2.7 (rho), +2.3
(T); live ~ -25 and -21 (T sign flipped).
```fortran
      conductive_dlnkap_dlnrho = (weight_h1*cond_h1*dlnkap_dlnrho_h1 + weight_he4*cond_he4*dlnkap_dlnrho_he4 &
                                  + weight_ox*cond_ox*dlnkap_dlnrho_ox) &
                                 / (weight_h1*cond_h1 + weight_he4*cond_he4 + weight_ox*cond_ox)
      conductive_dlnkap_dlnt   = (weight_h1*cond_h1*dlnkap_dlnt_h1 + weight_he4*cond_he4*dlnkap_dlnt_he4 &
                                  + weight_ox*cond_ox*dlnkap_dlnt_ox) &
                                 / (weight_h1*cond_h1 + weight_he4*cond_he4 + weight_ox*cond_ox)
```
Not covered: the extrapolation branch (130-158, logrho < -6 with
logR >= 0, effectively unreachable) has its own sign/pure-H issues.
Effect: kappa_c value unchanged; derivatives enter the blend with
weight kr/(kr+kc) (~5e-4 in solar/PMS -> pins numerically unchanged,
maybe iteration counts) and dominate in degenerate cores; PHYSICAL
consumers: rotation_stability_setup:474 qchit (rotating RGB
diffusion coefficients) and the pulse kap_T column. test_kap runs
without conductive opacity. Confidence: 0.9+.

## 9. ll95tbl slot collision -- kap/opacity_table_lib.f90:44,177

start_index = [0,10,..,100,109,118] assumes a packed layout but the
readers/consumers add the UNPACKED ix. Row 11 (Z=0.06) ix=10 -> slot
110 = (X=0, Z=0.08); row 12 (Z=0.08) ix=10 -> 119 = (X=0, Z=0.10).
The ix=10 tables are read last, so the X=0 tables at Z=0.08/0.10 are
overwritten by X=0.94/Z=0.06 and X=0.92/Z=0.08 (GS98.OP17 index:
#12, #13, #125, #126). Fix (option B, no consumer changes):
`n_opal95_xz = 130`, `start_index = [0,10,20,...,100,110,120]`; keep
ll95tbl's num_xz = 126 (table count). Unwritten slots 109/119/129/130
are unreachable by the stencils. Compile to confirm no hard-coded
126 dimensioning survives.
Effect: 4D path with X < 0.1 and Z-stencil reaching >= 0.08 (He cores
with C+O >~ 0.04) and the HB branch (Z > 0.1, logT > 7 -> getopal95 at
Z=0.1: returns table #126 instead of (X=0, Z=0.10), ~ +0.28 dex before
the (1-Z)/0.9 blend). Solar/MS/RGB untouched (Z <= 0.02). No pins
move; all TAHB/ZAHB/HeIgnite cases move physically. Confidence: 0.95.

## 9b. NEW: zahn_coupling_factor gate regression -- secular_transport.f90:196,239

F77 seculr.f:108,147 `IF(LVFC) CALL GETFC`; modernized (94c7f45)
`if(star%ctrl%use_diffusion_advection_transport) call
zahn_coupling_factor`. diffusion_velocity_scales.f90:74-93 selects on
lvfc and multiplies the mixing coefficient by star%vfc, which is
zeroed at secular_transport:184-186 -> velocity_form_correction=T
with use_diffusion_advection_transport=F gives ZERO rotational mixing
silently. Fix: `if(star%ctrl%lvfc)` at both sites. No pin sets lvfc.
Confidence: certain.

## 10. zahn_coupling_factor alpha/2 -- zahn_coupling_factor.f90:62-64

Zahn 1992: nu_h = (r/C_h)|2V - alpha U|, alpha = 1/2 dln(r^2 Omega)/
dln r; the code's circ_velocity already holds 2V (comment "DETERMINE
2V"), and the coded alpha is dln(r^2 Omega)/dln r = 2 alpha_Zahn, so
the consumer forms |2V - 2 alpha U|. No consistent reinterpretation
rescues it; the header states the 1/2. Fix: prefix `0.5D0*` at :62.
Effect: lvfc=T rotating runs only; O(1) in vfc (4x for the U=const
solid-body illustration, then capped at 1). No pins. Caveat: any
historically tuned FC (= C_h/30 in this mode) absorbed the wrong
alpha; LVFC mode is marked "Not recommended" in the templates.

## 11. compute_quadrupole drho/dr -- compute_quadrupole.f90:95

`-weight_minus*exp(ln10*log_density(zone_index))` should reference
`zone_index-1` (sibling stencil at 101-103 is correct). Current
expression collapses to 1/2 of a one-sided forward difference.
Fix: `log_density(zone_index-1)` at :95. Effect: instability-
transport rotating runs; 4piG rho'/g term is ~0.3-0.4 of the 6/r^2
term in the solar radiative interior -> quadrupole ~10-20%, ES
correction ratio a few-10%, Omega/Li at percent level. Pins:
cm_rot_base/disk/kawaler/settle_he/solid, difrotmix x3, GS_rot,
3Msolar I0/I3cz2. Header says 1/R^4 where the code uses 1/R^3 --
check separately. Confidence: high.

## 12. trapzd offset -- numerics_lib.f90:1728-1731

`rhot = rhop+drho*del` (and w2t, eta22t) use the constant del where
(y-b1) is needed (smt already uses y^2-b1^2). With jmax=2 the single
midpoint takes the endpoint values. Fix: `*(y-b1)` on the three
lines. Effect: rotating runs only; aint biased ~1-3% -> f_P by ~1e-8
(Sun) to 1e-6 (rapid rotators). All rotating pins move at bit level.
Confidence: high that it is a bug, low that it matters.

## 13. envelope_derivs dydx(3) -- core/envelope_derivs.f90:95

Kippenhahn-Thomas/Endal-Sofia: dP/dM_P = -G M f_P/(4 pi r^4), dr/dM =
1/(4 pi r^2 rho); chain rule with x = logP gives 1/f_P on BOTH
dlnm/dlnP and dlnr/dlnP. dydx(1) divides, dydx(3) multiplies ->
mass continuity violated by f_P^-2. Interior (henyey_coefficients
:221) and atmosphere (g_eff = g f_P) agree with 1/f_P. Fix:
`*pressure_rotation_factor` -> `/pressure_rotation_factor` at :95.
Effect: fp_rot = 1.0d0 exactly for non-rotating (read_starting_model
:466, henyey_iterate:164) -> non-rotating pins bit-identical;
rotating: 1-f_P ~ 1.4e-5 (Sun) -> R, Teff ~1e-6; PMS/rapid ~1e-4.
Confidence: high.

## 14. model_to_equal N15-as-Z -- io/model_to_equal.f90:125-127,233-236,255,316-319,329

composition(8) (N15) -> composition(i_metals) at all five sites;
equal_to_model applies the change back to slot 3. LATENT: on the old
settling path (use_diffusion_z + use_new_diffusion_routines=F) the
mid array is overwritten by lax_wendroff_step1 before its only read,
and metal_abundance_change is used only as EZ - EZ_orig (cancels to
round-off). All 13 pins with use_diffusion_z use the new routines.
Fix for hygiene; no pin moves. Confidence: high.

## Pin impact summary

| item | pins that move |
|---|---|
| 1 alpha-capture | none (logT >= 7.73 unreached) |
| 2 deuterium | ALL (early pre-MS), expected_test_net.out |
| 3 Itoh units | none |
| 4 liburn | all with use_extended_composition (testsuite solar, matrix, 0.3 Msun examples) |
| 5 sneut | expected_test_net.out |
| 6 SCV weight | probably none of the solar pins; SCV-only/cool decks |
| 7 rhoofp06 | ALL OPAL06 pins + expected_test_eos.out |
| 8 condopacpint | none numerically (maybe iters); rotating RGB physically |
| 9 ll95tbl | none; TAHB/ZAHB/HeIgnite physically |
| 9b LVFC gate | none |
| 10 alpha/2 | none |
| 11 quadrupole | all instability-transport rotating pins |
| 12 trapzd | all rotating (bit level) |
| 13 envelope fp | all rotating |
| 14 N15 | none |

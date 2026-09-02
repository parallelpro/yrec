# F: seculr + shape (secular AM/composition transport, rotational shape factors)

Second independent pass, 2026-09-01. Read-only; no files edited, nothing built or run.

## Files read in full

rotation/seculr/ (all 16):
- /Applications/YREC/src/rotation/seculr/am_advection_diffusion_coeffs.f90 (697 lines) vs dadcoeft.f
- /Applications/YREC/src/rotation/seculr/am_diffusion_coeffs.f90 vs dcoeft.f
- /Applications/YREC/src/rotation/seculr/am_transport_grid.f90 vs rotgrid.f
- /Applications/YREC/src/rotation/seculr/banded_solver.f90 vs bandw.f
- /Applications/YREC/src/rotation/seculr/check_angular_momentum.f90 vs checkj.f
- /Applications/YREC/src/rotation/seculr/check_composition.f90 vs checkc.f
- /Applications/YREC/src/rotation/seculr/circulation_velocities.f90 vs vcirc.f
- /Applications/YREC/src/rotation/seculr/composition_diffusion_coeffs.f90 vs ccoeft.f
- /Applications/YREC/src/rotation/seculr/composition_grid.f90 vs mixgrid.f
- /Applications/YREC/src/rotation/seculr/compute_quadrupole.f90 vs getqua.f
- /Applications/YREC/src/rotation/seculr/diffuse_composition.f90 vs mixcom.f
- /Applications/YREC/src/rotation/seculr/diffuse_composition_driver.f90 vs ndifcom.f
- /Applications/YREC/src/rotation/seculr/diffusion_velocity_scales.f90 vs codiff.f
- /Applications/YREC/src/rotation/seculr/equal_grid_to_model.f90 vs eq2mod.f
- /Applications/YREC/src/rotation/seculr/secular_transport.f90 vs seculr.f
- /Applications/YREC/src/rotation/seculr/zahn_coupling_factor.f90 vs getfc.f

rotation/shape/ (all 5):
- /Applications/YREC/src/rotation/shape/equipotential_integrand.f90 vs func.f
- /Applications/YREC/src/rotation/shape/rotation_shape_factors.f90 vs fpft.f
- /Applications/YREC/src/rotation/shape/shape.f90 vs shape.f
- /Applications/YREC/src/rotation/shape/shell_inertia_integral.f90 vs intmom.f
- /Applications/YREC/src/rotation/shape/zone_moments_of_inertia.f90 vs momi.f

Supporting code read for interface/semantic checks (not audited in full unless noted): numerics_lib.f90 tridia/ctridi (1188-1327), trapzd/qgauss/polint (819-872, 1660-1817, read in full); rotation/evolve_angular_momentum.f90 (60-95, 228-380); rotation/mid_timestep_model.f90 (350-466, zone_min/zone_max semantics); rotation/solid_body_omega.f90 (70-104); mixing/burn_settle_mix.f90 (115-140); rotation_scratch_lib.f90 declarations. F77 ancestors taken from `git show 6cd5673:src/<name>.f` for all of the above plus trapzd.f, getw.f, setupv.f, tridia.f, qgauss.f, model_to_equal.f, equal_to_model.f.

Interface check: every bare external in the assigned files was grepped with `grep -rn "call <name>" /Applications/YREC/src` (no head) and actual vs dummy lists compared: tridia, ctridi, osplin, equal_spaced_grid, homogenize_convection_zones, solid_body_omega, matt_wind, wind_spindown_matt, rotation_stability_setup, compute_quadrupole, circulation_velocities, zahn_coupling_factor, diffusion_velocity_scales, am_transport_grid, equal_grid_to_model, diffuse_composition_driver (both call sites incl. burn_settle_mix.f90:127), diffuse_composition, composition_grid, composition_diffusion_coeffs, am_diffusion_coeffs, am_advection_diffusion_coeffs, banded_solver, check_angular_momentum, check_composition, shape, zone_moments_of_inertia, shell_inertia_integral, trapzd, qgauss, polint, omega_from_j. All argument counts/orders/kinds match.

## Verified clean (transliteration matches F77 and the numerics/physics checked out)

- banded_solver.f90: identical to bandw.f; band convention A(i,j) = element (i, i+j-5), forward elimination fill-in stays within the 10-column band, back-substitution index arithmetic correct; residual-check loops in am_advection_diffusion_coeffs.f90:465-509 use the same convention.
- composition_diffusion_coeffs.f90: conservative form Δm_i(X^{n+1}-X^n) = (dt/dχ)[ECOD2_{i+1}(X_{i+1}-X_i) - ECOD2_i(X_i-X_{i-1})], zero-flux boundary rows, ECOD2 at interfaces (i, i-1). Identical to ccoeft.f.
- am_diffusion_coeffs.f90: AM analogue with I_i weights and Jacobian ρ r² (I/M) dχ/dr; wind-loss term seeded through tridia dj_n_seed consistently. Identical to dcoeft.f.
- am_transport_grid.f90 (rotgrid) and equal_grid_to_model.f90 (eq2mod): grid build, Jacobians (AM: ρr²(I/M)dχ/dr; mixing: ρr²dχ/dr), lumped-CZ end masses, fractional dJ/J interpolation back plus global ratio enforcing ΣΔJ = SUMDJ. Both directions conserve AM to the tolerance by construction; both branches of the I0/I1 CZ search present in eq2mod.
- composition_grid.f90, diffuse_composition_driver.f90: GOTO-elimination equivalent to mixgrid.f/ndifcom.f.
- circulation_velocities.f90, diffusion_velocity_scales.f90: identical to vcirc.f/codiff.f (ZM98 mu-inhibition, GSF, Endal-Sofia shear, Zahn coupling-factor branches, LNOJ zeroing).
- compute_quadrupole.f90: tridiagonal operator, Thomas solve, surface BC identical to getqua.f (but see finding 1 for the density derivative).
- check_composition.f90: identical to checkc.f; mean molecular weight expression consistent.
- shape.f90: Radau equation r dη/dr + 6(ρ/ρ̄)(η+1) + η(η-1) = 6 and the RK4/R0 iteration identical to shape.f; A = (5/3)ω²r0³/(GM(2+η)) reproduces the Roche limit for η = 3 (see finding 11 for the central seed).
- rotation_shape_factors.f90, equipotential_integrand.f90: identical to fpft.f/func.f. Checked that Φ_dis = (4π/3) aint/r³ with aint = ∫(ρ/m)ω²r0^7(5+η)/(2+η)dr0 is the interior Clairaut quadrupole potential with G cancelled by the ε ∝ 1/GM definition, that dΦ/dr and dΦ/dθ in func are the consistent derivatives, and that the area element 4π r sinθ r0 sqrt((1-aP2)²+(3a cs)²) over [0, π/2] is the correct surface measure. (Comments in func quote 12πG/5 and 4πG/5 factors that correspond to a different normalisation of AINT; the code is self-consistent.)
- zone_moments_of_inertia.f90, shell_inertia_integral.f90: identical to momi.f/intmom.f; rotation parameter matches shape's A.
- secular_transport.f90 region loop, wind/solid-body branches, composition passes (species 1-4 inside the iteration loop, 5-11/15 after), CZ-flag save/restore: equivalent to seculr.f.

## Findings

## /Applications/YREC/src/rotation/seculr/diffuse_composition.f90:158-166 -- i1 undefined when the unstable region's top is radiative (or zone_end == nz): He4 renormalisation silently skipped / out-of-bounds writes
- class: logical / numerical (uninitialized)
- severity: high
- confidence: high (defect), medium (frequency)
- provenance: inherited (mixcom.f:110-119 has the same missing ELSE; F77 relied on SAVEd stale I1), modernization changed the failure mode
- detail: The "CHECK FOR UPPER CONVECTION ZONE" block sets i1 only inside `if (convective_flag(zone_end).and.zone_end.lt.num_zones)`; the lower block has an `else i0 = zone_begin` (line 155) but there is no `else i1 = zone_end`. A region whose top zone is radiative (velocity fell below 1e-20, mu-barrier zeroing, no_am_transport_in_core splitting) or a region ending at nz leaves i1 undefined. i1 then bounds the species-mass sum (170), the outer-CZ copy `if (zone_end.lt.i1)` (198-201) and the He4 closure loop (212-216). F77 (SAVE) reused I1 from a previous call (e.g. M from a region that touched the surface CZ), so `IF(IEND.LT.I1)` copied HCOMP(ID,IEND) into IEND+1..M, overwriting the composition above the region; in the modern local (no SAVE, no init) GNU -finit-local-zero gives i1=0 so all three loops are empty and He4 is not recomputed as 1-X-Z-He3 after X/He3/Z were diffused, leaving the mass fractions not summing to 1 for that region (the pre-final pass whose result is kept when the iteration converges), while ifort/ifx gets garbage i1 and writes composition(:,zone_end+1..garbage). Observable: silent non-conservation of ΣX_i in rotating models with mu-barriered/radiative-topped unstable regions; Intel builds can corrupt memory.

## /Applications/YREC/src/rotation/seculr/compute_quadrupole.f90:92-95 -- interior dρ/dr uses ρ(i) twice instead of ρ(i-1): density-derivative term of the Sweet/Zahn quadrupole operator halved
- class: physical
- severity: medium
- confidence: high
- provenance: inherited (getqua.f:56-57 `DRHO = 0.5D0*(FPL*EXP(CLN*HD(I+1))+(FMI-FPL)*EXP(CLN*HD(I))-FMI*EXP(CLN*HD(I)))/DR`)
- detail: The centred derivative template used everywhere else in the routine (e.g. DRHOW2 at lines 62-64 of getqua.f and the modern equivalent) is 0.5[FPL f(i+1) + (FMI-FPL) f(i) - FMI f(i-1)]/DR. For the density it takes f(i-1) -> ρ(i) so the expression collapses to 0.5·FPL·(ρ(i+1)-ρ(i))/DR, i.e. half of a one-sided forward difference rather than the centred difference. This enters `diag(i) = ... - four_pi_g*drho_dr/gravity(i)`, the 4πG(dρ/dr)/g term of the quadrupole ODE, so the quadrupole (and hence the Eddington-Sweet velocity correction in circulation_velocities) is computed with a systematically underestimated density-gradient term throughout the interior. Not a crash; a biased ES circulation velocity.

## /Applications/YREC/src/rotation/seculr/zahn_coupling_factor.f90:61-66 -- alpha missing the factor 1/2 of Zahn (1992) / the header's own definition
- class: physical
- severity: medium
- confidence: high
- provenance: inherited (getfc.f:44-45)
- detail: Header (and getfc.f:30) defines f_c = C/30 · U/|2V - αU| with α = ½ d ln(r²Ω)/d ln r. The code computes half_dlnj_dlnr(i) = [Ω_i r_i² - Ω_{i-1} r_{i-1}²]/(Δr) /(Ω_mid r_mid) = d(r²Ω)/dr /(rΩ) = d ln(r²Ω)/d ln r, i.e. 2α, despite the variable being named half_dlnj_dlnr. The denominator |V - αU| at line ~103 (`abs(v - half_dlnj_dlnr*hv)`) therefore weights the advective term twice too strongly, reducing the variable coupling factor VFC (ratio of mixing to AM-transport diffusion coefficient) wherever the rotation profile is steep. Affects only the use_variable_coupling (LVFC) option.

## /Applications/YREC/src/rotation/seculr/am_advection_diffusion_coeffs.f90:552-559 -- "CORRECTIONS TOO LARGE" now exits only the innermost loop; F77 GOTO 950 aborted the whole substep sequence and cut the timestep
- class: logical (modernization control-flow divergence)
- severity: medium
- confidence: high
- provenance: modernization (dadcoeft.f:445-451 `LOKAD=.FALSE.; GOTO 950` where label 950 is outside all three DO loops; modern line 559 is an unlabelled `exit` inside `do coeff_iter_idx`)
- detail: In dadcoeft.f a >10% omega correction jumps straight to the timestep-cut block (NCUT+1, NSTEP doubled, restart from label 5). The modern `exit` leaves only the coeff_iter loop; execution continues at line 598 into the theta-iteration bookkeeping, the theta loop keeps iterating (re-solving from the unchanged omega_working and hitting the same rejection each time until max_diffusion_iters), then the substep is closed at 623-632 (omega_substep_start updated, wind_loss_delta added to total_angular_momentum_start) and the next substep begins. Because diffusion_converged is reset to .false. at the top of every theta iteration (186) and set .true. by any later converged substep (594), a rejected substep is forgotten if the final substep converges: that substep's AM transport is silently skipped (its wind torque is still applied through the conservation rescale) and no timestep cut happens. Only affects the optional use_diffusion_advection_transport (LDIFAD) path.

## /Applications/YREC/src/rotation/seculr/am_advection_diffusion_coeffs.f90:249-252, 281, 683-686 -- shear/GSF equal-grid diffusion coefficients are never assigned: secular-shear and GSF mixing silently dropped from the LDIFAD mixing coefficient
- class: logical (dead data path)
- severity: medium
- confidence: high
- provenance: inherited (dadcoeft.f:187-188, 215, 550: ESS/EGSF in COMMON/DIFMIX/ are only ever multiplied in place; `git grep -n "ESS(\|EGSF(" 6cd5673 -- src` shows no assignment anywhere) 
- detail: rot_scr%shear_diffusion_coeff_eqgrid / gsf_diffusion_coeff_eqgrid are declared in rotation_scratch_lib.f90:60-61 and referenced only here (`grep -rn` over src confirms). They are never filled from the model-grid rot_scr%shear_diffusion_coeff / gsf_diffusion_coeff that circulation_velocities.f90:524-527 computes, so they are identically zero. Consequences: (a) the in-place rescaling at 249-252 is a no-op; (b) line 281 adds zero to the diffusive coefficient (AM transport still gets SS/GSF via es_diffusive_velocity, so only partially affected); (c) mixing_diffusion_raw at 683-686 contains only the ES estimate, so with LDIFAD the composition mixing coefficient ignores secular shear and GSF entirely, regardless of secular_shear_mixing_scale / gsf_mixing_scale. A user enabling those scales under LDIFAD gets no effect and no warning. Even if the arrays were populated, the in-place multiplication by (domega_dr_prev_substep/domega_dr(i))² compounds across coeff_iter, theta_iter, substeps and cut retries because domega_dr is frozen at the start of the call (latent second defect).

## /Applications/YREC/src/rotation/evolve_angular_momentum.f90:88,327,335 + /Applications/YREC/src/rotation/seculr/secular_transport.f90:307 -- redo_flag (intent(out)) left unset on the "no unstable zone" early return; caller reads an uninitialized local
- class: interface / uninitialized
- severity: medium
- confidence: high
- provenance: inherited (seculr.f:197 `IF(.NOT.LDO) GOTO 9999` leaves LREDO untouched; getw.f:83 SAVE, LREDO never initialised before the CALL SECULR at getw.f:273), modernization changed stale-SAVE into uninitialized-local
- detail: secular_transport declares `logical, intent(out) :: redo_flag` (line 72) and `double precision, intent(out) :: cz_moment_of_inertia` (74) but returns at 307 before either is necessarily assigned (cz_moment_of_inertia is assigned at 274 only inside the wind-loss branch). evolve_angular_momentum.f90 never initialises redo_needed_flag (plain local, no SAVE) and tests it at 335 to decide whether to roll back the composition/J arrays and `cycle retry`. In F77 the SAVEd LREDO could retain .TRUE. from an earlier timestep-cut, producing a spurious second cut when the retried step found no instability; in the modern build GNU -finit-local-zero masks it (.false.) but ifort/ifx reads garbage, so a fully stable rotating model can enter the retry/rollback path at random. Also note that when secular_transport is entered with wind loss disabled the caller's moment_of_inertia_cz is only meaningful because mid_timestep_model happens to have filled it earlier.

## /Applications/YREC/src/rotation/seculr/secular_transport.f90:146-147,334 -- constant_diffusion_coeff_flag/constant_diffusion_coeff read but never assigned (uninitialized on Intel)
- class: uninitialized / dead branch
- severity: medium
- confidence: medium
- provenance: inherited (seculr.f uses LCODM/CODM without declaring COMMON/MAG/, which codiff.f does declare -- the F77 statics were zero-initialised implicit SAVE locals, so the branch was dead); modernization kept them as plain uninitialized locals
- detail: The header comment at 140-145 acknowledges they are "permanently at their (undefined/zero-initialized) default". On GNU with -finit-local-zero the `.and.` short-circuits to .false. and the branch is dead as in F77. On ifort/ifx a logical local is not zero-initialised, so `if(constant_diffusion_coeff_flag.and.constant_diffusion_coeff.gt.0.0D0)` can be true, forcing `zone_begin = zone_min - 1, zone_end = zone_max, scan_start_zone = zone_max+1`, i.e. treating the whole radiative interior as one unstable region with whatever velocities exist. The intended control (LCODM/CODM, now star%ctrl constant-D option in diffusion_velocity_scales) is simply not wired here, so the "constant background diffusion coefficient treats entire domain as unstable" feature documented in the MHP 8/13 comment does not work even when the user enables it.

## /Applications/YREC/src/numerics/numerics_lib.f90:1729-1732 (trapzd, called from rotation_shape_factors.f90:96) -- midpoint refinement interpolates ρ, ω², η+2 at b1+del instead of at the midpoint y
- class: numerical
- severity: low
- confidence: high
- provenance: inherited (trapzd.f:27-30 `RHOT = RHOP+DRHO*DEL`, `W2T = W2P+DW2*DEL`, `ETA22T = ETA22P+DETA2*DEL`)
- detail: In the n>1 branch the loop variable y advances (y = y+del) and the mass is correctly interpolated with (y²-b1²), but density, ω² and η+2 use the constant offset DEL, i.e. the value at the end of the first sub-interval for every midpoint. With rotation_shape_factors' jmax=2 there is exactly one refinement (it=1, del=dr), so the single midpoint sample uses ρ, ω², η of shell i (the upper end) rather than the shell-midpoint values; the subsequent 2-point Richardson extrapolation (Simpson weights 1:4:1) therefore carries a first-order error ∝ (ρ_i-ρ_{i-1})/ρ in each shell's contribution to AINT, the quadrupole potential integral that feeds <g>, <g⁻¹>, FP and FT. The term is small (rotational distortion) so the structural effect is modest, but the "Richardson-extrapolated" integral is not what the code claims.

## /Applications/YREC/src/rotation/seculr/check_angular_momentum.f90:207 -- omega-gradient-reversal guard uses threshold 1.0 rad/s: dead for any real star; header still advertises it
- class: logical / comment-code mismatch
- severity: low
- confidence: medium (that it is unintended; high that it is inert)
- provenance: inherited (checkj.f:81 `IF(OMEGA(I)-OMEGA(I-1).GT.1.0d0)` -- note the lower-case d0 differs from the 1.0D-15 used at checkj.f:110,124, suggesting a later hand-edit)
- detail: The first-pass reversal detector requires ΔΩ > 1 rad/s between adjacent zones; stellar Ω is 1e-6..1e-3 rad/s so it never fires, hence the solid-body repair (zone_scan block 205-288, including the 1e-15 refinement tests at 242/257) is unreachable. The header ("SECOND, IT GUARDS AGAINST REVERSAL OF ANGULAR VELOCITY GRADIENTS") and the OMEGA GRADIENT REVERSAL message are misleading. If this was a deliberate disable it should be a control; if not, spurious positive gradients from the diffusion solve are no longer repaired.

## /Applications/YREC/src/rotation/seculr/check_angular_momentum.f90:104-105,158 and am_advection_diffusion_coeffs.f90:106-107,590 -- per-iteration history arrays sized 16/50 indexed by user control max_diffusion_iters with no bound check
- class: numerical (array bound)
- severity: low
- confidence: high
- provenance: inherited (checkj.f DELJI(16)/IMAXI(16); dadcoeft.f DWIT(50)/IWIT(50))
- detail: max_delta_j_by_iter(iteration_number) is written for iteration_number = 1..star%ctrl%max_diffusion_iters (itdif2, namelist integer, default 1, no validation in read_controls.f90/inlist_new_read.inc). A user value > 16 (or > 50 for the LDIFAD history arrays) writes past the local array and corrupts the stack silently under the default build (no -fcheck=bounds).

## /Applications/YREC/src/rotation/shape/shape.f90:60 -- central seed η₂(1) = 6(1-ρ/ρ̄) is 7x the Radau small-r limit
- class: physical
- severity: low
- confidence: medium
- provenance: inherited (shape.f:33)
- detail: Near the centre with η ∝ r² the Radau equation r dη/dr + 6(ρ/ρ̄)(η+1) + η(η-1) = 6 gives 7η ≈ 6(1-ρ/ρ̄), i.e. η₂ → (6/7)(1-ρ/ρ̄); the code uses 6(1-ρ/ρ̄) (even the crude dη/dr = 0 assumption gives 6/5). The error is damped by the RK4 integration (perturbations decay roughly as r⁻⁵ in ln r) so only the first few zones' η₂, r0 and FP/FT carry it; for coarse central zoning or when shape() is called on a sub-range starting at zone 1 (solid_body_omega) it is visible in the innermost zones' moment of inertia.

## /Applications/YREC/src/rotation/shape/shell_inertia_integral.f90:51 + /Applications/YREC/src/rotation/solid_body_omega.f90:93-94 -- "dI/dω" is actually a·dI/da = (ω/2)dI/dω; Newton denominator for the solid-body ω solve is dimensionally inconsistent
- class: physical / comment-code mismatch
- severity: low
- confidence: high
- provenance: inherited (intmom.f:18-21; momi.f header "QIW = DI/D(OMEGA)")
- detail: series_sum_domega = Σ j c_j (jη+5) a^j is a·d(series)/da; with a ∝ ω² the true dI/dω = (2/ω)·di_domega_per_mass. solid_body_omega uses delta_omega = ΔJ/(I_tot + ω·ΣQIW), which should be ΔJ/(I_tot + ω dI/dω) = ΔJ/(I_tot + 2ΣQIW). Since ω ~ 1e-5 the correction term is effectively dropped rather than doubled, so the Newton step is only O(a) wrong and still converges within the 20-iteration cap for slow rotators; near break-up (a ~ 0.1) the iteration converges more slowly or hits the cap and the message-free fall-through at solid_body_omega.f90:91-100 accepts an unconverged ω.

## /Applications/YREC/src/rotation/seculr/diffuse_composition.f90:60-64 vs 208-216 -- He4 closure applied over i0..i1 but species diffused only over zone_begin..zone_end plus copied CZs; species-conservation correction is commented out in both versions
- class: numerical
- severity: low
- confidence: high
- provenance: inherited (mixcom.f:150-159 commented-out RATIO block)
- detail: The equal-grid solve conserves Σ EM·X exactly (dcomp), but the osculatory-spline map back to the model grid (dcomp2) does not, and the F77 authors' correction loop is commented out in both. Unlike the AM path (equal_grid_to_model applies the global RATIO), composition therefore drifts by the spline's interpolation error every diffusion substep; for species 1-4 the He4 closure hides the drift by absorbing it into Y, for species 5-11 it accumulates in the metals. Worth quantifying before trusting surface Li/Be depletion at the percent level.

## Weak/uncertain observations

- diffuse_composition_driver.f90:42 (and other headers) describe zone_max as "the outermost radiative zone"; mid_timestep_model.f90:436-466 actually sets it to the lowest convective shell of the surface CZ (i+1). Comment/code mismatch only, but it is the reason convective_flag(zone_end) is normally true and finding 1 is latent rather than universal.
- rot_scr%log_luminosity_mid is linear L/Lsun (HL convention), not a log; misnomer only.
- am_transport_grid.f90 lumps the lower CZ into the first equal-grid point only when convective_flag(zone_begin); equal_grid_to_model does the same, so consistent, but composition_grid.f90's surface lumping loop runs regardless of convective_flag(zone_end) (harmless because it stops at the first radiative zone).
- composition_grid.f90 ~ lines 60-70: `i0 = idx + 1`, `i1 = idx - 1` use the loop variable after loop exit (stale), exactly as mixgrid.f used stale I; i0/i1 are never read.
- circulation_velocities.f90:355,362,411-413 take log(omega(i)) without guarding omega <= 0; omega is protected upstream by check_angular_momentum's negative-J test, but a zero-omega initial model would produce -Inf here. `wmid` (403,408) is computed and never used.
- compute_quadrupole.f90:24 says the surface solution varies as 1/R**4 while line 109 and the BC coefficient FX=(R/R+dr)**3 implement 1/R**3; and the row-1 derivative (getqua.f:41) is a forward difference while the header claims centred derivatives. Cosmetic.
- am_advection_diffusion_coeffs.f90:191-204: lrossby/pmmsoltau are uninitialized locals (header admits it); the Rossby-scaled saturation branch is therefore build-dependent (dead on GNU, random on Intel). Same class as finding 7, LDIFAD only.
- am_advection_diffusion_coeffs.f90:553-555, 579-580, 589, 598-603, 653, 676-679, 695: unconditional write(*,*) diagnostics on every iteration of the production LDIFAD path (F77 had them too); floods stdout and slows the solve.
- am_advection_diffusion_coeffs.f90:634-635: empty `if (substep_idx > num_substeps) then; end if` left by the GOTO conversion.
- check_angular_momentum.f90:131 `if(cut_count.gt.0)` with message "3 ATTEMPTS AT CUTTING TIMESTEP FAILED": the first negative J/M stops the run (F77 identical; IREDO.GT.3 commented out). Message is misleading.
- check_composition.f90:78: species 5-11 are range-checked only when iteration_number == max_diffusion_iters, but secular_transport diffuses them after the loop exits early on convergence, so the CNO/Li pass is normally never checked for negative abundances (inherited).
- numerics_lib.f90:1318-1319 (tridia): the `fj` correction is a no-op (out of scope, noted for H sweep).
- rotation_shape_factors.f90 / equipotential_integrand.f90: only the interior Clairaut term (∫_0^r) of the quadrupole potential is included; the exterior-shell term r²∫_r^R ρ dε/dr' dr' of the full Clairaut potential is absent. I could not confirm from the papers whether Endal & Sofia (1976) intentionally drop it; flagged for someone with the reference.

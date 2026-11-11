# Seventh evaluation: report findings to bug_sweep numbers (round 0, 2026-11-10)

One row per report finding (split findings one row per part), from the "New findings (not in bug_sweep)"
sections of audit/module-evaluation-7-2026-11-10/*.md. "M196".."M274" are the new entries of
docs/bug_sweep.md §14 (branch r0-sweep); "plan item (not registered)" follows §10's rule (style or
structure only). Kind: A byte-identical, B changes results or pinned text, as the report says
(corrected where the check found otherwise; the entry says so).

| report | finding | M number (new or existing) | Kind | one-line title |
|---|---|---|---|---|
| io | IO-1 | M196 | A | The run plan's targets echo labels rows by position |
| io | IO-2 | M197 | A | Old group names refused without ERROR head, outside run_messages |
| io | IO-3 | plan item (not registered) | A | remove_file public, used only inside output_files (declarations) |
| io | IO-4 | plan item (not registered) | A | The interior GYRE digit count written three times (formula once) |
| io | IO-5 | M198 | A | legacy_controls' Sources note misdescribes 23 names as F77 |
| io | IO-6 | M199 | A | defaults/ holds two undocumented, stale rename maps |
| io | IO-7 | M200 | A | Stale text in io's tools and tests; coverage figure |
| state_setup | ST-1 | M201 | A | The star tools count an internal write as a read |
| state_setup | ST-2 | M202 | B | omega_substep missed by the per-run reset |
| state_setup | ST-3 | M203 | A | The lifetime check equates class with frequency |
| state_setup | ST-4 | plan item (not registered) | A | Star vectors without slot names (unnamed values) |
| state_setup | ST-5 | plan item (not registered) | A | Legacy names for the shared constants (naming) |
| state_setup | ST-6 | M204 | A | "Most are run copies": 7 of 31 |
| state_setup | ST-7 | M205 | A | Small stale statements (rotation_configuration's in M267; long cells: plan) |
| eos | E-1 | M206 | B | SCV downward T smoothing reads unbuilt columns |
| eos | E-2 | M207 | A | OPAL 1995 pressure check's bare-number message |
| eos | E-3 | M208 | A | Fermi-Dirac table extrapolated silently above log T 9 and past X = −3 |
| eos | E-4 | M209 | A | Second X quadratic's overflow check reads the first |
| eos | E-5 | M210 | A | "Add an OPAL vintage" step 5 does not fit a vintage |
| eos | E-6 | M211 | A | test_eos.py docstring names eos_get_gamma1 |
| eos | E-7 | plan item (not registered) | A | A value written twice and the unnamed 1e35 sentinel family (unnamed values) |
| eos | E-8 | plan item (not registered) | A | Declarations: pure, explicit shape, the test's slots |
| net | NET-1 | M212 | A | net's docs put the rate scales on the star |
| net | NET-2 | M213 | A | The light-element recipe names 2 files of a dozen |
| net | NET-3 | plan item (not registered) | A | burn_step's literal thresholds, F77 shapes, loose intents (its stale path is in M274) |
| net | NET-4 | plan item (not registered) | A | The intermediate-screening exponent written twice; Graboske's literals |
| net | NET-5 | M214 | A | nuclear_logT_cutoffs' description and "below" against "at or below" |
| numerics_math | NM-1 | M215 | A | Two checks without a mutation |
| numerics_math | NM-2 | M216 | A | The crmath bound unchecked between 23 and 99 |
| numerics_math | NM-3 | M217 | A | lir checks pin an intent(out) result |
| numerics_math | NM-4 | plan item (not registered) | A | osplin's and lir's optional ierr (declarations; cspline's with M99) |
| numerics_math | NM-5a | M218 | A | The spline's end column "linear" only for a natural end |
| numerics_math | NM-5b | M218 | A | The zero-width interval comment |
| numerics_math | NM-5c | M215 | A | mutate_numerics' docstring on a nonzero exit |
| numerics_math | NM-5d | M218 | A | Two recipes' missing files |
| numerics_math | NM-5e | plan item (not registered) | A | Style: tiny, one, double precision |
| kap_atm | KA-1 | M219 | A | Alexander 2006 cache set before a build that can fail; not reset by a re-read |
| kap_atm | KA-2 | M220 (corrects M170) | A | Five opens without iostat, two OPAL92 opens without status; reads without iostat |
| kap_atm | KA-3 | M221 | A | Allard [Fe/H]/[α/Fe] absent from the file: indices 0 and −1, format error |
| kap_atm | KA-4a | M222 | B | Conductive low-density branch: sign of the mixture term, extrapolation direction |
| kap_atm | KA-4b | M223 | A | OPAL95 low-X, low-T shortcut: unset density_shifted, stale extrapolation flags |
| kap_atm | KA-4c | M224 | A | OPAL95 fixed-Z table refuses 0.08 < Z ≤ 0.10 |
| kap_atm | KA-4d | M225 | A | Alexander 2006 fourth X node 1 − Z regardless of index_x for Z > 0.06 |
| kap_atm | KA-4e | M226 | A | Helium-core branch without else: undefined opacity, ierr 0 |
| kap_atm | KA-4f | M227 | A (doc) / B (code) | LAOL89 tolerance: ln in code, log10 in description |
| kap_atm | KA-4g | M228 | A | Second LAOL89 T grid converted over the first table's row count |
| kap_atm | KA-4h | M229 | A | Kurucz 1990 density count unbounded; malformed rows skipped silently |
| kap_atm | KA-4i | M230 | A | Alexander 2006 header check: "RUN STOPPED" and carries on; 16 of 19 log R columns |
| kap_atm | KA-5a | M231 | A | Castelli/Kurucz sentinels above ~48 kK |
| kap_atm | KA-5b | M232 | A | Kurucz atmosphere lookups: no upper range check |
| kap_atm | KA-5c | M233 | A | Allard extension four widths, headers say one |
| kap_atm | KA-5d | M234 | A | Kurucz general gravity search gap; all −999 row out of bounds |
| kap_atm | KA-5e | M235 | B | Allard per-lookup run-log notes |
| kap_atm | KA-6 | M236 | A | kap_init/atm_init install settings unchecked |
| kap_atm | KA-7 | O19 (widened) | A | Accretion surface's Allard lookup fails unless atm_choice is 4 |
| kap_atm | KA-8 | M237 | A (descriptions) / B (banner) | "Pure C and N" for the C and O table |
| kap_atm | KA-9 | plan item (not registered) | B | kap and atm write their own messages (messages layer; M112 covers the run-log writes) |
| kap_atm | KA-10 | plan item (not registered) | A | test_kap/test_atm pins without ok/FAIL checks; untested paths and coverage |
| kap_atm | KA-11 | M238 | A | "Add a table" recipes omit files |
| kap_atm | KA-12 | M239 (the unused Kurucz load); the rest plan item (not registered) | A | Kurucz table loaded for choice 4 and never read; OPAL95 scratch and atm results in module state (plan) |
| kap_atm | KA-13 | M240 | A | Stale or inexact text in kap and atm |
| kap_atm | KA-14 | plan item (not registered) | A (B for computing constants) | Code-style residue: commented-out code, history, operators, unnamed constants |
| convection_diffusion | CD-1 | M241 | B | Semiconvection's only_check_core: zone 1 whatever it is; stale bounds with no zone |
| convection_diffusion | CD-2 | M242 | B | Starspot Jacobian: unspotted ∇_rad and no d ln S/d ln T (beside O1) |
| convection_diffusion | CD-3a | M243 | B | Original settling: Noerdlinger dD1/dX sign; dD2/dX lacks the term |
| convection_diffusion | CD-3b | M244 | B | Original settling: lnΛ twice with choices 2/3 and the full Thoul solution |
| convection_diffusion | CD-3c | M245 | B | Original settling: midpoint change stored one slot off |
| convection_diffusion | CD-3d | M246 | B | Both solvers' convergence on the largest signed change |
| convection_diffusion | CD-4 | M247 | A | First-midpoint search over num_eq_points model zones |
| convection_diffusion | CD-5 | M248 | A | The twelfth zone slot unused; later zones dropped without ierr |
| convection_diffusion | CD-6 | M249 | B | Adiabatic overshoot region stale without zones; intermediate zone's |
| convection_diffusion | CD-7 | M250 | A (doc) / B (change) | settling_z_factor scales Li and Be |
| convection_diffusion | CD-8 | M251 | A | diffuse_helium always true; documented case never happens |
| convection_diffusion | CD-9 | M252 | A | light(3,json) literal |
| convection_diffusion | CD-10 | M253 | A | microdiff's iteration repeats one identical solve |
| convection_diffusion | CD-11 | M254 | A / B (the `[` format: 12 pinned logs) | Run-log formats that cannot hold their values |
| convection_diffusion | CD-12 | M255 | A | No test or pinned deck reaches semiconvection; README says it does |
| convection_diffusion | CD-13 | M256 | A | Recipes name too few files |
| convection_diffusion | CD-14 | plan item (not registered) | A / B | Formulas and constants more than once (the differing H_P and overshoot-extent copies relate to O10; 2.2's precision to M257) |
| convection_diffusion | CD-15 | plan item (not registered) | B | Messages layer (kernels print; run-log writes are M112's) |
| convection_diffusion | CD-16 | M257 | A | Stale or inexact text |
| convection_diffusion | CD-17 | M258 | A | The controls' reach mis-stated |
| rotation | ROT-1 | M259 | B | A widened transport region copies the end zone's X into the convection zone beyond it |
| rotation | ROT-2 | M260 | A | The sub-step controls are unchecked; a non-positive one hangs the step or drops the ω limit |
| rotation | ROT-3 | M261 | B | A redone sub-step burns deuterium from the failed attempt's end rate |
| rotation | ROT-4 | M262 | B | The equipotential's θ gravity drops the quadrupole term's r and reverses its sign |
| rotation | ROT-5 | M263 | B | The viscosity's Coulomb logarithm uses n_e = Z̄ρ/m_u |
| rotation | ROT-6 | plan item (not registered) | A | A star-blind kernel (Eddington–Sweet velocity) prints to standard output: the messages layer |
| rotation | ROT-7 | plan item (not registered) | A | burn_and_settle's "in place" comment and loose intents; corrected: true for the original settling routines, which still convert in place (original_settling.f90:83, :486-492; settling_grid.f90:423-427), so it waits on CD-3's fix-or-retire decision; the stale-text family is M257 |
| rotation | ROT-8 | plan item (not registered) | A | 86 lines of commented-out code in five live files; retired/ headers |
| rotation | ROT-9 | M264 | B (pinned test text) | Rotation's stop and cut messages say things that are not so |
| rotation | ROT-10 | M265 | A | rotation/README's STAR_BLIND_FILES rule is not the list's; five star-free files unguarded |
| rotation | ROT-11 | plan item (not registered) | A | substep_time_fraction on the star, wind_now twice, rotation_settings not protected |
| rotation | ROT-12 | plan item (not registered) | B | Each torque law written twice (its stale domega header sentence is in M267) |
| rotation | ROT-13 | M266 | A | The transport-mechanism recipe omits the controls' declarations and the finishing steps |
| rotation | ROT-14 | M267 | A | Stale rotation comments and docs (the README test-list bullet fixed by a8c9da6d) |
| rotation | ROT-15 | plan item (not registered) | A | Declarations and style rules across rotation/ |
| core | C-1 | M268 | A | The star calibration writes cards 51 and 52 on card 50 |
| core | C-2 | M269 | A | henyey_iterate leaves its intent(out) divergence flag undefined |
| core | C-3 | M270 | B | The starting model's warning prints an unset flag and calls itself an error |
| core | C-4 | M271 | A | The C entry cuts a long path silently; the entry and the prologue close the host's units |
| core | C-5 | plan item (not registered) | A/B | core has no messages module; terminal-only events, printing kernels |
| core | C-6 | plan item (not registered) | A | The envelope kernel, integrands and Henyey elimination have no unit test (README:89's claim with it) |
| core | C-7 | M272 | A | evolve_step's header claims cross-step state and an evolve_state that do not exist |
| core | C-8 | plan item (not registered) | A | Drivers take the star and its members; argument-less internal steps |
| core | C-9 | plan item (not registered) | A | Commented-out code in core (its stale burn_step path is in M274) |
| core | C-10 | M273 (the year's inverse, ln 10, the loop ends); rest plan item (not registered) | B/A | The inverse year, ln 10 as 2.3026, the calibration's literal loops; atm_choice literals, slots, weights tables are the plan's |
| core | C-11 | plan item (not registered) | A | core/README lacks the run flow, conventions, tests, recipes (the README standard) |
| core | C-12 | M274 | A | Stale core text (with burn_step's stale path) |

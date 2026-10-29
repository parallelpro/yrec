# net/ review (2026-10-29, read-only evaluation; paths under src/)

Scores vs rotation/: Readability 4, Organization 4, Maintainability 3, Extensibility 3, Docs accuracy 3.

In good shape: the table-driven network (rows carry rate law, fit, screening, Q, chain, change list;
burn_implicitly and the energy sum read them), the self-checking table, the closed-form trace integrator,
pure kernels, ierr not stop, a small facade (35 names, all used outside). Main obstacle: the index layout
around the table -- reaction slots 1-13 share cross_section_scale with literal slots 14-17 (pep, Be7, hep);
derivative-scale arrays sized 8 with slot 8 aliased to Be7+p; slots 9 and 13 "unused" but special-cased;
chain energies leave net/ as five named fields core re-indexes by hand. Second: docs and comments drifted
behind the 2026-10 restructures.

## Possible bugs
1. Output eps omits the alpha captures -- core/structure/structure_equations.f90:241-243 (feeds profile
   eps_nuc, io/profile_output.f90:77, and GYRE eps, stitched_model.f90:597). VERIFIED by the integrator:
   eps_total = energy_gen_component(1..5) (pp, ppII, ppIII, CNO, 3a); C12(a,g), C13(a,n), N14(a,g)
   (eps_alpha_capture) go only into luminosity_breakdown (structure_equations.f90:184-185). The pulse
   column eps_eps_T (structure_equations.f90:296) is the derivative of the full total, so the pair is
   inconsistent in He burning. Inherited (coefft.f SESUM = EG(1..5)). H-burning output stays byte-identical;
   He-burning output moves (B). Fix together with structural rec 2.
2. nuclear_energy.f90:156-161 with :207-208 -- the capped deuterium rate keeps its uncapped derivatives:
   in the surface CZ the cap sets eps_D = (N_A/m_d/tau) X_D Q (independent of rho, T), yet
   sum_energy_generation adds eps_D*1 to d eps/d ln rho and eps_D*(fit's dln/dlnT) to d eps/d ln T. Only the
   Henyey Jacobian is wrong (slower Newton in pre-MS D burning). Zero derivatives when the cap binds (B).
3. nuclear_energy.f90:193-197 -- dead code and a trap: zeroes the derivatives of r_ao16_unused /
   r_c12c12_unused, which no row fills; a newcomer who puts a reaction in slot 9 or 13 silently loses its
   Jacobian.
4. reaction_rates.f90:151-153 and README:157 -- logT_cutoffs(5) does NOT gate the solar neutrino fluxes;
   nuclear_energy.f90:103-104 gates them only on calc_neutrinos (as F77 engeb). Doc error.
5. reaction_rates.f90:158-161 with star_info cross_section_scale(17) -- slots 8-13 of the scale vector are
   read in and never used (filled from sstandard when LNEWNUC=.false., io/read_controls.f90:478); every
   Caughlan & Fowler row has scale_slot = 0. sstandard(10) to scale C12(a,g) silently does nothing.
6. reaction_rates.f90:525 -- frac_n15_pg computed and never read.
7. test_net.f90:226-230 -- sets star%xa(i_h2,1) (unused), calls deuterium_rate(..., 1.ge.star%jcz, ...): the
   printed value depends on star_info's default init; comment and label stale.

## Findings (ranked)
1. reaction_rates.f90:37-55, :158-162, state/star_info_lib.f90:248, core/run/monte_carlo.f90:90-93 -- the
   scale layout blocks the README's add-reaction recipe: raising num_reactions to 14 aliases xs_pep=14; a new
   S-factor row has no derivative-scale slot (qs0e_scale sized iq_be7p=8; slot 8 is Be7+p); sizes 17/8/8 are
   literals in star_info; monte_carlo writes slots 1,2,3,16 by number. Fix: per-slot arrays sized
   num_reactions (rate, S', S'' scales), a separate named "extra scales" vector for pep, Be7e, Be7p, hep,
   Be7p's derivative scales in their own fields; delete slots 9 and 13 (dropping summed zeros is exact).
   A. M-L (state/, setup/user_inputs, io/read_controls' legacy sstandard map, monte_carlo).
2. README.md stale/wrong: :98, :129 (Henyey calls compute_nuclear_energy from
   core/structure/shell_equations.f90:266, not via compute_zone_energy); :105 (burn_li_be_midpoint is
   rotation/retired/*.f90.txt, not compiled; no -D flag); :192 (the energy-source list exists:
   core/structure/energy_sources.f90, examples/other_energy); :157 (logT_cutoffs(5)); :42 (says
   ctrl%calc_neutrinos; net reads net_settings); :54 (table check also tests reactants subset of solved
   species); :183-188 recipes omit files (see trace). Rewrite to rotation's standard (errors/messages
   convention, what each test guards, coverage figure, recipes naming every file). A, S.
3. Two formulas in two copies: nuclear_energy.f90:148-161 vs reaction_rates.f90:200-203 and
   light_elements.f90:173-178 -- the d(p,g) fit coefficients (3.72, 0.112, 3.38, 2.65) retyped for the
   T-derivative; the overturn cap written twice with its 1e-6 and 1 s thresholds. Fix:
   deuterium_proton_capture_rate also returns dln/dlnT; one overturn_capped_rate function. A, S.
4. test_net.f90:421-428, :455-458, :185 -- coverage gaps: the BPS thermal-neutrino branch
   (nuclear_energy.f90:332-376; the default; 340 of 520 decks) never run or derivative-checked because
   install_production_settings forces Itoh; its comment (:449-451) wrongly says defaults have Itoh and
   fluxes on (physics.defaults:258,324 has both off). No energy point carries deuterium
   (compute_deuterium_burning and its cap untested). Untested: strong screening, burn_implicitly's
   non-convergence ierr, compute_nuclear_energy below the cutoff. :185 prints scale slots 1-16 of 17.
   Add an Itoh-off point at log T 8.1, a pre-MS surface-zone point with D and the cap binding (FD-checked);
   mutation-confirm each check. A, S-M.
5. Unnamed numerical thresholds: network.f90:286 (1e-5), :552 (1e-8); nuclear_energy.f90:122 (1e-30), :209
   (1e-12), :159/:260 and light_elements.f90:176 (1e-6), :156/light_elements:174 (1 s), :286 (literal 10 =
   n_nu_fluxes); light_elements.f90:101/:134 (700). Physical prefactors: reaction_rates.f90:298/299/316;
   :498 (8.584981d22 = N_A/m_Be7, bare mass); :579 etc. (1.25388d22 = N_A/(4.002603*12), neutral masses,
   unlike bare masses elsewhere). Named parameters with unit and source. A, S.
6. reaction_rates.f90:144-176 -- four modules in one: the run's settings (incl. Li/Be scales and neutrino
   switches), the index layouts, the screening plasma (:238-390), the rate laws. light_elements imports
   reaction_rates only for settings. Split into net_configuration (settings and layouts) and a screening
   module. A, M.
7. reaction_rates.f90:213-236, :480, :511, :637 -- nuclear_rates filled piecemeal: compute_zone_conditions
   does two jobs; three branch/hep routines take nr intent(inout). Store zone as a component of nr; make the
   branchings functions. A, S-M.
8. nuclear_energy.f90:183-217 -- sum_energy_generation also zeroes reactions_per_gram (silencing the
   neutrino fluxes when eps <= 1e-12) and the dead slots. Move the flux gate into
   compute_solar_neutrino_fluxes. A, S.
9. Misleading names: c12_alpha_scale/o16_gamma_scale (reaction_rates.f90:163) scale N15(p,a)C12 and
   N15(p,g)O16 -> n15_pa_scale/n15_pg_scale; net_lib.f90:30 calls iq_be7p a "Q-value" (it is the S'/S'' slot;
   net also has q_mev); qs0e/qqs0ee -> s_prime_scale/s_double_prime_scale; rate_row vs "rate slot" in
   docs; eg at nuclear_energy.f90:256. A, S (crosses star_info and configure_network).
10. Stale comments: net_lib.f90:11 ("equilibrium burns"); network.f90:6 (mixing/), :26-30 (b_j omits the Be7
   branch), :316-319 (check omits the energy test); light_elements.f90:4-12 ("dburn, sub-stepped");
   thermal_neutrinos.f90:230-233; test_net.f90:97, :150-153 (core/read_input.f90, parmin -- neither exists),
   :210; nuclear_energy.f90:387 (rates_times_reactants), :468 (star%neutrino_flux); reaction_rates.f90:887-888
   (first_zeroed_rxn). A, S.
11. reaction_rates.f90:680-1064 (36% of the file) and nuclear_energy.f90:381-517 -- the Sources text: 55
   lines duplicated between files; contradictions (hep scale x1.913 Wolfs at :782 vs x0.1625 Carlson at
   nuclear_energy:455; S0(N15 p,a) "675000" at :958 while the code uses 67500; J=8 differs); describes
   deleted code. Keep the glossary (:680-692), one verbatim block marked historical (or move to docs/),
   drop duplicates; ALL-CAPS prose out. A, S.
12. reaction_rates.f90:155-170 -- net_settings_t defaults differ from the control defaults
   (weak_screening_threshold 0 vs 0.03; cutoffs 0 vs 6.5...): a caller skipping configure_network gets
   intermediate screening everywhere and burning at any T. Default to the controls' values, or a
   configured flag. A, S.
13. thermal_neutrinos.f90:1528-1541, :1590-1610 -- the only thing keeping compute_nuclear_energy impure is
   four data statements (implicitly SAVEd) in ifermi12/zfermim12; make them parameter arrays and declare
   the chain pure. Merge the two wrappers (:31-105; ymass discarded). A, S.
14. network.f90:320-330 -- the table self-check misses invariants: chain in 1..n_chains (else an OOB write at
   :472); rate_row in range; every shares_rate row has an owning row; rate_law known (reaction_rates.f90:427
   silently gives zero otherwise); an S-factor row's scale_slot fits the scale arrays; each slot one owner
   (reaction_counts :308 overwrites). Also run the check in configure_network (with ierr), not only in
   test_net. A, S.
15. Small: nuclear_energy.f90:97-98 computes network_mass_fractions(xa) twice; network.f90:451 redundant
   use; light_elements.f90:205 unused use math_lib; s0_n15* public needlessly (reaction_rates.f90:65); ymass
   discarded; reaction_rates.f90:674-677 blank lines. Outside: core/run/run_driver.f90:48 bare use net_lib
   (uses nothing); burn_step.f90:6-8 names compute_reaction_rates (not exported) and an "equilibrium burn".

## Structural recommendations
1. One reaction layout (net_configuration): reaction slots, per-slot scale arrays sized num_reactions, a
   named extra-scale list; rv_*, extra slots and star_info's scale sizes derive from these; monte_carlo and
   user_inputs address scales by name; slots 9 and 13 go. "Add an S-factor reaction with input S0, S', S''"
   becomes one row, one slot and one control mapping. M-L. A.
2. Chain energies as an array: nuclear_energy_result%eps_by_chain(n_chains) + chain_name(:) instead of five
   named chain fields and eps_alpha_capture; shell_equations, structure_equations (1..5/1..6 loops, i_lum_*)
   and outputs loop over chains; eps_nuc = sum over all chains (fixes bug 1). M. B for He-burning
   eps_nuc/GYRE eps; A elsewhere. KEEP star%luminosity_breakdown's layout and i_lum_* names unchanged.
3. One species table: a parameter array of species_t (name, xa slot, A, mass excess, solved/explicit,
   screening isotope) from which s_*, n_solved/n_network, network_species_slot, mass_number,
   mass_excess_mev, screening_mass_fractions derive; solved-before-explicit ordering checked. S-M. A.

## Newcomer traces
Unscaled C&F-style reaction among tracked species: reaction_rates.f90 (r_* slot raising num_reactions, a
rate-law constant, a rate routine with exact derivatives, a case in compute_screened_rate), network.f90 (one
row), regenerate expected_test_net.out -- two files, good. Stuck points: slots 9/13 trap; the prefactor
convention N_A/(A1 A2) with neutral masses, abundance factors omitted, 1/3! for 3a; the four Graboske
screening factors with no helper from Z1, Z2; Q net of neutrinos (checked only by q_from_masses).
S-factor reaction with input S0/S'/S'': ~10 files (net_settings_t sizes, iq_be7p/xs_* renumbering,
star_info_lib.f90:248, setup/user_inputs.f90 three set_* routines and _ref/_a98 data, registry, defaults,
sync, monte_carlo.f90:90-93, read_controls' sstandard(17), test_net:185). New chain: chain_*/n_chains,
nuclear_energy_result field, store_chain_energies, shell_equations.f90:106/270-275,
structure_equations.f90:179-185/241-252, profile/history columns. Neutrinos: neutrino_flux_layout,
compute_solar_neutrino_fluxes.
Add an isotope (Ne20 solved): composition_layout (+ the literal 15 in 39 declarations across 17 files --
known Phase 5 debt); network.f90 s_* (solved before s_o18 renumbers and changes pivot order, B),
n_solved/n_network, network_species_slot, mass_number, mass_excess_mev; screening_mass_fractions; thermal
neutrino composition; io/model_files.f90:45-74 slots(15) and header (format change);
core/model/starting_model.f90; element_diffusion/species_table.f90; rotation num_tracked_species; outputs;
tests. ~12-15 files; the README lists 3 steps.

Fine as is: burn_implicitly, compute_abundance_rates, reaction_rate(_derivative), branch_fraction,
q_from_masses/baryon_imbalance/write_network_table; burn_exposure, surviving_fraction, accreted_survival,
burn_deuterium_region, compute_li_be_rates; the rate-law routines (reaction_rates.f90:445-632); Timmes'
sneut/azbar bodies (kept as authored by decision); the reaction table's design.

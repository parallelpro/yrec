# burn -- readability wave 1 (R1-R3)

Worktree `/Applications/YREC-wt/burn`, branch `rs/burn`, files `src/core/burn_lib.f90`, `src/net/net_lib.f90`. `src/deps.mk` unchanged (regenerated after every batch; no use/module change affected it).

| batch | commit | verification |
|---|---|---|
| R1 comments / dead code | `09f9d0d` | gate1 IDENTICAL; pins 37 passed (`logs/burn.pins.R1.log`, PINS EXIT 0) |
| R2 named indices / constants | `01ef296` | gate1 IDENTICAL; pins 37 passed (`logs/burn.pins.R2.log`, PINS EXIT 0) |
| R3 shared helpers | `b47a358` | gate1 IDENTICAL; pins 37 passed (`logs/burn.pins.R3.log`, PINS EXIT 0); aux 22 passed (`logs/burn.aux.log`, AUX EXIT 0) |

Build warnings in the two files: none new (the only one, `light_element_save(3,json)` moved to static storage in liburn, `-Wsurprising`, predates the sweep).

## R1 -- comments and dead code (09f9d0d)

- burn_lib.f90: headers/comments no longer reference dburn.f90/liburn.f90/liburn2.f90, `common/deuter/`, `common/masschg/`, `common/fluxes/`, `common/nuloss/`, read_input.f90 or "COMMON STATEMENT"; dburn header lists the four real dburn/dburnm differences (rate source, timestep units, skip threshold, accretion weighting) [item 9, comment half]; deutrate/engeb/lirate88 headers name the `star%` members they write [item 7]; PEP/PET comment states they are absolute d eps/d ln rho, d eps/d ln T [item 5, comment half]; DRATT/DRATRO comments drop the wrong log-base claim; the two neutrino-loss branches are labelled Itoh 1996 vs Beaudet-Petrosian-Salpeter 1967 [item 25, label option]; sentinel comments at `cz_base_radius_prev = 0` and `log_rate_li6_prev <= 0` [item 20]; comment on why the 1e-5 rate-floor sweep runs after the energy sums.
- burn_lib.f90 dead code: eqburn's four `else dx_dt/dy_dt/dc_dt/do_dt = 0` branches (all four outputs zeroed unconditionally at the top) and the CN block re-indented under the H-burning `if` [item 11]; engeb local `en` (assigned twice, never read) and its `use luout_lib`; deutrate's empty `if(i.eq.star%jcz)` block; three commented-out EPP lines in compute_neutrino_emission; liburn/liburn2 `li6/li7/be9_substep_depletion(json)` arrays -> scalars [item 19]; `if (refine_idx > 11)` -> `if (.not. converged)`; bare `continue`/`return` before `end`, dead `! 67 CONTINUE` labels [item 22]; unused imports (`i_grad_ad`, `i_grad_rad`, `luout_lib` in liburn2, `phys_const_lib` in eqburn/dburn); mis-indented `end do`.
- net_lib.f90: stale "Completed the merge" header paragraph and "this module's own liburn" removed [item 6]; nulosses/azbar/sneut intents, commented-out declaration remnants removed, the tfac2 bugsweep note moved out of the continued `parameter` statement [item 24, part]; rates: never-read `o16_gamma_frac`, the `rate(i).lt.1.d-30` guard inside the i=1..7 loop (subsumed by the later `.le.1.d-5` sweep over all 13; rate(1..7) not read in between), "FLUX COMMON BLOCK" comment -> names `star%ctrl%weak_screening_threshold`, the stale "unsuffixed literals" header sentence [item 8]; safedivexp `implicit none` + intents + one-line description; ifermi12/zfermim12 commented-out KC 2025 lines.

## R2 -- named indices and constants (01ef296)

- net_lib.f90 module header: `integer, parameter, public` blocks `r_pp=1 .. r_c12c12_unused=13, num_reactions=13`; `s_pep=14, s_be7e=15, s_be7p=16, s_hep=17` (cross_section_scale slots); `iq_be7p=8` (slot 8 of q1..q5 / qs0e_scale / qqs0ee_scale); `iso_n=1 .. iso_mg24=13, num_isotopes=13`; `nu_ionmax=4`; `double precision, parameter :: gyr_amu_per_sec_gram = 5.240358d-8` (the C21 of engeb/rates/deutrate); `integer, parameter :: s0_n15pg_kevb = 64, s0_n15pa_kevb = 67500` (kept integer = original literal kind) [items 2, 14, 15, 24].
- net_lib.f90 rates: local `num_isotopes/num_reactions/c21` DATA removed; arrays dimensioned by the named counts (`q6..q8(r_po16)`); `nz` -> `first_zeroed_rxn` with a comment [item 12, rename half]; all literal subscripts of `rate`, `screening_factor`, `mass_frac`, `q1..q5`, `qs0e_scale`, `qqs0ee_scale`, `cross_section_scale` named; `rate(k)*c21` -> `rate(r_*)*gyr_amu_per_sec_gram`; commented-out `O16GAMMA*64.`/`C12ALPHA*67500.` lines dropped. neutrino/nulosses: one `nu_ionmax`.
- burn_lib.f90 engeb: same names at 99 reaction, 15 isotope, 7 Be7+p-slot, 4 cross-section subscripts; private `nrxns/num_isotopes/years_per_sec_over_amu` DATA removed; work arrays dimensioned by `num_reactions/num_isotopes`, `eg(50)` -> `eg(num_reactions)` [item 17, size half]; `f1..f4` -> `frac_be7_ecap/frac_be7_pcap/frac_n15_pa/frac_n15_pg` [item 13, engeb half]; dummies `dlnepsilon_dlnrho/dlnepsilon_dlnt` -> `deps_dlnrho/deps_dlnt` [item 5]; `nz` -> `first_zeroed_rxn`.
- burn_lib.f90 deutrate: `dl,tl,x,i,itlvl` -> `log_density, log_temperature, hydrogen_fraction, zone_idx, iteration_level`; local `c21` -> `gyr_amu_per_sec_gram` [item 18].
- burn_lib.f90 eqburn: `zone_avg_abundance(1,2,3,4,5,6,7,9)` -> `i_h1, i_he4, i_metals, i_he3, i_c12, i_c13, i_n14, i_o16` (star_info_lib parameters, already public).
- burn_lib.f90 dburn/dburnm/liburn/liburn2/lirate88: `composition(12..15,:)` and `star%ctrl%accreted_composition(12..15)` -> `i_h2/i_li6/i_li7/i_be9` (111 sites); `only:` lists extended. Declarations stay `composition(15,json)`.
- All callers of the renamed routines pass arguments positionally (checked by grep of `src/`), so the dummy renames are interface-safe.

## R3 -- shared helpers (b47a358)

- burn_lib.f90: `locate_cz_base(radius, env_cz_zone, env_cz_zone_old, num_zones, cz_base_zone, cz_base_zone_old, cz_base_frac)` and `average_cz_abundances(composition, mass_coordinate, shell_mass, cz_base_zone, cz_base_zone_old, cz_base_frac, num_zones, li6/li7/be9_cz_start, log_rate_*_cz_start, li6/li7/be9_cz_end, cz_mass_end, log_rate_*_cz_end)`, private module subroutines ahead of liburn, replacing the 63-line and 76-line blocks that liburn and liburn2 carried byte-for-byte (diffed before extraction; the "SKIP IF WHOLE CZ IS BELOW THE BURNING THRESHOLD" early return stays in the callers). Same `use` set as the callers (math_lib overrides exp/log). Eight block-only locals removed from each caller. Net -57 lines [item 3, roadmap 3].

## Deferred (cross-domain)

1. **`star%reaction_rate_1..13(json)` -> `star%reaction_rate(13,json)`** [item 1]. `src/state/star_info_lib.f90`: replace the thirteen `reaction_rate_N(json)` members by one `reaction_rate(num_reactions,json)` (num_reactions from net_lib, or a copy in star_info_lib). Then burn_lib.f90 compute_neutrino_emission (the 13 `star%reaction_rate_k(shell_index) = reaction_rate(r_*)*gyr_amu_per_sec_gram` lines) -> `star%reaction_rate(1:num_reactions,shell_index) = reaction_rate(1:num_reactions)*gyr_amu_per_sec_gram` (same product per element, same order), and every reader: `grep -n "reaction_rate_[0-9]" src/` -- mixing/mix.f90, util/timestep_limit_hburn.f90, core/henyey*, io/*. eqburn's nine `rate_*` dummies and rates' thirteen `rate_*(json)` outputs could then take array sections.
2. **Names for `neutrino_flux(9)` and `(10)`**. `src/state/star_info_lib.f90` lines ~60-72: the comment says slots 9-10 are spares, but burn_lib.f90 compute_neutrino_emission stores fictional He3+He3 and He3+He4 "fluxes" there (used by the neutrino-flux table). Add `i_nu_he3he3 = 9, i_nu_he3he4 = 10` next to `i_nu_f17 = 8`, fix the comment, and replace the two literal subscripts in burn_lib.
3. **Move the reaction/scale/isotope index block to the owning state module** (structural rec. 1). If a `net_def`-style home is wanted, the `r_*`, `s_*`, `iq_be7p`, `iso_*`, `num_reactions`, `num_isotopes` parameters now at the top of net_lib.f90 can move verbatim; sstandard readers in `src/io/read_controls.f90` / `map_user_inputs.f90` (`cross_section_scale(14..17)`) should then use `s_pep..s_hep`.
4. **`gyr_amu_per_sec_gram` into `phys_const_lib`** [item 14 as proposed]: `src/state/phys_const_lib.f90` add `double precision, parameter :: gyr_amu_per_sec_gram = 5.240358d-8` and delete the net_lib copy (engeb, rates, deutrate all `use phys_const_lib` already; deutrate would need nothing else). Left in net_lib because phys_const_lib is outside the assignment.
5. **`star%cross_section_scale(14..17)`, `qs0e_scale(8)`, `qqs0ee_scale(8)` in `src/io`** (SSTANDARD namelist mapping): use `s_pep`, `s_be7e`, `s_be7p`, `s_hep`, `iq_be7p` from net_lib where the namelist values are copied.

## Reverted (changed numbers)

None. Every batch was byte-identical on the 37-case pin selection at the first attempt.

## Skipped (disagree with reviewer / too risky)

- item 4 (split compute_neutrino_emission into three): moves 250 lines inside engeb's internal-procedure family; no shared code, no reader benefit that the three section comments do not already give; risk of touching the bit-drift zone.
- item 10 (merge eqburn's adjacent `zone_begin.ne.zone_end` / `.eq.` if/else pairs, burn_lib 90 and 117): possible as pure block movement, but each pair handles one concern under its own comment (mass-averaged abundances, then mass-averaged rates); the merged form interleaves the two 30-line bodies. Judged restyling of clear code; left as is.
- item 12 (logical `has_hydrogen` + `rate(1:7)=0` rewrite): kept the `do i=first_zeroed_rxn,num_reactions` idiom (renamed, commented); the rewrite reorders stores and adds an array-section assignment for a two-line readability gain.
- item 16 (screening 8-liner helper): both engeb copies are inside the measured bit-drift zone; the `rates` copy alone would be a one-caller helper.
- item 17 (split `q1..q5(8)` slot 8 into `q1_be7p` scalars): changes the DATA tables that engeb and rates share; slot named `iq_be7p` instead.
- item 21 (make liburn/liburn2 safedivexp use uniform): class B by the reviewer's own note.
- item 23 (delete engeb's 150-line header, point at rates): the header carries the BP00/S0 provenance next to the DATA it documents; removing one copy makes the engeb tables undocumented in place. Left for the engeb/rates de-duplication (R5/R6).
- item 25 (move the Beaudet branch to its own internal subroutine): labelled instead (R1); moving it is an extraction inside engeb.
- SUMMARY 1.1 #1, #3, #4, #9 and all of R6: out of scope by the assignment.
- alfmlt/phmlt/cmxmlt, burn_settle_mix / burn_mix_extrapolated (BS extrapolation): not touched.

## Anything the audit missed

- `star_info_lib` already exports `i_h1..i_be9` for the 15 composition slots; burn_lib used bare 12..15 in 111 places (now named) -- worth the same sweep in mixing/ and rotation/ where `composition(13..15,:)` recurs.
- `rate(9)`/`rate(13)` are set to zero twice in `rates` (inside the `do i=first_zeroed_rxn,num_reactions` loop and again unconditionally) -- harmless, left as is (the second store is the documented "XEROING OUT OF REACTIONS 9 AND 13").
- `locate_cz_base` reads the loop index after a `do ... exit` search (`if (zone_idx .lt. 1)`); legal but easy to misread -- a found-flag would be clearer (not changed: same-numbers rule is safe here, but it is restyling).
- engeb's `q6/q7/q8(7)` vs `q1..q5(8)` mismatch [item 17] is now visible from the declarations (`q6(r_po16)` vs `q1(iq_be7p)`).

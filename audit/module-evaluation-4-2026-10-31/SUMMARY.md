# Module evaluation 4, 2026-10-31 -- summary

The same five read-only reviews as the three before (audit/module-evaluation-2026-10-29,
module-reevaluation-2026-10-29, module-evaluation-3-2026-10-30), at yrec-next e8c847ba,
scored this time out of 10 (10 = at or above rotation/'s bar with nothing material left;
the earlier 1-5 scores are doubled for comparison). Since the third evaluation the code
changed only by the controls split (s%control/s%physics replace s%ctrl/s%job) and docs.
The open defects of docs/bug_sweep.md section 9 (M89-M126) weighed in the scores but were
not re-reported; each report lists new findings and what would raise each score.

| Module | First | Second | Third | Now (/10) | Now: readability / organization / maintainability / extensibility / docs |
|---|---|---|---|---|---|
| math | 6.4 | 8.8 | 9.6 | 9.3 | 10 / 10 / 8.5 / 10 / 8 |
| numerics | 6.2 | 7.6 | 8.8 | 8.5 | 9 / 9 / 7.5 / 9 / 8 |
| net | 6.8 | 8.0 | 8.6 | 8.2 | 9 / 8.5 / 8.5 / 7.5 / 7.5 |
| setup | 6.6 | 7.4 | 8.2 | 8.2 | 9 / 8 / 8 / 8 / 8 |
| state | 4.8 | 6.8 | 7.6 | 7.8 | 8 / 7.5 / 8 / 8 / 7.5 |
| io | 4.8 | 6.6 | 7.6 | 7.5 | 7.5 / 7.5 / 7.5 / 7.5 / 7.5 |
| eos | 5.2 | 6.6 | 7.6 | 7.3 | 7.5 / 8 / 7 / 6.5 / 7.5 |

Where the code did not change (eos, net, numerics, math) the small drops come from new
findings and a stricter reading, not from regressions; state rose with the split, io's
gains from it (organization, maintainability) were offset by new findings.

New findings the reviews flag first (not yet in bug_sweep; checked by hand where noted).
Each one recorded is now an open item of docs/bug_sweep.md section 10, M127-M165; its
number is given in brackets:
- net: the Li6 rate scale is taken against a reference S-factor of 5.5 MeV b
  (io/legacy_controls.f90:226, li6_reference_xs = 5.5d3 keV b, checked), but the
  Caughlan & Fowler 1988 fit net computes has S(0) of about 3 MeV b (reviewer 2.69;
  a rough check gives about 3.1), so the 1,036 shipped decks that set
  li6_rate_value = 3100 burn Li6 at roughly half the rate they ask for. Inherited (F77).
  [M127; recomputed: the fit's S(0) is 3144 keV b (the reviewer's 2689 carries an extra
  factor of the reduced mass), and the 1,175 decks that set the override and follow the
  light elements burn Li6 at 0.57 of the rate they ask for, Be9 at 1.095 (0.94 by
  (p,g)); Li7's reference is its fit's.]
- net: burn_implicitly can return ierr = 0 with all-NaN abundances when the NaN is only
  in the Newton right-hand side (seen at -O0 and -O2 in a scratch build). Latent. [M131]
- numerics: the envelope integration drops bsstep's derivatives_ierr, so an EOS or
  opacity failure inside a midpoint substep ends the run with exit 0 instead of 1
  [M132]. splint on a decreasing table [M134], osplin on a decreasing table, simeqc with
  no right-hand side and polint at n = 1 [M135] return plausible wrong values with
  ierr 0. Latent.
- eos: README's source recipe has a new source compute its mean weights from a fresh
  record (R/mu = NaN) [M148]; compute_mean_weights ignores use_diffusion_z, so OPAL/MHD
  and Yale/SCV use two mixture rules (a 3.9e-5 step in 1/mu_ion at the OPAL edge; F77)
  [M129; the 3.9e-5 is spread across OPAL's edge ramp, and 186 inlists set
  use_diffusion_z]; SCV's H+ fraction counts H2 once in its denominator (up to 33% high
  where H2 and H+ coexist at log P <= 10; F77) [M128].
- io: following the "add a table path" recipe literally gives code that does not
  compile, and no generator catches it [M149]; the registry describes seven controls by
  what their run copies do [M152; eight rows]; FGONG writes 0 for the mixing length
  [M130, with globals 4 and 5].
- state/setup: atm_choice in three homes [M140]; one run copy that is really a per-step
  value [left out: the run does change it, as its comment says]; the Monte-Carlo
  record's misleading names [M151]; the EOS and kap halves of the surface mixture
  installed directly from five core call sites [left out: structure; setup/README's
  rule names the installer calls it covers, and M85, M86 and M96 hold the recompute].

The other findings recorded, by report: net M141 (check_boundaries misses pure
procedures), M146 (the chain recipe and L_H), M147 (the isotope recipe), M139 (neutrino
fluxes without calc_neutrinos), M161-M162 (stale pointers); eos M136-M137 (readers, the
Saha cutoff), M144 (the defaults check), M156-M157 (which source fills which field;
counts, comments, the MHD header); io M138 (unchecked writes), M153-M155 (the records,
the registry's columns, the disk fields), M165 (stale comments); numerics and math M133
(bsstep's NaN), M142 (the `**` checker), M145 (pow_i's zero base), M158-M160 (status text,
hunt's header, HOW_YREC_WORKS); state and setup M143 (contract coverage), M150 (the
run-copy recipe), M163-M164 (stale statements). Left out besides the two above: net's
network_checks gaps (nothing claims those checks), README:157's pins sentence (it holds)
and net_lib's dead exports; numerics/README's sentence on the searches (not false) and
mutate_numerics.py's dead strings (every mutation applies); io's interpret_kind_cards
argument and imports through legacy_controls (structure); the builders' purity and
max_runs (structure); eos's history comments (style). bug_sweep section 10's
introduction gives the reasons.

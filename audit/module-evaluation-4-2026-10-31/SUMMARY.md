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

New findings the reviews flag first (not yet in bug_sweep; checked by hand where noted):
- net: the Li6 rate scale is taken against a reference S-factor of 5.5 MeV b
  (io/legacy_controls.f90:226, li6_reference_xs = 5.5d3 keV b, checked), but the
  Caughlan & Fowler 1988 fit net computes has S(0) of about 3 MeV b (reviewer 2.69;
  a rough check gives about 3.1), so the 1,036 shipped decks that set
  li6_rate_value = 3100 burn Li6 at roughly half the rate they ask for. Inherited (F77).
- net: burn_implicitly can return ierr = 0 with all-NaN abundances when the NaN is only
  in the Newton right-hand side (seen at -O0 and -O2 in a scratch build). Latent.
- numerics: the envelope integration drops bsstep's derivatives_ierr, so an EOS or
  opacity failure inside a midpoint substep ends the run with exit 0 instead of 1.
  splint on a decreasing table, osplin on a decreasing table, simeqc with no right-hand
  side and polint at n = 1 return plausible wrong values with ierr 0. Latent.
- eos: README's source recipe has a new source compute its mean weights from a fresh
  record (R/mu = NaN); compute_mean_weights ignores use_diffusion_z, so OPAL/MHD and
  Yale/SCV use two mixture rules (a 3.9e-5 step in 1/mu_ion at the OPAL edge; F77);
  SCV's H+ fraction counts H2 once in its denominator (up to 33% high where H2 and H+
  coexist at log P <= 10; F77).
- io: following the "add a table path" recipe literally gives code that does not
  compile, and no generator catches it; the registry describes seven controls by what
  their run copies do; FGONG writes 0 for the mixing length.
- state/setup: atm_choice in three homes; one run copy that is really a per-step value;
  the Monte-Carlo record's misleading names; the EOS and kap halves of the surface
  mixture installed directly from five core call sites.

# Module re-evaluation, 2026-10-29 -- summary

The same five read-only reviews as audit/module-evaluation-2026-10-29/,
repeated at yrec-next 7aa138b2 after the module series (docs/ROADMAP.md,
"Module evaluation follow-up"; docs/bug_sweep.md section 7), on the same
1-5 scale relative to rotation/. Each report marks the first evaluation's
findings fixed, partly fixed or open, lists what the series introduced, and
says what would raise each score by one point.

| Module | Before (avg) | Now (avg) | Now: readability / organization / maintainability / extensibility / docs |
|---|---|---|---|
| math | 3.2 | 4.4 | 4.5 / 4.5 / 4 / 4.5 / 4.5 |
| net | 3.4 | 4.0 | 4.5 / 4.5 / 4 / 3.5 / 3.5 |
| numerics | 3.1 | 3.8 | 4 / 4 / 3.5 / 4 / 3.5 |
| setup | 3.3 | 3.7 | 4 / 3.5 / 3.5 / 3.5 / 4 |
| state | 2.4 | 3.4 | 3.5 / 3 / 3.5 / 3.5 / 3.5 |
| eos | 2.6 | 3.3 | 3.5 / 3.5 / 3 / 3 / 3.5 |
| io | 2.4 | 3.3 | 3.5 / 3 / 3 / 3.5 / 3.5 |

Checked against the source by hand when this summary was written:
- io/run_cards.f90:39 forces lfirst(1) after it has been copied into
  first_call_flag, so card 1 of ten shipped decks evolves an empty model
  (a one-line note in the ROADMAP since August; never fixed).
- numerics/quadratic_spline.f90: when the last point in range lies on node
  n-1 and points above the table follow, the points above use arcs that were
  never computed (introduced by M9's rewrite; latent).
- eos/README.md's "add an EOS source" recipe offers a step after OPAL that
  blend_eos's early returns skip; it also names the renamed read_controls.f90.
- state/star_info_lib.f90's Sources notes were damaged by the member-rename
  script ("zahn_coupling (zahn_coupling)", "controls.!").
- setup/user_inputs.f90's s0_p_o16_ref is the single-precision 9.4E0.
- state/physical_constants.f90's cgl has no initializer (0 until run_setup).

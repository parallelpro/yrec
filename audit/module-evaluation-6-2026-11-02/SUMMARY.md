# Module evaluation 6, 2026-11-02 -- io and state/setup after the cleanup round

Two read-only reviews at yrec-next 2632f06c, out of 10 against rotation/'s bar
(10 = at or above it with nothing material left; 9 = only minor items left).

| Module | Fourth | Fifth | Now | Now: readability / organization / maintainability / extensibility / docs |
|---|---|---|---|---|
| io | 7.5 | 8.3 | 8.9 | 9 / 9 / 9 / 9 / 8.5 |
| state | 7.8 | 8.3 | 8.9 | 9 / 9 / 9 / 9 / 8.5 |
| setup | 8.2 | 8.8 | 9.2 | 9.5 / 9 / 9 / 9.5 / 9 |

The cleanup round did what the fifth reports asked (io: 37 of 40 items done, 2 partly,
1 open -- CI's first GitHub run; state: all but item 9 and the Phase 5 type split).
Every criterion but docs is now at 9 for io and state; docs accuracy (8.5) holds the
averages just below. io has no open bug of its own in docs/bug_sweep.md.

New findings (not in bug_sweep when this summary was written):
- io: diffuse_lithium silently turns on metal settling, the Thoul routines without
  the fit and Coulomb log 4 (derived_controls.f90:294-307, from the F77), which its
  description does not say -- 169 shipped inlists set use_diffusion_z = .FALSE. with
  diffuse_lithium = .TRUE. and run with metal settling; descriptions cite other
  controls by code name in 5 places and inlist name in 10; gen_controls_state's
  refusal (the first make runs) prints no cells guide; open_run_log cuts its path at
  the first ".short" anywhere; GSM files are still written in place (M177's class);
  interior and profile GYRE names can collide from 10,000 profiles; a removed deck's
  stale hashes.txt row passes; test_io's run-log checks log deck_lines before filling
  them (pass only under -finit-local-zero); stale counts and docstrings.
- state/setup: disk_gate_age_gyr (legacy SAGE) is written only after a converged
  step (step_driver.f90:557), so a later Monte-Carlo run or a card that re-reads its
  model starts from the previous one's age -- with accretion and disk locking a
  stale age can switch the mass change off for the whole run (F77 main.f:791; no
  shipped deck reaches it); a contract overclaim in state/README; card_rescale_targets'
  rows by literal 1-4 at 58 places; check_max_runs blind to a literal (50) in
  controls_buffer; Monte Carlo overwriting cards 2-3's stop ages without the solar
  calibration (F77 main.f:243-244); a few stale sentences.

# Module evaluation 5, 2026-11-01 -- io and state/setup after round 3

Two read-only reviews at yrec-next 6a75b72f, after round 3 ("io and state to 9",
docs/ROADMAP.md), scored out of 10 against rotation/'s bar, as in
audit/module-evaluation-4-2026-10-31 (the first three rounds' 1-5 scores doubled).

| Module | First | Second | Third | Fourth | Now | Now: readability / organization / maintainability / extensibility / docs |
|---|---|---|---|---|---|---|
| io | 4.8 | 6.6 | 7.6 | 7.5 | 8.3 | 8.5 / 8.5 / 8 / 8 / 8.5 |
| state | 4.8 | 6.8 | 7.6 | 7.8 | 8.3 | 8.5 / 8 / 8.5 / 8.5 / 8 |
| setup | 6.6 | 7.4 | 8.2 | 8.2 | 8.8 | 9.5 / 9 / 8.5 / 9 / 8 |

Round 3 did what the fourth reports asked almost item for item (io: 31 of 39 done,
5 partly, 3 open), checked in the code. What keeps io and state below 9, per the reports:
- docs left behind by the merge: setup/README describing the read as it was before
  READ3; state/README not naming state/run_messages.f90; controls_buffer.f90's
  comments (now the only description of each control) naming 5 controls that do not
  exist and 24 F77 routines; generated headers telling the reader to edit the registry;
- tests and checks that check less than they say: test_controls' whole-dump hashes,
  the read-twice test reading one deck twice, check_max_runs comparing a value with
  itself; 57 of 149 star-taking procedures without a contract (outside the four
  required directories);
- open bugs, three of them in the round's own goal (M108, M111, M112), and M117,
  M169, M171-M173; a new one: the gray-atmosphere latch use_ttau_relation is not reset
  when the Monte-Carlo re-arm or a calibration cycle restores atm_choice_in_use;
- CI configured and replayed locally but never run on GitHub (needs a push);
- by design: the buffer's 406 public variables, a domain's control still 4 hand edits
  and 3 generator runs that make does not run, the flat 265-member star type.

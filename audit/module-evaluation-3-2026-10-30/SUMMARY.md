# Module evaluation 3, 2026-10-30 -- summary

The same five read-only reviews as audit/module-evaluation-2026-10-29/ and
audit/module-reevaluation-2026-10-29/, repeated a third time at yrec-next
6e058523 after round 2 (docs/ROADMAP.md, "Round 2: the re-evaluation's
items" and "Third evaluation"), on the same 1-5 scale relative to
rotation/. Each review scored its modules and listed the defects it found;
the 37 items are recorded as open in docs/bug_sweep.md section 9 (M89-M125),
recorded, not fixed.

| Module | First (avg) | After round 1 (avg) | Now (avg) | Now: readability / organization / maintainability / extensibility / docs |
|---|---|---|---|---|
| math | 3.2 | 4.4 | 4.8 | 5 / 5 / 4.5 / 5 / 4.5 |
| numerics | 3.1 | 3.8 | 4.4 | 4.5 / 4.5 / 4 / 4.5 / 4.5 |
| net | 3.4 | 4.0 | 4.3 | 4.5 / 4.5 / 4.5 / 4 / 4 |
| setup | 3.3 | 3.7 | 4.1 | 4.5 / 4 / 4 / 4 / 4 |
| eos | 2.6 | 3.3 | 3.8 | 4 / 4 / 3.5 / 3.5 / 4 |
| io | 2.4 | 3.3 | 3.8 | 4 / 3.5 / 3.5 / 4 / 4 |
| state | 2.4 | 3.4 | 3.8 | 4 / 3.5 / 4 / 4 / 3.5 |

Item 1 was checked in the code by hand when the list was drawn up; every
item was then checked again in the code and against the F77 import
(6cd56739) before its bug_sweep entry was written. All 37 hold, several
with corrections the entries give; none was dropped.

Physics and results (a fix changes results):
1. net: the deuterium burn starts from the start-of-step X and He3 and
   undoes the step's burning, settling and clipping of them (M89).
2. net: the screening's mixture leaves out the metals the Itoh et al.
   mixture counts; about 1.4% on the weak-screening exponent at the solar
   centre (M90).
3. net: the S-factor law's constant 70.62860 carries a CODATA 1973 N_A,
   3.4e-5 below 1/m_u (M91).
4. eos: SCV's general path before the envelope table is built gives NaN
   c_P and grad_ad (zero smix), also after a first-card Z rescale's X = 0
   build (M92).
5. eos: the fully ionized density iteration's failure returns ierr 0 with
   the caller's guess (M93).
6. eos: OPAL 1995 reads X column 6 and temperature row nt+1 past its
   arrays (M94).
7. state: a stop disarmed in one Monte Carlo run stays disarmed in later
   runs (M95).
8. state: after a later card's Z rescale, the convection zones' O16 floor
   keeps the first card's metals (M96).
9. state: a Monte Carlo run's L_sun is 3.844e33 x its target, not the
   deck's (M97).
10. io: the run log opens after the starting and final models, so a failed
    open leaves no run log (M98).

Latent defects and error paths:
11. numerics: a failed natural_spline leaves y2 unset, and none of the 23
    callers passes ierr (M99).
12. numerics: osplin's value on a repeated interior node depends on the
    run of points (M100).
13. numerics: lir gives NaN with ierr 0 when its stencil holds a repeated
    abscissa (M101).
14. numerics: ludcmp calls an all-NaN row singular; the network calls
    simeqc's bad input a singular Jacobian (M102).
15. numerics: hunt's answer for a NaN depends on the guess, against its
    header (M103).
16. kap: OPAL 1992's x_index and rho_search_index are read before they are
    set (M104).
17. eos: a second eos_init keeps each OPAL vintage's first table and leaves
    the new unit open (M105).
18. eos: the Fermi-Dirac reader writes at the file's own indices (M106).
19. eos: SCV's search start idtt, idp is read before it is set (M107).
20. io: a run_kind outside 1-3 is accepted silently (M108).
21. io: an alias wins over a new name set to the default; an array split
    across the two loses the alias's elements (M109).
22. io: write_mod_model has no iostat and no ierr (M110).
23. io: the interior GYRE file name and the isochrone record have five
    digits for the model number (M111).
24. io: direct writes to run_log_unit bypass log_line's guard (M112).
25. io: envelope_print_points (npenv) is read and used nowhere (M113).
26. io/net: net_report skips io_units' queue for messages before the run
    log opens (M114).
27. state: mc_run_start, mc_run_end and the Monte Carlo parameter file are
    not checked (M115).
28. state tools: check_star_contracts and gen_star_docs miss assignments
    after `case (k);` and `real(dp) function` headers (M116).
29. setup: the registry calls fpft_average_weight an averaging weight; it
    is the shape iteration's tolerance, and the copy is right (M117).
30. net: at log T exactly on a cutoff the comparisons disagree, and
    net_configuration documents the helium cutoff the other way (M118).

Tests that check less than they say:
31. eos: test_eos's "OPAL 2006's first derivatives" check tests the 2001
    table (M119).
32. io: test_controls dumps two fields the read does not set and misses six
    scales it does; the pristine test skips the buffer members (M120).

Documentation traps:
33. net: reaction_rates' Sources give q7 with -ln S0; the fits use +ln S0
    (M121).
34. eos: opal_vintage_conflict is produced by nothing (M122).
35. eos: adopt_table_result and blend_table_result drop the ionization
    fractions eos/README tells a new source to fill (M123).
36. math: pow "matches x**y exactly", "within a few ulp of y ln x", and
    "CI runs test_math.py" (M124).

Infrastructure:
37. CI has not run since 2026-09-04, and as configured would fail before
    any test (M125).

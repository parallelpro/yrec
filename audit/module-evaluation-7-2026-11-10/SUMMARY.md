# Module evaluation 7, 2026-11-10 -- every package of src/

Nine read-only reviews at yrec-next 823326a7 (the day it was first pushed), out of 10 against the
standard rotation/ set (10 = at or above it with nothing material left; 9 = only minor items left;
the bar is written out in each report's scope and was the same for every reviewer). For the first
time every package is scored: core, rotation, kap, atm, convection and element_diffusion had not
been scored on this scale before; rotation is scored against the standard as written, not against
itself.

| Package | Fourth (10-31) | Sixth (11-02) | Now | Now: readability / organization / maintainability / extensibility / docs |
|---|---|---|---|---|
| setup | 8.2 | 9.2 | 9.5 | 10 / 9 / 9.5 / 9.5 / 9.5 |
| io | 7.5 | 8.9 | 9.3 | 9.5 / 9 / 9.5 / 9.5 / 9 |
| math | 9.3 | -- | 9.3 | 10 / 10 / 8.5 / 10 / 8 |
| state | 7.8 | 8.9 | 9.1 | 9.5 / 9 / 9 / 9 / 9 |
| numerics | 8.5 | -- | 8.5 | 9 / 9 / 7.5 / 9 / 8 |
| net | 8.2 | -- | 8.1 | 9 / 8.5 / 8.5 / 7.5 / 7 |
| rotation | -- | -- | 8.0 | 7.5 / 8 / 8 / 8.5 / 8 |
| eos | 7.3 | -- | 7.2 | 7.5 / 8 / 7 / 6.5 / 7 |
| core | -- | -- | 6.7 | 7 / 7 / 6.5 / 6.5 / 6.5 |
| convection | -- | -- | 6.4 | 6 / 7 / 6 / 6.5 / 6.5 |
| atm | -- | -- | 6.3 | 6.5 / 6.5 / 6 / 6 / 6.5 |
| kap | -- | -- | 6.2 | 6.5 / 7 / 5.5 / 6 / 6 |
| element_diffusion | -- | -- | 5.9 | 5.5 / 6.5 / 5.5 / 6 / 6 |

core by sub-directory: structure 7.7, step 7.0, diagnostics 7.0, run 6.8, model 6.6, envelope 6.2.

## Reading the table

- The packages the rounds since 10-31 worked on (io, state, setup) are past 9. math holds 9.3.
- numerics, net and eos have not changed since 10-31: their scores hold, net's and eos's docs
  dropped half a point because the round-3 moves (rate scales, builders) left their READMEs stale.
- rotation set the standard and has fallen behind the packages raised to it since: F77 declarations
  (344 `double precision`, 117 `.lt.`-style operators, 245 explicit-shape dummies), a writable
  settings record where kap/eos/net/atm make theirs `protected`, 86 lines of commented-out code in
  live files, capitalized messages naming routines that no longer exist.
- The packages scored for the first time share the architecture (facade, settings record from
  setup/, star-free kernels, all enforced) and the same missing middle layer:
  - no messages module: kernels write to the terminal and the run log themselves (core 139 direct
    writes, atm 115, kap 58; M112);
  - tests that byte-compare printed numbers, with no ok/FAIL checks and no mutation runs (kap, atm,
    convection, element_diffusion), or no unit tests at all for big kernels (core's 1,052-line
    envelope integration and the Henyey elimination);
  - READMEs without the flow, conventions, coverage figure and traced recipes;
  - F77 shape: ALL-CAPS comments (485 lines in convection and element_diffusion), long routines,
    `double precision`, argument-less internal steps.

## New findings worth the authors' attention first (details and file:line in each report)

Bugs (no shipped deck is known to reach most of them; "F77" = inherited):
- core: a star calibration that reaches its 25th pair writes cards 51-52 past the 50-slot arrays
  (F77); `henyey_iterate` returns without setting its divergence flag on its commonest path (works
  only because gfortran keeps the caller's value); the start-model mismatch warning prints an unset
  logical (F77); the C entry cuts namelist paths at 256 characters and closes a host's units 7-99.
- rotation: the composition write-back gives a convection zone beside a radiative end zone the end
  zone's abundance instead of its change (species not conserved; 2026, O13's fix); non-positive
  sub-step controls hang the step (F77); a redone sub-step advances the deuterium rates twice (F77);
  the shape integrand's theta-gradient term lacks a factor r and has the wrong sign (F77; ~1e-3 near
  break-up); `viscos`' electron density makes the ion viscosity ~4% high in the solar core (F77).
- kap/atm: the Alexander 2006 cache can return another composition's opacity with ierr 0 after a
  failed build or a second `kap_init` (2026); 6 of kap's 9 table opens end in a runtime error on a
  missing file and the two OPAL92 opens create it (contradicts M170; F77); an [Fe/H] the Allard file
  lacks dies on a format error (F77); two-Z blending with equal table Zs gives NaN with ierr 0; an
  unknown `atm_choice` silently runs HSRA; every run log says "pure C and N tables" for the C and O
  table; fourteen defects found by the August/September sweeps are open but were never registered in
  bug_sweep (e.g. the conductive opacity's low-density branch, Castelli/Kurucz -999 sentinels above
  ~48 kK, an unset `density_shifted` in OPAL95).
- convection/element_diffusion: semiconvection treats the first mixed zone as the core whatever it
  is (F77 2025 change a8e60674), and nothing pins semiconvection (the cited deck is fully convective);
  starspots' Jacobian terms use the unspotted gradient (F77; convergence only); the original settling
  routines carry four inherited defects (Noerdlinger sign, lnΛ counted twice, midpoint off by one, no
  `abs` in the convergence test); zones past the eleventh convection zone dropped silently;
  `light(3,json)` declared with a literal.
- state: both star tools count an internal WRITE into a star member as a read; `omega_substep` is a
  second member the per-run reset misses (F77; ~1e-14).
- eos: SCV smoothing reads never-built table entries below short rows (F77; up to ~0.5%); OPAL 1995's
  pressure-check error prints two bare numbers; the Fermi-Dirac table's other edges are unchecked.
- io: the run plan labels the rescale targets by position; the old-group-name refusal skips
  `report_error`. numerics/math: two checks never mutation-tested; crmath's |n| < 100 bound untested.

## CI

CI's first GitHub run (42885208151, Linux) passed the build, the generators' checks, the library
link check and 86 of 88 tests; it failed test_rotation:133 (`sum dJ`, a J-conservation residual at
round-off) and test_eos:201 (a finite-difference relative error's fourth digit): both pins of
round-off that another libm moves, not code defects. Both are now ok/FAIL checks (merge 942cd521),
and the round-off columns CI had not yet reached (test_rotation's j/m, X and He4 changes,
test_structure's finite-difference errors) have absolute bounds measured under crmath and +-4 ulp
inputs, used in CI's tolerant mode only; local runs still compare byte for byte.

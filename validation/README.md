# Reference data of the physics validation tier

Small files that the tier's tests (`test_validation_*.py`,
docs/PHYSICS_VALIDATION.md) compare YREC with. Nothing here was made by
YREC. No test fetches anything from the network: these copies were
made offline, from the files MESA 26.04.1 ships, on 2026-11-17.

## solar/

| File | What | Source | How it was made |
|---|---|---|---|
| `bison_l0.txt` | The Sun's radial (l = 0) p-mode frequencies, n = 6 to 28, with their 1-sigma errors (uHz) | BiSON: Broomhall et al. 2009, MNRAS 396, L100 | Copied from `$MESA_DIR/astero/test_suite/surface_effects/inlist_astero_search_controls` (`freq_target(0,:)`, `freq_sigma(0,:)`); the radial order n = index + 5 was assigned by matching YREC's and MESA's solar models (3033.754 uHz is n = 21). MESA's file does not name its source; the values are those of Broomhall et al.'s table as far as could be checked offline. |
| `sound_speed_mdi.txt` | The Sun's sound speed c(r) at 79 radii, 0.09 to 0.99 R, with its 1-sigma error, from the inversion of SOHO/MDI mean-multiplet frequencies | Rhodes, Kosovichev, Scherrer, Schou & Reiter 1997, Solar Physics 175, 287 (their table 6b) | Columns 1, 3, 4 and 7 of `$MESA_DIR/data/star_data/solar_csound.data`, whose header gives the source. |

The scalar references (the convection zone's depth, the surface helium,
Z/X, the nominal L, R and Teff) are in `test_validation_solar.py`'s
header, each with its source.

The cross-code reference (MESA's track and frequencies for the 1 Msun
benchmark) is with its case, `sample_cases/solar_vs_MESA/mesa/`, and the
rotating solar deck's landmark rotation profiles are in
`rotation/` here; each has its own note on how it was made.

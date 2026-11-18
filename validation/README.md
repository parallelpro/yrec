# Reference data of the physics validation tier

Small files that the tier's tests (`test_validation_*.py`,
docs/PHYSICS_VALIDATION.md) compare YREC with. The files in `solar/`
are published results, not made by YREC; the file in `rotation/` is
YREC's own output, a regression pin by tolerance. No test fetches
anything from the network: the solar files were copied on 2026-11-17
from the files MESA 26.04.1 ships, and each was compared with its
published source the same day (the "Checked" column).

## solar/

| File | What | Source | How it was made | Checked against the source |
|---|---|---|---|---|
| `bison_l0.txt` | The Sun's radial (l = 0) p-mode frequencies, n = 6 to 28, with their 1-sigma errors (uHz) | BiSON, two tables joined at 2000 uHz: n = 6 to 13, Davies, Broomhall, Chaplin, Elsworth & Hale 2014, MNRAS 439, 2025; n = 14 to 28, Broomhall et al. 2009, MNRAS 396, L100, table 3 (corrected to the quiet Sun; table 1 is the raw fit) | Copied from `$MESA_DIR/astero/test_suite/surface_effects/inlist_astero_search_controls` (`freq_target(0,:)`, `freq_sigma(0,:)`; that test's `README.rst` names the two papers). The radial order n = index + 5 is the papers' (3033.754 uHz is n = 21). | n = 14 to 28: 15 of 15 frequencies and errors equal to table 3 of arXiv:0903.5219v1 (the journal's tables are images). n = 6 to 13: 8 of 8 equal to the table on the journal's page (doi:10.1093/mnras/stu080). |
| `sound_speed_mdi.txt` | The Sun's sound speed c(r) at 79 radii, 0.09 to 0.99 R, with its standard error, from the inversion of SOHO/MDI mean-multiplet frequencies (144 days) | Rhodes, Kosovichev, Scherrer, Schou & Reiter 1997, Solar Physics 175, 287: the electronic table 6b of the paper's SOHO CD-ROM version, https://soho.nascom.nasa.gov/publications/CDROM1/papers/rhodes/table6b.txt | Columns 1, 3, 4 and 7 of `$MESA_DIR/data/star_data/solar_csound.data`, whose header gives the source. | All 79 rows of the four columns equal to the file at that address. It is the inversion's published result (the inferred c, not the reference Model S's, which is the table's column 2). The printed article was not read (paywalled). |

The scalar references (the convection zone's depth, the surface helium,
Z/X, the nominal L, R and Teff) are in `test_validation_solar.py`'s
header, each with its source; docs/PHYSICS_VALIDATION.md, "Sources",
lists every published value the tier uses with where it was read.

## rotation/

| File | What | How it was made |
|---|---|---|
| `omega_landmarks.txt` | The angular velocity of `testsuite/Test_solarGS98_difrotmixper18p4d1p0` at five ages and six mass coordinates | **YREC's own values**, written by `test_validation_rotation.py` with `YREC_VALIDATION_WRITE_LANDMARKS=1` (its header names the commit). Not a reference from outside the code: the test compares a run with them at 1%, which says that the rotation profile has not moved, not that it is right. |

## Elsewhere

The cross-code reference (MESA's track and frequencies for the 1 Msun
benchmark) is with its case, `sample_cases/solar_vs_MESA/mesa/`, whose
README says how it was made.

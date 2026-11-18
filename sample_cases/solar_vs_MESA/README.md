# One 1 Msun track in YREC and in MESA

The cross-code benchmark of the physics validation tier
(docs/PHYSICS_VALIDATION.md, section 4). A 1 Msun star with the physics
of YREC's standard solar deck, as closely as MESA can be set to it:

| | YREC | MESA 26.04.1 |
|---|---|---|
| Composition | GS98, X = 0.719620, Z = 0.016492 | the same (`initial_zfracs = 3`) |
| Mixing length | 1.707575 | 1.707575, `MLT_option = 'Henyey'` |
| Element diffusion, rotation, overshoot | none | none |
| Atmosphere | Eddington gray T(tau) | `atm_option = 'T_tau'`, Eddington, varying opacity |
| Opacity | OP (`GS98.OP17`), Ferguson et al. 2005 below, Potekhin conduction | `OP_gs98`, `lowT_fa05_gs98`, MESA's conduction |
| Equation of state | OPAL 2006, SCV where OPAL ends | `use_OPAL_SCVH` (Skye and FreeEOS off) |
| Nuclear rates | YREC's (`use_new_nuclear_rates`) | MESA's defaults (`basic.net`) |
| Start | the deuterium birthline model (log R = 0.63) | a pre-main-sequence model (`pre_ms_T_c = 9.9d5`) |
| G, Msun | 6.67252e-8, a star of 1.98900e33 g | 6.67430e-8, 1.98841e33 g (G M equal to 3e-5) |

| File | |
|---|---|
| `one_msun_no_diffusion.inlist` | The YREC deck: the standard solar deck with its stop at central hydrogen 0.05, its opacity table at the model's Z, and the testsuite's solver tolerances (its header says why). |
| `mesa/inlist`, `mesa/inlist_project`, `mesa/*_columns.list` | The MESA side. `inlist_stop_solar_age` (run a, to 4.568 Gyr, writing the last model's GYRE file) or `inlist_stop_track` (run b, to central hydrogen 0.05) is copied to `inlist_stop`. |
| `mesa/history_track.data` | Run b's history, every model (247), the columns the comparison reads, 12 digits. |
| `mesa/history_solar_age.data` | Run a's last history row (4.568 Gyr). |
| `mesa/gyre_solar_age_modes.txt` | GYRE 9.1.1's l = 0 and 1 modes of run a's last model (`mesa/gyre_solar_age.in`). |
| `mesa/inlist_stop_opacity_experiment` | Run b with MESA's opacity raised by 2.5% between log T = 6.3 and 6.75 (below). |
| `compare.py` | Prints YREC minus MESA at matched central hydrogen and at the solar age. Standard library only. |

`test_validation_mesa.py` (repository root) runs YREC and makes the
comparison against the stored MESA files, with the bands below.

## Running it

```
cd sample_cases/solar_vs_MESA
mkdir -p yrec/output && cp one_msun_no_diffusion.inlist yrec/
(cd yrec && YREC_INPUT=$YREC/input YREC_START=$YREC/startmodels YREC_OUTPUT=output \
   $YREC/src/yrec one_msun_no_diffusion.inlist)            # 45 s, 515 models
python3 -I compare.py yrec/output/history.data             # against the stored MESA track

# to remake the MESA side (2.5 minutes a run):
cp -R $MESA_DIR/star/work mesa_run && cp mesa/inlist mesa/inlist_project mesa/*_columns.list mesa_run/
cp mesa/inlist_stop_track mesa_run/inlist_stop
(cd mesa_run && ./mk && ./rn)
python3 -I compare.py yrec/output/history.data mesa_run/LOGS/history.data
```

## Results (2026-11-17, YREC valid-tier at 1179afbf, MESA 26.04.1, GYRE 9.1.1)

YREC minus MESA (dex; the age relative, the central hydrogen absolute):

| At | log L | log Teff | log R | log T_c | log rho_c | age | X_c |
|---|---|---|---|---|---|---|---|
| X_c = 0.6 | -0.0094 | -0.0031 | +0.0016 | -0.0032 | -0.0046 | +2.6% | |
| X_c = 0.35 | -0.0076 | -0.0030 | +0.0022 | -0.0026 | -0.0009 | +4.0% | |
| X_c = 0.1 | +0.0111 | -0.0021 | +0.0098 | +0.0030 | +0.0336 | +7.7% | |
| 4.568 Gyr | -0.0141 | -0.0034 | -0.0003 | -0.0058 | -0.0145 | | +0.0148 |

At 4.568 Gyr YREC's model is the Sun's (log L = +0.00009, log R =
-0.000002: the deck's X and mixing length are a solar calibration);
MESA's is 3.3% brighter and 45 K hotter with the same X, Z and mixing
length.

Frequencies at 4.568 Gyr (GYRE on each code's model with its
atmosphere): YREC's are higher by 1.8 uHz at n = 6 to 5.2 uHz at
n = 21 (l = 0; 1.6e-3 to 2.3e-3, the same for l = 1); the radial large
separation over n = 15..25 is 135.756 against 135.715 uHz (+3.1e-4).
Both codes' sound-travel `delta_nu` agree to 0.54% (142.66 and 141.89).

Bands of `test_validation_mesa.py` (about twice the largest measured):
log L 0.030, log Teff 0.007, log R 0.020, log T_c 0.012, log rho_c
0.070, the age 15%, X_c at the solar age 0.030, delta_nu 1%, a
frequency 0.5%, the large separation 0.2%.

## What explains the differences

Measured, not assumed:

- **Opacity: more than half of the luminosity offset.** Along the two
  4.568 Gyr models, at equal temperature and corrected to equal density
  (kappa ~ rho^0.5), YREC's opacity is higher by 0.3% at log T = 7.0,
  0.8% at 6.8, 1.8% at 6.7, 2.5% at 6.6 and 3.1% at 6.35, just below the
  convection zone. Part of it is in the tables themselves: at the nodes
  the two codes' OP tables share (X = 0.7, Z = 0.02, log R = -1.5 and
  -2), YREC's `GS98.OP17` equals MESA's `OP_gs98` at log T >= 7.0 and is
  0.7 to 1.6% higher between log T = 6.2 and 6.8. With MESA's opacity
  raised by 2.5% between log T = 6.3 and 6.75
  (`inlist_stop_opacity_experiment`), MESA's log L falls by 0.0059 dex at
  X_c = 0.6, 0.0055 at 0.35 and 0.0076 at 4.568 Gyr (63%, 73% and 54% of
  the offsets above), its log Teff by 0.0012 (40%), its log T_c by
  0.0011 at X_c = 0.6 (33%), and its age at X_c = 0.35 rises by 1.3%
  (of 4.0%). Which table is the OP data as published was not determined.
- **OP or OPAL in MESA: a quarter of that size.** MESA with `gs98` (OPAL)
  in place of `OP_gs98` is fainter by 0.0030 dex at X_c = 0.6 and 0.0046
  at 4.568 Gyr; log Teff moves by 0.0004.
- **The equation of state: not a cause.** P/(rho T) at equal temperature
  along the two models agrees to 3e-4 through the chemically uniform
  envelope (log T = 6.9 to 5.5).
- **G and the mass: 3e-4 in L.** G M is equal to 3e-5; at equal G M the
  luminosity scales as 1/M (YREC's star is 3.0e-4 more massive).
- **Not separated:** the nuclear rates and screening, the opacity above
  log T = 6.9 (0.3 to 0.4% higher in YREC), the conductive opacity, the
  mixing-length formulations (they move Teff and R, not L) and the
  codes' difference schemes. Together they are the rest: 0.0035 dex in
  log L at X_c = 0.6, 0.0019 in log Teff.
- **X_c = 0.1 amplifies.** Near central exhaustion the central hydrogen
  is a poor clock (it falls ever more slowly while the shell around the
  centre burns), so the offsets at X_c = 0.1 (the age +7.7%, the central
  density +8%) are larger than at matched age and are dominated by how
  centrally concentrated each code's burning is, which follows from the
  0.7% difference in central temperature.

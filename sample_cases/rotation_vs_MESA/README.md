# Rotation against MESA: a solid body under J conservation

A 1 Msun star (X = 0.71221, Z = 0.018126, GS98, mixing length 1.825)
rotates as a solid body from the pre-main sequence to 1 Gyr, with no
wind, no disk locking, no mass loss and no settling. Its angular
momentum J is then conserved and its rotation rate is omega = J/I. The
two codes' rotation curves agree to the extent that their moments of
inertia I(t) do, so the comparison checks YREC's J-to-omega inversion
(`omega_from_j`, with the rotational distortion) and the structure it
reads, against MESA's.

| File | |
|---|---|
| `solid_body_no_wind.inlist` | The YREC deck: `force_solid_body_rotation`, no wind (`use_wind_torque`), no disk locking, no settling, `max_age(2) = 1.0D9`. |
| `mesa/inlist_project`, `mesa/history_columns.list` | The MESA side (26.04.1): a pre-main-sequence model (`pre_ms_T_c = 3d5`), rigid rotation through a uniform angular momentum viscosity of 1e20 cm^2/s, no element diffusion. |
| `compare.py` | Prints each code's J conservation and, against age, YREC over MESA for omega, I = J/omega, R and k^2 = J/(omega M R^2). Needs numpy. |

The codes start from different models: YREC at log R = 0.625, MESA at
log R = 0.284. MESA's `new_omega` gives it YREC's J there (k^2 is
about 0.2 in both on the Hayashi track): 1.804029e-5 rad/s for YREC's
`initial_omega` of 3.95e-6 (fast: omega/omega_crit about 0.25 on the
main sequence), 4.567162e-7 rad/s for 1.0e-7 (slow: about 0.006, where
the distortion is negligible).

## Running it

```
cd sample_cases/rotation_vs_MESA
mkdir -p yrec/output && cp solid_body_no_wind.inlist yrec/
(cd yrec && YREC_INPUT=$YREC/input YREC_START=$YREC/startmodels YREC_OUTPUT=output \
   $YREC/src/yrec solid_body_no_wind.inlist)
cp -R $MESA_DIR/star/work mesa_run && cp mesa/* mesa_run/
(cd mesa_run && ./mk && ./rn)
python compare.py yrec/output/history.data mesa_run/LOGS/history.data
```

## Results (2026-10-28, YREC yrec-next, MESA 26.04.1)

J conservation over the run: YREC within 2.3e-5 (J moves where the
mesh is rezoned, by up to 1e-6 per model near the main sequence),
MESA within 3e-12. MESA's centre and surface rotate
together to 4e-8.

YREC over MESA minus one, at the same J:

| Age (yr) | omega, slow | omega, fast | I, fast | R, fast | k^2, fast |
|---|---|---|---|---|---|
| 3e6 | -0.210 | -0.208 | +0.261 | +0.117 | +0.011 |
| 1e7 | -0.035 | -0.031 | +0.030 | -0.003 | +0.037 |
| 3e7 | -0.014 | -0.007 | +0.006 | -0.003 | +0.011 |
| 4e7 | +0.014 | +0.018 | -0.019 | -0.039 | +0.063 |
| 5e7 | +0.001 | +0.004 | -0.005 | -0.024 | +0.045 |
| 1e8 | -0.008 | -0.004 | +0.003 | -0.023 | +0.050 |
| 3e8 | -0.002 | +0.001 | -0.002 | -0.022 | +0.043 |
| 1e9 | -0.002 | +0.000 | -0.001 | -0.023 | +0.045 |

- **On the main sequence** (50 Myr to 1 Gyr) the surface omega agrees
  to 0.8% (slow) and 0.4% (fast), because the moments of inertia
  agree to 0.8% and 0.5%.
- **On the pre-main sequence** the differences follow the radius: the
  two stars contract on different schedules from their different
  starting models (YREC's is 12% larger at 3 Myr), and at 40 Myr they
  reach the main sequence at slightly different ages.
- **k^2 and R** differ by 4-5% and -2.3% on the main sequence, in
  opposite directions, so that I = k^2 M R^2 agrees to 0.5%: the stars
  differ mostly in their outer layers, which move R but hold little of
  the moment of inertia. (Why the radii differ, among the codes'
  opacities, equations of state and atmospheres, was not examined.)
- **The fast case agrees as well as the slow one**, so the rotational
  distortion of the moment of inertia (Law 1980 in YREC) is not a
  source of difference at omega/omega_crit = 0.25.

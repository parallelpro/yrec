#!/usr/bin/env python3
"""YREC against MESA on one 1 Msun track (README.md).

    python3 -I compare.py YREC_history.data [MESA_history.data]

Prints YREC minus MESA at matched central hydrogen (0.6, 0.35, 0.1) and
at the solar age (4.568 Gyr): log L, log Teff, log R, the central
temperature and density (dex), and the age at each central hydrogen
(relative). The MESA history defaults to the stored reference,
mesa/history_track.data. Standard library only; test_validation_mesa.py
imports the functions.
"""
import pathlib
import sys

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src" / "tools"))
import validate_model as vm  # noqa: E402

SOLAR_AGE = 4.568e9
MATCHED_XC = (0.6, 0.35, 0.1)
# quantity: (YREC's history column, MESA's, the difference printed)
QUANTITIES = {
    "log_L": ("log_L", "log_L", "dex"),
    "log_Teff": ("log_Teff", "log_Teff", "dex"),
    "log_R": ("log_R", "log_R", "dex"),
    "log_center_T": ("log_center_T", "log_center_T", "dex"),
    "log_center_Rho": ("log_center_Rho", "log_center_Rho", "dex"),
    "star_age": ("star_age", "star_age", "relative"),
    "center_h1": ("center_h1", "center_h1", "absolute"),
}
# the main sequence: the first crossing after this age (the central
# hydrogen is not monotonic before, while He3 and C12 come to equilibrium)
MAIN_SEQUENCE_FROM = 5.0e7


def at(history, key, value, name):
    """Column `name` where column `key` first crosses `value` on the main
    sequence, linear between the two models around it."""
    xs, ages = history[key], history["star_age"]
    for i in range(1, len(xs)):
        a, b = xs[i - 1], xs[i]
        if ages[i] < MAIN_SEQUENCE_FROM or a == b:
            continue
        if (a - value) * (b - value) <= 0.0:
            f = (value - a) / (b - a)
            return history[name][i - 1] + f * (history[name][i] - history[name][i - 1])
    raise ValueError(f"{key} never crosses {value}")


def differences(yrec, mesa, key, value, names):
    """YREC minus MESA for each named quantity where `key` = `value`
    (relative for the age)."""
    out = {}
    for name in names:
        y_col, m_col, kind = QUANTITIES[name]
        y, m = at(yrec, key, value, y_col), at(mesa, key, value, m_col)
        out[name] = y / m - 1.0 if kind == "relative" else y - m
    return out


def compare(yrec, mesa):
    """{point label: {quantity: YREC - MESA}} at the matched points."""
    structure = ["log_L", "log_Teff", "log_R", "log_center_T", "log_center_Rho"]
    out = {}
    for xc in MATCHED_XC:
        out[f"X_c = {xc}"] = differences(yrec, mesa, "center_h1", xc, structure + ["star_age"])
    out["4.568 Gyr"] = differences(yrec, mesa, "star_age", SOLAR_AGE, structure + ["center_h1"])
    return out


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    _, yrec = vm.read_mesa_table(argv[1])
    _, mesa = vm.read_mesa_table(argv[2] if len(argv) > 2 else HERE / "mesa" / "history_track.data")
    table = compare(yrec, mesa)
    names = ["log_L", "log_Teff", "log_R", "log_center_T", "log_center_Rho", "star_age", "center_h1"]
    print("YREC minus MESA (dex; the age relative; X_c absolute)")
    print(f"{'':12s}" + "".join(f"{n:>16s}" for n in names))
    for label, row in table.items():
        print(f"{label:12s}" + "".join(f"{row[n]:+16.5f}" if n in row else f"{'':16s}" for n in names))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

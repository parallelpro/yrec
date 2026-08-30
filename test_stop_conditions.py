"""End-to-end tests for the MESA-style structure-limit stop conditions
(2026): log_L / Teff / log_g / nu_max upper+lower limits checked by
check_stop_conditions after every converged model, and the seismic
scaling-relation observables behind them (nu_max, delta_nu_rho -- optional
history columns, off by default).

Each case runs the standard solar model with one limit set so the
starting structure already violates it: the run must end the kind
cards almost immediately (a clean exit-0 with a short history), and
the STOP diagnostic must appear in the run log.
"""
import math
import os
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
CONVERTER = REPO / "src" / "tools" / "upgrade_inlist.py"
CASE_DIR = REPO / "examples" / "run_standard_solar_model"
CASE = "Test_solar_noGS_norot"
ENV = {"YREC_INPUT": REPO / "input", "YREC_START": REPO / "startmodels"}


def _run_with_controls(tmp_path, extra_controls, extra_star_job=""):
    """Convert the solar pair, splice extra lines into the namelist
    groups, run."""
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    conv = tmp_path / "inlist"
    r = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE_DIR / f"{CASE}.nml1"),
         str(CASE_DIR / f"{CASE}.nml2"), "-o", str(conv)],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    text = conv.read_text()
    text = text.replace("&controls", "&controls\n" + extra_controls, 1)
    if extra_star_job:
        text = text.replace("&star_job", "&star_job\n" + extra_star_job, 1)
    conv.write_text(text)
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env.update({k: str(v) for k, v in ENV.items()})
    r = subprocess.run([str(YREC), "inlist"], cwd=tmp_path, env=env,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]
    logs = list((tmp_path / "output").glob("*.log"))
    assert logs, "no run log produced"
    log = logs[0].read_text()
    hist = (tmp_path / "output" / "history.data").read_text().splitlines()
    return log, hist


def _assert_stopped_early(log, hist, needle):
    assert needle in log, f"expected '{needle}' STOP diagnostic in run.log"
    # header block is 7 lines; a limit hit at the first checked model of
    # each card leaves only a handful of history rows (vs ~460 for the
    # full solar run)
    assert len(hist) < 40, f"run did not stop early: {len(hist)} history lines"


# per-kind-card arrays, like the other stopping criteria: the solar
# deck's card 2 is the evolving card, so the limits are set on (2)
@pytest.mark.parametrize("control,needle", [
    ("Teff_upper_limit(2) = 4.0d3", "STOP: Teff"),
    ("log_L_lower_limit(2) = 2.0d0", "STOP: log_L"),
    ("nu_max_lower_limit(2) = 5.0d3", "STOP: nu_max"),
    ("log_g_upper_limit(2) = 3.0d0", "STOP: log_g"),
])
def test_structure_limit_stops_run(tmp_path, control, needle):
    log, hist = _run_with_controls(tmp_path, f" {control}\n")
    _assert_stopped_early(log, hist, needle)


def test_seismic_columns_are_opt_in(tmp_path):
    cols = tmp_path / "hist_cols"
    cols.write_text(
        "model_number\nlog_Teff\nnu_max\ndelta_nu_rho\ndelta_nu\ndelta_Pg\n")
    log, hist = _run_with_controls(
        tmp_path, " nu_max_lower_limit(2) = 5.0d3\n",
        extra_star_job=f" history_columns_file = '{cols}'\n")
    _assert_stopped_early(log, hist, "STOP: nu_max")
    names = hist[5].split()
    assert names == ["model_number", "log_Teff", "nu_max", "delta_nu_rho",
                     "delta_nu", "delta_Pg"], names
    row = dict(zip(names, hist[-1].split()))
    nu_max = float(row["nu_max"])
    delta_nu_rho = float(row["delta_nu_rho"])
    delta_nu = float(row["delta_nu"])
    delta_pg = float(row["delta_Pg"])
    # a ~1 Msun pre-main-sequence/early model: positive, sub-solar
    # nu_max (large radius), physically plausible values
    assert 0.0 < nu_max < 5.0e3, nu_max
    assert 0.0 < delta_nu_rho < 2.0e2, delta_nu_rho
    # the sound-travel-time integral should land near the mean-density
    # scaling estimate (same star, two estimators)
    assert 0.0 < delta_nu < 2.0e2, delta_nu
    assert 0.5 < delta_nu / delta_nu_rho < 2.0, (delta_nu, delta_nu_rho)
    # g-mode cavity: zero when the early model has no radiative
    # region; a nearly-convective pre-MS star has a vanishing cavity
    # integral, so the spacing can be legitimately enormous -- only
    # sign and finiteness are guaranteed here (the solar-age value is
    # checked against the Sun in the commit's validation run)
    assert delta_pg >= 0.0 and math.isfinite(delta_pg), delta_pg


def test_seismic_columns_absent_by_default(tmp_path):
    log, hist = _run_with_controls(tmp_path, " Teff_upper_limit(2) = 4.0d3\n")
    _assert_stopped_early(log, hist, "STOP: Teff")
    names = hist[5].split()
    assert "nu_max" not in names and "delta_nu_rho" not in names, (
        "seismic columns must be opt-in via history_columns_file")

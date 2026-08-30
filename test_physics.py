"""Loose-tolerance physics assertions on the standard solar model.

The byte-identical pins (test_all.py, test_mesa_baselines.py) are the
refactoring gate: they catch ANY change, including intentional ones,
and are reseeded deliberately when output changes on purpose. This
test is their complement: it asserts that the solar model is still
physically a solar model -- surface luminosity/radius/Teff near solar
at the solar age, a plausibly depleted central hydrogen abundance --
with tolerances wide enough to survive legitimate development
(output-format changes, small physics refinements) and tight enough
to catch a real physics regression (wrong tables, broken solver,
mis-wired composition).

Runs the noGS_norot case through the new-style inlist (converted on
the fly, profiles off) exactly like test_mesa_baselines, then reads
the final history row.
"""
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


def test_solar_model_physics(tmp_path):
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    conv = tmp_path / "inlist"
    r = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE_DIR / f"{CASE}.nml1"),
         str(CASE_DIR / f"{CASE}.nml2"), "-o", str(conv)],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env.update({k: str(v) for k, v in ENV.items()})
    r = subprocess.run([str(YREC), "inlist"], cwd=tmp_path, env=env,
                       capture_output=True, text=True, timeout=1200)
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]

    hist = (tmp_path / "output" / "history.data").read_text().splitlines()
    names = hist[5].split()
    row = dict(zip(names, hist[-1].split()))

    age = float(row["star_age"])
    log_l = float(row["log_L"])
    log_r = float(row["log_R"])
    log_teff = float(row["log_Teff"])
    xc = float(row["center_h1"])
    nz = int(row["num_zones"])

    # the run must reach the configured solar-age stop
    assert age > 4.5e9, f"stopped early: age = {age:.3e} yr"
    # a 1 Msun GS98 model at solar age is the Sun to a few percent
    assert abs(log_l) < 0.02, f"log L/Lsun = {log_l}"
    assert abs(log_r) < 0.02, f"log R/Rsun = {log_r}"
    assert abs(log_teff - 3.7615) < 0.010, f"log Teff = {log_teff}"
    # central hydrogen depletion after 4.57 Gyr of pp burning
    assert 0.30 < xc < 0.42, f"center X = {xc}"
    # sane mesh
    assert 500 < nz < 3000, f"num_zones = {nz}"

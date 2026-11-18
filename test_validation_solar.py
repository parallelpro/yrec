"""Solar truth checks (physics validation tier, piece 3).

The calibrated standard solar model with the GS98 mixture and element
diffusion (testsuite/Test_solarGS98_dif) at the solar age, against the
Sun: observed and helioseismic values with their uncertainties. No
standard solar model matches all of them at 1 sigma, so each check
states a band and the doc gives the measured offset beside the
literature's typical one (docs/PHYSICS_VALIDATION.md section 3).

  quantity            observed                 source
  L                   3.828e33 erg/s           IAU 2015 Resolution B3 (nominal)
  R                   6.957e10 cm              IAU 2015 Resolution B3 (nominal)
  Teff                5772 K                   IAU 2015 Resolution B3 (nominal)
  age                 4.568 Gyr                the deck's (meteoritic, Bouvier & Wadhwa 2010: 4.5682)
  R_cz/R              0.713 +- 0.001           Basu & Antia 1997, MNRAS 287, 189
  surface Y           0.2485 +- 0.0035         Basu & Antia 2004, ApJ 606, L85
  surface Z/X         0.0229 (10%)             Grevesse & Sauval 1998, Space Sci. Rev. 85, 161
  c(r)                MDI inversion            Rhodes et al. 1997 (validation/solar/sound_speed_mdi.txt)
  nu(l = 0, n)        BiSON, n = 6..28         Broomhall et al. 2009 (validation/solar/bison_l0.txt)

The frequencies need GYRE (GYRE_DIR); those tests skip without it.
test_checks_reject_a_model_without_diffusion shows the checks can fail:
the standard solar deck without diffusion (examples/
run_standard_solar_model) has the Sun's L, R and Teff and fails the
convection zone's depth, the surface helium and the sound speed.

Local (two decks to the solar age, twice each: about 4 minutes).
"""
import bisect
import math
import os
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
sys.path.insert(0, str(REPO / "src" / "tools"))
import validate_model as vm  # noqa: E402
from conftest import INPUT_ROOT  # noqa: E402

# this repository's own build (to validate another build, run these tests
# from its tree); the tables and starting models are conftest.py's roots
YREC = REPO / "src" / "yrec"

# not in CI (minutes, the full input tree): deselect with `pytest -m "not local"`
pytestmark = pytest.mark.local
DATA = REPO / "validation" / "solar"

OBSERVED = {
    "L_cgs": 3.828e33, "R_cm": 6.957e10, "Teff_K": 5772.0,
    "R_cz": (0.713, 0.001), "Y_surface": (0.2485, 0.0035), "Z_over_X": 0.0229,
}

# The bands (measured value beside each, 2026-11-17) and what sets them.
BAND = {
    # the deck is a solar calibration (its X, Z and mixing length were
    # fitted to L and R to 3e-5); later changes of the code moved it
    "L": 1.0e-3,                # measured +6.5e-5
    "R": 1.0e-3,                # measured -1.6e-4
    "Teff_K": 1.0,              # measured +0.5 K
    # three times the helioseismic 1-sigma
    "R_cz_sigma": 3.0,          # measured -0.3 sigma (0.7127)
    "Y_surface_sigma": 3.0,     # measured -1.5 sigma (0.2434)
    # the calibration's target, 0.02292 for GS98
    "Z_over_X": 1.0e-2,         # measured -0.02% of 0.0229
    # the profile's own base of the convection zone (gradr = grada,
    # linear between two shells 9e-4 R apart) against the history's column
    "R_cz_rederived": 1.0e-4,   # measured below 1e-6
    # c(r) against the inversion, 0.1 <= r/R <= 0.95 (71 points): a
    # standard model with the GS98 mixture is at about 1e-3 rms
    "sound_speed_rms": 2.0e-3,  # measured 8.9e-4
    "sound_speed_max": 5.0e-3,  # measured 2.7e-3 (below the convection zone)
    # radial modes: the low orders, least affected by the surface layers
    "low_order_relative": 2.0e-3,   # measured +3.9e-4 to +9.4e-4, n = 6..10
    # the large separation, n = 15..25, model over observed: the surface
    # effect raises a model's high-order frequencies
    "large_separation": (1.000, 1.020),   # measured 1.0093
}
MEASURED = {}


def record(name, value):
    MEASURED[name] = value
    return value


def _data_rows(name):
    return [[float(x) for x in line.split()] for line in (DATA / name).read_text().splitlines()
            if line.strip() and not line.startswith("#")]


def _final_model(tmp, inlist_text):
    out = vm.run_with_final_profile(YREC, tmp, inlist_text, pulse_formats=("GYRE",))
    _, hist = vm.read_mesa_table(out["history"] / "history.data")
    header, cols = vm.read_gyre_text(out["GYRE"] / "profile1.data.GYRE")
    _, prof = vm.read_profile(out["GYRE"] / "profile1.data")
    return {"hist": {k: v[-1] for k, v in hist.items()}, "header": header, "cols": cols,
            "prof": prof, "pulse": out["GYRE"] / "profile1.data.GYRE", "tmp": tmp}


def _need_yrec(table):
    if not YREC.exists():
        pytest.skip("build src/yrec first")
    if not table.exists():
        pytest.skip(f"needs the full input/ tree ({table})")


@pytest.fixture(scope="module")
def sun(tmp_path_factory):
    _need_yrec(INPUT_ROOT / "eos" / "opal2006" / "EOSOPAL06Z0.018776")
    deck = REPO / "testsuite" / "Test_solarGS98_dif.inlist"
    return _final_model(tmp_path_factory.mktemp("solar"), deck.read_text())


@pytest.fixture(scope="module")
def sun_without_diffusion(tmp_path_factory):
    _need_yrec(INPUT_ROOT / "eos" / "opal2006" / "EOSOPAL06Z0.016492")
    tmp = tmp_path_factory.mktemp("solar_nodif")
    case = REPO / "examples" / "run_standard_solar_model"
    conv = subprocess.run(
        [sys.executable, str(REPO / "src" / "tools" / "upgrade_inlist.py"),
         str(case / "Test_solar_noGS_norot.nml1"), str(case / "Test_solar_noGS_norot.nml2"),
         "-o", str(tmp / "inlist_conv")], capture_output=True, text=True)
    assert conv.returncode == 0, conv.stderr
    return _final_model(tmp / "run", (tmp / "inlist_conv").read_text())


# ---------------------------------------------------------------------
# the checks, each of one model
# ---------------------------------------------------------------------
def check_global(m, tag=""):
    """L, R and Teff at the solar age."""
    row = m["hist"]
    assert row["star_age"] == pytest.approx(4.568e9, rel=1e-9)
    # the history's log_L and log_R are in the deck's units, the nominal ones
    assert abs(record(tag + "L/Lsun - 1", m["header"]["L_star"] / OBSERVED["L_cgs"] - 1.0)) <= BAND["L"]
    assert abs(record(tag + "R/Rsun - 1", m["header"]["R_star"] / OBSERVED["R_cm"] - 1.0)) <= BAND["R"]
    assert abs(record(tag + "Teff - 5772 K", 10.0 ** row["log_Teff"] - OBSERVED["Teff_K"])) <= BAND["Teff_K"]


def convection_zone_base(m):
    """r/R where gradr = grada at the base of the surface convection
    zone, from the profile (linear between the two shells around it)."""
    p = m["prof"]
    nz = m["hist"]["num_zones"]
    k = max(j for j in range(1, nz) if p["mixing_type"][j - 1] == 0 and p["mixing_type"][j] == 1)
    below = p["gradr"][k - 1] - p["grada"][k - 1]
    above = p["gradr"][k] - p["grada"][k]
    log_r = p["logR"][k - 1] - below / (above - below) * (p["logR"][k] - p["logR"][k - 1])
    return 10.0 ** log_r / m["header"]["R_star"]


def check_convection_zone_depth(m, tag=""):
    value, sigma = OBSERVED["R_cz"]
    # the history's cz_bot_radius is in Rsun
    r_cz = record(tag + "R_cz/R (history)", m["hist"]["cz_bot_radius"] / 10.0 ** m["hist"]["log_R"])
    assert abs(record(tag + "R_cz/R, profile's gradr = grada minus history",
                      convection_zone_base(m) - r_cz)) <= BAND["R_cz_rederived"]
    assert abs(record(tag + "R_cz: (model - 0.713)/0.001", (r_cz - value) / sigma)) <= BAND["R_cz_sigma"]


def check_surface_helium(m, tag=""):
    value, sigma = OBSERVED["Y_surface"]
    y = record(tag + "surface Y", m["hist"]["surface_he4"])
    assert abs(record(tag + "Y: (model - 0.2485)/0.0035", (y - value) / sigma)) <= BAND["Y_surface_sigma"]


def check_surface_metals(m, tag=""):
    row = m["hist"]
    assert row["surface_z_div_x"] == pytest.approx(row["surface_z"] / row["surface_h1"], rel=1e-12)
    assert abs(record(tag + "(Z/X)/0.0229 - 1", row["surface_z_div_x"] / OBSERVED["Z_over_X"] - 1.0)) \
        <= BAND["Z_over_X"]


def sound_speed_differences(m):
    """(c_model - c_sun)/c_sun at the inversion's radii, the model's
    sqrt(Gamma_1 P/rho) interpolated linearly in r/R."""
    x = [r / m["header"]["R_star"] for r in m["cols"]["r"]]
    c = vm.sound_speed(m["cols"])
    out = []
    for r0, c_sun, _sigma, _width in _data_rows("sound_speed_mdi.txt"):
        j = bisect.bisect_left(x, r0)
        c_model = c[j - 1] + (c[j] - c[j - 1]) * (r0 - x[j - 1]) / (x[j] - x[j - 1])
        out.append((r0, c_model / c_sun - 1.0))
    return out


def check_sound_speed(m, tag=""):
    inner = [d for r0, d in sound_speed_differences(m) if 0.1 <= r0 <= 0.95]
    assert len(inner) == 71
    rms = record(tag + "sound speed, rms dc/c, 0.1 to 0.95 R", math.sqrt(sum(d * d for d in inner) / len(inner)))
    worst = record(tag + "sound speed, max |dc/c|, 0.1 to 0.95 R", max(abs(d) for d in inner))
    assert rms <= BAND["sound_speed_rms"], rms
    assert worst <= BAND["sound_speed_max"], worst


CHECKS = [check_global, check_convection_zone_depth, check_surface_helium, check_surface_metals,
          check_sound_speed]


@pytest.mark.parametrize("check", CHECKS, ids=lambda f: f.__name__)
def test_solar_model(sun, check):
    check(sun)


def test_radial_modes(sun):
    """The model's radial modes (GYRE) against BiSON's."""
    if vm.gyre_program() is None:
        pytest.skip("GYRE_DIR does not name a GYRE installation (bin/gyre)")
    cols = vm.run_gyre(sun["pulse"], "GYRE", sun["tmp"] / "gyre_run",
                       [(0, "LINEAR", 800, 4500, 400)], tag="p")
    model = vm.mode_table(cols, 0)
    observed = {int(n): f for n, f, _sigma in _data_rows("bison_l0.txt")}
    assert sorted(observed) == list(range(6, 29))
    low = [model[n] / observed[n] - 1.0 for n in range(6, 11)]
    record("nu(model)/nu(BiSON) - 1, n = 6", low[0])
    record("nu(model)/nu(BiSON) - 1, n = 10", low[-1])
    record("nu(model) - nu(BiSON), n = 21 (uHz)", model[21] - observed[21])
    assert max(abs(d) for d in low) <= BAND["low_order_relative"], low
    dnu_model, _ = vm.large_separation(model, 15, 25)
    dnu_observed, _ = vm.large_separation(observed, 15, 25)
    record("large separation, n = 15..25, model (uHz)", dnu_model)
    record("large separation, n = 15..25, BiSON (uHz)", dnu_observed)
    lo, hi = BAND["large_separation"]
    assert lo <= record("large separation, model/BiSON", dnu_model / dnu_observed) <= hi


def test_checks_reject_a_model_without_diffusion(sun_without_diffusion):
    """The standard solar model without element diffusion is calibrated to
    the same L, R and Teff, and fails the three helioseismic checks: the
    convection zone is too shallow, the surface helium is the initial
    one, and the sound speed is off by several times the band."""
    m = sun_without_diffusion
    check_global(m, tag="no diffusion: ")
    for check in (check_convection_zone_depth, check_surface_helium, check_sound_speed):
        with pytest.raises(AssertionError):
            check(m, tag="no diffusion: ")


def test_report(sun):
    """Print the measured values (pytest -s), and write them to
    YREC_VALIDATION_REPORT when it names a file."""
    text = "".join(f"{value:14.6f}  {name}\n" for name, value in sorted(MEASURED.items()))
    print("\nmeasured:\n" + text)
    path = os.environ.get("YREC_VALIDATION_REPORT")
    if path:
        pathlib.Path(path).write_text(text)

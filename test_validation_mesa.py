"""Cross-code benchmark, YREC against MESA (physics validation tier, piece 4).

One 1 Msun track with the two codes' physics matched as closely as they
allow (sample_cases/solar_vs_MESA: GS98, X = 0.719620, Z = 0.016492,
mixing length 1.707575, no diffusion, no rotation, no overshoot, an
Eddington gray atmosphere, OP opacities, OPAL equations of state). YREC
is run here; MESA 26.04.1's side is the stored reference under
sample_cases/solar_vs_MESA/mesa/ (README.md there says how it was made
and how to remake it), so the test needs no MESA.

Compared: log L, log Teff, log R, the central temperature and density
at matched central hydrogen (0.6, 0.35, 0.1) and at 4.568 Gyr, the age
at each central hydrogen, both codes' sound-travel delta_nu, and (with
GYRE) the two models' frequencies at 4.568 Gyr.

The bands are about twice the measured differences. The two codes are
independent implementations with different opacity tables, nuclear
rates and difference schemes, so the bands say how far apart two
correct codes are here, and the README lists what each known difference
explains. A band is not a pin: YREC may move inside it.

Local (three YREC runs, about 2 minutes); the frequency comparison
skips without GYRE_DIR.
"""
import os
import pathlib
import sys

import pytest

from conftest import INPUT_ROOT

REPO = pathlib.Path(__file__).parent
CASE = REPO / "sample_cases" / "solar_vs_MESA"
sys.path.insert(0, str(REPO / "src" / "tools"))
sys.path.insert(0, str(CASE))
import validate_model as vm  # noqa: E402
import compare  # noqa: E402

# this repository's own build (to validate another build, run these tests
# from its tree); the tables and starting models are conftest.py's roots
YREC = REPO / "src" / "yrec"

# not in CI (minutes): deselect with `pytest -m "not local"`
pytestmark = pytest.mark.local

# |YREC - MESA| allowed (dex; the age relative), and the largest measured
# at the three central hydrogens and at 4.568 Gyr (2026-11-17).
BAND = {
    "log_L": 0.030,            # measured -0.0141 (4.568 Gyr), +0.0111 (X_c = 0.1)
    "log_Teff": 0.007,         # measured -0.0034
    "log_R": 0.020,            # measured +0.0098 (X_c = 0.1)
    "log_center_T": 0.012,     # measured -0.0058 (4.568 Gyr)
    "log_center_Rho": 0.070,   # measured +0.0336 (X_c = 0.1)
    "star_age": 0.15,          # measured +7.7% (X_c = 0.1)
    "center_h1": 0.030,        # measured +0.0148 (4.568 Gyr)
    # the same integral in both codes, [2 int_0^R dr/c]^-1, at 4.568 Gyr
    "delta_nu": 0.010,         # measured +5.4e-3
    # frequencies at 4.568 Gyr, l = 0 and 1, n = 5..31, relative
    "frequency": 0.005,        # measured 1.6e-3 to 2.3e-3 (+1.8 to +5.2 uHz)
    # the radial large separation over n = 15..25
    "large_separation": 0.002,  # measured +3.1e-4
}
MEASURED = {}


def record(name, value):
    MEASURED[name] = value
    return value


@pytest.fixture(scope="module")
def mesa():
    _, track = vm.read_mesa_table(CASE / "mesa" / "history_track.data")
    _, solar = vm.read_mesa_table(CASE / "mesa" / "history_solar_age.data")
    modes = {}
    for line in (CASE / "mesa" / "gyre_solar_age_modes.txt").read_text().splitlines():
        if not line.startswith("#"):
            ell, n_pg, _n_p, _n_g, freq = line.split()
            modes[(int(ell), int(n_pg))] = float(freq)
    return {"track": track, "solar": {k: v[-1] for k, v in solar.items()}, "modes": modes}


@pytest.fixture(scope="module")
def yrec(tmp_path_factory):
    if not YREC.exists():
        pytest.skip("build src/yrec first")
    if not (INPUT_ROOT / "opacity" / "atomic" / "legacy" / "opal95" / "GS98.OP17").exists():
        pytest.skip("needs the standard solar deck's tables under input/")
    tmp = tmp_path_factory.mktemp("mesa")
    deck = (CASE / "one_msun_no_diffusion.inlist").read_text()
    track_text = vm.edit_inlist(deck, "&control", ["history_columns_file = 'history_all.list'"])
    _, track = vm.read_mesa_table(vm.run_yrec(YREC, tmp / "track", track_text) / "history.data")
    stop = " max_age(2) = 1.2D10"
    assert stop in deck
    out = vm.run_with_final_profile(YREC, tmp / "solar_age", deck.replace(stop, " max_age(2) = 4.568D9"))
    _, solar = vm.read_mesa_table(out["history"] / "history.data")
    return {"track": track, "solar": {k: v[-1] for k, v in solar.items()},
            "pulse": out["GYRE"] / "profile1.data.GYRE", "tmp": tmp}


def test_tracks_at_matched_points(yrec, mesa):
    """log L, log Teff, log R, the central T and rho at matched central
    hydrogen and at the solar age; the age at each central hydrogen."""
    table = compare.compare(yrec["track"], mesa["track"])
    assert list(table) == ["X_c = 0.6", "X_c = 0.35", "X_c = 0.1", "4.568 Gyr"]
    for label, row in table.items():
        for name, diff in row.items():
            record(f"{label}: {name}", diff)
            assert abs(diff) <= BAND[name], (label, name, diff)


def test_both_runs_reach_the_same_states(yrec, mesa):
    """The comparison is of the same phase: both tracks pass central
    hydrogen 0.05 after the solar age, on the main sequence, and both
    stars keep their mass and (no diffusion) their surface helium."""
    for side in (yrec, mesa):
        h = side["track"]
        assert h["center_h1"][-1] < 0.05 < h["center_h1"][-2]
        assert h["star_age"][-1] > compare.SOLAR_AGE
        assert h["surface_he4"][-1] == pytest.approx(0.26386, abs=5e-5)
    assert yrec["solar"]["star_age"] == pytest.approx(mesa["solar"]["star_age"], rel=1e-9)


def test_delta_nu_at_the_solar_age(yrec, mesa):
    """Both codes print the same sound-travel integral."""
    ratio = yrec["solar"]["delta_nu"] / mesa["solar"]["delta_nu"] - 1.0
    assert abs(record("4.568 Gyr: delta_nu, YREC/MESA - 1", ratio)) <= BAND["delta_nu"]
    # MESA's own acoustic radius gives its delta_nu
    assert 1.0e6 / (2.0 * mesa["solar"]["acoustic_radius"]) == pytest.approx(mesa["solar"]["delta_nu"], rel=1e-9)


def test_frequencies_at_the_solar_age(yrec, mesa):
    """GYRE on YREC's model against GYRE on MESA's (stored), l = 0 and 1."""
    if vm.gyre_program() is None:
        pytest.skip("GYRE_DIR does not name a GYRE installation (bin/gyre)")
    cols = vm.run_gyre(yrec["pulse"], "GYRE", yrec["tmp"] / "gyre",
                       [(0, "LINEAR", 800, 4500, 400), (1, "LINEAR", 800, 4500, 400)], tag="yrec")
    for ell in (0, 1):
        mine = vm.mode_table(cols, ell)
        theirs = {n: f for (l_, n), f in mesa["modes"].items() if l_ == ell and n > 0}
        orders = sorted(set(mine) & set(theirs))
        assert orders == list(range(5, 32)), (ell, orders)
        rel = [mine[n] / theirs[n] - 1.0 for n in orders]
        record(f"4.568 Gyr: l = {ell}, largest |nu(YREC)/nu(MESA) - 1|", max(abs(d) for d in rel))
        record(f"4.568 Gyr: l = {ell}, n = 21, nu(YREC) - nu(MESA) (uHz)", mine[21] - theirs[21])
        assert max(abs(d) for d in rel) <= BAND["frequency"], (ell, rel)
        if ell == 0:
            dnu_y, _ = vm.large_separation(mine, 15, 25)
            dnu_m, _ = vm.large_separation(theirs, 15, 25)
            assert abs(record("4.568 Gyr: large separation n = 15..25, YREC/MESA - 1",
                              dnu_y / dnu_m - 1.0)) <= BAND["large_separation"]


def test_comparison_can_fail(yrec, mesa):
    """A track whose luminosity is 10% high (0.041 dex), or whose star is
    3% larger, leaves the bands; so does a model 1% denser on the mean
    (its frequencies scale as the square root of the mean density)."""
    bright = dict(yrec["track"])
    bright["log_L"] = [x + 0.0414 for x in bright["log_L"]]
    table = compare.compare(bright, mesa["track"])
    assert any(abs(row["log_L"]) > BAND["log_L"] for row in table.values())
    large = dict(yrec["track"])
    large["log_R"] = [x + 0.0128 for x in large["log_R"]]
    table = compare.compare(large, mesa["track"])
    assert any(abs(row["log_R"]) > BAND["log_R"] for row in table.values())
    theirs = {n: f for (l_, n), f in mesa["modes"].items() if l_ == 0 and n > 0}
    dense = {n: f * 1.01 ** 0.5 for n, f in theirs.items()}
    assert max(abs(dense[n] / theirs[n] - 1.0) for n in theirs) > BAND["frequency"] * 0.99
    assert abs(vm.large_separation(dense, 15, 25)[0] / vm.large_separation(theirs, 15, 25)[0] - 1.0) \
        > BAND["large_separation"]


def test_report(yrec, mesa):
    """Print the measured values (pytest -s), and write them to
    YREC_VALIDATION_REPORT when it names a file."""
    text = "".join(f"{value:+12.5f}  {name}\n" for name, value in MEASURED.items())
    print("\nmeasured (YREC minus MESA):\n" + text)
    path = os.environ.get("YREC_VALIDATION_REPORT")
    if path:
        pathlib.Path(path).write_text(text)

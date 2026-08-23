"""MESA-style output acceptance test (2026).

Runs the solar noGS/norot case twice from the same converted inlist:
once with the stamped `use_legacy_output = .true.` (legacy .track) and
once with the stamp removed (MESA mode -> CASE.history in MESA's
history.data layout). Asserts:
  - MESA mode writes the .history file with the documented layout
    (data starts at line 7; line 6 is the column-name row) and one row
    per converged model;
  - the legacy .track row count matches;
  - model_number / star_age / log_L / log_Teff agree numerically with
    the legacy .track columns (same physics, different format);
  - MESA mode suppresses the per-shell .store stream.

Local test: needs the full input/ tree.
"""
import os
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
CONVERTER = REPO / "src" / "tools" / "upgrade_inlist.py"
CASE = REPO / "examples" / "run_standard_solar_model"
NML1 = "Test_solar_noGS_norot.nml1"
NML2 = "Test_solar_noGS_norot.nml2"
ENV = {"YREC_INPUT": REPO / "input", "YREC_START": REPO / "startmodels"}


def _run(workdir, inlist_text):
    workdir.mkdir()
    (workdir / "output").mkdir()
    (workdir / "inlist").write_text(inlist_text)
    env = dict(os.environ)
    env.update({k: str(v) for k, v in ENV.items()})
    r = subprocess.run([str(YREC), "inlist"], cwd=workdir, env=env,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]
    return workdir / "output"


def test_mesa_history_matches_legacy_track(tmp_path):
    if not YREC.exists():
        pytest.skip("build src/yrec first")
    if not (REPO / "input").exists():
        pytest.skip("needs the full input/ tree (local only)")

    inlist = tmp_path / "inlist_conv"
    conv = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE / NML1), str(CASE / NML2),
         "-o", str(inlist)], capture_output=True, text=True)
    assert conv.returncode == 0, conv.stderr
    stamped = inlist.read_text()
    assert "use_legacy_output = .true." in stamped
    unstamped = "\n".join(l for l in stamped.splitlines()
                          if "use_legacy_output" not in l) + "\n"

    legacy_out = _run(tmp_path / "legacy", stamped)
    mesa_out = _run(tmp_path / "mesa", unstamped)

    base = NML1.rsplit(".", 1)[0]
    track = legacy_out / f"{base}.track"
    history = mesa_out / f"{base}.history"
    assert history.exists(), sorted(p.name for p in mesa_out.iterdir())

    # legacy .track rows (skip # comments and the header-name row)
    track_rows = [l.split() for l in track.read_text().splitlines()
                  if l.strip() and not l.lstrip().startswith("#")
                  and not l.lstrip().startswith("Step")]
    hist_lines = history.read_text().splitlines()
    name_row = hist_lines[5].split()
    assert name_row[0] == "model_number" and "log_L" in name_row
    hist_rows = [l.split() for l in hist_lines[6:] if l.strip()]
    assert len(hist_rows) == len(track_rows), (len(hist_rows), len(track_rows))

    icol = {n: i for i, n in enumerate(name_row)}
    for tr, hr in zip(track_rows, hist_rows):
        assert int(hr[icol["model_number"]]) == int(tr[0])
        # star_age is years in MESA convention; .track column 3 is Gyr
        # .track prints 8 significant digits; .history carries 16, so
        # compare at the track's own precision.
        assert abs(float(hr[icol["star_age"]]) / 1e9 - float(tr[2])) <= \
            1e-7 * max(1.0, abs(float(tr[2])))
        for name, j in (("log_L", 3), ("log_Teff", 6)):
            a, b = float(hr[icol[name]]), float(tr[j])
            assert abs(a - b) <= 1e-7 * max(1.0, abs(b)), (name, a, b)

    # MESA mode suppresses the legacy per-shell stream
    store = mesa_out / f"{base}.store"
    if store.exists():
        body = [l for l in store.read_text().splitlines()
                if l.strip() and not l.startswith("#")]
        assert not any(l.startswith("MOD") for l in body), \
            ".store still receiving per-model blocks in MESA mode"

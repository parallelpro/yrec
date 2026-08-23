"""MESA-style output acceptance test (2026).

Runs the solar noGS/norot case from one converted inlist in three
configurations and checks the whole MESA-mode output contract:

  legacy   (stamped use_legacy_output = .true.)  -- the oracle .track
  mesa     (stamp removed)                       -- default MESA mode
  custom   (mesa + history_columns_file)         -- column selection

Assertions: MESA mode writes exactly {history.data, profile*.data,
CASE.log}; the history file has MESA's layout (data from line 7,
names on line 6) with one row per converged model matching the
legacy .track numerically; profiles appear every profile_interval
models with zone 1 = the surface and num_zones rows; the history
profile_number column maps models to profiles (YREC's replacement
for profiles.index); a history_columns_file selects exactly the
requested columns in the requested order.

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
BASE = NML1.rsplit(".", 1)[0]


def _run(workdir, inlist_text, extra_files=None):
    workdir.mkdir()
    (workdir / "output").mkdir()
    (workdir / "inlist").write_text(inlist_text)
    for name, text in (extra_files or {}).items():
        (workdir / name).write_text(text)
    env = dict(os.environ)
    env.update({k: str(v) for k, v in ENV.items()})
    r = subprocess.run([str(YREC), "inlist"], cwd=workdir, env=env,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]
    return workdir / "output"


def _parse_mesa_file(path):
    lines = path.read_text().splitlines()
    names = lines[5].split()
    rows = [l.split() for l in lines[6:] if l.strip()]
    icol = {n: i for i, n in enumerate(names)}
    return names, rows, icol


def test_mesa_output_contract(tmp_path):
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

    # ---- exact file set ----
    produced = sorted(p.name for p in mesa_out.iterdir())
    profiles = [p for p in produced if p.startswith("profile")]
    assert produced == sorted([f"{BASE}.log", "history.data"] + profiles), \
        produced
    assert profiles, "no profile files written (profile_interval default)"

    # ---- history vs the legacy track ----
    track_rows = [l.split() for l in
                  (legacy_out / f"{BASE}.track").read_text().splitlines()
                  if l.strip() and not l.lstrip().startswith("#")
                  and not l.lstrip().startswith("Step")]
    names, hist_rows, icol = _parse_mesa_file(mesa_out / "history.data")
    assert names[0] == "model_number" and "profile_number" in names
    assert len(hist_rows) == len(track_rows)
    for tr, hr in zip(track_rows, hist_rows):
        assert int(float(hr[icol["model_number"]])) == int(tr[0])
        assert abs(float(hr[icol["star_age"]]) / 1e9 - float(tr[2])) <= \
            1e-7 * max(1.0, abs(float(tr[2])))
        for name, j in (("log_L", 3), ("log_Teff", 6)):
            a, b = float(hr[icol[name]]), float(tr[j])
            assert abs(a - b) <= 1e-7 * max(1.0, abs(b)), (name, a, b)

    # ---- profile_number column maps models to profiles ----
    expect_num = 0
    for hr in hist_rows:
        model = int(float(hr[icol["model_number"]]))
        pnum = int(float(hr[icol["profile_number"]]))
        if model % 50 == 0:   # profile_interval default
            expect_num += 1
            assert pnum == expect_num, (model, pnum, expect_num)
        else:
            assert pnum == 0, (model, pnum)
    assert expect_num == len(profiles)

    # ---- profile layout: zone 1 = surface, num_zones rows ----
    pnames, prows, pcol = _parse_mesa_file(mesa_out / "profile1.data")
    assert pnames[0] == "zone" and "logT" in pnames and "h1" in pnames
    plines = (mesa_out / "profile1.data").read_text().splitlines()
    glob_names = plines[1].split()
    glob_vals = plines[2].split()
    nz = int(glob_vals[glob_names.index("num_zones")])
    assert len(prows) == nz
    assert int(float(prows[0][pcol["zone"]])) == 1
    assert int(float(prows[-1][pcol["zone"]])) == nz
    # surface first: mass decreases from zone 1 toward the center row
    assert float(prows[0][pcol["mass"]]) > float(prows[-1][pcol["mass"]])

    # ---- custom history columns + FGONG pulse output ----
    columns = "! test subset\nstar_age\nlog_Teff\nlog_L\n"
    custom = unstamped.replace(
        "&star_job",
        "&star_job\n history_columns_file = 'my_columns.list'\n"
        " pulse_gyre_interval = 100\n", 1)
    custom = custom.replace(
        "&controls",
        "&controls\n pulse_format = 'FGONG'\n", 1)
    custom_out = _run(tmp_path / "custom", custom,
                      extra_files={"my_columns.list": columns})
    cnames, crows, _ = _parse_mesa_file(custom_out / "history.data")
    assert cnames == ["star_age", "log_Teff", "log_L"], cnames
    assert len(crows) == len(hist_rows)

    fgongs = sorted((tmp_path / "custom").glob("fgong_profile_*.fgong"))
    assert fgongs, "no FGONG files written"
    flines = fgongs[0].read_text().splitlines()
    nn, iconst, ivar, ivers = (int(x) for x in flines[4].split())
    assert (iconst, ivar, ivers) == (15, 40, 300)
    data = flines[5:]
    n_glob_lines = (iconst + 4) // 5
    n_var_lines = (ivar + 4) // 5
    assert len(data) == n_glob_lines + nn * n_var_lines,         (len(data), nn, n_glob_lines, n_var_lines)

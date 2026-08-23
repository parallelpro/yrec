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
    assert names[:3] == ["model_number", "profile_number", "num_zones"]
    # integer id columns are written as true integers
    assert "." not in hist_rows[0][0] and "." not in hist_rows[0][2]
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

    # ---- extended model: profiles reach the top of the atmosphere ----
    # The interior grid stops at the fitting point (T ~ 1e4-1e6 K);
    # profiles must continue through the envelope + atmosphere, so the
    # outermost point sits at/above the photosphere with T near (in
    # fact below) Teff.  Guards against the truncation bug where
    # zone 1 was the fitting point (T = 15145 K vs Teff = 4436 K).
    prof_model = next(int(float(hr[icol["model_number"]])) for hr in hist_rows
                      if int(float(hr[icol["profile_number"]])) == 1)
    log_teff_1 = next(float(hr[icol["log_Teff"]]) for hr in hist_rows
                      if int(float(hr[icol["model_number"]])) == prof_model)
    t_outer = 10.0 ** float(prows[0][pcol["logT"]])
    assert t_outer < 1.2 * 10.0 ** log_teff_1, (t_outer, 10 ** log_teff_1)

    # ---- custom history columns + FGONG pulse output ----
    columns = "! test subset\nstar_age\nlog_Teff\nlog_L\n"
    custom = unstamped.replace(
        "&star_job",
        "&star_job\n history_columns_file = 'my_columns.list'\n", 1)
    custom = custom.replace(
        "&controls",
        "&controls\n pulse_format = 'FGONG'\n"
        " write_profile_flag = .false.\n write_pulse_flag = .true.\n", 1)
    custom_out = _run(tmp_path / "custom", custom,
                      extra_files={"my_columns.list": columns})
    cnames, crows, _ = _parse_mesa_file(custom_out / "history.data")
    assert cnames == ["star_age", "log_Teff", "log_L"], cnames
    assert len(crows) == len(hist_rows)

    # write_profile_flag = .false. -> pulse-only mode: pulse files
    # carry the profile numbering, no profile<N>.data
    prof_data = [p.name for p in custom_out.glob("profile*.data")]
    assert not prof_data, prof_data
    fgongs = sorted(custom_out.glob("profile*.data.FGONG"))
    assert fgongs, "no FGONG files written"
    flines = fgongs[0].read_text().splitlines()
    nn, iconst, ivar, ivers = (int(x) for x in flines[4].split())
    assert (iconst, ivar, ivers) == (15, 40, 300)
    data = flines[5:]
    n_glob_lines = (iconst + 4) // 5
    n_var_lines = (ivar + 4) // 5
    assert len(data) == n_glob_lines + nn * n_var_lines, \
        (len(data), nn, n_glob_lines, n_var_lines)

    # extended model in the pulse file too: the first point (FGONG is
    # surface-to-center) is the top of the atmosphere -- radius at or
    # above the photospheric R_star (glob 2), temperature below Teff
    # (glob 14), not the fitting-point value.
    def fgong_vals(lines):
        vals = []
        for l in lines:
            l = l.rstrip()
            vals += [float(l[i:i + 16]) for i in range(0, len(l), 16)]
        return vals
    glob = fgong_vals(data[:n_glob_lines])[:iconst]
    surf = fgong_vals(data[n_glob_lines:n_glob_lines + n_var_lines])[:ivar]
    r_star, teff = glob[1], glob[13]
    assert surf[0] >= 0.999 * r_star, (surf[0], r_star)
    assert surf[2] < 1.2 * teff, (surf[2], teff)


def test_default_columns_lists_in_sync():
    """defaults/{history,profile}_columns.list must list exactly the
    writers' column tables (they are the user-facing documentation of
    what history_columns_file / profile_columns_file may contain)."""
    import re
    src = (REPO / "src" / "io" / "yrec_output.f90").read_text()

    def harvest(sub):
        m = re.search(r"subroutine " + sub + r"\(names\)(.*?)end subroutine " + sub,
                      src, re.S)
        return [n for _, n in
                re.findall(r"names\((\d+)\)\s*=\s*'([^']+)'", m.group(1))]

    def listed(fname):
        out = []
        for line in (REPO / "src" / "defaults" / fname).read_text().splitlines():
            line = line.split("!")[0].strip()
            if line:
                out.append(line)
        return out

    assert listed("history_columns.list") == harvest("history_column_names")
    assert listed("profile_columns.list") == harvest("profile_column_names")


def test_gsm_pulse_output(tmp_path):
    """GSM (GYRE-HDF5) pulse output -- opt-in: needs an HDF5-enabled
    build (`make USE_HDF5=1` in src/) and YREC_TEST_HDF5=1 in the
    environment, since the default build compiles the GSM writer as a
    reporting stub."""
    if os.environ.get("YREC_TEST_HDF5") != "1":
        pytest.skip("set YREC_TEST_HDF5=1 (and build with USE_HDF5=1)")
    if not YREC.exists() or not (REPO / "input").exists():
        pytest.skip("needs src/yrec and the full input/ tree")
    inlist = tmp_path / "inlist_conv"
    conv = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE / NML1), str(CASE / NML2),
         "-o", str(inlist)], capture_output=True, text=True)
    assert conv.returncode == 0, conv.stderr
    text = "\n".join(l for l in inlist.read_text().splitlines()
                     if "use_legacy_output" not in l) + "\n"
    text = text.replace(
        "&controls",
        "&controls\n pulse_format = 'GSM'\n write_pulse_flag = .true.\n"
        " write_profile_flag = .false.\n", 1)
    out = _run(tmp_path / "gsm", text)
    gsms = sorted(out.glob("profile*.data.GSM"))
    assert gsms, "no GSM files written"
    blob = gsms[0].read_bytes()[:4]
    assert blob == b"\x89HDF", blob
    try:
        import h5py
    except ImportError:
        return
    with h5py.File(gsms[0], "r") as f:
        n = int(f.attrs["n"][0]) if hasattr(f.attrs["n"], "__len__") \
            else int(f.attrs["n"])
        assert int(f.attrs["version"][0] if hasattr(f.attrs["version"],
                   "__len__") else f.attrs["version"]) == 101
        assert len(f["r"]) == n
        assert sorted(f.keys()) == sorted(
            ["r", "M_r", "L_r", "P", "T", "rho", "nabla", "N2", "Gamma_1",
             "nabla_ad", "delta", "kap", "kap_kap_T", "kap_kap_rho",
             "eps", "eps_eps_T", "eps_eps_rho", "Omega_rot"])

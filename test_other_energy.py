"""The extra energy source hook (src/core/structure/energy_sources.f90).

Runs examples/other_energy/yrec_extra_heat (built by `make other_energy`
in src/) on the config-matrix base deck, in a scratch directory:
  - with eps_extra = 0 its outputs are byte-identical to the plain yrec's
    (a registered source that adds nothing changes nothing);
  - with eps_extra > 0 the star is brighter and the history column
    lum_extra is positive.
"""
import os
import pathlib
import re
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
HEATED = REPO / "examples" / "other_energy" / "yrec_extra_heat"
DECK = REPO / "examples" / "run_config_matrix" / "cm_base.inlist"
VERSION = re.compile(r"^# YREC v\S* \([0-9a-f]+(-dirty)?\)$", re.M)


def _run(exe, workdir, eps):
    workdir.mkdir()
    (workdir / "output").mkdir()
    deck = DECK.read_text()
    deck = deck.replace("&control", '&control\n history_columns_file = "history_columns.list"', 1)
    (workdir / "deck.inlist").write_text(deck)
    (workdir / "history_columns.list").write_text("model_number\nlog_L\nlog_Teff\nlum_extra\n")
    (workdir / "extra_heat.nml").write_text(
        f"&extra_heat\n eps_extra = {eps}\n m_inner_msun = 0.2d0\n m_outer_msun = 0.5d0\n/\n")
    env = dict(os.environ)
    env["YREC_OUTPUT"] = "output"
    result = subprocess.run([str(exe), "deck.inlist"], cwd=workdir, env=env,
                            capture_output=True, text=True)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    return workdir / "output"


def _history(outdir):
    lines = (outdir / "history.data").read_text().splitlines()
    for i, line in enumerate(lines):
        if line.split()[:1] == ["model_number"]:
            names = line.split()
            rows = [l.split() for l in lines[i + 1:] if l.strip()]
            return {n: [float(r[k]) for r in rows] for k, n in enumerate(names)}
    raise AssertionError("no column header in history.data")


def test_extra_energy_source(tmp_path):
    if not HEATED.exists() or not YREC.exists():
        pytest.skip("build src/yrec and `make other_energy` in src/ first")
    plain = _run(YREC, tmp_path / "plain", "0.0d0")
    zero = _run(HEATED, tmp_path / "zero", "0.0d0")
    heated = _run(HEATED, tmp_path / "heated", "1.0d0")

    files = sorted(p.name for p in plain.iterdir())
    assert files == sorted(p.name for p in zero.iterdir())
    for name in files:
        a, b = (plain / name).read_bytes(), (zero / name).read_bytes()
        if a != b:
            a, b = VERSION.sub("", a.decode()), VERSION.sub("", b.decode())
        assert a == b, f"{name}: a zero source changed the output"

    h0, h1 = _history(plain), _history(heated)
    assert all(v == 0.0 for v in h0["lum_extra"])
    assert h1["lum_extra"][-1] > 0.0
    assert h1["log_L"][-1] > h0["log_L"][-1]

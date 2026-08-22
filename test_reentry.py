"""Phase-five step C acceptance test (2026 -- ROADMAP.md): re-entrancy.

The contract: two run_yrec calls in one process produce, on the second
call, outputs byte-identical to a fresh single-process run. src/
test_reentry runs the engine twice with the same namelists; this test
runs a fresh `yrec` once in one working directory and `test_reentry`
in another, then byte-compares every output file (version-hash line
excluded, as in the Stage-0 discipline).

Uses the solar noGS/norot case (~18 s), which terminates NORMALLY
(run count exhausted). The m0030 case cannot serve here: it
historically ends via the BSSTEP numerics stop, which kills the
process mid-call -- the documented residual (numerics gates are
opt-in) means such cases are not re-enterable until their callers
pass ierr. Three engine runs, ~1 minute total. Local test: needs
the full input/ tree (not in CI's sparse checkout).
"""
import os
import pathlib
import shutil
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
REENTRY = REPO / "src" / "test_reentry"
CASE = REPO / "examples" / "run_standard_solar_model"
NML1 = "Test_solar_noGS_norot.nml1"
NML2 = "Test_solar_noGS_norot.nml2"


def _run(binary, workdir):
    (workdir / "output").mkdir()
    for name in (NML1, NML2):
        shutil.copy(CASE / name, workdir / name)
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(binary), NML1, NML2], cwd=workdir,
        env=env, capture_output=True, text=True, timeout=600,
    )
    return result


def _strip(path):
    return [l for l in path.read_text(errors="replace").splitlines()
            if not l.startswith("# YREC v")]


def test_second_call_equals_fresh_process(tmp_path):
    if not YREC.exists() or not REENTRY.exists():
        pytest.skip("build src/yrec and src/test_reentry first")
    fresh = tmp_path / "fresh"
    twice = tmp_path / "twice"
    fresh.mkdir()
    twice.mkdir()
    r1 = _run(YREC, fresh)
    assert r1.returncode == 0, r1.stdout[-500:] + r1.stderr[-500:]
    r2 = _run(REENTRY, twice)
    assert r2.returncode == 0, r2.stdout[-500:] + r2.stderr[-500:]
    assert "second run done" in r2.stdout

    fresh_files = sorted(p.name for p in (fresh / "output").iterdir())
    twice_files = sorted(p.name for p in (twice / "output").iterdir())
    assert fresh_files == twice_files, (fresh_files, twice_files)
    diffs = []
    for name in fresh_files:
        if _strip(fresh / "output" / name) != _strip(twice / "output" / name):
            diffs.append(name)
    assert not diffs, (
        f"second in-process run differs from a fresh process in: {diffs} "
        "-- some SAVEd state is escaping the reset (see "
        "core/yrec_reset.f90 and the ROADMAP step-C notes)"
    )

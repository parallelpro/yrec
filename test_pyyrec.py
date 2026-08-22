"""libyrec/pyyrec acceptance test (2026 -- ROADMAP.md).

Runs the solar noGS/norot case three ways -- once through the yrec CLI
binary in a fresh process, then TWICE through pyyrec inside this very
Python process -- and byte-compares every output file of both library
runs against the CLI run (version-hash line excluded, the Stage-0
discipline). That proves the milestone's three claims at once: the
shared library runs the engine, the namelist paths reach parmin via
the override channel (not argv -- python's argv would be nonsense to
getarg), and the second in-process call is byte-identical
(re-entrancy through the C API).

The CLI run is the oracle rather than the case's standard/ because
the .short stream echoes the input paths verbatim: standard/ was
pinned with the harness's relative-path placeholder substitution,
while this test (like test_reentry.py) uses the YREC_INPUT/YREC_START
environment overrides -- same physics, different echoed strings.

Local test: needs the full input/ tree (not in CI's sparse checkout),
~1 min for the three runs. Build first: `make` and `make lib` in src/.
"""
import os
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
CASE = REPO / "examples" / "run_standard_solar_model"
NML1 = "Test_solar_noGS_norot.nml1"
NML2 = "Test_solar_noGS_norot.nml2"
ENV = {"YREC_INPUT": REPO / "input", "YREC_START": REPO / "startmodels"}

sys.path.insert(0, str(REPO))
import pyyrec  # noqa: E402


def _strip(path):
    return [l for l in path.read_text(errors="replace").splitlines()
            if not l.startswith("# YREC v")]


def _prepare(workdir):
    workdir.mkdir()
    (workdir / "output").mkdir()
    for name in (NML1, NML2):
        (workdir / name).write_text((CASE / name).read_text())


def _run_cli(workdir):
    _prepare(workdir)
    env = dict(os.environ)
    env.update({k: str(v) for k, v in ENV.items()})
    result = subprocess.run(
        [str(YREC), NML1, NML2], cwd=workdir,
        env=env, capture_output=True, text=True, timeout=600,
    )
    assert result.returncode == 0, result.stdout[-500:] + result.stderr[-500:]
    return workdir / "output"


def _run_lib(workdir):
    _prepare(workdir)
    ierr = pyyrec.run(NML1, NML2, cwd=workdir, env=ENV)
    assert ierr == 0
    return workdir / "output"


def _compare(label, out, oracle):
    oracle_files = sorted(p.name for p in oracle.iterdir())
    got_files = sorted(p.name for p in out.iterdir())
    assert got_files == oracle_files, (label, got_files, oracle_files)
    diffs = [n for n in oracle_files
             if _strip(oracle / n) != _strip(out / n)]
    assert not diffs, (
        f"{label} pyyrec run differs from the CLI run in: {diffs}"
    )


def test_two_inprocess_runs_match_cli(tmp_path):
    try:
        pyyrec.library_path()
    except OSError:
        pytest.skip("build libyrec first (`make lib` in src/)")
    if not YREC.exists():
        pytest.skip("build src/yrec first")
    if not (REPO / "input").exists():
        pytest.skip("needs the full input/ tree (local only)")

    oracle = _run_cli(tmp_path / "cli")
    _compare("first", _run_lib(tmp_path / "lib1"), oracle)
    _compare("second", _run_lib(tmp_path / "lib2"), oracle)

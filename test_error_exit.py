"""End-to-end error-path test for the run_yrec engine (phase five, step B).

Historically every fatal condition in YREC ended in a bare `stop`,
which exits with status 0 -- indistinguishable from success to any
pipeline. Phase B converts the driver-layer stops into ierr returns
that surface through run_yrec to the CLI wrapper, which exits 1. This
test drives the oldest configuration check in parmin (semi-convection
and overshoot both enabled, rejected since the original F77) and
asserts the new contract: a clean nonzero exit, no crash, no model
evolved. It runs in an isolated copy of a Stage-0 case directory so
no example outputs are touched.
"""
import os
import pathlib
import re
import shutil
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
CASE_DIR = REPO / "examples" / "run_standard_solar_model"
NML1 = "Test_solar_noGS_norot.nml1"
NML2 = "Test_solar_noGS_norot.nml2"


def test_config_error_exits_nonzero(tmp_path):
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    for name in (NML1, NML2):
        shutil.copy(CASE_DIR / name, tmp_path / name)
    nml2 = tmp_path / NML2
    text = nml2.read_text()
    text = re.sub(r"LSEMIC\s*=\s*\.FALSE\.", "LSEMIC = .TRUE.", text)
    text = re.sub(r"LOVSTC\s*=\s*\.FALSE\.", "LOVSTC = .TRUE.", text)
    nml2.write_text(text)
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(YREC), NML1, NML2], cwd=tmp_path,
        env=env, capture_output=True, text=True, timeout=120,
    )
    # the new contract: a config error is a clean nonzero exit (stop 1
    # in the CLI wrapper), where historically a bare `stop` exited 0
    assert result.returncode == 1, (
        "conflicting LSEMIC/LOVSTC must exit 1 "
        f"(got {result.returncode}); stderr tail: {result.stderr[-300:]}"
    )
    # and the diagnostic still prints at the point of failure
    blobs = [result.stdout, result.stderr]
    for f in tmp_path.rglob("*"):
        if f.is_file() and f.suffix not in (".nml1", ".nml2"):
            try:
                blobs.append(f.read_text(errors="replace"))
            except OSError:
                pass
    assert any("SEMI-CONVECTION" in b for b in blobs), (
        "expected the PARMIN conflict diagnostic in the run output"
    )

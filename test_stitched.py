"""The stitched model: its index map, N^2 and seismic integrals on
synthetic arrays (core/diagnostics/stitched_kernels.f90), and the star's
stitched model of a starting model (core/diagnostics/stitched_model.f90).

Runs src/test_stitched (built by `make tests` in src/) on
src/core/test/test_structure_tams.inlist, the 1 Msun TAMS start model
with the base solar physics. Every check line must end in ok; the checks
are listed in the header of src/core/test/test_stitched.f90, and each
fails under a mutation of the code it covers
(src/core/test/test_stitched.mutations.toml). The binary exits 1 when a
check fails. The run prints the ERROR lines of the two error paths it
checks (a column without a value, an envelope integration that fails).
"""
import os
import pathlib
import subprocess

import pytest

import conftest  # noqa: F401  (sets YREC_INPUT and YREC_START)

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_stitched"
DECK = REPO / "src" / "core" / "test" / "test_structure_tams.inlist"
NUM_CHECKS = 24


def test_stitched(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_stitched not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_OUTPUT"] = "output"
    (tmp_path / "output").mkdir()
    result = subprocess.run([str(BINARY), str(DECK)], cwd=tmp_path, env=env,
                            capture_output=True, text=True)
    verdicts = [line.rstrip() for line in result.stdout.splitlines()
                if line.rstrip().endswith(("  ok", "  FAIL"))]
    failed = [line for line in verdicts if line.endswith("FAIL")]
    assert not failed, "\n".join(failed)
    assert len(verdicts) == NUM_CHECKS, "\n".join(verdicts)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    assert "test_stitched: done" in result.stdout

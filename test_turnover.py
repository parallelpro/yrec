"""The turnover walk on synthetic profiles: the base of the surface
convection zone and its turnover time (core/diagnostics/turnover_walk.f90).

Runs src/test_turnover (built by `make tests` in src/). Every check line
must end in ok; the checks are listed in the header of
src/core/test/test_turnover.f90, and each fails under a mutation of the
walk (src/core/test/test_turnover.mutations.toml). The binary exits 1
when a check fails.
"""
import pathlib
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_turnover"
NUM_CHECKS = 12


def test_turnover(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_turnover not built (run `make tests` in src/)")
    result = subprocess.run([str(BINARY)], cwd=tmp_path, capture_output=True, text=True)
    verdicts = [line.rstrip() for line in result.stdout.splitlines()
                if line.rstrip().endswith(("  ok", "  FAIL"))]
    failed = [line for line in verdicts if line.endswith("FAIL")]
    assert not failed, "\n".join(failed)
    assert len(verdicts) == NUM_CHECKS, "\n".join(verdicts)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    assert "test_turnover: done" in result.stdout

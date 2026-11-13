"""The timestep limits' kernels against their formulas.

Runs src/test_timestep_limits (built by `make tests` in src/) and
byte-compares its report with the pinned
src/core/test/expected_test_timestep_limits.out: on synthetic models of
five points (no deck, no tables), each limit of
src/core/step/timestep_limits.f90 (the structure change, the change of
omega, hydrogen and helium burning, the HR-diagram motion) against its
formula, with its switches, clipping and floors, and the tightest of
them within the growth limit. The binary exits nonzero when a check
fails (see the header of src/core/test/test_timestep_limits.f90). See
test_eos.py for the conventions.

The values printed are the checks' relative errors, round-off of order
1e-16 that another platform moves, or the limit itself (1e20 when
unlimited); in CI's tolerant mode (YREC_TEST_TOLERANT=1, conftest.py)
numbers below ABS_FLOOR compare equal.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_timestep_limits"
EXPECTED = REPO / "src" / "core" / "test" / "expected_test_timestep_limits.out"
ABS_FLOOR = 1.0e-12


def test_timestep_limits(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_timestep_limits not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True, timeout=60
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_timestep_limits",
                   abs_floor=ABS_FLOOR)

"""The Li/Be burn's kernels against their rules.

Runs src/test_li_be_burn (built by `make tests` in src/) and
byte-compares its report with the pinned
src/core/test/expected_test_li_be_burn.out. On a synthetic model of
twelve shells (no deck, no tables) it checks
src/core/step/li_be_burn.f90: the mixed region's averages; Li6, Li7 and
Be9 in radiative shells against the closed-form survival
X exp(-Lambda); the mixed region as one unit; the accreted share (the
zone gains x_acc M_acc); a base that moves during the step. The binary
exits nonzero when a check fails (see the header of
src/core/test/test_li_be_burn.f90). See test_eos.py for the conventions.

The values printed are the checks' relative errors, round-off of order
1e-16 that another platform moves (the accreted gas's Simpson rule:
1e-8), or the value checked; in CI's tolerant mode (YREC_TEST_TOLERANT=1,
conftest.py) numbers below ABS_FLOOR compare equal.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_li_be_burn"
EXPECTED = REPO / "src" / "core" / "test" / "expected_test_li_be_burn.out"
ABS_FLOOR = 1.0e-7


def test_li_be_burn(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_li_be_burn not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True, timeout=60
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_li_be_burn",
                   abs_floor=ABS_FLOOR)

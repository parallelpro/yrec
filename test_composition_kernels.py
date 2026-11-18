"""The composition step's kernels against their rules.

Runs src/test_composition_kernels (built by `make tests` in src/) and
byte-compares its report with the pinned
src/core/test/expected_test_composition_kernels.out. On synthetic models
of a few shells (no deck, no tables) it checks
src/core/step/composition_kernels.f90: the deep-mixing flags, the
mass-weighted mean of the mixed regions (each species keeps its mass)
and the clip and closure of the mass fractions. The binary exits nonzero
when a check fails (see the header of
src/core/test/test_composition_kernels.f90). See test_eos.py for the
conventions.

The values printed are the checks' errors, round-off of order 1e-16 that
another platform moves, or the value checked; in CI's tolerant mode
(YREC_TEST_TOLERANT=1, conftest.py) numbers below ABS_FLOOR compare
equal.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_composition_kernels"
EXPECTED = REPO / "src" / "core" / "test" / "expected_test_composition_kernels.out"
ABS_FLOOR = 1.0e-12


def test_composition_kernels(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_composition_kernels not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True, timeout=60
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_composition_kernels",
                   abs_floor=ABS_FLOOR)

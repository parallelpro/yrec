"""The envelope triangle's kernels: geometry and fit.

Runs src/test_triangle_fit (built by `make tests` in src/) and
byte-compares its report with the pinned
src/core/test/expected_test_triangle_fit.out: the planes of
fit_triangle_coefficients against linear functions (forward and reverse
relations), its non-finite coefficients for degenerate fitting points,
the point-in-triangle test, the reflection of a vertex, and walks of
reflect_until_inside. The binary exits nonzero when a check fails (see
the header of src/core/test/test_triangle_fit.f90); the timeout guards
against a walk that does not end. See test_eos.py for the conventions.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_triangle_fit"
EXPECTED = REPO / "src" / "core" / "test" / "expected_test_triangle_fit.out"


def test_triangle_fit(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_triangle_fit not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True, timeout=60
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_triangle_fit")

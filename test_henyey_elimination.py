"""The linear algebra of the Henyey solve against a dense solve.

Runs src/test_henyey_elimination (built by `make tests` in src/) and
byte-compares its report with the pinned
src/core/test/expected_test_henyey_elimination.out: henyey_eliminate
and henyey_solve on systems of 2, 3 and 40 shells with pseudo-random
terms and partials give the corrections of the Newton system that
src/core/structure/henyey_elimination.f90's header writes out (solved
here densely by LU decomposition) to a relative 1e-10, and the residual
maxima exactly; a zero pivot in any of the elimination's three divisions
or the surface system's two returns numerics_singular. The binary exits
nonzero when a check fails (see the header of
src/core/test/test_henyey_elimination.f90). See test_eos.py for the
conventions.

The values printed are the checks' errors, round-off of order 1e-16 to
1e-12 that another platform moves; in CI's tolerant mode
(YREC_TEST_TOLERANT=1, conftest.py) numbers below ABS_FLOOR compare
equal. A check's verdict (ok or FAIL) is always compared.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_henyey_elimination"
EXPECTED = REPO / "src" / "core" / "test" / "expected_test_henyey_elimination.out"
ABS_FLOOR = 1.0e-9


def test_henyey_elimination(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_henyey_elimination not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True, timeout=60
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_henyey_elimination",
                   abs_floor=ABS_FLOOR)

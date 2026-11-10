"""Standalone rotation regression test.

Runs src/test_rotation (built by `make tests` in src/) and byte-compares
its stdout with the pinned src/rotation/test/expected_test_rotation.out.
The secular-transport cases also check that the angular momentum solve
conserves the region's J without a wind (|sum dJ| < 1e-12 J): the sum
is round-off, so the pin keeps the check's ok, not its digits, and the
binary exits nonzero when the check fails. See test_eos.py for the
conventions.

Their changes of j/m, X and He4 (the rows after each
"zone  d(j/m)/(j/m) dX dHe4") are small changes of order-one
quantities, X and He4 after the mixing minus before, and j/m after the
map back relative to before, so their low digits are those quantities'
round-off: in CI's tolerant mode (YREC_TEST_TOLERANT=1, conftest.py)
the three columns are compared to an absolute 1e-11 (CHANGE_ABS_TOL).
Measured over 22 runs, the most any of them moved is 1.8e-13 for X and
He4 and 2.4e-13 for j/m. The runs: a correctly-rounded libm (the crmath
build), and omega, the diffusion coefficients or both perturbed by up to
four ulp per zone. That leaves margins of 57 and 42. A defect stays far
above the bound: the composition diffusion coefficients off by 1e-3
move X by up to 6e-6 (by 6e-9 for 1e-6), and the angular momentum's off
by 1e-3 move j/m by up to 9e-5. A wrong or weakened boundary flux, or a
dropped term, fails too. The local byte-exact comparison is unchanged.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_rotation"
EXPECTED = REPO / "src" / "rotation" / "test" / "expected_test_rotation.out"
CHANGE_HEADER = "zone  d(j/m)/(j/m) dX dHe4"
CHANGE_ABS_TOL = 1.0e-11


def change_tolerances(expected):
    """CHANGE_ABS_TOL for the d(j/m)/(j/m), dX and dHe4 tokens (1-3) of
    each row (zone and three numbers) after a CHANGE_HEADER line, keyed
    by (line number from 1, token)."""
    lines = expected.splitlines()
    tol = {}
    for i, line in enumerate(lines):
        if line != CHANGE_HEADER:
            continue
        for k in range(i + 1, len(lines)):
            row = lines[k].split()
            if len(row) != 4 or not row[0].isdigit():
                break
            for tok in (1, 2, 3):
                tol[(k + 1, tok)] = CHANGE_ABS_TOL
    return tol


def test_rotation_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_rotation not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    # a few printed values are round-off around zero (|x| < 1e-12), which
    # another platform's libm may move: they count as equal; the changes
    # of j/m, X and He4 are compared to CHANGE_ABS_TOL (the docstring)
    expected = EXPECTED.read_text()
    tolerances = change_tolerances(expected)
    assert len(tolerances) == 4*8*3, "the pin's change rows were not found"
    compare_output(result.stdout, expected, "test_rotation", abs_floor=1.0e-12,
                   abs_tol=tolerances)

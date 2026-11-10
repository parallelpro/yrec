"""The rezone's regridding kernels on a synthetic model.

Runs src/test_regrid (built by `make tests` in src/) and byte-compares
its report with the pinned baseline src/core/test/expected_test_regrid.out:
the shell masses of a mesh sum to the total mass and follow the centre,
interior and outer rules; profiles linear in log m are reproduced at the
new points; and the changes of each species' mass, of X + Y + Z and of
the angular momentum and rotational kinetic energy across a regrid with
a composition step at a convective-core edge stay within 10% of their
pinned values. The binary exits nonzero when a check fails (see the
header of src/core/test/test_regrid.f90). See test_eos.py for the
conventions.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_regrid"


def test_regrid_conservation(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_regrid not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    report = result.stdout.split("test_regrid report\n", 1)[1]
    expected = REPO / "src" / "core" / "test" / "expected_test_regrid.out"
    # The checks of sums, rules and linear profiles print residuals that
    # are zero but for round-off (at most 1.4e-15, against bounds of n eps
    # to 1e-13), and their digits may move on another platform: below
    # 1e-12 two values count as equal (YREC_TEST_TOLERANT, conftest.py);
    # the binary's own bounds judge them. The relative changes of the
    # species' masses, J and E_rot (1e-5 and up) are compared as before.
    compare_output(report, expected.read_text(), "test_regrid", abs_floor=1.0e-12)

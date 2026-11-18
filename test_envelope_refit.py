"""The kernels that turn envelope points into interior shells.

Runs src/test_envelope_refit (built by `make tests` in src/), which needs
no deck and no tables, and byte-compares its report with the pinned
src/core/test/expected_test_envelope_refit.out: on a synthetic model
whose profiles are linear in log m, append_envelope_shells keeps the
points inside the fitting mass, puts the last zone at the fitting mass
and on the profiles in both of its branches, and returns ierr 1 for a
degenerate interval; shift_envelope_onto_shell joins the envelope to the
last shell; refit_shell_masses' masses sum to the star's; and
guess_new_shell_rotation follows its three laws. The binary exits
nonzero when a check fails (see the header of
src/core/test/test_envelope_refit.f90). See test_eos.py for the
conventions.

The values printed are the checks' residuals, zero but for round-off
(bounds of 1e-14 to 1e-11): in CI's tolerant mode (YREC_TEST_TOLERANT=1,
conftest.py) numbers below ABS_FLOOR compare equal. A check's verdict
(ok or FAIL) is always compared.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_envelope_refit"
EXPECTED = REPO / "src" / "core" / "test" / "expected_test_envelope_refit.out"
ABS_FLOOR = 1.0e-10


def test_envelope_refit(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_envelope_refit not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    report = result.stdout.split("test_envelope_refit report\n", 1)[1]
    compare_output(report, EXPECTED.read_text(), "test_envelope_refit",
                   abs_floor=ABS_FLOOR)

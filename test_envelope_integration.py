"""The envelope integration kernel and its integrands on a solar surface.

Runs src/test_envelope_integration (built by `make tests` in src/) with
the microphysics and integration controls of
src/core/test/test_structure_tams.inlist (no model is read) and
byte-compares its report with the pinned
src/core/test/expected_test_envelope_integration.out: the integrands
against their equations at a point (hydrostatic equilibrium in tau and
in P, mass continuity, the Eddington T(tau)); the integration's status,
events and end points (tau = 2/3 with T = Teff, the fitting mass); each
Bulirsch-Stoer step of the atmosphere and the envelope against a
fourth-order Runge-Kutta integration of the same equations; the vertex
as the interpolation to the fitting mass; the thin envelope's vertex at
the surface; the retried start; and the error contracts (maxstp in the
atmosphere and in the envelope, Bulirsch-Stoer limits out of range, a
step that underflows). The binary exits nonzero when a check fails (see
the header of src/core/test/test_envelope_integration.f90). See
test_eos.py for the conventions.

The values printed are the checks' errors, which another platform's
libm and tables move at round-off; in CI's tolerant mode
(YREC_TEST_TOLERANT=1, conftest.py) numbers below ABS_FLOOR compare
equal. A check's verdict (ok or FAIL) is always compared.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_envelope_integration"
DECK = REPO / "src" / "core" / "test" / "test_structure_tams.inlist"
EXPECTED = REPO / "src" / "core" / "test" / "expected_test_envelope_integration.out"
ABS_FLOOR = 1.0e-6


def test_envelope_integration(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_envelope_integration not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_OUTPUT"] = "output"
    (tmp_path / "output").mkdir()
    result = subprocess.run(
        [str(BINARY), str(DECK)], cwd=tmp_path, env=env,
        capture_output=True, text=True, timeout=600,
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    report = result.stdout.split("test_envelope_integration report\n", 1)[1]
    compare_output(report, EXPECTED.read_text(), "test_envelope_integration",
                   abs_floor=ABS_FLOOR)

"""The mass-change kernels against the conservation laws they keep.

Runs src/test_mass_change (built by `make tests` in src/) and
byte-compares its stdout with the pinned
src/core/test/expected_test_mass_change.out. On a synthetic 40-shell
model (no deck, no tables) it checks that after a loss or gain of mass
the shell masses sum to the new total and the convection zone's masses
scale about its base (with dm(1) from the center when the zone reaches
it), that a thermal wind leaves the zone J_cz + dJ for walpcz = 0, -1
and -2, and that accretion conserves every basic species. The binary
exits nonzero when a check fails (see the header of
src/core/test/test_mass_change.f90). See test_eos.py for the conventions.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_mass_change"
EXPECTED = REPO / "src" / "core" / "test" / "expected_test_mass_change.out"


def test_mass_change_conservation(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_mass_change not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    # Every check's value is the residual of an identity the kernels
    # keep (a sum, a scaling, a conservation law), zero but for
    # round-off (at most 2e-15), and its digits may move on another
    # platform: below 1e-12 two values count as equal (YREC_TEST_TOLERANT,
    # conftest.py); the binary's own bounds judge them.
    compare_output(result.stdout, EXPECTED.read_text(), "test_mass_change", abs_floor=1.0e-12)

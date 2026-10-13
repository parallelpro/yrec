"""Rotation's kernels against the laws they keep.

Runs src/test_rotation_kernels (built by `make tests` in src/) and
byte-compares its stdout with the pinned
src/rotation/test/expected_test_rotation_kernels.out. On a synthetic
60-zone model (no deck, no tables) it checks the spherical limits of
the moments of inertia and the shape factors and how fp and ft grow with
omega; that omega_from_j gives each radiative zone and the convective
envelope back their angular momentum, the envelope as a solid body; that
the angular momentum diffusion conserves J without a wind, loses exactly
the wind's surface term with one, and keeps a solid body; that the
composition mixing conserves each species; and the Kawaler wind law's
exponents below and above saturation. The binary exits nonzero when a
check fails (see the header of src/rotation/test/test_rotation_kernels.f90).
See test_eos.py for the conventions.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_rotation_kernels"
EXPECTED = REPO / "src" / "rotation" / "test" / "expected_test_rotation_kernels.out"


def test_rotation_kernels(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_rotation_kernels not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_rotation_kernels")

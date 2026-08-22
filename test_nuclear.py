"""Standalone nuclear-domain regression test (2026).

Runs src/test_nuclear (built by `make tests` in src/) with YREC_INPUT
pointing at the repository's input/ tables, and compares its stdout
against the pinned src/nuclear/test/expected_test_nuclear.out
(byte-exact locally; ULP-tolerant under YREC_TEST_TOLERANT=1 in CI,
via conftest.compare_output).

Exercises remap's cross-section scales (Solar Fusion II over SFI
references), the 13 reaction rates over a (logT, logRho) grid, the
Itoh et al. (1996) neutrino losses, azbar, and deutrate -- the
composition-in/rates-out core of the domain. engeb (the per-zone
energy-generation driver) is model-coupled and covered by Stage-0.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_nuclear"
EXPECTED = REPO / "src" / "nuclear" / "test" / "expected_test_nuclear.out"


def test_nuclear_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_nuclear not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_nuclear")

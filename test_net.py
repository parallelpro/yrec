"""net/'s unit test (2026).

Runs src/test_net (built by `make tests` in src/) with YREC_INPUT
pointing at the repository's input/ tables, and compares its stdout
against the pinned src/net/test/expected_test_net.out (byte-exact
locally; ULP-tolerant under YREC_TEST_TOLERANT=1 in CI, via
conftest.compare_output). The program also exits 1 when one of its
ok/FAIL checks fails.

test_net checks net/ piece by piece (its header lists the checks): the
settings and their install, the reaction table's self-check, the rates
and their derivatives, the implicit burn and its error codes, the
energy generation with the thermal neutrino losses and deuterium, the
screening, the light elements and the trace-element integrator, the
solar neutrino fluxes, and the legacy cross-section map.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_net"
EXPECTED = REPO / "src" / "net" / "test" / "expected_test_net.out"


def test_net_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_net not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_net")

"""kap's test program: pins and ok/FAIL checks.

Runs src/test_kap (built by `make tests` in src/) against the
repository's opacity tables and byte-compares its stdout with the pinned
src/kap/test/expected_test_kap.out: the opacity on a grid of points for
every family, two error paths, and the ok/FAIL checks (exactness at the
table nodes, derivatives against differences, continuity at the ramp's
ends and the tables' edges, the blends' and the conductive mixture's
formulas, the error contracts), each confirmed to fail under a mutation
in src/kap/test/test_kap.mutations.toml (python3 src/tools/mutate.py
kap). The binary exits nonzero when a check fails. See test_eos.py for
the conventions.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_kap"
EXPECTED = REPO / "src" / "kap" / "test" / "expected_test_kap.out"


def test_kap_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_kap not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_kap")

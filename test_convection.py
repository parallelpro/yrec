"""Standalone convection regression test.

Runs src/test_convection (built by `make tests` in src/) and byte-compares
its stdout with the pinned src/mixing/test/expected_test_convection.out.
See test_eos.py for the conventions.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_convection"
EXPECTED = REPO / "src" / "mixing" / "test" / "expected_test_convection.out"


def test_convection_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_convection not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_convection")

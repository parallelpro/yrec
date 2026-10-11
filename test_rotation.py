"""Standalone rotation regression test.

Runs src/test_rotation (built by `make tests` in src/) and byte-compares
its stdout with the pinned src/rotation/test/expected_test_rotation.out.
See test_eos.py for the conventions.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_rotation"
EXPECTED = REPO / "src" / "rotation" / "test" / "expected_test_rotation.out"


def test_rotation_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_rotation not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_rotation")

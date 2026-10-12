"""Standalone element-diffusion regression test.

Runs src/test_element_diffusion (built by `make tests` in src/) and byte-compares
its stdout with the pinned src/element_diffusion/test/expected_test_element_diffusion.out.
See test_eos.py for the conventions.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_element_diffusion"
EXPECTED = REPO / "src" / "element_diffusion" / "test" / "expected_test_element_diffusion.out"


def test_element_diffusion_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_element_diffusion not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_element_diffusion")

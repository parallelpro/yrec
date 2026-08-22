"""Standalone atm-domain regression test (phase three, ROADMAP.md stage 2).

Runs src/test_atm (built by `make tests` in src/) against the
repository's Kurucz, Castelli/Kurucz, and Allard atmosphere tables
and byte-compares its stdout with the pinned
src/atm/test/expected_test_atm.out. See test_eos.py for the
conventions.
"""
import os
import pathlib
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_atm"
EXPECTED = REPO / "src" / "atm" / "test" / "expected_test_atm.out"


def test_atm_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_atm not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert result.stdout == EXPECTED.read_text(), (
        "test_atm output differs from src/atm/test/expected_test_atm.out; "
        "if the change is deliberate, regenerate and commit the baseline."
    )

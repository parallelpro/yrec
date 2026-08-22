"""Standalone kap-domain regression test (phase three, ROADMAP.md stage 2).

Runs src/test_kap (built by `make tests` in src/) against the
repository's OPAL95/GS98.OP17 opacity table and byte-compares its
stdout with the pinned src/kap/test/expected_test_kap.out. See
test_eos.py for the conventions.
"""
import os
import pathlib
import subprocess

import pytest

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
    assert result.stdout == EXPECTED.read_text(), (
        "test_kap output differs from src/kap/test/expected_test_kap.out; "
        "if the change is deliberate, regenerate and commit the baseline."
    )

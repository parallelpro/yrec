"""Standalone convection test.

Runs src/test_convection (built by `make tests` in src/), which exits 1
when one of its ok/FAIL checks fails, and compares its stdout with the
pinned src/convection/test/expected_test_convection.out: byte for byte,
except that a check line's measured error is masked (only its label and
its verdict are compared), since round-off differs between platforms.
See test_eos.py for the conventions.
"""
import os
import pathlib
import re
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_convection"
EXPECTED = REPO / "src" / "convection" / "test" / "expected_test_convection.out"


def test_convection_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_convection not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    check = re.compile(r"^(.*\S)\s+[-+0-9.Ee]+\s+(ok|FAIL)$")

    def masked(text):
        out = []
        for line in text.splitlines():
            m = check.match(line)
            out.append(f"{m.group(1)} [{m.group(2)}]" if m else line)
        return "\n".join(out) + "\n"
    compare_output(masked(result.stdout), masked(EXPECTED.read_text()), "test_convection")

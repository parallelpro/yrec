"""Standalone eos-domain regression test (phase three, ROADMAP.md stage 2).

Runs src/test_eos (built by `make tests` in src/) with YREC_INPUT
pointing at the repository's input/ tables, and byte-compares its
stdout against the pinned src/eos/test/expected_test_eos.out.

Unlike the Stage-0 full-model regression, this exercises the eos
facade entries (eos_get, eos_get_gamma1) directly over a fixed grid,
without needing a stellar model to converge -- including branch
combinations no Stage-0 case reaches (e.g. eos_get_gamma1's Yale/SCV
branch). The MHD section self-reports SKIPPED unless YREC_MHD_TABLES
supplies tables (none ship with YREC); when tables are provided the
output will differ from the pinned baseline, which only covers the
shipped-data configuration.

To regenerate the baseline after a deliberate eos change: rebuild,
run src/test_eos with YREC_INPUT set, review the diff like any other,
and commit the new expected file.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_eos"
EXPECTED = REPO / "src" / "eos" / "test" / "expected_test_eos.out"


def test_eos_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_eos not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env.pop("YREC_MHD_TABLES", None)  # baseline covers shipped data only
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_eos")

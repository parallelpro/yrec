"""Every generated source file matches its generator.

Runs each generator in src/tools/ with --check: it writes nothing and fails
if a committed output differs from what the generator would write, so a hand
edit to a generated file, or a source change without regenerating, is caught.
"""
import pathlib
import subprocess
import sys

import pytest

SRC = pathlib.Path(__file__).parent / "src"
GENERATORS = ["gen_controls_state.py", "gen_inlist.py", "gen_controls_dump.py",
              "gen_makefile_deps.py", "gen_default_columns.py"]


@pytest.mark.parametrize("generator", GENERATORS)
def test_generated_files_current(generator):
    result = subprocess.run([sys.executable, str(SRC / "tools" / generator), "--check"],
                            cwd=SRC, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr

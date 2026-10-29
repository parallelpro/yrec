"""Every generated source file matches its generator, and the hand-kept
list of path controls covers the registry.

Runs each generator in src/tools/ with --check: it writes nothing and fails
if a committed output differs from what the generator would write, so a hand
edit to a generated file, or a source change without regenerating, is caught.
"""
import pathlib
import re
import subprocess
import sys

import pytest

SRC = pathlib.Path(__file__).parent / "src"
GENERATORS = ["gen_controls_state.py", "gen_inlist.py", "gen_controls_dump.py",
              "gen_makefile_deps.py", "gen_default_columns.py", "gen_star_docs.py"]


@pytest.mark.parametrize("generator", GENERATORS)
def test_generated_files_current(generator):
    result = subprocess.run([sys.executable, str(SRC / "tools" / generator), "--check"],
                            cwd=SRC, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


# String controls that are not file paths (names, codes, the output
# file's own name and prefix).
NOT_PATHS = {"amix", "aiso", "anewcp", "atmp", "awind", "descrip", "pulse_format",
             "profile_data_prefix", "star_history_name",
             # read by the output from star%ctrl (see expand_paths)
             "history_columns_file", "profile_columns_file"}


def test_every_path_control_is_expanded():
    """run_files' expand_paths expands the {YREC_...} placeholders in
    every string control of the registry that names a file."""
    registry = (SRC / "defaults" / "controls_registry.tsv").read_text().splitlines()
    paths = {row.split("\t")[1].lower() for row in registry[1:]
             if row.split("\t")[2].startswith("string")} - NOT_PATHS
    code = (SRC / "io" / "run_files.f90").read_text()
    body = re.search(r"subroutine expand_paths(.*?)end subroutine", code, re.S).group(1)
    expanded = {n.lower() for n in re.findall(r"expand_value\((\w+)", body)}
    assert paths <= expanded, sorted(paths - expanded)

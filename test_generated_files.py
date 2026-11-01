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


def test_every_path_a_deck_sets_is_a_path_control():
    """A shipped deck that sets a string control to a {YREC_...} path sets a
    control of type path, which the read expands (io/expand_paths.inc is
    generated from the type path rows, so a path control needs no hand
    edit): the deck's name for it (inlist name, alias or legacy name) is a
    path row's, or a control nothing uses (surface_composition_file)."""
    registry = [row.split("\t") for row in
                (SRC / "defaults" / "controls_registry.tsv").read_text().splitlines()]
    cols = registry[0]
    rows = [dict(zip(cols, r)) for r in registry[1:]]
    kind = {}
    for r in rows:
        for name in [r["proposed_name"], r["legacy"]] + r["aliases"].split(","):
            if name and name != "-":
                kind[name.lower()] = "path" if r["used"] == "none" else r["type"]
    repo = SRC.parent
    decks = [p for top in ("examples", "testsuite") for pat in ("*.inlist", "*.nml1", "*.nml2")
             for p in sorted((repo / top).rglob(pat))]
    entry = re.compile(r"(?im)^\s*&?\s*(\w+)\s*(?:\([^)]*\))?\s*=\s*['\"]\{YREC_")
    wrong = set()
    for d in decks:
        for name in entry.findall(d.read_text(errors="replace")):
            if not kind.get(name.lower(), "path").startswith("path"):
                wrong.add((name, d.name))
    assert decks and not wrong, sorted(wrong)[:10]

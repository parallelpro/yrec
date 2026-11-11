#!/usr/bin/env python3
"""Mutation check of test_numerics: `python3 src/tools/mutate.py numerics`,
whose manifest is numerics/test/test_numerics.mutations.toml (each
mutation's file, text, replacement and the check that must print FAIL).
This name stays for the commands that use it:

    python3 src/numerics/test/mutate_numerics.py [--only REGEX] [--check]

from the repository root, after `make tests` in src/.
"""
import pathlib
import subprocess
import sys

TOOL = pathlib.Path(__file__).resolve().parents[2] / "tools" / "mutate.py"
sys.exit(subprocess.call([sys.executable, str(TOOL), "numerics"] + sys.argv[1:]))

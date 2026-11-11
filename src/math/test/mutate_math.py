#!/usr/bin/env python3
"""Mutation check of test_math: `python3 src/tools/mutate.py math`, whose
manifest is math/test/test_math.mutations.toml (each mutation of
math_lib.f90 with the backend it is built under, and the check that must
print FAIL). This name stays for the commands that use it:

    python3 src/math/test/mutate_math.py [--only REGEX] [--check]

from the repository root, after `make tests` in src/.
"""
import pathlib
import subprocess
import sys

TOOL = pathlib.Path(__file__).resolve().parents[2] / "tools" / "mutate.py"
sys.exit(subprocess.call([sys.executable, str(TOOL), "math"] + sys.argv[1:]))

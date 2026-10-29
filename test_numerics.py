"""numerics_lib: pins and checks.

Runs src/test_numerics (built by `make tests` in src/) and compares its
stdout with the pinned src/numerics/test/expected_test_numerics.out.
The pins (every public routine on small polynomial and rational data,
error paths included) are compared byte for byte; they guard
refactors. The checks test behaviour: the splines' exactness on cubics
and lines and their orders, osplin's monotonicity, convexity and
consistency, the Lagrange weights' exactness, lir's and the searches'
contracts, the linear solvers against residuals and dense solves with
their singular cases, bsstep against sin and cos with its error codes,
and the status codes (see the header of
src/numerics/test/test_numerics.f90 and numerics/README.md section 6).
The binary exits nonzero when a check fails; a check's line is compared
by its label and status, since its error is round-off that differs
between platforms. src/numerics/test/mutate_numerics.py confirms that
each check fails under a mutation of the code it covers.
See test_eos.py for the conventions.
"""
import os
import pathlib
import re
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_numerics"
EXPECTED = REPO / "src" / "numerics" / "test" / "expected_test_numerics.out"


def test_numerics_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_numerics not built (run `make tests` in src/)")
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
    compare_output(masked(result.stdout), masked(EXPECTED.read_text()), "test_numerics")


def test_numerics_mutations_apply():
    """Every mutation of mutate_numerics.py still finds its text (a
    rewrite of a mutated routine must update its mutations)."""
    script = REPO / "src" / "numerics" / "test" / "mutate_numerics.py"
    result = subprocess.run(["python3", str(script), "--check"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr

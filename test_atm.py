"""atm/: pins and checks.

Runs src/test_atm (built by `make tests` in src/) against the
repository's Kurucz, Castelli/Kurucz and Allard atmosphere tables and
compares its stdout with the pinned src/atm/test/expected_test_atm.out.
The pins (the lookups of every table through atm_lib, the three T-tau
relations on a grid of optical depths, an error path) are compared
byte for byte; they guard refactors. The checks test behaviour: the
T-tau relations' photospheres and limits, the lookups' exactness at
the nodes and on a plane or a polynomial, their continuity, their
refusals and error contracts, and the choice descriptor (see the
header of src/atm/test/test_atm.f90 and atm/README.md section 4). The
binary exits nonzero when a check fails; a check's line is compared by
its label and status, since its error is round-off that differs
between platforms. `python3 src/tools/mutate.py atm` confirms that each
check fails under a mutation of the code it covers. See test_eos.py for
the conventions.
"""
import os
import pathlib
import re
import subprocess

import pytest

from conftest import compare_output

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
    check = re.compile(r"^(.*\S)\s+[-+0-9.Ee]+\s+(ok|FAIL)$")
    def masked(text):
        out = []
        for line in text.splitlines():
            m = check.match(line)
            out.append(f"{m.group(1)} [{m.group(2)}]" if m else line)
        return "\n".join(out) + "\n"
    compare_output(masked(result.stdout), masked(EXPECTED.read_text()), "test_atm")


def test_atm_mutations_apply():
    """Every mutation of test_atm.mutations.toml still finds its text (a
    rewrite of a mutated routine must update its mutations)."""
    result = subprocess.run(
        ["python3", str(REPO / "src" / "tools" / "mutate.py"), "atm", "--check"],
        capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr

"""math_lib's backend contract.

Runs src/test_math (math_lib under the build's backend) and, when it
was built (the MESA SDK present), src/test_math_crmath (under the
correctly-rounded crmath backend); both are built by `make tests` in
src/. Each prints its checks with ok or FAIL and exits nonzero on a
FAIL: the integer paths, the real paths within the composition's error,
a single-precision exponent, negative bases, pow(0, y) under each
backend, and an exponent beyond the integer range (see the header of
src/math/test/test_math.f90). The output, which holds no numbers, is
compared with the pinned list of checks for the backend.
See test_eos.py for the conventions.
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
TESTS = REPO / "src" / "math" / "test"


@pytest.mark.parametrize("binary", ["test_math", "test_math_crmath"])
def test_math(binary, tmp_path):
    path = REPO / "src" / binary
    if not path.exists():
        pytest.skip(f"src/{binary} not built (run `make tests` in src/)")
    result = subprocess.run([str(path)], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    backend = result.stdout.splitlines()[0].split()[-1]
    expected = TESTS / f"expected_test_math_{backend}.out"
    compare_output(result.stdout, expected.read_text(), binary)

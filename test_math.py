"""math_lib's backend contract.

Runs src/test_math (math_lib under the build's backend) and, when it
was built (the MESA SDK present), src/test_math_crmath (under the
correctly-rounded crmath backend); both are built by `make tests` in
src/. Each prints its checks with ok or FAIL and exits nonzero on a
FAIL: the integer paths, the real paths within the composition's error,
a single-precision exponent, negative bases, pow(0, y) under each
backend, and an exponent beyond the integer range (see the header of
src/math/test/test_math.f90). The output, which holds no numbers, is
compared with the pinned list of checks for the backend. Without the
MESA SDK the crmath leg is skipped with a reason that says so (CI has
no SDK: there only the intrinsic backend is tested).
src/math/test/mutate_math.py confirms that each check fails under a
mutation of the math_lib code it covers.
See test_eos.py for the conventions.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
TESTS = REPO / "src" / "math" / "test"


def mesa_sdk():
    """The MESA SDK as the Makefile finds it (MESASDK_DIR, else
    MESASDK_ROOT, else /Applications/mesasdk), or None without crmath."""
    sdk = pathlib.Path(os.environ.get("MESASDK_DIR") or os.environ.get("MESASDK_ROOT")
                       or "/Applications/mesasdk")
    return sdk if (sdk / "include" / "crmath.mod").exists() else None


@pytest.mark.parametrize("binary", ["test_math", "test_math_crmath"])
def test_math(binary, tmp_path):
    path = REPO / "src" / binary
    if binary == "test_math_crmath" and mesa_sdk() is None:
        pytest.skip("crmath backend NOT TESTED: no MESA SDK (include/crmath.mod under "
                    "MESASDK_DIR, MESASDK_ROOT or /Applications/mesasdk), so make tests "
                    "does not build test_math_crmath")
    if not path.exists():
        pytest.skip(f"src/{binary} not built (run `make tests` in src/)")
    result = subprocess.run([str(path)], cwd=tmp_path, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr
    backend = result.stdout.splitlines()[0].split()[-1]
    expected = TESTS / f"expected_test_math_{backend}.out"
    compare_output(result.stdout, expected.read_text(), binary)


def test_math_mutations_apply():
    """Every mutation of mutate_math.py still finds its text in
    math_lib.f90."""
    script = TESTS / "mutate_math.py"
    result = subprocess.run(["python3", str(script), "--check"], capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr

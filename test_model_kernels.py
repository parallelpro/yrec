"""model/'s star-free kernels: find_shell_boundaries on synthetic profiles.

Runs src/test_model_kernels (built by `make tests` in src/). Every check
line must end in ok; the checks are listed in the header of
src/core/test/test_model_kernels.f90, and each fails under a mutation of
the kernel (src/core/test/test_model_kernels.mutations.toml). The binary
exits 1 when a check fails.
"""
import pathlib
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_model_kernels"
NUM_CHECKS = 9


def test_model_kernels(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_model_kernels not built (run `make tests` in src/)")
    result = subprocess.run([str(BINARY)], cwd=tmp_path, capture_output=True, text=True)
    verdicts = [line.rstrip() for line in result.stdout.splitlines()
                if line.rstrip().endswith(("  ok", "  FAIL"))]
    failed = [line for line in verdicts if line.endswith("FAIL")]
    assert not failed, "\n".join(failed)
    assert len(verdicts) == NUM_CHECKS, "\n".join(verdicts)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    assert "test_model_kernels: done" in result.stdout

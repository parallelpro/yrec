"""run/'s contracts: the star calibration's last card and verdicts, and the
C entry's path limit and the units it closes.

Runs src/test_run_contracts (built by `make tests` in src/) on the
configuration matrix's cm_env_overshoot deck, in a scratch directory with
an empty output/. Every check line must end in ok; the checks and what
each covers are in the header of src/core/test/test_run_contracts.f90,
and each fails under a mutation of the code it covers
(src/core/test/test_run_contracts.mutations.toml). The binary exits 1
when a check fails. Only the verdicts are compared: the deck's own output
is the deck pins' business.
"""
import os
import pathlib
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_run_contracts"
DECK = REPO / "examples" / "run_config_matrix" / "cm_env_overshoot.inlist"
NUM_CHECKS = 17


def test_run_contracts(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_run_contracts not built (run `make tests` in src/)")
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    env.pop("YREC_OUTPUT", None)
    result = subprocess.run(
        [str(BINARY), str(DECK.resolve())], cwd=tmp_path, env=env,
        capture_output=True, text=True, timeout=600,
    )
    verdicts = [line.rstrip() for line in result.stdout.splitlines()
                if line.rstrip().endswith(("  ok", "  FAIL"))]
    failed = [line for line in verdicts if line.endswith("FAIL")]
    assert not failed, "\n".join(failed)
    assert len(verdicts) == NUM_CHECKS, "\n".join(verdicts)
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr[-3000:]
    assert "test_run_contracts: done" in result.stdout

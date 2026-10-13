"""The structure equations' partials against finite differences.

Runs src/test_structure (built by `make tests` in src/) on two decks and
byte-compares its report with the pinned baselines in src/core/test/:
  - src/core/test/test_structure_tams.inlist: the 1 Msun TAMS start model
    with the base solar physics (OPAL 2006 + SCV EOS, OPAL95 opacities),
    tolerance 3e-3 (the tables' own derivatives);
  - testsuite/Test_solarGS98_yaleeos.inlist: the 1 Msun birthline model
    with the Yale EOS, tolerance 1e-4, which also checks the eps_grav
    partials.
The binary exits nonzero when an entry outside the documented
approximations fails (see the header of src/core/test/test_structure.f90).
See test_eos.py for the conventions.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_structure"
CASES = [
    ("tams", REPO / "src" / "core" / "test" / "test_structure_tams.inlist", "3e-3"),
    ("yale", REPO / "testsuite" / "Test_solarGS98_yaleeos.inlist", "1e-4"),
]


@pytest.mark.parametrize("label,deck,tolerance", CASES, ids=[c[0] for c in CASES])
def test_structure_jacobian(label, deck, tolerance, tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_structure not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_OUTPUT"] = "output"
    (tmp_path / "output").mkdir()
    result = subprocess.run(
        [str(BINARY), str(deck), tolerance], cwd=tmp_path, env=env,
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    report = result.stdout.split("test_structure report\n", 1)[1]
    expected = REPO / "src" / "core" / "test" / f"expected_test_structure_{label}.out"
    compare_output(report, expected.read_text(), f"test_structure[{label}]")

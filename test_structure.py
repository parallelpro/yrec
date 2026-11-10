"""The structure equations' partials against finite differences.

Runs src/test_structure (built by `make tests` in src/) on two decks and
byte-compares its report with the pinned baselines in src/core/test/:
  - src/core/test/test_structure_tams.inlist: the 1 Msun TAMS start model
    with the base solar physics (OPAL 2006 + SCV EOS, OPAL95 opacities),
    tolerance 3e-3 (the tables' own derivatives);
  - testsuite/Test_solarGS98_yaleeos.inlist: the 1 Msun birthline model
    with the Yale EOS, tolerance 1e-4, which also checks the eps_grav
    partials;
  - extra: the TAMS deck again with a test energy source registered
    (src/core/structure/energy_sources.f90), which checks its partials
    and the size of the term it adds.
The binary exits nonzero when an entry outside the documented
approximations fails (see the header of src/core/test/test_structure.f90).
See test_eos.py for the conventions.

The rel.diff column (and the extra case's eps_extra rel.diff) is the
finite difference's relative error, printed to two digits. Where the
partial is analytic it is round-off amplified by 1/h: the pins' values
of that kind are at most 4.1e-9, and other platforms move them. So in
CI's tolerant mode (YREC_TEST_TOLERANT=1, conftest.py) that column is
compared to an absolute 1e-8 (REL_DIFF_ABS_TOL). Measured, the most one
moved is 7.4e-10 (8.0e-10 to 6.4e-11; margin 13), under a
correctly-rounded libm (the crmath build) and with each shell's state
perturbed by up to two ulp (6 runs). The bound lies below every value
the tables or the truncation set (3.7e-7 and up, which do not move). A
partial off by 1e-6 moves its rel.diff by 1e-6, and the comparison
fails. The local byte-exact comparison is unchanged.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_structure"
REL_DIFF_ABS_TOL = 1.0e-8


def rel_diff_tolerances(expected):
    """REL_DIFF_ABS_TOL for the rel.diff token of each partial's row
    (name, analytic, f.d., rel.diff, status) and of each eps_extra line,
    keyed by (line number from 1, token)."""
    tol = {}
    for n, line in enumerate(expected.splitlines(), 1):
        t = line.split()
        if len(t) == 5 and t[4] in ("ok", "known", "FAIL"):
            tol[(n, 3)] = REL_DIFF_ABS_TOL
        elif len(t) == 7 and t[0] == "eps_extra" and t[5] == "rel.diff":
            tol[(n, 6)] = REL_DIFF_ABS_TOL
    return tol


CASES = [
    ("tams", REPO / "src" / "core" / "test" / "test_structure_tams.inlist", "3e-3"),
    ("yale", REPO / "testsuite" / "Test_solarGS98_yaleeos.inlist", "1e-4"),
    ("extra", REPO / "src" / "core" / "test" / "test_structure_tams.inlist", "3e-3"),
]


@pytest.mark.parametrize("label,deck,tolerance", CASES, ids=[c[0] for c in CASES])
def test_structure_jacobian(label, deck, tolerance, tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_structure not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_OUTPUT"] = "output"
    (tmp_path / "output").mkdir()
    result = subprocess.run(
        [str(BINARY), str(deck), tolerance] + (["extra"] if label == "extra" else []),
        cwd=tmp_path, env=env,
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stdout[-3000:] + result.stderr
    report = result.stdout.split("test_structure report\n", 1)[1]
    expected = REPO / "src" / "core" / "test" / f"expected_test_structure_{label}.out"
    expected_text = expected.read_text()
    tolerances = rel_diff_tolerances(expected_text)
    assert tolerances, "the pin's rel.diff tokens were not found"
    compare_output(report, expected_text, f"test_structure[{label}]", abs_tol=tolerances)

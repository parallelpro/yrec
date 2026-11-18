"""Shared helpers for the per-module standalone tests.

compare_output: byte-exact by default (the local discipline). When
YREC_TEST_TOLERANT=1 (set only by CI), numeric tokens may differ by a
tiny relative tolerance -- the standalone baselines are pinned on
macOS/arm64, and x86_64 gfortran rounds the last printed digit
differently in a handful of values (FMA/codegen ULP differences, seen
on the first CI run). Everything non-numeric, and the line structure,
must still match exactly. The tolerance is 1e-9: the eos grid prints
second-derivative quantities whose differencing amplifies platform ULP
noise to ~1e-11 relative (worst observed across CI runs 4-5: 9.6e-12),
so 1e-9 keeps two orders of margin over the noise while remaining far
below anything a genuine physics or refactor defect would produce.
Two values both below abs_floor count as equal (round-off around zero),
and a test may give single tokens an absolute tolerance (abs_tol) where
a printed number is a small change of an order-one quantity, so that
its low digits are that quantity's round-off (test_rotation.py's
composition and j/m changes).
"""
import os
import pathlib
import sys

# The table and starting-model roots of every program a test starts.
# Decks name their files as {YREC_INPUT}/... and {YREC_START}/...; the
# in-code fallback ('../../input') assumes the examples/<case>/ depth and
# silently breaks for testsuite/ (one level deep), so the harness sets
# the two variables once, absolutely. They are this repository's input/
# and startmodels/, whatever the shell exports: a test run tests its own
# tree. YREC_TEST_EXTERNAL_ROOTS=1 is the one opt-in that takes them
# from the environment (src/tools/harness.py has the rule).
# A test file names the roots as INPUT_ROOT and START_ROOT (imported
# from here), never as its own tree's directories, so that the opt-in
# and pytest's header hold for every test.
_REPO = pathlib.Path(__file__).parent
sys.path.insert(0, str(_REPO / "src" / "tools"))
import harness  # noqa: E402

_ROOT_LINES = harness.describe(_REPO)   # before the two are set below
_ROOTS = harness.apply(_REPO, os.environ)
INPUT_ROOT = _ROOTS["YREC_INPUT"]
START_ROOT = _ROOTS["YREC_START"]

REL_TOL = 1.0e-9


def pytest_report_header(config):
    """What this run tests: the binaries' directory and the two roots."""
    return [f"yrec binaries: {_REPO / 'src'}"] + _ROOT_LINES


def pytest_configure(config):
    """The `local` marker: the physics validation tier's tests that CI does
    not run (docs/PHYSICS_VALIDATION.md); `pytest -m "not local"` leaves
    them out."""
    config.addinivalue_line(
        "markers", "local: a validation test not run in CI (it needs GYRE, the full "
        "input tree, or minutes)")


def _to_float(tok):
    try:
        return float(tok)
    except ValueError:
        return None


def compare_output(actual: str, expected: str, label: str, abs_floor: float = 0.0,
                   abs_tol=None):
    """Byte-compare, or (YREC_TEST_TOLERANT=1) compare numbers to REL_TOL;
    two numbers both below abs_floor in magnitude (round-off around zero)
    count as equal, and so do two within abs_tol[(line, token)] of each
    other (line numbers from 1, tokens from 0, as split() gives them)."""
    if actual == expected:
        return
    if os.environ.get("YREC_TEST_TOLERANT") != "1":
        raise AssertionError(
            f"{label} output differs from its pinned baseline; "
            "if the change is deliberate, regenerate and commit the baseline."
        )
    a_lines = actual.splitlines()
    e_lines = expected.splitlines()
    assert len(a_lines) == len(e_lines), (
        f"{label}: line count {len(a_lines)} != baseline {len(e_lines)}"
    )
    for ln, (al, el) in enumerate(zip(a_lines, e_lines), 1):
        if al == el:
            continue
        at, et = al.split(), el.split()
        assert len(at) == len(et), f"{label}:{ln}: token count differs\n{al}\n{el}"
        for tok, (atok, etok) in enumerate(zip(at, et)):
            if atok == etok:
                continue
            af, ef = _to_float(atok), _to_float(etok)
            assert af is not None and ef is not None, (
                f"{label}:{ln}: non-numeric difference {atok!r} vs {etok!r}"
            )
            if abs(af) < abs_floor and abs(ef) < abs_floor:
                continue
            if abs_tol and abs(af - ef) <= abs_tol.get((ln, tok), -1.0):
                continue
            denom = max(abs(af), abs(ef), 1.0e-300)
            assert abs(af - ef) / denom <= REL_TOL, (
                f"{label}:{ln}: {atok} vs {etok} beyond rel tol {REL_TOL}"
            )

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
"""
import os
import pathlib

# Absolute table/startmodel roots for every test-driven yrec run.
# Decks reference {YREC_INPUT}/{YREC_START}; the in-code fallback
# ('../../input') assumes the examples/<case>/ depth and silently
# breaks for testsuite/ (one level deep), so the harness pins the
# real locations once, absolutely, unless the caller already did.
_REPO = pathlib.Path(__file__).parent
os.environ.setdefault("YREC_INPUT", str(_REPO / "input"))
os.environ.setdefault("YREC_START", str(_REPO / "startmodels"))

REL_TOL = 1.0e-9


def _to_float(tok):
    try:
        return float(tok)
    except ValueError:
        return None


def compare_output(actual: str, expected: str, label: str):
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
        for atok, etok in zip(at, et):
            if atok == etok:
                continue
            af, ef = _to_float(atok), _to_float(etok)
            assert af is not None and ef is not None, (
                f"{label}:{ln}: non-numeric difference {atok!r} vs {etok!r}"
            )
            denom = max(abs(af), abs(ef), 1.0e-300)
            assert abs(af - ef) / denom <= REL_TOL, (
                f"{label}:{ln}: {atok} vs {etok} beyond rel tol {REL_TOL}"
            )

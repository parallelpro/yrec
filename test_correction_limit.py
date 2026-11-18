"""The limit on the Henyey correction, solver_max_correction, on the deck
that needs it.

examples/run_from_seed_to_start/m3000gs98z040000_prems rescales a 2.8
Msun seed to 3 Msun and from Z = 0.0176 to 0.04, and evolves it for five
models. In its second model Newton asks for a change of about one dex in
the central temperature; with most of it applied, the next linearization
has no physical state and the run stops with eos's "beta = P_gas/P is not
positive" (docs/bug_sweep.md M316). The deck sets solver_max_correction
= 0.1 and runs to its stop.

Two cases: the deck as shipped reaches its five models; with the limit
off (0, the default) it stops with exit status 1 and the error. Nothing
is byte-compared. Each run takes about a second.
"""
import os
import pathlib
import re
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
DECK = REPO / "examples" / "run_from_seed_to_start" / "m3000gs98z040000_prems.inlist"
TABLE = REPO / "input" / "eos" / "opal2006" / "EOSOPAL06Z0.040000000"


def _run(tmp_path, text):
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    if not TABLE.exists():
        pytest.skip("the Z = 0.04 OPAL 2006 EOS table is not in this checkout")
    (tmp_path / DECK.name).write_text(text)
    out = tmp_path / "out"
    out.mkdir()
    env = dict(os.environ)
    env.update({"YREC_INPUT": str(REPO / "input"),
                "YREC_START": str(REPO / "startmodels"),
                "YREC_OUTPUT": str(out)})
    return subprocess.run([str(YREC), DECK.name], cwd=tmp_path, env=env,
                          capture_output=True, text=True, timeout=600)


def test_deck_runs_with_its_limit(tmp_path):
    text = DECK.read_text()
    assert re.search(r"(?m)^\s*solver_max_correction\s*=\s*0\.1\b", text), \
        "the deck no longer sets solver_max_correction = 0.1"
    result = _run(tmp_path, text)
    tail = result.stdout[-1500:] + result.stderr[-500:]
    assert result.returncode == 0, tail
    assert "run finished: model budget exhausted" in result.stdout, tail
    final = re.findall(r"final model\s+(\d+)", result.stdout)
    assert final and int(final[-1]) == 5, tail


def test_deck_stops_without_the_limit(tmp_path):
    text, n = re.subn(r"(?m)^(\s*solver_max_correction\s*=\s*)0\.1\b", r"\g<1>0.0",
                      DECK.read_text())
    assert n == 1
    result = _run(tmp_path, text)
    assert result.returncode == 1, result.stdout[-1500:]
    assert "beta = P_gas/P is not positive" in result.stdout + result.stderr

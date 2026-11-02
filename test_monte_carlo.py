"""Two Monte-Carlo runs of the same samples are the same run (docs/bug_sweep.md M169).

Each Monte-Carlo run is the deck's run again with other samples: begin_monte_carlo_run
(core/run/monte_carlo.f90) re-arms every run copy of a control and resets the state a run
leaves on the star that the next run would read before writing it (the start-of-step
deuterium rates). This test runs the shipped deck examples/run_config_matrix/cm_monte_carlo
with two runs of identical samples and checks that the second run's history rows are the
first's, character for character. Before the reset the second run differed from the first
at the 1e-11 level here (and by much more on decks that burn deuterium for longer).
Local test: needs the full input/ tree.
"""
import os
import pathlib
import re
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
DECK = REPO / "examples" / "run_config_matrix"
STEM = "cm_monte_carlo"


def history_runs(path):
    """the history file's data rows, split into Monte-Carlo runs: the model number restarts at
    every run's first card (this deck reads its starting model on card 1 only)"""
    lines = path.read_text().splitlines()
    runs, last = [], None
    for line in lines[6:]:
        if not line.strip():
            continue
        number = int(line.split()[0])
        if last is None or number < last:
            runs.append([])
        runs[-1].append(line)
        last = number
    return runs


def test_identical_samples_give_identical_runs(tmp_path):
    if not YREC.exists():
        pytest.skip("src/yrec not built (run `make` in src/)")
    if not (REPO / "input" / "eos").exists():
        pytest.skip("needs the full input/ tree")
    text = (DECK / f"{STEM}.inlist").read_text()
    assert re.search(r"(?m)^ mc_run_end = 1\s*$", text), "cm_monte_carlo no longer makes one run"
    text = re.sub(r"(?m)^ mc_run_end = 1\s*$", " mc_run_end = 2", text).replace(STEM, "mc_twice")
    (tmp_path / "mc_twice.inlist").write_text(text)
    record = (DECK / "mc_input.dat").read_text().splitlines()[:2]
    assert record[0].startswith("run  1"), record
    (tmp_path / "mc_input.dat").write_text(
        "\n".join(record + [record[0].replace("run  1", "run  2", 1), record[1]]) + "\n")
    (tmp_path / "output").mkdir()
    env = dict(os.environ, YREC_INPUT=str(REPO / "input"), YREC_START=str(REPO / "startmodels"))
    result = subprocess.run([str(YREC), "mc_twice.inlist"], cwd=tmp_path, env=env,
                            capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, result.stdout[-2000:] + result.stderr[-2000:]
    runs = history_runs(tmp_path / "output" / "history.data")
    assert len(runs) == 2 and runs[0], f"expected two runs, got {[len(r) for r in runs]}"
    first_difference = next((i for i, (a, b) in enumerate(zip(*runs)) if a != b), None)
    assert len(runs[0]) == len(runs[1]) and first_difference is None, (
        f"run 2 differs from run 1 ({len(runs[0])} and {len(runs[1])} rows; first differing row "
        f"{first_difference}): state from run 1 carried into run 2")

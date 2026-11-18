"""The shipped low-mass decks run to their stops.

Each of these decks used to end with "the previous model failed to
converge" (docs/bug_sweep.md O3, M294): the blended equation of state
jumped along a fully convective model at the OPAL 2006 table's ragged
edge (M356), and where the OPAL density search gave up (M357). Each
case runs the deck as shipped, from its own directory, and requires a
clean exit, the run's last line naming the stop the deck asks for, and
no density-search fallback in the run log. Nothing is byte-compared,
so the test is the same on every platform.

Run times (release build, macOS/arm64, 2026-11-17): 10 to 30 s a deck,
2 minutes in all.
"""
import os
import pathlib
import re
import shutil
import subprocess

import pytest

from conftest import INPUT_ROOT, START_ROOT

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"

# (deck, the reason the run must end with, the fewest models it takes)
CASES = [
    ("examples/run_from_dbl_to_zams/Test_m0040_feh0_noGS_ZAMS.inlist",
     "central X below stop value", 200),
    ("examples/run_from_dbl_to_zams/Test_m0050_feh0_noGS_ZAMS.inlist",
     "central X below stop value", 200),
    ("examples/run_from_dbl_to_zams/Test_m0060_feh0_noGS_ZAMS.inlist",
     "central X below stop value", 150),
    ("testsuite/Test_m0p05feh+0p0GS98_allard_15Gyr.inlist",
     "target age reached ( 1.500E+10 yr)", 300),
    ("testsuite/Test_m0p1feh+0p0GS98_allard_ZAMS.inlist",
     "central X below stop value", 150),
    ("sample_cases/lowerMS/SamplelowerMS_m0p15feh+0p0GS98_kurucz_ZAMS.inlist",
     "central X below stop value", 150),
]


@pytest.mark.parametrize("deck,reason,min_models", CASES,
                         ids=[pathlib.Path(c[0]).stem for c in CASES])
def test_low_mass_deck_reaches_its_stop(deck, reason, min_models, tmp_path):
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    source = REPO / deck
    if not source.exists():
        pytest.skip(f"{deck} is not in this checkout")
    text = source.read_text()
    # the deck's files beside it (a parameter file, an .nml2), and the
    # output subdirectories it names
    work = tmp_path / "deck"
    work.mkdir()
    for other in source.parent.iterdir():
        if other.is_file() and other.suffix not in (".inlist", ".nml1", ".log") \
                and other.stat().st_size < 5_000_000:
            shutil.copy(other, work / other.name)
    (work / source.name).write_text(text)
    out = tmp_path / "out"
    out.mkdir()
    for sub in set(re.findall(r'\{YREC_OUTPUT\}/([^"]*)/[^/"]*"', text)):
        (out / sub).mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env.update({"YREC_INPUT": str(INPUT_ROOT),
                "YREC_START": str(START_ROOT),
                "YREC_OUTPUT": str(out)})
    result = subprocess.run([str(YREC), source.name], cwd=work, env=env,
                            capture_output=True, text=True, timeout=1200)
    tail = result.stdout[-1500:] + result.stderr[-500:]
    assert result.returncode == 0, tail
    finished = re.findall(r"run finished: (.*)", result.stdout)
    assert finished and finished[-1].strip() == reason, tail
    final = re.findall(r"final model\s+(\d+)", result.stdout)
    assert final and int(final[-1]) >= min_models, tail
    logs = "".join(p.read_text(errors="replace")
                   for p in out.rglob("*") if p.is_file()
                   and p.suffix in (".log", ".short"))
    assert "falling back" not in logs, "an OPAL density search fell back"

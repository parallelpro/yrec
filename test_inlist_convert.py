"""Inlist-revamp acceptance test (2026).

For each covered example case: convert the legacy .nml1/.nml2 pair to a
single new-style inlist (&star_job + &controls, readable names) with
tools/upgrade_inlist.py, run yrec on BOTH, and require byte-identical
outputs (version-hash line excluded). This is the whole revamp's
correctness contract: renaming the user-facing controls must not change
one byte of physics.

Local test: needs the full input/ tree (not in CI's sparse checkout).
"""
import os
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
CONVERTER = REPO / "src" / "tools" / "upgrade_inlist.py"
CASE = REPO / "examples" / "run_standard_solar_model"
PAIRS = [
    ("Test_solar_noGS_norot.nml1", "Test_solar_noGS_norot.nml2"),
    ("Test_solar_GS_rot_gray_OPALSCV_SF2_GS98_OP_CF10.nml1",
     "Test_solar_GS_rot_gray_OPALSCV_SF2_GS98_OP_CF10.nml2"),
    ("Test_solarcal_noGS_norot.nml1", "Test_solarcal_noGS_norot.nml2"),
]
ENV = {"YREC_INPUT": REPO / "input", "YREC_START": REPO / "startmodels"}


def _strip(path):
    return [l for l in path.read_text(errors="replace").splitlines()
            if not l.startswith("# YREC v")]


def _run(workdir, args):
    workdir.mkdir()
    (workdir / "output").mkdir()
    env = dict(os.environ)
    env.update({k: str(v) for k, v in ENV.items()})
    r = subprocess.run([str(YREC)] + args, cwd=workdir, env=env,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]
    return workdir / "output"


@pytest.mark.parametrize("nml1,nml2", PAIRS)
def test_converted_inlist_is_byte_identical(tmp_path, nml1, nml2):
    if not YREC.exists():
        pytest.skip("build src/yrec first")
    if not (REPO / "input").exists():
        pytest.skip("needs the full input/ tree (local only)")

    newstyle = tmp_path / "newstyle"

    # legacy run
    legacy = tmp_path / "legacy"
    legacy.mkdir()
    (legacy / "output").mkdir()
    for n in (nml1, nml2):
        (legacy / n).write_text((CASE / n).read_text())
    env = dict(os.environ)
    env.update({k: str(v) for k, v in ENV.items()})
    r = subprocess.run([str(YREC), nml1, nml2], cwd=legacy, env=env,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]

    # convert + new-style run (single file)
    newstyle.mkdir()
    (newstyle / "output").mkdir()
    inlist = newstyle / "inlist"
    conv = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE / nml1), str(CASE / nml2),
         "-o", str(inlist)], capture_output=True, text=True)
    assert conv.returncode == 0, conv.stderr
    r = subprocess.run([str(YREC), "inlist"], cwd=newstyle, env=env,
                       capture_output=True, text=True, timeout=600)
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]

    legacy_files = sorted(p.name for p in (legacy / "output").iterdir())
    new_files = sorted(p.name for p in (newstyle / "output").iterdir())
    assert legacy_files == new_files, (legacy_files, new_files)
    # inlist_used is the verbatim copy of the INPUT inlists -- it
    # legitimately differs between the legacy pair and the converted
    # single inlist; everything the run COMPUTES must still match.
    legacy_files = [n for n in legacy_files if n != "inlist_used"]
    diffs = [n for n in legacy_files
             if _strip(legacy / "output" / n) != _strip(newstyle / "output" / n)]
    assert not diffs, (
        f"converted inlist changed outputs: {diffs} -- a rename in "
        "defaults/controls_registry.tsv or the generated reader is not "
        "value-preserving"
    )

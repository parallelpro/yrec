"""Byte-pinned baselines for the MESA-format outputs.

Runs the standard solar cases in MESA mode (converted single inlist,
profiles + pulse on) and byte-compares every produced data file --
history.data, profile{N}.data, pulse files, and the .mod model --
against examples/run_standard_solar_model/standard/mesa/<case>/.

This is the MESA-format counterpart of the legacy .short pin (the
.track/.store pins it superseded are retired): the retire-legacy
campaign only removes a legacy file once its content is pinned here
in the format that replaces it. Unlike the
legacy files, the MESA files carry no git hash (version_number is the
bare release string), so the comparison is byte-exact with no masking.

Reference handling matches test_all.py: the first run seeds the
reference directory (reported as a skip so seeding is visible), and
baselines live under the gitignored **/standard/ tree -- they are
machine-local and regenerated deliberately, never committed.
"""
import os
import pathlib
import shutil
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
CONVERTER = REPO / "src" / "tools" / "upgrade_inlist.py"
CASE_DIR = REPO / "examples" / "run_standard_solar_model"
REF_ROOT = CASE_DIR / "standard" / "mesa"
ENV = {"YREC_INPUT": REPO / "input", "YREC_START": REPO / "startmodels"}

# noGS_norot pins the baseline physics; the GS_rot case additionally
# exercises the rotation history/profile columns (omega, j_rot, D_omega,
# D_mix_rot) and the wind/turnover path that feeds them.
CASES = [
    "Test_solar_noGS_norot",
    "Test_solar_GS_rot_gray_OPALSCV_SF2_GS98_OP_CF10",
]


def _mesa_inlist(case, tmp):
    conv = tmp / f"{case}.inlist"
    r = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE_DIR / f"{case}.nml1"),
         str(CASE_DIR / f"{case}.nml2"), "-o", str(conv)],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    text = conv.read_text()
    assert "use_legacy_output" not in text   # retired; nothing stamped
    text = text.replace(
        "&controls",
        "&controls\n write_profile_flag = .true.\n write_pulse_flag = .true.\n",
        1)
    return text


def _pinned_names(outdir):
    """Every data file the run produced that this test pins: history,
    profiles (and their pulse companions), and the .mod model. The .log
    run log carries the git-hash banner and is not pinned."""
    names = sorted(p.name for p in outdir.glob("profile*"))
    names += sorted(p.name for p in outdir.glob("*.mod"))
    if (outdir / "history.data").is_file():
        names.insert(0, "history.data")
    return names


@pytest.mark.parametrize("case", CASES)
def test_mesa_baseline(case, tmp_path):
    work = tmp_path / case
    work.mkdir()
    (work / "output").mkdir()
    (work / "inlist").write_text(_mesa_inlist(case, tmp_path))
    env = dict(os.environ)
    env.update({k: str(v) for k, v in ENV.items()})
    r = subprocess.run([str(YREC), "inlist"], cwd=work, env=env,
                       capture_output=True, text=True, timeout=1200)
    assert r.returncode == 0, r.stdout[-500:] + r.stderr[-500:]

    outdir = work / "output"
    names = _pinned_names(outdir)
    assert "history.data" in names, names
    assert any(n.startswith("profile") for n in names), names

    ref = REF_ROOT / case
    if not ref.is_dir():
        ref.mkdir(parents=True)
        for n in names:
            shutil.copyfile(outdir / n, ref / n)
        pytest.skip(f"seeded MESA baseline: {ref}")

    ref_names = sorted(p.name for p in ref.iterdir())
    assert ref_names == sorted(names), (
        f"produced file set changed: baseline {ref_names} vs run {sorted(names)}")
    diffs = [n for n in names
             if (ref / n).read_bytes() != (outdir / n).read_bytes()]
    assert not diffs, f"byte difference vs baseline in: {diffs}"

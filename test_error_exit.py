"""End-to-end error-path test for the run_yrec engine (phase five, step B).

Historically every fatal condition in YREC ended in a bare `stop`,
which exits with status 0 -- indistinguishable from success to any
pipeline. Phase B converts the driver-layer stops into ierr returns
that surface through run_yrec to the CLI wrapper, which exits 1. This
test drives the oldest configuration check in parmin (semi-convection
and overshoot both enabled, rejected since the original F77) and
asserts the new contract: a clean nonzero exit, no crash, no model
evolved. It runs in an isolated copy of a Stage-0 case directory so
no example outputs are touched.
"""
import os
import pathlib
import re
import shutil
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
CASE_DIR = REPO / "examples" / "run_standard_solar_model"
CONVERTER = REPO / "src" / "tools" / "upgrade_inlist.py"
NML1 = "Test_solar_noGS_norot.nml1"
NML2 = "Test_solar_noGS_norot.nml2"


def test_config_error_exits_nonzero(tmp_path):
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    for name in (NML1, NML2):
        shutil.copy(CASE_DIR / name, tmp_path / name)
    nml2 = tmp_path / NML2
    text = nml2.read_text()
    text = re.sub(r"LSEMIC\s*=\s*\.FALSE\.", "LSEMIC = .TRUE.", text)
    text = re.sub(r"LOVSTC\s*=\s*\.FALSE\.", "LOVSTC = .TRUE.", text)
    nml2.write_text(text)
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(YREC), NML1, NML2], cwd=tmp_path,
        env=env, capture_output=True, text=True, timeout=120,
    )
    # the new contract: a config error is a clean nonzero exit (stop 1
    # in the CLI wrapper), where historically a bare `stop` exited 0
    assert result.returncode == 1, (
        "conflicting LSEMIC/LOVSTC must exit 1 "
        f"(got {result.returncode}); stderr tail: {result.stderr[-300:]}"
    )
    # and the diagnostic still prints at the point of failure
    blobs = [result.stdout, result.stderr]
    for f in tmp_path.rglob("*"):
        if f.is_file() and f.suffix not in (".nml1", ".nml2"):
            try:
                blobs.append(f.read_text(errors="replace"))
            except OSError:
                pass
    assert any("SEMI-CONVECTION" in b for b in blobs), (
        "expected the PARMIN conflict diagnostic in the run output"
    )


def test_bad_opacity_table_exits_nonzero(tmp_path):
    """A malformed opacity table must surface as a clean nonzero exit
    (2026 ierr-not-stop: the table-reader stops became ierr returns
    threaded up through setupopac/kap_init/setups/star_setup)."""
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    for name in (NML1, NML2):
        shutil.copy(CASE_DIR / name, tmp_path / name)
    bad_table = tmp_path / "bad_opal95_table"
    # a TABLE header whose X/Z fields (cols 37-43 / 55-61) read as
    # 0.9 -- incompatible with the expected first grid table, failing
    # the reader's compatibility check: the shortest path into the
    # table-reader ierr plumbing
    bad_table.write_text("TABLE" + " " * 31 + " 0.9000"
                         + " " * 11 + " 0.9000\n")
    nml1 = tmp_path / NML1
    text = nml1.read_text()
    text = re.sub(r'FLIV95\s*=\s*"[^"]*"', f'FLIV95 = "{bad_table}"', text)
    nml1.write_text(text)
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(YREC), NML1, NML2], cwd=tmp_path,
        env=env, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 1, (
        "malformed OPAL95 table must exit 1 "
        f"(got {result.returncode}); stderr tail: {result.stderr[-300:]}"
    )
    blobs = [result.stdout, result.stderr]
    for f in tmp_path.rglob("*"):
        if f.is_file() and f.suffix not in (".nml1", ".nml2"):
            try:
                blobs.append(f.read_text(errors="replace"))
            except OSError:
                pass
    assert any("OPAL95" in b for b in blobs), (
        "expected the table-reader diagnostic in the run output"
    )


def test_advective_transport_rejected(tmp_path):
    """LDIFAD (use_diffusion_advection_transport) has never worked: its
    equal-grid advective coefficients are never computed. Where it would
    take effect -- a rotating run with instability transport -- the run
    must stop with a clean nonzero exit instead of dividing by zero."""
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    rot = "Test_solar_GS_rot_gray_OPALSCV_SF2_GS98_OP_CF10"
    for ext in (".nml1", ".nml2"):
        shutil.copy(CASE_DIR / (rot + ext), tmp_path / (rot + ext))
    nml2 = tmp_path / (rot + ".nml2")
    text, n = re.subn(r"LDIFAD\s*=\s*\.FALSE\.", "LDIFAD = .TRUE.",
                      nml2.read_text())
    assert n == 1, "expected one LDIFAD = .FALSE. in the rotating deck"
    nml2.write_text(text)
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(YREC), rot + ".nml1", rot + ".nml2"], cwd=tmp_path,
        env=env, capture_output=True, text=True, timeout=600,
    )
    assert result.returncode == 1, (
        "LDIFAD in a rotating run must exit 1 "
        f"(got {result.returncode}); stdout tail: {result.stdout[-300:]}"
    )
    assert "use_diffusion_advection_transport" in result.stdout, (
        "expected the secular_transport diagnostic on the terminal"
    )


def test_old_group_names_refused(tmp_path):
    """An inlist whose groups carry their 2026-08 names (&star_job and
    &controls, now &control and &physics) must stop the run with exit 1
    and the command that renames them -- and that command must give back
    the inlist exactly."""
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    inlist = tmp_path / "deck.inlist"
    conv = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE_DIR / NML1),
         str(CASE_DIR / NML2), "-o", str(inlist)],
        capture_output=True, text=True)
    assert conv.returncode == 0, conv.stderr
    current = inlist.read_text()
    old = re.sub(r"(?m)^&control\b", "&star_job", current, count=1)
    old = re.sub(r"(?m)^&physics\b", "&controls", old, count=1)
    assert old.count("&star_job") == 1 and old.count("&controls") == 1
    inlist.write_text(old)
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(YREC), "deck.inlist"], cwd=tmp_path,
        env=env, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 1, (
        "a &star_job inlist must exit 1 "
        f"(got {result.returncode}); stdout tail: {result.stdout[-300:]}"
    )
    assert "upgrade_inlist.py deck.inlist" in result.stdout, result.stdout[-300:]
    fix = subprocess.run([sys.executable, str(CONVERTER), str(inlist)],
                         capture_output=True, text=True)
    assert fix.returncode == 0, fix.stderr
    assert inlist.read_text() == current


def test_nonpositive_triangle_size_refused(tmp_path):
    """A zero envelope-triangle size (TRIDT = 0) once made the
    triangle's reflection walk loop forever; the controls read now
    refuses it with exit 1 and a message naming the control."""
    if not YREC.exists():
        pytest.skip("src/yrec not built")
    for name in (NML1, NML2):
        shutil.copy(CASE_DIR / name, tmp_path / name)
    nml2 = tmp_path / NML2
    text, n = re.subn(r"(?im)^(\s*TRIDT\s*=\s*)[-+.0-9dDeE]+", r"\g<1>0.0", nml2.read_text())
    assert n == 1, "expected one TRIDT in the solar deck"
    nml2.write_text(text)
    (tmp_path / "output").mkdir()
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(YREC), NML1, NML2], cwd=tmp_path,
        env=env, capture_output=True, text=True, timeout=120,
    )
    assert result.returncode == 1, (
        f"TRIDT = 0 must exit 1 (got {result.returncode}); "
        f"stdout tail: {result.stdout[-300:]}"
    )
    assert "envelope triangle's size must be positive" in result.stdout, result.stdout[-300:]

"""Phase-five step C acceptance test (2026 -- ROADMAP.md): re-entrancy.

The contract: two run_yrec calls in one process produce, on the second
call, outputs byte-identical to a fresh single-process run. src/
test_reentry runs the engine twice with the same namelists; this test
runs a fresh `yrec` once in one working directory and `test_reentry`
in another, then byte-compares every output file (version-hash line
excluded, as in the Stage-0 discipline).

Uses the solar noGS/norot case (~18 s), which terminates NORMALLY
(run count exhausted). The m0030 case cannot serve here: it
historically ends via the BSSTEP numerics stop, which kills the
process mid-call -- the documented residual (numerics gates are
opt-in) means such cases are not re-enterable until their callers
pass ierr. Three engine runs, ~1 minute total. Local test: needs
the full input/ tree (not in CI's sparse checkout).
"""
import os
import pathlib
import shutil
import subprocess

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
REENTRY = REPO / "src" / "test_reentry"
CASE = REPO / "examples" / "run_standard_solar_model"
NML1 = "Test_solar_noGS_norot.nml1"
NML2 = "Test_solar_noGS_norot.nml2"


def _run(binary, workdir, extra_args=()):
    (workdir / "output").mkdir()
    for name in (NML1, NML2):
        shutil.copy(CASE / name, workdir / name)
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(binary), *extra_args, NML1, NML2], cwd=workdir,
        env=env, capture_output=True, text=True, timeout=600,
    )
    return result


# A short variant of the solar case that sets a legacy-spelled PHYSICS
# control the solar case itself omits (S0_1_1, the p-p S-factor, doubled)
# and writes under its own file stem. Used as the FIRST run of a
# two-inlist process: before the read_controls fix (bugsweep sec-11
# batch 0) the DATA-initialised NAMELIST locals were SAVEd, so the
# second run inherited 8.02d-22 instead of the 4.01d-22 default.
A_STEM = "A_carryover"


def _write_case_a(workdir):
    nml1 = (CASE / NML1).read_text()
    nml1 = nml1.replace("NMODLS(2) = 1000", "NMODLS(2) = 5")
    nml1 = nml1.replace("Test_solar_noGS_norot", A_STEM)
    assert "NMODLS(2) = 5" in nml1 and A_STEM in nml1
    nml2 = (CASE / NML2).read_text()
    marker = " LNEWNUC = .TRUE."
    assert marker in nml2 and "S0_1_1" not in nml2
    nml2 = nml2.replace(marker, marker + "\n S0_1_1 = 8.02D-22", 1)
    (workdir / f"{A_STEM}.nml1").write_text(nml1)
    (workdir / f"{A_STEM}.nml2").write_text(nml2)
    return f"{A_STEM}.nml1", f"{A_STEM}.nml2"


def _strip(path):
    return [l for l in path.read_text(errors="replace").splitlines()
            if not l.startswith("# YREC v")]


def _assert_second_matches_fresh(fresh, twice, ignore_stem=None):
    fresh_files = sorted(p.name for p in (fresh / "output").iterdir())
    twice_files = sorted(p.name for p in (twice / "output").iterdir()
                         if ignore_stem is None or ignore_stem not in p.name)
    assert fresh_files == twice_files, (fresh_files, twice_files)
    diffs = []
    for name in fresh_files:
        if _strip(fresh / "output" / name) != _strip(twice / "output" / name):
            diffs.append(name)
    assert not diffs, (
        f"second in-process run differs from a fresh process in: {diffs} "
        "-- some SAVEd state is escaping the reset (see "
        "core/yrec_reset.f90 and the ROADMAP step-C notes)"
    )


def test_second_call_equals_fresh_process(tmp_path):
    if not YREC.exists() or not REENTRY.exists():
        pytest.skip("build src/yrec and src/test_reentry first")
    fresh = tmp_path / "fresh"
    twice = tmp_path / "twice"
    fresh.mkdir()
    twice.mkdir()
    r1 = _run(YREC, fresh)
    assert r1.returncode == 0, r1.stdout[-500:] + r1.stderr[-500:]
    r2 = _run(REENTRY, twice)
    assert r2.returncode == 0, r2.stdout[-500:] + r2.stderr[-500:]
    assert "second run done" in r2.stdout
    _assert_second_matches_fresh(fresh, twice)


def test_second_call_does_not_inherit_omitted_controls(tmp_path):
    """Run A (S0_1_1 doubled, 5 models) then B (the solar case, which
    omits S0_1_1) in one process; B must match a fresh B. Guards the
    read_controls reset block: a control set by an earlier inlist must
    not leak into a later one that leaves it at its default."""
    if not YREC.exists() or not REENTRY.exists():
        pytest.skip("build src/yrec and src/test_reentry first")
    fresh = tmp_path / "fresh"
    twice = tmp_path / "twice"
    fresh.mkdir()
    twice.mkdir()
    r1 = _run(YREC, fresh)
    assert r1.returncode == 0, r1.stdout[-500:] + r1.stderr[-500:]
    a_nml = _write_case_a(twice)
    r2 = _run(REENTRY, twice, extra_args=a_nml)
    assert r2.returncode == 0, r2.stdout[-500:] + r2.stderr[-500:]
    assert "second run done" in r2.stdout
    # run A must really have taken the doubled S-factor, else the test
    # proves nothing
    assert (twice / "output" / f"{A_STEM}.log").exists()
    _assert_second_matches_fresh(fresh, twice, ignore_stem=A_STEM)

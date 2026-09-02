"""Phase-five step C acceptance test (2026 -- ROADMAP.md): re-entrancy.

The contract: two run_yrec calls in one process produce, on the second
call, outputs byte-identical to a fresh single-process run. src/
test_reentry runs the engine twice (same inlist, or A then B); this
test runs a fresh `yrec` once in one working directory and
`test_reentry` in another, then byte-compares every output file
(version-hash line excluded, as in the Stage-0 discipline).

Uses the solar noGS/norot case (~18 s), which terminates NORMALLY
(run count exhausted), converted to a single new-style inlist at
test time (tools/upgrade_inlist.py, as test_mesa_baselines does) --
the legacy .nml1/.nml2 pair is kept only as the converter's fixture.
The m0030 case cannot serve here: it historically ends via the
BSSTEP numerics stop, which kills the process mid-call -- the
documented residual (numerics gates are opt-in) means such cases are
not re-enterable until their callers pass ierr. Local test: needs
the full input/ tree (not in CI's sparse checkout).
"""
import os
import pathlib
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
YREC = REPO / "src" / "yrec"
REENTRY = REPO / "src" / "test_reentry"
CONVERTER = REPO / "src" / "tools" / "upgrade_inlist.py"
CASE = REPO / "examples" / "run_standard_solar_model"
STEM = "Test_solar_noGS_norot"
INLIST = f"{STEM}.inlist"


def _write_solar_inlist(workdir):
    conv = workdir / INLIST
    r = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE / f"{STEM}.nml1"),
         str(CASE / f"{STEM}.nml2"), "-o", str(conv)],
        capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return conv.read_text()


def _run(binary, workdir, inlists):
    """Write the solar inlist into workdir and run `binary` on the given
    inlist names (any other named inlist must already exist there)."""
    (workdir / "output").mkdir()
    _write_solar_inlist(workdir)
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env["YREC_START"] = str(REPO / "startmodels")
    result = subprocess.run(
        [str(binary), *inlists], cwd=workdir,
        env=env, capture_output=True, text=True, timeout=600,
    )
    return result


# A short variant of the solar case that sets a control the solar case
# itself omits (s0_pp, the p-p S-factor, doubled) and writes under its
# own file stem. Used as the FIRST run of a two-inlist process: before
# the read_controls fix (bugsweep sec-11 batch 0) the DATA-initialised
# NAMELIST carriers were SAVEd -- and the new-style read seeds s0_pp
# from that carrier -- so the second run inherited 8.02d-22 instead of
# the 4.01d-22 default.
A_STEM = "A_carryover"


def _write_case_a(workdir):
    text = _write_solar_inlist(workdir)
    text = text.replace("max_model_number(2) = 1000", "max_model_number(2) = 5")
    text = text.replace(STEM, A_STEM)
    marker = " use_new_nuclear_rates = .TRUE."
    assert marker in text and "s0_pp" not in text and A_STEM in text
    assert "max_model_number(2) = 5" in text
    text = text.replace(marker, marker + "\n s0_pp = 8.02D-22", 1)
    (workdir / f"{A_STEM}.inlist").write_text(text)
    return f"{A_STEM}.inlist"


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
    r1 = _run(YREC, fresh, [INLIST])
    assert r1.returncode == 0, r1.stdout[-500:] + r1.stderr[-500:]
    r2 = _run(REENTRY, twice, [INLIST])
    assert r2.returncode == 0, r2.stdout[-500:] + r2.stderr[-500:]
    assert "second run done" in r2.stdout
    _assert_second_matches_fresh(fresh, twice)


def test_second_call_does_not_inherit_omitted_controls(tmp_path):
    """Run A (s0_pp doubled, 5 models) then B (the solar case, which
    omits s0_pp) in one process; B must match a fresh B. Guards the
    read_controls reset block: a control set by an earlier inlist must
    not leak into a later one that leaves it at its default."""
    if not YREC.exists() or not REENTRY.exists():
        pytest.skip("build src/yrec and src/test_reentry first")
    fresh = tmp_path / "fresh"
    twice = tmp_path / "twice"
    fresh.mkdir()
    twice.mkdir()
    r1 = _run(YREC, fresh, [INLIST])
    assert r1.returncode == 0, r1.stdout[-500:] + r1.stderr[-500:]
    a_inlist = _write_case_a(twice)
    r2 = _run(REENTRY, twice, [a_inlist, INLIST])
    assert r2.returncode == 0, r2.stdout[-500:] + r2.stderr[-500:]
    assert "second run done" in r2.stdout
    # run A must really have taken the doubled S-factor, else the test
    # proves nothing
    assert (twice / "output" / f"{A_STEM}.log").exists()
    _assert_second_matches_fresh(fresh, twice, ignore_stem=A_STEM)

#!/usr/bin/env python3
"""Line coverage of a source directory, measured with gcov.

  python3 src/tools/coverage.py --work DIR [--dir rotation] [--deck PATH ...]
                                [--test BINARY ...] [--min PCT]

Copies src/ to DIR/src (the working tree is not touched), builds yrec and
the test programs there with -O0 --coverage, runs the test programs and
the decks with that build (each deck in a fresh directory, with
YREC_INPUT and YREC_START pointing at the repository's), then prints, for
each file under src/<dir> (test/ excluded), gcov's line coverage: the
executable lines run at least once over the executable lines. Comments
and commented-out code are not executable lines, so they do not count;
nor does an `end module` line, where gfortran places the procedures it
generates for a module's derived types.
Exits 1 when the total is below --min.

Defaults: --dir rotation; the tests test_rotation and
test_rotation_kernels; the decks examples/run_config_matrix/cm_rot_*.inlist
(fully convective: they do not reach the secular transport, so pass
decks with a radiative core for a full measure); --jobs 4 decks at once.
"""
import argparse
import concurrent.futures
import glob
import os
import pathlib
import re
import shutil
import subprocess
import sys

REPO = pathlib.Path(__file__).resolve().parents[2]


def is_binary(path):
    """A compiled executable (Mach-O or ELF), not a source or a script."""
    if not path.is_file():
        return False
    with open(path, "rb") as f:
        head = f.read(4)
    return head in (b"\xcf\xfa\xed\xfe", b"\xce\xfa\xed\xfe", b"\x7fELF")


def build(work):
    src = work / "src"
    if src.exists():
        shutil.rmtree(src)
    def skip(directory, names):
        # build products: objects, modules, archives, coverage data and
        # the compiled executables at the top of src/
        top = pathlib.Path(directory) == REPO / "src"
        return [n for n in names
                if n.endswith((".o", ".mod", ".smod", ".a", ".gcda", ".gcno", ".gcov"))
                or n == "__pycache__"
                or (top and is_binary(pathlib.Path(directory) / n))]
    shutil.copytree(REPO / "src", src, ignore=skip)
    flags = ["GNU_REL_FFLAGS=-O0 --coverage", "LDFLAGS=--coverage"]
    for target in (["default"], ["tests"]):
        r = subprocess.run(["make", "-j8"] + flags + target, cwd=src,
                           capture_output=True, text=True)
        if r.returncode != 0:
            sys.exit("coverage build failed:\n" + r.stdout[-3000:] + r.stderr[-3000:])
    return src


def run_tests(src, tests):
    for t in tests:
        r = subprocess.run([str(src / t)], cwd=src, capture_output=True, text=True)
        print(f"  {t}: exit {r.returncode}")


def run_decks(src, work, decks, jobs):
    env = dict(os.environ, YREC_INPUT=str(REPO / "input"),
               YREC_START=str(REPO / "startmodels"), YREC_OUTPUT="output")
    def one(d):
        deck = pathlib.Path(d)
        if not deck.is_absolute():
            deck = REPO / deck
        rundir = work / "runs" / deck.stem
        if rundir.exists():
            shutil.rmtree(rundir)
        (rundir / "output").mkdir(parents=True)
        shutil.copy(deck, rundir)
        r = subprocess.run([str(src / "yrec"), deck.name], cwd=rundir, env=env,
                           capture_output=True, text=True)
        return deck.name, r.returncode
    # gcov's counters are merged safely by concurrent runs of one binary
    with concurrent.futures.ThreadPoolExecutor(jobs) as pool:
        for name, rc in pool.map(one, decks):
            print(f"  {name}: exit {rc}")


def line_coverage(src, directory):
    rows = []
    files = sorted(p for p in (src / directory).rglob("*.f90") if "/test/" not in str(p))
    for f in files:
        rel = f.relative_to(src)
        obj_dir = f.parent
        r = subprocess.run(["gcov", "-t", "-o", str(obj_dir), str(f)], cwd=src,
                           capture_output=True, text=True)
        executable = covered = 0
        for line in r.stdout.splitlines():
            m = re.match(r"^\s*([^:]+):\s*(\d+):(.*)$", line)
            if not m or m.group(2) == "0":
                continue
            count = m.group(1).strip()
            # gfortran puts a module's generated copy and finalization
            # procedures for its derived types on the end module line
            if count == "-" or re.match(r"\s*end\s+module\b", m.group(3), re.I):
                continue
            executable += 1
            if count not in ("#####", "====="):
                covered += 1
        rows.append((str(rel), executable, covered))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--work", required=True)
    ap.add_argument("--dir", default="rotation")
    ap.add_argument("--deck", action="append")
    ap.add_argument("--test", action="append")
    ap.add_argument("--min", type=float, default=0.0)
    ap.add_argument("--nobuild", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    work = pathlib.Path(a.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    tests = a.test or ["test_rotation", "test_rotation_kernels"]
    decks = a.deck or sorted(glob.glob(str(REPO / "examples/run_config_matrix/cm_rot_*.inlist")))
    src = work / "src" if a.nobuild else build(work)
    print("tests:"); run_tests(src, tests)
    print("decks:"); run_decks(src, work, decks, a.jobs)
    rows = line_coverage(src, a.dir)
    tot_e = sum(r[1] for r in rows); tot_c = sum(r[2] for r in rows)
    print(f"\n{'file':52s} {'lines':>6s} {'run':>6s} {'cover':>7s}")
    for name, e, c in rows:
        if e:
            print(f"{name:52s} {e:6d} {c:6d} {100.0*c/e:6.1f}%")
    total = 100.0*tot_c/max(tot_e, 1)
    print(f"{'total':52s} {tot_e:6d} {tot_c:6d} {total:6.1f}%")
    if total < a.min:
        sys.exit(f"coverage {total:.1f}% is below {a.min}%")


if __name__ == "__main__":
    main()

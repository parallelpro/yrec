#!/usr/bin/env python3
"""Deck comparison: proves a refactor byte-identical on the decks and pins.

Usage (from anywhere; the tree is the checkout or worktree under test):
    python3 -I src/tools/compare_decks.py TREE [LABEL] [--set quick|full|rot|all]
        [--jobs N] [--only REGEX] [--nobuild] [--snapshot] [--work DIR] [--pins DIR]

Before a refactor, take a baseline from the commit it starts from; after
it (and after every commit of it), compare:
    python3 -I src/tools/compare_decks.py . --set all --snapshot
    python3 -I src/tools/compare_decks.py . --set all
A refactor is done when the second prints RESULT: ALL IDENTICAL.

What it does, in order:
  1. Builds TREE/src (make -j8) unless --nobuild, and runs
     tools/check_boundaries.py there (its last line is printed; a
     failure does not stop the decks). A failed build stops it.
  2. Runs every deck of the set with TREE/src/yrec, N at a time
     (--jobs, default 6), each from its own directory, with
     YREC_INPUT=TREE/input, YREC_START=TREE/startmodels and YREC_OUTPUT
     unset, so its files go to <deck directory>/output (git-ignored).
     A deck's old files there ({deck}.*) are removed before it runs.
  3. With --snapshot, copies each deck's output files ({deck}.*) into
     the baseline, replacing that deck's earlier ones. Otherwise it
     byte-compares every one of them with the baseline's copy, and
     reports a file the baseline lacks (new:NAME) and a baseline file
     the run did not write (gone:NAME). The only text ignored is the
     version line, "# YREC v<version> (<hash>[-dirty])", which differs
     between commits. The files not named after a deck, output/history.data
     and output/inlist_used, which the decks of a directory share, are
     not compared.
  4. Compares each deck's .log with its pin, PINS/<deck dir>/standard/
     <deck>.log, where there is one, with the same rule (a snapshot too).

The per-deck lines:
    IDENT      12s same= 2 diff=0 pin=OK   examples/run_config_matrix/cm_base
    DIFF       12s same= 1 diff=1 pin=DIFF examples/run_config_matrix/cm_base cm_base.log
  The tag is IDENT when the deck ran (exit status 0), wrote at least one
  file, and every file, the file list and the pin (if any) agree; RC<n>
  when yrec exited with status n (999: it ran out its hour); DIFF
  otherwise. same/diff count the deck's files equal to and different
  from the baseline's; pin= is blank for a deck without a pin. The
  differing files follow, then new: and gone: names. A snapshot prints
  OK (or RC<n>), the time, the number of files taken and the pin.
  Then a pins line (how many pins were compared, and where), and last:
    RESULT: ALL IDENTICAL       (exit status 0)
    RESULT: <n> NOT IDENTICAL   (exit status 1)
  A snapshot ends with RESULT: SNAPSHOT OF <n> DECKS, or RESULT: <k>
  FAILED (exit status 1) when a deck did not run. Exit status 2: the
  build failed, src/yrec is missing (--nobuild), the tree has no
  src/tools/compare_decks, or there is no baseline to compare with.

The sets (DECKS_DIR is src/tools/compare_decks, whose README.md says why
each of its decks is there):
  quick  the configuration matrix (examples/run_config_matrix), cm_accretion
         and the noGS_norot standard solar model: 18 decks, half a minute
  full   quick, the 12 testsuite solar decks and 10 more examples: 40 decks
  rot    the full set's decks whose name contains "rot" (norot ones too), the
         15 rv_ decks and o24_rescale_rotating: 29 decks
  all    full and rot: 56 decks, about 10 minutes with 6 jobs (the longest
         deck, ZAHB_m1p0feh+0p0_TAHB, takes about 5)
--only REGEX keeps the decks of the set whose file name matches.

Where things go:
  The work directory (--work; default .compare_decks/ in the main checkout,
  which .gitignore lists) holds the baseline, base/<deck dir>/output/, with
  SNAPSHOTS.log naming the commit each snapshot came from, and each run's
  build log and yrec terminal output, runs/<LABEL>/ (LABEL defaults to the
  tree's directory name). The pins (--pins; default the main checkout, the
  first worktree `git worktree list` names) are the git-ignored
  <deck dir>/standard/*.log, which only the main checkout has; a tree
  without them gets no pin check (the pins line says so). The tool writes
  nothing in TREE but the build's files and the decks' output/ files.

Never run two of these, or a build, in the same tree at once: the decks
of one directory share its output/ directory.
"""
import argparse
import datetime
import glob
import os
import re
import shutil
import subprocess as sp
import sys
import time
from concurrent.futures import ThreadPoolExecutor

DECKS_DIR = "src/tools/compare_decks"
VER = re.compile(rb"^# YREC v\S* \([0-9a-f]+(-dirty)?\)$", re.M)


def same_file(a, b):
    """byte-equal, the version line aside"""
    x, y = open(a, "rb").read(), open(b, "rb").read()
    return x == y or VER.sub(b"", x) == VER.sub(b"", y)


def deck_name(f1):
    return f1.replace(".inlist", "").replace(".nml1", "")


def deck_sets(tree):
    """the four sets, each a list of (directory, deck file, .nml2 file or '')"""
    cm = "examples/run_config_matrix"
    ssm = "examples/run_standard_solar_model"
    quick = [(cm, f, "") for f in sorted(os.listdir(f"{tree}/{cm}")) if f.endswith(".inlist")]
    quick += [(DECKS_DIR, "cm_accretion.inlist", "")]
    quick += [(ssm, "Test_solar_noGS_norot.nml1", "Test_solar_noGS_norot.nml2")]
    full = quick + [("testsuite", f"Test_solar{s}.inlist", "") for s in (
        "AAG21_dif", "GS98_SF3", "GS98_allard", "GS98_base", "GS98_dif",
        "GS98_difrotmixper18p4d0p1", "GS98_difrotmixper18p4d1p0",
        "GS98_difrotmixsolidper18p4d1p0", "GS98_kurucz", "GS98_yaleeos",
        "M22M_dif", "M22M_difkur")]
    full += [("examples/run_from_dbl_to_zams", "Test_m0030_feh0_noGS_ZAMS.inlist", ""),
             ("examples/run_from_seed_to_start", "m0030gs98z018804_prems.inlist", ""),
             ("examples/run_from_start_to_dbl", "m0030gs98z018804_Dbl.inlist", "")]
    full += [(ssm, f"{c}.nml1", f"{c}.nml2") for c in (
        "Test_solar_GS_norot_gray_OPALSCV_SF2_GS98_OP_CF10",
        "Test_solar_GS_rot_gray_OPALSCV_SF2_GS98_OP_CF10", "Test_solarcal_noGS_norot")]
    full += [("examples/test_nuclear_reaction_rates", "Test_m10000_feh0_noGSstdcno_ZAMS.inlist", ""),
             ("examples/test_nuclear_reaction_rates", "Test_m1000_feh0_noGSstdpp_ZAMS.inlist", ""),
             ("examples/run_from_zahb_to_tahb", "ZAHB_m1p0feh+0p0_TAHB.inlist", ""),
             ("examples/run_rotation_to_RG",
              "Test_3Msolar_GS_rot_gray_OPALSCV_SF2_GS98_OP_CF10_TAMS.inlist", "")]
    variants = sorted(f for f in os.listdir(f"{tree}/{DECKS_DIR}") if f.startswith("rv_"))
    rot = [d for d in full if "rot" in d[1]] + [(DECKS_DIR, f, "") for f in variants]
    rot += [(DECKS_DIR, "o24_rescale_rotating.inlist", "")]
    return {"quick": quick, "full": full, "rot": rot,
            "all": full + [d for d in rot if d not in full]}


def fail(message):
    print(message, file=sys.stderr)
    sys.exit(2)


def git(tree, *args):
    r = sp.run(["git", "-C", tree, *args], capture_output=True, text=True)
    return r.stdout.strip() if r.returncode == 0 else ""


def main_checkout(tree):
    """the main worktree (the first `git worktree list` names), else tree"""
    for line in git(tree, "worktree", "list", "--porcelain").splitlines():
        if line.startswith("worktree "):
            return line[len("worktree "):]
    return tree


def build(tree, rundir):
    t0 = time.time()
    r = sp.run(["make", "-j8"], cwd=f"{tree}/src", capture_output=True, text=True)
    with open(f"{rundir}/build.log", "w") as f:
        f.write(r.stdout + r.stderr)
    errs = [l for l in (r.stdout + r.stderr).splitlines() if "Error" in l]
    print(f"build: rc={r.returncode} {time.time() - t0:.0f}s errors={len(errs)}")
    for l in errs[:10]:
        print("  ", l)
    if r.returncode:
        sys.exit(2)
    r = sp.run([sys.executable, "tools/check_boundaries.py"], cwd=f"{tree}/src",
               capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    print("boundaries:", out.splitlines()[-1] if out else "", "rc", r.returncode)


def run_deck(tree, rundir, deck):
    """runs one deck from its directory; (deck, exit status, seconds)"""
    ddir, f1, f2 = deck
    cwd, base = f"{tree}/{ddir}", deck_name(f1)
    os.makedirs(f"{cwd}/output", exist_ok=True)
    for old in glob.glob(f"{glob.escape(cwd)}/output/{glob.escape(base)}.*"):
        os.remove(old)
    env = dict(os.environ, YREC_INPUT=f"{tree}/input", YREC_START=f"{tree}/startmodels")
    env.pop("YREC_OUTPUT", None)
    args = [f"{tree}/src/yrec", f1] + ([f2] if f2 else [])
    t0 = time.time()
    try:
        r = sp.run(args, cwd=cwd, env=env, capture_output=True, text=True, timeout=3600)
        rc, out = r.returncode, r.stdout + r.stderr
    except sp.TimeoutExpired:
        rc, out = 999, "TIMEOUT"
    with open(f"{rundir}/{base}.out", "w") as f:
        f.write(out)
    return deck, rc, time.time() - t0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0],
                                 epilog="See the module docstring for the rules and output.")
    ap.add_argument("tree", help="the checkout or worktree under test")
    ap.add_argument("label", nargs="?", help="names the run logs' directory (default: the tree's name)")
    ap.add_argument("--set", default="quick", choices=["quick", "full", "rot", "all"])
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--only", default="", help="keep the decks whose file name matches this regex")
    ap.add_argument("--nobuild", action="store_true", help="use the tree's src/yrec as it is")
    ap.add_argument("--snapshot", action="store_true", help="write the baseline instead of comparing")
    ap.add_argument("--work", help="baseline and run logs (default: <main checkout>/.compare_decks)")
    ap.add_argument("--pins", help="root of the standard/ pins (default: the main checkout)")
    a = ap.parse_args()

    tree = os.path.abspath(a.tree)
    if not os.path.isdir(f"{tree}/{DECKS_DIR}"):
        fail(f"{tree} has no {DECKS_DIR} (the extra decks)")
    main_tree = main_checkout(tree)
    work = os.path.abspath(a.work or f"{main_tree}/.compare_decks")
    pins = os.path.abspath(a.pins or main_tree)
    baseline = f"{work}/base"
    rundir = f"{work}/runs/{a.label or os.path.basename(tree)}"
    os.makedirs(rundir, exist_ok=True)

    decks = deck_sets(tree)[a.set]
    if a.only:
        decks = [d for d in decks if re.search(a.only, d[1])]
    names = [deck_name(d[1]) for d in decks]
    assert len(set(names)) == len(names), "two decks of the set share a name"
    print(f"tree {tree}, {len(decks)} decks ({a.set}{', only ' + a.only if a.only else ''}), "
          f"baseline {baseline}, pins {pins}")
    if not a.snapshot:
        if not os.path.isdir(baseline):
            fail(f"no baseline at {baseline}: take one with --snapshot first")
        log = f"{baseline}/SNAPSHOTS.log"
        if os.path.exists(log):
            print("baseline:", open(log).read().splitlines()[-1])

    if not a.nobuild:
        build(tree, rundir)
    if not os.path.exists(f"{tree}/src/yrec"):
        fail(f"{tree}/src/yrec is not built")
    t0 = time.time()
    with ThreadPoolExecutor(a.jobs) as ex:
        results = list(ex.map(lambda d: run_deck(tree, rundir, d), decks))
    print(f"ran {len(decks)} decks in {time.time() - t0:.0f}s")

    bad, pin_dirs, pin_bad = 0, {}, 0
    for (ddir, f1, f2), rc, dt in results:
        base = deck_name(f1)
        outdir, refdir = f"{tree}/{ddir}/output", f"{baseline}/{ddir}/output"
        outs = sorted(glob.glob(f"{glob.escape(outdir)}/{glob.escape(base)}.*"))
        # the run log against the pin, where there is one
        pinfile, log = f"{pins}/{ddir}/standard/{base}.log", f"{outdir}/{base}.log"
        pin = ""
        if os.path.exists(pinfile):
            pin = "pin=OK" if os.path.exists(log) and same_file(log, pinfile) else "pin=DIFF"
            pin_dirs[ddir] = pin_dirs.get(ddir, 0) + 1
            pin_bad += pin == "pin=DIFF"
        if a.snapshot:
            os.makedirs(refdir, exist_ok=True)
            for old in glob.glob(f"{glob.escape(refdir)}/{glob.escape(base)}.*"):
                os.remove(old)
            for o in outs:
                shutil.copy(o, refdir)
            bad += rc != 0
            print(f"{'OK ' if rc == 0 else 'RC' + str(rc)} {dt:6.0f}s snap {len(outs):2d} files  "
                  f"{ddir}/{base} {pin}".rstrip())
            continue
        # every output file of the deck against the baseline's
        diffs, missing, same = [], [], 0
        for o in outs:
            name = os.path.basename(o)
            ref = f"{refdir}/{name}"
            if not os.path.exists(ref):
                missing.append(name)
            elif same_file(o, ref):
                same += 1
            else:
                diffs.append(name)
        # and every baseline file of the deck must have been written again
        gone = [os.path.basename(p)
                for p in sorted(glob.glob(f"{glob.escape(refdir)}/{glob.escape(base)}.*"))
                if not os.path.exists(f"{outdir}/{os.path.basename(p)}")]
        ok = rc == 0 and outs and not diffs and not missing and not gone and pin != "pin=DIFF"
        bad += not ok
        tag = "IDENT" if ok else (f"RC{rc}" if rc else "DIFF")
        notes = diffs + [f"new:{m}" for m in missing] + [f"gone:{g}" for g in gone]
        print(f"{tag:8s} {dt:4.0f}s same={same:2d} diff={len(diffs)} {pin:8s} {ddir}/{base} "
              f"{' '.join(notes)}".rstrip())

    where = ", ".join(f"{d} {n}" for d, n in sorted(pin_dirs.items()))
    n_pins = sum(pin_dirs.values())
    print(f"pins: {n_pins - pin_bad} OK, {pin_bad} DIFF ({where}) under {pins}" if n_pins
          else f"pins: none under {pins} for these decks")
    if a.snapshot:
        head = git(tree, "rev-parse", "--short=8", "HEAD") or "?"
        dirty = "-dirty" if git(tree, "status", "--porcelain", "--untracked-files=no") else ""
        stamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
        with open(f"{baseline}/SNAPSHOTS.log", "a") as f:
            f.write(f"{stamp} {head}{dirty} {tree} --set {a.set}"
                    f"{' --only ' + a.only if a.only else ''}: {len(decks)} decks, {bad} failed\n")
        print(f"RESULT: SNAPSHOT OF {len(decks)} DECKS" if bad == 0 else f"RESULT: {bad} FAILED")
    else:
        print("RESULT: ALL IDENTICAL" if bad == 0 else f"RESULT: {bad} NOT IDENTICAL")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()

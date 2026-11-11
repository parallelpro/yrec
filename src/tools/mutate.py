#!/usr/bin/env python3
"""Mutation check of a package's test programs: each mutation of the code
must make the check that covers it print FAIL.

    python3 src/tools/mutate.py <pkg> [--program NAME] [--only REGEX] [--check]

from the repository root, after `make tests` in src/. <pkg> is a directory
of src/ (math, numerics, core, ...); its manifests are
src/<pkg>/test/<program>.mutations.toml, one per test program (--program
picks one). A manifest is TOML (tomllib, Python 3.11 or later):

    program = "test_numerics"     # make target in src/, and the binary src/<program>
    timeout = 120                 # seconds (optional, default 120)
    args = []                     # (optional) the program's arguments; {repo} is the
                                  # repository root
    [env]                         # (optional) variables for the run
    [backends]                    # (optional) backend name = the program built and run
    crmath = "test_math_crmath"   #   for that backend's mutations

    [[mutation]]
    id = "cubic_spline-01"        # unique in the manifest
    file = "numerics/cubic_spline.f90"   # under src/
    old = '...'                   # text that occurs exactly once in the file
    new = '...'                   # its replacement
    check = "cspline/splint: clamped spline exact on a cubic"
    backend = "crmath"            # (optional) a key of [backends]

For each mutation it writes the mutant, runs `make <program>` in src/,
runs src/<program> in a fresh temporary directory (holding an empty
output/; YREC_INPUT and YREC_START are the repository's unless set), and
reads the lines of its standard output that start with `check`. The
mutation is caught when one of those lines ends in FAIL; when no such
line is printed and the program exits nonzero (a crash, an error stop);
or when the program runs longer than the timeout. Otherwise it survived.
A mutant that does not compile is reported as BUILD FAILED, and counts
as not caught. A mutation of a backend whose toolchain is missing (crmath
without the MESA SDK) is skipped. --only runs the mutations whose id or
check matches REGEX. --check only verifies the manifests: required
fields, unique ids, known backends, and each `old` occurring exactly once
in its file (run it after editing a mutated routine).

Every file is restored after its mutation, whatever happens (an error, a
timeout, Ctrl-C or SIGTERM), and the programs are rebuilt from the
restored sources at the end. make compares 1 s file times, so the tool
waits a second before writing a mutant and before restoring a file, and
touches the restored file. Exit status 1 if a mutation survived or did
not build, or (--check) if a manifest is wrong.
"""
import argparse
import os
import pathlib
import re
import signal
import subprocess
import sys
import tempfile
import time

try:
    import tomllib
except ImportError:  # Python < 3.11
    sys.exit("mutate.py needs Python 3.11 or later (tomllib)")

REPO = pathlib.Path(__file__).resolve().parents[2]
SRC = REPO / "src"
REQUIRED = ("id", "file", "old", "new", "check")


def crmath_available():
    """The Makefile's rule: MESASDK_DIR, else MESASDK_ROOT, else
    /Applications/mesasdk, holding include/crmath.mod."""
    sdk = os.environ.get("MESASDK_DIR") or os.environ.get("MESASDK_ROOT") or "/Applications/mesasdk"
    return (pathlib.Path(sdk) / "include" / "crmath.mod").exists()


# A backend whose toolchain may be missing: its mutations are skipped then.
BACKEND_AVAILABLE = {"crmath": crmath_available}


def manifests(pkg, program=None):
    paths = sorted((SRC / pkg / "test").glob("*.mutations.toml"))
    if program:
        paths = [p for p in paths if p.name == f"{program}.mutations.toml"]
    if not paths:
        sys.exit(f"mutate.py: no src/{pkg}/test/{program or '*'}.mutations.toml")
    return [(p, tomllib.loads(p.read_text())) for p in paths]


def problems_of(path, m):
    """What is wrong with one manifest, as messages."""
    out = []
    rel = path.relative_to(REPO).as_posix()
    if not isinstance(m.get("program"), str):
        out.append(f"{rel}: no `program`")
    backends = m.get("backends", {})
    seen = set()
    for k, mut in enumerate(m.get("mutation", []), 1):
        missing = [f for f in REQUIRED if not isinstance(mut.get(f), str)]
        if missing:
            out.append(f"{rel}: mutation {k}: missing {', '.join(missing)}")
            continue
        if mut["id"] in seen:
            out.append(f"{rel}: id {mut['id']} twice")
        seen.add(mut["id"])
        if "backend" in mut and mut["backend"] not in backends:
            out.append(f"{rel}: {mut['id']}: backend {mut['backend']} not in [backends]")
        f = SRC / mut["file"]
        n = f.read_text().count(mut["old"]) if f.is_file() else -1
        if n != 1:
            out.append(f"MISSING  {mut['id']}  {mut['check']}  ({mut['file']}: "
                       + ("no such file" if n < 0 else f"the old text occurs {n} times") + ")")
    if not m.get("mutation"):
        out.append(f"{rel}: no [[mutation]]")
    return out


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def run_program(m, binary):
    """(stdout, return code), or None after the timeout."""
    env = dict(os.environ)
    env.setdefault("YREC_INPUT", str(REPO / "input"))
    env.setdefault("YREC_START", str(REPO / "startmodels"))
    env.update({k: str(v) for k, v in m.get("env", {}).items()})
    args = [a.replace("{repo}", str(REPO)) for a in m.get("args", [])]
    with tempfile.TemporaryDirectory(prefix="mutate_") as work:
        (pathlib.Path(work) / "output").mkdir()
        try:
            r = run([str(SRC / binary)] + args, cwd=work, env=env, timeout=m.get("timeout", 120))
        except subprocess.TimeoutExpired:
            return None
    return r.stdout, r.returncode


def caught(check, stdout, returncode):
    lines = [l for l in stdout.splitlines() if l.startswith(check)]
    return any(l.rstrip().endswith("FAIL") for l in lines) or (not lines and returncode != 0)


def interrupt(signum, frame):
    raise KeyboardInterrupt(f"signal {signum}")


def mutate(m, mut, binary):
    """'caught', 'caught (timeout)', 'SURVIVED' or 'BUILD FAILED', with
    the file restored afterwards."""
    f = SRC / mut["file"]
    original = f.read_text()
    assert original.count(mut["old"]) == 1, (mut["id"], mut["file"])
    try:
        time.sleep(1)
        f.write_text(original.replace(mut["old"], mut["new"]))
        b = run(["make", binary], cwd=SRC)
        if b.returncode != 0:
            return "BUILD FAILED", (b.stdout[-800:] + b.stderr[-800:])
        result = run_program(m, binary)
        if result is None:
            return "caught", f"the mutant did not finish in {m.get('timeout', 120)} s"
        return ("caught" if caught(mut["check"], *result) else "SURVIVED"), ""
    finally:
        f.write_text(original)
        time.sleep(1)
        f.touch()


def main():
    ap = argparse.ArgumentParser(description="mutation check (see the module docstring)")
    ap.add_argument("pkg", help="a directory of src/ with test/*.mutations.toml")
    ap.add_argument("--program", help="one test program's manifest")
    ap.add_argument("--only", default="", help="the mutations whose id or check matches REGEX")
    ap.add_argument("--check", action="store_true", help="only verify the manifests")
    a = ap.parse_args()
    found = manifests(a.pkg, a.program)
    if a.check:
        bad = total = 0
        for path, m in found:
            problems = problems_of(path, m)
            total += len(m.get("mutation", []))
            bad += len(problems)
            for p in problems:
                print(p)
        print(f"{total - bad} of {total} mutations apply")
        return 1 if bad else 0
    for path, m in found:
        problems = problems_of(path, m)
        if problems:
            print("\n".join(problems))
            return 1
    for sig in (signal.SIGTERM, signal.SIGHUP):
        signal.signal(sig, interrupt)
    failed = skipped = 0
    built = set()
    try:
        for path, m in found:
            for mut in m["mutation"]:
                if a.only and not (re.search(a.only, mut["id"]) or re.search(a.only, mut["check"])):
                    continue
                backend = mut.get("backend")
                binary = m.get("backends", {}).get(backend, m["program"])
                where = mut["file"] + (f", {backend}" if backend else "")
                if backend in BACKEND_AVAILABLE and not BACKEND_AVAILABLE[backend]():
                    print(f"skipped   {mut['id']}  {mut['check']}  ({where}: no toolchain)")
                    skipped += 1
                    continue
                built.add(binary)
                verdict, detail = mutate(m, mut, binary)
                failed += verdict != "caught"
                print(f"{verdict:<9} {mut['id']}  {mut['check']}  ({where})"
                      + (f"\n{detail}" if detail and verdict != "caught" else
                         f"  [{detail}]" if detail else ""), flush=True)
    finally:
        for binary in sorted(built):
            r = run(["make", binary], cwd=SRC)
            if r.returncode != 0:
                print(f"rebuilding {binary} from the restored sources failed:\n{r.stderr[-800:]}")
    print(f"{failed} mutation(s) survived or did not build, {skipped} skipped")
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except KeyboardInterrupt as e:
        print(f"mutate.py: interrupted ({e or 'Ctrl-C'}); the mutated file is restored")
        sys.exit(130)

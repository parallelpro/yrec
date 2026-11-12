#!/usr/bin/env python3
"""Run the unit test programs built with the development flags.

    make -C src -j8 devel_tests && python3 src/tools/devel_tests.py [--all] [--only REGEX]

`make devel_tests` builds the test programs with DEV_FFLAGS (gfortran:
-fcheck=all -ffpe-trap=invalid,zero,overflow, and -g for the line of a
trap; src/Makefile's build modes). This runs each program the way its
root test runs it (a fresh temporary directory, YREC_INPUT and YREC_START
the repository's), with its arguments for each case. A case passes when
the program exits 0, prints no gfortran runtime error ("Fortran runtime
error", "Program received signal", "Error termination"), and prints no
line whose last word is FAIL. Nothing is byte-compared: the devel flags
move round-off, so the pinned outputs are the release build's (`make
tests` and pytest). Runtime warnings (array temporaries, which
-fcheck=all reports) are counted, not failures.

On macOS the program runs as a copy with its dSYM (dsymutil), so that a
trap's backtrace names the file and line, as it does on Linux.

A case in KNOWN_FAILURES is skipped and listed with its defect until the
defect is fixed (--all runs it too). Exit status 1 if a case fails.
CI's devel job runs this.
"""
import argparse
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile

REPO = pathlib.Path(__file__).resolve().parents[2]
SRC = REPO / "src"
SOLAR = "examples/run_standard_solar_model/Test_solar_noGS_norot"

# (case, program, arguments, extra environment, working directory or None
# for the temporary one). {repo} is the repository root, {work} the
# temporary directory.
CASES = [
    ("test_math", "test_math", [], {}, None),
    ("test_math_crmath", "test_math_crmath", [], {}, None),
    ("test_numerics", "test_numerics", [], {}, None),
    ("test_eos", "test_eos", [], {"YREC_MHD_TABLES": None}, None),
    ("test_kap", "test_kap", [], {}, None),
    ("test_atm", "test_atm", [], {}, None),
    ("test_net", "test_net", [], {}, None),
    ("test_convection", "test_convection", [], {}, None),
    ("test_element_diffusion", "test_element_diffusion", [], {}, None),
    ("test_rotation", "test_rotation", [], {}, None),
    ("test_rotation_kernels", "test_rotation_kernels", [], {}, None),
    ("test_triangle_fit", "test_triangle_fit", [], {}, None),
    ("test_mass_change", "test_mass_change", [], {}, None),
    ("test_regrid", "test_regrid", [], {}, None),
    ("test_structure[tams]", "test_structure",
     ["{repo}/src/core/test/test_structure_tams.inlist", "3e-3"], {"YREC_OUTPUT": "output"}, None),
    ("test_structure[yale]", "test_structure",
     ["{repo}/testsuite/Test_solarGS98_yaleeos.inlist", "1e-4"], {"YREC_OUTPUT": "output"}, None),
    ("test_structure[extra]", "test_structure",
     ["{repo}/src/core/test/test_structure_tams.inlist", "3e-3", "extra"],
     {"YREC_OUTPUT": "output"}, None),
    ("test_io", "test_io", [], {}, None),
    ("test_controls", "test_controls", [], {"YREC_OUTPUT": "{work}"}, None),
    ("test_controls[solar]", "test_controls", [f"{{repo}}/{SOLAR}.nml1", f"{{repo}}/{SOLAR}.nml2"],
     {"YREC_OUTPUT": "{work}"}, "{repo}/examples/run_standard_solar_model"),
]
# the configuration-matrix decks through the controls read, from their
# own directory as test_controls.py runs them
for _deck in sorted((REPO / "examples" / "run_config_matrix").glob("*.inlist")):
    CASES.append((f"test_controls[{_deck.stem}]", "test_controls", [str(_deck)],
                  {"YREC_OUTPUT": "{work}"}, str(_deck.parent)))

# Cases that fail under the devel flags today, each with its defect (the
# first the run meets; the program stops there). Remove an entry when its
# defect is fixed, and see what the run meets next.
_O4 = ("O4 (docs/bug_sweep.md): the starting model's stitched atmosphere reads "
       "env_log10_radius(0), core/diagnostics/stitched_model.f90:475, a bounds error")
KNOWN_FAILURES = {
    "test_eos": "the fully ionized density iteration's second-order step overflows far from "
                "its root, eos/yale/yale_eos.f90:925 (iterate_fully_ionized_density), an "
                "overflow trap at test_eos's check of a start 100 dex off (M93's)",
    "test_structure[tams]": _O4,
    "test_structure[yale]": _O4,
    "test_structure[extra]": _O4,
}

RUNTIME_ERROR = re.compile(r"Fortran runtime error|Program received signal|Error termination")


def run_case(case, program, args, env_extra, cwd, timeout):
    """(passed, summary lines)"""
    binary = SRC / program
    if not binary.exists():
        return None, [f"{program} not built"]
    env = dict(os.environ)
    env.setdefault("YREC_INPUT", str(REPO / "input"))
    env.setdefault("YREC_START", str(REPO / "startmodels"))
    with tempfile.TemporaryDirectory(prefix="devel_") as work:
        (pathlib.Path(work) / "output").mkdir()
        for k, v in env_extra.items():
            if v is None:
                env.pop(k, None)
            else:
                env[k] = v.replace("{work}", work)
        if sys.platform == "darwin" and shutil.which("dsymutil"):
            # macOS keeps the debug information in the objects: a copy of
            # the program with its dSYM beside it lets a trap's backtrace
            # name the file and line
            binary = pathlib.Path(work) / program
            shutil.copy2(SRC / program, binary)
            subprocess.run(["dsymutil", str(binary)], capture_output=True)
        argv = [str(binary)] + [a.replace("{repo}", str(REPO)) for a in args]
        where = (cwd or "{work}").replace("{repo}", str(REPO)).replace("{work}", work)
        try:
            r = subprocess.run(argv, cwd=where, env=env, capture_output=True, text=True,
                               errors="replace", timeout=timeout)
        except subprocess.TimeoutExpired:
            return False, [f"no end within {timeout} s"]
    text = r.stdout + "\n" + r.stderr
    errors = [l for l in text.splitlines() if RUNTIME_ERROR.search(l)]
    fails = [l for l in r.stdout.splitlines() if l.split() and l.split()[-1] == "FAIL"]
    warnings = text.count("Fortran runtime warning")
    out = []
    if r.returncode != 0:
        out.append(f"exit status {r.returncode}")
    if errors:
        lines = text.splitlines()
        k = lines.index(errors[0])
        out += [l for l in lines[max(0, k - 2):k + 30]
                if l.strip() and "in ???" not in l and not l.startswith("Fortran runtime warning")][:16]
    out += [f"FAIL line: {l.strip()}" for l in fails[:10]]
    if warnings:
        out.append(f"({warnings} runtime warnings)")
    return r.returncode == 0 and not errors and not fails, out


def main():
    ap = argparse.ArgumentParser(description="run the devel-built unit tests (see the docstring)")
    ap.add_argument("--all", action="store_true", help="run the known failures too")
    ap.add_argument("--only", default="", help="the cases whose name matches REGEX")
    ap.add_argument("--timeout", type=int, default=900, help="seconds per case")
    a = ap.parse_args()
    failed = 0
    for case, program, args, env_extra, cwd in CASES:
        if a.only and not re.search(a.only, case):
            continue
        if case in KNOWN_FAILURES and not a.all:
            print(f"skipped  {case}: known failure, {KNOWN_FAILURES[case]}", flush=True)
            continue
        ok, lines = run_case(case, program, args, env_extra, cwd, a.timeout)
        if ok is None:
            if program == "test_math_crmath":
                print(f"skipped  {case}: not built (no MESA SDK)", flush=True)
                continue
            ok = False
        failed += not ok
        print(f"{'ok      ' if ok else 'FAILED  '} {case}", flush=True)
        for l in lines:
            if not ok or l.startswith("("):
                print(f"    {l}")
    print(f"{failed} case(s) failed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())

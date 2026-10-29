#!/usr/bin/env python3
"""Mutation check of test_math: each check must fail under a mutation of
the math_lib code it covers, under the backend that code belongs to.

    python3 src/math/test/mutate_math.py [--only REGEX] [--check]

from the repository root, after `make tests` in src/. For each mutation
it edits math/math_lib.f90, rebuilds the test program of the mutated
backend (test_math for the intrinsic branch, test_math_crmath for the
crmath branch, which needs the MESA SDK: without it those mutations are
reported as skipped), runs it and expects the named check to print FAIL
(a crash, a nonzero exit without the check's line, or a run longer than
120 s also counts as caught). The file is restored afterwards, whatever
happens; make uses 1 s mtimes, so the script sleeps a second around each
write. Exit status 1 if any mutation survives. --check only verifies
that every mutation's text occurs exactly once in math_lib.f90.
"""
import argparse
import os
import pathlib
import re
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parents[3]
SRC = REPO / "src"
MATH_LIB = SRC / "math" / "math_lib.f90"
BINARY = {"intrinsic": "test_math", "crmath": "test_math_crmath"}

# (backend, the text replaced (it must occur exactly once), the
# replacement, the label of the check that must fail)
MUTATIONS = [
    ("crmath", "            pow_x = pow_x*x\n         end do\n         if (iy < 0) pow_x = 1.0d0/pow_x\n      else if",
     "            pow_x = pow_x*x*(1.0d0+epsilon(1.0d0))\n         end do\n         if (iy < 0) pow_x = 1.0d0/pow_x\n      else if",
     "pow: integer exponents are repeated multiplication"),
    ("crmath", "      if (iy < 0) pow_x = 1.0d0/pow_x\nend function pow_i",
     "      if (iy < 0) pow_x = 1.0d0/(pow_x*(1.0d0+epsilon(1.0d0)))\nend function pow_i",
     "pow: integer exponents are repeated multiplication"),
    ("crmath", "            exp10_x = exp10_x*10.0d0\n",
     "            exp10_x = exp10_x*10.0d0*(1.0d0+epsilon(1.0d0))\n",
     "exp10: exact for 0 <= n <= 22"),
    ("crmath", "      else\n         pow_x = exp(log(x)*y)\n      end if",
     "      else\n         pow_x = exp(log(x)*y*(1.0d0+1.0d-12))\n      end if",
     "pow, exp10: x**y and 10**x to a few ulp"),
    ("crmath", "      pow_x = pow_r(x, real(y, dp))",
     "      pow_x = pow_r(x, real(nint(y*1.0e7), dp)*1.0d-7)",
     "pow: a single-precision exponent is its binary value"),
    ("crmath", "         if (mod(y, 2.0d0) /= 0.0d0) pow_x = -pow_x",
     "         if (.false.) pow_x = -pow_x",
     "pow: (-2)**101 and (-2)**100"),
    ("crmath", "      if (x == 0.0d0) then\n         pow_x = 0.0d0\n         return\n      end if\n! the test on |y|",
     "! the test on |y|",
     "pow: 0**0 and 0**(-1) are 0 under crmath"),
    ("crmath", "      if (abs(y) < 100.0d0 .and. y == aint(y)) then",
     "      if (abs(y) < 1.0d10 .and. y == aint(y)) then",
     "pow: y = 3e9 takes the composition"),
    ("crmath", "         do i = 1, abs(iy)\n            pow_x = pow_x*x\n         end do\n         if (iy < 0) pow_x = 1.0d0/pow_x\n      else if",
     "         do i = 1, abs(iy)\n            pow_x = pow_x*abs(x)\n         end do\n         if (iy < 0) pow_x = 1.0d0/pow_x\n      else if",
     "pow: (-2)**3 = -8"),
    ("crmath", "      else\n         pow_x = exp(log(x)*y)",
     "      else\n         pow_x = exp(log(abs(x))*y)",
     "pow: (-2)**0.5 is NaN"),
    ("crmath", "      if (x == 0.0d0) then\n         pow_x = 0.0d0\n         return\n      end if\n! the test on |y|",
     "      if (x == 0.0d0) then\n         pow_x = 1.0d-300\n         return\n      end if\n! the test on |y|",
     "pow: 0**y = 0 for y > 0"),
    ("intrinsic", "      pow_x = x**y\nend function pow_r\n",
     "      pow_x = merge(1.0d-300, x**y, x == 0.0d0)\nend function pow_r\n",
     "pow: 0**y = 0 for y > 0"),
    ("intrinsic", "      pow_x = x**y\nend function pow_r\n",
     "      pow_x = merge(0.0d0, x**y, x == 0.0d0)\nend function pow_r\n",
     "pow: 0**0 = 1 and 0**(-1) = +Inf under the intrinsic backend"),
    ("intrinsic", "      pow_x = x**y\nend function pow_r\n",
     "      pow_x = x**y*(1.0d0+epsilon(1.0d0))\nend function pow_r\n",
     "pow, exp10: the intrinsic backend is x**y, bit for bit"),
    ("intrinsic", "      exp10_x = 10.0d0**x\n",
     "      exp10_x = 10.0d0**(x*(1.0d0+1.0d-12))\n",
     "pow, exp10: x**y and 10**x to a few ulp"),
    ("intrinsic", "      pow_x = x**iy\n",
     "      pow_x = x**iy*(1.0d0+epsilon(1.0d0))\n",
     "pow, exp10: integer exponents are x**y, x**n, 10**x"),
    ("intrinsic", "      pow_x = x**y\nend function pow_r_sp",
     "      pow_x = x**(real(nint(y*1.0e7), dp)*1.0d-7)\nend function pow_r_sp",
     "pow: a single-precision exponent is its binary value"),
]


def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)


def crmath_available():
    """The Makefile's rule: MESASDK_DIR, else MESASDK_ROOT, else
    /Applications/mesasdk, holding include/crmath.mod."""
    sdk = os.environ.get("MESASDK_DIR") or os.environ.get("MESASDK_ROOT") or "/Applications/mesasdk"
    return (pathlib.Path(sdk) / "include" / "crmath.mod").exists()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default="")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    original = MATH_LIB.read_text()
    if a.check:
        missing = [label for backend, old, new, label in MUTATIONS if original.count(old) != 1]
        for label in missing:
            print(f"MISSING  {label}")
        print(f"{len(MUTATIONS) - len(missing)} of {len(MUTATIONS)} mutations apply")
        return 1 if missing else 0
    survived = skipped = 0
    for backend, old, new, label in MUTATIONS:
        if a.only and not re.search(a.only, label):
            continue
        if backend == "crmath" and not crmath_available():
            print(f"skipped   {label}  (crmath: no MESA SDK)")
            skipped += 1
            continue
        assert original.count(old) == 1, (backend, old)
        binary = BINARY[backend]
        try:
            time.sleep(1)
            MATH_LIB.write_text(original.replace(old, new))
            b = run(["make", binary], cwd=SRC)
            if b.returncode != 0:
                print(f"BUILD FAILED  {label}\n{b.stdout[-800:]}{b.stderr[-800:]}")
                survived += 1
                continue
            try:
                r = run([str(SRC / binary)], cwd=SRC, timeout=120)
            except subprocess.TimeoutExpired:
                print(f"caught    {label}  ({backend}; the mutant did not finish in 120 s)")
                continue
            lines = [l for l in r.stdout.splitlines() if l.startswith(label)]
            caught = any(l.rstrip().endswith("FAIL") for l in lines) or \
                (not lines and r.returncode != 0)
            print(f"{'caught  ' if caught else 'SURVIVED'}  {label}  ({backend})")
            survived += 0 if caught else 1
        finally:
            MATH_LIB.write_text(original)
            time.sleep(1)
            MATH_LIB.touch()
    for binary in BINARY.values():
        if binary == "test_math_crmath" and not crmath_available():
            continue
        run(["make", binary], cwd=SRC)
    print(f"{survived} mutation(s) survived, {skipped} skipped")
    return 1 if survived else 0


if __name__ == "__main__":
    sys.exit(main())

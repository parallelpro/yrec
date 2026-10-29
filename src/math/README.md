# math/ — elementary math, for reproducibility

`math_lib` is the one place YREC gets its elementary functions (`exp`, `log`, `log10`, `sin`, ..., `pow`, `exp10`). It sits below every other directory (it uses no YREC module), and its concern is reproducibility: with the correctly-rounded backend, a run gives the same bits on every platform. It is separate from `numerics/` because it is a dependency root with its own usage rule, which `tools/check_boundaries.py` enforces on every file.

## 1. What lives where

| Module (file) | What it holds |
|---|---|
| `math_lib` (`math_lib.f90`) | `pow` (generic: real, integer and single-precision exponents), `exp10`, and `math_backend` (`'crmath'` or `'intrinsic'`). Under the crmath backend it also re-exports crmath's correctly-rounded `exp`, `expm1`, `log`, `log1p`, `log10`, `log2`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh` and `hypot`, which shadow the intrinsics. |

## 2. The two backends

The backend is chosen at build time:

- `make` builds the **intrinsic** backend: `pow(x, y)` is `x**y`, `pow(x, n)` is `x**n` and `exp10(x)` is `10.0d0**x`, verbatim, and the elementary functions are the compiler's. This is the historical behaviour bit for bit, so converting a call site to `pow` or `exp10` leaves results unchanged. The last bits depend on the platform's libm.
- `make USE_CRMATH=1` builds the **crmath** backend (Townsend's Fortran wrapper over crlibm, from the MESA SDK): every elementary function returns the correctly rounded result, so with `-ffp-contract=off` (always set in the Makefile) and the same gfortran major version, results are the same on every platform. `pow` and `exp10` follow MESA's definitions: an integer-valued exponent with |y| < 100 takes repeated multiplication, anything else the composition `exp(log(x)*y)` (`exp(x*ln10)`), which crmath determines completely.

Where the two backends differ (each is checked by `test_math`):

| Case | intrinsic | crmath |
|---|---|---|
| `pow(0, y)`, y > 0 | 0 | 0 |
| `pow(0, 0)` | 1 | 0 (MESA's definition) |
| `pow(0, y)`, y < 0 | +Inf | 0 (MESA's definition) |
| `pow(x, n)`: n an integer-valued real with \|n\| < 100, or an integer | libm's `pow` (real n) or the compiler's power by squaring (integer n) | repeated multiplication |
| `pow(x, y)`, general | libm's `pow` | `exp(log(x)*y)`, within a few ulp of y ln x |
| `pow(x, y)`, x < 0, integer-valued \|y\| ≥ 100 | libm's `pow` | `exp(log(-x)*y)` with the sign of (-1)^y |
| `exp10(n)`, integer n | libm's `pow(10, n)` | exact for 0 ≤ n ≤ 22, 1/10^\|n\| (correctly rounded) for -22 ≤ n < 0, repeated multiplication up to \|n\| < 100 |

The `pow(0, y ≤ 0)` difference is MESA's and kept (2026-10, a revisable default): no YREC call site is known to reach it.

## 3. The usage rule

`tools/check_boundaries.py` (`check_math_lib`) enforces, for every source file outside `math/` and the test programs:

- a procedure that calls a function math_lib shadows under crmath (`exp`, `log`, `log10`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`, `atan2`, `sinh`, `cosh`, `tanh`, `hypot`, and `expm1`, `log1p`, `log2`: the checker reads them from math_lib's `public ::` lines under `#ifdef YREC_CRMATH`) has `use math_lib` in scope, its own or its host's;
- no real exponent appears after `**`: write `pow(x, y)`, and `exp10(x)` for 10 to a real power. Integer exponents (`x**2`, `x**n`) are exact repeated multiplication and stay as `**`. Parameter initializers are exempt (the compiler folds them, correctly rounded).

`use math_lib` is the one `use` without an `only:` list (docs/GUIDELINES.md, "Declarations and interfaces", rule 2): under crmath its procedures shadow the intrinsics, which an `only:` list cannot express for a file that must also build without them.

## 4. Adding a shadowed function

To route another elementary function through math_lib (say `erf`):

1. Check that crmath provides it (`$MESASDK_DIR/include/crmath.mod`; the wrapper's source lists its functions).
2. Add it to the `public ::` lines under `#ifdef YREC_CRMATH` in `math_lib.f90`. Under the intrinsic backend nothing is needed: the intrinsic is used.
3. Run `python3 tools/check_boundaries.py` (in `src/`): it reads the shadowed names from those `public ::` lines, so its callers now must `use math_lib`; add it where the checker asks.
4. Add a check to `math/test/test_math.f90`, under both backends.

A function crmath lacks would need a definition in both branches of the `#ifdef`, as `pow` and `exp10` have.

## 5. Its test

`test_math` (`make tests` in `src/`, then `pytest test_math.py`) checks math_lib under the build's backend, and `test_math_crmath` under crmath whenever the MESA SDK is installed (the Makefile builds that copy of math_lib in `math/test/crmath/`, apart from the main build). The checks: the integer paths (repeated multiplication and exact powers of ten under crmath; `x**y`, `x**n`, `10**x` verbatim under the intrinsic backend), the real paths against `x**y` within the composition's error, a single-precision exponent, negative bases (`(-2)**3`, `(-2)**101`, `(-2)**100`, `(-2)**0.5`), `pow(0, y)` under each backend, and an exponent beyond the default integer range. Each prints ok or FAIL, and the program stops with status 1 on a FAIL. The negative-base check fails on the math_lib of 2c88d896, whose crmath `pow` gave NaN for `(-2)**101`.

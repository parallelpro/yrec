"""pyyrec -- run the YREC stellar-evolution engine in-process.

The 2026 modernization made the engine an embeddable, re-enterable
library (phase five: run_yrec + yrec_reset; see src/ROADMAP.md). This
package is the thin ctypes binding over `libyrec` (built by
`make lib` in src/): each `run()` call executes one full YREC job --
exactly what `yrec <nml1> <nml2>` does from a working directory --
and repeated calls in one Python process are byte-identical to fresh
processes (enforced by test_pyyrec.py).

Typical use, mirroring the examples/ layout::

    import pyyrec
    pyyrec.run("Test_solar_noGS_norot.nml1",
               "Test_solar_noGS_norot.nml2",
               cwd="examples/run_standard_solar_model",
               env={"YREC_INPUT": ".../input",
                    "YREC_START": ".../startmodels"})

Outputs land in the case directory's output/ tree as usual; read them
back with your own tooling (nothing is returned in-memory yet -- the
named-index structure API is a later milestone).

Caveats (documented residuals of the modernization, same as the CLI):
- Configurations that end via a legacy numerics-gate `stop` (e.g. the
  m0030 case's historical BSSTEP termination) end the *process*, not
  the call. Avoid embedding such cases until the numerics-gate ierr
  opt-in lands.
- The engine is one-instance-per-process (module state); use
  `multiprocessing` for parallel parameter scans, one engine per
  worker process.
"""
import ctypes
import os
import pathlib

__all__ = ["run", "YrecError", "library_path"]

_REPO = pathlib.Path(__file__).resolve().parent.parent
_LIB = None


class YrecError(RuntimeError):
    """The engine reported a nonzero ierr for this job."""


def library_path():
    """Path of the shared library run() will load (or raise if unbuilt)."""
    for name in ("libyrec.dylib", "libyrec.so"):
        p = _REPO / "src" / name
        if p.exists():
            return p
    raise OSError(
        f"libyrec not found under {_REPO / 'src'} -- build it with "
        "`make lib` in src/"
    )


def _lib():
    global _LIB
    if _LIB is None:
        lib = ctypes.CDLL(str(library_path()))
        lib.yrec_run.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        lib.yrec_run.restype = ctypes.c_int
        _LIB = lib
    return _LIB


def run(control_nml="", physics_nml="", cwd=None, env=None, check=True):
    """Run one YREC job in-process.

    control_nml / physics_nml: the .nml1/.nml2 paths, relative to cwd
    (empty strings mean the engine's historical defaults,
    yrec8.nml1/yrec8.nml2). cwd: working directory for the job (where
    the namelists live and output/ is written); default is the current
    directory. env: extra environment entries (e.g. YREC_INPUT /
    YREC_START placeholder roots used by the example namelists),
    applied for the duration of the call. check: raise YrecError on
    nonzero ierr instead of returning it.

    Returns the engine's ierr (0 on success).
    """
    lib = _lib()
    old_cwd = os.getcwd()
    saved = {}
    try:
        if env:
            for key, value in env.items():
                saved[key] = os.environ.get(key)
                os.environ[key] = str(value)
        if cwd is not None:
            os.chdir(cwd)
        ierr = lib.yrec_run(
            str(control_nml).encode(), str(physics_nml).encode()
        )
    finally:
        os.chdir(old_cwd)
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
    if check and ierr != 0:
        raise YrecError(f"yrec_run returned ierr={ierr}")
    return ierr

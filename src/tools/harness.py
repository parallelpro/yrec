"""What the test harnesses share: where a run reads the tables and the
starting models.

A test run tests its own tree. pytest (conftest.py), devel_tests.py and
mutate.py take `input/` and `startmodels/` from the repository they are
in, whatever YREC_INPUT and YREC_START name in the shell: a shell that
exports the two for one checkout would otherwise make every other
worktree read that checkout's tables, without a word.

The one opt-in is the environment variable

    YREC_TEST_EXTERNAL_ROOTS=1

With it YREC_INPUT and YREC_START are taken from the environment (each
one that is set and not empty; the other stays the repository's). It is
for a worktree without the 3 GB input/ tree, or for tables kept
elsewhere; the run then tests this tree's code on those files.

conftest.py imports this module, sets the two variables for every
program a test starts, and offers the two roots as INPUT_ROOT and
START_ROOT; pytest's header prints them. Standard library only.
"""
import os
import pathlib

OPT_IN = "YREC_TEST_EXTERNAL_ROOTS"
VARIABLES = {"YREC_INPUT": "input", "YREC_START": "startmodels"}


def resolve(repo, environ=None):
    """{variable: path} for YREC_INPUT and YREC_START: the repository's
    directories, or the environment's values under the opt-in."""
    environ = os.environ if environ is None else environ
    external = environ.get(OPT_IN) == "1"
    roots = {}
    for variable, directory in VARIABLES.items():
        roots[variable] = pathlib.Path(repo) / directory
        if external and environ.get(variable):
            roots[variable] = pathlib.Path(environ[variable])
    return roots


def describe(repo, environ=None):
    """The lines a harness prints: the two roots, whose they are, and
    each value of the environment that was not used."""
    environ = os.environ if environ is None else environ
    roots = resolve(repo, environ)
    own = resolve(repo, {})
    lines = []
    for variable, what in (("YREC_INPUT", "tables"), ("YREC_START", "starting models")):
        if roots[variable] == own[variable]:
            whose = "this repository's"
        else:
            whose = f"NOT this repository's: {OPT_IN}=1 takes {variable} from the environment"
        lines.append(f"yrec {what}: {roots[variable]} ({whose})")
        named = environ.get(variable)
        if named and pathlib.Path(named) != roots[variable]:
            lines.append(f"  the environment's {variable}={named} is not used "
                         f"({OPT_IN}=1 would use it)")
    return lines


def apply(repo, env):
    """Set YREC_INPUT and YREC_START in env (a dict, or os.environ) to the
    roots resolve() gives for it; returns them."""
    roots = resolve(repo, env)
    for variable, path in roots.items():
        env[variable] = str(path)
    return roots


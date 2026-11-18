"""What the test harnesses share: where a run reads the tables and the
starting models, and a stand-in for a deck's directory.

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
START_ROOT; pytest's header prints them.

A test writes only under its temporary directory. A deck that names
files relative to its own directory (cm_monte_carlo: mc_input.dat, and
output/cm_monte_carlo.mc1 and .mc2, which a read creates) is read from
stand_in_directory()'s copy of that directory, and the decks that write
their final model under {YREC_START} are read with stand_in_tree()'s
copy of the starting models. Standard library only.
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


def stand_in_directory(directory, scratch):
    """Make scratch (a path that does not exist yet) a stand-in for a
    deck's directory: a link to every file beside the deck, and an empty
    output/. A program run there finds what the deck names relative to its
    directory, and what it creates under such names stays in scratch; a
    fresh checkout, which has no output/ directories (git ignores them),
    behaves as a used one. Returns scratch."""
    scratch = pathlib.Path(scratch)
    scratch.mkdir(parents=True)
    for entry in sorted(pathlib.Path(directory).iterdir()):
        if entry.is_file():
            (scratch / entry.name).symlink_to(entry)
    (scratch / "output").mkdir()
    return scratch


def stand_in_tree(directory, scratch):
    """Make scratch (a path that does not exist yet) a stand-in for a whole
    tree: its directories made anew, its files linked. A program given
    scratch for directory reads every file of the tree, and a file it
    creates there stays in scratch (a deck whose last_model_file is under
    {YREC_START}: the read of its controls creates that file). Returns
    scratch."""
    directory, scratch = pathlib.Path(directory), pathlib.Path(scratch)
    scratch.mkdir(parents=True)
    for entry in sorted(directory.rglob("*")):
        target = scratch / entry.relative_to(directory)
        if entry.is_dir():
            target.mkdir()
        else:
            target.symlink_to(entry)
    return scratch

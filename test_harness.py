"""The test harness tests its own tree (bug_sweep M374).

pytest (conftest.py), src/tools/devel_tests.py and src/tools/mutate.py
give every program they start this repository's input/ and startmodels/
as YREC_INPUT and YREC_START, whatever the shell exports for another
checkout; YREC_TEST_EXTERNAL_ROOTS=1 is the one opt-in that takes the
two from the environment (src/tools/harness.py). These tests hold
that rule: the function's cases, conftest.py as a fresh interpreter
imports it under a shell that names another tree, the two tools' report
of the roots, and the convention that keeps the rule true for every
root test (none names its own tree's input/ or startmodels/: they
import INPUT_ROOT and START_ROOT from conftest.py).

A test writes only under its temporary directory (M373): the two
stand-ins that test_controls.py and devel_tests.py read the decks
through are checked here too.

Needs no build and no table.
"""
import os
import pathlib
import re
import subprocess
import sys

REPO = pathlib.Path(__file__).parent
sys.path.insert(0, str(REPO / "src" / "tools"))
import harness  # noqa: E402

ELSEWHERE = {"YREC_INPUT": "/another/checkout/input", "YREC_START": "/another/checkout/startmodels"}
OWN = {"YREC_INPUT": REPO / "input", "YREC_START": REPO / "startmodels"}

# test_element_diffusion.py is not to be edited (docs/GUIDELINES.md); it sets
# YREC_INPUT to its own tree's input/ itself, so the opt-in does not reach it.
NOT_CONVERTED = {"test_element_diffusion.py"}


def test_the_roots_are_the_repositorys_whatever_the_environment_names():
    assert harness.resolve(REPO, {}) == OWN
    assert harness.resolve(REPO, ELSEWHERE) == OWN
    # the opt-in needs its exact value
    for value in ("", "0", "yes"):
        assert harness.resolve(REPO, dict(ELSEWHERE, YREC_TEST_EXTERNAL_ROOTS=value)) == OWN


def test_the_opt_in_takes_each_root_the_environment_names():
    opted = dict(ELSEWHERE, YREC_TEST_EXTERNAL_ROOTS="1")
    assert harness.resolve(REPO, opted) == {k: pathlib.Path(v) for k, v in ELSEWHERE.items()}
    # one that is unset or empty stays the repository's
    assert harness.resolve(REPO, {"YREC_TEST_EXTERNAL_ROOTS": "1", "YREC_INPUT": "/tables",
                                        "YREC_START": ""}) == \
        {"YREC_INPUT": pathlib.Path("/tables"), "YREC_START": OWN["YREC_START"]}
    env = dict(opted)
    assert harness.apply(REPO, env) == harness.resolve(REPO, opted)
    assert (env["YREC_INPUT"], env["YREC_START"]) == (ELSEWHERE["YREC_INPUT"], ELSEWHERE["YREC_START"])


def test_the_report_names_the_roots_and_what_was_not_used():
    plain = harness.describe(REPO, {})
    assert plain == [f"yrec tables: {OWN['YREC_INPUT']} (this repository's)",
                     f"yrec starting models: {OWN['YREC_START']} (this repository's)"]
    ignored = "\n".join(harness.describe(REPO, ELSEWHERE))
    for variable, value in ELSEWHERE.items():
        assert f"the environment's {variable}={value} is not used" in ignored
    assert str(OWN["YREC_INPUT"]) in ignored and "NOT this repository's" not in ignored
    opted = "\n".join(harness.describe(REPO, dict(ELSEWHERE, YREC_TEST_EXTERNAL_ROOTS="1")))
    assert f"yrec tables: {ELSEWHERE['YREC_INPUT']} (NOT this repository's" in opted
    assert "is not used" not in opted


def _conftest_roots(extra):
    """What conftest.py leaves in the environment of the programs a test
    starts, and offers as INPUT_ROOT and START_ROOT, in a fresh interpreter
    whose environment names another checkout."""
    env = {k: v for k, v in os.environ.items() if k != "YREC_TEST_EXTERNAL_ROOTS"}
    env.update(ELSEWHERE, **extra)
    code = ("import os, conftest\n"
            "print(os.environ['YREC_INPUT']); print(os.environ['YREC_START'])\n"
            "print(conftest.INPUT_ROOT); print(conftest.START_ROOT)\n")
    r = subprocess.run([sys.executable, "-c", code], cwd=REPO, env=env, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return r.stdout.split()


def test_conftest_sets_this_trees_roots_under_a_shell_that_names_another():
    own = [str(OWN["YREC_INPUT"]), str(OWN["YREC_START"])]
    assert _conftest_roots({}) == own + own
    elsewhere = [ELSEWHERE["YREC_INPUT"], ELSEWHERE["YREC_START"]]
    assert _conftest_roots({"YREC_TEST_EXTERNAL_ROOTS": "1"}) == elsewhere + elsewhere


def test_the_two_tools_report_this_trees_roots():
    """devel_tests.py and mutate.py print the roots before their cases (a
    selection that matches nothing: the report alone)."""
    env = {k: v for k, v in os.environ.items() if k != "YREC_TEST_EXTERNAL_ROOTS"}
    env.update(ELSEWHERE)
    for command in (["devel_tests.py", "--only", "^no such case$"],
                    ["mutate.py", "numerics", "--only", "^no such mutation$"]):
        r = subprocess.run([sys.executable, "-I", str(REPO / "src" / "tools" / command[0])] + command[1:],
                           cwd=REPO, env=env, capture_output=True, text=True)
        assert r.returncode == 0, (command, r.stdout, r.stderr)
        assert f"yrec tables: {OWN['YREC_INPUT']} (this repository's)" in r.stdout, (command, r.stdout)
        assert f"the environment's YREC_START={ELSEWHERE['YREC_START']} is not used" in r.stdout


def test_no_test_file_or_tool_names_its_own_trees_roots():
    """The rule holds for every root test only while each takes the roots
    from conftest.py, and for the two tools while they take them from
    harness.py: none builds <its tree>/input or /startmodels, and
    none keeps a root the environment set (setdefault)."""
    own_tree = re.compile(r"""REPO\s*/\s*["'](input|startmodels)["']|setdefault\(\s*["']YREC_(INPUT|START)""")
    files = [p for p in sorted(REPO.glob("test_*.py")) if p.name not in NOT_CONVERTED | {"test_harness.py"}]
    files += [REPO / "src" / "tools" / name for name in ("devel_tests.py", "mutate.py", "validate_model.py")]
    assert len(files) > 40
    offenders = [f"{p.relative_to(REPO)}:{n}: {line.strip()}"
                 for p in files for n, line in enumerate(p.read_text().splitlines(), 1) if own_tree.search(line)]
    assert not offenders, "\n".join(offenders)


def test_a_stand_in_directory_reads_as_the_deck_directory_and_keeps_what_is_created(tmp_path):
    """harness.stand_in_directory: the files beside a deck are there (as
    links), output/ is there and empty whether or not the deck's directory
    has one, and a file created in the stand-in is not in the tree."""
    deck_directory = REPO / "examples" / "run_config_matrix"
    stand_in = harness.stand_in_directory(deck_directory, tmp_path / "stand_in")
    assert (stand_in / "mc_input.dat").read_bytes() == (deck_directory / "mc_input.dat").read_bytes()
    assert (stand_in / "cm_monte_carlo.inlist").is_symlink()
    assert (stand_in / "output").is_dir() and not list((stand_in / "output").iterdir())
    (stand_in / "output" / "made_by_a_read").write_text("")
    (stand_in / "made_beside_the_deck").write_text("")
    assert not (deck_directory / "output" / "made_by_a_read").exists()
    assert not (deck_directory / "made_beside_the_deck").exists()


def test_a_stand_in_tree_reads_as_the_tree_and_keeps_what_is_created(tmp_path):
    """harness.stand_in_tree: every file of the tree at its place (as a
    link), every directory made anew, so that a file created anywhere in
    the stand-in is not in the tree."""
    tree = tmp_path / "tree"
    (tree / "dbl" / "deep").mkdir(parents=True)
    (tree / "top.first").write_text("top")
    (tree / "dbl" / "deep" / "model.first").write_text("model")
    stand_in = harness.stand_in_tree(tree, tmp_path / "stand_in")
    assert (stand_in / "top.first").read_text() == "top"
    assert (stand_in / "dbl" / "deep" / "model.first").read_text() == "model"
    assert (stand_in / "dbl" / "deep" / "model.first").is_symlink() and not (stand_in / "dbl").is_symlink()
    (stand_in / "dbl" / "deep" / "final.mod").write_text("")
    assert sorted(p.name for p in (tree / "dbl" / "deep").iterdir()) == ["model.first"]

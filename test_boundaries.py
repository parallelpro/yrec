"""Domain-boundary regression test (phase three, docs/ROADMAP.md stage 2).

Runs src/tools/check_boundaries.py, which enforces that every
cross-domain call in src/ goes through the target domain's public
surface (its <domain>_lib facade entries, or the documented de-facto
public routines of the facade-less domains). See that script's header
for the maintenance rules. Also runs the star-contract check and the
style ratchet (src/tools/check_style.py with src/tools/style_counts.tsv).
"""
import pathlib
import subprocess
import sys


def test_domain_boundaries():
    script = pathlib.Path(__file__).parent / "src" / "tools" / "check_boundaries.py"
    result = subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_star_contracts():
    """The `star reads:` / `star writes:` lists declared in procedure
    headers match the code (src/tools/check_star_contracts.py)."""
    script = pathlib.Path(__file__).parent / "src" / "tools" / "check_star_contracts.py"
    result = subprocess.run(
        [sys.executable, str(script)], capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_style_ratchet():
    """No style count of a directory of src/ rises above its allowed value
    in src/tools/style_counts.tsv (src/tools/check_style.py, after its
    scanner's self-test). The failure lists the new occurrences."""
    script = pathlib.Path(__file__).parent / "src" / "tools" / "check_style.py"
    for args in (["--self-test"], []):
        result = subprocess.run(
            [sys.executable, str(script)] + args, capture_output=True, text=True
        )
        assert result.returncode == 0, result.stdout + result.stderr


def test_star_tools_self_test():
    """The star tools' statement scan sees the forms a contract or a member
    table depends on: an assignment after `case (k);`, a labelled
    statement, an if-statement's condition, a ';' or `s%x` inside a string,
    and `real(dp) function` / `type(t) function` headers
    (check_star_contracts.py and gen_star_docs.py --self-test)."""
    tools = pathlib.Path(__file__).parent / "src" / "tools"
    for tool in ("check_star_contracts.py", "gen_star_docs.py"):
        result = subprocess.run(
            [sys.executable, str(tools / tool), "--self-test"], capture_output=True, text=True
        )
        assert result.returncode == 0, tool + "\n" + result.stdout + result.stderr

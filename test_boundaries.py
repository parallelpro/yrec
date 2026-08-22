"""Domain-boundary regression test (phase three, ROADMAP.md stage 2).

Runs src/tools/check_boundaries.py, which enforces that every
cross-domain call in src/ goes through the target domain's public
surface (its <domain>_lib facade entries, or the documented de-facto
public routines of the facade-less domains). See that script's header
for the maintenance rules.
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

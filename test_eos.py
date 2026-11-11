"""Standalone eos-domain regression test (phase three, docs/ROADMAP.md stage 2).

Runs src/test_eos (built by `make tests` in src/) with YREC_INPUT
pointing at the repository's input/ tables, and byte-compares its
stdout against the pinned src/eos/test/expected_test_eos.out.

Unlike the Stage-0 full-model regression, this exercises the eos
facade entries (eos_get, eos_get_gamma1) directly over a fixed grid,
without needing a stellar model to converge -- including branch
combinations no Stage-0 case reaches (e.g. eos_get_gamma1's Yale/SCV
branch). The MHD section self-reports SKIPPED unless YREC_MHD_TABLES
supplies tables (none ship with YREC); when tables are provided the
output will differ from the pinned baseline, which only covers the
shipped-data configuration.

The derivative check prints each analytic derivative and its finite
difference; their relative difference, whose leading digits are
round-off where the two agree, is not printed but checked (ok/FAIL)
where they must agree, and the binary exits nonzero when a check fails.

To regenerate the baseline after a deliberate eos change: rebuild,
run src/test_eos with YREC_INPUT set, review the diff like any other,
and commit the new expected file.
"""
import os
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_eos"
EXPECTED = REPO / "src" / "eos" / "test" / "expected_test_eos.out"


def test_eos_standalone(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_eos not built (run `make tests` in src/)")
    env = dict(os.environ)
    env["YREC_INPUT"] = str(REPO / "input")
    env.pop("YREC_MHD_TABLES", None)  # baseline covers shipped data only
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, env=env, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    compare_output(result.stdout, EXPECTED.read_text(), "test_eos")


SRC = REPO / "src"


def _fortran_value(text):
    """A Fortran literal default (.true., .false., 6.0d0, -1, ...) as Python."""
    text = text.strip().lower()
    if text in (".true.", ".false."):
        return text == ".true."
    return float(text.replace("d", "e"))


def test_eos_settings_defaults_are_the_controls():
    """eos_settings_t's declared defaults are the defaults of the controls
    the star layer copies into it (eos_configuration's header says so; M144).

    Each field the registry maps a control to (controls_registry.tsv's
    settings column, eos_settings%<field>) must default to that control's
    default in defaults/physics.defaults (generated from the controls'
    declarations; `make check` keeps it current), and opal_vintage, which
    derive_controls takes from the three OPAL vintage switches, must
    default to opal_vintage_none while all three default to .false. (its
    own registry row is a derived control, with no default of its own).
    The Fortran test's "the record's defaults are the controls'" checks the
    record against literals; this one checks it against the controls.
    """
    config = (SRC / "eos" / "eos_configuration.f90").read_text()
    block = config[config.index("type, public :: eos_settings_t"):
                   config.index("end type eos_settings_t")]
    declared = {}
    for line in block.splitlines():
        code = line.split("!")[0]
        if "::" not in code:
            continue
        for item in code.split("::", 1)[1].split(","):
            if "=" in item:
                name, value = item.split("=", 1)
                declared[name.strip()] = value.strip()
    defaults = {}
    for line in (SRC / "defaults" / "physics.defaults").read_text().splitlines():
        code = line.split("!")[0]
        if "=" in code:
            name, value = code.split("=", 1)
            defaults[name.strip()] = value.strip()
    mapped = {}
    rows = (SRC / "defaults" / "controls_registry.tsv").read_text().splitlines()
    header = rows[0].split("\t")
    for row in rows[1:]:
        cells = dict(zip(header, row.split("\t")))
        for target in cells.get("settings", "").split(","):
            if target.startswith("eos_settings%"):
                mapped[target.split("%", 1)[1]] = cells["proposed_name"]
    assert mapped, "no control maps to an eos_settings_t field"
    vintage_control = mapped.pop("opal_vintage", None)
    for field, control in sorted(mapped.items()):
        assert field in declared, f"eos_settings_t%{field} has no declared default"
        assert _fortran_value(declared[field]) == _fortran_value(defaults[control]), (
            f"eos_settings_t%{field} defaults to {declared[field]}, "
            f"its control {control} to {defaults[control]}")
    assert declared["opal_vintage"] == "opal_vintage_none"
    assert vintage_control is None or vintage_control not in defaults, vintage_control
    for switch in ("use_opal95_eos", "use_opal2001_eos", "use_opal2006_eos"):
        assert _fortran_value(defaults[switch]) is False, switch
    unmapped = set(declared) - set(mapped) - {"opal_vintage"}
    assert not unmapped, f"eos_settings_t fields no control sets: {sorted(unmapped)}"

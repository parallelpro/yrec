"""io's own procedures against their contracts.

Runs src/test_io (built by `make tests` in src/) in a scratch directory,
where it writes every file it reads (no deck, no tables), and compares
its stdout with the pinned src/io/test/expected_test_io.out: the list
of checks, each `ok`. The binary exits nonzero when a check fails. It
checks the path placeholders' expansion (expand_value, with the
environment set through the C library), directory_of,
pulse_format_code, the columns-file selection (parse_columns) and its
errors, nml_file_has_group, a bit-exact .mod write and read with its
replace/append modes, an unknown global and the reader's errors, the
model files' format tag, synthetic YREC7 and MODEL2 starting models
(every value placed where the layout says, and the files cut or
damaged), the history and profile column tables (unique names,
valid kinds, the default selections), the writes' failures (a failed
write, a directory that cannot be made, a full disk through a file-size
limit, and two runs that write one path, which is not a full disk) and
FGONG's globals. Each check was confirmed to fail
under a mutation of the code it covers (see the header of
src/io/test/test_io.f90).
"""
import pathlib
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_io"
EXPECTED = REPO / "src" / "io" / "test" / "expected_test_io.out"


def test_io(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_io not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "FAIL" not in result.stdout, result.stdout
    compare_output(result.stdout, EXPECTED.read_text(), "test_io")

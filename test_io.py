"""io's own procedures against their contracts.

Runs src/test_io (built by `make tests` in src/) in a scratch directory,
where it writes every file it reads (no deck, no tables), and compares
its stdout with the pinned src/io/test/expected_test_io.out: the list
of checks, each `ok`. The binary exits nonzero when a check fails. It
checks the path placeholders' expansion (expand_value, with the
environment set through the C library), directory_of, the run log's
path (run_log_path: only a .short ending becomes .log), pulse_format_code, the columns-file selection (parse_columns) and its
errors, nml_file_has_group, a bit-exact .mod write and read with its
replace/append modes, an unknown global and the reader's errors, the
model files' format tag, synthetic YREC7 and MODEL2 starting models
(every value placed where the layout says, and the files cut or
damaged), the history and profile column tables (unique names,
valid kinds, the default selections), the writes' failures (a failed
write, a directory that cannot be made, a full disk through a file-size
limit, and two runs that write one path, which is not a full disk; a
file written in one go renamed to its path) and FGONG's globals. Each
check was confirmed to fail under a mutation of the code it covers (see
the header of src/io/test/test_io.f90).

test_two_runs_copy_their_decks_to_one_path runs two copies of the binary
at once, each copying its own deck to one inlist_used many times over,
as the decks of one directory all write output/inlist_used when they
start (docs/bug_sweep.md M177). test_two_runs_write_one_history_path
runs two at once on one history path, one writing many rows into one
open while the other opens the path anew many times, as the decks of one
directory all write output/history.data (M174).
"""
import pathlib
import subprocess
import time

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


# each run's copies of its deck (two runs of this many copies each met the
# old code's race several times per run: a false "is the disk full?" or a
# partial file at the path)
COPIES = 2000


def test_two_runs_copy_their_decks_to_one_path(tmp_path):
    """Two runs that copy their decks to one output/inlist_used at once
    (M177) both succeed, the path holds one run's whole copy whenever it is
    read while they run and when they end, and no file of either run's own
    is left beside it. Each run writes its copy under its own name and
    renames it over the path (output_files' close_replacement_file)."""
    if not BINARY.exists():
        pytest.skip("src/test_io not built (run `make tests` in src/)")
    decks = []
    for k, n_lines in ((1, 120), (2, 80)):
        deck = tmp_path / f"deck{k}.nml"
        deck.write_text("".join(f" line_{i} = {k}  ! deck {k}, a line of a deck of some kilobytes\n"
                                for i in range(n_lines)))
        decks.append(deck)
    whole = {f"! ==== copied verbatim from: {d.name}\n" + d.read_text() for d in decks}
    out = tmp_path / "output"
    out.mkdir()
    dest = out / "inlist_used"
    runs = [subprocess.Popen([str(BINARY), "copy", d.name, "output/inlist_used", str(COPIES)],
                             cwd=tmp_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             text=True) for d in decks]
    partial = 0
    deadline = time.monotonic() + 300
    while any(r.poll() is None for r in runs) and time.monotonic() < deadline:
        try:
            text = dest.read_text()
        except FileNotFoundError:
            continue
        if text not in whole:
            partial += 1
    results = [r.communicate(timeout=60) for r in runs]
    for k, (r, (stdout, stderr)) in enumerate(zip(runs, results), 1):
        assert r.returncode == 0, f"run {k}: {stdout}{stderr}"
    assert partial == 0, f"the path held a partial copy {partial} times while the runs wrote it"
    assert dest.read_text() in whole, "the copy left at the path is not one run's whole copy"
    assert sorted(p.name for p in out.iterdir()) == ["inlist_used"], sorted(p.name for p in out.iterdir())


# one run's rows into one open of the path, and the other run's opens of it
# (the code before the fix failed the first run with "is the disk full?" in
# five runs of five, within half a second)
ROWS, OPENS = 20000, 4000


def test_two_runs_write_one_history_path(tmp_path):
    """Two runs that keep one path open for writing at once, as the decks of
    one directory all write output/history.data (M174): one writes many rows
    into one open, each flushed and checked on disk as a history row is,
    while the other opens the path anew many times. Both succeed: each open
    makes a new file (output_files' open_output_file removes the path
    first), so the first run's checks see its own rows and never a file cut
    under it; and the path holds one open's whole file, nothing beside it."""
    if not BINARY.exists():
        pytest.skip("src/test_io not built (run `make tests` in src/)")
    out = tmp_path / "output"
    out.mkdir()
    runs = [subprocess.Popen([str(BINARY), "history", "output/history.data", times, rows],
                             cwd=tmp_path, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            for times, rows in (("1", str(ROWS)), (str(OPENS), "3"))]
    results = [r.communicate(timeout=300) for r in runs]
    for k, (r, (stdout, stderr)) in enumerate(zip(runs, results), 1):
        assert r.returncode == 0, f"run {k}: {stdout}{stderr}"
    lines = (out / "history.data").read_text().splitlines()
    opens = {l.split()[1] for l in lines}
    assert len(opens) == 1 and [l.split()[3] for l in lines] == [str(k) for k in range(1, len(lines) + 1)] \
        and len(lines) in (3, ROWS), lines[:5]
    assert sorted(p.name for p in out.iterdir()) == ["history.data"], sorted(p.name for p in out.iterdir())

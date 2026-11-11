"""Writing generated files, or checking them.

A generator calls emit(path, text) for each output (or writes through
`with Out(path) as f:`) and finish(name) at the end. Run plainly, it
writes the files whose text changes and leaves the others alone (their
times too, so that make recompiles only what uses a changed file: make
runs the controls generators, src/Makefile); run with --check, it writes
nothing and exits 1 if any committed output differs from what it would
write (`make check` in src/ and test_generated_files.py run every generator
this way).
"""
import io
import sys

CHECK = "--check" in sys.argv
_stale = []


def emit(path, text):
    current = path.is_file() and path.read_text() == text
    if CHECK:
        if not current:
            _stale.append(path)
    elif not current:
        path.write_text(text)


class Out(io.StringIO):
    """a text buffer that emit()s to path when its with-block ends"""

    def __init__(self, path):
        super().__init__()
        self.path = path

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        if exc[0] is None:
            emit(self.path, self.getvalue())
        return False


def finish(generator):
    if not CHECK:
        return
    if _stale:
        for p in _stale:
            print(f"{p.name} is stale: rerun tools/{generator}", file=sys.stderr)
        sys.exit(1)
    print(f"{generator}: outputs up to date")
    sys.exit(0)

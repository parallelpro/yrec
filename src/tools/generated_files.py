"""Writing generated files, or checking them.

A generator calls emit(path, text) for each output (or writes through
`with Out(path) as f:`) and finish(name) at the end. Run plainly, it
writes the files; run with --check, it writes nothing and exits 1 if any
committed output differs from what it would write (test_generated_files.py
runs every generator this way).
"""
import io
import sys

CHECK = "--check" in sys.argv
_stale = []


def emit(path, text):
    if CHECK:
        if not path.is_file() or path.read_text() != text:
            _stale.append(path)
    else:
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

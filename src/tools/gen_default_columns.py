#!/usr/bin/env python3
"""Generate src/io/default_columns.inc from defaults/*_columns.list, and
check the .list files against the writers' column tables.

defaults/history_columns.list and defaults/profile_columns.list are
the authority on the default output column selection: a blank
history_columns_file / profile_columns_file control selects exactly
the uncommented names in these files, in file order. Opt-in columns
ship commented out ('!' immediately followed by the name); remove
the '!' and rebuild (the Makefile reruns this generator when a .list
changes) to add one to the default output.

The selection is compiled in rather than read at runtime so a run
never depends on locating the source tree from its working
directory (the {YREC_INPUT}-style relative fallback has a history of
breaking) and a deployed bare binary keeps working.

The column tables are Fortran: history_columns in io/history_output.f90
(named indices ih_<name>) and profile_columns in
core/diagnostics/stitched_model.f90 (ip_<name>), one row per column,
`cols(ih_<name>) = output_column('<name>', <kind>)`. --check also
checks that each table names every column of its .list and no other,
that each row sits at the index named after its column, and that the
indices run 1..n. test_mesa_output.py::test_default_columns_lists_in_sync
uses column_tables() for the same checks and runs --check.

Usage:
  python3 tools/gen_default_columns.py           # (re)write the .inc
  python3 tools/gen_default_columns.py --check   # exit 1 if stale or out of sync (make check, CI)
"""
import pathlib
import re
import sys

SRC = pathlib.Path(__file__).resolve().parent.parent
OUT = SRC / "io" / "default_columns.inc"
WIDTH = 24  # output_columns' column_name_len

# (tag, .list file, Fortran file, table function, index prefix)
TABLES = [
    ("hist", "history_columns.list", "io/history_output.f90",
     "history_columns", "ih_"),
    ("prof", "profile_columns.list", "core/diagnostics/stitched_model.f90",
     "profile_columns", "ip_"),
]
ROW_RE = re.compile(
    r"^\s*cols\((\w+)\)\s*=\s*output_column\('([^']+)',\s*(\w+)\)", re.M)


def names_from(fname):
    """The default selection: the uncommented names, in file order."""
    out = []
    for line in (SRC / "defaults" / fname).read_text().splitlines():
        line = line.split("!")[0].strip()
        if line:
            if len(line) > WIDTH:
                sys.exit(f"{fname}: '{line}' longer than {WIDTH} chars")
            out.append(line)
    return out


def list_entries(fname):
    """Every column a .list documents: its names, and '!name' entries
    (opt-in columns); '! prose' comments are skipped."""
    out = []
    for line in (SRC / "defaults" / fname).read_text().splitlines():
        m = re.match(r"^(!?)([A-Za-z]\w*)", line)
        if m:
            out.append(m.group(2))
    return out


def column_table(path, function):
    """The rows of a Fortran column table: (index name, column name, kind)."""
    text = (SRC / path).read_text()
    m = re.search(r"function " + function + r"\(\)(.*?)end function " + function,
                  text, re.S)
    if not m:
        sys.exit(f"{path}: no function {function}")
    return ROW_RE.findall(m.group(1))


def index_values(path, prefix):
    """The named indices <prefix><name> = <n> declared in a Fortran file."""
    text = (SRC / path).read_text()
    return {k.lower(): int(v) for k, v in
            re.findall(r"\b(" + prefix + r"\w+)\s*=\s*(\d+)\b", text)}


def column_tables():
    """{tag: (rows, indices, list entries, list file)} for both outputs."""
    return {tag: (column_table(path, func), index_values(path, prefix),
                  list_entries(lst), lst, prefix)
            for tag, lst, path, func, prefix in TABLES}


def sync_errors():
    """What is out of sync between the tables, their indices and the .lists."""
    errors = []
    for tag, (rows, idx, entries, lst, prefix) in column_tables().items():
        names = [n for _, n, _ in rows]
        if sorted(names) != sorted(entries):
            errors.append(f"{lst} vs the {tag} table: only in the list "
                          f"{sorted(set(entries) - set(names))}, only in the "
                          f"table {sorted(set(names) - set(entries))}, "
                          f"duplicated {sorted({n for n in entries if entries.count(n) > 1})}")
        for index, name, _ in rows:
            if index.lower() != (prefix + name).lower():
                errors.append(f"{tag} row {name} sits at index {index}")
        values = sorted(idx.get((prefix + n).lower(), -1) for n in names)
        if values != list(range(1, len(names) + 1)):
            errors.append(f"{tag} indices are not 1..{len(names)}: {values}")
    return errors


def emit(tag, names):
    lines = [f"      integer, parameter :: n_{tag}_default = {len(names)}",
             f"      character(len=column_name_len), parameter :: "
             f"{tag}_default_names(n_{tag}_default) = [ &"]
    for i, n in enumerate(names):
        sep = ", &" if i < len(names) - 1 else " ]"
        lines.append(f"           '{n:<{WIDTH}}'{sep}")
    return lines


def main():
    lines = [
        "! GENERATED by tools/gen_default_columns.py from",
        "! defaults/history_columns.list and",
        "! defaults/profile_columns.list -- do not edit; edit the .list",
        "! files and rebuild (make reruns the generator). These are the",
        "! default output column selections used when",
        "! history_columns_file / profile_columns_file are blank.",
    ]
    hist = names_from("history_columns.list")
    prof = names_from("profile_columns.list")
    lines += emit("hist", hist)
    lines += emit("prof", prof)
    content = "\n".join(lines) + "\n"

    if "--check" in sys.argv:
        errors = sync_errors()
        for e in errors:
            print(e, file=sys.stderr)
        if not OUT.is_file() or OUT.read_text() != content:
            print("default_columns.inc is stale: rerun "
                  "tools/gen_default_columns.py", file=sys.stderr)
            sys.exit(1)
        if errors:
            sys.exit(1)
        print("default_columns.inc up to date; column tables in sync")
        return
    OUT.write_text(content)
    print(f"{len(hist)} history + {len(prof)} profile default columns "
          "-> io/default_columns.inc")


if __name__ == "__main__":
    main()

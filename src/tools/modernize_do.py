#!/usr/bin/env python3
"""Wave 1 of the goto/save eradication (2026): numbered DO -> block DO.

Transform, uniformly semantics-preserving:
  - `do <label> [,] var = ...`  ->  `do var = ...`
  - the labeled terminator LINE IS KEPT AS-IS (so any `goto <label>`
    from inside the loop still lands on it and then falls into the
    new `end do` -- exactly F77's jump-to-terminator = cycle
    semantics), and one `end do` is appended per DO that shared the
    terminator (F77 shared-terminator nesting).
  - a terminator that is a real statement (F77 allowed `100 x=y` as a
    loop end) is likewise kept, with `end do` appended after it.

Usage: python3 tools/modernize_do.py <dir-or-file> [...]
"""
import pathlib
import re
import sys


DO_RE = re.compile(r"^(\s*)do\s+(\d+)\s*,?\s*(\S.*)$", re.I)
LABEL_RE = re.compile(r"^(\s*)(\d+)\s+(.*)$")


def code(line):
    return line.split("!")[0]


def convert(path):
    lines = path.read_text().splitlines(keepends=True)
    out = []
    stack = []  # open numbered-do labels, innermost last
    changed = 0
    for line in lines:
        c = code(line)
        m = DO_RE.match(c)
        if m and not line.lstrip().startswith("!"):
            indent, label, rest = m.groups()
            # keep any trailing comment
            comment = line[len(c):] if len(line) > len(c) else "\n"
            out.append(f"{indent}do {rest.strip()}{comment.rstrip()}\n"
                       if comment.strip() else f"{indent}do {rest.strip()}\n")
            stack.append((label, indent))
            changed += 1
            continue
        lm = LABEL_RE.match(c)
        if lm and stack and lm.group(2) == stack[-1][0]:
            label = lm.group(2)
            out.append(line)
            while stack and stack[-1][0] == label:
                _, indent = stack.pop()
                out.append(f"{indent}end do\n")
            continue
        out.append(line)
    if stack:
        raise RuntimeError(f"{path}: unterminated numbered do {stack}")
    if changed:
        path.write_text("".join(out))
    return changed


def main():
    total = 0
    for arg in sys.argv[1:]:
        p = pathlib.Path(arg)
        files = [p] if p.is_file() else sorted(p.rglob("*.f90"))
        for f in files:
            if "test" in f.parts:
                continue
            n = convert(f)
            if n:
                print(f"{f}: {n} numbered do -> block do")
                total += n
    print(f"TOTAL converted: {total}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""goto campaign pass 1: loop-control gotos -> cycle / exit.

Transforms, with full do-nesting verification (innermost loop only,
since unlabeled cycle/exit bind to the innermost construct):
  - `go to L` targeting the innermost open loop's kept labeled
    terminator (`L continue` + `end do`)          -> cycle
  - `go to L` targeting the first statement after the innermost open
    loop's `end do`                               -> exit
Only the goto text changes; labels stay (cleaned up separately once
unreferenced).
"""
import pathlib
import re
import sys


def code(l):
    return l.split("!")[0]


GOTO_RE = re.compile(r"\bgo\s*to\s+(\d+)\b", re.I)
DO_RE = re.compile(r"^\s*(\d+\s+)?do\b(?!\s*\w*\s*=\s*\()", re.I)
ENDDO_RE = re.compile(r"^\s*(\d+\s+)?end\s*do\b", re.I)
LABEL_RE = re.compile(r"^\s*(\d+)\s")


def transform(path):
    lines = path.read_text().splitlines(keepends=True)
    n = len(lines)
    labels = {}
    for i, l in enumerate(lines):
        m = LABEL_RE.match(code(l))
        if m:
            labels[m.group(1)] = i

    # build do-nesting: for each line, the stack of (do_line, end_line)
    stack, pairs, open_at = [], {}, [[] for _ in range(n)]
    for i, l in enumerate(lines):
        c = code(l)
        if DO_RE.match(c) and not re.match(r"^\s*(\d+\s+)?do\s*while", c, re.I) \
           and not re.match(r"^\s*(\d+\s+)?do\s*$", c) is None or True:
            pass
    # simpler second pass: treat any `do` (incl. do while, bare do) as opener
    stack = []
    enddo_of = {}
    for i, l in enumerate(lines):
        c = code(l)
        if re.match(r"^\s*(\d+\s+)?do\b", c, re.I) and "end do" not in c.lower():
            stack.append(i)
        elif ENDDO_RE.match(c):
            if stack:
                start = stack.pop()
                enddo_of[start] = i
        open_at[i] = list(stack)
    if stack:
        raise RuntimeError(f"{path}: unbalanced do nesting")

    changed = 0
    out = list(lines)
    for i, l in enumerate(lines):
        c = code(l)
        m = GOTO_RE.search(c)
        if not m:
            continue
        lab = m.group(1)
        tgt = labels.get(lab)
        if tgt is None or not open_at[i]:
            continue
        inner_do = open_at[i][-1]
        inner_end = enddo_of.get(inner_do)
        if inner_end is None:
            continue
        # cycle: target is the labeled terminator line immediately
        # before this loop's `end do` (wave-1 kept terminators), or any
        # labeled `continue` directly inside this loop just before end
        tgt_stmt = re.sub(r"^\s*\d+\s*", "", code(lines[tgt])).strip().lower()
        if tgt == inner_end - 1 and tgt_stmt == "continue" and tgt > i:
            out[i] = GOTO_RE.sub("cycle", l, count=1)
            changed += 1
            continue
        # also: target IS the labeled end do of the innermost loop
        if tgt == inner_end and tgt > i:
            out[i] = GOTO_RE.sub("cycle", l, count=1)
            changed += 1
            continue
        # exit: target is the first statement after the innermost end do
        j = inner_end + 1
        while j < n and (not code(lines[j]).strip()):
            j += 1
        if tgt == j and tgt > i:
            out[i] = GOTO_RE.sub("exit", l, count=1)
            changed += 1
            continue
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
            k = transform(f)
            if k:
                print(f"{f}: {k}")
                total += k
    print(f"TOTAL transformed: {total}")


if __name__ == "__main__":
    main()

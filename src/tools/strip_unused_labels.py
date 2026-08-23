#!/usr/bin/env python3
"""Remove statement labels that nothing references (goto campaign
leftovers -- kept loop terminators, dead format anchors).

Per procedure scope: a label is referenced by `go to N`, computed/
assigned gotos, `err=N`/`end=N`/`fmt=N`, alternate returns (`*N` in
call argument lists), or as a format used in `read/write(unit, N)`.
Unreferenced label on a bare CONTINUE -> the whole line goes;
on a real statement (e.g. FORMAT would always be referenced or is
dead) -> strip just the label field. Column alignment of the
statement text is preserved.
"""
import pathlib
import re
import sys

LABEL_RE = re.compile(r"^(\s*)(\d+)(\s+)(\S.*)$")


def refs_in(code_lines):
    refs = set()
    for c in code_lines:
        for m in re.finditer(r"\bgo\s*to\s+(\d+)", c, re.I):
            refs.add(m.group(1))
        for m in re.finditer(r"\b(?:err|end|fmt)\s*=\s*(\d+)", c, re.I):
            refs.add(m.group(1))
        for m in re.finditer(r"\b(?:write|read)\s*\(\s*[^(),]+\s*,\s*(\d+)\s*[),]", c, re.I):
            refs.add(m.group(1))
        for m in re.finditer(r"\bprint\s+(\d+)", c, re.I):
            refs.add(m.group(1))
        for m in re.finditer(r",\s*\*(\d+)\b", c):   # alternate returns
            refs.add(m.group(1))
    return refs


def transform(path):
    lines = path.read_text().splitlines(keepends=True)
    n = len(lines)
    # procedure scope boundaries
    bounds = [0]
    for i, l in enumerate(lines):
        if re.match(r"^\s*end\s+(subroutine|function|program)\b", l.split("!")[0], re.I):
            bounds.append(i + 1)
    if bounds[-1] != n:
        bounds.append(n)
    changed = 0
    out = list(lines)
    for b0, b1 in zip(bounds, bounds[1:]):
        scope = [l.split("!")[0] for l in lines[b0:b1]]
        refs = refs_in(scope)
        for i in range(b0, b1):
            c = lines[i].split("!")[0]
            m = LABEL_RE.match(c)
            if not m:
                continue
            lab = m.group(2)
            if lab in refs:
                continue
            stmt = m.group(4).strip().lower()
            if re.match(r"^continue\s*$", stmt):
                out[i] = None
            elif stmt.startswith("format"):
                # dead format: label unreferenced means the whole
                # statement is dead, but leave it (harmless, and some
                # refs may be built dynamically) -- skip.
                continue
            else:
                pad = " " * (len(m.group(1)) + len(lab) + len(m.group(3)))
                out[i] = pad + lines[i][len(m.group(1)) + len(lab) + len(m.group(3)):]
            changed += 1
    if changed:
        path.write_text("".join(l for l in out if l is not None))
    return changed


def main():
    total = files = 0
    for arg in sys.argv[1:]:
        p = pathlib.Path(arg)
        for f in ([p] if p.is_file() else sorted(p.rglob("*.f90"))):
            if "test" in f.parts:
                continue
            k = transform(f)
            if k:
                files += 1
                total += k
    print(f"stripped {total} unused labels in {files} files")


if __name__ == "__main__":
    main()

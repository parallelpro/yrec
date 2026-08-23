#!/usr/bin/env python3
"""goto campaign pass 3: the F77 if/else idiom -> real if/else.

    if (cond) go to L1        if (.not. cond) then
    <A>                          <A>
    go to L2          ->      else
 L1 <first stmt of B>            <B...>
    <B...>                    end if
 L2 <after>                 <after>

Requires: L1 and L2 each single-reference, both spans block-balanced,
no external jumps into either span. Conservative: skip on any doubt.
"""
import pathlib
import re
import sys


def code(l):
    return l.split("!")[0]


def span_ok(lines, a, b, label_refs):
    depth = 0
    for j in range(a, b):
        cj = code(lines[j]).strip().lower()
        if not cj:
            continue
        cjs = re.sub(r"^\d+\s*", "", cj)
        if re.match(r"do\b", cjs) and not cjs.startswith("double"):
            depth += 1
        elif cjs.endswith("then") and cjs.startswith("if"):
            depth += 1
        elif re.match(r"end\s*(do|if)\b", cjs):
            depth -= 1
        elif re.match(r"(else\b|case\b|contains\b|end\s+(subroutine|function)|entry\b)", cjs):
            if depth == 0:
                return False
        elif cjs.startswith("if") and cjs.endswith("&"):
            return False
        if depth < 0:
            return False
    if depth != 0:
        return False
    for j in range(a, b):
        mm = re.match(r"^\s*(\d+)\s", code(lines[j]))
        if mm:
            for r in label_refs.get(mm.group(1), []):
                if r < a or r >= b:
                    return False
    return True


def transform(path):
    lines = path.read_text().splitlines(keepends=True)
    labels, label_refs = {}, {}
    for i, l in enumerate(lines):
        c = code(l)
        m = re.match(r"^\s*(\d+)\s", c)
        if m:
            labels[m.group(1)] = i
        for g in re.finditer(r"\bgo\s*to\s+(\d+)", c, re.I):
            label_refs.setdefault(g.group(1), []).append(i)
        for g in re.finditer(r"\b(?:err|end)\s*=\s*(\d+)", c, re.I):
            label_refs.setdefault(g.group(1), []).append(i)
    edits = []
    for i, l in enumerate(lines):
        c = code(l)
        m = re.match(r"^(\s*)if\s*(\(.*\))\s*go\s*to\s+(\d+)\s*$", c, re.I)
        if not m:
            continue
        indent, cond, lab1 = m.groups()
        if cond.count("(") != cond.count(")"):
            continue
        t1 = labels.get(lab1)
        if t1 is None or t1 <= i or len(label_refs.get(lab1, [])) != 1:
            continue
        # the line before L1 must be an unconditional goto L2
        j = t1 - 1
        while j > i and not code(lines[j]).strip():
            j -= 1
        m2 = re.match(r"^\s*go\s*to\s+(\d+)\s*$", code(lines[j]), re.I)
        if not m2:
            continue
        lab2 = m2.group(1)
        t2 = labels.get(lab2)
        if t2 is None or t2 <= t1:
            continue  # L2 may be multi-ref: other jumps to it exit the new block, which is legal
        if not span_ok(lines, i + 1, j, label_refs):
            continue
        if not span_ok(lines, t1, t2, label_refs):
            continue
        edits.append((i, indent, cond, j, t1, t2))
    # non-crossing check + apply by t2 descending
    spans = [(e[0], e[5]) for e in edits]
    def crosses(a, b):
        return (a[0] < b[0] < a[1] < b[1]) or (b[0] < a[0] < b[1] < a[1])
    keep = [e for e in edits
            if not any(crosses((e[0], e[5]), s) for s in spans if s != (e[0], e[5]))]
    changed = 0
    for i, indent, cond, jgoto, t1, t2 in sorted(keep, key=lambda e: -e[5]):
        lines.insert(t2, f"{indent}end if\n")
        # strip the label from L1's line (it becomes plain first stmt of else)
        lm = re.match(r"^(\s*)(\d+)\s+(.*)$", lines[t1], flags=re.S)
        body = lm.group(3) if lm else lines[t1]
        lines[t1] = f"{indent}   {body}"
        lines[jgoto] = f"{indent}else\n"
        lines[i] = f"{indent}if (.not. {cond}) then\n"
        changed += 1
    if changed:
        path.write_text("".join(lines))
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
    print(f"TOTAL if/else: {total}")


if __name__ == "__main__":
    main()

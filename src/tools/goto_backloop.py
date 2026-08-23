#!/usr/bin/env python3
"""goto campaign pass 4: single-reference backward conditional gotos
-> do / exit loops.

 <L> <first stmt>              do
     <body>            ->         <first stmt>
     if (cond) go to L             <body>
                                   if (.not. (cond)) exit
                                end do
Requires: L single-reference, backward jump, span block-balanced, no
external jumps into the span. Conservative: skip on any doubt.
"""
import pathlib
import re
import sys


def code(l):
    return l.split("!")[0]


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
        indent, cond, lab = m.groups()
        if cond.count("(") != cond.count(")"):
            continue
        tgt = labels.get(lab)
        if tgt is None or tgt >= i or len(label_refs.get(lab, [])) != 1:
            continue
        # span balance from label line to goto
        depth = 0
        ok = True
        for j in range(tgt, i):
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
                    ok = False
                    break
            elif cjs.startswith("if") and cjs.endswith("&"):
                ok = False
                break
            if depth < 0:
                ok = False
                break
        if not ok or depth != 0:
            continue
        bad = False
        for j in range(tgt + 1, i):
            mm = re.match(r"^\s*(\d+)\s", code(lines[j]))
            if mm:
                for r in label_refs.get(mm.group(1), []):
                    if r < tgt or r > i:
                        bad = True
        if bad:
            continue
        edits.append((tgt, i, indent, cond))
    spans = [(e[0], e[1]) for e in edits]
    def crosses(a, b):
        return (a[0] < b[0] < a[1] < b[1]) or (b[0] < a[0] < b[1] < a[1])
    keep = [e for e in edits
            if not any(crosses((e[0], e[1]), s) for s in spans if s != (e[0], e[1]))]
    changed = 0
    for tgt, i, indent, cond in sorted(keep, key=lambda e: -e[1]):
        lines[i] = (f"{indent}if (.not. {cond}) exit\n{indent}end do\n")
        lm = re.match(r"^(\s*)(\d+)\s+(.*)$", lines[tgt], flags=re.S)
        body = lm.group(3) if lm else lines[tgt].lstrip()
        lines[tgt] = f"{indent}do\n{indent}   {body}"
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
    print(f"TOTAL back-loops: {total}")


if __name__ == "__main__":
    main()

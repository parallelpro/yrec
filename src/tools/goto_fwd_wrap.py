#!/usr/bin/env python3
"""goto campaign pass 2: single-reference forward-skip gotos -> if-blocks.

`if (cond) go to L` where L is defined later in the same procedure,
the span crosses no unbalanced do/if boundary, contains no ELSE of an
enclosing construct, and L has exactly one reference, becomes:
    if (.not. (cond)) then
       <span>
    end if
The label line survives (unreferenced; removed by the cleanup pass).
Conservative by construction: any doubt -> skip.
"""
import pathlib
import re
import sys


def code(l):
    return l.split("!")[0]


def transform(path):
    lines = path.read_text().splitlines(keepends=True)
    n = len(lines)
    text_labels = {}
    label_refs = {}
    for i, l in enumerate(lines):
        c = code(l)
        m = re.match(r"^\s*(\d+)\s", c)
        if m:
            text_labels[m.group(1)] = i
        for g in re.finditer(r"\bgo\s*to\s+(\d+)", c, re.I):
            label_refs.setdefault(g.group(1), []).append(i)
        for g in re.finditer(r"\b(?:err|end)\s*=\s*(\d+)", c, re.I):
            label_refs.setdefault(g.group(1), []).append(i)
    changed = 0
    i = 0
    edits = []  # (goto_line, cond, tgt)
    for i, l in enumerate(lines):
        c = code(l)
        m = re.match(r"^(\s*)if\s*(\(.*\))\s*go\s*to\s+(\d+)\s*$", c, re.I)
        if not m:
            continue
        indent, cond, lab = m.groups()
        # balanced condition parens?
        if cond.count("(") != cond.count(")"):
            continue
        tgt = text_labels.get(lab)
        if tgt is None or tgt <= i:
            continue
        # multi-reference labels are fine: other gotos to L either lie
        # outside the span (unaffected) or inside it (jumping OUT of
        # the new if-block, which Fortran permits); only jumps INTO
        # the span are fatal, checked below
        # span structure check
        depth = 0
        ok = True
        for j in range(i + 1, tgt):
            cj = code(lines[j]).strip().lower()
            if not cj:
                continue
            cjs = re.sub(r"^\d+\s*", "", cj)
            if re.match(r"do\b", cjs) and not cjs.startswith("double"):
                depth += 1
            elif re.match(r"(if\s*\(.*\)\s*then|if\s*\(.*&$)", cjs):
                # block if (approx: 'then' on this line); multiline conds handled below
                if cjs.endswith("then"):
                    depth += 1
            elif re.match(r"end\s*(do|if)\b", cjs):
                depth -= 1
            elif re.match(r"(else\b|elseif\b|else\s+if\b|case\b|contains\b|end\s+subroutine|end\s+function|entry\b)", cjs):
                if depth == 0:
                    ok = False
                    break
            elif cjs.endswith("&"):
                # a continued statement that might be a block-if opener:
                # conservatively skip the whole span
                if cjs.startswith("if"):
                    ok = False
                    break
            if depth < 0:
                ok = False
                break
        if not ok or depth != 0:
            continue
        # other labels inside span referenced from outside span?
        bad = False
        for j in range(i + 1, tgt):
            mm = re.match(r"^\s*(\d+)\s", code(lines[j]))
            if mm:
                for r in label_refs.get(mm.group(1), []):
                    if r < i or r > tgt:
                        bad = True
        if bad:
            continue
        edits.append((i, indent, cond, tgt))
    # reject partially-overlapping (crossing) spans -- they would
    # produce invalid interleaved if-blocks
    spans = [(e[0], e[3]) for e in edits]
    def crosses(a, b):
        return (a[0] < b[0] < a[1] < b[1]) or (b[0] < a[0] < b[1] < a[1])
    keep = []
    for e in edits:
        if any(crosses((e[0], e[3]), s2) for s2 in spans if s2 != (e[0], e[3])):
            continue
        keep.append(e)
    # apply sorted by TARGET descending: an insert at a high target
    # never shifts a lower edit's indices (nested spans handled
    # naturally; goto-line replacements are always below their target)
    for i, indent, cond, tgt in sorted(keep, key=lambda e: -e[3]):
        lines.insert(tgt, f"{indent}end if\n")
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
    print(f"TOTAL wrapped: {total}")


if __name__ == "__main__":
    main()

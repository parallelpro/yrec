#!/usr/bin/env python3
"""goto campaign pass 6: search-loop-with-fallthrough -> exit + guard.

    do v = e1, e2                    do v = e1, e2
       ... goto L ...        ->         ... exit ...
    end do                           end do
    <S: not-found code>              if (v > e2) then
 L  <after>                             <S>
                                     end if
                                  L  <after>
Exactness: with unit step, natural completion leaves v == e2+1 and
exit leaves v <= e2, so `v > e2` reproduces the fall-through/skip
split precisely. Conditions: all refs to L are gotos inside this
loop at any depth; e2 is a bare name or literal not assigned in the
loop body or S; no other labels in S referenced elsewhere; S is
block-balanced with no else at depth 0.
"""
import pathlib
import re
import sys


def code(l):
    return l.split("!")[0]


def transform(path):
    lines = path.read_text().splitlines(keepends=True)
    n = len(lines)
    labels, refs = {}, {}
    for i, l in enumerate(lines):
        c = code(l)
        m = re.match(r"^\s*(\d+)\s", c)
        if m:
            labels[m.group(1)] = i
        for g in re.finditer(r"\bgo\s*to\s+(\d+)", c, re.I):
            refs.setdefault(g.group(1), []).append(i)
        for g in re.finditer(r"\b(?:err|end)\s*=\s*(\d+)", c, re.I):
            refs.setdefault(g.group(1), []).append(-1)  # poison
    # find simple unit-step do loops
    changed = 0
    i = 0
    while i < n:
        c = code(lines[i])
        m = re.match(r"^(\s*)do\s+(\w+)\s*=\s*([\w\d+-]+)\s*,\s*(\w+)\s*$", c)
        if not m or lines[i].lstrip().startswith("!"):
            i += 1
            continue
        indent, var, e1, e2 = m.groups()
        # find matching end do
        depth = 1
        j = i + 1
        while j < n and depth:
            cj = code(lines[j]).strip().lower()
            cjs = re.sub(r"^\d+\s*", "", cj)
            if re.match(r"do\b", cjs) and not cjs.startswith("double"):
                depth += 1
            elif re.match(r"end\s*do\b", cjs):
                depth -= 1
            j += 1
        end_do = j - 1
        # gotos inside body all to same forward label?
        body_gotos = []
        okflag = True
        for k in range(i + 1, end_do):
            for g in re.finditer(r"\bgo\s*to\s+(\d+)", code(lines[k]), re.I):
                body_gotos.append((k, g.group(1)))
        if not body_gotos:
            i = end_do + 1
            continue
        labs = {g[1] for g in body_gotos}
        if len(labs) != 1:
            i = end_do + 1
            continue
        lab = labs.pop()
        tgt = labels.get(lab)
        if tgt is None or tgt <= end_do:
            i = end_do + 1
            continue
        rl = refs.get(lab, [])
        if -1 in rl or any(r < i or r > end_do for r in rl):
            i = end_do + 1
            continue
        # e2 must not be assigned in body or S; var not assigned in S
        span_text = "".join(code(lines[k]) for k in range(i + 1, tgt))
        if not re.match(r"^\d+$", e2):
            if re.search(r"(?<![%\w])" + re.escape(e2) + r"\s*=[^=]", span_text):
                i = end_do + 1
                continue
        if re.search(r"(?<![%\w])" + re.escape(var) + r"\s*=[^=]",
                     "".join(code(lines[k]) for k in range(end_do + 1, tgt))):
            pass  # var assigned in S is fine (guard evaluated first)
        # S structure: balanced, no labels referenced elsewhere, no else@0
        depth = 0
        sok = True
        for k in range(end_do + 1, tgt):
            ck = code(lines[k]).strip().lower()
            if not ck:
                continue
            cks = re.sub(r"^\d+\s*", "", ck)
            if re.match(r"do\b", cks) and not cks.startswith("double"):
                depth += 1
            elif cks.endswith("then") and cks.startswith("if"):
                depth += 1
            elif re.match(r"end\s*(do|if)\b", cks):
                depth -= 1
            elif re.match(r"(else\b|case\b|contains\b|end\s+(subroutine|function))", cks):
                if depth == 0:
                    sok = False
                    break
            if depth < 0:
                sok = False
                break
            mm = re.match(r"^\s*(\d+)\s", code(lines[k]))
            if mm:
                for r in refs.get(mm.group(1), []):
                    if r < end_do or r > tgt:
                        sok = False
        if not sok or depth != 0:
            i = end_do + 1
            continue
        # apply: gotos -> exit; wrap S
        for k, _ in body_gotos:
            lines[k] = re.sub(r"\bgo\s*to\s+" + lab + r"\b", "exit", lines[k],
                              flags=re.I)
        lines.insert(tgt, f"{indent}end if\n")
        lines.insert(end_do + 1, f"{indent}if ({var} > {e2}) then\n")
        changed += len(body_gotos)
        # indices shifted; write and let the caller re-invoke
        path.write_text("".join(lines))
        return changed
    return 0


def main():
    total = 0
    for arg in sys.argv[1:]:
        p = pathlib.Path(arg)
        files = [p] if p.is_file() else sorted(p.rglob("*.f90"))
        for f in files:
            if "test" in f.parts:
                continue
            # iterate per file until no change
            while True:
                k = transform(f)
                if not k:
                    break
                print(f"{f}: {k}")
                total += k
    print(f"TOTAL search-loops: {total}")


if __name__ == "__main__":
    main()

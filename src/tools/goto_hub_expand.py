#!/usr/bin/env python3
"""goto campaign pass 5: expand small terminal hubs at their jump sites.

A hub is a labeled block that (a) every reference jumps FORWARD to,
(b) is at most MAXLEN statements, (c) ends in an unconditional
transfer (return / stop / return 1), and (d) contains no internal
labels or block constructs. Each `go to HUB` (bare or as a one-line
logical if) is replaced by the hub's statements inline (wrapped in an
if-block for the conditional form). The hub block itself stays if it
is reachable by fall-through; if the statement before its label is an
unconditional transfer, the now-unreferenced block is deleted.
"""
import pathlib
import re
import sys

MAXLEN = 8


def code(l):
    return l.split("!")[0]


def transform(path):
    src = path.read_text().splitlines(keepends=True)
    labels, refs = {}, {}
    for i, l in enumerate(src):
        c = code(l)
        m = re.match(r"^\s*(\d+)\s", c)
        if m:
            labels[m.group(1)] = i
        for g in re.finditer(r"\bgo\s*to\s+(\d+)", c, re.I):
            refs.setdefault(g.group(1), []).append(i)
        for g in re.finditer(r"\b(?:err|end)\s*=\s*(\d+)", c, re.I):
            refs.setdefault(g.group(1), []).append(("io", i))
    changed = 0
    for lab, tgt in labels.items():
        rlist = refs.get(lab, [])
        if not rlist or any(isinstance(r, tuple) for r in rlist):
            continue
        if any(r > tgt for r in rlist):
            continue  # only forward hubs
        # collect hub block: from label line until unconditional transfer
        block = []
        j = tgt
        ok = False
        while j < len(src) and len(block) < MAXLEN + 1:
            cj = code(src[j]).rstrip()
            stripped = re.sub(r"^\s*\d+\s*", "", cj).strip()
            if j > tgt and re.match(r"^\s*\d+\s", cj):
                break  # another label -> not a simple hub
            low = stripped.lower()
            if re.match(r"(do\b|if\s*\(.*then|else|end\s*(do|if)|case|contains|end\s+(subroutine|function))", low):
                break
            block.append(stripped)
            if re.match(r"(return\b|stop\b)", low) or low.startswith("return"):
                ok = True
                break
            if low.endswith("&"):
                j += 1
                # absorb continuation lines
                while j < len(src) and (code(src[j]).rstrip().endswith("&")):
                    block.append(code(src[j]).strip())
                    j += 1
                if j < len(src):
                    block.append(code(src[j]).strip())
            j += 1
        if not ok or not block:
            continue
        hub_end = j
        # apply at each ref
        out = []
        local_changed = 0
        for i, l in enumerate(src):
            c = code(l)
            m1 = re.match(r"^(\s*)go\s*to\s+" + lab + r"\s*$", c, re.I)
            m2 = re.match(r"^(\s*)if\s*(\(.*\))\s*go\s*to\s+" + lab + r"\s*$", c, re.I)
            if m1:
                ind = m1.group(1)
                for b in block:
                    out.append(f"{ind}{b}\n")
                local_changed += 1
                continue
            if m2 and m2.group(2).count("(") == m2.group(2).count(")"):
                ind, cond = m2.groups()
                out.append(f"{ind}if {cond} then\n")
                for b in block:
                    out.append(f"{ind}   {b}\n")
                out.append(f"{ind}end if\n")
                local_changed += 1
                continue
            out.append(l)
        if local_changed != len(rlist):
            continue  # some ref form unhandled -> skip whole hub
        src = out
        changed += local_changed
        # rebuild indices for next hub
        labels.clear(); refs.clear()
        for i, l in enumerate(src):
            c = code(l)
            m = re.match(r"^\s*(\d+)\s", c)
            if m:
                labels[m.group(1)] = i
            for g in re.finditer(r"\bgo\s*to\s+(\d+)", c, re.I):
                refs.setdefault(g.group(1), []).append(i)
            for g in re.finditer(r"\b(?:err|end)\s*=\s*(\d+)", c, re.I):
                refs.setdefault(g.group(1), []).append(("io", i))
    if changed:
        path.write_text("".join(src))
    return changed


def main():
    total = 0
    for arg in sys.argv[1:]:
        p = pathlib.Path(arg)
        files = [p] if p.is_file() else sorted(p.rglob("*.f90"))
        for f in files:
            if "test" in f.parts:
                continue
            try:
                k = transform(f)
            except Exception as e:
                print(f"{f}: SKIP ({e})")
                continue
            if k:
                print(f"{f}: {k}")
                total += k
    print(f"TOTAL hub-expanded: {total}")


if __name__ == "__main__":
    main()

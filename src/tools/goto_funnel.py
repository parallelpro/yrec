#!/usr/bin/env python3
"""Expand the stage-3 ierr funnels (go to 900 -> inline block) and
delete the funnel when it becomes unreachable."""
import pathlib
import re
import sys


def code(l):
    return l.split("!")[0]


def transform(path):
    lines = path.read_text().splitlines(keepends=True)
    text = "".join(lines)
    changed_total = 0
    # find each funnel: `  900 continue` ... stop  (allow 900-style labels)
    while True:
        m = re.search(
            r"^(\s*)(\d+)\s+continue\s*\n"
            r"\s*if \(present\(ierr\)\) then\s*\n"
            r"\s*ierr = jerr\s*\n"
            r"\s*return\s*\n"
            r"\s*end if\s*\n"
            r"\s*stop\s*\n", text, re.M)
        if not m:
            break
        lab = m.group(2)
        block = ("if (present(ierr)) then\n"
                 "   ierr = jerr\n"
                 "   return\n"
                 "end if\n"
                 "stop\n")
        # expand refs
        def repl_cond(mm):
            ind, cond = mm.group(1), mm.group(2)
            b = "".join(f"{ind}   {x}\n" for x in block.splitlines())
            return f"{ind}if {cond} then\n{b}{ind}end if\n"
        def repl_bare(mm):
            ind = mm.group(1)
            return "".join(f"{ind}{x}\n" for x in block.splitlines())
        text2, n1 = re.subn(
            r"^(\s*)if\s*(\([^\n]*\))\s*go\s*to\s+" + lab + r"\s*$\n",
            repl_cond, text, flags=re.M)
        text2, n2 = re.subn(
            r"^(\s*)go\s*to\s+" + lab + r"\s*$\n", repl_bare, text2, flags=re.M)
        # delete the funnel block itself (preceded by return -> unreachable)
        text2 = text2.replace(m.group(0), "", 1)
        # remove the now-dangling comment header lines just above? leave.
        if re.search(r"\bgo\s*to\s+" + lab + r"\b", text2):
            print(f"{path}: label {lab} still referenced -- aborting this funnel")
            break
        text = text2
        changed_total += n1 + n2
    if changed_total:
        path.write_text(text)
    return changed_total


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
    print(f"TOTAL funnels expanded: {total}")


if __name__ == "__main__":
    main()

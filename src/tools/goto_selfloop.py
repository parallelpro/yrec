#!/usr/bin/env python3
"""goto campaign pass 7: label if (cond) then ... goto label; end if
-> do while (cond) ... end do.  Also accepts `label continue` followed
immediately by the if-block. Body must contain no other labels and no
other gotos to this label elsewhere in the file."""
import pathlib
import re
import sys


def code(l):
    return l.split("!")[0]


def scope_bounds(lines, i):
    """Return (lo, hi) line range of the procedure containing line i."""
    lo, hi = 0, len(lines)
    for k in range(i, -1, -1):
        if re.match(r"^\s*(pure\s+|elemental\s+|recursive\s+)*(subroutine|(double\s+precision\s+|real\S*\s+|integer\s+|logical\s+)?function)\s+\w", code(lines[k]), re.I):
            lo = k
            break
    for k in range(i, len(lines)):
        if re.match(r"^\s*end\s+(subroutine|function)\b", code(lines[k]), re.I):
            hi = k
            break
    return lo, hi


def transform(path):
    lines = path.read_text().splitlines(keepends=True)
    n = len(lines)
    text = "".join(code(l) for l in lines)
    changed = 0
    i = 0
    while i < n:
        c = code(lines[i])
        m = re.match(r"^(\s*)(\d+)\s+continue\s*$", c)
        start = i
        if m:
            lab, indent = m.group(2), m.group(1)
            j = i + 1
            mif = re.match(r"^(\s*)if\s*(\(.*\))\s*then\s*$", code(lines[j]) if j < n else "")
        else:
            mif = re.match(r"^(\s*)(\d+)\s+if\s*(\(.*\))\s*then\s*$", c)
            if mif:
                indent, lab = mif.group(1), mif.group(2)
                mif = re.match(r"^(\s*)(\d+)\s+(if\s*\(.*\)\s*then)\s*$", c)
                j = i
            else:
                i += 1
                continue
        if not mif:
            i += 1
            continue
        cond = re.search(r"if\s*(\(.*\))\s*then", code(lines[j])).group(1)
        # find matching end if with goto lab immediately before
        depth = 1
        k = j + 1
        goto_line = None
        ok = True
        while k < n and depth:
            ck = code(lines[k]).strip().lower()
            if re.match(r"^\d+", ck):
                ok = False  # no labels inside
                break
            if re.match(r"if\s*\(.*\)\s*then", ck) or re.match(r"(else\s*)?if\s*\(.*\)then", ck):
                if ck.startswith("if"):
                    depth += 1
            if re.match(r"end\s*if\b", ck):
                depth -= 1
                if depth == 0:
                    break
            k += 1
        if not ok or depth != 0:
            i += 1
            continue
        end_if = k
        # last executable stmt before end_if must be `goto lab` at depth1
        last = end_if - 1
        while last > j and not code(lines[last]).strip():
            last -= 1
        if not re.match(r"^\s*go\s*to\s+" + lab + r"\s*$", code(lines[last])):
            i += 1
            continue
        # no other refs to lab in the enclosing procedure
        lo, hi = scope_bounds(lines, i)
        scope_text = "".join(code(l) for l in lines[lo:hi])
        refs = len(re.findall(r"\bgo\s*to\s+" + lab + r"\b", scope_text)) + \
               len(re.findall(r"\b(?:err|end)\s*=\s*" + lab + r"\b", scope_text))
        if refs != 1:
            i += 1
            continue
        # body must not contain other gotos to lab (covered) -- apply
        newblock = [f"{indent}do while {cond}\n"]
        newblock += lines[j + 1:last]
        newblock.append(f"{indent}end do\n")
        lines[start:end_if + 1] = newblock
        n = len(lines)
        text = "".join(code(l) for l in lines)
        changed += 1
        i = start + len(newblock)
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
    print(f"TOTAL self-loops: {total}")


if __name__ == "__main__":
    main()

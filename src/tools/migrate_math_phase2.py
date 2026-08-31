#!/usr/bin/env python3
"""Phase-2 transforms, v2.

powers <files> [--apply]
    Convert ** with a non-integer-literal exponent to pow()/exp10().
    Strings and comments are masked; whitespace around ** handled;
    integer-typed exponents are fine because math_lib's pow is generic
    over (dp,dp) and (dp,integer). Multi-line (continued) exponents
    are reported for hand conversion.

uses <files> [--apply]
    Insert `use math_lib` into every program unit that references an
    elementary math name (exp/log/log10/sin/cos/tan/pow/exp10),
    right before the unit's `implicit none` (YREC's uniform style) or
    after the unit header when there is none. Module-level insertion
    covers contained procedures.

Run from /Applications/YREC/src.
"""
import pathlib
import re
import sys

IDENT = r"[a-zA-Z][a-zA-Z0-9_]*"
MATH_REF = re.compile(r"(?i)(?<![a-z0-9_])(exp|log|log10|sin|cos|tan|pow|exp10)\s*\(")


def mask_line(line):
    """Return code part with strings blanked and comment removed."""
    out = []
    q = None
    for ch in line:
        if q:
            out.append(" ")
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            out.append(" ")
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


def grab_exponent(s, i):
    j = i
    while j < len(s) and s[j] == " ":
        j += 1
    start = j
    if j < len(s) and s[j] in "+-":
        j += 1
    if j < len(s) and s[j] == "(":
        depth = 0
        while j < len(s):
            if s[j] == "(":
                depth += 1
            elif s[j] == ")":
                depth -= 1
                if depth == 0:
                    j += 1
                    break
            j += 1
        return start, j
    m = re.match(rf"({IDENT}|\d+\.?\d*(?:[deDE][+-]?\d+)?)", s[j:])
    if not m:
        return start, start
    j += m.end()
    while j < len(s) and s[j] in "(%":
        if s[j] == "%":
            mm = re.match(IDENT, s[j+1:])
            if not mm:
                break
            j += 1 + mm.end()
        else:
            depth = 0
            while j < len(s):
                if s[j] == "(":
                    depth += 1
                elif s[j] == ")":
                    depth -= 1
                    if depth == 0:
                        j += 1
                        break
                j += 1
    return start, j


def grab_base(s, i):
    j = i
    while j > 0 and s[j-1] == " ":
        j -= 1
    end = j
    while True:
        if j > 0 and s[j-1] == ")":
            depth = 0
            while j > 0:
                j -= 1
                if s[j] == ")":
                    depth += 1
                elif s[j] == "(":
                    depth -= 1
                    if depth == 0:
                        break
        m = re.search(rf"({IDENT}|\d+\.?\d*(?:[deDE][+-]?\d+)?|\.\d+)$", s[:j])
        if m:
            j = m.start()
        if j > 0 and s[j-1] == "%":
            j -= 1
            continue
        break
    return j, end


def powers(paths, apply=False):
    n_conv = n_manual = 0
    for path in paths:
        p = pathlib.Path(path)
        text = p.read_text(errors="replace")
        lines = text.splitlines()
        changed = False
        for li in range(len(lines)):
            line = lines[li]
            code = mask_line(line)
            positions = [m.start() for m in re.finditer(r"\*\*", code)]
            for pos in reversed(positions):
                es, ee = grab_exponent(code, pos + 2)
                etok = line[es:ee]
                if ee == es or line[es:].lstrip().startswith("&") \
                        or code[pos+2:].strip() == "":
                    print(f"{p}:{li+1}: MANUAL (continued/empty exponent): "
                          f"{line.strip()[:70]}")
                    n_manual += 1
                    continue
                if re.fullmatch(r"\(?\s*\d+\s*\)?", etok):
                    continue                       # integer literal: exact
                bs, be = grab_base(code, pos)
                btok = line[bs:be]
                if not btok.strip():
                    print(f"{p}:{li+1}: MANUAL (no base): {line.strip()[:70]}")
                    n_manual += 1
                    continue
                if re.fullmatch(r"10(\.0?[dD]0|\.[dD]0|\.0|\.)?", btok.strip()):
                    repl = f"exp10({etok})"
                else:
                    repl = f"pow({btok}, {etok})"
                line = line[:bs] + repl + line[ee:]
                code = mask_line(line)
                n_conv += 1
                changed = True
            lines[li] = line
        if changed and apply:
            p.write_text("\n".join(lines) + ("\n" if text.endswith("\n") else ""))
    print(f"powers: {n_conv} converted, {n_manual} manual"
          f"{' (applied)' if apply else ' (dry run)'}")


UNIT_RE = re.compile(
    r"(?im)^\s*(?:module|program|(?:pure\s+|elemental\s+|impure\s+|recursive\s+"
    r"|double\s+precision\s+|real(?:\(\w+\))?\s+|integer\s+|logical\s+"
    r"|character(?:\(.*?\))?\s+)*(?:subroutine|function))\s+" + IDENT)
END_UNIT_RE = re.compile(r"(?im)^\s*end\s*(module|subroutine|function|program)?\s*"
                         + f"(?:{IDENT})?\\s*$")


def uses(paths, apply=False):
    for path in paths:
        p = pathlib.Path(path)
        text = p.read_text(errors="replace")
        if "use math_lib" in text:
            continue
        lines = text.splitlines()
        # locate top-level program units (depth 0 tracking via end counts)
        units = []          # (header_idx, kind)
        depth = 0
        for i, l in enumerate(lines):
            code = mask_line(l)
            if UNIT_RE.match(code) and "end" != code.split()[0].lower() \
                    and not re.match(r"(?i)\s*module\s+procedure", code):
                if depth == 0:
                    units.append(i)
                depth += 1
            elif END_UNIT_RE.match(code) and code.strip().lower() != "end if":
                depth = max(0, depth - 1)
        if not units:
            continue
        # which top-level units reference math names?
        spans = [(u, (units[k+1] if k+1 < len(units) else len(lines)))
                 for k, u in enumerate(units)]
        inserts = []
        for (a, b) in spans:
            body = "\n".join(mask_line(l) for l in lines[a:b])
            if MATH_REF.search(body):
                # insert before the unit's first `implicit none`, else
                # right after the header line
                for i in range(a, b):
                    if re.match(r"(?i)\s*implicit\s+none", lines[i]):
                        inserts.append(i)
                        break
                else:
                    j = a
                    while j < b - 1 and mask_line(lines[j]).rstrip().endswith("&"):
                        j += 1
                    inserts.append(j + 1)
        if not inserts:
            continue
        for i in sorted(set(inserts), reverse=True):
            lines.insert(i, "      use math_lib")
        print(f"{p}: use math_lib x{len(set(inserts))}")
        if apply:
            p.write_text("\n".join(lines) + ("\n" if text.endswith("\n") else ""))


if __name__ == "__main__":
    mode = sys.argv[1]
    apply = "--apply" in sys.argv
    args = [a for a in sys.argv[2:] if a != "--apply"]
    if mode == "powers":
        powers(args, apply)
    elif mode == "uses":
        uses(args, apply)

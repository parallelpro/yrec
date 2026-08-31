#!/usr/bin/env python3
"""Reproducibility-campaign source transforms for YREC.

phase1: normalize specific-name intrinsics to generics (dexp -> exp,
        dlog10 -> log10, dlog -> log) so crmath's shadowing can
        intercept them. Token-boundary replace; byte-identical output
        (generics resolve to the same double-precision intrinsics).

phase2 <file>...: report (and with --apply rewrite) real-exponent **
        expressions as pow()/exp10() calls, and say whether the file
        needs `use math_lib`. The exponent classifier: real literal or
        a name declared double precision/real in the same file ->
        convert; integer literal/name -> leave; unknown -> flag for
        hand review. The base extractor walks back over one balanced
        token (identifier with %/() chains, literal, or paren group).

Run from /Applications/YREC/src.
"""
import pathlib
import re
import sys

RENAMES = [("dexp", "exp"), ("dlog10", "log10"), ("dlog", "log")]


def phase1(apply=False):
    total = 0
    for p in sorted(pathlib.Path(".").rglob("*.f90")):
        t = orig = p.read_text(errors="replace")
        for old, new in RENAMES:
            t = re.sub(rf"(?i)(?<![a-z0-9_]){old}(?![a-z0-9_])", new, t)
        if t != orig:
            n = sum(len(re.findall(rf"(?i)(?<![a-z0-9_]){o}(?![a-z0-9_])", orig))
                    for o, _ in RENAMES)
            total += n
            print(f"{p}: {n} renames")
            if apply:
                p.write_text(t)
    print(f"total: {total} renames{' (applied)' if apply else ' (dry run)'}")


IDENT = r"[a-z][a-z0-9_]*"


def decl_types(code):
    """name -> 'real'|'int' from declarations in the file (crude but
    effective for YREC's flat style)."""
    types = {}
    for m in re.finditer(
            r"(?im)^\s*(double\s+precision|real(?:\*8|\(8\)|\(dp\))?|integer(?:\*\d)?)"
            r"[^:!\n]*::\s*([^!\n]+)$", code):
        kind = "int" if m.group(1).lower().startswith("integer") else "real"
        for piece in re.split(r",(?![^()]*\))", m.group(2)):
            mm = re.match(rf"\s*({IDENT})", piece, re.I)
            if mm:
                types[mm.group(1).lower()] = kind
    # old-style: DOUBLE PRECISION a,b(10),c  /  INTEGER i,j
    for m in re.finditer(
            r"(?im)^\s*(double\s+precision|integer(?:\*\d)?|real\*8)\s+([^:!\n]+)$",
            code):
        if "::" in m.group(0) or "function" in m.group(0).lower():
            continue
        kind = "int" if m.group(1).lower().startswith("integer") else "real"
        for piece in re.split(r",(?![^()]*\))", m.group(2)):
            mm = re.match(rf"\s*({IDENT})", piece, re.I)
            if mm:
                types[mm.group(1).lower()] = kind
    return types


def grab_base(s, i):
    """Balanced token ending at s[i-1] (exclusive end i). Returns start
    index, walking back over )..( groups and identifier chains with
    % components and array refs."""
    j = i
    while True:
        j0 = j
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
        m = re.search(rf"({IDENT}|\d[\d.]*(?:[de][+-]?\d+)?)$", s[:j], re.I)
        if m:
            j = m.start()
        if j > 0 and s[j-1] == "%":
            j -= 1
            continue
        if j == j0:
            break
        if not (j > 0 and s[j-1] in ")%"):
            break
    return j


def grab_exponent(s, i):
    """Token starting at s[i]: sign? then literal/identifier(+refs) or
    paren group. Returns (end_index, text)."""
    j = i
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
        return j, s[i:j]
    m = re.match(rf"({IDENT}|\d+\.?\d*(?:[de][+-]?\d+)?(?:_dp)?)", s[j:], re.I)
    if m:
        j += m.end()
        # trailing array ref / component chain
        while j < len(s) and s[j] in "(%":
            if s[j] == "%":
                mm = re.match(IDENT, s[j+1:], re.I)
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
    return j, s[i:j]


def classify_exp(tok, types):
    tok = tok.strip().lstrip("+-").strip()
    if tok.startswith("("):
        inner = tok[1:-1].strip()
        return classify_exp(inner, types) if inner else "unknown"
    if re.fullmatch(r"\d+", tok):
        return "int"
    if re.fullmatch(r"\d*\.\d*([de][+-]?\d+)?|\d+[de][+-]?\d+", tok, re.I):
        return "real"
    m = re.match(rf"({IDENT})", tok, re.I)
    if m:
        return types.get(m.group(1).lower(), "unknown")
    return "unknown"


def phase2(paths, apply=False):
    for path in paths:
        p = pathlib.Path(path)
        text = p.read_text(errors="replace")
        types = decl_types(text)
        lines = text.splitlines()
        changed = False
        for li, line in enumerate(lines):
            if "!" in line:
                code, comment = line.split("!", 1)
                comment = "!" + comment
            else:
                code, comment = line, ""
            out = code
            # right-to-left so indices stay valid
            positions = [m.start() for m in re.finditer(r"\*\*", out)]
            for pos in reversed(positions):
                ei, etok = grab_exponent(out, pos + 2)
                kind = classify_exp(etok, types)
                if kind == "int":
                    continue
                bi = grab_base(out, pos)
                btok = out[bi:pos].strip()
                if kind == "unknown":
                    print(f"{p}:{li+1}: REVIEW  {btok} ** {etok}")
                    continue
                if re.fullmatch(r"10(\.0?d0|\.d0|\.0)?", btok, re.I):
                    repl = f"exp10({etok.strip().lstrip('+')})" \
                        if not etok.strip().startswith("-") \
                        else f"exp10({etok.strip()})"
                else:
                    repl = f"pow({btok}, {etok.strip()})"
                print(f"{p}:{li+1}: {btok}**{etok}  ->  {repl}")
                out = out[:bi] + repl + out[ei:]
                changed = True
            lines[li] = out + comment
        if changed and apply:
            p.write_text("\n".join(lines) + ("\n" if text.endswith("\n") else ""))
            print(f"{p}: written")


if __name__ == "__main__":
    mode = sys.argv[1]
    apply = "--apply" in sys.argv
    args = [a for a in sys.argv[2:] if a != "--apply"]
    if mode == "phase1":
        phase1(apply)
    elif mode == "phase2":
        phase2(args, apply)

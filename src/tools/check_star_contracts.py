#!/usr/bin/env python3
"""Star-contract checker for YREC (2026-10).

A procedure may declare, in the comment block directly above its
`subroutine`/`function` line, which members of the global `star` it
reads and writes:

    !   star reads: jcz, convective_turnover_timescale
    !   star writes: none

Long lists continue on further comment lines indented past the label.
Members of the two input bundles are written ctrl%<name> and job%<name>.
The lists cover the procedure's own body (with its internal procedures,
which share its scope); calls to other procedures are visible in the
code and carry their own contracts.

A write is an assignment to star%<member>, or star%<member> passed as an
actual argument to a dummy declared intent(out) or intent(inout). A
member that is both read and written is listed under writes only.

The check is opt-in: only procedures that carry the two lines are
checked, so contracts can be added driver by driver (docs/GUIDELINES.md,
"the star object").

usage: check_star_contracts.py          check; exit 1 on any mismatch
       check_star_contracts.py --fix    rewrite the declared lists
"""
import collections
import pathlib
import re
import sys

SRC = pathlib.Path(__file__).resolve().parent.parent

PROC_RE = re.compile(
    r"^\s*(?:(?:recursive|pure|elemental|double precision|integer|logical|"
    r"real(?:\*8)?|character(?:\*\d+)?)\s+)*(subroutine|function)\s+(\w+)"
    r"\s*(?:\(([^)]*)\))?", re.I)
END_RE = re.compile(r"^\s*end\s*(subroutine|function)\b", re.I)
STAR_RE = re.compile(r"(?<![\w%])star%(ctrl%\w+|job%\w+|\w+)", re.I)
ASSIGN_RE = re.compile(r"^\s*(?:if\s*\(.*\)\s*)?star%(ctrl%\w+|job%\w+|\w+)"
                       r"(?:\s*\([^=]*\))?(?:%\w+(?:\s*\([^=]*\))?)*\s*=(?!=)", re.I)
CALL_RE = re.compile(r"\bcall\s+(\w+)\s*\((.*)\)\s*$", re.I)
INTENT_RE = re.compile(r"^\s*(?:double precision|integer|logical|real|character|"
                       r"type\s*\(\s*\w+\s*\))[^:]*intent\s*\(\s*(in|out|inout)\s*\)"
                       r"[^:]*::\s*(.*)$", re.I)
READS_RE = re.compile(r"^(\s*!\s*)star reads:\s*(.*)$", re.I)
WRITES_RE = re.compile(r"^(\s*!\s*)star writes:\s*(.*)$", re.I)


def strip_comment(line):
    out, q = [], None
    for ch in line:
        if q:
            out.append(ch)
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out).rstrip()


def logical_lines(lines):
    """(first physical line index, joined code) for each statement."""
    buf, start = "", None
    for i, raw in enumerate(lines):
        code = strip_comment(raw)
        if not code.strip():
            continue
        s = code.strip()
        if buf and s.startswith("&"):
            s = s[1:]
        if start is None:
            start = i
        if s.endswith("&"):
            buf += s[:-1] + " "
            continue
        yield start, buf + s
        buf, start = "", None


def split_args(s):
    depth, cur, parts = 0, "", []
    for ch in s:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch == "," and depth == 0:
            parts.append(cur.strip())
            cur = ""
        else:
            cur += ch
    parts.append(cur.strip())
    return parts


SPELLING = collections.defaultdict(collections.Counter)


def note_spelling(name):
    SPELLING[name.lower()][name] += 1
    return name.lower()


def spelled(name):
    """The code's most common spelling of a member (Fortran ignores
    case; the lists keep the spelling readers see in the code)."""
    forms = SPELLING.get(name)
    return forms.most_common(1)[0][0] if forms else name


def scan(files):
    """Every procedure: file, statement line, dummies, intents, star
    refs/writes in its own scope, its parent (for internal procedures),
    and its calls."""
    procs = {}
    for path in files:
        lines = path.read_text(errors="replace").splitlines()
        stack = []
        for idx, code in logical_lines(lines):
            m = PROC_RE.match(code)
            if m and not re.match(r"^\s*end\b", code, re.I):
                name = m.group(2).lower()
                p = dict(name=name, path=path, line=idx,
                         dummies=[a.strip().lower() for a in (m.group(3) or "").split(",") if a.strip()],
                         intents={}, reads=set(), writes=set(), calls=[],
                         parent=stack[-1]["name"] if stack else None)
                procs.setdefault(name, p)
                stack.append(p)
                continue
            if END_RE.match(code):
                if stack:
                    stack.pop()
                continue
            if not stack:
                continue
            p = stack[-1]
            mi = INTENT_RE.match(code)
            if mi:
                for n in re.sub(r"\([^)]*\)", "", mi.group(2)).split(","):
                    n = n.split("=")[0].strip().lower()
                    if n:
                        p["intents"][n] = mi.group(1).lower()
                continue
            ma = ASSIGN_RE.match(code)
            if ma:
                p["writes"].add(note_spelling(ma.group(1)))
                rhs = code[code.index("=", ma.end() - 1) + 1:]
                for f in STAR_RE.findall(rhs):
                    p["reads"].add(note_spelling(f))
                lhs = code[:ma.end() - 1]
                # star refs inside the target's subscripts are reads
                for f in STAR_RE.findall(lhs)[1:]:
                    p["reads"].add(note_spelling(f))
                continue
            mc = CALL_RE.search(code)
            if mc:
                p["calls"].append((mc.group(1).lower(), split_args(mc.group(2))))
            for f in STAR_RE.findall(code):
                p["reads"].add(note_spelling(f))
    return procs


def footprint(procs, name):
    """reads, writes of a procedure's scope (itself and its internal
    procedures), with star% actual arguments to intent(out/inout)
    dummies counted as writes."""
    scope = [p for p in procs.values() if p["name"] == name or _ancestor(procs, p, name)]
    reads, writes = set(), set()
    for p in scope:
        reads |= p["reads"]
        writes |= p["writes"]
        for callee, args in p["calls"]:
            q = procs.get(callee)
            for i, a in enumerate(args):
                kw = re.match(r"^(\w+)\s*=\s*(.*)$", a)
                dummy = None
                if kw and not a.lower().startswith("star%"):
                    dummy, a = kw.group(1).lower(), kw.group(2)
                ma = re.match(r"^star%(ctrl%\w+|job%\w+|\w+)", a, re.I)
                if not ma or q is None:
                    continue
                if dummy is None:
                    dummy = q["dummies"][i] if i < len(q["dummies"]) else None
                if q["intents"].get(dummy) in ("out", "inout"):
                    writes.add(ma.group(1).lower())
    reads -= writes
    return reads, writes


def _ancestor(procs, p, name):
    while p["parent"]:
        if p["parent"] == name:
            return True
        p = procs.get(p["parent"])
        if p is None:
            return False
    return False


def declared(lines, stmt_line):
    """The declared contract above a statement: (reads, writes, span)
    or None. span is the (first, last) line of the two entries."""
    i = stmt_line - 1
    block = []
    while i >= 0 and lines[i].lstrip().startswith("!"):
        block.append(i)
        i -= 1
    block.reverse()
    r_at = next((j for j in block if READS_RE.match(lines[j])), None)
    w_at = next((j for j in block if WRITES_RE.match(lines[j])), None)
    if r_at is None or w_at is None:
        return None

    def collect(at):
        m = (READS_RE.match(lines[at]) or WRITES_RE.match(lines[at]))
        indent = len(m.group(1))
        text = [m.group(2)]
        j = at + 1
        while j < stmt_line and j in block and not READS_RE.match(lines[j]) \
                and not WRITES_RE.match(lines[j]):
            cm = re.match(r"^(\s*!\s*)(\S.*)$", lines[j])
            if not cm or len(cm.group(1)) <= indent:
                break
            text.append(cm.group(2))
            j += 1
        names = {n.strip().lower() for n in ",".join(text).split(",") if n.strip()}
        names.discard("none")
        return names, j - 1

    reads, r_end = collect(r_at)
    writes, w_end = collect(w_at)
    return reads, writes, (min(r_at, w_at), max(r_end, w_end))


def render(label, names, prefix="!   ", width=72):
    items = sorted((spelled(n) for n in names), key=str.lower) or ["none"]
    out, cur = [], prefix + label + ":"
    for k, n in enumerate(items):
        piece = " " + n + ("," if k < len(items) - 1 else "")
        if len(cur) + len(piece) > width:
            out.append(cur)
            cur = prefix + "  " + piece.lstrip()
        else:
            cur += piece
    out.append(cur)
    return out


def main():
    fix = "--fix" in sys.argv[1:]
    files = sorted(f for f in SRC.rglob("*.f90") if "/test/" not in f.as_posix())
    procs = scan(files)
    problems, fixes = [], collections.defaultdict(list)
    for p in procs.values():
        lines = p["path"].read_text(errors="replace").splitlines()
        d = declared(lines, p["line"])
        if d is None:
            continue
        reads, writes = footprint(procs, p["name"])
        dr, dw, span = d
        rendered = render("star reads", reads) + render("star writes", writes)
        if (dr, dw) != (reads, writes) or (fix and lines[span[0]:span[1] + 1] != rendered):
            rel = p["path"].relative_to(SRC)
            msg = [f"{rel}: {p['name']}"]
            for label, have, want in (("reads", dr, reads), ("writes", dw, writes)):
                if want - have:
                    msg.append(f"    {label} missing: {', '.join(spelled(n) for n in sorted(want - have))}")
                if have - want:
                    msg.append(f"    {label} not used: {', '.join(spelled(n) for n in sorted(have - want))}")
            problems.append("\n".join(msg))
            fixes[p["path"]].append((span, rendered))
    if fix:
        for path, edits in fixes.items():
            lines = path.read_text(errors="replace").splitlines()
            for (a, b), new in sorted(edits, reverse=True):
                lines[a:b + 1] = new
            path.write_text("\n".join(lines) + "\n")
        print(f"rewrote {sum(len(e) for e in fixes.values())} contract(s)")
        return 0
    if problems:
        print("STAR CONTRACT MISMATCHES (update the header, or run "
              "tools/check_star_contracts.py --fix):")
        print("\n".join(problems))
        return 1
    n = sum(1 for p in procs.values()
            if declared(p["path"].read_text(errors="replace").splitlines(), p["line"]))
    print(f"Star contracts OK ({n} declared).")
    return 0


if __name__ == "__main__":
    sys.exit(main())

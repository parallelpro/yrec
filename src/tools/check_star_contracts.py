#!/usr/bin/env python3
"""Star-contract checker for YREC (2026-10).

A procedure may declare, in the comment block directly above its
`subroutine`/`function` line, which members of the star it reads and
writes -- through the global `star`, or through the argument `s`
(type(star_info)) that core's procedures receive it as:

    !   star reads: envelope_cz_bottom_index_start, convective_turnover_timescale
    !   star writes: none

Long lists continue on further comment lines indented past the label.
Members of the two input bundles are written ctrl%<name> and job%<name>.
The lists cover the procedure's own body (with its internal procedures,
which share its scope); calls to other procedures are visible in the
code and carry their own contracts.

A write is an assignment to star%<member>, a read statement's or an
allocate's target, or star%<member> passed as an actual argument to a
dummy declared intent(out) or intent(inout). A member that is both read
and written is listed under writes only.

Procedures are told apart by module (and host): two modules may each
have a wind_now, and a call resolves as the compiler does: an internal
procedure of the caller or its hosts, the caller's module, a module the
scope uses (by name in an only-list, or whole), then the one procedure
of that name in the tree.

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
MODULE_RE = re.compile(r"^\s*module\s+(\w+)\s*$", re.I)
USE_RE = re.compile(r"^\s*use\s*(?:,\s*(?:non_)?intrinsic\s*::\s*|::\s*)?(\w+)\s*(?:,\s*only\s*:\s*(.*))?$",
                    re.I)
READ_STMT_RE = re.compile(r"^\s*(?:if\s*\(.*\)\s*)?read\s*\(", re.I)
ALLOC_RE = re.compile(r"^\s*(?:if\s*\(.*\)\s*)?(?:allocate|deallocate)\s*\(", re.I)
INTENT_RE = re.compile(r"^\s*(?:double precision|integer|logical|real|character|"
                       r"type\s*\(\s*\w+\s*\))[^:]*intent\s*\(\s*(in|out|inout)\s*\)"
                       r"[^:]*::\s*(.*)$", re.I)
READS_RE = re.compile(r"^(\s*!\s*)star reads:\s*(.*)$", re.I)
# The star as an argument: a dummy s of type(star_info), whose s%member
# references count as star%member in its scope (and its internal
# procedures').
STAR_DUMMY_RE = re.compile(r"^\s*type\s*\(\s*star_info\s*\)[^:]*::\s*s\s*$", re.I)
S_RE = re.compile(r"(?<![\w%])s%", re.I)
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


def matching_paren(s, i):
    depth, q = 0, None
    for j in range(i, len(s)):
        ch = s[j]
        if q:
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
            if depth == 0:
                return j
    return len(s) - 1


def io_targets(items, name_re):
    """the names a read statement's io-list assigns: its items, and the
    items of its implied-do loops (not their bounds)"""
    out = []
    for part in split_args(items):
        part = part.strip()
        if part.startswith("(") and matching_paren(part, 0) == len(part) - 1:
            inner = split_args(part[1:-1])
            j = next((k for k, x in enumerate(inner) if re.match(r"^\w+\s*=", x)), len(inner))
            out += io_targets(",".join(inner[:j]), name_re)
            continue
        m = name_re.match(part)
        if m:
            out.append(m.group(1))
    return out


def parse_use(code):
    """(module, {local name: remote name} or None for a whole use)"""
    m = USE_RE.match(code)
    if not m:
        return None
    if m.group(2) is None:
        return m.group(1).lower(), None
    names = {}
    for item in m.group(2).split(","):
        item = item.strip().lower()
        if "=>" in item:
            local, remote = (x.strip() for x in item.split("=>", 1))
            names[local] = remote
        elif item:
            names[item] = item
    return m.group(1).lower(), names


def scan(files):
    """Every procedure, keyed by (module, name) -- an internal procedure
    by its host's key and its name: file, statement line, dummies,
    intents, star refs/writes in its own scope, its host (parent), the
    modules its scope uses, and its calls."""
    procs = {}
    for path in files:
        lines = path.read_text(errors="replace").splitlines()
        stack, module, module_uses = [], path.stem.lower(), []
        for idx, code in logical_lines(lines):
            mm = MODULE_RE.match(code)
            if mm and not stack:
                module, module_uses = mm.group(1).lower(), []
                continue
            m = PROC_RE.match(code)
            if m and not re.match(r"^\s*end\b", code, re.I):
                name = m.group(2).lower()
                key = (stack[-1]["key"] if stack else (module,)) + (name,)
                p = dict(name=name, key=key, module=module, path=path, line=idx,
                         dummies=[a.strip().lower() for a in (m.group(3) or "").split(",") if a.strip()],
                         intents={}, reads=set(), writes=set(), calls=[], uses=[],
                         module_uses=module_uses, parent=stack[-1]["key"] if stack else None)
                procs.setdefault(key, p)
                stack.append(p)
                continue
            if END_RE.match(code):
                if stack:
                    stack.pop()
                continue
            u = parse_use(code)
            if u:
                (stack[-1]["uses"] if stack else module_uses).append(u)
                continue
            if not stack:
                continue
            p = stack[-1]
            if STAR_DUMMY_RE.match(code):
                p["star_s"] = True
            if any(q.get("star_s") for q in stack):
                code = S_RE.sub("star%", code)
            mi = INTENT_RE.match(code)
            if mi:
                for n in re.sub(r"\([^)]*\)", "", mi.group(2)).split(","):
                    n = n.split("=")[0].strip().lower()
                    if n:
                        p["intents"][n] = mi.group(1).lower()
                continue
            mr = READ_STMT_RE.match(code)
            if mr:
                j = matching_paren(code, mr.end() - 1)
                targets = io_targets(code[j + 1:], re.compile(r"^star%(ctrl%\w+|job%\w+|\w+)", re.I))
                for f in targets:
                    p["writes"].add(note_spelling(f))
                for f in STAR_RE.findall(code):
                    if f.lower() not in {x.lower() for x in targets}:
                        p["reads"].add(note_spelling(f))
                continue
            if ALLOC_RE.match(code):
                for f in STAR_RE.findall(code):
                    p["writes"].add(note_spelling(f))
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


_BY_NAME = {}


def resolve(procs, caller, name):
    """the procedure a call from caller names, as the compiler finds it:
    scope by scope from the caller out to its module, a procedure the
    scope contains, then one a module it uses provides (by name in an
    only-list, then a whole use); else the only top-level procedure of
    that name in the tree (a facade's re-export)"""
    def from_uses(uses):
        for mod, only in uses:
            if only is not None and name in only:
                r = procs.get((mod, only[name]))
                if r is not None:
                    return r
        for mod, only in uses:
            if only is None:
                r = procs.get((mod, name))
                if r is not None:
                    return r
        return None
    q = caller
    while q is not None:
        r = procs.get(q["key"] + (name,)) or from_uses(q["uses"])
        if r is not None:
            return r
        q = procs.get(q["parent"]) if q["parent"] else None
    r = procs.get((caller["module"], name)) or from_uses(caller["module_uses"])
    if r is not None:
        return r
    if id(procs) not in _BY_NAME:
        by = collections.defaultdict(list)
        for p in procs.values():
            if p["parent"] is None:
                by[p["name"]].append(p)
        _BY_NAME[id(procs)] = by
    top = _BY_NAME[id(procs)].get(name, [])
    return top[0] if len(top) == 1 else None


def footprint(procs, key):
    """reads, writes of a procedure's scope (itself and its internal
    procedures), with star% actual arguments to intent(out/inout)
    dummies counted as writes. key is the procedure's key in procs."""
    scope = [p for p in procs.values() if p["key"][:len(key)] == key]
    reads, writes = set(), set()
    for p in scope:
        reads |= p["reads"]
        writes |= p["writes"]
        for callee, args in p["calls"]:
            q = resolve(procs, p, callee)
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
        reads, writes = footprint(procs, p["key"])
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

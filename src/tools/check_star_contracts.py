#!/usr/bin/env python3
"""Star-contract checker for YREC (2026-10).

A procedure that takes the star as its argument `s` (type(star_info);
there is no global star) declares, in the comment block directly above
its `subroutine`/`function` line, which members of the star it reads and
writes:

    !   star reads: envelope_cz_bottom_index_start, convective_turnover_timescale
    !   star writes: none

Long lists continue on further comment lines indented past the label.
Members of the two controls records are written control%<name> and
physics%<name>.
The lists cover the procedure's own body (with its internal procedures,
which share its scope); calls to other procedures are visible in the
code and carry their own contracts.

A write is an assignment to s%<member>, a read statement's or an
allocate's target, or s%<member> passed as an actual argument to a
dummy declared intent(out) or intent(inout), whether or not the callee
assigns it. (gen_star_docs.py, which builds state/README's member
tables, counts an intent(inout) argument as a write only when the
callee assigns the dummy, directly or through its own callees, and as a
read only when the callee reads it; so a member passed to an inout dummy
the callee only reads is a write here and a read there.) A member that
is both read and written is listed under writes only.

Statements are read as the compiler reads them: continuation lines are
joined, a line is split at each ';' outside strings (so `case (1);
s%x = v` is an assignment), statement labels are dropped and string
literals blanked. A procedure's first line may carry any prefixes and a
type-spec: `pure real(dp) function`, `type(t) function`. An
if-statement's condition is read, as are a target's subscripts.

Procedures are told apart by module (and host): two modules may each
have a wind_now, and a call resolves as the compiler does, scope by
scope from the caller out to its module: a procedure the scope
contains, then one a module the scope uses provides (by name in an
only-list, then a whole use); failing those, the one top-level
procedure of that name in the tree (a facade's re-export).

Where they are required: every procedure that takes the star in the
directories of REQUIRED_DIRS (all of core/, setup, io and rotation:
every directory whose code takes the star) declares the two lines, and
the check fails on one that does not; a file in CONTRACT_EXEMPT is left
out, with the reason. Elsewhere a contract is optional, and checked when
present (docs/GUIDELINES.md, "The star object: six rules").

usage: check_star_contracts.py              check; exit 1 on any mismatch
                                            (make check and test_boundaries.py
                                            run it)
       check_star_contracts.py --fix        rewrite the declared lists
       check_star_contracts.py --self-test  the scan on synthetic cases
                                            (test_boundaries.py runs it)
"""
import collections
import pathlib
import re
import sys

SRC = pathlib.Path(__file__).resolve().parent.parent

# A procedure's first statement: any of the prefixes (recursive, pure,
# elemental, impure, module) and a type-spec -- double precision, real(dp),
# real*8, integer, logical, complex, character(...), type(...), class(...)
# -- before subroutine or function.
TYPE_SPEC = (r"(?:double\s+precision|integer|logical|complex|real|character|type|class)"
             r"(?:\s*\*\s*\d+|\s*\([^()]*\))?")
PROC_RE = re.compile(
    r"^\s*(?:(?:recursive|pure|elemental|impure|module|" + TYPE_SPEC + r")\s+)*"
    r"(subroutine|function)\s+(\w+)\s*(?:\(([^)]*)\))?", re.I)
END_RE = re.compile(r"^\s*end\s*(subroutine|function)\b", re.I)
STAR_RE = re.compile(r"(?<![\w%])star%(control%\w+|physics%\w+|\w+)", re.I)
ASSIGN_RE = re.compile(r"^\s*(?:if\s*\(.*\)\s*)?star%(control%\w+|physics%\w+|\w+)"
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


# A statement label (a number before the statement), which the patterns
# below would otherwise take for the start of the statement.
LABEL_RE = re.compile(r"^\s*\d+\s+(?=\S)")


def mask_strings(code):
    """code with the characters inside string literals blanked, so that a
    message mentioning s%x or a ';' inside a string is not code"""
    out, q = [], None
    for ch in code:
        if q:
            out.append(ch if ch == q else " ")
            if ch == q:
                q = None
        else:
            if ch in "'\"":
                q = ch
            out.append(ch)
    return "".join(out)


def split_statements(code):
    """the statements of one line: split on ';' outside strings, each
    without its statement label (`case (1);  s%x = v` is two statements,
    `case (1)` and `s%x = v`)"""
    parts, cur, q = [], "", None
    for ch in code:
        if q:
            cur += ch
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            cur += ch
        elif ch == ";":
            parts.append(cur)
            cur = ""
        else:
            cur += ch
    parts.append(cur)
    return [LABEL_RE.sub("", p).strip() for p in parts if p.strip()]


def logical_lines(lines):
    """(first physical line index, code) for each statement: continuation
    lines joined, then split into statements (split_statements)."""
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
        for stmt in split_statements(mask_strings(buf + s)):
            yield start, stmt
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
                targets = io_targets(code[j + 1:], re.compile(r"^star%(control%\w+|physics%\w+|\w+)", re.I))
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
                # star refs in an if-statement's condition and inside the
                # target's subscripts are reads
                for f in STAR_RE.findall(code[:ma.start(1)] + " " + code[ma.end(1):ma.end() - 1]):
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
                ma = re.match(r"^star%(control%\w+|physics%\w+|\w+)", a, re.I)
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


# The directories where every procedure that takes the star declares its
# contract: all of core/ (the run and step drivers, the model, the
# structure, the envelope and the diagnostics), the settings builders and
# installs, io's readers and writers, and rotation's drivers, which keep
# the star for its state (relative to src/, tests aside). The physics
# domains never take the star (check_boundaries' STAR_FREE_DIRS), and
# state/ declares it without procedures that take it.
REQUIRED_DIRS = ("core", "setup", "io", "rotation")
# Files of those directories whose procedures declare none, and why.
CONTRACT_EXEMPT = {
    "io/controls_sync.f90": "generated by tools/gen_controls_state.py; its procedures copy every "
                            "control between the buffer and the two records, and every run copy",
}


def missing_contracts(procs, root):
    """the procedures in REQUIRED_DIRS (paths relative to root) that take
    the star and declare no contract"""
    out = []
    for p in procs.values():
        rel = p["path"].relative_to(root).as_posix()
        if not p.get("star_s") or rel in CONTRACT_EXEMPT \
                or not any(rel.startswith(d + "/") for d in REQUIRED_DIRS):
            continue
        if declared(p["path"].read_text(errors="replace").splitlines(), p["line"]) is None:
            out.append(f"{rel}: {'/'.join(p['key'][1:])} takes the star and declares no contract "
                       "(add '!   star reads: none' and '!   star writes: none' above it, then --fix)")
    return sorted(out)


# The self-test (--self-test, run by test_boundaries.py): the scan on a
# synthetic module whose procedures write the forms the patterns must see.
# Each case is (procedure, the reads and writes the scan must find).
SELF_TEST_SOURCE = """\
module star_tool_cases
      use star_info_lib, only: star_info
      use, intrinsic :: iso_fortran_env, only: dp => real64
contains
subroutine case_labels(s, k, val)
      type(star_info), intent(inout) :: s
      integer, intent(in) :: k
      real(dp), intent(in) :: val
      select case (k)
      case (1);      s%case_one = val
      case default;  s%case_other = val
      end select
end subroutine case_labels
subroutine statement_parts(s)
      type(star_info), intent(inout) :: s
      if (s%flag_in) s%part_a = 1; s%part_b = s%part_c
   10 s%labelled = 2
      write(*,*) 'a; s%not_code = 3'
end subroutine statement_parts
real(dp) function typed_real(s)
      type(star_info), intent(in) :: s
      typed_real = s%real_read
end function typed_real
pure type(probe_t) function typed_record(s) result(r)
      type(star_info), intent(in) :: s
      r%x = s%record_read
end function typed_record
elemental real(dp) function no_star(x)
      real(dp), intent(in) :: x
      no_star = x
end function no_star
end module star_tool_cases
"""
SELF_TEST_CASES = {
    "case_labels": (set(), {"case_one", "case_other"}),
    "statement_parts": ({"flag_in", "part_c"}, {"part_a", "part_b", "labelled"}),
    "typed_real": ({"real_read"}, set()),
    "typed_record": ({"record_read"}, set()),
    "no_star": (set(), set()),
}


def self_test():
    """ok/FAIL per case; returns the number of failures"""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        path = pathlib.Path(d) / "star_tool_cases.f90"
        path.write_text(SELF_TEST_SOURCE)
        procs = scan([path])
        # the coverage rule: a star-taking procedure without a contract is
        # reported in a required directory (core/run, core/model: all of core/;
        # rotation/), not elsewhere (state/) or in an exempt file
        bare = SELF_TEST_SOURCE.replace("star_tool_cases", "bare_cases")
        for rel in ("core/run/bare_cases.f90", "core/model/model_cases.f90", "rotation/rotation_cases.f90",
                    "state/state_cases.f90", "io/controls_sync.f90"):
            (pathlib.Path(d) / rel).parent.mkdir(parents=True, exist_ok=True)
            (pathlib.Path(d) / rel).write_text(bare.replace("bare_cases", pathlib.Path(rel).stem))
        (pathlib.Path(d) / "core/run/declared_cases.f90").write_text(
            bare.replace("bare_cases", "declared_cases").replace(
                "subroutine case_labels(", "!   star reads: none\n!   star writes: case_one, case_other\n"
                "subroutine case_labels("))
        found = missing_contracts(scan(sorted(f for f in pathlib.Path(d).rglob("*.f90") if f != path)),
                                  pathlib.Path(d))
    want = sorted(f"{f}: {n}" for f in ("core/run/bare_cases.f90", "core/run/declared_cases.f90",
                                         "core/model/model_cases.f90", "rotation/rotation_cases.f90")
                  for n in ("case_labels", "statement_parts", "typed_real", "typed_record")
                  if not (f == "core/run/declared_cases.f90" and n == "case_labels"))
    got = sorted(m.split(" takes")[0] for m in found)
    failures = int(got != want)
    print(f"{'ok  ' if got == want else 'FAIL'} required contracts: {len(got)} reported"
          + ("" if got == want else f" {got} (want {want})"))
    for name, (reads, writes) in SELF_TEST_CASES.items():
        key = ("star_tool_cases", name)
        got = footprint(procs, key) if key in procs else None
        good = got == (reads, writes)
        failures += not good
        print(f"{'ok  ' if good else 'FAIL'} {name}: "
              + (f"reads {sorted(got[0])}, writes {sorted(got[1])}" if got else "not found as a procedure")
              + ("" if good else f" (want reads {sorted(reads)}, writes {sorted(writes)})"))
    return failures


def main():
    if "--self-test" in sys.argv[1:]:
        return 1 if self_test() else 0
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
    missing = missing_contracts(procs, SRC)
    if missing:
        print(f"STAR CONTRACTS MISSING (required in {', '.join(REQUIRED_DIRS)}):")
        print("\n".join("  " + m for m in missing))
        return 1
    n = sum(1 for p in procs.values()
            if declared(p["path"].read_text(errors="replace").splitlines(), p["line"]))
    print(f"Star contracts OK ({n} declared).")
    return 0


if __name__ == "__main__":
    sys.exit(main())

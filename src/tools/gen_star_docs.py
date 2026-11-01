#!/usr/bin/env python3
"""Generate state/README.md's lifetime-class table and member tables.

Source of truth: the star_info type in state/star_info_lib.f90. Its
members are grouped into lifetime classes, each opened by a banner

    ! ==== <class> ====
    ! <what the class holds> ... Valid: <when its values can be trusted>

and each member has a comment directly above its declaration (one
member per declaration):

    ! <meaning> [<unit>]; <writer>, <writer>. <notes>

The writers named there are modules (or procedures) that write the
member; the README lists every writer and reader, from the code:

- written by: the procedures in which s%member is assigned (or a
  component of it), is a read statement's target or is allocated, or is
  passed to a dummy argument that is intent(out), or intent(inout) (or
  without an intent) and that the callee assigns, directly, in its
  internal procedures or through its own callees. A member written
  inside a callee therefore appears under its caller, when the caller
  passes s%member explicitly. More than six procedures: the modules only.
- read by: the modules where the value of s%member is used: in an
  expression, a subscript or a loop bound, or passed to a dummy that the
  callee reads (an intent(out) dummy, or an intent(inout) one the callee
  never reads, is not a read). More than seven modules: a count by
  directory.

The analysis reuses check_star_contracts.py's statement scan (joined
continuation lines, statements split at ';' outside strings, labels
dropped, string literals blanked, any type-spec before `function`, the
star argument s, the declared intents). Tests (*/test/) and retired code
are not scanned.

Fails, writing nothing, when a member has no comment in that form, names
a writer that does not write it, or is never written or never read
outside the tests; when the first column of README section 4's table
(the members removed) names a live member; or when a member breaks a
class rule:
- an earlier value, named <member>_start, <member>_prev or
  <member>_iteration_start, is in the Step snapshot class, and that
  class holds nothing else;
- any other member that one module alone writes and reads is in the
  class "Private to one module", and that class holds nothing else.

usage: gen_star_docs.py              rewrite state/README.md's generated parts
       gen_star_docs.py --check      exit 1 if they are stale
                                     (test_generated_files.py runs this)
       gen_star_docs.py --self-test  the writer and reader analysis on
                                     synthetic cases (test_boundaries.py)
"""
import collections
import pathlib
import re
import sys

import check_star_contracts as csc
from generated_files import emit, finish

SRC = pathlib.Path(__file__).resolve().parent.parent
LIB = SRC / "state" / "star_info_lib.f90"
README = SRC / "state" / "README.md"

# ---------------------------------------------------------------------
# star_info's classes and member comments

CLASS_RE = re.compile(r"^!\s*====\s*(.+?)\s*=+\s*$")
MEMBER_DECL_RE = re.compile(r"^\s*(double precision|integer|logical|character\s*\(\s*len\s*=\s*\d+\s*\)|"
                            r"type\s*\(\s*(\w+)\s*\))\s*::\s*(.*)$", re.I)
DOC_RE = re.compile(r"^(?P<meaning>[^\[]+?)\s*\[(?P<unit>[^\]]*)\];\s*(?P<writers>\w+(?:\s*,\s*\w+)*)"
                    r"(?:\.\s*(?P<notes>.*?))?\.?$", re.S)


def parse_star_info(path):
    """(classes, members, errors): each class with its name, description
    and members; each member with its declaration and its comment's
    meaning, unit, writers and notes."""
    lines = path.read_text().splitlines()
    i0 = next(i for i, l in enumerate(lines) if re.match(r"\s*type, public :: star_info\b", l))
    i1 = next(i for i, l in enumerate(lines) if i > i0 and re.match(r"\s*end type star_info\b", l))
    classes, members, errors = [], [], []
    comment, cls, in_banner = [], None, False
    for k in range(i0 + 1, i1):
        s = lines[k].strip()
        m = CLASS_RE.match(s)
        if m:
            cls = dict(name=m.group(1), text=[], members=[])
            classes.append(cls)
            in_banner, comment = True, []
            continue
        if not s:
            if comment and not in_banner:
                errors.append(f"star_info_lib.f90:{k}: a comment that is not directly above a declaration")
            in_banner, comment = False, []
            continue
        if s.startswith("!"):
            (cls["text"] if in_banner and cls else comment).append(s[1:].strip())
            continue
        d = MEMBER_DECL_RE.match(s)
        if not d or s.endswith("&"):
            errors.append(f"star_info_lib.f90:{k + 1}: not a one-line member declaration: {s}")
            comment = []
            continue
        decl = d.group(3)
        if "," in re.sub(r"\([^()]*\)", "", decl.split("=")[0]):
            errors.append(f"star_info_lib.f90:{k + 1}: declare one member per line: {s}")
        name = re.match(r"\w+", decl).group(0)
        doc = " ".join(comment)
        md = DOC_RE.match(doc)
        e = dict(name=name, shape=re.sub(r"\s*=.*$", "", decl).strip(), record=d.group(2),
                 cls=cls["name"] if cls else None, meaning=doc, unit="", writers=[], notes="")
        if cls is None:
            errors.append(f"{name}: declared before the first class banner")
        if md:
            e.update(meaning=md.group("meaning").strip(), unit=md.group("unit").strip(),
                     writers=[w.strip() for w in md.group("writers").split(",")],
                     notes=(md.group("notes") or "").strip())
        else:
            errors.append(f"{name}: its comment must read '<meaning> [<unit>]; <writers>. <notes>'"
                          f" (found {doc!r})")
        members.append(e)
        if cls:
            cls["members"].append(e)
        comment = []
    for c in classes:
        text = " ".join(c["text"])
        what, _, valid = text.partition("Valid:")
        c["what"], c["valid"] = what.strip(), valid.strip()
        if not valid:
            errors.append(f"class {c['name']}: its banner needs a 'Valid:' sentence")
    return classes, members, errors


# the class rules
SNAPSHOT_CLASS = "Step snapshot"
PRIVATE_CLASS = "Private to one module"
SNAPSHOT_NAME_RE = re.compile(r"_(start|prev|iteration_start)$", re.I)
# members named like earlier values that are not (none since the
# 2026-10-29 renames; a name listed here is exempt from the rule)
NAME_RULE_PENDING = set()


def class_rule_errors(members, writers, readers):
    out = []
    for e in members:
        key = e["name"].lower()
        mods = sorted(set(writers[key]) | readers[key])
        snapshot = bool(SNAPSHOT_NAME_RE.search(e["name"])) and key not in NAME_RULE_PENDING
        if snapshot and e["cls"] != SNAPSHOT_CLASS:
            out.append(f"{e['name']}: an earlier value (<member>_start, _prev, _iteration_start) "
                       f"belongs in the class '{SNAPSHOT_CLASS}', not '{e['cls']}'")
        elif e["cls"] == SNAPSHOT_CLASS and not snapshot:
            out.append(f"{e['name']}: the class '{SNAPSHOT_CLASS}' holds earlier values only "
                       f"(<member>_start, _prev, _iteration_start)")
        elif not snapshot and len(mods) == 1 and e["cls"] != PRIVATE_CLASS:
            out.append(f"{e['name']}: only {mods[0]} uses it, so it belongs in the class "
                       f"'{PRIVATE_CLASS}' (or is a local of {mods[0]}), not '{e['cls']}'")
        elif e["cls"] == PRIVATE_CLASS and len(mods) > 1:
            out.append(f"{e['name']}: in the class '{PRIVATE_CLASS}', but {', '.join(mods)} use it")
    return out


# ---------------------------------------------------------------------
# who writes and who reads each member

MODULE_RE = re.compile(r"^\s*module\s+(\w+)\s*$", re.I)
READ_STMT_RE = re.compile(r"^\s*(?:if\s*\(.*\)\s*)?read\s*\(", re.I)
ALLOC_RE = re.compile(r"^\s*(?:if\s*\(.*\)\s*)?(?:allocate|deallocate)\s*\(", re.I)
DECL_RE = re.compile(r"^\s*(double precision|integer|logical|real|character|complex|"
                     r"type\s*\(\s*\w+\s*\))\b[^:]*::\s*(.*)$", re.I)
SKIP_RE = re.compile(r"^\s*(use\b|implicit\b|save\b|data\b|parameter\b|external\b|"
                     r"intrinsic\b|public\b|private\b|contains\b|\d*\s*format\b|interface\b|"
                     r"end\s*interface\b|procedure\b|module\s+procedure\b)", re.I)
STAR_PAT = r"star%(\w+)"
# an actual argument that is a variable (with subscripts or components),
# not an expression
WHOLE_REF = r"(\s*\([^()]*(\([^()]*\)[^()]*)*\))?(%\w+(\s*\([^()]*\))?)*$"


matching_paren, io_targets = csc.matching_paren, csc.io_targets


def declared_names(text):
    return [m.group(1).lower() for m in
            (re.match(r"\s*(\w+)", p) for p in csc.split_args(text)) if m]


def scan(files):
    """every procedure, keyed as check_star_contracts keys them ((module,
    name), an internal procedure by its host's key and its name):
    module, file, dummies, intents, locals, the modules its scope uses,
    and its executable statements (the star argument's s% written
    star%)"""
    procs, order = {}, []
    for path in files:
        lines = path.read_text(errors="replace").splitlines()
        stack, module, module_uses = [], path.stem.lower(), []
        for _, code in csc.logical_lines(lines):
            mm = MODULE_RE.match(code)
            if mm and not stack:
                module, module_uses = mm.group(1).lower(), []
                continue
            m = csc.PROC_RE.match(code)
            if m and not re.match(r"^\s*end\b", code, re.I):
                name = m.group(2).lower()
                key = (stack[-1]["key"] if stack else (module,)) + (name,)
                p = dict(name=name, key=key, module=module, path=path,
                         dummies=[a.strip().lower() for a in (m.group(3) or "").split(",") if a.strip()],
                         intents={}, local=set(), stmts=[], uses=[], module_uses=module_uses,
                         parent=stack[-1]["key"] if stack else None, star_s=False)
                if key not in procs:
                    procs[key] = p
                    order.append(p)
                stack.append(p)
                continue
            if csc.END_RE.match(code):
                if stack:
                    stack.pop()
                continue
            u = csc.parse_use(code)
            if u:
                (stack[-1]["uses"] if stack else module_uses).append(u)
                continue
            if not stack:
                continue
            p = stack[-1]
            if csc.STAR_DUMMY_RE.match(code):
                p["star_s"] = True
                continue
            if any(q["star_s"] for q in stack):
                code = csc.S_RE.sub("star%", code)
            mi = csc.INTENT_RE.match(code)
            if mi:
                for n in declared_names(mi.group(2)):
                    p["intents"][n] = mi.group(1).lower()
                continue
            md = DECL_RE.match(code)
            if md:
                p["local"] |= set(declared_names(md.group(2))) - set(p["dummies"])
                continue
            if not SKIP_RE.match(code):
                p["stmts"].append(code)
    return procs, order


resolve = csc.resolve


def visible_dummies(procs, p):
    """name -> the procedure whose dummy it is: p's own and its hosts'
    (host association), unless p declares the name itself"""
    vis, q, hidden = {}, p, set()
    while q is not None:
        for d in q["dummies"]:
            if d not in vis and d not in hidden:
                vis[d] = q
        hidden |= q["local"] | set(q["dummies"])
        q = procs.get(q["parent"]) if q["parent"] else None
    return vis


def classify(stmt, pat):
    """(written, read): the names that pat (one group) finds in a
    statement other than a call"""
    find = re.compile(r"(?<![\w%])" + pat, re.I)
    m = re.match(r"^\s*(?:if\s*\(.*\)\s*)?" + pat + r"(?:\s*\([^=]*\))?(?:%\w+(?:\s*\([^=]*\))?)*\s*=(?!=)",
                 stmt, re.I)
    if m:
        eq = stmt.index("=", m.end() - 1)
        read = {f.lower() for f in find.findall(stmt[eq + 1:])}
        # an if-statement's condition and the target's subscripts are reads
        read |= {f.lower() for f in find.findall(stmt[:m.start(1)] + " " + stmt[m.end(1):eq])}
        return {m.group(1).lower()}, read
    mr = READ_STMT_RE.match(stmt)
    if mr:
        j = matching_paren(stmt, mr.end() - 1)
        targets = {n.lower() for n in io_targets(stmt[j + 1:], re.compile(r"^" + pat, re.I))}
        return targets, {f.lower() for f in find.findall(stmt)} - targets
    if ALLOC_RE.match(stmt):
        return {f.lower() for f in find.findall(stmt)}, set()
    return set(), {f.lower() for f in find.findall(stmt)}


def analyse(procs, order):
    """each procedure's star writes and reads (s_w, s_r) and the dummies
    it assigns and reads (d_w, d_r), the call arguments refined by the
    callees' intents and use, to a fixed point"""
    for p in order:
        p.update(s_w=set(), s_r=set(), d_w=set(), d_r=set(), args=[])
    for p in order:
        vis = visible_dummies(procs, p)
        dpat = ("(" + "|".join(re.escape(d) for d in sorted(vis, key=len, reverse=True)) + r")\b") if vis else None

        def note(kind, written, read):
            if kind == "s":
                p["s_w"].update(written)
                p["s_r"].update(read)
            else:
                for n in written:
                    vis[n]["d_w"].add(n)
                for n in read:
                    vis[n]["d_r"].add(n)

        pats = [("s", STAR_PAT)] + ([("d", dpat)] if dpat else [])
        for stmt in p["stmts"]:
            mc = csc.CALL_RE.search(stmt)
            if not mc:
                for kind, pat in pats:
                    note(kind, *classify(stmt, pat))
                continue
            for kind, pat in pats:
                note(kind, set(), classify(stmt[:mc.start()], pat)[1])
            q = resolve(procs, p, mc.group(1).lower())
            for i, a in enumerate(csc.split_args(mc.group(2))):
                kw = re.match(r"^(\w+)\s*=\s*(.*)$", a)
                dummy = None
                if kw:
                    dummy, a = kw.group(1).lower(), kw.group(2)
                elif q is not None and i < len(q["dummies"]):
                    dummy = q["dummies"][i]
                for kind, pat in pats:
                    if re.match(r"^" + pat + WHOLE_REF, a, re.I):
                        m0 = re.match(r"^" + pat, a, re.I)
                        name = m0.group(1).lower()
                        note(kind, set(), classify(a[m0.end():], pat)[1])
                        p["args"].append((q, dummy, kind, name, p if kind == "s" else vis[name]))
                    else:
                        note(kind, set(), classify(a, pat)[1])
    changed = True
    while changed:
        changed = False
        for p in order:
            for q, dummy, kind, name, owner in p["args"]:
                intent = q["intents"].get(dummy) if q is not None and dummy is not None else "unknown"
                if intent in ("in", "unknown"):
                    wr, rd = False, True
                elif intent == "out":
                    wr, rd = True, False
                else:
                    wr = dummy in q["d_w"]
                    rd = dummy in q["d_r"] or (intent is None and not wr)
                tw, tr = ("s_w", "s_r") if kind == "s" else ("d_w", "d_r")
                for flag, key in ((wr, tw), (rd, tr)):
                    if flag and name not in owner[key]:
                        owner[key].add(name)
                        changed = True


def star_access():
    """writers[member] = {module: {procedures}}, readers[member] =
    {modules}, and each module's file"""
    files = sorted(f for f in SRC.rglob("*.f90")
                   if "/test/" not in f.as_posix() and "/retired/" not in f.as_posix())
    procs, order = scan(files)
    analyse(procs, order)
    writers = collections.defaultdict(lambda: collections.defaultdict(set))
    readers = collections.defaultdict(set)
    where = {}
    for p in order:
        where[p["module"]] = p["path"].relative_to(SRC).as_posix()
        for n in p["s_w"]:
            writers[n][p["module"]].add(p["name"])
        for n in p["s_r"]:
            readers[n].add(p["module"])
    return writers, readers, where


# ---------------------------------------------------------------------
# rendering

IDENT_RE = re.compile(r"(?<![\w`*%.])([A-Za-z*][\w%*]*(?:\([\w,:+*-]*\))?)(?![\w`])")
# file names and paths: a known directory first, or a file extension last
PATH_RE = re.compile(r"(?<![\w`/])((?:(?:src|core|io|state|setup|rotation|net|kap|eos|atm|numerics|math|"
                     r"convection|element_diffusion|defaults|tools|docs)/[\w*./-]*\w)|"
                     r"(?:[\w*-]+/)*[\w*-]+\.(?:f90|inc|md|py|list|defaults|tsv|mod))(?![\w`])")
# member names that read as plain words in prose
PLAIN_WORDS = {"control", "physics", "stitched"}


def escape(text):
    """|, <, > and * outside code spans"""
    out = re.sub(r"(`[^`]*`)|([<>*])", lambda m: m.group(1) or {"<": "&lt;", ">": "&gt;",
                                                                 "*": "\\*"}[m.group(2)], text)
    return out.replace("|", "\\|")


def markdown(text, member_names):
    """a table cell: paths and code names in backticks, the rest escaped"""
    def code(m):
        w = m.group(1)
        base = re.match(r"[\w%*]+", w).group(0)
        if re.search(r"[_%*]", base) or (base in member_names and base not in PLAIN_WORDS):
            return f"`{w}`"
        return w
    out = []
    for part in re.split(r"(`[^`]*`)", PATH_RE.sub(lambda m: f"`{m.group(1)}`", text)):
        out.append(part if part.startswith("`") else IDENT_RE.sub(code, part))
    return escape("".join(out))


def sentence(text):
    return text[:1].upper() + text[1:]


def written_by(w, where):
    mods = sorted(w, key=lambda m: where[m])
    n = sum(len(v) for v in w.values())
    if n > 6:
        return ", ".join(f"`{m}`" for m in mods) + f" ({n} procedures)"
    return "; ".join(f"`{m}` ({', '.join(sorted(w[m]))})" for m in mods)


def read_by(r, where):
    mods = sorted(r, key=lambda m: where[m])
    if len(mods) > 7:
        by_dir = collections.Counter(where[m].split("/")[0] for m in mods)
        dirs = ", ".join(f"{d}/ {c}" for d, c in sorted(by_dir.items(), key=lambda x: (-x[1], x[0])))
        return f"{len(mods)} modules ({dirs})"
    return ", ".join(f"`{m}`" for m in mods)


def anchor(i, name):
    return f"3.{i} {name}"


def render(classes, members, writers, readers, where):
    names = {e["name"] for e in members}
    cls_rows = ["| Class | Members | What it is | When its value is valid |", "|---|---|---|---|"]
    for i, c in enumerate(classes, 1):
        cls_rows.append(f"| [{c['name']}](#{anchor(i, c['name']).replace(' ', '-').replace('.', '').lower()}) "
                        f"| {len(c['members'])} | {markdown(c['what'], names)} | "
                        f"{markdown(sentence(c['valid']), names)} |")
    n_rec = sum(1 for e in members if e["record"])
    cls_rows.append("")
    cls_rows.append(f"{len(members)} members in all, {n_rec} of them records "
                    f"(the controls among them).")
    tables = []
    for i, c in enumerate(classes, 1):
        tables += [f"### {anchor(i, c['name'])}", "",
                   "| member | meaning | unit | written by | read by | notes |",
                   "|---|---|---|---|---|---|"]
        for e in c["members"]:
            key = e["name"].lower()
            shape = f"`{e['shape']}`" + (f" (type `{e['record']}`)" if e["record"] else "")
            tables.append(f"| {shape} | {markdown(e['meaning'], names)} | {escape(e['unit'])} | "
                          f"{written_by(writers[key], where)} | {read_by(readers[key], where)} | "
                          f"{markdown(e['notes'], names)} |")
        tables.append("")
    return "\n".join(cls_rows) + "\n", "\n".join(tables)


REMOVED_HEADING = "## 4. Members removed"


def removed_names(text):
    """the names in backticks in the first column of README section 4's
    table, with their line numbers"""
    i = text.find(REMOVED_HEADING)
    if i < 0:
        sys.exit(f"state/README.md: the section '{REMOVED_HEADING}' is missing")
    j = text.find("\n## ", i + len(REMOVED_HEADING))
    first = text[:i].count("\n") + 1
    out = []
    for k, line in enumerate(text[i:j if j > 0 else len(text)].splitlines()):
        cells = line.split("|")
        if not line.startswith("|") or len(cells) < 3 or set(cells[1].strip()) <= set("-"):
            continue
        out += [(first + k, n) for n in re.findall(r"`([A-Za-z_]\w*)`", cells[1])]
    return out


def splice(text, tag, body):
    begin = f"<!-- BEGIN GENERATED {tag}: tools/gen_star_docs.py from state/star_info_lib.f90 -->"
    end = f"<!-- END GENERATED {tag} -->"
    i, j = text.find(begin), text.find(end)
    if i < 0 or j < i:
        sys.exit(f"state/README.md: the markers of the generated {tag} are missing:\n  {begin}\n  {end}")
    return text[:i + len(begin)] + "\n" + body + text[j:]


# The self-test (--self-test, run by test_boundaries.py): the writer and
# reader analysis on check_star_contracts' synthetic module, plus a call
# that hands members to dummies the callee assigns after a case label or
# only reads: (procedure, the members it must be found to write and read).
SELF_TEST_SOURCE = csc.SELF_TEST_SOURCE.replace("end module star_tool_cases", """\
subroutine assign_after_label(k, x, y)
      integer, intent(in) :: k
      real(dp), intent(inout) :: x, y
      select case (k)
      case (1);  x = 2*y
      end select
end subroutine assign_after_label
subroutine pass_members(s)
      type(star_info), intent(inout) :: s
      call assign_after_label(1, s%passed_written, s%passed_read)
end subroutine pass_members
end module star_tool_cases""")
SELF_TEST_CASES = dict(csc.SELF_TEST_CASES, pass_members=({"passed_read"}, {"passed_written"}))


def self_test():
    """ok/FAIL per case; returns the number of failures"""
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        path = pathlib.Path(d) / "star_tool_cases.f90"
        path.write_text(SELF_TEST_SOURCE)
        procs, order = scan([path])
        analyse(procs, order)
    failures = 0
    for name, (reads, writes) in SELF_TEST_CASES.items():
        p = procs.get(("star_tool_cases", name))
        got = (p["s_r"] - p["s_w"], p["s_w"]) if p else None
        good = got == (reads, writes)
        failures += not good
        print(f"{'ok  ' if good else 'FAIL'} {name}: "
              + (f"reads {sorted(got[0])}, writes {sorted(got[1])}" if got else "not found as a procedure")
              + ("" if good else f" (want reads {sorted(reads)}, writes {sorted(writes)})"))
    return failures


def main():
    if "--self-test" in sys.argv[1:]:
        sys.exit(1 if self_test() else 0)
    classes, members, errors = parse_star_info(LIB)
    writers, readers, where = star_access()
    for e in members:
        key = e["name"].lower()
        procs = {p for ps in writers[key].values() for p in ps}
        for w in e["writers"]:
            if w.lower() not in writers[key] and w.lower() not in procs:
                errors.append(f"{e['name']}: its comment names {w} as a writer; the code's writers are "
                              f"{', '.join(sorted(writers[key])) or 'none'}")
        if not writers[key]:
            errors.append(f"{e['name']}: nothing outside the tests writes it")
        if not readers[key]:
            errors.append(f"{e['name']}: nothing outside the tests reads it")
    errors += class_rule_errors(members, writers, readers)
    text = README.read_text()
    live = {e["name"].lower() for e in members}
    for line, n in removed_names(text):
        if n.lower() in live:
            errors.append(f"state/README.md:{line}: section 4 lists {n} as removed, but it is a member")
    if errors:
        print("gen_star_docs.py: state/star_info_lib.f90 or state/README.md needs attention:", file=sys.stderr)
        print("\n".join("  " + e for e in errors), file=sys.stderr)
        sys.exit(1)
    cls_table, member_tables = render(classes, members, writers, readers, where)
    text = splice(text, "classes", cls_table)
    text = splice(text, "members", member_tables)
    emit(README, text)
    finish("gen_star_docs.py")
    print(f"state/README.md: {len(members)} members in {len(classes)} classes")


if __name__ == "__main__":
    main()

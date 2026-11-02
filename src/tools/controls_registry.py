"""The controls registry, defaults/controls_registry.tsv: one row per
control, read by the controls generators (gen_controls_state, gen_inlist,
gen_controls_dump), tools/upgrade_inlist.py and tools/check_boundaries.py;
and the one parser of io/controls_buffer.f90's declarations
(buffer_members, buffer_constants), which the generators share.

Every control has a code name, a controls_buffer member that holds its one
default as its initializer and its description as the comment directly
above its declaration (buffer_members); both kinds of deck are read into it
(tools/gen_inlist.py).

Columns:
  group          the namelist group, control or physics (the inlist's
                 &control and &physics); a derived member, which no deck
                 sets, takes the group of its topic (the controls it is
                 derived from). The star holds each control in the record
                 of its group, star%control or star%physics.
  legacy         the name in the legacy CONTROL/PHYSICS namelists
                 (io/legacy_controls.f90, frozen); '-': the legacy decks do
                 not read it (every control added since the freeze, and the
                 derived members)
  type           real, integer, logical, string or path, with (dims) for an
                 array, a size by its name in state/control_slots.f90 or
                 state/reaction_layout.f90 (real(max_runs)); path is a string that names a file the run reads or
                 writes, whose {YREC_INPUT}, {YREC_START} or {YREC_OUTPUT}
                 placeholder the read expands (io/expand_paths.inc, generated)
  proposed_name  the name in a single-file inlist (for a derived member,
                 its code name)
  status         live: read by the inlist format under proposed_name (and
                     its aliases) and, unless legacy is '-', by the legacy
                     decks under legacy;
                 legacy: read from the legacy decks only, and used;
                 unused: read from the legacy decks only, and used nowhere
                     (dropped from the inlist format);
                 derived: a member no deck sets, which the read computes
                     from the controls (or, for the deck's file names, from
                     the command line)
  mesa_note      the MESA control it corresponds to, if any
  canonical      the code name: the controls_buffer member the value lands
                 in, for every row
  home           control or physics, the group: the star's record of that
                 group (star%control, star%physics), read after the controls
                 read and never changed after it; buffer: the controls
                 read's buffer only (an input-only or unused control, used
                 read or none). tools/check_boundaries.py checks it.
  run_copy       for a control the run changes (the run cards the
                 calibrations rewrite, the stops they disarm, the latches,
                 the Monte-Carlo samples): the star member the run keeps as
                 its working copy, declared and documented in star_info_lib
                 and set from the control by the controls read
                 (controls_sync's arm_run_copies); the control keeps the
                 deck's value. Blank for the others.
                 tools/check_boundaries.py checks it.
  used           run: read after the controls read; read: only during the
                 read (map_user_inputs, derive_controls: an input-only
                 control); none: nowhere. tools/check_boundaries.py checks it.
  aliases        older inlist names the inlist format still accepts (comma
                 separated; a live control whose proposed_name replaced a
                 misleading one). An alias counts, element by element, where
                 a deck sets it and not the inlist name.
  settings       the settings fields a domain's builder in setup/ sets from
                 the control: each field a builder assigns in a statement
                 that names the control (on its right-hand side or in its
                 condition), as <record>%<field> (the record the builder
                 returns: rotation_settings%max_diffusion_iters,
                 kap_table_paths%conductive), comma separated; blank for a
                 control no builder statement names (core reads it, or a
                 builder derives from it through a local, as the wind law
                 from its inputs). tools/check_boundaries.py checks it
                 against the builders (builder_copies) and prints the cell
                 to write.
"""
import pathlib
import re

SRC = pathlib.Path(__file__).resolve().parent.parent
REG = SRC / "defaults" / "controls_registry.tsv"

COLUMNS = ["group", "legacy", "type", "proposed_name", "status", "mesa_note",
           "canonical", "home", "run_copy", "used", "aliases", "settings"]
STATUSES = {"live", "legacy", "unused", "derived"}
DEAD = {"unused", "legacy"}       # read from the legacy decks only
GROUPS = {"control", "physics"}
HOMES = GROUPS | {"buffer"}
USES = {"run", "read", "none"}


# What gen_inlist and the checks print when a row is filled wrongly: the
# cells of each kind of row.
CELLS_GUIDE = """How to fill a registry row (tools/controls_registry.py describes the columns):
  a new control, its inlist name its code name (the usual case): its declaration
      in controls_buffer, in its group's section, with its default as the
      initializer and its description as the comment above it; then the row:
      group control or physics; legacy -; type as controls_buffer declares it;
      proposed_name and canonical both the code name (the controls_buffer member);
      status live; home the group (buffer if only the read uses it); used run
      (read or none with home buffer); run_copy and aliases empty.
  an inlist name other than the code name (my_table_file for my_table_path):
      proposed_name the inlist name, canonical the code name; the rest as above.
  a control a domain reads: as above, and settings its settings field as
      record%field (rotation_settings%my_scale; check_boundaries prints the
      cell once the builder copies it).
  a control the run changes: as above, and run_copy the star member the run
      changes (state/README.md, "Add a value the run changes").
  a table path (or any file the run reads or writes): as above, type path (the read
      expands its {YREC_...} placeholder), group and home control, used run.
  a control the legacy decks read: legacy its name in io/legacy_controls.f90's
      frozen namelists; no new name goes there."""


def rows():
    """the registry's rows as dicts, checked for their columns and values"""
    header, *body = REG.read_text().splitlines()
    cols = header.split("\t")
    if cols != COLUMNS:
        raise SystemExit(f"controls_registry.tsv: columns {cols}, expected {COLUMNS}")
    out = []
    for n, line in enumerate(body, 2):
        fields = line.split("\t")
        if len(fields) != len(cols):
            raise SystemExit(f"controls_registry.tsv:{n}: {len(fields)} fields, expected {len(cols)}")
        r = dict(zip(cols, fields))
        if r["group"] not in GROUPS:
            raise SystemExit(f"controls_registry.tsv:{n}: group {r['group']!r}, expected control or physics")
        if r["home"] in GROUPS and r["home"] != r["group"]:
            raise SystemExit(f"controls_registry.tsv:{n}: home {r['home']}, but group {r['group']}: "
                             "the star holds a control in the record of its group")
        if r["status"] not in STATUSES or r["home"] not in HOMES or r["used"] not in USES:
            raise SystemExit(f"controls_registry.tsv:{n}: status/home/used "
                             f"{r['status']!r}/{r['home']!r}/{r['used']!r}")
        if not r["canonical"] or not r["home"]:
            raise SystemExit(f"controls_registry.tsv:{n}: every row has a code name (canonical, "
                             "a controls_buffer member) and its home\n" + CELLS_GUIDE)
        if (r["home"] == "buffer") != (r["used"] in ("read", "none")):
            raise SystemExit(f"controls_registry.tsv:{n}: home buffer is for an input-only "
                             "or unused control, and only for it")
        if r["run_copy"] and (r["home"] in ("buffer", "") or not re.fullmatch(r"[a-z]\w*", r["run_copy"])):
            raise SystemExit(f"controls_registry.tsv:{n}: a run copy is a star member's name, "
                             "for a control the star holds")
        if r["aliases"] and r["status"] != "live":
            raise SystemExit(f"controls_registry.tsv:{n}: only a live control has aliases")
        out.append(r)
    copies = [r["run_copy"] for r in out if r["run_copy"]]
    if len(copies) != len(set(copies)):
        raise SystemExit("controls_registry.tsv: two controls share a run copy")
    return out


def live(rs):
    """the rows the inlist format reads"""
    return [r for r in rs if r["status"] not in DEAD and r["status"] != "derived"]


# The layout modules whose named sizes the controls' arrays use (max_runs,
# num_reactions, ...): each integer constant given by a literal.
LAYOUTS = [SRC / "state" / "control_slots.f90", SRC / "state" / "reaction_layout.f90"]


def named_sizes():
    """{name: value} of the layout modules' literal integer constants"""
    sizes = {}
    for f in LAYOUTS:
        for line in _joined(f.read_text()):
            m = re.match(r"(?i)^\s*integer\s*,\s*parameter\b[^:]*::\s*(.*)$", line)
            if m:
                for ent in _entities(m.group(1)):
                    em = re.match(r"^(\w+)\s*=\s*(\d+)$", ent.strip())
                    if em:
                        sizes[em.group(1).lower()] = em.group(2)
    return sizes


def size_imports(dims_list):
    """the use statements that bring the named sizes in dims_list (max_runs,
    num_sstandard, ...) into scope: ['use control_slots, only: max_runs', ...]"""
    wanted = {w.lower() for d in dims_list for w in re.findall(r"[a-z]\w*", d, re.I)}
    out = []
    for f in LAYOUTS:
        module = re.search(r"(?im)^\s*module\s+(\w+)", f.read_text()).group(1)
        names = []
        for line in _joined(f.read_text()):
            m = re.match(r"(?i)^\s*integer\s*,\s*parameter\b[^:]*::\s*(.*)$", line)
            if m:
                names += [e.split("=")[0].strip().lower() for e in _entities(m.group(1))]
        mine = [n for n in names if n in wanted]
        if mine:
            out.append(f"      use {module}, only: {', '.join(mine)}")
    return out


def type_of(typespec, dims):
    """a declaration's type in the registry's spelling: real(50), string, ...
    (typespec as declared, attributes after a comma ignored; a named size,
    max_runs or num_reactions, written as its value)"""
    base = re.sub(r"\s+", " ", typespec.lower().split(",")[0].strip())
    if base.startswith("character"):
        base = "character"
    kind = ("real" if base in ("double precision", "real(dp)") else
            "string" if base.startswith("character") else base)
    sizes = named_sizes()
    dims = re.sub(r"\w+", lambda m: sizes.get(m.group(0).lower(), m.group(0)),
                  dims.strip("()").replace(" ", ""))
    return kind + (f"({dims})" if dims else "")


BUFFER = SRC / "io" / "controls_buffer.f90"
_DECL_RE = re.compile(
    r"^\s*(double precision|real\s*\(\s*dp\s*\)|integer|logical|character\s*\(\s*len\s*=\s*\d+\s*\))"
    r"\s*(,\s*(?:public|parameter)\b[^:]*)?::\s*(.*)$", re.IGNORECASE)


def _joined(text):
    """declaration lines with their &-continuations joined, comments dropped"""
    out, buf = [], ""
    for raw in text.splitlines():
        code = raw.split("!")[0]
        if buf:
            code = code.lstrip().lstrip("&")
        buf += code
        if buf.rstrip().endswith("&"):
            buf = buf.rstrip()[:-1]
            continue
        if buf.strip():
            out.append(buf)
        buf = ""
    return out


def _entities(text):
    """a declaration's entity list split on its top-level commas"""
    ents, depth, cur = [], 0, ""
    for ch in text:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            ents.append(cur.strip())
            cur = ""
        else:
            cur += ch
    return ents + [cur.strip()] if cur.strip() else ents


def _descriptions(text):
    """{member: the comment lines directly above its declaration, joined}
    (a topic line, '! ----', or a section line, '! ====', ends a comment);
    a declaration of several members, or a member without a comment, is
    an error"""
    docs, block, problems = {}, [], []
    lines = text.splitlines()
    in_module = False
    for i, raw in enumerate(lines):
        s = raw.strip()
        if re.match(r"(?i)^implicit none", s):
            in_module, block = True, []
            continue
        if not in_module:
            continue
        if "/*" in s:
            problems.append(f"controls_buffer.f90:{i + 1}: '/*' (make compiles with -cpp, which "
                            "reads it as the start of a C comment)")
        if s.startswith("!") and not s.startswith(("! ----", "! ====")):
            block.append(s[1:].strip())
            continue
        m = _DECL_RE.match(raw.split("!")[0])
        if m and "parameter" not in (m.group(2) or "").lower():
            stmt, k = raw.split("!")[0].rstrip(), i
            while stmt.endswith("&"):
                k += 1
                stmt = stmt[:-1] + lines[k].split("!")[0].strip().lstrip("&")
            names = [re.match(r"\s*(\w+)", e).group(1).lower()
                     for e in _entities(_DECL_RE.match(stmt).group(3))]
            if len(names) > 1:
                problems.append(f"controls_buffer.f90:{i + 1}: declare one member per line ({', '.join(names)})")
            elif not block:
                problems.append(f"controls_buffer.f90:{i + 1}: {names[0]} has no description (the "
                                "comment directly above its declaration)")
            else:
                docs[names[0]] = " ".join(block)
        block = []
    return docs, problems


# The value a member without an initializer starts at (static storage's
# zero, blank or .false.), which the records' components declare.
ZERO = {"real": "0.0d0", "integer": "0", "logical": ".false.", "string": "' '"}


def _declarations():
    """(members, constants) of controls_buffer.f90, in declaration order:
    each a dict (see buffer_members); a constant's has name, spelled and
    init only"""
    docs, problems = _descriptions(BUFFER.read_text())
    if problems:
        raise SystemExit("\n".join(problems))
    members, constants, params = [], [], {}
    for line in _joined(BUFFER.read_text()):
        m = _DECL_RE.match(line)
        if not m:
            continue
        decl = re.sub(r"\s+", "", m.group(1).lower()).replace("doubleprecision", "real(dp)")
        constant = "parameter" in (m.group(2) or "").lower()
        for ent in _entities(m.group(3)):
            em = re.match(r"^(\w+)\s*(\([^)]*\))?\s*(?:=\s*(.*))?$", ent, re.S)
            if not em:
                raise SystemExit(f"controls_buffer.f90: cannot parse {ent!r}")
            spelled, dims, init = em.group(1), (em.group(2) or "").replace(" ", ""), em.group(3)
            if constant:
                params[spelled.lower()] = init.strip()
                constants.append(dict(name=spelled.lower(), spelled=spelled, init=init.strip()))
                continue
            if init is not None:
                init = params.get(init.strip().lower(), init.strip())
            kind = type_of(decl, dims)
            members.append(dict(name=spelled.lower(), spelled=spelled, decl=decl, dims=dims, init=init,
                                default=init if init is not None else ZERO[kind.split("(")[0]],
                                type=kind, doc=docs.get(spelled.lower(), "")))
    if not members:
        raise SystemExit("controls_buffer.f90: no members parsed")
    return members, constants


def buffer_members():
    """controls_buffer's members, in declaration order: dicts with
    name (lower case), spelled (as declared: G_cgs), decl (the type spec,
    as the records declare it: real(dp), integer, logical,
    character(len=256)), dims (as written, blanks dropped: '(max_runs)',
    or ''), init (the initializer as written, a named constant replaced by
    its value; None without one), default (init, or the zero, blank or
    .false. a member without one starts at), type (the registry's
    spelling: real(50), string, ...) and doc (the comment directly above
    the declaration). Named constants are not members. Stops when a
    declaration declares several members or a member has no
    description."""
    return _declarations()[0]


def buffer_constants():
    """controls_buffer's named constants, in declaration order, as declared
    (li6_reference_xs, ...): public by name with the members"""
    return [c["spelled"] for c in _declarations()[1]]


# The F77 routines' names (defaults/file_rename_map.tsv's subroutine
# column), less the ones that are English words: a description names one
# only as "legacy X", since in lower case it reads as a name of today's code.
RENAME_MAP = SRC / "defaults" / "file_rename_map.tsv"
ENGLISH = {"wind", "solid", "reduce", "func", "mu", "remap", "stitch", "physic"}


def _code_words():
    """every name the Fortran sources use outside comments: identifiers,
    and the words of string literals (column names, file names)"""
    words = set()
    for p in sorted(SRC.rglob("*.f90")) + sorted(SRC.rglob("*.inc")):
        for line in p.read_text(errors="replace").splitlines():
            code, quote = "", None
            for ch in line:
                if quote:
                    quote = None if ch == quote else quote
                elif ch in "'\"":
                    quote = ch
                elif ch == "!":
                    break
                code += ch
            words |= {w.lower() for w in re.findall(r"[A-Za-z_]\w*", code)}
    return words


def description_problems(rs, members):
    """the descriptions that name what the code does not: a lower-case
    snake_case name that is no control's (code, inlist, alias, legacy or
    run-copy name, or a MESA name of mesa_note) and nothing the Fortran
    sources use; a prefix pattern (tol_*) that starts no such name; an F77
    routine or a legacy name in lower case, unless written "legacy X" (a
    path, with a / and an extension or a placeholder, is not checked)"""
    names, legacy = set(), set()
    for r in rs:
        names |= {r[k].lower() for k in ("proposed_name", "canonical", "run_copy") if r[k]}
        names |= {a.lower() for a in r["aliases"].split(",") if a}
        names |= {w.lower() for w in re.findall(r"\w+", r["mesa_note"])}
        if r["legacy"] != "-":
            legacy.add(r["legacy"].lower())
    legacy -= names
    f77 = set()
    for line in RENAME_MAP.read_text().splitlines():
        cells = line.split("\t")
        if len(cells) > 2 and "->" in cells[2]:
            f77.add(cells[2].split("->")[0].strip().lower())
    old = (f77 - ENGLISH) | legacy
    known = names | _code_words()
    problems = []
    for m in members:
        doc = m["doc"]
        words = []
        for token in re.findall(r"[\w.{}*/-]+", doc):
            if "/" in token and (re.search(r"\.\w+$|[{}]", token) or token.endswith("/")):
                continue
            words += re.findall(r"([A-Za-z_]\w*)(\*?)", token)
        for w, star in dict.fromkeys(words):
            lw = w.lower()
            if star:
                if not any(k.startswith(lw) for k in known):
                    problems.append(f"{m['name']}: its description's {w}* matches no name")
            elif w == lw and "_" in w and lw not in known and lw not in legacy:
                problems.append(f"{m['name']}: its description names {w}, which is no control "
                                "and nothing the code uses")
            elif w == lw and lw in old and not re.search(r"legacy\s+" + re.escape(w) + r"\b", doc):
                problems.append(f"{m['name']}: its description names {w}, an F77 or legacy name: "
                                "say what the value means, in today's names (\"(legacy "
                                f"{w.upper()})\" may follow)")
    return problems


def declared_type(cell):
    """a registry type cell as type_of spells a declaration's type: a path
    is a string, and a named size (real(max_runs)) is its value"""
    base, _, dims = re.sub(r"^path", "string", cell.replace(" ", "")).partition("(")
    return type_of(base, dims.rstrip(")"))


def check_rows_against_buffer(rs, members):
    """the problems between the registry and controls_buffer: a row whose
    code name (canonical) is not a buffer member, or has another type; a
    member without a row"""
    by_name = {m["name"]: m for m in members}
    problems = []
    for r in rs:
        m = by_name.get(r["canonical"].lower())
        if m is None:
            problems.append(f"{r['proposed_name']}: its code name (canonical) "
                            f"{r['canonical'] or '(blank)'} is not declared in controls_buffer.f90")
        elif m["type"] != declared_type(r["type"]):
            problems.append(f"{r['proposed_name']}: registry type {r['type']}, "
                            f"declared {m['type']} ({m['name']})")
    named = {r["canonical"].lower() for r in rs}
    problems += [f"{m['name']}: a controls_buffer member without a registry row"
                 for m in members if m["name"] not in named]
    return problems


SETUP = SRC / "setup"


def _builder_statements(text):
    """a builder file's statements: comments dropped, &-continuations
    joined"""
    out, buf = [], ""
    for raw in text.splitlines():
        code = raw.split("!")[0].rstrip()
        if buf:
            code = code.lstrip().lstrip("&")
        buf += code
        if buf.endswith("&"):
            buf = buf[:-1]
            continue
        if buf.strip():
            out.append(buf)
        buf = ""
    return out


def builder_copies():
    """{control (lower case): sorted [record%field, ...]}: each field of a
    record a builder in setup/ returns (result(settings), result(paths),
    declared type(<record>_t)) that the builder assigns in a statement
    naming the control, s%physics%<c> or s%control%<c>, on its right-hand
    side or in its one-line condition; a subscript is kept
    (rotation_settings%mixing_scale(i_mech_gsf))"""
    copies = {}
    for f in sorted(SETUP.glob("*.f90")):
        statements = _builder_statements(f.read_text(errors="replace"))
        record, result = None, None
        for st in statements:
            m = re.match(r"(?i)^\s*(?:(?:pure|elemental|impure|recursive)\s+)*function\s+\w+\s*\([^)]*\)"
                         r"\s*result\s*\(\s*(\w+)\s*\)", st)
            if m:
                record, result = None, m.group(1).lower()
                continue
            if re.match(r"(?i)^\s*end\s*function\b", st):
                record, result = None, None
                continue
            if result is None:
                continue
            m = re.match(r"(?i)^\s*type\s*\(\s*(\w+?)_t\s*\)\s*::\s*(\w+)\s*$", st)
            if m and m.group(2).lower() == result:
                record = m.group(1).lower()
                continue
            m = re.match(r"(?i)^\s*(?:if\s*\((?P<cond>.*)\)\s*)?(?P<r>\w+)%(?P<field>\w+(?:\s*\([^=]*?\))?)"
                         r"\s*=(?!=)(?P<rhs>.*)$", st)
            if record and m and m.group("r").lower() == result:
                name = re.sub(r"\s+", "", m.group("field")).lower()
                field = f"{record}%{name}"
                for c in re.findall(r"(?i)\bs%(?:physics|control)%(\w+)", (m.group("cond") or "") + m.group("rhs")):
                    copies.setdefault(c.lower(), set()).add(field)
    return {c: sorted(fs) for c, fs in copies.items()}


def settings_cell(fields):
    """a registry settings cell for these fields"""
    return ",".join(sorted(fields))


# the run cards' capacity, state/control_slots.f90's max_runs, read from
# there (its one source). check_boundaries' check_max_runs holds the rest
# to it: every per-run control's registry type cell names max_runs at the
# run cards' dimension, and every run-card size in the four generated
# files equals its value.
MAX_RUNS = int(named_sizes()["max_runs"])

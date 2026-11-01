"""The controls registry, defaults/controls_registry.tsv: one row per
control, read by the controls generators (gen_controls_state, gen_inlist,
gen_controls_dump), tools/upgrade_inlist.py and tools/check_boundaries.py.

Every control has a code name, a controls_buffer member that holds its one
default as its initializer; both kinds of deck are read into it
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
  type           real, integer, logical or string, with (dims) for an array
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
  doc            what the control holds (what the run does with its run
                 copy is the copy's comment in state/star_info_lib.f90)
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
"""
import pathlib
import re

SRC = pathlib.Path(__file__).resolve().parent.parent
REG = SRC / "defaults" / "controls_registry.tsv"

COLUMNS = ["group", "legacy", "type", "proposed_name", "status", "mesa_note", "doc",
           "canonical", "home", "run_copy", "used", "aliases"]
STATUSES = {"live", "legacy", "unused", "derived"}
DEAD = {"unused", "legacy"}       # read from the legacy decks only
GROUPS = {"control", "physics"}
HOMES = GROUPS | {"buffer"}
USES = {"run", "read", "none"}


# What gen_inlist and the checks print when a row is filled wrongly: the
# cells of each kind of row.
CELLS_GUIDE = """How to fill a registry row (tools/controls_registry.py describes the columns):
  a new control, its inlist name its code name (the usual case):
      group control or physics; legacy -; type as controls_buffer declares it;
      proposed_name and canonical both the code name (the controls_buffer member);
      status live; home the group (buffer if only the read uses it); used run
      (read or none with home buffer); run_copy and aliases empty.
  an inlist name other than the code name (my_table_file for my_table_path):
      proposed_name the inlist name, canonical the code name; the rest as above.
  a control the run changes: as above, and run_copy the star member the run
      changes (state/README.md, "Add a value the run changes").
  a table path: as above, with group and home control and used run.
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


# controls_buffer's named dimension
MAX_RUNS = 50


def type_of(typespec, dims):
    """a declaration's type in the registry's spelling: real(50), string, ...
    (typespec as declared, attributes after a comma ignored)"""
    base = re.sub(r"\s+", " ", typespec.lower().split(",")[0].strip())
    if base.startswith("character"):
        base = "character"
    kind = ("real" if base in ("double precision", "real(dp)") else
            "string" if base.startswith("character") else base)
    dims = dims.strip("()").replace(" ", "").replace("max_runs", str(MAX_RUNS))
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


def buffer_members():
    """controls_buffer's members, in declaration order: dicts with
    name (lower case), decl (the type spec as written: real(dp),
    character(len=256), ...), dims (as written: '(max_runs)', or ''),
    init (the initializer as written, a named constant replaced by its
    value; None without one) and type (the registry's spelling:
    real(50), string, ...). Named constants are not members."""
    members, params = [], {}
    for line in _joined(BUFFER.read_text()):
        m = _DECL_RE.match(line)
        if not m:
            continue
        decl = re.sub(r"\s+", "", m.group(1).lower()).replace("doubleprecision", "double precision")
        constant = "parameter" in (m.group(2) or "").lower()
        for ent in _entities(m.group(3)):
            em = re.match(r"^(\w+)\s*(\([^)]*\))?\s*(?:=\s*(.*))?$", ent, re.S)
            if not em:
                raise SystemExit(f"controls_buffer.f90: cannot parse {ent!r}")
            name, dims, init = em.group(1).lower(), em.group(2) or "", em.group(3)
            if constant:
                params[name] = init.strip()
                continue
            if init is not None:
                init = params.get(init.strip().lower(), init.strip())
            members.append(dict(name=name, decl=decl, dims=dims.replace(" ", ""), init=init,
                                type=type_of(decl, dims)))
    return members


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
        elif m["type"] != r["type"].replace(" ", "").replace("path", "string"):
            problems.append(f"{r['proposed_name']}: registry type {r['type']}, "
                            f"declared {m['type']} ({m['name']})")
    named = {r["canonical"].lower() for r in rs}
    problems += [f"{m['name']}: a controls_buffer member without a registry row"
                 for m in members if m["name"] not in named]
    return problems

"""The controls registry, defaults/controls_registry.tsv: one row per
control, read by the controls generators (gen_controls_state, gen_inlist),
tools/upgrade_inlist.py and tools/check_boundaries.py.

Columns:
  group          the namelist group, control or physics (the inlist's
                 &control and &physics); a derived member, which no deck
                 sets, takes the group of its topic (the controls it is
                 derived from). The star holds each control in the record
                 of its group, star%control or star%physics.
  legacy         the name in the legacy CONTROL/PHYSICS decks ('-': none)
  type           real, integer, logical or string, with (dims) for an array
  proposed_name  the name in a single-file inlist (for a derived member,
                 the member's name)
  status         live: read by both formats, the inlist under
                     proposed_name (and its aliases), the legacy decks under
                     legacy (tools/gen_inlist.py works out how: the inlist
                     name may be the legacy name, the canonical name, or a
                     third name it copies);
                 legacy: read from the legacy decks only, and used;
                 unused: read from the legacy decks only, and used nowhere
                     (dropped from the inlist format);
                 derived: a controls_buffer member that the read computes
                     (a flag derived from other controls) or that the run
                     keeps (the Newton factor fcorr)
  mesa_note      the MESA control it corresponds to, if any
  doc            what it does
  canonical      the controls_buffer member the value lands in (star%<x>
                 for a star member; blank for a legacy-only variable)
  home           control or physics, the group: the star's record of that
                 group (star%control, star%physics), read after the controls
                 read and never changed after it; buffer: the controls
                 read's buffer only (an input-only or unused control, used
                 read or none); blank without a canonical member.
                 tools/check_boundaries.py checks it.
  run_copy       for a control the run changes (the run cards the
                 calibrations rewrite, the stops they disarm, the latches,
                 the working values): the star member the run keeps as its
                 working copy, declared and documented in star_info_lib and
                 set from the control by the controls read
                 (controls_sync's arm_run_copies); the control keeps the
                 deck's value. Blank for the others.
                 tools/check_boundaries.py checks it.
  used           run: read after the controls read; read: only during the
                 read (map_user_inputs, derive_controls: an input-only
                 control); none: nowhere. tools/check_boundaries.py checks it.
  aliases        older inlist names the inlist format still accepts (comma
                 separated; a live control whose proposed_name replaced a
                 misleading one). A deck that sets the proposed name wins
                 over one that sets an alias.
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
HOMES = GROUPS | {"buffer", ""}
USES = {"run", "read", "none", ""}


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
        if bool(r["home"]) != (bool(r["canonical"]) and not r["canonical"].startswith("star%")):
            raise SystemExit(f"controls_registry.tsv:{n}: a buffer member needs a home, and only it")
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

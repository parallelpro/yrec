"""The controls registry, defaults/controls_registry.tsv: one row per
control, read by the controls generators (gen_controls_state, gen_inlist),
tools/upgrade_inlist.py and tools/check_boundaries.py.

Columns:
  group          the namelist group, control or physics; '-' for a
                 derived member, which no deck sets
  legacy         the name in the legacy CONTROL/PHYSICS decks ('-': none)
  type           real, integer, logical or string, with (dims) for an array
  proposed_name  the name in a single-file inlist (for a derived member,
                 the member's name)
  status         verified, canonical: a live control (canonical: the
                     inlist name is the code's name);
                 keep: a live control whose inlist name is its legacy name;
                 legacy: read from the legacy decks only, and used;
                 unused: read from the legacy decks only, and used nowhere
                     (dropped from the inlist format);
                 derived: a controls_buffer member that the read computes
                     (a unit number, a flag derived from other controls) or
                     that the run keeps (the Newton factor fcorr)
  mesa_note      the MESA control it corresponds to, if any
  doc            what it does
  canonical      the controls_buffer member the value lands in (star%<x>
                 for a star member; blank for a legacy-only variable)
  home           ctrl: star%ctrl, never changed after the read; job:
                 star%job, which code outside the read assigns (the run
                 cards, the calibration, the latches); blank without a
                 canonical member. tools/check_boundaries.py checks it.
  used           run: read after the controls read; read: only during the
                 read (map_user_inputs, derive_controls: an input-only
                 control); none: nowhere. tools/check_boundaries.py checks it.
"""
import pathlib
import re

SRC = pathlib.Path(__file__).resolve().parent.parent
REG = SRC / "defaults" / "controls_registry.tsv"

COLUMNS = ["group", "legacy", "type", "proposed_name", "status", "mesa_note", "doc",
           "canonical", "home", "used"]
STATUSES = {"verified", "canonical", "keep", "legacy", "unused", "derived"}
DEAD = {"unused", "legacy"}       # read from the legacy decks only
HOMES = {"ctrl", "job", ""}
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
        if r["status"] not in STATUSES or r["home"] not in HOMES or r["used"] not in USES:
            raise SystemExit(f"controls_registry.tsv:{n}: status/home/used "
                             f"{r['status']!r}/{r['home']!r}/{r['used']!r}")
        if bool(r["home"]) != (bool(r["canonical"]) and not r["canonical"].startswith("star%")):
            raise SystemExit(f"controls_registry.tsv:{n}: a buffer member needs a home, and only it")
        out.append(r)
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

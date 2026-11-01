#!/usr/bin/env python3
"""Generate the new-style inlist reader + .defaults docs from the registry.

Source of truth: defaults/controls_registry.tsv. Emits:
  io/inlist_new_decl.inc    -- read_inlist's declarations: the new-name
                               variables that do not already exist in the
                               reader's scope, and the &control / &physics
                               namelist groups (named like the legacy
                               groups, which they hide in read_inlist)
  io/legacy_mirror.inc      -- set_legacy_defaults' seeding of every legacy
                               name from its canonical twin's default
  io/adopt_canonical.inc    -- adopt_canonical_names' copies of the legacy
                               names' values into their canonical homes
  io/inlist_new_read.inc    -- read_inlist's body: seed new names from
                               the legacy defaults already assigned by
                               read_input, read the two groups from the
                               one file, then map back onto the legacy
                               names so ALL downstream code is untouched
  defaults/control.defaults, defaults/physics.defaults -- documentation
  defaults/controls_registry.md -- the registry as a table (every row's
                               names, home, use, type, default, description)

Entry classes (per live row): keep (new == legacy: the legacy variable
itself joins the new group), direct (new == the canonical internal name
read_input already copies into post-read: that variable joins the group),
newlocal (a fresh local is declared). A row's aliases (older inlist names)
are locals in the group too; one the deck sets counts when the deck does
not set the inlist name. Rows with status 'unused' are dead legacy
controls, accepted by the legacy path only and read nowhere; rows with
status 'legacy' are accepted by the legacy path only but still used
(documented separately in the .defaults files).

Defaults strategy: the generated code NEVER re-states default values.
Pre-read it seeds every new-style variable from its legacy twin (which
set_legacy_defaults has already defaulted); post-read it copies back.
Defaults thus stay defined in one place and the two paths agree by
construction. The .defaults files document the same values: each
control's default is read from set_legacy_defaults (io/legacy_controls.f90),
else from its initializer in io/controls_buffer.f90. The registry holds
no defaults; a live control without one in the code is an error.
"""
import pathlib

from generated_files import emit, finish
import controls_registry
import re
import sys

SRC = pathlib.Path(__file__).resolve().parent.parent
REG = SRC / "defaults" / "controls_registry.tsv"


def joined_decl_lines(text):
    out, buf = [], ""
    for l in text.splitlines():
        c = l.split("!")[0]
        buf = (buf + " " + c.strip().lstrip("&")) if buf else c
        if c.rstrip().endswith("&"):
            buf = buf.rstrip().rstrip("&")
            continue
        out.append(buf)
        buf = ""
    return out


def legacy_decl(name, decl_lines):
    """Return (type_spec, dims) for a legacy namelist variable."""
    for c in decl_lines:
        lc = c.lower()
        if "::" not in lc:
            continue
        typ, rest = lc.split("::", 1)
        for g in re.split(r",(?![^()]*\))", rest):
            m = re.match(r"\s*([a-z][a-z0-9_]*)\s*(\(([^)]*)\))?", g)
            if m and m.group(1) == name:
                return typ.strip(), (m.group(3) or "")
    return None, None


def reader_defaults(joined, rows=()):
    """{name: default as written in the code}: set_legacy_defaults'
    assignments, then controls_buffer's declaration initializers, a
    legacy name taking its canonical twin's (the generated mirrors)"""
    lc = (SRC / "io" / "legacy_controls.f90").read_text()
    body = lc[lc.index("subroutine set_legacy_defaults"):
              lc.index("end subroutine set_legacy_defaults")]
    buffer, params = {}, {}
    for c in joined((SRC / "io" / "controls_buffer.f90").read_text()):
        if "::" not in c:
            continue
        into = params if "parameter" in c.split("::")[0].lower() else buffer
        for g in re.split(r",(?![^()]*\))", c.split("::", 1)[1]):
            m = re.match(r"\s*([a-z]\w*)(\([^)]*\))?\s*=\s*(.+?)\s*$", g, re.I)
            if m:
                into[m.group(1).lower()] = m.group(3)
    # an initializer that names a parameter (li6_reference_xs) is its value
    buffer = {k: params.get(v.lower(), v) for k, v in buffer.items()}
    vals = {}
    for c in joined(body):
        m = re.match(r"\s*([a-z]\w*)\s*=\s*(.+?)\s*$", c, re.I)
        if m:
            value = m.group(2)
            # a mirror of a canonical control takes that control's default
            vals[m.group(1).lower()] = buffer.get(value.lower(), value)
    for r in rows:
        if r["home"] and r["canonical"] in buffer:
            vals.setdefault(r["legacy"].lower(), buffer[r["canonical"]])
    for name, value in buffer.items():
        vals.setdefault(name, value)
    return vals


def shown(value):
    """a default as an inlist would write it"""
    if value.startswith("repeat(achar(0)"):
        return "''"
    m = re.match(r"reshape\(\(/(.*?)/\)", value) or re.match(r"\(/(.*)/\)$", value)
    if m:
        value = m.group(1)
    return ", ".join(x.strip() for x in value.split(","))


def main():
    rows = [dict(group=r["group"], legacy=r["legacy"], typ=r["type"],
                 new=r["proposed_name"], status=r["status"], mesa=r["mesa_note"],
                 doc=r["doc"], canonical=r["canonical"], home=r["home"], used=r["used"],
                 run_copy=r["run_copy"],
                 aliases=[a for a in r["aliases"].split(",") if a])
            for r in controls_registry.rows()]
    derived = [r for r in rows if r["status"] == "derived"]
    rows = [r for r in rows if r["status"] != "derived"]
    # the controls reader's declarations
    reader = "".join((SRC / "io" / f).read_text() for f in (
        "controls_reader.f90", "legacy_controls.f90", "run_files.f90",
        "derived_controls.f90", "run_cards.f90"))
    decls = joined_decl_lines(reader)
    for extra in ("io/controls_buffer.f90", "state/physical_constants.f90",
                  "kap/opacity_tables.f90"):
        decls += joined_decl_lines((SRC / extra).read_text())
    code_default = reader_defaults(joined_decl_lines, rows)

    live = [r for r in rows if r["status"] not in controls_registry.DEAD]
    for r in live:
        if r["new"] == r["legacy"]:
            r["cls"] = "keep"
        elif r["canonical"] == r["new"].lower():
            r["cls"] = "direct"
        else:
            r["cls"] = "newlocal"
        if r["aliases"] and r["cls"] == "keep":
            sys.exit(f"gen_inlist: {r['legacy']}: an alias needs an inlist name other than "
                     "the legacy one")
    # each row's type is its legacy variable's (which the inlist name shares)
    wrong = []
    for r in rows:
        typ, dims = legacy_decl(r["legacy"], decls)
        if typ is not None and controls_registry.type_of(typ, dims) != r["typ"].replace(" ", ""):
            wrong.append(f"{r['legacy']}: registry {r['typ']}, declared "
                         f"{controls_registry.type_of(typ, dims)}")
    if wrong:
        sys.exit("gen_inlist: registry types disagree with the declarations:\n  " + "\n  ".join(wrong))

    # ---- declarations ----
    decl = ["! GENERATED by tools/gen_inlist.py from defaults/controls_registry.tsv",
            "! -- do not edit by hand. read_inlist's declarations: the inlist",
            "! variables that have no existing symbol in the reader's scope, and",
            "! the &control and &physics groups under the new names (they hide",
            "! legacy_controls' groups of the same names in read_inlist). See",
            "! defaults/*.defaults for docs."]
    for r in live:
        if r["cls"] != "newlocal" and not r["aliases"]:
            continue
        typ, dims = legacy_decl(r["legacy"], decls)
        if typ is None:
            print("WARNING: no legacy declaration found for", r["legacy"])
            typ, dims = "double precision", ""
        d = f"({dims})" if dims else ""
        if r["cls"] == "newlocal":
            decl.append(f"      {typ} :: {r['new']}{d}")
        for alias in r["aliases"]:
            decl.append(f"      {typ} :: {alias}{d}")
    for gname in ("control", "physics"):
        members = [n for r in live if r["group"] == gname
                   for n in [r["new"] if r["cls"] != "keep" else r["legacy"]] + r["aliases"]]
        decl.append(f"      namelist /{gname}/ &")
        line = "           & "
        outl = []
        for i, mname in enumerate(members):
            piece = mname + ("," if i < len(members) - 1 else "")
            if len(line) + len(piece) > 72:
                outl.append(line.rstrip() + " &")
                line = "           & "
            line += piece + " "
        outl.append(line.rstrip())
        decl.extend(outl)
    emit(SRC / "io" / "inlist_new_decl.inc", "\n".join(decl) + "\n")

    # ---- the legacy <-> canonical translation ----
    # every legacy variable the read copies into a controls_buffer member:
    # set_legacy_defaults seeds it from the member's default, and
    # adopt_canonical_names copies the value read back into the member
    copies = [r for r in rows if r["home"] and r["legacy"] != "-"
              and r["legacy"].lower() != r["canonical"]]
    head = ["! GENERATED by tools/gen_inlist.py from defaults/controls_registry.tsv",
            "! -- do not edit by hand."]
    emit(SRC / "io" / "legacy_mirror.inc", "\n".join(
        head + ["! set_legacy_defaults: each legacy name takes its canonical twin's",
                "! default (the controls_buffer initializer, in the buffer's pristine",
                "! state at this point)."]
        + [f"      {r['legacy']} = {r['canonical']}" for r in copies]) + "\n")
    emit(SRC / "io" / "adopt_canonical.inc", "\n".join(
        head + ["! adopt_canonical_names: each value read under a legacy name, copied",
                "! into its canonical home in controls_buffer."]
        + [f"      {r['canonical']} = {r['legacy']}" for r in copies]) + "\n")

    # ---- read block ----
    rd = ["! GENERATED by tools/gen_inlist.py -- do not edit by hand.",
          "! read_inlist's body: the inlist read. Runs instead of the legacy",
          "! CONTROL/PHYSICS read; set_legacy_defaults has already run, so",
          "! seeding the new names from their legacy twins gives the same",
          "! defaults, and the copies after the read hand the values back to",
          "! the legacy names, which the rest of the read (adopt_canonical_names",
          "! on) reads.",
          "! Seed new-style names from the legacy defaults:"]
    for r in live:
        if r["cls"] == "keep":
            continue
        rd.append(f"      {r['new']} = {r['legacy']}")
        for alias in r["aliases"]:
            rd.append(f"      {alias} = {r['legacy']}")
    rd += [
        "! Read &control and &physics, in either order, from the one file.",
        "      open(newunit=namelist_unit, file=control_nml_file, status='OLD', iostat=ios, iomsg=msg)",
        "      call check_io(ios, msg, 'open the inlist', control_nml_file, ierr)",
        "      if (ierr /= 0) return",
        "      read(unit=namelist_unit, nml=control, iostat=ios, iomsg=msg)",
        "      call check_io(ios, msg, 'read &control from', control_nml_file, ierr, namelist_unit)",
        "      if (ierr /= 0) return",
        "      rewind(namelist_unit)",
        "      read(unit=namelist_unit, nml=physics, iostat=ios, iomsg=msg)",
        "      call check_io(ios, msg, 'read &physics from', control_nml_file, ierr, namelist_unit)",
        "      if (ierr /= 0) return",
        "      close(namelist_unit)",
    ]
    aliased = [r for r in live if r["aliases"]]
    if aliased:
        rd.append("! An older inlist name (an alias) counts when the deck sets it and not")
        rd.append("! the inlist name:")
    for r in aliased:
        ne = ".neqv." if r["typ"].startswith("logical") else "/="
        for alias in r["aliases"]:
            rd.append(f"      if (any([{alias} {ne} {r['legacy']}]) .and. "
                      f".not. any([{r['new']} {ne} {r['legacy']}])) &")
            rd.append(f"           {r['new']} = {alias}")
    rd.append("! Map the new-style values back onto the legacy names:")
    for r in live:
        if r["cls"] == "keep":
            continue
        rd.append(f"      {r['legacy']} = {r['new']}")
    emit(SRC / "io" / "inlist_new_read.inc", "\n".join(rd) + "\n")

    # ---- defaults docs ----
    for gname in ("control", "physics"):
        out = [f"! {gname}.defaults -- every &{gname} variable with its default.",
               "! GENERATED by tools/gen_inlist.py from",
               "! defaults/controls_registry.tsv (edit that, then regenerate).",
               "! The defaults are the reader's: set_legacy_defaults, else the",
               "! controls_buffer initializer.",
               "", f"&{gname}", ""]
        for r in live:
            if r["group"] != gname:
                continue
            name = r["new"] if r["cls"] != "keep" else r["legacy"]
            if r["doc"]:
                out.append(f"      ! {r['doc']}")
            if r["legacy"] != name:
                out.append(f"      ! (legacy name: {r['legacy']})")
            if r["aliases"]:
                out.append(f"      ! (also accepted: {', '.join(r['aliases'])})")
            if r["legacy"].lower() not in code_default:
                sys.exit(f"gen_inlist: no default for {r['legacy']} in "
                         "set_legacy_defaults or controls_buffer")
            out.append(f"      {name} = {shown(code_default[r['legacy'].lower()])}")
            out.append("")
        out += ["/", ""]
        dead = [r for r in rows if r["status"] == "unused" and r["group"] == gname]
        if dead:
            out.append("! Unused legacy controls (accepted by the legacy reader only,")
            out.append("! referenced nowhere in the code; dropped from the new format):")
            out.append("!   " + ", ".join(r["legacy"] for r in dead))
        for r in rows:
            if r["status"] == "legacy" and r["group"] == gname:
                out.append(f"! Legacy-reader only, dropped from the new format: {r['legacy']}, the")
                out.append(f"! {r['doc']}")
        emit(SRC / "defaults" / f"{gname}.defaults", "\n".join(out) + "\n")

    # ---- the registry as a document ----
    def cell(text):
        return text.replace("|", "\\|")

    md = ["# The controls registry",
          "",
          "GENERATED by tools/gen_inlist.py from defaults/controls_registry.tsv",
          "(edit that, then regenerate). Every control: its name in a single-file",
          "inlist, in the legacy CONTROL/PHYSICS decks and in the code",
          "(controls_buffer, star%control or star%physics), its default and what it does.",
          "",
          "- **home**: `control` or `physics`, the group: the star holds the control",
          "  in the record of its namelist group, star%control or star%physics, never",
          "  changed after the controls read; `buffer` is the read's buffer only (an",
          "  input-only or unused control). `tools/check_boundaries.py` checks it.",
          "- **run copy**: for a control the run changes (the run cards the",
          "  calibrations rewrite, the stops they disarm, the latches, the working",
          "  values), the star member the run keeps as its working copy, set from",
          "  the control by the controls read (star_info_lib documents it); the",
          "  control keeps the deck's value.",
          "- **used**: `run`, read after the controls read; `read`, only during the",
          "  read (an input-only control: `map_user_inputs` or `derive_controls`",
          "  turns it into others); `none`, nowhere. `tools/check_boundaries.py`",
          "  checks it.",
          "- **status**: `live`, read by both formats (an inlist name in",
          "  parentheses is an older one the inlist format still accepts);",
          "  `legacy`, read from the legacy decks only; `unused`, read from the",
          "  legacy decks only and used nowhere; `derived`, a member no deck sets.",
          "  The columns are described in `tools/controls_registry.py`.",
          ""]
    for gname in ("control", "physics"):
        md += [f"## &{gname}", "",
               "| inlist name | legacy name | code name | home | run copy | used | type | default | status | description |",
               "|---|---|---|---|---|---|---|---|---|---|"]
        for r in rows:
            if r["group"] != gname:
                continue
            name = r["new"] if r["status"] not in controls_registry.DEAD else "(legacy only)"
            name = f"`{name}`" + "".join(f" (`{a}`)" for a in r["aliases"])
            default = shown(code_default.get(r["legacy"].lower(), "")) or " "
            copy = f"`{r['run_copy']}`" if r["run_copy"] else "-"
            md.append(f"| {name} | `{r['legacy']}` | `{r['canonical'] or '-'}` | {r['home'] or '-'} | "
                      f"{copy} | {r['used'] or '-'} | {r['typ']} | `{cell(default)}` | {r['status']} | "
                      f"{cell(r['doc'])} |")
        md.append("")
    md += ["## Members no deck sets", "",
           "The read computes them (unit numbers, flags derived from other controls)",
           "or the run keeps them. Each is in the group of its topic.", "",
           "| code name | group | home | run copy | used | type | default | description |",
           "|---|---|---|---|---|---|---|---|"]
    for r in derived:
        default = shown(code_default.get(r["canonical"], "0")) or " "
        copy = f"`{r['run_copy']}`" if r["run_copy"] else "-"
        md.append(f"| `{r['canonical']}` | {r['group']} | {r['home']} | {copy} | {r['used']} | {r['typ']} | "
                  f"`{cell(default)}` | {cell(r['doc'])} |")
    emit(SRC / "defaults" / "controls_registry.md", "\n".join(md) + "\n")

    finish("gen_inlist.py")
    print(f"live: {len(live)} (keep {sum(r['cls']=='keep' for r in live)}, "
          f"direct {sum(r['cls']=='direct' for r in live)}, "
          f"newlocal {sum(r['cls']=='newlocal' for r in live)}); "
          f"dead: {len(rows)-len(live)}")


if __name__ == "__main__":
    main()

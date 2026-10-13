#!/usr/bin/env python3
"""Generate the new-style inlist reader + .defaults docs from the registry.

Source of truth: defaults/controls_registry.tsv. Emits:
  io/inlist_new_decl.inc    -- read_inlist's declarations: the new-name
                               variables that do not already exist in the
                               reader's scope, and the &control / &physics
                               namelist groups (named like the legacy
                               groups, which they hide in read_inlist)
  io/inlist_new_read.inc    -- read_inlist's body: seed new names from
                               the legacy defaults already assigned by
                               read_input, read the two groups from the
                               one file, then map back onto the legacy
                               names so ALL downstream code is untouched
  defaults/control.defaults, defaults/physics.defaults -- documentation

Entry classes (per row): keep (new == legacy: the legacy variable itself
joins the new group), direct (new == the canonical internal name read_input
already copies into post-read: that variable joins the group), newlocal
(a fresh local is declared). Rows with status 'unused' are dead legacy
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


def reader_defaults(joined):
    """{name: default as written in the code}: set_legacy_defaults'
    assignments, then controls_buffer's declaration initializers"""
    lc = (SRC / "io" / "legacy_controls.f90").read_text()
    body = lc[lc.index("subroutine set_legacy_defaults"):
              lc.index("end subroutine set_legacy_defaults")]
    buffer = {}
    for c in joined((SRC / "io" / "controls_buffer.f90").read_text()):
        if "::" not in c or "parameter" in c.split("::")[0].lower():
            continue
        for g in re.split(r",(?![^()]*\))", c.split("::", 1)[1]):
            m = re.match(r"\s*([a-z]\w*)(\([^)]*\))?\s*=\s*(.+?)\s*$", g, re.I)
            if m:
                buffer[m.group(1).lower()] = m.group(3)
    vals = {}
    for c in joined(body):
        m = re.match(r"\s*([a-z]\w*)\s*=\s*(.+?)\s*$", c, re.I)
        if m:
            value = m.group(2)
            # a mirror of a canonical control takes that control's default
            vals[m.group(1).lower()] = buffer.get(value.lower(), value)
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
    rows = []
    header, *body = REG.read_text().splitlines()
    for l in body:
        r = dict(zip(header.split("\t"), l.split("\t")))
        rows.append(dict(group=r["group"], legacy=r["legacy"], typ=r["type"],
                         new=r["proposed_name"], status=r["status"],
                         mesa=r.get("mesa_note", ""), doc=r.get("doc", "")))
    # the controls reader: its declarations and its post-read copies
    parmin = "".join((SRC / "io" / f).read_text() for f in (
        "read_controls.f90", "legacy_controls.f90", "run_files.f90",
        "derived_controls.f90", "run_cards.f90"))
    decls = joined_decl_lines(parmin)
    for extra in ("io/controls_buffer.f90", "state/physical_constants.f90",
                  "state/integration_parameters.f90",
                  "kap/opacity_tables.f90"):
        decls += joined_decl_lines((SRC / extra).read_text())
    code_default = reader_defaults(joined_decl_lines)
    # canonical names = simple post-read copies canon = legacy
    canon = {}
    for c in joined_decl_lines(parmin):
        m = re.match(r"^\s*([a-z]\w*)\s*=\s*([a-z]\w*)\s*$", c, re.I)
        if m:
            canon.setdefault(m.group(2).lower(), m.group(1).lower())

    live = [r for r in rows if r["status"] not in ("unused", "legacy")]
    for r in live:
        if r["new"] == r["legacy"]:
            r["cls"] = "keep"
        elif canon.get(r["legacy"]) == r["new"].lower():
            r["cls"] = "direct"
        else:
            r["cls"] = "newlocal"

    # ---- declarations ----
    decl = ["! GENERATED by tools/gen_inlist.py from defaults/controls_registry.tsv",
            "! -- do not edit by hand. read_inlist's declarations: the inlist",
            "! variables that have no existing symbol in the reader's scope, and",
            "! the &control and &physics groups under the new names (they hide",
            "! legacy_controls' groups of the same names in read_inlist). See",
            "! defaults/*.defaults for docs."]
    for r in live:
        if r["cls"] != "newlocal":
            continue
        typ, dims = legacy_decl(r["legacy"], decls)
        if typ is None:
            print("WARNING: no legacy declaration found for", r["legacy"])
            typ, dims = "double precision", ""
        d = f"({dims})" if dims else ""
        decl.append(f"      {typ} :: {r['new']}{d}")
    for gname in ("control", "physics"):
        members = [r["new"] if r["cls"] != "keep" else r["legacy"]
                   for r in live if r["group"] == gname]
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

    # ---- read block ----
    rd = ["! GENERATED by tools/gen_inlist.py -- do not edit by hand.",
          "! read_inlist's body: the inlist read. Runs INSTEAD of the legacy",
          "! CONTROL/PHYSICS read; parmin's own default assignments have",
          "! already run, so seeding new names from legacy twins reproduces",
          "! the exact legacy defaults, and the post-read copies hand the",
          "! user's values back to the legacy names all downstream code",
          "! (including parmin's own post-processing) reads.",
          "! Seed new-style names from the legacy defaults:"]
    for r in live:
        if r["cls"] == "keep":
            continue
        rd.append(f"      {r['new']} = {r['legacy']}")
    rd += [
        "! Read &control and &physics, in either order, from the one file.",
        "      open(unit=standard_unit, file=control_nml_file, status='OLD', iostat=ios, iomsg=msg)",
        "      call check_io(ios, msg, 'open the inlist', control_nml_file, ierr)",
        "      if (ierr /= 0) return",
        "      read(unit=standard_unit, nml=control, iostat=ios, iomsg=msg)",
        "      call check_io(ios, msg, 'read &control from', control_nml_file, ierr, standard_unit)",
        "      if (ierr /= 0) return",
        "      rewind(standard_unit)",
        "      read(unit=standard_unit, nml=physics, iostat=ios, iomsg=msg)",
        "      call check_io(ios, msg, 'read &physics from', control_nml_file, ierr, standard_unit)",
        "      if (ierr /= 0) return",
        "      close(standard_unit)",
        "! Map the new-style values back onto the legacy names:",
    ]
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

    finish("gen_inlist.py")
    print(f"live: {len(live)} (keep {sum(r['cls']=='keep' for r in live)}, "
          f"direct {sum(r['cls']=='direct' for r in live)}, "
          f"newlocal {sum(r['cls']=='newlocal' for r in live)}); "
          f"dead: {len(rows)-len(live)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate the star's two controls records and the buffer sync module.

Source of truth: io/controls_buffer.f90's declarations (names, types,
dims, and declaration-time defaults: each control's one default) and the
registry (defaults/controls_registry.tsv): its home column (control or
physics: the star's record of the control's namelist group, star%control
or star%physics; buffer: the read's buffer only) and its run_copy column
(the star member the run keeps as its working copy of a control it
changes). Every buffer member needs a registry row of its declared type.
Emits:

  state/control_group_def.inc, state/physics_group_def.inc -- the
      component lists of the types control_group_controls and
      physics_group_controls (star_info_lib defines each around an
      include of its file): the members whose registry home is control,
      or physics. Every component gets a default initializer: the
      buffer's own where it has one, otherwise the explicit
      zero/blank/.false. that static module storage gave it -- so
      `star%control = control_group_controls()` reproduces the pristine
      pre-read_input state exactly.

  io/controls_sync.f90      -- module controls_sync (part of the controls
      read, so in io/; the only code that assigns the two records,
      tools/check_boundaries.py checks) with
      seed_controls_buffer  (the two records reset to their declared
      defaults; buffer member = star%<group>%member, and each member
      whose home is buffer -- an input-only or unused control, which the
      star does not hold -- set to its declared default),
      store_controls_to_star (the star's records from the buffer), one
      assignment per member, in declaration order, and arm_run_copies
      (each run copy, a star member star_info_lib declares, set from
      its control: star%<run_copy> = star%<group>%member).

Regenerate (and commit the results) whenever controls_buffer's member
list or the registry's home or run_copy column changes:
python3 tools/gen_controls_state.py
The build breaks loudly on drift (unknown/missing components).
"""
import pathlib
import re

from generated_files import Out, finish
import controls_registry

SRC = pathlib.Path(__file__).resolve().parent.parent
LIB = SRC / "io" / "controls_buffer.f90"
SLOTS = SRC / "state" / "control_slots.f90"
# the star's records, by registry home: the namelist group's record and its
# generated component list
RECORDS = {"control": SRC / "state" / "control_group_def.inc",
           "physics": SRC / "state" / "physics_group_def.inc"}
SYNC = SRC / "io" / "controls_sync.f90"

# Each member's home, the record of its group (star%control, star%physics)
# or the read's buffer, is the registry's home column
# (tools/controls_registry.py; check_boundaries checks it against the code);
# a control the run changes has a run copy, a star member.
TYPE_RE = re.compile(
    r"^\s*(double precision|real\s*\(\s*dp\s*\)|integer|logical|character\s*\(\s*len\s*=\s*\d+\s*\))"
    r"\s*(,\s*public)?\s*::\s*(.*)$",
    re.IGNORECASE,
)

# Named integer constants declared in controls_buffer or control_slots
# (max_runs) that the buffer's array dims use. star_info_lib includes the generated type
# bodies without `use controls_buffer`, so the emitted dims carry the
# literal value instead of the name.
PARAM_RE = re.compile(
    r"^\s*integer\s*,\s*parameter\s*(,\s*public)?\s*::\s*(\w+)\s*=\s*(\d+)\s*$",
    re.IGNORECASE,
)

ZERO = {
    "double precision": "0.0d0",
    "integer": "0",
    "logical": ".false.",
    "character": "' '",
}


def joined_decl_lines(text):
    """Join &-continued lines, strip comments."""
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


def split_entities(s):
    """Split a declaration's entity list on top-level commas."""
    ents, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            ents.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        ents.append(cur.strip())
    return ents


# A real named constant of controls_buffer (a default other members take,
# e.g. li6_reference_xs): an initializer that names one is written as its
# value, since the records and controls_sync do not import it.
REAL_PARAM_RE = re.compile(r"^\s*real\s*\(\s*dp\s*\)\s*,\s*parameter\s*(,\s*public)?\s*::\s*(.*)$",
                           re.IGNORECASE)


def parse():
    members = []  # (typespec, name, dims, init)
    params = {}   # name -> literal value
    real_params = {}
    for line in joined_decl_lines(SLOTS.read_text()):
        pm = PARAM_RE.match(line)
        if pm:
            params[pm.group(2)] = pm.group(3)
    for line in joined_decl_lines(LIB.read_text()):
        pm = PARAM_RE.match(line)
        if pm:
            params[pm.group(2)] = pm.group(3)
            continue
        rm = REAL_PARAM_RE.match(line)
        if rm:
            for ent in split_entities(rm.group(2)):
                name, value = (x.strip() for x in ent.split("=", 1))
                real_params[name.lower()] = value
            continue
        m = TYPE_RE.match(line)
        if not m:
            continue
        typespec = re.sub(r"\s+", " ", m.group(1).lower()).replace(" (", "(")
        if typespec.replace(" ", "") == "real(dp)":
            typespec = "double precision"   # star_info_lib's spelling
        for ent in split_entities(m.group(3)):
            em = re.match(r"^(\w+)\s*(\([^)]*\))?\s*(=\s*(.*))?$", ent, re.S)
            if not em:
                raise SystemExit(f"cannot parse entity: {ent!r}")
            name, dims, init = em.group(1), em.group(2) or "", em.group(4)
            for pname, pval in params.items():
                dims = re.sub(r"\b" + pname + r"\b", pval, dims)
            if init is None:
                base = "character" if typespec.startswith("character") else typespec
                init = ZERO[base]
            init = real_params.get(init.strip().lower(), init)
            members.append((typespec, name, dims, init.strip()))
    if not members:
        raise SystemExit("no members parsed from controls_buffer.f90")
    return members


def homes(members):
    """({member: 'control', 'physics' or 'buffer'}, {member: its run copy}) from the
    registry, which must give every buffer member a row of its declared
    type"""
    rows = [r for r in controls_registry.rows() if r["home"]]
    home, copy, problems = {}, {}, []
    for typespec, name, dims, _ in members:
        mine = [r for r in rows if r["canonical"] == name.lower()]
        if not mine:
            problems.append(f"{name}: no registry row (canonical column)")
            continue
        declared = controls_registry.type_of(typespec, dims)
        for r in mine:
            if r["type"].replace(" ", "") != declared:
                problems.append(f"{name}: registry type {r['type']} ({r['legacy']}), declared {declared}")
        if len({r["home"] for r in mine}) > 1 or len({r["run_copy"] for r in mine}) > 1:
            problems.append(f"{name}: two homes or run copies in the registry")
        home[name] = mine[0]["home"]
        if mine[0]["run_copy"]:
            copy[name] = mine[0]["run_copy"]
    known = {n.lower() for _, n, _, _ in members}
    for r in rows:
        if r["canonical"] not in known:
            problems.append(f"registry row {r['legacy']}: canonical {r['canonical']} is not in controls_buffer")
    if problems:
        raise SystemExit("gen_controls_state: the registry and controls_buffer disagree:\n  "
                         + "\n  ".join(problems))
    return home, copy


def write_def(path, members, home, which, what):
    with Out(path) as f:
        f.write(
            "! GENERATED by tools/gen_controls_state.py from\n"
            "! io/controls_buffer.f90 and defaults/controls_registry.tsv --\n"
            "! do not edit; regenerate instead.\n"
            f"! {what}: one component per buffer member, declaration\n"
            "! order, every one default-initialized (see the generator's header).\n"
        )
        for typespec, name, dims, init in members:
            if home[name] in which:
                f.write(f"            {typespec} :: {name}{dims} = {init}\n")


def main():
    members = parse()
    home, copy = homes(members)
    write_def(RECORDS["control"], members, home, {"control"},
              "The &control group's controls (star%control)")
    write_def(RECORDS["physics"], members, home, {"physics"},
              "The &physics group's controls (star%physics)")
    with Out(SYNC) as f:
        f.write(
            "!----------------------------------------------------------------------\n"
            "! controls_sync\n"
            "!----------------------------------------------------------------------\n"
            "! GENERATED by tools/gen_controls_state.py from\n"
            "! io/controls_buffer.f90 and defaults/controls_registry.tsv --\n"
            "! do not edit; regenerate instead.\n"
            "! The copies between controls_buffer, the buffer the decks are\n"
            "! read into, and the star. read_controls runs:\n"
            "!   call seed_controls_buffer(s)   (s%control, s%physics reset to\n"
            "!                                   their defaults; buffer <- them,\n"
            "!                                   and the defaults of the members\n"
            "!                                   the star does not hold)\n"
            "!   call read_input(s, ierr)       (the deck read over the buffer)\n"
            "!   call store_controls_to_star(s) (s%control, s%physics <- buffer)\n"
            "!   call arm_run_copies(s)         (each run copy <- its control)\n"
            "! All take the star as an argument s.\n"
            "module controls_sync\n"
        )
        names = [name for _, name, _, _ in members]
        line = "      use controls_buffer, only: "
        for i, name in enumerate(names):
            piece = name + (", " if i < len(names) - 1 else "")
            if len(line) + len(piece) > 76:
                f.write(line.rstrip() + " &\n")
                line = "           "
            line += piece
        f.write(line.rstrip() + "\n")
        f.write(
            "      use star_info_lib, only: star_info, control_group_controls, physics_group_controls\n"
            "      implicit none\n"
            "      private\n"
            "      public :: seed_controls_buffer, store_controls_to_star, arm_run_copies\n"
            "\n"
            "contains\n"
            "\n"
            "subroutine seed_controls_buffer(s)\n"
            "      type(star_info), intent(inout) :: s\n"
            "! the declared defaults, whatever an earlier read left on s\n"
            "      s%control = control_group_controls()\n"
            "      s%physics = physics_group_controls()\n"
        )
        for _, name, _, init in members:
            if home[name] == "buffer":
                f.write(f"      {name} = {init}\n")
            else:
                f.write(f"      {name} = s%{home[name]}%{name}\n")
        f.write(
            "end subroutine seed_controls_buffer\n"
            "\n"
            "subroutine store_controls_to_star(s)\n"
            "      type(star_info), intent(inout) :: s\n"
        )
        for _, name, _, _ in members:
            if home[name] != "buffer":
                f.write(f"      s%{home[name]}%{name} = {name}\n")
        f.write(
            "end subroutine store_controls_to_star\n"
            "\n"
            "! The run's working copies of the controls it changes (the registry's\n"
            "! run_copy column), from the controls just stored: star_info_lib\n"
            "! documents each copy and who changes it.\n"
            "subroutine arm_run_copies(s)\n"
            "      type(star_info), intent(inout) :: s\n"
        )
        for _, name, _, _ in members:
            if name in copy:
                f.write(f"      s%{copy[name]} = s%{home[name]}%{name}\n")
        f.write(
            "end subroutine arm_run_copies\n"
            "\n"
            "end module controls_sync\n"
        )
    finish("gen_controls_state.py")
    print(f"{len(members)} members -> {', '.join(p.name for p in RECORDS.values())}, "
          f"{SYNC.name}; {len(copy)} run copies")


if __name__ == "__main__":
    main()

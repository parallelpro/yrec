#!/usr/bin/env python3
"""Convert a legacy YREC deck to a single-file inlist, or rename the
groups of an inlist written before 2026-10.

Usage:
    python3 tools/upgrade_inlist.py CASE.nml1 CASE.nml2 [-o inlist_CASE]
    python3 tools/upgrade_inlist.py DECK.inlist [DECK2.inlist ...]

A legacy pair (the CONTROL namelist in .nml1, PHYSICS in .nml2) becomes
one file holding a &control and a &physics group, with every variable
renamed per defaults/controls_registry.tsv. Values, array subscripts,
comments and line layout are preserved verbatim -- only the names
change. Unused legacy controls (status 'unused' in the registry) and
legacy-reader-only ones (status 'legacy') are dropped with a warning.
The converted file must produce byte-identical outputs to the legacy
pair (test_inlist_convert.py enforces this on the examples).

An inlist whose groups are named &star_job and &controls (their names
from 2026-08 to 2026-10) has them renamed &control and &physics in
place (or written to -o, for a single file); nothing else changes.
"""
import argparse
import pathlib
import re
import sys

SRC = pathlib.Path(__file__).resolve().parent.parent
REG = SRC / "defaults" / "controls_registry.tsv"


def load_registry():
    ren = {"control": {}, "physics": {}}
    dead = {"control": set(), "physics": set()}
    header, *body = REG.read_text().splitlines()
    for l in body:
        r = dict(zip(header.split("\t"), l.split("\t")))
        grp, legacy, new, status = r["group"], r["legacy"], r["proposed_name"], r["status"]
        if status in ("unused", "legacy"):
            dead[grp].add(legacy)
        else:
            ren[grp][legacy] = new
    return ren, dead


def convert_group(text, group_in, group_out, renames, dead, other_renames):
    """Rename variables inside one namelist group's record text.

    A line whose assignments all resolve in the OTHER group's rename
    map belongs to a control whose new-style home moved across groups
    (e.g. legacy CONTROL FSHORT -> &physics log_output_file); such
    lines are renamed and returned separately for the caller to splice
    into the other output group."""
    # capture the record: &group ... / (terminator on its own or line end)
    m = re.search(r"[$&]" + group_in + r"\b(.*?)(^\s*/\s*$|[$&]end)",
                  text, re.S | re.I | re.M)
    if not m:
        sys.exit(f"no &{group_in} group found")
    body = m.group(1)
    out_lines = []
    dropped = []
    moved = []
    # a variable assignment starts a line or follows a comma:
    # rename NAME when followed by '=' or '(' (array element)
    pat = re.compile(r"(^|[,\s])([A-Za-z]\w*)(\s*(?:\([^)]*\))?\s*=)")

    def rename(mm):
        name = mm.group(2)
        low = name.lower()
        if low in dead:
            dropped.append(name)
            return mm.group(0)  # marked after; line dropped below
        return mm.group(1) + renames.get(low, name) + mm.group(3)

    for line in body.splitlines():
        code, sep, comment = line.partition("!")
        # drop whole lines that only set dead OR retired controls. A name
        # absent from the registry altogether is a RETIRED control (the
        # 2026 retire-legacy campaign deletes registry rows outright);
        # passing it through would make the run's namelist read fail, so
        # it is dropped with a marker instead.
        assigns = re.findall(r"([A-Za-z]\w*)\s*(?:\([^)]*\))?\s*=", code)
        gone = lambda a: a.lower() in dead or a.lower() not in renames
        migrated = lambda a: (a.lower() in other_renames
                              and a.lower() not in renames
                              and a.lower() not in dead)
        if assigns and all(migrated(a) for a in assigns):
            def rename_other(mm):
                return (mm.group(1) + other_renames[mm.group(2).lower()]
                        + mm.group(3))
            moved.append(pat.sub(rename_other, code)
                         + (sep + comment if sep else ""))
            continue
        if assigns and all(gone(a) for a in assigns):
            dropped.extend(assigns)
            tag = ("dropped unused legacy control"
                   if all(a.lower() in dead for a in assigns)
                   else "dropped retired control")
            out_lines.append(f"! ({tag}) " + line.strip())
            continue
        new_code = pat.sub(rename, code)
        out_lines.append(new_code + (sep + comment if sep else ""))
    return "&" + group_out + "\n".join(out_lines) + "\n/\n", dropped, moved


# the group headers of an inlist from 2026-08 to 2026-10
OLD_GROUPS = [(re.compile(r"^(\s*)&star_job\b", re.I | re.M), r"\1&control"),
              (re.compile(r"^(\s*)&controls\b", re.I | re.M), r"\1&physics")]


def rename_groups(path, output):
    text = pathlib.Path(path).read_text()
    for (pattern, new), old in zip(OLD_GROUPS, ("&star_job", "&controls")):
        text, n = pattern.subn(new, text)
        if n != 1:
            sys.exit(f"{path}: expected one {old} group, found {n}")
    pathlib.Path(output or path).write_text(text)
    print(f"wrote {output or path}")


def convert_pair(nml1, nml2, output):
    ren, dead = load_registry()
    t1 = pathlib.Path(nml1).read_text(errors="replace")
    t2 = pathlib.Path(nml2).read_text(errors="replace")
    g1, d1, m1 = convert_group(t1, "control", "control", ren["control"],
                               dead["control"], ren["physics"])
    g2, d2, m2 = convert_group(t2, "physics", "physics", ren["physics"],
                               dead["physics"], ren["control"])
    # splice cross-group migrations into their new home (before the '/')
    if m1:
        g2 = g2[:g2.rfind("\n/\n")] + "\n" + "\n".join(m1) + "\n/\n"
    if m2:
        g1 = g1[:g1.rfind("\n/\n")] + "\n" + "\n".join(m2) + "\n/\n"
    out = pathlib.Path(output or ("inlist_" + pathlib.Path(nml1).stem))
    # (use_legacy_output is retired: every run produces the unified
    # MESA-style output set, so nothing is stamped any more.)
    out.write_text("! Converted from " + nml1 + " + " + nml2 +
                   " by tools/upgrade_inlist.py\n" + g1 + "\n" + g2)
    for d in d1 + d2:
        print(f"WARNING: dropped unused/retired legacy control {d}", file=sys.stderr)
    print(f"wrote {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+",
                    help="a legacy CASE.nml1 CASE.nml2 pair, or inlists "
                         "with &star_job/&controls groups")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()
    old_style = [any(p.search(pathlib.Path(f).read_text(errors="replace"))
                     for p, _ in OLD_GROUPS) for f in args.files]
    if all(old_style):
        if args.output and len(args.files) > 1:
            sys.exit("-o takes a single inlist")
        for f in args.files:
            rename_groups(f, args.output)
    elif len(args.files) == 2 and not any(old_style):
        convert_pair(*args.files, args.output)
    else:
        sys.exit("give a legacy .nml1 .nml2 pair, or inlists whose groups "
                 "are &star_job and &controls")


if __name__ == "__main__":
    main()

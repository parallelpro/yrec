#!/usr/bin/env python3
"""Convert a legacy YREC inlist pair (.nml1 + .nml2) to new-style.

Usage:
    python3 tools/upgrade_inlist.py CASE.nml1 CASE.nml2 [-o inlist_CASE]

Produces one new-style file containing &star_job (from &control) and
&controls (from &physics), with every variable renamed per
defaults/controls_registry.tsv. Values, array subscripts, comments and
line layout are preserved verbatim -- only the names change. Unused
legacy controls (status 'unused' in the registry) are dropped with a
warning. The converted file must produce byte-identical outputs to the
legacy pair (test_inlist_convert.py enforces this on the examples).
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
    for l in REG.read_text().splitlines()[1:]:
        parts = l.split("\t")
        grp, legacy, new, status = parts[0], parts[1], parts[4], parts[5]
        if status == "unused":
            dead[grp].add(legacy)
        else:
            ren[grp][legacy] = new
    return ren, dead


def convert_group(text, group_in, group_out, renames, dead):
    """Rename variables inside one namelist group's record text."""
    # capture the record: &group ... / (terminator on its own or line end)
    m = re.search(r"[$&]" + group_in + r"\b(.*?)(^\s*/\s*$|[$&]end)",
                  text, re.S | re.I | re.M)
    if not m:
        sys.exit(f"no &{group_in} group found")
    body = m.group(1)
    out_lines = []
    dropped = []
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
        if assigns and all(gone(a) for a in assigns):
            dropped.extend(assigns)
            tag = ("dropped unused legacy control"
                   if all(a.lower() in dead for a in assigns)
                   else "dropped retired control")
            out_lines.append(f"! ({tag}) " + line.strip())
            continue
        new_code = pat.sub(rename, code)
        out_lines.append(new_code + (sep + comment if sep else ""))
    return "&" + group_out + "\n".join(out_lines) + "\n/\n", dropped


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("nml1")
    ap.add_argument("nml2")
    ap.add_argument("-o", "--output")
    args = ap.parse_args()
    ren, dead = load_registry()
    t1 = pathlib.Path(args.nml1).read_text(errors="replace")
    t2 = pathlib.Path(args.nml2).read_text(errors="replace")
    g1, d1 = convert_group(t1, "control", "star_job", ren["control"],
                           dead["control"])
    g2, d2 = convert_group(t2, "physics", "controls", ren["physics"],
                           dead["physics"])
    out = pathlib.Path(args.output or
                       ("inlist_" + pathlib.Path(args.nml1).stem))
    # (use_legacy_output is retired: every run produces the unified
    # MESA-style output set, so nothing is stamped any more.)
    out.write_text("! Converted from " + args.nml1 + " + " + args.nml2 +
                   " by tools/upgrade_inlist.py\n" + g1 + "\n" + g2)
    for d in d1 + d2:
        print(f"WARNING: dropped unused/retired legacy control {d}", file=sys.stderr)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()

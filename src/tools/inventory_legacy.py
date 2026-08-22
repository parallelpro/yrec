#!/usr/bin/env python3
"""Inventory of legacy constructs: blanket SAVE and goto (2026).

Position (see ROADMAP.md "Legacy constructs: position and policy"):
- Blanket `save` as hidden static state is the enemy of re-entrancy;
  every occurrence is either (a) already covered by the phase-five
  reset machinery, (b) a benign first-call table cache, or (c) unaudited
  debt. This tool counts them so (c) shrinks measurably.
- `goto` is tolerated in hot F77 kernels (replacing it risks byte
  drift for zero functional gain) and eliminated opportunistically
  everywhere else; the ierr funnels and evolve_step's status contract
  already removed the structurally worst ones.

Run from src/: python3 tools/inventory_legacy.py
"""
import pathlib
import re
import sys

root = pathlib.Path(__file__).resolve().parent.parent


def code(line):
    return line.split("!")[0]


def main():
    rows = []
    for p in sorted(root.rglob("*.f90")):
        rel = p.relative_to(root)
        if "test" in rel.parts:
            continue
        text = p.read_text(errors="replace").splitlines()
        blanket_save = sum(1 for l in text if re.match(r"\s*save\s*$", code(l)))
        named_save = sum(1 for l in text if re.match(r"\s*save\s*::|\s*save\s+\w", code(l)))
        gotos = sum(len(re.findall(r"\bgo\s*to\s+\d", code(l), re.I)) for l in text)
        if blanket_save or named_save or gotos:
            rows.append((str(rel), blanket_save, named_save, gotos))
    dom = {}
    for rel, bs, ns, g in rows:
        d = rel.split("/")[0]
        a = dom.setdefault(d, [0, 0, 0])
        a[0] += bs
        a[1] += ns
        a[2] += g
    print(f"{'domain':12s} {'blanket-save':>12s} {'named-save':>10s} {'goto':>6s}")
    for d in sorted(dom):
        bs, ns, g = dom[d]
        print(f"{d:12s} {bs:12d} {ns:10d} {g:6d}")
    tb = sum(v[0] for v in dom.values())
    tn = sum(v[1] for v in dom.values())
    tg = sum(v[2] for v in dom.values())
    print(f"{'TOTAL':12s} {tb:12d} {tn:10d} {tg:6d}")
    if "-v" in sys.argv:
        print("\nper file (blanket-save desc):")
        for rel, bs, ns, g in sorted(rows, key=lambda r: -r[1]):
            print(f"  {rel:44s} save={bs:2d} named={ns:2d} goto={g:3d}")


if __name__ == "__main__":
    main()

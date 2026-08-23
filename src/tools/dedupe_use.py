#!/usr/bin/env python3
"""Merge duplicate USE statements within each specification part.

Within one contiguous use-block (comments/blanks allowed between):
  - several `use M, only: ...` -> one, only-list unioned in first-
    appearance order (case-insensitive dedupe, first spelling kept);
  - a blanket `use M` absorbs plain only-lists of the same module,
    but `=>` renames are kept as their own statement (a blanket
    import does not provide the local rename).
Statements keep the position of their first occurrence.
"""
import pathlib
import re
import sys

USE_RE = re.compile(r"^(\s*)use\s+(\w+)\s*(,\s*only\s*:\s*(.*))?$", re.I)


def parse_stmt(text):
    body = ""
    for raw in text:
        c = raw.split("!")[0].rstrip()
        if c.endswith("&"):
            c = c[:-1]
        body += " " + c.lstrip("& ")
    m = USE_RE.match(body.strip().replace("  ", " ") if False else body.lstrip())
    m = re.match(r"^use\s+(\w+)\s*(?:,\s*only\s*:\s*(.*))?$", body.strip(), re.I)
    if not m:
        return None
    mod = m.group(1)
    only = m.group(2)
    items = None
    if only is not None:
        items = [x.strip() for x in only.split(",") if x.strip()]
    return mod, items


def emit(indent, mod, items):
    if items is None:
        return [f"{indent}use {mod}\n"]
    line = f"{indent}use {mod}, only: "
    out = []
    cur = line
    for i, it in enumerate(items):
        piece = it + (", " if i < len(items) - 1 else "")
        if len(cur) + len(piece) > 78 and cur.strip() != f"use {mod}, only:":
            out.append(cur.rstrip() + " &\n")
            cur = indent + "     " + piece
        else:
            cur += piece
    out.append(cur.rstrip() + "\n")
    return out


def transform(path):
    lines = path.read_text().splitlines(keepends=True)
    n = len(lines)
    out = []
    i = 0
    changed = 0
    while i < n:
        if not re.match(r"^\s*use\s+\w+", lines[i], re.I):
            out.append(lines[i])
            i += 1
            continue
        # collect the block
        entries = []   # (indent, mod, items, passthrough_lines)
        inter = []     # (position_index, comment/blank line)
        j = i
        while j < n:
            l = lines[j]
            if re.match(r"^\s*use\s+\w+", l, re.I):
                stmt = [l]
                while stmt[-1].split("!")[0].rstrip().endswith("&"):
                    j += 1
                    stmt.append(lines[j])
                parsed = parse_stmt(stmt)
                indent = re.match(r"^(\s*)", stmt[0]).group(1)
                mod, items = parsed
                entries.append([indent, mod, items, stmt])
                j += 1
            elif l.strip() == "" or l.lstrip().startswith("!"):
                inter.append((len(entries), l))
                j += 1
            else:
                break
        # merge per module
        merged = []  # list of dicts in first-appearance order
        bykey = {}
        renames_kept = []
        for indent, mod, items, stmt in entries:
            key = mod.lower()
            if key not in bykey:
                rec = {"indent": indent, "mod": mod, "items": items,
                       "blanket": items is None}
                bykey[key] = rec
                merged.append(rec)
                continue
            changed += 1
            rec = bykey[key]
            if items is None:
                rec["blanket"] = True
                rec["items"] = None
                continue
            if rec["blanket"]:
                ren = [x for x in items if "=>" in x]
                if ren:
                    renames_kept.append((rec, ren))
                continue
            seen = {re.sub(r"\s+", "", x).lower() for x in rec["items"]}
            for x in items:
                if re.sub(r"\s+", "", x).lower() not in seen:
                    rec["items"].append(x)
                    seen.add(re.sub(r"\s+", "", x).lower())
        # re-emit: merged uses first (original module order), then the
        # interleaved comments/blanks in their original relative spots
        # (simplest faithful option: comments that sat between use
        # lines are appended after the block, preserving content).
        if changed == 0 or all(
                len([e for e in entries if e[1].lower() == k]) == 1
                for k in bykey):
            # nothing merged in this block; emit verbatim
            for _, _, _, stmt in entries:
                out.extend(stmt)
            for _, l in inter:
                out.append(l)
            i = j
            continue
        for rec in merged:
            out.extend(emit(rec["indent"], rec["mod"],
                            None if rec["blanket"] else rec["items"]))
            for r2, ren in renames_kept:
                if r2 is rec:
                    out.extend(emit(rec["indent"], rec["mod"], ren))
        for _, l in inter:
            out.append(l)
        i = j
    if changed:
        path.write_text("".join(out))
    return changed


def main():
    total = files = 0
    for arg in sys.argv[1:]:
        p = pathlib.Path(arg)
        for f in ([p] if p.is_file() else sorted(p.rglob("*.f90"))):
            k = transform(f)
            if k:
                files += 1
                total += k
    print(f"merged {total} duplicate use statements in {files} files")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Domain-boundary checker for YREC (phase three, ROADMAP.md stage 2).

Enforces the public/private rule established in stage 1: a file outside
a physics domain may only call that domain's *public* entry points.
The public surface below was derived empirically on 2026-08-21, right
after stage 1 closed the last bypass, and is maintained by hand -- the
YREC analogue of MESA's public/ directories, expressed as data until
per-module library builds enforce it structurally.

Adding a name here is an API decision: it should accompany a deliberate
new public entry (ideally on a <domain>_lib module), not silence a
failure. Removing a name (surface shrinkage) needs no action.

Comment text is stripped before scanning, so historical `! CALL FOO`
lines don't trip the check.

Exit status: 0 if the observed cross-domain surface is a subset of the
allowlist; 1 otherwise, listing each violation.
"""
import re
import sys
import pathlib
import collections

SRC = pathlib.Path(__file__).resolve().parent.parent

DOMAINS = ["eos", "kap", "atm", "nuclear", "wind", "mixing",
           "rotation", "numerics"]

# domain -> names callable from outside that domain.
PUBLIC = {
    # The three eos_lib facade entries. Everything else in eos/
    # (eqstat/eqstat2, the OPAL/MHD/SCV/Yale internals, mu) is private.
    "eos": {"eos_get", "eos_get_gamma1", "eos_init"},
    # The three kap_lib facade entries.
    "kap": {"kap_get", "kap_init", "kap_update_surface_tables"},
    # atm_lib's three entries, plus surfbc (the solver's boundary-
    # condition wrapper, sole caller core/crrect.f90) and the turnover/
    # diagnostics consumed by core/io/rotation (calcad, gettau).
    "atm": {"atm_get", "atm_get_surface_pt", "atm_init",
            "surfbc", "calcad", "gettau"},
    # Every nuclear_lib module procedure is public by construction.
    "nuclear": {"neutrino", "nulosses", "azbar", "sneut", "rates",
                "eqburn", "dburn", "dburnm", "deutrate", "engeb",
                "liburn", "liburn2", "lirate88", "safedivexp",
                "ifermi12", "zfermim12"},
    # mixing has no facade module (phase-two finding: mix.f90 is
    # already the orchestrator); these are its de-facto public surface.
    "mixing": {"mix", "mixcz", "convec", "bursmix", "rotmix", "hsubp",
               "oversh", "sconvec"},
    # rotation deliberately has no facade (multi-primitive surface,
    # user decision during the phase-two sweep). NOTE: "func" is only
    # here because numerics_lib's qgauss hard-codes a call to it (the
    # F77 fixed-name-integrand idiom) -- a pre-existing backwards
    # numerics->rotation dependency, recorded in ROADMAP.md.
    "rotation": {"getw", "getrot", "fpft", "momi", "ovrot", "viscos",
                 "wczimp", "grsett", "microdiff", "ndifcom", "mixcom",
                 "mixgrid", "func"},
    "wind": {"massloss", "mwind", "mcowind"},
    # numerics is the shared numerics library: all-public by design.
    "numerics": None,
}

CALL_RE = re.compile(r"\bcall\s+([a-z0-9_]+)\s*\(", re.I)
DEF_RE = re.compile(
    r"^\s*(?:recursive\s+)?(?:double precision\s+function|"
    r"integer\s+function|logical\s+function|function|subroutine)"
    r"\s+([a-z0-9_]+)",
    re.M | re.I,
)


def strip_comments(text):
    return "\n".join(line.split("!")[0] for line in text.split("\n"))


def main():
    defs = {}
    for dom in DOMAINS:
        for path in (SRC / dom).rglob("*.f90"):
            for m in DEF_RE.finditer(strip_comments(path.read_text())):
                defs[m.group(1).lower()] = dom

    violations = collections.defaultdict(set)
    for path in SRC.rglob("*.f90"):
        rel = path.relative_to(SRC)
        caller_dom = rel.parts[0] if rel.parts[0] in DOMAINS else "app"
        code = strip_comments(path.read_text())
        for m in CALL_RE.finditer(code):
            name = m.group(1).lower()
            dom = defs.get(name)
            if dom is None or dom == caller_dom:
                continue
            allowed = PUBLIC[dom]
            if allowed is not None and name not in allowed:
                violations[(dom, name)].add(str(rel))

    if violations:
        print("Domain-boundary violations (cross-domain calls to "
              "non-public entries):")
        for (dom, name) in sorted(violations):
            callers = ", ".join(sorted(violations[(dom, name)]))
            print(f"  {dom}/{name}  called from: {callers}")
        print("\nRoute these through the domain's <domain>_lib facade, "
              "or (deliberately) add the name to PUBLIC in "
              "tools/check_boundaries.py.")
        return 1

    print("Domain boundaries OK: every cross-domain call goes through "
          "a public entry.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

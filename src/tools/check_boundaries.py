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

DOMAINS = ["eos", "kap", "atm", "net", "wind", "mixing",
           "rotation", "numerics"]

# domain -> names callable from outside that domain.
PUBLIC = {
    # The three eos_lib facade entries. Everything else in eos/
    # (eqstat/eqstat2, the OPAL/MHD/SCV/Yale internals, mu) is private.
    "eos": {"eos_get", "eos_get_r", "eos_get_gamma1", "eos_init", "eos_set_mixture"},
    # The kap_lib facade entries (kap_get_r is the named-index
    # result-array variant of kap_get).
    "kap": {"kap_get", "kap_get_r", "kap_init", "kap_update_surface_tables"},
    # atm_lib's three entries, plus surfbc (the solver's boundary-
    # condition wrapper, sole caller core/crrect.f90) and the turnover/
    # diagnostics consumed by core/io/rotation (calcad, gettau).
    "atm": {"atm_get_surface_pt", "atm_init",
            # pure surface-pressure lookups, called by the star
            # layer's envelope integrator (core/envint_lib) since the
            # atm split -- clean physics services
            "surfp", "kcsurfp", "alsurfp"},
    # Every net_lib module procedure is public by construction.
    "net": {"neutrino", "nulosses", "azbar", "sneut", "rates",
                "eqburn", "dburn", "dburnm", "deutrate", "engeb",
                "liburn", "liburn2", "lirate88", "safedivexp",
                "ifermi12", "zfermim12"},
    # mixing has no facade module (phase-two finding: mix.f90 is
    # already the orchestrator); these are its de-facto public surface.
    # tpgrad joined when phase four's step 1 moved it here from misc/
    # (it is MLT convection physics): it is legitimately called from
    # core (physic/coefft/starin) and atm (qenv) pending the
    # star-layer consolidation.
    "mixing": {"mix", "homogenize_convection_zones", "find_convection_zones",
               "burn_settle_mix", "rotmix", "compute_scale_height",
               "overshoot_boundaries", "semiconvection",
               "temperature_gradients", "temperature_gradients_r"},
    # rotation deliberately has no facade (multi-primitive surface,
    # user decision during the phase-two sweep). "func" was here
    # because numerics' qgauss hard-coded a call to it; phase four's
    # step 2 made it a procedure argument, but func stays public since
    # fpft passes it across the module boundary into numerics. solid
    # joined when step 1 moved it here from misc/ (rotation geometry,
    # legitimately called from setup/midmod and wind/wcz).
    "rotation": {"evolve_angular_momentum", "omega_from_j",
                 "rotation_shape_factors", "zone_moments_of_inertia",
                 "am_convective_regions", "viscos",
                 "enforce_rotation_profile", "gravitational_settling",
                 "microdiff", "diffuse_composition_driver",
                 "diffuse_composition", "composition_grid",
                 "equipotential_integrand", "solid_body_omega"},
    "wind": {"massloss", "matt_wind", "wind_spindown_matt"},
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

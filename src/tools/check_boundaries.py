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

# The controls-read BUFFER (io/controls_lib.f90) is reader-internal:
# star%ctrl / star%job are the authoritative post-read home, and every
# consumer reads them. Only the read pipeline itself may import the
# buffer's bare names -- read_controls (the reader), controls_sync
# (the generated buffer<->star copies), map_user_inputs (runs inside
# the read, before the store), and net's test harness (which seeds
# the buffer to simulate that pipeline). A new name here is a design
# regression, not an allowlist entry.
CONTROLS_BUFFER_IMPORTERS = {
    "io/read_controls.f90",
    "state/controls_sync_lib.f90",
    "setup/map_user_inputs.f90",
    "net/test/test_net.f90",
}

# Named constants in controls_lib that are not buffer state: any file
# may `use controls_lib, only: <these>` (an only-list drawn entirely
# from these is not a buffer import). max_runs dimensions the per-run
# arrays; ichi_*/itime_* name the chi_grid_scale/atime slots.
CONTROLS_LIB_CONSTANT_RE = re.compile(r"^(max_runs|ichi_\w+|itime_\w+)$")

# domain -> names callable from outside that domain.
PUBLIC = {
    # The three eos_lib facade entries. Everything else in eos/
    # (eqstat/eqstat2, the OPAL/MHD/SCV/Yale internals, mu) is private.
    # eos_get is the named-index result-array query (the former
    # eos_get_r); the long-form engine eos_eval is domain-internal.
    # eos_set_debye_huckel_z: kap/setupopac hands the eos domain the
    # LAOL89 table's 18-element metal mixture (2026 wave 2).
    # eos_get_gamma1 is test-only: its sole caller is
    # eos/test/test_eos.f90 (the calcad diagnostic it served is gone).
    "eos": {"eos_get", "eos_get_gamma1", "eos_init", "eos_set_mixture",
            "eos_set_debye_huckel_z"},
    # The kap_lib facade entries (kap_get_r is the named-index
    # result-array variant of kap_get).
    # kap_get is the named-index result-array query (the former
    # kap_get_r); the long-form engine kap_eval is domain-internal.
    "kap": {"kap_get", "kap_init", "kap_update_surface_tables"},
    # atm_lib's three entries, plus surfbc (the solver's boundary-
    # condition wrapper, sole caller core/crrect.f90) and the turnover/
    # diagnostics consumed by core/io/rotation (calcad, gettau).
    "atm": {"atm_get_surface_pt", "atm_init",
            "ttau_log10_temperature", "ttau_start_log10_temperature",
            "ttau_photosphere_x_limit", "hsra_t_tau_offset",
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
               "temperature_gradients"},
    # rotation deliberately has no facade (multi-primitive surface,
    # user decision during the phase-two sweep). "func" was here
    # because numerics' qgauss hard-coded a call to it; phase four's
    # step 2 made it a procedure argument, but func stays public since
    # fpft passes it across the module boundary into numerics.
    # rotation_shape_factors (the former solid) joined when step 1
    # moved it here from misc/ (rotation geometry, legitimately called
    # from core and setup/rezone).
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



MATH_FUNCS = r"(?:exp|log|log10|sin|cos|tan|asin|acos|atan|atan2|sinh|cosh|tanh)"

def mask_strings(line):
    out, q = [], None
    for ch in line:
        if q:
            out.append(" ")
            if ch == q:
                q = None
        elif ch in "'\"":
            q = ch
            out.append(" ")
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out)


STAR_BLIND_FILES = [
    # envint purity split (2026): the integration kernel takes its
    # configuration explicitly and must never touch the star layer
    "core/envint_kernel.f90",
    "atm/ttau_lib.f90",
]

def check_star_blind(src_root):
    problems = []
    for rel in STAR_BLIND_FILES:
        text = (src_root / rel).read_text(errors="replace")
        code = "\n".join(l.split("!")[0] for l in text.splitlines())
        if "star_info_lib" in code or "star%" in code:
            problems.append(f"{rel}: declared star-blind but references star_info")
    return problems


def check_math_lib(src_root):
    """Reproducibility contract (2026): every file calling an elementary
    transcendental must `use math_lib` (so USE_CRMATH builds shadow the
    intrinsics), and no real-exponent ** may exist (a hidden libm pow;
    write pow()/exp10()). Integer-literal exponents are exact and fine."""
    import re as _re
    ref = _re.compile(r"(?i)(?<![a-z0-9_])" + MATH_FUNCS + r"\s*\(")
    problems = []
    for f in sorted(src_root.rglob("*.f90")):
        rel = f.relative_to(src_root).as_posix()
        if "/test/" in rel or rel.startswith("math/"):
            continue
        raw = f.read_text(errors="replace")
        code_lines = [mask_strings(l) for l in raw.splitlines()]
        code = "\n".join(code_lines)
        if ref.search(code) and "use math_lib" not in raw:
            problems.append(f"{rel}: calls elementary math without `use math_lib`")
        for i, cl in enumerate(code_lines, 1):
            for m in _re.finditer(r"\*\*\s*(\(?\s*[+-]?\s*)([A-Za-z0-9_.]+)", cl):
                tok = m.group(2)
                if _re.fullmatch(r"\d+", tok):
                    continue
                if _re.fullmatch(r"\d+\.(?![0-9dDeE])", tok) and \
                        cl[m.end(2):].lstrip().startswith(("lt.", "gt.", "le.",
                                                           "ge.", "eq.", "ne.")):
                    continue   # maximal-munch: 2.lt. is integer 2 + .lt.
                problems.append(f"{rel}:{i}: real-exponent ** "
                                f"(use pow()/exp10()): ...{cl.strip()[:60]}")
    return problems


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

    buffer_violations = []
    for path in SRC.rglob("*.f90"):
        rel = path.relative_to(SRC).as_posix()
        if rel in CONTROLS_BUFFER_IMPORTERS or rel == "io/controls_lib.f90":
            continue
        code_lines = strip_comments(path.read_text()).splitlines()
        for i, line in enumerate(code_lines):
            if line.strip().lower().startswith("use controls_lib"):
                stmt = line
                while stmt.rstrip().endswith("&") and i + 1 < len(code_lines):
                    i += 1
                    stmt = stmt.rstrip()[:-1] + code_lines[i].lstrip().lstrip("&")
                m = re.match(r"^\s*use\s+controls_lib\s*,\s*only\s*:\s*(.*)$",
                             stmt, re.IGNORECASE)
                if m and all(CONTROLS_LIB_CONSTANT_RE.match(x.strip())
                             for x in m.group(1).split(",")):
                    continue
                buffer_violations.append(rel)
                break
    if buffer_violations:
        print("controls_lib buffer imported outside the read pipeline "
              "(read star%ctrl / star%job instead):")
        for rel in sorted(buffer_violations):
            print(f"  {rel}")
        return 1

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

    blind_problems = check_star_blind(SRC)
    if blind_problems:
        print("STAR-BLIND CONTRACT VIOLATIONS:")
        for b in blind_problems:
            print("  " + b)
        return 1

    math_problems = check_math_lib(SRC)
    if math_problems:
        print("MATH-LIB CONTRACT VIOLATIONS (reproducibility campaign):")
        for m in math_problems:
            print("  " + m)
        return 1

    print("Domain boundaries OK: every cross-domain call goes through "
          "a public entry; math-lib contract holds.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

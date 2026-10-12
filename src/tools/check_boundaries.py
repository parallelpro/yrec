#!/usr/bin/env python3
"""Domain-boundary checker for YREC (phase three, docs/ROADMAP.md stage 2).

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

DOMAINS = ["eos", "kap", "atm", "net", "convection",
           "rotation", "element_diffusion", "numerics"]

# The controls-read BUFFER (io/controls_buffer.f90) is reader-internal:
# star%ctrl / star%job are the authoritative post-read home, and every
# consumer reads them. Only the read pipeline itself may import the
# buffer's bare names -- read_controls (the reader), controls_sync
# (the generated buffer<->star copies), map_user_inputs (runs inside
# the read, before the store), and net's test harness (which seeds
# the buffer to simulate that pipeline). A new name here is a design
# regression, not an allowlist entry.
CONTROLS_BUFFER_IMPORTERS = {
    "io/read_controls.f90",
    "io/legacy_controls.f90",
    "io/run_files.f90",
    "io/derived_controls.f90",
    "io/run_cards.f90",
    "state/controls_sync.f90",
    "setup/user_inputs.f90",
    "net/test/test_net.f90",
}

# Named constants in controls_buffer that are not buffer state: any file
# may `use controls_buffer, only: <these>` (an only-list drawn entirely
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
            "eos_set_debye_huckel_z", "eos_set_debye_huckel_composition",
            "eos_update_envelope_tables", "set_eos_settings"},
    # The kap_lib facade entries (kap_get_r is the named-index
    # result-array variant of kap_get).
    # kap_get is the named-index result-array query (the former
    # kap_get_r); the long-form engine kap_eval is domain-internal.
    "kap": {"kap_get", "kap_init", "kap_update_surface_tables"},
    # atm_lib's entries: the table loads, the tabulated surface lookup,
    # the Allard Teff range, and the analytic T-tau relations the star
    # layer's envelope integrator uses.
    "atm": {"atm_init", "atm_get_surface_pt", "atm_allard_teff_max",
            "set_atm_settings", "ttau_log10_temperature",
            "ttau_start_log10_temperature", "ttau_photosphere_x_limit",
            "hsra_t_tau_offset"},
    # net/'s public procedures: what net_lib re-exports (2026-10 names;
    # Timmes' sneut/azbar/ifermi12/zfermim12 keep theirs).
    "net": {"compute_network_rates", "set_net_settings",
            "burn_implicitly", "compute_abundance_rates",
            "compute_nuclear_energy",
            "deuterium_rate", "burn_deuterium_region", "compute_li_be_rates"},
    # convection_lib's entries (2026-10: mixing/ became convection/; its
    # burn-settle-mix sub-stepping moved to rotation).
    "convection": {"compute_temperature_gradients", "find_convection_zones",
                   "homogenize_convection_zones", "extend_by_semiconvection",
                   "set_convection_settings", "set_convection_envelope_metals"},
    # rotation deliberately has no facade (multi-primitive surface, a
    # user decision, confirmed 2026-10-12): its modules are used by name,
    # and these are the procedures core and setup may call
    # (set_rotation_settings installs its settings record each run card).
    "rotation": {"evolve_angular_momentum", "omega_from_j",
                 "rotation_shape_factors", "zone_moments_of_inertia",
                 "am_convective_regions", "viscos", "set_rotation_settings"},
    # element diffusion (gravitational settling and microscopic
    # diffusion of the elements), split out of rotation in 2026-10: the
    # original routines and the newer ones, chosen by
    # use_new_diffusion_routines; diffuse_elements, the facade's one
    # entry, chooses.
    "element_diffusion": {"diffuse_elements", "reset_element_diffusion"},
    # numerics is the shared numerics library: everything numerics_lib
    # re-exports is public (FACADE_ONLY_DOMAINS keeps its five topic
    # modules behind it).
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
    "core/envelope/envelope_integration_kernel.f90",
    "atm/t_tau_relations.f90",
]

def check_star_blind(src_root):
    problems = []
    for rel in STAR_BLIND_FILES:
        text = (src_root / rel).read_text(errors="replace")
        code = "\n".join(l.split("!")[0] for l in text.splitlines())
        if "star_info_lib" in code or "star%" in code:
            problems.append(f"{rel}: declared star-blind but references star_info")
    return problems


# Physics domains whose code may not touch the star object (2026-10):
# their state comes in as arguments and their configuration as a
# settings record (net: net_lib's net_settings, filled by the star
# layer's configure_network). Index constants from star_info_lib (the
# composition slots i_*, json, ...) are allowed; `star` itself is not.
STAR_FREE_DIRS = ["net", "kap", "eos", "atm", "numerics", "convection"]

def check_star_free(src_root):
    problems = []
    use_star = re.compile(r"(?im)^\s*use\s+star_info_lib\b")
    for d in STAR_FREE_DIRS:
        for f in sorted((src_root / d).rglob("*.f90")):
            rel = f.relative_to(src_root).as_posix()
            if "/test/" in rel:
                continue
            code = strip_comments(f.read_text(errors="replace"))
            for m in use_star.finditer(code):
                problems.append(f"{rel}: uses star_info_lib (STAR_FREE_DIRS; the "
                                "slot constants are in composition_layout)")
            if re.search(r"(?i)\bstar%", code):
                problems.append(f"{rel}: references star% (STAR_FREE_DIRS)")
    return problems


# Domains that use the star object for their state but read their
# configuration from a settings record (2026-10): their code may not
# read star%ctrl or star%job, except the listed members, which are
# state shared with core rather than settings (rotation: the settling
# latch, which settling turns off when X or Y runs out).
CONTROLS_FREE_DIRS = {"rotation": {"job%diffuse_helium_active"}}

def check_controls_free(src_root):
    problems = []
    read = re.compile(r"(?i)\bstar%((?:ctrl|job)%\w+)")
    for d, allowed in CONTROLS_FREE_DIRS.items():
        for f in sorted((src_root / d).rglob("*.f90")):
            rel = f.relative_to(src_root).as_posix()
            if "/test/" in rel:
                continue
            for name in sorted(set(m.lower() for m in read.findall(
                    strip_comments(f.read_text(errors="replace"))))):
                if name not in allowed:
                    problems.append(f"{rel}: reads star%{name} (CONTROLS_FREE_DIRS; "
                                    f"add it to {d}'s settings record)")
    return problems


# Domains whose modules only their facade may expose (docs/GUIDELINES.md,
# naming rule 2): outside the domain's directory, `use` may name only
# <domain>_lib. Add a domain once its users go through its facade.
FACADE_ONLY_DOMAINS = ["net", "kap", "eos", "atm", "numerics", "convection",
                       "element_diffusion"]

def check_facade_imports(src_root):
    problems = []
    use_re = re.compile(r"(?im)^\s*use\s+(\w+)")
    mod_re = re.compile(r"(?im)^\s*module\s+(?!procedure)(\w+)")
    for d in FACADE_ONLY_DOMAINS:
        internal = set()
        for f in (src_root / d).rglob("*.f90"):
            internal |= {m.lower() for m in mod_re.findall(f.read_text(errors="replace"))}
        internal.discard(f"{d}_lib")
        for f in sorted(src_root.rglob("*.f90")):
            rel = f.relative_to(src_root).as_posix()
            if rel.startswith(d + "/"):
                continue
            for m in use_re.findall(strip_comments(f.read_text(errors="replace"))):
                if m.lower() in internal:
                    problems.append(f"{rel}: uses {m}, internal to {d}/ (use {d}_lib)")
    return problems


def check_math_lib(src_root):
    """Reproducibility contract (2026): every file calling an elementary
    transcendental must `use math_lib` (so USE_CRMATH builds shadow the
    intrinsics), and no real-exponent ** may exist (a hidden libm pow;
    write pow()/exp10()). Integer-literal exponents are exact and fine.
    Parameter initializers are exempt: the compiler folds them, correctly
    rounded, and a math_lib procedure may not appear in a constant
    expression."""
    import re as _re
    ref = _re.compile(r"(?i)(?<![a-z0-9_])" + MATH_FUNCS + r"\s*\(")
    param_decl = _re.compile(r"(?i)^\s*(double\s+precision|real|integer)\b[^:]*\bparameter\b[^:]*::")

    def without_parameter_decls(lines):
        kept, stmt = [], []
        for line in lines:
            stmt.append(line)
            if line.rstrip().endswith("&"):
                continue
            if not param_decl.search(" ".join(stmt)):
                kept.extend(stmt)
            stmt = []
        return kept + stmt

    problems = []
    for f in sorted(src_root.rglob("*.f90")):
        rel = f.relative_to(src_root).as_posix()
        if "/test/" in rel or rel.startswith("math/"):
            continue
        raw = f.read_text(errors="replace")
        code_lines = [mask_strings(l) for l in raw.splitlines()]
        code_kept = without_parameter_decls(code_lines)
        # Scope by scope: a procedure that calls elementary math must see
        # math_lib through its own use or a host's (module, host
        # procedure). A sibling internal procedure does not see another
        # sibling's use, so a file-level test would miss it.
        stack = []
        for i, cl in enumerate(code_kept, 1):
            line = cl.split("!")[0]
            low = line.strip().lower()
            if not low:
                continue
            if _re.match(r"end(\s*(subroutine|function|module|program)\b.*)?$", low):
                if stack:
                    stack.pop()
                continue
            if _re.match(r"(module|program)\s+\w+\s*$", low) and \
                    not low.startswith("module procedure"):
                stack.append([low.split()[1], False])
                continue
            m = _re.match(r"(?:(?:recursive|pure|elemental|impure)\s+)*"
                          r"(?:(?:double\s+precision|integer|logical|real|"
                          r"complex|character)(?:\s*\([^)]*\))?\s+)?"
                          r"(subroutine|function)\s+(\w+)", low)
            if m:
                stack.append([m.group(2), False])
                continue
            if _re.match(r"use\s+math_lib\b", low):
                if stack:
                    stack[-1][1] = True
                continue
            if ref.search(line) and stack and not any(s[1] for s in stack):
                problems.append(f"{rel}:{i}: elementary math in {stack[-1][0]} "
                                "without `use math_lib` in scope")
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
            if "test" in path.relative_to(SRC).parts:
                continue   # a test program's own procedures are not domain API
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
        if rel in CONTROLS_BUFFER_IMPORTERS or rel == "io/controls_buffer.f90":
            continue
        code_lines = strip_comments(path.read_text()).splitlines()
        for i, line in enumerate(code_lines):
# legacy_controls re-exports the buffer (its namelists bind to it)
            if re.match(r"(?i)^\s*use\s+legacy_controls\b", line):
                buffer_violations.append(rel)
                break
            if line.strip().lower().startswith("use controls_buffer"):
                stmt = line
                while stmt.rstrip().endswith("&") and i + 1 < len(code_lines):
                    i += 1
                    stmt = stmt.rstrip()[:-1] + code_lines[i].lstrip().lstrip("&")
                m = re.match(r"^\s*use\s+controls_buffer\s*,\s*only\s*:\s*(.*)$",
                             stmt, re.IGNORECASE)
                if m and all(CONTROLS_LIB_CONSTANT_RE.match(x.strip())
                             for x in m.group(1).split(",")):
                    continue
                buffer_violations.append(rel)
                break
    if buffer_violations:
        print("controls_buffer buffer imported outside the read pipeline "
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

    free_problems = check_star_free(SRC)
    if free_problems:
        print("STAR-FREE CONTRACT VIOLATIONS:")
        for b in free_problems:
            print("  " + b)
        return 1

    controls_problems = check_controls_free(SRC)
    if controls_problems:
        print("CONTROLS-FREE CONTRACT VIOLATIONS:")
        for b in controls_problems:
            print("  " + b)
        return 1

    facade_problems = check_facade_imports(SRC)
    if facade_problems:
        print("FACADE CONTRACT VIOLATIONS:")
        for b in facade_problems:
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

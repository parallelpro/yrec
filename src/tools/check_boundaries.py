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
# star%control and star%physics are the authoritative post-read home, and
# every consumer reads them (or the run's working copy of a control the
# run changes). Only the read pipeline itself may import the
# buffer's bare names -- read_controls (the reader), controls_sync
# (the generated buffer<->star copies), map_user_inputs (runs inside
# the read, before the store), and net's test harness (which seeds
# the buffer to simulate that pipeline). A new name here is a design
# regression, not an allowlist entry.
CONTROLS_BUFFER_IMPORTERS = {
    "io/controls_reader.f90",
    "io/legacy_controls.f90",
    "io/derived_controls.f90",
    "io/run_cards.f90",
    "io/controls_check.f90",
    "io/controls_sync.f90",
    "io/user_inputs.f90",
    "net/test/test_net.f90",
    "io/test/test_controls.f90",
}

# Named constants in controls_buffer that are not buffer state: any file
# may `use controls_buffer, only: <these>` (an only-list drawn entirely
# from these is not a buffer import). max_runs dimensions the per-run
# arrays; ichi_*/itime_* name the chi_grid_scale/atime slots (declared in
# state/control_slots.f90, which new code uses; the buffer re-exports
# them).
CONTROLS_LIB_CONSTANT_RE = re.compile(r"^(max_runs|ichi_\w+|itime_\w+)$")

# domain -> names callable from outside that domain.
PUBLIC = {
    # eos_lib's facade: eos_init installs the settings and reads the
    # tables, eos_get evaluates a zone into an eos_result, eos_set_mixture
    # takes the envelope composition (and builds the SCV envelope table).
    # The blended, Yale, SCV, OPAL and MHD modules are eos-internal (the
    # test-only gamma1_at_density and set_eos_settings are used by
    # eos/test/test_eos.f90 from their own modules).
    "eos": {"eos_get", "eos_init", "eos_set_mixture"},
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
    # cross_section_scales_from_sstandard: the controls read maps the
    # legacy sstandard vector onto the reaction scales with it.
    "net": {"compute_network_rates", "set_net_settings",
            "cross_section_scales_from_sstandard",
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
    "element_diffusion": {"diffuse_elements", "set_diffusion_settings"},
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



# The elementary functions math_lib shadows under the crmath backend: the
# names on its public lines inside `#ifdef YREC_CRMATH` (math/README.md,
# section 4). A caller of one must have `use math_lib` in scope.
def shadowed_math_functions(src_root):
    names, in_crmath = [], False
    for line in (src_root / "math" / "math_lib.f90").read_text().splitlines():
        directive = line.strip().lower()
        if directive.startswith("#ifdef") and "yrec_crmath" in directive:
            in_crmath = True
        elif directive.startswith(("#else", "#endif")):
            in_crmath = False
        elif in_crmath:
            m = re.match(r"(?i)^\s*public\s*::(.*)$", line)
            if m:
                names += [n.strip().lower() for n in m.group(1).split(",") if n.strip()]
    if not names:
        raise SystemExit("check_boundaries: no shadowed functions found in math/math_lib.f90")
    return names


def math_funcs_pattern(src_root):
    names = sorted(shadowed_math_functions(src_root), key=len, reverse=True)
    return r"(?:" + "|".join(names) + r")"

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
    # the structure equations at one shell (2026-10, GUIDELINES star rule 6)
    "core/structure/shell_equations.f90",
    # the atmosphere and envelope integrands, called back by the kernel
    "core/envelope/envelope_integrands.f90",
    # the envelope triangle's geometry and fit (2026-10, star rule 6)
    "core/structure/triangle_fit.f90",
    # the mass change's masses, wind angular momentum and accreted
    # composition (2026-10, GUIDELINES star rule 6)
    "core/step/mass_change_kernel.f90",
    # the rezone's regridding arithmetic (2026-10, GUIDELINES star rule 6)
    "core/model/regridding.f90",
    "atm/t_tau_relations.f90",
    # rotation's kernels (2026-10, star rule 6; test_rotation_kernels):
    # the rotation curve from j, the moments of inertia, the shape and
    # its quadratures, the angular momentum diffusion system and its
    # solve, the wind laws, the viscosity, the convection zones for the
    # transport
    "rotation/omega_profile.f90",
    "rotation/shape/moments_of_inertia.f90",
    "rotation/shape/rotating_shape.f90",
    "rotation/shape/shell_quadrature.f90",
    "rotation/seculr/am_transport_coefficients.f90",
    "rotation/magnetic_braking.f90",
    "rotation/viscosity.f90",
    "rotation/am_regions.f90",
    # the instability velocities and the transport-region scan (2026-10)
    "rotation/seculr/instability_velocities.f90",
    "rotation/seculr/transport_regions.f90",
    # one region's angular momentum diffusion: the system, its solve and
    # the map back (2026-10-28)
    "rotation/seculr/region_transport.f90",
    # the mechanisms' velocities, slot by slot (2026-10-28)
    "rotation/seculr/transport_mechanisms.f90",
    # the equation of state's solvers and interpolators (2026-10-29;
    # test_eos): the blend, the mixture, the Yale, SCV, OPAL and MHD EOS
    "eos/blended_eos.f90",
    "eos/eos_mixture.f90",
    "eos/yale/yale_eos.f90",
    "eos/scv/scv_eos.f90",
    "eos/opal/opal_eos.f90",
    "eos/opal/opal95_eos.f90",
    "eos/mhd/mhd_eos.f90",
]

# Scratch modules that only their own domain may use (2026-10): their
# records are working arrays with a one-step lifetime, and what the star
# layer needs from them goes through arguments or the star.
SCRATCH_OWNERS = {"rotation_scratch_lib": "rotation",
                  "element_diffusion_scratch": "element_diffusion"}

def check_scratch_private(src_root):
    problems = []
    use_re = re.compile(r"(?im)^\s*use\s+(\w+)")
    for f in sorted(src_root.rglob("*.f90")):
        rel = f.relative_to(src_root).as_posix()
        for m in use_re.finditer(strip_comments(f.read_text(errors="replace"))):
            owner = SCRATCH_OWNERS.get(m.group(1).lower())
            if owner and not rel.startswith(owner + "/"):
                problems.append(f"{rel}: uses {m.group(1)} ({owner}'s scratch)")
    return problems


# Scratch modules that declare types only (2026-10): their owner's driver
# keeps the instances for the step and passes them down as arguments, so
# the module holds no state at all. Any declaration outside a type
# block, other than a parameter, is a module variable.
SCRATCH_TYPES_ONLY = ["rotation/rotation_scratch_lib.f90",
                      "element_diffusion/element_diffusion_scratch.f90"]

def check_scratch_types_only(src_root):
    problems = []
    for rel in SCRATCH_TYPES_ONLY:
        code = strip_comments((src_root / rel).read_text(errors="replace"))
        depth = 0
        for line in join_continuations(code):
            s = line.strip().lower()
            if re.match(r"type\b\s*(,.*)?::", s) or re.match(r"type\s+\w+\s*$", s):
                depth += 1
                continue
            if re.match(r"end\s*type\b", s):
                depth -= 1
                continue
            if s.startswith("contains"):
                break
            if depth == 0 and "::" in s and not re.match(r"(public|private)\b", s) \
                    and "parameter" not in s:
                problems.append(f"{rel}: module variable ({line.strip()})")
    return problems


def check_star_blind(src_root):
    problems = []
    for rel in STAR_BLIND_FILES:
        text = (src_root / rel).read_text(errors="replace")
        code = "\n".join(l.split("!")[0] for l in text.splitlines())
        if "star_info_lib" in code or "star%" in code:
            problems.append(f"{rel}: declared star-blind but references star_info")
    return problems


# log10 G (2026-10-29): physical_constants' cgl is the star layer's
# (run_setup sets it from the control G_cgs). The physics domains, rotation
# and the star-blind kernels take it from their settings records
# (log10_gravitational_constant), so they give the right answer without
# run_setup having run.
GLOBAL_G_RE = re.compile(r"(?i)(?<![\w%])cgl\b")

def check_no_global_g(src_root):
    problems = []
    files = set(STAR_BLIND_FILES)
    for d in STAR_FREE_DIRS + sorted(CONTROLS_FREE_DIRS):
        files |= {f.relative_to(src_root).as_posix() for f in (src_root / d).rglob("*.f90")
                  if "/test/" not in f.relative_to(src_root).as_posix()}
    for rel in sorted(files):
        code = strip_comments((src_root / rel).read_text(errors="replace"))
        if GLOBAL_G_RE.search(code):
            problems.append(f"{rel}: reads physical_constants' cgl (take log10 G from "
                            "the settings record)")
    return problems


# The settings install (2026-10-29): every star-layer call of a domain's
# installer (set_<domain>_settings, the kap_init, eos_init and atm_init of
# the domains that also load tables, and the setters of the state that
# follows the surface mixture: set_convection_envelope_metals,
# eos_set_mixture, kap_update_surface_tables) is in
# setup/settings_install.f90, one procedure per lifetime (per run, per
# Monte-Carlo run, per card, per surface mixture), so that file is the
# list of what is installed when. A domain's
# calls inside its own directory are its own business, and a set_*_settings
# the star layer defines for itself (core's set_structure_settings, the
# shell equations' record filled at every call) is not a domain's. The
# test programs install records themselves.
INSTALL_FILES = ["setup/settings_install.f90"]
INSTALLER_CALL_RE = re.compile(
    r"(?i)\bcall\s+(set_\w+_settings|set_convection_envelope_metals|kap_init|eos_init|atm_init|"
    r"eos_set_mixture|kap_update_surface_tables)\s*\(")

def check_install_sites(src_root):
    layer = [f for f in sorted(src_root.rglob("*.f90"))
             if f.relative_to(src_root).parts[0] not in DOMAINS + ["math"]
             and "/test/" not in f.relative_to(src_root).as_posix()]
    own = set()
    for f in layer:
        own |= {m.group(1).lower() for m in re.finditer(
            r"(?im)^\s*subroutine\s+(set_\w+_settings)\b", strip_comments(f.read_text(errors="replace")))}
    problems = []
    for f in layer:
        rel = f.relative_to(src_root).as_posix()
        if rel in INSTALL_FILES:
            continue
        for m in INSTALLER_CALL_RE.finditer(strip_comments(f.read_text(errors="replace"))):
            if m.group(1).lower() not in own:
                problems.append(f"{rel}: calls {m.group(1)} (install settings in "
                                "setup/settings_install.f90)")
    return problems


# Every field of a domain's settings record is set by its builder
# (2026-10-29): a field the builder does not set keeps its declared
# default for every run, which is how a newly added field is forgotten.
# The builders are setup/'s functions that return a record of a
# *_settings_t type (`type(<t>) :: settings` in the function); a field is
# set when the function assigns settings%<field> (or an element of it).
# A field that is deliberately left at its default is listed here, with
# the reason.
SETTINGS_DEFAULT_ONLY = {}

def settings_type_fields(src_root):
    """{type name: [field, ...]} for the *_settings_t types the domains
    declare"""
    out = {}
    for d in DOMAINS:
        for f in sorted((src_root / d).rglob("*.f90")):
            if "/test/" in f.relative_to(src_root).as_posix():
                continue
            code = "\n".join(join_continuations(strip_comments(f.read_text(errors="replace"))))
            for m in re.finditer(r"(?ims)^\s*type\s*(?:,\s*public\s*)?(?:::)?\s*(\w+_settings_t)\s*$(.*?)^\s*end\s*type\b",
                                 code):
                fields = []
                for line in m.group(2).splitlines():
                    if "::" not in line:
                        continue
                    decl = line.split("::", 1)[1]
                    while re.search(r"\([^()]*\)", decl):
                        decl = re.sub(r"\([^()]*\)", "", decl)
                    decl = re.sub(r"\[[^\]]*\]", "", decl)
                    for item in decl.split(","):
                        name = item.split("=")[0].strip().lower()
                        if re.match(r"^\w+$", name):
                            fields.append(name)
                out[m.group(1).lower()] = fields
    return out

def check_settings_fields(src_root):
    types = settings_type_fields(src_root)
    built = {}
    for f in sorted((src_root / "setup").glob("*.f90")):
        code = "\n".join(join_continuations(strip_comments(f.read_text(errors="replace"))))
        for m in re.finditer(r"(?ims)^\s*(?:(?:pure|elemental|impure|recursive)\s+)*function\s+(\w+)\s*\([^)]*\)\s*result\s*\(\s*settings\s*\)(.*?)^\s*end\s*function\b",
                             code):
            t = re.search(r"(?im)^\s*type\s*\(\s*(\w+_settings_t)\s*\)\s*::\s*settings\s*$", m.group(2))
            if t:
                built[t.group(1).lower()] = (m.group(1), {a.lower() for a in re.findall(
                    r"(?im)^\s*(?:if\s*\(.*\)\s*)?settings%(\w+)(?:\s*\([^=]*\))?\s*=(?!=)", m.group(2))})
    problems = []
    for t, fields in sorted(types.items()):
        if t not in built:
            problems.append(f"{t}: no builder in setup/ returns it")
            continue
        name, assigned = built[t]
        for fld in fields:
            if fld not in assigned and fld not in SETTINGS_DEFAULT_ONLY.get(t, {}):
                problems.append(f"{t}%{fld}: {name} does not set it (copy its control, or list it "
                                "in SETTINGS_DEFAULT_ONLY with the reason)")
    return problems


# Physics domains whose code may not touch the star object (2026-10):
# their state comes in as arguments and their configuration as a
# settings record (net: net_lib's net_settings, filled by the star
# layer's net_settings_from_star, setup/). Index constants from star_info_lib (the
# composition slots i_*, json, ...) are allowed; `star` itself is not.
STAR_FREE_DIRS = ["net", "kap", "eos", "atm", "numerics", "convection",
                  "element_diffusion"]

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
# read star%control or star%physics, except the listed members, which are state shared with
# core rather than settings (none left: the settling latch is
# star%settling_active since 2026-10-29).
CONTROLS_FREE_DIRS = {"rotation": set()}

def check_controls_free(src_root):
    problems = []
    # the global star, or the star as an argument s
    read = re.compile(r"(?i)\b(?:star|s)%((?:control|physics)%\w+)")
    for d, allowed in CONTROLS_FREE_DIRS.items():
        for f in sorted((src_root / d).rglob("*.f90")):
            rel = f.relative_to(src_root).as_posix()
            if "/test/" in rel:
                continue
            for name in sorted(set(m.lower() for m in read.findall(
                    strip_comments(f.read_text(errors="replace"))))):
                if name not in allowed:
                    problems.append(f"{rel}: reads s%{name} (CONTROLS_FREE_DIRS; "
                                    f"add it to {d}'s settings record)")
    return problems


# The controls' homes and uses (2026-10-29, run copies 2026-10-30):
# defaults/controls_registry.tsv says, for each controls_buffer member,
# whether it lives in the star's record of its group (star%control,
# star%physics) or only in the read's buffer, whether the run keeps a
# working copy of it (run_copy, a star member), and where it is used. Code
# outside the controls read never assigns a control: a control the run
# changes has a run copy, which the run assigns instead. The records hold
# the controls that code outside the read reads, directly or through
# their run copies; `used` is run when it does, read when only the read
# does (an input-only control: home buffer), none when nothing does (home
# buffer). A run copy is assigned outside the read (otherwise the code
# reads the control). Tests are not code here. A write is an assignment to
# s%<x> (or star%...), or the member passed to an intent(out) or
# intent(inout) dummy (tools/check_star_contracts.py's analysis).
CONTROLS_READ_FILES = {
    "io/controls_reader.f90", "io/legacy_controls.f90",
    "io/derived_controls.f90", "io/run_cards.f90", "io/controls_buffer.f90",
    "io/controls_sync.f90", "io/user_inputs.f90",
}


def check_control_homes(src_root):
    import controls_registry
    import check_star_contracts as contracts
    rows = [r for r in controls_registry.rows() if r["home"]]
    outside = [p for p in sorted(src_root.rglob("*.f90"))
               if "/test/" not in p.relative_to(src_root).as_posix()
               and p.relative_to(src_root).as_posix() not in CONTROLS_READ_FILES]
    code = "\n".join(strip_comments(p.read_text(errors="replace")) for p in outside).lower()
    written = set(m.group(1) for m in re.finditer(
        r"\b(?:s|star)%(?:control|physics)%(\w+)\s*(?:\([^=\n]*\))?\s*=(?!=)", code))
    star_written = set(m.group(1) for m in re.finditer(
        r"\b(?:s|star)%(\w+)\s*(?:\([^=\n]*\))?\s*=(?!=)", code))
    procs = contracts.scan(outside)
    for name in procs:
        _, writes = contracts.footprint(procs, name)
        written |= {w.split("%", 1)[1] for w in writes if w.startswith(("control%", "physics%"))}
        star_written |= {w for w in writes if "%" not in w}
    referenced = set(re.findall(r"\b(?:s|star)%(?:control|physics)%(\w+)", code))
    read_files = [src_root / f for f in sorted(CONTROLS_READ_FILES)
                  if f not in ("io/controls_buffer.f90", "io/controls_sync.f90")
                  and (src_root / f).exists()]
    # a namelist group's list of names, or a use statement's only-list, is
    # not a use
    read_code = "\n".join(l for p in read_files
                          for l in join_continuations(strip_comments(p.read_text(errors="replace")))
                          if not re.match(r"(?i)^\s*(namelist|use)\b", l)).lower()
    problems = []
    for r in rows:
        m = r["canonical"]
        if m in written:
            problems.append(f"{m}: code outside the read assigns the control; give it a run "
                            "copy (registry run_copy, a star member) and assign that")
        copy = r["run_copy"]
        if copy and copy not in star_written:
            problems.append(f"{m}: its run copy {copy} is never assigned outside the read "
                            "(read the control instead, and drop run_copy)")
        home = r["group"] if m in referenced or copy else "buffer"
        if home != r["home"]:
            problems.append(f"{m}: registry home {r['home']}, but code outside the read "
                            + ("never reads it (buffer)" if home == "buffer" else
                               f"reads it (the record of its group, {home})"))
        # the read's own copies (m = legacy, m = s%<group>%m) are not uses
        rest = re.sub(r"(?m)^\s*" + m + r"\s*(\([^=\n]*\))?\s*=.*$", "", read_code)
        rest = re.sub(r"(?m)^\s*" + re.escape(r["legacy"].lower()) + r"\s*=\s*" + m + r"\s*$", "", rest)
        used = ("run" if m in referenced or copy else
                "read" if re.search(r"(?<![%\w])" + m + r"\b|%(?:control|physics)%" + m + r"\b", rest)
                else "none")
        if used != r["used"]:
            problems.append(f"{m}: registry used {r['used']}, but the code says {used}")
    return sorted(set(problems))

# The controls' records (2026-10-30). The star holds the deck's controls in
# two records named for the inlist's namelist groups, star%control (the
# &control group) and star%physics (the &physics group), each control in the
# record of its registry group; the controls read alone writes them
# (io/controls_sync.f90, generated by tools/gen_controls_state.py). The
# values the run changes are not controls: each is a star member of its own,
# the control's run copy (the registry's run_copy column), which the read
# arms from the control. The records' names before 2026-10-30 (ctrl, and job
# for the run-changed part) must not come back, in code or in prose.
OLD_RECORDS = ("ctrl", "job")
OLD_RECORD_RE = re.compile(r"(?i)(?<![\w%])(?:s|star)%(?:" + "|".join(OLD_RECORDS) + r")\b")
# the bare form of a contract line or a comment (<old record>%<member>), and
# an old record under any other name of the star (a test's probe)
OLD_MEMBER_RE = re.compile(r"(?i)(?<![\w%])(?:" + "|".join(OLD_RECORDS) + r")%\w"
                           r"|%(?:" + "|".join(OLD_RECORDS) + r")%\w")
TEXT_SUFFIXES = {".f90", ".inc", ".txt", ".py", ".md", ".tsv", ".defaults", ".list"}
CONTROLS_WRITER = "io/controls_sync.f90"
CONTROL_RECORDS = {"control": "state/control_group_def.inc",
                   "physics": "state/physics_group_def.inc"}


def check_old_records(src_root):
    """The old record names (OLD_RECORDS) after the star (s or
    star) appear nowhere in src/, and neither the bare <old record>%<name>
    nor <anything>%<old record>%<name> in the Fortran sources or the
    READMEs"""
    problems = []
    for f in sorted(src_root.rglob("*")):
        if not f.is_file() or f.suffix not in TEXT_SUFFIXES:
            continue
        rel = f.relative_to(src_root).as_posix()
        bare = f.suffix in (".f90", ".inc", ".md")
        for i, line in enumerate(f.read_text(errors="replace").splitlines(), 1):
            if OLD_RECORD_RE.search(line) or (bare and OLD_MEMBER_RE.search(line)):
                problems.append(f"{rel}:{i}: names an old controls record: {line.strip()[:70]}")
    return problems


def check_controls_writer(src_root):
    """star%control and star%physics are assigned only in
    io/controls_sync.f90 (tests aside): no assignment, read target,
    allocation or intent(out/inout) actual argument elsewhere"""
    import check_star_contracts as contracts
    files = [p for p in sorted(src_root.rglob("*.f90"))
             if "/test/" not in p.relative_to(src_root).as_posix()
             and p.relative_to(src_root).as_posix() != CONTROLS_WRITER]
    problems = []
    assign = re.compile(r"(?i)^\s*(?:if\s*\(.*\)\s*)?(?:s|star)%(control|physics)\b"
                        r"(?:%\w+)?(?:\s*\([^=]*\))?\s*=(?!=)")
    for p in files:
        rel = p.relative_to(src_root).as_posix()
        for stmt in join_continuations(strip_comments(p.read_text(errors="replace"))):
            if assign.match(stmt):
                problems.append(f"{rel}: assigns s%{assign.match(stmt).group(1)} "
                                f"({stmt.strip()[:60]})")
    procs = contracts.scan(files)
    for key in procs:
        _, writes = contracts.footprint(procs, key)
        for w in sorted(writes):
            if w.split("%")[0] in CONTROL_RECORDS:
                problems.append(f"{procs[key]['path'].relative_to(src_root).as_posix()}: "
                                f"{key[-1]} writes {w} (an intent(out/inout) argument or a read)")
    return sorted(set(problems))


# The run cards' capacity, max_runs, is state/control_slots.f90's
# parameter; tools/controls_registry.py repeats it (MAX_RUNS) to read the
# registry's per-run types, and the two must agree.
def check_max_runs(src_root):
    import controls_registry
    m = re.search(r"(?im)^\s*integer\s*,\s*parameter\s*,\s*public\s*::\s*max_runs\s*=\s*(\d+)",
                  (src_root / "state" / "control_slots.f90").read_text())
    if not m:
        return ["state/control_slots.f90: no max_runs parameter"]
    if int(m.group(1)) != controls_registry.MAX_RUNS:
        return [f"state/control_slots.f90 has max_runs = {m.group(1)}, tools/controls_registry.py "
                f"MAX_RUNS = {controls_registry.MAX_RUNS}"]
    return []


def check_control_records(src_root):
    """Every control is in the record of its registry group: the
    registry's home is the group (or buffer), and the generated component
    lists hold exactly the controls of their group"""
    import controls_registry
    rows = [r for r in controls_registry.rows() if r["home"]]
    problems = [f"{r['canonical']}: registry home {r['home']}, group {r['group']}"
                for r in rows if r["home"] not in ("buffer", r["group"])]
    want = {g: {r["canonical"].lower() for r in rows if r["home"] == g} for g in CONTROL_RECORDS}
    decl = re.compile(r"^\s*(?:double precision|real\s*\(\s*dp\s*\)|integer|logical|character\(len=\d+\))"
                      r"\s*::\s*(\w+)", re.I)
    for g, rel in CONTROL_RECORDS.items():
        have = {m.group(1).lower() for m in map(decl.match, (src_root / rel).read_text().splitlines()) if m}
        for n in sorted(have - want[g]):
            problems.append(f"{rel}: {n} is in star%{g}, but its registry home is not {g}")
        for n in sorted(want[g] - have):
            problems.append(f"{rel}: {n} (registry home {g}) is missing from star%{g}")
    return problems


def star_member_types(src_root):
    """{member: registry-style type (real(50), logical, ...)} of star_info's
    one-line declarations"""
    import controls_registry
    decl = re.compile(r"(?i)^\s*(double precision|real\s*\(\s*dp\s*\)|integer|logical|"
                      r"character\s*\(\s*len\s*=\s*\d+\s*\))"
                      r"\s*::\s*(\w+)\s*(\([^)]*\))?")
    out = {}
    for line in (src_root / "state" / "star_info_lib.f90").read_text().splitlines():
        m = decl.match(line)
        if m:
            out[m.group(2).lower()] = controls_registry.type_of(m.group(1), m.group(3) or "")
    return out


def check_control_copies(src_root):
    """A star member that copies a control says so, where it is set.
    A run copy (registry run_copy) is declared in star_info_lib with the
    control's type and shape, its comment names the control
    (<group>%<name>) and controls_sync among its writers, and
    controls_sync's arm_run_copies sets it from the control, once and
    nowhere else (a writer the comment names may restore it from its own
    control). Any other copy of a control into a whole star member is
    in a module the member's comment names as a writer, and the comment
    names the control: outside the read and the tests s%m = s%<group>%x;
    in the read's own files (CONTROLS_READ_FILES but controls_sync)
    s%m = x, with x the control's buffer name or its legacy name (READ_COPY
    below). A run copy is set from its control by arm_run_copies alone, in
    the read as outside it (but for that restore by a named writer)."""
    import controls_registry
    import gen_star_docs
    _, members, _ = gen_star_docs.parse_star_info(src_root / "state" / "star_info_lib.f90")
    doc = {e["name"].lower(): e for e in members}
    types = star_member_types(src_root)
    problems = []
    arm = (src_root / CONTROLS_WRITER).read_text().lower()
    arm = arm[arm.find("subroutine arm_run_copies"):arm.find("end subroutine arm_run_copies")]

    def names_control(e, ref):
        text = (e["meaning"] + " " + e["notes"]).lower()
        return re.search(r"(?<![\w%])" + re.escape(ref.lower()) + r"\b", text) is not None

    copies = {}
    for r in controls_registry.rows():
        c = r["run_copy"].lower()
        if not c:
            continue
        ref = f"{r['home']}%{r['canonical']}".lower()
        copies[c] = ref
        e = doc.get(c)
        if e is None:
            problems.append(f"{r['canonical']}: its run copy {c} is not a star member")
            continue
        if types.get(c) != r["type"].replace(" ", ""):
            problems.append(f"{c}: declared {types.get(c)}, but its control {ref} is {r['type']}")
        if not names_control(e, ref):
            problems.append(f"{c}: its comment does not name its control {ref}")
        if "controls_sync" not in [w.lower() for w in e["writers"]]:
            problems.append(f"{c}: its comment does not name controls_sync (arm_run_copies) "
                            "as a writer")
        n = len(re.findall(r"(?m)^\s*s%" + re.escape(c) + r"\s*=\s*s%" + re.escape(ref) + r"\s*$", arm))
        if n != 1:
            problems.append(f"{c}: arm_run_copies sets it from {ref} {n} times, not once")
    copy_stmt = re.compile(r"(?i)^\s*(?:if\s*\(.*\)\s*)?s%(\w+)\s*=\s*s%((?:control|physics)%\w+)\s*$")
    module_re = re.compile(r"(?im)^\s*module\s+(?!procedure)(\w+)")
    for p in sorted(src_root.rglob("*.f90")):
        rel = p.relative_to(src_root).as_posix()
        if "/test/" in rel or rel in CONTROLS_READ_FILES:
            continue
        text = p.read_text(errors="replace")
        mods = module_re.findall(text)
        module = mods[0].lower() if mods else p.stem.lower()
        for stmt in join_continuations(strip_comments(text)):
            m = copy_stmt.match(stmt)
            if not m:
                continue
            member, ref = m.group(1).lower(), m.group(2).lower()
            e = doc.get(member)
            if e is None:
                continue
            if member in copies:
                # restoring a run copy to its own control (the surface boundary's
                # switch back to the deck's atmosphere) is a writer's business
                if ref == copies[member] and module in [w.lower() for w in e["writers"]]:
                    continue
                problems.append(f"{rel}: sets the run copy {member} from {ref}; only the "
                                "controls read arms it (arm_run_copies), and a writer its "
                                "comment names may restore it from its own control")
                continue
            if module not in [w.lower() for w in e["writers"]]:
                problems.append(f"{rel}: copies {ref} into {member}, but {member}'s comment does "
                                f"not name {module} as a writer")
            if not names_control(e, ref):
                problems.append(f"{rel}: copies {ref} into {member}, but {member}'s comment does "
                                f"not name {ref}")
    # the read's own copies: s%m = <the control's buffer or legacy name>
    by_name = {}
    for r in controls_registry.rows():
        if not r["canonical"]:
            continue
        ref = f"{r['group']}%{r['canonical']}".lower()
        by_name[r["canonical"].lower()] = ref
        if r["legacy"] and r["legacy"] != "-":
            by_name.setdefault(r["legacy"].lower(), ref)
    for rel in sorted(CONTROLS_READ_FILES - {CONTROLS_WRITER}):
        p = src_root / rel
        if not p.exists():
            continue
        text = p.read_text(errors="replace")
        mods = module_re.findall(text)
        module = mods[0].lower() if mods else p.stem.lower()
        for stmt in join_continuations(strip_comments(text)):
            m = READ_COPY.match(stmt)
            if not m or m.group(2).lower() not in by_name or m.group(1).lower() not in doc:
                continue
            member, ref = m.group(1).lower(), by_name[m.group(2).lower()]
            e = doc[member]
            if member in copies:
                problems.append(f"{rel}: sets the run copy {member} from {m.group(2)} ({ref}); "
                                "only arm_run_copies arms it, from the control")
                continue
            if module not in [w.lower() for w in e["writers"]]:
                problems.append(f"{rel}: copies {ref} ({m.group(2)}) into {member}, but {member}'s "
                                f"comment does not name {module} as a writer")
            if not names_control(e, ref):
                problems.append(f"{rel}: copies {ref} ({m.group(2)}) into {member}, but {member}'s "
                                f"comment does not name {ref}")
    return sorted(set(problems))


# a copy of one name into a whole star member, s%m = x, inside the read
READ_COPY = re.compile(r"(?i)^\s*(?:if\s*\(.*\)\s*)?s%(\w+)\s*=\s*(\w+)\s*$")

# The star layer's access to the star object (2026-10). Routines come in
# three kinds, never mixed in one routine:
#   orchestration (core's drivers, io's writers): reads and writes the
#     star, which it takes as an argument s (run_yrec owns the run's
#     star; star_info_lib declares no instance since 2026-10);
#   kernels (the physics and numerics inside core and rotation): strict
#     intent(in)/intent(out) arguments, never the star -- listed in
#     STAR_BLIND_FILES above once extracted;
#   domains (eos, kap, atm, net, convection, numerics): star-free
#     (STAR_FREE_DIRS).
# STAR_ACCESS_BUDGET is a ratchet on the rest: per directory, the number
# of star% references and of files that import a module-level star may
# not grow. Every count reached zero in 2026-10 (the star is an argument
# s everywhere), so the budget now keeps a global star from coming back.
STAR_ACCESS_BUDGET = {
    # directory: (star% references, files importing the global star)
    # core (2026-10): run_yrec owns the star; everything takes it as s
    "core": (0, 0),
    # io (2026-10): the readers and writers take the star as an argument s
    "io": (0, 0),
    # rotation (2026-10): the drivers take the star as an argument s; the
    # rest is star-free
    "rotation": (0, 0),
    # setup (2026-10): the settings builders and map_user_inputs take the
    # star as an argument s
    "setup": (0, 0),
    # state (2026-10): controls_sync takes the star as an argument s
    "state": (0, 0),
}


def join_continuations(code):
    """Join free-form continuation lines (trailing &, optional leading &)."""
    out, cur = [], ""
    for line in code.splitlines():
        s = line.rstrip()
        if cur:
            s = s.lstrip()
            if s.startswith("&"):
                s = s[1:]
        if s.endswith("&"):
            cur += s[:-1] + " "
            continue
        out.append(cur + s)
        cur = ""
    if cur:
        out.append(cur)
    return out


def imports_global_star(code):
    """True if a use statement brings star_info_lib's star into scope."""
    for stmt in join_continuations(code):
        m = re.match(r"(?i)^\s*use\s+star_info_lib\b\s*(,\s*only\s*:(.*))?$", stmt)
        if not m:
            continue
        if m.group(1) is None:
            return True
        names = [x.split("=>")[-1].strip().lower() for x in m.group(2).split(",")]
        if "star" in names:
            return True
    return False


def star_access_counts(src_root):
    counts = {}
    for d in STAR_ACCESS_BUDGET:
        refs = files = 0
        for f in sorted((src_root / d).rglob("*.f90")):
            if "/test/" in f.relative_to(src_root).as_posix():
                continue
            code = strip_comments(f.read_text(errors="replace"))
            refs += len(re.findall(r"(?i)\bstar%", code))
            files += imports_global_star(code)
        counts[d] = (refs, files)
    return counts


def check_star_access_budget(src_root):
    problems, notes = [], []
    for d, (refs, files) in star_access_counts(src_root).items():
        max_refs, max_files = STAR_ACCESS_BUDGET[d]
        if refs > max_refs:
            problems.append(f"{d}: {refs} star% references, budget {max_refs} "
                            "(pass what a kernel needs as arguments instead)")
        if files > max_files:
            problems.append(f"{d}: {files} files import the global star, "
                            f"budget {max_files}")
        if (refs, files) < (max_refs, max_files) and refs <= max_refs and files <= max_files:
            notes.append(f"{d}: now ({refs}, {files}), budget ({max_refs}, {max_files}) "
                         "-- lower STAR_ACCESS_BUDGET to match")
    return problems, notes


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


# Directories whose code does no I/O and never stops the run (2026-10):
# their routines return a status code and their callers report
# (numerics/README.md, section 3).
NO_IO_DIRS = ["numerics", "math"]
IO_STATEMENT = re.compile(r"(?i)(?<![\w%])(?:(?:write|read|open|close|inquire)\s*\(|print\s*[*'0-9\x22])")
STOP_STATEMENT = re.compile(r"(?i)(?<![\w%])(?:error\s+)?stop(?![\w%])")

def check_no_io(src_root):
    problems = []
    for d in NO_IO_DIRS:
        for f in sorted((src_root / d).rglob("*.f90")):
            rel = f.relative_to(src_root).as_posix()
            if "/test/" in rel:
                continue
            for i, line in enumerate(f.read_text(errors="replace").splitlines(), 1):
                code = mask_strings(line)
                if IO_STATEMENT.search(code):
                    problems.append(f"{rel}:{i}: I/O statement in a no-I/O directory: {line.strip()[:60]}")
                if STOP_STATEMENT.search(code):
                    problems.append(f"{rel}:{i}: stop in a no-I/O directory (return a status): {line.strip()[:60]}")
    return problems


# Facades whose users must name what they take (docs/GUIDELINES.md,
# "Declarations and interfaces", rule 2: every `use` has an only: list;
# math_lib is the exception, since under USE_CRMATH its procedures
# shadow the intrinsics). Add a module once its users comply.
ONLY_LIST_MODULES = ["numerics_lib"]

def check_only_lists(src_root):
    problems = []
    use_re = re.compile(r"(?i)^\s*use\s+(\w+)\s*(,\s*only\s*:)?")
    for f in sorted(src_root.rglob("*.f90")):
        rel = f.relative_to(src_root).as_posix()
        for i, line in enumerate(strip_comments(f.read_text(errors="replace")).splitlines(), 1):
            m = use_re.match(line)
            if m and m.group(1).lower() in ONLY_LIST_MODULES and not m.group(2):
                problems.append(f"{rel}:{i}: `use {m.group(1)}` without an only: list")
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
    ref = _re.compile(r"(?i)(?<![a-z0-9_])" + math_funcs_pattern(src_root) + r"\s*\(")
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
# legacy_controls is the read's (read_legacy_deck reads a legacy pair
# into the buffer)
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
              "(read star%control / star%physics instead):")
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

    scratch_problems = check_scratch_private(SRC)
    if scratch_problems:
        print("Scratch used outside its domain (SCRATCH_OWNERS; pass the "
              "values as arguments or keep them on the star):")
        for b in scratch_problems:
            print(f"  {b}")
        return 1

    types_problems = check_scratch_types_only(SRC)
    if types_problems:
        print("Scratch module holds state (SCRATCH_TYPES_ONLY; the driver "
              "owns the instances and passes them down):")
        for b in types_problems:
            print(f"  {b}")
        return 1

    blind_problems = check_star_blind(SRC)
    if blind_problems:
        print("STAR-BLIND CONTRACT VIOLATIONS:")
        for b in blind_problems:
            print("  " + b)
        return 1

    g_problems = check_no_global_g(SRC)
    if g_problems:
        print("GLOBAL-G CONTRACT VIOLATIONS (physical_constants' cgl):")
        for b in g_problems:
            print("  " + b)
        return 1

    install_problems = check_install_sites(SRC)
    if install_problems:
        print("SETTINGS-INSTALL CONTRACT VIOLATIONS (INSTALL_FILES):")
        for b in install_problems:
            print("  " + b)
        return 1

    field_problems = check_settings_fields(SRC)
    if field_problems:
        print("SETTINGS-FIELD CONTRACT VIOLATIONS (setup/'s builders):")
        for b in field_problems:
            print("  " + b)
        return 1

    free_problems = check_star_free(SRC)
    if free_problems:
        print("STAR-FREE CONTRACT VIOLATIONS:")
        for b in free_problems:
            print("  " + b)
        return 1

    budget_problems, budget_notes = check_star_access_budget(SRC)
    if budget_problems:
        print("STAR-ACCESS BUDGET EXCEEDED (tools/check_boundaries.py):")
        for b in budget_problems:
            print("  " + b)
        return 1
    for n in budget_notes:
        print("star-access budget: " + n)

    controls_problems = check_controls_free(SRC)
    if controls_problems:
        print("CONTROLS-FREE CONTRACT VIOLATIONS:")
        for b in controls_problems:
            print("  " + b)
        return 1

    for check, header in ((check_old_records, "OLD CONTROLS RECORDS (ctrl and job are star%control, "
                           "star%physics and the run copies since 2026-10-30):"),
                          (check_controls_writer, "CONTROLS-RECORD WRITERS (only io/controls_sync.f90 "
                           "assigns star%control and star%physics):"),
                          (check_control_records, "CONTROLS-RECORD GROUPS (a control is in the record "
                           "of its registry group):"),
                          (check_control_copies, "CONTROL COPIES (a star member that copies a control "
                           "says so, where it is set):"),
                          (check_max_runs, "MAX_RUNS (state/control_slots.f90 is the one source):")):
        found = check(SRC)
        if found:
            print(header)
            for b in found:
                print("  " + b)
            return 1

    home_problems = check_control_homes(SRC)
    if home_problems:
        print("CONTROLS REGISTRY HOME/USE VIOLATIONS (defaults/controls_registry.tsv):")
        for b in home_problems:
            print("  " + b)
        return 1

    facade_problems = check_facade_imports(SRC)
    if facade_problems:
        print("FACADE CONTRACT VIOLATIONS:")
        for b in facade_problems:
            print("  " + b)
        return 1

    io_problems = check_no_io(SRC)
    if io_problems:
        print("NO-I/O CONTRACT VIOLATIONS (numerics/, math/):")
        for m in io_problems:
            print("  " + m)
        return 1

    only_problems = check_only_lists(SRC)
    if only_problems:
        print("ONLY-LIST CONTRACT VIOLATIONS (use <facade>, only: ...):")
        for m in only_problems:
            print("  " + m)
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

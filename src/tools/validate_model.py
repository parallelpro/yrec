#!/usr/bin/env python3
"""Read YREC's model files and re-derive the structure from primitive columns.

The physics validation tier's shared library (docs/PHYSICS_VALIDATION.md)
and a command-line audit of one model:

    python3 -I src/tools/validate_model.py output/profile6.data.GYRE \\
        [output/profile6.data.FGONG] [output/profile6.data.GSM]

It needs only the standard library (CI installs pytest and nothing else).

Readers
    read_mesa_table   history.data and profile<N>.data, YREC's or MESA's
    read_gyre_text    the GYRE text model (schema 101, 18 columns)
    read_fgong        FGONG (ivers 1300)
    read_gsm          GSM (HDF5): with h5py, else through the h5dump program
    read_gyre_summary a GYRE summary file written with summary_file_format
                      = 'TXT'

Re-derivations, each from primitive columns only (r, m, P, T, rho,
Gamma_1, L, eps, the gradients), never from the column under test
    hydrostatic_residual   d ln P/dr against -G m rho/(r^2 P), point by point
    hydrostatic_defect     the same law integrated between two points:
                           ln P(b) - ln P(a) + int G m rho/(r^2 P) dr
    continuity_residual    dm/dr against 4 pi r^2 rho
    luminosity_residual    dL/dm against eps
    brunt_n2_density       g [(1/Gamma_1) dlnP/dr - dlnrho/dr]
    brunt_n2_thermal       g^2 rho delta (grad_ad - grad)/P
    lamb_s2, gravity, sound_speed, acoustic_delta_nu, period_spacing

A model is a dict of equal-length lists ordered from the centre outward,
whatever order its file uses.
"""
import math
import os
import pathlib
import shutil
import subprocess
import sys

# a sibling module: `python3 -I` leaves the script's directory off sys.path
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import harness  # noqa: E402

# YREC's default G: log10 G = -7.17571 (state/physical_constants.f90). The
# GYRE text and GSM files carry no G, so a pulsation code must be told
# this value (GYRE: &constants G_GRAVITY); FGONG carries it as global 15.
YREC_G = 10.0 ** -7.17571
FOUR_PI = 4.0 * math.pi

GYRE_COLUMNS = ["r", "M_r", "L_r", "P", "T", "rho", "nabla", "N2", "Gamma_1",
                "nabla_ad", "delta", "kap", "kap_kap_T", "kap_kap_rho", "eps",
                "eps_eps_T", "eps_eps_rho", "Omega_rot"]


# ---------------------------------------------------------------------
# readers
# ---------------------------------------------------------------------
def _number(token):
    token = token.strip('"')
    try:
        return int(token)
    except ValueError:
        pass
    try:
        return float(token)
    except ValueError:
        # Fortran drops the E of a three-digit exponent: 1.0-100
        for sign in ("-", "+"):
            head, sep, tail = token[1:].rpartition(sign)
            if sep and head and tail.isdigit():
                try:
                    return float(token[0] + head + "E" + sign + tail)
                except ValueError:
                    pass
        return token


def read_mesa_table(path):
    """A MESA-layout table: (header dict, column dict of lists in file order)."""
    lines = pathlib.Path(path).read_text().splitlines()
    header = dict(zip(lines[1].split(), (_number(t) for t in lines[2].split())))
    names = lines[5].split()
    cols = {n: [] for n in names}
    for line in lines[6:]:
        tokens = line.split()
        if len(tokens) != len(names):
            continue
        for n, t in zip(names, tokens):
            cols[n].append(_number(t))
    return header, cols


def read_profile(path):
    """A profile as (header, columns), the columns' lists reversed to run
    from the centre outward (the file's zone 1 is the surface), in the
    file's own units (profile_model converts them)."""
    header, cols = read_mesa_table(path)
    model = {n: v[::-1] for n, v in cols.items()}
    return header, model


def read_gyre_text(path):
    """The GYRE text model: (header dict n, M_star, R_star, L_star, version
    and rows, the number of data rows, each checked to carry its own
    number; model dict of the 18 columns, centre outward)."""
    lines = pathlib.Path(path).read_text().splitlines()
    h = lines[0].split()
    header = {"n": int(h[0]), "M_star": float(h[1]), "R_star": float(h[2]),
              "L_star": float(h[3]), "version": int(h[4])}
    model = {n: [] for n in GYRE_COLUMNS}
    rows = 0
    for line in lines[1:]:
        tokens = line.split()
        if not tokens:
            continue
        rows += 1
        if int(tokens[0]) != rows:
            raise ValueError(f"{path}: row {rows} is numbered {tokens[0]}")
        for n, t in zip(GYRE_COLUMNS, tokens[1:]):
            model[n].append(_number(t))
    header["rows"] = rows
    return header, model


def read_fgong(path):
    """FGONG: (header dict with text, nn, iconst, ivar, ivers, glob; var,
    a list of ivar lists in FILE order, surface to centre)."""
    lines = pathlib.Path(path).read_text().splitlines()
    nn, iconst, ivar, ivers = (int(x) for x in lines[4].split())
    width = 27 if ivers >= 1000 else 16
    values = []
    for line in lines[5:]:
        line = line.rstrip("\n")
        for i in range(0, len(line), width):
            field = line[i:i + width].strip()
            if field:
                values.append(_number(field))
    glob = values[:iconst]
    body = values[iconst:]
    var = [body[k::ivar][:nn] for k in range(ivar)]
    header = {"text": lines[:4], "nn": nn, "iconst": iconst, "ivar": ivar,
              "ivers": ivers, "glob": glob, "n_values": len(values)}
    return header, var


def fgong_model(header, var):
    """The FGONG columns as a centre-outward model in the GYRE column
    names (those FGONG holds), with N2 rebuilt from A4 = r N^2/g."""
    g_newton = header["glob"][14]
    m_star = header["glob"][0]
    n = header["nn"]
    rev = [v[::-1] for v in var]
    model = {
        "r": rev[0],
        "M_r": [m_star * math.exp(x) if x > -1.0e90 else 0.0 for x in rev[1]],
        "T": rev[2], "P": rev[3], "rho": rev[4], "X": rev[5], "L_r": rev[6],
        "kap": rev[7], "eps": rev[8], "Gamma_1": rev[9], "nabla_ad": rev[10],
        "delta": rev[11], "cp": rev[12], "mu_e_inv": rev[13], "A4": rev[14],
        "Z": rev[16], "R_minus_r": rev[17], "eps_grav": rev[18],
    }
    n2 = []
    for i in range(n):
        r, m = model["r"][i], model["M_r"][i]
        n2.append(model["A4"][i] * g_newton * m / r ** 3 if r > 0.0 else 0.0)
    model["N2"] = n2
    return model


def _h5dump_program():
    for candidate in (shutil.which("h5dump"),
                      os.path.join(os.environ.get("MESASDK_ROOT", ""), "bin", "h5dump")):
        if candidate and os.path.exists(candidate):
            return candidate
    return None


def gsm_readable():
    """Can read_gsm work here (h5py, or an h5dump program)?"""
    try:
        import h5py  # noqa: F401
        return True
    except ImportError:
        return _h5dump_program() is not None


def read_gsm(path):
    """GSM: (header dict of the root attributes; model dict of datasets,
    centre outward)."""
    try:
        import h5py
    except ImportError:
        h5py = None
    if h5py is not None:
        with h5py.File(path, "r") as f:
            header = {}
            for k, v in f.attrs.items():
                header[k] = v[0].item() if hasattr(v, "__len__") else v
            model = {k: [float(x) for x in f[k][...]] for k in f.keys()}
        return header, model
    prog = _h5dump_program()
    if prog is None:
        raise RuntimeError("read_gsm needs h5py or an h5dump program")
    text = subprocess.run([prog, "-m", "%.17g", "-w", "0", str(path)],
                          capture_output=True, text=True, check=True).stdout
    header, model = {}, {}
    name, kind = None, None
    in_data = False
    for line in text.splitlines():
        s = line.strip()
        if s.startswith("ATTRIBUTE "):
            name, kind = s.split('"')[1], "attr"
        elif s.startswith("DATASET "):
            name, kind = s.split('"')[1], "data"
            model[name] = []
        elif s.startswith("DATA {"):
            in_data = True
        elif in_data and s.startswith("}"):
            in_data = False
        elif in_data and s.startswith("("):
            values = [_number(t) for t in s.split(":", 1)[1].replace(",", " ").split()]
            if kind == "attr":
                header[name] = values[0]
            else:
                model[name].extend(values)
    return header, model


def read_gyre_summary(path):
    """A GYRE summary file in its TXT format: (header dict, column dict)."""
    lines = pathlib.Path(path).read_text().splitlines()
    header = dict(zip(lines[2].split(), (_number(t) for t in lines[3].split())))
    names = lines[5].split()
    cols = {n: [] for n in names}
    for line in lines[6:]:
        tokens = line.split()
        if len(tokens) != len(names):
            continue
        for n, t in zip(names, tokens):
            cols[n].append(_number(t))
    return header, cols


def profile_model(header, prof, msun_g, lsun_cgs):
    """A profile's columns as a centre-outward model in the GYRE column
    names and cgs units. msun_g and lsun_cgs convert the profile's mass
    (Msun) and luminosity (Lsun) columns; they are the run's own
    constants, which the caller takes from the run (a pulse file's M_star
    over the history's star_mass), not from this module."""
    model = {
        "r": [10.0 ** x for x in prof["logR"]],
        "M_r": [x * msun_g for x in prof["mass"]],
        "P": [10.0 ** x for x in prof["logP"]],
        "T": [10.0 ** x for x in prof["logT"]],
        "rho": [10.0 ** x for x in prof["logRho"]],
        "L_r": [x * lsun_cgs for x in prof["luminosity"]],
    }
    for gyre_name, prof_name in (("Gamma_1", "gamma1"), ("nabla", "gradT"),
                                 ("nabla_ad", "grada"), ("delta", "delta"),
                                 ("kap", "opacity"), ("eps", "eps_nuc"),
                                 ("N2", "brunt_N2"), ("Omega_rot", "omega")):
        if prof_name in prof:
            model[gyre_name] = list(prof[prof_name])
    return model


# ---------------------------------------------------------------------
# re-derivations
# ---------------------------------------------------------------------
def centered(y, x):
    """dy/dx by the centred difference (y[i+1]-y[i-1])/(x[i+1]-x[i-1]);
    None at the two ends and where x does not advance."""
    n = len(y)
    out = [None] * n
    for i in range(1, n - 1):
        dx = x[i + 1] - x[i - 1]
        if dx != 0.0:
            out[i] = (y[i + 1] - y[i - 1]) / dx
    return out


def gravity(model, g_newton=YREC_G):
    return [g_newton * m / (r * r) if r > 0.0 else 0.0
            for m, r in zip(model["M_r"], model["r"])]


def sound_speed(model):
    return [math.sqrt(g1 * p / rho)
            for g1, p, rho in zip(model["Gamma_1"], model["P"], model["rho"])]


def lamb_s2(model, ell=1):
    """The Lamb frequency squared, l(l+1) c^2/r^2."""
    return [ell * (ell + 1) * g1 * p / (rho * r * r) if r > 0.0 else 0.0
            for g1, p, rho, r in zip(model["Gamma_1"], model["P"], model["rho"], model["r"])]


def hydrostatic_residual(model, g_newton=YREC_G):
    """(d ln P/dr) / (-G m rho/(r^2 P)) - 1 at interior points of the
    grid, by centred differences; None at the ends and the centre."""
    lnp = [math.log(p) for p in model["P"]]
    dlnp = centered(lnp, model["r"])
    out = [None] * len(lnp)
    for i, d in enumerate(dlnp):
        r = model["r"][i]
        if d is None or r <= 0.0:
            continue
        rhs = -g_newton * model["M_r"][i] * model["rho"][i] / (r * r * model["P"][i])
        if rhs != 0.0:
            out[i] = d / rhs - 1.0
    return out


def hydrostatic_defect(model, a, b, g_newton=YREC_G):
    """ln P(b) - ln P(a) + int_a^b G m rho/(r^2 P) dr by the trapezoidal
    rule over the model's own points a..b: zero for a hydrostatic model,
    up to the rule's error. Also returns the integral, for scale."""
    total = 0.0
    prev = None
    for i in range(a, b + 1):
        r = model["r"][i]
        f = g_newton * model["M_r"][i] * model["rho"][i] / (r * r * model["P"][i]) if r > 0.0 else 0.0
        if prev is not None:
            total += 0.5 * (f + prev[1]) * (r - prev[0])
        prev = (r, f)
    return math.log(model["P"][b]) - math.log(model["P"][a]) + total, total


def continuity_residual(model):
    """(dm/dr) / (4 pi r^2 rho) - 1 by centred differences (None where m
    does not change, as in the massless atmosphere)."""
    dm = centered(model["M_r"], model["r"])
    out = [None] * len(dm)
    for i, d in enumerate(dm):
        r = model["r"][i]
        if d is None or r <= 0.0 or d == 0.0:
            continue
        out[i] = d / (FOUR_PI * r * r * model["rho"][i]) - 1.0
    return out


def continuity_defect(model, a, b):
    """m(b) - m(a) - int_a^b 4 pi r^2 rho dr (trapezoidal in r^3, exact
    for constant density), relative to m(b) - m(a)."""
    total = 0.0
    for i in range(a + 1, b + 1):
        r0, r1 = model["r"][i - 1], model["r"][i]
        total += 0.5 * (model["rho"][i - 1] + model["rho"][i]) * FOUR_PI / 3.0 * (r1 ** 3 - r0 ** 3)
    span = model["M_r"][b] - model["M_r"][a]
    return (span - total) / span


def luminosity_defect(model, eps_total, a, b):
    """L(b) - L(a) - int_a^b eps dm (trapezoidal in m), relative to
    L(b) - L(a); eps_total is the whole source per gram (nuclear minus
    neutrino plus gravitational)."""
    total = 0.0
    for i in range(a + 1, b + 1):
        total += 0.5 * (eps_total[i - 1] + eps_total[i]) * (model["M_r"][i] - model["M_r"][i - 1])
    span = model["L_r"][b] - model["L_r"][a]
    return (span - total) / span


def structure_terms(model, g_newton=YREC_G):
    """The right-hand sides of the structure equations at every point,
    from r, m, P and rho only, and the logs they are differenced against:
      dlnP_dlnm = -G m^2/(4 pi r^4 P)      hydrostatic, mass coordinate
      dlnr_dlnm = m/(4 pi r^3 rho)         continuity, mass coordinate
      dr_dlnP   = -P r^2/(G m rho)         hydrostatic, pressure coordinate
      dm_dlnP   = -4 pi r^4 P/(G m)        continuity over hydrostatic
    (the mass coordinate suits the interior, where YREC's own difference
    scheme is the trapezoidal rule in log m; the pressure coordinate the
    envelope and atmosphere, which are integrated in log P and hold
    almost no mass)."""
    r, m, p, rho = model["r"], model["M_r"], model["P"], model["rho"]
    n = len(r)
    return {
        "lnm": [math.log(x) if x > 0.0 else None for x in m],
        "lnP": [math.log(x) for x in p],
        "lnr": [math.log(x) if x > 0.0 else None for x in r],
        "lnT": [math.log(x) for x in model["T"]],
        "dlnP_dlnm": [-g_newton * m[k] ** 2 / (FOUR_PI * r[k] ** 4 * p[k]) for k in range(n)],
        "dlnr_dlnm": [m[k] / (FOUR_PI * r[k] ** 3 * rho[k]) for k in range(n)],
        "dr_dlnP": [-p[k] * r[k] ** 2 / (g_newton * m[k] * rho[k]) for k in range(n)],
        "dm_dlnP": [-FOUR_PI * r[k] ** 4 * p[k] / (g_newton * m[k]) for k in range(n)],
    }


def quad3(x0, x1, x2, f0, f1, f2):
    """The integral over [x0, x2] of the parabola through three points
    (Simpson's rule on an uneven mesh)."""
    h0, h1 = x1 - x0, x2 - x1
    h = h0 + h1
    return h / 6.0 * ((2.0 - h1 / h0) * f0 + h * h / (h0 * h1) * f1 + (2.0 - h0 / h1) * f2)


def pair_defects(y, x, f, lo, hi, scale=None):
    """For each pair of neighbours k, k+1 in lo..hi: y[k+1] - y[k] against
    the trapezoidal integral of f = dy/dx, as (difference - integral)
    over the integral, or over scale when one is given. None where x does
    not advance."""
    out = []
    for k in range(lo, hi):
        dx = x[k + 1] - x[k]
        if dx == 0.0:
            out.append(None)
            continue
        q = 0.5 * (f[k] + f[k + 1]) * dx
        d = y[k + 1] - y[k]
        out.append((d - q) / (scale if scale is not None else q))
    return out


def triple_defects(y, x, f, lo, hi, scale=None):
    """The same over each three neighbours k, k+1, k+2 in lo..hi with the
    three-point rule quad3: the law at one order higher than the pair
    form, so it shows the difference scheme's own truncation error where
    the pair form shows only how well the scheme's equations are solved."""
    out = []
    for k in range(lo, hi - 1):
        if x[k + 1] == x[k] or x[k + 2] == x[k + 1]:
            out.append(None)
            continue
        q = quad3(x[k], x[k + 1], x[k + 2], f[k], f[k + 1], f[k + 2])
        d = y[k + 2] - y[k]
        out.append((d - q) / (scale if scale is not None else q))
    return out


def photosphere_index(model, r_star):
    """The last point at or below r_star (the pulse header's R_star)."""
    return max(k for k, r in enumerate(model["r"]) if r <= r_star * (1.0 + 1.0e-12))


def brunt_n2_file_rule(model, g_newton=YREC_G):
    """N^2 by the rule YREC documents for its files: the thermal form
    where that is not positive (convective), else the centred density
    form; each end copies its neighbour."""
    thermal = brunt_n2_thermal(model, g_newton)
    density = brunt_n2_density(model, g_newton)
    n = len(thermal)
    out = list(thermal)
    for j in range(1, n - 1):
        if thermal[j] > 0.0 and density[j] is not None and model["r"][j] > 0.0:
            out[j] = density[j]
    if n >= 2:
        out[0], out[-1] = out[1], out[-2]
    return out


def brunt_brackets(model):
    """The two dimensionless forms of N^2 P/(g^2 rho): the structure's
    dlnrho/dlnP - 1/Gamma_1 (centred differences) and the thermal
    delta (grad_ad - grad). They are equal where the composition is
    uniform; their difference elsewhere is the composition term."""
    dlnrho = centered([math.log(x) for x in model["rho"]], [math.log(x) for x in model["P"]])
    structure = [None if d is None else d - 1.0 / g1 for d, g1 in zip(dlnrho, model["Gamma_1"])]
    thermal = [d * (na - nb) for d, na, nb in zip(model["delta"], model["nabla_ad"], model["nabla"])]
    return structure, thermal


def brunt_n2_density(model, g_newton=YREC_G):
    """N^2 = g [(1/Gamma_1) dlnP/dr - dlnrho/dr], centred differences."""
    g = gravity(model, g_newton)
    dlnp = centered([math.log(p) for p in model["P"]], model["r"])
    dlnrho = centered([math.log(x) for x in model["rho"]], model["r"])
    return [None if a is None or b is None else gi * (a / g1 - b)
            for gi, a, b, g1 in zip(g, dlnp, dlnrho, model["Gamma_1"])]


def brunt_n2_thermal(model, g_newton=YREC_G):
    """N^2 = g^2 rho delta (grad_ad - grad)/P: the thermal part, the whole
    of N^2 where the composition is uniform."""
    g = gravity(model, g_newton)
    return [gi * gi * rho * d * (na - nb) / p
            for gi, rho, d, na, nb, p in zip(g, model["rho"], model["delta"],
                                             model["nabla_ad"], model["nabla"], model["P"])]


def acoustic_delta_nu(model, r_outer):
    """[2 int_0^R dr/c]^-1 in uHz, trapezoidal over the points at or
    below r_outer (the photosphere)."""
    c = sound_speed(model)
    total = 0.0
    for i in range(1, len(c)):
        if model["r"][i] > r_outer * (1.0 + 1.0e-12):
            break
        total += 0.5 * (1.0 / c[i] + 1.0 / c[i - 1]) * (model["r"][i] - model["r"][i - 1])
    return 1.0e6 / (2.0 * total)


def period_spacing(model, n2, ell=1):
    """The asymptotic g-mode period spacing 2 pi^2/sqrt(l(l+1)) /
    int N dr/r over the innermost cavity (N^2 > 0, from the centre
    outward to the first point where it is not), trapezoidal in ln r. In
    seconds; 0 without a cavity."""
    total = 0.0
    inside = False
    prev = None
    for i in range(1, len(n2)):
        r = model["r"][i]
        if n2[i] is not None and n2[i] > 0.0 and r > 0.0:
            inside = True
            point = (math.log(r), math.sqrt(n2[i]))
            if prev is not None:
                total += 0.5 * (point[1] + prev[1]) * (point[0] - prev[0])
            prev = point
        elif inside:
            break
    return 2.0 * math.pi ** 2 / math.sqrt(ell * (ell + 1)) / total if total > 0.0 else 0.0


# ---------------------------------------------------------------------
# running YREC and GYRE
# ---------------------------------------------------------------------
REPO = pathlib.Path(__file__).resolve().parents[2]


def edit_inlist(text, group, lines):
    """inlist text with the lines added at the head of &control or
    &physics. A name the deck also sets must be replaced in the deck's
    text instead (the deck's later line would win)."""
    if group not in text:
        raise ValueError(f"no {group} group")
    return text.replace(group, group + "\n" + "".join(f" {line}\n" for line in lines), 1)


def all_columns_file(list_name):
    """The text of a columns file with every column of
    src/defaults/<list_name> switched on ('!name' marks an opt-in one)."""
    import re
    text = (REPO / "src" / "defaults" / list_name).read_text()
    return re.sub(r"^!(?=[A-Za-z])", "", text, flags=re.M)


def run_yrec(yrec, workdir, inlist_text, timeout=3600):
    """Run yrec on inlist_text in workdir (made here, with output/ and the
    two all-columns files profile_all.list and history_all.list), with
    this repository's input/ and startmodels/; returns the output
    directory. Raises RuntimeError when the run fails."""
    workdir = pathlib.Path(workdir)
    (workdir / "output").mkdir(parents=True)
    (workdir / "inlist").write_text(inlist_text)
    (workdir / "profile_all.list").write_text(all_columns_file("profile_columns.list"))
    (workdir / "history_all.list").write_text(all_columns_file("history_columns.list"))
    # this repository's tables and starting models, whatever YREC_INPUT and
    # YREC_START name in the caller's environment (harness.py has the
    # rule and its one opt-in)
    env = dict(os.environ)
    harness.apply(REPO, env)
    env["YREC_OUTPUT"] = "output"
    r = subprocess.run([str(yrec), "inlist"], cwd=workdir, env=env,
                       capture_output=True, text=True, timeout=timeout)
    if r.returncode != 0:
        raise RuntimeError(f"yrec failed in {workdir}:\n" + r.stdout[-800:] + r.stderr[-800:])
    return workdir / "output"


def run_with_final_profile(yrec, workdir, inlist_text, pulse_formats=("GYRE",), timeout=3600):
    """Run a deck to its end and write the profile and pulse files of its
    LAST model. YREC writes profiles every profile_interval models, so the
    deck is run once to learn the last model's number N (the run is
    deterministic) and once per pulse format with profile_interval = N.
    Returns {"history": dir of the first run, <format>: output dir}; the
    last model's files are profile1.data and profile1.data.<format>."""
    workdir = pathlib.Path(workdir)
    text = edit_inlist(inlist_text, "&control", ["profile_columns_file = 'profile_all.list'",
                                                 "history_columns_file = 'history_all.list'"])
    out = {"history": run_yrec(yrec, workdir / "count", text, timeout)}
    _, hist = read_mesa_table(out["history"] / "history.data")
    last = hist["model_number"][-1]
    for fmt in pulse_formats:
        with_files = edit_inlist(text, "&physics", [
            f"profile_interval = {last}", "write_profile_flag = .true.",
            "write_pulse_flag = .true.", f"pulse_format = '{fmt}'"])
        out[fmt] = run_yrec(yrec, workdir / fmt.lower(), with_files, timeout)
        if not (out[fmt] / "profile1.data").exists():
            raise RuntimeError(f"no profile at the last model {last} in {out[fmt]}")
    return out


GYRE_INLIST = """&constants
  G_GRAVITY = {g_newton:.16e}
/
&model
  model_type = 'EVOL'
  file = '{model_file}'
  file_format = '{file_format}'
{model_extra}/
{modes}
&osc
  outer_bound = '{outer_bound}'
/
&rot
/
&num
  diff_scheme = 'COLLOC_GL4'
/
{scans}
&grid
  w_osc = 10
  w_exp = 2
  w_ctr = 10
/
&ad_output
  summary_file = '{summary}'
  summary_file_format = 'TXT'
  summary_item_list = 'l,n_pg,n_p,n_g,freq,E_norm'
  freq_units = 'UHZ'
/
&nad_output
/
"""

# GYRE's name for each file YREC writes, and what its reader needs besides
GYRE_FILE_FORMATS = {"GYRE": ("MESA", ""), "GSM": ("GSM", ""),
                     "FGONG": ("FGONG", "  data_format = '(1P,5(X,E26.18E3))'\n")}


def gyre_program():
    """The gyre executable under $GYRE_DIR, or None."""
    root = os.environ.get("GYRE_DIR")
    if not root:
        return None
    prog = pathlib.Path(root) / "bin" / "gyre"
    return prog if prog.exists() else None


def run_gyre(model_file, pulse_format, workdir, scans, g_newton=YREC_G,
             outer_bound="VACUUM", tag="modes", timeout=1800):
    """Adiabatic modes of a YREC pulse file with GYRE. scans is a list of
    (l, grid_type, freq_min_uHz, freq_max_uHz, n_freq), grid_type 'LINEAR'
    for p modes or 'INVERSE' for g modes. G_GRAVITY is set to g_newton:
    the GYRE text and GSM files carry no G, and GYRE's own (6.67430e-8)
    is 2.7e-4 from YREC's. Returns the summary's columns."""
    prog = gyre_program()
    if prog is None:
        raise RuntimeError("GYRE_DIR does not name a GYRE installation")
    workdir = pathlib.Path(workdir)
    workdir.mkdir(parents=True, exist_ok=True)
    file_format, model_extra = GYRE_FILE_FORMATS[pulse_format]
    modes = "".join(f"&mode\n  l = {ell}\n  tag = 'l{ell}_{i}'\n/\n"
                    for i, (ell, *_rest) in enumerate(scans))
    scan_text = "".join(
        f"&scan\n  grid_type = '{grid}'\n  freq_min = {fmin}\n  freq_max = {fmax}\n"
        f"  n_freq = {nfreq}\n  freq_min_units = 'UHZ'\n  freq_max_units = 'UHZ'\n"
        f"  tag_list = 'l{ell}_{i}'\n/\n"
        for i, (ell, grid, fmin, fmax, nfreq) in enumerate(scans))
    summary = f"{tag}_summary.txt"
    (workdir / f"{tag}.in").write_text(GYRE_INLIST.format(
        g_newton=g_newton, model_file=model_file, file_format=file_format,
        model_extra=model_extra, modes=modes, outer_bound=outer_bound,
        scans=scan_text, summary=summary))
    env = dict(os.environ)
    env.setdefault("OMP_NUM_THREADS", "4")
    r = subprocess.run([str(prog), f"{tag}.in"], cwd=workdir, env=env,
                       capture_output=True, text=True, timeout=timeout)
    (workdir / f"{tag}.log").write_text(r.stdout + r.stderr)
    if r.returncode != 0 or not (workdir / summary).exists():
        raise RuntimeError(f"gyre failed in {workdir} ({tag}):\n" + r.stdout[-1500:] + r.stderr[-800:])
    _, cols = read_gyre_summary(workdir / summary)
    return cols


def mode_table(cols, ell):
    """{n_pg: frequency in uHz} of the degree-l modes of a GYRE summary."""
    return {n: f for l_, n, f in zip(cols["l"], cols["n_pg"], cols["Re(freq)"]) if l_ == ell}


def large_separation(freqs_by_order, n_lo, n_hi):
    """The mean large separation over the radial orders n_lo..n_hi: the
    slope of a least-squares line through (n, nu). Also the orders used."""
    pts = [(n, f) for n, f in sorted(freqs_by_order.items()) if n_lo <= n <= n_hi]
    count = len(pts)
    mean_n = sum(n for n, _ in pts) / count
    mean_f = sum(f for _, f in pts) / count
    slope = sum((n - mean_n) * (f - mean_f) for n, f in pts) / sum((n - mean_n) ** 2 for n, _ in pts)
    return slope, count


def mean_period_spacing(freqs_by_order, n_lo, n_hi):
    """The mean period spacing (s) of the g modes with n_lo <= |n_pg| <=
    n_hi (n_pg is negative for g modes): the slope of a least-squares line
    through (|n|, period). Also the orders used."""
    pts = [(-n, 1.0e6 / f) for n, f in sorted(freqs_by_order.items()) if n_lo <= -n <= n_hi]
    count = len(pts)
    mean_n = sum(n for n, _ in pts) / count
    mean_p = sum(per for _, per in pts) / count
    slope = sum((n - mean_n) * (per - mean_p) for n, per in pts) / sum((n - mean_n) ** 2 for n, _ in pts)
    return slope, count


def is_strictly_increasing(values):
    return all(b > a for a, b in zip(values, values[1:]))


def worst(values, lo=None, hi=None):
    """(largest |value|, its index) over values[lo:hi], skipping None."""
    best, where = 0.0, None
    for i, v in enumerate(values):
        if v is None or (lo is not None and i < lo) or (hi is not None and i >= hi):
            continue
        if abs(v) >= best:
            best, where = abs(v), i
    return best, where


def rel_diff(a, b):
    scale = max(abs(a), abs(b))
    return abs(a - b) / scale if scale > 0.0 else 0.0


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        return 2
    for path in argv[1:]:
        name = pathlib.Path(path).name
        if name.endswith(".FGONG"):
            header, var = read_fgong(path)
            model = fgong_model(header, var)
            g_newton = header["glob"][14]
            r_star = header["glob"][1]
        elif name.endswith(".GYRE"):
            header, model = read_gyre_text(path)
            g_newton, r_star = YREC_G, header["R_star"]
        elif name.endswith(".GSM"):
            header, model = read_gsm(path)
            g_newton, r_star = YREC_G, header["R_star"]
        else:
            print(f"{name}: give a pulse file (.GYRE, .FGONG or .GSM); a profile needs "
                  "the run's solar constants (see test_validation_contract.py)")
            continue
        n = len(model["r"])
        iph = photosphere_index(model, r_star)
        terms = structure_terms(model, g_newton)
        print(f"{name}: {n} points ({n - 1 - iph} above the photosphere), G = {g_newton:.6e}")
        print(f"  r strictly increasing: {is_strictly_increasing(model['r'])}")
        # below the photosphere: the three-point rule in ln P (in the
        # interior it shows the difference scheme's truncation, a few
        # 1e-4; in the envelope the output mesh's, about 1e-3; at the
        # fitting point the surface boundary condition's mismatch)
        # (a pair of points closer than 2e-5 in ln r is the envelope's
        # repeat of the fitting point, V2 of docs/PHYSICS_VALIDATION.md:
        # the rule means nothing across it)
        close = [k for k in range(iph) if terms["lnr"][k + 1] - terms["lnr"][k] < 2.0e-5]
        for k in close:
            print(f"  a near-duplicate pair at points {k}, {k + 1} (r/R = {model['r'][k] / r_star:.5f}): "
                  f"ln rho steps by {math.log(model['rho'][k + 1] / model['rho'][k]):+.1e} outward")
        for label, y, f in (("hydrostatic equilibrium", model["r"], terms["dr_dlnP"]),
                            ("mass continuity", model["M_r"], terms["dm_dlnP"])):
            defects = triple_defects(y, terms["lnP"], f, 0, iph)
            for k in close:
                for j in (k - 1, k):
                    if 0 <= j < len(defects):
                        defects[j] = None
            v, i = worst(defects)
            ordered = sorted(abs(d) for d in defects if d is not None)
            print(f"  {label}, three-point rule in ln P, below the photosphere: median "
                  f"{ordered[len(ordered) // 2]:.1e}, worst {v:.1e} at point {i} "
                  f"(r/R = {model['r'][i] / r_star:.5f})")
        if n - 1 > iph + 1:
            v, i = worst(pair_defects(model["r"], terms["lnP"], terms["dr_dlnP"], iph + 1, n - 1))
            print(f"  hydrostatic equilibrium in the atmosphere, pair rule: worst {v:.1e}")
        if "nabla" in model:
            rule = brunt_n2_file_rule(model, g_newton)
            scale = [g * g * rho / pr for g, rho, pr in zip(gravity(model, g_newton), model["rho"], model["P"])]
            v, i = worst([(a - b) / s for a, b, s in zip(rule, model["N2"], scale)])
            print(f"  N2 against its rule (thermal form where not positive, else the density form), "
                  f"in g^2 rho/P: worst {v:.1e}")
        print(f"  [2 int dr/c]^-1 to R_star: {acoustic_delta_nu(model, r_star):.4f} uHz; "
              f"to the top: {acoustic_delta_nu(model, model['r'][-1]):.4f} uHz")
        print(f"  sqrt(2) pi^2/int N dr/r from the innermost point, l = 1: "
              f"{period_spacing(model, model['N2']):.2f} s (the centre's part is left out)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))

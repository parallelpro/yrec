"""Structure and pulse-file contract (physics validation tier, piece 1).

What a produced model must satisfy, checked from the files alone: the
structure equations re-derived in Python from primitive columns (r, m,
P, T, rho, L, eps, Gamma_1, the gradients), never from the column under
test. Nothing here is a pin of YREC against itself; every expected
value is a law (docs/PHYSICS_VALIDATION.md section 1 has the measured
values and what sets each tolerance).

The model is the standard solar deck (examples/run_standard_solar_model,
Test_solar_noGS_norot), run twice:

  gyre   to the solar age with every profile column, the GYRE text file
         of each profile (pulse_format = 'GYRE') and the interior-only
         GYRE file every 100 models: six models, 2 Myr to 3.4 Gyr
  fgong  to model 100 with pulse_format = 'FGONG': the same two first
         models in FGONG
  gsm    as fgong with pulse_format = 'GSM', only with YREC_TEST_HDF5=1
         (an HDF5 build, as test_mesa_output.py's GSM test)

test_mesa_output.py checks the output's layout (the file set, the
headers' shape, FGONG's record count); this file checks its physics.

Each check is a function of one model that raises AssertionError.
test_each_check_can_fail applies a named mutation to a copy of the model
(most are bugs this code once had: docs/bug_sweep.md) and requires the
check to fail, so a check that cannot fail is itself a test failure.

Known defects are expected failures (xfail, strict), each named: the
check states the law at the tolerance a correct code would meet, and
turns into a failure of this file (XPASS) when the defect is fixed, at
which point the xfail mark is removed.

YREC_EXE names another yrec to validate (default src/yrec). Run time:
about 35 s. Needs only the tables CI's sparse checkout has.
"""
import copy
import math
import os
import pathlib
import re
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).parent
sys.path.insert(0, str(REPO / "src" / "tools"))
import validate_model as vm  # noqa: E402

YREC = pathlib.Path(os.environ.get("YREC_EXE", REPO / "src" / "yrec"))
CONVERTER = REPO / "src" / "tools" / "upgrade_inlist.py"
CASE = REPO / "examples" / "run_standard_solar_model"
DECK = "Test_solar_noGS_norot"
DEFAULTS = REPO / "src" / "defaults"

# YREC's Stefan-Boltzmann constant (state/physical_constants.f90, CODATA 1986)
SIGMA_SB = 5.67051e-5

# ---------------------------------------------------------------------
# Tolerances. Each is about three times the largest value measured over
# the six models (2026-11-17, macOS/arm64, gfortran 15.2; the measured
# value is in the comment), and each names what sets it.
# ---------------------------------------------------------------------
TOL = {
    # header numbers against the history row and the laws that tie them:
    # round-off of exp10/log10 chains
    "header": 1.0e-12,
    # interior, YREC's own scheme (trapezoid in log m): the Henyey
    # solver's convergence (solver_tolerances 3e-7 to 9e-7 on this deck)
    "interior_hydrostatic_pair": 6.0e-6,      # measured 1.3e-6
    "interior_continuity_pair": 6.0e-6,       # measured 2.0e-6
    # interior at one order higher: the scheme's truncation error, h^2/12
    # with h = the shell's step in ln P (0.046 from tol_dp_env_max = 0.02)
    "interior_hydrostatic_triple": 1.1e-3,    # measured 3.5e-4
    "interior_continuity_triple": 8.0e-4,     # measured 2.6e-4
    # dlnT/dlnP = nabla in the interior: the scheme centres nabla on the
    # pair, the file has it at the points; largest at a convective boundary
    "interior_gradient_pair": 1.3e-3,         # measured 4.2e-4
    # envelope (integrated in log P to 1e-9, written every 0.05 dex): the
    # three-point rule's error on that output mesh
    "envelope_hydrostatic_triple": 1.8e-3,    # measured 5.7e-4
    "envelope_continuity_triple": 1.2e-3,     # measured 3.8e-4
    "envelope_gradient_triple": 2.0e-2,       # measured 6.2e-3 (H ionization)
    # atmosphere heights: plane-parallel gravity (2 (r - R)/R = 2.5e-3 at
    # the top) and the depth quadrature dtau/<rho kappa> (h^2/12 = 1e-3)
    "atmosphere_hydrostatic_pair": 1.8e-2,    # measured 6.1e-3
    "atmosphere_gradient_pair": 1.0e-2,       # measured 3.1e-3 (radiative layers)
    # the luminosity equation: this deck's solver tolerances leave it
    # open in proportion to the time step (defect V1); a converged model
    # closes to 4e-6 (measured with the tolerances 1000 times tighter)
    "luminosity_closure": 1.5e-3,             # measured 4.6e-4
    "luminosity_centre": 2.5e-3,              # measured 8.5e-4 (L = eps m at shell 1)
    "luminosity_closure_converged": 2.0e-5,   # what V1's fix would reach
    # N^2 against the rule the file documents: round-off of a difference
    # of nearly equal terms
    "brunt_rule": 1.0e-9,                     # measured 6.3e-11
    # N^2's two forms where the composition is uniform and the layer
    # radiative, as brackets of N^2 P/(g^2 rho): the centred difference
    # of ln rho over a shell (h^2/6 of terms of order 0.6)
    "brunt_forms": 3.0e-3,                    # measured 9.0e-4
    # derived profile columns against their definitions: round-off
    "derived_columns": 1.0e-12,
    # the same quantity in the profile and in a pulse file: round-off of
    # the profile's log columns
    "profile_vs_pulse": 1.0e-12,
    # the fitting point: the surface boundary condition is interpolated
    # in the envelope triangle (tri_delta_teffl = 0.005, tri_delta_logl =
    # 0.04), second order: a triangle 10 times smaller gave jumps 100
    # times smaller
    "fit_jump_lnr": 2.0e-5,                   # measured 5.1e-6
    "fit_jump_lnT": 6.0e-4,                   # measured 2.0e-4
    "fit_jump_lnrho": 1.0e-3,                 # measured 3.3e-4
    "fit_jump_mass": 1.0e-7,                  # measured 2.3e-8 (one shell's dm/M)
    # the photosphere: the hydrostatic and gradient defects of the pair
    # that straddles it, as in the atmosphere
    "photosphere_hydrostatic_pair": 1.0e-2,   # measured 2.8e-3
}

# measured values, filled by the checks (pytest -s prints them at the end;
# YREC_VALIDATION_REPORT=<path> also writes them there)
MEASURED = {}
RECORDING = [True]   # off while a mutated model is checked


def record(name, value):
    if RECORDING[0]:
        MEASURED[name] = max(MEASURED.get(name, 0.0), abs(value))
    return abs(value)


def worst(values):
    return max((abs(v) for v in values if v is not None), default=0.0)


# ---------------------------------------------------------------------
# the runs
# ---------------------------------------------------------------------
def _deck_constant(inlist_text, name, default_file="physics.defaults"):
    """A solar constant of the run: the deck's value, else the default's."""
    for text in (inlist_text, (DEFAULTS / default_file).read_text()):
        m = re.search(r"^\s*" + name + r"\s*=\s*([-+0-9.dDeE]+)", text, re.M)
        if m:
            return float(m.group(1).lower().replace("d", "e"))
    raise KeyError(name)


def _run(workdir, inlist_text):
    return vm.run_yrec(YREC, workdir, inlist_text, timeout=900)


_edit = vm.edit_inlist


@pytest.fixture(scope="module")
def runs(tmp_path_factory):
    if not YREC.exists():
        pytest.skip("build src/yrec first")
    if not (REPO / "input" / "eos" / "opal2006" / "EOSOPAL06Z0.016492").exists():
        pytest.skip("needs the standard solar deck's tables under input/")
    tmp = tmp_path_factory.mktemp("contract")
    inlist = tmp / "inlist_conv"
    conv = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE / f"{DECK}.nml1"),
         str(CASE / f"{DECK}.nml2"), "-o", str(inlist)], capture_output=True, text=True)
    assert conv.returncode == 0, conv.stderr
    base = inlist.read_text()
    base = _edit(base, "&control", ["profile_columns_file = 'profile_all.list'",
                                    "history_columns_file = 'history_all.list'"])
    gyre = _edit(base, "&control", ["pulse_gyre_interval = 100"])
    gyre = _edit(gyre, "&physics", ["write_profile_flag = .true.", "write_pulse_flag = .true.",
                                    "pulse_format = 'GYRE'"])
    short = base.replace("max_model_number(2) = 1000", "max_model_number(2) = 98", 1)
    assert short != base
    fgong = _edit(short, "&physics", ["write_profile_flag = .true.", "write_pulse_flag = .true.",
                                      "pulse_format = 'FGONG'"])
    out = {"gyre": _run(tmp / "gyre", gyre), "fgong": _run(tmp / "fgong", fgong),
           "inlist": gyre}
    if os.environ.get("YREC_TEST_HDF5") == "1":
        gsm = _edit(short, "&physics", ["write_profile_flag = .true.", "write_pulse_flag = .true.",
                                        "pulse_format = 'GSM'"])
        out["gsm"] = _run(tmp / "gsm", gsm)
    return out


class Model:
    """One model of the gyre run: its GYRE file (header, cols), its
    profile (prof_header, prof), its history row, the run's constants,
    and the indices of the fitting point (nz - 1, the last interior
    point) and of the photosphere (iph)."""

    def __init__(self, outdir, number, inlist_text):
        self.number = number
        self.header, self.cols = vm.read_gyre_text(outdir / f"profile{number}.data.GYRE")
        self.prof_header, self.prof = vm.read_profile(outdir / f"profile{number}.data")
        self.hist_header, hist = vm.read_mesa_table(outdir / "history.data")
        row = hist["model_number"].index(self.prof_header["model_number"])
        self.hist = {k: v[row] for k, v in hist.items()}
        self.nz = self.hist["num_zones"]
        self.iph = vm.photosphere_index(self.cols, self.header["R_star"])
        self.n = len(self.cols["r"])
        self.lsun = _deck_constant(inlist_text, "solar_luminosity_cgs")
        self.rsun = _deck_constant(inlist_text, "solar_radius_cgs")
        self.msun = _deck_constant(inlist_text, "solar_mass_cgs")
        self.label = f"model {self.prof_header['model_number']}"


@pytest.fixture(scope="module")
def models(runs):
    numbers = sorted(int(re.match(r"profile(\d+)\.data$", p.name).group(1))
                     for p in runs["gyre"].glob("profile*.data")
                     if re.match(r"profile\d+\.data$", p.name))
    assert len(numbers) >= 6, numbers
    return [Model(runs["gyre"], k, runs["inlist"]) for k in numbers]


# ---------------------------------------------------------------------
# the checks: each takes one Model and raises AssertionError
# ---------------------------------------------------------------------
def check_header(b):
    """The pulse header, the profile header and the history row describe
    one star, in cgs, tied by L = 4 pi R^2 sigma Teff^4 and g = G M/R^2."""
    h, ph, row = b.header, b.prof_header, b.hist
    assert h["version"] == 101, h["version"]
    assert h["n"] == h["rows"] == ph["num_zones"] == len(b.prof["zone"]) == b.n
    assert b.prof["zone"][-1] == 1 and b.prof["zone"][0] == b.n   # zone 1 is the surface
    assert ph["model_number"] == row["model_number"]
    tol = TOL["header"]
    assert record("header: profile age vs history", vm.rel_diff(ph["star_age"], row["star_age"])) <= tol
    assert record("header: profile log_Teff vs history",
                  abs(ph["log_Teff"] - row["log_Teff"])) <= tol
    assert record("header: L_star vs Lsun 10^log_L",
                  vm.rel_diff(h["L_star"], b.lsun * 10.0 ** row["log_L"])) <= tol
    assert record("header: R_star vs Rsun 10^log_R",
                  vm.rel_diff(h["R_star"], b.rsun * 10.0 ** row["log_R"])) <= tol
    teff = 10.0 ** row["log_Teff"]
    assert record("header: L = 4 pi R^2 sigma Teff^4",
                  h["L_star"] / (vm.FOUR_PI * h["R_star"] ** 2 * SIGMA_SB * teff ** 4) - 1.0) <= tol
    assert record("header: g = G M/R^2",
                  vm.YREC_G * h["M_star"] / h["R_star"] ** 2 / 10.0 ** row["log_g"] - 1.0) <= tol
    # the photosphere is a mesh point, and the star's mass and luminosity
    # are the outermost point's
    assert record("header: a mesh point at R_star", b.cols["r"][b.iph] / h["R_star"] - 1.0) <= tol
    assert b.cols["M_r"][-1] == h["M_star"]
    assert record("header: L_r at the top vs L_star", b.cols["L_r"][-1] / h["L_star"] - 1.0) <= tol
    assert b.iph > b.nz - 1 and b.n > b.iph + 1, (b.nz, b.iph, b.n)   # all three regions


def check_star_mass_units(b):
    """The star's mass in grams is its mass in Msun times the run's solar
    mass (defect M317: the solar starting models hold 1.98900e33 g for a
    star_mass of 1 while solar_mass_cgs is 1.9891e33)."""
    assert record("M317: M_star vs star_mass x Msun",
                  vm.rel_diff(b.header["M_star"], b.hist["star_mass"] * b.msun)) <= 1.0e-9


def check_profile_mass_continuous(b):
    """The profile's mass column does not jump across the photosphere
    (defect M317 again: the atmosphere rows are star_mass, the rows below
    m/solar_mass_cgs)."""
    mass = b.prof["mass"]
    assert record("M317: profile mass jump at the photosphere",
                  mass[b.iph + 1] - mass[b.iph]) <= 1.0e-9


def check_monotonic(b):
    """r increases strictly from the centre to the top of the atmosphere,
    m never decreases and never exceeds the star's mass, P decreases
    strictly, T never rises, and everything is positive."""
    c = b.cols
    assert vm.is_strictly_increasing(c["r"]), b.label
    assert all(b2 >= a for a, b2 in zip(c["M_r"], c["M_r"][1:])), b.label
    assert c["M_r"][-1] <= b.header["M_star"] and c["M_r"][0] > 0.0
    assert vm.is_strictly_increasing([-p for p in c["P"]]), b.label
    assert all(b2 <= a * (1.0 + 1.0e-12) for a, b2 in zip(c["T"], c["T"][1:])), b.label
    for name in ("r", "M_r", "P", "T", "rho", "Gamma_1", "nabla_ad", "delta", "kap"):
        assert min(c[name]) > 0.0, (b.label, name)
    assert min(c["L_r"]) > 0.0 and min(c["eps"]) >= 0.0


def check_interior_equations(b):
    """Hydrostatic equilibrium and mass continuity in the interior, in the
    mass coordinate: to the solver's tolerance in YREC's own scheme (the
    pair form), and to the scheme's truncation error one order higher (the
    three-point form). dlnT/dlnP is the file's nabla."""
    t = vm.structure_terms(b.cols)
    lo, hi = 0, b.nz - 1
    hyd = worst(vm.pair_defects(t["lnP"], t["lnm"], t["dlnP_dlnm"], lo, hi))
    con = worst(vm.pair_defects(t["lnr"], t["lnm"], t["dlnr_dlnm"], lo, hi))
    hyd3 = worst(vm.triple_defects(t["lnP"], t["lnm"], t["dlnP_dlnm"], lo, hi))
    con3 = worst(vm.triple_defects(t["lnr"], t["lnm"], t["dlnr_dlnm"], lo, hi))
    grad = worst(vm.pair_defects(t["lnT"], t["lnP"], b.cols["nabla"], lo, hi))
    assert record("interior: hydrostatic, pair (ln m)", hyd) <= TOL["interior_hydrostatic_pair"], b.label
    assert record("interior: continuity, pair (ln m)", con) <= TOL["interior_continuity_pair"], b.label
    assert record("interior: hydrostatic, three-point", hyd3) <= TOL["interior_hydrostatic_triple"], b.label
    assert record("interior: continuity, three-point", con3) <= TOL["interior_continuity_triple"], b.label
    assert record("interior: dlnT/dlnP = nabla, pair", grad) <= TOL["interior_gradient_pair"], b.label


def check_envelope_equations(b):
    """The same laws in the envelope, in the pressure coordinate, by the
    three-point rule on the output mesh."""
    t = vm.structure_terms(b.cols)
    lo, hi = b.nz, b.iph
    hyd = worst(vm.triple_defects(b.cols["r"], t["lnP"], t["dr_dlnP"], lo, hi))
    con = worst(vm.triple_defects(b.cols["M_r"], t["lnP"], t["dm_dlnP"], lo, hi))
    grad = worst(vm.triple_defects(t["lnT"], t["lnP"], b.cols["nabla"], lo, hi))
    assert record("envelope: hydrostatic, three-point (ln P)", hyd) <= TOL["envelope_hydrostatic_triple"], b.label
    assert record("envelope: continuity, three-point (ln P)", con) <= TOL["envelope_continuity_triple"], b.label
    assert record("envelope: dlnT/dlnP = nabla, three-point", grad) <= TOL["envelope_gradient_triple"], b.label


def _radiative_atmosphere_pairs(b):
    c = b.cols
    return [k for k in range(b.iph, b.n - 1)
            if c["nabla"][k] < c["nabla_ad"][k] and c["nabla"][k + 1] < c["nabla_ad"][k + 1]]


def check_atmosphere_equations(b):
    """The atmosphere: its heights follow hydrostatic equilibrium, its
    mass is the star's (it adds none), its luminosity is constant, and
    where it is radiative dlnT/dlnP is the file's nabla."""
    c = b.cols
    t = vm.structure_terms(c)
    lo, hi = b.iph + 1, b.n - 1
    hyd = worst(vm.pair_defects(c["r"], t["lnP"], t["dr_dlnP"], lo, hi))
    assert record("atmosphere: hydrostatic, pair (ln P)", hyd) <= TOL["atmosphere_hydrostatic_pair"], b.label
    assert all(m == b.header["M_star"] for m in c["M_r"][b.iph + 1:]), b.label
    assert all(lum == c["L_r"][-1] for lum in c["L_r"][b.iph + 1:]), b.label
    grad = vm.pair_defects(t["lnT"], t["lnP"], c["nabla"], b.iph, hi)
    radiative = _radiative_atmosphere_pairs(b)
    assert radiative, b.label
    worst_grad = max(abs(grad[k - b.iph]) for k in radiative if k >= lo)
    assert record("atmosphere: dlnT/dlnP = nabla, pair, radiative", worst_grad) <= TOL["atmosphere_gradient_pair"], b.label
    across = worst(vm.pair_defects(c["r"], t["lnP"], t["dr_dlnP"], b.iph, b.iph + 1))
    assert record("photosphere: hydrostatic across it, pair", across) <= TOL["photosphere_hydrostatic_pair"], b.label


def check_atmosphere_gradient_everywhere(b):
    """dlnT/dlnP is the file's nabla in every atmosphere layer (defect V3:
    where the atmosphere is convectively unstable the file holds grad_ad,
    while the temperature follows the T(tau) relation, grad_rad)."""
    t = vm.structure_terms(b.cols)
    grad = worst(vm.pair_defects(t["lnT"], t["lnP"], b.cols["nabla"], b.iph, b.n - 1))
    assert record("V3: atmosphere dlnT/dlnP = nabla, every layer", grad) <= TOL["atmosphere_gradient_pair"], b.label


def _energy_sources(b):
    return [a + n + g for a, n, g in zip(b.prof["eps_nuc"], b.prof["eps_neu"], b.prof["eps_grav"])]


def _luminosity_closure(b):
    """(L(fit) - L(1) - int eps dm)/L_star and L(1)/(eps(1) m(1)) - 1, with
    eps = eps_nuc + eps_neu + eps_grav of the profile and the trapezoidal
    rule in log m, YREC's scheme (dL/dlog m = ln10 m eps)."""
    c = b.cols
    eps = _energy_sources(b)
    t = vm.structure_terms(c)
    rhs = [e * m for e, m in zip(eps, c["M_r"])]
    total = sum(0.5 * (rhs[k] + rhs[k + 1]) * (t["lnm"][k + 1] - t["lnm"][k]) for k in range(b.nz - 1))
    closure = (c["L_r"][b.nz - 1] - c["L_r"][0] - total) / b.header["L_star"]
    centre = c["L_r"][0] / (eps[0] * c["M_r"][0]) - 1.0
    return closure, centre


def check_luminosity_equation(b):
    """dL/dm = eps_nuc + eps_neu + eps_grav through the interior, and
    L = eps m at the innermost point; no source above the fitting point."""
    closure, centre = _luminosity_closure(b)
    assert record("luminosity: (L - int eps dm)/L_star", closure) <= TOL["luminosity_closure"], b.label
    assert record("luminosity: L(1)/(eps(1) m(1)) - 1", centre) <= TOL["luminosity_centre"], b.label
    c = b.cols
    assert record("luminosity: L(fit) vs L_star", c["L_r"][b.nz - 1] / b.header["L_star"] - 1.0) <= 1.0e-6
    assert all(e == 0.0 for e in c["eps"][b.nz:]), b.label
    # the pulse file's eps is the nuclear rate of the profile
    assert c["eps"][:b.nz] == pytest.approx(b.prof["eps_nuc"][:b.nz], rel=TOL["profile_vs_pulse"], abs=1e-300)


def check_luminosity_equation_converged(b):
    """The luminosity equation closes as the hydrostatic one does (defect
    V1: this deck's solver tolerances leave it open by a few 1e-4 of L on
    the main sequence, in proportion to the time step)."""
    closure, centre = _luminosity_closure(b)
    assert record("V1: (L - int eps dm)/L_star, main sequence", closure) <= TOL["luminosity_closure_converged"], b.label
    assert abs(centre) <= TOL["luminosity_closure_converged"], b.label


def check_brunt(b):
    """N^2: the file's column is the rule it documents, re-derived here
    from rho, P, Gamma_1, r, m, delta and the gradients; and where the
    composition is uniform and the layer radiative the rule's two forms,
    the structure's and the thermal one, agree."""
    c = b.cols
    rule = vm.brunt_n2_file_rule(c)
    scale = [g * g * rho / p for g, rho, p in zip(vm.gravity(c), c["rho"], c["P"])]
    diff = worst([(a - f) / s for a, f, s in zip(rule, c["N2"], scale)])
    assert record("N2: file vs its rule, in g^2 rho/P", diff) <= TOL["brunt_rule"], b.label
    structure, thermal = vm.brunt_brackets(c)
    x = b.prof["h1"]
    mixing = b.prof["mixing_type"]
    compared = 0
    forms = 0.0
    for j in range(3, b.nz - 3):
        if any(mixing[j + d] != 0 for d in (-2, -1, 0, 1, 2)):
            continue
        if max(abs(x[j + d] - x[j]) for d in (-1, 1)) > 1.0e-7:
            continue
        compared += 1
        forms = max(forms, abs(structure[j] - thermal[j]))
    b.homogeneous_points = compared
    if compared:
        assert record("N2: structure vs thermal bracket, uniform radiative", forms) <= TOL["brunt_forms"], b.label
    # a radiative layer has N^2 > 0 and a superadiabatic one N^2 < 0, in
    # the interior away from the boundaries
    for j in range(3, b.nz - 3):
        if all(mixing[j + d] == 0 for d in (-2, -1, 0, 1, 2)):
            assert c["N2"][j] > 0.0, (b.label, j)
        if c["nabla"][j] > c["nabla_ad"][j]:
            assert c["N2"][j] < 0.0, (b.label, j)


def check_brunt_sign_in_convection(b):
    """N^2 is not positive in a convective layer (defect V4: where YREC
    sets nabla = nabla_ad exactly, as in a convective core, the thermal
    form is exactly zero, the rule's test `>= 0` takes the centred
    density difference, and N^2 is that difference's error, of either
    sign)."""
    c = b.cols
    mixing = b.prof["mixing_type"]
    scale = [g * g * rho / p for g, rho, p in zip(vm.gravity(c), c["rho"], c["P"])]
    positive = [c["N2"][j] / scale[j] for j in range(3, b.nz - 3)
                if all(mixing[j + d] == 1 for d in (-2, -1, 0, 1, 2)) and c["N2"][j] > 0.0]
    record("V4: N^2 > 0 in convective layers, largest, in g^2 rho/P", max(positive, default=0.0))
    record("V4: N^2 > 0 in convective layers, points", len(positive))
    assert not positive, (b.label, len(positive), max(positive))


def check_derived_columns(b):
    """The profile's derived columns against their definitions from
    primitive columns: the Lamb frequency, the sound speed, gravity, the
    Schwarzschild ratio, and N^2 (the pulse file's)."""
    c, p = b.cols, b.prof
    tol = TOL["derived_columns"]
    lamb = vm.lamb_s2(c, ell=1)
    cs = vm.sound_speed(c)
    grav = vm.gravity(c)
    assert record("profile: lamb_S2 = 2 Gamma_1 P/(rho r^2)",
                  max(vm.rel_diff(a, e) for a, e in zip(p["lamb_S2"], lamb))) <= tol
    assert record("profile: csound = sqrt(Gamma_1 P/rho)",
                  max(vm.rel_diff(a, e) for a, e in zip(p["csound"], cs))) <= tol
    assert record("profile: grav = G m/r^2 (interior)",
                  max(vm.rel_diff(a, e) for a, e in zip(p["grav"][:b.nz], grav[:b.nz]))) <= tol
    assert record("profile: gradr_div_grada",
                  max(vm.rel_diff(a, r / g) for a, r, g in zip(p["gradr_div_grada"], p["gradr"], p["grada"]))) <= tol
    assert p["brunt_N2"] == c["N2"]
    # the gradient is the radiative one in a radiative shell and lies
    # between the adiabatic and the radiative in a convective one
    for j in range(b.nz):
        if p["mixing_type"][j] == 0:
            assert p["gradT"][j] == pytest.approx(p["gradr"][j], rel=1e-9), (b.label, j)
        else:
            assert p["grada"][j] * (1 - 1e-9) <= p["gradT"][j] <= p["gradr"][j] * (1 + 1e-9), (b.label, j)


PROFILE_TO_PULSE = (("Gamma_1", "gamma1"), ("nabla", "gradT"), ("nabla_ad", "grada"),
                    ("delta", "delta"), ("kap", "opacity"), ("Omega_rot", "omega"))


def check_profile_vs_pulse(b):
    """The same quantity in the profile and in the pulse file."""
    c, p = b.cols, b.prof
    tol = TOL["profile_vs_pulse"]
    worst_rel = 0.0
    for name, values in (("r", [10.0 ** x for x in p["logR"]]), ("P", [10.0 ** x for x in p["logP"]]),
                         ("T", [10.0 ** x for x in p["logT"]]), ("rho", [10.0 ** x for x in p["logRho"]]),
                         ("L_r", [x * b.lsun for x in p["luminosity"]])):
        worst_rel = max(worst_rel, max(vm.rel_diff(a, e) for a, e in zip(c[name], values)))
    for gyre_name, prof_name in PROFILE_TO_PULSE:
        worst_rel = max(worst_rel, max(vm.rel_diff(a, e) for a, e in zip(c[gyre_name], p[prof_name])))
    assert record("profile vs pulse: r P T rho L Gamma_1 gradients delta kap omega", worst_rel) <= tol, b.label
    # the mass, in the interior and envelope rows' unit (solar_mass_cgs)
    worst_m = max(vm.rel_diff(a, e * b.msun) for a, e in zip(c["M_r"][:b.iph + 1], p["mass"][:b.iph + 1]))
    assert record("profile vs pulse: m (below the photosphere)", worst_m) <= tol, b.label
    # the composition: h1 is uniform above the fitting point
    assert all(x == p["h1"][b.nz - 1] for x in p["h1"][b.nz:]), b.label


def _fit_jumps(b):
    """The steps from the last interior point to the first envelope point
    beyond what the laws give for that pair: in ln r against hydrostatic
    equilibrium, in ln T against nabla, in ln rho against the interior's
    own dlnrho/dlnP just below."""
    c = b.cols
    t = vm.structure_terms(c)
    k = b.nz - 1
    dlnp = t["lnP"][k + 1] - t["lnP"][k]
    v = -1.0 / (0.5 * (t["dr_dlnP"][k] / c["r"][k] + t["dr_dlnP"][k + 1] / c["r"][k + 1]))
    jump_r = (t["lnr"][k + 1] - t["lnr"][k]) + dlnp / v
    jump_t = (t["lnT"][k + 1] - t["lnT"][k]) - 0.5 * (c["nabla"][k] + c["nabla"][k + 1]) * dlnp
    slope = math.log(c["rho"][k] / c["rho"][k - 1]) / (t["lnP"][k] - t["lnP"][k - 1])
    jump_rho = math.log(c["rho"][k + 1] / c["rho"][k]) - slope * dlnp
    jump_m = (c["M_r"][k + 1] - c["M_r"][k]) / b.header["M_star"]
    return jump_r, jump_t, jump_rho, jump_m


def check_fit_join(b):
    """The interior and the envelope meet at the fitting point: r, T, rho
    and m continue across it."""
    jump_r, jump_t, jump_rho, jump_m = _fit_jumps(b)
    assert record("fit: jump in ln r", jump_r) <= TOL["fit_jump_lnr"], b.label
    assert record("fit: jump in ln T", jump_t) <= TOL["fit_jump_lnT"], b.label
    assert record("fit: jump in ln rho", jump_rho) <= TOL["fit_jump_lnrho"], b.label
    assert jump_m >= 0.0 and record("fit: step in m/M", jump_m) <= TOL["fit_jump_mass"], b.label


def check_fit_join_density(b):
    """The density does not rise outward across the fitting point (defect
    V2: the envelope's repeat of the fitting point is kept when its radius
    lands above the interior's, a pair about 1e-6 R apart whose density
    rises by the boundary condition's mismatch)."""
    c = b.cols
    k = b.nz - 1
    step = math.log(c["rho"][k + 1] / c["rho"][k])
    record("V2: ln rho step outward at the fit (positive = inversion)", max(step, 0.0))
    assert step <= 0.0, (b.label, step)


CHECKS = [check_header, check_monotonic, check_interior_equations, check_envelope_equations,
          check_atmosphere_equations, check_luminosity_equation, check_brunt,
          check_derived_columns, check_profile_vs_pulse, check_fit_join]


@pytest.mark.parametrize("check", CHECKS, ids=lambda f: f.__name__)
def test_contract(models, check):
    for b in models:
        check(b)


def test_brunt_forms_were_compared(models):
    """check_brunt's comparison of the two forms needs uniform radiative
    layers; at least the pre-main-sequence model with a radiative core
    and the zero-age models have them."""
    for b in models:
        check_brunt(b)
    assert sum(b.homogeneous_points for b in models) >= 100, [b.homogeneous_points for b in models]


@pytest.mark.xfail(strict=True, reason="M317 (docs/bug_sweep.md, open): the solar starting models' "
                   "total mass is 1.98900e33 g, solar_mass_cgs is 1.9891e33; 5.03e-5")
def test_star_mass_is_its_mass_in_grams(models):
    for b in models:
        check_star_mass_units(b)


@pytest.mark.xfail(strict=True, reason="M317 (docs/bug_sweep.md, open): the profile's mass column "
                   "jumps by 5.03e-5 at the photosphere (atmosphere rows are star_mass)")
def test_profile_mass_is_continuous(models):
    for b in models:
        check_profile_mass_continuous(b)


@pytest.mark.xfail(strict=True, reason="V1 (docs/PHYSICS_VALIDATION.md): the standard solar deck's "
                   "solver tolerances leave the luminosity equation open by 4.5e-4 L")
def test_luminosity_equation_is_converged(models):
    for b in models:
        if b.hist["star_age"] > 1.0e8:
            check_luminosity_equation_converged(b)


@pytest.mark.xfail(strict=True, reason="V2 (docs/PHYSICS_VALIDATION.md): a near-duplicate point with "
                   "a density inversion at the fitting point")
def test_no_density_inversion_at_the_fit(models):
    for b in models:
        check_fit_join_density(b)


@pytest.mark.xfail(strict=True, reason="V4 (docs/PHYSICS_VALIDATION.md): N^2 is the centred "
                   "difference's error, of either sign, where nabla = nabla_ad exactly")
def test_brunt_is_not_positive_in_convection(models):
    for b in models:
        check_brunt_sign_in_convection(b)


@pytest.mark.xfail(strict=True, reason="V3 (docs/PHYSICS_VALIDATION.md): the atmosphere's nabla is "
                   "grad_ad where it is convectively unstable, its T(tau) follows grad_rad")
def test_atmosphere_gradient_everywhere(models):
    for b in models:
        check_atmosphere_gradient_everywhere(b)


# ---------------------------------------------------------------------
# the other files of the same models
# ---------------------------------------------------------------------
def test_interior_gyre_file(runs, models):
    """The interior-only GYRE file (pulse_gyre_interval) is the first nz
    points of the model's full file, with the fitting point's m, r and L
    as M_star, R_star and L_star."""
    found = 0
    for b in models:
        path = runs["gyre"] / f"profile{b.prof_header['model_number']:05d}.data.GYRE"
        if not path.exists():
            continue
        found += 1
        header, cols = vm.read_gyre_text(path)
        assert header["version"] == 101 and header["n"] == header["rows"] == b.nz
        for name in vm.GYRE_COLUMNS:
            assert cols[name] == b.cols[name][:b.nz], (b.label, name)
        k = b.nz - 1
        assert (header["M_star"], header["R_star"], header["L_star"]) == \
            (b.cols["M_r"][k], b.cols["r"][k], b.cols["L_r"][k])
    assert found >= 3, found


FGONG_TO_GYRE = ("r", "T", "P", "rho", "L_r", "kap", "eps", "Gamma_1", "nabla_ad", "delta")


def test_fgong_file(runs, models):
    """FGONG holds the same model as the GYRE file: the globals (M, R, L,
    the age, Teff, G), each column it shares, m from ln(m/M), and A4 =
    r N^2/g. The fgong run stops at model 100; its first two profiles are
    the gyre run's."""
    compared = 0
    for b in models[:2]:
        header, var = vm.read_fgong(runs["fgong"] / f"profile{b.number}.data.FGONG")
        assert (header["iconst"], header["ivar"], header["ivers"]) == (15, 40, 1300)
        assert header["nn"] == b.n and header["n_values"] == 15 + 40 * b.n
        glob = header["glob"]
        assert glob[:3] == [b.header["M_star"], b.header["R_star"], b.header["L_star"]]
        assert glob[12] == pytest.approx(b.hist["star_age"], rel=1e-15)
        assert glob[13] == pytest.approx(10.0 ** b.hist["log_Teff"], rel=1e-14)
        assert glob[14] == pytest.approx(vm.YREC_G, rel=1e-15)
        f = vm.fgong_model(header, var)
        tol = TOL["profile_vs_pulse"]
        for name in FGONG_TO_GYRE:
            assert f[name] == pytest.approx(b.cols[name], rel=1e-15, abs=1e-300), (b.label, name)
        assert record("FGONG: m from ln(m/M) vs the GYRE file",
                      max(vm.rel_diff(a, e) for a, e in zip(f["M_r"], b.cols["M_r"]))) <= tol
        assert record("FGONG: N^2 from A4 vs the GYRE file, in g^2 rho/P",
                      max(abs(a - e) * p / (g * g * rho) for a, e, p, g, rho in
                          zip(f["N2"], b.cols["N2"], b.cols["P"], vm.gravity(b.cols), b.cols["rho"]))) <= tol
        assert f["X"] == pytest.approx(b.prof["h1"], rel=1e-15)
        assert f["Z"] == pytest.approx(b.prof["z"], rel=1e-15)
        assert f["R_minus_r"] == pytest.approx([b.header["R_star"] - r for r in b.cols["r"]], rel=1e-9, abs=1.0)
        assert f["eps_grav"][:b.nz] == pytest.approx(b.prof["eps_grav"][:b.nz], rel=1e-15, abs=1e-300)
        compared += 1
    assert compared == 2


def test_gsm_file(runs, models):
    """GSM holds the GYRE text file's header and columns (an HDF5 build:
    YREC_TEST_HDF5=1)."""
    if "gsm" not in runs:
        pytest.skip("set YREC_TEST_HDF5=1 (and build with USE_HDF5=1)")
    if not vm.gsm_readable():
        pytest.skip("reading GSM needs h5py or the h5dump program")
    for b in models[:2]:
        header, cols = vm.read_gsm(runs["gsm"] / f"profile{b.number}.data.GSM")
        assert header["version"] == 110 and header["n"] == b.n
        assert (header["M_star"], header["R_star"], header["L_star"]) == \
            (b.header["M_star"], b.header["R_star"], b.header["L_star"])
        assert sorted(cols) == sorted(vm.GYRE_COLUMNS)
        for name in vm.GYRE_COLUMNS:
            assert cols[name] == pytest.approx(b.cols[name], rel=1e-15, abs=1e-300), (b.label, name)


# ---------------------------------------------------------------------
# each check can fail
# ---------------------------------------------------------------------
def _scaled(b, name, factor, lo=None, hi=None):
    c = b.cols[name]
    for k in range(lo or 0, len(c) if hi is None else hi):
        c[k] *= factor


def mutate_gravitational_constant(b):
    """Another G (GYRE's default, 6.67430e-8): 2.7e-4 in g."""
    _scaled(b, "M_r", 6.67430e-8 / vm.YREC_G)
    b.header["M_star"] *= 6.67430e-8 / vm.YREC_G


def mutate_density_units(b):
    """rho low by 1e-4 everywhere (an eos or units slip)."""
    _scaled(b, "rho", 1.0 - 1.0e-4)


def mutate_envelope_gradient_slots(b):
    """The envelope's actual and adiabatic gradients swapped (stitched_model
    Sources 5, the pre-2026 state)."""
    c = b.cols
    for k in range(b.nz, b.iph + 1):
        c["nabla"][k], c["nabla_ad"][k] = c["nabla_ad"][k], c["nabla"][k]


def mutate_atmosphere_heights(b):
    """The atmosphere points at the photosphere plus one step's length
    (stitched_model Sources 1, the legacy stitch): heights not accumulated."""
    c = b.cols
    r_phot = c["r"][b.iph]
    steps = [c["r"][k] - c["r"][k - 1] for k in range(b.iph + 1, b.n)]
    for i, k in enumerate(range(b.iph + 1, b.n)):
        c["r"][k] = r_phot + steps[i]


def mutate_brunt_thermal_only(b):
    """N^2 from the thermal form everywhere (the composition term dropped:
    the N^2 bug of docs/ROADMAP.md section 9)."""
    b.cols["N2"] = vm.brunt_n2_thermal(b.cols)
    b.prof["brunt_N2"] = list(b.cols["N2"])


def mutate_brunt_gamma1(b):
    """N^2 differenced with Gamma_1 off by 1e-3 (another Gamma_1 than the
    file's column)."""
    saved = list(b.cols["Gamma_1"])
    _scaled(b, "Gamma_1", 1.0 + 1.0e-3)
    b.cols["N2"] = vm.brunt_n2_file_rule(b.cols)
    b.cols["Gamma_1"] = saved
    b.prof["brunt_N2"] = list(b.cols["N2"])


def mutate_delta_sign_convention(b):
    """delta written as dlnrho/dlnT (negative) instead of -dlnrho/dlnT."""
    _scaled(b, "delta", -1.0)


def mutate_luminosity_units(b):
    """L in another solar luminosity (3.846e33 for 3.828e33)."""
    _scaled(b, "L_r", 3.846e33 / 3.828e33)


def mutate_eps_neutrinos_dropped(b):
    """eps_grav left out of the profile (an energy term dropped)."""
    b.prof["eps_grav"] = [0.0] * len(b.prof["eps_grav"])


def mutate_radius_header(b):
    """R_star the fitting point's radius, not the photosphere's (the
    truncation test_mesa_output.py guards against)."""
    b.header["R_star"] = b.cols["r"][b.nz - 1]
    b.iph = vm.photosphere_index(b.cols, b.header["R_star"])


def mutate_lamb_degree(b):
    """lamb_S2 for l = 2 under the l = 1 name."""
    b.prof["lamb_S2"] = [3.0 * x for x in b.prof["lamb_S2"]]


def mutate_fit_temperature(b):
    """The envelope 0.2% hotter than the interior at the fitting point (a
    surface boundary condition ten times worse)."""
    _scaled(b, "T", 1.002, lo=b.nz)


def mutate_profile_pressure(b):
    """The profile's logP is another model's (1e-6 dex off the pulse file)."""
    b.prof["logP"] = [x + 1.0e-6 for x in b.prof["logP"]]


def mutate_radius_order(b):
    """Two neighbouring points exchanged (a non-monotonic radius)."""
    for name in b.cols:
        b.cols[name][b.iph + 2], b.cols[name][b.iph + 3] = b.cols[name][b.iph + 3], b.cols[name][b.iph + 2]


MUTATIONS = [
    (mutate_gravitational_constant, check_interior_equations),
    (mutate_density_units, check_interior_equations),
    (mutate_envelope_gradient_slots, check_envelope_equations),
    (mutate_envelope_gradient_slots, check_brunt),
    (mutate_atmosphere_heights, check_atmosphere_equations),
    (mutate_brunt_thermal_only, check_brunt),
    (mutate_brunt_gamma1, check_brunt),
    (mutate_delta_sign_convention, check_brunt),
    (mutate_luminosity_units, check_luminosity_equation),
    (mutate_luminosity_units, check_header),
    (mutate_eps_neutrinos_dropped, check_luminosity_equation),
    (mutate_radius_header, check_header),
    (mutate_lamb_degree, check_derived_columns),
    (mutate_fit_temperature, check_fit_join),
    (mutate_profile_pressure, check_profile_vs_pulse),
    (mutate_radius_order, check_monotonic),
]


@pytest.mark.parametrize("mutation,check", MUTATIONS,
                         ids=[f"{m.__name__[7:]}->{c.__name__[6:]}" for m, c in MUTATIONS])
def test_each_check_can_fail(models, mutation, check):
    """Every check fails on a model with the named mutation, on at least
    one of the six models (and passes without it: test_contract)."""
    failed = 0
    RECORDING[0] = False
    try:
        for b in models:
            mutated = copy.deepcopy(b)
            mutation(mutated)
            try:
                check(mutated)
            except AssertionError:
                failed += 1
    finally:
        RECORDING[0] = True
    assert failed >= 1, f"{check.__name__} did not notice {mutation.__name__}"
    # the checks named in CHECKS are all exercised by some mutation
    assert {c for _, c in MUTATIONS} == set(CHECKS)


def test_report(models):
    """Print the measured values (pytest -s), and write them to
    YREC_VALIDATION_REPORT when it names a file."""
    lines = [f"{value:10.3e}  {name}" for name, value in sorted(MEASURED.items())]
    text = "\n".join(lines) + "\n"
    print("\nmeasured (largest over the models):\n" + text)
    path = os.environ.get("YREC_VALIDATION_REPORT")
    if path:
        pathlib.Path(path).write_text(text)

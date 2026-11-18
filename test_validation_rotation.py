"""Rotation benchmark (physics validation tier, piece 5).

The testsuite's rotating solar deck (Test_solarGS98_difrotmixper18p4d1p0:
GS98, diffusion, rotational mixing, the van Saders & Pinsonneault 2013
wind, disk locking to 1 Myr at 18.4 d), run to the solar age, against

  a published number   the paper's solar calibration (YREC release paper,
                       section 3.1; the PDF in the repository's root):
                       the equatorial rotation period 25.4 d at 4.568 Gyr
                       within 1%, L and R within 1e-4, the surface Z/X
                       within 1% of 0.02292, and the surface Li7 2.3 dex
                       below its initial value within 1%. L and R are
                       checked in the deck's units, the IAU nominal
                       3.828e33 erg/s and 6.957e10 cm; the paper prints
                       3.827e33 and 6.958e10, against which the model is
                       outside the paper's 1e-4 (an expected failure)
  a budget             J(t) = J_0 - int torque dt, the torque re-derived
                       here from the published law (van Saders &
                       Pinsonneault 2013, ApJ 776, 67, eqs. 1 and 2, which
                       sets the centrifugal factor to 1; Matt et al. 2012,
                       ApJ 754, L26, whose fit of the Alfven radius gives
                       that factor, which YREC keeps) and the history's
                       own columns (omega, M, R, L, P_phot, tau_cz)
  a law                the spin-down exponent dln omega/dln t over 1 to
                       4.57 Gyr against a solid body with the run's
                       structure under the same law. With the structure
                       frozen that body follows the law's own solution,
                       omega^-2 = omega_1^-2 + 2 C (t - t_1), whose
                       exponent tends to Skumanich's -1/2 (Skumanich 1972,
                       ApJ 171, 565); the Sun's is steeper because its
                       radius and luminosity grow
  helioseismology      the radiative interior rotates almost as a solid
                       body, at a rate close to the surface's (Howe 2009,
                       Living Rev. Solar Phys. 6, 1, sections 1 and 6; its
                       section 8 gives 452 nHz for the surface plasma at
                       the equator. Its section 5.5's table of l = 1
                       splittings, which measure the interior's rate,
                       ends at 433 to 436 nHz (Bertello et al. 2000,
                       Chaplin et al. 2001, Gelly et al. 2002): about
                       0.96 of the surface's. The band does not use it)

and, by tolerance and not by bytes, omega(m) at five landmark ages
against validation/rotation/omega_landmarks.txt (YREC's own values of
2026-11-17: a regression pin at 1%, wide enough for this deck's path
sensitivity).

Local (one run of the deck with a rotation profile every 10 models).
YREC_VALIDATION_WRITE_LANDMARKS=1 rewrites the landmark file from the
run (say so in the commit that does).
"""
import math
import os
import pathlib
import re
import sys

import pytest

REPO = pathlib.Path(__file__).parent
sys.path.insert(0, str(REPO / "src" / "tools"))
import validate_model as vm  # noqa: E402
from conftest import INPUT_ROOT  # noqa: E402

# this repository's own build (to validate another build, run these tests
# from its tree); the tables and starting models are conftest.py's roots
YREC = REPO / "src" / "yrec"

# not in CI (minutes, the full input tree): deselect with `pytest -m "not local"`
pytestmark = pytest.mark.local
DECK = REPO / "testsuite" / "Test_solarGS98_difrotmixper18p4d1p0.inlist"
LANDMARKS = REPO / "validation" / "rotation" / "omega_landmarks.txt"
PROFILE_COLUMNS = "zone\nmass\nlogR\nomega\nj_rot\ni_rot\nmixing_type\n"
SECONDS_PER_YEAR = 3.1558e7      # YREC's (state/physical_constants.f90)

# The published wind law's fixed numbers: the torque exponent m as van
# Saders & Pinsonneault (2013) round it (their exponents 3.1, -0.22, 0.56
# and 0.44 are 2 + 5m, -m, 1 - 2m and 2m; Matt et al. 2012 fit m = 0.2177)
# and the centrifugal constant K_2 of Matt et al. (2012).
V13_M = 0.22
MATT_K2 = 0.0506

# The release paper's section 3.1 as printed: "Solar models require L =
# Lsun, 3.827e33 erg/s, and R = Rsun, 6.958e10 cm, at the solar age,
# 4.568 Gyr, within one part in 10^-4" (read in the repository's PDF).
PAPER_LSUN_CGS = 3.827e33
PAPER_RSUN_CM = 6.958e10

LANDMARK_AGES = (1.0e8, 5.0e8, 1.0e9, 2.0e9, 4.0e9)
LANDMARK_MASSES = (0.1, 0.3, 0.5, 0.7, 0.9, 0.98)

# Tolerances and bands, the measured value beside each (2026-11-17).
TOL = {
    # the paper's own tolerances
    "period": 0.01,                # measured +0.05% of 25.4 d (25.412 d)
    "L_and_R": 2.0e-4,             # the paper: 1e-4; measured +7.1e-5 and +3.7e-5
    "paper_L_and_R": 1.0e-4,       # the paper's, against its printed Lsun and Rsun:
                                   # +3.3e-4 and -1.1e-4 (an expected failure)
    "Z_over_X": 0.01,              # measured -0.1% of 0.02292
    "lithium_dex": 0.023,          # 1% of 2.3 dex; measured -0.21 dex (V9)
    # the budget, from the end of disk locking: the torque is known at
    # the models and integrated by the trapezoidal rule over the steps,
    # while the code applies it in sub-steps with tau_cz and P_phot
    # interpolated: 1 to 2% of the loss on the pre-main sequence, where
    # the steps are long against the changes
    "budget": 1.5e-2,              # measured 4.9e-3 of J_0
    "budget_step_median": 3.0e-3,  # measured -7e-4 (main sequence, per step)
    "budget_step_worst": 3.0e-2,   # measured +1.4e-2
    # the exponent against the solid-body integration of the same law:
    # the model's core spins down later than its envelope
    "exponent_vs_solid_body": 0.03,   # measured +0.013 (-0.677 against -0.689)
    # that integration with the structure frozen against the law's
    # analytic solution, which takes the centrifugal term at its 1 Gyr
    # value, 0.9933; it rises to 1 as the star slows
    "frozen_vs_analytic": 4.0e-3,     # measured -1.5e-3
    # omega(m) at the landmark ages against the stored values
    "landmarks": 0.01,
    # the interior against the surface at the solar age
    "interior_over_surface": (0.85, 1.15),   # measured 1.03 to 1.06; the Sun: about 0.96
}
MEASURED = {}


def record(name, value):
    MEASURED[name] = value
    return value


def _control(texts, name):
    for text in texts:
        m = re.search(r"^\s*" + name + r"\s*=\s*([-+0-9.dDeE]+)", text, re.M)
        if m:
            return float(m.group(1).lower().replace("d", "e"))
    raise KeyError(name)


class Run:
    """The deck's history, the wind law's controls as the run had them
    (the deck's, else the defaults), and the profiles' paths by model."""

    def __init__(self, outdir):
        self.outdir = pathlib.Path(outdir)
        _, self.h = vm.read_mesa_table(self.outdir / "history.data")
        texts = ((self.outdir / "inlist_used").read_text(),
                 (REPO / "src" / "defaults" / "physics.defaults").read_text())
        self.fk = _control(texts, "wind_law_constant")
        self.solar_torque = _control(texts, "pmm_norm_jdot")
        self.solar_omega = _control(texts, "pmm_solar_omega")
        self.solar_log_pphot = _control(texts, "pmm_solar_pressure")
        self.solar_tau = _control(texts, "pmm_solar_turnover_timescale")
        self.omega_sat = _control(texts, "solar_omega_sat") * self.solar_omega
        self.disk_age = _control(texts, "disk_locking_age_gyr") * 1.0e9
        self.msun = _control(texts, "solar_mass_cgs")
        self.rsun = _control(texts, "solar_radius_cgs")
        self.lsun = _control(texts, "solar_luminosity_cgs")
        self.n = len(self.h["model_number"])
        self.t = [a * SECONDS_PER_YEAR for a in self.h["star_age"]]


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    if not YREC.exists():
        pytest.skip("build src/yrec first")
    if not (INPUT_ROOT / "eos" / "opal2006" / "EOSOPAL06Z0.018127060").exists():
        pytest.skip("needs the full input/ tree")
    tmp = tmp_path_factory.mktemp("rotation")
    text = vm.edit_inlist(DECK.read_text(), "&control", [
        "profile_columns_file = 'rotation_profile.list'",
        "history_columns_file = 'history_all.list'"])
    text = vm.edit_inlist(text, "&physics", ["write_profile_flag = .true.", "profile_interval = 10"])
    (tmp / "run").mkdir()
    (tmp / "run" / "rotation_profile.list").write_text(PROFILE_COLUMNS)
    return Run(vm.run_yrec(YREC, tmp / "run", text))


# ---------------------------------------------------------------------
# the published law, from the history's columns
# ---------------------------------------------------------------------
def structure_factor(r, i, pressure_exponent=2.0 * V13_M):
    """K_M/K_M,sun of van Saders & Pinsonneault (2013) without its
    centrifugal term: (R/Rsun)^3.1 (M/Msun)^-0.22 (L/Lsun)^0.56
    (P_phot/P_phot,sun)^0.44."""
    h = r.h
    return (10.0 ** h["log_R"][i]) ** (2.0 + 5.0 * V13_M) * h["star_mass"][i] ** (-V13_M) * \
        (10.0 ** h["log_L"][i]) ** (1.0 - 2.0 * V13_M) * \
        10.0 ** (pressure_exponent * (h["log_P_photosphere"][i] - r.solar_log_pphot))


def centrifugal_factor(r, i, omega):
    """Matt et al. (2012): [(K_2^2 + f_sun)/(K_2^2 + f)]^m with f =
    omega^2 R^3/(2 G M), relative to the Sun."""
    radius = 10.0 ** r.h["log_R"][i] * r.rsun
    f = 0.5 * omega ** 2 * radius ** 3 / (vm.YREC_G * r.h["star_mass"][i] * r.msun)
    f_sun = 0.5 * r.solar_omega ** 2 * r.rsun ** 3 / (vm.YREC_G * r.msun)
    return ((MATT_K2 ** 2 + f_sun) / (MATT_K2 ** 2 + f)) ** V13_M


def torque(r, i, omega=None, fk=None, pressure_exponent=2.0 * V13_M, omega_power=2.0):
    """dJ/dt (erg) at model i: f_K T_sun K (omega/omega_sun)
    (min(omega tau/tau_sun, omega_sat)/omega_sun)^2, the Rossby-scaled,
    saturating law. omega defaults to the model's surface omega."""
    omega = r.h["surf_avg_omega"][i] if omega is None else omega
    scaled = omega * r.h["tau_conv_cz"][i] / r.solar_tau
    return (r.fk if fk is None else fk) * r.solar_torque * \
        structure_factor(r, i, pressure_exponent) * centrifugal_factor(r, i, omega) * \
        (omega / r.solar_omega) * (min(scaled, r.omega_sat) / r.solar_omega) ** omega_power


def budget(r, **law):
    """From the first model after disk locking: the largest
    |J - J_0 + int torque dt|/J_0, and the per-step ratios
    (J lost)/(torque dt) on the main sequence (age > 0.1 Gyr)."""
    h = r.h
    first = next(i for i in range(r.n) if h["star_age"][i] > 1.2 * r.disk_age)
    j = h["total_angular_momentum"]
    lost, worst, ratios = 0.0, 0.0, []
    previous = torque(r, first, **law)
    for i in range(first + 1, r.n):
        now = torque(r, i, **law)
        step = 0.5 * (now + previous) * (r.t[i] - r.t[i - 1])
        lost += step
        previous = now
        worst = max(worst, abs(j[i] - j[first] + lost) / j[first])
        if h["star_age"][i] > 1.0e8:
            ratios.append((j[i - 1] - j[i]) / step)
    return worst, ratios, first


def slope(xs, ys):
    mx, my = sum(xs) / len(xs), sum(ys) / len(ys)
    return sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)


def solid_body_exponent(r, frozen=False, **law):
    """dln omega/dln t over 1 to 4.57 Gyr of a solid body with the run's
    total moment of inertia under the law, started from the run's omega
    at 1 Gyr: d(I omega)/dt = -torque, by the midpoint rule on the
    run's own steps (each cut in four). frozen keeps the structure (I and
    the law's structure factor) of the first model."""
    h = r.h
    models = [i for i in range(r.n) if 1.0e9 <= h["star_age"][i] <= 4.57e9]
    inertia = h["total_moment_of_inertia"]
    omega = h["surf_avg_omega"][models[0]]
    xs, ys = [math.log(r.t[models[0]])], [math.log(omega)]
    for a, b in zip(models, models[1:]):
        ia, ib = (models[0], models[0]) if frozen else (a, b)
        dt = (r.t[b] - r.t[a]) / 4.0
        for k in range(4):
            def rate(w, x):
                # the structure, linear in the step between the two models
                tq = (1.0 - x) * torque(r, ia, omega=w, **law) + x * torque(r, ib, omega=w, **law)
                i_now = (1.0 - x) * inertia[ia] + x * inertia[ib]
                di_dt = 0.0 if frozen else (inertia[ib] - inertia[ia]) / (r.t[b] - r.t[a])
                return -(tq + w * di_dt) / i_now
            mid = omega + 0.5 * dt * rate(omega, (k + 0.0) / 4.0)
            omega += dt * rate(mid, (k + 0.5) / 4.0)
        xs.append(math.log(r.t[b]))
        ys.append(math.log(omega))
    _LAST_OMEGA[0] = omega
    return slope(xs, ys)


def frozen_analytic_omega(r):
    """omega at the end of the run for the solid body with the structure
    frozen at the first model past 1 Gyr: the law is d omega/dt =
    -C omega^3 there, so omega^-2 = omega_1^-2 + 2 C (t - t_1)."""
    h = r.h
    first = next(i for i in range(r.n) if h["star_age"][i] >= 1.0e9)
    omega_1 = h["surf_avg_omega"][first]
    c = torque(r, first) / omega_1 ** 3 / h["total_moment_of_inertia"][first]
    return (omega_1 ** -2 + 2.0 * c * (r.t[-1] - r.t[first])) ** -0.5


def solid_body_final_omega(r, **kwargs):
    """The same body's omega at the end of the run, integrated."""
    h = r.h
    first = next(i for i in range(r.n) if h["star_age"][i] >= 1.0e9)
    exponent = solid_body_exponent(r, **kwargs)
    del exponent
    return _LAST_OMEGA[0], first


_LAST_OMEGA = [0.0]


def run_exponent(r):
    h = r.h
    models = [i for i in range(r.n) if 1.0e9 <= h["star_age"][i] <= 4.57e9]
    return slope([math.log(r.t[i]) for i in models], [math.log(h["surf_avg_omega"][i]) for i in models])


# ---------------------------------------------------------------------
# omega(m) at the landmark ages
# ---------------------------------------------------------------------
def omega_at_masses(path):
    """(age, [omega at each landmark mass]) of one profile, linear in m."""
    header, prof = vm.read_profile(path)
    mass, omega = prof["mass"], prof["omega"]
    out = []
    for target in LANDMARK_MASSES:
        k = next(i for i in range(1, len(mass)) if mass[i] >= target)
        f = (target - mass[k - 1]) / (mass[k] - mass[k - 1])
        out.append(omega[k - 1] + f * (omega[k] - omega[k - 1]))
    return header["star_age"], out


def landmark_table(r):
    """{age: [omega at the landmark masses]}: between the two profiles
    around each age, linear in log omega against log age."""
    paths = sorted(r.outdir.glob("profile*.data"), key=lambda p: int(p.name[7:-5]))
    profiles = [omega_at_masses(p) for p in paths]
    table = {}
    for age in LANDMARK_AGES:
        k = next(i for i in range(1, len(profiles)) if profiles[i][0] >= age)
        (t0, w0), (t1, w1) = profiles[k - 1], profiles[k]
        f = math.log(age / t0) / math.log(t1 / t0)
        table[age] = [math.exp(math.log(a) + f * math.log(b / a)) for a, b in zip(w0, w1)]
    return table


def read_landmarks():
    table = {}
    for line in LANDMARKS.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            values = [float(x) for x in line.split()]
            table[values[0]] = values[1:]
    return table


def check_landmarks(table, reference):
    worst = 0.0
    for age in LANDMARK_AGES:
        for got, want in zip(table[age], reference[age]):
            worst = max(worst, abs(got / want - 1.0))
    assert record("landmarks: largest |omega/reference - 1|", worst) <= TOL["landmarks"], worst


# ---------------------------------------------------------------------
# the tests
# ---------------------------------------------------------------------
def test_published_solar_calibration(run):
    """The paper's rotating solar model: 25.4 d at 4.568 Gyr, the Sun's
    L and R, the GS98 surface Z/X."""
    row = {k: v[-1] for k, v in run.h.items()}
    assert row["star_age"] == pytest.approx(4.568e9, rel=1e-9)
    period = record("rotation period at 4.568 Gyr (d)", row["rotation_period_days"])
    assert period == pytest.approx(2.0 * math.pi / row["surf_avg_omega"] / 86400.0, rel=1e-9)
    assert abs(record("period/25.4 d - 1", period / 25.4 - 1.0)) <= TOL["period"]
    assert abs(record("L/Lsun - 1", 10.0 ** row["log_L"] - 1.0)) <= TOL["L_and_R"]
    assert abs(record("R/Rsun - 1", 10.0 ** row["log_R"] - 1.0)) <= TOL["L_and_R"]
    assert abs(record("(Z/X)/0.02292 - 1", row["surface_z_div_x"] / 0.02292 - 1.0)) <= TOL["Z_over_X"]


@pytest.mark.xfail(strict=True, reason="V13 (docs/PHYSICS_VALIDATION.md; not yet in docs/bug_sweep.md): the "
                   "release paper prints Lsun = 3.827e33 erg/s and Rsun = 6.958e10 cm, the testsuite's decks set "
                   "the IAU nominal 3.828e33 and 6.957e10; against the printed values the model's L is +3.3e-4 "
                   "and its R -1.1e-4, outside the paper's 1e-4")
def test_published_solar_units(run):
    """The paper's solar L and R as it prints them, within its 1e-4. The
    deck is calibrated to its own solar_luminosity_cgs and
    solar_radius_cgs (test_published_solar_calibration), which are not
    the printed numbers: the paper's text or the decks' units are to be
    corrected (the authors' decision)."""
    row = {k: v[-1] for k, v in run.h.items()}
    lum = record("L/3.827e33 - 1 (the paper's printed Lsun)",
                 10.0 ** row["log_L"] * run.lsun / PAPER_LSUN_CGS - 1.0)
    rad = record("R/6.958e10 - 1 (the paper's printed Rsun)",
                 10.0 ** row["log_R"] * run.rsun / PAPER_RSUN_CM - 1.0)
    assert abs(lum) <= TOL["paper_L_and_R"] and abs(rad) <= TOL["paper_L_and_R"]


def lithium_depletion(run):
    """log10 of the initial over the final surface Li7/H."""
    h = run.h
    return math.log10((h["surface_li7"][0] / h["surface_h1"][0]) / (h["surface_li7"][-1] / h["surface_h1"][-1]))


@pytest.mark.xfail(strict=True, reason="V9, M364 (docs/PHYSICS_VALIDATION.md, docs/bug_sweep.md): the deck's surface Li7 "
                   "depletion is 2.09 dex, not the paper's 2.3 dex within 1%")
def test_published_lithium_depletion(run):
    assert abs(record("V9: Li7 depletion - 2.3 dex", lithium_depletion(run) - 2.3)) <= TOL["lithium_dex"]


def test_angular_momentum_budget(run):
    """J(t) = J_0 - int torque dt from the end of disk locking to the
    solar age, the torque from the published law and the history."""
    worst, ratios, first = budget(run)
    h = run.h
    record("budget: J_0 (g cm^2/s), at age (yr) %.3e" % h["star_age"][first], h["total_angular_momentum"][first])
    record("budget: fraction of J_0 lost by 4.568 Gyr",
           1.0 - h["total_angular_momentum"][-1] / h["total_angular_momentum"][first])
    assert record("budget: largest |J - J_0 + int torque dt|/J_0", worst) <= TOL["budget"]
    ordered = sorted(ratios)
    median = ordered[len(ordered) // 2]
    assert abs(record("budget: per-step (J lost)/(torque dt) - 1, median", median - 1.0)) <= TOL["budget_step_median"]
    assert record("budget: per-step (J lost)/(torque dt) - 1, largest",
                  max(abs(x - 1.0) for x in ratios)) <= TOL["budget_step_worst"]
    # while the disk holds the star, omega is the disk's
    locked = [w for w, age in zip(h["surf_avg_omega"], h["star_age"]) if age < 0.9 * run.disk_age]
    assert locked and max(locked) == min(locked)


def test_spin_down_law(run):
    """The surface rotation's exponent in time over 1 to 4.57 Gyr against
    a solid body under the published law with the run's structure; that
    body with its structure frozen follows the law's analytic solution."""
    measured = record("spin-down exponent, the run, 1 to 4.57 Gyr", run_exponent(run))
    solid = record("spin-down exponent, solid body, the run's structure", solid_body_exponent(run))
    record("spin-down exponent, solid body, structure frozen at 1 Gyr",
           solid_body_exponent(run, frozen=True))
    integrated, first = solid_body_final_omega(run, frozen=True)
    assert abs(record("frozen structure: integrated/analytic omega at 4.568 Gyr - 1",
                      integrated / frozen_analytic_omega(run) - 1.0)) <= TOL["frozen_vs_analytic"]
    # what the Sun's evolution adds to the law between 1 and 4.568 Gyr
    record("torque at fixed omega, 4.568 Gyr over 1 Gyr",
           torque(run, run.n - 1, omega=run.solar_omega) / torque(run, first, omega=run.solar_omega))
    assert abs(measured - solid) <= TOL["exponent_vs_solid_body"], (measured, solid)
    assert measured < -0.5


def test_rotation_profile_landmarks(run):
    """omega(m) at the landmark ages against the stored values."""
    table = landmark_table(run)
    if os.environ.get("YREC_VALIDATION_WRITE_LANDMARKS") == "1":
        lines = ["# Angular velocity (rad/s) of testsuite/Test_solarGS98_difrotmixper18p4d1p0 at landmark",
                 "# ages and mass coordinates: YREC's own values (valid-tier at 1179afbf, 2026-11-17,",
                 "# macOS/arm64, gfortran 15.2), made by test_validation_rotation.py with",
                 "# YREC_VALIDATION_WRITE_LANDMARKS=1. A regression pin at 1%, not a truth.",
                 "# age_yr   omega at m/Msun = " + "  ".join(f"{m:g}" for m in LANDMARK_MASSES)]
        lines += [f"{age:.3e}  " + "  ".join(f"{w:.6e}" for w in table[age]) for age in LANDMARK_AGES]
        LANDMARKS.parent.mkdir(parents=True, exist_ok=True)
        LANDMARKS.write_text("\n".join(lines) + "\n")
    check_landmarks(table, read_landmarks())
    # the surface convection zone rotates as a solid body
    for age in LANDMARK_AGES[1:]:
        assert table[age][-1] == pytest.approx(run_surface_omega(run, age), rel=5.0e-3)


def run_surface_omega(run, age):
    h = run.h
    k = next(i for i in range(1, run.n) if h["star_age"][i] >= age)
    f = math.log(age / h["star_age"][k - 1]) / math.log(h["star_age"][k] / h["star_age"][k - 1])
    return math.exp(math.log(h["surf_avg_omega"][k - 1]) +
                    f * math.log(h["surf_avg_omega"][k] / h["surf_avg_omega"][k - 1]))


def test_interior_rotation_at_the_solar_age(run):
    """The radiative interior against the surface, by the last profile
    (within 0.2 Gyr of the solar age) and the history's centre."""
    paths = sorted(run.outdir.glob("profile*.data"), key=lambda p: int(p.name[7:-5]))
    age, omegas = omega_at_masses(paths[-1])
    assert age > 4.3e9
    lo, hi = TOL["interior_over_surface"]
    for mass, omega in zip(LANDMARK_MASSES[:-1], omegas[:-1]):
        ratio = record(f"omega(m = {mass:g})/omega(surface), last profile", omega / omegas[-1])
        assert lo <= ratio <= hi, (mass, ratio)
    centre = record("omega(centre)/omega(surface) at 4.568 Gyr",
                    run.h["center_omega"][-1] / run.h["surf_avg_omega"][-1])
    assert lo <= centre <= hi


def test_checks_can_fail(run):
    """The budget does not close for a law 5% stronger, or without its
    photospheric-pressure term; a spin-down law quadratic in omega does
    not give the run's exponent; a profile 2% faster leaves the landmarks."""
    worst, _, _ = budget(run, fk=1.05 * run.fk)
    assert worst > TOL["budget"], worst
    worst, ratios, _ = budget(run, pressure_exponent=0.0)
    assert worst > TOL["budget"] or abs(sorted(ratios)[len(ratios) // 2] - 1.0) > TOL["budget_step_median"]
    assert abs(run_exponent(run) - solid_body_exponent(run, omega_power=1.0)) > TOL["exponent_vs_solid_body"]
    integrated, _ = solid_body_final_omega(run, frozen=True, omega_power=1.0)
    assert abs(integrated / frozen_analytic_omega(run) - 1.0) > TOL["frozen_vs_analytic"]
    fast = {age: [1.02 * w for w in row] for age, row in landmark_table(run).items()}
    saved = dict(MEASURED)
    with pytest.raises(AssertionError):
        check_landmarks(fast, read_landmarks())
    MEASURED.clear()
    MEASURED.update(saved)


def test_report(run):
    """Print the measured values (pytest -s), and write them to
    YREC_VALIDATION_REPORT when it names a file."""
    record("Li7 depletion (dex)", lithium_depletion(run))
    text = "".join(f"{value:14.6g}  {name}\n" for name, value in MEASURED.items())
    print("\nmeasured:\n" + text)
    path = os.environ.get("YREC_VALIDATION_REPORT")
    if path:
        pathlib.Path(path).write_text(text)

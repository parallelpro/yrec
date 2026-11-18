"""Seismic observables against GYRE (physics validation tier, piece 2).

The standard solar model at the solar age (examples/run_standard_solar_
model, Test_solar_noGS_norot), its pulse files given to GYRE (adiabatic
modes), and YREC's own seismic scalars (the history columns delta_nu,
delta_nu_rho, nu_max and delta_Pg) compared with what GYRE finds:

  - GYRE reads each format YREC writes (the GYRE text file, FGONG, and
    GSM with an HDF5 build) and finds the same frequencies from each;
  - delta_nu and delta_Pg are the integrals they are documented to be,
    re-derived here from the pulse file's r, P, rho, Gamma_1 and N^2;
  - the sound-travel integral against GYRE's radial large separation,
    and the buoyancy integral against GYRE's l = 1 period spacing;
  - the scaling relations nu_max and delta_nu_rho against their
    definitions, and delta_nu_rho against GYRE.

Local: it needs GYRE (GYRE_DIR, as GYRE's own setup has it; skipped
cleanly without) and src/yrec. GYRE is told YREC's G (&constants
G_GRAVITY): the GYRE text and GSM files carry none. Measured numbers
are in docs/PHYSICS_VALIDATION.md section 2. Run time: about 80 s
(three YREC runs of 25 s, GYRE under 10 s).

Expected failures (xfail, strict) name two defects, V5 and V6
(docs/PHYSICS_VALIDATION.md).
"""
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

# p modes: radial and dipole, 800 to 4500 uHz (the solar acoustic cutoff
# is about 5000); g modes: dipole, 3 to 150 uHz (orders 3 to about 220)
P_SCANS = [(0, "LINEAR", 800, 4500, 400), (1, "LINEAR", 800, 4500, 400)]
G_SCANS = [(1, "INVERSE", 3, 150, 1500)]

# Tolerances: measured value beside each (2026-11-17, GYRE 9.1.1), and
# what sets it.
TOL = {
    # the same model in two formats: round-off of FGONG's ln(m/M) and A4
    "formats_uhz": 1.0e-9,                 # measured 1.4e-12 uHz
    # a history scalar against its integral re-derived from the pulse
    # file: round-off
    "scalar_rederived": 1.0e-9,
    # delta_nu (sound travel, centre to photosphere) over GYRE's radial
    # large separation at n = 15..25: the integral stops at the
    # photosphere (+3.5%) and those orders are not asymptotic (+1.5%)
    "delta_nu_over_gyre": (1.03, 1.07),    # measured 1.0509
    # the same integral through the atmosphere against GYRE's separation
    # at n = 25..31, the asymptotic law itself
    "asymptotic_delta_nu": 6.0e-3,         # measured -1.7e-3
    # the mean-density scaling against GYRE at n = 15..25 (a scaling
    # relation with the reference 135.1 uHz)
    "delta_nu_rho_vs_gyre": 1.5e-2,        # measured -4.8e-3
    # the buoyancy integral with the centre's part against GYRE's l = 1
    # period spacing at n_g = 120..200 on the centre-extended file
    "asymptotic_period_spacing": 3.0e-3,   # measured 3.5e-4
    # what the two defects would meet when fixed
    "delta_pg_vs_gyre": 1.0e-2,            # V5: measured +6.1e-2
    "gyre_file_vs_extended": 5.0e-3,       # V6: measured -2.4e-2
}
MEASURED = {}


def record(name, value):
    MEASURED[name] = value
    return value


@pytest.fixture(scope="module")
def sun(tmp_path_factory):
    if vm.gyre_program() is None:
        pytest.skip("GYRE_DIR does not name a GYRE installation (bin/gyre)")
    if not YREC.exists():
        pytest.skip("build src/yrec first")
    if not (REPO / "input" / "eos" / "opal2006" / "EOSOPAL06Z0.016492").exists():
        pytest.skip("needs the standard solar deck's tables under input/")
    tmp = tmp_path_factory.mktemp("gyre")
    inlist = tmp / "inlist_conv"
    conv = subprocess.run(
        [sys.executable, str(CONVERTER), str(CASE / f"{DECK}.nml1"),
         str(CASE / f"{DECK}.nml2"), "-o", str(inlist)], capture_output=True, text=True)
    assert conv.returncode == 0, conv.stderr
    formats = ["GYRE", "FGONG"] + (["GSM"] if os.environ.get("YREC_TEST_HDF5") == "1" else [])
    out = vm.run_with_final_profile(YREC, tmp / "yrec", inlist.read_text(), pulse_formats=formats)
    _, hist = vm.read_mesa_table(out["history"] / "history.data")
    header, cols = vm.read_gyre_text(out["GYRE"] / "profile1.data.GYRE")
    result = {"tmp": tmp, "out": out, "formats": formats, "header": header, "cols": cols,
              "inlist": inlist.read_text(),
              "hist": {k: v[-1] for k, v in hist.items()}}
    # GYRE on each format (p modes), and the g modes of the text file
    result["p"] = {fmt: vm.run_gyre(out[fmt] / f"profile1.data.{fmt}", fmt, tmp / "gyre",
                                    P_SCANS, tag=f"p_{fmt}") for fmt in formats}
    result["g"] = vm.run_gyre(out["GYRE"] / "profile1.data.GYRE", "GYRE", tmp / "gyre",
                              G_SCANS, tag="g_GYRE")
    extended = tmp / "gyre" / "centre_extended.GYRE"
    write_gyre_text(extended, header, centre_extended(cols))
    result["g_extended"] = vm.run_gyre(extended, "GYRE", tmp / "gyre", G_SCANS, tag="g_extended")
    return result


def write_gyre_text(path, header, cols):
    n = len(cols["r"])
    lines = ["%6d %26.16E %26.16E %26.16E %6d" % (n, header["M_star"], header["R_star"],
                                                 header["L_star"], header["version"])]
    for k in range(n):
        lines.append("%6d" % (k + 1) + "".join(" %26.16E" % cols[name][k] for name in vm.GYRE_COLUMNS))
    pathlib.Path(path).write_text("\n".join(lines) + "\n")


def centre_extended(cols, points=30):
    """The model with points added between the centre and its innermost
    mesh point, by the expansions about the centre: rho constant,
    m ~ r^3, P = P_c - (2 pi/3) G rho^2 r^2, N^2 ~ r^2 (the other columns
    constant, L ~ m). The innermost mesh point's own N^2, which the file
    copies from its neighbour, follows the same law."""
    r1, m1, p1, rho1 = cols["r"][0], cols["M_r"][0], cols["P"][0], cols["rho"][0]
    k_p = 2.0 * math.pi / 3.0 * vm.YREC_G * rho1 ** 2
    p_c = p1 + k_p * r1 ** 2
    n2_per_r2 = cols["N2"][1] / cols["r"][1] ** 2
    new = {name: [] for name in vm.GYRE_COLUMNS}
    for i in range(points):
        x = (i + 1.0) / (points + 1.0)
        for name in vm.GYRE_COLUMNS:
            new[name].append(cols[name][0])
        new["r"][-1] = x * r1
        new["M_r"][-1] = m1 * x ** 3
        new["L_r"][-1] = cols["L_r"][0] * x ** 3
        new["P"][-1] = p_c - k_p * (x * r1) ** 2
        new["N2"][-1] = n2_per_r2 * (x * r1) ** 2
    out = {name: new[name] + list(cols[name]) for name in vm.GYRE_COLUMNS}
    out["N2"][points] = n2_per_r2 * r1 ** 2
    return out


def yrec_delta_pg_rule(cols, n_ie):
    """delta_Pg as observables.f90 documents it: sqrt(2) pi^2 over the sum
    of N/r times half the span to the neighbours, from the second point
    outward while N^2 > 0, to the first point where it is not."""
    total, inside = 0.0, False
    for k in range(1, n_ie - 1):
        n2 = cols["N2"][k]
        if n2 > 0.0:
            inside = True
            total += math.sqrt(n2) / cols["r"][k] * 0.5 * (cols["r"][k + 1] - cols["r"][k - 1])
        elif inside:
            break
    return math.sqrt(2.0) * math.pi ** 2 / total


def buoyancy_period_spacing(cols):
    """sqrt(2) pi^2 / int_0^r2 N dr/r over the cavity from the centre:
    the trapezoid in ln r from the innermost point plus the centre's part,
    N(r_1) (N is proportional to r there)."""
    total, k = math.sqrt(cols["N2"][1]) * cols["r"][0] / cols["r"][1], 0
    while cols["N2"][k + 1] > 0.0:
        n_a = math.sqrt(cols["N2"][k]) if k > 0 else math.sqrt(cols["N2"][1]) * cols["r"][0] / cols["r"][1]
        total += 0.5 * (n_a + math.sqrt(cols["N2"][k + 1])) * math.log(cols["r"][k + 1] / cols["r"][k])
        k += 1
    return math.sqrt(2.0) * math.pi ** 2 / total


def test_gyre_reads_every_format(sun):
    """GYRE accepts each pulse file and finds the same p modes from each."""
    ref0, ref1 = vm.mode_table(sun["p"]["GYRE"], 0), vm.mode_table(sun["p"]["GYRE"], 1)
    assert len(ref0) >= 25 and len(ref1) >= 25, (len(ref0), len(ref1))
    assert sorted(ref0) == list(range(min(ref0), max(ref0) + 1))   # no radial order missed
    for fmt in sun["formats"][1:]:
        for ell, ref in ((0, ref0), (1, ref1)):
            got = vm.mode_table(sun["p"][fmt], ell)
            assert sorted(got) == sorted(ref), fmt
            worst = max(abs(got[n] - ref[n]) for n in ref)
            assert record(f"{fmt} vs GYRE text, l = {ell}, max |dnu| (uHz)", worst) <= TOL["formats_uhz"]


def test_delta_nu_is_the_sound_travel_integral(sun):
    """The history's delta_nu is [2 int_0^R dr/c]^-1 over the pulse
    file's points up to the photosphere."""
    rederived = vm.acoustic_delta_nu(sun["cols"], sun["header"]["R_star"])
    record("delta_nu, history (uHz)", sun["hist"]["delta_nu"])
    assert abs(rederived / sun["hist"]["delta_nu"] - 1.0) <= TOL["scalar_rederived"]


def test_delta_nu_against_gyre(sun):
    """The sound-travel integral against GYRE's radial large separation."""
    radial = vm.mode_table(sun["p"]["GYRE"], 0)
    dnu_mid, count = vm.large_separation(radial, 15, 25)
    assert count == 11
    record("GYRE large separation, l = 0, n = 15..25 (uHz)", dnu_mid)
    record("GYRE nu(l = 0, n = 21) (uHz)", radial[21])
    ratio = record("delta_nu / GYRE(15..25)", sun["hist"]["delta_nu"] / dnu_mid)
    lo, hi = TOL["delta_nu_over_gyre"]
    assert lo <= ratio <= hi, ratio
    # the asymptotic law: the integral through the atmosphere against the
    # highest orders below the cutoff
    dnu_high, count = vm.large_separation(radial, 25, 31)
    assert count == 7
    record("GYRE large separation, l = 0, n = 25..31 (uHz)", dnu_high)
    through = record("[2 int dr/c]^-1 to the top of the atmosphere (uHz)",
                     vm.acoustic_delta_nu(sun["cols"], sun["cols"]["r"][-1]))
    assert abs(record("integral through the atmosphere / GYRE(25..31) - 1",
                      through / dnu_high - 1.0)) <= TOL["asymptotic_delta_nu"]


def test_scaling_relations(sun):
    """nu_max and delta_nu_rho are their definitions (the controls'
    solar references, the run's own solar constants), and delta_nu_rho
    is near GYRE's large separation."""
    row = sun["hist"]
    # a control's value on this run: the deck's, else the default
    texts = (sun["inlist"], (REPO / "src" / "defaults" / "physics.defaults").read_text())

    def control(name):
        for text in texts:
            m = re.search(r"^\s*" + name + r"\s*=\s*([-+0-9.dDeE]+)", text, re.M)
            if m:
                return float(m.group(1).lower().replace("d", "e"))
        raise KeyError(name)
    g_sun = vm.YREC_G * control("solar_mass_cgs") / control("solar_radius_cgs") ** 2
    nu_max = control("nu_max_sun") * 10.0 ** row["log_g"] / g_sun * \
        math.sqrt(control("Teff_sun") / 10.0 ** row["log_Teff"])
    assert nu_max == pytest.approx(row["nu_max"], rel=1e-9)
    dnu_rho = control("delta_nu_sun") * math.sqrt(row["star_mass"]) * 10.0 ** (-1.5 * row["log_R"])
    assert dnu_rho == pytest.approx(row["delta_nu_rho"], rel=1e-9)
    dnu_mid, _ = vm.large_separation(vm.mode_table(sun["p"]["GYRE"], 0), 15, 25)
    record("delta_nu_rho, history (uHz)", row["delta_nu_rho"])
    assert abs(record("delta_nu_rho / GYRE(15..25) - 1", row["delta_nu_rho"] / dnu_mid - 1.0)) \
        <= TOL["delta_nu_rho_vs_gyre"]


def test_delta_pg_is_the_buoyancy_sum(sun):
    """The history's delta_Pg is the sum observables.f90 documents, over
    the pulse file's N^2 (to the photosphere)."""
    n_ie = vm.photosphere_index(sun["cols"], sun["header"]["R_star"]) + 1
    record("delta_Pg, history (s)", sun["hist"]["delta_Pg"])
    assert yrec_delta_pg_rule(sun["cols"], n_ie) == pytest.approx(sun["hist"]["delta_Pg"],
                                                                  rel=TOL["scalar_rederived"])


def test_buoyancy_integral_against_gyre(sun):
    """The asymptotic law for g modes: sqrt(2) pi^2/int N dr/r, with the
    centre's part, is GYRE's l = 1 period spacing at high order, on the
    file with the centre resolved."""
    spacing, count = vm.mean_period_spacing(vm.mode_table(sun["g_extended"], 1), 120, 200)
    assert count == 81
    record("GYRE period spacing, l = 1, n_g = 120..200, centre-extended file (s)", spacing)
    integral = record("sqrt(2) pi^2 / int N dr/r with the centre (s)",
                      buoyancy_period_spacing(sun["cols"]))
    assert abs(record("buoyancy integral / GYRE(centre-extended) - 1",
                      integral / spacing - 1.0)) <= TOL["asymptotic_period_spacing"]


@pytest.mark.xfail(strict=True, reason="V5 (docs/PHYSICS_VALIDATION.md): delta_Pg's sum starts at "
                   "the second mesh point and leaves out the centre's part, 6% of the integral")
def test_delta_pg_against_gyre(sun):
    spacing, _ = vm.mean_period_spacing(vm.mode_table(sun["g_extended"], 1), 120, 200)
    assert abs(record("V5: delta_Pg / GYRE(centre-extended) - 1",
                      sun["hist"]["delta_Pg"] / spacing - 1.0)) <= TOL["delta_pg_vs_gyre"]


@pytest.mark.xfail(strict=True, reason="V6 (docs/PHYSICS_VALIDATION.md): the pulse files' innermost "
                   "point is at r = 0.01 R; GYRE's g-mode periods depend on how it fills the centre")
def test_g_modes_do_not_depend_on_the_centre_fill(sun):
    as_written, _ = vm.mean_period_spacing(vm.mode_table(sun["g"], 1), 120, 200)
    extended, _ = vm.mean_period_spacing(vm.mode_table(sun["g_extended"], 1), 120, 200)
    record("GYRE period spacing, l = 1, n_g = 120..200, file as written (s)", as_written)
    assert abs(record("V6: GYRE(as written) / GYRE(centre-extended) - 1",
                      as_written / extended - 1.0)) <= TOL["gyre_file_vs_extended"]


def test_p_modes_do_not_depend_on_the_centre_fill(sun):
    """The p modes do not: the centre-extended file gives the same radial
    frequencies (so V6 is a g-mode matter)."""
    extended = sun["tmp"] / "gyre" / "centre_extended.GYRE"
    radial = vm.mode_table(vm.run_gyre(extended, "GYRE", sun["tmp"] / "gyre", P_SCANS[:1],
                                       tag="p_extended"), 0)
    ref = vm.mode_table(sun["p"]["GYRE"], 0)
    worst = max(abs(radial[n] - ref[n]) for n in ref)
    assert record("p modes, centre-extended vs as written, max |dnu| (uHz)", worst) <= 5.0e-3


# ---------------------------------------------------------------------
# each comparison can fail: GYRE is given a file that no longer is the
# model YREC's scalars describe
# ---------------------------------------------------------------------
def test_checks_can_fail(sun):
    tmp = sun["tmp"] / "gyre"
    # Gamma_1 high by 2%: the sound speed by 1%, so GYRE's separation
    # leaves the asymptotic integral of the model as YREC wrote it
    cols = {k: list(v) for k, v in sun["cols"].items()}
    cols["Gamma_1"] = [1.02 * x for x in cols["Gamma_1"]]
    write_gyre_text(tmp / "mutant_gamma1.GYRE", sun["header"], cols)
    radial = vm.mode_table(vm.run_gyre(tmp / "mutant_gamma1.GYRE", "GYRE", tmp, P_SCANS[:1],
                                       tag="p_mutant"), 0)
    dnu_high, _ = vm.large_separation(radial, 25, 31)
    through = vm.acoustic_delta_nu(sun["cols"], sun["cols"]["r"][-1])
    assert abs(through / dnu_high - 1.0) > TOL["asymptotic_delta_nu"]
    # ... and a delta_nu computed with that Gamma_1 is not the history's
    assert abs(vm.acoustic_delta_nu(cols, sun["header"]["R_star"]) / sun["hist"]["delta_nu"] - 1.0) \
        > TOL["scalar_rederived"]
    # N^2 high by 5% in the file GYRE reads: its period spacing leaves
    # the integral of YREC's N^2 by 2.5%
    ext = centre_extended(sun["cols"])
    ext["N2"] = [1.05 * x for x in ext["N2"]]
    write_gyre_text(tmp / "mutant_n2.GYRE", sun["header"], ext)
    g = vm.mode_table(vm.run_gyre(tmp / "mutant_n2.GYRE", "GYRE", tmp, G_SCANS, tag="g_mutant"), 1)
    spacing, _ = vm.mean_period_spacing(g, 120, 200)
    assert abs(buoyancy_period_spacing(sun["cols"]) / spacing - 1.0) > TOL["asymptotic_period_spacing"]
    # ... and the sum over that N^2 is not the history's delta_Pg
    cols = {k: list(v) for k, v in sun["cols"].items()}
    cols["N2"] = [1.05 * x for x in cols["N2"]]
    n_ie = vm.photosphere_index(cols, sun["header"]["R_star"]) + 1
    assert abs(yrec_delta_pg_rule(cols, n_ie) / sun["hist"]["delta_Pg"] - 1.0) > TOL["scalar_rederived"]


def test_report(sun):
    """Print the measured values (pytest -s), and write them to
    YREC_VALIDATION_REPORT when it names a file."""
    text = "".join(f"{value:14.6f}  {name}\n" for name, value in sorted(MEASURED.items()))
    print("\nmeasured:\n" + text)
    path = os.environ.get("YREC_VALIDATION_REPORT")
    if path:
        pathlib.Path(path).write_text(text)

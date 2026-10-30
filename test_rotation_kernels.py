"""Rotation's kernels against the laws they keep.

Runs src/test_rotation_kernels (built by `make tests` in src/) and
byte-compares its stdout with the pinned
src/rotation/test/expected_test_rotation_kernels.out. On a synthetic
60-zone model (no deck, no tables) it checks the spherical limits of
the moments of inertia and the shape factors and how fp and ft grow with
omega; that omega_from_j gives each radiative zone and the convective
envelope back their angular momentum, the envelope as a solid body; that
the angular momentum diffusion conserves J without a wind, loses exactly
the wind's surface term with one, and keeps a solid body; that the
composition mixing conserves each species; the Kawaler wind law's
exponents below and above saturation; the wind's spin-down against the
analytic laws, at second order in the step; the angular momentum
diffusion's eigenmodes and its convergence to the continuum's decay (h^2
in space, dt in time); the instability velocities' laws
(omega^2 scalings, no shear term or GSF for a solid body, the mu
factor's bounds, the shear's onset and stabilization); the region scan;
the transport grid and the map back's conservation; one region's
transport and the wind's terms for it; the stability
coefficients against their formulas; and check_angular_momentum's and
check_composition's contracts. The binary exits nonzero when a
check fails (see the header of src/rotation/test/test_rotation_kernels.f90).
See test_eos.py for the conventions.
"""
import pathlib
import re
import subprocess

import pytest

from conftest import compare_output

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_rotation_kernels"
EXPECTED = REPO / "src" / "rotation" / "test" / "expected_test_rotation_kernels.out"


def test_mechanism_controls():
    """Each transport mechanism's two controls, as
    src/state/transport_mechanism_layout.f90 names them, are registry
    controls, and rotation_settings_from_star copies them into the
    mechanism's slot of velocity_scale and mixing_scale."""
    layout = (REPO / "src" / "state" / "transport_mechanism_layout.f90").read_text()
    slots = {int(k): name for name, k in
             re.findall(r"(i_mech_\w+)\s*=\s*(\d+)", layout)}

    def names(array):
        m = re.search(array + r"\(n_transport_mechanisms\)\s*=\s*&?\s*"
                      r"\[character\(len=32\) ::([^\]]*)\]", layout)
        assert m, array
        return re.findall(r"'([^']+)'", m.group(1))

    velocity = names("transport_mechanism_velocity_scale_control")
    mixing = names("transport_mechanism_mixing_scale_control")
    assert len(velocity) == len(mixing) == len(slots)
    registry = (REPO / "src" / "defaults" / "controls_registry.tsv").read_text()
    canonical = {line.split("\t")[3] for line in registry.splitlines()
                 if line.count("\t") >= 3 and not line.startswith("#")}
    setup = (REPO / "src" / "setup" / "rotation_settings_setup.f90").read_text()
    for k, slot in sorted(slots.items()):
        for field, control in (("velocity_scale", velocity[k - 1]),
                               ("mixing_scale", mixing[k - 1])):
            assert control in canonical, f"{control} is not a registry control"
            copy = rf"settings%{field}\({slot}\)\s*=\s*s%physics%{control}\b"
            assert re.search(copy, setup), \
                f"rotation_settings_from_star does not copy {control} into {field}({slot})"


def test_rotation_kernels(tmp_path):
    if not BINARY.exists():
        pytest.skip("src/test_rotation_kernels not built (run `make tests` in src/)")
    result = subprocess.run(
        [str(BINARY)], cwd=tmp_path, capture_output=True, text=True
    )
    assert result.returncode == 0, result.stdout + result.stderr
    # The checks' own tolerances decide (the binary exits nonzero on a
    # FAIL); the pin then keeps the list of checks and their statuses.
    # Their relative errors are round-off and differ between platforms,
    # so a check's line is compared by its label and status only.
    check = re.compile(r"^(.*\S)\s+[-+0-9.Ee]+\s+(ok|FAIL)$")
    def masked(text):
        out = []
        for line in text.splitlines():
            m = check.match(line)
            out.append(f"{m.group(1)} [{m.group(2)}]" if m else line)
        return "\n".join(out) + "\n"
    compare_output(masked(result.stdout), masked(EXPECTED.read_text()),
                   "test_rotation_kernels", abs_floor=1.0e-12)

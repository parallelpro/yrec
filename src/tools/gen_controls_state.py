#!/usr/bin/env python3
"""Generate star%ctrl's type body and the buffer sync module.

Source of truth: io/controls_lib.f90's declarations (names, types,
dims, and declaration-time defaults -- the same single place the
legacy read has always taken its pristine state from). Emits:

  state/controls_state_def.inc  -- the component list of the
      controls_state type (star_info_lib defines
      `type, public :: controls_state` around an include of this
      file). Every component gets a default initializer: the
      buffer's own where it has one, otherwise the explicit
      zero/blank/.false. that static module storage gave it -- so
      `star%ctrl = controls_state()` reproduces the pristine
      pre-read_input state exactly.

  state/controls_sync_lib.f90   -- module controls_sync_lib with
      seed_controls_buffer  (buffer member = star%ctrl%member) and
      store_controls_to_star (star%ctrl%member = buffer member),
      one assignment per member, in declaration order.

Regenerate (and commit the results) whenever controls_lib's member
list changes:  python3 tools/gen_controls_state.py
The build breaks loudly on drift (unknown/missing components).
"""
import pathlib
import re

SRC = pathlib.Path(__file__).resolve().parent.parent
LIB = SRC / "io" / "controls_lib.f90"
DEF = SRC / "state" / "controls_state_def.inc"
JOBDEF = SRC / "state" / "job_controls_def.inc"
SYNC = SRC / "state" / "controls_sync_lib.f90"

# Members that live on star%job instead of star%ctrl: the run-list /
# calibration-protocol card arrays and latches, which the calibration
# protocol (setcal/chkcal/setscal/chkscal), the stop-disarm pass, and
# the MC loop WRITE after the read -- mutable job state, excluded
# from ctrl's immutable-after-read contract. Still namelist/card
# targets: read_input reads them into the buffer; the sync store routes
# them to star%job.
JOB_MEMBERS = {
    "mixing_length_array", "rescale_params", "rescale_kind",
    "num_models", "num_runs", "initial_x_array", "initial_z_array",
    "has_senv0_array", "senv0_array", "rsclzc", "rsclzm1", "rsclzm2",
    "target_end_age", "end_age_stop_active", "central_deuterium_stop",
    "central_hydrogen_stop", "central_helium_stop",
    "timestep_override", "timestep_override_active", "first_call_flag",
    # batch 3: namelist-read members the run mutates -- the
    # model-restore set (getyrec7/getmodel2 write them back from
    # stored models), driver/output toggles, per-card working copies,
    # and config the physics adjusts in place.
    "rotation_active", "use_extended_composition", "core_overshoot_active", "lovstm",
    "envelope_overshoot_active", "use_semiconvection",
    "instability_transport_active", "disk_locking_active",
    "disk_omega_rad_s", "disk_locking_age_gyr", "use_diffusion_z",
    "wind_saturation_omega", "use_wind_torque", "atm_choice",
    "use_structure_dt_limits", "atm_step_begin", "atm_step_min",
    "atm_step_max", "env_step_begin", "env_step_min", "env_step_max",
    "calc_envelope_flag", "initial_envelope_x",
    "initial_envelope_z", "requested_envelope_mass",
    "change_envelope_mass_flag", "core_mass_reduction_factor",
    "num_core_shells_added", "fcorr", "num_rotation_structure_iters", "diffusion_timestep_factor", "mc_run_start", "mc_run_end",
    "log_L_upper_limit", "log_L_lower_limit", "Teff_upper_limit",
    "Teff_lower_limit", "log_g_upper_limit", "log_g_lower_limit",
    "nu_max_upper_limit", "nu_max_lower_limit",
    "acfpft", "max_domega_global", "structfactor", "fgrz",
    "diffuse_helium_active", "use_mass_accretion", "new_species_value",
    "target_radius_rsun", "target_teff", "pulse_format",
}

TYPE_RE = re.compile(
    r"^\s*(double precision|integer|logical|character\s*\(\s*len\s*=\s*\d+\s*\))"
    r"\s*(,\s*public)?\s*::\s*(.*)$",
    re.IGNORECASE,
)

ZERO = {
    "double precision": "0.0d0",
    "integer": "0",
    "logical": ".false.",
    "character": "' '",
}


def joined_decl_lines(text):
    """Join &-continued lines, strip comments."""
    out, buf = [], ""
    for raw in text.splitlines():
        code = raw.split("!")[0]
        if buf:
            code = code.lstrip().lstrip("&")
        buf += code
        if buf.rstrip().endswith("&"):
            buf = buf.rstrip()[:-1]
            continue
        if buf.strip():
            out.append(buf)
        buf = ""
    return out


def split_entities(s):
    """Split a declaration's entity list on top-level commas."""
    ents, depth, cur = [], 0, ""
    for ch in s:
        if ch in "([":
            depth += 1
        elif ch in ")]":
            depth -= 1
        if ch == "," and depth == 0:
            ents.append(cur.strip())
            cur = ""
        else:
            cur += ch
    if cur.strip():
        ents.append(cur.strip())
    return ents


def parse():
    members = []  # (typespec, name, dims, init)
    for line in joined_decl_lines(LIB.read_text()):
        m = TYPE_RE.match(line)
        if not m:
            continue
        typespec = re.sub(r"\s+", " ", m.group(1).lower()).replace(" (", "(")
        for ent in split_entities(m.group(3)):
            em = re.match(r"^(\w+)\s*(\([^)]*\))?\s*(=\s*(.*))?$", ent, re.S)
            if not em:
                raise SystemExit(f"cannot parse entity: {ent!r}")
            name, dims, init = em.group(1), em.group(2) or "", em.group(4)
            if init is None:
                base = "character" if typespec.startswith("character") else typespec
                init = ZERO[base]
            members.append((typespec, name, dims, init.strip()))
    if not members:
        raise SystemExit("no members parsed from controls_lib.f90")
    return members


def main():
    members = parse()
    found_job = {n for _, n, _, _ in members} & JOB_MEMBERS
    missing = JOB_MEMBERS - found_job
    if missing:
        raise SystemExit(f"JOB_MEMBERS not found in controls_lib: {missing}")
    with DEF.open("w") as f:
        f.write(
            "! GENERATED by tools/gen_controls_state.py from\n"
            "! io/controls_lib.f90 -- do not edit; regenerate instead.\n"
            "! One component per buffer member, declaration order, every\n"
            "! one default-initialized (see the generator's header).\n"
        )
        for typespec, name, dims, init in members:
            if name in JOB_MEMBERS:
                continue
            f.write(f"            {typespec} :: {name}{dims} = {init}\n")
    with JOBDEF.open("w") as f:
        f.write(
            "! GENERATED by tools/gen_controls_state.py from\n"
            "! io/controls_lib.f90 -- do not edit; regenerate instead.\n"
            "! The JOB_MEMBERS subset (run-list / calibration-protocol\n"
            "! card arrays and latches -- mutable job state), included\n"
            "! into star_job. Default-initialized like controls_state, so\n"
            "! declaration/star0-restore state is pristine.\n"
        )
        for typespec, name, dims, init in members:
            if name not in JOB_MEMBERS:
                continue
            f.write(f"            {typespec} :: {name}{dims} = {init}\n")
    with SYNC.open("w") as f:
        f.write(
            "!----------------------------------------------------------------------\n"
            "! controls_sync_lib\n"
            "!----------------------------------------------------------------------\n"
            "! GENERATED by tools/gen_controls_state.py from\n"
            "! io/controls_lib.f90 -- do not edit; regenerate instead.\n"
            "! The two copies between read_input's namelist read BUFFER\n"
            "! (controls_lib's module variables) and the authoritative\n"
            "! star%ctrl. read_controls drives the sequence:\n"
            "!   star%ctrl = controls_state()   (pristine defaults)\n"
            "!   call seed_controls_buffer      (buffer <- star%ctrl)\n"
            "!   call read_input(...)               (namelist reads over buffer)\n"
            "!   call store_controls_to_star    (star%ctrl <- buffer)\n"
            "module controls_sync_lib\n"
            "      use controls_lib\n"
            "      use star_info_lib, only: star\n"
            "      implicit none\n"
            "\n"
            "contains\n"
            "\n"
            "subroutine seed_controls_buffer\n"
        )
        for _, name, _, _ in members:
            home = "job" if name in JOB_MEMBERS else "ctrl"
            f.write(f"      {name} = star%{home}%{name}\n")
        f.write(
            "end subroutine seed_controls_buffer\n"
            "\n"
            "subroutine store_controls_to_star\n"
        )
        for _, name, _, _ in members:
            home = "job" if name in JOB_MEMBERS else "ctrl"
            f.write(f"      star%{home}%{name} = {name}\n")
        f.write(
            "end subroutine store_controls_to_star\n"
            "\n"
            "end module controls_sync_lib\n"
        )
    print(f"{len(members)} members -> {DEF.name}, {SYNC.name}")


if __name__ == "__main__":
    main()

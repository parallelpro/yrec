"""Controls regression test: every deck's controls, exactly as read_controls
produces them.

Runs src/test_controls (built by `make tests` in src/) on the pristine
defaults and on every deck in examples/ and testsuite/: the single-file
inlists and the legacy CONTROL/PHYSICS (.nml1/.nml2) pairs. Each deck runs
from its own directory (some name inputs relative to it) with YREC_OUTPUT a
scratch directory. Each dump is normalized: the YREC_OUTPUT, YREC_INPUT and
YREC_START directories become placeholders and NUL bytes are escaped. A deck
whose starting model is another run's output (the testsuite's HeIgnite,
ZAHB and TAHB chains) stops inside read_controls; its hash records that. Then it compares
  * the pristine dump with src/io/test/expected_controls/pristine.txt;
  * for the representative decks (the legacy pairs and the configuration
    matrix), the lines that differ from pristine with
    src/io/test/expected_controls/<deck>.txt;
  * a hash of every deck's dump with src/io/test/expected_controls/hashes.txt.
Regenerate the expected files with YREC_REGEN_CONTROLS=1 after a deliberate
change.
"""
import concurrent.futures
import hashlib
import os
import pathlib
import re
import subprocess
import tempfile

import pytest

REPO = pathlib.Path(__file__).parent
BINARY = REPO / "src" / "test_controls"
EXPECTED = REPO / "src" / "io" / "test" / "expected_controls"
REGEN = os.environ.get("YREC_REGEN_CONTROLS") == "1"


def decks():
    """(name, [files]) for every deck: inlists, and .nml1 with its .nml2"""
    out = []
    for top in ("examples", "testsuite"):
        for p in sorted((REPO / top).rglob("*.inlist")):
            out.append((str(p.relative_to(REPO)), [p]))
        for p in sorted((REPO / top).rglob("*.nml1")):
            q = p.with_suffix(".nml2")
            if q.exists():
                out.append((str(p.relative_to(REPO)) + "+nml2", [p, q]))
    return out


def representative(name):
    return name.startswith("examples/run_config_matrix/") or name.endswith("+nml2")


def run(files):
    """the normalized dump; decks run from their own directory (some name
    inputs relative to it), with their output in a scratch directory"""
    env = dict(os.environ)
    with tempfile.TemporaryDirectory() as tmp:
        env["YREC_OUTPUT"] = tmp
        cwd = files[0].parent if files else tmp
        result = subprocess.run([str(BINARY)] + [str(f) for f in files], cwd=cwd, env=env,
                                capture_output=True)
    text = result.stdout.decode("latin-1").replace("\x00", "\\0")
    for var in ("YREC_OUTPUT", "YREC_INPUT", "YREC_START"):
        if env.get(var):
            text = text.replace(env[var], "{" + var + "}")
    lines = text.splitlines()
    marker = next((i for i, l in enumerate(lines) if l.startswith("# controls")), None)
    if marker is None:
        return [f"status = {result.returncode}", "no dump"]
    return lines[marker:] + ([f"status = {result.returncode}"] if result.returncode else [])


def differing(dump, pristine):
    base = dict(l.split(" = ", 1) for l in pristine[1:] if " = " in l)
    return [l for l in dump if " = " not in l or base.get(l.split(" = ", 1)[0]) != l.split(" = ", 1)[1]]


def expect(path, lines):
    text = "\n".join(lines) + "\n"
    if REGEN:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
        return True
    return path.exists() and path.read_text() == text


def test_controls_every_deck():
    if not BINARY.exists():
        pytest.skip("src/test_controls not built (run `make tests` in src/)")
    pristine = run([])
    failures = []
    if not expect(EXPECTED / "pristine.txt", pristine):
        failures.append("pristine")
    all_decks = decks()
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as pool:
        dumps = list(pool.map(lambda d: run(d[1]), all_decks))
    hashes = []
    for (name, _), dump in zip(all_decks, dumps):
        hashes.append(f"{hashlib.sha256(chr(10).join(dump).encode()).hexdigest()}  {name}")
        if representative(name):
            fname = name.replace("/", "__") + ".txt"
            if not expect(EXPECTED / fname, differing(dump, pristine)):
                failures.append(name)
    if not expect(EXPECTED / "hashes.txt", hashes):
        old = dict(l.split("  ", 1)[::-1] for l in (EXPECTED / "hashes.txt").read_text().splitlines()) \
            if (EXPECTED / "hashes.txt").exists() else {}
        for line in hashes:
            h, name = line.split("  ", 1)
            if old.get(name) != h:
                failures.append(name)
    assert not failures, f"controls differ for {len(failures)} deck(s): {failures[:10]}"


def test_zoning_tolerances_without_new_rate_variables():
    """use_new_rate_variables = .false.: the zoning tolerances are the deck's
    zoning_tolerances (legacy hpttol), at their documented defaults when the
    deck sets none, and the timestep factors its timestep_factors (atime).
    Every shipped deck sets the switch, so only this test reads that path."""
    if not BINARY.exists():
        pytest.skip("src/test_controls not built (run `make tests` in src/)")
    deck = REPO / "examples" / "run_config_matrix" / "cm_base.inlist"
    text = deck.read_text()
    assert "use_new_rate_variables = .TRUE." in text
    text = text.replace("use_new_rate_variables = .TRUE.", "use_new_rate_variables = .FALSE.")
    defaults = (REPO / "src" / "defaults" / "physics.defaults").read_text()

    def documented(name):
        line = next(l for l in defaults.splitlines() if l.strip().startswith(name + " = "))
        return [float(v.replace("d", "e")) for v in line.split(" = ", 1)[1].split(",")]

    def dumped(dump, name):
        vals = {}
        for l in dump:
            if l.startswith(f"ctrl%{name}("):
                vals[int(l.split("(")[1].split(")")[0])] = float(l.split(" = ")[1])
        return [vals[i] for i in sorted(vals)]

    with tempfile.TemporaryDirectory() as tmp:
        off = pathlib.Path(tmp) / "cm_base_lnewvars_off.inlist"
        off.write_text(text)
        env = dict(os.environ, YREC_OUTPUT=tmp)
        result = subprocess.run([str(BINARY), str(off)], cwd=deck.parent, env=env,
                                capture_output=True)
    dump = result.stdout.decode("latin-1").splitlines()
    assert "ierr = 0" in dump
    zoning = dumped(dump, "chi_grid_scale")
    assert zoning == documented("zoning_tolerances"), zoning
    assert all(v > 0 for i, v in enumerate(zoning) if i != 6), "a zero zoning tolerance"
    assert dumped(dump, "atime") == documented("timestep_factors")


def test_pristine_controls_are_the_defaults(tmp_path):
    """controls_state(), the pristine star%ctrl and star%job, is what a deck
    that sets nothing reads: each default is written once, as its
    controls_buffer initializer, and the legacy names take theirs from it
    (defaults/*.defaults document the same initializers). The exceptions are
    not defaults: the unit numbers the read assigns, the expanded
    inlist_used_file, and what the read derives (npenv's 0 meaning every
    point, the run kinds, the OPAL EOS table path)."""
    if not BINARY.exists():
        pytest.skip("src/test_controls not built (run `make tests` in src/)")
    derived = re.compile(r"^(ctrl%(\w*_unit\d?|unit_\w+|iopale|icondopacp|"
                         r"inlist_used_file|npenv)|job%(opal_eos_table_path|rescale_kind))(\(|$)")
    deck = tmp_path / "empty.inlist"
    deck.write_text("&control\n first_model_file = 'empty.inlist'\n/\n&physics\n/\n")
    env = dict(os.environ, YREC_OUTPUT=str(tmp_path))

    def dump(args):
        out = subprocess.run([str(BINARY)] + args, cwd=tmp_path, env=env,
                             capture_output=True).stdout.decode("latin-1").splitlines()
        k = next(i for i, l in enumerate(out) if l.startswith("# controls"))
        return dict(l.split(" = ", 1) for l in out[k + 1:] if " = " in l)

    pristine, empty = dump([]), dump([str(deck)])
    assert empty.get("ierr") == "0"
    diff = [k for k in pristine if k.startswith(("ctrl%", "job%"))
            and pristine[k] != empty.get(k) and not derived.match(k)]
    assert not diff, [(k, pristine[k], empty.get(k)) for k in diff[:10]]


def test_inert_controls_warn_on_the_terminal(tmp_path):
    """A deck that sets a control nothing uses gets a warning on the
    terminal (controls_check's warn_inert_controls), and the run log stays
    as it was: a giant deck sets use_rotation_scaled_solar_wind = .TRUE.
    (its wind is commented out), cm_base sets it .FALSE.; both set
    surface_composition_file, which nothing writes."""
    if not BINARY.exists():
        pytest.skip("src/test_controls not built (run `make tests` in src/)")
    giant = next((REPO / "examples" / "giant_differential_rotation").glob("*RGp2diffrotI0.inlist"))
    base = REPO / "examples" / "run_config_matrix" / "cm_base.inlist"
    out = {}
    for deck in (giant, base):
        env = dict(os.environ, YREC_OUTPUT=str(tmp_path / deck.stem))
        r = subprocess.run([str(BINARY), str(deck)], cwd=deck.parent, env=env, capture_output=True)
        out[deck] = r.stdout.decode("latin-1")
        log = "".join(f.read_text(errors="replace") for f in (tmp_path / deck.stem).glob("*")
                      if f.is_file() and f.suffix in (".log", ".short"))
        assert "WARNING (controls)" not in log
    wind = "WARNING (controls): use_rotation_scaled_solar_wind"
    assert wind in out[giant] and wind not in out[base]
    assert all("WARNING (controls): surface_composition_file" in o for o in out.values())


def test_card_one_starts_from_the_starting_model(tmp_path):
    """Card 1 has no previous card's model, so it always loads the starting
    model: a deck that sets first_call_flag(1) = .FALSE. (ten shipped decks
    did) reads job%first_call_flag(1) = T, and the run plan says so; card
    2's flag stays as the deck sets it."""
    if not BINARY.exists():
        pytest.skip("src/test_controls not built (run `make tests` in src/)")
    deck = REPO / "examples" / "run_config_matrix" / "cm_base.inlist"
    text = deck.read_text()
    assert "first_call_flag(1) = .TRUE." in text
    off = tmp_path / "cm_base_card1_false.inlist"
    off.write_text(text.replace("first_call_flag(1) = .TRUE.", "first_call_flag(1) = .FALSE."))
    env = dict(os.environ, YREC_OUTPUT=str(tmp_path))
    result = subprocess.run([str(BINARY), str(off)], cwd=deck.parent, env=env, capture_output=True)
    out = result.stdout.decode("latin-1")
    dump = out.splitlines()
    assert "ierr = 0" in dump
    assert "job%first_call_flag(1) = T" in dump
    assert "job%first_call_flag(2) = F" in dump
    assert "card  1: rescale and evolve the starting model" in out, out[-500:]

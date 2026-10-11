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

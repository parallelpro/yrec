#!/usr/bin/env python3
"""Convergence of a rotating run in its rotation sub-step.

  python3 src/tools/substep_convergence.py --deck PATH --work DIR
        [--factors 0.04,0.02,0.01,0.005,0.0025] [--yrec src/yrec] [--jobs 6]

The rotation step is cut into sub-steps that each change omega by at
most time_dw_mix (diffusion_timestep_factor); nothing else reads it. The
deck is run once per factor, with every time_dw_mix line rewritten,
each in its own directory under DIR (YREC_INPUT and YREC_START point at
the repository's). For the last model of each run the tool prints the
age, the surface omega, the total angular momentum and log L, and, for
each quantity q, the differences between successive factors, d_k =
q(f_k) - q(f_{k+1}), and the observed order p = log2(d_k/d_{k+1}) (the
factors halving). A first-order splitting gives p near 1. The runs
must end at the same age (a max_age stop), or the comparison mixes in
the age difference.
"""
import argparse
import concurrent.futures
import math
import os
import pathlib
import re
import shutil
import subprocess

REPO = pathlib.Path(__file__).resolve().parents[2]
COLUMNS = ("star_age", "surf_avg_omega", "total_angular_momentum", "log_L")


def last_model(history):
    lines = history.read_text().splitlines()
    names = lines[5].split()
    row = [l for l in lines[6:] if l.strip()][-1].split()
    return {c: float(row[names.index(c)]) for c in COLUMNS}


def run(deck_text, factor, work, yrec, stem):
    rundir = work / f"{stem}_f{factor:g}"
    if rundir.exists():
        shutil.rmtree(rundir)
    (rundir / "output").mkdir(parents=True)
    text, n = re.subn(r"(?im)^(\s*time_dw_mix\s*=\s*)[-+0-9.eEdD]+", rf"\g<1>{factor:.6e}", deck_text)
    if n == 0:
        raise SystemExit("the deck sets no time_dw_mix")
    deck = rundir / f"{stem}.inlist"
    deck.write_text(text)
    env = dict(os.environ, YREC_INPUT=str(REPO / "input"),
               YREC_START=str(REPO / "startmodels"), YREC_OUTPUT="output")
    r = subprocess.run([str(yrec), deck.name], cwd=rundir, env=env,
                       capture_output=True, text=True)
    if r.returncode != 0:
        raise SystemExit(f"factor {factor}: yrec exit {r.returncode}\n{r.stdout[-2000:]}")
    return last_model(rundir / "output" / "history.data")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--deck", required=True)
    ap.add_argument("--work", required=True)
    ap.add_argument("--factors", default="0.04,0.02,0.01,0.005,0.0025")
    ap.add_argument("--yrec", default=str(REPO / "src" / "yrec"))
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    factors = [float(f) for f in a.factors.split(",")]
    deck = pathlib.Path(a.deck).resolve()
    work = pathlib.Path(a.work).resolve()
    work.mkdir(parents=True, exist_ok=True)
    text = deck.read_text()
    with concurrent.futures.ThreadPoolExecutor(a.jobs) as pool:
        results = list(pool.map(lambda f: run(text, f, work, pathlib.Path(a.yrec).resolve(), deck.stem), factors))
    print(f"{'time_dw_mix':>12s} " + " ".join(f"{c:>24s}" for c in COLUMNS))
    for f, r in zip(factors, results):
        print(f"{f:12.4e} " + " ".join(f"{r[c]:24.16e}" for c in COLUMNS))
    for c in COLUMNS[1:]:
        d = [results[k][c] - results[k + 1][c] for k in range(len(results) - 1)]
        orders = [math.log2(abs(d[k] / d[k + 1])) if d[k + 1] != 0 and d[k] != 0 else float("nan")
                  for k in range(len(d) - 1)]
        rel = [x / abs(results[-1][c]) for x in d]
        print(f"{c}: successive differences (relative) " + " ".join(f"{x:.2e}" for x in rel))
        print(f"{'':{len(c)}s}  observed order " + " ".join(f"{p:.2f}" for p in orders))


if __name__ == "__main__":
    main()

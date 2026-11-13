#!/usr/bin/env python3
"""
run_all.py -- reproduce every figure and every number in the project.

    python run_all.py            # everything (~10 minutes)
    python run_all.py e1 e5      # selected experiments
    python run_all.py --tests    # only the correctness checks

Each experiment writes
    figures/figNN_*.pdf + .png     vector figure for the report, raster for slides
    results/*.txt                  the console output of the run
    results/*.json / *.csv         the numbers quoted in the report

Author: Aryan Bandyopadhyay (2411014) -- PHY341/745 DIY project.
"""

import os
import subprocess
import sys
import time

ROOT = os.path.dirname(os.path.abspath(__file__))

EXPERIMENTS = [
    ("e1", "experiments/e1_validation.py", "Convergence, cost and stability"),
    ("e2", "experiments/e2_shock_formation.py", "Shock formation and thickness"),
    ("e3", "experiments/e3_spectra.py", "Energy spectra and the k^-2 law"),
    ("e4", "experiments/e4_structure_functions.py", "Structure functions and exponents"),
    ("e5", "experiments/e5_grid_viscosity.py", "Numerical viscosity of the grid"),
]


def run(cmd, label):
    print("=" * 78)
    print(f"  {label}")
    print("=" * 78)
    t0 = time.time()
    result = subprocess.run([sys.executable] + cmd, cwd=ROOT)
    dt = time.time() - t0
    print(f"  -> {label}: {'ok' if result.returncode == 0 else 'FAILED'} in {dt:.1f} s\n")
    return result.returncode, dt


def main(argv):
    if "--tests" in argv:
        sys.exit(run(["tests/test_core.py"], "correctness checks")[0])

    wanted = [a for a in argv if not a.startswith("-")] or [k for k, _, _ in EXPERIMENTS]
    total, failed = 0.0, []
    for key, path, label in EXPERIMENTS:
        if key not in wanted:
            continue
        code, dt = run([path], f"{key.upper()} -- {label}")
        total += dt
        if code != 0:
            failed.append(key)

    print("=" * 78)
    print(f"  total wall time: {total/60:.1f} min")
    if failed:
        print(f"  FAILED: {', '.join(failed)}")
        sys.exit(1)
    print("  all experiments completed; see figures/ and results/")
    print("=" * 78)


if __name__ == "__main__":
    main(sys.argv[1:])

"""
Run the three experiments in order and time them.

    python run_all.py              # everything (about 45 s)
    python run_all.py e3           # just one
    python run_all.py --tests      # the eleven checks instead
"""

import os
import subprocess
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable or "python3"

EXPERIMENTS = {
    "e1": ("experiments/e1_sanity.py", "Sanity checks: diffusion, order, stability"),
    "e2": ("experiments/e2_shock.py", "The shock: profile, strain, width"),
    "e3": ("experiments/e3_spectrum.py", "Fourier spectrum of the shock"),
}


def run(script, label):
    t0 = time.time()
    print(f"\n{'=' * 78}\n  {label}\n  {script}\n{'=' * 78}")
    result = subprocess.run([PY, os.path.join(HERE, script)], cwd=HERE)
    status = 'ok' if result.returncode == 0 else 'FAILED'
    print(f"\n  -> {status} in {time.time() - t0:.1f} s")
    return result.returncode


def main():
    args = [a for a in sys.argv[1:]]
    if "--tests" in args:
        run("tests/test_burgers.py", "Tests")
        return

    chosen = [a for a in args if a in EXPERIMENTS] or list(EXPERIMENTS)
    t0 = time.time()
    for key in chosen:
        script, label = EXPERIMENTS[key]
        if run(script, label) != 0:
            sys.exit(1)
    print(f"\n{'=' * 78}")
    print(f"  all done in {(time.time() - t0) / 60:.1f} min -- see figures/ and results/")
    print(f"{'=' * 78}")


if __name__ == "__main__":
    main()

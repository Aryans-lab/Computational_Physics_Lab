#!/usr/bin/env python3
"""
run_all.py
----------
Regenerate every assignment: bundle the library, then run each
``assignments/*/question*.py`` with its default paths (inputs from ``data/``,
results into ``output/`` and ``figures/`` next to the script).

Any script is runnable from anywhere with no arguments — data paths are
resolved relative to the script itself.  This is also the smoke test used in
CI, and the way to reproduce the committed ``output/`` files and plots.

Usage:
    python scripts/run_all.py
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Headless plot generation, wherever this runs (CI, laptop, lab machine).
ENV = {**os.environ, "MPLBACKEND": "Agg", "MATPLOTLIBRC": "/dev/null"}


def main() -> int:
    # 1. Make sure every folder has a current copy of the library.
    proc = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "bundle.py")], env=ENV
    )
    if proc.returncode != 0:
        return proc.returncode

    # 2. Run every question driver with default paths.
    scripts = sorted((ROOT / "assignments").glob("*/question*.py"))
    if not scripts:
        print("error: no question scripts found", file=sys.stderr)
        return 1

    print(f"\nRunning {len(scripts)} scripts...\n")
    failures: list[tuple[Path, subprocess.CompletedProcess]] = []
    for script in scripts:
        started = time.perf_counter()
        proc = subprocess.run(
            [sys.executable, str(script)],
            cwd=ROOT,
            env=ENV,
            capture_output=True,
            text=True,
        )
        elapsed = time.perf_counter() - started
        rel = script.relative_to(ROOT)
        if proc.returncode == 0:
            print(f"  [ok]   {rel}  ({elapsed:5.2f}s)")
        else:
            print(f"  [FAIL] {rel}  ({elapsed:5.2f}s)")
            failures.append((script, proc))

    if failures:
        print("\n" + "=" * 72)
        print(f"{len(failures)} script(s) failed:\n")
        for script, proc in failures:
            print(f"--- {script.relative_to(ROOT)} ---")
            print(proc.stdout[-2000:])
            print(proc.stderr[-2000:])
        return 1

    total = len(scripts)
    print(f"\nAll {total} scripts finished successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

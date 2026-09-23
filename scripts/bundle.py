#!/usr/bin/env python3
"""
bundle.py
---------
Copy the canonical ``mylib.py`` (at the repository root) into every folder
under ``assignments/`` so that each assignment is a self-contained
submission unit (front code + library in the same directory, as required by
the course convention).

The per-assignment copies are *generated artifacts*: the single source of
truth is ``mylib.py`` at the root, and the copies are git-ignored.  Run this
script any time the library changes — it is idempotent, and it reports any
out-of-date copies it overwrites (so library drift can never go unnoticed).

Usage:
    python scripts/bundle.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
CANONICAL = ROOT / "mylib.py"
TARGETS = sorted((ROOT / "assignments").glob("*/"))


def main() -> int:
    if not CANONICAL.is_file():
        print(f"error: canonical library not found at {CANONICAL}", file=sys.stderr)
        return 1
    if not TARGETS:
        print(f"error: no assignment folders found under {ROOT / 'assignments'}",
              file=sys.stderr)
        return 1

    src = CANONICAL.read_bytes()
    stale = 0
    for folder in TARGETS:
        dst = folder / "mylib.py"
        if dst.is_file() and dst.read_bytes() != src:
            stale += 1
            print(f"  out of date -> re-bundled {dst.relative_to(ROOT)}")
        elif not dst.is_file():
            print(f"  bundled       -> {dst.relative_to(ROOT)}")
        dst.write_bytes(src)

    print(f"Bundle complete: {len(TARGETS)} folders, {stale} out-of-date copy(ies) "
          f"overwritten. Source: mylib.py ({len(src)} bytes).")
    return 0


if __name__ == "__main__":
    sys.exit(main())

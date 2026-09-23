"""Shared pytest fixtures: put the repository root on sys.path so that
``import mylib`` resolves to the canonical library regardless of where
pytest is invoked from."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

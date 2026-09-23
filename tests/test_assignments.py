"""Golden-output regression tests for every assignment driver.

Each ``assignments/*/questionN.py`` is executed inside a sandboxed *copy* of
the assignments tree (with the canonical library bundled in), using its
default paths, and the regenerated ``output/qN_output.txt`` must be
byte-identical to the committed one.  This guarantees that the numbers in
the repository's output files are exactly what the current code produces.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
ENV = {**os.environ, "MPLBACKEND": "Agg", "MATPLOTLIBRC": "/dev/null"}


def question_scripts() -> list[Path]:
    return sorted((ROOT / "assignments").glob("*/question*.py"))


def sandbox_script(sandbox: Path, script: Path) -> Path:
    """The sandboxed copy of an assignment script (paths relative to assignments/)."""
    return sandbox / script.relative_to(ROOT / "assignments")


@pytest.fixture(scope="module")
def sandbox(tmp_path_factory) -> Path:
    """A throwaway copy of the assignments tree with the library bundled."""
    target = tmp_path_factory.mktemp("assignments")
    for folder in sorted((ROOT / "assignments").iterdir()):
        if folder.is_dir():
            shutil.copytree(folder, target / folder.name)
    for folder in sorted(target.iterdir()):
        shutil.copyfile(ROOT / "mylib.py", folder / "mylib.py")
    return target


def test_every_question_has_a_committed_output():
    for script in question_scripts():
        m = re.match(r"question(\d+)\.py", script.name)
        assert m, f"unexpected script name {script.name}"
        committed = script.parent / "output" / f"q{m.group(1)}_output.txt"
        assert committed.is_file(), f"missing committed output: {committed}"


@pytest.mark.parametrize(
    "script", question_scripts(), ids=lambda p: str(p.relative_to(ROOT))
)
def test_golden_output(sandbox: Path, script: Path) -> None:
    m = re.match(r"question(\d+)\.py", script.name)
    name = script.name
    sb = sandbox_script(sandbox, script)
    assert sb.is_file()

    proc = subprocess.run(
        [sys.executable, str(sb)],
        cwd=sandbox,
        env=ENV,
        capture_output=True,
        text=True,
        timeout=600,
    )
    assert proc.returncode == 0, (
        f"{name} failed:\nSTDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}"
    )

    produced = sb.parent / "output" / f"q{m.group(1)}_output.txt"
    committed = script.parent / "output" / f"q{m.group(1)}_output.txt"
    assert produced.is_file(), f"{name} produced no output file"
    assert produced.read_text() == committed.read_text(), (
        f"{name}: regenerated output differs from the committed golden file"
    )


def test_repository_outputs_are_self_consistent(sandbox: Path) -> None:
    """Re-running each script twice must give identical outputs (determinism)."""
    for script in question_scripts():
        m = re.match(r"question(\d+)\.py", script.name)
        sb = sandbox_script(sandbox, script)
        run1 = subprocess.run(
            [sys.executable, str(sb)],
            cwd=sandbox, env=ENV, capture_output=True, text=True, timeout=600,
        )
        produced = sb.parent / "output" / f"q{m.group(1)}_output.txt"
        first = produced.read_text()
        run2 = subprocess.run(
            [sys.executable, str(sb)],
            cwd=sandbox, env=ENV, capture_output=True, text=True, timeout=600,
        )
        assert run1.returncode == 0 and run2.returncode == 0
        assert produced.read_text() == first, f"{script.name} is not deterministic"

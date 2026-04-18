from __future__ import annotations

import runpy
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_DIR = ROOT / "examples"

pytestmark = pytest.mark.integration


@pytest.mark.parametrize(
    "script_name",
    [
        "dataset_builder_demo.py",
        "simple_workflow_demo.py",
        "signalome_workflow_demo.py",
    ],
)
def test_public_example_script_runs(script_name: str) -> None:
    runpy.run_path(str(EXAMPLES_DIR / script_name), run_name="__main__")

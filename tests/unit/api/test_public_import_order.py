from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
SRC = ROOT / "src"


def _run_fresh(source: str) -> None:
    environment = dict(os.environ)
    environment["PYTHONPATH"] = os.pathsep.join(
        (str(SRC), environment.get("PYTHONPATH", ""))
    )
    completed = subprocess.run(
        [sys.executable, "-c", source],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr


def test_ordinary_then_specialist_imports_preserve_owner_identities() -> None:
    _run_fresh(
        """
import phospy.api as ordinary
import phospy.api.advanced as advanced
from phospy.datasets.builder import AnalysisReadyDatasetBuilder
from phospy.signalome.config import SignalomeConfig

assert len(ordinary.__all__) == 66
assert len(advanced.__all__) == 126
assert ordinary.AnalysisReadyDatasetBuilder is AnalysisReadyDatasetBuilder
assert advanced.SignalomeConfig is SignalomeConfig
assert set(ordinary.__all__).isdisjoint(set(advanced.__all__))
"""
    )


def test_direct_specialist_import_initialises_parent_without_broken_identity() -> None:
    _run_fresh(
        """
import sys
import phospy.api.advanced as advanced
from phospy.api import KinaseWorkflow
from phospy.kinase.public import KinaseWorkflow as OwnerKinaseWorkflow
from phospy.kinase.scoring.config import KinaseScoringConfig

assert "phospy.api" in sys.modules
assert len(advanced.__all__) == 126
assert KinaseWorkflow is OwnerKinaseWorkflow
assert advanced.KinaseScoringConfig is KinaseScoringConfig
"""
    )


def test_root_version_import_is_lightweight_and_metadata_only() -> None:
    _run_fresh(
        """
import sys
import phospy

assert phospy.__version__
assert phospy.__all__ == []
assert not any(name.startswith("phospy.") for name in sys.modules)
for name in ("AnalysisReadyDatasetBuilder", "KinaseWorkflow", "SignalomeWorkflow"):
    assert not hasattr(phospy, name)
"""
    )


def test_retired_facades_and_root_conveniences_fail_in_fresh_interpreter() -> None:
    _run_fresh(
        """
import importlib

for module_name in ("phospy.advanced", "phospy._api_inventory"):
    try:
        importlib.import_module(module_name)
    except ModuleNotFoundError:
        pass
    else:
        raise AssertionError(f"retired module still resolves: {module_name}")

for statement in (
    "from phospy import AnalysisReadyDatasetBuilder",
    "from phospy import KinaseWorkflow",
    "from phospy import SignalomeWorkflow",
):
    try:
        exec(statement, {})
    except ImportError:
        pass
    else:
        raise AssertionError(f"retired root convenience still resolves: {statement}")
"""
    )

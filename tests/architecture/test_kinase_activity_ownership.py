from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from phospy.activities.contracts import MethodQuantitativeInputContract
from phospy.activities.semantics import ActivityInputMatrix
from phospy.advanced import KinaseActivityConfig as FacadeKinaseActivityConfig
from phospy.api import KinaseActivityResult as FacadeKinaseActivityResult
from phospy.kinase.activity.config import KinaseActivityConfig
from phospy.kinase.activity.results import KinaseActivityResult

PROJECT_ROOT = Path(__file__).resolve().parents[2]
RETIRED_ACTIVITY_ROUTES = (
    "phospy.science.activities",
    "phospy.science.activities.membership",
    "phospy.science.activities.methods.ksea_zscore",
    "phospy.science.activities.methods.ssgsea_substrate_enrichment",
    "phospy.science.activities.methods.weighted_substrate_activity",
    "phospy.science.activities.semantics",
    "phospy.science.tables.activity",
    "phospy.validation.workflows.activity",
    "phospy.workflows.kinase.activity_runner",
    "phospy.workflows.kinase.membership",
    "phospy.contracts.configs.kinase",
    "phospy.science.configs.kinase",
)


def _fresh_process(source: str) -> subprocess.CompletedProcess[str]:
    environment = os.environ.copy()
    environment["PYTHONPATH"] = str(PROJECT_ROOT / "src")
    return subprocess.run(
        [sys.executable, "-c", source],
        cwd=PROJECT_ROOT,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
    )


@pytest.mark.parametrize("module_name", RETIRED_ACTIVITY_ROUTES)
def test_retired_activity_routes_fail_in_fresh_process(module_name: str) -> None:
    completed = _fresh_process(
        "import importlib\n"
        f"name = {module_name!r}\n"
        "try:\n"
        "    importlib.import_module(name)\n"
        "except ModuleNotFoundError as error:\n"
        "    if error.name is None or not name.startswith(error.name):\n"
        "        raise\n"
        "else:\n"
        "    raise AssertionError(f'retired route imported: {name}')\n"
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize("owner_first", (False, True))
def test_activity_facade_and_owner_import_orders_preserve_identity(
    owner_first: bool,
) -> None:
    imports = (
        "from phospy.kinase.activity.config import KinaseActivityConfig as OwnerConfig\n"
        "from phospy.kinase.activity.results import KinaseActivityResult as OwnerResult\n"
        "from phospy.advanced import KinaseActivityConfig as PublicConfig\n"
        "from phospy.api import KinaseActivityResult as PublicResult\n"
        if owner_first
        else "from phospy.advanced import KinaseActivityConfig as PublicConfig\n"
        "from phospy.api import KinaseActivityResult as PublicResult\n"
        "from phospy.kinase.activity.config import KinaseActivityConfig as OwnerConfig\n"
        "from phospy.kinase.activity.results import KinaseActivityResult as OwnerResult\n"
    )
    completed = _fresh_process(
        imports
        + "assert OwnerConfig is PublicConfig\n"
        + "assert OwnerResult is PublicResult\n"
    )
    assert completed.returncode == 0, completed.stderr


def test_supported_activity_facades_share_the_concrete_owner_definitions() -> None:
    assert FacadeKinaseActivityConfig is KinaseActivityConfig
    assert FacadeKinaseActivityResult is KinaseActivityResult
    assert KinaseActivityConfig.__module__ == "phospy.kinase.activity.config"
    assert KinaseActivityResult.__module__ == "phospy.kinase.activity.results"
    assert ActivityInputMatrix.__module__ == "phospy.activities.semantics"
    assert MethodQuantitativeInputContract.__module__ == "phospy.activities.contracts"

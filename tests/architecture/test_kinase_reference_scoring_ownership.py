from __future__ import annotations

import ast
import importlib
import os
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

import phospy.advanced as advanced
import phospy.api as api
import phospy.kinase.results as kinase_result_contracts
from phospy.errors.validation import PhosPyValidationError
from phospy.kinase.contracts import ResolvedKinaseWorkflowRequest
from phospy.kinase.prediction.results import KinasePredictionResult
from phospy.kinase.references.models import (
    KinaseLibraryResource,
    KinaseLibraryResourceLoadRequest,
)
from phospy.kinase.scoring.config import KinaseAttritionPolicy, KinaseScoringConfig
from phospy.kinase.scoring.policies import ProfileSelfInclusionPolicy
from phospy.kinase.scoring.results import KinaseScoringResult

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"

_RETIRED_ROUTES = (
    "phospy.contracts.results.kinase",
    "phospy.contracts.results.kinase_boundary",
    "phospy.contracts.configs.prediction",
    "phospy.contracts.kinase_reference_projection",
    "phospy.science.configs.prediction",
    "phospy.science.prediction.motif_scoring",
    "phospy.science.prediction.scoring",
    "phospy.science.prediction",
    "phospy.science.prediction.candidates",
    "phospy.science.prediction.execution",
    "phospy.science.prediction.internal_view",
    "phospy.science.prediction.models",
    "phospy.science.prediction.policies",
    "phospy.science.prediction.sampling_core",
    "phospy.science.prediction.sampling_runtime",
    "phospy.science.prediction.scientific_policies",
    "phospy.science.prediction.sequence_validation",
    "phospy.science.prediction.svm",
    "phospy.science.scoring",
    "phospy.science.scoring.policy_models",
    "phospy.science.references.kinase_library",
    "phospy.science.references.kinase_library_loading",
    "phospy.science.references.kinase_library_models",
    "phospy.science.references.validation.kinase_library",
    "phospy.science.tables.kinase",
    "phospy.validation.references.kinase_library",
    "phospy.workflows.kinase.contributions",
    "phospy.workflows.kinase.kinase_library_scoring",
    "phospy.workflows.kinase.component_models",
    "phospy.workflows.kinase.prediction_runner",
    "phospy.workflows.kinase.reference_projection",
    "phospy.workflows.kinase.scoring_mode_contracts",
    "phospy.workflows.kinase.scoring_runner",
    "phospy.workflows.kinase.scoring_transforms",
    "phospy.workflows.kinase.science",
    "phospy.workflows.kinase",
)


def _absolute_imports(path: Path) -> set[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    imports: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imports.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            imports.add(node.module)
    return imports


def test_supported_facades_preserve_the_single_new_owner_identities() -> None:
    assert KinaseScoringConfig.__module__ == "phospy.kinase.scoring.config"
    assert KinaseAttritionPolicy.__module__ == "phospy.kinase.scoring.config"
    assert ProfileSelfInclusionPolicy.__module__ == "phospy.kinase.scoring.policies"
    assert KinaseScoringResult.__module__ == "phospy.kinase.scoring.results"
    assert KinaseLibraryResource.__module__ == "phospy.kinase.references.models"
    assert KinaseLibraryResourceLoadRequest.__module__ == (
        "phospy.kinase.references.models"
    )

    assert advanced.KinaseScoringConfig is KinaseScoringConfig
    assert advanced.KinaseAttritionPolicy is KinaseAttritionPolicy
    assert advanced.ProfileSelfInclusionPolicy is ProfileSelfInclusionPolicy
    assert advanced.KinaseLibraryResource is KinaseLibraryResource
    assert advanced.KinaseLibraryResourceLoadRequest is KinaseLibraryResourceLoadRequest
    assert api.KinaseScoringResult is KinaseScoringResult


def test_prediction_and_workflow_lifecycle_objects_have_approved_owners() -> None:
    assert KinasePredictionResult.__module__ == "phospy.kinase.prediction.results"
    assert api.KinasePredictionResult is KinasePredictionResult
    assert kinase_result_contracts.KinasePredictionResult is KinasePredictionResult
    assert ResolvedKinaseWorkflowRequest.__module__ == "phospy.kinase.contracts"
    assert hasattr(
        importlib.import_module("phospy.kinase.activity.results"),
        "KinaseActivityResult",
    )
    assert hasattr(importlib.import_module("phospy.kinase.public"), "KinaseWorkflow")


def test_retired_prediction_package_is_absent() -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("phospy.science.prediction")


def test_prediction_owned_frame_guards_preserve_former_behavior() -> None:
    site_key = (
        "phospy:v1|organism=rat|protein_namespace=protein_id|"
        "protein_identifier=P1|residue=S|position=1"
    )
    pred_mat = pd.DataFrame({"K1": [0.5]}, index=[site_key])
    substrate_list = pd.DataFrame(
        {"kinase": ["K1"], "substrate_site": [site_key], "score": [0.5], "rank": [1]}
    )

    with pytest.raises(
        PhosPyValidationError,
        match=r"^prediction_result\.pred_mat must be a pandas DataFrame$",
    ):
        KinasePredictionResult.from_trusted_owned(
            pred_mat="not-a-frame",  # type: ignore[arg-type]
            substrate_list="also-not-a-frame",  # type: ignore[arg-type]
        )
    with pytest.raises(
        PhosPyValidationError,
        match=r"^prediction_result\.substrate_list must be a pandas DataFrame when provided$",
    ):
        KinasePredictionResult.from_trusted_owned(
            pred_mat=pred_mat,
            substrate_list="not-a-frame",  # type: ignore[arg-type]
        )

    trusted = KinasePredictionResult.from_trusted_owned(
        pred_mat=pred_mat,
        substrate_list=substrate_list,
    )
    assert trusted._pred_mat is pred_mat
    assert trusted._substrate_list is substrate_list
    assert trusted.pred_mat is not pred_mat
    assert trusted.pred_mat.equals(pred_mat)
    assert trusted.substrate_list is not substrate_list
    trusted_substrate_list = trusted.substrate_list
    assert trusted_substrate_list is not None
    assert trusted_substrate_list.equals(substrate_list)

    optional = KinasePredictionResult.from_trusted_owned(
        pred_mat=pred_mat,
        substrate_list=None,
    )
    assert optional.substrate_list is None

    ordinary = KinasePredictionResult(pred_mat, substrate_list)
    assert ordinary._pred_mat is not pred_mat
    assert ordinary._substrate_list is not substrate_list


@pytest.mark.parametrize(
    "routes",
    [
        (
            "phospy.kinase.prediction.results",
            "phospy.kinase.results",
            "phospy.api",
        ),
        (
            "phospy.api",
            "phospy.kinase.results",
            "phospy.kinase.prediction.results",
        ),
        (
            "phospy.api",
            "phospy.kinase.prediction.results",
            "phospy.kinase.results",
        ),
    ],
)
def test_prediction_result_supported_routes_are_import_order_independent(
    routes: tuple[str, ...],
) -> None:
    code = "\n".join(
        (
            "import importlib",
            f"routes = {routes!r}",
            "modules = {route: importlib.import_module(route) for route in routes}",
            "owner_module = importlib.import_module('phospy.kinase.prediction.results')",
            "owner = owner_module.KinasePredictionResult",
            "assert importlib.import_module('phospy.api').KinasePredictionResult is owner",
            "kinase = importlib.import_module('phospy.kinase.results')",
            "assert kinase.KinasePredictionResult is owner",
            "package = importlib.import_module('phospy.kinase.prediction')",
            "assert '__getattr__' not in vars(package)",
        )
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(SRC_ROOT)
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_every_retired_route_fails_in_a_fresh_process() -> None:
    code = "\n".join(
        (
            "import importlib",
            f"routes = {_RETIRED_ROUTES!r}",
            "for route in routes:",
            "    try:",
            "        importlib.import_module(route)",
            "    except ModuleNotFoundError as error:",
            "        assert error.name == route or route.startswith(error.name + '.')",
            "    else:",
            "        raise AssertionError(route)",
        )
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(SRC_ROOT)
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


@pytest.mark.parametrize("order", ["supported-first", "retired-first", "alternating"])
def test_supported_and_retired_routes_are_import_order_independent(order: str) -> None:
    code = "\n".join(
        (
            "import importlib",
            f"order = {order!r}",
            f"retired = {_RETIRED_ROUTES!r}",
            "supported = ('phospy.kinase.prediction.results', "
            "'phospy.kinase.prediction.tables', 'phospy.api')",
            "def load_supported():",
            "    modules = [importlib.import_module(route) for route in supported]",
            "    owner = modules[0].KinasePredictionResult",
            "    assert modules[2].KinasePredictionResult is owner",
            "def reject(route):",
            "    try:",
            "        importlib.import_module(route)",
            "    except ModuleNotFoundError as error:",
            "        assert error.name == route or route.startswith(error.name + '.')",
            "    else:",
            "        raise AssertionError(route)",
            "if order == 'supported-first':",
            "    load_supported()",
            "    for route in retired: reject(route)",
            "elif order == 'retired-first':",
            "    for route in retired: reject(route)",
            "    load_supported()",
            "else:",
            "    for route in retired:",
            "        load_supported()",
            "        reject(route)",
        )
    )
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(SRC_ROOT)
    completed = subprocess.run(
        [sys.executable, "-c", code],
        cwd=PROJECT_ROOT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stderr


def test_general_references_do_not_import_kinase_implementation() -> None:
    offenders = sorted(
        f"{path.relative_to(PROJECT_ROOT).as_posix()} -> {target}"
        for path in (SRC_ROOT / "phospy" / "references").rglob("*.py")
        for target in _absolute_imports(path)
        if target == "phospy.kinase" or target.startswith("phospy.kinase.")
    )
    assert offenders == []


def test_signalome_consumes_only_published_kinase_contracts() -> None:
    allowed = {
        "phospy.kinase.result_boundary",
        "phospy.kinase.results",
        "phospy.kinase.prediction.tables",
        "phospy.kinase.scoring.results",
    }
    offenders = sorted(
        f"{path.relative_to(PROJECT_ROOT).as_posix()} -> {target}"
        for root in (SRC_ROOT / "phospy" / "signalome",)
        for path in root.rglob("*.py")
        for target in _absolute_imports(path)
        if (
            target == "phospy.kinase"
            or target.startswith("phospy.kinase.")
            or target.startswith("phospy.workflows.kinase")
        )
        if target not in allowed
    )
    assert offenders == []


def test_scoring_reference_definitions_remain_references_owned_once() -> None:
    references = importlib.import_module("phospy.references.tables")
    assert references.KinaseSubstrateReference.__module__ == "phospy.references.tables"
    source = (SRC_ROOT / "phospy" / "references" / "tables.py").read_text(
        encoding="utf-8"
    )
    tree = ast.parse(source)
    names = [
        node.name
        for node in ast.walk(tree)
        if isinstance(node, (ast.ClassDef, ast.FunctionDef))
    ]
    for name in (
        "KinaseSubstrateReference",
        "_classify_duplicate_and_conflicting_pair_records",
        "_build_pair_classification_records",
    ):
        assert names.count(name) == 1

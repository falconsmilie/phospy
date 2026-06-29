from __future__ import annotations

import ast
import inspect
from pathlib import Path

import phospy.api as public_api
from phospy.api.requests import (
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    KinaseWorkflowRequest,
    PhosphositeImportRequest,
    SignalomeWorkflowRequest,
)
from phospy.workflows.differential.validator import DifferentialAnalysisValidator
from phospy.workflows.kinase.validator import KinaseWorkflowValidator
from phospy.workflows.signalome.validator import SignalomeWorkflowValidator

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SRC_ROOT = PROJECT_ROOT / "src"
API_ROOT = SRC_ROOT / "phospy" / "api"
DOCS_ROOT = PROJECT_ROOT / "docs"
REQUEST_CONTRACTS_PATH = SRC_ROOT / "phospy" / "contracts" / "requests.py"
DATASET_VALIDATION_PREFIX = "phospy.validation.datasets"

DATASET_VALIDATOR_EXPORT_NAMES = frozenset(
    {
        "AnalysisReadyDatasetModelBoundaryValidator",
        "BatchCorrectionAdequacyValidator",
        "DISPLAY_SITE_CONTEXT_COLUMNS",
        "DatasetInputSourceValidator",
        "DatasetPreprocessingConfigValidator",
        "PhosphositeImportRequestValidator",
        "enforce_display_id_column",
        "enforce_site_key_column",
        "enforce_site_key_index",
        "enforce_site_key_matches_metadata",
        "enforce_unique_display_site_identity_rows",
        "enforce_unique_site_key_identity",
    }
)

DATASET_REQUEST_DTOS = (
    DatasetBuildRequest,
    DifferentialAnalysisRequest,
    KinaseWorkflowRequest,
    PhosphositeImportRequest,
    SignalomeWorkflowRequest,
)


def _imported_modules(path: Path) -> list[str]:
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    modules: list[str] = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            modules.extend(alias.name for alias in node.names)
        if isinstance(node, ast.ImportFrom) and node.module is not None:
            modules.append(node.module)
    return modules


def test_api_package_does_not_import_private_dataset_validation_modules() -> None:
    offenders = sorted(
        str(path.relative_to(PROJECT_ROOT))
        for path in API_ROOT.rglob("*.py")
        if any(
            module == DATASET_VALIDATION_PREFIX
            or module.startswith(f"{DATASET_VALIDATION_PREFIX}.")
            for module in _imported_modules(path)
        )
    )

    assert offenders == []


def test_api_all_has_no_dataset_validator_exports() -> None:
    explicit_name_leaks = sorted(
        DATASET_VALIDATOR_EXPORT_NAMES.intersection(public_api.__all__)
    )
    validator_name_leaks = sorted(
        name for name in public_api.__all__ if name.endswith("Validator")
    )
    dataset_module_leaks = sorted(
        name
        for name in public_api.__all__
        if getattr(getattr(public_api, name), "__module__", "").startswith(
            DATASET_VALIDATION_PREFIX
        )
    )

    assert explicit_name_leaks == []
    assert validator_name_leaks == []
    assert dataset_module_leaks == []


def test_workflow_validators_compose_shared_and_domain_validation() -> None:
    differential_source = inspect.getsource(DifferentialAnalysisValidator.run)
    kinase_source = inspect.getsource(KinaseWorkflowValidator.run)
    signalome_source = inspect.getsource(
        SignalomeWorkflowValidator._require_site_identity_and_protein_grouping_metadata
    )
    signalome_grouping_source = inspect.getsource(
        SignalomeWorkflowValidator._require_signalome_protein_grouping_metadata
    )

    for source in (differential_source, kinase_source, signalome_source):
        assert "enforce_workflow_site_identity_contract(" in source

    assert "self._dataset_eligibility_validator.run(" in differential_source
    assert "self._design_validator.run(" in differential_source
    assert "self._technical_replicate_planner.run(" in differential_source
    assert "enforce_localisation_requirement(" in kinase_source
    assert "enforce_localisation_requirement(" in signalome_source
    assert "enforce_required_non_empty_string_column(" in signalome_grouping_source


def test_dataset_request_dtos_do_not_call_private_dataset_validators() -> None:
    request_source = REQUEST_CONTRACTS_PATH.read_text(encoding="utf-8")

    for request_type in DATASET_REQUEST_DTOS:
        assert "__post_init__" not in request_type.__dict__

    assert DATASET_VALIDATION_PREFIX not in request_source
    for symbol_name in DATASET_VALIDATOR_EXPORT_NAMES:
        assert symbol_name not in request_source


def test_docs_do_not_present_dataset_validators_as_user_import_path() -> None:
    docs_with_private_imports = sorted(
        str(path.relative_to(PROJECT_ROOT))
        for path in DOCS_ROOT.rglob("*.md")
        if "from phospy.validation.datasets" in path.read_text(encoding="utf-8")
        or "import phospy.validation.datasets" in path.read_text(encoding="utf-8")
    )
    ownership_doc = (DOCS_ROOT / "validation-ownership.md").read_text(encoding="utf-8")

    assert docs_with_private_imports == []
    assert "Users validate dataset inputs by constructing datasets" in ownership_doc
    assert "not supported user entrypoints" in ownership_doc

from __future__ import annotations

import importlib

import pytest

import phospy
import phospy.api as public_api

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

PUBLIC_API_SMOKE_IMPORTS = frozenset(
    {
        "AnalysisReadyDatasetBuilder",
        "AnalysisReadyPhosphoDataset",
        "DatasetBuildRequest",
        "DifferentialAnalysisWorkflow",
        "KinaseWorkflow",
        "PhosPyValidationError",
        "ReferenceBundle",
        "SignalomeWorkflow",
    }
)


def _assert_from_import_fails(module_name: str, symbol_name: str) -> None:
    with pytest.raises(ImportError):
        exec(f"from {module_name} import {symbol_name}", {})


def test_phospy_api_does_not_export_dataset_validator_names() -> None:
    for symbol_name in DATASET_VALIDATOR_EXPORT_NAMES:
        assert symbol_name not in public_api.__all__
        assert not hasattr(public_api, symbol_name)
        _assert_from_import_fails("phospy.api", symbol_name)

        assert symbol_name not in phospy.__all__
        assert not hasattr(phospy, symbol_name)
        _assert_from_import_fails("phospy", symbol_name)


def test_phospy_api_all_contains_no_validator_implementations() -> None:
    validator_names = sorted(
        name for name in public_api.__all__ if name.endswith("Validator")
    )
    dataset_validation_exports = sorted(
        name
        for name in public_api.__all__
        if getattr(getattr(public_api, name), "__module__", "").startswith(
            "phospy.validation.datasets"
        )
    )

    assert validator_names == []
    assert dataset_validation_exports == []


def test_retired_dataset_validation_package_is_absent() -> None:
    with pytest.raises(ModuleNotFoundError) as exc_info:
        importlib.import_module("phospy.validation.datasets")
    assert exc_info.value.name == "phospy.validation.datasets"


def test_internal_validation_and_evidence_owner_imports_still_work() -> None:
    from phospy.datasets.validation.model_boundary import (
        AnalysisReadyDatasetModelBoundaryValidator,
    )
    from phospy.evidence import site_metadata

    assert site_metadata.__name__ == "phospy.evidence.site_metadata"
    assert AnalysisReadyDatasetModelBoundaryValidator.__name__ == (
        "AnalysisReadyDatasetModelBoundaryValidator"
    )


def test_public_api_import_smoke_still_passes() -> None:
    namespace: dict[str, object] = {}

    for symbol_name in PUBLIC_API_SMOKE_IMPORTS:
        exec(f"from phospy.api import {symbol_name}", namespace)

    assert PUBLIC_API_SMOKE_IMPORTS.issubset(namespace)

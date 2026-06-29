from __future__ import annotations

import pytest

import phospy
import phospy.api as public_api
import phospy.validation.datasets as dataset_validation

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


def test_dataset_validation_package_does_not_advertise_public_exports() -> None:
    wildcard_namespace: dict[str, object] = {}

    exec("from phospy.validation.datasets import *", wildcard_namespace)

    assert dataset_validation.__all__ == []
    assert DATASET_VALIDATOR_EXPORT_NAMES.isdisjoint(wildcard_namespace)
    assert all(
        not name.endswith("Validator")
        for name in wildcard_namespace
        if not name.startswith("__")
    )


def test_internal_dataset_validation_submodule_imports_still_work() -> None:
    from phospy.validation.datasets import site_metadata
    from phospy.validation.datasets.analysis_ready import (
        AnalysisReadyDatasetModelBoundaryValidator,
    )

    assert site_metadata.__name__ == "phospy.validation.datasets.site_metadata"
    assert AnalysisReadyDatasetModelBoundaryValidator.__name__ == (
        "AnalysisReadyDatasetModelBoundaryValidator"
    )


def test_public_api_import_smoke_still_passes() -> None:
    namespace: dict[str, object] = {}

    for symbol_name in PUBLIC_API_SMOKE_IMPORTS:
        exec(f"from phospy.api import {symbol_name}", namespace)

    assert PUBLIC_API_SMOKE_IMPORTS.issubset(namespace)

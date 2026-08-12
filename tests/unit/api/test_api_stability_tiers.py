from __future__ import annotations

import importlib.util
from pathlib import Path
from types import ModuleType

import pytest

import phospy.advanced as advanced_api
import phospy.api as public_api
from phospy._api_inventory import (
    ADVANCED_API_STABILITY_JUSTIFICATIONS,
    ADVANCED_PUBLIC_API_BASELINE_COUNT,
    STABLE_PUBLIC_API_BASELINE_COUNT,
)

ROOT = Path(__file__).resolve().parents[3]
API_GUIDE = ROOT / "docs" / "api" / "guide.md"
POLICY_ADR = ROOT / "docs" / "adr" / "adr_0031_public_api_stability_tiers.md"

DATASET_INTERNAL_DIAGNOSTIC_NAMES = frozenset(
    {
        "DatasetProcessingState",
        "MissingDataState",
        "NormalisationState",
        "SiteMatrixState",
        "SiteSequenceResolutionState",
        "TotalProteinCorrectionState",
        "BatchCorrectionDiagnostics",
        "BatchCorrectionReport",
        "RUVReadinessState",
        "RuvReadinessState",
        "IntensityScaleState",
    }
)


def test_api_tiers_are_explicit_disjoint_and_drive_all() -> None:
    stable_inventory = public_api._STABLE_PUBLIC_API
    advanced_inventory = public_api._ADVANCED_SUPPORTED_API
    internal_inventory = public_api._INTERNAL_EXPERIMENTAL_API
    stable = set(stable_inventory)
    advanced = set(advanced_inventory)
    internal = set(internal_inventory)

    assert stable
    assert advanced
    assert internal
    assert len(stable_inventory) == len(stable)
    assert len(advanced_inventory) == len(advanced)
    assert len(internal_inventory) == len(internal)
    assert stable.isdisjoint(advanced)
    assert stable.isdisjoint(internal)
    assert advanced.isdisjoint(internal)
    assert set(public_api.__all__) == stable
    assert set(advanced_api.__all__) == advanced


def test_api_surface_counts_do_not_increase_without_contract_review() -> None:
    assert len(public_api.__all__) == STABLE_PUBLIC_API_BASELINE_COUNT
    assert len(advanced_api.__all__) == ADVANCED_PUBLIC_API_BASELINE_COUNT


def test_advanced_api_names_have_stability_justifications() -> None:
    justifications = ADVANCED_API_STABILITY_JUSTIFICATIONS

    assert set(justifications) == set(advanced_api.__all__)
    assert all(justification.strip() for justification in justifications.values())


def test_inventory_names_are_not_orphaned_from_public_facades() -> None:
    for name in public_api.__all__:
        assert getattr(public_api, name) is not None

    for name in advanced_api.__all__:
        assert getattr(advanced_api, name) is not None


def test_advanced_exports_are_grouped_and_documented() -> None:
    advanced = set(public_api._ADVANCED_SUPPORTED_API)

    assert {
        "DatasetBatchCorrectionConfig",
        "SpsRuvBatchCorrectionConfig",
        "KinaseLibraryResourceLoader",
        "filter_differential_results",
        "rank_differential_results",
    } <= advanced
    assert {
        "DifferentialModelDiagnostics",
        "KinaseEligibilityReport",
        "KinaseWorkflowAttritionProvenance",
    } <= advanced
    assert {
        "publish_dataset",
        "publish_kinase_workflow",
        "publish_signalome_workflow",
    } <= advanced

    guide = API_GUIDE.read_text(encoding="utf-8")
    adr = POLICY_ADR.read_text(encoding="utf-8")
    for required in (
        "Stable public API",
        "Advanced supported API",
        "Internal / experimental API",
    ):
        assert required in public_api.__doc__
        assert required in guide
        assert required in adr
    assert "phospy.advanced" in guide
    assert "phospy.advanced" in adr


def test_advanced_exports_are_not_stable_api_exports() -> None:
    advanced = set(public_api._ADVANCED_SUPPORTED_API)

    assert advanced.isdisjoint(public_api.__all__)
    assert not any(name.endswith("Validator") for name in advanced_api.__all__)


def test_internal_inventory_documents_removed_previous_exports() -> None:
    internal = set(public_api._INTERNAL_EXPERIMENTAL_API)

    assert {
        "DatasetProcessingState",
        "ReferenceBundleValidationReport",
        "BatchCorrectionReport",
        "RuvReadinessState",
        "IMPORTER_QUALITY_STATUS_REPORTED",
        "DifferentialPolicyProvenance",
    } <= internal
    assert internal.isdisjoint(public_api.__all__)


def test_api_datasets_does_not_export_internal_processing_state_models() -> None:
    import phospy.api.datasets as dataset_api

    assert DATASET_INTERNAL_DIAGNOSTIC_NAMES.isdisjoint(set(dataset_api.__all__))


def test_internal_processing_state_models_are_not_importable_from_api_datasets() -> (
    None
):
    import phospy.api.datasets as dataset_api

    for symbol_name in DATASET_INTERNAL_DIAGNOSTIC_NAMES:
        assert not hasattr(dataset_api, symbol_name)
        with pytest.raises(ImportError):
            exec(f"from phospy.api.datasets import {symbol_name}", {})


def test_analysis_ready_dataset_is_public_from_api_datasets() -> None:
    import phospy.api.datasets as dataset_api
    from phospy.science.datasets.models import AnalysisReadyPhosphoDataset

    namespace: dict[str, object] = {}

    assert set(dataset_api.__all__) == {"AnalysisReadyPhosphoDataset"}
    assert hasattr(dataset_api, "AnalysisReadyPhosphoDataset")
    assert dataset_api.AnalysisReadyPhosphoDataset is AnalysisReadyPhosphoDataset
    exec("from phospy.api.datasets import AnalysisReadyPhosphoDataset", namespace)
    assert namespace["AnalysisReadyPhosphoDataset"] is AnalysisReadyPhosphoDataset


def test_api_dataset_submodule_does_not_export_internal_experimental_names() -> None:
    import phospy.api.datasets as dataset_api

    assert set(dataset_api.__all__).isdisjoint(public_api._INTERNAL_EXPERIMENTAL_API)


def test_dataset_diagnostic_public_reexport_modules_do_not_exist() -> None:
    assert importlib.util.find_spec("phospy.api.diagnostics") is None
    assert importlib.util.find_spec("phospy.api.advanced_datasets") is None


def test_api_submodules_respect_stability_tier_inventory() -> None:
    import phospy.api.configs as configs_api
    import phospy.api.datasets as dataset_api
    import phospy.api.requests as requests_api
    import phospy.api.results as results_api

    internal = set(public_api._INTERNAL_EXPERIMENTAL_API)
    advanced = set(public_api._ADVANCED_SUPPORTED_API)
    submodules: tuple[ModuleType, ...] = (
        dataset_api,
        results_api,
        configs_api,
        requests_api,
    )

    for submodule in submodules:
        leaked_names = set(submodule.__all__) & internal
        assert not leaked_names, (
            f"{submodule.__name__} exports internal/experimental names: "
            f"{sorted(leaked_names)}"
        )
        leaked_advanced_names = set(submodule.__all__) & advanced
        assert not leaked_advanced_names, (
            f"{submodule.__name__} exports advanced names from the stable API route"
        )

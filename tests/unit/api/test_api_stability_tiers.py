from __future__ import annotations

from pathlib import Path

import phospy.api as public_api

ROOT = Path(__file__).resolve().parents[3]
API_GUIDE = ROOT / "docs" / "api" / "guide.md"
POLICY_ADR = ROOT / "docs" / "adr" / "adr_0031_public_api_stability_tiers.md"


def test_api_tiers_are_explicit_disjoint_and_drive_all() -> None:
    stable = set(public_api._STABLE_PUBLIC_API)
    advanced = set(public_api._ADVANCED_SUPPORTED_API)
    internal = set(public_api._INTERNAL_EXPERIMENTAL_API)

    assert stable
    assert advanced
    assert internal
    assert stable.isdisjoint(advanced)
    assert stable.isdisjoint(internal)
    assert advanced.isdisjoint(internal)
    assert set(public_api.__all__) == stable | advanced


def test_advanced_exports_are_grouped_and_documented() -> None:
    advanced = set(public_api._ADVANCED_SUPPORTED_API)

    assert {
        "DatasetBatchCorrectionConfig",
        "SpsRuvBatchCorrectionConfig",
        "KinaseLibraryResourceLoader",
        "filter_differential_results",
        "rank_differential_results",
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


def test_internal_inventory_documents_removed_previous_exports() -> None:
    internal = set(public_api._INTERNAL_EXPERIMENTAL_API)

    assert {
        "DatasetProcessingState",
        "ReferenceBundleValidationReport",
        "BatchCorrectionReport",
        "IMPORTER_QUALITY_STATUS_REPORTED",
        "DifferentialPolicyProvenance",
    } <= internal
    assert internal.isdisjoint(public_api.__all__)

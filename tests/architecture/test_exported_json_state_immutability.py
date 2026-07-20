from __future__ import annotations

import dataclasses
import inspect
from collections.abc import Mapping
from typing import cast

import pandas as pd

import phospy.api as public_api
import phospy.api.results as public_results
import phospy.contracts.results as contract_results
from phospy.api.configs import ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL
from phospy.contracts.configs import EnrichmentConfig
from phospy.contracts.results import (
    EnrichmentWorkflowResult,
    ImporterMissingIntensitySummary,
    ImporterQualityReport,
    PhosphositeImportResult,
    ResultCaveat,
)
from phospy.provenance.hashing import hash_json_payload
from phospy.provenance.immutability import FrozenJsonMapping
from phospy.provenance.models import JsonValue
from phospy.science.enrichment.models import GeneSetCollection

_IMMUTABLE_JSON_FIELDS = {
    "phospy.contracts.results.base.ImporterMissingIntensitySummary."
    "missing_values_by_sample_id": (
        "typed importer count map; recursively frozen after count validation"
    ),
    "phospy.contracts.results.base.ImporterMissingIntensitySummary."
    "missing_values_by_source_column": (
        "typed importer count map; recursively frozen after count validation"
    ),
    "phospy.contracts.results.base.ImporterQualityReport.format_specific": (
        "importer-owned JSON metadata; recursively frozen at construction"
    ),
    "phospy.contracts.results.base.PhosphositeImportResult.diagnostics": (
        "importer-owned JSON diagnostics; recursively frozen at construction"
    ),
    "phospy.contracts.results.enrichment.EnrichmentWorkflowResult.diagnostics": (
        "contract-owned enrichment diagnostics; recursively frozen"
    ),
    "phospy.contracts.results.enrichment.EnrichmentWorkflowResult."
    "method_metadata": "contract-owned enrichment method metadata; recursively frozen",
    "phospy.contracts.results.enrichment.EnrichmentWorkflowResult."
    "background_summary": "contract-owned enrichment background summary; recursively frozen",
    "phospy.contracts.results.enrichment.EnrichmentWorkflowResult."
    "set_collection_summary": (
        "contract-owned enrichment set-collection summary; recursively frozen"
    ),
    "phospy.science.result_caveats.ResultCaveat.details": (
        "common result caveat JSON details; recursively frozen"
    ),
}

_REVIEWED_JSON_FIELD_ALLOWLIST = {
    "phospy.contracts.results.base.PhosphositeImportResult._sample_column_mapping": (
        "private scalar str->str map; public sample_column_mapping returns a copy"
    ),
    "phospy.contracts.results.kinase.KinaseWorkflowAttritionProvenance.metrics": (
        "PHOSPY-REV-004 owns kinase attrition JSON-state closure"
    ),
    "phospy.contracts.results.kinase.KinaseWorkflowAttritionProvenance.policy": (
        "PHOSPY-REV-004 owns kinase attrition JSON-state closure"
    ),
    "phospy.contracts.results.kinase.KinaseWorkflowAttritionProvenance."
    "policy_violations": "PHOSPY-REV-004 owns kinase attrition JSON-state closure",
    "phospy.science.datasets.preprocessing.protein_aware_preparation."
    "ProteinAwarePreparationReport.policy_parameters": (
        "PHOSPY-REV-003 domain inventory entry for preprocessing diagnostics"
    ),
    "phospy.science.differential.models.results.DifferentialAnalysisResult."
    "workflow_provenance": (
        "PHOSPY-REV-003 domain inventory entry for differential provenance"
    ),
    "phospy.science.differential.models.results.DifferentialAnalysisResult."
    "_contrast_tables": "private DataFrame mapping; DataFrame ownership is PHOSPY-REV-002",
    "phospy.science.prediction.models.KinaseScoringResult.score_scale_metadata": (
        "PHOSPY-REV-004 domain inventory entry for kinase score-scale metadata"
    ),
}


def _exported_result_classes() -> tuple[type[object], ...]:
    classes: set[type[object]] = set()
    for module in (public_results, contract_results):
        for name in getattr(module, "__all__", ()):
            value = getattr(module, name, None)
            if inspect.isclass(value):
                classes.add(value)
    return tuple(sorted(classes, key=lambda cls: f"{cls.__module__}.{cls.__name__}"))


def _json_like_field_ids() -> set[str]:
    field_ids: set[str] = set()
    for owner in _exported_result_classes():
        if not dataclasses.is_dataclass(owner):
            continue
        for field in dataclasses.fields(owner):
            annotation = str(field.type)
            if "Mapping[" not in annotation and "dict[" not in annotation:
                continue
            field_ids.add(f"{owner.__module__}.{owner.__qualname__}.{field.name}")
    return field_ids


def _minimal_import_result(
    diagnostics: Mapping[str, object] | None = None,
) -> PhosphositeImportResult:
    index = pd.Index(["MAPK1;S10;"], name="site_id")
    return PhosphositeImportResult(
        phospho_matrix_candidate=pd.DataFrame({"sample_a": [1.0]}, index=index),
        site_metadata_candidate=pd.DataFrame({"site": ["S10"]}, index=index),
        sample_column_mapping={"Intensity A": "sample_a"},
        diagnostics=diagnostics,
    )


def _minimal_enrichment_result() -> EnrichmentWorkflowResult:
    return EnrichmentWorkflowResult(
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=GeneSetCollection(
            sets={"mapk_pathway": ("AKT1", "MAPK1")},
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        ),
        config=EnrichmentConfig(),
        diagnostics={"nested": {"items": [1]}},
        method_metadata={"nested": {"items": [1]}},
        background_summary={"nested": {"items": [1]}},
        set_collection_summary={"nested": {"items": [1]}},
    )


def test_exported_frozen_json_fields_are_registered_and_protected() -> None:
    expected = set(_IMMUTABLE_JSON_FIELDS) | set(_REVIEWED_JSON_FIELD_ALLOWLIST)

    assert _json_like_field_ids() == expected

    quality_report = ImporterQualityReport(format_specific={"nested": {"items": [1]}})
    missing_summary = ImporterMissingIntensitySummary(
        missing_values_by_sample_id={"sample_a": 1},
        missing_values_by_source_column={"Intensity A": 1},
    )
    import_result = _minimal_import_result(diagnostics={"nested": {"items": [1]}})
    enrichment_result = _minimal_enrichment_result()
    caveat = ResultCaveat(
        code="json_state",
        severity="info",
        message="JSON state is frozen.",
        details={"nested": {"items": [1]}},
    )

    protected_values = {
        "phospy.contracts.results.base.ImporterMissingIntensitySummary."
        "missing_values_by_sample_id": missing_summary.missing_values_by_sample_id,
        "phospy.contracts.results.base.ImporterMissingIntensitySummary."
        "missing_values_by_source_column": (
            missing_summary.missing_values_by_source_column
        ),
        "phospy.contracts.results.base.ImporterQualityReport.format_specific": (
            quality_report.format_specific
        ),
        "phospy.contracts.results.base.PhosphositeImportResult.diagnostics": (
            import_result.diagnostics
        ),
        "phospy.contracts.results.enrichment.EnrichmentWorkflowResult."
        "diagnostics": enrichment_result.diagnostics,
        "phospy.contracts.results.enrichment.EnrichmentWorkflowResult."
        "method_metadata": enrichment_result.method_metadata,
        "phospy.contracts.results.enrichment.EnrichmentWorkflowResult."
        "background_summary": enrichment_result.background_summary,
        "phospy.contracts.results.enrichment.EnrichmentWorkflowResult."
        "set_collection_summary": enrichment_result.set_collection_summary,
        "phospy.science.result_caveats.ResultCaveat.details": caveat.details,
    }

    assert set(protected_values) == set(_IMMUTABLE_JSON_FIELDS)
    for field_id, value in protected_values.items():
        assert isinstance(value, FrozenJsonMapping), field_id
        assert not isinstance(value, dict), field_id

    nested = caveat.details["nested"]
    assert isinstance(nested, FrozenJsonMapping)
    assert nested["items"] == (1,)
    thawed = caveat.details.copy()
    thawed["nested"]["items"].append(2)  # type: ignore[union-attr]
    assert caveat.details.copy()["nested"]["items"] == [1]

    first_payload = caveat.to_payload()
    second_payload = caveat.to_payload()
    assert first_payload == second_payload
    assert hash_json_payload(cast(JsonValue, first_payload)) == hash_json_payload(
        cast(JsonValue, second_payload)
    )
    assert hash_json_payload(
        cast(JsonValue, enrichment_result.diagnostics.copy())
    ) == hash_json_payload(cast(JsonValue, enrichment_result.diagnostics.copy()))


def test_no_result_api_exports_new_immutable_container() -> None:
    for module in (public_api, public_results, contract_results):
        assert "FrozenJsonMapping" not in getattr(module, "__all__", ())

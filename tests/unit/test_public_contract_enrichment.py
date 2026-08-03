from __future__ import annotations

from collections.abc import Iterator, Mapping, Sequence
from dataclasses import fields
from typing import get_args, get_origin, get_type_hints

import numpy as np
import pandas as pd
import pytest

import phospy
import phospy.api as public_api
import phospy.api.requests as request_models
import phospy.api.results as result_models
import phospy.api.workflows as workflow_models
import phospy.workflows as native_workflows
from phospy.api import (
    ContractValidationError,
    EnrichmentConfig,
    EnrichmentIdentifierKind,
    EnrichmentOutsideBackgroundPolicy,
    EnrichmentResultRecord,
    EnrichmentSetCollection,
    EnrichmentWorkflow,
    EnrichmentWorkflowRequest,
    EnrichmentWorkflowResult,
    GeneSetCollection,
    MultipleTestingCorrection,
    PtmSetCollection,
    WorkflowValidationError,
)
from phospy.api.configs import (
    ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
    ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
    ENRICHMENT_METHOD_OVER_REPRESENTATION,
    ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP,
    ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_ERROR,
    MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG,
    MULTIPLE_TESTING_CORRECTION_NONE,
)


class _DuplicateKeyMapping(Mapping[object, object]):
    def __iter__(self) -> Iterator[object]:
        return iter(("duplicate", "duplicate"))

    def __len__(self) -> int:
        return 2

    def __getitem__(self, key: object) -> object:
        if key == "duplicate":
            return "value"
        raise KeyError(key)


_ENRICHMENT_JSON_MAPPING_FIELDS = (
    "diagnostics",
    "method_metadata",
    "background_summary",
    "set_collection_summary",
)


def _gene_collection() -> GeneSetCollection:
    return GeneSetCollection(
        sets={"mapk_pathway": ("AKT1", "MAPK1")},
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        term_names={"mapk_pathway": "MAPK pathway"},
        source_name="unit_test",
    )


def _ptm_collection() -> PtmSetCollection:
    return PtmSetCollection(
        sets={"motif_sites": ("rat|P12345|S10", "rat|P12345|T20")},
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        source_name="unit_test",
    )


def test_enrichment_request_constructs_from_selected_identifiers() -> None:
    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1", " MAPK1 "),
        background_universe=("AKT1", "MAPK1", "MTOR"),
        config=EnrichmentConfig(),
    )

    assert request.selected_identifiers == ("AKT1", " MAPK1 ")
    assert request.input_table is None
    assert request.background_universe == ("AKT1", "MAPK1", "MTOR")
    assert request.config.method == ENRICHMENT_METHOD_OVER_REPRESENTATION


def test_enrichment_request_constructs_from_input_table() -> None:
    table = pd.DataFrame({"site_key": ["rat|P12345|S10"]})

    request = EnrichmentWorkflowRequest(
        identifier_column="site_key",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        set_collection=_ptm_collection(),
        input_table=table,
        background_universe=("rat|P12345|S10", "rat|P12345|T20"),
    )

    assert request.input_table is table
    assert request.selected_identifiers is None
    assert request.set_collection.collection_kind == "ptm_set"


def test_enrichment_config_defaults_are_explicit() -> None:
    config = EnrichmentConfig()

    assert config.method == ENRICHMENT_METHOD_OVER_REPRESENTATION
    assert config.multiple_testing_correction == (
        MULTIPLE_TESTING_CORRECTION_BENJAMINI_HOCHBERG
    )
    assert config.min_set_size is None
    assert config.max_set_size is None
    assert (
        config.selected_outside_background_policy
        == ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_ERROR
    )
    assert (
        config.set_member_outside_background_policy
        == ENRICHMENT_OUTSIDE_BACKGROUND_POLICY_DROP
    )
    assert config.minimum_retained_foreground_fraction is None

    uncorrected = EnrichmentConfig(
        multiple_testing_correction=MULTIPLE_TESTING_CORRECTION_NONE
    )
    assert uncorrected.multiple_testing_correction == "none"


def test_enrichment_config_publishing_preset_has_configurable_set_size_bounds() -> None:
    default_publishing = EnrichmentConfig.publishing()
    custom_publishing = EnrichmentConfig.publishing(min_set_size=3, max_set_size=250)

    assert default_publishing.min_set_size == 5
    assert default_publishing.max_set_size == 500
    assert default_publishing.method == ENRICHMENT_METHOD_OVER_REPRESENTATION
    assert custom_publishing.min_set_size == 3
    assert custom_publishing.max_set_size == 250


def test_enrichment_unsupported_method_rejected() -> None:
    with pytest.raises(ContractValidationError, match="enrichment.method"):
        EnrichmentConfig(method="competitive")  # type: ignore[arg-type]


def test_enrichment_request_construction_is_passive() -> None:
    table = pd.DataFrame({"gene_symbol": ["AKT1"]})
    config = object()

    request = EnrichmentWorkflowRequest(
        identifier_column=" gene_symbol ",
        identifier_kind="accession",  # type: ignore[arg-type]
        set_collection=_ptm_collection(),
        input_table=table,
        selected_identifiers=(),
        background_universe=(),
        config=config,  # type: ignore[arg-type]
    )

    assert request.identifier_column == " gene_symbol "
    assert request.identifier_kind == "accession"
    assert request.set_collection is not None
    assert request.input_table is table
    assert request.selected_identifiers == ()
    assert request.background_universe == ()
    assert request.config is config


def test_enrichment_background_universe_is_required_and_non_empty() -> None:
    with pytest.raises(TypeError):
        EnrichmentWorkflowRequest(
            identifier_column="gene_symbol",
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            set_collection=_gene_collection(),
            selected_identifiers=("AKT1",),
        )

    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        selected_identifiers=("AKT1",),
        background_universe=(),
    )

    with pytest.raises(WorkflowValidationError, match="background_universe"):
        EnrichmentWorkflow().run(request)


def test_enrichment_request_requires_exactly_one_identifier_source() -> None:
    table = pd.DataFrame({"gene_symbol": ["AKT1"]})

    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        input_table=table,
        selected_identifiers=("AKT1",),
        background_universe=("AKT1", "MAPK1"),
    )
    with pytest.raises(WorkflowValidationError, match="exactly one"):
        EnrichmentWorkflow().run(request)

    missing_source_request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        background_universe=("AKT1", "MAPK1"),
    )
    with pytest.raises(WorkflowValidationError, match="exactly one"):
        EnrichmentWorkflow().run(missing_source_request)


def test_enrichment_gene_and_ptm_semantics_do_not_mix() -> None:
    with pytest.raises(ContractValidationError, match="gene_set_collection"):
        GeneSetCollection(
            sets={"bad": ("rat|P12345|S10",)},
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_SITE_KEY,
        )

    request = EnrichmentWorkflowRequest(
        identifier_column="gene_symbol",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_ptm_collection(),
        selected_identifiers=("AKT1",),
        background_universe=("AKT1", "MAPK1"),
    )

    with pytest.raises(WorkflowValidationError, match="must match"):
        EnrichmentWorkflow().run(request)


def test_enrichment_result_contract_is_shape_only() -> None:
    record = EnrichmentResultRecord(
        term_id="mapk_pathway",
        collection_kind="gene_set",
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        input_overlap_count=1,
        background_overlap_count=2,
        set_size=2,
        overlap_identifiers=("AKT1",),
    )

    result = EnrichmentWorkflowResult(
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        config=EnrichmentConfig(),
        records=(record,),
        diagnostics={"status": "not_computed_by_contract"},
    )

    assert result.records == (record,)
    assert result.records[0].p_value is None
    assert result.diagnostics == {"status": "not_computed_by_contract"}


def test_enrichment_summary_mappings_are_recursively_immutable() -> None:
    source_payloads = {
        field_name: {
            "nested": {
                "items": [{"score": 1.0}],
                "labels": ["kept"],
            }
        }
        for field_name in _ENRICHMENT_JSON_MAPPING_FIELDS
    }
    result = EnrichmentWorkflowResult(
        identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
        set_collection=_gene_collection(),
        config=EnrichmentConfig(),
        diagnostics=source_payloads["diagnostics"],
        method_metadata=source_payloads["method_metadata"],
        background_summary=source_payloads["background_summary"],
        set_collection_summary=source_payloads["set_collection_summary"],
    )

    for source in source_payloads.values():
        source["nested"]["items"][0]["score"] = 9.0
        source["nested"]["items"].append({"score": 10.0})
        source["nested"]["labels"].append("source-only")

    for field_name in _ENRICHMENT_JSON_MAPPING_FIELDS:
        mapping = getattr(result, field_name)
        nested = mapping["nested"]
        assert isinstance(nested, Mapping)
        assert nested["items"][0]["score"] == 1.0
        assert nested["labels"] == ("kept",)

        with pytest.raises(TypeError):
            nested["items"][0]["score"] = 2.0  # type: ignore[index]

        payload = mapping.copy()
        payload["nested"]["items"][0]["score"] = 3.0  # type: ignore[index]
        payload["nested"]["labels"].append("payload-only")  # type: ignore[union-attr]

        assert mapping.copy()["nested"]["items"] == [{"score": 1.0}]
        assert mapping.copy()["nested"]["labels"] == ["kept"]


@pytest.mark.parametrize("field_name", _ENRICHMENT_JSON_MAPPING_FIELDS)
@pytest.mark.parametrize(
    ("invalid_value", "expected"),
    (
        ({1: "numeric-key"}, "keys must be strings"),
        (_DuplicateKeyMapping(), "duplicate JSON object key"),
        ({"bad": float("nan")}, "finite JSON number"),
        ({"bad": float("inf")}, "finite JSON number"),
        ({"bad": {"unsupported"}}, "got set"),
        ({"bad": np.array([1])}, "got ndarray"),
        ({"bad": bytearray(b"x")}, "got bytearray"),
        ({"bad": object()}, "got object"),
    ),
)
def test_enrichment_summary_mappings_reject_invalid_json_values(
    field_name: str,
    invalid_value: Mapping[object, object],
    expected: str,
) -> None:
    kwargs: dict[str, object] = {name: {} for name in _ENRICHMENT_JSON_MAPPING_FIELDS}
    kwargs[field_name] = invalid_value

    with pytest.raises(ContractValidationError) as exc_info:
        EnrichmentWorkflowResult(
            identifier_kind=ENRICHMENT_IDENTIFIER_KIND_GENE_SYMBOL,
            set_collection=_gene_collection(),
            config=EnrichmentConfig(),
            **kwargs,
        )

    message = str(exc_info.value)
    assert f"enrichment_result.{field_name}" in message
    assert expected in message


def test_enrichment_public_contract_remains_typed_and_narrow() -> None:
    request_hints = get_type_hints(EnrichmentWorkflowRequest)
    result_hints = get_type_hints(EnrichmentWorkflowResult)
    config_hints = get_type_hints(EnrichmentConfig)

    assert request_hints["identifier_kind"] == EnrichmentIdentifierKind
    assert request_hints["set_collection"] == EnrichmentSetCollection
    assert get_origin(request_hints["background_universe"]) is Sequence
    assert result_hints["identifier_kind"] == EnrichmentIdentifierKind
    assert result_hints["set_collection"] == EnrichmentSetCollection
    assert get_args(result_hints["records"]) == (EnrichmentResultRecord, Ellipsis)
    assert config_hints["method"] == public_api.EnrichmentMethod
    assert config_hints["multiple_testing_correction"] == MultipleTestingCorrection
    assert (
        config_hints["selected_outside_background_policy"]
        == EnrichmentOutsideBackgroundPolicy
    )
    assert (
        config_hints["set_member_outside_background_policy"]
        == EnrichmentOutsideBackgroundPolicy
    )

    request_field_names = {field.name for field in fields(EnrichmentWorkflowRequest)}
    assert {"identifier_kind", "set_collection", "background_universe"} <= (
        request_field_names
    )
    assert "background_policy" not in request_field_names
    assert "analysis_level" not in request_field_names


def test_enrichment_public_imports_include_native_workflow() -> None:
    assert "EnrichmentWorkflowRequest" in request_models.__all__
    assert "EnrichmentWorkflowResult" in result_models.__all__
    assert "EnrichmentConfig" in public_api.__all__
    assert "EnrichmentWorkflowRequest" in public_api.__all__
    assert "EnrichmentWorkflowResult" in public_api.__all__
    assert "EnrichmentWorkflow" in public_api.__all__
    assert "EnrichmentWorkflowRequest" not in phospy.__all__
    assert "EnrichmentWorkflowResult" not in phospy.__all__
    assert "EnrichmentWorkflow" not in phospy.__all__
    assert "EnrichmentWorkflow" in workflow_models.__all__
    assert "EnrichmentWorkflow" in native_workflows.__all__
    assert public_api.EnrichmentWorkflow is workflow_models.EnrichmentWorkflow
    assert public_api.EnrichmentWorkflow is native_workflows.EnrichmentWorkflow
